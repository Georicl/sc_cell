# 1. 读取原来的输入，设置参考重抽样的参数
# 这一步检查同一个供者内部的参考细胞变化，会不会改变空间预测。
# 每类有放回抽取原来的细胞数，作者标签、空间输入与模型参数保持一致。
# 抽样会出现重复细胞，不代表增加了新的独立细胞或供者。

suppressPackageStartupMessages({
  library(Matrix)
  library(spacexr)
  library(SummarizedExperiment)
  library(SpatialExperiment)
})

project_dir <- "/Users/georicl/Documents/sc_sncell_github"
input_dir <- file.path(project_dir, "data/processed/CID4535_combine_v1/rctd_input")
output_dir <- file.path(project_dir, "data/processed/CID4535_reference_bootstrap_v1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 默认完成10次；终端也可以传入部分种子，分批运行。
# 例如 Rscript src/sc_sncell_github/rctd_reference_bootstrap.R 1001 1002
seed_args <- commandArgs(trailingOnly = TRUE)
bootstrap_seeds <- if (length(seed_args)) as.integer(seed_args) else 1001:1010
fit_seed <- 9L
max_cores <- 2L

genes <- read.delim(file.path(input_dir, "genes.tsv"), header = FALSE,
                    colClasses = "character")[[1]]
reference_metadata <- read.csv(file.path(input_dir, "reference_metadata.csv"),
                               row.names = 1, check.names = FALSE)
spatial_metadata <- read.csv(file.path(input_dir, "spatial_metadata.csv"),
                             row.names = 1, check.names = FALSE)
reference_counts <- as(readMM(file.path(input_dir, "reference_counts.mtx")), "CsparseMatrix")
spatial_counts <- as(readMM(file.path(input_dir, "spatial_counts.mtx")), "CsparseMatrix")
stopifnot(nrow(reference_counts) == length(genes),
          nrow(spatial_counts) == length(genes),
          ncol(reference_counts) == nrow(reference_metadata),
          ncol(spatial_counts) == nrow(spatial_metadata))
rownames(reference_counts) <- rownames(spatial_counts) <- genes
colnames(reference_counts) <- rownames(reference_metadata)
colnames(spatial_counts) <- rownames(spatial_metadata)
reference_metadata$reference_label <- factor(reference_metadata$reference_label)
spatial_spe <- SpatialExperiment(
  assays = list(counts = spatial_counts),
  colData = S4Vectors::DataFrame(spatial_metadata),
  spatialCoords = as.matrix(spatial_metadata[, c("x", "y")])
)
type_indices <- split(seq_len(ncol(reference_counts)), reference_metadata$reference_label)
baseline <- readRDS(file.path(project_dir, "data/processed/CID4535_rctd_v1/rctd_preprocessed.rds"))

# 2. 在每类内部抽样，并记录实际抽到了哪些原始细胞
# 保持每类数量相同，避免同时改变参考细胞组成与抽样因素。
for (bootstrap_seed in bootstrap_seeds) {
  run_dir <- file.path(output_dir, sprintf("seed_%d", bootstrap_seed))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(file.path(run_dir, "completed.rds"))) {
    message("已有完整结果，跳过：", bootstrap_seed)
    next
  }
  log_connection <- file(file.path(run_dir, "execution.log"), open = "wt")
  sink(log_connection, split = TRUE)
  sink(log_connection, type = "message")
  started <- Sys.time()
  cat("开始种子", bootstrap_seed, format(started), "\n")
  set.seed(bootstrap_seed)
  sampled <- unlist(lapply(type_indices, function(indices) {
    indices[sample.int(length(indices), size = length(indices), replace = TRUE)]
  }), use.names = FALSE)

  sampled_metadata <- reference_metadata[sampled, , drop = FALSE]
  sampled_metadata$source_barcode <- rownames(reference_metadata)[sampled]
  rownames(sampled_metadata) <- sprintf("boot_%d_%05d", bootstrap_seed, seq_along(sampled))
  sampled_counts <- reference_counts[, sampled, drop = FALSE]
  colnames(sampled_counts) <- rownames(sampled_metadata)
  stopifnot(all(table(sampled_metadata$reference_label) == table(reference_metadata$reference_label)))
  write.csv(sampled_metadata, file.path(run_dir, "sampled_cells.csv"))
  unique_counts <- aggregate(source_barcode ~ reference_label, sampled_metadata,
                             function(x) length(unique(x)))
  names(unique_counts)[2] <- "n_unique_source_cells"
  write.csv(unique_counts, file.path(run_dir, "unique_cells_by_type.csv"), row.names = FALSE)
  reference_se <- SummarizedExperiment(
    assays = list(counts = sampled_counts),
    colData = S4Vectors::DataFrame(sampled_metadata)
  )

  # 3. 重建参考表达特征，使用与基线一致的预处理和full模式
  # 每轮输入共同基因相同；方法内部筛选基因可能因参考变化而改变。
  # 这属于整条参考到空间映射流程的敏感性，另存基因列表用于核对。
  set.seed(fit_seed)
  rctd_data <- createRctd(
    spatial_experiment = spatial_spe,
    reference_experiment = reference_se,
    cell_type_col = "reference_label",
    ref_UMI_min = 100,
    ref_n_cells_min = 15
  )
  saveRDS(rctd_data$config, file.path(run_dir, "preprocessing_config.rds"))
  saveRDS(rctd_data$internal_vars[c("gene_list_reg", "gene_list_bulk")],
          file.path(run_dir, "selected_genes.rds"))
  writeLines(rctd_data$internal_vars$gene_list_reg, file.path(run_dir, "genes_reg.txt"))
  cat("预处理后spot数", ncol(rctd_data$spatial_experiment), "\n")
  result <- runRctd(rctd_data, rctd_mode = "full", max_cores = max_cores)

  # 4. 保存原生结果和逐spot的原始权重，不覆盖原来的基线
  w <- as.matrix(t(assay(result, "weights")))
  stopifnot(all(is.finite(w)), all(w >= 0), all(rowSums(w) > 0))
  saveRDS(result, file.path(run_dir, "rctd_full_result.rds"))
  write.csv(w, file.path(run_dir, "weights.csv"))
  retained_genes <- rctd_data$internal_vars$gene_list_reg
  baseline_genes <- baseline$internal_vars$gene_list_reg
  run_summary <- data.frame(
    seed = bootstrap_seed, n_spots = nrow(w), n_types = ncol(w),
    n_genes_reg = length(retained_genes),
    gene_jaccard_vs_baseline = length(intersect(retained_genes, baseline_genes)) /
      length(union(retained_genes, baseline_genes)),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  )
  write.csv(run_summary, file.path(run_dir, "run_summary.csv"), row.names = FALSE)
  capture.output(sessionInfo(), file = file.path(run_dir, "sessionInfo.txt"))
  record <- list(bootstrap_seed = bootstrap_seed, fit_seed = fit_seed,
                 rctd_mode = "full", max_cores = max_cores,
                 ref_UMI_min = 100, ref_n_cells_min = 15,
                 source_md5 = tools::md5sum(file.path(input_dir, c(
                   "genes.tsv", "reference_counts.mtx", "spatial_counts.mtx",
                   "reference_metadata.csv", "spatial_metadata.csv"))),
                 finished = Sys.time())
  saveRDS(record, file.path(run_dir, "completed.rds"))
  print(run_summary)
  sink(type = "message")
  sink()
  close(log_connection)
  message("完成种子 ", bootstrap_seed)
}
