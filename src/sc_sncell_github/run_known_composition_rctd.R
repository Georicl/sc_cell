# 1. 读取固定供者划分后的计数与标签
suppressPackageStartupMessages({
  library(Matrix)
  library(spacexr)
  library(SummarizedExperiment)
  library(SpatialExperiment)
})
project_dir <- "/Users/georicl/Documents/sc_sncell_github"
input_dir <- file.path(project_dir, "data/processed/Wu2021_known_composition_v1")
output_dir <- file.path(project_dir, "results/Wu2021_known_composition_v1")
genes <- read.delim(file.path(input_dir, "genes.tsv"), header = FALSE, colClasses = "character")[[1]]
ref_meta <- read.csv(file.path(input_dir, "reference_metadata.csv"), row.names = 1, check.names = FALSE)
sp_meta <- read.csv(file.path(input_dir, "spatial_metadata.csv"), row.names = 1, check.names = FALSE)
ref_counts <- as(readMM(file.path(input_dir, "reference_counts.mtx")), "CsparseMatrix")
sp_counts <- as(readMM(file.path(input_dir, "spatial_counts.mtx")), "CsparseMatrix")
rownames(ref_counts) <- rownames(sp_counts) <- genes
colnames(ref_counts) <- rownames(ref_meta)
colnames(sp_counts) <- rownames(sp_meta)
ref_meta$reference_label <- factor(ref_meta$reference_label)
stopifnot(all(ref_meta$donor_id %in% c("CID4471", "CID4535")),
          all(table(ref_meta$reference_label[ref_meta$nUMI >= 100]) >= 25),
          all(colSums(ref_counts) == ref_meta$nUMI), all(colSums(sp_counts) == sp_meta$nUMI))
reference_se <- SummarizedExperiment(assays = list(counts = ref_counts), colData = S4Vectors::DataFrame(ref_meta))
spatial_spe <- SpatialExperiment(assays = list(counts = sp_counts),
                               colData = S4Vectors::DataFrame(sp_meta),
                               spatialCoords = as.matrix(sp_meta[, c("x", "y")]))

# 2. 使用参考标签构建表达特征，再对模拟spot进行full模式拟合
set.seed(9)
started <- Sys.time()
rctd_data <- createRctd(spatial_spe, reference_se, cell_type_col = "reference_label",
                       ref_UMI_min = 100, ref_n_cells_min = 25)
saveRDS(rctd_data, file.path(input_dir, "rctd_preprocessed.rds"))
writeLines(rctd_data$internal_vars$gene_list_reg, file.path(output_dir, "rctd_selected_genes.txt"))
result <- runRctd(rctd_data, rctd_mode = "full", max_cores = 2)
saveRDS(result, file.path(input_dir, "rctd_result.rds"))
weights <- as.matrix(t(assay(result, "weights")))
stopifnot(all(is.finite(weights)), all(weights >= 0), all(rowSums(weights) > 0))
write.csv(weights, file.path(output_dir, "rctd_raw_weights.csv"))
write.csv(weights / rowSums(weights), file.path(output_dir, "rctd_relative_weights.csv"))
capture.output(sessionInfo(), file = file.path(output_dir, "rctd_sessionInfo.txt"))
saveRDS(list(seed = 9, ref_n_cells_min = 25, ref_UMI_min = 100,
             mode = "full", max_cores = 2, config = rctd_data$config,
             elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))),
        file.path(input_dir, "rctd_run_record.rds"))
cat("完成RCTD：", nrow(weights), "个模拟spot；", ncol(weights), "种细胞类型\n")
