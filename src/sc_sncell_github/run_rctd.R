library(Matrix)
library(spacexr)
library(SummarizedExperiment)
library(SpatialExperiment)

# 输入目录就是之前通过CID4535测试集中得出的数据集
input_dir <- paste0(
  "/Users/georicl/Documents/sc_sncell_github/",
  "data/processed/CID4535_combine_v1/rctd_input"
)


genes <- read.delim(file.path(input_dir, "genes.tsv"),
  header = FALSE, colClasses = "character"
)[[1]]

reference_metadata <- read.csv(
  file.path(input_dir, "reference_metadata.csv"),
  row.names = 1,
  check.names = FALSE
)

spatial_metadata <- read.csv(
  file.path(input_dir, "spatial_metadata.csv"),
  row.names = 1,
  check.names = FALSE
)

reference_counts <- as(
  readMM(file.path(input_dir, "reference_counts.mtx")),
  "CsparseMatrix"
)

spatial_counts <- as(
  readMM(file.path(input_dir, "spatial_counts.mtx")),
  "CsparseMatrix"
)


# 把genes和barcode的名字赋予矩阵
rownames(reference_counts) <- genes
colnames(reference_counts) <- rownames(reference_metadata)

rownames(spatial_counts) <- genes
colnames(spatial_counts) <- rownames(spatial_metadata)

#############
# 构建模型对象
#############
# 通过label标签转换为分类
reference_metadata$reference_label <- factor(
  reference_metadata$reference_label
)
# 创建单细胞参考数据
reference_se <- SummarizedExperiment(
  assays = list(counts = reference_counts),
  colData = S4Vectors::DataFrame(reference_metadata)
)
# 创建空间转录组数据
spatial_spe <- SpatialExperiment(
  assays = list(counts = spatial_counts),
  colData = S4Vectors::DataFrame(spatial_metadata),
  spatialCoords = as.matrix(spatial_metadata[, c("x", "y")])
)


print(reference_se)
print(spatial_spe)

print(sort(table(reference_se$reference_label)))


# 当前输入包含公共基因，并在这个范围里对UMI求和

reference_umi <- Matrix::colSums(assay(reference_se, "counts"))

# 本轮使用RCTD的参考细胞UMI：100
keep_reference <- reference_umi >= 100
# 获取UMI
ref_n <- table(reference_se$reference_label[keep_reference])

print(sort(ref_n))

###############################
# RCTD预处理： 构建表达特征、筛选基因和spot
###############################
set.seed(9)
# 创建RCTD对象
rctd_data <- createRctd(
  spatial_experiment = spatial_spe,
  reference_experiment = reference_se,
  cell_type_col = "reference_label",
  ref_UMI_min = 100,
  ref_n_cells_min = 15
)

# 找出被移除的空间spot
removed_spots <- setdiff(
  colnames(spatial_spe),
  colnames(rctd_data$spatial_experiment)
)

print(removed_spots)


output_dir <- paste0(
  "/Users/georicl/Documents/sc_sncell_github/",
  "data/processed/CID4535_rctd_v1"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 储存预处理后的rctd数据
saveRDS(
  rctd_data,
  file.path(output_dir, "rctd_preprocessed.rds")
)

rctd_result <- runRctd(
  rctd_data,
  rctd_mode = "full",
  max_cores = 2
)

saveRDS(
  rctd_result,
  file.path(output_dir, "rctd_full_result.rds")
)

weights <- assay(rctd_result, "weights")

print(dim(weights))
print(rownames(weights))

weights <- assay(rctd_result, "weights")

# 行是细胞类型，列是 spot
print(dim(weights))
print(rownames(weights))

# 查看前 5 个 spot 的 CAF、CD8 和癌上皮预测权重
print(
  weights[
    c("CAFs", "T cells CD8+", "Cancer Epithelial"),
    1:5
  ]
)

# 查看每个 spot 的权重之和
print(summary(colSums(weights)))

stopifnot(
  all(is.finite(weights)),
  all(weights >= 0)
)


# 转置后：每行一个 spot，每列一种细胞类型
write.csv(
  as.matrix(t(weights)),
  file.path(output_dir, "rctd_weights_spot_by_celltype.csv"),
  row.names = TRUE
)

saveRDS(
  list(
    seed = 9,
    rctd_mode = "full",
    max_cores = 2,
    spacexr_version = as.character(packageVersion("spacexr"))
  ),
  file.path(output_dir, "rctd_run_parameters.rds")
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)
