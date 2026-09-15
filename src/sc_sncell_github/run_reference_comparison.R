suppressPackageStartupMessages({
  library(Matrix)
  library(spacexr)
  library(SummarizedExperiment)
  library(SpatialExperiment)
})
root <- "/Users/georicl/Documents/sc_sncell_github"
base <- file.path(root, "data/processed/Wu2021_known_composition_v1")
sp_meta <- read.csv(file.path(base, "spatial_metadata.csv"), row.names=1, check.names=FALSE)
sp_counts <- as(readMM(file.path(base, "spatial_counts.mtx")), "CsparseMatrix")
genes <- read.delim(file.path(base, "genes.tsv"), header=FALSE)[[1]]
rownames(sp_counts) <- genes
colnames(sp_counts) <- rownames(sp_meta)
spatial <- SpatialExperiment(assays=list(counts=sp_counts), colData=S4Vectors::DataFrame(sp_meta), spatialCoords=as.matrix(sp_meta[,c("x","y")]))
for (variant in c("donor_balanced", "epithelial_split", "balanced_split")) {
  inp <- file.path(root, "data/processed/Wu2021_reference_comparison_v1", variant)
  out <- file.path(root, "results/Wu2021_reference_comparison_v1", variant)
  cat("START", variant, format(Sys.time()), "\n")
  meta <- read.csv(file.path(inp,"reference_metadata.csv"), row.names=1, check.names=FALSE)
  counts <- as(readMM(file.path(inp,"reference_counts.mtx")), "CsparseMatrix")
  stopifnot(identical(genes, read.delim(file.path(inp,"genes.tsv"),header=FALSE)[[1]]),
            all(meta$donor_id %in% c("CID4471", "CID4535")),
            all(table(meta$reference_label[meta$nUMI >= 100]) >= 25),
            all(colSums(counts) == meta$nUMI), all(colSums(sp_counts) == sp_meta$nUMI))
  rownames(counts) <- genes; colnames(counts) <- rownames(meta)
  meta$reference_label <- factor(meta$reference_label)
  reference <- SummarizedExperiment(assays=list(counts=counts),colData=S4Vectors::DataFrame(meta))
  set.seed(9)
  start <- Sys.time()
  prepared <- createRctd(spatial, reference, cell_type_col="reference_label", ref_UMI_min=100, ref_n_cells_min=25)
  saveRDS(prepared, file.path(inp,"rctd_preprocessed.rds"))
  writeLines(prepared$internal_vars$gene_list_reg, file.path(out,"rctd_selected_genes.txt"))
  result <- runRctd(prepared,rctd_mode="full",max_cores=2)
  saveRDS(result,file.path(inp,"rctd_result.rds"))
  weights <- as.matrix(t(assay(result,"weights")))
  stopifnot(nrow(weights)==488, setequal(rownames(weights),rownames(sp_meta)), all(is.finite(weights)), all(weights>=0), all(rowSums(weights)>0))
  write.csv(weights,file.path(out,"rctd_raw_weights.csv"))
  write.csv(weights/rowSums(weights),file.path(out,"rctd_relative_weights.csv"))
  saveRDS(list(seed=9,config=prepared$config,mode="full",max_cores=2,elapsed_seconds=as.numeric(difftime(Sys.time(),start,units="secs"))),file.path(inp,"rctd_run_record.rds"))
  capture.output(sessionInfo(),file=file.path(out,"rctd_sessionInfo.txt"))
  cat("DONE",variant,format(Sys.time()),"\n")
}
