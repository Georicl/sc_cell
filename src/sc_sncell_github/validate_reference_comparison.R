suppressPackageStartupMessages(library(SummarizedExperiment))
root <- "/Users/georicl/Documents/sc_sncell_github"
checks <- list()
for (variant in c("donor_balanced","epithelial_split","balanced_split")) {
  inp <- file.path(root,"data/processed/Wu2021_reference_comparison_v1",variant)
  out <- file.path(root,"results/Wu2021_reference_comparison_v1",variant)
  result <- readRDS(file.path(inp,"rctd_result.rds"))
  raw <- as.matrix(read.csv(file.path(out,"rctd_raw_weights.csv"),row.names=1,check.names=FALSE))
  relative <- as.matrix(read.csv(file.path(out,"rctd_relative_weights.csv"),row.names=1,check.names=FALSE))
  actual <- as.matrix(t(assay(result,"weights")))
  checks[[paste0(variant,"_rds_matches_csv")]] <- isTRUE(all.equal(actual,raw,tolerance=1e-12))
  checks[[paste0(variant,"_normalization")]] <- isTRUE(all.equal(raw/rowSums(raw),relative,tolerance=1e-12))
  record <- readRDS(file.path(inp,"rctd_run_record.rds"))
  checks[[paste0(variant,"_run_settings")]] <- record$seed==9 && record$mode=="full" && record$max_cores==2
}
capture.output(checks,file=file.path(root,"results/Wu2021_reference_comparison_v1/rds_validation.txt"))
print(checks)
stopifnot(all(unlist(checks)))
