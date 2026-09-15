suppressPackageStartupMessages(library(SummarizedExperiment))
root <- "/Users/georicl/Documents/sc_sncell_github"
checks <- list()
for(donor in c("CID3941","CID3948","CID4461")) for(variant in c("baseline","donor_balanced","epithelial_split","balanced_split")) {
  inp <- file.path(root,"data/processed/Wu2021_tcell_validation_v1",donor,variant)
  out <- file.path(root,"results/Wu2021_tcell_validation_v1",donor,variant)
  actual <- as.matrix(t(assay(readRDS(file.path(inp,"rctd_result.rds")),"weights")))
  raw <- as.matrix(read.csv(file.path(out,"rctd_raw_weights.csv"),row.names=1,check.names=FALSE))
  relative <- as.matrix(read.csv(file.path(out,"rctd_relative_weights.csv"),row.names=1,check.names=FALSE))
  tag <- paste(donor,variant,sep="_")
  checks[[paste0(tag,"_rds_csv")]] <- isTRUE(all.equal(actual,raw,tolerance=1e-12))
  checks[[paste0(tag,"_normalization")]] <- isTRUE(all.equal(raw/rowSums(raw),relative,tolerance=1e-12))
  record <- readRDS(file.path(inp,"rctd_run_record.rds"))
  checks[[paste0(tag,"_settings")]] <- record$seed==9 && record$mode=="full" && record$max_cores==2
}
capture.output(checks,file=file.path(root,"results/Wu2021_tcell_validation_v1/rds_validation.txt"))
print(checks)
stopifnot(all(unlist(checks)))
