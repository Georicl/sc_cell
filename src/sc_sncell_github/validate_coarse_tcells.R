suppressPackageStartupMessages(library(SummarizedExperiment))
root<-"/Users/georicl/Documents/sc_sncell_github"
checks<-list()
for(donor in c("CID4067","CID4463","CID3941","CID3948","CID4461")) {
  inp<-file.path(root,"data/processed/Wu2021_coarse_T_v1",donor)
  out<-file.path(root,"results/Wu2021_coarse_T_v1",donor)
  actual<-as.matrix(t(assay(readRDS(file.path(inp,"rctd_result.rds")),"weights")))
  raw<-as.matrix(read.csv(file.path(out,"rctd_raw_weights.csv"),row.names=1,check.names=FALSE))
  relative<-as.matrix(read.csv(file.path(out,"rctd_relative_weights.csv"),row.names=1,check.names=FALSE))
  checks[[paste0(donor,"_rds_csv")]]<-isTRUE(all.equal(actual,raw,tolerance=1e-12))
  checks[[paste0(donor,"_normalization")]]<-isTRUE(all.equal(raw/rowSums(raw),relative,tolerance=1e-12))
  record<-readRDS(file.path(inp,"rctd_run_record.rds"))
  checks[[paste0(donor,"_settings")]]<-record$seed==9 && record$mode=="full" && record$max_cores==2
}
capture.output(checks,file=file.path(root,"results/Wu2021_coarse_T_v1/rds_validation.txt"))
print(checks);stopifnot(all(unlist(checks)))
