suppressPackageStartupMessages(library(SummarizedExperiment))
root<-"/Users/georicl/Documents/sc_sncell_github"
samples<-commandArgs(trailingOnly=TRUE)
if(!length(samples)) samples<-"CID4535"
checks<-list()
for(sample in samples) {
  inp<-file.path(root,"data/processed/Wu2021_multislice_CAF_T_v1",sample)
  out<-file.path(root,"results/Wu2021_multislice_CAF_T_v1",sample)
  actual<-as.matrix(t(assay(readRDS(file.path(inp,"rctd_result.rds")),"weights")))
  raw<-as.matrix(read.csv(file.path(out,"rctd_raw_weights.csv"),row.names=1,check.names=FALSE))
  relative<-as.matrix(read.csv(file.path(out,"rctd_relative_weights.csv"),row.names=1,check.names=FALSE))
  checks[[paste0(sample,"_rds_csv")]]<-isTRUE(all.equal(actual,raw,tolerance=1e-12))
  checks[[paste0(sample,"_normalization")]]<-isTRUE(all.equal(raw/rowSums(raw),relative,tolerance=1e-12))
  record<-readRDS(file.path(inp,"rctd_run_record.rds"))
  checks[[paste0(sample,"_parameters")]]<-record$seed==9 && record$mode=="full" && record$max_cores==2 && record$config$gene_cutoff_reg==0.0002
}
name<-if(length(samples)==6) "rds_validation_all.txt" else paste0("rds_validation_",paste(samples,collapse="_"),".txt")
capture.output(checks,file=file.path(root,"results/Wu2021_multislice_CAF_T_v1",name))
print(checks);stopifnot(all(unlist(checks)))
