# 从已完成的模型输入中导出参考特征，不重新拟合。
root <- "/Users/georicl/Documents/sc_sncell_github"
out <- file.path(root,"results/Wu2021_tcell_label_review_v1")
dir.create(out,showWarnings=FALSE,recursive=TRUE)
for(v in c("baseline","donor_balanced","epithelial_split","balanced_split")) {
  x <- readRDS(file.path(root,"data/processed/Wu2021_tcell_validation_v1/CID3941",v,"rctd_preprocessed.rds"))
  write.csv(x$cell_type_info$info[[1]],file.path(out,paste0(v,"_rctd_profiles.csv")))
}
selection_function <- get("getDeGenes", envir=asNamespace("spacexr"))
writeLines(deparse(selection_function),file.path(out,"installed_gene_selection_function.txt"))
