args <- commandArgs(trailingOnly = TRUE)




library(Rcpp,lib.loc="/data/msb/tools/virfinder/r_packages")
library(qvalue,lib.loc="/data/msb/tools/virfinder/r_packages")
library(glmnet,lib.loc="/data/msb/tools/virfinder/r_packages")
library(VirFinder,lib.loc="/data/msb/tools/virfinder/r_packages")


setwd(args[1])

infafile <- args[2]

predResult <- VF.pred(infafile)

predResult[order(predResult$pvalue),]

predResult$qvalue <- VF.qvalue(predResult$pvalue)

write.table(predResult, args[3], sep="\t")




