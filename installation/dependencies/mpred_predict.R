library(base)
load("mpred_models.RData")

args = commandArgs(trailingOnly=TRUE)
if (! file.exists(args[1])) {
    write("Usage: mpred_predict.R <unique_kmer33_55.txt>", stderr())
    quit(status=1)
}


file_unique_kmers = args[1]
data <- read.csv(file_unique_kmers, header = F, sep = "\t")
colnames(data) <- c("dataset","33","55")

source("mpred_function_predict_memory.R")

#data <- read.delim(file = "unique_kmer33_55.txt", sep="\t", header=F)
#data <- read.csv("unique_kmer33_55_EXAMPLE.txt", header = F, sep = "\t")
results <- as.data.frame(t(apply(data, 1, predict_memory)))
write.table(results, file = "metaspades_prediction.tsv", sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
