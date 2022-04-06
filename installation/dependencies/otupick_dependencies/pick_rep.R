#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))

## Usage: Rscript pick_rep.R 

df <- read.delim(args[1])

spl <- split(df, df$new_group)

best_bin <- function(x) {
  # Select what bins have the best quality score per cluster
  max_qual <- which(x$quality_score == max(x$quality_score))
  max_qual_df <- x[max_qual,]
  
  if (nrow(max_qual_df) > 1){
    # From those, select the ones with smallest number of contigs
    min_cont <- which(max_qual_df$contig_number == min(max_qual_df$contig_number))
    min_cont_df <- x[min_cont,]
    
    # From those, pick the greatest n50 value
    max_n50 <- which(min_cont_df$n50 == max(min_cont_df$n50))
    max_50_df <- x[max_n50,]
    
    # Pick the lowest strain_heterogeneity
    x[which.max(max_50_df$strain_heterogeneity),]
    
  } else {
    max_qual_df
  }
}

best_bins_list <- (lapply(spl, best_bin))

write.table((bind_rows(best_bins_list)), file='bestbins_metadata.tsv', quote=FALSE, sep='\t', row.names=F)

print("Best bins are picked")
