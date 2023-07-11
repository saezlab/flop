library(tidyverse)
library(stats)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
path_file <- args[grep("--file=", args)] %>%
  sub("rank_analysis.R", "", .) %>%
  sub("--file=", "", .)
source(paste0(path_file, "rank_analysis_helper.R"))

args <- commandArgs(trailingOnly = FALSE)
datafile <- args[grep("--file", args) + 1][2]
dataset_id <- args[grep("--dataset",args)+1]

merged_data <- read_tsv(datafile)
cor_results <- corr_analysis(merged_data, "pipeline")
write_tsv(cor_results, file = paste0(dataset_id, "__total__rank.tsv"))
# #Correlation split
# #By norm
# norm_cor_results <- corr_analysis(merged_data, "norm")
# write_tsv(norm_cor_results, file = paste0(dataset_id, "__norm__rank.tsv"))
# #By diffexp
# diffexp_cor_results <- corr_analysis(merged_data, "diffexp")
# write_tsv(diffexp_cor_results, file = paste0(dataset_id, "__diffexp__rank.tsv"))

