library(tidyverse)

path_file <- args[grep("--file=", args)] %>%
  sub("jaccard_analysis.R", "", .) %>%
  sub("--file=", "", .)
source(paste0(path_file, "jaccard_analysis_helper.R"))

args <- commandArgs(trailingOnly = FALSE)
datafile <- args[grep("--file", args) + 1][2]
dataset_id <- args[grep("--dataset", args) + 1]

merged_data <- read_tsv(datafile)
jaccard_results <- jaccard_analysis(merged_data)
write_tsv(jaccard_results, file = paste0(dataset_id, "__jaccard.tsv"))

