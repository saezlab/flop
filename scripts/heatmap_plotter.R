library(tidyverse)
library(pheatmap)
library(decoupleR)
library(ggplot2)

args <- commandArgs(trailingOnly = FALSE)
input_file <- args[grep("--input",args)+1]
res_decouple <- read.table(file = input_file, header = TRUE, sep = "\t")

run_info <- sub('__decoupleroutput.tsv', '', input_file)


mat_consensus <- res_decouple %>%
  filter(statistic=='consensus') %>%
  pivot_wider_profile(id_cols = source, names_from = condition, 
                      values_from = score) %>%
  as.matrix()

heatmap <- pheatmap(mat_consensus, cluster_rows = F, cluster_cols = F, cellwidth = 30, cellheight = 15, main = run_info)
ggsave(filename=paste(run_info, "__heatmap.png"), plot=heatmap)
