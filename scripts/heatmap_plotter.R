library(tidyverse)
library(pheatmap)
library(decoupleR)
library(ggplot2)

args <- commandArgs(trailingOnly = FALSE)
dataID <- args[grep("--ID",args)+1]
input_files <- list.files('.\\results', pattern=paste(dataID, '.*\\.tsv', sep=''))
for (input_file in input_files){
  res_decouple <- read.table(file = paste('.\\results\\',input_file, sep=''), header = TRUE, sep = "\t")
  annotated_tab <- res_decouple %>%
    separate(X, c('sampleID', 'normalization', 'differential_analysis', 'functional_analysis' ), sep='__', remove=FALSE) %>%
    select(-sampleID)

  numeric_info <- annotated_tab %>% select(-(2:4)) %>% column_to_rownames(var='X')
  annotations <- annotated_tab %>% select(1:4) %>% column_to_rownames(var='X')
  stat <- strsplit(input_file, split='__')[[1]][2]
  heatmap <- pheatmap(numeric_info, annotation_row = annotations, cellwidth=20, cellheight=10, main=stat)
  

  run_info <- sub('__decoupleroutput.tsv', '', input_file)
  ggsave(filename=paste('.\\results\\heatmaps\\', run_info, "__heatmap.png", sep=''), plot=heatmap, width=10, height=5)
}

