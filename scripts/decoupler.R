library(decoupleR)
library(tidyverse)
library(pheatmap)


args <- commandArgs(trailingOnly = FALSE)
input_file <- args[grep("--input",args)+1]
input_data <- read.table(file = input_file, header = TRUE, sep = "\t") %>%
  na.omit() %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "ID") %>% 
  as.matrix()

network <- get_progeny() %>% dplyr::rename(mor=weight)

res_decouple <- decouple(input_data, 
                         network, 
                         .source='source', 
                         .target='target',
                         minsize=0)
output_file <- paste(sub('__decouplerinput.tsv', '', input_file), 'decoupleroutput.tsv', sep='__')
write.table(res_decouple, output_file, sep='\t', quote=FALSE, row.names=FALSE)
