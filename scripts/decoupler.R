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
#pipeline_files <- list.files(pattern = paste(dataset_id, '__.*__de\\.', sep=''))

network <- get_progeny() %>% dplyr::rename(mor=weight)

res_decouple <- decouple(input_data, 
                         network, 
                         .source='source', 
                         .target='target',
                         minsize=0)
output_file <- paste(sub('__decouplerinput.tsv', '', input_file), 'decoupleroutput.tsv', sep='__')
write.table(res_decouple, output_file, sep='\t', quote=FALSE, row.names=FALSE)

# Transform to matrix
mat_consensus <- res_decouple %>%
  filter(statistic=='consensus') %>%
  pivot_wider_profile(id_cols = source, names_from = condition, 
                      values_from = score) %>%
  as.matrix()

pheatmap(mat_consensus, cluster_rows = F, cluster_cols = F, cellwidth = 20, cellheight = 10)



