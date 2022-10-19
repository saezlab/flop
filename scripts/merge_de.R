library(tidyverse)

###Main###
#get DE files from given dataset ID
args <- commandArgs(trailingOnly = FALSE)
dataset_id <- args[grep("--dataset",args)+1]
param_id <- args[grep("--param",args)+1]
pipeline_files <- list.files(pattern = paste(dataset_id, '__.*__de\\.', sep=''))
out_table <- read.table(file = pipeline_files[1], header = TRUE, sep = "\t") %>%
  select(ID)

#read files and write output
for(filename in pipeline_files){
  proc_data <- read.table(file = filename, header = TRUE, sep = "\t")
  selected_data <- proc_data %>% select(ID, !!param_id)
  iter_name <- sub('__de.tsv', '', filename)
  out_table <- left_join(out_table, selected_data, by='ID') %>% 
    rename(!!iter_name := !!param_id)
}
output_filename <- paste(dataset_id, param_id, 'decouplerinput.tsv', sep='__')
write.table(out_table, output_filename, sep='\t', quote=FALSE, row.names=FALSE)

