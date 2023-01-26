library(tidyverse)

###Main###
#get DE files from given dataset ID
args <- commandArgs(trailingOnly = FALSE)
dataset_id <- args[grep("--dataset",args)+1]
biocontext <- args[grep("--bio",args)+1]
param_id <- args[grep("--param",args)+1]
pipeline_files <- args[grep("--files",args)+1] %>% strsplit(., split=' ') %>% .[[1]]
out_table <- readRDS(pipeline_files[1]) %>%
  select(ID)

#read files and write output
for(filename in pipeline_files){
  proc_data <- readRDS(filename)
  selected_data <- proc_data %>% select(ID, !!param_id)
  iter_name <- sub('__de.rds', '', filename) %>% sub(pattern = paste0(dataset_id, "__"), replacement = "", .)
  out_table <- full_join(out_table, selected_data, by='ID') %>% 
    rename(!!iter_name := !!param_id)
}
output_filename <- paste(dataset_id, biocontext, param_id, 'decouplerinput.tsv', sep='__')
write.table(out_table, output_filename, sep='\t', quote=FALSE, row.names=FALSE)

