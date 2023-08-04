library(tidyverse)
library(qs)

# get DE files from given dataset ID, biological context and statistcal parameter
args <- commandArgs(trailingOnly = FALSE)
dataset_id <- args[grep("--dataset",args)+1]
biocontext <- args[grep("--bio",args)+1]
param_id <- args[grep("--param",args)+1]
threshold <- as.numeric(args[grep("--threshold",args)+1])
pipeline_files <- args[grep("--files",args)+1] %>% strsplit(., split=' ') %>% .[[1]]
out_table <- qread(pipeline_files[1]) %>%
  select(ID)

#read files and write output
failing_thresh <- c()
for(filename in pipeline_files){
  proc_data <- qread(filename)
  selected_data <- proc_data %>% select(ID, !!param_id)
  num_selected_genes <- proc_data %>% filter(padj < 0.05) %>% nrow()
  if(num_selected_genes < threshold){
    failing_thresh <- c(failing_thresh, filename)
  }
  iter_name <- sub('__de.qs', '', filename) %>% sub(pattern = paste0(dataset_id, "__"), replacement = "", .)
  out_table <- full_join(out_table, selected_data, by='ID') %>% 
    rename(!!iter_name := !!param_id)
}
output_filename <- paste(dataset_id, biocontext, param_id, 'decouplerinput.tsv', sep='__')

# If at least one pipeline passes the threshold, write the output of all pipelines
if(length(failing_thresh) < 6){
  write.table(out_table, output_filename, sep='\t', quote=FALSE, row.names=FALSE)
}


