library(tidyverse)

args <- commandArgs(trailingOnly = FALSE)
dataset_id <- args[grep("--dataset",args)+1]
pipeline_files <- args[grep("--files",args)+1] %>% strsplit(., split=' ') %>% .[[1]] %>% as_tibble()
run_info <- pipeline_files %>% separate(sep='__', col=value, into=c('datasetID', 'biocontext', 'statparam', 'dcmethod', 'resource', 'ext')) %>% select(-ext)
bio_contexts <- run_info %>% select(biocontext) %>% distinct() %>% pull()
resources <- run_info %>% select(resource) %>% distinct() %>% pull()
statparams <- run_info %>% select(statparam) %>% distinct() %>% pull()
merged_data <- tibble()

for(statparam in statparams){
  for(bio_context in bio_contexts){
    for(resource in resources){
      print(paste(bio_context, '__', statparam,'__cons__', resource, '__decoupleroutput.tsv', sep=''))
      file_list <- pipeline_files %>% pull() %>% print() %>% str_subset(.,pattern=paste(bio_context,'__', statparam,'__cons__', resource, '__decoupleroutput.tsv', sep='')) %>% print()
      for(file in file_list) {
        file_data <- read_tsv(file, show_col_types=FALSE) %>% 
            separate(...1, into=c('bio_context','norm', 'diffexp', 'dcmethod'), sep='__', remove=FALSE) %>% 
            dplyr::rename(obs_id = `...1`) %>%
            select(-dcmethod) %>% 
            pivot_longer(cols=c(5:ncol(.)), names_to = 'items', values_to = 'scores') %>% 
            mutate(statparam=!!statparam, resource=!!resource, pipeline=paste(norm, diffexp, sep = '+'))
        merged_data <- rbind(merged_data, file_data)
      }
    }
  }
}

filename <- paste(dataset_id, "result.tsv", sep='__')
write_tsv(merged_data, filename)