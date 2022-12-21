library(tidyverse)

args <- commandArgs(trailingOnly = FALSE)
dataset_id <- args[grep("--dataset",args)+1]
pipeline_files <- args[grep("--files",args)+1] %>% strsplit(., split=' ') %>% .[[1]] %>% as_tibble()
bio_contexts <- pipeline_files %>% pull() %>% sub(pattern='__.*__.*__.*__decoupleroutput.tsv', replacement='') %>% unique()
resources <- pipeline_files %>% filter(str_detect(.,pattern='__decoupleroutput.tsv')) %>% pull() %>% sub(pattern='__decoupleroutput.tsv', replacement='', .)%>% sub(pattern='.*__.*__.*__', replacement='', .) %>% unique()
merged_data <- tibble()
statparams <- c('stat', 'logFC')

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