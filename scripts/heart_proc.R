library(tidyverse)

data <- readRDS('./unproc_data/MetaHeart/MetaHeart__counts.rds')



for(i in 1:length(data)){
  dir.create(paste0('./data/', names(data[i])))
  int_data <- data[[i]]
  count_data <- int_data$gex %>%
    as.data.frame() %>%
    rownames_to_column(var = 'gene_symbol') %>%
    mutate(across(!contains('gene_symbol'), as.integer))
  write_tsv(count_data, paste0('./data/', names(data[i]), '/', names(data[i]), '__countdata.tsv'))

  metadata <- int_data$target %>%
    dplyr::rename(sample_ID = Sample, group = HeartFailure) %>%
    .[!duplicated(as.list(.))] %>%
    relocate(sample_ID, group) %>% #select(-SampleType, -TechnicalTime, -batch_sign) %>% select(-DCM)
    group_by(group) %>% select(-HTx)
    #%>% select(-Disease)
    # mutate_at(vars(Age), ~replace_na(., mean(., na.rm = TRUE))) %>% select(-Disease)
  


  write_tsv(metadata, paste0('./data/', names(data[i]), '/', names(data[i]), '__metadata.tsv'))

  contrast <- data.frame(group1 = 'yes', group2 = 'no')
  write_tsv(contrast, paste0('./data/', names(data[i]), '/', names(data[i]), '__contrast.tsv'))
}

