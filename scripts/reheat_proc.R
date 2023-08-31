library(tidyverse)


data <- readRDS('./unproc_data/MetaHeart/MetaHeart__counts.rds')

# Preprocess the reheat datasets.
# please check the manuscript for more information about this. Where the age column contained NAs, 
# we filled those values with the average of the column. 

selected_datasets <- c('Spurell19', 'Liu15_R', 'Pepin19', 'Schiano17', 'Yang14')
dir.create('./flop_data/')

for(i in 1:length(data)){
  
  if(!names(data[i]) %in% selected_datasets){
    next
  }

  if(names(data[i]) == 'Liu15_R'){
    names(data)[i] <- 'Liu15'
  }

  dir.create(paste0('./flop_data/reheat_', names(data[i])))

  int_data <- data[[i]]
  count_data <- int_data$gex %>%
    as.data.frame() %>%
    rownames_to_column(var = 'gene_symbol') %>%
    mutate(across(!contains('gene_symbol'), as.integer))
  write_tsv(count_data, paste0('./flop_data/reheat_', names(data[i]), '/reheat_', names(data[i]),'__countdata.tsv'))

  metadata <- int_data$target %>%
      dplyr::rename(sample_ID = Sample, group = HeartFailure) %>%
      filter(!duplicated(.)) %>%
      relocate(sample_ID, group) %>%
      group_by(group) %>%
      {
        if (names(data[i]) == 'Yang14') {
          mutate(., Age = replace_na(Age, mean(Age, na.rm = TRUE))) %>%
          select(., -HTx)
        } else if (names(data[i]) == 'Spurell19') {
          select(., -DCM)
        } else if (names(data[i]) == 'Schiano17') {
          select(., -Disease, -DCM)
        } else if (names(data[i]) == 'Liu15') {
          mutate(., Age = replace_na(Age, mean(Age, na.rm = TRUE))) %>%
          select(., -Disease, -DCM)
        } else if (names(data[i]) == 'Pepin19') {
          mutate(., Age = replace_na(Age, mean(Age, na.rm = TRUE))) %>%
          select(., -HTx)
        } else {
          .
        }
      }

  


  write_tsv(metadata, paste0('./flop_data/reheat_', names(data[i]), '/reheat_', names(data[i]), '__metadata.tsv'))

  contrast <- data.frame(group1 = 'yes', group2 = 'no')
  write_tsv(contrast, paste0('./flop_data/reheat_', names(data[i]), '/reheat_', names(data[i]), '__contrast.tsv'))
}

