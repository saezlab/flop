#Packages
library(tidyverse)

ID_gen <- function(n = 10) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

n=5
random_ident <- ID_gen(n)
set.seed(1)

countfile <- list.files(path = "./unproc_data/GTex/", pattern = '*.gct', full.names = T)
metafiles <- list.files(path = "./unproc_data/GTex/", pattern = '*SampleAttributesDS.txt', full.names = T)
megacounts <- tibble()

data <- read_tsv(countfile, name_repair = make.names)  %>% dplyr::rename('gene_symbol'='Description') %>% dplyr::select(-Name) 
data_samples <- colnames(data)

metadata <- read_tsv(metafiles) %>% 
  dplyr::select(SAMPID, SMTSD) %>% 
  mutate(SAMPID = make.names(SAMPID)) %>% 
  mutate(SMTSD = gsub(' - ', '_', SMTSD)) %>%
  mutate(SMTSD = gsub(' ', '_', SMTSD)) %>%
  mutate(SMTSD = gsub('-', '_', SMTSD)) %>%
  mutate(SMTSD = gsub('\\(', '', SMTSD)) %>%
  mutate(SMTSD = tolower(gsub('\\)', '', SMTSD))) %>%
  dplyr::rename('sample_ID'='SAMPID', group='SMTSD') %>% 
  filter(sample_ID %in% data_samples)

thresh_meta <- metadata %>%
  group_by(group) %>%
  filter(n()>=150)

tissues <- thresh_meta %>% 
  dplyr::select(group) %>%
  distinct() %>%
  pull()

n=5

for (i in random_ident){
  dir.create(paste("./data/GTex_", i, "/", sep=''))
  megacounts <- tibble()
  megametadata <- tibble()
  
  sampled_meta <- thresh_meta %>% 
    group_by(group) %>%
    print() %>% 
    slice_sample(n=150)
  
  sampled_data <- data %>% 
    dplyr::select(gene_symbol, sampled_meta$sample_ID)
  
  write_tsv(sampled_meta, paste("./data/GTex_", i, "/metadata.tsv", sep=''))
  write_tsv(sampled_data, paste("./data/GTex_", i, "/countdata.tsv", sep=''))
  
}
  

