#Packages
library(tidyverse)

ID_gen <- function(n = 10) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

set.seed(1)
n=1
random_ident <- ID_gen(n)


countfile <- list.files(path = "./unproc_data/CCLE/", pattern = '*.gct', full.names = T)

data <- read_tsv(countfile, name_repair = make.names)  %>% dplyr::rename('gene_symbol'='Description') %>% dplyr::select(-Name) 
data_samples <- colnames(data)[-1]

metadata <- data_samples %>% 
  as_tibble() %>%
  separate(value, c('cell_line', 'group'), sep='_', extra='merge', remove=F) %>%
  dplyr::rename('sample_ID'='value')

thresh_meta <- metadata %>%
  group_by(group) %>%
  filter(n()>=20)

tissues <- thresh_meta %>% 
  dplyr::select(group) %>%
  distinct() %>%
  pull()


for (i in random_ident){
  dir.create(paste("./flop_data/CCLE_", i, "/", sep=''))
  
  sampled_meta <- thresh_meta %>% 
    group_by(group) %>%
    print() %>% 
    slice_sample(n=20)
  
  sampled_data <- data %>% 
    dplyr::select(gene_symbol, sampled_meta$sample_ID)

  count_name <- paste0("./flop_data/CCLE_", i, "/CCLE_", i, '__countdata.tsv')
  meta_name <- paste0("./flop_data/CCLE_", i, "/CCLE_", i, '__metadata.tsv')
  
  write_tsv(sampled_meta, meta_name)
  write_tsv(sampled_data, count_name)

}
  
