library(tidyverse)
files <- list.files(pattern = 'Raw_data\\.tsv')

for(filename in files){
  datamat <- read.table(file = filename, header = TRUE, sep = "\t") %>% select(-X) %>% as.data.frame()
  new_data <- aggregate(. ~ gene_symbol, datamat, sum, na.rm = TRUE) %>% 
    mutate_if(is.numeric, round) %>%
    write.table(., filename, sep='\t', quote=FALSE, row.names=FALSE)
}

