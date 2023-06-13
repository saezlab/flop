#Packages
library(tidyverse)
library(org.Hs.eg.db)
source('https://raw.githubusercontent.com/martingarridorc/biokit/f5298742edaa3118e22c3372e5ac1387024c492c/R/translate_matrix.R')
source('https://raw.githubusercontent.com/martingarridorc/biokit/f5298742edaa3118e22c3372e5ac1387024c492c/R/validators.R')
#Directory settings
db <- org.Hs.eg.db
filenames <- list.files(path = "../unproc_data/GSE186341/", pattern = '*.csv', , full.names = T)
megacounts <- tibble()
megacontrast <- tibble()
megametadata <- tibble()
for(filename in filenames){
  data <- read.csv(filename) %>% 
    column_to_rownames(.,'X') %>% 
    as.matrix()
  cellID <- strsplit(filename, split='/')[[1]] %>% tail(1) %>% strsplit(., split='_') %>% unlist() %>% .[2]
  
  translated_data <- translateMatrixWithDb(data, db, 'ENTREZID', 'SYMBOL', sum) %>% 
    as.data.frame() %>% 
    rownames_to_column(var='gene_symbol') %>% 
    as_tibble() %>% 
    rename_with( ~ paste(cellID, .x, sep='_')) %>%
    dplyr::rename(gene_symbol=paste(cellID, "gene_symbol", sep='_'))
  
  samplenames <- translated_data %>% 
    colnames() %>% 
    as_tibble() %>% 
    filter(value!='gene_symbol') %>%
    dplyr::rename(sample_ID=value) %>% 
    separate(.,col='sample_ID', sep=c('_'), into=c('cell', 'drug', 'dose', 'time'), remove=FALSE) %>%
    mutate(time = str_replace_all(time, "\\.*\\..*", ""),
    group = paste(cell, drug, sep="_"))

  if(cellID=='ASPC'){
    megacounts <- translated_data %>% 
      dplyr::select(gene_symbol)
  }

  metadata <- samplenames %>% 
    filter(drug!='UNTREATED')
  
  genecounts <- translated_data %>% 
    dplyr::select(gene_symbol, as.vector(metadata$sample_ID)) 
  
  megacounts <- megacounts %>% inner_join(genecounts, by='gene_symbol')

  treatments <- metadata %>% 
    filter(drug != "DMSO") %>%
    distinct(group) %>% 
    pull()

  contrast_mat <- tibble(group1 = treatments, group2 = paste(cellID, 'DMSO', sep = "_"))

  megacontrast <- rbind(megacontrast, contrast_mat)
  megametadata <- rbind(megametadata, metadata)
}

write_tsv(megacounts, file = '../data/GSE186341/GSE186341__countdata.tsv')
write_tsv(megacontrast, file = '../data/GSE186341/GSE186341__contrast.tsv')
write_tsv(megametadata, file = '../data/GSE186341/GSE186341__metadata.tsv')
