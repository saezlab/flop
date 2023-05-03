#Packages
library(tidyverse)
library(org.Hs.eg.db)
source('https://raw.githubusercontent.com/martingarridorc/biokit/f5298742edaa3118e22c3372e5ac1387024c492c/R/translate_matrix.R')
source('https://raw.githubusercontent.com/martingarridorc/biokit/f5298742edaa3118e22c3372e5ac1387024c492c/R/validators.R')
#Directory settings
db <- org.Hs.eg.db
filename <- "./data/GSE186341_ASPC_dream_counts.csv.gz"
megacounts <- tibble()
megacontrast <- tibble()
megametadata <- tibble()
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
samplenames %>% distinct(drug) %>% pull()
sel_drugs <- c("KW2449", "MK2206", "NILOTINIB", "DMSO")

metadata <- samplenames %>% 
  filter(drug!='UNTREATED') %>%
  filter(drug %in% sel_drugs) %>% 
  group_by(drug) %>% 
  slice_sample(n = 4)

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


write_tsv(megacounts, file = 'test__countdata.tsv')
write_tsv(megacontrast, file = 'test__contrast.tsv')
write_tsv(megametadata, file = 'test__metadata.tsv')
