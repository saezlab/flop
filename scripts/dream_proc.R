#Packages
library(tidyverse)
library(org.Hs.eg.db)
source('https://raw.githubusercontent.com/martingarridorc/biokit/f5298742edaa3118e22c3372e5ac1387024c492c/R/translate_matrix.R')
source('https://raw.githubusercontent.com/martingarridorc/biokit/f5298742edaa3118e22c3372e5ac1387024c492c/R/validators.R')
#Directory settings
db <- org.Hs.eg.db
filenames <- list.files(path = "./unproc_data/GSE186341/", pattern = '*.csv', , full.names = T)
megacounts <- tibble()
for(filename in filenames){
  data <- read.csv(filename) %>% 
    column_to_rownames(.,'X') %>% 
    as.matrix()
  cellID <- strsplit(filename, split='_')[[1]][3]
  
  translated_data <- translateMatrixWithDb(data, db, 'ENTREZID', 'SYMBOL', sum) %>% 
    as.data.frame() %>% 
    rownames_to_column(var='gene_symbol') %>% 
    as_tibble()
  
  samplenames <- translated_data %>% 
    colnames() %>% 
    as_tibble() %>% 
    filter(value!='gene_symbol') %>%
    dplyr::rename(sample_ID=value) %>% 
    separate(.,col='sample_ID', sep=c('_'), into=c('group', 'dose', 'time'), remove=FALSE) %>%
    mutate(time = str_replace_all(time, "\\.*\\..*", ""))
  
  drugs <- (samplenames %>% 
              filter(group!='DMSO') %>% 
              filter(group!='UNTREATED') %>% 
              distinct(group) %>% 
              as.vector())$group
  if(cellID=='ASPC'){
    megacounts <- translated_data %>% 
      dplyr::select(gene_symbol)
  }
  
  for(drug in drugs){
    metadata <- samplenames %>% 
      filter(group=='DMSO' | group==as.symbol(drug))
    check_treat <- (metadata %>% 
      dplyr::select(group) %>%
      filter(group==as.symbol(drug)) %>% 
      count() %>% 
      as.numeric())[1]
    
    if(check_treat<2){
      print('Not enough treatment replicates')
      break
    } 
    genecounts <- translated_data %>% 
      dplyr::select(as.vector(metadata$sample_ID)) %>% 
      rename_with( ~ paste(cellID, drug, .x, sep='__')) %>%
      mutate(gene_symbol=translated_data$gene_symbol) %>%
      relocate(gene_symbol)
    
    count_name <- paste("./data/GSE186341", cellID, drug, '_countdata.tsv', sep = '_')
    meta_name <- paste("./data/GSE186341", cellID, drug, '_metadata.tsv', sep='_')
    
    write.table(metadata, meta_name, sep='\t',quote=FALSE, row.names=FALSE)
    write.table(genecounts, count_name, sep='\t',quote=FALSE, row.names=FALSE)
  }
}
