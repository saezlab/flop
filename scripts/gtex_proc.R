#Packages
library(tidyverse)
library(rstudioapi)
library(org.Hs.eg.db)

set.seed(1)

#Directory settings
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

db <- org.Hs.eg.db
countfile <- list.files(path = "./data/GTex/", pattern = '*.gct', full.names = T)
metafiles <- list.files(path = "./data/GTex/", pattern = '*SampleAttributesDS.txt', full.names = T)
megacounts <- tibble()

data <- read_tsv(countfile) %>% dplyr::rename('gene_symbol'='Description') %>% dplyr::select(-Name)

data_samples <- colnames(data)

metadata <- read_tsv(metafiles) %>% 
  dplyr::select(SAMPID, SMTSD) %>% 
  mutate(SMTSD = gsub(' - ', '-', tolower(SMTSD))) %>% 
  mutate(SMTSD = gsub(' ', '-', SMTSD)) %>% 
  dplyr::rename('sample_ID'='SAMPID', group='SMTSD') %>% 
  filter(sample_ID %in% data_samples)

thresh_meta <- metadata %>%
  group_by(group) %>%
  filter(n()>=150)

tissues <- thresh_meta %>% 
  dplyr::select(group) %>%
  distinct() %>%
  pull() %>%
  tolower()

n=5

for (i in c(1:n)){
  dir.create(paste("./data/GTex_", i, "/", sep=''))
  
  sampled_meta <- thresh_meta %>% 
    group_by(group) %>%
    slice_sample(n=150)
  
  sampled_data <- data %>% 
    dplyr::select(gene_symbol, sampled_meta$sample_ID)
  
  comparisons <- character()
  for (tissue1 in tissues){
    for (tissue2 in tissues){
      if (tissue1==tissue2){
        break
      }
      else {
        contrast <- paste0(sort(c(tissue1, tissue2))[1],'_v_', sort(c(tissue1, tissue2))[2])
        if (contrast %in% comparisons){
          break
        }else{
          comparisons <- c(comparisons, contrast)
          sel_meta <- sampled_meta %>% filter(group == tissue1 | group == tissue2)
          contr_samples <- sel_meta %>% dplyr::select(sample_ID) %>% pull()
          sel_counts <- sampled_data %>% dplyr::select(gene_symbol,!!contr_samples) %>% group_by(gene_symbol) %>% summarise_all(sum)
          
          count_name <- paste("./data/GTex_", i, "/", paste(contrast, 'countdata.tsv', sep = '__'), sep='')
          meta_name <- paste("./data/GTex_", i, "/", paste(contrast, 'metadata.tsv', sep='__'), sep='')
          
          write_tsv(sel_meta, meta_name)
          write_tsv(sel_counts, count_name)
        }
      }
    }
  }
}
  


#metadata %>% 
#  group_by(group) %>% 
#  count() %>% 
#  ggplot(.,aes(x=group,y=n, fill=n<100)) + 
#  geom_col() +
#  scale_fill_discrete('blue', 'green') +
#  ggtitle('Threshold = 100') +
#  theme(axis.text.x = element_text(angle=-60, vjust = 0, hjust=0), plot.margin = margin(r=30)) +
#  geom_hline(aes(yintercept=100, color='red'))
