#Packages
library(tidyverse)
library(vsn)
library(limma)
library(rstudioapi)
library(edgeR)

#Directory settings
source(this.dir <- dirname(parent.frame(2)$ofile))
source(setwd(this.dir))

#File imports
source('./functions/normalization.R')
source('./functions/diffexp.R')

###MAIN###
#Read files
transcriptdata <- read_tsv(file='./data/GSE85214-expression.txt')
metadata <- read_tsv(file='./data/GSE85214-metadata.txt') %>% 
  rename_all(make.names) %>%
  as_tibble()

#Prefiltering to remove low expressed genes. Threshold 500 counts
filt_data <- transcriptdata %>%
  filter(rowSums(select_if(., is.numeric))>500)

#Normalization
norm_data <- vsn_norm(filt_data)

#Differential expression analysis
inter_contr <- limma_anal(norm_data, metadata, TRUE, TRUE) %>% 
  select(gene_symbol, adj.P.Val) %>% 
  rename(adj_pvalue_inter_contr = adj.P.Val)
no_inter_no_contr <- limma_anal(norm_data, metadata, FALSE, FALSE) %>% 
  select(gene_symbol, adj.P.Val) %>% 
  rename(adj_pvalue_none = adj.P.Val)

#Results formatting
res <- left_join(inter_contr, no_inter_no_contr, by = 'gene_symbol')

#Statistical analysis
pearsoncorr <- res %>% select(-gene_symbol) %>% cor(.)
pearsoncorr