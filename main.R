#Packages
library(tidyverse)
library(vsn)
library(limma)
library(rstudioapi)
library(edgeR)

#Directory settings
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#File imports
source('./modules/normalization.R')
source('./modules/diffexp.R')

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
vsn_results <- vsn_norm(filt_data)
log2quant_results <- log2quant_norm(filt_data)
tmm_results <- tmm_norm(filt_data)

#Differential expression analysis
limma_results <- limma_anal(vsn_results, metadata)

edger_results <- edger_anal(vsn_results, metadata)

#Results formatting
res <- left_join(inter_contr, no_inter_no_contr, by = 'gene_symbol')

#Statistical analysis
pearsoncorr <- res %>% select(-gene_symbol) %>% cor(.)
pearsoncorr