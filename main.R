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

diffanal_types = c('limma', 'edger')
norm_types <- c('vsn', 'log2quant', 'tmm')

#Normalization
norm_func <- function(filt_data, norm_types) {
  vsn_results <- vsn_norm(filt_data)
  log2quant_results <- log2quant_norm(filt_data)
  tmm_results <- tmm_norm(filt_data)
  norm_results <- tibble(norm = norm_types, 
                         data=list(vsn_results, log2quant_results, tmm_results))
  return(norm_results)
}

#Differential expression analysis
diff_anal <- function(norm_data, metadata, diffanal_types) {
  limma_results <- limma_anal(norm_data, metadata)
  edger_results <- edger_anal(norm_data, metadata)
  diffanal_results <- tibble(diffanal = diffanal_types, 
                             data=list(limma_results, edger_results))
  return(diffanal_results)
}

#benchmarking
norm_results <- norm_func(filt_data, norm_types)
diffanal_results <- list()
for (i in norm_types){
  norm_data <- norm_results %>% filter(norm==i) %>% unnest() %>% select(-norm)
  diffanal_results[[i]]<- diff_anal(norm_data, metadata, diffanal_types)
    }
diffanal_results
res <- filt_data %>% select(gene_symbol)
for (i in norm_types){
  for (j in diffanal_types){
    old_name <- paste('padj', j, sep='_')
    new_name <- paste('padj', i, j, sep='_')
    adj_val <- diffanal_results[[i]] %>% 
      filter(diffanal==j) %>% 
      unnest() %>% 
      select(gene_symbol, old_name) %>%
      dplyr::rename(., !!new_name := !!old_name)
    res <- left_join(res, adj_val, by='gene_symbol')
  }
}
deseq2_results <- deseq2_anal(filt_data, metadata) 
res <- deseq2_results %>% 
  select(gene_symbol, padj_deseq2) %>%
  left_join(res, ., by='gene_symbol')

#Statistical analysis
pearsoncor<- res %>% select(-gene_symbol) %>% cor(., use='complete.obs')

#Heatmap
pearsoncor %>% 
  as_tibble(rownames = "var1") %>% 
  gather(var2, value, -var1) %>% 
  
  ggplot(aes(x = var1, y = var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, digits = 2))) +
  labs(x = "", y = "", fill = "Corr", title = "Correlation Matrix") +
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_viridis_c(
    limits = c(0,1)
  )
          