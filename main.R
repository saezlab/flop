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

res_parser <- function(diffanal_results, 
                       stat_param,
                       norm_met=norm_types, 
                       diffanal_met=diffanal_types){
  res <- diffanal_results$tmm$data[[1]] %>% 
    select(gene_symbol)
  deseq_log <- FALSE
  for (i in norm_met){
    for (j in diffanal_met){
      if(j=='deseq2'){
        if(deseq_log){
          next
        } else {
          new_name <- paste(stat_param, j, sep='_')
          param <- diffanal_results[[i]] %>% 
            filter(diffanal==j) %>% 
            unnest() %>% 
            select(gene_symbol, new_name)
          res <- left_join(res, param, by='gene_symbol')
          deseq_log <- TRUE
        }
      } else {
        old_name <- paste(stat_param, j, sep='_')
        new_name <- paste(stat_param, i, j, sep='_')
        param <- diffanal_results[[i]] %>% 
          filter(diffanal==j) %>% 
          unnest() %>% 
          select(gene_symbol, old_name) %>%
          dplyr::rename(., !!new_name := !!old_name)
        res <- left_join(res, param, by='gene_symbol')
      }
    }
  }
  return(res)
}

corr_plotter <- function(corr_res, stat_param, cor_type){
  corr_res %>% 
    as_tibble(rownames = "var1") %>% 
    gather(var2, value, -var1) %>% 
    
    ggplot(aes(x = var1, y = var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, digits = 2))) +
    labs(x = "", y = "", fill = "Corr", title = paste("Correlation Matrix:", stat_param, cor_type, sep=' ')) +
    theme(axis.text.x = element_text(angle = 90))+
    scale_fill_viridis_c(
      limits = c(0,1)
    )
}

correlation_heatmap <- function(diffanal_results, stats_param, cor_type = 'pearson'){
  res <- res_parser(diffanal_results, stats_param)
  corr_res <- res %>% 
    select(-gene_symbol) %>% 
    cor(., use='complete.obs', method=cor_type)
  corr_plotter(corr_res, stats_param, cor_type)
}

###MAIN###
#Read files
transcriptdata <- read_tsv(file='./data/GSE85214-expression.txt')
metadata <- read_tsv(file='./data/GSE85214-metadata.txt') %>% 
  rename_all(make.names) %>%
  as_tibble()

#Prefiltering to remove low expressed genes. Threshold 500 counts
filt_data <- transcriptdata %>%
  filter(rowSums(select_if(., is.numeric))>500)

diffanal_types = c('limma', 'edger', 'deseq2')
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
deseq2_results <- deseq2_anal(filt_data, metadata) 

#benchmarking
norm_results <- norm_func(filt_data, norm_types)
diffanal_results <- list()
for (i in norm_types){
  norm_data <- norm_results %>% 
    filter(norm==i) %>% 
    unnest() %>% 
    select(-norm)
  diffanal_results[[i]]<- diff_anal(norm_data, metadata, diffanal_types)
}
diffanal_results[['vsn']] <- diffanal_results[['vsn']] %>% 
  add_row(diffanal = "deseq2", data= list(deseq2_results))

#Heatmap
correlation_heatmap(diffanal_results, 'padj', 'spearman')
correlation_heatmap(diffanal_results, 'padj', 'pearson')
correlation_heatmap(diffanal_results, 'logFC', 'spearman')
correlation_heatmap(diffanal_results, 'logFC', 'pearson')


  
