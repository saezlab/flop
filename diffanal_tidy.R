this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

###Packages
library(tidyverse)
library(vsn)
library(limma)
library(ggfortify)
library(rstudioapi)    

###FUNCTIONS###
#vsn normalization
#This function performs vsn normalization of the data. It takes a wide-format 
#tibble and returns a matrix. 
vsn_norm <- function(data, plots=FALSE){
  vsn_matrix <- data %>% select(-gene_symbol) %>% as.matrix(.) %>% justvsn(.)
  if(plots == TRUE){
    meanSdPlot(vsn_matrix, ranks=TRUE)
  }
  return(vsn_matrix)
}


#Differential expression analysis
#This function takes as input a long-format tibble with all the gene expression 
#levels plus the metadata. It allows limma analysis with or without intersection
#and contrast. Returns a tibble containing the full results of the analysis
#indexed by gene symbol.
#' Title
#'
#' @param total_data 
#' @param intersection 
#' @param contrast 
#'
#' @return
#' @export
#'
#' @examples
diff_anal <- function(total_data, intersection=FALSE, contrast=FALSE){
  data <- total_data %>% 
    select(gene_symbol, expr_lev, sample) %>% 
    pivot_wider(., names_from = sample, values_from = expr_lev)
  genenames <- data %>% select(gene_symbol)
  vsnres <- vsn_norm(data)

  designmat <- total_data %>%
    select(sample, disease.status, fraction) %>%
    distinct(sample, .keep_all = TRUE) %>%
    { if(intersection) model.matrix(~ 0 + disease.status + fraction, data=.)
      else model.matrix(~disease.status + fraction, data=.)
        }
  
  fit <- lmFit(vsnres, design = designmat)
  
  {if(contrast)
    fit <- "disease.statusSA-disease.statusHC" %>%
    makeContrasts(levels = designmat, contrasts = .) %>%
    contrasts.fit(fit = fit, contrasts = .) %>%
    eBayes(.)
    else fit <- eBayes(fit)
  }
  
  if(contrast) coefs="disease.statusSA-disease.statusHC" else coefs=NULL
  de_table <- topTable(fit, number = Inf, coef = coefs, sort.by='none') %>%
    as_tibble() %>%
    mutate(gene_symbol = genenames)
  
  saveRDS(designmat, filename)
  
  return(de_table)
  
}

###MAIN###
#read files
transcriptdata <- read_tsv(file='./Datasets/GSE85214-expression.txt') %>%
  pivot_longer(., cols= -gene_symbol, names_to = 'sample', values_to='expr_lev')
total_data <- read_tsv(file='./Datasets/GSE85214-metadata.txt') %>%
  left_join(transcriptdata, ., by=c('sample'='Sample_geo_accession'))%>%
  rename(., 'disease.status' = 'disease status')
#Prefiltering to remove low expressed genes. Threshold 500 counts
filt_data <- total_data %>% 
  group_by(gene_symbol) %>% 
  filter(sum(expr_lev) >500) %>%
  ungroup()

inter_contr <- diff_anal(filt_data, TRUE, TRUE) %>% 
  select(gene_symbol, adj.P.Val) %>% 
  rename(adj_pvalue_inter_contr = adj.P.Val)
inter <- diff_anal(filt_data, TRUE, FALSE) %>% 
  select(gene_symbol, adj.P.Val) %>% 
  rename(adj_pvalue_inter = adj.P.Val)
no_inter_no_contr <- diff_anal(filt_data, FALSE, FALSE) %>% 
  select(gene_symbol, adj.P.Val) %>% 
  rename(adj_pvalue_none = adj.P.Val)

res <- reduce(list(inter_contr, inter, no_inter_no_contr), dplyr::left_join, by = 'gene_symbol')
pearsoncorr <- res %>% select(-gene_symbol) %>% cor(.)
pearsoncorr