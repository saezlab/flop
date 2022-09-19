this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

###Packages
library(tidyverse)
library(vsn)
library(limma)
library(ggfortify)
library(rstudioapi)
library(edgeR)

###NORMALIZATION###
#vsn normalization
#This function performs vsn normalization of the data.
vsn_norm <- function(data){
  genenames <- data %>% select(gene_symbol)
  vsn_matrix <- data %>% 
    select(-gene_symbol) %>% 
    as.matrix(.) %>% 
    justvsn(.) %>%
    as_tibble(.) %>%
    add_column(gene_symbol = genenames)
  return(vsn_matrix)
}

#Log2 quantile normalization
#This function performs a log2 transformation followed by a quantile 
#normalization of the data.
log2quant_norm <- function(data) {
  genenames <- data %>% select(gene_symbol)
  log2quant_matrix <- data %>% 
    select_if(is.numeric)  %>% 
    as.matrix() %>% 
    +1 %>% 
    log2() %>% 
    normalizeQuantiles() %>%
    as_tibble() %>%
    add_column(gene_symbol = genenames) %>%
    relocate(gene_symbol)
  return(log2quant_matrix)
}

#TMM normalization
#This function performs TMM normalization of the data.
tmm_norm <- function(data) {
  genenames <- data %>% select(gene_symbol)
  tmm_matrix <- data %>% 
    select_if(is.numeric) %>% 
    as.matrix() %>% 
    DGEList(count=.) %>% 
    calcNormFactors(., method="TMM") %>%
    .$counts %>%
    as_tibble(.) %>%
    add_column(gene_symbol = genenames) %>%
    relocate(gene_symbol)
  return(tmm_matrix)
}

###dIFFERENTIAL EXPRESSION ANALYSIS###
#This function takes as input a long-format tibble with all the gene expression 
#levels plus the metadata. It allows limma analysis with or without intersection
#and contrast. Returns a tibble containing the full results of the analysis
#indexed by gene symbol.
diff_anal <- function(data_norm, metadata, intersection=FALSE, contrast=FALSE){
  genenames <- data_norm %>% select(gene_symbol)
  vsn_matrix <- data_norm %>% select(-gene_symbol)
  
  designmat <- metadata %>%
    select(Sample.Title, disease.status, fraction) %>%
    { if(intersection) model.matrix(~ 0 + disease.status + fraction, data=.)
      else model.matrix(~disease.status + fraction, data=.)
    }
  
  fit <- lmFit(vsn_matrix, design = designmat)
  
  {if(contrast)
    fit <- "disease.statusSA-disease.statusHC" %>%
    makeContrasts(levels = designmat, contrasts = .) %>%
    contrasts.fit(fit = fit, contrasts = .) %>%
    eBayes(.)
    else fit <- eBayes(fit)
  }
  
  if(contrast) coefs="disease.statusSA-disease.statusHC" 
  else coefs="disease.statusSA"
  de_table <- topTable(fit, number = Inf, coef = coefs, sort.by='none') %>%
    as_tibble() %>%
    add_column(gene_symbol = genenames)
  
  return(de_table)
  
}

###MAIN###
#read files
transcriptdata <- read_tsv(file='./data/GSE85214-expression.txt')
metadata <- read_tsv(file='./data/GSE85214-metadata.txt') %>% 
  rename_all(make.names) %>%
  as_tibble()
#Prefiltering to remove low expressed genes. Threshold 500 counts
filt_data <- transcriptdata %>%
  filter(rowSums(select_if(., is.numeric))>500)
#Normalization
norm_data <- vsn_norm(filt_data)

inter_contr <- diff_anal(norm_data, metadata, TRUE, TRUE) %>% 
  select(gene_symbol, adj.P.Val) %>% 
  rename(adj_pvalue_inter_contr = adj.P.Val)
no_inter_no_contr <- diff_anal(norm_data, metadata, FALSE, FALSE) %>% 
  select(gene_symbol, adj.P.Val) %>% 
  rename(adj_pvalue_none = adj.P.Val)

res <- reduce(list(inter_contr, inter, no_inter_no_contr), dplyr::left_join, by = 'gene_symbol')
pearsoncorr <- res %>% select(-gene_symbol) %>% cor(.)
pearsoncorr