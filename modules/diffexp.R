###Packages
library(tidyverse)
library(limma)
library(edgeR)
library(DESeq2)


#' Limma differential expression analysis
#' This function performs a limma differnetial expression analysis with 
#' intersection and contrast. 
#' @param data_norm a tibble containing normalised values of gene counts; genes 
#' as rows, samples as columns.
#' @param metadata a tibble containing the metadata related to the normalised
#' data.
#'
#' @return a tibble containing statistical parameters for each gene.
#' @importFrom limma lmFit makeContrasts contrasts.fit topTable
#' @importFrom dplyr %>% select add_column as_tibble 
#' @export
#'
#' @examples
#' limma_anal(vsn_counts_exp1, metadata_exp1)
limma_anal <- function(data_norm, metadata){
  norm_counts <- data_norm %>% select_if(is.numeric)
  
  designmat <- metadata %>%
    select(Sample.Title, disease.status, fraction) %>%
    model.matrix(~ 0 + disease.status + fraction, data=.)
  
  fit <- lmFit(norm_counts, design = designmat)
  coefs="disease.statusSA-disease.statusHC" 
  fit <- coefs %>%
    makeContrasts(levels = designmat, contrasts = .) %>%
    contrasts.fit(fit = fit, contrasts = .) %>%
    eBayes(.)
  
  results <- topTable(fit, number = Inf, coef = coefs, sort.by='none') %>%
    as_tibble() %>%
    mutate(gene_symbol = data_norm$gene_symbol) %>%
    relocate(gene_symbol) %>%
    dplyr::rename(logFC_limma = logFC, padj_limma = adj.P.Val)
  
  return(results)
}

#' EdgeR differential expression analysis
#'
#' @param data_norm a tibble containing normalised values of gene counts; genes 
#' as rows, samples as columns.
#' @param metadata a tibble containing the metadata related to the normalised
#' data.
#'
#' @return a tibble containing statistical parameters for each gene.
#' @importFrom edgeR estimateDisp glmFit glmLRT topTags DGElist
#' @importFrom dplyr %>% select as_tibble add_column
#' @export
#'
#' @examples
#' edger_anal(vsn_results, metadata_cancer)
edger_anal <- function(data_norm, metadata) {
  norm_counts <- data_norm %>% select_if(is.numeric) %>% DGEList(count=.)
  
  metadata$disease.status <- factor(metadata$disease.status, levels=c('HC', 'SA'))
  metadata$fraction <- factor(metadata$fraction, levels=c('Total', 'Polysome'))
  
  designmat <- metadata %>%
    select(Sample.Title, disease.status, fraction) %>%
    model.matrix(~ 0 + disease.status + fraction, data=.)
  
  y <- estimateDisp(norm_counts, designmat)
  fit <- glmQLFit(y, designmat)
  coefs<- makeContrasts(levels = designmat, 
                            contrasts = "disease.statusSA-disease.statusHC") 
  
  results <- glmQLFTest(fit, contrast=coefs) %>% 
    topTags(., n= Inf, sort.by='none') %>% 
    as.data.frame() %>% 
    as_tibble() %>%
    mutate(gene_symbol = data_norm$gene_symbol) %>%
    relocate(gene_symbol) %>%
    dplyr::rename(logFC_edger = logFC, padj_edger = PValue)
  
  return(results)
}

#' DESeq2 differential expression analysis
#'
#' @param raw_data a tibble containing raw gene counts; genes as rows, samples
#' as columns
#' @param metadata a tibble containing the metadata related to the normalised
#' data.
#'
#' @return a tibble containing statistical parameters for each gene.
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom dplyr %>% select as_tibble add_column
#' @export
#'
#' @examples
#' deseq2_anal(raw_counts, metadata_cancer)
deseq2_anal <- function(raw_data, metadata) {
  gene_counts <- raw_data %>% 
    as.data.frame(.)
  rownames(gene_counts) <- raw_data$gene_symbol
  gene_counts <- gene_counts %>% 
    select(-gene_symbol)
  metadata$disease.status <- factor(metadata$disease.status, levels=c('HC', 'SA'))
  metadata$fraction <- factor(metadata$fraction, levels=c('Total', 'Polysome'))
  
  formatted_data <- DESeqDataSetFromMatrix(countData=gene_counts, 
                                           colData = metadata, 
                                           design = ~ 0+ disease.status + fraction)
  
  results <- DESeq(formatted_data) %>% 
    results(., contrast=c('disease.status', 'SA', 'HC')) %>%
    as_tibble() %>%
    mutate(gene_symbol = raw_data$gene_symbol)  %>%
    relocate(gene_symbol) %>%
    dplyr::rename(logFC_deseq2 = log2FoldChange, padj_deseq2= padj)
  return(results)
}
