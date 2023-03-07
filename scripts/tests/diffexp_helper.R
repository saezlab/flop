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
#' data containing the following columns: sample_ID, group (2 values only) 
#'
#' @return a tibble containing statistical parameters for each gene.
#' @importFrom limma lmFit makeContrasts contrasts.fit topTable
#' @importFrom dplyr %>% select add_column as_tibble 
#' @export
#'
#' @examples
#' limma_analysis(vsn_counts_exp1, metadata_exp1)
limma_analysis <- function(data_norm, metadata){
  norm_counts <- data_norm %>% select_if(is.numeric)
  
  metadata <- metadata %>%
    select(sample_ID, group)
  
  designmat <- metadata %>%
    select(sample_ID, group) %>%
    model.matrix(~ 0 + group, data=.)
  studygroups <- levels(metadata$group)
  
  coefs <- paste("group", studygroups[1], "-", "group", studygroups[2], sep='')
  
  fit <- lmFit(norm_counts, design = designmat)
  fit <- coefs %>%
    makeContrasts(levels = designmat, contrasts = .) %>%
    contrasts.fit(fit = fit, contrasts = .) %>%
    eBayes(.)
  
  results <- topTable(fit, number = Inf, coef = coefs, sort.by='none') %>%
    as_tibble() %>%
    mutate(ID = data_norm$ID) %>%
    dplyr::rename(padj = adj.P.Val, stat = t) %>%
    select(ID, logFC, stat, padj)
  
  return(results)
}

#' DESeq2 differential expression analysis
#'
#' @param counts a tibble containing non-normalised gene counts; genes as rows, samples
#' as columns
#' @param metadata a tibble containing the metadata related to the normalised
#' data containing the following columns: sample_ID, group (2 values only) 
#'
#' @return a tibble containing statistical parameters for each gene.
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom dplyr %>% select as_tibble add_column
#' @export
#'
#' @examples
#' deseq2_analysis(raw_counts, metadata_cancer)
deseq2_analysis <- function(counts, metadata) {
  gene_counts <- counts %>% 
    as.data.frame(.)
  rownames(gene_counts) <- counts$gene_symbol
  gene_counts <- gene_counts %>% 
    select(-gene_symbol)
  
  formatted_data <- DESeqDataSetFromMatrix(countData = gene_counts, 
                                           colData = metadata, 
                                           design = ~ 0+ group)
  studygroups <- levels(metadata$group)
  coefs <- paste("group", studygroups[1], "-", "group", studygroups[2], sep='')
  
  results <- DESeq(formatted_data) %>% 
    results(., contrast=c('group', studygroups[1], studygroups[2])) %>%
    as_tibble() %>%
    mutate(ID = counts$gene_symbol)  %>%
    dplyr::rename(logFC = log2FoldChange) %>%
    select(ID, logFC, stat, padj)
  return(results)
}

#' EdgeR differential expression analysis
#'
#' @param counts a tibble containing non-normalised values of gene counts; genes 
#' as rows, samples as columns.
#' @param metadata a tibble containing the metadata related to the normalised
#' data containing the following columns: sample_ID, group (2 values only) 
#'
#' @return a tibble containing statistical parameters for each gene.
#' @importFrom edgeR estimateDisp glmFit glmLRT topTags DGElist
#' @importFrom dplyr %>% select as_tibble add_column
#' @export
#'
#' @examples
#' edger_analysis(vsn_results, metadata_cancer)
edger_analysis <- function(counts, metadata) {
  dge_obj <- counts %>% select_if(is.numeric) %>% DGEList(count=.)
  
  designmat <- metadata %>%
    select(sample_ID, group) %>%
    model.matrix(~ 0 + group, data=.)
  
  studygroups <- levels(metadata$group)
  coefs <- paste("group", studygroups[1], "-", "group", studygroups[2], sep = '')
  
  dge_obj <- calcNormFactors(dge_obj)

  contrast_mat<- makeContrasts(levels = designmat, 
                        contrasts = coefs) 
  dge_obj <- estimateDisp(dge_obj, designmat)
  
  fit <- glmQLFit(dge_obj, designmat)
  
  results <- glmQLFTest(fit, contrast = contrast_mat) %>% 
    topTags(., n= Inf, sort.by = "none", adjust.method="BH") %>% 
    as.data.frame() %>% 
    as_tibble() %>%
    mutate(ID = counts$gene_symbol) %>%
    dplyr::rename(padj = FDR) %>%
    mutate(stat = sign(logFC) * sqrt(F)) %>%
    select(ID, logFC, stat, padj)
  
  return(results)
}

