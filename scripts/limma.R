library(tidyverse)
library(limma)
###Functions###

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
#' limma_anal(vsn_counts_exp1, metadata_exp1)
limma_anal <- function(data_norm, metadata){
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
    mutate(ID = data_norm$gene_symbol) %>%
    dplyr::rename(padj = adj.P.Val) %>%
    select(ID, logFC, B, padj)
  
  return(results)
}

###Main###
args <- commandArgs(trailingOnly = FALSE)
norm_file <- args[grep("--norm",args)+1]
meta_file <- args[grep("--meta",args)+1]
norm_counts <- read.table(file = norm_file, header = TRUE, sep = "\t")
metadata <- read.table(file = meta_file, header = TRUE, sep = "\t", stringsAsFactors=TRUE)
results <- limma_anal(norm_counts, metadata)

write.table(results, 'limma_results.tsv', sep='\t', quote=FALSE, row.names=FALSE)
print('Done!')
