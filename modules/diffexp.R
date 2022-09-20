###Packages
library(tidyverse)
library(limma)

#' Limma differential expression analysis
#' This function performs a limma differnetial expression analysis with 
#' intersection and contrast. 
#' @param data_norm a tibble containing normalised values of gene counts, genes 
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
  genenames <- data_norm %>% select(gene_symbol)
  vsn_matrix <- data_norm %>% select(-gene_symbol)
  
  designmat <- metadata %>%
    select(Sample.Title, disease.status, fraction) %>%
    model.matrix(~ 0 + disease.status + fraction, data=.)
  
  fit <- lmFit(vsn_matrix, design = designmat)
  coefs="disease.statusSA-disease.statusHC" 
  fit <- coefs %>%
    makeContrasts(levels = designmat, contrasts = .) %>%
    contrasts.fit(fit = fit, contrasts = .) %>%
    eBayes(.)
  
  de_table <- topTable(fit, number = Inf, coef = coefs, sort.by='none') %>%
    as_tibble() %>%
    add_column(gene_symbol = genenames)
  
  return(de_table)
}