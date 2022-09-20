###Packages
library(tidyverse)
library(limma)

#' Limma differential expression analysis
#' This function performs a limma differnetial expression analysis with or 
#' without intersection and contrast. 
#' @param data_norm a tibble containing normalised values of gene counts, genes 
#' as rows, samples as columns.
#' @param metadata a tibble containing the metadata related to the normalised
#' data.
#' @param intersection boolean value, indicates if an intersection must be taken
#' into account. Default is FALSE
#' @param contrast boolean value, indicates if contrast analysis is performed. 
#' Default is FALSE
#'
#' @return a tibble containing statistical parameters for each gene.
#' @export
#'
#' @examples
#' limma_anal(vsn_counts_exp1, metadata_exp1, intersection=TRUE, contrast=TRUE)
limma_anal <- function(data_norm, metadata, intersection=FALSE, contrast=FALSE){
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