library(tidyverse)
library(limma)
library(edgeR)
###Functions###

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
  coefs <- paste("group", studygroups[1], "-", "group", studygroups[2], sep='')
  
  dge_obj <- calcNormFactors(dge_obj)

  contrast_mat<- makeContrasts(levels = designmat, 
                        contrasts = coefs) 
  dge_obj <- estimateDisp(dge_obj, designmat)
  
  fit <- glmQLFit(dge_obj, designmat)
  
  results <- glmQLFTest(fit, contrast=contrast_mat) %>% 
    topTags(., n= Inf, sort.by='none', adjust.method='BH') %>% 
    as.data.frame() %>% 
    as_tibble() %>%
    mutate(ID = counts$gene_symbol) %>%
    dplyr::rename(padj = FDR, stat = F) %>%
    select(ID, logFC, stat, padj)
  
  return(results)
}

###Main###
args <- commandArgs(trailingOnly = FALSE)
counts_file <- args[grep("--counts",args)+1]
meta_file <- args[grep("--meta",args)+1]
counts <- read.table(file = counts_file, header = TRUE, sep = "\t")
metadata <- read.table(file = meta_file, header = TRUE, sep = "\t", stringsAsFactors=TRUE)
dataset_id <- strsplit(counts_file, split='_')[[1]][1]
results <- edger_analysis(counts, metadata)

output_filename <- paste(dataset_id, 'NA', 'edger', 'de', sep = '__') %>%
  paste(., '.tsv', sep='')
write.table(results, output_filename, sep='\t', quote=FALSE, row.names=FALSE)
print('Done!')
