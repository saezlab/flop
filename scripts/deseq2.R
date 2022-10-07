library(tidyverse)
library(limma)
library(DESeq2)
###Functions###

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
#' deseq2_anal(raw_counts, metadata_cancer)
deseq2_anal <- function(counts, metadata) {
  gene_counts <- counts %>% 
    as.data.frame(.)
  rownames(gene_counts) <- counts$gene_symbol
  gene_counts <- gene_counts %>% 
    select(-gene_symbol)
  
  formatted_data <- DESeqDataSetFromMatrix(countData=gene_counts, 
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


###Main###
args <- commandArgs(trailingOnly = FALSE)
counts_file <- args[grep("--counts",args)+1]
meta_file <- args[grep("--meta",args)+1]
counts <- read.table(file = counts_file, header = TRUE, sep = "\t")
metadata <- read.table(file = meta_file, header = TRUE, sep = "\t", stringsAsFactors=TRUE)
results <- deseq2_anal(counts, metadata)

write.table(results, 'deseq2_results.tsv', sep='\t', quote=FALSE, row.names=FALSE)
print('Done!')
