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

###Main###
args <- commandArgs(trailingOnly = FALSE)
counts_file <- args[grep("--counts",args)+1]
meta_file <- args[grep("--meta",args)+1]
dataset_id <- strsplit(counts_file, split='_')[[1]][1]
counts <- read.table(file = counts_file, header = TRUE, sep = "\t", row.names=NULL)
metadata <- read.table(file = meta_file, header = TRUE, sep = "\t", stringsAsFactors=TRUE, row.names=NULL)
results <- deseq2_analysis(counts, metadata)

output_filename <- paste(dataset_id, 'NA', 'deseq2', 'de', sep = '__') %>%
  paste(., '.tsv', sep='')
write.table(results, output_filename, sep='\t', quote=FALSE, row.names=FALSE)
print('Done!')
