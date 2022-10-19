library(tidyverse)
library(limma)
library(edgeR)

###Functions###

#' TMM normalization
#' This function performs TMM normalization of the data using edgeR package
#' @param data tibble containing gene counts, genes as rows and
#' samples as columns
#'
#' @return a tibble containing the TMM normalised values of the data
#' @importFrom dplyr %>% select select_if as_tibble add_column
#' @importFrom edgeR DGEList calcNormFactors
#' @export
#'
#' @examples
#' tmm_norm(gene_counts)
tmm_norm <- function(data) {
  tmm_matrix <- data %>% 
    select_if(is.numeric) %>% 
    as.matrix() %>% 
    DGEList(count=.) %>% 
    calcNormFactors(., method="TMM") %>%
    cpm(.) %>%
    as_tibble(.) %>%
    mutate(ID = data$gene_symbol) %>%
    relocate(ID)
  return(tmm_matrix)
}

###Main###
args <- commandArgs(trailingOnly = FALSE)
counts_file <- args[grep("--counts",args)+1]
counts <- read.table(file = counts_file, header = TRUE, sep = "\t")
dataset_id <- strsplit(counts_file, split='_')[[1]][1]
results <- tmm_norm(counts)

output_filename <- paste(dataset_id, 'tmm', 'norm', sep = '__') %>%
  paste(., '.tsv', sep='')
write.table(results, output_filename, sep='\t', quote=FALSE, row.names=FALSE)
print('Done!')
