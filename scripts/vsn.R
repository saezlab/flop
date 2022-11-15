library(tidyverse)
library(vsn)

###Functions###

#' Vsn normalization
#' This function performs vsn normalization of the data.
#' @param data a tibble containing gene counts, genes as rows and
#' samples as columns
#'
#' @return a tibble containing the vsn normalised values of the data
#' @importFrom dplyr %>% select as_tibble add_column
#' @importFrom limma justvsn
#' @export
#'
#' @examples
#' vsn_norm(gene_counts)
vsn_norm <- function(data){
  vsn_matrix <- data %>% 
    select_if(is.numeric) %>% 
    as.matrix(.) %>% 
    justvsn(.) %>%
    as_tibble(.) %>%
    mutate(ID = data$gene_symbol) %>%
    relocate(ID)
  return(vsn_matrix)
}

###Main###
args <- commandArgs(trailingOnly = FALSE)
counts_file <- args[grep("--counts",args)+1]
counts <- read.table(file = counts_file, header = TRUE, sep = "\t")
dataset_id <- strsplit(counts_file, split='__')[[1]][1]
results <- vsn_norm(counts)

output_filename <- paste(dataset_id, 'vsn', 'norm', sep = '__') %>%
  paste(., '.tsv', sep='')
write.table(results, output_filename, sep='\t', quote=FALSE, row.names=FALSE)
print('Done!')
