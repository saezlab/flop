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
    mutate(gene_symbol = data$gene_symbol) %>%
    relocate(gene_symbol)
  return(vsn_matrix)
}

###Main###
args <- commandArgs(trailingOnly = FALSE)
counts_file <- args[grep("--counts",args)+1]
counts <- read.table(file = counts_file, header = TRUE, sep = "\t")

results <- vsn_norm(counts)

write.table(results, 'vsn_norm.tsv', sep='\t', quote=FALSE, row.names=FALSE)
print('Done!')
