library(tidyverse)
library(limma)

###Functions###

#' Log2 quantile normalization
#' This function performs a log2 transformation followed by a quantile
#' normalization of the data
#' @param data a tibble containing gene counts, genes as rows and
#' samples as columns
#'
#' @return a tibble containing the log2 + quantile normalised values of the data
#' @importFrom dplyr %>% select select_if as_tibble add_column relocate
#' @importFrom limma normalizeQuantiles
#' @export
#'
#' @examples
#' log2quant_norm(gene_counts)
log2quant_norm <- function(data) {
  log2quant_matrix <- data %>% 
    select_if(is.numeric)  %>% 
    as.matrix() %>% 
    +1 %>% 
    log2() %>% 
    normalizeQuantiles() %>%
    as_tibble() %>%
    mutate(gene_symbol = data$gene_symbol) %>%
    relocate(gene_symbol)
  return(log2quant_matrix)
}

###Main###
args <- commandArgs(trailingOnly = FALSE)
counts_file <- args[grep("--counts",args)+1]
counts <- read.table(file = counts_file, header = TRUE, sep = "\t")

results <- log2quant_norm(counts)

write.table(results, 'log2quant_norm.tsv', sep='\t', quote=FALSE, row.names=FALSE)
print('Done!')
