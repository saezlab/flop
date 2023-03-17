library(tidyverse)
library(edgeR)

# Insert roxygen skeleton:

#' @title Filtering helper
#' @description This function filters the genes based on the expression level
#' @param counts A dataframe containing the gene counts
#' @param metadata A dataframe containing the sample metadata
#' @return A dataframe containing the filtered gene counts
#' @importFrom dplyr %>% select arrange
#' @importFrom edgeR DGEList filterByExpr
#' 
filtering <- function(counts, metadata) {
    samples <- metadata %>%
        select(sample_ID, group) %>%
        arrange(group)

    dge_obj <- counts %>%
        column_to_rownames("gene_symbol") %>%
        DGEList(count = ., samples = samples)

    filtered_genes <- filterByExpr(dge_obj)
    filtered_counts <- counts[filtered_genes,]
    return (filtered_counts)
}