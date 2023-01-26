library(tidyverse)
library(edgeR)


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