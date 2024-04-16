library(tidyverse)
library(edgeR)



#' @title Filtering helper
#' @description This function filters the genes based on the expression level
#' @param counts A dataframe containing the gene counts
#' @param metadata A dataframe containing the sample metadata
#' @return A dataframe containing the filtered gene counts
#' @importFrom dplyr %>% select arrange
#' @importFrom edgeR DGEList filterByExpr
#' 
filtering <- function(counts, metadata, additional_args) {

    libsize <- eval(parse(text = additional_args$filterbyexpr_libsize))
    mincount <- as.numeric(additional_args$filterbyexpr_mincount)
    mintotalcount <- as.numeric(additional_args$filterbyexpr_mintotalcount)
    large_number <- as.numeric(additional_args$filterbyexpr_largen)
    min_prop <- as.numeric(additional_args$filterbyexpr_minprop)


    covariates <- metadata %>%
        dplyr::select(-sample_ID) %>% 
        .[, sapply(., Negate(anyNA)), drop = FALSE] %>%
        sapply(., function(x) n_distinct(x)) %>% as.data.frame() %>% filter(. > 1) %>% rownames()

    designmat <- metadata %>%
        model.matrix(as.formula(rlang::parse_expr(paste(c("~ 0", covariates), collapse = " + "))), data=.)

    samples <- metadata %>%
        select(sample_ID, group) %>%
        arrange(group)

    dge_obj <- counts %>%
        column_to_rownames("gene_symbol") %>%
        DGEList(count = ., samples = samples)

    filtered_genes <- filterByExpr(dge_obj, 
                                    design = designmat,
                                    lib.size = libsize,
                                    min.count = mincount,
                                    min.total.count = mintotalcount,
                                    large.n = large_number,
                                    min.prop = min_prop)
    filtered_counts <- counts[filtered_genes,]
    return (filtered_counts)
}