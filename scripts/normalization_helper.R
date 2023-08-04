library(tidyverse)
library(vsn)
library(limma)
library(edgeR)

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
vsn_norm <- function(data, ...){
  vsn_matrix <- data %>%
    select_if(is.numeric) %>%
    as.matrix(.) %>%
    justvsn(.) %>%
    as_tibble(.) %>%
    mutate(ID = data$gene_symbol) %>%
    relocate(ID)
  return(vsn_matrix)
}

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
tmm_norm <- function(data, ...) {
  tmm_matrix <- data %>%
    select_if(is.numeric) %>%
    as.matrix() %>%
    DGEList(count = .) %>%
    calcNormFactors(., method = "TMM") %>%
    cpm(., log = T, prior.count = 3 ) %>%
    as_tibble(.) %>%
    mutate(ID = data$gene_symbol) %>%
    relocate(ID)
  return(tmm_matrix)
}

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
log2quant_norm <- function(data, ...) {
  log2quant_matrix <- data %>%
    select_if(is.numeric)  %>%
    as.matrix() %>%
    +1 %>%
    log2() %>%
    normalizeQuantiles() %>%
    as_tibble() %>%
    mutate(ID = data$gene_symbol) %>%
    relocate(ID)
  return(log2quant_matrix)
}

#' voom normalization
#' This function performs a voom transformation followed by a
#' normalization of the data
#' @param data a tibble containing gene counts, genes as rows and
#' samples as columns
#' @param metadata a tibble containing the sample_ID and group#'
#' @return a tibble containing the voom normalised values of the data
#' @importFrom dplyr %>% select select_if as_tibble add_column relocate
#' @importFrom limma voom
#' @export
#'
#' @examples
#' voom_norm(gene_counts, metadata)
voom_norm <- function(data, metadata) {
  covariates <- metadata %>%
    dplyr::select(-sample_ID) %>%
    sapply(., function(x) n_distinct(x)) %>%
    as.data.frame() %>%
    filter(. > 1) %>%
    rownames()

  designmat <- metadata %>%
    model.matrix(as.formula(rlang::parse_expr(paste(c("~ 0", covariates), collapse = " + "))), data=.)

  voom_mat <- data %>%
    select_if(is.numeric) %>%
    as.matrix() %>%
    voom(., designmat) %>%
    .$E %>%
    as_tibble(.) %>%
    mutate(ID = data$gene_symbol) %>%
    relocate(ID)
  return(voom_mat)
}
