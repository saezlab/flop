library(tidyverse)
library(stats)
library(cowplot)
library(egg)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Rank Analysis: correlation
corr_analysis <- function(data, corrparam) {
  cor_results <- merged_data %>%
    group_by(statparam, resource, bio_context, status, main_dataset) %>%
    group_split() %>% 
    purrr::map(., function(x) {
      to_cor <- x %>%
        select(!!corrparam, scores, items) %>%
        pivot_wider(
          names_from = !!corrparam,
          values_from = scores,
          values_fn = {mean}
        ) %>%
        column_to_rownames(var = "items") %>%
        select(order(colnames(.))) %>%
        as.matrix()

      cor_results <- cor(to_cor, method = "spearman") %>%
        as.data.frame() %>%
        rownames_to_column(var = "feature_1") %>%
        pivot_longer(-feature_1) %>%
        mutate(
          statparam = unique(x$statparam),
          bio_context = unique(x$bio_context),
          resource = unique(x$resource),
          status = unique(x$status),
          main_dataset = unique(x$main_dataset)
        )

      return(cor_results)

    }) %>%
    bind_rows() %>%
    subset(feature_1 != name) %>%
    rowwise() %>%
    mutate(
      id = paste0(
        sort(c(feature_1, name))[1],
        " - ",
        sort(c(feature_1, name))[2]
      )
    ) %>%
    distinct(
      id,
      statparam,
      bio_context,
      resource,
      status,
      main_dataset,
      .keep_all = TRUE)
  return(cor_results)
}

args <- commandArgs(trailingOnly = FALSE)
datafile <- args[grep("--file", args) + 1][2]
dataset_id <- args[grep("--dataset",args)+1]

merged_data <- read_tsv(datafile)
cor_results <- corr_analysis(merged_data, "pipeline")
write_tsv(cor_results, file = paste0(dataset_id, "__total__rank.tsv"))
#Correlation split
#By norm
norm_cor_results <- corr_analysis(merged_data, "norm")
write_tsv(norm_cor_results, file = paste0(dataset_id, "__norm__rank.tsv"))
#By diffexp
diffexp_cor_results <- corr_analysis(merged_data, "diffexp")
write_tsv(diffexp_cor_results, file = paste0(dataset_id, "__diffexp__rank.tsv"))



