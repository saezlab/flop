library(tidyverse)
library(qs)

args <- commandArgs(trailingOnly = FALSE)
file_dir <- args[grep("--file_dir", args) + 1]

counts_file <- list.files(path = file_dir, pattern='countdata.tsv', full.names = TRUE)
meta_file <- list.files(path = file_dir, pattern='metadata.tsv', full.names = TRUE)
# the contrast file is optional. If not supplied, all possible contrasts will be created
contrast_file <- list.files(path = file_dir, pattern='contrast.tsv', full.names = TRUE)

subset_id <- file_dir %>% strsplit('/') %>% unlist() %>% tail(1)
# The counts should include a column named gene_symbol,
# and the samples should be in the columns
counts <- read_tsv(counts_file)
# The metadata should include a column named sample_ID
# that refer to the samples found in the counts file and a column named group
metadata <- read_tsv(meta_file)
# if the contrast file is provided, read it in. Otherwise, create all possible contrasts
if (length(contrast_file) != 0){
  contrast <- read_tsv(contrast_file) %>%
    mutate(samples = map2(group1, group2, ~ metadata %>% filter(group == .x | group == .y) %>%
    dplyr::select(sample_ID) %>%
    pull()))

  for (i in 1:nrow(contrast)){
    sel_groups <- contrast[i,] %>%
      dplyr::select(group1, group2) %>%
      c(., recursive = TRUE) %>% unname()

    sel_samples <- contrast[i,][3] %>%
      dplyr::select(samples) %>%
      c(., recursive = TRUE) %>% unname()

    sel_counts <- counts %>%
      dplyr::select(gene_symbol,!!sel_samples) %>%
      group_by(gene_symbol) %>%
      summarise_all(sum)

    sel_metadata <- metadata %>%
      filter(sample_ID %in% sel_samples)

    contrast_name <- paste(sel_groups[1],'_v_', sel_groups[2], sep = '')
    count_name <- paste(subset_id, contrast_name, 'countdata.qs', sep = '__')
    meta_name <- paste(subset_id, contrast_name, 'metadata.qs', sep='__')
    qsave(sel_metadata, meta_name)
    qsave(sel_counts, count_name)
  }
} else if (length(contrast_file) == 0) {
  groups <- metadata %>%
    dplyr::select(group) %>%
    distinct() %>%
    pull()

  comparisons <- c()
  for (group1 in groups){
    for (group2 in groups){
      if (group1 == group2){
        break
      }
      else {
        contrast_name <- paste0(sort(c(group1, group2))[1],'_v_', sort(c(group1, group2))[2])
        if (contrast_name %in% comparisons){
          break
        } else {
          comparisons <- c(comparisons, contrast_name)

          sel_metadata <- metadata %>%
            filter(group == group1 | group == group2)

          contr_samples <- sel_metadata %>%
            dplyr::select(sample_ID) %>% pull()

          sel_counts <- counts %>%
            dplyr::select(gene_symbol,!!contr_samples) %>%
            group_by(gene_symbol) %>%
            summarise_all(sum)

          count_name <- paste(subset_id, contrast_name, 'countdata.qs', sep = '__')
          meta_name <- paste(subset_id, contrast_name, 'metadata.qs', sep='__')
          qsave(sel_metadata, meta_name)
          qsave(sel_counts, count_name)
        }
      }
    }
  }
}


