library(tidyverse)

args <- commandArgs(trailingOnly = FALSE)
path_file <- args[grep("--file=", args)] %>%
  sub("top_bottom_overlap_analysis.R", "", .) %>%
  sub("--file=", "", .)
source(paste0(path_file, "top_bottom_overlap_analysis_helper.R"))
func_datafile <- args[grep("--func_file", args) + 1]
de_datafile <- args[grep("--de_file", args) + 1]
dataset_id <- args[grep("--dataset",args)+1]
pval_cutoff <- as.numeric(args[grep("--pval_thresh",args)+1])

# Top-bottom overlap analysis for gene set values
func_merged_data <- read_tsv(func_datafile)
func_jaccard_results <- jaccard_analysis(func_merged_data, pval_cutoff)

statparams <- func_merged_data %>% distinct(statparam) %>% pull()
status <- func_merged_data %>% distinct(status) %>% pull()
pipelines <- func_merged_data %>% distinct(pipeline) %>% pull()
bio_contexts <- func_merged_data %>% distinct(bio_context) %>% pull()


de_merged_data <- read_tsv(de_datafile)

# Top-bottom overlap analysis for the DE space
de_jaccard_results <- tibble()
for(bio_context in bio_contexts){
  de_subset_act <- de_merged_data %>%
          select(ID, contains(bio_context), -contains('padj')) %>%
          pivot_longer(-c(ID), names_to = "runID", values_to = "act") %>%
          separate(runID, into = c("statparam", "status", "pipeline", "bio_context", "main_dataset", "subset"), sep = "__")  %>% 
          mutate(resource = 'DE') %>% 
          dplyr::rename('items' = 'ID')
  de_subset_padj <- de_merged_data %>%
        select(ID, contains(bio_context) & contains('padj')) %>%
        pivot_longer(-c(ID), names_to = "runID", values_to = "padj") %>%
        separate(runID, into = c("statparam", "status", "pipeline", "bio_context", "main_dataset", "subset"), sep = "__")  %>% 
        select(-statparam) %>%
        mutate(resource = 'DE') %>% 
        dplyr::rename('items' = 'ID')
  de_subset <- full_join(de_subset_act, de_subset_padj, by = c("status", "pipeline", "bio_context", "main_dataset", "subset", "items", "resource")) %>%
  filter(!is.na(act) & !is.na(padj))
  de_subset_jaccard_results <- jaccard_analysis(de_subset, pval_cutoff)
  de_jaccard_results <- bind_rows(de_jaccard_results, de_subset_jaccard_results)
}

jaccard_results <- bind_rows(func_jaccard_results, de_jaccard_results)

write_tsv(jaccard_results, file = paste0(dataset_id, "__overlap.tsv"))