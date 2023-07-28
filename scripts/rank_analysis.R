library(tidyverse)
library(stats)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
args <- commandArgs(trailingOnly = FALSE)
path_file <- args[grep("--file=", args)] %>%
  sub("rank_analysis.R", "", .) %>%
  sub("--file=", "", .)
source(paste0(path_file, "rank_analysis_helper.R"))


args <- commandArgs(trailingOnly = FALSE)
func_datafile <- args[grep("--func_file", args) + 1]
de_datafile <- args[grep("--de_file", args) + 1]
dataset_id <- args[grep("--dataset",args)+1]


# func_datafile <- "./data/flop_results/funcomics/fullmerged/Spurell19__fullmerge.tsv"
# de_datafile <- "./data/flop_results/diffexp/Spurell19__deresults.tsv"

func_merged_data <- read_tsv(func_datafile)
func_cor_results <- corr_analysis(func_merged_data)


statparams <- func_merged_data %>% distinct(statparam) %>% pull()
status <- func_merged_data %>% distinct(status) %>% pull()
pipelines <- func_merged_data %>% distinct(pipeline) %>% pull()
bio_contexts <- func_merged_data %>% distinct(bio_context) %>% pull()


de_merged_data <- read_tsv(de_datafile)

de_cor_results <- tibble()
for(bio_context in bio_contexts){
  de_subset <- de_merged_data %>%
          select(ID, contains(bio_context), -contains('padj')) %>%
          pivot_longer(-ID, names_to = "runID", values_to = "scores") %>%
          separate(runID, into = c("statparam", "status", "pipeline", "bio_context", "main_dataset", "subset"), sep = "__")  %>%
          mutate(resource = 'DE') %>% 
          dplyr::rename('items' = 'ID') %>%
    filter(!is.na(scores))
  de_subset_cor_results <- corr_analysis(de_subset)
  de_cor_results <- bind_rows(de_cor_results, de_subset_cor_results)
}

cor_results <- bind_rows(func_cor_results, de_cor_results)

write_tsv(cor_results, file = paste0(dataset_id, "__total__rank.tsv"))




# DE data
