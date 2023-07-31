library(tidyverse)

args <- commandArgs(trailingOnly = FALSE)
path_file <- args[grep("--file=", args)] %>%
  sub("jaccard_analysis.R", "", .) %>%
  sub("--file=", "", .)
source(paste0(path_file, "jaccard_analysis_helper.R"))

args <- commandArgs(trailingOnly = FALSE)
func_datafile <- args[grep("--func_file", args) + 1]
de_datafile <- args[grep("--de_file", args) + 1]
dataset_id <- args[grep("--dataset",args)+1]

func_merged_data <- read_tsv(func_datafile)
func_jaccard_results <- jaccard_analysis(func_merged_data)

write_tsv(func_jaccard_results, file = paste0("GSE186341__jaccard__func.tsv"))

statparams <- func_merged_data %>% distinct(statparam) %>% pull()
status <- func_merged_data %>% distinct(status) %>% pull()
pipelines <- func_merged_data %>% distinct(pipeline) %>% pull()
bio_contexts <- func_merged_data %>% distinct(bio_context) %>% pull()



de_merged_data <- read_tsv(de_datafile)

de_jaccard_results <- tibble()
for(bio_context in bio_contexts){
  de_subset <- de_merged_data %>%
          select(ID, contains(bio_context), -contains('padj')) %>%
          pivot_longer(-ID, names_to = "runID", values_to = "act") %>%
          separate(runID, into = c("statparam", "status", "pipeline", "bio_context", "main_dataset", "subset"), sep = "__")  %>% 
          mutate(resource = 'DE') %>% 
          dplyr::rename('items' = 'ID') %>%
    filter(!is.na(act))
  de_subset_jaccard_results <- jaccard_analysis(de_subset)
  de_jaccard_results <- bind_rows(de_jaccard_results, de_subset_jaccard_results)
}

jaccard_results <- bind_rows(func_jaccard_results, de_jaccard_results)

write_tsv(jaccard_results, file = paste0(dataset_id, "__jaccard.tsv"))