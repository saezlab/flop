library(tidyverse)
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

args <- commandArgs(trailingOnly = FALSE)
dataset_id <- args[grep("--dataset",args)+1]
status <- args[grep("--status",args)+1]
file_directory <- args[grep("--file",args)+1][2]
pipeline_files <- file_directory %>%
  read_lines() %>%
  as_tibble()
run_info <- pipeline_files %>%
  separate(
    sep="__",
    col=value,
    into=c("datasetID", "bio_context",
    "statparam", "dcmethod",
    "resource", "ext")
  ) %>%
  select(-ext)
bio_contexts <- run_info %>%
  select(bio_context) %>%
  distinct() %>%
  pull()
resources <- run_info %>%
  select(resource) %>%
  distinct() %>%
  pull()
statparams <- run_info %>%
  select(statparam) %>%
  distinct() %>%
  pull()
merged_data <- tibble()

for(statparam in statparams){
  for(bio_context in bio_contexts){
    for(resource in resources){
      file_list <- pipeline_files %>%
        pull() %>%
          str_subset(.,pattern=paste(bio_context,"__", statparam,"__cons__", resource, "__decoupleroutput.tsv", sep=""))
      for(file in file_list) {
        file_data <- read_tsv(file, show_col_types=FALSE) %>%
            separate(...1, into=c("bio_context","status", "norm", "diffexp", "dcmethod"), sep="__", remove=FALSE) %>%
            dplyr::rename(obs_id = `...1`) %>%
            select(-dcmethod) %>%
            pivot_longer(cols=c(6:ncol(.)), names_to = "items", values_to = "scores") %>%
            mutate(statparam=!!statparam, resource=!!resource, pipeline=paste(norm, diffexp, sep = "+"))
        merged_data <- rbind(merged_data, file_data)
      }
    }
  }
}

filename <- paste(dataset_id, status, "result.tsv", sep = "__")
write_tsv(merged_data, filename)
