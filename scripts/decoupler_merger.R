library(tidyverse)
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

args <- commandArgs(trailingOnly = FALSE)
dataset_id <- args[grep("--dataset",args)+1]
status <- args[grep("--status",args)+1]
file_directory <- args[grep("--file",args)+1][2]

# Read all the DecoupleR output file paths written in the summary file from Nextflow 
# and merges them into a single file (datasetID__status__result.tsv)
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
dcmethods <- run_info %>%
  select(dcmethod) %>%
  distinct() %>%
  pull()
merged_data <- tibble()

for(statparam in statparams){
  for(bio_context in bio_contexts){
    for(resource in resources){
      for(dcmethod in dcmethods){
        file_list <- pipeline_files %>%
          pull() %>%
            str_subset(.,pattern=paste(bio_context,"__", statparam,"__", dcmethod, "__", resource, "__decoupleroutput.tsv", sep=""))
        for(file in file_list) {
          file_data <- read_tsv(file, show_col_types=FALSE) %>%
              separate(...1, into=c("bio_context","status", "norm", "diffexp", "dcmethod", 'metric'), sep="__", remove=FALSE) %>%
              dplyr::rename(obs_id = `...1`) %>%
              pivot_longer(cols=c(8:ncol(.)), names_to = "items", values_to = "scores") %>%
              select(-obs_id) %>%
              pivot_wider(names_from = metric, values_from = scores) %>%
              mutate(statparam=!!statparam, resource=!!resource, pipeline=paste(norm, diffexp, sep = "+"),
              obs_id = paste(bio_context,status, norm, diffexp, dcmethod, sep = "__")) %>%
              relocate(obs_id)
          merged_data <- rbind(merged_data, file_data)
        }
      }
    }
  }
}

filename <- paste(dataset_id, status, "result.tsv", sep = "__")
write_tsv(merged_data, filename)

