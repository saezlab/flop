library(tidyverse)
library(qs)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

args <- commandArgs(trailingOnly = FALSE)
dataset_id <- args[grep("--dataset",args)+1]
file_list <- args[grep("--files",args)+1] %>% 
    strsplit(., split=' ') %>% 
    unlist() %>% 
    read_tsv(col_names = FALSE) %>% 
    bind_rows() %>% 
    pull(X1)

megadata <- tibble()
for(file in file_list){
    status <- file %>% strsplit(., split = "__") %>% unlist() %>% .[3]
    pipeline <- file %>% strsplit(., split = "__") %>% unlist() %>% .[4:5] %>% paste(., collapse = "+")
    bio_context <- file %>% strsplit(., split = "__") %>% unlist() %>% .[2]
    subset_id <- file %>% strsplit(., split = "__") %>% unlist() %>% .[1]
    data <- qread(file) %>%
    mutate(status = status,
           bio_context = bio_context,
           subset = subset_id,
           pipeline = pipeline,
           subset_id = basename(subset_id)) %>%
    separate(subset_id, into = c("main_dataset", "subset"), sep = "_") %>%
    separate(pipeline, into = c("norm", "diffexp"), sep = "\\+", remove = FALSE) %>%
    replace_na(list(subset = "-"))
    megadata <- rbind(megadata, data)
}

write_tsv(megadata, file = paste0(dataset_id, "__deresults.tsv"))
