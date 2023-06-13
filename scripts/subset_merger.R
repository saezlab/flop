library(tidyverse)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

status <- c("filtered", "unfiltered")
args <- commandArgs(trailingOnly = FALSE)
dataset_id <- args[grep("--dataset", args) + 1]
file_list <- args[grep("--files", args) + 1] %>% strsplit(., split = " ") %>% unlist() %>% as_tibble()
files_info <- file_list %>%
    dplyr::rename(filename = value) %>%
    separate(
        sep = "__",
        col = filename,
        into = c("datasetID", "status", "ext"),
        remove = FALSE
    ) %>%
    select(- ext) %>%
    separate(
        sep = "_",
        col = datasetID,
        into = c("dataset", "subset")
    ) %>%
    replace_na(list(subset = "-"))

database_merged <- tibble()
for(status_i in status){
    selected_files <- files_info %>%
        filter(
            dataset == !!dataset_id,
            status == !!status_i
        )
    subsets <- selected_files %>%
        select(subset) %>%
        distinct() %>%
        pull()
    for(subset in subsets){
        single_file <- selected_files %>%
            filter(subset == !!subset)
        data <- single_file %>%
            select(filename) %>%
            pull() %>%
            read_tsv(.) %>%
            mutate(
                status = !!status_i,
                subset = !!subset,
                main_dataset = !!dataset_id
            )
        database_merged <- rbind(database_merged, data)
    }
}

fullmerge_name <- paste(dataset_id, "fullmerge.tsv", sep = "__")
write_tsv(database_merged, file = fullmerge_name)






