library(tidyverse)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

status <- c("filtered", "unfiltered")

dataset_id <- args[grep("--dataset", args) + 1]
file_directory <- args[grep("--files", args) + 1]

colnames(files_info) <- c("path", "filename", "dataset", "subset", "status")

for(status_i in status){
    file_info <- list.files(
        path = paste("./results/files", dataset, status_i, sep = "/"),
        pattern = "__result.tsv",
        full.names = TRUE
        ) %>%
        as_tibble() %>%
        mutate(
            filename = sub("./results/files/.*/.*/", "", value)
        ) %>%
        separate(
            sep = "__",
            col = filename,
            into = c("datasetID", "ext"),
            remove = FALSE
        ) %>%
        select(- ext) %>%
        separate(
            sep = "_",
            col = datasetID,
            into = c("dataset", "subset")
        ) %>%
        dplyr::rename(path = value) %>%
        mutate(status = !!status_i) %>%
        replace_na(list(subset = "-"))

    files_info <- rbind(files_info, file_info)
}


main_datasets <- files_info %>% select(dataset) %>% distinct() %>% pull()

for(main_dataset in main_datasets){
    database_merged <- tibble()
    for(status_i in status){
        selected_files <- files_info %>%
            filter(
                dataset == !!main_dataset,
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
                select(path) %>%
                pull() %>%
                read_tsv(.) %>%
                mutate(
                    status = !!status_i,
                    subset = !!subset,
                    main_dataset = !!main_dataset
                )
            database_merged <- rbind(database_merged, data)
        }
    }
    fullmerge_name <- paste(main_dataset, "fullmerge.tsv", sep = "__")
    write_tsv(database_merged, file = fullmerge_name)
}





