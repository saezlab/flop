library(tidyverse)

dir.create('benchmark_data', showWarnings = FALSE)

main_dataset <- 'Simonson2023'

files_list <- list.files(path = './benchmark', pattern=paste0(main_dataset, '.*\\.csv$'), recursive=TRUE, full.names = TRUE)

old_metadata <- files_list %>% 
    grep(pattern='coldata', ., value = TRUE) %>% 
    read_csv() %>% 
    .[,-1]
counts <- files_list %>% grep(pattern='pbulk', ., value = TRUE) %>% read_csv() %>% column_to_rownames('...1') %>% t(.) %>% as.data.frame(.) %>% rownames_to_column('gene_symbol') %>% as_tibble()
# count number of columns -1
length(counts) - 1
old_metadata %>% 
    nrow()

cell_types <- old_metadata %>% 
    pull(cell_type) %>% 
    unique()

for(cell_type in cell_types){
    metadata <- old_metadata %>% 
        mutate(group = case_when(
            cell_type == !!cell_type ~ !!cell_type,
            cell_type != !!cell_type ~ paste0('No', !!cell_type)
        )) %>% 
        mutate(sample_id = paste0(sample_id, '_', cell_type)) %>%
        select(sample_id, group, disease, sex, age) %>%
        dplyr::rename('sample_ID' = 'sample_id')

    contrasts <- tibble(
        group1 = cell_type,
        group2 = paste0('No', cell_type)
    )

    # file output

    dataset_dir = paste0('benchmark_data/', main_dataset, '_', cell_type, '/')

    dir.create(dataset_dir, showWarnings = FALSE)

    counts %>% 
        write_tsv(paste0(dataset_dir, main_dataset, '_', cell_type, '__countdata.tsv'))

    metadata %>%
        write_tsv(paste0(dataset_dir, main_dataset, '_', cell_type, '__metadata.tsv'))

    contrasts %>%
        write_tsv(paste0(dataset_dir, main_dataset, '_', cell_type, '__contrast.tsv'))

    contrasts

}
