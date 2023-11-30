library(tidyverse)

dir.create('benchmark_data', showWarnings = FALSE)

main_datasets <- c('Simonson2023', 'Armute2023', 'Chaffin2022', 'Koenig2022', 'Kuppe2022', 'Reichart2022')

covariates <- list('Simonson2023' = c('heart_failure', 'sex', 'age'),
'Armute2023' = c('heart_failure'),
'Chaffin2023' = c('heart_failure', 'sex', 'age'),
'Koenig2022' = c('heart_failure'),
'Kuppe2022' = c('sex', 'patient_group'),
'Reichart' = c('sex', 'heart_failure'))

for(main_dataset in main_datasets){

    files_list <- list.files(path = './benchmark', pattern=paste0(main_dataset, '.*\\.csv$'), recursive=TRUE, full.names = TRUE)

    old_metadata <- files_list %>% 
        grep(pattern='coldata', ., value = TRUE) %>% 
        read_csv() %>% 
        .[,-1]
    counts <- files_list %>% grep(pattern='pbulk', ., value = TRUE) %>% read_csv() %>% {if('sample_id' %in% colnames(.)) column_to_rownames(., 'sample_id') else column_to_rownames(., '...1')} %>%
        t(.) %>% as.data.frame(.) %>% rownames_to_column('gene_symbol') %>% as_tibble()

    cell_types <- old_metadata %>% 
        pull(cell_type) %>% 
        unique()

    for(cell_type in cell_types){
        metadata <- old_metadata %>% 
            mutate(group = case_when(
                cell_type == !!cell_type ~ !!cell_type,
                cell_type != !!cell_type ~ paste0('No', !!cell_type)
            )) %>% 
            select(colname, group, !!covariates[[main_dataset]]) %>%
            dplyr::rename('sample_ID' = 'colname')

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
}



