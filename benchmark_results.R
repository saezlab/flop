library(tidyverse)

results <- read_tsv('flop_results/funcomics/fullmerged/Simonson2023__fullmerge.tsv')
pk_cell_links <- read_delim(file='pk_data_celltypelinks.csv', delim = ';') %>% select(PK, Data)

biomarkers <- results %>%
    filter(resource == 'biomarkers') %>%
    group_by(bio_context, pipeline, status) %>%
    arrange(desc(act), .by_group = TRUE)

# map names from the data column in pk_cell_links to the items column in biomarkers using the common column PK
biomarkers_data <- biomarkers %>%
    left_join(pk_cell_links, by = c('items' = 'PK'))

# create a new column with 1 if the name in Data is present in Biocontext
biomarkers_data_filtered <- biomarkers_data %>%
    mutate(present = case_when(
        str_detect(pattern = Data, string = bio_context) & (str_detect(pattern = Data, string = items)) ~ 1,

        TRUE ~ 0
    ))



