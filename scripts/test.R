library(tidyverse)

stat <- 'logFC'
dc_method <- 'cons'

input_files <- list.files('.\\results', pattern=stat)
merged_file <- data.frame()
for(input_file in input_files){
    filedata <- read.csv(paste('.\\results\\',input_file,sep=''), sep='\t', header = TRUE) %>%
        as_tibble()
    cell_treat_info <- strsplit(input_file, split='__')[[1]][1]

    merged_file <- rbind(merged_file, filedata)
}

regex <- paste('.*', dc_method, '.*', sep = '')
filtered_file <- merged_file %>% 
    filter(str_detect(X, !!regex))
write.table(filtered_file, paste(stat, dc_method, 'test.tsv', sep='__'),sep='\t', quote=FALSE, row.names=FALSE)