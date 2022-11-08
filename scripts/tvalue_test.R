library(tidyverse)
library(ggplot2)
library(nlme)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


cell_lines <- list.files('./results/',pattern='_.*\\.tsv') %>% sub(pattern='_.*\\.tsv', replacement='') %>% unique()
merged_data <- tibble()

for(cell_line in cell_lines){
  file_list <- list.files('./results/', pattern=paste(cell_line, '_.*__stat__cons__decoupleroutput.tsv', sep='')) %>% paste('./results/', ., sep='')
  for(file in file_list) {
    file_data <- read_tsv(file) %>% separate(...1, into=c('cell_line_treatment','norm', 'diffexp', 'dcmethod'), sep='__', remove=FALSE) %>% separate(cell_line_treatment, into=c('cell_line', 'treatment'), sep='_') %>%
      pivot_longer(cols=c('Androgen':ncol(.)), names_to = 'pathway', values_to = 'scores') %>% select(-dcmethod)
    merged_data <- rbind(merged_data, file_data)
  }
}

aspc_test <- merged_data %>% filter(cell_line=='ASPC')
pathways <- aspc_test %>% distinct(pathway) %>% pull()
result_df <- tibble()
for(pathway in pathways){
  selected_data <- aspc_test %>% filter(pathway==!!pathway)
  result <- lme(scores ~ 0 + norm, data = selected_data, random = ~ 1 | treatment) %>% summary()
  result_pathway <- result$tTable %>% as.data.frame() %>% rownames_to_column(var='norm') %>% mutate(norm = sub('norm','',norm)) %>% as_tibble() %>% mutate(pathway=!!pathway)
  result_df <- rbind(result_df, result_pathway)
}



  

mergeddata <- right_join(tsquare, martin) %>% pivot_longer(-ID, names_to = 'method', values_to = 'scores')

ggplot(mergeddata, aes(scores)) +
  geom_histogram(alpha=0.5, bins=50) +
  facet_grid(~method)

ggplot(mergeddata, aes(scores, fill=method)) +
  geom_histogram(alpha=0.5, bins=100,position = 'identity')

new <- mergeddata %>% pivot_wider(names_from = method, values_from = scores) %>% column_to_rownames(var = 'ID') %>%as.matrix()

right_join(tsquare, martin) %>% summarize(cor=cor(stat_martin, stat_tsquare))
S