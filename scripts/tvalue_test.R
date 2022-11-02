library(tidyverse)
library(ggplot2)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

tsquare <- read_tsv('GSE151251__NA__edger__de__t2equalsF.tsv') %>% select (ID, stat) %>%
  dplyr::rename(stat_tsquare = stat)
martin <- read_tsv('GSE151251__NA__edger__de__martinsform.tsv') %>% select (ID, stat) %>%
  dplyr::rename(stat_martin = stat)
mergeddata <- left_join(tsquare, martin) %>% pivot_longer(-ID, names_to = 'method', values_to = 'scores') %>% filter(!is.infinite(scores)) %>% filter(!is.na(scores))

ggplot(mergeddata, aes(x=method, y=scores)) +
  geom_boxplot()

new <- mergeddata %>% pivot_wider(names_from = method, values_from = scores) %>% column_to_rownames(var = 'ID') %>%as.matrix() %>% 
cor(new)
