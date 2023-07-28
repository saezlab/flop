library(tidyverse)


unfiltered_de <- read_tsv('./data/flop_results/diffexp/GSE186341__deresults.tsv')
remaining_biocontexts <- unfiltered_de %>% 
    group_by(status, bio_context, pipeline) %>% 
    filter(padj<0.05) %>% 
    count() %>% 
    filter(n>=30) %>%
    ungroup() %>% 
    group_by(status, bio_context) %>%
    count() %>%
    filter(all(n==6)) %>%
    distinct(bio_context, status) %>% 
    separate(bio_context, into=c('cont1', 'cont2', sep = '_v_')) %>%
    select(-cont2) %>%
    separate(cont1, into=c('cell_line', 'treatment'))
