library(tidyverse)
path_file = './scripts/'
source(paste0(path_file, "rank_analysis_helper.R"))
source(paste0(path_file, "jaccard_analysis_helper.R"))

datasets <- list.files('./flop_results/diffexp/') %>% gsub('__deresults.tsv', '', .)

for(dataset in datasets){
    diffexp_data <- list.files('./flop_results/diffexp/', pattern = dataset, full.names = TRUE) %>% read_tsv(.) %>%
        pivot_longer(c(logFC, stat), names_to = 'statparam', values_to = 'scores') %>%
        mutate(resource = 'DE') %>% dplyr::rename('items' = 'ID')
    diffexp_rank <- corr_analysis(diffexp_data, "pipeline")
    diffexp_jaccard <- jaccard_analysis(diffexp_data)

    funcomics_rank <- list.files('./flop_results/funcomics/rank/', pattern = paste0(dataset, '__total__rank.tsv'), full.names = TRUE) %>% read_tsv(.)
    funcomics_jaccard <- list.files('./flop_results/funcomics/jaccard/', pattern = paste0(dataset, '__jaccard.tsv'), full.names = TRUE) %>% read_tsv(.)
    results_rank <- bind_rows(diffexp_rank, funcomics_rank)
    results_jaccard <- bind_rows(diffexp_jaccard, funcomics_jaccard)
    write_tsv(results_rank, file = paste0("./flop_results/funcomics/rank/", dataset, "__merged__rank.tsv"))
    write_tsv(results_jaccard, file = paste0("./flop_results/funcomics/jaccard/", dataset, "__merged__jaccard.tsv"))

}
diffexp_data <- read_tsv('./flop_results/diffexp/') %>% 
    pivot_longer(c(logFC, stat), names_to = 'statparam', values_to = 'scores') %>%
    mutate(resource = 'DE') %>% dplyr::rename('items' = 'ID')
diffexp_rank <- corr_analysis(diffexp_data, "pipeline")
diffexp_jaccard <- jaccard_analysis(diffexp_data)

funcomics_rank <- read_tsv('./flop_results/funcomics/rank/MetaHeart__total__rank.tsv')
funcomics_jaccard <- read_tsv('./flop_results/funcomics/jaccard/MetaHeart__jaccard.tsv')
results_rank <- bind_rows(diffexp_rank, funcomics_rank)
results_jaccard <- bind_rows(diffexp_jaccard, funcomics_jaccard)
write_tsv(total_results, file = "./flop_results/funcomics/rank/MetaHeart__merged__rank.tsv")
write_tsv(jaccard_total, file = "./flop_results/funcomics/jaccard/MetaHeart__merged__jaccard.tsv")
