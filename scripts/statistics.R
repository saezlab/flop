library(tidyverse)
library(stats)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Data import
files_info <- tibble()
rank_info <- list.files(
    path = "../results/rank/",
    pattern = "*__rank.tsv",
    full.names = TRUE
) %>%
    as_tibble() %>%
    separate(
        col = value,
        into = c("dataset", "type", "ext"),
        sep = "__", remove = FALSE
    ) %>%
    select(-ext) %>%
    mutate(
        dataset = sub("../results/rank/", "", dataset),
        analysis = "rank"
    ) %>%
    dplyr::rename(path = value)

randindex_info <- list.files(
    path = "../results/rand_index/",
    pattern = "*__randindex.tsv",
    full.names = TRUE
) %>%
    as_tibble() %>%
    separate(
        col = value,
        into = c("dataset", "type", "ext"),
        sep = "__", remove = FALSE
    ) %>%
    select(-ext) %>%
    mutate(
        dataset = sub("../results/rand_index/", "", dataset),
        analysis = "randindex"
    ) %>%
    dplyr::rename(path = value)

jaccard_info <- list.files(
    path = "../results/jaccard/",
    pattern = "*__jaccard_index.tsv",
    full.names = TRUE
) %>%
    as_tibble() %>%
    separate(
        col = value,
        into = c("dataset", "type", "ext"),
        sep = "__", remove = FALSE
    ) %>%
    select(-ext) %>%
    mutate(
        dataset = sub("../results/jaccard/", "", dataset),
        analysis = "jaccard"
    ) %>%
    dplyr::rename(path = value)

files_info <- rbind(rank_info, randindex_info, jaccard_info)
datasets <- files_info %>%
    distinct(dataset) %>%
    pull()
types <- files_info %>%
    distinct(type) %>%
    pull()
analysis_types <- files_info %>%
    distinct(analysis) %>%
    pull()

results_rank <- files_info %>%
    filter(type == "total", analysis == "rank") %>%
    pull(path) %>%
    lapply(., read_tsv) %>%
    bind_rows() %>%
    filter(statparam == "stat")

results_jaccard <- files_info %>%
    filter(analysis == "jaccard") %>%
    pull(path) %>%
    lapply(., read_tsv) %>%
    bind_rows() %>%
    filter(statparam == "stat")

randindex_reader <- function(path) {
    main_dataset <- strsplit(path, split = "/")[[1]][4] %>% sub("__randindex.tsv", "", .)
    curated_data <- read_tsv(path) %>% mutate(main_dataset = !!main_dataset, id = paste(pipeline1, "-", pipeline2))
    return(curated_data)
}

results_randindex <- files_info %>%
    filter(analysis == "randindex", type == "randindex.tsv", dataset == "GSE186341") %>%
    pull(path) %>%
    lapply(., read_tsv) %>%
    bind_rows() %>%
    filter(k == "11" | k == "32") %>%
    dplyr::rename(value = scores)

# Rank statistical analysis
# Filtered-Unfiltered
wilcox_rank_results <- list()
for (dataset in datasets) {
    wilcox_rank_data <- results_rank %>%
        filter(main_dataset == !!dataset)
    wilcox_rank_results[[dataset]] <- wilcox.test(value ~ status, data = wilcox_data, alternative = "greater", paired = TRUE, p.adjust.methods = "BH")
}

wilcox_rand_results <- list()
resources <- results_randindex %>%
    distinct(resource) %>%
    pull()
k_vals <- c("22", "32", "2")
for (resource in resources) {
    for (k in k_vals) {
        wilcox_rand_data <- results_randindex %>%
            filter(resource == !!resource, k == !!k)
        wilcox_rand_results[[resource]][[k]] <- wilcox.test(value ~ status, data = wilcox_rand_data, alternative = "greater", paired = TRUE, p.adjust.methods = "BH")
    }
}


rank_corr_mean <- results_rank %>%
    filter(statparam == "stat") %>%
    group_by(main_dataset, id, status, resource) %>%
    summarise(mean = mean(value), sd = sd(value), above_0.8 = mean>0.8, above_0.5 = mean>0.5)

rank_corr_mean_perc <- rank_corr_mean %>%
    group_by(main_dataset, resource) %>%
    group_split %>%
    purrr::map(function(x){
        tibble(main_dataset = unique(x$main_dataset), resource = unique(x$resource), 
        perc_0.8 = nrow(subset(x, above_0.8 == TRUE))/nrow(x), 
        perc_0.5 = nrow(subset(x, above_0.5 == TRUE))/nrow(x))
    }) %>%
    bind_rows()

jaccard_mean <- results_jaccard %>%
    filter(statparam == "stat") %>%
    group_by(main_dataset, id, status, resource) %>%
    summarise(mean = mean(value), sd = sd(value), above_0.5 = mean>0.5, above_0.8 = mean>0.8)

jaccard_corr_mean_perc <- jaccard_mean %>%
    group_by(main_dataset, resource) %>%
    group_split %>%
    purrr::map(function(x){
        tibble(main_dataset = unique(x$main_dataset), resource = unique(x$resource), 
        perc_0.5 = nrow(subset(x, above_0.5 == TRUE))/nrow(x),
        perc_0.8 = nrow(subset(x, above_0.8 == TRUE))/nrow(x))
    }) %>%
    bind_rows()

rand_mean <- results_randindex %>%
    group_by(main_dataset, id, status, resource, k) %>%
    summarise(mean = mean(value), sd = sd(value), above_0.8 = mean>0.8, above_0.5 = mean>0.5)

rand_mean_perc <- rand_mean %>%
    group_by(main_dataset, resource, k) %>%
    group_split %>%
    purrr::map(function(x){
        tibble(main_dataset = unique(x$main_dataset), k = unique(x$k), resource = unique(x$resource),
        perc_0.8 = nrow(subset(x, above_0.8 == TRUE))/nrow(x),
        perc_0.5 = nrow(subset(x, above_0.5 == TRUE))/nrow(x))
    }) %>%
    bind_rows()

# Barplots:
# Rank
rank_corr_mean_perc %>%
    ggplot(aes(x = main_dataset, y = perc_0.8, fill = resource)) +
    geom_bar(stat = "identity", position = "dodge") +
    cowplot::theme_cowplot() +
    theme()
