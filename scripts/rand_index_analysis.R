library(tidyverse)
library(nlme)
library(fossil)

args <- commandArgs(trailingOnly = FALSE)
path_file <- args[grep("--file=", args)] %>%
  sub("rand_index_analysis.R", "", .) %>%
  sub("--file=", "", .)
dataset_id <- args[grep("--dataset",args)+1]
datafile <- args[grep("--file", args) + 1][2]
k_val <- args[grep("--k_val", args) + 1]
k_type <- args[grep("--k_type", args) + 1]
source(paste0(path_file, "rand_index_helper.R"))

merged_data <- read_tsv(datafile)
bio_context <- merged_data %>% distinct(bio_context) %>% pull()
pipelines <- merged_data %>% distinct(pipeline) %>% pull()
resources <- merged_data %>% distinct(resource) %>% pull()
status <- merged_data %>% distinct(status) %>% pull()
statparam <- "stat"


if (k_type == 'range'){
  k_values <- seq(from = 1, to = max(length(bio_context), as.numeric(k_val)), by = 1)
  print(k_values)
} else if (k_type == 'discrete'){
  k_values <- max(k_val, length(bio_context)) %>% str_split(., " ") %>% unlist() %>% as.numeric()
}


#K rand index variation
rand_results_long <- tibble()
for (status_i in status) {
  for (resource in resources) {
    for (i in k_values) {
      print(paste0("k = ", i, " resource = ", resource, " status = ", status_i))
      rand_results <- matrix(
        nrow = length(pipelines),
        ncol = length(pipelines),
        0,
        dimnames = list(pipelines, pipelines)
      ) %>%
      as.data.frame()
      cluster_results <- clustering_k(merged_data, i, resource, status_i)
      for (pipeline1 in pipelines) {
        for (pipeline2 in pipelines) {
          rand_result <- adj.rand.index(
            cluster_results[[status_i]][[resource]][[pipeline1]]$segments$clust,
            cluster_results[[status_i]][[resource]][[pipeline2]]$segments$clust
          )
          rand_results[pipeline1, pipeline2] <- rand_result
        }
      }

      rand_results_long <- rand_results %>%
        rownames_to_column(var = "pipeline1") %>%
        pivot_longer(.,
          cols = -pipeline1,
          names_to = "pipeline2",
          values_to = "scores"
        ) %>%
        mutate(k = i) %>%
        mutate(k = i, main_dataset = !!dataset_id) %>%
        rowwise() %>%
        mutate(
          id = paste(
            sort(c(pipeline1, pipeline2))[1],
            sort(c(pipeline1, pipeline2))[2], sep = " - "
          ),
          resource = !!resource,
          status = !!status_i
        ) %>%
        filter(pipeline1 != pipeline2) %>%
        distinct(id, .keep_all = TRUE) %>%
        rbind(., rand_results_long)
    }
  }
}

summary_sd <- rand_results_long %>%
  group_by(k, resource) %>%
  summarise(main_dataset = dataset_id, sd = sd(scores))

write_tsv(rand_results_long, file = paste0(dataset_id, "__randindex.tsv"))
write_tsv(summary_sd, file = paste0(dataset_id, "__summary__randindex.tsv"))
