library(tidyverse)
library(nlme)
library(fossil)

datafile <- "./results/fullmerged/GSE186341__fullmerge.tsv"
dataset_id <- "GSE186341"
path_file <- "./scripts/"
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

args <- commandArgs(trailingOnly = FALSE)
path_file <- args[grep("--file=", args)] %>%
  sub("rand_index_analysis.R", "", .) %>%
  sub("--file=", "", .)
dataset_id <- args[grep("--dataset",args)+1]
datafile <- args[grep("--file", args) + 1][2]
source(paste0(path_file, "rand_index_helper.R"))

merged_data <- read_tsv(datafile)
bio_context <- merged_data %>% distinct(bio_context) %>% pull()
pipelines <- merged_data %>% distinct(pipeline) %>% pull()
resources <- merged_data %>% distinct(resource) %>% pull()
status <- merged_data %>% distinct(status) %>% pull()
statparam <- "stat"

#K rand index variation
sep_bio_context <- bio_context %>% as_tibble() %>% separate(col = value, into =c('cell_line', 'treatment'))
cell_lines <- sep_bio_context %>% distinct(cell_line) %>% mutate(num_clust = seq(1,11))
treatments <- sep_bio_context %>% distinct(treatment) %>% mutate(num_clust = seq(1,32))
cell_line_clust <- tibble(bio_context) %>% separate(col = bio_context, into =c('cell_line', 'treatment'), remove = FALSE) %>% inner_join(.,cell_lines, by = "cell_line") %>% select(-treatment)
treatment_clust <- tibble(bio_context) %>% separate(col = bio_context, into =c('cell_line', 'treatment'), remove = FALSE) %>% inner_join(.,treatments, by = "treatment") %>% select(-cell_line)
k_values <- seq(1, 32, 1)
ground_truths <- c(cell_line_clust, treatment_clust)
rand_results_long <- tibble()


for (status_i in status) {
  for (resource in resources) {
    for (i in k_values) {
      rand_results <- matrix(
        nrow = length(pipelines),
        ncol = 2,
        0,
        dimnames = list(pipelines, c("cell_line", "treatment"))
      ) %>%
      as.data.frame()
      cluster_results <- clustering_k(merged_data, i, resource, status_i)
      for (pipeline1 in pipelines) {
        if (i == 11) {
          sorted_cell_line <- cluster_results[[status_i]][[resource]][[pipeline1]]$labels %>% left_join(cell_line_clust, by=c('label' = 'bio_context')) %>% pull(num_clust)
          rand_result <- adj.rand.index(
            cluster_results[[status_i]][[resource]][[pipeline1]]$labels$clust,
            sorted_cell_line
          )
          rand_results[pipeline1, "cell_line"] <- rand_result

        } else if (i == 32) {
          sorted_treatment <- cluster_results[[status_i]][[resource]][[pipeline1]]$labels %>% left_join(treatment_clust, by=c('label' = 'bio_context')) %>% pull(num_clust)
          rand_result <- adj.rand.index(
            cluster_results[[status_i]][[resource]][[pipeline1]]$labels$clust,
            sorted_treatment
          )
        rand_results[pipeline1, "treatment"] <- rand_result
        }
      }

      rand_results_long <- rand_results %>%
        rownames_to_column(var = "pipeline1") %>%
        pivot_longer(.,
          cols = -pipeline1,
          names_to = "ground_truth",
          values_to = "scores"
        ) %>%
        mutate(k = i, main_dataset = !!dataset_id) %>%
        filter(!(ground_truth == "cell_line" & k == 32)) %>%
        filter(!(ground_truth == "treatment" & k == 11)) %>%
        rowwise() %>%
        mutate(
          resource = !!resource,
          status = !!status_i
        ) %>%
        rbind(., rand_results_long)
    }
  }
}

write_tsv(rand_results_long, file = paste0("./results/rand_index/", dataset_id, "__bio__adj__randindex.tsv"))

