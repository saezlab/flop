library(tidyverse)
library(ggplot2)
library(nlme)
library(fossil)
library(cowplot)
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

args <- commandArgs(trailingOnly = FALSE)
path_file <- args[grep("--file=", args)] %>%
  sub("rand_index_analysis.R", "", .) %>%
  sub("--file=", "", .)
dataset_id <- args[grep("--dataset",args)+1]
status <- args[grep("--status", args) + 1]
datafile <- args[grep("--file", args) + 1][2]

#dataset_id <- "GSE186341"
#path_file <- ""
#status <- "filtered"
#datafile <- "GSE186341__result.tsv"
source(paste0(path_file, "dendro_helpers.R"))
#source("dendro_helpers.R")

merged_data <- read_tsv(datafile)

bio_context <- merged_data %>% distinct(bio_context) %>% pull()
pipelines <- merged_data %>% distinct(pipeline) %>% pull()
items <- merged_data %>% distinct(items) %>% pull()
resources <- merged_data %>% distinct(resource) %>% pull()
statparam <- "stat"

clustering_k <- function(merged_data, k, resource) {
  cluster_results <- list()
  for (pipeline in pipelines) {
    filt_data <- merged_data %>%
      filter(
        pipeline == !!pipeline,
        statparam == !!statparam,
        resource == !!resource
        ) %>%
      mutate(id = sub(pattern = "__.*", "", obs_id))
    numeric_data <- filt_data %>%
      pivot_wider(
        names_from = items,
        values_from = scores,
        id_cols = id
      ) %>%
      column_to_rownames("id")
    cluster <- numeric_data %>%
      dist() %>%
      hclust(., "ave")
    hcdata <- dendro_data_k(cluster, k)
    cluster_results[[resource]][[pipeline]] <- hcdata
  }
  return(cluster_results)
}

#K rand index variation
k_values <- seq(from = 1, to = min(length(bio_context), 32), by = 1)
print(k_values)
rand_results_long <- tibble()
for (resource in resources) {
  for (i in k_values) {
    rand_results <- matrix(
      nrow = 5,
      ncol = 5,
      0,
      dimnames = list(pipelines, pipelines)
    ) %>%
    as.data.frame()
    row.names(rand_results)
    cluster_results <- clustering_k(merged_data, i, resource)
    for (pipeline1 in pipelines) {
      for (pipeline2 in pipelines) {
        rand_result <- rand.index(
          cluster_results[[resource]][[pipeline1]]$segments$clust,
          cluster_results[[resource]][[pipeline2]]$segments$clust
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
      rowwise() %>%
      mutate(
        diff = paste(
          sort(c(pipeline1, pipeline2))[1],
          sort(c(pipeline1, pipeline2))[2], sep = " - "
        ),
        resource = resource
      ) %>%
      filter(pipeline1 != pipeline2) %>%
      distinct(diff, .keep_all = TRUE) %>%
      rbind(., rand_results_long)
  }
}

boxplot_k <- ggplot(
  rand_results_long,
  aes(x = factor(k), y = scores)
  ) +

  ggdist::stat_halfeye(
    width = 0.8,
    justification = -.2,
    point_colour = NA
  ) +

  geom_boxplot(
    width = .15,
    outlier.shape = NA
  ) +

  geom_point(
    position = ggpp::position_jitternudge(
      x = -0.2,
      width = 0.07,
      nudge.from = "jittered.x"
    ),
    size = 0.6
  ) +

  xlab("Number of K clusters") +
  ylab("Rand Index score") +
  scale_x_discrete(breaks = c(1, levels(factor(rand_results_long$k))[5:9==9])) +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(vjust = 0, hjust = 0.5),
    plot.margin = margin(10, 10, 10, 20),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(2, "lines")
  ) +

  facet_grid(
    rows = vars(resource),
    labeller = labeller(
      resource = c(
        dorothea = "DoRothEA",
        msigdb_hallmarks = "MSigDB Hallmarks",
        progeny = "PROGENy"
      )
    )
  )

boxplot_pipelines <- ggplot(
  rand_results_long,
  aes(x = diff, y = scores, fill = diff)
  ) +

  ggdist::stat_halfeye(
    adjust = .5,
    width = .5,
    .width = 0,
    justification = -.2,
    point_colour = NA
  ) +

  geom_boxplot(
    width = .15,
    outlier.shape = NA
  ) +

  geom_point(
    aes(color = factor(diff)),
    position = ggpp::position_jitternudge(
      x = -0.3,
      width = 0.1,
      nudge.from = "jittered.x"
    ),
    size = 0.6
  ) +

  ylab("Rand Index score") +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
    axis.title.x = element_blank(),
    plot.margin = margin(10, 10, 10, 20),
    axis.text.y.left = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.spacing = unit(2, "lines")
  ) +

  facet_grid(
    rows = vars(resource),
    labeller = labeller(
      resource = c(
        dorothea = "DoRothEA",
        msigdb_hallmarks = "MSigDB Hallmarks",
        progeny = "PROGENy"
      )
    )
  )

aligned_plots <- egg::ggarrange(
  boxplot_k,
  boxplot_pipelines,
  ncol = 2,
  widths = c(2, 1)
)

randind_plot <- ggdraw() + draw_plot(aligned_plots, 0, 0, 1, 1) +
  draw_plot_label("A", 0, 1) +
  draw_plot_label("B", 0.65, 1)

randind_filename <- paste(
  dataset_id,
  status,
  statparam,
  "rand_index_plots.png",
  sep = "__"
)

save_plot(
  filename = randind_filename,
  plot = randind_plot,
  device = "png",
  dpi = 300,
  base_height = 10,
  base_width = 10
)
