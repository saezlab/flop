library(tidyverse)
library(stats)
library(cowplot)
library(egg)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

files_info <- list.files(
        path = "./results/rank/",
        pattern = "*__rank.tsv",
        full.names = TRUE
    ) %>%
    as_tibble() %>%
    separate(
        col = value,
        into = c("dataset", "type", "ext"),
        sep = "__", remove = FALSE
    ) %>%
    select(- ext) %>%
    mutate(dataset = sub("./results/rank/", "", dataset)) %>%
    dplyr::rename(path = value)
datasets <- files_info %>% distinct(dataset) %>% pull()
types <- files_info %>% distinct(type) %>% pull()

for (dataset in datasets) {
    for (type in types) {
        if (type == 'total') {
            datafile <- files_info %>%
                filter(
                    dataset == !!dataset,
                    type == !!type
                ) %>%
                pull(path)
            rank_results <- read_tsv(datafile)
            total_plot <- ggplot(
                rank_results,
                aes(x = id, y = value, fill = status)
                ) +
                facet_grid(
                    cols = vars(statparam),
                    rows = vars(resource),
                    labeller = labeller(
                    resource = c(
                        dorothea = "DoRothEA",
                        msigdb_hallmarks = "MSigDB Hallmarks",
                        progeny = "PROGENy"
                    ),
                    statparam = c(logFC = "LogFC", stat = "T values")
                    )
                ) +
                guides(x = guide_axis(angle = 60)) +
                geom_violin(
                    width = 2,
                    position = position_dodge(0.5),
                    alpha = 0.5,
                    color = NA
                ) +
                geom_boxplot(
                    width = 0.3,
                    position = position_dodge(0.5),
                    outlier.alpha = 0.3,
                    outlier.size = 0.6
                ) +
                # geom_point(
                    # position = position_jitterdodge(
                        # jitter.width = 0.1,
                        # dodge.width = 0.5
                    # ),
                    # size = 0.6,
                    # alpha = 0.01,
                # ) +
                # geom_boxplot() +
                theme_cowplot() +
                ylab("Score") +
                theme(
                    legend.position = "top",
                    axis.title.x = element_blank(),
                    plot.margin = margin(10, 10, 10, 24)
                )
            rankplot_filename <- paste(
                dataset,
                "pipeline__rank.png",
                sep = "__"
            )
            save_plot(
                filename = rankplot_filename,
                plot = total_plot, device = "png",
                dpi = 300,
                base_height = 10,
                base_width = 10
            )
        } else if (type == "diffexp" | type == "norm") {
            datafile <- files_info %>%
                filter(
                    dataset == !!dataset,
                    type == !!type
                ) %>%
                pull(path)
            
            rank_split_results <- read_tsv(datafile)

            split_corplot <- ggplot(
                rank_split_results,
                aes(x = id, y = value, fill = status)
                ) +
                facet_grid(
                    cols = vars(statparam),
                    rows = vars(resource),
                    labeller = labeller(
                    resource = c(
                        dorothea = "DoRothEA",
                        msigdb_hallmarks = "MSigDB Hallmarks",
                        progeny = "PROGENy"
                    ),
                    statparam = c(logFC = "LogFC", stat = "T values")
                    )
                ) +
                guides(x = guide_axis(angle = 60)) +
                geom_violin(
                    width = 2,
                    position = position_dodge(0.5),
                    alpha = 0.5,
                    color = NA
                ) +
                geom_boxplot(
                    width = 0.3,
                    position = position_dodge(0.5),
                    outlier.alpha = 0.3,
                    outlier.size = 0.6
                ) +
                theme_cowplot() +
                theme(
                    legend.position = "top",
                    axis.title.x = element_blank(),
                    strip.background.y = element_blank(),
                    strip.text.y = element_blank()
                )
            rankplot_split_filename <- paste(
                dataset,
                type,
                "split__rank.png",
                sep = "__"
            )
            save_plot(
                filename = rankplot_split_filename,
                plot = split_corplot, device = "png",
                dpi = 300,
                base_height = 10,
                base_width = 10
            )
        }
    }
}

