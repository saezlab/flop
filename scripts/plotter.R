library(tidyverse)
library(stats)
library(cowplot)
library(egg)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

merged_files <- list.files(path = "./results/full_merge/", pattern = "*__fullmerge.tsv", full.names = TRUE) %>% as_tibble() %>%
    mutate(dataset = sub("./results/full_merge/", "", sub("__fullmerge.tsv", "", value)))
datasets <- merged_files %>% select(dataset) %>% pull()

for(dataset_id in datasets){
    merged_data <- merged_files %>%
        filter(dataset == !!dataset_id) %>%
        select(value) %>%
        pull() %>%
        read_tsv()

    bio_context <- merged_data %>% distinct(bio_context) %>% pull()
    pipelines <- merged_data %>% distinct(pipeline) %>% pull()
    items <- merged_data %>% distinct(items) %>% pull()
    resources <- merged_data %>% distinct(resource) %>% pull()
    status <- merged_data %>% distinct(status) %>% pull()

    cor_results <- merged_data %>%
    group_by(statparam, resource, bio_context, status, main_dataset) %>%
    group_split() %>% 
    purrr::map(., function(x) {
        to_cor <- x %>%
        select(pipeline, scores, items) %>%
        pivot_wider(
            names_from = pipeline,
            values_from = scores,
            values_fn = {mean}
            ) %>%
        column_to_rownames(var = "items")

        cor_results <- cor(to_cor, method = "spearman") %>%
        as.data.frame() %>%
        rownames_to_column(var = "feature_1") %>%
        pivot_longer(-feature_1) %>%
        mutate(
            statparam = unique(x$statparam),
            bio_context = unique(x$bio_context),
            resource = unique(x$resource),
            status = unique(x$status),
            main_dataset = unique(x$main_dataset)

        )

        return(cor_results)

    }) %>%
    bind_rows()

    cor_filt_results <- cor_results %>%
        subset(feature_1 != name) %>%
        rowwise() %>%
        mutate(
            id = paste0(
            sort(c(feature_1, name))[1],
            " - ",
            sort(c(feature_1, name))[2]
            )
        ) %>%
        distinct(
            id,
            statparam,
            bio_context,
            resource,
            status,
            main_dataset,
            .keep_all = TRUE)

    total_plot <- ggplot(
        cor_filt_results,
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
        dataset_id,
        "merged__rank.png",
        sep = "__"
    )
    save_plot(
        filename = rankplot_filename,
        plot = total_plot, device = "png",
        dpi = 300,
        base_height = 10,
        base_width = 10
    )
}
