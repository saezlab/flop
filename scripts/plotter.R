library(tidyverse)
library(stats)
library(cowplot)
library(egg)
library(grid)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
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
    filter(analysis == "randindex", dataset == "GSE186341") %>%
    pull(path) %>%
    lapply(., randindex_reader) %>%
    bind_rows() %>%
    filter(k == "11" | k == "32") %>%
    dplyr::rename(value = scores)

plotter <- function(results_df, category) {
    if(grepl("Rand", category)){
        datasets <- "GSE186341"
    } else {
        datasets <- c("GTex", "CCLE", "GSE186341")
    }
    resources <- results_df %>%
        distinct(resource) %>%
        pull()
    for (dataset in datasets) {
        plots_list <- list()
        for (resource in resources) {
            toplot <- results_df %>%
                dplyr::filter(main_dataset == !!dataset, resource == !!resource) %>%
                mutate(id_2 = paste0(status, " - ", id))

            id_order <- toplot %>%
                group_by(id_2) %>%
                summarise(median_score = median(value, na.rm = TRUE)) %>%
                arrange(., median_score) %>%
                mutate(id_2 = as_factor(id_2))

            sorted_data <- toplot %>%
                arrange(., factor(id_2, levels = id_order$id_2)) %>%
                mutate(
                    id_2 = as_factor(id_2),
                    category = !!category
                )

            id_heatmap_data <- id_order %>%
                separate(id_2, into = c("is_filtered", "pipeline_a", "pipeline_b"), sep = " - ", remove = FALSE) %>%
                separate(pipeline_a, into = c("norm_method__a", "diffexp_method__a"), sep = "\\+") %>%
                separate(pipeline_b, into = c("norm_method__b", "diffexp_method__b"), sep = "\\+") %>%
                dplyr::select(-median_score) %>%
                pivot_longer(cols = -c(id_2)) %>%
                mutate(
                    category = case_when(
                        grepl("filtered", name) ~ "Filtering",
                        grepl("norm", name) ~ "Normalization",
                        grepl("diffexp", name) ~ "DE"
                    ),
                    method = value,
                    val = 1,
                    pipeline_num = case_when(
                        grepl("norm_method__a", name) ~ "Pipeline 1",
                        grepl("norm_method__b", name) ~ "Pipeline 2",
                        grepl("diffexp_method__a", name) ~ "Pipeline 1",
                        grepl("diffexp_method__b", name) ~ "Pipeline 2",
                        grepl("is_filtered", name) ~ value
                    )
                ) %>%
                mutate(category = fct_relevel(category, c("Filtering", "Normalization", "DE")))
            hm_plot <- id_heatmap_data %>%
                ggplot(aes(y = id_2, x = method, fill = as_factor(pipeline_num), alpha = 0.5)) +
                geom_tile() +
                guides(x = guide_axis(angle = 60)) +
                scale_fill_discrete() +
                facet_grid(cols = vars(category), scales = "free", space = "free",
                    labeller = as_labeller(c(
                        `Filtering` = "F",
                        `Normalization` = "N",
                        `DE` = "DE"))) +
                theme_cowplot() +
                theme(
                    legend.position = "none",
                    axis.title = element_blank()
                )

            x_min <- ifelse(grepl("Rank", category), -1, 0)

        
            box_plot <- ggplot(sorted_data) +
                annotate(geom = "rect", xmin = 0.8, xmax = 1, ymin = -Inf, ymax = Inf,
                    fill = "#439425", colour = NA, alpha = 0.2) +
                annotate(geom = "rect", xmin = 0.5, xmax = 0.8, ymin = -Inf, ymax = Inf,
                    fill = "#b7ba2e", colour = NA, alpha = 0.2) +
                annotate(geom = "rect", xmin = x_min, xmax = 0.5, ymin = -Inf, ymax = Inf,
                    fill = "#d04a35", colour = NA, alpha = 0.2) +
                {if (category == "Rank correlation") list(
                    geom_boxplot(aes(y = id_2, x = value)),
                    geom_violin(aes(y = id_2, x = value)))
                } +
                {if (category == "Rand Index") list(
                    geom_point(aes(y = id_2, x = value, shape = factor(k))),
                    labs(shape = "K"))
                } +
                {if (category == "Jaccard Index") list(
                    geom_boxplot(aes(y = id_2, x = value)),
                    geom_violin(aes(y = id_2, x = value)))
                } +
                xlim(c(x_min, 1)) +
                facet_grid(cols = vars(category), rows = vars(resource), scales = "free", space = "free",
                    labeller = as_labeller(c(
                        `dorothea` = "DoRothEA",
                        `msigdb_hallmarks` = "MSigDB Hallmarks",
                        `progeny` = "PROGENy",
                        `Rank correlation`= "Rank correlation"))) +
                theme_cowplot() +
                theme(
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.title = element_blank(),
                    axis.line.y = element_blank(), 
                )

            p1 <- cowplot::plot_grid(hm_plot, box_plot, nrow = 1, rel_widths = c(0.7, 1), align = "h")
            plots_list[[resource]] <- p1
            
        }

        aligned_plots <- egg::ggarrange(
            plots = plots_list,
            nrow = length(plots_list)
        )

        merged_plot <- ggdraw() + draw_plot(aligned_plots, 0, 0, 1, 1)

        merged_filename <- paste0(dataset, "__", category)

        save_plot(
            filename = paste0("../results/plots/", merged_filename, ".png"),
            plot = merged_plot,
            device = "png",
            dpi = 300,
            base_height = 18,
            base_width = 13
        )

        save_plot(
            filename = paste0("../results/plots/", merged_filename, ".svg"),
            plot = merged_plot,
            device = "svg",
            dpi = 300,
            base_height = 18,
            base_width = 13
        )
    }
}

rank_plotter <- function(results_df, category){
    plots_list <- list()
    datasets <- c("GTex", "CCLE", "GSE186341")
    for(dataset in datasets){
        toplot <- results_df %>%
            dplyr::filter(main_dataset == !!dataset) %>%
            mutate(id_2 = paste0(status, " - ", id),
                    resource = case_when(grepl("msigdb_hallmarks", name) ~ "MSigDB Hallmarks",
                                        grepl("dorothea", name) ~ "DoRothEA",
                                        grepl("progeny", name) ~ "PROGENy")) 

        id_order <- toplot %>%
            group_by(id_2) %>%
            summarise(median_score = median(value, na.rm = TRUE)) %>%
            arrange(.,median_score) %>%
            mutate(id_2 = as_factor(id_2))
        
        sorted_data <- toplot %>%
            arrange(., factor(id_2, levels = id_order$id_2)) %>%
            mutate(id_2 = as_factor(id_2),
            category = !!category)
        
        id_heatmap_data <- id_order %>%
            separate(id_2, into = c("is_filtered", "pipeline_a", "pipeline_b"), sep = " - ", remove = FALSE) %>%
            separate(pipeline_a, into = c("norm_method__a", "diffexp_method__a"), sep = "\\+") %>%
            separate(pipeline_b, into = c("norm_method__b", "diffexp_method__b"), sep = "\\+") %>%
            dplyr::select(-median_score) %>%
            pivot_longer(cols = -c(id_2)) %>%
            mutate(category = case_when(grepl("filtered", name) ~ "Filtering",
                                        grepl("norm", name) ~ "Normalization",
                                        grepl("diffexp", name) ~ "DE"),
                    method = value,
                    val = 1,
                    pipeline_num = case_when(grepl("norm_method__a", name) ~ "Pipeline 1",
                                            grepl("norm_method__b", name) ~ "Pipeline 2",
                                            grepl("diffexp_method__a", name) ~ "Pipeline 1",
                                            grepl("diffexp_method__b", name) ~ "Pipeline 2",
                                            grepl("is_filtered", name) ~ value
                                            )
            ) %>%
            mutate(category = fct_relevel(category, c("Filtering", "Normalization", "DE")))
        hm_plot <- id_heatmap_data %>%
            ggplot(aes(y = id_2, x = method, fill = as_factor(pipeline_num), alpha=0.5)) +
            geom_tile() +
            guides(x = guide_axis(angle = 60)) +
            scale_fill_discrete() +
            facet_grid(cols = vars(category), scales = "free", space = "free",
                    labeller = as_labeller(c(
                        `Filtering` = "F",
                        `Normalization` = "N",
                        `DE` = "DE"))) +
            theme_cowplot() +
            theme(legend.position = "none",
                axis.title = element_blank()
            )
        
        x_min <- ifelse(grepl("Rank",category), -1, 0)
        
        box_plot <- ggplot(sorted_data) + 
            annotate(geom = "rect", xmin = 0.8, xmax = 1, ymin = -Inf, ymax = Inf,
               fill = "#439425", colour = NA, alpha = 0.2) +
            annotate(geom = "rect", xmin = 0.5, xmax = 0.8, ymin = -Inf, ymax = Inf,
               fill = "#b7ba2e", colour = NA, alpha = 0.2) +
            annotate(geom = "rect", xmin = -1, xmax = 0.5, ymin = -Inf, ymax = Inf,
               fill = "#d04a35", colour = NA, alpha = 0.2) +
            geom_boxplot(aes(y = id_2, x = value)) +
            geom_violin(aes(y = id_2, x = value)) +
            xlim(c(x_min,1)) +
            facet_grid(cols = vars(category), scales = "free", space = "free",
            labeller = as_labeller(c(
                `dorothea` = "DoRothEA",
                `msigdb_hallmarks` = "MSigDB Hallmarks",
                `progeny` = "PROGENy",
                `Rank correlation` = "Rank correlation"))) +
            theme_cowplot() +
            theme(
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title = element_blank(),
                axis.line.y = element_blank()
                )

        p1 <- cowplot::plot_grid(hm_plot, box_plot, nrow = 1, rel_widths = c(0.7, 1), align = "h")
        plots_list[[dataset]] <- p1
    }
    aligned_plots <- egg::ggarrange(
            plots = plots_list,
            nrow = 3
            )

    merged_plot <- ggdraw() + draw_plot(aligned_plots, 0, 0, 1, 1) +
        draw_plot_label("A", 0, 1) +
        draw_plot_label("B", 0, 0.667) +
        draw_plot_label("C", 0, 0.333)

    merged_filename <- "rank_plot"

    save_plot(
        filename = paste0("../results/plots/", merged_filename, ".png"),
        plot = merged_plot,
        device = "png",
        dpi = 300,
        base_height = 20,
        base_width = 16
    )

    save_plot(
        filename = paste0("../results/plots/", merged_filename, ".svg"),
        plot = merged_plot,
        device = "svg",
        dpi = 300,
        base_height = 20,
        base_width = 16
    )
}

rank_plotter(results_rank, "Rank correlation")

plotter(results_rank, "Rank correlation")
plotter(results_randindex, "Rand Index")
plotter(results_jaccard, "Jaccard Index")
