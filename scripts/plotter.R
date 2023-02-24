library(tidyverse)
library(stats)
library(cowplot)
library(egg)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
files_info <- tibble()
rank_info <- list.files(
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
    mutate(
        dataset = sub("./results/rank/", "", dataset),
        analysis = "rank") %>%
    dplyr::rename(path = value)

randindex_info <- list.files(
        path = "./results/rand_index/",
        pattern = "*__randindex.tsv",
        full.names = TRUE
    ) %>%
    as_tibble() %>%
    separate(
        col = value,
        into = c("dataset", "type", "ext"),
        sep = "__", remove = FALSE
    ) %>%
    select(- ext) %>%
    mutate(
        dataset = sub("./results/rand_index/", "", dataset),
        analysis = "randindex") %>%
    dplyr::rename(path = value)

files_info <- rbind(rank_info, randindex_info)
datasets <- files_info %>% distinct(dataset) %>% pull()
types <- files_info %>% distinct(type) %>% pull()
analysis_types <- files_info %>% distinct(analysis) %>% pull()

# for (dataset in datasets) {
    for (analysis_type in analysis_types) {
        if (analysis_type == "rank") {
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
        } else if (analysis_type == "randindex") {
            datafile <- files_info %>%
                filter(dataset == !!dataset) %>%
                pull(path)

            rand_results_long <- read_tsv(datafile) %>% filter(status == "filtered")
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

        }
    }
# }

results_rank <- files_info %>%
    filter(type=='total', analysis=='rank') %>%
    pull(path) %>%
    lapply(., read_tsv) %>%
    bind_rows()

randindex_reader <- function(path) {
    main_dataset <- strsplit(path, split = "/")[[1]][4] %>% sub("__randindex.tsv", "", .)
    curated_data <- read_tsv(path) %>% mutate(main_dataset = !!main_dataset, id = paste(pipeline1, "-", pipeline2))
    return(curated_data)
}

results_randindex <- files_info %>%
    filter(analysis=='randindex') %>%
    pull(path) %>%
    lapply(., randindex_reader) %>%
    bind_rows() %>%
    dplyr::rename(value = scores)

plotter <- function(results_df, category){
    for(dataset in datasets){
        toplot <- results_df %>%
            dplyr::filter(main_dataset == !!dataset) %>%
            mutate(id_2 = paste0(status, " - ", id)) 

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
                                            )) %>%
            mutate(category = fct_relevel(category, c("Filtering", "Normalization", "DE")))
        hm_plot <- id_heatmap_data %>%
            ggplot(aes(y = id_2, x = method, fill = as_factor(pipeline_num), alpha=0.5)) +
            geom_tile() +
            guides(x = guide_axis(angle = 60)) +
            scale_fill_discrete() +
            facet_grid(cols = vars(category), scales = "free", space = "free") +
            theme_cowplot() +
            theme(legend.position = "none",
            axis.title = element_blank())
        
        x_min <- ifelse(grepl("Rank",category), -1, 0)

        g <- rasterGrob(matrix(adjustcolor(c(rep("#d04a35", 10),"#ba892e", "#899425", "#449105"), alpha=0.2), nrow=1),
            width=unit(1,"npc"), height = unit(1,"npc"), 
            interpolate = TRUE, just=c(0.5,0.5)) 
        
        box_plot <- ggplot(sorted_data) + 
            # annotate(geom = "rect", xmin = 0.9, xmax = 1, ymin = -Inf, ymax = Inf,
            #    fill = "#439425", colour = NA, alpha = 0.2) +
            # annotate(geom = "rect", xmin = 0.5, xmax = 0.9, ymin = -Inf, ymax = Inf,
            #    fill = "#b7ba2e", colour = NA, alpha = 0.2) +
            # annotate(geom = "rect", xmin = -1, xmax = 0.5, ymin = -Inf, ymax = Inf,
            #    fill = "#d04a35", colour = NA, alpha = 0.2)
            annotation_custom(g, xmin=x_min, xmax=1, ymin=-Inf, ymax=Inf) +
            geom_boxplot(aes(y = id_2, x = value)) +
            geom_violin(aes(y = id_2, x = value)) +
            xlim(c(x_min,1)) +
            facet_grid(cols = vars(category), scales = "free", space = "free") +
            theme_cowplot() +
            theme(
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title = element_blank(),
                axis.line.y = element_blank()
                )

        p1 <- cowplot::plot_grid(hm_plot, box_plot, nrow = 1, rel_widths = c(0.5, 1), align = "h")
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

    merged_filename <- paste0("./results/plots/"category, ".png")

    save_plot(
        filename = merged_filename,
        plot = merged_plot,
        device = "png",
        dpi = 300,
        base_height = 20,
        base_width = 18
    )
}


plotter(results_rank, "Rank correlation")
plotter(results_randindex, "Rand Index across K(1 to 32)")






