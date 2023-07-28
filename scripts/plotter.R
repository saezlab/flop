library(tidyverse)
library(stats)
library(cowplot)
library(egg)
library(grid)
library(ComplexHeatmap)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
files_info <- tibble()
rank_info <- list.files(
    path = "./flop_results/funcomics/rank/",
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
        dataset = sub("./flop_results/funcomics/rank/", "", dataset),
        analysis = "rank"
    ) %>%
    dplyr::rename(path = value)

randindex_info <- list.files(
    path = "./flop_results/funcomics/rand_index/",
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
        dataset = sub("./flop_results/funcomics/rand_index/", "", dataset),
        analysis = "randindex"
    ) %>%
    dplyr::rename(path = value)

jaccard_info <- list.files(
    path = "./flop_results/funcomics/jaccard/",
    pattern = "*__jaccard.tsv",
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
        dataset = sub("./flop_results/funcomics/jaccard/", "", dataset),
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
    filter(statparam == "stat") #%>% 
    # mutate(status = case_when(
    #         grepl("\\bfiltered\\b", status) ~ "Filtered",
    #         grepl("\\bunfiltered\\b", status) ~ "Unfiltered"
    #     ))

results_rank_norm <- files_info %>%
    filter(type == "norm", analysis == "rank") %>%
    pull(path) %>%
    lapply(., read_tsv) %>%
    bind_rows() %>%
    filter(statparam == "stat")

results_rank_diffexp <- files_info %>%
    filter(type == "diffexp", analysis == "rank") %>%
    pull(path) %>%
    lapply(., read_tsv) %>%
    bind_rows() %>%
    filter(statparam == "stat")

results_jaccard <- files_info %>%
    filter(type == 'jaccard.tsv', analysis == "jaccard") %>%
    pull(path) %>% 
    lapply(., read_tsv) %>%
    bind_rows() %>%
    filter(statparam == "stat") %>%
    # mutate(status = case_when(
    #         grepl("\\bfiltered\\b", status) ~ "Filtered",
    #         grepl("\\bunfiltered\\b", status) ~ "Unfiltered"
    #     )) %>%
    select(-feature_1, -name) %>%
    separate(id, into = c('feature_1', 'name'), sep = " - ", remove = FALSE)

# Dotplots 

results_rank_points <- results_rank %>%
    filter(type == "correlation", statparam == 'stat') %>%
    group_by(id, resource) %>%
    summarise(feature_1 = feature_1, name = name, sd = sd(value), value = mean(value))

results_occurrence <- results_rank %>%
    filter(type == "occurrence", statparam == 'stat') %>%
    group_by(id, resource) %>%
    summarise(feature_1 = feature_1, name = name, sd = sd(value), value = mean(value))

occurrence_boxplots <- ggplot(results_occurrence, aes(x=id, y=value)) +
    geom_boxplot() +
    facet_grid(~resource) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

ocurrence_dotplot <- ggplot(results_occurrence, aes(x=feature_1, y=name)) +
    geom_point(aes(color = log(value), size = sd)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))+
    facet_grid(~resource)

ggsave("./flop_results/paper_plots/correlation_dotplots.png", ocurrence_dotplot, width = 15, height = 5)



correlation_dotplot <- ggplot(results_rank_points, aes(x=feature_1, y=name)) +
    geom_point(aes(color = value, size = sd)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))

ggsave("./flop_results/paper_plots/correlation_dotplots.png", correlation_dotplot, width = 15, height = 5)

library(ggdendro)

results_rank_points <- results_rank_points %>%
    filter(resource == 'dorothea')


# Figure 2
# A

resources <- c('DE', 'DoRothEA')
heatmaps_list <- c()
col_fun = circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
lgd = Legend(col_fun = col_fun, title = 'Spearman rank correlation', legend_height = unit(8, "cm"), title_position = "leftcenter-rot")
lgd_anno = Legend(labels = c('Filtered', 'Unfiltered'), title = "Status", legend_gp = gpar(fill = c('#47d3b4', '#de4ce3')))
for(resource in resources){
    correlation_matrix <- results_rank %>%
        filter(type == "correlation", statparam == 'stat', main_dataset %in% c('Spurell19', 'Sweet18', 'vanHeesch19', 'Yang14'), resource == ifelse(!!resource == 'DE', !!resource, tolower(!!resource))) %>%
        select(feature_1,name,value) %>%
        arrange(feature_1) %>%
        bind_rows(., tibble(
            feature_1 = .$name,
            name = .$feature_1,
            value = .$value
        )) %>%
        bind_rows(., tibble(
            feature_1 = unique(c(.$feature_1, .$name)),
            name = unique(c(.$feature_1, .$name)),
            value = 1)

        ) %>%
        pivot_wider(names_from=name, values_from=value, values_fn = mean) %>%
        column_to_rownames('feature_1') %>% select(rownames(.)) %>% as.matrix()

    pipelines_status <- rownames(correlation_matrix) %>% gsub('-.*', '', .)
    
    hc = HeatmapAnnotation(Filtering= pipelines_status,  col = list(Filtering = c("filtered" = "#47d3b4", "unfiltered" = "#de4ce3")), show_annotation_name = FALSE, show_legend = FALSE)
    hr = HeatmapAnnotation(Filtering= pipelines_status, col = list(Filtering = c("filtered" = "#47d3b4", "unfiltered" = "#de4ce3")), which = "row", show_annotation_name = FALSE, show_legend = FALSE)

    heatmap <- ComplexHeatmap::Heatmap(correlation_matrix, cluster_rows = TRUE, cluster_columns = TRUE, col = col_fun,
        column_title = resource, show_heatmap_legend = FALSE, width = unit(8, "cm"), height = unit(8, "cm"), bottom_annotation = hc, right_annotation = hr)

    heatmaps_list <- heatmaps_list + heatmap
}

heatmap_grob <- grid::grid.grabExpr(draw(heatmaps_list, ht_gap = unit(2, "cm")))
legend_grob <- grid::grid.grabExpr(draw(lgd))
legend_anno_grob <- grid::grid.grabExpr(draw(lgd_anno))


results_rank_correlation <- results_rank %>% filter(type == 'correlation', statparam == 'stat', main_dataset != 'GSE186341') %>% mutate(Space = ifelse(resource=='DE', 'Gene space', 'Functional space')) %>%
    mutate(resource = case_when(grepl('dorothea', resource) ~ 'DoRothEA',
                                grepl('msigdb_hallmarks', resource) ~ 'MSigDB Hallmarks',
                                grepl('progeny', resource) ~ 'PROGENy',
                                grepl('DE', resource) ~ 'DE'))

boxplot_rank <- ggplot(results_rank_correlation, aes(x = resource, y = value, fill = Space)) +
    geom_boxplot(outlier.alpha = 0) +
    ggpubr::stat_compare_means(aes(x=resource, y=value), comparisons = list(c('DE', 'DoRothEA'), c('DE', 'MSigDB Hallmarks'), c('DE', 'PROGENy')), size = 4, method = 'wilcox.test', method.args = list(alternative = 'greater', paired = TRUE, p.adjust.methods = "BH")) +
    ylab('Spearman rank correlation') +
    scale_fill_hue() +
    ggforce::geom_sina(size = 0.5, alpha = 0.2) +
    theme_cowplot() +
    ylim(0, 1.3) +
    theme(
        plot.background = element_rect(fill = "white",colour = NA),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 12),
        axis.title.x = element_blank(),
        text = element_text(size = 12))

results_jaccard_boxplot <- results_jaccard %>% filter(statparam == 'stat', main_dataset != 'GSE186341') %>% mutate(Space = ifelse(resource=='DE', 'Gene space', 'Functional space')) %>%
    mutate(resource = case_when(grepl('dorothea', resource) ~ 'DoRothEA',
                                grepl('msigdb_hallmarks', resource) ~ 'MSigDB Hallmarks',
                                grepl('progeny', resource) ~ 'PROGENy',
                                grepl('DE', resource) ~ 'DE'))


boxplot_jaccard <- ggplot(results_jaccard_boxplot, aes(x = resource, y = value, fill = Space)) +
    geom_boxplot(outlier.alpha = 0) +
    ggpubr::stat_compare_means(aes(x=resource, y=value), comparisons = list(c('DE', 'DoRothEA'), c('DE', 'MSigDB Hallmarks'), c('DE', 'PROGENy')), size = 4, method = 'wilcox.test', method.args = list(alternative = 'greater', paired = TRUE, p.adjust.methods = "BH")) +
    ylab('Jaccard index') +
    scale_fill_hue() +
    ggforce::geom_sina(size = 0.5, alpha = 0.2) +
    theme_cowplot() +
    ylim(0, 1.3) +
    theme(
        plot.background = element_rect(fill = "white",colour = NA),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 12),
        axis.title.x = element_blank(),
        text = element_text(size = 12), 
        axis.text = element_text(size = 12))

boxplots_merged <- egg::ggarrange(plots=list(boxplot_rank, boxplot_jaccard), nrow = 1, align = "h", padding = unit(10, "cm"))



figure_2 <- ggdraw() +
    draw_plot(heatmap_grob, 0.04, 0.17, 0.96, 1) +
    draw_plot(legend_grob, -0.44, 0.24, 1, 1) +
    draw_plot(legend_anno_grob, -0.44, 0.05, 1, 1) +
    draw_plot(boxplot_rank, 0.08, 0, 0.32, 0.39) +
    draw_plot(boxplot_jaccard, 0.44, 0, 0.445, 0.39) +
    draw_label("A.", size = 20, x = 0.05, y = 0.95) +
    draw_label("C.", size = 20, x = 0.05, y = 0.4) +
    draw_label("B.", size = 20, x = 0.45, y = 0.95) +
    draw_label("D.", size = 20, x = 0.45, y = 0.4)


ggsave(plot = figure_2, filename = "./flop_results/paper_plots/figure_2.png", width = 30, height = 30, units = 'cm', dpi = 300)




# Figure 3
run_info <- read_tsv('./flop_results/traces/trace-20230719-36429954__ccle.txt') %>% mutate(vmem = as.numeric(gsub(' GB', '', vmem))) %>%
    rowwise() %>%
    filter(!process %in% c('jaccard_analysis', 'rank_analysis', 'output_merge_de', 'subset_merger', 'downstream_merge_de', 'decoupler_merger')) %>%
    mutate(realtime = case_when(
        length(unlist(strsplit(realtime, ' '))) == 1 ~ paste0("0h 0m ", realtime),
        length(unlist(strsplit(realtime, ' '))) == 2 ~ paste0("0h ", realtime),
        length(unlist(strsplit(realtime, ' '))) == 3 ~ realtime)) %>%
    ungroup() %>%
    mutate(realtime = hms(realtime)) %>% 
    mutate(realtime = hour(realtime)*60 + minute(realtime) + second(realtime)/60) %>%
    filter(!process %in% c('jaccard_analysis', 'rank_analysis', 'output_merge_de', 'subset_merger', 'downstream_merge_de', 'decoupler_merger'))


memory_plot <- ggplot(run_info, aes(x = process, y = vmem)) +
    geom_boxplot() + 
    theme_cowplot() + 
    ylab('Memory usage (GB)') +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    axis.title.x = element_blank()) 

time_plot <- ggplot(run_info, aes(x = process, y = realtime)) +
    geom_boxplot() + 
    theme_cowplot() + 
    ylab('Time (minutes)') +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    axis.title.x = element_blank()) 

boxplots_fig3 <- egg::ggarrange(plots=list(memory_plot, time_plot), nrow = 2, align = 'v')

figure_3 <- ggdraw() +
    draw_plot(boxplots_fig3, 0.05, 0.05, 0.9, 0.9) +
    draw_label('A.', size = 20, x = 0.05, y = 0.95) +
    draw_label('B.', size = 20, x = 0.05, y = 0.5)

ggsave(plot = figure_3, filename = "./flop_results/paper_plots/figure_3.png", width = 30, height = 25, units = 'cm', dpi = 300)






#Figure 4

results_rank_figure4 <- results_rank %>% filter(statparam == 'stat', main_dataset == 'GSE186341', type == 'correlation', resource == 'dorothea') %>% mutate(Space = ifelse(resource=='DE', 'Gene space', 'Functional space')) %>%
    mutate(resource = case_when(grepl('dorothea', resource) ~ 'DoRothEA')
    ) %>%
    arrange(feature_1) %>%
    group_by(statparam, bio_context, resource, status, main_dataset) %>%
    group_split() %>%
    purrr::map(., function(x){
        pipelines <- unique(c(x$feature_1, x$name))
        corr_1_df <- x %>% select(statparam, bio_context, resource, status, main_dataset) %>% distinct() %>%
            bind_cols(., tibble(
                feature_1 = pipelines,
                name = pipelines,
                value = 1))
        inverted_corrs <- x %>% rename(feature_1 = name, name = feature_1)
        updated_df <- bind_rows(corr_1_df, inverted_corrs, x)
        return(updated_df)
    }) %>% bind_rows() %>%
    select(feature_1, name, value) %>%
    pivot_wider(names_from=name, values_from=value, values_fn = mean) %>% column_to_rownames('feature_1') %>% select(rownames(.)) %>% as.matrix()

results_jaccard_figure4 <- results_jaccard %>% filter(statparam == 'stat', main_dataset == 'GSE186341', resource == 'dorothea') %>% mutate(Space = ifelse(resource=='DE', 'Gene space', 'Functional space')) %>%
    mutate(resource = case_when(grepl('dorothea', resource) ~ 'DoRothEA')
    ) %>% 
    arrange(feature_1) %>%
    group_by(statparam, bio_context, resource, status, main_dataset) %>%
    group_split() %>%
    purrr::map(., function(x){
        pipelines <- unique(c(x$feature_1, x$name))
        corr_1_df <- x %>% select(statparam, bio_context, resource, main_dataset) %>% distinct() %>%
            bind_cols(., tibble(
                feature_1 = pipelines,
                name = pipelines,
                value = 1))
        inverted_corrs <- x %>% rename(feature_1 = name, name = feature_1)
        updated_df <- bind_rows(corr_1_df, inverted_corrs, x)
        return(updated_df)
    }) %>% bind_rows() %>%
    select(feature_1, name, value) %>%
    group_by(feature_1, name) %>%
    pivot_wider(names_from=name, values_from=value, values_fn = mean) %>% column_to_rownames('feature_1') %>% select(rownames(.)) %>% as.matrix()



resource <- c('DoRothEA')
heatmaps_list <- c()
col_fun = circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
lgd = Legend(col_fun = col_fun, title = NULL, legend_height = unit(8, "cm"))

pipelines_status <- rownames(results_rank_figure4) %>% gsub('-.*', '', .)

hc = HeatmapAnnotation(Filtering= pipelines_status,  col = list(Filtering = c("filtered" = "#47d3b4", "unfiltered" = "#de4ce3")), show_annotation_name = FALSE, show_legend = FALSE)
hr = HeatmapAnnotation(Filtering= pipelines_status, col = list(Filtering = c("filtered" = "#47d3b4", "unfiltered" = "#de4ce3")), which = "row", show_annotation_name = FALSE, show_legend = FALSE)

heatmap <- ComplexHeatmap::Heatmap(results_rank_figure4, cluster_rows = TRUE, cluster_columns = TRUE,
    column_title = paste0("Rank Correlation, ", resource), show_heatmap_legend = FALSE, width = unit(8, "cm"), height = unit(8, "cm"), bottom_annotation = hc, right_annotation = hr, col = col_fun)
heatmap_jaccard <- ComplexHeatmap::Heatmap(results_jaccard_figure4, cluster_rows = TRUE, cluster_columns = TRUE,
    column_title = paste0("Jaccard Index, ", resource), show_heatmap_legend = FALSE, width = unit(8, "cm"), height = unit(8, "cm"), bottom_annotation = hc, right_annotation = hr, col = col_fun)

heatmaps_list <- heatmap + heatmap_jaccard


heatmap_grob <- grid::grid.grabExpr(draw(heatmaps_list, ht_gap = unit(2, "cm")))
legend_grob <- grid::grid.grabExpr(draw(lgd))

figure_4 <- ggdraw() +
    draw_plot(heatmap_grob) +
    draw_plot(legend_grob, 0.05, 0.39, 0.1, 0.5) +
    draw_plot(legend_anno_grob, 0.05, 0.05, 0.1, 0.5) +
    draw_label("A.", size = 20, x = 0.1, y = 0.98) +
    draw_label("B.", size = 20, x = 0.45, y = 0.98)

ggsave(plot = figure_4, filename = "./flop_results/paper_plots/figure_4.png", width = 35, height = 16, units = 'cm', dpi = 300)
# cor_dendro <- as.dendrogram(hclust(d = dist(x = correlation_matrix)))
# dendro_plot <- ggdendrogram(data = cor_dendro, rotate = TRUE)
# dendro_order <- order.dendrogram(cor_dendro)

# results_rank_points$feature_1 <- factor(x = results_rank_points$feature_1,
#                                levels = rownames(correlation_matrix)[dendro_order], 
#                                ordered = TRUE)
# results_rank_points$name <- factor(x = results_rank_points$name,
#                                levels = rownames(correlation_matrix)[dendro_order], 
#                                ordered = TRUE)


# heatmap_plot <- ggplot(data = results_rank_points, aes(x = feature_1, y = name)) +
#   geom_tile(aes(fill = value)) +
#   scale_fill_gradient2(midpoint = 0.5) +
#   theme(axis.text.y = element_text(size = 6))

heatmap <- ComplexHeatmap::Heatmap(correlation_matrix, cluster_rows = TRUE, cluster_columns = TRUE)


# results_randindex <- files_info %>%
#     filter(analysis == "randindex", dataset == "GSE186341", type == "randindex.tsv") %>%
#     pull(path) %>%
#     lapply(., read_tsv) %>%
#     bind_rows() %>%
#     mutate(status = case_when(
#             grepl("\\bfiltered\\b", status) ~ "Filtered",
#             grepl("\\bunfiltered\\b", status) ~ "Unfiltered"
#         )) %>%
#     dplyr::rename(value = scores)

# results_randindex_bio <- files_info %>%
#     filter(analysis == "randindex", dataset == "GSE186341", type == "bio") %>%
#     pull(path) %>%
#     lapply(., read_tsv) %>%
#     bind_rows() %>%
#     filter(k == "10" | k == "22") %>%
#     dplyr::rename(value = scores)

# results_randindex_summary <- files_info %>%
#     filter(analysis == "randindex", dataset == "GSE186341", type == "summary") %>%
#     pull(path) %>%
#     lapply(., read_tsv) %>%
#     bind_rows()

jaccard_matrix <- results_jaccard %>%
    mutate(feature_1 = paste(status, feature_1, sep='+'),
        name = paste(status, name, sep='+')) %>%
    filter(resource == 'DE') %>%
    select(feature_1,name,value) %>%
    pivot_wider(names_from=name, values_from=value, values_fn = mean) %>%
    column_to_rownames('feature_1') %>% as.matrix()

Heatmap(jaccard_matrix, cluster_rows = FALSE, cluster_columns = FALSE)





plotter <- function(results_df, category) {
    if(grepl("Rand", category)){
        datasets <- "GSE186341"
        results_df <- results_df %>% filter(k == 10 | k == 22)
    } else {
        datasets <- results_df %>% distinct(main_dataset) %>% pull()
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
                summarise(mean_score = mean(value, na.rm = TRUE)) %>%
                arrange(., mean_score) %>%
                mutate(id_2 = as_factor(id_2))

            sorted_data <- toplot %>%
                arrange(., factor(id_2, levels = id_order$id_2)) %>%
                mutate(
                    id_2 = as_factor(id_2),
                    category = !!category
                )

            id_heatmap_data <- id_order %>%
                separate(id_2, into = c("is_filtered", "pipeline_a", "pipeline_b"), sep = " - ", remove = FALSE) %>%
                # separate(id_2, into = c("is_filtered", "pipelines"), sep = " - ", remove = FALSE) %>%
                # separate(pipelines, into = c("pipeline_a", "pipeline_b"), sep = "-", remove = TRUE) %>%
                separate(pipeline_a, into = c("norm_method__a", "diffexp_method__a"), sep = "\\+") %>%
                separate(pipeline_b, into = c("norm_method__b", "diffexp_method__b"), sep = "\\+") %>%
                dplyr::select(-mean_score) %>%
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
                scale_fill_manual(values = c(
                    "Filtered" = "purple",
                    "Unfiltered" = "orange",
                    "Pipeline 1" = "#439425",
                    "Pipeline 2" = "darkred"
                )) +
                guides(x = guide_axis(angle = 60)) +
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
                    geom_boxplot(aes(y = id_2, x = value)))
                } +
                xlim(c(x_min, 1)) +
                facet_grid(cols = vars(category), rows = vars(resource), scales = "free", space = "free",
                    labeller = as_labeller(c(
                        `dorothea` = "DoRothEA",
                        `msigdb_hallmarks` = "MSigDB Hallmarks",
                        `progeny` = "PROGENy",
                        `Rank correlation`= "Rank correlation",
                        `Rand Index`= "Rand Index",
                        `DE` = "DE",
                        `Jaccard Index`= "Jaccard Index"))) +
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
            filename = paste0("./flop_results/paper_plots/", merged_filename, ".png"),
            plot = merged_plot,
            device = "png",
            dpi = 300,
            base_height = 18,
            base_width = 13
        )

        save_plot(
            filename = paste0("./flop_results/paper_plots/", merged_filename, ".svg"),
            plot = merged_plot,
            device = "svg",
            dpi = 300,
            base_height = 18,
            base_width = 13
        )
    }
}
sub_jac <- results_jaccard %>% filter(main_dataset == 'GSE186341', statparam == 'stat')
plotter(sub_jac, 'Jaccard Index')
merged_plotter <- function(results_df, category) {
    datasets <- results_df %>% distinct(main_dataset) %>% pull()
    resources <- results_df %>%
        distinct(resource) %>%
        pull()

    plots_list <- list()
    for (resource in resources) {
        toplot <- results_df %>%
            dplyr::filter(resource == !!resource) %>%
            mutate(id_2 = paste0(status, " - ", id))

        id_order <- toplot %>%
            group_by(id_2) %>%
            summarise(mean_score = mean(value, na.rm = TRUE)) %>%
            arrange(., mean_score) %>%
            mutate(id_2 = as_factor(id_2))

        sorted_data <- toplot %>%
            arrange(., factor(id_2, levels = id_order$id_2)) %>%
            mutate(
                id_2 = as_factor(id_2),
                category = !!category
            )

        id_heatmap_data <- id_order %>%
            separate(id_2, into = c("is_filtered", "pipeline_a", "pipeline_b"), sep = " - ", remove = FALSE) %>%
            # separate(id_2, into = c("is_filtered", "pipelines"), sep = " - ", remove = FALSE) %>%
            # separate(pipelines, into = c("pipeline_a", "pipeline_b"), sep = "-", remove = TRUE) %>%
            separate(pipeline_a, into = c("norm_method__a", "diffexp_method__a"), sep = "\\+") %>%
            separate(pipeline_b, into = c("norm_method__b", "diffexp_method__b"), sep = "\\+") %>%
            dplyr::select(-mean_score) %>%
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
            scale_fill_manual(values = c(
                "Filtered" = "purple",
                "Unfiltered" = "orange",
                "Pipeline 1" = "#439425",
                "Pipeline 2" = "darkred"
            )) +
            guides(x = guide_axis(angle = 60)) +
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
                geom_boxplot(aes(y = id_2, x = value)))
            } +
            xlim(c(x_min, 1)) +
            facet_grid(cols = vars(category), rows = vars(resource), scales = "free", space = "free",
                labeller = as_labeller(c(
                    `dorothea` = "DoRothEA",
                    `msigdb_hallmarks` = "MSigDB Hallmarks",
                    `progeny` = "PROGENy",
                    `Rank correlation`= "Rank correlation",
                    `Rand Index`= "Rand Index",
                    `DE` = "DE",
                    `random` = 'Random DoRothEA',
                    `Jaccard Index`= "Jaccard Index",
                    `DE progeny` = "DE PROGENy",
                    `DE msigdb_hallmarks` = "DE MSigDB Hallmarks",
                    `DE dorothea` = "DE DoRothEA",
                    `DE TFs` = 'DE TFs'))) +
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
        ncol = ifelse(length(plots_list) > 3, 2, 1)
        
        nrow = ceiling(length(plots_list)/ncol)

        aligned_plots <- egg::ggarrange(
            plots = plots_list,
            nrow = nrow,
            ncol = ncol
        )

        merged_plot <- ggdraw() + draw_plot(aligned_plots, 0, 0, 1, 1)

        merged_filename <- paste0('merged', "__", category)

        save_plot(
            filename = paste0("./flop_results/paper_plots/", merged_filename, ".png"),
            plot = merged_plot,
            device = "png",
            dpi = 300,
            base_height = 6*nrow,
            base_width = 13*ncol
        )

        save_plot(
            filename = paste0("./flop_results/paper_plots/", merged_filename, ".svg"),
            plot = merged_plot,
            device = "svg",
            dpi = 300,
            base_height = 6*nrow,
            base_width = 13*ncol
        )
    
}

merged_plotter(results_rank, "Rank correlation")
merged_plotter(results_jaccard, "Jaccard Index")




plotter(results_rank, "Rank correlation")
# plotter(results_randindex, "Rand Index")
plotter(results_jaccard, "Jaccard Index")

randindex_plotter <- function(results_df, category) {
    datasets <- "GSE186341"
    results_df <- results_df %>% dplyr::filter(k == 10 | k == 22)
    resources <- results_df %>%
        distinct(resource) %>%
        pull()
    for (dataset in datasets) {
        plots_list <- list()
        for (resource in resources) {
            toplot <- results_df %>%
                dplyr::filter(main_dataset == !!dataset, resource == !!resource) %>%
                mutate(id_2 = paste0(status, " - ", pipeline1))

            id_order <- toplot %>%
                group_by(id_2) %>%
                summarise(mean_score = mean(value, na.rm = TRUE)) %>%
                arrange(., mean_score) %>%
                mutate(id_2 = as_factor(id_2))

            sorted_data <- toplot %>%
                arrange(., factor(id_2, levels = id_order$id_2)) %>%
                mutate(
                    id_2 = as_factor(id_2),
                    category = !!category,
                    ground_truth = case_when(
                        grepl("cell_line", ground_truth) ~ "k=10 vs Cell line",
                        grepl("treatment", ground_truth) ~ "k=22 vs Treatment"
                    )
                )

            id_heatmap_data <- id_order %>%
                separate(id_2, into = c("is_filtered", "pipeline_a"), sep = " - ", remove = FALSE) %>%
                separate(pipeline_a, into = c("norm_method__a", "diffexp_method__a"), sep = "\\+") %>%
                dplyr::select(-mean_score) %>%
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
                        grepl("norm_method__a", name) ~ "Pipeline",
                        grepl("diffexp_method__a", name) ~ "Pipeline",
                        grepl("is_filtered", name) ~ value
                    )
                ) %>%
                mutate(category = fct_relevel(category, c("Filtering", "Normalization", "DE")))
            hm_plot <- id_heatmap_data %>%
                ggplot(aes(y = id_2, x = method, fill = as_factor(pipeline_num), alpha = 0.5)) +
                geom_tile() +
                scale_fill_manual(values = c(
                    "filtered" = "purple",
                    "unfiltered" = "orange",
                    "Pipeline" = "#439425"
                )) +
                guides(x = guide_axis(angle = 60)) +
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
                geom_point(aes(y = id_2, x = value, shape = factor(ground_truth))) +
                labs(shape = "Hierarchical vs custom") +
                xlim(c(x_min, 1)) +
                facet_grid(cols = vars(category), rows = vars(resource), scales = "free", space = "free",
                    labeller = as_labeller(c(
                        `dorothea` = "DoRothEA",
                        `msigdb_hallmarks` = "MSigDB Hallmarks",
                        `progeny` = "PROGENy",
                        `Rank correlation`= "Rank correlation",
                        `Rand Index`= "Rand Index",
                        `Jaccard Index`= "Jaccard Index"))) +
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

        merged_filename <- paste0(dataset, "__bio__", category)

        save_plot(
            filename = paste0("./flop_results/paper_plots/", merged_filename, ".png"),
            plot = merged_plot,
            device = "png",
            dpi = 300,
            base_height = 13,
            base_width = 13
        )

        save_plot(
            filename = paste0("./flop_results/paper_plots/", merged_filename, ".svg"),
            plot = merged_plot,
            device = "svg",
            dpi = 300,
            base_height = 13,
            base_width = 13
        )
    }
}

randindex_plotter(results_randindex_bio, "Rand Index")



# get the highest sd and the k it belongs to, per resource

summary_plotter <- function(results_df, summary) {
    plots_list <- list()
    dataset <- "GSE186341"

    resources <- results_df %>%
        distinct(resource) %>%
        pull()

    category <- 'Rand Index'

    for (resource in resources) {
        summary <- results_randindex_summary %>% 
            group_by(resource) %>%
            summarise(
                max_sd = max(sd),
                k = filter(., sd == max_sd) %>% pull(k)
            )

        k_value <- summary %>%
            filter(resource == !!resource) %>%
            pull(k)
        results_subset <- results_df %>%
            filter(resource == !!resource & k == !!k_value)
        toplot <- results_subset %>%
            mutate(id_2 = paste0(status, " - ", id)) 

        id_order <- toplot %>%
            group_by(id_2) %>%
            summarise(mean_score = mean(value, na.rm = TRUE)) %>%
            arrange(.,mean_score) %>%
            mutate(id_2 = as_factor(id_2))
        
        sorted_data <- toplot %>%
            arrange(., factor(id_2, levels = id_order$id_2)) %>%
            mutate(id_2 = as_factor(id_2),
            category = !!category)
        
        id_heatmap_data <- id_order %>%
            separate(id_2, into = c("is_filtered", "pipelines"), sep = " - ", remove = FALSE) %>%
            separate(pipelines, into = c("pipeline_a", "pipeline_b"), sep = "-", remove = TRUE) %>%
            separate(pipeline_a, into = c("norm_method__a", "diffexp_method__a"), sep = "\\+") %>%
            separate(pipeline_b, into = c("norm_method__b", "diffexp_method__b"), sep = "\\+") %>%
            dplyr::select(-mean_score) %>%
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
            scale_fill_manual(values = c(
                "Filtered" = "purple",
                "Unfiltered" = "orange",
                "Pipeline 1" = "#439425",
                "Pipeline 2" = "darkred"
            )) +
            guides(x = guide_axis(angle = 60)) +
            facet_grid(cols = vars(category), scales = "free", space = "free",
                    labeller = as_labeller(c(
                        `Filtering` = "F",
                        `Normalization` = "N",
                        `DE` = "DE"))) +
            theme_cowplot() +
            theme(legend.position = "none",
                axis.title = element_blank()
            )
        
        x_min <- 0
        box_plot <- ggplot(sorted_data) +
            annotate(geom = "rect", xmin = 0.8, xmax = 1, ymin = -Inf, ymax = Inf,
                fill = "#439425", colour = NA, alpha = 0.2) +
            annotate(geom = "rect", xmin = 0.5, xmax = 0.8, ymin = -Inf, ymax = Inf,
                fill = "#b7ba2e", colour = NA, alpha = 0.2) +
            annotate(geom = "rect", xmin = x_min, xmax = 0.5, ymin = -Inf, ymax = Inf,
                fill = "#d04a35", colour = NA, alpha = 0.2) +
            geom_point(aes(y = id_2, x = value))+
            xlim(c(x_min, 1)) +
            facet_grid(cols = vars(category), rows = vars(resource), scales = "free", space = "free",
                labeller = as_labeller(c(
                    `dorothea` = paste0("DoRothEA, k = ", k_value),
                    `msigdb_hallmarks` = paste0("MSigDB Hallmarks, k = ", k_value),
                    `progeny` = paste0("PROGENy, k = ", k_value),
                    `Rank correlation`= "Rank correlation",
                    `Rand Index`= "Rand Index",
                    `Jaccard Index`= "Jaccard Index"))) +
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
            nrow = 3
            )

    merged_plot <- ggdraw() + draw_plot(aligned_plots, 0, 0, 1, 1) +
        draw_plot_label("A", 0, 1) +
        draw_plot_label("B", 0, 0.667) +
        draw_plot_label("C", 0, 0.333)

    merged_filename <- paste0(dataset, "__maxsd__", category)

    save_plot(
        filename = paste0("./flop_results/paper_plots/", merged_filename, ".png"),
        plot = merged_plot,
        device = "png",
        dpi = 300,
        base_height = 20,
        base_width = 16
    )

    save_plot(
        filename = paste0("./flop_results/paper_plots/", merged_filename, ".svg"),
        plot = merged_plot,
        device = "svg",
        dpi = 300,
        base_height = 20,
        base_width = 16
    )
}


summary_plotter(results_randindex, results_randindex_summary)

summarised_results_rank <- results_rank %>%
    group_by(feature_1, name, status) %>%
    summarise(meanval = mean(value), stdev = sd(value)) %>%
    ungroup() %>%
    # mutate(value = ifelse(status == "filtered", value, value*-1),
    mutate(comparison = paste0(feature_1, " - ", name),
    ranking = rank(-meanval)) %>%
    mutate(angle = 90 - 360 * (.$ranking-0.5)/nrow(.),
    hjust = ifelse(angle < -90, 1, 0)) %>%
    mutate(angle = ifelse(angle < -90, angle+180, angle))


sorted_results_rank <- results_rank %>%
    mutate(feature_1 = case_when(
            grepl("deseq2", feature_1) ~ "DESeq2",
            grepl("edger", feature_1) ~ "edgeR",
            grepl("tmm", feature_1) ~ "TMM+l",
            grepl("voom", feature_1) ~ "voom+l",
            grepl("vsn", feature_1) ~ "vsn+l",
            grepl("log2quant", feature_1) ~ "logQ+l"),
            name = case_when(
            grepl("deseq2", name) ~ "DESeq2",
            grepl("edger", name) ~ "edgeR",
            grepl("tmm", name) ~ "TMM+l",
            grepl("voom", name) ~ "voom+l",
            grepl("vsn", name) ~ "vsn+l",
            grepl("log2quant", name) ~ "logQ+l")
        ) %>%
    left_join(., (summarised_results_rank %>% dplyr::select(comparison, ranking, status)), by = c("id" = "comparison", "status" = "status")) %>%
    mutate(ranking = ranking)

circular_plot_data <- results_rank %>%
    group_by(feature_1, name, status) %>%
    summarise(meanval = mean(value), stdev = sd(value)) %>%
    ungroup() %>%
    mutate(feature_1 = case_when(
            grepl("deseq2", feature_1) ~ "DESeq2",
            grepl("edger", feature_1) ~ "edgeR",
            grepl("tmm", feature_1) ~ "TMM+l",
            grepl("voom", feature_1) ~ "voom+l",
            grepl("vsn", feature_1) ~ "vsn+l",
            grepl("log2quant", feature_1) ~ "logQ+l"),
            name = case_when(
            grepl("deseq2", name) ~ "DESeq2",
            grepl("edger", name) ~ "edgeR",
            grepl("tmm", name) ~ "TMM+l",
            grepl("voom", name) ~ "voom+l",
            grepl("vsn", name) ~ "vsn+l",
            grepl("log2quant", name) ~ "logQ+l")
        ) %>%
    mutate(comparison = paste0(feature_1, " - ", name),
        ranking = rank(-meanval)) %>%
    rbind(tibble(feature_1 = rep("", ceiling(nrow(.))), 
        name = rep("", ceiling(nrow(.))), 
        status = rep("", ceiling(nrow(.))), 
        meanval = rep(NA, ceiling(nrow(.))), 
        stdev = rep(NA, ceiling(nrow(.))),
        comparison = rep("", ceiling(nrow(.))),
        ranking = seq(nrow(.)+1, nrow(.) + ceiling(nrow(.)))
    )) %>%
    mutate(angle = - 360 * (.$ranking-0.5)/nrow(.),
    hjust = ifelse(angle < -90, 0, 1)) %>%
    mutate(angle = ifelse(angle < -90, angle+180, angle))


# summarised_results_rank_norm <- results_rank_norm %>%
#     group_by(feature_1, name, status) %>%
#     summarise(meanval = mean(value), stdev = sd(value)) %>%
#     ungroup() %>%
#     # mutate(value = ifelse(status == "filtered", value, value*-1),
#     mutate(comparison = paste0(feature_1, " - ", name),
#     ranking = rank(-meanval)) %>%
#     mutate(angle = 90 - 360 * (.$ranking-0.5)/nrow(.),
#     hjust = ifelse(angle < -90, 1, 0)) %>%
#     mutate(angle = ifelse(angle < -90, angle+180, angle))

# summarised_results_rank_diffexp <- results_rank_diffexp %>%
#     group_by(feature_1, name, status) %>%
#     summarise(meanval = mean(value), stdev = sd(value)) %>%
#     ungroup() %>%
#     # mutate(value = ifelse(status == "filtered", value, value*-1),
#     mutate(comparison = paste0(feature_1, " - ", name),
#     ranking = rank(-meanval)) %>%
#     mutate(angle = 90 - 360 * (.$ranking-0.5)/nrow(.),
#     hjust = ifelse(angle < -90, 1, 0)) %>%
#     mutate(angle = ifelse(angle < -90, angle+180, angle))
    

rank_circular_barplot <- ggplot(circular_plot_data) +
    geom_segment(x = 1, xend = 30.5, y = 1, yend = 1, color = "grey", linewidth = 0.5) +
    geom_col(aes(x = as.factor(ranking), y = meanval, fill = status),width = 0.5) +
    geom_boxplot(data = sorted_results_rank, aes(x = as.factor(ranking), y = value), outlier.alpha = 0, inherit.aes = FALSE) +
    coord_polar(clip = 'off', start = -pi/2) +
    theme_minimal() +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(9,4), "cm"),
        legend.position = c(0.55,0.55),
        legend.title = element_blank(),
        text = element_text(size = 25)
    ) +
    scale_y_continuous(limits = c(-2, 1), oob = scales::rescale_none) +
    scale_fill_manual(values =  c(
                "Filtered" = "#5ab4ac",
                "Unfiltered" = "#d8b365"
            )) +
    # geom_errorbar(aes(x = ranking, ymin=meanval-stdev, ymax=meanval+stdev), position = position_dodge())  +
    annotate(geom = "rect", ymin = 1.08, ymax = 1.35, xmin = -Inf, xmax = nrow(circular_plot_data[circular_plot_data$meanval>0.8 & !is.na(circular_plot_data$meanval),]) + 0.5,
                fill = "#439425", colour = NA, alpha = 0.2) +
    annotate(geom = "rect", ymin = 1.08, ymax = 1.35, xmin = nrow(circular_plot_data[circular_plot_data$meanval>0.8 & !is.na(circular_plot_data$meanval),]) + 0.5, xmax = nrow(circular_plot_data[circular_plot_data$meanval>0.5 & !is.na(circular_plot_data$meanval),]) + 0.5,
                fill = "#b7ba2e", colour = NA, alpha = 0.2) +
    annotate(geom = "rect", ymin = 1.08, ymax = 1.35, xmin = nrow(circular_plot_data[circular_plot_data$meanval>0.5 & !is.na(circular_plot_data$meanval),]) + 0.5, xmax = nrow(circular_plot_data[!is.na(circular_plot_data$meanval),]) + 0.5,
                fill = "#d04a35", colour = NA, alpha = 0.2) +
    geom_text(aes(x = ranking, y = 1.45, label = comparison, hjust = hjust, angle = angle), color="black", size=7) +
    geom_text(aes(x = ranking, y = 1.12, label = round(meanval, digits = 2), hjust = hjust, angle = angle), color="black", size=6)


filtering_boxplot <- ggplot(results_rank) +
    geom_boxplot(aes(x = status, y = value, fill = status)) +
    theme_cowplot() +
    scale_fill_manual(values = c(
                "Filtered" = "#5ab4ac",
                "Unfiltered" = "#d8b365"
            )) +
    scale_y_continuous(limits = c(-1,1)) +
    theme(axis.title = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        legend.position = "none",
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()
        # plot.margin = unit(c(15,3,15,3), "cm")
        )

# heatmap_plot_diffexp <- ggplot(summarised_results_rank_diffexp) +
#     geom_tile(aes(x = feature_1, y = name, fill = meanval)) +
#     theme_cowplot() +
#     theme(
#         axis.text.x = element_text(angle = 90, hjust = 1),
#         axis.title = element_blank()
#     ) +
#     facet_grid(cols = vars(status), scales = "free", space = "free",
#         labeller = as_labeller(c(
#             `filtered` = "F",
#             `unfiltered` = "U"))) +
#     scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = mean(summarised_results_rank_diffexp$meanval)) +
#     guides(fill = guide_colorbar(barwidth = 10, barheight = 1, title.position = "top", title.hjust = 0.5)) +
#     labs(fill = "Rank correlation")

# filtering_plots <- egg::ggarrange(circular_barplot, filtering_boxplot, widths = c(2, 1))

# filtering_plots <- cowplot::plot_grid(circular_barplot, filtering_boxplot, rel_widths = c(1, 0.4))

filtering_plots <- ggdraw() + 
    draw_plot(circular_barplot, 0, 0, 1, 1) +
    draw_plot(filtering_boxplot, 0.1, 0.55, 0.3, 0.38) 


# save_plot(
#     filename = paste0("./flop_results/plots/", "filtering_plots", ".png"),
#     plot = filtering_plots,
#     device = "png",
#     dpi = 300,
#     base_height = 14,
#     base_width = 14
# )

# figure 2: rank

summarised_results_rank <- results_rank %>%
    group_by(feature_1, name, status) %>%
    summarise(meanval = mean(value), stdev = sd(value)) %>%
    ungroup() %>%
    # mutate(value = ifelse(status == "filtered", value, value*-1),
    mutate(comparison = paste0(feature_1, " - ", name),
    ranking = rank(-meanval)) %>%
    mutate(angle = 90 - 360 * (.$ranking-0.5)/nrow(.),
    hjust = ifelse(angle < -90, 1, 0)) %>%
    mutate(angle = ifelse(angle < -90, angle+180, angle))


sorted_results_rank <- results_rank %>%
    mutate(feature_1 = case_when(
            grepl("deseq2", feature_1) ~ "DESeq2",
            grepl("edger", feature_1) ~ "edgeR",
            grepl("tmm", feature_1) ~ "TMM+l",
            grepl("voom", feature_1) ~ "voom+l",
            grepl("vsn", feature_1) ~ "vsn+l",
            grepl("log2quant", feature_1) ~ "logQ+l"),
            name = case_when(
            grepl("deseq2", name) ~ "DESeq2",
            grepl("edger", name) ~ "edgeR",
            grepl("tmm", name) ~ "TMM+l",
            grepl("voom", name) ~ "voom+l",
            grepl("vsn", name) ~ "vsn+l",
            grepl("log2quant", name) ~ "logQ+l")
        ) %>%
    left_join(., (summarised_results_rank %>% dplyr::select(comparison, ranking, status)), by = c("id" = "comparison", "status" = "status")) %>%
    mutate(ranking = ranking)

circular_plot_data <- results_rank %>%
    group_by(feature_1, name, status) %>%
    summarise(meanval = mean(value), stdev = sd(value)) %>%
    ungroup() %>%
    mutate(feature_1 = case_when(
            grepl("deseq2", feature_1) ~ "DESeq2",
            grepl("edger", feature_1) ~ "edgeR",
            grepl("tmm", feature_1) ~ "TMM+l",
            grepl("voom", feature_1) ~ "voom+l",
            grepl("vsn", feature_1) ~ "vsn+l",
            grepl("log2quant", feature_1) ~ "logQ+l"),
            name = case_when(
            grepl("deseq2", name) ~ "DESeq2",
            grepl("edger", name) ~ "edgeR",
            grepl("tmm", name) ~ "TMM+l",
            grepl("voom", name) ~ "voom+l",
            grepl("vsn", name) ~ "vsn+l",
            grepl("log2quant", name) ~ "logQ+l")
        ) %>%
    mutate(comparison = paste0(feature_1, " - ", name),
        ranking = rank(-meanval)) %>%
    rbind(tibble(feature_1 = rep("", ceiling(nrow(.))), 
        name = rep("", ceiling(nrow(.))), 
        status = rep("", ceiling(nrow(.))), 
        meanval = rep(NA, ceiling(nrow(.))), 
        stdev = rep(NA, ceiling(nrow(.))),
        comparison = rep("", ceiling(nrow(.))),
        ranking = seq(nrow(.)+1, nrow(.) + ceiling(nrow(.)))
    )) %>%
    mutate(angle = - 360 * (.$ranking-0.5)/nrow(.),
    hjust = ifelse(angle < -90, 0, 1)) %>%
    mutate(angle = ifelse(angle < -90, angle+180, angle))
    

rank_circular_barplot <- ggplot(circular_plot_data) +
    geom_segment(x = 1, xend = 30.5, y = 1, yend = 1, color = "grey", linewidth = 0.5) +
    geom_col(aes(x = as.factor(ranking), y = meanval, fill = status),width = 0.5) +
    geom_boxplot(data = sorted_results_rank, aes(x = as.factor(ranking), y = value), outlier.alpha = 0, inherit.aes = FALSE) +
    coord_polar(clip = 'off', start = -pi/2) +
    theme_minimal() +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(9,4), "cm"),
        legend.position = c(0.55,0.55),
        legend.title = element_blank(),
        text = element_text(size = 25)
    ) +
    scale_y_continuous(limits = c(-2, 1), oob = scales::rescale_none) +
    scale_fill_manual(values =  c(
                "Filtered" = "#5ab4ac",
                "Unfiltered" = "#d8b365"
            )) +
    # geom_errorbar(aes(x = ranking, ymin=meanval-stdev, ymax=meanval+stdev), position = position_dodge())  +
    annotate(geom = "rect", ymin = 1.08, ymax = 1.5, xmin = -Inf, xmax = nrow(circular_plot_data[circular_plot_data$meanval>0.8 & !is.na(circular_plot_data$meanval),]) + 0.5,
                fill = "#439425", colour = NA, alpha = 0.2) +
    annotate(geom = "rect", ymin = 1.08, ymax = 1.5, xmin = nrow(circular_plot_data[circular_plot_data$meanval>0.8 & !is.na(circular_plot_data$meanval),]) + 0.5, xmax = nrow(circular_plot_data[circular_plot_data$meanval>0.5 & !is.na(circular_plot_data$meanval),]) + 0.5,
                fill = "#b7ba2e", colour = NA, alpha = 0.2) +
    annotate(geom = "rect", ymin = 1.08, ymax = 1.5, xmin = nrow(circular_plot_data[circular_plot_data$meanval>0.5 & !is.na(circular_plot_data$meanval),]) + 0.5, xmax = nrow(circular_plot_data[!is.na(circular_plot_data$meanval),]) + 0.5,
                fill = "#d04a35", colour = NA, alpha = 0.2) +
    geom_text(aes(x = ranking, y = 1.7, label = comparison, hjust = hjust, angle = angle), color="black", size=7) +
    geom_text(aes(x = ranking, y = 1.12, label = round(meanval, digits = 2), hjust = hjust, angle = angle), color="black", size=6)



rank_filtering_boxplot <- ggplot(results_rank, aes(x = status, y = value)) +
    ggforce::geom_sina(aes(color = status), jitter_y = FALSE, size = 2, alpha = 0.5, pch = 21) +
    # ggpubr::ggboxplot(results_rank, x = 'status', y = 'value', fill = 'status') +
    ggpubr::stat_compare_means(aes(x='status', y='value'), comparisons = list(c('Filtered', 'Unfiltered')), size = 7, method = 'wilcox.test', method.args = list(alternative = 'greater', paired = TRUE, p.adjust.methods = "BH")) +
    theme_cowplot() +
    scale_color_manual(values = c(
                "Filtered" = "#5ab4ac",
                "Unfiltered" = "#d8b365"
            )) +
    scale_y_continuous(limits = c(-1,1.1)) +
    theme(axis.title = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 20),
        legend.position = 'none',
        text = element_text(size = 20),
        plot.margin = unit(c(1,1,1,1), "cm"))

rank_dataset_sina_plot <- results_rank %>%
    group_by(id, main_dataset, status, resource) %>%
    summarise(meanval = mean(value), stdev = sd(value)) %>%
    ggplot(., aes(x = main_dataset, y = meanval)) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey", linewidth = 1) +
    ggforce::geom_sina(aes(color = status, group = main_dataset), jitter_y = FALSE, size = 2) +
    # ggpubr::stat_compare_means(aes(x='main_dataset', y='value'), comparisons = list(c('CCLE', 'GTex'), c('CCLE', 'GSE186341'), c('GTex', 'GSE186341')), size = 7) +
    scale_colour_manual(values = c(
                "Filtered" = "#5ab4ac",
                "Unfiltered" = "#d8b365"
            )) +
    theme_cowplot() +
    # scale_y_continuous(limits = c(0,1)) +
    theme(axis.title = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 20, angle = 60, hjust = 1),
        legend.position = 'none',
        plot.margin = unit(c(1,1,1,1), "cm"))

viewings_rank <- egg::ggarrange(rank_dataset_sina_plot, rank_filtering_boxplot, ncol=2, widths = c(1,1))

figure_2 <- ggdraw() +
    draw_plot(viewings_rank, 0, 0, 2/3, 0.45) +
    draw_plot(rank_circular_barplot, -0.1, -0.1, 1.2, 1.2) +
    draw_label('A.', x = 0.03, y = 0.99, size = 25) +
    draw_label('B.', x = 0.03, y = 0.44, size = 25) +
    draw_label('C.', x = 0.34, y = 0.44, size = 25)

save_plot(
    filename = paste0("./flop_results/paper_plots/", "figure_2", ".png"),
    plot = figure_2,
    device = "png",
    dpi = 300,
    base_height = 15,
    base_width = 18
)



# plot2
jaccard_dataset_sina_plot <- results_jaccard %>%
    group_by(id, main_dataset, status, resource) %>%
    summarise(meanval = mean(value), stdev = sd(value)) %>%
    ggplot(., aes(x = main_dataset, y = meanval)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey", linewidth = 1) +
    ggforce::geom_sina(aes(color = status, group = main_dataset), jitter_y = FALSE, size = 2) +
    scale_colour_manual(values = c(
                "Filtered" = "#5ab4ac",
                "Unfiltered" = "#d8b365"
            )) +
    theme_cowplot() +
    scale_y_continuous(limits = c(0,1)) +
    theme(axis.title = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 20),
        legend.position = 'none',
        plot.margin = unit(c(1,1,1,1), "cm")
        )

jaccard_filtering_boxplot <- ggplot(results_jaccard, aes(x = status, y = value)) +
    ggforce::geom_sina(aes(color = status), jitter_y = FALSE, size = 2, alpha = 0.5, pch = 21) +
    ggpubr::stat_compare_means(aes(x='status', y='value'), comparisons = list(c('Filtered', 'Unfiltered')), size = 7, method = 'wilcox.test', method.args = list(alternative = 'greater', paired = TRUE, p.adjust.methods = "BH")) +
    theme_cowplot() +
    scale_color_manual(values = c(
                "Filtered" = "#5ab4ac",
                "Unfiltered" = "#d8b365"
            )) +
    # scale_y_continuous(limits = c(-1,1)) +
    theme(axis.title = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 20),
        legend.position = 'none',
        plot.margin = unit(c(1,1,1,1), "cm")
        )


plots_list <- c()

library(ComplexUpset)
upset_jaccard <- results_jaccard %>%
    group_by(id, main_dataset, resource) %>%
    mutate(resource = case_when(
            grepl("dorothea", resource) ~ "DoRothEA",
            grepl("msigdb_hallmarks", resource) ~ "MSigDB",
            grepl("progeny", resource) ~ "PROGENy",
            TRUE ~ resource)
        ) %>%
    summarise(meanval = mean(value), stdev = sd(value)) %>%
    group_by(main_dataset, resource) %>%
    slice_max(order_by = meanval, n = 5) %>%
    mutate(value = TRUE) %>%
    distinct(id, main_dataset, resource, value) %>%
    pivot_wider(names_from  = resource,
              values_from = value,
              values_fill = FALSE)


resources <- upset_jaccard %>% ungroup() %>% select(-id, -main_dataset) %>% colnames()
plots_list <- list()
for(dataset in datasets){
    upset_jaccard_test <- upset_jaccard %>%
        filter(main_dataset==!!dataset) %>%
        ungroup() %>%
        select(-main_dataset) %>%
        column_to_rownames('id')
    plots_list[[dataset]] <- upset(upset_jaccard_test, intersect = resources,
        intersections = list(c("DoRothEA", "MSigDB"),
                                    c("DoRothEA", "MSigDB", "PROGENy"),
                                    c("MSigDB", "PROGENy"),
                                    c("DoRothEA", "PROGENy")),
        base_annotations=list('Intersection size'=intersection_size(counts=TRUE, text = list(size=7, nudge_y = 0.7, color = 'black')) + ylim(0, 5) + ylab(dataset)),
        name = "resource",
        set_sizes=FALSE,
        sort_intersections=FALSE,
        height_ratio = 1,
        themes=upset_modify_themes(
            list(
                'Intersection size'=theme(text = element_text(size=20),
                axis.title.y = element_text(size = 20)
            )))
        )[[1]]

}
plots_list[['inters']] <- upset(upset_jaccard_test, intersect = resources,
        intersections = list(c("DoRothEA", "MSigDB"),
                                    c("MSigDB", "PROGENy"),
                                    c("DoRothEA", "PROGENy"),
                                    c("DoRothEA", "MSigDB", "PROGENy")),
        name = "resource",
        set_sizes=FALSE,
        sort_intersections=FALSE,
        themes=upset_modify_themes(
            list(
                'intersections_matrix'=theme(axis.text = element_text(size=20),
                axis.title.x = element_blank()
            )))
        )[[2]]



upset_plots <- egg::ggarrange(plots = plots_list, nrow = 8, ncol = 1)

summarised_results_jaccard <- results_jaccard %>%
    group_by(feature_1, name, status) %>%
    mutate(feature_1 = case_when(
            grepl("deseq2", feature_1) ~ "DESeq2",
            grepl("edger", feature_1) ~ "edgeR",
            grepl("tmm", feature_1) ~ "TMM+l",
            grepl("voom", feature_1) ~ "voom+l",
            grepl("vsn", feature_1) ~ "vsn+l",
            grepl("log2quant", feature_1) ~ "logQ+l"),
            name = case_when(
            grepl("deseq2", name) ~ "DESeq2",
            grepl("edger", name) ~ "edgeR",
            grepl("tmm", name) ~ "TMM+l",
            grepl("voom", name) ~ "voom+l",
            grepl("vsn", name) ~ "vsn+l",
            grepl("log2quant", name) ~ "logQ+l")
        ) %>%
    summarise(meanval = mean(value), stdev = sd(value)) %>%
    ungroup() %>%
    # mutate(value = ifelse(status == "filtered", value, value*-1),
    mutate(comparison = paste0(feature_1, " - ", name)) %>%
    mutate(ranking = rank(meanval)) %>%
    mutate(angle = 90 - 360 * (.$ranking-0.5)/nrow(.),
    hjust = ifelse(angle < -90, 1, 0)) %>%
    mutate(angle = ifelse(angle < -90, angle+180, angle))


sorted_results_jaccard <- results_jaccard %>%
    mutate(feature_1 = case_when(
            grepl("deseq2", feature_1) ~ "DESeq2",
            grepl("edger", feature_1) ~ "edgeR",
            grepl("tmm", feature_1) ~ "TMM+l",
            grepl("voom", feature_1) ~ "voom+l",
            grepl("vsn", feature_1) ~ "vsn+l",
            grepl("log2quant", feature_1) ~ "logQ+l"),
            name = case_when(
            grepl("deseq2", name) ~ "DESeq2",
            grepl("edger", name) ~ "edgeR",
            grepl("tmm", name) ~ "TMM+l",
            grepl("voom", name) ~ "voom+l",
            grepl("vsn", name) ~ "vsn+l",
            grepl("log2quant", name) ~ "logQ+l")
        ) %>%
    mutate(id = paste0(feature_1, " - ", name)) %>%
    left_join(., (summarised_results_jaccard %>% dplyr::select(comparison, ranking, status)), by = c("id" = "comparison", "status" = "status")) %>%
    mutate(ranking = ranking)

jaccard_circular_plot_data <- results_jaccard %>%
    group_by(feature_1, name, status) %>%
    summarise(meanval = mean(value), stdev = sd(value)) %>%
    ungroup() %>%
    mutate(feature_1 = case_when(
            grepl("deseq2", feature_1) ~ "DESeq2",
            grepl("edger", feature_1) ~ "edgeR",
            grepl("tmm", feature_1) ~ "TMM+l",
            grepl("voom", feature_1) ~ "voom+l",
            grepl("vsn", feature_1) ~ "vsn+l",
            grepl("log2quant", feature_1) ~ "logQ+l"),
            name = case_when(
            grepl("deseq2", name) ~ "DESeq2",
            grepl("edger", name) ~ "edgeR",
            grepl("tmm", name) ~ "TMM+l",
            grepl("voom", name) ~ "voom+l",
            grepl("vsn", name) ~ "vsn+l",
            grepl("log2quant", name) ~ "logQ+l")
        ) %>%
    mutate(comparison = paste0(feature_1, " - ", name),
    ranking = rank(meanval)) %>%
    rbind(tibble(feature_1 = rep("", ceiling(nrow(.))), 
        name = rep("", ceiling(nrow(.))), 
        status = rep("", ceiling(nrow(.))), 
        meanval = rep(NA, ceiling(nrow(.))), 
        stdev = rep(NA, ceiling(nrow(.))),
        comparison = rep("", ceiling(nrow(.))),
        ranking = seq(nrow(.)+1, nrow(.) + ceiling(nrow(.)))
    )) %>%
    mutate(angle = 90 - 360 * (.$ranking-0.5)/nrow(.),
    hjust = ifelse(angle < -90, 0, 1)) %>%
    mutate(angle = ifelse(angle < -90, angle+180, angle))
    

jaccard_circular_barplot <- ggplot(jaccard_circular_plot_data) +
    geom_segment(x = 1, xend = 30.5, y = 1, yend = 1, color = "grey", linewidth = 0.5) +
    geom_col(aes(x = as.factor(ranking), y = meanval, fill = status),width = 0.5) +
    geom_boxplot(data = sorted_results_jaccard, aes(x = as.factor(ranking), y = value), outlier.alpha = 0, inherit.aes = FALSE) +
    coord_polar(clip = 'off', start = pi) +
    theme_minimal() +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(9,4), "cm"),
        legend.position = c(0.4,0.5),
        legend.title = element_blank(),
        text = element_text(size = 25)
    ) +
    scale_y_continuous(limits = c(-2, 1), oob = scales::rescale_none) +
    scale_fill_manual(values =  c(
                "Filtered" = "#5ab4ac",
                "Unfiltered" = "#d8b365"
            )) +
    # geom_errorbar(aes(x = ranking, ymin=meanval-stdev, ymax=meanval+stdev), position = position_dodge())  +
    annotate(geom = "rect", ymin = 1.08, ymax = 1.5, xmin = -Inf, xmax = (nrow(jaccard_circular_plot_data[!is.na(jaccard_circular_plot_data$meanval),]) + 0.5) - (nrow(jaccard_circular_plot_data[jaccard_circular_plot_data$meanval>0.5 & !is.na(jaccard_circular_plot_data$meanval),]) + 0.5),
                fill = "#d04a35", colour = NA, alpha = 0.2) +
    annotate(geom = "rect", ymin = 1.08, ymax = 1.5, xmin = (nrow(jaccard_circular_plot_data[!is.na(jaccard_circular_plot_data$meanval),]) + 0.5) - (nrow(jaccard_circular_plot_data[jaccard_circular_plot_data$meanval>0.5 & !is.na(jaccard_circular_plot_data$meanval),]) + 0.5), xmax = (nrow(jaccard_circular_plot_data[!is.na(jaccard_circular_plot_data$meanval),]) + 0.5) - (nrow(jaccard_circular_plot_data[jaccard_circular_plot_data$meanval>0.8 & !is.na(jaccard_circular_plot_data$meanval),]) + 0.5),
                fill = "#b7ba2e", colour = NA, alpha = 0.2) +
    annotate(geom = "rect", ymin = 1.08, ymax = 1.5, xmin = (nrow(jaccard_circular_plot_data[!is.na(jaccard_circular_plot_data$meanval),]) + 0.5) - (nrow(jaccard_circular_plot_data[jaccard_circular_plot_data$meanval>0.8 & !is.na(jaccard_circular_plot_data$meanval),]) + 0.5), xmax = nrow(jaccard_circular_plot_data[!is.na(jaccard_circular_plot_data$meanval),]) + 0.5,
                fill = "#439425", colour = NA, alpha = 0.2) +
    geom_text(aes(x = ranking, y = 1.7, label = comparison, hjust = hjust, angle = angle), color="black", size=7) +
    geom_text(aes(x = ranking, y = 1.12, label = round(meanval, digits = 2), hjust = hjust, angle = angle), color="black", size=6)


viewings_jaccard <- egg::ggarrange(jaccard_filtering_boxplot, jaccard_dataset_sina_plot, nrow=2)



options(repr.plot.width=18, repr.plot.height=18)
figure_3 <- ggdraw() +
    draw_plot(viewings_jaccard, 0.5, 0, 0.3, 1) +
    draw_plot(jaccard_circular_barplot, -0.05, 0, 1, 1) +
    draw_plot(upset_plots, 0.8, 0, 0.2, 1) +
    draw_label('A.', x = 0.03, y = 0.99, size = 25) +
    draw_label('B.', x = 0.5, y = 0.99, size = 25) +
    draw_label('C.', x = 0.5, y = 0.49, size = 25) +
    draw_label('D.', x = 0.8, y = 0.99, size = 25)

save_plot(
    filename = paste0("./flop_results/paper_plots/", "figure_3", ".png"),
    plot = figure_3,
    device = "png",
    dpi = 300,
    base_height = 18,
    base_width = 18
)

# Plot 3

rand_bio_plot <- ggplot(results_randindex_bio) +
    geom_boxplot(aes(x = factor(k), y = value)) +
    theme_cowplot() +
    # scale_fill_manual(values = c(
    #             "11" = "#5ab4ac",
    #             "32" = "#d8b365"
    #         )) +
    scale_y_continuous(limits = c(0,1)) +
    theme(axis.title = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = 'none',
        plot.margin = unit(c(1,1,1,1), "cm")
        )

get_k_results <- function(summary){
    resource <- summary$resource
    k_value <- summary$k
    results_subset <- results_randindex %>%
    filter(k == !!k_value, resource == !!resource)
    return(results_subset)
}

results_maxvar_summary <- results_randindex_summary %>% 
            group_by(resource) %>%
            summarise(
                max_sd = max(sd),
                k = filter(., sd == max_sd) %>% pull(k)
            ) %>% group_by(resource) %>% group_split() %>%
            purrr::map(., get_k_results) %>% 
            bind_rows()


rand_dataset_sina_plot <- results_maxvar_summary %>%
    mutate(resource = case_when(
            grepl("dorothea", resource) ~ "DoRothEA",
            grepl("msigdb_hallmarks", resource) ~ "MSigDB",
            grepl("progeny", resource) ~ "PROGENy")
        ) %>%
    ggplot(., aes(x = resource, y = value)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey", linewidth = 1) +
    ggforce::geom_sina(aes(color = status), jitter_y = FALSE, size = 2) +
    ggpubr::stat_compare_means(aes(group = status), size = 7, method = 'wilcox.test', method.args = list(alternative = 'greater', paired = TRUE, p.adjust.methods = "BH"), label = 'p.format') +
    scale_colour_manual(values = c(
                "filtered" = "#5ab4ac",
                "unfiltered" = "#d8b365"
            )) +
    theme_cowplot() +
    scale_y_continuous(limits = c(0,1)) +
    theme(axis.title = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 20),
        legend.position = 'none',
        plot.margin = unit(c(1,1,1,1), "cm")
        )

# rand_filtering_plot <- results_maxvar_summary %>%
#     mutate(resource = case_when(
#             grepl("dorothea", resource) ~ "DoRothEA",
#             grepl("msigdb_hallmarks", resource) ~ "MSigDB",
#             grepl("progeny", resource) ~ "PROGENy")
#         ) %>%
# ggplot(., aes(x = resource, y = value, fill = status)) +
#     geom_boxplot() +
#     scale_fill_manual(values = c(
#                 "filtered" = "#5ab4ac",
#                 "unfiltered" = "#d8b365"
#             )) +
#     geom_text(aes(x = 3, y = 1, label = "**"), color="black", size=10) +
#     geom_text(aes(x = 1, y = 1, label = "***"), color="black", size=10) +
#     theme_cowplot() +
#     scale_y_continuous(limits = c(0,1)) +
#     theme(axis.title = element_blank(),
#         legend.title = element_blank(),
#         axis.text = element_text(size = 20),
#         legend.position = 'none'
#         )

figure_4 <- ggdraw() +
    draw_plot(egg::ggarrange(rand_bio_plot, rand_dataset_sina_plot, widths = c(1,3)), 0, 0, 1, 1) +
    draw_label('A.', x = 0.03, y = 0.97, size = 25) +
    draw_label('B.', x = 0.3, y = 0.97, size = 25)

save_plot(
    filename = paste0("./flop_results/paper_plots/", "figure_4", ".png"),
    plot = figure_4,
    device = "png",
    dpi = 300,
    base_height = 6,
    base_width = 18
)
#


# full plot
left_column_plots <- egg::ggarrange(rank_filtering_boxplot, rank_dataset_sina_plot, jaccard_filtering_boxplot, jaccard_dataset_sina_plot, nrow = 4, padding = 1)

below_circle_left <- egg::ggarrange(plots_list$CCLE, plots_list$GSE186341, plots_list$GTex, plots_list$inters, rand_bio_plot, nrow = 5)

below_circle_right <- egg::ggarrange(rand_dataset_sina_plot, rand_filtering_plot, nrow = 2)

options(repr.plot.width=18, repr.plot.height=21)




rank_jaccard_plots <- ggdraw() + 
    draw_plot(circular_barplot, -0.25, 0.05, 1.4, 1) +
    draw_plot(left_column_plots, 0.03, 0, 0.28, 1) +
    draw_plot(below_circle_left, 0.35, 0, 0.3, 0.5) +
    draw_plot(below_circle_right, 0.7, 0, 0.3, 0.5) +
    draw_label('A.', x = 0.03, y = 0.99, size = 25) +
    draw_label('B.', x = 0.03, y = 0.74, size = 25) +
    draw_label('D.', x = 0.03, y = 0.49, size = 25) +
    draw_label('E.', x = 0.03, y = 0.24, size = 25) +
    draw_label('C.', x = 0.38, y = 0.99, size = 25) +
    draw_label('F.', x = 0.38, y = 0.49, size = 25) +
    draw_label('G.', x = 0.38, y = 0.1, size = 25) +
    draw_label('H.', x = 0.69, y = 0.49, size = 25) +
    draw_label('I.', x = 0.69, y = 0.24, size = 25)



save_plot(
    filename = paste0("./flop_results/paper_plots/", "filtering_plots", ".png"),
    plot = rank_jaccard_plots,
    device = "png",
    dpi = 300,
    base_height = 21,
    base_width = 18
)
