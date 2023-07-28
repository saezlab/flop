library(tidyverse)
path_file = './scripts/'
source(paste0(path_file, "rank_analysis_helper.R"))
source(paste0(path_file, "jaccard_analysis_helper.R"))

datasets <- list.files('./flop_results/diffexp/') %>% gsub('__deresults.tsv', '', .)
resources <- list.files('./scripts/dc_resources/') %>% gsub('__source.tsv', '', .)

for(dataset in datasets){
    diffexp_data <- list.files('./flop_results/diffexp/', pattern = dataset, full.names = TRUE) %>% read_tsv(.) %>%
        pivot_longer(stat, names_to = 'statparam', values_to = 'scores') %>%
        mutate(resource = 'DE') %>% dplyr::rename('items' = 'ID') %>%
        filter(!is.na(scores))
    diffexp_rank <- corr_analysis(diffexp_data, "pipeline")
    diffexp_jaccard <- jaccard_analysis(diffexp_data)

    subset_rank <- tibble()
    subset_jaccard <- tibble()
    for(resource in resources){
        resource_df <- read_tsv(paste0('./scripts/dc_resources/', resource, '__source.tsv')) 
        regulon <- resource_df %>% pull(target) %>% unique()
        subset_diffexp_data <- diffexp_data %>% filter(items %in% regulon)
        subset_diffexp_rank <- corr_analysis(subset_diffexp_data, "pipeline") %>% mutate(resource = paste('DE', !!resource, sep = ' '))
        subset_diffexp_jaccard <- jaccard_analysis(subset_diffexp_data) %>% mutate(resource = paste('DE', !!resource, sep = ' '))
        subset_rank <- bind_rows(subset_rank, subset_diffexp_rank)
        subset_jaccard <- bind_rows(subset_jaccard, subset_diffexp_jaccard)
    }

    funcomics_fullmerge <- list.files('./flop_results/funcomics/fullmerged/', pattern = paste0(dataset, '__fullmerge.tsv'), full.names = TRUE) %>% read_tsv(.)
    tfs <- funcomics_fullmerge %>% filter(resource == 'dorothea') %>% pull(items) %>% unique()
    subset_diffexp_data_tfs <- diffexp_data %>% filter(items %in% tfs)
    subset_diffexp_rank_tfs <- corr_analysis(subset_diffexp_data_tfs, "pipeline") %>% mutate(resource = 'DE TFs')
    subset_diffexp_jaccard_tfs <- jaccard_analysis(subset_diffexp_data_tfs) %>% mutate(resource = 'DE TFs')


    funcomics_rank <- list.files('./flop_results/funcomics/rank/', pattern = paste0(dataset, '__total__rank.tsv'), full.names = TRUE) %>% read_tsv(.)
    funcomics_jaccard <- list.files('./flop_results/funcomics/jaccard/', pattern = paste0(dataset, '__jaccard.tsv'), full.names = TRUE) %>% read_tsv(.)
    results_rank <- bind_rows(subset_diffexp_rank_tfs, diffexp_rank, subset_rank, funcomics_rank)
    results_jaccard <- bind_rows(subset_diffexp_jaccard_tfs, diffexp_jaccard, subset_jaccard, funcomics_jaccard)
    write_tsv(results_rank, file = paste0("./flop_results/funcomics/rank/", dataset, "__merged__rank.tsv"))
    write_tsv(results_jaccard, file = paste0("./flop_results/funcomics/jaccard/", dataset, "__merged__jaccard.tsv"))

}




# check which datasets pass the filtering criteria
sign_res <- tibble()
for(dataset in datasets){
    diffexp_data <- list.files('./flop_results/diffexp/', pattern = dataset, full.names = TRUE) %>% read_tsv(.) %>%
        group_by(status, bio_context, pipeline, main_dataset, subset) %>%
        count(padj <= 0.05) %>%
        ungroup() %>%
        complete(status, bio_context, pipeline, main_dataset, subset, `padj <= 0.05`, fill = list(n = 0)) %>%
        filter(`padj <= 0.05` == TRUE) %>%
        select(-`padj <= 0.05`)
    sign_res <- bind_rows(sign_res, diffexp_data)
}

heatmap_matrix <- sign_res %>% mutate(full_pipeline = paste(status, pipeline, sep = '+')) %>% 
    arrange(main_dataset) %>%
    select(-status, -pipeline, -bio_context, -subset) %>%
    pivot_wider(names_from = main_dataset, values_from = n) %>%
    column_to_rownames('full_pipeline') %>% as.matrix() %>% +1 %>% log()


ComplexHeatmap::Heatmap(heatmap_matrix, cluster_rows = TRUE, cluster_columns = TRUE,
row_split = factor(c(rep("Filtered", 6), rep("Unfiltered", 6)), levels = c("Filtered", "Unfiltered")),  
column_split = factor(colnames(heatmap_matrix), levels = colnames(heatmap_matrix)),
cluster_row_slices = TRUE,
cluster_column_slices = TRUE,
layer_fun = function(j, i, x, y, width, height, fill) {
    v = pindex(heatmap_matrix, i, j)
    if(min(v) >= log(2)) {
            grid.rect(gp = gpar(lwd = 10, fill = "transparent", col = "green"))
        }
}
)

# diffexp_data <- read_tsv('./flop_results/diffexp/') %>% 
#     pivot_longer(c(logFC, stat), names_to = 'statparam', values_to = 'scores') %>%
#     mutate(resource = 'DE') %>% dplyr::rename('items' = 'ID')
# diffexp_rank <- corr_analysis(diffexp_data, "pipeline")
# diffexp_jaccard <- jaccard_analysis(diffexp_data)

# funcomics_rank <- read_tsv('./flop_results/funcomics/rank/MetaHeart__total__rank.tsv')
# funcomics_jaccard <- read_tsv('./flop_results/funcomics/jaccard/MetaHeart__jaccard.tsv')
# results_rank <- bind_rows(diffexp_rank, funcomics_rank)
# results_jaccard <- bind_rows(diffexp_jaccard, funcomics_jaccard)
# write_tsv(total_results, file = "./flop_results/funcomics/rank/MetaHeart__merged__rank.tsv")
# write_tsv(jaccard_total, file = "./flop_results/funcomics/jaccard/MetaHeart__merged__jaccard.tsv")
