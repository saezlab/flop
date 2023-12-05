library(tidyverse)
library(ROCR)
library(cowplot)

results <- list.files(path='./flop_results/funcomics/fullmerged/', pattern='__fullmerge.tsv', full.names=TRUE) %>%
    lapply(., read_tsv) %>% bind_rows()
pk_cell_links <- read_delim(file='pk_data_celltypelinks.csv', delim = ';') %>% select(PK, Data)

biomarkers <- results %>%
    filter(resource == 'biomarkers', statparam == 'stat') %>%
    group_by(bio_context, pipeline, status) %>%
    arrange(desc(act), .by_group = TRUE)

# map names from the data column in pk_cell_links to the items column in biomarkers using the common column PK
biomarkers_data <- biomarkers %>%
    left_join(pk_cell_links, by = c('items' = 'PK'))

# if the name in the groundtruth column is not NA, then substitute that value in the items column
biomarkers_data <- biomarkers_data %>%
    mutate(items = ifelse(is.na(Data), items, Data)) %>%
    select(-Data)

celltypes <- biomarkers_data %>% pull(items) %>% unique()
selected_celltypes <- c('CM', 'Fib', 'Endo', 'Lymphoid', 'Myeloid', 'vSMCs', 'PC')

non_selected_celltypes <- celltypes %>% setdiff(selected_celltypes)

# random sample 7 celltypes from the non-selected celltypes
full_auc_tibble <- tibble()
full_auprc_tibble <- tibble()



# biomarkers_data_subset %>% group_by(main_dataset, bio_context, pipeline) %>% arrange(desc(abs(act)), .by_group = TRUE) %>% slice(1) %>% print(n=50)

# Just one iteration, only using the celltypes of the dataset and separating per dataset
biomarkers_data_subset <- biomarkers_data %>%
    mutate(pipeline=paste0(status, '_', pipeline)) %>%
    group_by(main_dataset, bio_context, pipeline, status) %>%
    group_split() %>%
    purrr::map(., function(x){
        biocontext <- x %>% pull(bio_context) %>% unique()
        celltypes <- selected_celltypes
        # detect and retrieve cell type that is contained in the biocontext name
        real_celltype <- strsplit(biocontext, '_') %>% unlist() %>% .[1]
        # remove the real celltype from the rest of celltypes
        false_celltypes <- celltypes %>% setdiff(real_celltype)
        cell_types <- c(real_celltype, false_celltypes)

        subset <- x %>% filter(items %in% cell_types)
        return(subset)
    }) %>% bind_rows()

auc_tibble <- biomarkers_data_subset %>%
    mutate(is_tp = ifelse(str_detect(bio_context, items), 1, 0)) %>%
    group_by(main_dataset, pipeline) %>%
    arrange(desc(act), .by_group = TRUE) %>%
    group_split() %>%
    purrr::map(., function(x){
        PR_object <- prediction(x$act, x$is_tp) #Evaluate classification
        auc_obj <- ROCR::performance(PR_object, measure = "auc")
        pipeline <- x$pipeline %>% unique()
        # add row to auc_tibble
        auc_tibble <- tibble(main_dataset = unique(x$main_dataset), random_iter = 1, pipeline = pipeline, auc = auc_obj@y.values[[1]])
        return(auc_tibble)
    }) %>% bind_rows()

auprc_tibble <- biomarkers_data_subset %>%
    mutate(is_tp = ifelse(str_detect(bio_context, items), 1, 0)) %>%
    group_by(main_dataset, pipeline) %>%
    arrange(desc(act), .by_group = TRUE) %>%
    group_split() %>%
    purrr::map(., function(x){
        PR_object <- prediction(x$act, x$is_tp) #Evaluate classification
        auprc_obj <- ROCR::performance(PR_object, measure = "aucpr")
        pipeline <- x$pipeline %>% unique()
        # add row to auc_tibble
        auprc_tibble <- tibble(main_dataset = unique(x$main_dataset), pipeline = pipeline, auprc = auprc_obj@y.values[[1]])
        return(auprc_tibble)
    }) %>% bind_rows()

biomarkers_results <- left_join(auc_tibble, auprc_tibble, by = c('main_dataset', 'pipeline'))

# scatter plot of AUC and AUPRC
biomarkers_results %>%
    ggplot(aes(x = auc, y = auprc, color = pipeline)) +
    geom_point() +
    theme_cowplot() +
    labs(x = 'AUC', y = 'AUPRC', title = 'AUC vs AUPRC for the different pipelines') +
    theme(plot.title = element_text(hjust = 0.5))

# boxplots of AUC and AUPRC
biomarkers_results %>%
    mutate(pipeline = fct_reorder(pipeline, auc)) %>%
    ggplot(aes(x = pipeline, y = auc, fill = pipeline)) +
    geom_boxplot() +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = 0.5, linetype = 'dashed') +
    ylim(0.5, 1) +
    labs(x = 'Pipeline', y = 'AUC', title = 'AUC of the different pipelines for the different celltypes') +
    theme(plot.title = element_text(hjust = 0.5))

biomarkers_results %>%
    mutate(pipeline = fct_reorder(pipeline, auprc)) %>%
    ggplot(aes(x = pipeline, y = auprc, fill = pipeline)) +
    geom_boxplot() +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    # geom_hline(yintercept = 0.5, linetype = 'dashed') +
    ylim(0, 1) +
    labs(x = 'Pipeline', y = 'AUPRC', title = 'AUPRC of the different pipelines for the different celltypes and datasets') +
    theme(plot.title = element_text(hjust = 0.5))


# BENCHMARK 2: tfs to cell type
library(stringr)

atac_tf_activity <- read_csv('atac_tf_activity.csv') %>% 
    dplyr::rename('TF' = '...1') %>%
        separate(TF, into = c("part1", "part2", "TF"), sep = "\\.") %>%
        select(-part1, -part2)

results_collectri <- results %>% filter(resource == 'collectri')
tfs_collectri <- results_collectri %>% pull(items) %>% unique()

# filter the TF activity only for TFs in collectri
atac_tf_activity_filtered <- atac_tf_activity %>% filter(TF %in% tfs_collectri)

tfs_groundtruth <- atac_tf_activity_filtered %>% pivot_longer(-TF, names_to='cell_type', values_to='act_atac') %>%
    group_by(cell_type) %>%
    # for tfs that are mapped to more that one cell line, keep the one with the highest activity
    filter(TF %in% tfs_collectri) %>%
    group_by(TF) %>%
    arrange(desc(act_atac), .by_group = TRUE) %>%
    slice(1) %>%
    group_by(cell_type) %>%
    arrange(desc(act_atac), .by_group = TRUE) %>%
    slice(1:10)

# count how many cell types per TF
tfs_groundtruth %>% group_by(cell_type) %>% count() %>% arrange(desc(n)) %>% print(n=50)

tfs_groundtruth_results_auc <- left_join(results_collectri, tfs_groundtruth, by = c('items' = 'TF')) %>%
    mutate(pipeline=paste0(status, '_', pipeline)) %>%
    mutate(is_tp = ifelse(str_detect(bio_context, cell_type), 1, 0)) %>%
    mutate(is_tp = ifelse(is.na(is_tp), 0, is_tp)) %>%
    group_by(main_dataset, pipeline) %>%
    arrange(desc(act), .by_group = TRUE) %>%
    group_split() %>%
    purrr::map(., function(x){
        PR_object <- prediction(x$act, x$is_tp) #Evaluate classification
        auc_obj <- ROCR::performance(PR_object, measure = "auc")
        pipeline <- x$pipeline %>% unique()
        # add row to auc_tibble
        auc_tibble <- tibble(main_dataset = unique(x$main_dataset), pipeline = pipeline, auc = auc_obj@y.values[[1]])
        return(auc_tibble)
    }) %>% bind_rows() 

# using only tfs that are markers of the cell type
tfs_groundtruth_results_auprc <- left_join(results_collectri, tfs_groundtruth, by = c('items' = 'TF')) %>%
    mutate(pipeline=paste0(status, '_', pipeline)) %>%
    mutate(is_tp = ifelse(str_detect(bio_context, cell_type), 1, 0)) %>%
    filter(is.na(is_tp) == FALSE) %>%
    group_by(main_dataset, pipeline) %>%
    arrange(desc(act), .by_group = TRUE) %>%
    group_split() %>%
    purrr::map(., function(x){
        PR_object <- prediction(x$act, x$is_tp) #Evaluate classification
        auprc_obj <- ROCR::performance(PR_object, measure = "aucpr")
        pipeline <- x$pipeline %>% unique()
        # add row to auc_tibble
        auc_tibble <- tibble(main_dataset = unique(x$main_dataset), pipeline = pipeline, auprc = auprc_obj@y.values[[1]])
        return(auc_tibble)
    }) %>% bind_rows() 

tfs_results <- left_join(tfs_groundtruth_results_auc, tfs_groundtruth_results_auprc, by = c('main_dataset', 'pipeline'))    

tfs_results %>%
    ggplot(aes(x = auc, y = auprc, color = pipeline)) +
    geom_point() +
    theme_cowplot() +
    labs(x = 'AUC', y = 'AUPRC', title = 'AUC vs AUPRC for the different pipelines') +
    theme(plot.title = element_text(hjust = 0.5))
    
tfs_results %>%
    mutate(pipeline = fct_reorder(pipeline, auc)) %>%
    ggplot(aes(x = pipeline, y = auc, fill = pipeline)) +
    geom_boxplot() +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = 0.5, linetype = 'dashed') +
    ylim(0.5, 1) +
    labs(x = 'Pipeline', y = 'AUC', title = 'AUC of the different pipelines for the different celltypes') +
    theme(plot.title = element_text(hjust = 0.5))

tfs_results %>%
    mutate(pipeline = fct_reorder(pipeline, auprc)) %>%
    ggplot(aes(x = pipeline, y = auprc, fill = pipeline)) +
    geom_boxplot() +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylim(0, 1) +
    labs(x = 'Pipeline', y = 'AUPRC', title = 'AUPRC of the different pipelines for the different celltypes') +
    theme(plot.title = element_text(hjust = 0.5))

tfs_groundtruth_results_cell <- inner_join(results_collectri, tfs_groundtruth, by = c('items' = 'TF')) %>%
    mutate(pipeline=paste0(status, '_', pipeline)) %>%
    mutate(is_tp = ifelse(str_detect(bio_context, cell_type), 1, 0)) %>%
    group_by(main_dataset, bio_context) %>%
    arrange(desc(act), .by_group = TRUE) %>%
    group_split() %>%
    purrr::map(., function(x){
        print(unique(x$is_tp))
        if(length(unique(x$is_tp)) == 1){
            return()
        }
        PR_object <- prediction(x$act, x$is_tp) #Evaluate classification
        auc_obj <- ROCR::performance(PR_object, measure = "auc")
        bio_context <- x$bio_context %>% unique()
        # add row to auc_tibble
        auc_tibble <- tibble(main_dataset = unique(x$main_dataset), random_iter = 1, bio_context = bio_context, auc = auc_obj@y.values[[1]])
        return(auc_tibble)
    }) %>% bind_rows()

    
tfs_groundtruth_results_cell %>%
    mutate(pipeline = fct_reorder(bio_context, auc)) %>%
    ggplot(aes(x = bio_context, y = auc, fill = bio_context)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = 0.5, linetype = 'dashed') +
    ylim(0, 1) +
    labs(x = 'Pipeline', y = 'AUC', title = 'AUC of the different pipelines for the different celltypes') +
    theme(plot.title = element_text(hjust = 0.5))
