library(tidyverse)
library(ROCR)

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

for(i in 1:100){
    biomarkers_data_subset <- biomarkers_data %>%
        mutate(pipeline=paste0(status, '_', pipeline)) %>%
        group_by(main_dataset, bio_context, pipeline, status) %>%
        group_split() %>%
        purrr::map(., function(x){
            biocontext <- x %>% pull(bio_context) %>% unique()
            # detect and retrieve cell type that is contained in the biocontext name
            real_celltype <- strsplit(biocontext, '_') %>% unlist() %>% .[1]
            false_celltype <- non_selected_celltypes %>% sample(1)
            cell_types <- c(real_celltype, false_celltype)

            subset <- x %>% filter(items %in% cell_types)
            return(subset)
        }) %>% bind_rows()

    auc_tibble <- biomarkers_data_subset %>%
        mutate(is_tp = ifelse(str_detect(bio_context, items), 1, 0)) %>%
        group_by(pipeline) %>%
        arrange(desc(act), .by_group = TRUE) %>%
        group_split() %>%
        purrr::map2(., i, function(x, iter){
            PR_object <- prediction(x$act, x$is_tp) #Evaluate classification
            auc_obj <- ROCR::performance(PR_object, measure = "auc")
            pipeline <- x$pipeline %>% unique()
            # add row to auc_tibble
            auc_tibble <- tibble(random_iter = iter, pipeline = pipeline, auc = auc_obj@y.values[[1]])
            return(auc_tibble)
        }) %>% bind_rows()

    full_auc_tibble <- bind_rows(full_auc_tibble, auc_tibble)

    auprc_tibble <- biomarkers_data_subset %>%
        mutate(is_tp = ifelse(str_detect(bio_context, items), 1, 0)) %>%
        group_by(pipeline) %>%
        arrange(desc(act), .by_group = TRUE) %>%
        group_split() %>%
        purrr::map2(., i, function(x, iter){
            PR_object <- prediction(x$act, x$is_tp) #Evaluate classification
            auprc_obj <- ROCR::performance(PR_object, measure = "aucpr")
            pipeline <- x$pipeline %>% unique()
            # add row to auc_tibble
            auprc_tibble <- tibble(random_iter = iter, pipeline = pipeline, auprc = auprc_obj@y.values[[1]])
            return(auprc_tibble)
        }) %>% bind_rows()

    full_auprc_tibble <- bind_rows(full_auprc_tibble, auprc_tibble)
}

write_tsv(full_auc_tibble, '100iter_auc_tibble.tsv')
write_tsv(full_auprc_tibble, '100iter_auprc_tibble.tsv')


full_auc_tibble %>%
    mutate(pipeline = fct_reorder(pipeline, auc)) %>%
    ggplot(aes(x = pipeline, y = auc, fill = pipeline)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = 0.5, linetype = 'dashed') +
    ylim(0, 1) +
    labs(x = 'Pipeline', y = 'AUC', title = 'AUC of the different pipelines for the different celltypes') +
    theme(plot.title = element_text(hjust = 0.5))

full_auprc_tibble %>%
    mutate(pipeline = fct_reorder(pipeline, auprc)) %>%
    ggplot(aes(x = pipeline, y = auprc, fill = pipeline)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    # geom_hline(yintercept = 0.5, linetype = 'dashed') +
    ylim(0, 1) +
    labs(x = 'Pipeline', y = 'AUPRC', title = 'AUPRC of the different pipelines for the different celltypes and datasets') +
    theme(plot.title = element_text(hjust = 0.5))
# plot boxplots, pipeline in x axis and auc in y axis. Colours and sorted from highest to lowest. put ticks every 0.1



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
        false_celltype <- celltypes %>% setdiff(real_celltype)%>% sample(1)
        cell_types <- c(real_celltype, false_celltype)

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
        auc_tibble <- tibble(main_dataset = x$main_dataset, random_iter = 1, pipeline = pipeline, auc = auc_obj@y.values[[1]])
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

auc_tibble %>%
    mutate(pipeline = fct_reorder(pipeline, auc)) %>%
    ggplot(aes(x = pipeline, y = auc, fill = pipeline)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = 0.5, linetype = 'dashed') +
    ylim(0, 1) +
    labs(x = 'Pipeline', y = 'AUC', title = 'AUC of the different pipelines for the different celltypes') +
    theme(plot.title = element_text(hjust = 0.5))

auprc_tibble %>%
    mutate(pipeline = fct_reorder(pipeline, auprc)) %>%
    ggplot(aes(x = pipeline, y = auprc, fill = pipeline)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    # geom_hline(yintercept = 0.5, linetype = 'dashed') +
    ylim(0, 1) +
    labs(x = 'Pipeline', y = 'AUPRC', title = 'AUPRC of the different pipelines for the different celltypes and datasets') +
    theme(plot.title = element_text(hjust = 0.5))


# Per celltype

auc_tibble_cell <- biomarkers_data_subset %>%
    mutate(is_tp = ifelse(str_detect(bio_context, items), 1, 0)) %>%
    group_by(main_dataset, bio_context) %>%
    arrange(desc(act), .by_group = TRUE) %>%
    group_split() %>%
    purrr::map(., function(x){
        PR_object <- prediction(x$act, x$is_tp) #Evaluate classification
        auc_obj <- ROCR::performance(PR_object, measure = "auc")
        bio_context <- x$bio_context %>% unique()
        # add row to auc_tibble
        auc_tibble <- tibble(main_dataset = unique(x$main_dataset), random_iter = 1, bio_context = bio_context, auc = auc_obj@y.values[[1]])
        return(auc_tibble)
    }) %>% bind_rows()

auc_tibble_cell %>%
    mutate(pipeline = fct_reorder(bio_context, auc)) %>%
    ggplot(aes(x = pipeline, y = auc, fill = bio_context)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    # geom_hline(yintercept = 0.5, linetype = 'dashed') +
    ylim(0, 1) +
    labs(x = 'Pipeline', y = 'AUROC', title = 'AUROC of the different pipelines for the different celltypes and datasets') +
    theme(plot.title = element_text(hjust = 0.5))

auprc_tibble_cell <- biomarkers_data_subset %>%
    mutate(is_tp = ifelse(str_detect(bio_context, items), 1, 0)) %>%
    group_by(main_dataset, bio_context) %>%
    arrange(desc(act), .by_group = TRUE) %>%
    group_split() %>%
    purrr::map(., function(x){
        PR_object <- prediction(x$act, x$is_tp) #Evaluate classification
        auprc_obj <- ROCR::performance(PR_object, measure = "aucpr")
        bio_context <- x$bio_context %>% unique()
        # add row to auc_tibble
        auprc_tibble <- tibble(main_dataset = unique(x$main_dataset), bio_context = bio_context, auprc = auprc_obj@y.values[[1]])
        return(auprc_tibble)
    }) %>% bind_rows()

auprc_tibble_cell %>%
    mutate(pipeline = fct_reorder(bio_context, auprc)) %>%
    ggplot(aes(x = pipeline, y = auprc, fill = bio_context)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    # geom_hline(yintercept = 0.5, linetype = 'dashed') +
    ylim(0, 1) +
    labs(x = 'Pipeline', y = 'AUROC', title = 'AUROC of the different pipelines for the different celltypes and datasets') +
    theme(plot.title = element_text(hjust = 0.5))






