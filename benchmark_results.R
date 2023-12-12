library(tidyverse)
library(ROCR)
library(cowplot)

# data exploration 
raw_metadata <- list.files(path='./benchmark/coldata', pattern='coldata.csv', full.names=TRUE) %>%
    purrr::map(., function(x){
        dataset <- x %>% strsplit(split = '/') %>% unlist() %>% .[4] %>% strsplit(split = '_') %>% unlist() %>% .[1]
        df <- read_csv(x) %>% mutate(dataset = dataset)
        return(df)
    }) %>% bind_rows()

num_samples <- raw_metadata %>% distinct(sample_id, dataset) %>% group_by(dataset) %>% summarise(n = n()) %>% arrange(desc(n)) %>% mutate(dataset = fct_reorder(dataset, n))

ggplot(num_samples, aes(x = dataset, y = n, fill=dataset)) +
    geom_bar(stat = 'identity') +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = 'Dataset', y = 'Number of samples', title = 'Number of samples per dataset') +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(family='Calibri'), legend.position = 'none')

# check how many genes per gene set
deresults <- read_tsv('flop_results/diffexp/Simonson2023__deresults.tsv')
hgenes <- deresults %>% pull(ID) %>% unique()
biomarkers <- read_tsv('./scripts/dc_resources/biomarkers__source.tsv') %>%
    left_join(pk_cell_links, by = c('source' = 'PK')) %>%
    filter(!is.na(Data))

# plot the percentage of genes in the biomarkers that are also in the DE results, per biomarker
biomarkers %>%
    mutate(overlap = ifelse(target %in% hgenes, 1, 0)) %>%
    group_by(source) %>%
    summarise(perc_overlap = sum(overlap) / n()) %>%
    mutate(source = fct_reorder(source, perc_overlap)) %>%
    ggplot(aes(x = source, y = perc_overlap, fill = source)) +
    geom_bar(stat = 'identity') +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = 'Gene set', y = 'Share of covered genes') +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(family='Calibri'), legend.position = 'none')

# plot the percentage of genes in the biomarkers that are also in the DE results

# benchmark 1: biomarkers
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

full_auc_tibble <- tibble()
full_auprc_tibble <- tibble()

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

# check the activity distributions per bio_context.
means <- biomarkers_data_subset %>%
  group_by(bio_context, items) %>%
  summarize(mean_act = mean(act, na.rm = TRUE)) %>%
  ungroup()

p <- biomarkers_data_subset %>%
  ggplot(aes(x = act, fill = items)) +
  geom_density(alpha = 0.2) +
  theme_cowplot() +
  theme(text = element_text(family='Calibri')) +
  facet_wrap(bio_context~items)

mean_labels <- means %>%
  mutate(label = paste("Mean =", round(mean_act, 2)),
         y_pos = 0.45) 
         
p + geom_text(data = mean_labels, aes(x = mean_act, y = y_pos, label = label), 
              size = 3, vjust = "top", color = "black")

# compute ROC curves, plot
roc_tibble <- biomarkers_data_subset %>%
    mutate(is_tp = ifelse(str_detect(bio_context, items), 1, 0)) %>%
    group_by(pipeline) %>%
    arrange(desc(act), .by_group = TRUE) %>%
    group_split() %>%
    purrr::map(., function(x){
        PR_object <- prediction(x$act, x$is_tp) #Evaluate classification
        pipeline <- x$pipeline %>% unique()
        roc_datapoints <- ROCR::performance(PR_object, measure = "tpr", x.measure = "fpr")
        roc_tibble <- tibble(pipeline = pipeline, fpr = roc_datapoints@x.values[[1]], tpr = roc_datapoints@y.values[[1]])
        return(roc_tibble)
    }) %>% bind_rows()

roc_tibble %>%
    ggplot(aes(x = fpr, y = tpr, color = pipeline, label = pipeline)) +
    geom_line() +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
    theme_cowplot() +
    labs(x = 'False positive rate', y = 'True positive rate', title = 'ROC curves') +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(family='Calibri')) 

# compute PR curves, plot
baseline_auprc_biomarkers <- biomarkers_data_subset %>% 
    ungroup() %>%
    group_split() %>%
    purrr::map(., function(x){
        x <- x %>% mutate(is_tp = ifelse(str_detect(bio_context, items), 1, 0))
        number_positives <- x %>% filter(is_tp == 1) %>% nrow()
        total_number <- x %>% nrow()
        baseline <- number_positives / total_number
        return(baseline)
    }) %>% unlist()

prc_tibble <- biomarkers_data_subset %>%
    mutate(is_tp = ifelse(str_detect(bio_context, items), 1, 0)) %>%
    group_by(pipeline) %>%
    arrange(desc(act), .by_group = TRUE) %>%
    group_split() %>%
    purrr::map(., function(x){
        PR_object <- prediction(x$act, x$is_tp) #Evaluate classification
        pipeline <- x$pipeline %>% unique()
        roc_datapoints <- ROCR::performance(PR_object, measure = "prec", x.measure = "rec")
        roc_tibble <- tibble(pipeline = pipeline, rec = roc_datapoints@x.values[[1]], prec = roc_datapoints@y.values[[1]])
        return(roc_tibble)
    }) %>% bind_rows()

prc_tibble %>%
    ggplot(aes(x = rec, y = prec, color = pipeline, label = pipeline)) +
    geom_line() +
    geom_hline(yintercept = baseline_auprc_biomarkers, linetype = 'dashed') +
    theme_cowplot() +
    ylim(0, 1) +
    xlim(0, 1) +
    labs(x = 'False positive rate', y = 'True positive rate', title = 'Precision-recall curves') +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(family='Calibri')) 


# compute AUC and AUPRC, plot
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
        auc_tibble <- tibble(main_dataset = unique(x$main_dataset), pipeline = pipeline, auc = auc_obj@y.values[[1]])
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
    ggplot(aes(x = auc, y = auprc, color = pipeline, label = pipeline)) +
    geom_point() +
    geom_text_repel(size = 5, family = 'Calibri') +
    theme_cowplot() +
    labs(x = 'AUC', y = 'AUPRC', title = 'AUC vs AUPRC for the different pipelines') +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(family='Calibri'), legend.position = 'none')
    #add labels to the points

# boxplots of AUC and AUPRC
biomarkers_results %>%
    mutate(pipeline = fct_reorder(pipeline, auc)) %>%
    ggplot(aes(x = pipeline, y = auc, fill = pipeline)) +
    geom_boxplot() +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = 0.5, linetype = 'dashed') +
    ylim(0.5, 1) +
    labs(x = 'Pipeline', y = 'AUC', title = 'AUROC per dataset and pipeline') +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(family='Calibri'), legend.position = 'none')

biomarkers_results %>%
    mutate(pipeline = fct_reorder(pipeline, auc)) %>%
    ggplot(aes(x = pipeline, y = auprc, fill = pipeline)) +
    geom_boxplot() +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = baseline_auprc_biomarkers, linetype = 'dashed') +
    labs(x = 'Pipeline', y = 'AUPRC', title = 'AUPRC per dataset and pipeline') +
    theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(family='Calibri'), 
        legend.position = 'none') +
    scale_y_continuous(breaks=seq(0, 1, 0.2), limits=c(0,1))
    # add ticks from 0 to 1 every 0.2 points


# BENCHMARK 2: tfs to cell type
library(stringr)

atac_tf_activity <- read_csv('atac_tf_activity.csv') %>% 
    dplyr::rename('TF' = '...1') %>%
        separate(TF, into = c("part1", "part2", "TF"), sep = "\\.") %>%
        select(-part1, -part2)

results_collectri <- results %>% filter(resource == 'collectri')
tfs_collectri <- results_collectri %>% pull(items) %>% unique()


tfs_collectri_geneset <- read_tsv('./scripts/dc_resources/collectri__source.tsv') %>%
    left_join(tfs_groundtruth, by = c('source' = 'TF')) %>%
    filter(!is.na(cell_type))

# plot the percentage of genes in the biomarkers that are also in the DE results, per biomarker
tfs_collectri_geneset %>%
    mutate(overlap = ifelse(target %in% hgenes, 1, 0)) %>%
    group_by(source, cell_type) %>%
    summarise(perc_overlap = sum(overlap) / n()) %>%
    mutate(source = fct_reorder(cell_type, mean(perc_overlap))) %>%
    ggplot(aes(x = perc_overlap, fill = source)) +
    geom_density(alpha = 0.2) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = 'Share of covered genes') +
    xlim(0, 1) +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(family='Calibri'), legend.position = 'none') +
    facet_wrap(~cell_type)

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
    slice(1:10) %>%
    mutate(cell_type = case_when(
        cell_type == 'Pericyte' ~ 'PC',
        TRUE ~ cell_type
    ))

# count how many cell types per TF
tfs_groundtruth %>% group_by(cell_type) %>% count() %>% arrange(desc(n)) %>% print(n=50)


tfs_data_subset <- left_join(results_collectri, tfs_groundtruth, by = c('items' = 'TF')) %>%
    mutate(pipeline=paste0(status, '_', pipeline)) %>%
    mutate(is_tp = ifelse(str_detect(bio_context, cell_type), 1, 0)) %>%
    filter(is.na(is_tp) == FALSE)

# check the activity distributions per bio_context.
means <- tfs_data_subset %>%
  group_by(bio_context, cell_type) %>%
  summarize(mean_act = mean(act, na.rm = TRUE)) %>%
  ungroup()

p <- tfs_data_subset %>%
  ggplot(aes(x = act, fill = bio_context)) +
  geom_density(alpha = 0.2) +
  theme_cowplot() +
  theme(text = element_text(family='Calibri')) +
  facet_wrap(bio_context~cell_type)

mean_labels <- means %>%
  mutate(label = paste("Mean =", round(mean_act, 2)),
         y_pos = 0.45) 
         
p + geom_text(data = mean_labels, aes(x = mean_act, y = y_pos, label = label), 
              size = 3, vjust = "top", color = "black")



# compute ROC curves, plot
roc_tibble <- tfs_data_subset %>%
    group_by(pipeline) %>%
    arrange(desc(act), .by_group = TRUE) %>%
    group_split() %>%
    purrr::map(., function(x){
        PR_object <- prediction(x$act, x$is_tp) #Evaluate classification
        pipeline <- x$pipeline %>% unique()
        roc_datapoints <- ROCR::performance(PR_object, measure = "tpr", x.measure = "fpr")
        roc_tibble <- tibble(pipeline = pipeline, fpr = roc_datapoints@x.values[[1]], tpr = roc_datapoints@y.values[[1]])
        return(roc_tibble)
    }) %>% bind_rows()

roc_tibble %>%
    ggplot(aes(x = fpr, y = tpr, color = pipeline, label = pipeline)) +
    geom_line() +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
    theme_cowplot() +
    labs(x = 'False positive rate', y = 'True positive rate', title = 'ROC curves') +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(family='Calibri')) 

# compute PR curves, plot
baseline_auprc_tfs <- tfs_data_subset %>% 
    ungroup() %>%
    group_split() %>%
    purrr::map(., function(x){
        number_positives <- x %>% filter(is_tp == 1) %>% nrow()
        total_number <- x %>% nrow()
        baseline <- number_positives / total_number
        return(baseline)
    }) %>% unlist()


prc_tibble <- tfs_data_subset %>%
    group_by(pipeline) %>%
    arrange(desc(act), .by_group = TRUE) %>%
    group_split() %>%
    purrr::map(., function(x){
        PR_object <- prediction(x$act, x$is_tp) #Evaluate classification
        pipeline <- x$pipeline %>% unique()
        roc_datapoints <- ROCR::performance(PR_object, measure = "prec", x.measure = "rec")
        roc_tibble <- tibble(pipeline = pipeline, rec = roc_datapoints@x.values[[1]], prec = roc_datapoints@y.values[[1]])
        return(roc_tibble)
    }) %>% bind_rows()

prc_tibble %>%
    ggplot(aes(x = rec, y = prec, color = pipeline, label = pipeline)) +
    geom_line() +
    geom_hline(yintercept = baseline_auprc_tfs, linetype = 'dashed') +
    theme_cowplot() +
    ylim(0, 1) +
    xlim(0, 1) +
    labs(x = 'False positive rate', y = 'True positive rate', title = 'Precision-recall curves') +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(family='Calibri')) 

tfs_groundtruth_results_auc <- tfs_data_subset %>%
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

tfs_groundtruth_results_auprc <- tfs_data_subset %>%
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

#scatter
tfs_results %>%
    mutate(pipeline = fct_reorder(pipeline, auc)) %>%
    ggplot(aes(x = auc, y = auprc, color = pipeline, label = pipeline)) +
    geom_point() +
    geom_text_repel(size = 5, family = 'Calibri') +
    theme_cowplot() +
    labs(x = 'AUC', y = 'AUPRC', title = 'AUC vs AUPRC for the different pipelines') +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(family='Calibri'), legend.position = 'none')
    
#boxplots
tfs_results %>%
    mutate(pipeline = fct_reorder(pipeline, auc)) %>%
    ggplot(aes(x = pipeline, y = auc, fill = pipeline)) +
    geom_boxplot() +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = 0.5, linetype = 'dashed') +
    ylim(0.5, 1) +
    labs(x = 'Pipeline', y = 'AUC', title = 'AUROC per dataset and pipeline') +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(family='Calibri'), legend.position = 'none')

tfs_results %>%
    mutate(pipeline = fct_reorder(pipeline, auc)) %>%
    ggplot(aes(x = pipeline, y = auprc, fill = pipeline)) +
    geom_boxplot() +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = baseline_auprc_tfs, linetype = 'dashed') +
    labs(x = 'Pipeline', y = 'AUPRC', title = 'AUPRC per dataset and pipeline') +
    theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(family='Calibri'), 
        legend.position = 'none') +
    scale_y_continuous(breaks=seq(0, 1, 0.2), limits=c(0,1))





# BENCHMARK 3: cytokines associated to hallmarks
# read in the cytokine activity data
cytokine_results <- read_tsv('./flop_results/funcomics/fullmerged/Atlascyt__fullmerge.tsv') %>%
    filter(resource == 'hallmarks_mouse', statparam == 'stat') %>%
    separate(bio_context, into = c('cell_type', 'rm'), sep = '_v_', remove=FALSE) %>%
    separate(cell_type, into = c('cell_type', 'true_cytokine'), sep = '_') %>%
    select(-rm)


celltypes_cytokines <- cytokine_results %>% pull(bio_context) %>% unique() %>% strsplit(., split = '_v_') %>% unlist() %>% tibble('bio_context'=.) %>% separate(bio_context, into=c('cell_type', 'pred_cytokine'), remove=TRUE) %>% distinct()

cytokine_hallmark_association <- read_tsv('unproc_data/hallmarks/GSE202186_map-scRNAseq-cytokines-dictionary.txt') %>% select(cytokine, hallmark) %>% distinct() %>% filter(!is.na(hallmark)) %>% dplyr::rename('true_hallmark' = 'hallmark') %>% distinct()
    # substitute the alpha hallmark
    # mutate(true_hallmark = ifelse(str_detect(true_hallmark, 'GAMMA'), 'HALLMARK_INTERFERON_ALPHA_RESPONSE', true_hallmark)) %>% distinct()

# plot number of cytokines per hallmark
cytokine_hallmark_association %>%
    group_by(true_hallmark) %>%
    count() %>%
    arrange(desc(n)) %>%
    ggplot(aes(x = true_hallmark, y = n, fill = true_hallmark)) +
    geom_bar(stat = 'identity') +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = 'Hallmark', y = 'Number of cytokines', title = 'Number of cytokines per hallmark') +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(family='Calibri'), legend.position = 'none')

cytokine_results <- full_join(cytokine_results, cytokine_hallmark_association, by = c('true_cytokine' = 'cytokine'))

selected_hallmarks <- cytokine_hallmark_association %>% pull(true_hallmark) %>% unique()

cytokines_data_subset <- cytokine_results %>%
    mutate(pipeline=paste0(status, '_', pipeline)) %>%
    group_by(main_dataset, bio_context, pipeline, status) %>%
    filter(items %in% selected_hallmarks) %>%
    filter(!is.na(true_hallmark)) %>%
    mutate(is_tp = ifelse(true_hallmark == items, 1, 0))

means <- cytokines_data_subset %>%
    group_by(true_hallmark, items) %>%
    summarize(mean_act = mean(act, na.rm = TRUE)) %>%
    ungroup()

p <- cytokines_data_subset %>%
    ggplot(aes(x = act, fill = items)) +
    geom_density(alpha = 0.2) +
    theme_cowplot() +
    theme(text = element_text(family='Calibri'), legend.position = 'none') +
    facet_wrap(true_hallmark~items)

mean_labels <- means %>%
    mutate(label = paste("Mean =", round(mean_act, 2)),
            y_pos = 0.45) 
         
p + geom_text(data = mean_labels, aes(x = mean_act, y = y_pos, label = label), 
              size = 3, vjust = "top", color = "black")


roc_tibble <- cytokines_data_subset %>%
    group_by(pipeline) %>%
    arrange(desc(act), .by_group = TRUE) %>%
    group_split() %>%
    purrr::map(., function(x){
        PR_object <- prediction(x$act, x$is_tp) #Evaluate classification
        pipeline <- x$pipeline %>% unique()
        roc_datapoints <- ROCR::performance(PR_object, measure = "tpr", x.measure = "fpr")
        roc_tibble <- tibble(pipeline = pipeline, fpr = roc_datapoints@x.values[[1]], tpr = roc_datapoints@y.values[[1]])
        return(roc_tibble)
    }) %>% bind_rows()

roc_tibble %>%
    ggplot(aes(x = fpr, y = tpr, color = pipeline, label = pipeline)) +
    geom_line() +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
    theme_cowplot() +
    labs(x = 'False positive rate', y = 'True positive rate', title = 'ROC curves') +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(family='Calibri')) 

# compute PR curves, plot
baseline_auprc_hallmarks <- cytokines_data_subset %>% 
    ungroup() %>%
    group_split() %>%
    purrr::map(., function(x){
        number_positives <- x %>% filter(is_tp == 1) %>% nrow()
        total_number <- x %>% nrow()
        baseline <- number_positives / total_number
        return(baseline)
    }) %>% unlist()


prc_tibble <- cytokines_data_subset %>%
    group_by(pipeline) %>%
    arrange(desc(act), .by_group = TRUE) %>%
    group_split() %>%
    purrr::map(., function(x){
        PR_object <- prediction(x$act, x$is_tp) #Evaluate classification
        pipeline <- x$pipeline %>% unique()
        roc_datapoints <- ROCR::performance(PR_object, measure = "prec", x.measure = "rec")
        roc_tibble <- tibble(pipeline = pipeline, rec = roc_datapoints@x.values[[1]], prec = roc_datapoints@y.values[[1]])
        return(roc_tibble)
    }) %>% bind_rows()

prc_tibble %>%
    ggplot(aes(x = rec, y = prec, color = pipeline, label = pipeline)) +
    geom_line() +
    geom_hline(yintercept = baseline_auprc_hallmarks, linetype = 'dashed') +
    theme_cowplot() +
    ylim(0, 1) +
    xlim(0, 1) +
    labs(x = 'Recall', y = 'Precision', title = 'Precision-recall curves') +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(family='Calibri')) 

hallmarks_groundtruth_results_auc <- cytokines_data_subset %>%
    group_by(pipeline) %>%
    arrange(desc(act), .by_group = TRUE) %>%
    group_split() %>%
    purrr::map(., function(x){
        PR_object <- prediction(x$act, x$is_tp) #Evaluate classification
        auc_obj <- ROCR::performance(PR_object, measure = "auc")
        pipeline <- x$pipeline %>% unique()
        # add row to auc_tibble
        auc_tibble <- tibble(pipeline = pipeline, auc = auc_obj@y.values[[1]])
        return(auc_tibble)
    }) %>% bind_rows() 

hallmarks_groundtruth_results_auprc <- cytokines_data_subset %>%
    group_by(pipeline) %>%
    arrange(desc(act), .by_group = TRUE) %>%
    group_split() %>%
    purrr::map(., function(x){
        PR_object <- prediction(x$act, x$is_tp) #Evaluate classification
        auprc_obj <- ROCR::performance(PR_object, measure = "aucpr")
        pipeline <- x$pipeline %>% unique()
        # add row to auc_tibble
        auc_tibble <- tibble(pipeline = pipeline, auprc = auprc_obj@y.values[[1]])
        return(auc_tibble)
    }) %>% bind_rows() 

hallmarks_results <- left_join(hallmarks_groundtruth_results_auc, hallmarks_groundtruth_results_auprc, by = c('pipeline'))   

order_pipelines_auc <- hallmarks_results %>% arrange(desc(auc)) %>% pull(pipeline)

#scatter
hallmarks_results %>%
    mutate(pipeline = factor(pipeline, levels = order_pipelines_auc)) %>%
    ggplot(aes(x = auc, y = auprc, color = pipeline, label = pipeline)) +
    geom_point() +
    ggrepel::geom_text_repel(size = 5, family = 'Calibri') +
    theme_cowplot() +
    labs(x = 'AUC', y = 'AUPRC', title = 'AUC vs AUPRC for the different pipelines') +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(family='Calibri'), legend.position = 'none')


hallmarks_groundtruth_results_auc <- cytokines_data_subset %>%
    group_by(cell_type, pipeline) %>%
    arrange(desc(act), .by_group = TRUE) %>%
    group_split() %>%
    purrr::map(., function(x){
        PR_object <- prediction(x$act, x$is_tp) #Evaluate classification
        auc_obj <- ROCR::performance(PR_object, measure = "auc")
        pipeline <- x$pipeline %>% unique()
        # add row to auc_tibble
        auc_tibble <- tibble(cell_type = unique(x$cell_type), pipeline = pipeline, auc = auc_obj@y.values[[1]])
        return(auc_tibble)
    }) %>% bind_rows() 

hallmarks_groundtruth_results_auprc <- cytokines_data_subset %>%
    group_by(cell_type, pipeline) %>%
    arrange(desc(act), .by_group = TRUE) %>%
    group_split() %>%
    purrr::map(., function(x){
        PR_object <- prediction(x$act, x$is_tp) #Evaluate classification
        auprc_obj <- ROCR::performance(PR_object, measure = "aucpr")
        pipeline <- x$pipeline %>% unique()
        # add row to auc_tibble
        auc_tibble <- tibble(cell_type = unique(x$cell_type), pipeline = pipeline, auprc = auprc_obj@y.values[[1]])
        return(auc_tibble)
    }) %>% bind_rows() 

hallmarks_results <- left_join(hallmarks_groundtruth_results_auc, hallmarks_groundtruth_results_auprc, by = c('cell_type', 'pipeline'))   
#boxplots
hallmarks_results %>%
    mutate(pipeline = fct_reorder(pipeline, auc)) %>%
    ggplot(aes(x = pipeline, y = auc, fill = pipeline)) +
    geom_boxplot() +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = 0.5, linetype = 'dashed') +
    ylim(0.5, 1) +
    labs(x = 'Pipeline', y = 'AUC', title = 'AUROC per dataset and pipeline') +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(family='Calibri'), legend.position = 'none')

hallmarks_results %>%
    mutate(pipeline = fct_reorder(pipeline, auc)) %>%
    ggplot(aes(x = pipeline, y = auprc, fill = pipeline)) +
    geom_boxplot() +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = baseline_auprc_hallmarks, linetype = 'dashed') +
    labs(x = 'Pipeline', y = 'AUPRC', title = 'AUPRC per dataset and pipeline') +
    theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(family='Calibri'), 
        legend.position = 'none') +
    scale_y_continuous(breaks=seq(0, 1, 0.2), limits=c(0,1))


