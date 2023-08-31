library(tidyverse)
library(ggforce)
library(ggpubr)
library(cowplot)
library(ggrepel)

# reheat
reheat_de <- read_tsv('./flop_results/diffexp/reheat__deresults.tsv') %>%
  pivot_longer(-ID) %>%
  separate(name, into = c('colid', 'filtering', 'de', 'contrast', 'dataset', 'study'), sep = '__') %>%
  pivot_wider(names_from = colid, values_from = value)
full_df <- list.files('./flop_results/funcomics/fullmerged/', full.names = TRUE) %>% lapply( read_tsv) %>% bind_rows() %>%
  dplyr::filter(main_dataset %in% c("Liu15", "Pepin19", "Schiano17", "Spurell19", "Yang14"))
similarity_df <- list.files('flop_results/funcomics/jaccard/', full.names = TRUE) %>% lapply( read_tsv) %>% bind_rows() %>%
  dplyr::filter(main_dataset %in% c("Liu15", "Pepin19", "Schiano17", "Spurell19", "Yang14"))
correlation_df <-  list.files('flop_results/funcomics/rank/', full.names = TRUE) %>% lapply( read_tsv) %>% bind_rows() %>%
  dplyr::filter(main_dataset %in% c("Liu15", "Pepin19", "Schiano17", "Spurell19", "Yang14"))

# ccle and panacea
panacea_spearman <- read_tsv('flop_results/funcomics/rank/GSE186341__total__rank.tsv')
ccle_spearman <- read_tsv('flop_results/funcomics/rank/CCLE__total__rank.tsv')
panacea_similarity <- read_tsv('flop_results/funcomics/overlap/GSE186341__overlap.tsv')
ccle_similarity <- read_tsv('flop_results/funcomics/overlap/CCLE__overlap.tsv')


p_cutoff <- 0.05
lfc_cutoff <- log2(2)

#### fig 2 #### 
plot_reheat_detailed <- function(int_pipelines, int_dataset) {

  # volcanos
  toplot <- reheat_de %>%
    mutate(pipeline = paste0(filtering, '+', str_replace(de , 'NA\\+', ''))) %>%
    dplyr::filter(pipeline %in% int_pipelines & study == int_dataset) %>%
    mutate(status = case_when(padj <= p_cutoff & logFC >= lfc_cutoff ~ 'Up',
                              padj <= p_cutoff & logFC <= -lfc_cutoff ~ 'Down',
                              TRUE ~ 'Not')) %>%
    drop_na(padj, logFC)
  
  diff_feat_table <- toplot %>%
    group_by(pipeline, status) %>%
    dplyr::count()
  
  set.seed(149)
  vplot <- toplot %>%
    ggplot(aes(x = logFC, y = -log10(padj), fill = status)) +
    scale_fill_manual(values = c('Up' = 'red', 'Down' = 'blue', 'Not' = 'black')) +
    scale_color_manual(values = c('Up' = 'red', 'Down' = 'blue', 'Not' = 'black')) +
    geom_point(pch = 21, alpha = 0.3) +
    facet_wrap(facets = vars(pipeline), scales = 'free', nrow = 1) +
    geom_label(data = diff_feat_table %>% dplyr::filter(status == 'Up'), 
               aes(x = Inf, y = Inf, label = n), hjust = 1, vjust = 1, color = 'red', alpha = 0.5, fill = 'white') +
    geom_label(data = diff_feat_table %>% dplyr::filter(status == 'Down'), 
               aes(x = -Inf, y = Inf, label = n), hjust = 0, vjust = 1, color = 'blue', alpha = 0.5, fill = 'white') +
    geom_label(data = diff_feat_table %>% dplyr::filter(status == 'Not'), 
               aes(x = 0, y = Inf, label = n),  vjust = 1, color = 'black', alpha = 0.5, fill = 'white') +
    theme_cowplot()
  
  # gene corr
  toplot <- reheat_de %>%
    mutate(pipeline = paste0(filtering, '+', str_replace(de , 'NA\\+', ''))) %>%
    dplyr::filter(pipeline %in% int_pipelines & study == int_dataset) %>%
    dplyr::select(pipeline, stat, ID)  %>%
    pivot_wider(names_from = pipeline, values_from = stat) 
  
  genecorrplot <- ggplot(toplot, aes(x = !!sym(int_pipelines[1]), y = !!sym(int_pipelines[2]), fill = !!sym(int_pipelines[1]))) +
    scale_fill_gradient2(high =  'red',low = 'blue', mid = 'black', midpoint = 0) +
    geom_point(alpha = 0.5, pch = 21) +
    stat_cor(method = 'spearman') +
    ggtitle('Gene space') +
    theme_cowplot() +
    theme(legend.position = 'none')
  
  # genespace corr
  hallmarks_df <- full_df %>%
    mutate(pipeline = str_replace_all(obs_id, 'yes_v_no__|__ulm', '') %>% str_replace_all('__', '\\+') %>% str_replace_all('\\+NA', '')) %>%
    dplyr::filter(pipeline %in%int_pipelines & main_dataset == int_dataset & statparam == 'stat') %>%
    dplyr::filter(resource == 'msigdb_hallmarks')
  
  toplot <- hallmarks_df %>%
    dplyr::select(items, pipeline, act)  %>%
    pivot_wider(names_from = pipeline, values_from = act) 
  
  genesetcorrplot <- ggplot(toplot, aes(x = !!sym(int_pipelines[1]), y = !!sym(int_pipelines[2]), fill = !!sym(int_pipelines[1]))) +
    scale_fill_gradient2(high =  'red',low = 'blue', mid = 'black', midpoint = 0) +
    geom_point(alpha = 0.5, pch = 21) +
    stat_cor(method = 'spearman') +
    ggtitle('Gene set space') +
    theme_cowplot() +
    theme(legend.position = 'none')
  
  # top gene set barplot
  toplot <- hallmarks_df %>%
    mutate(status = case_when(act > 0 ~ 'Up', act <= 0 ~ 'Down', TRUE ~ 'Not')) %>%
    group_by(pipeline, status) %>%
    arrange(desc(abs(act))) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    arrange(act) %>%
    mutate(items = fct_inorder(items))
  
  barplot <- toplot %>%
    ggplot(aes( x = act, y = items, fill = status)) +
    geom_col() +
    xlab("Activity score") + ylab('Gene set') +
    scale_fill_manual(values = c('Up' = 'red', 'Down' = 'blue')) +
    facet_grid(cols = vars(pipeline)) +
    theme_cowplot() 
  
  p <- plot_grid(
    plot_grid(vplot, genecorrplot, rel_widths = c(1, 0.55), nrow = 1, labels = c('A', 'B')),
    plot_grid(genesetcorrplot, barplot, rel_widths = c( 0.55,1), nrow = 1, labels = c('C', 'D')),
    nrow = 2)

  
  return(p)
  
}

p <- plot_reheat_detailed( c('unfiltered+edger', 'filtered+deseq2'), 'Spurell19')
ggsave('flop_results/plots/fig2.png', p,  width = 15, height = 10, dpi = 300)
supp1 <- plot_reheat_detailed(c('unfiltered+tmm+limma', 'filtered+tmm+limma'), 'Spurell19')
ggsave('flop_results/plots/supp1.png', supp1,  width = 15, height = 10, dpi = 300)
supp2 <- plot_reheat_detailed(c('unfiltered+edger', 'filtered+deseq2'), 'Yang14')
ggsave('flop_results/plots/supp2.png', supp2,  width = 15, height = 10, dpi = 300)

#### fig 3 ####
correlation_toplot <- correlation_df %>%
  mutate(id = paste0(pmin(feature_1, name), ' - ', pmax(feature_1, name))) %>%
  dplyr::filter(type == 'correlation' & statparam == 'stat') %>%
  group_by(id, resource) %>%
  summarise(spearman_average = mean(value), spearman_sem =  sd(value)/sqrt(n())) %>%
  ungroup() %>%
  separate(id, into = c('pipeline_a', 'pipeline_b'), sep = ' - ') %>%
  
  dplyr::filter(pipeline_a != pipeline_b) %>%
  mutate(pipeline_a = factor(pipeline_a, levels = sort(unique(pipeline_a))),
         pipeline_b = factor(pipeline_b, levels = sort(unique(pipeline_b))))

p <- ggplot(correlation_toplot, aes(x = pipeline_a, y = pipeline_b, size = -spearman_sem, fill = spearman_average)) +
  geom_point(pch = 21) +
  scale_fill_viridis_c(option = "plasma", limits = c(0, 1)) +
  facet_wrap( vars(resource)) +
  xlab('Pipeline A') + ylab('Pipeline B') +
  guides(x = guide_axis(angle = 60)) +
  theme_cowplot() 

ggsave('flop_results/plots/fig3a.png', p, width = 12, height = 10, dpi = 300)

# values for text
correlation_toplot %>% 
  mutate(unfil_limma_a = ifelse(grepl('unfiltered', pipeline_a) & grepl('limma', pipeline_a), 'Yes', 'No')) %>%
  mutate(unfil_limma_b = ifelse(grepl('unfiltered', pipeline_b) & grepl('limma', pipeline_b), 'Yes', 'No')) %>%
  mutate(unfil_limma = ifelse(unfil_limma_a == 'Yes' | unfil_limma_b == 'Yes', 'Yes', 'No')) %>%
  group_by(unfil_limma, resource) %>%
  summarise(mean(spearman_average))

correlation_toplot %>% 
  mutate(unfil_limma_a = ifelse(grepl('unfiltered', pipeline_a) & grepl('limma', pipeline_a), 'Yes', 'No')) %>%
  mutate(unfil_limma_b = ifelse(grepl('unfiltered', pipeline_b) & grepl('limma', pipeline_b), 'Yes', 'No')) %>%
  mutate(unfil_limma = ifelse(unfil_limma_a == 'Yes' | unfil_limma_b == 'Yes', 'Yes', 'No')) %>%
  filter(unfil_limma == 'Yes') %>%
  mutate(is_de = ifelse(resource == 'DE', 'Yes', 'No')) %>%
wilcox.test(spearman_average ~ is_de, data = ., alternative = "less", p.adjust.methods = "BH")


correlation_toplot %>% 
  mutate(filtered_a = ifelse(!grepl('unfiltered', pipeline_a), 'Yes', 'No')) %>%
  mutate(filtered_b = ifelse(!grepl('unfiltered', pipeline_b), 'Yes', 'No')) %>%
  mutate(filtered = ifelse(filtered_a == 'Yes' & filtered_b == 'Yes', 'Yes', 'No')) %>%
  group_by(filtered, resource) %>%
  summarise(mean(spearman_average))

correlation_toplot %>% 
  mutate(filtered_a = ifelse(!grepl('unfiltered', pipeline_a), 'Yes', 'No')) %>%
  mutate(filtered_b = ifelse(!grepl('unfiltered', pipeline_b), 'Yes', 'No')) %>%
  mutate(filtered = ifelse(filtered_a == 'Yes' & filtered_b == 'Yes', 'Yes', 'No')) %>%
  filter(filtered == 'Yes') %>%
  mutate(is_de = ifelse(resource == 'DE', 'Yes', 'No')) %>%
wilcox.test(spearman_average ~ is_de, data = ., alternative = "less", p.adjust.methods = "BH")

correlation_toplot %>% 
  mutate(unfil_deseq_edger_a = ifelse(grepl('unfiltered', pipeline_a) & (grepl('edger', pipeline_a) | grepl('deseq2', pipeline_a)), 'Yes', 'No')) %>%
  mutate(unfil_deseq_edger_b =  ifelse(grepl('unfiltered', pipeline_b) & (grepl('edger', pipeline_b) | grepl('deseq2', pipeline_b)), 'Yes', 'No')) %>%
  mutate(unfil_deseq_edger = ifelse(unfil_deseq_edger_a == 'Yes' & unfil_deseq_edger_b == 'Yes', 'Yes', 'No')) %>%
  group_by(unfil_deseq_edger, resource) %>%
  summarise(mean(spearman_average))

correlation_toplot %>% 
  mutate(unfil_deseq_edger_a = ifelse(grepl('unfiltered', pipeline_a) & (grepl('edger', pipeline_a) | grepl('deseq2', pipeline_a)), 'Yes', 'No')) %>%
  mutate(unfil_deseq_edger_b =  ifelse(grepl('unfiltered', pipeline_b) & (grepl('edger', pipeline_b) | grepl('deseq2', pipeline_b)), 'Yes', 'No')) %>%
  mutate(unfil_deseq_edger = ifelse(unfil_deseq_edger_a == 'Yes' & unfil_deseq_edger_b == 'Yes', 'Yes', 'No')) %>%
  filter(unfil_deseq_edger == 'Yes') %>%
  mutate(is_de = ifelse(resource == 'DE', 'Yes', 'No')) %>%
wilcox.test(spearman_average ~ is_de, data = ., alternative = "less", p.adjust.methods = "BH")

#### fig 4 ####
similarity_toplot <- similarity_df %>%
  mutate(id = paste0(pmin(feature_1, name), ' - ', pmax(feature_1, name))) %>%
  dplyr::filter(type == 'agreement' & statparam == 'stat') %>%
  dplyr::filter(main_dataset %in% c("Liu15", "Pepin19", "Schiano17", "Spurell19", "Yang14")) %>%
  group_by(id, resource) %>%
  summarise(similarity_average = mean(value), similarity_sem =  sd(value)/sqrt(n())) %>%
  ungroup() %>%
  separate(id, into = c('pipeline_a', 'pipeline_b'), sep = ' - ') %>%
  dplyr::filter(pipeline_a != pipeline_b) %>%
  
  mutate(pipeline_a = factor(pipeline_a, levels = sort(unique(pipeline_a))),
         pipeline_b = factor(pipeline_b, levels = sort(unique(pipeline_b))))

p <- ggplot(similarity_toplot, aes(x = pipeline_a, y = pipeline_b, size = -similarity_sem, fill = similarity_average)) +
  geom_point(pch = 21) +
  scale_fill_viridis_c(option = "plasma", limits = c(0, 1)) +
  facet_wrap( vars(resource)) +
  xlab('Pipeline A') + ylab('Pipeline B') +
  guides(x = guide_axis(angle = 60)) +
  theme_cowplot() 

ggsave('flop_results/plots/fig4a.png', p, width = 12, height = 10, dpi = 300)

# values for text
similarity_toplot %>% 
  mutate(unfil_limma_a = ifelse(pipeline_a == 'unfiltered-vsn+limma', 'Yes', 'No')) %>%
  mutate(unfil_limma_b = ifelse(pipeline_b == 'unfiltered-vsn+limma', 'Yes', 'No')) %>%
  mutate(unfil_limma = ifelse(unfil_limma_a == 'Yes' | unfil_limma_b == 'Yes', 'Yes', 'No')) %>%
  group_by(unfil_limma, resource) %>%
  summarise(mean(similarity_average))

similarity_toplot %>% 
  mutate(unfil_limma_a = ifelse(pipeline_a == 'unfiltered-vsn+limma', 'Yes', 'No')) %>%
  mutate(unfil_limma_b = ifelse(pipeline_b == 'unfiltered-vsn+limma', 'Yes', 'No')) %>%
  mutate(unfil_limma = ifelse(unfil_limma_a == 'Yes' | unfil_limma_b == 'Yes', 'Yes', 'No')) %>%
  filter(unfil_limma == 'Yes') %>%
  mutate(is_de = ifelse(resource == 'DE', 'Yes', 'No')) %>%
wilcox.test(similarity_average ~ is_de, data = ., alternative = "less", p.adjust.methods = "BH")

similarity_toplot %>% 
  mutate(unfiltered_limma_a = ifelse(grepl('unfiltered', pipeline_a) & grepl('limma', pipeline_a), 'Yes', 'No')) %>%
  mutate(unfiltered_limma_b = ifelse(grepl('unfiltered', pipeline_b) & grepl('limma', pipeline_a), 'Yes', 'No')) %>%
  mutate(unfiltered_limma = ifelse(unfiltered_limma_a == 'Yes' | unfiltered_limma_b == 'Yes', 'Yes', 'No')) %>%
  group_by(unfiltered_limma, resource) %>%
  summarise(mean(similarity_average))

similarity_toplot %>% 
  mutate(unfil_deseq_edger_a = ifelse(grepl('unfiltered', pipeline_a) & (grepl('edger', pipeline_a) | grepl('deseq2', pipeline_a)), 'Yes', 'No')) %>%
  mutate(unfil_deseq_edger_b =  ifelse(grepl('unfiltered', pipeline_b) & (grepl('edger', pipeline_b) | grepl('deseq2', pipeline_b)), 'Yes', 'No')) %>%
  mutate(unfil_deseq_edger = ifelse(unfil_deseq_edger_a == 'Yes' & unfil_deseq_edger_b == 'Yes', 'Yes', 'No')) %>%
  group_by(unfil_deseq_edger, resource) %>%
  summarise(mean(spearman_average))

# sup about variance
sim_values <- similarity_df %>%
  mutate(id = paste0(pmin(feature_1, name), ' - ', pmax(feature_1, name))) %>%
  dplyr::filter(type == 'agreement' & statparam == 'stat') %>%
  dplyr::filter(!main_dataset %in% c('vanHeesch19', 'Sweet18')) %>%
  transmute(id = paste(id, main_dataset, resource, sep = '___'), similarity = value)
  

cor_values <- correlation_df %>%
  mutate(id = paste0(pmin(feature_1, name), ' - ', pmax(feature_1, name))) %>%
  dplyr::filter(feature_1 != name) %>%
  dplyr::filter(type == 'correlation' & statparam == 'stat') %>%
  dplyr::filter(!main_dataset %in% c('vanHeesch19', 'Sweet18')) %>%
  transmute(id = paste(id, main_dataset, resource, sep = '___'), correlation = value)

toplot <- full_join(cor_values, sim_values, by = 'id') %>%
  pivot_longer(-id)

p <- ggplot(toplot, aes(x = name, y = value)) +
  geom_line(aes(group = id), alpha = 0.05) +
  geom_sina(alpha = 0.5, pch = 21) +
  geom_boxplot(alpha = 0.5, pch = 21) +
  stat_compare_means(method = 'wilcox.test', paired = TRUE, comparisons = list(c('correlation', 'similarity'))) +
  theme_cowplot()

ggsave('flop_results/plots/suppfig3.png', p, width = 5, height = 5, dpi = 300)

#### fig 5 ####
spearman_df <- bind_rows(ccle = ccle_spearman, panacea = panacea_spearman, .id = 'dataset') %>%
  mutate(id = paste0(pmin(feature_1, name), ' - ', pmax(feature_1, name))) %>%
  dplyr::filter(type == 'correlation' & statparam == 'stat') %>%
  group_by(dataset, id, resource) %>%
  summarise(spearman_average = mean(value), spearman_sem =  sd(value)/sqrt(n())) %>%
  ungroup() %>%
  separate(id, into = c('pipeline_a', 'pipeline_b'), sep = ' - ') %>%
  dplyr::filter(pipeline_a != pipeline_b) %>%
  mutate(pipeline_a = factor(pipeline_a, levels = sort(unique(pipeline_a))),
         pipeline_b = factor(pipeline_b, levels = sort(unique(pipeline_b)))) 

p <- ggplot(spearman_df, aes(x = pipeline_a, y = pipeline_b, size = -spearman_sem, fill = spearman_average)) +
  geom_point(pch = 21) +
  scale_fill_viridis_c(option = "plasma", limits = c(0, 1)) +
  facet_grid( cols = vars(resource), rows = vars(dataset)) +
  xlab('Pipeline A') + ylab('Pipeline B') +
  guides(x = guide_axis(angle = 60)) +
  theme_cowplot() 

ggsave('flop_results/plots/fig5.png', p, width = 14, height = 8, dpi = 300)

similarity_df <- bind_rows(ccle = ccle_similarity, panacea = panacea_similarity, .id = 'dataset') %>%
  mutate(id = paste0(pmin(feature_1, name), ' - ', pmax(feature_1, name))) %>%
  dplyr::filter(type == 'agreement' & statparam == 'stat') %>%
  group_by(dataset, id, resource) %>%
  summarise(similarity_average = mean(value), similarity_sem =  sd(value)/sqrt(n())) %>%
  ungroup() %>%
  separate(id, into = c('pipeline_a', 'pipeline_b'), sep = ' - ') %>%
  dplyr::filter(pipeline_a != pipeline_b) %>%
  mutate(pipeline_a = factor(pipeline_a, levels = sort(unique(pipeline_a))),
         pipeline_b = factor(pipeline_b, levels = sort(unique(pipeline_b))))

p <- ggplot(similarity_df, aes(x = pipeline_a, y = pipeline_b, size = -similarity_sem, fill = similarity_average)) +
  geom_point(pch = 21) +
  scale_fill_viridis_c(option = "plasma", limits = c(0, 1)) +
  facet_grid( cols = vars(resource), rows = vars(dataset)) +
  xlab('Pipeline A') + ylab('Pipeline B') +
  guides(x = guide_axis(angle = 60)) +
  theme_cowplot() 

# values for text
spearman_df %>% 
  mutate(unfil_limma_a = ifelse(pipeline_a == 'unfiltered-vsn+limma', 'Yes', 'No')) %>%
  mutate(unfil_limma_b = ifelse(pipeline_b == 'unfiltered-vsn+limma', 'Yes', 'No')) %>%
  mutate(unfil_limma = ifelse(unfil_limma_a == 'Yes' | unfil_limma_b == 'Yes', 'Yes', 'No')) %>%
  group_by(unfil_limma, resource) %>%
  summarise(mean(similarity_average))

similarity_toplot %>% 
  mutate(unfiltered_limma_a = ifelse(grepl('unfiltered', pipeline_a) & grepl('limma', pipeline_a), 'Yes', 'No')) %>%
  mutate(unfiltered_limma_b = ifelse(grepl('unfiltered', pipeline_b) & grepl('limma', pipeline_a), 'Yes', 'No')) %>%
  mutate(unfiltered_limma = ifelse(unfiltered_limma_a == 'Yes' | unfiltered_limma_b == 'Yes', 'Yes', 'No')) %>%
  group_by(unfiltered_limma, resource) %>%
  summarise(mean(similarity_average))

similarity_toplot %>% 
  mutate(unfil_deseq_edger_a = ifelse(grepl('unfiltered', pipeline_a) & (grepl('edger', pipeline_a) | grepl('deseq2', pipeline_a)), 'Yes', 'No')) %>%
  mutate(unfil_deseq_edger_b =  ifelse(grepl('unfiltered', pipeline_b) & (grepl('edger', pipeline_b) | grepl('deseq2', pipeline_b)), 'Yes', 'No')) %>%
  mutate(unfil_deseq_edger = ifelse(unfil_deseq_edger_a == 'Yes' & unfil_deseq_edger_b == 'Yes', 'Yes', 'No')) %>%
  group_by(unfil_deseq_edger, resource) %>%
  summarise(mean(spearman_average))

