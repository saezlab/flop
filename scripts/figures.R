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
full_df <- list.files('./flop_results/funcomics/fullmerged/', full.names = TRUE) %>% lapply( read_tsv) %>% bind_rows() #%>%
  # dplyr::filter(subset %in% c("Liu15", "Pepin19", "Schiano17", "Spurell19", "Yang14") & main_dataset == 'reheat')
similarity_df <- list.files('flop_results/funcomics/overlap/', full.names = TRUE) %>% lapply(read_tsv) %>% bind_rows() # %>%
  # dplyr::filter(subset %in% c("Liu15", "Pepin19", "Schiano17", "Spurell19", "Yang14") & main_dataset == 'reheat')
correlation_df <-  list.files('flop_results/funcomics/rank/', full.names = TRUE) %>% lapply( read_tsv) %>% bind_rows() # %>%
  # dplyr::filter(subset %in% c("Liu15", "Pepin19", "Schiano17", "Spurell19", "Yang14") & main_dataset == 'reheat')


p_cutoff <- 0.05
lfc_cutoff <- log2(2)

#### fig 2 #### 
plot_reheat_detailed <- function(int_pipelines, int_dataset) {

  # volcanos
  toplot <- reheat_de %>%
    mutate(pipeline = paste0(filtering, '+', str_replace(de , 'NA\\+', ''))) %>% 
    dplyr::filter(pipeline %in% !!int_pipelines & study == !!int_dataset) %>%
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
    dplyr::filter(pipeline %in% int_pipelines & main_dataset == 'reheat' & subset == int_dataset & statparam == 'stat') %>%
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
ggsave('flop_results/plots/supp2.png', supp1,  width = 15, height = 10, dpi = 300)
supp2 <- plot_reheat_detailed(c('unfiltered+edger', 'filtered+deseq2'), 'Yang14')
ggsave('flop_results/plots/supp3.png', supp2,  width = 15, height = 10, dpi = 300)

custom_sort <- function(x) {
  
  # Split the vector based on "filtered" and "unfiltered". In case none match, get all
  filtered <- x[grep("^filtered", x)]
  unfiltered <- x[grep("unfiltered", x)]

  ifelse(length(filtered) == 0, filtered <- x, filtered <- filtered)
  
  # Within the filtered group, sort strings containing "NA" first
  filtered_nas <- sort(filtered[grep("NA+", filtered)])
  filtered_non_nas <- sort(filtered[!filtered %in% filtered_nas])
  
  # Within the unfiltered group, sort strings containing "NA" first
  unfiltered_nas <- sort(unfiltered[grep("NA+", unfiltered)])
  unfiltered_non_nas <- sort(unfiltered[!unfiltered %in% unfiltered_nas])


  # Combine all sorted vectors
  c(filtered_nas, filtered_non_nas, unfiltered_nas, unfiltered_non_nas)
}


custom_sort(c('filtered-NA+', 'unfiltered-log1quant+limma'))

#### fig 3 ####
correlation_toplot <- correlation_df %>%
  rowwise() %>%
  mutate(id = paste0(custom_sort(c(feature_1, name))[1], ' - ', custom_sort(c(feature_1, name))[2])) %>%
  dplyr::filter(type == 'correlation' & statparam == 'stat') %>%
  group_by(id, resource, main_dataset) %>%
  summarise(spearman_average = mean(value), spearman_sem =  sd(value)/sqrt(n())) %>%
  ungroup() %>%
  separate(id, into = c('pipeline_a', 'pipeline_b'), sep = ' - ') %>%
  dplyr::filter(pipeline_a != pipeline_b) %>%
  mutate(pipeline_a = factor(pipeline_a, levels = custom_sort(unique(pipeline_a))),
         pipeline_b = factor(pipeline_b, levels = custom_sort(unique(pipeline_b))))


# dotplot
p <- ggplot(correlation_toplot, aes(x = pipeline_a, y = pipeline_b, size = -spearman_sem, fill = spearman_average)) +
  geom_point(pch = 21) +
  scale_fill_viridis_c(option = "plasma", limits = c(0, 1)) +
  facet_wrap(vars(resource)) +
  xlab('Pipeline A') + ylab('Pipeline B') +
  guides(x = guide_axis(angle = 60)) +
  theme_cowplot() 

ggsave('flop_results/plots/fig3.png', p, width = 12, height = 10, dpi = 300)



# boxplot
correlation_toplot_boxplot <- correlation_df %>%
  rowwise() %>%
  mutate(id = paste0(custom_sort(c(feature_1, name))[1], ' - ', custom_sort(c(feature_1, name))[2])) %>%
  dplyr::filter(type == 'correlation' & statparam == 'stat') %>%
  separate(id, into = c('pipeline_a', 'pipeline_b'), sep = ' - ') %>%
  dplyr::filter(pipeline_a != pipeline_b) %>%
  mutate(pipeline_a = factor(pipeline_a, levels = custom_sort(unique(pipeline_a))),
         pipeline_b = factor(pipeline_b, levels = custom_sort(unique(pipeline_b)))) %>%
  select(-feature_1, -name)

correlation_toplot_dupl <- correlation_toplot_boxplot %>%
  rename(pipeline_a = pipeline_b, pipeline_b = pipeline_a) %>%
  bind_rows(correlation_toplot_boxplot)

dataset.labs <- c('CCLE', 'PANACEA', 'ReHeaT')
names(dataset.labs) <- c('CCLE', 'GSE186341', 'reheat')

resource_comparisons <- list( c('DE', 'dorothea'), c('DE', 'msigdb_hallmarks'), c('DE', 'progeny'))
boxplot_fig2 <- correlation_toplot_dupl %>%
  drop_na(value) %>%
  ggplot(aes(y = pipeline_b, x = value, fill = resource)) +
  geom_boxplot(outlier.alpha = 0.1, outlier.size=0.5, lwd = 0.2) +
  facet_wrap(~main_dataset, ncol=1, labeller = labeller(main_dataset = dataset.labs)) +
  # stat_compare_means(aes(group = resource), method = 'wilcox.test', ref.group = 'DE', label = "p.signif", label.y = 0) +
  ylab('Pipeline') + xlab('Spearman correlation') + 
  theme_cowplot() +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(size=10),
        text = element_text(family='Calibri', size=10),
        axis.text.y = element_text(size = 10))
        
ggsave('plots/fig2_new.png', boxplot_fig2, width = 17, height = 20, dpi = 300, units = 'cm')



fig2_theme <- theme_cowplot() + theme(
  axis.text.x = element_text(size=8, family='Arial'),
  axis.text.y = element_text(size=8, family='Arial'),
  axis.title = element_text(family='Arial', size=8),
  legend.text = element_text(size=8, family='Arial'),
  legend.title = element_text(size=8, family='Arial')
)

facet_b_fig2 <- correlation_toplot_dupl %>% 
  filter(resource != 'DE') %>%
  mutate(Filtered = ifelse(!grepl('unfiltered', pipeline_a), 'Yes', 'No')) %>%
  ggplot(aes(y = spearman_average, x = main_dataset, fill = Filtered)) +
  geom_boxplot(outlier.alpha = 0.1, outlier.size=0.5, lwd = 0.2) +
  stat_compare_means(aes(group = Filtered), method = 'wilcox.test', label = "p.signif", label.y=1.05, size=3) +
  scale_fill_brewer(palette="Dark2") +
  ylab('') + xlab('Dataset') +
  scale_x_discrete(labels = c('CCLE', 'PANACEA', 'reheat')) +
  ylim(0.5, 1.1) +
  fig2_theme +
  theme(axis.text.y = element_blank())

facet_a_fig2 <- correlation_toplot_dupl %>% 
  mutate(Space = ifelse(resource == 'DE', 'DGE sp.', 'Functional sp.')) %>%
  ggplot(aes(y = spearman_average, x = main_dataset, fill = Space)) +
  geom_boxplot(outlier.alpha = 0.1, outlier.size=0.5, lwd = 0.2) +
  stat_compare_means(aes(group = Space), method = 'wilcox.test', label = "p.signif", label.y = 1.05, size=3) +
  ylab('Spearman correlation') + xlab('Dataset') +
  scale_x_discrete(labels = c('CCLE', 'PANACEA', 'reheat')) +
  ylim(0.5, 1.1) +
  fig2_theme 


facet_c_fig2_data <- correlation_df %>%
# filter rows that have filtered in name
  dplyr::filter(!grepl('unfiltered', feature_1) & !grepl('unfiltered', name)) %>%
  rowwise() %>%
  mutate(id = paste0(custom_sort(c(feature_1, name))[1], ' - ', custom_sort(c(feature_1, name))[2])) %>%
  dplyr::filter(type == 'correlation' & statparam == 'stat') %>%
  group_by(id) %>%
  summarise(spearman_average = mean(value), spearman_sd =  sd(value)) %>%
  ungroup() %>%
  separate(id, into = c('pipeline_a', 'pipeline_b'), sep = ' - ') %>%
  dplyr::filter(pipeline_a != pipeline_b) %>%
  mutate(pipeline_a = str_replace(pipeline_a, '.*-', ''),
         pipeline_b = str_replace(pipeline_b, '.*-', ''))  %>%
  mutate(pipeline_a = factor(pipeline_a, levels = custom_sort(unique(pipeline_a))),
         pipeline_b = factor(pipeline_b, levels = custom_sort(unique(pipeline_b))))


facet_c_fig2 <- facet_c_fig2_data %>%
  ggplot(aes(y = pipeline_b, x = pipeline_a, size = spearman_sd, fill = spearman_average)) +
  geom_point(pch = 21) +
  scale_fill_viridis_c(option = "plasma", limits = c(0.4, 1), name='Avg value') +
  scale_size(trans = 'reverse', name='Sd', limits = c(0.2,0)) +
  xlab('') + ylab('') +
  guides(x = guide_axis(angle = 60)) +
  fig2_theme

facet_d_fig2_data <- similarity_df %>%
# filter rows that have filtered in name
  dplyr::filter(!grepl('unfiltered', feature_1) & !grepl('unfiltered', name)) %>%
  rowwise() %>%
  mutate(id = paste0(custom_sort(c(feature_1, name))[1], ' - ', custom_sort(c(feature_1, name))[2])) %>%
  dplyr::filter(type == 'agreement' & statparam == 'stat') %>%
  group_by(id) %>%
  summarise(similarity_average = mean(value), similarity_sd =  sd(value)) %>%
  ungroup() %>%
  separate(id, into = c('pipeline_a', 'pipeline_b'), sep = ' - ') %>%
  dplyr::filter(pipeline_a != pipeline_b) %>%
  mutate(pipeline_a = str_replace(pipeline_a, '.*-', ''),
         pipeline_b = str_replace(pipeline_b, '.*-', ''))  %>%
  mutate(pipeline_a = factor(pipeline_a, levels = custom_sort(unique(pipeline_a))),
         pipeline_b = factor(pipeline_b, levels = custom_sort(unique(pipeline_b))))


facet_d_fig2 <- facet_d_fig2_data %>%
  ggplot(aes(y = pipeline_b, x = pipeline_a, size = similarity_sd, fill = similarity_average)) +
  geom_point(pch = 21) +
  scale_fill_viridis_c(option = "plasma", limits = c(0.4, 1)) +
  scale_size(trans = 'reverse', limits = c(0.2, 0)) +
  xlab('') + ylab('') +
  guides(x = guide_axis(angle = 60)) +
  fig2_theme +
  theme(axis.text.y = element_blank())

panel_cd__fig2_legend <- cowplot::get_legend(facet_c_fig2) %>% ggpubr::as_ggplot()
panel_a__fig2_legend <- cowplot::get_legend(facet_a_fig2) %>% ggpubr::as_ggplot()
panel_b__fig2_legend <- cowplot::get_legend(facet_b_fig2) %>% ggpubr::as_ggplot()

facet_a_fig2 <- facet_a_fig2 + theme(legend.position = 'none')
facet_b_fig2 <- facet_b_fig2 + theme(legend.position = 'none')
facet_c_fig2 <- facet_c_fig2 + theme(legend.position = 'none')
facet_d_fig2 <- facet_d_fig2 + theme(legend.position = 'none')

fig_2_new <- egg::ggarrange(facet_a_fig2, facet_b_fig2, facet_c_fig2, facet_d_fig2, ncol = 2, nrow = 2, padding = unit(4, "cm"), labels=c('A', 'B', 'C', 'D'))
fig_2_legends <- egg::ggarrange(panel_a__fig2_legend, panel_b__fig2_legend, panel_cd__fig2_legend, ncol = 1, nrow = 3)
ggsave('plots/supfig2_new_facets_filt.svg', fig_2_new, width = 16, height = 16, dpi = 300, units = 'cm')
ggsave('plots/supfig2_new_facets_legend_filt.svg', fig_2_legends, height = 25, dpi = 300, units = 'cm')

t_test_adj <- correlation_toplot_dupl %>%
  drop_na(value) %>%
  mutate(is_de = ifelse(resource == 'DE', 'Yes', 'No')) %>%
  group_by(main_dataset, pipeline_a) %>%
  group_split() %>%
  purrr::map(., function(x){
    padj_corr <- x %>% compare_means(value ~ is_de, p.adjust.method = "BH", method='wilcox.test', data = .,  alternative = 'greater') %>% mutate(main_dataset = unique(x$main_dataset), pipeline_a = unique(x$pipeline_a))
    return(padj_corr)
  }) %>% bind_rows() %>% select(main_dataset, pipeline_a, p.format) %>%
  pivot_wider(names_from = main_dataset, values_from = p.format)


# values for text
correlation_toplot %>% 
  mutate(unfil_limma_a = ifelse(grepl('unfiltered', pipeline_a) & grepl('limma', pipeline_a), 'Yes', 'No')) %>%
  mutate(unfil_limma_b = ifelse(grepl('unfiltered', pipeline_b) & grepl('limma', pipeline_b), 'Yes', 'No')) %>%
  mutate(unfil_limma = ifelse(unfil_limma_a == 'Yes' | unfil_limma_b == 'Yes', 'Yes', 'No')) %>%
  group_by(unfil_limma, resource, main_dataset) %>%
  summarise(mean(spearman_average))

correlation_toplot %>% 
  mutate(unfil_limma_a = ifelse(grepl('unfiltered', pipeline_a) & grepl('limma', pipeline_a), 'Yes', 'No')) %>%
  mutate(unfil_limma_b = ifelse(grepl('unfiltered', pipeline_b) & grepl('limma', pipeline_b), 'Yes', 'No')) %>%
  mutate(unfil_limma = ifelse(unfil_limma_a == 'Yes' | unfil_limma_b == 'Yes', 'Yes', 'No')) %>%
  filter(unfil_limma == 'Yes') %>%
  mutate(is_de = ifelse(resource == 'DE', 'Yes', 'No')) %>%
  group_by(main_dataset) %>%
  group_split() %>%
  purrr::map(., function(x){
    x %>% distinct(main_dataset) %>% pull(main_dataset) %>% print()
    x %>% wilcox.test(spearman_average ~ is_de, data = ., alternative = "less", p.adjust.methods = "BH")
  })

correlation_toplot_dupl %>%
  filter(resource != 'DE', main_dataset=='GSE186341') %>%
  group_by(pipeline_a) %>%
  mutate(filtered = ifelse(!grepl('unfiltered', pipeline_a), 'Yes', 'No')) %>%
  group_by(filtered) %>%
  summarise(mean = mean(value), sd = sd(value), max = max(value), min = min(value)) %>% arrange(desc(sd)) %>% print(n=50)

correlation_toplot_dupl %>%
  filter(resource != 'DE', main_dataset=='reheat') %>%
  group_by(pipeline_a) %>%
  summarise(mean = mean(value), sd = sd(value), max = max(value), min = min(value)) %>% arrange(desc(sd)) %>% print(n=50)

correlation_toplot_dupl %>%
  filter(resource != 'DE', main_dataset=='CCLE') %>%
  group_by(pipeline_a) %>%
  summarise(mean = mean(value), sd = sd(value), max = max(value), min = min(value)) %>% arrange(desc(mean)) %>% print(n=50)

correlation_toplot_dupl %>%
  filter(resource != 'DE', main_dataset=='GSE186341') %>%
  group_by(pipeline_a) %>%
  summarise(mean = mean(value), sd = sd(value), max = max(value), min = min(value)) %>% arrange(desc(mean)) %>% print(n=50)

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
  group_by(main_dataset) %>%
  group_split() %>%
  purrr::map(., function(x){
    x %>% distinct(main_dataset) %>% pull(main_dataset) %>% print()
    x %>% wilcox.test(spearman_average ~ is_de, data = ., alternative = "less", p.adjust.methods = "BH")
  })

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
  rowwise() %>%
  mutate(id = paste0(custom_sort(c(feature_1, name))[1], ' - ', custom_sort(c(feature_1, name))[2])) %>%
  dplyr::filter(type == 'agreement' & statparam == 'stat') %>%
  dplyr::filter(subset %in% c("Liu15", "Pepin19", "Schiano17", "Spurell19", "Yang14") & main_dataset == 'reheat') %>%
  group_by(id, resource) %>%
  summarise(similarity_average = mean(value), similarity_sem =  sd(value)/sqrt(n())) %>%
  ungroup() %>%
  separate(id, into = c('pipeline_a', 'pipeline_b'), sep = ' - ') %>%
  dplyr::filter(pipeline_a != pipeline_b) %>%
  mutate(pipeline_a = factor(pipeline_a, levels = custom_sort(unique(pipeline_a))),
         pipeline_b = factor(pipeline_b, levels = custom_sort(unique(pipeline_b))))

p <- ggplot(similarity_toplot, aes(x = pipeline_a, y = pipeline_b, size = -similarity_sem, fill = similarity_average)) +
  geom_point(pch = 21) +
  scale_fill_viridis_c(option = "plasma", limits = c(0, 1)) +
  facet_wrap( vars(resource)) +
  xlab('Pipeline A') + ylab('Pipeline B') +
  guides(x = guide_axis(angle = 60)) +
  theme_cowplot() 

ggsave('flop_results/plots/fig4.png', p, width = 12, height = 10, dpi = 300)


similarity_toplot_boxplot <- similarity_df %>%
  rowwise() %>%
  mutate(id = paste0(custom_sort(c(feature_1, name))[1], ' - ', custom_sort(c(feature_1, name))[2])) %>%
  dplyr::filter(type == 'agreement' & statparam == 'stat') %>%
  separate(id, into = c('pipeline_a', 'pipeline_b'), sep = ' - ') %>%
  dplyr::filter(pipeline_a != pipeline_b) %>%
  mutate(pipeline_a = factor(pipeline_a, levels = custom_sort(unique(pipeline_a))),
         pipeline_b = factor(pipeline_b, levels = custom_sort(unique(pipeline_b))))


similarity_toplot_dupl <- similarity_toplot_boxplot %>%
  rename(pipeline_a = pipeline_b, pipeline_b = pipeline_a) %>%
  bind_rows(similarity_toplot_boxplot)


boxplot_fig_similarity <- similarity_toplot_dupl %>%
  drop_na(value) %>%
  ggplot(aes(y = pipeline_b, x = value, fill = resource)) +
  geom_boxplot(outlier.alpha = 0.1, outlier.size=0.5, lwd = 0.2) +
  facet_wrap(~main_dataset, ncol=1, labeller = labeller(main_dataset = dataset.labs)) +
  # stat_compare_means(aes(group = resource), method = 'wilcox.test', ref.group = 'DE', label = "p.signif", label.y = 0) +
  ylab('') + xlab('Similarity score') + 
  theme_cowplot() +
  theme(legend.position = 'none',
        axis.text.x = element_text(size=10),
        text = element_text(family='Calibri', size=10),
        axis.text.y = element_blank())

fig_2_new <- egg::ggarrange(boxplot_fig2, boxplot_fig_similarity, ncol = 2, nrow = 1)
        
ggsave('plots/fig2_new.png', fig_2_new, width = 17, height = 20, dpi = 300, units = 'cm')

ggsave('plots/fig2_new_legend.svg', fig_2_new, width = 17, height = 20, dpi = 300, units = 'cm')

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
  summarise(mean(similarity_average))

similarity_toplot_dupl %>%
  filter(resource != 'DE', main_dataset=='reheat') %>%
  group_by(pipeline_a) %>%
  summarise(mean = mean(value), sd = sd(value), max = max(value), min = min(value)) %>% arrange(desc(mean)) %>% print(n=50)

similarity_toplot_dupl %>%
  filter(resource != 'DE', main_dataset=='CCLE') %>%
  group_by(pipeline_a) %>%
  summarise(mean = mean(value), sd = sd(value), max = max(value), min = min(value)) %>% arrange(desc(mean)) %>% print(n=50)

similarity_toplot_dupl %>%
  filter(resource != 'DE', main_dataset=='GSE186341') %>%
  group_by(pipeline_a) %>%
  summarise(mean = mean(value), sd = sd(value), max = max(value), min = min(value)) %>% arrange(desc(mean)) %>% print(n=50)

# sup about variance
sim_values <- similarity_df %>%
  mutate(id = paste0(pmin(feature_1, name), ' - ', pmax(feature_1, name))) %>%
  dplyr::filter(type == 'agreement' & statparam == 'stat') %>%
  dplyr::filter(!subset %in% c('vanHeesch19', 'Sweet18')) %>%
  transmute(id = paste(id, main_dataset, resource, sep = '___'), similarity = value)
  

cor_values <- correlation_df %>%
  mutate(id = paste0(pmin(feature_1, name), ' - ', pmax(feature_1, name))) %>%
  dplyr::filter(feature_1 != name) %>%
  dplyr::filter(type == 'correlation' & statparam == 'stat') %>%
  dplyr::filter(!subset %in% c('vanHeesch19', 'Sweet18')) %>%
  transmute(id = paste(id, main_dataset, resource, sep = '___'), correlation = value)

toplot <- full_join(cor_values, sim_values, by = 'id') %>%
  pivot_longer(-id)

p <- ggplot(toplot, aes(x = name, y = value)) +
  geom_line(aes(group = id), alpha = 0.05) +
  geom_sina(alpha = 0.5, pch = 21) +
  geom_boxplot(alpha = 0.5, pch = 21) +
  stat_compare_means(method = 'wilcox.test', paired = TRUE, comparisons = list(c('correlation', 'similarity'))) +
  theme_cowplot()

ggsave('flop_results/plots/supp4.png', p, width = 5, height = 5, dpi = 300)

#### fig 5 ####
spearman_df <- bind_rows(ccle = ccle_spearman, panacea = panacea_spearman, .id = 'dataset') %>%
  rowwise() %>%
  mutate(id = paste0(custom_sort(c(feature_1, name))[1], ' - ', custom_sort(c(feature_1, name))[2])) %>%
  dplyr::filter(type == 'correlation' & statparam == 'stat') %>%
  group_by(dataset, id, resource) %>%
  summarise(spearman_average = mean(value), spearman_sem =  sd(value)/sqrt(n())) %>%
  ungroup() %>%
  separate(id, into = c('pipeline_a', 'pipeline_b'), sep = ' - ') %>%
  dplyr::filter(pipeline_a != pipeline_b) %>%
  mutate(pipeline_a = factor(pipeline_a, levels = custom_sort(unique(pipeline_a))),
         pipeline_b = factor(pipeline_b, levels = custom_sort(unique(pipeline_b)))) 

p <- ggplot(spearman_df, aes(x = pipeline_a, y = pipeline_b, size = -spearman_sem, fill = spearman_average)) +
  geom_point(pch = 21) +
  scale_fill_viridis_c(option = "plasma", limits = c(0, 1)) +
  facet_grid( cols = vars(resource), rows = vars(dataset)) +
  xlab('Pipeline A') + ylab('Pipeline B') +
  guides(x = guide_axis(angle = 60)) +
  theme_cowplot() 

ggsave('flop_results/plots/fig5.png', p, width = 14, height = 8, dpi = 300)

similarity_df <- bind_rows(ccle = ccle_similarity, panacea = panacea_similarity, .id = 'dataset') %>%
  rowwise() %>%
  mutate(id = paste0(custom_sort(c(feature_1, name))[1], ' - ', custom_sort(c(feature_1, name))[2])) %>%
  dplyr::filter(type == 'agreement' & statparam == 'stat') %>%
  group_by(dataset, id, resource) %>%
  summarise(similarity_average = mean(value), similarity_sem =  sd(value)/sqrt(n())) %>%
  ungroup() %>%
  separate(id, into = c('pipeline_a', 'pipeline_b'), sep = ' - ') %>%
  dplyr::filter(pipeline_a != pipeline_b) %>%
  mutate(pipeline_a = factor(pipeline_a, levels = custom_sort(unique(pipeline_a))),
         pipeline_b = factor(pipeline_b, levels = custom_sort(unique(pipeline_b))))

p <- ggplot(similarity_df, aes(x = pipeline_a, y = pipeline_b, size = -similarity_sem, fill = similarity_average)) +
  geom_point(pch = 21) +
  scale_fill_viridis_c(option = "plasma", limits = c(0, 1)) +
  facet_grid( cols = vars(resource), rows = vars(dataset)) +
  xlab('Pipeline A') + ylab('Pipeline B') +
  guides(x = guide_axis(angle = 60)) +
  theme_cowplot() 

# # values for text
spearman_df %>% 
  mutate(unfil_limma_a = ifelse(pipeline_a == 'unfiltered-vsn+limma', 'Yes', 'No')) %>%
  mutate(unfil_limma_b = ifelse(pipeline_b == 'unfiltered-vsn+limma', 'Yes', 'No')) %>%
  mutate(unfil_limma = ifelse(unfil_limma_a == 'Yes' | unfil_limma_b == 'Yes', 'Yes', 'No')) %>%
  group_by(unfil_limma, resource) %>%
  summarise(mean(spearman_average))

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
  summarise(mean(similarity_average))

## Supp figure

plot_reheat_tscores <- function(int_pipelines, int_dataset) {
  hallmark <- read_tsv("./scripts/dc_resources/msigdb_hallmarks__source.tsv")
  pick_HM <- hallmark %>% filter(source %in% c("HALLMARK_ADIPOGENESIS",
                                             "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                                             "HALLMARK_HEDGEHOG_SIGNALING",
                                             "HALLMARK_INTERFERON_GAMMA_RESPONSE"))

  toplot <- reheat_de %>%
      mutate(pipeline = paste0(filtering, '+', str_replace(de , 'NA\\+', ''))) %>%
      dplyr::filter(pipeline %in% int_pipelines & study == int_dataset) %>%
      dplyr::select(pipeline, stat, ID)  %>%
      pivot_wider(names_from = pipeline, values_from = stat) 

  toplot2 <- toplot %>%
      mutate(filter_status = ifelse(is.na(`filtered+deseq2`),"filtered","unfiltered"),
            `filtered+deseq2` = ifelse(is.na(`filtered+deseq2`),0,`filtered+deseq2`),
            `unfiltered+edger` = ifelse(is.na(`unfiltered+edger`),0,`unfiltered+edger`))%>%
      left_join(pick_HM, by=join_by(ID==target)) %>%
      filter(ID %in% pick_HM$target)%>%
      mutate(source = gsub("HALLMARK_","",source)) 

  pw_const_genes_plot <- toplot2 %>% 
    
      ggplot( aes(x = !!sym(int_pipelines[1]), y = !!sym(int_pipelines[2]))) +
      geom_point(alpha = 0.5, pch = 21,aes(fill = !!sym(int_pipelines[1]))) +
      geom_point(data = filter(toplot2,filter_status=="filtered"), col="orange") +
      scale_fill_gradient2(high =  'red',low = 'blue', mid = 'grey', midpoint = 0) +
      geom_hline(yintercept = 0, linetype = 2, color="grey")+
      geom_vline(xintercept = 0, linetype = 2, color="grey")+
      stat_cor(method = 'spearman',label.x = -13,label.y = 16) +
      ggtitle('Spurell19, PW-constituting genes') +
      theme_cowplot() +
      theme(legend.position = 'none')+facet_wrap(~ source)



  box_comp_plot <- toplot2 %>% 
      ggplot(aes(filter_status,!!sym(int_pipelines[1]),)) + geom_point() +
      geom_boxplot() + facet_wrap(~source) + stat_compare_means() + theme_cowplot()+ 
      ggtitle("T-score, unfiltered + edger pipeline") + xlab("") + ylab(" T-statistics")


  p <- plot_grid(
      plot_grid(pw_const_genes_plot, box_comp_plot, rel_widths = c(1, 1), nrow = 1, labels = c('A', 'B')),
      nrow = 1
  )
}

p <- plot_reheat_tscores(c('unfiltered+edger', 'filtered+deseq2'), 'Spurell19')

ggsave('flop_results/plots/supp1.png', p,  width = 15, height = 10, dpi = 300)