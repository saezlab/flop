library(tidyverse)
library(stats)
library(cowplot)
library(egg)

args <- commandArgs(trailingOnly = FALSE)
datafile <- args[grep("--file",args)+1][2]
merged_data <- read_tsv(datafile)

bio_context <- merged_data %>% distinct(bio_context) %>% pull()
pipelines <- merged_data %>% distinct(pipeline) %>% pull()
items <- merged_data %>% distinct(items) %>% pull()

## Rank Analysis: correlation
cor_results <- merged_data %>%
  group_by(statparam, resource, bio_context) %>%
  group_split() %>%
  print() %>% 
  purrr::map(., function(x) {
    to_cor <- x %>%
      select(pipeline, scores, items) %>%
      pivot_wider(names_from = pipeline, values_from = scores) %>%
      column_to_rownames(var = "items")
    print(head(to_cor))
    
    cor_results <- cor(to_cor, method = "spearman") %>%
      as.data.frame() %>%
      rownames_to_column(var = "feature_1") %>%
      pivot_longer(-feature_1) %>%
      mutate(statparam = unique(x$statparam), bio_context = unique(x$bio_context),
             resource = unique(x$resource))
    
    return(cor_results)
  }) %>%
  bind_rows()

cor_filt_results <- cor_results %>%
  subset(feature_1 != name) %>%
  rowwise() %>%
  mutate(id = paste0(sort(c(feature_1, name))[1], " - ", sort(c(feature_1, name))[2])) %>%
  distinct(id, statparam, bio_context, resource, .keep_all = T)

#Plot figure rank preservation
total_plot <- ggplot(cor_filt_results,aes(x = id, y = value, fill=id)) +
  facet_grid(cols = vars(statparam), rows=vars(resource), labeller=labeller(
    resource = c(dorothea = 'DoRothEA', msigdb_hallmarks = 'MSigDB Hallmarks', progeny='PROGENy'),
    statparam = c(logFC = 'LogFC', stat = 'T values')
  )) +
  guides(x = guide_axis(angle = 60)) +
  geom_boxplot() +
  theme_cowplot()+
  ylab('Score')+
  theme(
    legend.position = 'null',
    axis.title.x = element_blank(),
    plot.margin=margin(10,10,10,20)
        )
total_plot
save_plot(filename = "total.png", plot = total_plot, device='png', dpi=300, base_height=10, base_width =10)

#Correlation split
#By norm
norm_cor_results <- merged_data %>%
  group_by(statparam, resource, bio_context) %>%
  group_split() %>%
  purrr::map(., function(x) {
    
    to_cor <- x %>%
      select(norm, scores, items) %>%
      pivot_wider(names_from = norm, values_from = scores, values_fn = {mean}) %>%
      column_to_rownames(var = "items") %>%
      as.matrix()
    
    cor_results <- cor(to_cor, method = "spearman") %>%
      as.data.frame() %>%
      rownames_to_column(var = "feature_1") %>%
      pivot_longer(-feature_1) %>%
      mutate(statparam = unique(x$statparam), bio_context = unique(x$bio_context),
             resource = unique(x$resource))
    
    return(cor_results)
    
  }) %>%
  bind_rows()

norm_cor_filt_results <- norm_cor_results %>%
  subset(feature_1 != name) %>%
  rowwise() %>%
  mutate(id = paste0(sort(c(feature_1, name))[1], " - ", sort(c(feature_1, name))[2])) %>%
  distinct(id, statparam, bio_context, resource, .keep_all = T)

norm_corplot <- norm_cor_filt_results %>%
  subset(feature_1 != name) %>%
  mutate(id = paste0(feature_1, " - ", name)) %>%
  ggplot(aes(x = id, y = value, fill = name))+
  geom_boxplot() +
  ylim(-0.7,1.2)+
  scale_fill_manual(values=c( "grey", "grey", "#69b3a2"))+
  scale_alpha_manual(values=c(0.1,0.1,1))+
  facet_grid(cols = vars(statparam), rows=vars(resource)) +
  guides(x = guide_axis(angle = 60)) +
  theme_cowplot()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        strip.background.y = element_blank(),
        strip.text.y = element_blank()
        )

#By diffexp
diffexp_cor_results <- merged_data %>%
  group_by(statparam, resource, bio_context) %>%
  group_split() %>%
  purrr::map(., function(x) {
    
    to_cor <- x %>%
      select(diffexp, scores, items) %>%
      pivot_wider(names_from = diffexp, values_from = scores, values_fn = {mean}) %>%
      column_to_rownames(var = "items") %>%
      as.matrix()
    
    cor_results <- cor(to_cor, method = "spearman") %>%
      as.data.frame() %>%
      rownames_to_column(var = "feature_1") %>%
      pivot_longer(-feature_1) %>%
      mutate(statparam = unique(x$statparam), bio_context = unique(x$bio_context),
             resource = unique(x$resource))
    
    return(cor_results)
    
  }) %>%
  bind_rows()

  diffexp_cor_filt_results <- diffexp_cor_results %>%
  subset(feature_1 != name) %>%
  rowwise() %>%
  mutate(id = paste0(sort(c(feature_1, name))[1], "__", sort(c(feature_1, name))[2])) %>%
  distinct(id, statparam, bio_context, resource, .keep_all = T)

diffexp_corplot <- diffexp_cor_filt_results %>%
  subset(feature_1 != name) %>%
  mutate(id = paste0(feature_1, " - ", name)) %>%
  ggplot(aes(x = id, y = value, fill=name))+
  scale_fill_manual(values=rep('grey',2))+
  geom_boxplot() +
  ylim(-0.7,1.2)+
  facet_grid(cols = vars(statparam), rows=vars(resource)) +
  guides(x = guide_axis(angle = 60)) +
  theme_cowplot()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()
        )

aligned_plots <- egg::ggarrange(norm_corplot, diffexp_corplot, ncol = 2, widths = c(2,1))
final <- ggdraw() + draw_plot(aligned_plots, 0, 0, 1, 1)
save_plot(filename = 'split_methods_plots.png', plot = final, device='png', dpi=300, base_height=9, base_width =12)


#Comparison of outliers
max_corr <- cor_filt_results %>% 
  filter(statparam=='stat') %>%
  filter(resource=='msigdb_hallmarks') %>%
  .[which.max(.$value),] 

min_corr <- cor_filt_results %>% 
  filter(statparam=='stat') %>%
  filter(resource=='msigdb_hallmarks') %>%
  .[which.min(.$value),] 
  
max_data <- merged_data %>%
  filter(resource=='msigdb_hallmarks') %>%
  filter(pipeline== max_corr$feature_1|pipeline==max_corr$name, bio_context== max_corr$bio_context) %>%
  mutate('plot'='max')

min_data <- merged_data %>%
  filter(resource=='msigdb_hallmarks') %>%
  filter(pipeline== min_corr$feature_1|pipeline==min_corr$name, bio_context== max_corr$bio_context) %>%
  mutate('plot'='min')

plot_data <- rbind(max_data, min_data)

rank_comp<-ggplot(plot_data) +
  geom_boxplot(aes(x = pipeline, y = scores, fill=pipeline))+
  geom_point(aes(x = pipeline, y = scores)) +
  geom_line(aes(x  = pipeline, y = scores, group = items), alpha=0.3) +
  theme_cowplot() +
  theme(axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2,'lines'),
        axis.text.x = element_text(angle=60, vjust=0.9, hjust=1),
        legend.position = 'none'
  ) +
  facet_grid(~plot, scales = 'free')
labeled_rank_comp <- ggdraw() + draw_plot(rank_comp, 0, 0, 1, 1) + draw_plot_label('A', 0, 1) +
    draw_plot_label('B', 0.52, 1)
save_plot('rank_comp.png', labeled_rank_comp, device='png', dpi=300, base_height=6, base_width =6)
