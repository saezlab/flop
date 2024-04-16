library(tidyverse)
library(limma)
library(ggrepel)

args <- commandArgs(trailingOnly = FALSE)
input_config <- args[grep("-m",args)+1] %>% strsplit(., split = ' ') %>% unlist()
input_f <- input_config[1]
biocontext <- input_config[2]
resource <- input_config[3]
statparam <- input_config[4] # logFC or t-value
print(input_config)


# try-except block
try(subset <- input_config[5], silent=TRUE)


# read input file
df <- read_tsv(input_f)

# pivot to create MDS-compatible table
filtered_df <- df %>% 
  dplyr::filter(bio_context == !!biocontext, resource == !!resource, statparam == !!statparam) 

distinct_subsets <- unique(filtered_df$subset)

if(length(distinct_subsets) > 1 & exists('subset')){
  filtered_df <- filtered_df %>% 
    dplyr::filter(subset == subset)
} else if(length(distinct_subsets) > 1 & !exists('subset')){
  print('Multiple subsets found in the input file. The resulting plot will combine signal from different subsets.')
  filtered_df <- filtered_df %>%
  group_by(obs_id) %>%
  select(-subset) %>%
  summarise(across(where(is.numeric), mean))
}
  
toplot <- filtered_df %>%
  mutate(id = paste0(status, '__', pipeline)) %>%
  dplyr::select(id, items, act) %>%
  pivot_wider(names_from = id, values_from = act) %>%
  column_to_rownames(var = 'items')

# compute MDS using limma
mds <- plotMDS(toplot,plot = FALSE)
mds_toplot <- tibble(
  axis_1 = mds$x,
  axis_2 = mds$y,
  label = colnames(mds$distance.matrix.squared)
)
x_axis_label <- paste0('Leading ', statparam, ' dim 1 (', round(mds$var.explained[1]*100,2),'%)')
y_axis_label <- paste0('Leading ', statparam, ' dim 2 (', round(mds$var.explained[2]*100,2),'%)')

# create plot
mds_plot <- ggplot(mds_toplot, aes(x = axis_1, y = axis_2, label = label)) +
  geom_point() +
  geom_text_repel() +
  labs(x = x_axis_label, y = y_axis_label) +
  ggtitle(paste0('Selected contrast:\n', biocontext)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# save plot
output_filename <- paste0(biocontext, '_', resource, '_', statparam, '_mds_plot.png')
ggsave(output_filename, mds_plot, width = 7, height = 7, dpi = 300)