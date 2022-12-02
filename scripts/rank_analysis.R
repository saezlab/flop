library(tidyverse)
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data <- readRDS("./results/merged_data.rds") %>%
  rename(obs_id = `...1`)

cor_results <- data %>%
  group_by(statparam, resource, cell_line, treatment) %>%
  group_split() %>%
  purrr::map(., function(x) {
    
    to_cor <- x %>%
      select(pipeline, scores, items) %>%
      pivot_wider(names_from = pipeline, values_from = scores) %>%
      column_to_rownames(var = "items") %>%
      as.matrix()
    
    cor_results <- cor(to_cor, method = "spearman") %>%
      as.data.frame() %>%
      rownames_to_column(var = "feature_1") %>%
      pivot_longer(-feature_1) %>%
      mutate(statparam = unique(x$statparam), cell_line = unique(x$cell_line),
             resource = unique(x$resource), treatment = unique(x$treatment))
    
    return(cor_results)
    
  }) %>%
  bind_rows()


pdf("test.pdf", width = 32, height = 10)
cor_results %>%
  subset(feature_1 != name) %>%
  mutate(id = paste0(feature_1, "__", name)) %>%
  ggplot(aes(x = id, y = value, fill = cell_line)) +
  facet_grid(cols = vars(resource), rows=vars(statparam)) +
  guides(x = guide_axis(angle = 60)) +
  geom_boxplot()
dev.off()