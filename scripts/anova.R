library(tidyverse)
setwd('C:/Users/victo/OneDrive - Universidad Politécnica de Madrid/Documentos/1º Master/Internship/flop_benchmark/scripts/')
data <- read_tsv(file = 'logFC__cons__test.tsv')

long_data <- data %>%
  pivot_longer(-X, names_to = 'pathway', values_to = 'scores') %>% 
  separate(X, into = c("cl_treat", "norm", "diffexp", "decouplermethod"), sep = "__") %>%
  separate(cl_treat, into = c("cell_line", "treatment")) %>%
  write.table(., file='long__logFC__cons__test.tsv')

# discuss later the parametric assumption in consensus scores
#res <- kruskal.test(scores ~ norm, data = long_data)

ggplot(data = long_data, aes(x = scores)) +
  geom_histogram() +
  facet_wrap(facets =vars(pathway))


data %>%
  column_to_rownames("X") %>%
  as.matrix() %>%
  
  dist() %>%
  hclust() %>%
  plot()


library(nlme)


aov_res <- aov(scores ~ norm + pathway + cell_line + treatment, data = long_data) 


nlm_res <- 


summary(aov_res)
TukeyHSD(aov_res)

  

