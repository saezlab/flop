#Packages
library(tidyverse)

ID_gen <- function(n = 10) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

n=5
random_ident <- ID_gen(n)

set.seed(1)

countfile <- list.files(path = "./unproc_data/CCLE/", pattern = '*.gct', full.names = T)

data <- read_tsv(countfile, name_repair = make.names)  %>% dplyr::rename('gene_symbol'='Description') %>% dplyr::select(-Name) 
data_samples <- colnames(data)[-1]

metadata <- data_samples %>% 
  as_tibble() %>%
  separate(value, c('cell_line', 'group'), sep='_', extra='merge', remove=F) %>%
  dplyr::rename('sample_ID'='value')

thresh_meta <- metadata %>%
  group_by(group) %>%
  filter(n()>=20)

tissues <- thresh_meta %>% 
  dplyr::select(group) %>%
  distinct() %>%
  pull()


for (i in random_ident){
  dir.create(paste("./data/CCLE_", i, "/", sep=''))
  
  sampled_meta <- thresh_meta %>% 
    group_by(group) %>%
    print() %>% 
    slice_sample(n=20)
  
  sampled_data <- data %>% 
    dplyr::select(gene_symbol, sampled_meta$sample_ID)
  
  comparisons <- character()
  for (tissue1 in tissues){
    for (tissue2 in tissues){
      if (tissue1==tissue2){
        break
      }
      else {
        contrast <- paste0(sort(c(tissue1, tissue2))[1],'_v_', sort(c(tissue1, tissue2))[2])
        if (contrast %in% comparisons){
          break
        }else{
          comparisons <- c(comparisons, contrast)
          sel_meta <- sampled_meta %>% filter(group == tissue1 | group == tissue2)
          contr_samples <- sel_meta %>% dplyr::select(sample_ID) %>% pull()
          sel_counts <- sampled_data %>% dplyr::select(gene_symbol,!!contr_samples) %>% group_by(gene_symbol) %>% summarise_all(sum)
          
          count_name <- paste("./data/CCLE_", i, "/", paste(contrast, 'countdata.tsv', sep = '__'), sep='')
          meta_name <- paste("./data/CCLE_", i, "/", paste(contrast, 'metadata.tsv', sep='__'), sep='')
          
          write_tsv(sel_meta, meta_name)
          write_tsv(sel_counts, count_name)
        }
      }
    }
  }
}
  


#metadata %>% 
#  group_by(group) %>% 
#  count() %>% 
#  ggplot(.,aes(x=group,y=n, fill=n<20)) + 
#  geom_col() +
#  scale_fill_discrete('blue', 'green') +
#  ggtitle('Threshold = 20') +
#  theme(axis.text.x = element_text(angle=-60, vjust = 0, hjust=0), plot.margin = margin(r=30)) +
#  geom_hline(aes(yintercept=20, color='red'))