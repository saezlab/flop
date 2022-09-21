#benchmarking
norm_results <- norm_func(filt_data, norm_types)
diffanal_results <- list()
for (i in norm_types){
  norm_data <- norm_results %>% filter(norm==i) %>% unnest() %>% select(-norm)
  diffanal_results[[i]]<- diff_anal(norm_data, metadata, diffanal_types)
}
diffanal_results
res <- tibble()
for (i in norm_types){
  for (j in diffanal_types){
    adj_val <- diffanal_results[[i]] %>% 
      filter(diffanal==j) %>% 
      unnest() %>% 
      select(gene_symbol, padj) %>%
      mutate(norm := !!i, diffanal := !!j, .after=gene_symbol)
    res <- bind_rows(res, adj_val)
  }
}

deseq2_results <- deseq2_anal(filt_data, metadata) %>% 
  select(gene_symbol, padj) %>%
  mutate(norm = 'deseq2', diffanal = 'deseq2', .after=gene_symbol)

res <- bind_rows(res, deseq2_results)