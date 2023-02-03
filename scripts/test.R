library(tidyverse)
library(cowplot)
library(ggplot2)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

filtered_files <- list.files("./results/GSE186341/filtered/", pattern = "*de.rds", full.names = T)
unfiltered_files <- list.files("./results/GSE186341/unfiltered/", pattern = "*de.rds", full.names = T)

bio_contexts <- filtered_files %>% 
    as_tibble() %>% 
    separate(., col=value, into = c("path", "bio_context", "norm", "diffexpr"), sep = "__") %>% 
    select(bio_context) %>% 
    distinct() %>% 
    pull()

pipelines <- filtered_files %>% 
    as_tibble() %>% 
    separate(., col=value, into = c("path", "bio_context", "norm", "diffexpr"), sep = "__") %>%
    select(norm, diffexpr) %>%
    mutate(pipeline = paste(norm, "__", diffexpr, sep="")) %>%
    select(pipeline) %>%
    distinct() %>% 
    pull()

for (pipeline in pipelines){
    pipeline_files <- grep(pipeline, filtered_files, value = TRUE)
    tables <- lapply(pipeline_files, readRDS)
    names(tables) <- bio_contexts
    imap(tables[1:5], function(x, tit){
    x %>%
        ggplot(aes(x = logFC, y = -log10(padj))) +
        geom_point() +
        ggtitle(tit)
    }) %>%
    plot_grid(plotlist = .) %>%
    save_plot(., device = "png", filename = paste0("./results/GSE186341/filtered/", pipeline, ".png"), base_width=10, base_height=7)

}

for (pipeline in pipelines){
    pipeline_files <- grep(pipeline, unfiltered_files, value = TRUE)
    tables <- lapply(pipeline_files, readRDS)
    names(tables) <- bio_contexts
    imap(tables[1:5], function(x, tit){
    x %>%
        ggplot(aes(x = logFC, y = -log10(padj))) +
        geom_point() +
        ggtitle(tit)
    }) %>%
    plot_grid(plotlist = .) %>%
    save_plot(., device = "png", filename = paste0("./results/GSE186341/unfiltered/", pipeline, ".png"), base_width=10, base_height=7)

}
