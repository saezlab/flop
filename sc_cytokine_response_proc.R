library(Seurat)
library(tidyverse)

files <- list.files('./unproc_data/hallmarks', pattern = 'RDS$', full.names = TRUE)
dir.create('./hallmarks_benchmark/')

full_counts <- tibble(gene_symbol= character())
full_meta <- tibble()
full_contrast <- tibble()

# MSIGDB hallmarks mouse gene sets:


for(file in files){
    seu.s <- readRDS(file)
    seu.s@meta.data <- seu.s@meta.data %>% mutate(celltype = gsub('_', '', celltype),
    sample = gsub('-', '', sample),
    sample_ID = paste0(celltype, '_', sample, '_', rep))

    
    metadata <- seu.s@meta.data %>% rownames_to_column('barcodes') %>% as_tibble()
    celltype <- metadata %>% pull(celltype) %>% unique()

    initial_cells <- metadata %>% pull(barcodes) %>% unique()

    qc_data <- metadata %>%
        group_by(sample_ID, sample, celltype, rep) %>%
        summarise(ncells = n(), counts = sum(nCount_RNA)) %>%
        filter(ncells > 10, counts > 1000)

    seu.s_filt <- subset(x = seu.s, subset =  sample_ID %in% qc_data$sample_ID)

    pseudobulk_counts <- AggregateExpression(seu.s_filt, group.by = c("sample_ID"), return.seurat = FALSE) %>% .$RNA %>%
        as.data.frame(.)

    colnames(pseudobulk_counts) <- colnames(pseudobulk_counts) %>% gsub('-', '_', .)

    setdiff(colnames(pseudobulk_counts), qc_data$sample_ID)

    pseudobulk_counts <- pseudobulk_counts  %>%
        rownames_to_column('gene_symbol') %>% as_tibble()

    pseudobulk_meta <- qc_data %>% 
        ungroup() %>%
        select(sample_ID, celltype, sample, rep) %>% 
        mutate(
            'group' = paste0(celltype, '_', sample)) %>% 
        distinct() %>%
        select(sample_ID, group)

        
    # get names that are in the pseudobulk colnames but not in the sample_ID column in pseudobulk_meta
    control_group <- pseudobulk_meta %>% filter(str_detect(group, 'PBS')) %>% pull(group) %>% unique()
    treatment_groups <- pseudobulk_meta %>% filter(!str_detect(group, 'PBS')) %>% pull(group) %>% unique()

    if(length(treatment_groups) == 0 | length(control_group) == 0){
        next
    }

    pseudobulk_contrast <- tibble(group1 = treatment_groups, group2 = control_group)

    full_counts <- full_join(full_counts, pseudobulk_counts, by = 'gene_symbol')
    full_meta <- bind_rows(full_meta, pseudobulk_meta)
    full_contrast <- bind_rows(full_contrast, pseudobulk_contrast)

}

dir.create('./hallmarks_benchmark/Atlascyt')

write_tsv(full_counts, './hallmarks_benchmark/Atlascyt/Atlascyt__countdata.tsv')
write_tsv(full_meta, './hallmarks_benchmark/Atlascyt/Atlascyt__metadata.tsv')
write_tsv(full_contrast, './hallmarks_benchmark/Atlascyt/Atlascyt__contrast.tsv')


# we lost 4 celltypes due to pseudobulk filtering
full_meta %>% separate(group, c('celltype', 'cyt'), sep = '_') %>% 
    distinct(celltype)

metadata %>% 
  	ggplot(aes(color=sample, x=nCount_RNA, fill=sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	cowplot::theme_cowplot() +
  	ylab("Cell density") +
  	geom_vline(xintercept = 500, linetype = 'dashed') +
    theme(text = element_text(family='Calibri'), legend.position = 'none')

# 87 treatments
# metadata %>%
#   	ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
#   	geom_density(alpha = 0.2) +
#   	theme_classic() +
#   	geom_vline(xintercept = 0.8)

# exploratory analysis

full_qc_data <- tibble()

for(file in files){
    seu.s <- readRDS(file)
    metadata <- seu.s@meta.data %>% rownames_to_column('barcodes') %>% as_tibble()
    celltype <- metadata %>% pull(celltype) %>% unique() %>% gsub('_', '', .)

    initial_cells <- metadata %>% pull(barcodes) %>% unique()

    qc_data <- metadata %>%
        group_by(sample, celltype, rep) %>%
        summarise(ncells = n(), counts = sum(nCount_RNA))
    
    full_qc_data <- bind_rows(full_qc_data, qc_data)
}


full_qc_data %>%
    mutate(biocontext = paste0(celltype, '_', sample)) %>%
    ggplot(aes(x = ncells, y = counts, color = celltype)) +
    geom_point() +
    cowplot::theme_cowplot() +
    scale_y_log10() +
    theme(text = element_text(family='Calibri')) +
    geom_hline(yintercept = 1000) +
    geom_vline(xintercept = 10)

full_qc_data %>%
    mutate(biocontext = paste0(celltype, '_', sample)) %>%
    ggplot(aes(x = ncells, y = counts, color = celltype)) +
    geom_point() +
    cowplot::theme_cowplot() +
    theme(text = element_text(family='Calibri')) +
    geom_hline(yintercept = 1000) +
    geom_vline(xintercept = 10)

full_qc_data %>% group_by(celltype) %>% filter(sample=='PBS') %>% summarise(n=n())


# per celltype, if there are no pbs samples, remove all rows containing cell type
full_qc_data %>% filter(ncells > 10, counts > 1000) %>%
    group_by(celltype) %>%
    mutate(pbs = ifelse(sample == 'PBS', 1, 0)) %>%
    summarise(pbs = sum(pbs)) %>%
    filter(pbs != 0) %>%
    pull(celltype) %>%
    unique()


