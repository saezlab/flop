library(tidyverse)
library(qs)
#library(rstudioapi)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# counts_file <- "./data/CCLE_ABCDEF/BONE_v_BREAST__countdata.tsv"
# meta_file <- "./data/CCLE_ABCDEF/BONE_v_BREAST__metadata.tsv"
# pipeline <- c("voom_norm", "limma_analysis")
# status <- "filtered"
# path_file <- ""

###Main###
args <- commandArgs(trailingOnly = FALSE)
path_file <- args[grep("--file=", args)] %>%
    sub("diffexp_analysis.R", "", .) %>%
    sub("--file=", "", .)
counts_file <- args[grep("--counts", args) + 1]
meta_file <- args[grep("--meta", args) + 1]
status <- args[grep("--status", args) + 1]
subset_id <- args[grep("--dataset", args) + 1]
biocontext <- args[grep("--bio", args) + 1]
pipeline <- args[grep("--pipeline", args) + 1] %>%
    strsplit(., split = " ") %>%
    .[[1]]

source(paste0(path_file, "normalization_helper.R"))
source(paste0(path_file, "diffexp_helper.R"))
source(paste0(path_file, "filtering_helper.R"))

counts <- qread(counts_file)
metadata <- qread(meta_file) %>% mutate(group = factor(group))

if (status == 'filtered') {
    filt_func <- get('filtering')
    newcounts <- filt_func(counts, metadata)
} else {
    newcounts <- counts
}

if (length(pipeline) == 2) {
    norm_func <- get(pipeline[1])
    diffexp_func <- get(pipeline[2])
    results <- norm_func(newcounts, metadata) %>% diffexp_func(., metadata)
} else {
    diffexp_func <- get(pipeline[1])
    results <- diffexp_func(newcounts, metadata)
}

pipeline <- sub(pattern = "_.*", replacement =  "", pipeline)

if (length(pipeline) == 1) {
    pipeline <- c('NA', pipeline)
    print(pipeline)
}

output_filename <- paste(subset_id, biocontext, status, pipeline[1], pipeline[2], "de", sep = "__") %>%
  paste(., ".qs", sep = "")

qsave(results, output_filename)
