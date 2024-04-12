library(tidyverse)
library(qs)

# Performs the filtering, normalization and differential expression analysis. 
# It sources three different helper files containing the functions for each of these processes.
# The functions are called based on the pipeline argument.
# The output consists of a table consisting on log2 fold change, t-value and adjusted p-value for each gene.
###Main###
# Flop args
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

print(args[grep("--diffexp-limma-ndups", args) + 1])

# Diffexp args
additional_args <- list(
    filterbyexpr_libsize = args[grep("--filterbyexpr-libsize", args) + 1],
    filterbyexpr_mincount = args[grep("--filterbyexpr-mincount", args) + 1],
    filterbyexpr_mintotalcount = args[grep("--filterbyexpr-mintotalcount", args) + 1],
    filterbyexpr_largen = args[grep("--filterbyexpr-largen", args) + 1],
    filterbyexpr_minprop = args[grep("--filterbyexpr-minprop", args) + 1],
    diffexp_limma_ndups = args[grep("--diffexp-limma-ndups", args) + 1],
    diffexp_limma_spacing = args[grep("--diffexp-limma-spacing", args) + 1],
    diffexp_limma_block = args[grep("--diffexp-limma-block", args) + 1],
    diffexp_limma_weights = args[grep("--diffexp-limma-weights", args) + 1],
    diffexp_limma_method = args[grep("--diffexp-limma-method", args) + 1],
    diffexp_deseq2_test = args[grep("--diffexp-deseq2-test", args) + 1],
    diffexp_deseq2_fitType = args[grep("--diffexp-deseq2-fitType", args) + 1],
    diffexp_deseq2_quiet = args[grep("--diffexp-deseq2-quiet", args) + 1],
    diffexp_deseq2_minReplicatesForReplace = args[grep("--diffexp-deseq2-minReplicatesForReplace", args) + 1],
    diffexp_deseq2_parallel = args[grep("--diffexp-deseq2-parallel", args) + 1],
    diffexp_deseq2_betaprior = args[grep("--diffexp-deseq2-betaprior", args) + 1],
    diffexp_edger_calcnormfactors_method = args[grep("--diffexp-edger-calcnormfactors-method", args) + 1],
    diffexp_edger_calcnormfactors_refColumn = args[grep("--diffexp-edger-calcnormfactors-refColumn", args) + 1],
    diffexp_edger_calcnormfactors_logratiotrim = args[grep("--diffexp-edger-calcnormfactors-logratiotrim", args) + 1],
    diffexp_edger_calcnormfactors_sumtrim = args[grep("--diffexp-edger-calcnormfactors-sumtrim", args) + 1],
    diffexp_edger_calcnormfactors_doweighting = args[grep("--diffexp-edger-calcnormfactors-doweighting", args) + 1],
    diffexp_edger_calcnormfactors_acutoff = args[grep("--diffexp-edger-calcnormfactors-acutoff", args) + 1],
    diffexp_edger_calcnormfactors_p = args[grep("--diffexp-edger-calcnormfactors-p", args) + 1],
    diffexp_edger_estimatedisp_priordf = args[grep("--diffexp-edger-estimatedisp-priordf", args) + 1],
    diffexp_edger_estimatedisp_trendmethod = args[grep("--diffexp-edger-estimatedisp-trendmethod", args) + 1],
    diffexp_edger_estimatedisp_tagwise = args[grep("--diffexp-edger-estimatedisp-tagwise", args) + 1],
    diffexp_edger_estimatedisp_mixeddf = args[grep("--diffexp-edger-estimatedisp-mixeddf", args) + 1],
    diffexp_edger_estimatedisp_span = args[grep("--diffexp-edger-estimatedisp-span", args) + 1],
    diffexp_edger_estimatedisp_minrowsum = args[grep("--diffexp-edger-estimatedisp-minrowsum", args) + 1],
    diffexp_edger_estimatedisp_gridlength = args[grep("--diffexp-edger-estimatedisp-gridlength", args) + 1],
    diffexp_edger_estimatedisp_gridrange = args[grep("--diffexp-edger-estimatedisp-gridrange", args) + 1],
    diffexp_edger_estimatedisp_robust = args[grep("--diffexp-edger-estimatedisp-robust", args) + 1],
    diffexp_edger_estimatedisp_winsortailp = args[grep("--diffexp-edger-estimatedisp-winsortailp", args) + 1],
    diffexp_edger_estimatedisp_tol = args[grep("--diffexp-edger-estimatedisp-tol", args) + 1],
    diffexp_edger_glmqlfit_dispersion = args[grep("--diffexp-edger-glmqlfit-dispersion", args) + 1],
    diffexp_edger_glmqlfit_libsize = args[grep("--diffexp-edger-glmqlfit-libsize", args) + 1],
    diffexp_edger_glmqlfit_offset = args[grep("--diffexp-edger-glmqlfit-offset", args) + 1],
    diffexp_edger_glmqlfit_weights = args[grep("--diffexp-edger-glmqlfit-weights", args) + 1],
    diffexp_edger_glmqlfit_abundancetrend = args[grep("--diffexp-edger-glmqlfit-abundancetrend", args) + 1],
    diffexp_edger_glmqlfit_avelogcpm = args[grep("--diffexp-edger-glmqlfit-avelogcpm", args) + 1],
    diffexp_edger_glmqlfit_robust = args[grep("--diffexp-edger-glmqlfit-robust", args) + 1]
)

source(paste0(path_file, "normalization_helper.R"))
source(paste0(path_file, "diffexp_helper.R"))
source(paste0(path_file, "filtering_helper.R"))

counts <- qread(counts_file)
groups_sorting <- biocontext %>% strsplit("_v_") %>% unlist()
metadata <- qread(meta_file) %>% mutate(group = factor(group, levels = groups_sorting))

if (status == 'filtered') {
    filt_func <- get('filtering')
    newcounts <- filt_func(counts, metadata, additional_args)
} else {
    newcounts <- counts
}

print('filtering_ok')

if (length(pipeline) == 2) {
    norm_func <- get(pipeline[1])
    diffexp_func <- get(pipeline[2])
    results <- norm_func(newcounts, metadata) %>% diffexp_func(., metadata, additional_args)
} else {
    diffexp_func <- get(pipeline[1])
    results <- diffexp_func(newcounts, metadata, additional_args)
}

pipeline <- sub(pattern = "_.*", replacement =  "", pipeline)

if (length(pipeline) == 1) {
    pipeline <- c('NA', pipeline)
    print(pipeline)
}

output_filename <- paste(subset_id, biocontext, status, pipeline[1], pipeline[2], "de", sep = "__") %>%
  paste(., ".qs", sep = "")

qsave(results, output_filename)
