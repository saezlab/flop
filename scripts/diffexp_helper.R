library(tidyverse)
library(limma)
library(edgeR)
library(DESeq2)

#' Limma differential expression analysis
#' This function performs a limma differnetial expression analysis with 
#' intersection and contrast. 
#' @param data_norm a tibble containing normalised values of gene counts; genes 
#' as rows, samples as columns.
#' @param metadata a tibble containing the metadata related to the normalised
#' data containing the following columns: sample_ID, group (2 values only) 
#'
#' @return a tibble containing statistical parameters for each gene.
#' @importFrom limma lmFit makeContrasts contrasts.fit topTable
#' @importFrom dplyr %>% select add_column as_tibble 
#' @export
#'
#' @examples
#' limma_analysis(vsn_counts_exp1, metadata_exp1)
limma_analysis <- function(data_norm, metadata, additional_args){
  norm_counts <- data_norm %>% dplyr::select_if(is.numeric)
  
  # This removes covariates that have only one level or NA values
  covariates <- metadata %>%
    dplyr::select(-sample_ID) %>% 
    .[, sapply(., Negate(anyNA)), drop = FALSE] %>%
    sapply(., function(x) n_distinct(x)) %>% as.data.frame() %>% filter(. > 1) %>% rownames()

  designmat <- metadata %>%
    model.matrix(as.formula(rlang::parse_expr(paste(c("~ 0", covariates), collapse = " + "))), data=.)
  
  studygroups <- levels(metadata$group)
  
  coefs <- paste("group", studygroups[1], "-", "group", studygroups[2], sep='')
  
  ndups <- as.integer(additional_args$diffexp_limma_ndups)
  print('nduspok')
  spacing <- as.integer(additional_args$diffexp_limma_spacing)
  print('spacingok')
  block <- additional_args$diffexp_limma_block
  if (block == "NULL") {
      block <- NULL
  }
  print("blockok")
  weights <- additional_args$diffexp_limma_weights
  if (weights == "NULL") {
      weights <- NULL
  }
  method <- additional_args$diffexp_limma_method
  print('methodok')
  
  fit <- lmFit(norm_counts, 
               design = designmat, 
               ndups = ndups, 
               spacing = spacing, 
               block = block, 
               weights = weights, 
               method = method)

  fit <- coefs %>%
    makeContrasts(levels = designmat, contrasts = .) %>%
    contrasts.fit(fit = fit, contrasts = .) %>%
    eBayes(.)
  
  results <- topTable(fit, number = Inf, coef = coefs, sort.by='none') %>%
    as_tibble() %>%
    mutate(ID = data_norm$ID) %>%
    dplyr::rename(padj = adj.P.Val, stat = t) %>%
    dplyr::select(ID, logFC, stat, padj)
  
  return(results)
}

#' DESeq2 differential expression analysis
#'
#' @param counts a tibble containing non-normalised gene counts; genes as rows, samples
#' as columns
#' @param metadata a tibble containing the metadata related to the normalised
#' data containing the following columns: sample_ID, group (2 values only) 
#'
#' @return a tibble containing statistical parameters for each gene.
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom dplyr %>% select as_tibble add_column
#' @export
#'
#' @examples
#' deseq2_analysis(raw_counts, metadata_cancer)
deseq2_analysis <- function(counts, metadata, additional_args) {
  gene_counts <- counts %>% 
    as.data.frame(.)
  rownames(gene_counts) <- counts$gene_symbol
  gene_counts <- gene_counts %>% 
    dplyr::select(-gene_symbol)

  # This removes covariates that have only one level or NA values
  covariates <- metadata %>%
    dplyr::select(-sample_ID) %>% 
    .[, sapply(., Negate(anyNA)), drop = FALSE] %>%
    sapply(., function(x) n_distinct(x)) %>% as.data.frame() %>% filter(. > 1) %>% rownames()
  
  formatted_data <- DESeqDataSetFromMatrix(countData = gene_counts, 
                                           colData = metadata, 
                                           design = as.formula(rlang::parse_expr(paste(c("~ 0", covariates), collapse = " + "))))
  studygroups <- levels(metadata$group)
  coefs <- paste("group", studygroups[1], "-", "group", studygroups[2], sep='')

  deseq2_test <- additional_args$diffexp_deseq2_test
  deseq2_fitType <- additional_args$diffexp_deseq2_fitType
  deseq2_betaprior <- eval(parse(text=additional_args$diffexp_deseq2_betaprior))
  deseq2_quiet <- eval(parse(text=additional_args$diffexp_deseq2_quiet))
  deseq2_minReplicatesForReplace <- as.numeric(additional_args$diffexp_deseq2_minReplicatesForReplace)
  deseq2_parallel <- eval(parse(text=additional_args$diffexp_deseq2_parallel))
  
  results <- DESeq(formatted_data,
                    test = deseq2_test,
                    fitType = deseq2_fitType,
                    betaPrior = deseq2_betaprior,
                    quiet = deseq2_quiet,
                    minReplicatesForReplace = deseq2_minReplicatesForReplace,
                    parallel = deseq2_parallel) %>% 
    results(., contrast=c('group', studygroups[1], studygroups[2])) %>%
    as_tibble() %>%
    mutate(ID = counts$gene_symbol)  %>%
    dplyr::rename(logFC = log2FoldChange) %>%
    dplyr::select(ID, logFC, stat, padj)
  return(results)
}

#' EdgeR differential expression analysis
#'
#' @param counts a tibble containing non-normalised values of gene counts; genes 
#' as rows, samples as columns.
#' @param metadata a tibble containing the metadata related to the normalised
#' data containing the following columns: sample_ID, group (2 values only) 
#'
#' @return a tibble containing statistical parameters for each gene.
#' @importFrom edgeR estimateDisp glmFit glmLRT topTags DGElist
#' @importFrom dplyr %>% select as_tibble add_column
#' @export
#'
#' @examples
#' edger_analysis(vsn_results, metadata_cancer)

edger_analysis <- function(counts, metadata, additional_args) {
  dge_obj <- counts %>% dplyr::select_if(is.numeric) %>% DGEList(count=., samples = metadata)
  
  # This removes covariates that have only one level or NA values
  covariates <- metadata %>%
    dplyr::select(-sample_ID) %>% 
    .[, sapply(., Negate(anyNA)), drop = FALSE] %>%
    sapply(., function(x) n_distinct(x)) %>% as.data.frame() %>% filter(. > 1) %>% rownames()

  designmat <- metadata %>%
    model.matrix(as.formula(rlang::parse_expr(paste(c("~ 0", covariates), collapse = " + "))), data=.)
  
  studygroups <- levels(metadata$group)
  coefs <- paste("group", studygroups[1], "-", "group", studygroups[2], sep = '')
  
  # Extract additional arguments for calcNormFactors
  calcnormfactors_method <- additional_args$diffexp_edger_calcnormfactors_method
  print(calcnormfactors_method)
  calcnormfactors_refColumn <- additional_args$diffexp_edger_calcnormfactors_refColumn
  if(calcnormfactors_refColumn == "NULL") {
    calcnormfactors_refColumn <- NULL
  }
  print(calcnormfactors_refColumn)
  calcnormfactors_logratiotrim <- as.numeric(additional_args$diffexp_edger_calcnormfactors_logratiotrim)
  print(calcnormfactors_logratiotrim)
  calcnormfactors_sumtrim <- as.numeric(additional_args$diffexp_edger_calcnormfactors_sumtrim)
  print(calcnormfactors_sumtrim)
  calcnormfactors_doweighting <- eval(parse(text=additional_args$diffexp_edger_calcnormfactors_doweighting))
  print(calcnormfactors_doweighting)
  calcnormfactors_acutoff <- as.numeric(additional_args$diffexp_edger_calcnormfactors_acutoff)
  print(calcnormfactors_acutoff)
  calcnormfactors_p <- as.numeric(additional_args$diffexp_edger_calcnormfactors_p)
  print(calcnormfactors_p)
  
  # Normalize counts
  dge_obj <- calcNormFactors(dge_obj,
                              method = calcnormfactors_method,
                              refColumn = calcnormfactors_refColumn,
                              logratioTrim = calcnormfactors_logratiotrim,
                              sumTrim = calcnormfactors_sumtrim,
                              doWeighting = calcnormfactors_doweighting,
                              Acutoff = calcnormfactors_acutoff,
                              p = calcnormfactors_p)
  
  # Extract additional arguments for estimateDisp
  estimatedisp_priordf <- additional_args$diffexp_edger_estimatedisp_priordf
  if(estimatedisp_priordf == "NULL") {
    estimatedisp_priordf <- NULL
  }
  estimatedisp_trendmethod <- additional_args$diffexp_edger_estimatedisp_trendmethod
  estimatedisp_tagwise <- eval(parse(text=additional_args$diffexp_edger_estimatedisp_tagwise))
  estimatedisp_mixeddf <- eval(parse(text=additional_args$diffexp_edger_estimatedisp_mixeddf))
  estimatedisp_span <- additional_args$diffexp_edger_estimatedisp_span
  if(estimatedisp_span == "NULL") {
    estimatedisp_span <- NULL
  }
  estimatedisp_minrowsum <- as.numeric(additional_args$diffexp_edger_estimatedisp_minrowsum)
  estimatedisp_gridlength <- as.numeric(additional_args$diffexp_edger_estimatedisp_gridlength)
  estimatedisp_gridrange <- eval(parse(text = additional_args$diffexp_edger_estimatedisp_gridrange))
  estimatedisp_robust <- eval(parse(text=additional_args$diffexp_edger_estimatedisp_robust))
  estimatedisp_winsortailp <- eval(parse(text=additional_args$diffexp_edger_estimatedisp_winsortailp))
  estimatedisp_tol <- as.numeric(additional_args$diffexp_edger_estimatedisp_tol)
  
  # Estimate dispersion
  dge_obj <- estimateDisp(dge_obj,
                          design = designmat,
                          prior.df = estimatedisp_priordf,
                          trend.method = estimatedisp_trendmethod,
                          tagwise = estimatedisp_tagwise,
                          mixed.df = estimatedisp_mixeddf,
                          span = estimatedisp_span,
                          min.row.sum = estimatedisp_minrowsum,
                          grid.length = estimatedisp_gridlength,
                          grid.range = estimatedisp_gridrange,
                          robust = estimatedisp_robust,
                          winsor.tail.p = estimatedisp_winsortailp,
                          tol = estimatedisp_tol)
  
  # Extract additional arguments for glmQLFit
  glmqlfit_dispersion <- additional_args$diffexp_edger_glmqlfit_dispersion
  if(glmqlfit_dispersion == "NULL") {
    glmqlfit_dispersion <- NULL
  }
  
  glmqlfit_libsize <- additional_args$diffexp_edger_glmqlfit_libsize
  if(glmqlfit_libsize == "NULL") {
    glmqlfit_libsize <- NULL
  }
  glmqlfit_offset <- additional_args$diffexp_edger_glmqlfit_offset
  if(glmqlfit_offset == "NULL") {
    glmqlfit_offset <- NULL
  }
  glmqlfit_weights <- additional_args$diffexp_edger_glmqlfit_weights
  if(glmqlfit_weights == "NULL") {
    glmqlfit_weights <- NULL
  }
  glmqlfit_abundancetrend <- eval(parse(text=additional_args$diffexp_edger_glmqlfit_abundancetrend))
  glmqlfit_avelogcpm <- additional_args$diffexp_edger_glmqlfit_avelogcpm
  if(glmqlfit_avelogcpm == "NULL") {
    glmqlfit_avelogcpm <- NULL
  }
  glmqlfit_robust <- eval(parse(text=additional_args$diffexp_edger_glmqlfit_robust))
  
  # Perform glmQLFit
  fit <- glmQLFit(dge_obj,
                  design = designmat,
                  dispersion = glmqlfit_dispersion,
                  abundance.trend = glmqlfit_abundancetrend,
                  robust = glmqlfit_robust,
                  winsor.tail.p = estimatedisp_winsortailp,)
  
  # Perform glmQLFTest
  contrast_mat <- makeContrasts(levels = designmat, 
                                 contrasts = coefs) 
  
  results <- glmQLFTest(fit, contrast = contrast_mat) %>% 
    topTags(n = Inf, sort.by = "none", adjust.method = "BH") %>% 
    as.data.frame() %>% 
    as_tibble() %>%
    mutate(ID = counts$gene_symbol) %>%
    dplyr::rename(padj = FDR) %>%
    mutate(stat = sign(logFC) * sqrt(F)) %>%
    dplyr::select(ID, logFC, stat, padj)
  
  return(results)
}

