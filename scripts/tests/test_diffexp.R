library(tidyverse)
library(testthat)
library(usethis)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../diffexp_helper.R", chdir = TRUE)
source("../normalization_helper.R", chdir = TRUE)


test_vsn <- function() {
    gene_symbol <- c("gene1", "gene2", "gene3", "gene4", "gene5")
    s1 <- c(1, 4, 3, 6, 2)
    s2 <- c(2, 5, 4, 70, 3)
    s3 <- c(2000, 6, 50, 800, 4000)
    s4 <- c(4000, 7, 60, 900, 5000)
    gene_counts <- tibble(gene_symbol, s1, s2, s3, s4)
    results <- vsn_norm(gene_counts, minDataPointsPerStratum = 4)
    return(results)
}

# test that vsn_norm returns a non-empty tibble with the correct columns. Check also if the columns except for gene_symbol are numeric
test_that("vsn_norm returns a non-empty tibble with the correct columns", {
    results <- test_vsn()
    expect_true(nrow(results) > 0)
    expect_true(is.numeric(results$s1))
    expect_true(is.numeric(results$s2))
    expect_true(is.numeric(results$s3))
    expect_true(is.numeric(results$s4))
})

test_log2quant <- function() {
    gene_symbol <- c("gene1", "gene2", "gene3", "gene4", "gene5")
    s1 <- c(1, 4, 3, 6, 2)
    s2 <- c(2, 5, 4, 70, 3)
    s3 <- c(2000, 6, 50, 800, 4000)
    s4 <- c(4000, 7, 60, 900, 5000)
    gene_counts <- tibble(gene_symbol, s1, s2, s3, s4)
    results <- log2quant_norm(gene_counts)
    return(results)
}

# test that log2quant_norm returns a non-empty tibble with the correct columns. Check also if the columns except for gene_symbol are numeric
test_that("log2quant_norm returns a non-empty tibble with the correct columns", {
    results <- test_log2quant()
    expect_true(nrow(results) > 0)
    expect_true(is.numeric(results$s1))
    expect_true(is.numeric(results$s2))
    expect_true(is.numeric(results$s3))
    expect_true(is.numeric(results$s4))
})

test_tmm <- function() {
    gene_symbol <- c("gene1", "gene2", "gene3", "gene4", "gene5")
    s1 <- c(1, 4, 3, 6, 2)
    s2 <- c(2, 5, 4, 70, 3)
    s3 <- c(2000, 6, 50, 800, 4000)
    s4 <- c(4000, 7, 60, 900, 5000)
    gene_counts <- tibble(gene_symbol, s1, s2, s3, s4)
    results <- tmm_norm(gene_counts)
    return(results)
}

# test that tmm_norm returns a non-empty tibble with the correct columns. Check also if the columns except for gene_symbol are numeric
test_that("tmm_norm returns a non-empty tibble with the correct columns", {
    results <- test_tmm()
    expect_true(nrow(results) > 0)
    expect_true(is.numeric(results$s1))
    expect_true(is.numeric(results$s2))
    expect_true(is.numeric(results$s3))
    expect_true(is.numeric(results$s4))
})

test_voom <- function() {
    gene_symbol <- c("gene1", "gene2", "gene3", "gene4", "gene5")
    s1 <- c(1, 4, 3, 6, 2)
    s2 <- c(2, 5, 4, 70, 3)
    s3 <- c(2000, 6, 50, 800, 4000)
    s4 <- c(4000, 7, 60, 900, 5000)
    gene_counts <- tibble(gene_symbol, s1, s2, s3, s4)
    metadata <- tibble(sample_ID = c("s1", "s2", "s3", "s4"), group = as.factor(c("A", "A", "B", "B")))
    results <- voom_norm(gene_counts, metadata)
    return(results)
}

# test that voom_norm returns a non-empty tibble with the correct columns. Check also if the columns except for gene_symbol are numeric
test_that("voom_norm returns a non-empty tibble with the correct columns", {
    results <- test_voom()
    expect_true(nrow(results) > 0)
    expect_true(is.numeric(results$s1))
    expect_true(is.numeric(results$s2))
    expect_true(is.numeric(results$s3))
    expect_true(is.numeric(results$s4))
})

test_deseq2 <- function() {
    gene_symbol <- c("gene1", "gene2", "gene3", "gene4", "gene5")
    s1 <- c(1, 4, 3, 6, 2)
    s2 <- c(2, 5, 4, 70, 3)
    s3 <- c(2000, 6, 50, 800, 4000)
    s4 <- c(4000, 7, 60, 900, 5000)
    gene_counts <- tibble(gene_symbol, s1, s2, s3, s4)
    metadata <- tibble(sample_ID = c("s1", "s2", "s3", "s4"), group = as.factor(c("A", "A", "B", "B")))
    results <- deseq2_analysis(gene_counts, metadata)
    return(results)
}

# test if the function returns a non-empty tibble with the correct columns. Check also if the columns logFC, stat and padj are numeric
test_that("deseq2_analysis returns a non-empty tibble with the correct columns", {
    results <- test_deseq2()
    expect_true(nrow(results) > 0)
    expect_true(ncol(results) == 4)
    expect_true(is.numeric(results$logFC))
    expect_true(is.numeric(results$stat))
    expect_true(is.numeric(results$padj))
})

test_edger <- function() {
    gene_symbol <- c("gene1", "gene2", "gene3", "gene4", "gene5")
    s1 <- c(1, 4, 3, 6, 2)
    s2 <- c(2, 5, 4, 70, 3)
    s3 <- c(2000, 6, 50, 800, 4000)
    s4 <- c(4000, 7, 60, 900, 5000)
    gene_counts <- tibble(gene_symbol, s1, s2, s3, s4)
    metadata <- tibble(sample_ID = c("s1", "s2", "s3", "s4"), group = as.factor(c("A", "A", "B", "B")))
    results <- edger_analysis(gene_counts, metadata)
    return(results)
}

# test if the function returns a non-empty tibble with the correct columns. Check also if the columns logFC, stat and padj are numeric
test_that("edger_analysis returns a non-empty tibble with the correct columns", {
    results <- test_edger()
    expect_true(nrow(results) > 0)
    expect_true(ncol(results) == 4)
    expect_true(is.numeric(results$logFC))
    expect_true(is.numeric(results$stat))
    expect_true(is.numeric(results$padj))
})

devtools::test()
install.packages('htmltools')

unlink("C:/Users/victo/AppData/Local/R/win-library/4.2/00LOCK", recursive = TRUE)
