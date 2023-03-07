library(tidyverse)
library(testthat)
library(usethis)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../diffexp_helper.R", chdir = TRUE)
source("../normalization_helper.R", chdir = TRUE)

test_vsn <- function

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
test_that("deseq2_analysis returns a non-empty tibble with the correct columns", {
    results <- test_deseq2()
    expect_true(nrow(results) > 0)
    expect_true(ncol(results) == 4)
    expect_true(is.numeric(results$logFC))
    expect_true(is.numeric(results$stat))
    expect_true(is.numeric(results$padj))
})

test_that("edger_analysis returns a non-empty tibble with the correct columns", {
    results <- test_edger()
    expect_true(nrow(results) > 0)
    expect_true(ncol(results) == 4)
    expect_true(is.numeric(results$logFC))
    expect_true(is.numeric(results$stat))
    expect_true(is.numeric(results$padj))
})
