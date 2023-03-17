# Add instructions to readme

Created time: March 16, 2023 5:52 PM
Status: In progress

# Description

FLOP benchmarking is a unified workflow that analyses bulk RNA-seq counts using multiple combinations of filtering, normalisation and differential expression methods. It then evaluates the differences in the functional space between the different combinations of methods. 

## Installation

To install FLOP, you first will need to download the files from our GitHub repository:

```bash
git clone https://github.com/saezlab/flop_benchmark
```

Once downloaded, you can just execute the launcher command and it will install the necessary dependencies:

```bash
bash flop_launcher.sh
```

You are set up! The launcher will guide you through the rest of the analysis. You can reuse this command every time you launch a new analysis. 

# Input

FLOP works with two or three different files. Each different dataset folder should contain at least two of the three files specified below. FLOP support subseting of large datasets (in this case, the results are averaged for the evaluation modules), you can specify a subset using this naming format: {Dataset ID}_{Subset identifier} in the folder name and in the files names. One example data directory could be:

```bash
/data/
	./GSE186341/
		./GSE186341__countdata.tsv
		./GSE186341__metadata.tsv
		./GSE186341__contrast.tsv
	./GTEx_subset1/
		./GTEx_subset1__countdata.tsv
		./GTEx_subset1__metadata.tsv
		./GTEx_subset1__contrast.tsv
	./GTEX_subset2/
		./GTEx_subset2__countdata.tsv
		./GTEx_subset2__metadata.tsv
		./GTEx_subset2__contrast.tsv
	./CCLE/
		./CCLE__countdata.tsv
		./CCLE__metadata.tsv
```

## Count table

This tab-separated table should contain the gene counts of the dataset. It should have the following name pattern: {Dataset ID}__countdata.tsv 

The data should contain a column named gene_symbol, which contains the gene symbols, and the following columns should be the samples that will be included in the study. Be aware that the sample names should respect R’s naming regulations (see the command make.names). 

## Metadata table

This tab-separated table should at least contain the the group to which the samples belong. It should have the following name pattern: {Dataset ID}__metadata.tsv 

This table contains a column named sample_ID, with the names of the samples (again, they should respect R’s naming guidelines) and a column named group, which contains the group every sample belongs to. Other columns can be present, but they will not be considered in the analysis.

## Contrast table

This tab-separated table should contain the contrasts that should be calculated in the analysis. It should have the following name pattern: {Dataset ID}__contrast.tsv 

The table should contain two columns, group1 and group2, with the desired contrasts.

This table is optional. If included, FLOP will only calculate the specified contrast. If this table is not included, it will calculate every pairwise comparison between the different groups (be careful with choosing this option when analysing very large datasets, since it will greatly increase the execution time!).

# Output

The output consists in three or different files:

## General results file:

A long format table that contains the functional scores of all different prior knowledge sources and per pipeline, biological contrast, filtering status, subset (if applicable) and t-value or log fold change.

## Rank correlation results:

A long format table that contains the spearman rank correlation scores per comparison, parameter, filtering status, biological context and prior knowledge source. The comparisons are done at the pipeline level, normalisation level and differential expression method level.

## Jaccard index results

A long format table that contains the Jaccard index scores per pipeline comparison, parameter, filtering status, biological context and prior knowledge source. The number of top and bottom functional categories included vary between the different prior knowledge source:

- Dorothea: 15 top and 15 bottom
- MSigDB hallmarks: top 5 and bottom 5
- PROGENy: top 3 and bottom 3

## Rand index analysis

A long format table that contains the Rand index analysis for 11 and 32 K’s, per pipeline comparison, filtering status, biological context and prior knowledge source.

Since this analysis is only informative when there are biological factors (such a specific number of cell lines, treatments, etc.) that makes the clustering of the samples important for a study, we implemented this analysis as optimal inside the FLOP architecture. You can select the datasets that will be included in this analysis during the initial configuration.