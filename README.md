# README

# Description

FLOP is a unified workflow that analyses bulk RNA-seq counts using multiple combinations of filtering, normalisation and differential expression methods. It then evaluates the differences in the functional space between the different combinations of methods. 

## Installation

To install FLOP, you first will need to download the files from our GitHub repository:

```bash
git clone https://github.com/saezlab/flop
```

To run FLOP, you need to have conda installed in your computer. Please check [this link](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to learn how to install conda.

This command will install the necessary dependencies inside a dedicated environment.

```bash
cd flop/
conda env create -f ./scripts/flop_env.yaml
conda activate flop
```

Once installed, you are ready to run FLOP!

## Quick Start

It is possible to run an example version of FLOP with a test dataset containing three contrasts from the PANACEA study via:

```bash
bash flop_launcher.sh -t
```

## Run

Mode of usage:

```bash
bash flop_launcher.sh [-d data folder] [-e config set][-p pvalue threshold] [-f n DE genes threshold] [-t] [-s] [-h]
```

FLOP has several ways of personalization. These are all possible input parameters:

- -d: data folder, containing the subfolders with the datasets to be analyzed
- -e: config set, either 'desktop' or 'cluster'. If none specified, it defaults to desktop
- -t: test mode, runs the pipeline with the test dataset and default parameters. Bear in mind that you still need to specify a config set with -e
- -p: pvalue threshold that the genes or functional terms need to pass in order to be considered significant for the Jaccard Index module. Default is 1 (no filtering).
- -f: Minimum number of significant genes per contrast. Only contrasts that have a minimum of n genes with a pvalue below 0.05 will be considered for enrichment analysis. Default is 0 (no filtering).
- -s: Launches FLOP to run the analysis detailed in the accompanying study. It runs FLOP with CCLE, PANACEA and Reheat datasets, and then outputs the figures also shown in the study. It sets the pvalue threshold to 1 and the number of significant genes threshold to 30. For more information, please check the study. It is strongly advised to run this setup within a HPC environment (config set = cluster).
- -h: shows this help message

You can run FLOP with the minimal settings by using this command:

```bash
bash flop_launcher.sh [-d data_folder]
```

To run the analysis showed in the study using PANACEA, CCLE and ReHeat (it is recommended to run it in a HPC environment):

```bash
bash flop_launcher.sh -s -e cluster
```

# Input

FLOP works with two or three different files. Each different dataset folder should contain at least two of the three files specified below. FLOP support subseting of large datasets (in this case, the results are averaged for the evaluation modules), you can specify a subset using this naming format: {Dataset ID}_{Subset identifier} in the folder name and in the files names. One example data directory could be:

```bash
/data/
	./GSE186341/
		./GSE186341__countdata.tsv
		./GSE186341__metadata.tsv
		./GSE186341__contrast.tsv
	./Reheat_subset1/
		./Reheat_study1__countdata.tsv
		./Reheat_study1__metadata.tsv
		./Reheat_study1__contrast.tsv
	./Reheat_study2/
		./Reheat_study2__countdata.tsv
		./Reheat_study2__metadata.tsv
		./Reheat_study2__contrast.tsv
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

The output consists in three or four different files:

## General results file:

A long format table that contains the functional scores of all different prior knowledge sources and per pipeline, biological contrast, filtering status, subset (if applicable) and t-value or log fold change.

## Rank correlation results:

A long format table that contains the spearman rank correlation and coocurrence scores per comparison, parameter, filtering status, biological context and prior knowledge source.

## Top/bottom features overlap results

A long format table that contains the Top and Bottom features overlap index and similarity scores per pipeline comparison, parameter, filtering status, biological context and prior knowledge source. The number of top and bottom functional categories included vary between the different prior knowledge source:

- Dorothea: 15 top and 15 bottom
- MSigDB hallmarks: top 5 and bottom 5
- PROGENy: top 3 and bottom 3
