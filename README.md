# README

# Description

FLOP is a unified workflow that analyses bulk RNA-seq counts using multiple combinations of filtering, normalisation and differential expression methods. It then evaluates the differences in the functional space between the different combinations of methods. 

## Installation

To install FLOP, you first will need to download the files from our GitHub repository:

```bash
git clone https://github.com/saezlab/flop
```

To run FLOP, you need to have conda installed in your computer. Please check [this link](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to learn how to install conda.

This command will install the necessary dependencies inside an environment.

```bash
conda env create -f ./scripts/config_env.yaml
conda activate flop
```

Once installed, you are ready to run FLOP!

## Run

Mode of usage:

```bash
bash flop_launcher.sh [-d data folder] [-e config set][-p pvalue threshold] [-f n DE genes threshold] [-t] [-h]
```

FLOP has several ways of personalization. These are all possible input parameters:

- -d: data folder, containing the subfolders with the datasets to be analyzed
- -e: config set, either 'desktop' or 'cluster'
- -t: test mode, runs the pipeline with the test dataset and default parameters. Bear in mind that you still need to specify a config set with -e
- -p: pvalue threshold that the genes or functional terms need to pass in order to be considered significant for the Jaccard Index module. Default is 1 (no filtering).
- -f: Minimum number of significant genes per contrast. Only contrasts that have a minimum of n genes with a pvalue below 0.05 will be considered for enrichment analysis.
- -h: shows this help message

You can run FLOP with the minimal settings by using this command:

```bash
bash flop_launcher.sh [-d data_folder] [-e config_set]
```

Also, it is possible to run an example version of FLOP with a test dataset containing three contrasts from the PANACEA dataset via:

```bash
bash flop_launcher.sh -t [-e config_set]
```

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

The output consists in three or four different files:

## General results file:

A long format table that contains the functional scores of all different prior knowledge sources and per pipeline, biological contrast, filtering status, subset (if applicable) and t-value or log fold change.

## Rank correlation results:

A long format table that contains the spearman rank correlation and coocurrence scores per comparison, parameter, filtering status, biological context and prior knowledge source.

## Top/bottom features overlap results

A long format table that contains the Jaccard index and similarity scores per pipeline comparison, parameter, filtering status, biological context and prior knowledge source. The number of top and bottom functional categories included vary between the different prior knowledge source:

- Dorothea: 15 top and 15 bottom
- MSigDB hallmarks: top 5 and bottom 5
- PROGENy: top 3 and bottom 3
