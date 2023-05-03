<img src="https://github.com/saezlab/flop/blob/main/logo.png?raw=1" align="right" width="250" class="no-scaled-link" />

# FLOP: FunctionaL Omics Preprocessing platform 

## Description

FLOP is a unified workflow that analyses bulk RNA-seq counts using multiple combinations of filtering, normalisation and differential expression methods. It then evaluates the differences in the functional space between the different combinations of methods. 
FLOP is composed of six modules. The first three perform filtering, normalisation, differential expression analysis and functional analysis using different combination of methods or "pipelines". The remaining modules assess the differences in the results of the different pipelines. 

<p align="center" width="100%">
<img src="https://github.com/saezlab/flop/blob/main/man/images/graphicalabstractFLOP.svg" align="center" width="750">
</p>

## Installation

To install FLOP, you first will need to download the files from our GitHub repository:

```bash
git clone https://github.com/saezlab/flop
```

To run FLOP, you need to have conda installed in your computer. Please check [this link](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to learn how to install conda.

This command will install the necessary dependencies inside an environment.

```bash
conda env create -f ./scripts/flop_env.yaml
conda activate flop
```

Once installed, you are ready to run FLOP!

## Run

Mode of usage:

```bash
bash flop_launcher.sh [-d data_folder] [-e config_set] [-r perturbation_array] [-k k_val] [-b k_type] [-t] [-h]
```

FLOP has several ways of personalization. These are all possible input parameters:

- -d: data folder, containing the subfolders with the datasets to be analysed. Mandatory
- -e: config set, either 'desktop' or 'cluster’. Mandatory
- -r: perturbation array, a list of perturbational datasets to be included in the Rand Index analysis. Must be used jointly with -k and -b
- -t: test mode, runs the pipeline with the test dataset and the desktop config set. 
- -k: k value, the number of clusters to be used in the Rand Index analysis. Must be used jointly with -r
- -b: k value calculation, either 'range' or 'single'. Must be used jointly with -r
- -h: shows a help message

You can run FLOP with the minimal settings by using this command:

```bash
bash flop_launcher.sh [-d data_folder] [-e config_set]
```

Also, it is possible to run an example version of FLOP with a test dataset via:

```bash
bash flop_launcher.sh -t
```

## Input

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

### Count table

This tab-separated table should contain the gene counts of the dataset. It should have the following name pattern: {Dataset ID}__countdata.tsv 

The data should contain a column named gene_symbol, which contains the gene symbols, and the following columns should be the samples that will be included in the study. Be aware that the sample names should respect R’s naming regulations (see the command make.names). 

### Metadata table

This tab-separated table should at least contain the the group to which the samples belong. It should have the following name pattern: {Dataset ID}__metadata.tsv 

This table contains a column named sample_ID, with the names of the samples (again, they should respect R’s naming guidelines) and a column named group, which contains the group every sample belongs to. Other columns can be present, but they will not be considered in the analysis.

### Contrast table

This tab-separated table should contain the contrasts that should be calculated in the analysis. It should have the following name pattern: {Dataset ID}__contrast.tsv 

The table should contain two columns, group1 and group2, with the desired contrasts.

This table is optional. If included, FLOP will only calculate the specified contrast. If this table is not included, it will calculate every pairwise comparison between the different groups (be careful with choosing this option when analysing very large datasets, since it will greatly increase the execution time!).

## Output

The output consists in three or four different files:

### General results file:

A long format table that contains the functional scores of all different prior knowledge sources and per pipeline, biological contrast, filtering status, subset (if applicable) and t-value or log fold change.

### Rank correlation results:

A long format table that contains the spearman rank correlation scores per comparison, parameter, filtering status, biological context and prior knowledge source. The comparisons are done at the pipeline level, normalisation level and differential expression method level.

### Jaccard index results

A long format table that contains the Jaccard index scores per pipeline comparison, parameter, filtering status, biological context and prior knowledge source. The number of top and bottom functional categories included vary between the different prior knowledge source:

- Dorothea: 15 top and 15 bottom
- MSigDB hallmarks: top 5 and bottom 5
- PROGENy: top 3 and bottom 3

### Rand index analysis

A long format table that contains the Rand index analysis for a specified number of K values, per pipeline comparison, filtering status, biological context and prior knowledge source.

Since this analysis is only informative when there might be a ground truth (such a specific number of cell lines, treatments, etc.) that makes the clustering of the samples important for a study, we implemented this analysis as optional inside the FLOP architecture. You can select the datasets, the K value, and if this value is unique or the maximum of a range, that will be included in this analysis during the initial configuration.
