#!/bin/bash
banner=$"
##########################################################################
 Welcome to
.------..------..------..------.
|F.--. ||L.--. ||O.--. ||P.--. |
| :(): || :/\: || :/\: || :/\: |
| ()() || (__) || :\/: || (__) |
| '--'F|| '--'L|| '--'O|| '--'P|
\`------'\`------'\`------'\`------'

 The FunctionaL Omics Preprocessing platform is
 a workflow meant to evaluate the impact of different 
 normalization and differential expression tools on the 
 resulting functional space, in the context of bulk RNA-seq data.
 
 ##########################################################################
 "

# Change working directory to the directory where the script is located
cd -P -- "$(dirname -- "$0")"

help_txt="
Usage: flop_launcher.sh [-d data_folder] [-e config_set] [-t] [-h]

Argument list:
        -d: data folder, containing the subfolders with the datasets to be analyzed
        -e: config set, either 'desktop' or 'cluster'
        -t: test mode, runs the pipeline with the test dataset and default parameters. Bear in mind that you still need to specify a config set with -e
        -p: pvalue threshold that the genes or functional terms need to pass in order to be considered significant for the Jaccard Index module. Default is 1 (no filtering).
        -f: Minimum number of significant genes per contrast. Only contrasts that have a minimum of n genes with a pvalue below 0.05 will be considered for enrichment analysis.
        -h: shows this help message

For more information, please refer to the README file.
"

error_func () {
  echo "script usage: flop_launcher.sh [-d data_folder] [-e] [-a somevalue]" >&2
  echo "Try 'flop_launcher.sh -h' for more information." >&2
  exit 1
}

testdata_downloader () {
  curl -C - -O https://zenodo.org/record/8272401/files/test_data.zip?download=1
  mkdir ./test_data/
  mkdir ./test_data/test/
  unzip -n test_data.zip -d ./test_data/test/ 
}

while getopts 'd:e:f:p:ht' OPTION; do
  case "$OPTION" in
    d)
      data_folder="$OPTARG"
      ;;
    e)
      config_set="$OPTARG"
      ;;
    h)
      echo "$help_txt"
      exit
      ;;
    t)
      testdata_downloader
      data_folder="./test_data/"
      n_thresh=30
      p_thresh=0.05
      echo "##TEST MODE##"
      ;;
    f)
      n_thresh="$OPTARG"
      ;;
    p)
      p_thresh="$OPTARG"
      ;;
    ?)
      error_func
      ;;
  esac 
done
if [ $OPTIND -eq 1 ]; then error_func; fi
if [ -z $data_folder ]; then error_func; fi
if [ -z $n_thresh ]; then n_thresh=0; fi
if [ -z $config_set ]; then config_set='desktop'; fi
if [ -z $p_thresh ]; then p_thresh=1; fi
shift "$(($OPTIND -1))"

suffix="/" # add a slash to the end of the data folder if it is not already there
if [[ $data_folder == *$suffix ]]; then
  data_folder+=""
else
  data_folder+="/"
fi

parent_folder=$(dirname $data_folder)
num_dirs=$(ls -ld "$data_folder"* | awk '{print $NF}' | rev | cut -d "/" -f1 | rev | cut -d "/" -f3 | wc -l)
name_datasets=$(ls -ld "$data_folder"* | awk '{print $NF}' | rev | cut -d "/" -f1 | rev | sort | uniq | tr '\n' ' ')

# Ask if config is correct, if not, exit


echo "${banner}"
echo "
##SETTINGS##
Running option: $config_set
Data folder: $data_folder
Number of subsets found: $num_dirs
Datasets found: $name_datasets
Minimum number of significant genes per contrast: $n_thresh
Pvalue cutoff for the top-bottom overlapping module: $p_thresh"

# Run flop_benchmark


if [ $config_set = "desktop" ]; then
        echo "Running FLOP on a desktop computer..."
        nextflow -C flop.config run flop_desktop.nf -profile standard -resume --data_folder "$data_folder" --parent_folder "$parent_folder" --ngenes_threshold "$n_thresh" --pval_threshold "$p_thresh"
elif [ $config_set = "cluster" ]; then
        echo "Running FLOP on a slurm-controlled cluster..."
        nextflow -C flop.config run flop.nf -profile cluster --data_folder "$data_folder" --parent_folder "$parent_folder" --ngenes_threshold "$n_thresh" --pval_threshold "$p_thresh"
else
        echo "Valid options: desktop, cluster"
        error_func
fi


if [ $? -eq 0 ] 
then 
  echo "Analysis completed successfully! Your results are in $parent_folder/flop_results"
  exit 0
else 
  echo "The execution finished unsuccessfully. Please check the log file for more information." >&2 
  exit 1
fi



