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
Usage: flop_launcher.sh [-d data_folder] [-e config_set] [-r perturbation_array]

Argument list:
        -d: data folder, containing the subfolders with the datasets to be analyzed
        -e: config set, either 'desktop' or 'cluster'
        -r: perturbation array, a list of perturbational datasets to be included in the Rand Index analysis. Optional
        -h: shows this help message

For more information, please refer to the README file.
"

error_func () {
  echo "script usage: flop_launcher.sh [-d data_folder] [-e] [-a somevalue]" >&2
  echo "Try 'flop_launcher.sh -h' for more information." >&2
  exit 1
}

while getopts 'd:e:r:h' OPTION; do
  case "$OPTION" in
    d)
      data_folder="$OPTARG"
      ;;
    e)
      config_set="$OPTARG"
      echo "The value provided is $OPTARG"
      ;;
    r)
      perturbation_array="$OPTARG"
      ;;
    h)
      echo "$help_txt"
      exit
      ;;
    ?)
      error_func
      ;;
  esac
done
if [ $OPTIND -eq 1 ]; then error_func; fi
shift "$(($OPTIND -1))"

suffix="/" # add a slash to the end of the data folder if it is not already there
if [[ $data_folder == *$suffix ]]; then
  data_folder+=""
else
  data_folder+="/"
fi

# Description: Run flop_benchmark
echo "${banner}"

parent_folder=$(dirname $data_folder)
num_dirs=$(ls -ld "$data_folder"* | awk '{print $NF}' | rev | cut -d "/" -f1 | rev | cut -d "/" -f3 | wc -l)
name_datasets=$(ls -ld "$data_folder"* | awk '{print $NF}' | rev | cut -d "/" -f1 | rev | sort | uniq | tr '\n' ' ')

# Ask if config is correct, if not, exit
echo "
##SETTINGS##
Running option: $config_set
Data folder: $data_folder
Number of subsets found: $num_dirs
Datasets found: $name_datasets
Perturbational datasets included in the Rand Index analysis: $perturbation_array

Proceed? (y/n): "
read answer

if [ $answer != "y" ]; then
        echo "Aborting..."
        exit
fi

# Run flop_benchmark
if [ $config_set == "desktop" ]; then
        echo "Running FLOP on a desktop computer..."
        nextflow -C bq_slurm.config run flop.nf -profile standard -resume --data_folder "$data_folder" --parent_folder "$parent_folder" --perturbation "$perturbation_array"
elif [ $config_set == "cluster" ]; then
        echo "Running FLOP on a slurm-controlled cluster..."
        nextflow -C bq_slurm.config run flop.nf -profile cluster -resume --data_folder "$data_folder" --parent_folder "$parent_folder" --perturbation "$perturbation_array"
else
        echo "Valid options: desktop, cluster"
        exit 1
fi

echo "Analysis completed! Your results are in $parent_folder/flop_results"

