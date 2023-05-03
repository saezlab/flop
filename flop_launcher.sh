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
Usage: flop_launcher.sh [-d data_folder] [-e config_set] [-r perturbation_array] [-k k_val] [-b k_type] [-t] [-h]

Argument list:
        -d: data folder, containing the subfolders with the datasets to be analyzed
        -e: config set, either 'desktop' or 'cluster'
        -r: perturbation array, a list of perturbational datasets to be included in the Rand Index analysis.
        -k: k value, the number of clusters to be used in the Rand Index analysis.
        -b: k value calculation, either 'range' or 'single'.
        -t: test mode, runs the pipeline with the test dataset and the desktop config set
        -h: shows this help message

For more information, please refer to the README file.
"

error_func () {
  echo "script usage: flop_launcher.sh [-d data_folder] [-e] [-a somevalue]" >&2
  echo "Try 'flop_launcher.sh -h' for more information." >&2
  exit 1
}

testdata_downloader () {
  curl -C - -O https://filedn.eu/ld7S7VEWtgOf5uN0V7fbp84/test_data.zip
  unzip -n test_data.zip -d ./test_data
}

while getopts 'd:e:r:k:b:ht' OPTION; do
  case "$OPTION" in
    d)
      data_folder="$OPTARG"
      ;;
    e)
      config_set="$OPTARG"
      ;;
    r)
      perturbation_array="$OPTARG"
      ;;
    h)
      echo "$help_txt"
      exit
      ;;
    t)
      testdata_downloader
      data_folder="./test_data/"
      config_set="desktop"
      perturbation_array="test"
      k_val=3
      k_type="range"
      echo "##TEST MODE##"
      ;;
    k)
      k_val="$OPTARG"
      ;;
    b)
      k_type="$OPTARG"
      ;;
    ?)
      error_func
      ;;
  esac
done
if [ $OPTIND -eq 1 ]; then error_func; fi
if [ -z $k_val ] && [ ! -z $k_type ]; then error_func; fi
if [ ! -z $k_val ] && [ -z $k_type ]; then error_func; fi
if [ -z $data_folder ] || [ -z $config_set ]; then error_func; fi
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
Perturbational datasets included in the Rand Index analysis: $perturbation_array"
if [ ! -z $k_val ] && [ ! -z $k_type ]; then
  echo "k value(s): $k_val"
  echo "k value(s) calculation: $k_type"
fi

echo "Proceed? (y/n): "
read answer

if [ $answer != "y" ]; then
        echo "Aborting..."
        exit
fi

# check that k_

# Run flop_benchmark
if [ $config_set == "desktop" ]; then
        echo "Running FLOP on a desktop computer..."
        nextflow -C flop.config run flop.nf -profile standard -resume --data_folder "$data_folder" --parent_folder "$parent_folder" --perturbation "$perturbation_array" --k_val "$k_val" --k_type $k_type
elif [ $config_set == "cluster" ]; then
        echo "Running FLOP on a slurm-controlled cluster..."
        nextflow -C flop.config run flop.nf -profile cluster -resume --data_folder "$data_folder" --parent_folder "$parent_folder" --perturbation "$perturbation_array" --k_val "$k_val" --k_type $k_type
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


