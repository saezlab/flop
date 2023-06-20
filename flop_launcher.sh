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
        -t: test mode, runs the pipeline with the test dataset and default parameters. Bear in mind that you still need to specify a config set with -e
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
  curl -C - -O https://filedn.eu/ld7S7VEWtgOf5uN0V7fbp84/test_data.zip
  unzip -n test_data.zip -d ./test_data
}

while getopts 'd:e:r:k:b:f:ht' OPTION; do
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
      perturbation_array="test"
      k_val=3
      k_type="range"
      n_thresh=50
      echo "##TEST MODE##"
      ;;
    k)
      k_val="$OPTARG"
      ;;
    b)
      k_type="$OPTARG"
      ;;
    f)
      n_thresh="$OPTARG"
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
if [ -z $n_thresh ]; then n_thresh=0; fi
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
Perturbational datasets included in the Rand Index analysis: $perturbation_array
Minimum number of significant genes per contrast: $n_thresh"

if [ ! -z $k_val ] && [ ! -z $k_type ]; then
  echo "k value(s): $k_val"
  echo "k value(s) calculation: $k_type"
fi

# echo "Proceed? (y/n): "
# read answer

# if [ $answer != "y" ]; then
#         echo "Aborting..."
#         exit
# fi


# Run flop_benchmark

if [ -z $perturbation_array ]; then
  if [ $config_set = "desktop" ]; then
          echo "Running FLOP on a desktop computer... No cluster"
          nextflow -C flop.config run flop_desktop.nf -profile standard -resume --data_folder "$data_folder" --parent_folder "$parent_folder" --ngenes_threshold "$n_thresh"
  elif [ $config_set = "cluster" ]; then
          echo "Running FLOP on a slurm-controlled cluster... No cluster"
          nextflow -C flop.config run flop.nf -profile cluster --data_folder "$data_folder" --parent_folder "$parent_folder" --ngenes_threshold "$n_thresh"
  else
          echo "Valid options: desktop, cluster"
          error_func
  fi
else 
  if [ $config_set = "desktop" ]; then
          echo "Running FLOP on a desktop computer..."
          nextflow -C flop.config run flop_desktop.nf -profile standard -resume --data_folder "$data_folder" --parent_folder "$parent_folder" --perturbation "$perturbation_array" --k_val "$k_val" --k_type $k_type --ngenes_threshold "$n_thresh"
  elif [ $config_set = "cluster" ]; then
          echo "Running FLOP on a slurm-controlled cluster..."
          nextflow -C flop.config run flop.nf -profile cluster --data_folder "$data_folder" --parent_folder "$parent_folder" --perturbation "$perturbation_array" --k_val "$k_val" --k_type $k_type --ngenes_threshold "$n_thresh"
  else
          echo "Valid options: desktop, cluster"
          error_func
  fi
fi

if [ $? -eq 0 ] 
then 
  echo "Analysis completed successfully! Your results are in $parent_folder/flop_results"
  exit 0
else 
  echo "The execution finished unsuccessfully. Please check the log file for more information." >&2 
  exit 1
fi



