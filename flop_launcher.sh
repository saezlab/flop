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
Usage: flop_launcher.sh [-d data_folder] [-e config_set] [-t] [-h] [-m "plot_config"]

Argument list:
        -e: config set, either 'desktop' or 'cluster'. If none specified, it defaults to desktop
        -t: test mode, runs the pipeline with the test dataset and default parameters. Bear in mind that you still need to specify a config set with -e
        -s: Launches FLOP to run the analysis detailed in the accompanying study. It runs FLOP with CCLE, PANACEA and Reheat datasets. 
            It sets the pvalue threshold to 1 and the number of significant genes threshold to 30. For more information, please check the study. 
            It is strongly advised to run this setup within a HPC environment (config set = cluster).
        -m: plotter mode, generates a multidimensional scaling plot. The configuration string must be inputed as follows:

            -m '[file path for FLOP result (fullmerged)] [biological context] [resource] [DE metric (logFC or stat (t-value))]'

            The plot will be saved in the current directory.
        -h: shows this help message

For more information, please refer to the README file.
"

error_func () {
  echo "script usage: flop_launcher.sh [-d data_folder] [-e] [...]" >&2
  echo "Try 'flop_launcher.sh -h' for more information." >&2
  exit 1
}

testdata_downloader () {
  curl -C - -L -O https://zenodo.org/record/8272401/files/test_data.zip
  mkdir -p ./test_data/test/
  unzip -n test_data.zip -d ./test_data/test/ 
}

studydata_downloader () {
  curl -C - -L -O https://zenodo.org/record/8306225/files/flop_unproc_data.zip
  unzip -n flop_unproc_data.zip -d ./
}

pksource_downloader () {
  curl -C - -L -O https://zenodo.org/record/10980463/files/flop_pkresources_12042024.zip
  unzip -o flop_pkresources_12042024.zip -d ./scripts
}

plotter () {
  Rscript ./scripts/multid_plot.R -m "$1"
}

while getopts 'e:hstm:' OPTION; do
  case "$OPTION" in
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
    s)
      studydata_downloader
      pksource_downloader
      paper_mode=true
      data_folder="./flop_data/"
      n_thresh=30
      p_thresh=1
      ;;
    m)
      plot_config="$OPTARG"
      echo $plot_config
      plotter "$plot_config"
      echo "#Plot generated, file saved in $PWD#"
      exit
      ;;
    ?)
      error_func
      ;;
  esac 
done
if [ $OPTIND -eq 1 ]; then error_func; fi
if [ -z $paper_mode ]; then paper_mode=false; fi
if [ -z $config_set ]; then config_set='desktop'; fi
shift "$(($OPTIND -1))"

# Run flop_benchmark

if [ $config_set == "desktop" ] || [ $config_set == "cluster" ]; then
        nextflow -C flop.config run flop.nf -profile $config_set --data_folder "$data_folder" -params-file ./params_flop.json
else
        echo "Valid options: desktop, cluster"
        error_func
fi

if [ $paper_mode == true ]; then
  rm -r flop_pkresources_31082023.zip
  rm -r ./scripts/dc_resources
fi

