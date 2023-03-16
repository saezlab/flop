#!/bin/bash
source "$HOME/.bashrc"
# Change working directory to the directory where the script is located
cd -P -- "$(dirname -- "$0")"
echo $PWD


# Description: Install Nextflow
# Check if the directory Nextflow is already installed; if it is, skip step, if not, install it
if [ -f "nextflow" ]; then
        echo "Nextflow already installed, skipping installation"
else
        echo "Nextflow not installed, installing it"
        curl -s https://get.nextflow.io | bash
        chmod +x nextflow
        echo "Nextflow installed successfully"
fi


# Description: Install conda environment
conda init bash
# source $CONDA_PREFIX/etc/profile.d/conda.s./Data
# Check if conda env named flop_benchmark exists; if it does, skip step, if not, create it
if [ -d "$CONDA_PREFIX/envs/flop_benchmark" ]; then
        echo "Conda environment flop_benchmark already exists, skipping creation"
else
        echo "Conda environment flop_benchmark does not exist, creating it"
        conda env create -f scripts/config_env.yaml
        echo "Dependencies installed successfully" 
fi
conda activate flop_benchmark

# Description: Run flop_benchmark
echo '
##########################################################################
 Welcome to
 _______ _______ _______ _______ _______ _______ _______ _______
|  _______ ___     _______ _______                              |
| |   _   |   |   |   _   |   _   |                             |
| |.  1___|.  |   |.  |   |.  1   |                             |
| |.  __) |.  |___|.  |   |.  ____|                             |
| |:  |   |:  1   |:  1   |:  |                                 |
| |::.|   |::.. . |::.. . |::.|                                 |
| `---'\''   `-------`-------`---'\''                                 |
|  __                     __                        __          |
| |  |--.-----.-----.----|  |--.--------.---.-.----|  |--.      |
| |  _  |  -__|     |  __|     |        |  _  |   _|    <       |
| |_____|_____|__|__|____|__|__|__|__|__|___._|__| |__|__|      |
|_______ _______ _______ _______ _______ _______ _______ _______|

 The FunctionaL Omics Preprocessing benchmarking platform is
 a workflow meant to benchmark the impact of different 
 normalization and differential expression tools on the 
 resulting functional space, in the context of bulk RNA-seq data.
 
 ##########################################################################
 '

read -p "Please specify your folder containing the data to be analysed: " -e data_folder

parent_folder=$(dirname $data_folder)

# num_dirs=$(ls -l "$data_folder" | grep -c ^d)
num_dirs=$(find "$data_folder" -mindepth 1 -maxdepth 1 -not -empty -type d -printf '%f\n' | wc -l)
echo "Number of subsets found: $num_dirs"

name_datasets=$(find "$data_folder" -mindepth 1 -maxdepth 1 -not -empty -type d -printf '%f\n' | cut -d"_" -f1 | sort | uniq | tr '\n' ' ')
echo "Datasets found: $name_datasets"

echo "
Current options to run flop_benchmark are:
    1 - Run flop_benchmark on a desktop computer
    2 - Run flop_benchmark on a slurm-controlled cluster
    
Please select your option: "
read option

datasets_array=($name_datasets)
echo "Please select which of the following datasets will be included in the Rand Index analysis.
Leave blank if none:"
for index in ${!datasets_array[@]}; do
    echo $((index)) - "${datasets_array[index]}"
done
read perturbation

for index in $perturbation; do
    perturbation_array+=${datasets_array[$index]}" "
done

# Ask if config is correct, if not, exit
echo "
You have selected option $option.
Data folder is $data_folder.
Number of subsets found: $num_dirs
Datasets found: $name_datasets
Perturbational atasets included in the Rand Index analysis: $perturbation_array

Proceed? (y/n): "
read answer

if [ $answer != "y" ]; then
        echo "Aborting, exiting"
        exit
fi

# Run flop_benchmark
if [ $option -eq 1 ]; then
        echo "Running flop_benchmark on a desktop computer"
        ./nextflow -C bq_slurm.config run flop_benchmark.nf -profile standard --data_folder "$data_folder" --parent_folder "$parent_folder" --perturbation "$perturbation_array"
elif [ $option -eq 2 ]; then
        echo "Running flop_benchmark on a slurm-controlled cluster"
        ./nextflow -C bq_slurm.config run flop_benchmark.nf -profile cluster --data_folder "$data_folder" --parent_folder "$parent_folder" --perturbation "$perturbation_array"
else
        echo "Invalid option, aborting"
        exit
fi

echo "Analysis completed! Your results are in $parent_folder/flop_benchmark_results"

