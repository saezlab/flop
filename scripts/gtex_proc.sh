#!/bin/bash
#SBATCH --job-name GTEX_proc
#SBATCH --nodes 1
#SBATCH --ntasks 16
#SBATCH --mem 128GB
#SBATCH --time 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=victor.paton@embl.de 

source $HOME/.bashrc

conda activate flop

Rscript /net/data.isilon/ag-saez/bq_vpaton/flop/scripts/gtex_proc.R