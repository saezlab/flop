#!/bin/bash
#SBATCH --job-name FLOP_CCLE
#SBATCH --time 20-24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=victor.paton@embl.de 

echo "I, ${USER}, am running on host ${HOSTNAME}"
 
echo "STARTING JOB"
 
# Load R module
source ${HOME}/.bashrc
conda activate flop

# Executing an R script
bash /net/data.isilon/ag-saez/bq_vpaton/flop/flop_launcher.sh -d /net/data.isilon/ag-saez/bq_vpaton/flop/data/data_ccle/ -e cluster -f 30 -p 0.05

echo "FINISHED JOB"