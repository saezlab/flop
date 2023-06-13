#!/bin/bash
#SBATCH --job-name FLOP
#SBATCH --mem 150GB
#SBATCH --time 20-24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=victor.paton@embl.de 

echo "I, ${USER}, am running on host ${HOSTNAME}"
 
echo "STARTING JOB"
 
# Load R module
source ${HOME}/.bashrc
conda activate flop

# Executing an R script
bash /net/data.isilon/ag-saez/bq_vpaton/flop/flop_launcher.sh -d /net/data.isilon/ag-saez/bq_vpaton/flop/data/ -e cluster -r GSE186341 -k 35 -b range -f 50

echo "FINISHED JOB"