#!/bin/sh

#SBATCH --mem=100000 
#SBATCH --cpu-freq=High 

source $HOME/.bashrc

conda activate flop2

Rscript ccle_proc.R
