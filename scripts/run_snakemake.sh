#!/bin/bash
#SBATCH --output=logs/snakemake_%j.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00

source "$HOME/emilis.sh"
conda activate snakemake
cd ~/HW3
snakemake --use-conda --cores 8
