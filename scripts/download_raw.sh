#!/bin/bash
#SBATCH --job-name=download_raw
#SBATCH --output=logs/download_raw_%j.log
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=02:00:00

source "$HOME/emilis.sh"

OUT="$HOME/HW3/data/raw"
mkdir -p "$OUT"

while read -r sample; do
    echo "sample: $sample"
    fasterq-dump "$sample" -O "$OUT" --split-files --threads 2
    gzip -f "$OUT/$sample"*.fastq
done < "$1"
