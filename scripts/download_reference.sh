#!/bin/bash
#SBATCH --output=logs/download_reference_%j.log
#SBATCH --mem=1G
#SBATCH --time=02:00:00

OUT="$HOME/HW3/data/reference"
BASE=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50
mkdir -p "$OUT"

wget -c -nv -P "$OUT" $BASE/GRCh38.primary_assembly.genome.fa.gz
wget -c -nv -P "$OUT" $BASE/gencode.v50.primary_assembly.annotation.gtf.gz
wget -c -nv -P "$OUT" $BASE/gencode.v50.transcripts.fa.gz

gunzip "$OUT"/*.gz
