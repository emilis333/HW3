#!/bin/bash
#SBATCH --output=logs/download_raw_ena_%j.log
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=02:00:00

OUT="$HOME/HW3/data/raw"
mkdir -p "$OUT"

while read gsm; do
    curl -s "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=sample_alias=%22$gsm%22&fields=fastq_ftp&format=tsv" \
        | grep -o 'ftp[^;[:space:]]*fastq.gz' | xargs -I{} wget -c -P "$OUT" https://{}
done < "$HOME/HW3/data/assigned.txt"

