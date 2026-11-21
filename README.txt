Project folder must be ~/HW3. Commands should be run from project folder.
Conda environment `snakemake` with snakemake and pandas is required.

1. Download the reads from ENA (sample IDs must be in data/assigned.txt):
	sbatch scripts/download_raw_ena.sh

2. Download the reference and unpack it:
	sbatch scripts/download_reference.sh

3. Activate snakemake:
	conda activate snakemake

4. Run the workflow:
	sbatch scripts/run_snakemake.sh
