import pandas as pd

configfile: "config/config.yaml"

samples = pd.read_csv(config["samples"], sep="\t", index_col="sample")
SAMPLES = list(samples.index)

wildcard_constraints:
    read = "[12]"

rule all:
    input:
        expand("results/fastqc_raw/{sample}_{read}_fastqc.html",
               sample=SAMPLES, read=["1", "2"])

rule fastqc_raw:
    input:
        lambda wc: samples.loc[wc.sample, "fq" + wc.read]
    output:
        html = "results/fastqc_raw/{sample}_{read}_fastqc.html",
        zip  = "results/fastqc_raw/{sample}_{read}_fastqc.zip"
    log:
        "logs/fastqc_raw/{sample}_{read}.log"
    threads: 2
    resources:
        mem_mb = 2000,
        runtime = 60
    conda:
        "env/fastqc.yaml"
    shell:
        """
        outdir=$(dirname {output.html})
        base=$(basename {input} .fastq.gz)
        fastqc -t {threads} -o $outdir {input} > {log} 2>&1
        mv $outdir/${{base}}_fastqc.html {output.html}
        mv $outdir/${{base}}_fastqc.zip {output.zip}
        """
