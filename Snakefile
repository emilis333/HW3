import pandas as pd

configfile: "config/config.yaml"

samples = pd.read_csv(config["samples"], sep="\t", index_col="sample")
SAMPLES = list(samples.index)

wildcard_constraints:
    read = "[12]"

rule all:
    input:
        expand("results/fastqc_raw/{sample}_{read}_fastqc.html",
               sample=SAMPLES, read=["1", "2"]),
        expand("results/fastqc_trimmed/{sample}_{read}_fastqc.html",
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

rule trim:
    input:
        r1 = lambda wc: samples.loc[wc.sample, "fq1"],
        r2 = lambda wc: samples.loc[wc.sample, "fq2"]
    output:
        r1   = "results/trimmed/{sample}_1.fastq.gz",
        r2   = "results/trimmed/{sample}_2.fastq.gz",
        json = "results/trimmed/{sample}_fastp.json",
        html = "results/trimmed/{sample}_fastp.html"
    log:
        "logs/trim/{sample}.log"
    threads: 4
    resources:
        mem_mb = 4000,
        runtime = 120
    conda:
        "env/fastp.yaml"
    shell:
        """
        fastp --in1 {input.r1} --in2 {input.r2} \
              --out1 {output.r1} --out2 {output.r2} \
              --thread {threads} --detect_adapter_for_pe \
              --qualified_quality_phred 20 --length_required 30 \
              --json {output.json} --html {output.html} 2> {log}
        """

rule fastqc_trimmed:
    input:
        "results/trimmed/{sample}_{read}.fastq.gz"
    output:
        html = "results/fastqc_trimmed/{sample}_{read}_fastqc.html",
        zip  = "results/fastqc_trimmed/{sample}_{read}_fastqc.zip"
    log:
        "logs/fastqc_trimmed/{sample}_{read}.log"
    threads: 2
    resources:
        mem_mb = 2000,
        runtime = 60
    conda:
        "env/fastqc.yaml"
    shell:
        "fastqc -t {threads} -o results/fastqc_trimmed {input} > {log} 2>&1"
