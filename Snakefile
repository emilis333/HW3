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
               sample=SAMPLES, read=["1", "2"]),
        "results/featurecounts/counts.txt",
        expand("results/kallisto/{sample}/abundance.tsv", sample=SAMPLES),
        expand("results/rseqc/{sample}.bam_stat.txt", sample=SAMPLES),
        expand("results/rseqc/{sample}.infer_experiment.txt", sample=SAMPLES),

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

rule hisat2_index:
    input:
        config["genome_fasta"]
    output:
        multiext(config["hisat2_index"], ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2",
                 ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2")
    params:
        prefix = config["hisat2_index"]
    log:
        "logs/hisat2_index/index.log"
    threads: 8
    resources:
        mem_mb = 32000,
        runtime = 240
    conda:
        "env/align.yaml"
    shell:
        "hisat2-build -p {threads} {input} {params.prefix} > {log} 2>&1"

rule hisat2_align:
    input:
        r1 = "results/trimmed/{sample}_1.fastq.gz",
        r2 = "results/trimmed/{sample}_2.fastq.gz",
        idx = multiext(config["hisat2_index"], ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2",
                       ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2")
    output:
        bam = "results/hisat2/{sample}.sorted.bam",
        bai = "results/hisat2/{sample}.sorted.bam.bai"
    params:
        prefix = config["hisat2_index"]
    log:
        "logs/hisat2/{sample}.log"
    threads: 8
    resources:
        mem_mb = 16000,
        runtime = 240
    conda:
        "env/align.yaml"
    shell:
        """
        hisat2 -p {threads} -x {params.prefix} -1 {input.r1} -2 {input.r2} 2> {log} \
          | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """

rule featurecounts:
    input:
        bams = expand("results/hisat2/{sample}.sorted.bam", sample=SAMPLES),
        gtf = config["gtf"]
    output:
        "results/featurecounts/counts.txt"
    log:
        "logs/featurecounts/featurecounts.log"
    threads: 8
    resources:
        mem_mb = 8000,
        runtime = 120
    conda:
        "env/align.yaml"
    shell:
        """
        featureCounts -T {threads} -p --countReadPairs \
          -t exon -g gene_id -a {input.gtf} -o {output} {input.bams} 2> {log}
        """

rule kallisto_index:
    input:
        config["transcripts"]
    output:
        config["kallisto_index"]
    log:
        "logs/kallisto_index/index.log"
    threads: 1
    resources:
        mem_mb = 32000,
        runtime = 120
    conda:
        "env/kallisto.yaml"
    shell:
        "kallisto index -i {output} {input} > {log} 2>&1"

rule kallisto_quant:
    input:
        r1 = "results/trimmed/{sample}_1.fastq.gz",
        r2 = "results/trimmed/{sample}_2.fastq.gz",
        idx = config["kallisto_index"]
    output:
        tsv = "results/kallisto/{sample}/abundance.tsv",
        info = "results/kallisto/{sample}/run_info.json"
    params:
        outdir = "results/kallisto/{sample}"
    log:
        "logs/kallisto/{sample}.log"
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 120
    conda:
        "env/kallisto.yaml"
    shell:
        "kallisto quant -i {input.idx} -o {params.outdir} -t {threads} {input.r1} {input.r2} > {log} 2>&1"

rule gtf_to_bed:
    input:
        config["gtf"]
    output:
        config["bed"]
    conda:
        "env/rseqc.yaml"
    log:
        "logs/gtf_to_bed/gtf_to_bed.log"
    shell:
        """
        gtfToGenePred {input} {output}.genePred 2> {log}
        genePredToBed {output}.genePred {output} 2>> {log}
        rm {output}.genePred
        """

rule bam_stat:
    input:
        "results/hisat2/{sample}.sorted.bam"
    output:
        "results/rseqc/{sample}.bam_stat.txt"
    log:
        "logs/bam_stat/{sample}.log"
    conda:
        "env/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"

rule infer_experiment:
    input:
        bam = "results/hisat2/{sample}.sorted.bam",
        bed = config["bed"]
    output:
        "results/rseqc/{sample}.infer_experiment.txt"
    log:
        "logs/infer_experiment/{sample}.log"
    conda:
        "env/rseqc.yaml"
    shell:
        "infer_experiment.py -i {input.bam} -r {input.bed} > {output} 2> {log}"
