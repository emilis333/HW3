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
        expand("results/rseqc/{sample}.read_distribution.txt", sample=SAMPLES),
        expand("results/rseqc/{sample}.geneBodyCoverage.txt", sample=SAMPLES),
        expand("results/rseqc/{sample}.inner_distance_freq.txt", sample=SAMPLES),
        "results/deeptools/correlation_heatmap.png",
        "results/deeptools/pca.png",
        "results/multiqc/multiqc_report.html",

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
    log:
        "logs/gtf_to_bed/gtf_to_bed.log"
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 30
    conda:
        "env/rseqc.yaml"
    shell:
        """
        gtfToGenePred {input} {output}.genePred 2> {log}
        genePredToBed {output}.genePred {output} 2>> {log}
        rm {output}.genePred
        """

rule bam_stat: #a. alignment statistics
    input:
        "results/hisat2/{sample}.sorted.bam"
    output:
        "results/rseqc/{sample}.bam_stat.txt"
    log:
        "logs/bam_stat/{sample}.log"
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 60
    conda:
        "env/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"

rule infer_experiment: #b. strandedness
    input:
        bam = "results/hisat2/{sample}.sorted.bam",
        bed = config["bed"]
    output:
        "results/rseqc/{sample}.infer_experiment.txt"
    log:
        "logs/infer_experiment/{sample}.log"
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 60
    conda:
        "env/rseqc.yaml"
    shell:
        "infer_experiment.py -i {input.bam} -r {input.bed} > {output} 2> {log}"

rule read_distribution: #c. read distribution
    input:
        bam = "results/hisat2/{sample}.sorted.bam",
        bed = config["bed"]
    output:
        "results/rseqc/{sample}.read_distribution.txt"
    log:
        "logs/read_distribution/{sample}.log"
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 120
    conda:
        "env/rseqc.yaml"
    shell:
        "read_distribution.py -i {input.bam} -r {input.bed} > {output} 2> {log}"

rule gene_body_coverage: #d. gene body coverage
    input:
        bam = "results/hisat2/{sample}.sorted.bam",
        bai = "results/hisat2/{sample}.sorted.bam.bai",
        bed = config["bed_subset"]
    output:
        "results/rseqc/{sample}.geneBodyCoverage.txt"
    params:
        prefix = "results/rseqc/{sample}"
    log:
        "logs/gene_body_coverage/{sample}.log"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 240
    conda:
        "env/rseqc.yaml"
    shell:
        "geneBody_coverage.py -i {input.bam} -r {input.bed} -o {params.prefix} > {log} 2>&1"

rule bed_subset: #every 100th
    input:
        config["bed"]
    output:
        config["bed_subset"]
    threads: 1
    resources:
        mem_mb = 1000,
        runtime = 10
    shell:
        "awk 'NR % 100 == 0' {input} > {output}"

rule inner_distance: #e. insert size/fragment length distrubution
    input:
        bam = "results/hisat2/{sample}.sorted.bam",
        bed = config["bed"]
    output:
        "results/rseqc/{sample}.inner_distance_freq.txt"
    params:
        prefix = "results/rseqc/{sample}"
    log:
        "logs/inner_distance/{sample}.log"
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 120
    conda:
        "env/rseqc.yaml"
    shell:
        "inner_distance.py -i {input.bam} -r {input.bed} -o {params.prefix} > {log} 2>&1"

rule multibamsummary:
    input:
        bams = expand("results/hisat2/{sample}.sorted.bam", sample=SAMPLES),
        bais = expand("results/hisat2/{sample}.sorted.bam.bai", sample=SAMPLES)
    output:
        "results/deeptools/readCounts.npz"
    log:
        "logs/deeptools/multibamsummary.log"
    threads: 8
    resources:
        mem_mb = 8000,
        runtime = 120
    conda:
        "env/deeptools.yaml"
    shell:
        "multiBamSummary bins --bamfiles {input.bams} -p {threads} -o {output} > {log} 2>&1"

rule plot_correlation:
    input:
        "results/deeptools/readCounts.npz"
    output:
        plot = "results/deeptools/correlation_heatmap.png",
        tab  = "results/deeptools/correlation_matrix.tab"
    log:
        "logs/deeptools/plot_correlation.log"
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 30
    conda:
        "env/deeptools.yaml"
    shell:
        """
        plotCorrelation -in {input} --corMethod spearman --whatToPlot heatmap \
          --skipZeros --plotNumbers \
          -o {output.plot} --outFileCorMatrix {output.tab} > {log} 2>&1
        """

rule plot_pca:
    input:
        "results/deeptools/readCounts.npz"
    output:
        "results/deeptools/pca.png"
    log:
        "logs/deeptools/plot_pca.log"
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 30
    conda:
        "env/deeptools.yaml"
    shell:
        "plotPCA -in {input} -o {output} > {log} 2>&1"

rule multiqc:
    input:
        expand("results/fastqc_raw/{sample}_{read}_fastqc.zip",
               sample=SAMPLES, read=["1", "2"]),
        expand("results/fastqc_trimmed/{sample}_{read}_fastqc.zip",
               sample=SAMPLES, read=["1", "2"]),
        "results/featurecounts/counts.txt",
        expand("results/kallisto/{sample}/run_info.json", sample=SAMPLES),
        expand("results/rseqc/{sample}.{suffix}", sample=SAMPLES,
               suffix=["bam_stat.txt", "infer_experiment.txt", "read_distribution.txt",
                       "geneBodyCoverage.txt", "inner_distance_freq.txt"]),
        "results/deeptools/correlation_matrix.tab",
        "results/deeptools/pca.png"
    output:
        "results/multiqc/multiqc_report.html"
    log:
        "logs/multiqc/multiqc.log"
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 30
    conda:
        "env/multiqc.yaml"
    shell:
        "multiqc results logs -o results/multiqc -n multiqc_report --force > {log} 2>&1"
