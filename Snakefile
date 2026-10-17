configfile: "config.yaml"

import csv
import os
import re

# Reject unknown configuration keys
KNOWN = {
    "samplesheet", "ref", "region", "outdir",
    "platform", "executor", "containers"
}

unknown = sorted(set(config) - KNOWN)
if unknown:
    raise ValueError(f"unknown config key(s): {', '.join(unknown)}")

# Pipeline metadata
PIPELINE = {"name": "variant-call", "version": "1.0.0"}

OUT = config["outdir"]
BIN = os.path.join(workflow.basedir, "bin")

IMAGES = {
    name: f"docker://{image}"
    for name, image in config["containers"].items()
}

# Reference genome and associated indexes
REF = config["ref"]
REF_FAI = REF + ".fai"
REF_BWA = [
    REF + ext
    for ext in (".amb", ".ann", ".bwt", ".pac", ".sa")
]
REF_DICT = os.path.splitext(REF)[0] + ".dict"

# Read samplesheet
SHEET = config["samplesheet"]

with open(SHEET) as fh:
    ROWS = list(csv.DictReader(fh))

SAMPLES = [row["sample_id"] for row in ROWS]

PAIRED = [
    row["sample_id"]
    for row in ROWS
    if row["library_type"] == "paired"
]

SINGLE = [
    row["sample_id"]
    for row in ROWS
    if row["library_type"] == "single"
]

# Resolve FASTQ paths relative to samplesheet directory
def resolve(path):
    return path if os.path.isabs(path) else os.path.join(os.path.dirname(SHEET), path)

READS = {
    row["sample_id"]: [resolve(row["r1_fastq"])]
    + ([resolve(row["r2_fastq"])] if row["library_type"] == "paired" else [])
    for row in ROWS
}

# Match only samples of the specified library type
def one_of(names):
    return "|".join(map(re.escape, names)) if names else "(?!)"

# Return the appropriate trimmed FASTQs for each sample
def trimmed(wildcards):
    s = wildcards.sample
    if s in PAIRED:
        return [
            f"work/trim/{s}_R1.trim.fastq.gz",
            f"work/trim/{s}_R2.trim.fastq.gz"
        ]
    return [f"work/trim/{s}_R1.trim.fastq.gz"]


# Final outputs required from the pipeline
rule all:
    input:
        f"{OUT}/cohort.filtered.vcf.gz",
        f"{OUT}/variants.tsv",
        f"{OUT}/multiqc_report.html",
        f"{OUT}/samples.tsv",
        f"{OUT}/manifest.json"
rule validate:
    input:
        sheet=SHEET,
        fastqs=[fq for reads in READS.values() for fq in reads],
        ref=REF,
        ref_index=[REF_FAI] + REF_BWA,
        ref_dict=REF_DICT
    output:
        "work/validate/samplesheet.csv"
    log:
        "logs/validate.log"
    container:
        IMAGES["tools"]
    params:
        region=config["region"]
    shell:
        """
        ({BIN}/validate_samplesheet.sh {input.sheet} {input.ref} '{params.region}' &&
         cp {input.sheet} {output}) > {log} 2>&1
        """
rule fastqc:
    input:
        checked=rules.validate.output,
        read=lambda wc: READS[wc.sample][int(wc.read) - 1]
    output:
        html="work/fastqc/{sample}_R{read}_fastqc.html",
        zip="work/fastqc/{sample}_R{read}_fastqc.zip"
    log:
        "logs/fastqc/{sample}_R{read}.log"
    wildcard_constraints:
        read="[12]"
    params:
        outdir=lambda wc: f"work/fastqc/raw_{wc.sample}_R{wc.read}"
    container:
        IMAGES["fastqc"]
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -q -o {params.outdir} {input.read} > {log} 2>&1

        mv {params.outdir}/*_fastqc.html {output.html}
        mv {params.outdir}/*_fastqc.zip {output.zip}
        rmdir {params.outdir}
        """
rule fastp_paired:
    input:
        checked=rules.validate.output,
        r1=lambda wc: READS[wc.sample][0],
        r2=lambda wc: READS[wc.sample][1]
    output:
        r1="work/trim/{sample}_R1.trim.fastq.gz",
        r2="work/trim/{sample}_R2.trim.fastq.gz",
        html="work/trim/{sample}.fastp.html",
        json="work/trim/{sample}.fastp.json"
    log:
        "logs/fastp/{sample}.log"
    wildcard_constraints:
        sample=one_of(PAIRED)
    threads: 4
    container:
        IMAGES["fastp"]
    shell:
        """
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            --html {output.html} \
            --json {output.json} \
            --thread {threads} \
            > {log} 2>&1
        """
rule fastp_single:
    input:
        checked=rules.validate.output,
        r1=lambda wc: READS[wc.sample][0]
    output:
        r1="work/trim/{sample}_R1.trim.fastq.gz",
        html="work/trim/{sample}.fastp.html",
        json="work/trim/{sample}.fastp.json"
    log:
        "logs/fastp/{sample}.log"
    wildcard_constraints:
        sample=one_of(SINGLE)
    threads: 4
    container:
        IMAGES["fastp"]
    shell:
        """
        fastp \
            -i {input.r1} \
            -o {output.r1} \
            --html {output.html} \
            --json {output.json} \
            --thread {threads} \
            > {log} 2>&1
        """
rule bwa_mem:
    input:
        reads=trimmed,
        ref=REF,
        fai=REF_FAI,
        bwa_index=REF_BWA,
        ref_dict=REF_DICT
    output:
        bam="work/alignment/{sample}.sorted.bam",
        bai="work/alignment/{sample}.sorted.bam.bai"
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 8
    container:
        IMAGES["bwa"]
    params:
        rg=lambda wc: f"@RG\\tID:{wc.sample}\\tSM:{wc.sample}\\tPL:ILLUMINA"
    shell:
        """
        bwa mem -t {threads} -R '{params.rg}' {input.ref} {input.reads} 2> {log} |
            samtools sort -@ {threads} -o {output.bam} - 2>> {log}

        samtools index {output.bam} {output.bai} 2>> {log}
        """
rule markduplicates:
    input:
        bam="work/alignment/{sample}.sorted.bam",
        bai="work/alignment/{sample}.sorted.bam.bai"
    output:
        bam="work/dedup/{sample}.dedup.bam",
        bai="work/dedup/{sample}.dedup.bam.bai",
        metrics="work/dedup/{sample}.metrics.txt"
    log:
        "logs/markduplicates/{sample}.log"
    threads: 4
    container:
        IMAGES["gatk"]
    shell:
        """
        gatk MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --CREATE_INDEX true \
            > {log} 2>&1

        mv work/dedup/{wildcards.sample}.dedup.bai {output.bai}
        """
rule haplotypecaller:
    input:
        bam="work/dedup/{sample}.dedup.bam",
        bai="work/dedup/{sample}.dedup.bam.bai",
        ref=REF,
        fai=REF_FAI,
        ref_dict=REF_DICT
    output:
        gvcf="work/gvcf/{sample}.g.vcf.gz",
        tbi="work/gvcf/{sample}.g.vcf.gz.tbi"
    log:
        "logs/haplotypecaller/{sample}.log"
    threads: 8
    container:
        IMAGES["gatk"]
    params:
        region=config["region"]
    shell:
        """
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.gvcf} \
            -ERC GVCF \
            -L {params.region} \
            --native-pair-hmm-threads {threads} \
            > {log} 2>&1
        """
rule joint_genotype:
    input:
        gvcfs=expand("work/gvcf/{sample}.g.vcf.gz", sample=SAMPLES),
        indexes=expand("work/gvcf/{sample}.g.vcf.gz.tbi", sample=SAMPLES),
        ref=REF,
        fai=REF_FAI,
        ref_dict=REF_DICT
    output:
        vcf="work/joint/cohort.raw.vcf.gz",
        tbi="work/joint/cohort.raw.vcf.gz.tbi"
    log:
        "logs/joint_genotype.log"
    threads: 4
    container:
        IMAGES["gatk"]
    params:
        variants=lambda wc, input: " ".join(
            f"-V {gvcf}" for gvcf in input.gvcfs
        ),
        region=config["region"]
    shell:
        """
        gatk CombineGVCFs \
            -R {input.ref} \
            {params.variants} \
            -O work/joint/cohort.combined.g.vcf.gz \
            -L {params.region} \
            > {log} 2>&1

        gatk GenotypeGVCFs \
            -R {input.ref} \
            -V work/joint/cohort.combined.g.vcf.gz \
            -O {output.vcf} \
            -L {params.region} \
            >> {log} 2>&1
        """
rule filter:
    input:
        vcf=rules.joint_genotype.output.vcf,
        tbi=rules.joint_genotype.output.tbi
    output:
        vcf=f"{OUT}/cohort.filtered.vcf.gz",
        tbi=f"{OUT}/cohort.filtered.vcf.gz.tbi",
        table=f"{OUT}/variants.tsv"
    log:
        "logs/filter.log"
    threads: 2
    container:
        IMAGES["gatk"]
    shell:
        r"""
        gatk VariantFiltration \
            -V {input.vcf} \
            -O {output.vcf} \
            --filter-name "LowQD" \
            --filter-expression "QD < 2.0" \
            --filter-name "HighFS" \
            --filter-expression "FS > 60.0" \
            --filter-name "LowMQ" \
            --filter-expression "MQ < 40.0" \
            > {log} 2>&1

   bcftools index -t -f {output.vcf} 2>> {log}

        printf 'chrom\tpos\tref\talt\tqual\tfilter\n' > {output.table}

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\n' \
            {output.vcf} >> {output.table} 2>> {log}
        """
rule multiqc:
    input:
        fastqc=(
            expand("work/fastqc/{sample}_R1_fastqc.zip", sample=SAMPLES)
            +
            expand("work/fastqc/{sample}_R2_fastqc.zip", sample=PAIRED)
        ),
        fastp=expand("work/trim/{sample}.fastp.json", sample=SAMPLES)
    output:
        html=f"{OUT}/multiqc_report.html"
    log:
        "logs/multiqc.log"
    threads: 1
    container:
        IMAGES["multiqc"]
    shell:
        """
        multiqc {input.fastqc} {input.fastp} \
            --filename multiqc_report.html \
            --outdir {OUT} \
            --force \
            > {log} 2>&1
        """
# Stage 9: Publish final results

PUBLISH_SHEET = rules.validate.output

PUBLISH_RESULTS = [
    rules.filter.output.vcf,
    rules.filter.output.tbi,
    rules.filter.output.table,
    rules.multiqc.output.html
]

include: "rules/publish.smk"