# Snakefile

import os
import pandas as pd

# Configuration
configfile: "config.yaml"

# Load sample information
samples = pd.read_csv("data/samples.tsv", sep="\t").set_index("sample", drop=False)

# Define output directories
SNIPPY_DIR = "results/snippy"
BED_DIR = "results/bed_files"
INTERSECT_DIR = "results/intersect_230"

# Rule all
rule all:
    input:
        expand(f"{SNIPPY_DIR}/{{sample}}", sample=samples.index),
        BED_DIR,
        expand(f"{INTERSECT_DIR}/{{bed}}/{{sample}}_{{bed}}.vcf", 
               sample=samples.index, 
               bed=[b.replace('.bed', '') for b in os.listdir(BED_DIR) if b.endswith('.bed')])

# Snippy variant calling rule
rule run_snippy:
    input:
        ref = config["reference_genome"],
        contigs = "data/strains/{sample}.fasta"
    output:
        outdir = directory(f"{SNIPPY_DIR}/{{sample}}")
    params:
        ref_type = config["reference_type"]
    threads: config["snippy_threads"]
    conda:
        "envs/snippy_env.yaml"
    shell:
        "snippy --cpus {threads} --outdir {output.outdir} --{params.ref_type} {input.ref} "
        "--ctgs {input.contigs} --report"

# Rule to extract chromosome from a Snippy VCF file
rule extract_chromosome:
    input:
        vcf = f"{SNIPPY_DIR}/{config['representative_sample']}/snps.vcf"
    output:
        chrom = temp("chromosome.txt")
    conda:
        "envs/bcftools_env.yaml"
    shell:
        "bcftools view -H {input.vcf} | head -n 1 | cut -f1 > {output.chrom}"

# Rule for creating bed files
rule create_bed_files:
    input:
        gene_positions = "data/gene_positions.txt",
        chrom = "chromosome.txt"
    output:
        bed_dir = directory(BED_DIR)
    conda:
        "envs/python_env.yaml"
    shell:
        "python scripts/create_bed_files.py {input.gene_positions} {output.bed_dir} --chromosome $(cat {input.chrom})"

# Rule for intersecting VCF with BED files
rule intersect_vcf_bed:
    input:
        vcf = f"{SNIPPY_DIR}/{{sample}}/snps.vcf",
        bed = f"{BED_DIR}/{{bed}}.bed"
    output:
        intersected_vcf = f"{INTERSECT_DIR}/{{bed}}/{{sample}}_{{bed}}.vcf"
    conda:
        "envs/bedtools_env.yaml"
    shell:
        "intersectBed -a {input.vcf} -b {input.bed} -header > {output.intersected_vcf}"