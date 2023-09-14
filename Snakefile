#!/bin/env snakemake -s

# Rules to generates the build bowtie2 index and split fastq for alignement.

GENOME = "/home/sardine/Documents/genomes/Sc288_2m/SC288_with_micron.fa"
R1 = "/home/sardine/Documents/reads/AC1/AC1.end1_sub.fastq"
R2 = "/home/sardine/Documents/reads/AC1/AC1.end2_sub.fastq"

OUTPUT = "/home/sardine/Bureau/"
BANK_NAME = "hic_test"
THREADS = 16
RATE = 0.2
MODE = "full"
# START = "fastq"
# EXIT = "bam"
ROUNDS = 10
MAGNITUDE = 0.8
SENSITIVITY = "very-fast"
ENZYMES_0 = "DpnII"
ENZYMES_1  = "HinfI"

# Set parameters.
shell.prefix("set -euo pipefail;")

# Pipeline sub-workflows
include: 'rules/00_hicberg_step_0.smk'
include: 'rules/01_hicberg_step_1.smk'
include: 'rules/02_hicberg_step_2.smk'
include: 'rules/03_hicberg_step_3.smk'
include: 'rules/04_hicberg_step_4.smk'
include: 'rules/05_hicberg_step_5.smk'

# Specify at least one output of each rule (or one of the last output)  to ensure rule execution.
rule all:
    input: 
        "/home/sardine/Bureau/hic_test/index/",
        # "/home/sardine/Bureau/hic_test/1.sorted.bam",
        "/home/sardine/Bureau/hic_test/group1.1.bam",
        "/home/sardine/Bureau/hic_test/unrescued_map.cool",
        "/home/sardine/Bureau/hic_test/restriction_map.npy",
        "/home/sardine/Bureau/hic_test/contacts/matrices/rescued_map.cool",


