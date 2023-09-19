#!/bin/env snakemake -s

import numpy as np
import pandas as pd
from os.path import join
from snakemake.utils import validate

# Set parameters.
shell.prefix("set -euo pipefail;")

# LOAD CONFIG FILES
configfile: 'config/config.yaml'

samples = pd.read_csv(
    config['samples'], 
    sep=';', 
    dtype=str,
    comment='#',
).set_index(['library'], drop=False)

# Set paths for both input and output files.

OUT_DIR =  config['out_dir']
TMP = join(config['base_dir'], config['tmp_dir'])
REF_DIR = join(config['base_dir'], config['ref_dir'])
FASTQ_DIR = join(config['base_dir'], config['fastq_dir'])


SAMPLES = samples

libraries = np.unique(samples.library)
species = np.unique(samples.species)
sampling_rates = np.unique(samples.sampling_rate)
enzymes = samples.enzymes


print(f"OUT_DIR: {OUT_DIR}")
print(f"TMP: {TMP}")
print(f"REF_DIR: {REF_DIR}")
print(f"FASTQ_DIR: {FASTQ_DIR}")
print(f"SAMPLES: {SAMPLES}")

wildcard_constraints:
    libraries = "|".join(libraries),
    species = "|".join(species),
    sampling_rates = "|".join(sampling_rates),
    enzymes = "|".join(enzymes)

print(f"enzymes: {enzymes}")
print(f"species: {species}")
print(f"libraries: {libraries}")


# SAMPLING_RATES = expand(sampling_rate)
# ENZYMES = expand(enzymes)

# print(f"SAMPLING_RATES: {SAMPLING_RATES}")
# print(f"ENZYMES: {ENZYMES}")


# print(f"species: {species}")
# print('sampling_rate: {sampling_rate}')

# test_cmd = expand("hicberg pipeline -r {sampling_rate} -e {enzyme}", sampling_rate=sampling_rates, enzyme=enzymes)
# print(f"test_cmd: {test_cmd}")

# # test = lambda w: join(REF_DIR, config[w.species]['ref'])
# # test = [join(REF_DIR, config[w.species]['ref']) for w in species]

# for w in species:
#     print(f"{join(REF_DIR, config[w.species]['ref'])}")

# print(f"test: {test}")



########################################################################
########################### TEST PART ##################################
########################################################################

# Rules to generates the build bowtie2 index and split fastq for alignement.

GENOME = "/home/sardine/Documents/genomes/Sc288_2m/SC288_with_micron.fa"
R1 = "/home/sardine/Documents/reads/AC1/AC1.end1_sub.fastq"
R2 = "/home/sardine/Documents/reads/AC1/AC1.end2_sub.fastq"

# OUTPUT = "/home/sardine/Bureau/"
BANK_NAME = species
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

# print(f"test : {expand(join(OUT_DIR, '{libraries}', 'index/'), libraries=libraries)}")

# Specify at least one output of each rule (or one of the last output)  to ensure rule execution.
rule all:
    input: 
        # expand(join(OUT_DIR, '{libraries}', 'index/'), libraries=libraries), 
        expand(join(OUT_DIR, '{libraries}', "fragments_fixed_sizes.txt"), libraries=libraries),
        expand(join(OUT_DIR, '{libraries}', "1.sorted.bam"), libraries=libraries),
        # "/home/sardine/Bureau/hic_test/1.sorted.bam",
        # "/home/sardine/Bureau/hic_test/group1.1.bam",
        expand(join(OUT_DIR, '{libraries}', "group1.1.bam"), libraries=libraries),
        # "/home/sardine/Bureau/hic_test/unrescued_map.cool",
        expand(join(OUT_DIR, '{libraries}', "unrescued_map.cool"), libraries=libraries),
        # "/home/sardine/Bureau/hic_test/restriction_map.npy",
        expand(join(OUT_DIR, '{libraries}', "restriction_map.npy"), libraries=libraries),
        expand(join(OUT_DIR, '{libraries}', "contacts", "matrices", "rescued_map.cool"), libraries = libraries)
        # "/home/sardine/Bureau/hic_test/contacts/matrices/rescued_map.cool",


