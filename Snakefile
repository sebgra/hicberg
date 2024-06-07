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
REF_DIR = join(config['base_dir'], config['ref_dir'])
FASTQ_DIR = join(config['base_dir'], config['fastq_dir'])


libraries = np.unique(samples.library)
species = np.unique(samples.species)
sampling_rates = samples.sampling_rates
enzymes = samples.enzymes
modes = samples.modes
resolutions = samples.resolutions
kernel_sizes = samples.kernel_sizes
deviations = samples.deviations
max_reports = samples.max_reports
circularity = samples.circularity
distances = samples.distances
blacklists = samples.blacklists


wildcard_constraints:
    libraries = "|".join(libraries),
    species = "|".join(species),
    sampling_rates = "|".join(sampling_rates),
    enzymes = "|".join(enzymes),
    modes = "|".join(modes),
    resolutions = "|".join(resolutions),
    kernel_sizes = "|".join(kernel_sizes),
    deviations = "|".join(deviations),
    max_reports = "|".join(max_reports),
    circularity = "|".join(circularity),
    distances = "|".join(distances),
    blacklists = "|".join(blacklists)


########################################################################
########################### RULES PART #################################
########################################################################


# OUTPUT = "/home/sardine/Bureau/"
THREADS = 16
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
        expand(join(OUT_DIR, '{libraries}', "fragments_fixed_sizes.txt"), libraries=libraries),
        expand(join(OUT_DIR, '{libraries}', "1.sorted.bam"), libraries=libraries),
        expand(join(OUT_DIR, '{libraries}', "group1.1.bam"), libraries=libraries),
        expand(join(OUT_DIR, '{libraries}', "unrescued_map.cool"), libraries=libraries),
        expand(join(OUT_DIR, '{libraries}', "restriction_map.npy"), libraries=libraries),
        expand(join(OUT_DIR, '{libraries}', "contacts", "matrices", "rescued_map.cool"), libraries = libraries), 


