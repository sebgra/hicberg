import time
import glob, sys
from shutil import which
from os.path import join
from pathlib import Path
import subprocess as sp

import multiprocessing
from functools import partial
import logging

import numpy as np

import hicberg.align as hal
import hicberg.io as hio
import hicberg.utils as hut
import hicberg.plot as hpl
import hicberg.statistics as hst

from hicberg import logger


def check_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    return which(name) is not None

def pipeline(name :str = "sample",start_stage : str = "fastq", exit_stage : str = "None", genome : str = None,
            fq_for : str = None, fq_rev : str = None, sensitivity : str = "very-sensitive",
            max_alignment : int = None, mapq : int = 35, enzyme  : list[str] = ["DpnII", "HinfI"],
            circular : str = "", rate : float = 1.0, bins : int = 2000, nb_chunks : int = 1,
            mode : str = "full",  verbose : bool = False, cpus : int = 1, output_dir : str = None) -> None :


    logger.info('This is an info message from the main function')
    
    args = locals()

    if not check_tool("bowtie2"):
        logger.error("bowtie2 is not available on your system.")
        raise ValueError("bowtie2 is not available on your system.")

    if not check_tool("samtools"):
        logger.error("samtools is not available on your system.")
        raise ValueError("samtools is not available on your system.")
    
    stages = {"fastq": 0, "bam": 1, "stats": 2, "pairs": 3, "rescue": 4, "cool": 5}

    out_stage = {"None": None, "bam": 1, "stats": 2, "pairs": 3, "rescue": 4, "cool": 5}

    start_stage = stages[
        start_stage
    ]  # start_stage as variable of command line - default to "fastq" --> 0

    exit_stage = out_stage[exit_stage]

    logger.info("Start HiCBERG pipeline")

    # Keep track of the arguments used
    for arg in args:

        logger.info("%s: %s", arg, args[arg])


    if start_stage < 1 : 

        hio.create_folder(sample_name = name, output_dir = output_dir)

        hut.get_chromosomes_sizes(genome = genome, output_dir = output_dir)

        index = hal.hic_build_index(genome = genome, output_dir = output_dir, verbose = verbose)

        hal.hic_align(genome = genome, index = index, fq_for = fq_for, fq_rev = fq_rev, output_dir = output_dir, cpus = cpus, verbose = verbose)
    
    logger.info("Ending HiCBERG pipeline")


if __name__ == "__main__":

    pipeline(genome = "/home/sardine/Bureau/codes/hicberg/data_test/SC288_with_micron.fa", fq_for = "/home/sardine/Bureau/codes/hicberg/data_test/forward_reads_test.fq.gz", fq_rev = "/home/sardine/Bureau/codes/hicberg/data_test/reverse_reads_test.fq.gz", output_dir = "/home/sardine/Bureau/", cpus = 8)

