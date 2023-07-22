import time
import glob, sys
from glob import glob
from shutil import which, rmtree
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

        output_folder = hio.create_folder(sample_name = name, output_dir = output_dir)

        hut.get_chromosomes_sizes(genome = genome, output_dir = output_folder)
        hut.get_bin_table(bins = bins, output_dir = output_folder)

        index = hal.hic_build_index(genome = genome, output_dir = output_folder, cpus = cpus, verbose = verbose)

        hal.hic_align(genome = genome, index = index, fq_for = fq_for, fq_rev = fq_rev, output_dir = output_folder, cpus = cpus, verbose = verbose)

    if exit_stage == 1:

        logger.info(f"Ending HiCBERG pipeline at {exit_stage}")
        return
    

    if start_stage < 2: 

        hal.hic_view(cpus = cpus, output_dir = output_folder, verbose = True)
        hal.hic_sort(cpus = cpus, output_dir = output_folder, verbose = True)
        hut.classify_reads(output_dir = output_folder)

    if exit_stage == 2:

        logger.info(f"Ending HiCBERG pipeline at {exit_stage}")
        return
    
    if start_stage < 3:

        restriction_map = hst.get_restriction_map(genome = genome, enzyme = enzyme)
        hst.get_dist_frags(genome = genome, restriction_map = restriction_map, circular = circular, rate = rate, output_dir = output_folder)
        hio.build_pairs(output_dir = output_folder)
        hio.build_matrix(output_dir = output_folder)

        hst.log_bin_genome(genome = genome, output_dir = output_folder)
        hst.get_patterns(circular = circular, output_dir = output_folder)
        hst.generate_trans_ps(restriction_map = restriction_map, output_dir = output_folder)
        hst.generate_coverages(genome = genome, bins = bins, output_dir = output_folder)
        hst.generate_d1d2(output_dir = output_folder)
        

    if exit_stage == 3:
        logger.info(f"Ending HiCBERG pipeline at {exit_stage}")
        return
    
    if start_stage < 4:

        hut.chunk_bam(nb_chunks = nb_chunks, output_dir = output_folder)
        
        # Get chunks as lists
        forward_chunks, reverse_chunks = hut.get_chunks(output_folder)
   
        # Reattribute reads
        with multiprocessing.Pool(processes = cpus) as pool:

            results = pool.map_async(partial(hst.reattribute_reads, restriction_map = restriction_map, output_dir = output_folder),
            zip(forward_chunks, reverse_chunks))
            pool.close()
            pool.join()

        hio.merge_predictions(output_dir = output_folder)

        # Delete chunks
        folder_to_delete = Path(output_folder) / 'chunks'
        rmtree(folder_to_delete)


        print(f"output_folder : {output_folder}")


        hio.build_pairs(mode = True, output_dir = output_folder)
        hio.build_matrix(mode = True, output_dir = output_folder)

    if start_stage <= 5:

        
        hpl.plot_laws(output_dir = output_folder)
        hpl.plot_trans_ps(output_dir = output_folder)
        hpl.plot_coverages(bins = bins, output_dir = output_folder)


    
        
        
        # couple = ("/home/sardine/Bureau/sample/group2.1.bam", "/home/sardine/Bureau/sample/group2.2.bam")

        # hst.reattribute_reads(reads_couple = couple, restriction_map = restriction_map, output_dir = output_folder,)

    
    logger.info("Ending HiCBERG pipeline")


if __name__ == "__main__":

    pipeline(genome = "/home/sardine/Bureau/codes/hicberg/data_test/SC288_with_micron.fa", fq_for = "/home/sardine/Bureau/codes/hicberg/data_test/forward_reads_test.fq.gz", fq_rev = "/home/sardine/Bureau/codes/hicberg/data_test/reverse_reads_test.fq.gz", output_dir = "/home/sardine/Bureau/", cpus = 8, rate = 0.1, nb_chunks = 8)

