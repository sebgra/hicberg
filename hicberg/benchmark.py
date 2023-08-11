from os import path
from pathlib import Path
from multiprocessing import Process


from glob import glob
import tempfile as tmpf
import multiprocessing as mp
import subprocess as sp
from datetime import datetime


from itertools import product

from shutil import rmtree
import click
import logging

import numpy as np
import cooler

import hicberg.io as hio
import hicberg.utils as hut
# import hicberg.align as hal
import hicberg.statistics as hst
import hicberg.eval as hev
import hicberg.plot as hpl

BASE_MATRIX = "original_map.cool"
UNRESCUED_MATRIX = "unrescued_map.cool"
RESCUED_MATRIX = "rescued_map.cool"
RESTRICTION_MAP = "restriction_map.npy"
FORWARD_IN_FILE = "group1.1.in.bam"
REVERSE_IN_FILE = "group1.2.in.bam"
FORWARD_OUT_FILE = "group1.1.out.bam"
REVERSE_OUT_FILE = "group1.2.out.bam"

def benchmark(output_dir : str = None, chromosome : str = "", position : int = 0, trans_chromosome : str = None, trans_position : int = None, strides : list[int] = [], mode : str = "full", auto : int = None, bins : int = None, circular :str = "", genome : str = None):
    """
    AI is creating summary for benchmark

    Parameters
    ----------
    output_dir : str, optional
        [description], by default None
    chromosome : str, optional
        [description], by default ""
    position : int, optional
        [description], by default 0
    trans_chromosome : str, optional
        [description], by default None
    trans_position : int, optional
        [description], by default None
    strides : list[int], optional
        [description], by default []
    mode : str, optional
        [description], by default "full"
    auto : int, optional
        [description], by default None
    bins : int, optional
        [description], by default None
    """

    learning_status = False #False bckp
    picking_status = False # False bckp

    output_path = Path(output_dir)

    if not output_path.exists():
        raise ValueError(f"Output directory {output_dir} does not exist.")
    
    genome_path = output_path / genome

    if not genome_path.exists():
        raise ValueError(f"Genome file {genome} does not exist. PLease provide an existing genome file.")
    
    restriction_map_path = output_path / RESTRICTION_MAP

    if not restriction_map_path.exists():
        raise ValueError(f"Restriction map file {restriction_map_path} does not exist. PLease provide an existing restriction map file.")
    
    # Define file to store results

    header = f"date\tchrom\tpos\tstride\ttrans_chrom\ttrans_pos\tmode\tnb_reads\tscore\n"

    results = output_path / "benchmark.csv"

    base_matrix_path = output_path / BASE_MATRIX
    base_matrix = hio.load_cooler(base_matrix_path)

    restriction_map = hio.load_dictionary(restriction_map_path)

    OTHER = [[chromosome, position, trans_chromosome, trans_position, strides, auto, bins]]


    CHROMOSOME = list(chromosome)

    if trans_chromosome is not None:
        TRANS_CHROMOSOME = [str(t) for t in trans_chromosome.split(",")]
    else : 
        TRANS_CHROMOSOME = None

    if trans_position is not None:
        TRANS_POSITION = [int(p) for p in trans_position.split(",")]

    else :
        TRANS_POSITION = None


    NB_BINS = bins
    # POSITION = [int(p) for p in str(position).split(",")]
    POSITION = position
    BIN_SIZE = bins

    STRIDES = [int(s) for s in strides.split(",")]

    for mode in mode.split(","):

        # Pick reads

        if not picking_status:

            intervals_dictionary = hev.select_reads(matrix_file = output_path / BASE_MATRIX, position = POSITION, chromosome = CHROMOSOME, strides = STRIDES, trans_chromosome = TRANS_CHROMOSOME, trans_position = TRANS_POSITION, auto = auto, nb_bins = NB_BINS, output_dir = output_dir)
            indexes = hev.get_bin_indexes(matrix = base_matrix, dictionary = intervals_dictionary, )

            picking_status = True

        forward_in_path = output_path / FORWARD_IN_FILE
        reverse_in_path = output_path / REVERSE_IN_FILE
        forward_out_path = output_path / FORWARD_OUT_FILE
        reverse_out_path = output_path / REVERSE_OUT_FILE

        # Get corresponding indexes to the duplicated reads coordinates.

        print(f"indexes : {indexes}")

        # Re-build pairs and cooler matrix 

        hio.build_pairs(bam_for = forward_out_path, bam_rev = reverse_out_path, output_dir = output_path)
        hio.build_matrix(output_dir = output_path)
        
        # Reattribute reads from inner group

        if not learning_status : 
            ## Compute statistics

            p1 = Process(target = hst.get_patterns(forward_bam_file  = forward_out_path, reverse_bam_file = reverse_out_path, circular = circular, output_dir = output_path))
            p2 = Process(target = hst.generate_trans_ps(restriction_map = restriction_map, output_dir = output_path))
            p3 = Process(target = hst.generate_coverages(forward_bam_file = forward_out_path, reverse_bam_file  = reverse_out_path, genome = genome, bins = bins, output_dir = output_path))
            p4 = Process(target = hst.generate_d1d2(forward_bam_file = forward_out_path, reverse_bam_file  = reverse_out_path, output_dir = output_path))

            # Launch processes
            for process in [p1, p2, p3, p4]:
                process.start()
                process.join()

            learning_status = True

        # Reattribute reads
        hst.reattribute_reads(reads_couple = (forward_in_path, reverse_in_path), mode = mode, output_dir = output_path)
        hio.merge_predictions(output_dir = output_path, clean = True)
        hio.build_pairs(bam_for  = "group1.1.out.bam", bam_rev = "group1.2.out.bam", bam_for_rescued  = "group2.1.rescued.bam", bam_rev_rescued = "group2.2.rescued.bam", mode = True, output_dir = output_path)
        hio.build_matrix(mode = True, output_dir = output_path)

        rescued_matrix = hio.load_cooler(output_path / RESCUED_MATRIX)
        rescued_matrix_array = rescued_matrix.matrix(balance = False)

        rescued_matrix_path = output_path / RESCUED_MATRIX

        pearson = hst.pearson_score(original_matrix = base_matrix, rescued_matrix = rescued_matrix , markers = indexes)

        print(f"CHROMOSOME : {CHROMOSOME}")
        print(f"dtype : {type(CHROMOSOME)}")
        print(f"TRANS_CHROMOSOME : {TRANS_CHROMOSOME}")
        print(f"dtype : {type(TRANS_CHROMOSOME)}")


        chromosome_set = CHROMOSOME + TRANS_CHROMOSOME if TRANS_CHROMOSOME is not None else CHROMOSOME

        print(f"Chromosome set : {chromosome_set}")

        hpl.plot_benchmark(original_matrix = BASE_MATRIX, depleted_matrix = UNRESCUED_MATRIX, rescued_matrix = RESCUED_MATRIX, chromosomes = chromosome_set, output_dir = output_path)

        #TODO implement tidying plot function

        # number_reads = np.sum(rescued_matrix_array[indexes])
        number_reads = 10

        print(f"Pearson score : {pearson:9.4f} in mode {mode}")

        if not results.exists():
            with open(results, "w") as f_out:
                f_out.write(header)

                date = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
                f_out.write(f"{date}\t{CHROMOSOME}\t{POSITION}\t{STRIDES}\t{TRANS_CHROMOSOME}\t{TRANS_POSITION}\t{mode}\t{number_reads}\t{pearson:9.4f}\n")
                f_out.close()

        else :
            with open(results, "a") as f_out:
                date = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
                f_out.write(f"{date}\t{CHROMOSOME}\t{POSITION}\t{STRIDES}\t{TRANS_CHROMOSOME}\t{TRANS_POSITION}\t{mode}\t{number_reads}\t{pearson:9.4f}\n")
                f_out.close()
                
                
    return

