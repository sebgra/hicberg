import uuid
from os import mkdir #path
from shutil import rmtree
from pathlib import Path
from multiprocessing import Process


from glob import glob
# import tempfile as tmpf
# import multiprocessing as mp
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

from hicberg import benchmark_logger

BAM_FOR = "group1.1_subsampled.bam" 
BAM_REV = 'group1.2_subsampled.bam'
BASE_MATRIX = "original_map.cool"
UNRESCUED_MATRIX = "unrescued_map.cool"
RESCUED_MATRIX = "rescued_map.cool"
RESTRICTION_MAP = "restriction_map.npy"
FORWARD_IN_FILE = "group1.1.in.bam"
REVERSE_IN_FILE = "group1.2.in.bam"
FORWARD_OUT_FILE = "group1.1.out.bam"
REVERSE_OUT_FILE = "group1.2.out.bam"

# TODO : Complete docstring
def benchmark(output_dir : str = None, chromosome : str = "", position : int = 0, trans_chromosome : str = None, trans_position : int = None, strides : list[int] = [], mode : str = "full", auto : int = None, rounds : int = 1, magnitude : float = 1.0, bins : int = None, circular :str = "", genome : str = None, pattern : str = None, threshold : float = 0.0, jitter : int = 0, trend : bool = True, top : int = 100,  force : bool = False):
    """AI is creating summary for benchmark

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
    rounds : int, optional
        [description], by default 1
    magnitude : float, optional
        [description], by default 1.0
    bins : int, optional
        [description], by default None
    circular : str, optional
        [description], by default ""
    genome : str, optional
        [description], by default None
    force : bool, optional
        [description], by default False

    Raises
    ------
    ValueError
        [description]
    ValueError
        [description]
    ValueError
        [description]
    """    
    # logger.addHandler('hicberg_benchmark.log')

    args = locals()

    # Keep track of the arguments used
    for arg in args:

        benchmark_logger.info("%s: %s", arg, args[arg])

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
    header = f"id\tdate\tchrom\tpos\tstride\ttrans_chrom\ttrans_pos\tmode\tnb_reads\tscore\n"

    results = output_path / "benchmark.csv"

    base_matrix_path = output_path / BASE_MATRIX
    base_matrix = hio.load_cooler(base_matrix_path)
    bin_size = base_matrix.binsize

    restriction_map = hio.load_dictionary(restriction_map_path)

    flag_file = output_path / "flag.txt"
    # TODO : insert pattern based benchmarking mode

    ## Chromosomsight pre-call
    if pattern is not None:

        pre_recall_cmd = hev.chromosight_cmd_generator(file = base_matrix_path, pattern = pattern, untrend = trend, output_dir = output_path)

        print(f"Pre-call command : {pre_recall_cmd}")
        benchmark_logger.info("Starting Chromosight pre-call")
        sp.run(pre_recall_cmd, shell = True)

    # reformat inputs
    if pattern is None:
    

        chromosome = [str(c) for c in chromosome.split(",")]

        if trans_chromosome is not None and trans_chromosome != "-1":
            trans_chromosome = [str(t) for t in trans_chromosome.split(",")]

        else : 
            trans_chromosome = None

        if trans_position is not None and trans_position != "-1": # -1 usefull when not using trans cases with benchmark calling through Snakemake
            trans_position = [int(p) for p in trans_position.split(",")]

        else :
            trans_position = None



        nb_bins = bins
        # POSITION = [int(p) for p in str(position).split(",")]
        # POSITION = position


        strides = [int(s) for s in strides.split(",")]

    elif pattern is not None : # Pattern based benchmarking

        # TODO : add pattern based benchmarking args reformating
        chromosome = chromosome
        df = hev.get_top_pattern(file = output_path / "original.tsv", top = top, threshold = threshold, chromosome = chromosome).sort_values(by='start1', ascending=True)


        strides = [int(df.iloc[i].start1 - df['start1'].min()) for i in range(0, df.shape[0])]
        nb_bins = bins

        if trans_chromosome is not None and trans_chromosome != "-1":
            trans_chromosome = [str(t) for t in trans_chromosome.split(",")]

        if trans_position is not None and trans_position != "-1": # -1 usefull when not using trans cases with benchmark calling through Snakemake
            trans_position = [int(p) for p in trans_position.split(",")]
        

    for sub_mode in mode.split(","):

        # Pick reads

        if not picking_status:
            
            benchmark_logger.info("Picking reads")
            intervals_dictionary = hev.select_reads(bam_for = BAM_FOR, bam_rev = BAM_REV, matrix_file = output_path / BASE_MATRIX, position = position, chromosome = chromosome, strides = strides, trans_chromosome = trans_chromosome, trans_position = trans_position, auto = auto, nb_bins = nb_bins, output_dir = output_dir)
            
            benchmark_logger.info(f"intervals_dictionary : {intervals_dictionary}")
            indexes = hev.get_bin_indexes(matrix = base_matrix, dictionary = intervals_dictionary, )

            picking_status = True

        forward_in_path = output_path / FORWARD_IN_FILE
        reverse_in_path = output_path / REVERSE_IN_FILE
        forward_out_path = output_path / FORWARD_OUT_FILE
        reverse_out_path = output_path / REVERSE_OUT_FILE

        # Get corresponding indexes to the duplicated reads coordinates.

        # Re-build pairs and cooler matrix 

        hio.build_pairs(bam_for = forward_out_path, bam_rev = reverse_out_path, output_dir = output_path)
        hio.build_matrix(output_dir = output_path)

        unrescued_map_path = output_path / UNRESCUED_MATRIX
        
        # Reattribute reads from inner group

        if not learning_status : 
            ## Compute statistics

            p1 = Process(target = hst.get_patterns, kwargs = dict(forward_bam_file  = forward_out_path, reverse_bam_file = reverse_out_path, circular = circular, output_dir = output_path))
            p2 = Process(target = hst.generate_trans_ps, kwargs = dict(restriction_map = restriction_map, output_dir = output_path))
            p3 = Process(target = hst.generate_coverages, kwargs = dict(forward_bam_file = forward_out_path, reverse_bam_file  = reverse_out_path, genome = genome, bins = bin_size, output_dir = output_path))
            p4 = Process(target = hst.generate_d1d2, kwargs = dict(forward_bam_file = forward_out_path, reverse_bam_file  = reverse_out_path, output_dir = output_path))
            p5 = Process(target = hst.generate_density_map, kwargs = dict(matrix = unrescued_map_path, rounds = rounds, magnitude = magnitude, output_dir = output_path))
            
            # if  "full" in mode.split(","):

            benchmark_logger.info("Full mode selected. Learning step will be performed.")

            # Launch processes
            for process in [p1, p2, p3, p4, p5]:
                process.start()

            for process in [p1, p2, p3, p4, p5]:
                process.join()

            benchmark_logger.info("Learning step completed")

            # elif "ps_only" in mode.split(",") and len(mode.split(",")) == 1:

            #     benchmark_logger.info("PS only mode selected. Learning step will be performed.")

            #     # Launch processes
            #     for process in [p1, p2]:
            #         process.start()

            #     for process in [p1, p2]:
            #         process.join()

            #     benchmark_logger.info("Learning step completed")


            # elif "cover_only" in mode.split(",") and len(mode.split(",")) == 1:

            #     # Launch processes
            #     for process in [p3]:
            #         process.start()

            #     for process in [p3]:
            #         process.join()

            # elif "d1d2_only" in mode.split(",") and len(mode.split(",")) == 1:

            #     benchmark_logger.info("D1D2 only mode selected. Learning step will be performed.")

            #     # Launch processes
            #     for process in [p4]:
            #         process.start()

            #     for process in [p4]:
            #         process.join()

            #     benchmark_logger.info("Learning step completed")


            # elif "density_only" in mode.split(",") and len(mode.split(",")) == 1:

            #     benchmark_logger.info("Density only mode selected. Learning step will be performed.")

            #     # Launch processes
            #     for process in [p5]:
            #         process.start()

            #     for process in [p5]:
            #         process.join()

            #     benchmark_logger.info("Learning step completed")


            # elif "no_ps" in mode.split(",") and len(mode.split(",")) == 1:

            #     benchmark_logger.info("No PS mode selected. Learning step will be performed.")

            #     # Launch processes
            #     for process in [p3, p4, p5]:
            #         process.start()

            #     for process in [p3, p4, p5]:
            #         process.join()

            #     benchmark_logger.info("Learning step completed")


            # elif "no_cover" in mode.split(",") and len(mode.split(",")) == 1:

            #     benchmark_logger.info("No cover mode selected. Learning step will be performed.")

            #     # Launch processes
            #     for process in [p1, p2, p4, p5]:
            #         process.start()

            #     for process in [p1, p2, p4, p5]:
            #         process.join()

            #     benchmark_logger.info("Learning step completed")


            # elif "no_d1d2" in mode.split(",") and len(mode.split(",")) == 1:

            #     benchmark_logger.info("No D1D2 mode selected. Learning step will be performed.")
                

            #     # Launch processes
            #     for process in [p1, p2, p3, p5]:
            #         process.start()

            #     for process in [p1, p2, p3, p5]:
            #         process.join()

            #     benchmark_logger.info("Learning step completed")

            # elif "no_density" in mode.split(",") and len(mode.split(",")) == 1:

            #     benchmark_logger.info("No density mode selected. Learning step will be performed.")

            #     # Launch processes
            #     for process in [p1, p2, p3, p4]:
            #         process.start()

            #     for process in [p1, p2, p3, p4]:
            #         process.join()

            #     benchmark_logger.info("Learning step completed")

            # elif "d1d2_only" not in mode.split(",") and len(mode.split(",")) > 1:

            #     benchmark_logger.info("D1D2 only mode selected. Learning step will be performed.")

            #     # Launch processes
            #     for process in [p1, p2, p3]:
            #         process.start()

            #     for process in [p1, p2, p3]:
            #         process.join()

            #     benchmark_logger.info("Learning step completed")

            # elif  "random" in mode.split(",") and len(mode.split(",")) == 1:

            #     benchmark_logger.info("Random mode selected. No learning step will be performed.")

        learning_status = True

        # Reattribute reads

        benchmark_logger.info("Re-attributing reads")
        hst.reattribute_reads(reads_couple = (forward_in_path, reverse_in_path), mode = sub_mode, output_dir = output_path)
        hio.merge_predictions(output_dir = output_path, clean = True)
        hio.build_pairs(bam_for  = "group1.1.out.bam", bam_rev = "group1.2.out.bam", bam_for_rescued  = "group2.1.rescued.bam", bam_rev_rescued = "group2.2.rescued.bam", mode = True, output_dir = output_path)
        hio.build_matrix(mode = True, output_dir = output_path)

        rescued_matrix = hio.load_cooler(output_path / RESCUED_MATRIX)
        # rescued_matrix_array = rescued_matrix.matrix(balance = False)

        rescued_matrix_path = output_path / RESCUED_MATRIX

        pearson = hst.pearson_score(original_matrix = base_matrix, rescued_matrix = rescued_matrix , markers = indexes)

        if pattern is None:
            chromosome_set = [*chromosome, *trans_chromosome] if trans_chromosome is not None else chromosome
        
        else : 
            chromosome_set = [chromosome, *trans_chromosome] if trans_chromosome is not None else chromosome

        hpl.plot_benchmark(original_matrix = BASE_MATRIX, depleted_matrix = UNRESCUED_MATRIX, rescued_matrix = RESCUED_MATRIX, chromosomes = chromosome_set, output_dir = output_path)

        # Define unique id to keep track of the experiments
        id_tag = str(uuid.uuid4())[:8]

        number_reads = 10

        benchmark_logger.info(f"Pearson score : {pearson:9.4f} in mode {sub_mode}")

        if not results.exists():
            with open(results, "w") as f_out:
                f_out.write(header)

                date = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
                f_out.write(f"{id_tag}\t{date}\t{chromosome}\t{position}\t{strides}\t{trans_chromosome}\t{trans_position}\t{sub_mode}\t{number_reads}\t{pearson:9.4f}\n")
                f_out.close()

        else :
            with open(results, "a") as f_out:
                date = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
                f_out.write(f"{id_tag}\t{date}\t{chromosome}\t{position}\t{strides}\t{trans_chromosome}\t{trans_position}\t{sub_mode}\t{number_reads}\t{pearson:9.4f}\n")
                f_out.close()

        # tidy plots
        folder_path = Path(output_dir, id_tag)

        if folder_path.exists() and force: 

            rmtree(folder_path)

        if not folder_path.exists() : 
            
            mkdir(folder_path)
                
        files = [p for  p in output_path.glob("*")]

        for file in files :

            if Path(file).suffix == ".pdf" or Path(file).suffix == ".svg":

                Path(file).rename(folder_path / Path(file).name)

        benchmark_logger.info(f"Ending benchmark. Results stored in {folder_path}")

    forward_in_path.unlink()
    reverse_in_path.unlink()
    forward_out_path.unlink()
    reverse_out_path.unlink()

    open(flag_file, 'a').close()

    return

