import uuid
from os import mkdir #path
from shutil import rmtree, copy
from pathlib import Path
from functools import partial

from multiprocessing import Process, Pool


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

from hicberg import logger

# BAM_FOR = "group1.1_subsampled.bam" 
# BAM_REV = 'group1.2_subsampled.bam'
BAM_FOR = "group1.1.bam" 
BAM_REV = 'group1.2.bam'
BASE_MATRIX = "original_map.cool"
FRAGMENTS = "fragments_fixed_sizes.txt"
DIST_FRAGS = "dist.frag.npy"
CHROMOSOME_SIZES = "chromosome_sizes.npy"
XS = "xs.npy"
UNRESCUED_MATRIX = "unrescued_map.cool"
RESCUED_MATRIX = "rescued_map.cool"
RESTRICTION_MAP = "restriction_map.npy"
FORWARD_IN_FILE = "group1.1.in.bam"
REVERSE_IN_FILE = "group1.2.in.bam"
FORWARD_OUT_FILE = "group1.1.out.bam"
REVERSE_OUT_FILE = "group1.2.out.bam"

# TODO : Complete docstring
def benchmark(output_dir : str = None, chromosome : str = "", position : int = 0, trans_chromosome : str = None, trans_position : int = None, strides : list[int] = [], mode : str = "full", auto : int = None, kernel_size : int = 11, deviation : float = 0.5, bins : int = None, circular :str = "", genome : str = None, pattern : str = None, threshold : float = 0.0, jitter : int = 0, trend : bool = True, top : int = 100,  force : bool = False, iterations : int = 3, cpus : int = 8):
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
    kernel_size : int, optional
        [description], by default 11
    deviation : float, optional
        [description], by default 0.5
    bins : int, optional
        [description], by default None
    circular : str, optional
        [description], by default ""
    genome : str, optional
        [description], by default None
    pattern : str, optional
        [description], by default None
    threshold : float, optional
        [description], by default 0.0
    jitter : int, optional
        [description], by default 0
    trend : bool, optional
        [description], by default True
    top : int, optional
        [description], by default 100
    force : bool, optional
        [description], by default False
    iterations : int, optional
        [description], by default 3
    cpus : int, optional
        [description], by default 8

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
        logger.info("%s: %s", arg, args[arg])

    learning_status = False #False bckp
    picking_status = False # False bckp

    # Setting files paths
    output_path = Path(output_dir)
    output_path_chunks = Path(output_dir, "chunks")
    restriction_map_path = output_path / RESTRICTION_MAP
    bam_for_path = output_path / BAM_FOR
    bam_rev_path = output_path / BAM_REV

    fragments_path = output_path / FRAGMENTS 
    chromosome_size_path = output_path / CHROMOSOME_SIZES
    dist_frags_path = output_path / DIST_FRAGS
    xs_path = output_path / XS


    # Define unique id to keep track of the experiments
    id_tag = str(uuid.uuid4())[:8]

    # Define output child directory
    output_uniq_path = Path(output_dir, id_tag)
    output_data_path = Path(output_dir, id_tag, "data")
    output_plot_path = Path(output_dir, id_tag, "plots")

    if not output_uniq_path.exists():
        mkdir(output_uniq_path)

    if not output_data_path.exists():
        mkdir(output_data_path)
    if not output_plot_path.exists():
        mkdir(output_plot_path)


    if not output_path.exists():
        raise ValueError(f"Output directory {output_dir} does not exist.")

    if not restriction_map_path.exists():
        raise ValueError(f"Restriction map file {restriction_map_path} does not exist. PLease provide an existing restriction map file.")
    
    # Define file to store results
    header = f"id\tdate\tchrom\tpos\tstride\ttrans_chrom\ttrans_pos\tauto\tbins\tmode\tnb_reads\tpattern\tprecision\trecall\tf1_score\tscore\n"
    results = output_path / "benchmark.csv"

    # Copy files to data_path
    copy(output_path /BASE_MATRIX, output_data_path / BASE_MATRIX)
    copy(bam_for_path, output_data_path / BAM_FOR)
    copy(bam_rev_path, output_data_path / BAM_REV)
    copy(restriction_map_path, output_data_path / RESTRICTION_MAP)
    copy(fragments_path, output_data_path / FRAGMENTS)
    copy(chromosome_size_path, output_data_path / CHROMOSOME_SIZES)
    copy(dist_frags_path, output_data_path / DIST_FRAGS)
    copy(xs_path, output_data_path / XS)

    base_matrix = hio.load_cooler(output_data_path / BASE_MATRIX)
    bin_size = base_matrix.binsize
    flag_file = output_path / "flag.txt"



    forward_in_path = output_data_path / FORWARD_IN_FILE
    reverse_in_path = output_data_path / REVERSE_IN_FILE
    forward_out_path = output_data_path / FORWARD_OUT_FILE
    reverse_out_path = output_data_path / REVERSE_OUT_FILE

    ## Chromosomsight pre-call       
    # reformat inputs
    if pattern is None or pattern == "-1":
    
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

        strides = [int(s) for s in strides.split(",")]

    elif pattern is not None and pattern != "-1": # Pattern based benchmarking
        pre_recall_cmd = hev.chromosight_cmd_generator(file = output_data_path / BASE_MATRIX, pattern = pattern, untrend = trend, output_dir = output_data_path)

        logger.info("Starting Chromosight pre-call")
        sp.run(pre_recall_cmd, shell = True)

        chromosome = chromosome

        if trans_chromosome is not None and trans_chromosome != "-1":
            trans_chromosome = [str(t) for t in trans_chromosome.split(",")]

        else : 
            trans_chromosome = None

        if trans_position is not None and trans_position != "-1": # -1 usefull when not using trans cases with benchmark calling through Snakemake
            trans_position = [int(p) for p in trans_position.split(",")]

        else :
            trans_position = None

        df = hev.get_top_pattern(file = output_data_path / "original.tsv", top = top, threshold = threshold, chromosome = chromosome).sort_values(by='start1', ascending=True)
        print(df)

        position = df.iloc[0].start1 # select top score pattern

        if strides is None or strides == "-1":

            strides = [int(df.iloc[i].start1 - df['start1'].min()) for i in range(1, df.shape[0])]
        
        else : 

            strides = [int(s) for s in strides.split(",")]
        
        nb_bins = bins

    for sub_mode in mode.split(","):

        precision = None
        recall = None
        f1_score = None

        # Pick reads
        if not picking_status:
            
            # # chunk bam files
            forward_chunks, reverse_chunks = hut.get_chunks(output_dir = output_path.as_posix())

            # Multithread part 

            coordinates = hev.generate_dict_coordinates(matrix_file = output_path / BASE_MATRIX, position = position, chromosome = chromosome, strides = strides, trans_chromosome = trans_chromosome, trans_position = trans_position, auto = auto, nb_bins = nb_bins, output_dir = output_data_path)

            # TODO : pass coordinates to select_reads_multithreads to avoid seeding problems through multithreading

            # Pick reads reads
            with Pool(processes = cpus) as pool: # cpus

                res = pool.map(partial(hev.select_reads_multithreads, interval_dictionary =  coordinates, output_dir = output_data_path),
                zip(forward_chunks, reverse_chunks))
                pool.close()
                pool.join()

            # Get dictionrary of intervals from Pool
            intervals_dictionary = coordinates

            hio.merge_predictions(output_dir  = output_data_path, clean = True, stage = "benchmark", cpus = cpus)

            # TODO : put aside in function
            indexes = hev.get_bin_indexes(matrix = base_matrix, dictionary = intervals_dictionary, )
            picking_status = True


            # Get corresponding indexes to the duplicated reads coordinates.
            # Re-build pairs and cooler matrix 
            hio.build_pairs(bam_for = forward_out_path, bam_rev = reverse_out_path, output_dir = output_data_path)
            hio.build_matrix(balance = True, output_dir = output_data_path)

        unrescued_map_path = output_path / UNRESCUED_MATRIX
        
        # Reattribute reads from inner group
        if not learning_status : 
            ## Compute statistics
            p1 = Process(target = hst.get_patterns, kwargs = dict(forward_bam_file  = forward_out_path, reverse_bam_file = reverse_out_path, circular = circular, output_dir = output_data_path))
            p2 = Process(target = hst.generate_trans_ps, kwargs = dict(output_dir = output_data_path))
            p3 = Process(target = hst.generate_coverages, kwargs = dict(forward_bam_file = forward_out_path, reverse_bam_file  = reverse_out_path, genome = genome, bins = bin_size, output_dir = output_data_path))
            p4 = Process(target = hst.generate_d1d2, kwargs = dict(forward_bam_file = forward_out_path, reverse_bam_file  = reverse_out_path, output_dir = output_data_path))
            
            logger.info("Full mode selected. Learning step will be performed.")

            # Launch processes
            for process in [p1, p2, p3, p4]:
                process.start()

            # for process in [p1, p2, p3, p4, p5]:
            for process in [p1, p2, p3, p4]:
                process.join()

            logger.info("Learning step completed")

            # TODO : restore
            hst.compute_density(cooler_file = UNRESCUED_MATRIX, kernel_size = kernel_size, deviation = deviation, threads = cpus, output_dir  = output_data_path)

        learning_status = True


        for _ in range(iterations):
            # Reattribute reads
            logger.info("Re-attributing reads")
            
            # Get chunk_for_*.in.bam/chunk_rev_*.in.bam

            forward_chunks = sorted(glob(str(output_data_path / "chunk_for_*.in.bam")))
            reverse_chunks = sorted(glob(str(output_data_path / "chunk_rev_*.in.bam")))

            # Check if chunks are empty
            for forward_chunk, reverse_chunk in zip(forward_chunks, reverse_chunks):
                if hut.is_empty_alignment(forward_chunk) or hut.is_empty_alignment(reverse_chunk):
                    forward_chunks.remove(forward_chunk)
                    reverse_chunks.remove(reverse_chunk)

            # Reattribute reads
            with Pool(processes = cpus) as pool: # cpus

                res = pool.map(partial(hst.reattribute_reads, mode = sub_mode, restriction_map = output_data_path / RESTRICTION_MAP, output_dir = output_data_path),
                zip(forward_chunks, reverse_chunks))
                pool.close()
                pool.join()


            # hst.reattribute_reads(reads_couple = (forward_in_path, reverse_in_path), mode = sub_mode, output_dir = output_data_path)
            hio.merge_predictions(output_dir = output_data_path, clean = True, cpus = cpus)

            

            hio.build_pairs(bam_for  = "group1.1.out.bam", bam_rev = "group1.2.out.bam", bam_for_rescued  = "group2.1.rescued.bam", bam_rev_rescued = "group2.2.rescued.bam", mode = True, output_dir = output_data_path)
            hio.build_matrix(mode = True, balance = False, output_dir = output_data_path)

            rescued_matrix = hio.load_cooler(output_data_path / RESCUED_MATRIX)

            rescued_matrix_path = output_data_path / RESCUED_MATRIX

            pearson = hst.pearson_score(original_matrix = base_matrix, rescued_matrix = rescued_matrix , markers = indexes)

            if pattern is None or pattern == "-1":
                chromosome_set = [*chromosome, *trans_chromosome] if trans_chromosome is not None else chromosome
            
            else : 
                chromosome_set = [chromosome, *trans_chromosome] if trans_chromosome is not None else chromosome

                rescued_matrix_path = output_data_path / RESCUED_MATRIX
                post_recall_cmd = hev.chromosight_cmd_generator(file = rescued_matrix_path, pattern = pattern, untrend = trend, mode = True, output_dir = output_data_path)

                logger.info("Starting Chromosight post-call")
                sp.run(post_recall_cmd, shell = True)
                
                # TODO : move to top

                df_original = (output_path / "original.tsv").as_posix()
                df_rescued = (output_path / "rescued.tsv").as_posix()

                true_positives = hev.get_TP_table(df_pattern = output_data_path / "original.tsv", df_pattern_recall = output_data_path / "rescued.tsv", chromosome = chromosome, bin_size = bin_size, jitter = jitter, threshold = threshold)
                false_positives = hev.get_FP_table(df_pattern = output_data_path / "original.tsv", df_pattern_recall = output_data_path / "rescued.tsv", chromosome = chromosome, bin_size = bin_size, jitter = jitter, threshold = threshold)
                false_negatives = hev.get_FN_table(df_pattern = output_data_path / "original.tsv", df_pattern_recall = output_data_path / "rescued.tsv", chromosome = chromosome, bin_size = bin_size, jitter = jitter, threshold = threshold)

                # Get scores
                precision = hev.get_precision(df_pattern = df_original, df_pattern_recall = df_rescued, chromosome = chromosome, bin_size = bin_size, jitter = jitter, threshold = threshold)
                recall = hev.get_recall(df_pattern = df_original, df_pattern_recall = df_rescued, chromosome = chromosome, bin_size = bin_size, jitter = jitter, threshold = threshold)
                f1_score = hev.get_f1_score(df_pattern = df_original, df_pattern_recall = df_rescued, chromosome = chromosome, bin_size = bin_size, jitter = jitter, threshold = threshold)

                # Get plots
                hpl.plot_pattern_reconstruction(table = true_positives, original_cool = output_data_path / BASE_MATRIX, rescued_cool = rescued_matrix_path, chromosome = chromosome, threshold = threshold, case  = "true_positives",  output_dir = output_data_path)
                hpl.plot_pattern_reconstruction(table = false_positives, original_cool = output_data_path / BASE_MATRIX, rescued_cool = rescued_matrix_path, chromosome = chromosome, threshold = threshold, case  = "false_positives",  output_dir = output_data_path)
                hpl.plot_pattern_reconstruction(table = false_negatives, original_cool = output_data_path / BASE_MATRIX, rescued_cool = rescued_matrix_path, chromosome = chromosome, threshold = threshold, case  = "false_negatives",  output_dir = output_data_path)

            hpl.plot_benchmark(original_matrix = BASE_MATRIX, depleted_matrix = UNRESCUED_MATRIX, rescued_matrix = RESCUED_MATRIX, chromosomes = chromosome_set, output_dir = output_data_path)

            number_reads = 10

            logger.info(f"Pearson score : {pearson:9.4f} in mode {sub_mode}")

            if not results.exists():
                with open(results, "w") as f_out:
                    f_out.write(header)

                    date = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
                    f_out.write(f"{id_tag}\t{date}\t{chromosome}\t{position}\t{strides}\t{trans_chromosome}\t{trans_position}\t{auto}\t{bins}\t{sub_mode}\t{number_reads}\t{pattern}\t{precision}\t{recall}\t{f1_score}\t{pearson:9.4f}\n")
                    f_out.close()

            else :
                with open(results, "a") as f_out:
                    date = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
                    f_out.write(f"{id_tag}\t{date}\t{chromosome}\t{position}\t{strides}\t{trans_chromosome}\t{trans_position}\t{auto}\t{bins}\t{sub_mode}\t{number_reads}\t{pattern}\t{precision}\t{recall}\t{f1_score}\t{pearson:9.4f}\n")
                    f_out.close()

            logger.info(f"Ending benchmark")

    # # Clean up
    # forward_in_path.unlink()
    # reverse_in_path.unlink()
    # forward_out_path.unlink()
    # reverse_out_path.unlink()
    # (output_data_path / BASE_MATRIX).unlink()
    # (output_data_path / BAM_FOR).unlink()
    # (output_data_path / BAM_REV).unlink()
    # (output_data_path / RESTRICTION_MAP).unlink()
    # (output_data_path / FRAGMENTS).unlink()
    # (output_data_path / CHROMOSOME_SIZES).unlink()
    # (output_data_path / DIST_FRAGS).unlink()
    # (output_data_path / XS).unlink()

    files = [p for p in output_data_path.glob("*")]

    for file in files: 

        if Path(file).suffix in [ ".npy", ".tsv", ".bam", ".pairs"]:

            Path(file).unlink()

    open(flag_file, 'a').close()
    return

