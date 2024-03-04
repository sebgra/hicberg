# import os
from os import getcwd
import subprocess as sp
# import glob
from pathlib import Path
# import tempfile as tmpf
from itertools import combinations

import numpy as np
from numpy.random import choice
import pandas as pd
# from scipy import spatial
from scipy.stats import median_abs_deviation, pearsonr
import bioframe as bf
import pysam as ps

import matplotlib.pyplot as plt
import matplotlib.colors as plc
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

import cooler

import hicberg.io as hio
import hicberg.utils as hut

def get_interval_index(chromosome : str = "", value : int = None, intervals_dict : dict[str, list[(int, int)]] = None, chrom_sizes_dict : dict[str, int] = None) -> dict[str, list[(int, int)]]:
    """

    Parameters
    ----------
    chromosome : [str]
        chromosome linked to value and intervals to search in.
    value : [int]
        Value to search in intervals.
    intervals_dict : [dict]
        Dictionary of intervals as key : chromosome, value : [(low_limit_0, up_limit_0), (low_limit_1, up_limit_1), ..., (low_limit_n, up_limit_n)].
    chrom_sizes_dict : str
        Path to a dictionary containing chromosome sizes as {chromosome : size} saved in .npy format. By default chromosome_sizes.npy

    Returns
    -------
    [dict]
        Lists of intervals (sets) where the value is not in between as element 0 and where the value is  in between as element 1.
        One list per chromosome in a dictionary.
        If the value in not in any interval, return None.
    """

    # chrom_sizes_dict = hio.load_dictionary(chrom_sizes_dict)

    out_dict = {
        k: ([(None, None)], (None, None))
        for k in chrom_sizes_dict.keys()
    }

    for chrom in intervals_dict.keys():

        out_dict[chrom] = intervals_dict[chrom], (None, None)

        for idx, bounds in enumerate(intervals_dict.get(chrom)):

            if bounds is None:

                continue

            if bounds[0] <= value <= bounds[1]:

                out_dict[chrom] = (
                    intervals_dict[chrom][:idx] + intervals_dict[chrom][idx + 1 :],
                    intervals_dict[chrom][idx],
                )

                if len(out_dict[chrom][0]) == 0:
                    out_dict[chrom] = ([(None, None)], intervals_dict[chrom][idx])

    return out_dict


def select_reads(bam_for : str = "group1.1.bam", bam_rev : str = "group1.2.bam", matrix_file: str = "unrescued_map.cool", position : int = 0, chromosome : str | list[str] = "", bin_size : int  = 2000, chrom_sizes_dict : str =  "chromosome_sizes.npy", strides : list[int] = [0],
trans_chromosome :  str = None, output_dir : str = None, trans_position : list[int] = None, nb_bins : int = 1, random : bool = False, auto : int = None) -> dict[str, list[(int, int)]]:
    """
    Select reads from a given matrix and alignment files. Two groups of alignment files are produced (.bam):
    - group1.1.in.bam and group1.2.in.bam : Forward and reverse reads from the source and target genomic interval. These reads are going to be duplicated un each intervals.
    - group1.1.out.bam and group1.2.out.bam : Forward and reverse reads not included the source and target genomic interval. These reads are going to be conserved and not duplicated.

    Parameters
    ----------
    bam_for : str, optional
        Forward alignment file to select reads from, by default "group1.1.bam"
    bam_rev : str, optional
        Reverse alignment file to select reads from, by default "group1.2.bam"
    matrix_file : str, optional
        Name of the matrix from which the selection will be based (.cool format), by default "unrescued_map.cool"
    position : int, optional
        Genomic coordinates (0-based) used to defined the source genomic interval, by default 0
    chromosome : str or list[str], optional
        Chromosome associated to the 0-based coordinates given by "position', by default ""
    bin_size : int, optional
        Resolution of the given matrix, by default 2000
    chrom_sizes_dict : str, optional
        Path to a dictionary containing chromosome sizes as {chromosome : size} saved in .npy format, by default "chromosome_sizes.npy"
    strides : list[int], optional
        List of strides to apply from "position" to define target genomics intervals (in bp.), by default [0]
    trans_chromosome : str, optional
        Chromosomes to consider as trans-chromosomes to define genomic target intervals, by default None
    trans_position : list[int], optional
        Genomic coordinates (0-based) used to defined the trans-chromosomal target genomic interval, by default None
    nb_bins : int, optional
        Number of bins to consider on each side of the defined source and target intervals, by default 1
    random : bool, optional
        Set wether or not the source and target intervals are defined at random, by default False
    auto : int, optional
        Set the number of intervals to get while picking at random, by default None
    output_dir : str, optional
        Path to the folder where to save alignments, by default None

    Returns
    -------
    dict[str, list[(int, int)]]
        Dictionary of intervals as key : chromosome, value : [(low_limit_0, up_limit_0), (low_limit_1, up_limit_1), ..., (low_limit_n, up_limit_n)].

    """   
    
    output_path = Path(output_dir)

    if not output_path.exists():
        
        raise ValueError(f"Output path {output_path} does not exist. Please provide existing output path.")
    
    chrom_sizes_path = output_path / Path(chrom_sizes_dict)

    if not chrom_sizes_path.exists():

        raise ValueError(f"Chromosome sizes file {chrom_sizes_path} does not exist. Please provide existing chromosome sizes file.")
    
    if type(chromosome) == "list":

        chromosome = chromosome[0]

    chs = hio.load_dictionary(chrom_sizes_path)

    if trans_position is not None:
        random = False

    forward_file_path = output_path / Path(bam_for)
    reverse_file_path = output_path / Path(bam_rev)

    if not forward_file_path.exists():
            
        raise ValueError(f"Forward bam file {forward_file_path} does not exist. Please provide existing bam file.")
    
    if not reverse_file_path.exists():

        raise ValueError(f"Reverse bam file {reverse_file_path} does not exist. Please provide existing bam file.")
    
    matrix_file_path = output_path / Path(matrix_file)
    
    if not matrix_file_path.exists():

        raise ValueError(f"Matrix file {matrix_file_path} does not exist. Please provide existing matrix file.")
    
    matrix = hio.load_cooler(matrix_file_path)
    
    # alignment file handlers
    # Create handler for files to parse
    forward_file_handler = ps.AlignmentFile(forward_file_path, "rb")
    reverse_file_handler = ps.AlignmentFile(reverse_file_path, "rb")

    # Create handlers for files to write
    ## Files where selected reads and duplicates are going to be written
    selected_reads_forward = ps.AlignmentFile(
        output_path / "group1.1.in.bam", "wb", template = forward_file_handler
    )
    selected_reads_reverse = ps.AlignmentFile(
        output_path / "group1.2.in.bam", "wb", template = reverse_file_handler
    )

    ## Files where non selected reads are going to be written
    depleted_reads_forward = ps.AlignmentFile(
        output_path / "group1.1.out.bam", "wb", template = forward_file_handler
    )
    depleted_reads_reverse = ps.AlignmentFile(
        output_path / "group1.2.out.bam", "wb", template = reverse_file_handler
    )

    # get acces to dictionary containing chromosomes sizes to pick random position for trans-chromosomal duplication
    cs_disctionary = hio.load_dictionary(chrom_sizes_path)

    if auto is None:

        ## set areas and boundaries for intra-chromosomal duplications
        forward_intervals = [
            get_boundaries(
                position = position + stride,
                bins = bin_size,
                chromosome = chromosome,
                chrom_sizes_dict = chrom_sizes_path,
            )
            for stride in strides
        ]
        reverse_intervals = [
            get_boundaries(
                position = position + stride,
                bins = bin_size,
                chromosome = chromosome,
                chrom_sizes_dict = chrom_sizes_path,
            )
            for stride in strides
        ]

        # Define list of chromosome to target/duplicate read on.
        if type(chromosome) == list and len(chromosome) == 1:

            chromosome = chromosome[0]

        if trans_chromosome is not None:

            list_selected_chromosomes = list(chromosome.split()) + list(trans_chromosome)

        else:

            list_selected_chromosomes = chromosome #.split()

        # adjust intervals width
        if nb_bins > 1:

            forward_intervals = [
                (
                    interval[0] - (nb_bins - 1) * bin_size,
                    interval[1] + (nb_bins - 1) * bin_size,
                )
                for interval in forward_intervals
            ]
            reverse_intervals = [
                (
                    interval[0] - (nb_bins - 1) * bin_size,
                    interval[1] + (nb_bins - 1) * bin_size,
                )
                for interval in reverse_intervals
            ]

        # build dictionary key : chr , value : list of intervals to perform selection on
        dictionary_of_intervals = dict()  # Backup
        dictionary_of_intervals[chromosome] = forward_intervals

        if trans_chromosome is not None:

            for chrom, pos in zip(trans_chromosome, trans_position):

                if random:

                    trans_target_interval = [
                        get_boundaries(
                            position = np.random.randint(
                                low=1, high=cs_disctionary.item().get(chrom)
                            ),
                            bins = bin_size,
                            chromosome = chrom,
                            chrom_sizes_dict=chrom_sizes_path,
                        )
                    ]

                # set areas and boundaries for inter-chromosomal duplications
                else:

                    trans_target_interval = [
                        get_boundaries(
                            position=pos,
                            bins=bin_size,
                            chromosome=chrom,
                            chrom_sizes_dict=chrom_sizes_path,
                        )
                    ]

                if chrom not in dictionary_of_intervals.keys():

                    dictionary_of_intervals[chrom] = trans_target_interval

    if auto is not None:

        dictionary_of_intervals = draw_intervals(chrom_sizes_dict  = chrom_sizes_path, nb_intervals = auto, bins = bin_size)

        # If a randomly selected interval is empty, draw another set of intervals
        while check_emptiness(intervals = dictionary_of_intervals, matrix = matrix):
                
            dictionary_of_intervals = draw_intervals(chrom_sizes_dict  = chrom_sizes_path, nb_intervals = auto, bins = bin_size)

        list_selected_chromosomes = list(dictionary_of_intervals.keys())

    # parse both alignment files to eventually duplicate reads
    for forward_read, reverse_read in zip(forward_file_handler, reverse_file_handler):

        # Default save status
        save = True

        ## Order reads by coordinates
        ordered_forward_read, ordered_reverse_read = hut.get_ordered_reads(
            forward_read, reverse_read
        )

        # Avoid cases where chromosomes are not concerned for duplication
        if (
            forward_read.reference_name not in list_selected_chromosomes
            and reverse_read.reference_name not in list_selected_chromosomes
        ):  # or backp
            depleted_reads_forward.write(forward_read)
            depleted_reads_reverse.write(reverse_read)
            continue

        # Search for intervals where reads potentially belong
        forward_interval_search = get_interval_index(
            chromosome = ordered_forward_read.reference_name,
            value = ordered_forward_read.pos,
            intervals_dict = dictionary_of_intervals,
            chrom_sizes_dict = chs,
        )
        reverse_interval_search = get_interval_index(
            chromosome = ordered_reverse_read.reference_name,
            value = ordered_reverse_read.pos,
            intervals_dict = dictionary_of_intervals,
            chrom_sizes_dict = chs,
        )

        # Case where both forward and reverse reads belong to any of the selected of target interval
        if (
            forward_interval_search[ordered_forward_read.reference_name][1][0]
            is not None
        ) and (
            reverse_interval_search[ordered_reverse_read.reference_name][1][0]
            is not None
        ):

            # Do duplication and write + set save to false
            ## original forward and reverse reads save
            selected_reads_forward.write(ordered_forward_read)
            selected_reads_reverse.write(ordered_reverse_read)

            # Compute forward and reverse shifts
            forward_shift = (
                ordered_forward_read.pos
                - forward_interval_search[ordered_forward_read.reference_name][1][0]
            )  # [1] : interval where the read is initially - [0] : select lower bound of this interval
            reverse_shift = (
                reverse_interval_search[ordered_reverse_read.reference_name][1][1]
                - ordered_reverse_read.pos
            )  # [1] : interval where the read is initially - [1] : select upper bound of this interval

            # Duplicate and save both reads to corresponding files
            for chromosome in forward_interval_search.keys():

                # Check if forward read reference has valid intervals
                if forward_interval_search[chromosome][0][0][0] is None:
                    continue

                else:

                    for forward_interval in forward_interval_search[chromosome][
                        0
                    ]:  # parse all intervals where the read is not initially

                        ordered_forward_read.pos = forward_interval[0] + forward_shift
                        ordered_forward_read.reference_name = chromosome
                        ordered_forward_read.set_tag("XF", "Fake")
                        ordered_forward_read.set_tag("XC", "Case_0")
                        selected_reads_forward.write(ordered_forward_read)

            for chromosome in reverse_interval_search.keys():

                # Check if forward read reference has valid intervals
                if reverse_interval_search[chromosome][0][0][0] is None:
                    continue

                else:

                    for reverse_interval in reverse_interval_search[chromosome][
                        0
                    ]:  # parse all intervals where the read is not initially

                        ordered_reverse_read.pos = reverse_interval[1] - reverse_shift
                        ordered_reverse_read.reference_name = chromosome
                        ordered_reverse_read.set_tag("XF", "Fake")
                        ordered_reverse_read.set_tag("XC", "Case_0")
                        selected_reads_reverse.write(ordered_reverse_read)

            # set save at false
            save = False

        # Case where only forward read belong to any of selected or target interval
        elif (
            forward_interval_search[ordered_forward_read.reference_name][1][0]
            is not None
        ) and (
            reverse_interval_search[ordered_reverse_read.reference_name][1][0] is None
        ):

            ## original forward and reverse reads save
            selected_reads_forward.write(ordered_forward_read)
            selected_reads_reverse.write(ordered_reverse_read)

            # Compute forward  shifts
            forward_shift = (
                ordered_forward_read.pos
                - forward_interval_search[ordered_forward_read.reference_name][1][0]
            )  # [1] : interval where the read is initially - [0] : select lower bound of this interval

            # Do duplication of forward read on all potential interval and save all duplicates
            for chromosome in forward_interval_search.keys():

                # Check if forward read reference has valid intervals
                if forward_interval_search[chromosome][0][0][0] is None:
                    continue

                else:

                    for forward_interval in forward_interval_search[chromosome][
                        0
                    ]:  # parse all intervals where the read is not initially

                        ordered_forward_read.pos = forward_interval[0] + forward_shift
                        ordered_forward_read.reference_name = chromosome
                        ordered_forward_read.set_tag("XF", "Fake")
                        selected_reads_forward.write(ordered_forward_read)

            # set save to false
            save = False

        # Case where only reverse belong to selected or target intervals
        elif (
            forward_interval_search[ordered_forward_read.reference_name][1][0] is None
        ) and (
            reverse_interval_search[ordered_reverse_read.reference_name][1][0]
            is not None
        ):

            ## original forward and reverse reads save
            selected_reads_forward.write(ordered_forward_read)
            selected_reads_reverse.write(ordered_reverse_read)

            # Compute reverse  shifts
            reverse_shift = (
                reverse_interval_search[ordered_reverse_read.reference_name][1][1]
                - ordered_reverse_read.pos
            )  # [1] : interval where the read is initially - [1] : select upper bound of this interval

            # Do duplication of reverse read on all potential interval and save all duplicates
            for chromosome in reverse_interval_search.keys():

                # Check if forward read reference has valid intervals
                if reverse_interval_search[chromosome][0][0][0] is None:
                    continue

                else:

                    for reverse_interval in reverse_interval_search[chromosome][
                        0
                    ]:  # parse all intervals where the read is not initially

                        ordered_reverse_read.pos = reverse_interval[1] - reverse_shift
                        ordered_reverse_read.reference_name = chromosome
                        ordered_reverse_read.set_tag("XF", "Fake")
                        selected_reads_reverse.write(ordered_reverse_read)

            # set save at false
            save = False

        # Write unselected reads to their corresponding files
        if save:

            depleted_reads_forward.write(ordered_forward_read)
            depleted_reads_reverse.write(ordered_reverse_read)

    # close all files
    selected_reads_forward.close()
    selected_reads_reverse.close()

    depleted_reads_forward.close()
    depleted_reads_reverse.close()

    forward_file_handler.close()
    reverse_file_handler.close()

    return dictionary_of_intervals


def select_reads_multithreads(bam_couple : tuple[str, str], matrix_file: str = "unrescued_map.cool", position : int = 0, chromosome : str | list[str] = "", bin_size : int  = 2000, chrom_sizes_dict : str =  "chromosome_sizes.npy", strides : list[int] = [0],
trans_chromosome :  str = None, output_dir : str = None, trans_position : list[int] = None, nb_bins : int = 1, random : bool = False, auto : int = None) -> dict[str, list[(int, int)]]:
    """
    Select reads from a given matrix and alignment files. Two groups of alignment files are produced (.bam):
    - group1.1.in.bam and group1.2.in.bam : Forward and reverse reads from the source and target genomic interval. These reads are going to be duplicated un each intervals.
    - group1.1.out.bam and group1.2.out.bam : Forward and reverse reads not included the source and target genomic interval. These reads are going to be conserved and not duplicated.

    Parameters
    ----------
    bam_for : str, optional
        Forward alignment file to select reads from, by default "group1.1.bam"
    bam_rev : str, optional
        Reverse alignment file to select reads from, by default "group1.2.bam"
    matrix_file : str, optional
        Name of the matrix from which the selection will be based (.cool format), by default "unrescued_map.cool"
    position : int, optional
        Genomic coordinates (0-based) used to defined the source genomic interval, by default 0
    chromosome : str or list[str], optional
        Chromosome associated to the 0-based coordinates given by "position', by default ""
    bin_size : int, optional
        Resolution of the given matrix, by default 2000
    chrom_sizes_dict : str, optional
        Path to a dictionary containing chromosome sizes as {chromosome : size} saved in .npy format, by default "chromosome_sizes.npy"
    strides : list[int], optional
        List of strides to apply from "position" to define target genomics intervals (in bp.), by default [0]
    trans_chromosome : str, optional
        Chromosomes to consider as trans-chromosomes to define genomic target intervals, by default None
    trans_position : list[int], optional
        Genomic coordinates (0-based) used to defined the trans-chromosomal target genomic interval, by default None
    nb_bins : int, optional
        Number of bins to consider on each side of the defined source and target intervals, by default 1
    random : bool, optional
        Set wether or not the source and target intervals are defined at random, by default False
    auto : int, optional
        Set the number of intervals to get while picking at random, by default None
    output_dir : str, optional
        Path to the folder where to save alignments, by default None

    Returns
    -------
    dict[str, list[(int, int)]]
        Dictionary of intervals as key : chromosome, value : [(low_limit_0, up_limit_0), (low_limit_1, up_limit_1), ..., (low_limit_n, up_limit_n)].

    """   
    
    output_path = Path(output_dir)

    if not output_path.exists():
        
        raise ValueError(f"Output path {output_path} does not exist. Please provide existing output path.")
    
    chrom_sizes_path = output_path / Path(chrom_sizes_dict)

    if not chrom_sizes_path.exists():

        raise ValueError(f"Chromosome sizes file {chrom_sizes_path} does not exist. Please provide existing chromosome sizes file.")
    
    if type(chromosome) == "list":

        chromosome = chromosome[0]

    chs = hio.load_dictionary(chrom_sizes_path)

    if trans_position is not None:
        random = False

    # print(f"bam for : {bam_for}")
    # print(f"bam rev : {bam_rev}")
    forward_file_path = output_path / Path(bam_couple[0])
    reverse_file_path = output_path / Path(bam_couple[1])

    if not forward_file_path.exists():
            
        raise ValueError(f"Forward bam file {forward_file_path} does not exist. Please provide existing bam file.")
    
    if not reverse_file_path.exists():

        raise ValueError(f"Reverse bam file {reverse_file_path} does not exist. Please provide existing bam file.")
    
    matrix_file_path = output_path / Path(matrix_file)
    
    if not matrix_file_path.exists():

        raise ValueError(f"Matrix file {matrix_file_path} does not exist. Please provide existing matrix file.")
    
    matrix = hio.load_cooler(matrix_file_path)
    
    # alignment file handlers
    # Create handler for files to parse
    forward_file_handler = ps.AlignmentFile(forward_file_path, "rb")
    reverse_file_handler = ps.AlignmentFile(reverse_file_path, "rb")

    # Create handlers for files to write
    ## Files where selected reads and duplicates are going to be written
    selected_reads_forward = ps.AlignmentFile(
        output_path / f"{forward_file_path.stem}.in.bam", "wb", template = forward_file_handler
    )
    selected_reads_reverse = ps.AlignmentFile(
        output_path / f"{reverse_file_path.stem}.in.bam", "wb", template = reverse_file_handler
    )

    ## Files where non selected reads are going to be written
    depleted_reads_forward = ps.AlignmentFile(
        output_path / f"{forward_file_path.stem}.out.bam", "wb", template = forward_file_handler
    )
    depleted_reads_reverse = ps.AlignmentFile(
        output_path / f"{reverse_file_path.stem}.out.bam", "wb", template = reverse_file_handler
    )

    # get acces to dictionary containing chromosomes sizes to pick random position for trans-chromosomal duplication
    cs_disctionary = hio.load_dictionary(chrom_sizes_path)

    if auto is None:

        ## set areas and boundaries for intra-chromosomal duplications
        forward_intervals = [
            get_boundaries(
                position = position + stride,
                bins = bin_size,
                chromosome = chromosome,
                chrom_sizes_dict = chrom_sizes_path,
            )
            for stride in strides
        ]
        reverse_intervals = [
            get_boundaries(
                position = position + stride,
                bins = bin_size,
                chromosome = chromosome,
                chrom_sizes_dict = chrom_sizes_path,
            )
            for stride in strides
        ]

        # Define list of chromosome to target/duplicate read on.
        if type(chromosome) == list and len(chromosome) == 1:

            chromosome = chromosome[0]

        if trans_chromosome is not None:

            list_selected_chromosomes = list(chromosome.split()) + list(trans_chromosome)

        else:

            list_selected_chromosomes = chromosome #.split()

        # adjust intervals width
        if nb_bins > 1:

            forward_intervals = [
                (
                    interval[0] - (nb_bins - 1) * bin_size,
                    interval[1] + (nb_bins - 1) * bin_size,
                )
                for interval in forward_intervals
            ]
            reverse_intervals = [
                (
                    interval[0] - (nb_bins - 1) * bin_size,
                    interval[1] + (nb_bins - 1) * bin_size,
                )
                for interval in reverse_intervals
            ]

        # build dictionary key : chr , value : list of intervals to perform selection on
        dictionary_of_intervals = dict()  # Backup
        dictionary_of_intervals[chromosome] = forward_intervals

        if trans_chromosome is not None:

            for chrom, pos in zip(trans_chromosome, trans_position):

                if random:

                    trans_target_interval = [
                        get_boundaries(
                            position = np.random.randint(
                                low=1, high=cs_disctionary.item().get(chrom)
                            ),
                            bins = bin_size,
                            chromosome = chrom,
                            chrom_sizes_dict=chrom_sizes_path,
                        )
                    ]

                # set areas and boundaries for inter-chromosomal duplications
                else:

                    trans_target_interval = [
                        get_boundaries(
                            position=pos,
                            bins=bin_size,
                            chromosome=chrom,
                            chrom_sizes_dict=chrom_sizes_path,
                        )
                    ]

                if chrom not in dictionary_of_intervals.keys():

                    dictionary_of_intervals[chrom] = trans_target_interval

    if auto is not None:

        dictionary_of_intervals = draw_intervals(chrom_sizes_dict  = chrom_sizes_path, nb_intervals = auto, bins = bin_size)

        # If a randomly selected interval is empty, draw another set of intervals
        while check_emptiness(intervals = dictionary_of_intervals, matrix = matrix):
                
            dictionary_of_intervals = draw_intervals(chrom_sizes_dict  = chrom_sizes_path, nb_intervals = auto, bins = bin_size)

        list_selected_chromosomes = list(dictionary_of_intervals.keys())

    # parse both alignment files to eventually duplicate reads
    for forward_read, reverse_read in zip(forward_file_handler, reverse_file_handler):

        # Default save status
        save = True

        ## Order reads by coordinates
        ordered_forward_read, ordered_reverse_read = hut.get_ordered_reads(
            forward_read, reverse_read
        )

        # Avoid cases where chromosomes are not concerned for duplication
        if (
            forward_read.reference_name not in list_selected_chromosomes
            and reverse_read.reference_name not in list_selected_chromosomes
        ):  # or backp
            depleted_reads_forward.write(forward_read)
            depleted_reads_reverse.write(reverse_read)
            continue

        # Search for intervals where reads potentially belong
        forward_interval_search = get_interval_index(
            chromosome = ordered_forward_read.reference_name,
            value = ordered_forward_read.pos,
            intervals_dict = dictionary_of_intervals,
            chrom_sizes_dict = chs,
        )
        reverse_interval_search = get_interval_index(
            chromosome = ordered_reverse_read.reference_name,
            value = ordered_reverse_read.pos,
            intervals_dict = dictionary_of_intervals,
            chrom_sizes_dict = chs,
        )

        # Case where both forward and reverse reads belong to any of the selected of target interval
        if (
            forward_interval_search[ordered_forward_read.reference_name][1][0]
            is not None
        ) and (
            reverse_interval_search[ordered_reverse_read.reference_name][1][0]
            is not None
        ):

            # Do duplication and write + set save to false
            ## original forward and reverse reads save
            selected_reads_forward.write(ordered_forward_read)
            selected_reads_reverse.write(ordered_reverse_read)

            # Compute forward and reverse shifts
            forward_shift = (
                ordered_forward_read.pos
                - forward_interval_search[ordered_forward_read.reference_name][1][0]
            )  # [1] : interval where the read is initially - [0] : select lower bound of this interval
            reverse_shift = (
                reverse_interval_search[ordered_reverse_read.reference_name][1][1]
                - ordered_reverse_read.pos
            )  # [1] : interval where the read is initially - [1] : select upper bound of this interval

            # Duplicate and save both reads to corresponding files
            for chromosome in forward_interval_search.keys():

                # Check if forward read reference has valid intervals
                if forward_interval_search[chromosome][0][0][0] is None:
                    continue

                else:

                    for forward_interval in forward_interval_search[chromosome][
                        0
                    ]:  # parse all intervals where the read is not initially

                        ordered_forward_read.pos = forward_interval[0] + forward_shift
                        ordered_forward_read.reference_name = chromosome
                        ordered_forward_read.set_tag("XF", "Fake")
                        ordered_forward_read.set_tag("XC", "Case_0")
                        selected_reads_forward.write(ordered_forward_read)

            for chromosome in reverse_interval_search.keys():

                # Check if forward read reference has valid intervals
                if reverse_interval_search[chromosome][0][0][0] is None:
                    continue

                else:

                    for reverse_interval in reverse_interval_search[chromosome][
                        0
                    ]:  # parse all intervals where the read is not initially

                        ordered_reverse_read.pos = reverse_interval[1] - reverse_shift
                        ordered_reverse_read.reference_name = chromosome
                        ordered_reverse_read.set_tag("XF", "Fake")
                        ordered_reverse_read.set_tag("XC", "Case_0")
                        selected_reads_reverse.write(ordered_reverse_read)

            # set save at false
            save = False

        # Case where only forward read belong to any of selected or target interval
        elif (
            forward_interval_search[ordered_forward_read.reference_name][1][0]
            is not None
        ) and (
            reverse_interval_search[ordered_reverse_read.reference_name][1][0] is None
        ):

            ## original forward and reverse reads save
            selected_reads_forward.write(ordered_forward_read)
            selected_reads_reverse.write(ordered_reverse_read)

            # Compute forward  shifts
            forward_shift = (
                ordered_forward_read.pos
                - forward_interval_search[ordered_forward_read.reference_name][1][0]
            )  # [1] : interval where the read is initially - [0] : select lower bound of this interval

            # Do duplication of forward read on all potential interval and save all duplicates
            for chromosome in forward_interval_search.keys():

                # Check if forward read reference has valid intervals
                if forward_interval_search[chromosome][0][0][0] is None:
                    continue

                else:

                    for forward_interval in forward_interval_search[chromosome][
                        0
                    ]:  # parse all intervals where the read is not initially

                        ordered_forward_read.pos = forward_interval[0] + forward_shift
                        ordered_forward_read.reference_name = chromosome
                        ordered_forward_read.set_tag("XF", "Fake")
                        selected_reads_forward.write(ordered_forward_read)

            # set save to false
            save = False

        # Case where only reverse belong to selected or target intervals
        elif (
            forward_interval_search[ordered_forward_read.reference_name][1][0] is None
        ) and (
            reverse_interval_search[ordered_reverse_read.reference_name][1][0]
            is not None
        ):

            ## original forward and reverse reads save
            selected_reads_forward.write(ordered_forward_read)
            selected_reads_reverse.write(ordered_reverse_read)

            # Compute reverse  shifts
            reverse_shift = (
                reverse_interval_search[ordered_reverse_read.reference_name][1][1]
                - ordered_reverse_read.pos
            )  # [1] : interval where the read is initially - [1] : select upper bound of this interval

            # Do duplication of reverse read on all potential interval and save all duplicates
            for chromosome in reverse_interval_search.keys():

                # Check if forward read reference has valid intervals
                if reverse_interval_search[chromosome][0][0][0] is None:
                    continue

                else:

                    for reverse_interval in reverse_interval_search[chromosome][
                        0
                    ]:  # parse all intervals where the read is not initially

                        ordered_reverse_read.pos = reverse_interval[1] - reverse_shift
                        ordered_reverse_read.reference_name = chromosome
                        ordered_reverse_read.set_tag("XF", "Fake")
                        selected_reads_reverse.write(ordered_reverse_read)

            # set save at false
            save = False

        # Write unselected reads to their corresponding files
        if save:

            depleted_reads_forward.write(ordered_forward_read)
            depleted_reads_reverse.write(ordered_reverse_read)

    # close all files
    selected_reads_forward.close()
    selected_reads_reverse.close()

    depleted_reads_forward.close()
    depleted_reads_reverse.close()

    forward_file_handler.close()
    reverse_file_handler.close()

    return dictionary_of_intervals


def get_intervals_proportions(chrom_sizes_dict : str = "chromosome_sizes.npy", nb_intervals : int = 1) -> dict[str, float]:
    """
    Extract a dictionary containing the number of intervals to draw considering the size of each chromosome.

    Parameters
    ----------
    chrom_sizes_dict : str
        Path to a dictionary containing chromosome sizes as {chromosome : size} saved in .npy format. By default chromosome_sizes.npy
    nb_intervals : int, optional
        Number of intervals to draw, by default 1

    Returns
    -------
    dict[str, int]
        Dictionary containing proportion by intervals as {chromosome : proportion}.
    """ 
    
    chrom_sizes_path = Path(chrom_sizes_dict)
    chrom_sizes = hio.load_dictionary(chrom_sizes_path)

    # Get genome global length
    tot_genome_length = np.sum(list(chrom_sizes.values()))

    # Compute relative proportion of each chromosome through genome
    chr_proportions  = {k : (v / tot_genome_length)  for (k, v) in zip(chrom_sizes.keys(), chrom_sizes.values())}

    # Draw chromosomes considering their relative proportion in genome
    proportions_choice = choice(a = list(chr_proportions.keys()), size = nb_intervals, replace = True, p = list(chr_proportions.values()))

    # Get and return chromosomes picked and counts
    unique, counts = np.unique(proportions_choice, return_counts=True)

    return dict(zip(unique, counts))


def get_chromosomes_intervals(bins : int = 2000, chrom_sizes_dict : str = "chromosome_sizes.npy", chromosome : str = "") -> list[(int, int)]:
    """
    Get all possible intervals from a given chromosome considering a chromosome sizes dictionary.

    Parameters
    ----------
    bins : int, optional
        Size of the desired bin, by default 2000., by default 2000
    chrom_sizes_dict : str
        Path to a dictionary containing chromosome sizes as {chromosome : size} saved in .npy format. By default chromosome_sizes.npy
    chromosome : str, optional
        Chromosome to get intervals from, by default ""

    Returns
    -------
    list[(int, int)]
        List containing the intervals as [(start, end), ...]
    """
    chromosome_size_dictionary = hio.load_dictionary(chrom_sizes_dict)

    nb_bins = chromosome_size_dictionary.get(chromosome) // bins + 1 

    intervals = [(i * bins, (i + 1) * bins) for i in range(2, nb_bins - 1)] # -1 to exclude last bin to avoid anchor position error while constructing maps with cooler

    return intervals
    
def draw_intervals(nb_intervals : int = 1, bins : int = 2000, chrom_sizes_dict : str = "chromosome_sizes.npy") -> dict[str, list[tuple[int, int]]]:
    """
    Draw intervals from a given chromosome sizes dictionary.

    Parameters
    ----------
    nb_intervals : int, optional
        Number of intervals to draw, by default 1
    bins : int, optional
        Size of the desired bin, by default 2000., by default 2000
    chrom_sizes_dict : str
        Path to a dictionary containing chromosome sizes as {chromosome : size} saved in .npy format. By default chromosome_sizes.npy

    Returns
    -------
    dict[str, list[tuple[int, int]]]
        Dictionary containing the intervals as {chromosome : [(start, end), ...]}.
    """

    # get porportions by chromosome
    intervals_proportion = get_intervals_proportions(chrom_sizes_dict = chrom_sizes_dict, nb_intervals = nb_intervals)

    # Define placeholder for intervals
    selected_intervals = dict()

    for chrom in intervals_proportion.keys():

        candidate_intervals = get_chromosomes_intervals(chrom_sizes_dict = chrom_sizes_dict, chromosome = chrom, bins = bins)
        picked_interval_indexes = np.random.choice(np.arange(len(candidate_intervals)), size = intervals_proportion[chrom], replace = False)

        selected_intervals[chrom] = [candidate_intervals[idx] for idx in picked_interval_indexes]
    
    return selected_intervals

def get_boundaries(position : int = None, bins : int = 2000, chromosome : str | list[str] = None, chrom_sizes_dict : str = "chromosome_sizes.npy") -> tuple[int, int]:
    """
    Return boundaries surrounding position considering regular binning of bins.

    Parameters
    ----------
    position : int, optional
        Position to get surrounding boundaries from, by default None
    bins : int, optional
        Size of the desired bin, by default 2000., by default 2000
    chromosome : str or list[str], optional
        Chromosomes associated to positions to get surrounding boundaries from, by default None
    chrom_sizes_dict : str
        Path to a dictionary containing chromosome sizes as {chromosome : size} saved in .npy format. By default chromosome_sizes.npy

    Returns
    -------
    tuple[int, int]
        Tuple containing the boundaries of the position : (lower_bound, upper_bound).
    """

    chrom_sizes_path = Path(chrom_sizes_dict)
    chrom_sizes = hio.load_dictionary(chrom_sizes_path)

    if type(chromosome) == list and len(chromosome) == 1:  

        chromosome = chromosome[0]

    area_to_search = np.append(
        np.arange(2 * bins, chrom_sizes.get(chromosome) - 3 * bins, bins), # np.arange(0, cs_disctionary.item().get(chromosome) - bin_size, bin_size)
        chrom_sizes.get(chromosome),
    )

    if position % bins == 0:

        before_index = np.searchsorted(area_to_search, position, side="left")
        after_index = np.searchsorted(area_to_search, position, side="right")

    else:

        before_index = np.searchsorted(area_to_search, position, side="left") - 1
        after_index = np.searchsorted(area_to_search, position, side="right")

    lower_bound, upper_bound = area_to_search[before_index], area_to_search[after_index]

    return (lower_bound, upper_bound)

def check_emptiness(intervals : dict[str, list[(int, int)]], matrix : cooler.Cooler = None) -> bool:
    """
    Check if intervals are empty in a given matrix.

    Parameters
    ----------
    intervals : dict[str, list[(int, int)]
        Lists of intervals (sets) where the value is not in between as element 0 and where the value is  in between as element 1.
        One list per chromosome in a dictionary.
    matrix : cooler.Cooler, optional
        Cooler Hi-C matrix to be checked for emptiness, by default None

    Returns
    -------
    bool
        Returns True if one of the intervals are empty in the matrix, False otherwise.
    """

    for chrom in intervals.keys():
        for interval in intervals[chrom]:
            if matrix.matrix(balance = False).fetch((chrom, interval[0], interval[1])).sum() == 0:
                return True

    return False

def get_bin_indexes(matrix : str = None, dictionary : dict = None):
    """
    Return the bins indexes (in the cooler matrix) corresponding to the genomic coordinates in the dictionary.

    Parameters
    ----------
    matrix : str, optional
        Path to a cooler matrix file, by default None
    dictionary : dict, optional
        Dictionary of genomic coordinates where reads are going to be selected for duplication, by default None
    """    

    bin_list = []
    for chrm in dictionary.keys():
        for start, end in dictionary[chrm]:
            bin_list.append(int(matrix.bins().fetch(chrm + ':' + str(start) + '-' + str(end)).index[0]))

    return bin_list

def arange_multi(starts, stops=None, lengths=None):
    """
    Create concatenated ranges of integers for multiple start/length.
    Parameters
    ----------
    starts : numpy.ndarray
        Starts for each range
    stops : numpy.ndarray
        Stops for each range
    lengths : numpy.ndarray
        Lengths for each range. Either stops or lengths must be provided.
    Returns
    -------
    concat_ranges : numpy.ndarray
        Concatenated ranges.
    Notes
    -----
    See the following illustrative example:
    starts = np.array([1, 3, 4, 6])
    stops = np.array([1, 5, 7, 6])
    print arange_multi(starts, lengths)
    >>> [3 4 4 5 6]
    From: https://codereview.stackexchange.com/questions/83018/vectorized-numpy-version-of-arange-with-multiple-start-stop
    """

    if (stops is None) == (lengths is None):
        raise ValueError("Either stops or lengths must be provided!")

    if lengths is None:
        lengths = stops - starts

    if np.isscalar(starts):
        starts = np.full(len(stops), starts)

    # Repeat start position index length times and concatenate
    cat_start = np.repeat(starts, lengths)

    # Create group counter that resets for each start/length
    cat_counter = np.arange(lengths.sum()) - np.repeat(
        lengths.cumsum() - lengths, lengths
    )

    # Add group counter to group specific starts
    cat_range = cat_start + cat_counter

    return cat_range

def _overlap_intervals_legacy(starts1, ends1, starts2, ends2, closed=False, sort=False):
    """
    Take two sets of intervals and return the indices of pairs of overlapping intervals.
    Parameters
    ----------
    starts1, ends1, starts2, ends2 : numpy.ndarray
        Interval coordinates. Warning: if provided as pandas.Series, indices
        will be ignored.
    closed : bool
        If True, then treat intervals as closed and report single-point overlaps.
    Returns
    -------
    overlap_ids : numpy.ndarray
        An Nx2 array containing the indices of pairs of overlapping intervals.
        The 1st column contains ids from the 1st set, the 2nd column has ids
        from the 2nd set.
    """

    # for vec in [starts1, ends1, starts2, ends2]:
    #     if issubclass(type(vec), pd.core.series.Series):
    #         warnings.warn(
    #             "One of the inputs is provided as pandas.Series and its index "
    #             "will be ignored.",
    #             SyntaxWarning,
    #         )

    starts1 = np.asarray(starts1)
    ends1 = np.asarray(ends1)
    starts2 = np.asarray(starts2)
    ends2 = np.asarray(ends2)

    # Concatenate intervals lists
    n1 = len(starts1)
    n2 = len(starts2)
    starts = np.concatenate([starts1, starts2])
    ends = np.concatenate([ends1, ends2])

    # Encode interval ids as 1-based,
    # negative ids for the 1st set, positive ids for 2nd set
    ids = np.concatenate([-np.arange(1, n1 + 1), np.arange(1, n2 + 1)])

    # Sort all intervals together
    order = np.lexsort([ends, starts])
    starts, ends, ids = starts[order], ends[order], ids[order]

    # Find interval overlaps
    match_starts = np.arange(0, n1 + n2)
    match_ends = np.searchsorted(starts, ends, "right" if closed else "left")

    # Ignore self-overlaps
    match_mask = match_ends > match_starts + 1
    match_starts, match_ends = match_starts[match_mask], match_ends[match_mask]

    # Restore
    overlap_ids = np.vstack(
        [
            np.repeat(ids[match_starts], match_ends - match_starts - 1),
            ids[arange_multi(match_starts + 1, match_ends)],
        ]
    ).T

    # Drop same-set overlaps
    overlap_ids = overlap_ids[overlap_ids[:, 0] * overlap_ids[:, 1] <= 0]

    # Flip overlaps, such that the 1st column contains ids from the 1st set,
    # the 2nd column contains ids from the 2nd set.
    overlap_ids.sort(axis=-1)

    # Restore original indexes,
    overlap_ids[:, 0] = overlap_ids[:, 0] * (-1) - 1
    overlap_ids[:, 1] = overlap_ids[:, 1] - 1

    # Sort overlaps according to the 1st
    if sort:
        overlap_ids = overlap_ids[np.lexsort([overlap_ids[:, 1], overlap_ids[:, 0]])]

    return overlap_ids

def overlap_intervals(starts1, ends1, starts2, ends2, closed=False, sort=False):
    """
    Take two sets of intervals and return the indices of pairs of overlapping intervals.

    Parameters
    ----------
    starts1, ends1, starts2, ends2 : numpy.ndarray
        Interval coordinates. Warning: if provided as pandas.Series, indices
        will be ignored.

    closed : bool
        If True, then treat intervals as closed and report single-point overlaps.
    Returns
    -------
    overlap_ids : numpy.ndarray
        An Nx2 array containing the indices of pairs of overlapping intervals.
        The 1st column contains ids from the 1st set, the 2nd column has ids
        from the 2nd set.

    """

    # for vec in [starts1, ends1, starts2, ends2]:
    #     if issubclass(type(vec), pd.core.series.Series):
    #         warnings.warn(
    #             "One of the inputs is provided as pandas.Series and its index "
    #             "will be ignored.",
    #             SyntaxWarning,
    #         )

    starts1 = np.asarray(starts1)
    ends1 = np.asarray(ends1)
    starts2 = np.asarray(starts2)
    ends2 = np.asarray(ends2)

    # Concatenate intervals lists
    n1 = len(starts1)
    n2 = len(starts2)
    ids1 = np.arange(0, n1)
    ids2 = np.arange(0, n2)

    # Sort all intervals together
    order1 = np.lexsort([ends1, starts1])
    order2 = np.lexsort([ends2, starts2])
    starts1, ends1, ids1 = starts1[order1], ends1[order1], ids1[order1]
    starts2, ends2, ids2 = starts2[order2], ends2[order2], ids2[order2]

    # Find interval overlaps
    match_2in1_starts = np.searchsorted(starts2, starts1, "left")
    match_2in1_ends = np.searchsorted(starts2, ends1, "right" if closed else "left")
    # "right" is intentional here to avoid duplication
    match_1in2_starts = np.searchsorted(starts1, starts2, "right")
    match_1in2_ends = np.searchsorted(starts1, ends2, "right" if closed else "left")

    # Ignore self-overlaps
    match_2in1_mask = match_2in1_ends > match_2in1_starts
    match_1in2_mask = match_1in2_ends > match_1in2_starts
    match_2in1_starts, match_2in1_ends = (
        match_2in1_starts[match_2in1_mask],
        match_2in1_ends[match_2in1_mask],
    )
    match_1in2_starts, match_1in2_ends = (
        match_1in2_starts[match_1in2_mask],
        match_1in2_ends[match_1in2_mask],
    )

    # Generate IDs of pairs of overlapping intervals
    overlap_ids = np.block(
        [
            [
                np.repeat(ids1[match_2in1_mask], match_2in1_ends - match_2in1_starts)[
                    :, None
                ],
                ids2[arange_multi(match_2in1_starts, match_2in1_ends)][:, None],
            ],
            [
                ids1[arange_multi(match_1in2_starts, match_1in2_ends)][:, None],
                np.repeat(ids2[match_1in2_mask], match_1in2_ends - match_1in2_starts)[
                    :, None
                ],
            ],
        ]
    )

    if sort:
        # Sort overlaps according to the 1st
        overlap_ids = overlap_ids[np.lexsort([overlap_ids[:, 1], overlap_ids[:, 0]])]

    return overlap_ids

def overlap_intervals_outer(starts1, ends1, starts2, ends2, closed=False):
    """
    Take two sets of intervals and return the indices of pairs of overlapping intervals,
    as well as the indices of the intervals that do not overlap any other interval.

    Parameters
    ----------
    starts1, ends1, starts2, ends2 : numpy.ndarray
        Interval coordinates. Warning: if provided as pandas.Series, indices
        will be ignored.

    closed : bool
        If True, then treat intervals as closed and report single-point overlaps.

    Returns
    -------
    overlap_ids : numpy.ndarray
        An Nx2 array containing the indices of pairs of overlapping intervals.
        The 1st column contains ids from the 1st set, the 2nd column has ids
        from the 2nd set.

    no_overlap_ids1, no_overlap_ids2 : numpy.ndarray
        Two 1D arrays containing the indices of intervals in sets 1 and 2
        respectively that do not overlap with any interval in the other set.

    """

    ovids = overlap_intervals(starts1, ends1, starts2, ends2, closed=closed)
    no_overlap_ids1 = np.where(
        np.bincount(ovids[:, 0], minlength=starts1.shape[0]) == 0
    )[0]
    no_overlap_ids2 = np.where(
        np.bincount(ovids[:, 1], minlength=starts2.shape[0]) == 0
    )[0]
    return ovids, no_overlap_ids1, no_overlap_ids2

def intersect2D(a, b):
    """
    Find row intersection between 2D numpy arrays, a and b.
    Returns another numpy array with shared rows
    """
    
    return np.array([x for x in set(tuple(x) for x in a) & set(tuple(x) for x in b)])

def get_TP_table(df_pattern : pd.DataFrame, df_pattern_recall : pd.DataFrame, chromosome : str, bin_size : int, jitter : int  = 0, threshold : float = None) -> pd.DataFrame:
    """
    Take dataframe of pattern call (Chromosight) before and after reconstruction
    and return table of retieved patterns 

    Parameters
    ----------
    df_pattern : [dataframe]
        Pandas DataFrame given by chromosight before HiC map reconstruction.
    df_pattern_recall : [dataframe]
        Pandas DataFrame given by chromosight after HiC map reconstruction.
    chromosome : [str]
        chromosome to select position within.
    bin_size : [int]
        bin size used to construct the HiC map.
    jitter : [int]
        jitter to apply to the pattern recall table to allow overlapping with the pattern table. By default 0.
    threshold : [float]
        threshold to apply to the pattern table to select patterns to consider. By default None.

    Return
    ----------

    dataframe : [pandas DataFrame]
        Table containing components before and after pattern recall (True positives)
    """

    df_pattern = pd.read_csv(df_pattern, sep = "\t", header = 0)
    df_pattern_recall = pd.read_csv(df_pattern_recall, sep = "\t", header = 0)

    if threshold is None:
        # Selection of chromosomes of interest
        df_1 = df_pattern[df_pattern["chrom1"] == chromosome]
        df_2 = df_pattern_recall[df_pattern_recall["chrom1"] == chromosome]

    else :

        df_1 = df_pattern.query("chrom1 == @chromosome and score > @threshold")
        df_2 = df_pattern_recall.query("chrom1 == @chromosome and score > @threshold")

    if jitter != 0:

        df_2["start1"] = df_2["start1"].sub(jitter * bin_size) 
        df_2["end1"] = df_2["end1"].add(jitter * bin_size)

    # Getting overlaps indexes for left side (start1; end1) for both before and afetr pattern recall
    before_after_left = _overlap_intervals_legacy(starts1 = df_1["start1"], ends1 = df_1["end1"], starts2 = df_2["start1"], ends2 = df_2["end1"], closed=False, sort=False)
    # Getting overlaps indexes for left side (start2; end2) for both before and afetr pattern recall
    before_after_right = _overlap_intervals_legacy(starts1 = df_1["start2"], ends1 = df_1["end2"], starts2 = df_2["start2"], ends2 = df_2["end2"], closed=False, sort=False)

    # Getting intersection of indeces (left/right) common to  before and after tables
    selection = intersect2D(a = before_after_left , b = before_after_right)

    # Pick lines through indexes to construct final table
    lines_before = list()
    lines_after = list()
    
    for i in range(selection.shape[0]):
        lines_before.append(df_pattern[df_pattern['chrom1'] == chromosome].iloc[selection[i][0]])
        lines_after.append(df_pattern_recall[df_pattern_recall['chrom1'] == chromosome].iloc[selection[i][1]])

    final_before = pd.DataFrame(lines_before)
    final_before.rename(columns = {"chrom1" : "chrom1_before", "start1" : "start1_before", "end1": "end1_before" , "chrom2": "chrom2_before" , "start2": "start2_before" , "end2" : "end2_before" , "score" : "score_before"}, inplace = True)
    final_before.reset_index(drop=True, inplace=True)

    final_after = pd.DataFrame(lines_after)
    final_after.rename(columns = {"chrom1" : "chrom1_after", "start1" : "start1_after", "end1": "end1_after" , "chrom2": "chrom2_after" , "start2": "start2_after" , "end2" : "end2_after" , "score" : "score_after"}, inplace = True)
    final_after.reset_index(drop=True, inplace=True)

    true_positives_table = pd.concat([final_before,final_after], axis=1, ignore_index=False)
    true_positives_table["score"] = true_positives_table["score_after"]
    true_positives_table["start1"] = true_positives_table["start1_after"]
    true_positives_table["start2"] = true_positives_table["start2_after"]

    return true_positives_table

def get_FN_table(df_pattern : pd.DataFrame, df_pattern_recall : pd.DataFrame, chromosome : str, bin_size : int, jitter : int = 0, threshold : float = None) -> pd.DataFrame:
    """
    Take dataframe of pattern call (Chromosight) before and after reconstruction
    and return table of retieved patterns 

    Parameters
    ----------
    df_pattern : [dataframe]
        Pandas DataFrame given by chromosight before HiC map reconstruction.
    df_pattern_recall : [dataframe]
        Pandas DataFrame given by chromosight after HiC map reconstruction.
    chromosome : [str]
        chromosome to select position within.
    bin_size : [int]
        bin size used to construct the HiC map.
    jitter : [int]
        jitter to apply to the pattern recall table to allow overlapping with the pattern table. By default 0.
    threshold : [float]
        threshold to apply to the pattern table to select patterns to consider. By default None.

    Return
    ----------

    dataframe : [pandas DataFrame]
        Table containing components before and after pattern recall (True positives)
    """

    df_pattern = pd.read_csv(df_pattern, sep = "\t", header = 0)
    df_pattern_recall = pd.read_csv(df_pattern_recall, sep = "\t", header = 0)

    if threshold is None:
        # Selection of chromosomes of interest
        df_1 = df_pattern[df_pattern["chrom1"] == chromosome]
        df_2 = df_pattern_recall[df_pattern_recall["chrom1"] == chromosome]

    else :

        df_1 = df_pattern.query("chrom1 == @chromosome and score > @threshold")
        df_2 = df_pattern_recall.query("chrom1 == @chromosome and score > @threshold")
    
    if jitter != 0:

        df_2["start1"] = df_2["start1"].sub(jitter * bin_size) 
        df_2["end1"] = df_2["end1"].add(jitter * bin_size)

    _left, FN_left, FP_left = overlap_intervals_outer(starts1 = df_1["start1"], ends1 = df_1["end1"], starts2 = df_2["start1"], ends2 = df_2["end1"], closed=False)
    _right, FN_right, FP_right = overlap_intervals_outer(starts1 = df_1["start2"], ends1 = df_1["end2"], starts2 = df_2["start2"], ends2 = df_2["end2"], closed=False)

    FN_intersection = np.intersect1d(FN_left, FN_right)

    # Pick lines through indexes to construct final table
    lines_FN = list()

    for i in FN_intersection:
        lines_FN.append(df_1[df_1['chrom1'] == chromosome].iloc[i])

    final_FN = pd.DataFrame(lines_FN)

    return final_FN

def get_FP_table(df_pattern: pd.DataFrame, df_pattern_recall : pd.DataFrame, chromosome : str, bin_size : int, jitter : int = 0, threshold : float = None):
    """
    Take dataframe of pattern call (Chromosight) before and after reconstruction
    and return table of retieved patterns 

    Parameters
    ----------
    df_pattern : [dataframe]
        Pandas DataFrame given by chromosight before HiC map reconstruction.
    df_pattern_recall : [dataframe]
        Pandas DataFrame given by chromosight after HiC map reconstruction.
    chromosome : [str]
        chromosome to select position within.
    bin_size : [int]
        bin size used to construct the HiC map.
    jitter : [int]
        jitter to apply to the pattern recall table to allow overlapping with the pattern table. By default 0.
    threshold : [float]
        threshold to apply to the pattern table to select patterns to consider. By default None.

    Return
    ----------

    dataframe : [pandas DataFrame]
        Table containing components before and after pattern recall (True positives)
    """

    df_pattern = pd.read_csv(df_pattern, sep = "\t", header = 0)
    df_pattern_recall = pd.read_csv(df_pattern_recall, sep = "\t", header = 0)

    if threshold is None:
        # Selection of chromosomes of interest
        df_1 = df_pattern[df_pattern["chrom1"] == chromosome]
        df_2 = df_pattern_recall[df_pattern_recall["chrom1"] == chromosome]

    else :

        df_1 = df_pattern.query("chrom1 == @chromosome and score > @threshold")
        df_2 = df_pattern_recall.query("chrom1 == @chromosome and score > @threshold")
    
    if jitter != 0:

        df_2["start1"] = df_2["start1"].sub(jitter * bin_size) 
        df_2["end1"] = df_2["end1"].add(jitter * bin_size)

    _left, FN_left, FP_left = overlap_intervals_outer(starts1 = df_1["start1"], ends1 = df_1["end1"], starts2 = df_2["start1"], ends2 = df_2["end1"], closed=False)
    _right, FN_right, FP_right = overlap_intervals_outer(starts1 = df_1["start2"], ends1 = df_1["end2"], starts2 = df_2["start2"], ends2 = df_2["end2"], closed=False)

    FP_intersection = np.intersect1d(FP_left, FP_right)

    # Pick lines through indexes to construct final table
    lines_FP = list()

    if len(lines_FP) == 0 :

        return None

    else :

        for i in FP_intersection:
            lines_FP.append(df_1[df_1['chrom1'] == chromosome].iloc[i])

        final_FP = pd.DataFrame(lines_FP)

        return final_FP

def get_recall(df_pattern : pd.DataFrame, df_pattern_recall : pd.DataFrame, chromosome : str, bin_size : int, jitter : int = 0, threshold : float = None) -> float:
    """
    Take dataframe of pattern call (Chromosight) before and after reconstruction
    and return recall score corresponding to the number of true positives over the number of true positives and false negatives.

    Parameters
    ----------
    df_pattern : [dataframe]
        Pandas DataFrame given by chromosight before HiC map reconstruction.
    df_pattern_recall : [dataframe]
        Pandas DataFrame given by chromosight after HiC map reconstruction.
    chromosome : [str]
        chromosome to select position within.
    bin_size : [int]
        bin size used to construct the HiC map.
    jitter : [int]
        jitter to apply to the pattern recall table to allow overlapping with the pattern table. By default 0.
    threshold : [float]
        threshold to apply to the pattern table to select patterns to consider. By default None.


    Returns
    -------
    recall : float
        number of true positives over the number of true positives and false negatives
    """

    TP = get_TP_table(df_pattern = df_pattern , df_pattern_recall = df_pattern_recall , chromosome = chromosome, bin_size = bin_size, jitter = jitter, threshold = threshold).shape[0]
    FN = get_FN_table(df_pattern = df_pattern , df_pattern_recall = df_pattern_recall , chromosome = chromosome, bin_size = bin_size, jitter = jitter, threshold = threshold).shape[0]

    return  TP / (TP + FN)

def get_precision(df_pattern : pd.DataFrame, df_pattern_recall : pd.DataFrame, chromosome : str, bin_size : int, jitter : int = 0, threshold : float = None) -> float:
    """
    Take dataframe of pattern call (Chromosight) before and after reconstruction
    and return precision score corresponding to the number of true positives over the number of true positives and false positives.

    Parameters
    ----------
    df_pattern : [dataframe]
        Pandas DataFrame given by chromosight before HiC map reconstruction.
    df_pattern_recall : [dataframe]
        Pandas DataFrame given by chromosight after HiC map reconstruction.
    chromosome : [str]
        chromosome to select position within.
    bin_size : [int]
        bin size used to construct the HiC map.
    jitter : [int]
        jitter to apply to the pattern recall table to allow overlapping with the pattern table. By default 0.
    threshold : [float]
        threshold to apply to the pattern table to select patterns to consider. By default None.


    Returns
    -------
    precision : float
        number of true positives over the number of true positives and false positives
    """

    TP = get_TP_table(df_pattern = df_pattern , df_pattern_recall = df_pattern_recall , chromosome = chromosome, bin_size = bin_size, jitter = jitter, threshold = threshold).shape[0]

    if get_FP_table(df_pattern = df_pattern , df_pattern_recall = df_pattern_recall , chromosome = chromosome, bin_size = bin_size, jitter = jitter, threshold = threshold) is None:
        FP = 0
        return  TP / (TP + FP)
        
    FP = get_FP_table(df_pattern = df_pattern , df_pattern_recall = df_pattern_recall , chromosome = chromosome, bin_size = bin_size, jitter = jitter, threshold = threshold).shape[0]

    return  TP / (TP + FP)


def get_f1_score(df_pattern : pd.DataFrame, df_pattern_recall : pd.DataFrame, chromosome : str, bin_size : int, jitter : int = 0, threshold : float = None) -> float:
    """
    Take dataframe of pattern call (Chromosight) before and after reconstruction
    and return f1 score i.e. 2 * ((precision * recall) / (precision + recall)).

    Parameters
    ----------
    df_pattern : [dataframe]
        Pandas DataFrame given by chromosight before HiC map reconstruction.
    df_pattern_recall : [dataframe]
        Pandas DataFrame given by chromosight after HiC map reconstruction.
    chromosome : [str]
        chromosome to select position within.
    bin_size : [int]
        bin size used to construct the HiC map.
    jitter : [int]
        jitter to apply to the pattern recall table to allow overlapping with the pattern table. By default 0.
    threshold : [float]
        threshold to apply to the pattern table to select patterns to consider. By default None.


    Returns
    -------
    f1_score : float
        2 * ((precision * recall) / (precision + recall))
    """
    recall = get_recall(df_pattern = df_pattern , df_pattern_recall = df_pattern_recall , chromosome = chromosome, bin_size = bin_size, jitter = jitter, threshold = threshold)
    precision = get_precision(df_pattern = df_pattern , df_pattern_recall = df_pattern_recall , chromosome = chromosome, bin_size = bin_size, jitter = jitter, threshold = threshold)

    return 2 * ((precision * recall) / (precision + recall))

def get_top_pattern(file : str = None, top : int = 10, threshold :float = 0.0, chromosome : str = None) -> pd.DataFrame:
    """
    Get top patterns from a dataframe

    Parameters
    ----------
    df : pd.DataFrame, optional
        Dataframe containing patterns given by Chromosight, by default None
    top : int, optional
        Percentage of top patterns to get, by default 10
    threshold : float, optional
        Pattern Pearson score to consider to select pattern, by default 0.0
    chromosome : str, optional
        Chromosome to consider, by default None

    Returns
    -------
    pd.DataFrame
        Dataframe containing top percentage patterns.
    """
    df = pd.read_csv(file, sep = "\t", header = 0)
    df = df.query(f"score > {threshold}")


    if chromosome is not None:

        df = df.query(f"chrom1 == '{chromosome}' and chrom2 == '{chromosome}'")
    top_factor = (df.shape[0] * top) // 100
    df_top = df.sort_values(by='score', ascending=False).head(top_factor).reset_index(drop=True)

    return df_top

def hicberg_benchmark_cmd_generator(file : str = None, top : int = 10, threshold :float = 0.0, chromosome : str = None, mode : str = "full",  genome : str = None, bins : int = 0, output : str = None) -> str:
    """
    Get top patterns from a dataframe

    Parameters
    ----------
    df : pd.DataFrame, optional
        Dataframe containing patterns given by Chromosight, by default None
    top : int, optional
        Percentage of top patterns to get, by default 10
    threshold : float, optional
        Pattern Pearson score to consider to select pattern, by default 0.0
    chromosome : str, optional
        Chromosome to consider, by default None
    mode : str, optional
        Mode to consider for hicberg benchmark, by default "full"
    genome : str, optional
        Path to the genome to consider for hicberg benchmark, by default None
    bins : int, optional
        Number of bins to consider for hicberg benchmark, by default 0
    output : str, optional
        Path to the output folder to consider for hicberg benchmark, by default None

    Returns
    -------
    str
        HiC-BERG command line to run to evaluate reconstruction after pattern discarding
    """
    

    df = get_top_pattern(file = file, top = top, threshold = threshold, chromosome = chromosome).sort_values(by='start1', ascending=True)

    stride = [str(df.iloc[i].start1 - df['start1'].min()) for i in range(0, df.shape[0])]

    cmd = f"hicberg benchmark -c {chromosome} -p {df['start1'].min()} -s {','.join(list(stride))} -m {mode} -g {genome} -b {bins} -o {output}"

    return cmd


def chromosight_cmd_generator(file : str = None, pattern : str = "loops", untrend : bool = True, mode : bool = False,  output_dir : str = None) -> str:
    """
    Generate chromosight command line to run pattern detection.

    Parameters
    ----------
    file : str, optional
        Path to Hi-C balanced contact matrix in .cool format , by default None
    pattern : str, optional
        Pattern to detect, by default "loops"
    untrend : bool, optional
        Set if map has to be detrended, by default True
    mode : bool, optional
        Set either the detection has to be performed before or after map reconstruction, by default False
    output_dir : str, optional
        Path to the folder where to save alignments, by default None

    Returns
    -------
    str
        Chromosight command line to run to run pattern detection
    """    

    output_path = Path(output_dir)

    if not output_path.exists():
        raise ValueError(f"Output path {output_path} does not exist. Please provide existing output path.")
    
    matrix_path = Path(file)

    if matrix_path.suffix != ".cool":
        raise ValueError(f"Matrix path {matrix_path} is not a .cool file. Please provide a .cool file.")
    
    if not matrix_path.is_file():
        raise ValueError(f"Matrix path {matrix_path} does not exist. Please provide existing matrix path.")
    
    prefix = "original" if not mode else "rescued"

    if untrend :
        cmd = f"chromosight detect -P {pattern} -T {file} {output_path / prefix}"
    
    else :
        cmd = f"chromosight detect -P {pattern} {file} {output_path / prefix}"

    return cmd




