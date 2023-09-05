# import os
import subprocess as sp
# import glob
from pathlib import Path
# import tempfile as tmpf
import shutil
from itertools import combinations

import numpy as np
from numpy.random import choice

# from scipy import spatial
from scipy.stats import median_abs_deviation, pearsonr
import bioframe as bf
# import pysam
import pysam as ps

import matplotlib.pyplot as plt
# import matplotlib.colors
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

import cooler

import hicberg.io as hio
import hicberg.utils as hut

# TODO : Implement associated test

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


def select_reads(bam_for :str = "group1.1.bam", bam_rev : str = "group1.2.bam", matrix_file: str = "unrescued_map.cool", position : int = 0, chromosome : str | list[str] = "", bin_size : int  = 2000, chrom_sizes_dict : str =  "chromosome_sizes.npy", strides : list[int] = [0],
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
        
        raise ValueError(f"Output path {output_path} does not exist. Please provide existing ouput path.")
    
    chrom_sizes_path = output_path / Path(chrom_sizes_dict)

    print(f"Chromosome sizes path: {chrom_sizes_path}")

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
        # check intervals overlapping
        # print("----------------")
        # print(f"position : {type(position)}")
        # print(f"chromosome : {type(chromosome)}")
        # print(f"stride : {type(strides)}")

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

        if type(chromosome) == list:

            chromosome = chromosome[0]


        if trans_chromosome is not None:

            list_selected_chromosomes = list(chromosome.split()) + list(trans_chromosome)

        else:


            list_selected_chromosomes = list(chromosome.split())

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

        # print(f"trans_chromosome : {trans_chromosome}")

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

            # print(f"-- forward_intervals : {forward_intervals} --")


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

    
# TODO : Implement associated test

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


# TODO : Implement associated test
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
    
    
# TODO : implement associated test
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

# TODO : implement associated test
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

    if type(chromosome) == list:

        chromosome = chromosome[0]

    area_to_search = np.append(
        np.arange(2 * bins, chrom_sizes.get(chromosome) - 3 * bins, bins), # np.arange(0, cs_disctionary.item().get(chromosome) - bin_size, bin_size)
        chrom_sizes.get(chromosome),
    )

    # print(f"position : {position}")
    # print(f"bins : {bins}")

    if position % bins == 0:

        before_index = np.searchsorted(area_to_search, position, side="left")
        after_index = np.searchsorted(area_to_search, position, side="right")

    else:

        before_index = np.searchsorted(area_to_search, position, side="left") - 1
        after_index = np.searchsorted(area_to_search, position, side="right")

    lower_bound, upper_bound = area_to_search[before_index], area_to_search[after_index]

    return (lower_bound, upper_bound)

# TODO : implement test function
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

# TODO :  Implement associated test
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




