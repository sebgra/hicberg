import re
import glob
import subprocess as sp
from os import getcwd
from os.path import join
from pathlib import Path, PurePath
import multiprocessing
from functools import partial
import itertools

from typing import Iterator, Tuple

import numpy as np
from numpy.random import choice
import pandas as pd
import scipy.stats as stats
from scipy.stats import median_abs_deviation

import hicstuff.io as hico
from hicstuff.log import logger
import hicstuff.digest as hd

import pysam
from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, Analysis

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import hicberg.io as hio


DEFAULT_FRAGMENTS_LIST_FILE_NAME = "fragments_list.txt"
DEFAULT_INFO_CONTIGS_FILE_NAME = "info_contigs.txt"
DEFAULT_SPARSE_MATRIX_FILE_NAME = "abs_fragments_contacts_weighted.txt"

DEFAULT_KB_BINNING = 1
DEFAULT_THRESHOLD_SIZE = 0
# Most used enzyme for eukaryotes
DEFAULT_ENZYME = "DpnII"
# If using evenly-sized chunks instead of restriction
# enzymes, they shouldn't be too short
DEFAULT_MIN_CHUNK_SIZE = 50


def get_num_lines_alignment_file():
    pass

def get_restriction_table():
    pass

def write_frag_info():
    pass

def get_chromosomes_sizes(genome : str = None, output_dir : str = None) -> None:
    """
    Generate a dictionnary save in .npy format where keys are chromosome name and value are size in bp.

    Parameters
    ----------
    genome : str, optional
        Path to the genome, by default None

    output_dir : str, optional
        Path to the folder where to save the dictionary, by default None
    """

    genome_path = Path(genome)

    if not genome_path.is_file():

        raise IOError(f"Genome file {genome_path.name} not found. Please provide a valid path.")

    if output_dir is None:

        folder_path = Path(getcwd())

    else:

        folder_path = Path(output_dir)

    chrom_sizes = {}

    output_file = folder_path / "chromosome_sizes.npy"

    for rec in SeqIO.parse(genome_path,"fasta"):

        chrom_sizes[rec.id] = len(rec.seq)

    np.save(output_file, chrom_sizes)


def get_bin_table(chrom_sizes_dict : str, bins : int, output_dir : str = None) -> None:
    """
    Create bin table containing start and end position for fixed size bin per chromosome.

    Parameters
    ----------
    chrom_sizes_dict : str
        Path to a dictionary containing chromosome sizes as {chromosome : size} saved in .npy format
    bins : int
        Size of the desired bin.
    output_dir : str, optional
        Path to the folder where to save the dictionary, by default None

    """    

    chrom_sizes_dict_path = Path(chrom_sizes_dict)

    if not chrom_sizes_dict_path.is_file():

        raise IOError(f"Genome file {chrom_sizes_dict_path.name} not found. Please provide a valid path.")

    if output_dir is None:

        folder_path = Path(getcwd())

    else:

        folder_path = Path(output_dir)

    output_file = folder_path / "fragments_fixed_sizes.txt"

    chrom_size_dic = np.load(chrom_sizes_dict_path, allow_pickle=True).item()
    chr_count = 0

    with open (output_file, "w") as f_out:
        
        for chrom, length in zip(chrom_size_dic.keys(), chrom_size_dic.values()):

            curr_chr, curr_length = chrom, length
            chr_count += 1

            if (curr_length % bins) == 0:
                    interval_end = curr_length
            else:
                interval_end = (int((curr_length + bins) / bins)) * bins

                for val in range(0, interval_end, bins):
                    curr_start = val

                    if val + bins > curr_length:

                        curr_end = curr_length
                    else:
                        curr_end = val + bins
                    if (chr_count > 1) or (val > 0):
                        f_out.write("\n")
                    f_out.write(
                        str(curr_chr)
                        + "\t"
                        + str(curr_start)
                        + "\t"
                        + str(int(curr_end))
                        + "\t"
                    )

        # close the output fragment file
        f_out.close()

def is_duplicated(read : pysam.AlignedSegment) -> bool:
    """
    Check if read from pysam AlignmentFile is mapping more than once along the genome.

    Parameters
    ----------
    read : pysam.AlignedSegment
        pysam AlignedSegment object.

    Returns
    -------
    bool
        True if the read is duplicated i.e. mapping to more than one position.
    """    

    if "XS" in [x[0] for x in read.get_tags()]:
        return True

    else:
        return False

def is_poor_quality(read : pysam.AlignedSegment, mapq : int) -> bool:
    """
    Check if read from pysam AlignmentFile is under mapping quality threshold

    Parameters
    ----------
    read : pysam.AlignedSegment
        pysam AlignedSegment object.
    mapq : int
        Mapping quality threshold.

    Returns
    -------
    bool
        True if the read quality is below mapq threshold.
    """    
    if 0 < read.mapping_quality < mapq:
        return True

    else:
        return False

def is_unqualitative(read : pysam.AlignedSegment) -> bool:
    """
    Check if the read is unqualitative.

    Parameters
    ----------
    read : pysam.AlignedSegment
        pysam AlignedSegment object.

    Returns
    -------
    bool
        True if the read is unqualitative, False otherwise.
    """

    if read.mapping_quality == 0:

        return True

    else:

        return False

def is_unmapped(read : pysam.AlignedSegment) -> bool:
    """
    Check if read from pysam AlignmentFile is unmapped

    Parameters
    ----------
    read : pysam.AlignedSegment
        pysam AlignedSegment object.

    Returns
    -------
    bool
        True if the read is unmapped.
    """    
    if read.flag == 4:
        return True

    else:
        return False

def is_reverse(read : pysam.AlignedSegment) -> bool:
    """
    Check if read from pysam AlignmentFile is reverse

    Parameters
    ----------
    read : pysam.AlignedSegment
        pysam AlignedSegment object.

    Returns
    -------
    bool
        True if the read is reverse.
    """ 

    if read.flag == 16 or read.flag == 272:
        return True

    else:
        return False

def classify_reads(forward_bam_file : str = None, reverse_bam_file : str = None, chromosome_sizes : str = None, mapq : int = 30, output_dir : str = None) -> None:
    """
    Classification of pairs of reads in 2 different groups:
        Group 0) --> (Unmappable) - files :group0.1.bam and group0.2.bam
        Group 1) --> (Uniquely Mapped  Uniquely Mapped) - files :group1.1.bam and group1.2.bam
        Group 2) --> (Uniquely Mapped Multi Mapped) or (Multi Mapped  Multi Mapped).- files :group2.1.bam and group2.2.bam

    Parameters
    ----------
    forward_bam_file : str, optional
        Path to forward .bam alignment file, by default None
    reverse_bam_file : str, optional
        Path to reverse .bam alignment file, by default None
    chromosome_sizes : str, optional
        Path to a chromosome size dictionary save in .npy format, by default None
    mapq : int, optional
        Minimal read quality under which a Hi-C read pair will not be kept, by default 30
    output_dir : str, optional
        Path to the folder where to save the classified alignment files, by default None

    """

    forward_bam_file_path, reverse_bam_file_path = Path(forward_bam_file), Path(reverse_bam_file)

    chromosome_sizes_path = Path(chromosome_sizes)

    if not forward_bam_file_path.is_file():

        

        raise IOError(f"Forward alignment file {forward_bam_file_path.name} not found. Please provide a valid path.")

    if not reverse_bam_file_path.is_file():
            
            raise IOError(f"Reverse alignment file {reverse_bam_file_path.name} not found. Please provide a valid path.")
    
    if not chromosome_sizes_path.is_file():
            
            raise IOError(f"Chromosome sizes file {chromosome_sizes_path.name} not found. Please provide a valid path.")
    
    if output_dir is None:
            
            output_dir = Path(getcwd())
    else:

        output_dir = Path(output_dir)

    chromosome_sizes_dic = hio.load_dictionary(chromosome_sizes_path)

    #opening files to parse
    forward_bam_file = pysam.AlignmentFile(forward_bam_file_path, "rb")
    reverse_bam_file = pysam.AlignmentFile(reverse_bam_file_path, "rb")


    # create iterators

    forward_bam_file_iter = bam_iterator(forward_bam_file_path)
    reverse_bam_file_iter = bam_iterator(reverse_bam_file_path)

    # open the output files handlers
    unmapped_bam_file_foward = pysam.AlignmentFile(output_dir / "group0.1.bam", "wb", template = forward_bam_file)
    unmapped_bam_file_reverse = pysam.AlignmentFile(output_dir / "group0.2.bam", "wb", template = reverse_bam_file)

    uniquely_mapped_bam_file_foward = pysam.AlignmentFile(output_dir / "group1.1.bam", "wb", template = forward_bam_file)
    uniquely_mapped_bam_file_reverse = pysam.AlignmentFile(output_dir / "group1.2.bam", "wb", template = reverse_bam_file)

    multi_mapped_bam_file_foward = pysam.AlignmentFile(output_dir / "group2.1.bam", "wb", template = forward_bam_file)
    multi_mapped_bam_file_reverse = pysam.AlignmentFile(output_dir / "group2.2.bam", "wb", template = reverse_bam_file)

    for forward_block, reverse_block in zip(forward_bam_file_iter, reverse_bam_file_iter):
        


        unmapped_couple, multi_mapped_couple = False, False

        forward_reverse_combinations = list(itertools.product(tuple(forward_block), tuple(reverse_block)))

        for combination in forward_reverse_combinations:

            if is_unqualitative(combination[0]) or is_unmapped(combination[0]) or is_unqualitative(combination[1]) or is_unmapped(combination[1]):

                unmapped_couple = True

                break

            elif is_duplicated(combination[0]) or is_poor_quality(combination[0], mapq) or is_duplicated(combination[1]) or is_poor_quality(combination[1], mapq):

                multi_mapped_couple = True

                break

            for forward_read in forward_block:

                if unmapped_couple :

                    unmapped_bam_file_foward.write(forward_read)

                elif multi_mapped_couple:

                    forward_read.set_tag("XG", chromosome_sizes_dic[forward_read.reference_name])
                    multi_mapped_bam_file_foward.write(forward_read)

                else : 

                    forward_read.set_tag("XG", chromosome_sizes_dic[forward_read.reference_name])
                    uniquely_mapped_bam_file_foward.write(forward_read)

            for reverse_read in reverse_block:

                if unmapped_couple :

                    unmapped_bam_file_reverse.write(reverse_read)

                elif multi_mapped_couple:

                    reverse_read.set_tag("XG", chromosome_sizes_dic[reverse_read.reference_name])
                    multi_mapped_bam_file_reverse.write(reverse_read)

                else : 

                    reverse_read.set_tag("XG", chromosome_sizes_dic[reverse_read.reference_name])
                    uniquely_mapped_bam_file_reverse.write(reverse_read)

    #closing files
    forward_bam_file.close()
    reverse_bam_file.close()
    unmapped_bam_file_foward.close()
    unmapped_bam_file_reverse.close()
    uniquely_mapped_bam_file_foward.close()
    uniquely_mapped_bam_file_reverse.close()
    multi_mapped_bam_file_foward.close()
    multi_mapped_bam_file_reverse.close()

    print(f"Files for the different groups have been saved in {output_dir}")


def classify_reads_multi():
    pass

def classify_reads_multi():
    pass

def is_intra_chromosome(read_forward : pysam.AlignedSegment, read_reverse : pysam.AlignedSegment) -> bool:
    """
    Return True if two reads of a pair came from the same chromosome.

    Parameters
    ----------
    read_forward : pysam.AlignedSegment
        Forward read of the pair.
    read_reverse : pysam.AlignedSegment
        Reverse read of the pair.

    Returns
    -------
    bool
        True if the pair is intrachromosomic, False otherwise.
    """    

    if read_forward.query_name != read_reverse.query_name:
        raise ValueError("Reads are not comming from the same pair")

    if read_forward.reference_name == read_reverse.reference_name:
        return True
    else:
        return False 

def get_ordered_reads(read_forward : pysam.AlignedSegment, read_reverse : pysam.AlignedSegment) -> Tuple[pysam.AlignedSegment, pysam.AlignedSegment]:
    """
    Returns the ordered pair of reads in the same chromosome as the two reads .

    Parameters
    ----------
    read_forward : pysam.AlignedSegment
        Forward read to compare with the reverse read.
    read_reverse : pysam.AlignedSegment
        Reverse read to compare with the forward read.

    Returns
    -------
    Tuple[pysam.AlignedSegment, pysam.AlignedSegment]
        The ordered pair of reads in the same chromosome as the two reads.
    """    
    if read_forward.reference_name != read_reverse.reference_name:
            
        raise ValueError("The two reads must be mapped on the same chromosome.")
        
    if is_reverse(read_forward):

        forward_start = read_forward.reference_end
    
    elif not is_reverse(read_forward):
                
        forward_start = read_forward.reference_start

    if is_reverse(read_reverse):
            
        reverse_start = read_reverse.reference_end

    elif not is_reverse(read_reverse):
    
        reverse_start = read_reverse.reference_start

    
    if forward_start < reverse_start:
        
        return (read_forward, read_reverse)
    
    elif forward_start > reverse_start:

        return (read_reverse, read_forward)

def is_weird(read_forward : pysam.AlignedSegment, read_reverse : pysam.AlignedSegment) -> bool:
    """
    Check if two reads are forming a weird pattern .

    Parameters
    ----------
    read_forward : pysam.AlignedSegment
        Forward read of the pair
    read_reverse : pysam.AlignedSegment
        Reverse read of the pair

    Returns
    -------
    bool
        True if the two reads are forming a weird pattern, False otherwise.
    """
    
    if read_forward.query_name != read_reverse.query_name:

        raise ValueError("The two reads must be mapped on the same chromosome.")
    
    read_forward, read_reverse = get_ordered_reads(read_forward, read_reverse)

    if (
        (read_forward.flag == read_reverse.flag == 0)
        or (read_forward.flag == read_reverse.flag == 16)
        or (read_forward.flag == read_reverse.flag == 272)
        or (read_forward.flag == read_reverse.flag == 256)
        or (read_forward.flag == 256 and read_reverse.flag == 0)
        or (read_forward.flag == 0 and read_reverse.flag == 256)
        or (read_forward.flag == 16 and read_reverse.flag == 272)
        or (read_forward.flag == 272 and read_reverse.flag == 16)
    ):
        return True
    
    else:
        return False

def is_uncut(read_forward : pysam.AlignedSegment, read_reverse : pysam.AlignedSegment) -> bool:
    """
    Check if two reads are forming an uncut pattern .

    Parameters
    ----------
    read_forward : pysam.AlignedSegment
        Forward read of the pair
    read_reverse : pysam.AlignedSegment
        Reverse read of the pair

    Returns
    -------
    bool
        True if the two reads are forming an uncut pattern, False otherwise.
    """
        
    if read_forward.query_name != read_reverse.query_name:

        raise ValueError("The two reads must be mapped on the same chromosome.")
    
    read_forward, read_reverse = get_ordered_reads(read_forward, read_reverse)

    if (
        (read_forward.flag == 0 and read_reverse.flag == 16)
        or (read_forward.flag == 256 and read_reverse.flag == 16)
        or (read_forward.flag == 0 and read_reverse.flag == 272)
        or (read_forward.flag == 256 and read_reverse.flag == 272)
    ):
        return True
    else:
        return False

    

def is_circle(read_forward : pysam.AlignedSegment, read_reverse : pysam.AlignedSegment) -> bool:
    """
    Check if two reads are forming a loop pattern .

    Parameters
    ----------
    read_forward : pysam.AlignedSegment
        Forward read of the pair
    read_reverse : pysam.AlignedSegment
        Reverse read of the pair

    Returns
    -------
    bool
        True if the two reads are forming a loop pattern, False otherwise.
    """

    if read_forward.query_name != read_reverse.query_name:
        raise ValueError("Reads are not comming from the same pair")

    read_forward, read_reverse = get_ordered_reads(read_forward, read_reverse)

    if (
        (read_forward.flag == 16 and read_reverse.flag == 0)
        or (read_forward.flag == 272 and read_reverse.flag == 0)
        or (read_forward.flag == 16 and read_reverse.flag == 256)
        or (read_forward.flag == 272 and read_reverse.flag == 256)
    ):
        return True
    else:
        return False

def get_cis_distance(read_forward : pysam.AlignedSegment, read_reverse : pysam.AlignedSegment, circular : str = "") -> int:
    """
    Calculate the distance between two reads in the same pairwise alignment .

    Parameters
    ----------
    read_forward : pysam.aligned_segment
        Forward read of the pair
    read_reverse : pysam.AlignedSegment
        Reverse read of the pair
    circular : str, optional
        Name of the chromosomes to consider as circular, by default None, by default "".

    Returns
    -------
    int
        Genomic distance separating the two reads (bp).

    """    
    if read_forward.query_name != read_reverse.query_name:
        raise ValueError("Reads are not comming from the same pair")

    if is_intra_chromosome(read_forward, read_reverse):

        read_forward, read_reverse = get_ordered_reads(read_forward, read_reverse)

        if is_weird(read_forward, read_reverse):
            distance = np.abs(
                np.subtract(read_forward.reference_start, read_reverse.reference_start)
            )

        elif is_uncut(read_forward, read_reverse):
            distance = np.abs(np.subtract(read_forward.reference_start, read_reverse.reference_end))

        elif is_circle(read_forward, read_reverse):
            distance = np.abs(np.subtract(read_forward.reference_end, read_reverse.reference_start))

        # circular mode
        if read_forward.reference_name in circular:

            clockwise_distance = distance
            anti_clockwise_distance = np.subtract(read_forward.get_tag("XG"), distance)
            return np.min([clockwise_distance, anti_clockwise_distance])

        # linear mode
        else:
            return distance


def get_event_stats():
    pass

def bam_iterator(bam_file : str = None) -> Iterator[pysam.AlignedSegment]:
    """
    Returns an iterator for the given SAM/BAM file (must be query-sorted).
    In each call, the alignments of a single read are yielded as a 3-tuple: (list of primary pysam.AlignedSegment, list of supplementary pysam.AlignedSegment, list of secondary pysam.AlignedSegment).

    Parameters
    ----------
    bam : [str]
        Path to alignment file in .sam or .bam format.

    Yields
    -------
    Iterator[pysam.AlignedSegment]
        Yields a list containing pysam AlignementSegment objects, within which all the reads have the same id.
    """

    bam_path = Path(bam_file)

    if not bam_path.is_file():

        raise IOError(f"BAM file {bam_path.name} not found. Please provide a valid path.")

    bam_handler = pysam.AlignmentFile(bam_path, "rb")
    
    alignments = bam_handler.fetch(until_eof=True)
    current_aln = next(alignments)
    current_read_name = current_aln.query_name

    block = []
    block.append(current_aln)

    while True:
        try:
            next_aln = next(alignments)
            next_read_name = next_aln.query_name
            if next_read_name != current_read_name:
                yield (block)
                current_read_name = next_read_name
                block = []
                block.append(next_aln)

            else:
                block.append(next_aln)
        except StopIteration:
            break

    yield (block)

def block_counter():
    pass

def chunk_bam():
    pass

def get_pair_cover():
    pass

def get_pair_ps():
    pass

def get_trans_ps():
    pass

def get_d1d2():
    pass

def get_d1d2_distance():
    pass

def compute_propentsity():
    pass

def decompose_propentsity():
    pass

def check_propensity():
    pass

def draw_read_couple():
    pass

def reattribute_reads():
    pass

def reattribute_reads_multiprocess():
    pass

def inspect_reads():
    pass

def get_bin_rsites():
    pass

def get_full_bin_rsites():
    pass

def subsample_restriction_map(restriction_map : dict = None, rate : float = 1.0) -> dict[str, np.ndarray[int]]:
    """
    Subsample a restriction map by a given rate.

    Parameters
    ----------
    restriction_map : dict, optional
        Restriction map saved as a dictionary like chrom_name : list of restriction sites' position, by default None
    rate : float, optional
        Set the proportion of restriction sites to condiser. Avoid memory overflow when restriction maps are very dense, by default 1.0

    Returns
    -------
    dict[str, np.ndarray[int]]
        Dictionary of subsampled restriction map with keys as chromosome names and values as lists of restriction sites' position.

    """

    if (0.0 > rate) or (rate > 1.0):
        raise ValueError("Subsampling rate must be between 0.0 and 1.0.")
    

    subsampled_restriction_map = {}

    for chromosome in restriction_map:
            
        if int(len(restriction_map.get(str(chromosome))) * rate) < 5:

            subsampled_restriction_map[str(chromosome)] = restriction_map[str(chromosome)]

            continue

        size_sample = int(len(restriction_map.get(str(chromosome))) * rate)

        subsampled_restriction_map[str(chromosome)] = np.random.choice(
            restriction_map.get(str(chromosome)), size_sample, replace=False
        )
        
        subsampled_restriction_map[str(chromosome)] = np.sort(subsampled_restriction_map[str(chromosome)])

        if subsampled_restriction_map[str(chromosome)][0] != 0:
            subsampled_restriction_map[str(chromosome)][0] = 0

        if (
            subsampled_restriction_map[str(chromosome)][-1]
            != restriction_map.get(str(chromosome))[-1]
        ):
            subsampled_restriction_map[str(chromosome)][-1] = restriction_map.get(str(chromosome))[-1]        

    return subsampled_restriction_map

def fill_zeros_with_last():
    pass

def max_consecutive_nans(vector : np.ndarray) -> int:
    """
    Return the maximum number of consecutiver NaN values in a vecotr.

    Parameters
    ----------
    vector : np.ndarray
        Vector to get the maximum number of consecutive NaN values from.

    Returns
    -------
    int
        Number of maximum consecutive NaN values.
    """

    mask = np.concatenate(([False], np.isnan(vector), [False]))
    if ~mask.any():
        return 0
    else:
        idx = np.nonzero(mask[1:] != mask[:-1])[0]
        return (idx[1::2] - idx[::2]).max()

def mad_smoothing(vector : np.ndarray[int] = None, window_size : int | str = "auto", nmads :int = 1) -> np.ndarray[int]:
    """
    Apply MAD smoothing to an vector .

    Parameters
    ----------
    vector : np.ndarray[int], optional
        Data to smooth, by default None
    window_size : int or str, optional
        Size of the window to perform mean sliding average in. Window is center on current value as [current_value - window_size/2] U [current_value + window_size/2], by default "auto"
    nmads : int, optional
        number of median absolute deviation tu use, by default 1

    Returns
    -------
    np.ndarray[int]
        MAD smoothed vector.
    """

    mad = median_abs_deviation(vector)
    threshold = np.median(vector) - nmads * mad
    imputed_nan_data = np.where(vector < threshold, np.nan, vector)

    if window_size == "auto":
        # due to centered window, selected windows for rolling mean is :
        # [window_size / 2 <-- center_value --> window_size / 2]
        window_size = (max_consecutive_nans(imputed_nan_data) * 2) + 1

    averaged_data = (
        pd.Series(imputed_nan_data)
        .rolling(window=window_size, min_periods=1, center=True)
        .apply(lambda x: np.nanmean(x))
        .to_numpy()
    )

    return averaged_data
