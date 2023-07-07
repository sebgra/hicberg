import re
import glob
import subprocess as sp
from os import getcwd
from os.path import join
from pathlib import Path, PurePath
import multiprocessing
from functools import partial
import itertools

from typing import Iterator

import numpy as np
from numpy.random import choice
import pandas as pd
import scipy.stats as stats
from scipy.stats import median_abs_deviation

import hicstuff.io as hio
from hicstuff.log import logger
import hicstuff.digest as hd

import pysam
from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, Analysis

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


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




def attribute_xs():
    pass

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
        Path to the folder wher eto save the dictionary, by default None

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

    if read.flag == 16:
        return True

    else:
        return False

def classify_reads(forward_bam_file : str = None, reverse_bam_file : str = None, chromosome_sizes : str = None, mapq : int = 30, output_dir : str = None) -> None:
    """
    Classification of pairs of reads in 2 different groups:
        Group 0) --> (Unmappable)
        Group 1) --> (Uniquely Mapped  Uniquely Mapped)
        Group 2) --> (Uniquely Mapped Multi Mapped) or (Multi Mapped  Multi Mapped).

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

    chromosome_sizes_dic = hio.load_dictionary(chromosome_sizes_path)

    ## TODO : to complete and check for default file names as 1.sorted.bam and 2.sorted.bam


    return 0

def classify_reads_multi():
    pass

def classify_reads_multi():
    pass

def is_intra_chromosome():
    pass

def get_ordered_reads():
    pass

def is_weird():
    pass

def is_uncut():
    pass

def is_circle():
    pass

def get_cis_distance():
    pass

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

def subsample_restriction_map():
    pass

def fill_zeros_with_last():
    pass

def max_consecutive_nans():
    pass

def mad_smoothing():
    pass
