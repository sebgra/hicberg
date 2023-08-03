import pytest
import tempfile
import re
import glob
import subprocess as sp
from os.path import join
from pathlib import Path, PurePath
import multiprocessing
from functools import partial
import itertools

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

import hicberg.utils as hut
import hicberg.io as hio

from .conftest import temporary_folder

from .test_align import test_hic_build_index, test_hic_align, test_hic_view, test_hic_sort

GENOME = "data_test/SC288_with_micron.fa"
CHROM_SIZES_DIC = "chromosome_sizes.npy"
FRAG_FILE = "fragments_fixed_sizes.txt"

FORWARD_SORTED_BAM = "data_test/alignments/1.sorted.bam"

FORWARD_UNMAPPED_BAM = "group0.1.bam"
REVERSE_UNMAPPED_BAM = "group0.2.bam"
FORWARD_UNIQUE_BAM = "group1.1.bam"
REVERSE_UNIQUE_BAM = "group1.2.bam"
FORWARD_MULTI_BAM = "group2.1.bam"
REVERSE_MULTI_BAM = "group2.2.bam"



MIN_READ_MAPQ = 30
DICT_FIRST_KEY = "chr10"
DICT_FIRST_CHR_LAST_POS = 745751

FIRST_QUERY_NAME = "NS500150:487:HNLLNBGXC:1:11101:1047:14348"

HEADER = pysam.AlignmentHeader().from_dict({
            "HD": {"VN": "1.0", "SO": "unsorted"},
            "SQ": [
                {"SN": DICT_FIRST_KEY, "LN": DICT_FIRST_CHR_LAST_POS},
            ],
        })


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

BINS = 2000
MODE = "full"




def test_get_num_lines_alignment_file():
    pass

def test_get_restriction_table():
    pass

def test_write_frag_info():
    pass

@pytest.fixture(scope = "session")
def test_get_chromosomes_sizes(temporary_folder):
    """
    Test if the function creating a dictionary of chromosomes sizes from a genome is correctly performed.
    """

    temp_dir_path = Path(temporary_folder)
    hut.get_chromosomes_sizes(genome = GENOME, output_dir = temp_dir_path)

    chrom_sizes_dictionary_path = temp_dir_path / CHROM_SIZES_DIC

    yield temp_dir_path / CHROM_SIZES_DIC

    # Check if the sorted alignement files are created
    assert chrom_sizes_dictionary_path.is_file()

@pytest.fixture(scope = "session")
def test_get_bin_table(temporary_folder, test_get_chromosomes_sizes):
    """
    Test if the function creating a bin table from a genome is correctly performed.
    """
    temp_dir_path = Path(temporary_folder)
    hut.get_bin_table(chrom_sizes_dict = test_get_chromosomes_sizes, bins = 2000,  output_dir = temp_dir_path)

    bin_table_path = temp_dir_path / FRAG_FILE

    assert bin_table_path.is_file()

def test_is_duplicated():
    """
    Test if a read is duplicated
    """
    read_duplicated = pysam.AlignedSegment(header = HEADER)
    read_duplicated.query_name = "DUPLICATED"
    read_duplicated.query_sequence = "ATCG"
    read_duplicated.set_tag("XS", -1)

    duplicateness = hut.is_duplicated(read_duplicated)

    assert duplicateness 

def test_is_poor_quality():
    """
    Test if a read is poor quality.
    """
    read_poor_quality = pysam.AlignedSegment(header = HEADER)
    read_poor_quality.query_name = "POOR QUALITY"
    read_poor_quality.query_sequence = "ATCG"
    read_poor_quality.mapping_quality = 3

    qualitiveness = hut.is_poor_quality(read_poor_quality, MIN_READ_MAPQ)

    assert qualitiveness

def test_is_unqualitative():
    """
    Test if a read is unqualitative.
    """

    read_unqualitative = pysam.AlignedSegment(header = HEADER)
    read_unqualitative.query_name = "UNQUALITATIVE"
    read_unqualitative.query_sequence = "ATCG"
    read_unqualitative.mapping_quality = 0
    read_unqualitative.flag = 0

    qualitativeness = hut.is_unqualitative(read_unqualitative)

    assert qualitativeness 

def test_is_unmapped():
    """
    Test if a read is unmapped.
    """

    read_unmapped = pysam.AlignedSegment(header = HEADER)
    read_unmapped.query_name = "UNMAPPED"
    read_unmapped.query_sequence = "ATCG"
    read_unmapped.mapping_quality = 0
    read_unmapped.flag = 4
    mapping = hut.is_unmapped(read_unmapped)

    assert mapping


def test_is_reverse():
    """
    Test if a read is reverse.
    """

    read_reverse = pysam.AlignedSegment(header = HEADER)
    read_reverse.query_name = "REVERSE"
    read_reverse.query_sequence = "ATCG"
    read_reverse.mapping_quality = 0
    read_reverse.flag = 16
    reverse = hut.is_reverse(read_reverse)

    assert reverse

@pytest.fixture(scope = "session")
def test_classify_reads(temporary_folder, test_hic_sort, test_get_chromosomes_sizes):
    """
    Test if the creation of the different bam files per group is correctly performed.
    """

    temp_dir_path = Path(temporary_folder)
    hut.classify_reads(forward_bam_file = test_hic_sort[0], reverse_bam_file = test_hic_sort[1] , chromosome_sizes = test_get_chromosomes_sizes, mapq = MIN_READ_MAPQ,  output_dir = temp_dir_path)
    
    forward_unmapped_path = temp_dir_path / FORWARD_UNMAPPED_BAM
    reverse_unmapped_path = temp_dir_path / REVERSE_UNMAPPED_BAM
    forward_unique_path = temp_dir_path / FORWARD_UNIQUE_BAM
    reverse_unique_path = temp_dir_path / REVERSE_UNIQUE_BAM
    forward_multi_path = temp_dir_path / FORWARD_MULTI_BAM
    reverse_multi_path = temp_dir_path / REVERSE_MULTI_BAM

    yield forward_unique_path, reverse_unique_path, forward_multi_path, reverse_multi_path

    assert forward_unmapped_path.is_file()
    assert reverse_unmapped_path.is_file()
    assert forward_unique_path.is_file()
    assert reverse_unique_path.is_file()
    assert forward_multi_path.is_file()
    assert reverse_multi_path.is_file()

def test_classify_reads_multi():
    pass

def test_is_intra_chromosome():
    
    read_forward = pysam.AlignedSegment(header = HEADER)
    read_forward.query_name = "TEST"
    read_forward.reference_name = DICT_FIRST_KEY
    read_forward.reference_start = 500
    read_forward.cigarstring = "100M"
    read_forward.flag = 0

    read_reverse = pysam.AlignedSegment(header = HEADER)
    read_reverse.query_name = "TEST"
    read_reverse.reference_name = DICT_FIRST_KEY
    read_reverse.reference_start = 100
    read_reverse.cigarstring = "100M"
    read_reverse.flag = 16

    intra_chomosomic = hut.is_intra_chromosome(read_forward, read_reverse)

    assert intra_chomosomic

def test_get_ordered_reads():
    """
    Test if the function ordering the reads in a bam file is ordered reads correctly.
    """
    read_forward = pysam.AlignedSegment(header = HEADER)
    read_forward.query_name = "QUERY_1"
    read_forward.reference_name = DICT_FIRST_KEY
    read_forward.reference_start = 500
    read_forward.cigarstring = "100M"
    read_forward.flag = 0

    read_reverse = pysam.AlignedSegment(header = HEADER)
    read_reverse.query_name = "QUERY_1"
    read_reverse.reference_name = DICT_FIRST_KEY
    read_reverse.reference_start = 100
    read_reverse.cigarstring = "100M"
    read_reverse.flag = 16

    read_forward_2 = pysam.AlignedSegment(header = HEADER)
    read_forward_2.query_name = "QUERY_2"
    read_forward_2.reference_name = DICT_FIRST_KEY
    read_forward_2.reference_start = 500
    read_forward_2.cigarstring = "100M"
    read_forward_2.flag = 272

    read_reverse_2 = pysam.AlignedSegment(header = HEADER)
    read_reverse_2.query_name = "QUERY_2"
    read_reverse_2.reference_name = DICT_FIRST_KEY
    read_reverse_2.reference_start = 100
    read_reverse_2.cigarstring = "100M"
    read_reverse_2.flag = 0

    read_forward_3 = pysam.AlignedSegment(header = HEADER)
    read_forward_3.query_name = "QUERY_3"
    read_forward_3.reference_name = DICT_FIRST_KEY
    read_forward_3.reference_start = 180
    read_forward_3.cigarstring = "100M"
    read_forward_3.flag = 16

    read_reverse_3 = pysam.AlignedSegment(header = HEADER)
    read_reverse_3.query_name = "QUERY_3"
    read_reverse_3.reference_name = DICT_FIRST_KEY
    read_reverse_3.reference_start = 300
    read_reverse_3.cigarstring = "100M"
    read_reverse_3.flag = 0

    # Check case of direct and indirect reads that have to be reordered
    ordered_read_forward, ordered_read_reverse = hut.get_ordered_reads(read_forward, read_reverse)

    assert ordered_read_forward == read_reverse and ordered_read_reverse == read_forward

    # Check case of direct and indirect reads that have to be reordered
    ordered_read_forward, ordered_read_reverse = hut.get_ordered_reads(read_forward_2, read_reverse_2)

    assert ordered_read_forward == read_reverse_2 and ordered_read_reverse == read_forward_2

    # Check case of direct and indirect reads that do not have to be reordered
    ordered_read_forward, ordered_read_reverse = hut.get_ordered_reads(read_forward_3, read_reverse_3)

    
def test_is_weird():
    """
    Test if two reads are forming weird pattern.
    """

    read_forward = pysam.AlignedSegment(header = HEADER)

    read_forward.query_name = "TEST1"
    read_forward.reference_name = DICT_FIRST_KEY
    read_forward.reference_start = 500
    read_forward.cigarstring = "100M"
    read_forward.flag = 0

    read_reverse = pysam.AlignedSegment(header = HEADER)
    read_reverse.query_name = "TEST1"
    read_reverse.reference_name = DICT_FIRST_KEY
    read_reverse.reference_start = 100
    read_reverse.cigarstring = "100M"
    read_reverse.flag = 0

    read_forward_2 = pysam.AlignedSegment(header = HEADER)
    read_forward_2.query_name = "TEST2"
    read_forward_2.reference_name = DICT_FIRST_KEY
    read_forward_2.reference_start = 500
    read_forward_2.cigarstring = "100M"
    read_forward_2.flag = 16

    read_reverse_2 = pysam.AlignedSegment(header = HEADER)
    read_reverse_2.query_name = "TEST2"
    read_reverse_2.reference_name = DICT_FIRST_KEY
    read_reverse_2.reference_start = 100
    read_reverse_2.cigarstring = "100M"
    read_reverse_2.flag = 16

    read_forward_3 = pysam.AlignedSegment(header = HEADER)
    read_forward_3.query_name = "TEST3"
    read_forward_3.reference_name = DICT_FIRST_KEY
    read_forward_3.reference_start = 500
    read_forward_3.cigarstring = "100M"
    read_forward_3.flag = 272

    read_reverse_3 = pysam.AlignedSegment(header = HEADER)
    read_reverse_3.query_name = "TEST3"
    read_reverse_3.reference_name = DICT_FIRST_KEY
    read_reverse_3.reference_start = 100
    read_reverse_3.cigarstring = "100M"
    read_reverse_3.flag = 256

    weird = hut.is_weird(read_forward, read_reverse)

    assert weird 

    weird = hut.is_weird(read_forward_2, read_reverse_2)

    assert weird

    not_weird = hut.is_weird(read_forward_3, read_reverse_3)

    assert not not_weird


def test_is_uncut():
    """
    Test if two reads are forming uncut pattern.
    """
    
    read_forward = pysam.AlignedSegment(header = HEADER)
    read_forward.query_name = "TEST"
    read_forward.reference_name = DICT_FIRST_KEY
    read_forward.reference_start = 100
    read_forward.cigarstring = "100M"
    read_forward.flag = 0

    read_reverse = pysam.AlignedSegment(header = HEADER)
    read_reverse.query_name = "TEST"
    read_reverse.reference_name = DICT_FIRST_KEY
    read_reverse.reference_start = 500
    read_reverse.cigarstring = "100M"
    read_reverse.flag = 16

    uncut = hut.is_uncut(read_forward, read_reverse)

    assert uncut


def test_is_circle():
    """
    Test if two reads are forming circle pattern.
    """
    
    read_forward = pysam.AlignedSegment(header = HEADER)
    read_forward.query_name = "TEST"
    read_forward.reference_name = DICT_FIRST_KEY
    read_forward.reference_start = 100
    read_forward.cigarstring = "100M"
    read_forward.flag = 16

    read_reverse = pysam.AlignedSegment(header = HEADER)
    read_reverse.query_name = "TEST"
    read_reverse.reference_name = DICT_FIRST_KEY
    read_reverse.reference_start = 500
    read_reverse.cigarstring = "100M"
    read_reverse.flag = 0


    circle = hut.is_circle(read_forward, read_reverse)

    assert circle


def test_get_cis_distance():
    """
    Test if the distance between a intrachromosomal read pair is correctly computed.
    """
    read_forward = pysam.AlignedSegment(header = HEADER)
    read_forward.query_name = "TEST"
    read_forward.reference_name = DICT_FIRST_KEY
    read_forward.reference_start = 100
    read_forward.cigarstring = "100M"
    read_forward.flag = 0

    read_reverse = pysam.AlignedSegment(header = HEADER)
    read_reverse.query_name = "TEST"
    read_reverse.reference_name = DICT_FIRST_KEY
    read_reverse.reference_start = 500
    read_reverse.cigarstring = "100M"
    read_reverse.flag = 0

    distance = hut.get_cis_distance(read_forward, read_reverse)

    assert distance == 400

def test_get_event_stats():
    pass

def test_bam_iterator(bam_file = FORWARD_SORTED_BAM):
    """
    Test if the bam_iterator function is correctly built by returning the first iteration.
    """


    bam_path = Path(bam_file)
    bam_iterator = hut.bam_iterator(bam_path)
    first_iteration = next(bam_iterator)

    # Check if the first iteration return a list of 1
    assert(len(first_iteration) == 1)
    # Check if first iteration has the right query name
    assert first_iteration[0].query_name  == FIRST_QUERY_NAME

def test_block_counter(test_classify_reads):
    """
    Test if the block counter function is correctly counting the number of blocks in a bam file.
    """


    forward_bam_file, reverse_bam_file = str(test_classify_reads[2]), str(test_classify_reads[3])

    nb_forward_block, nb_reverse_block = hut.block_counter(forward_bam_file = forward_bam_file, reverse_bam_file = reverse_bam_file)
    
    assert nb_forward_block != 0 and nb_reverse_block != 0

@pytest.fixture(scope = "module")
def test_chunk_bam(temporary_folder, test_classify_reads):
    """
    Test if the chunk bam function is correctly chunking a couple of bam file.
    """

    temp_dir_path = Path(temporary_folder)
    forward_bam_file, reverse_bam_file = str(test_classify_reads[2]), str(test_classify_reads[3])


    hut.chunk_bam(forward_bam_file = forward_bam_file, reverse_bam_file = reverse_bam_file, nb_chunks = 12, output_dir = temp_dir_path)

    chunks_path = temp_dir_path / 'chunks'
    
    is_full =  any(chunks_path.iterdir())


    yield chunks_path

    assert is_full

def test_inspect_reads():
    pass

def test_get_bin_rsites():
    pass

def test_get_full_bin_rsites():
    pass

def test_subsample_restriction_map():
    pass

def test_fill_zeros_with_last():
    pass

def test_max_consecutive_nans():
    """
    Test if the function correctly return the maximum number of consecutive NaNs in a vector.
    """

    vector_test = np.array([np.nan, 2, 3, 4, np.nan, np.nan, 7, 8, 9])

    max_consecutive_nans = hut.max_consecutive_nans(vector_test)

    assert max_consecutive_nans == 2

# def test_mad_smoothing():
#     """
#     Test if the function correctly smooth a vector using MAD smoothing.
#     """

#     # TODO :  to implement
#     assert True


def test_get_chunks(test_chunk_bam):

    ret = test_chunk_bam

    assert True


