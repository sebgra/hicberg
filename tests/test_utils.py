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

import hicstuff.io as hio
from hicstuff.log import logger
import hicstuff.digest as hd

import pysam
from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, Analysis

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import hicberg.utils as hut

GENOME = "data_test/SC288_with_micron.fa"
CHROM_SIZES_DIC = "chromosome_sizes.npy"
FRAG_FILE = "fragments_fixed_sizes.txt"

FORWARD_SORTED_BAM = "data_test/alignments/1.sorted.bam"


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

@pytest.fixture(scope="session")
def temporary_folder():
    fn = tempfile.TemporaryDirectory()
    # yiled temporary folder name to be used as Path object
    yield fn.name


def test_attribute_xs():
    pass

def test_get_num_lines_alignment_file():
    pass

def test_get_restriction_table():
    pass

def test_write_frag_info():
    pass

@pytest.fixture(name = "get_chromosomes_sizes")
def test_get_chromosomes_sizes(temporary_folder):

    temp_dir_path = Path(temporary_folder)
    hut.get_chromosomes_sizes(genome = GENOME, output_dir = temp_dir_path)

    chrom_sizes_dictionary_path = temp_dir_path / CHROM_SIZES_DIC

    yield temp_dir_path / CHROM_SIZES_DIC

    # Check if the sorted alignement files are created
    assert chrom_sizes_dictionary_path.is_file()

def test_get_bin_table(temporary_folder, get_chromosomes_sizes):
    
    temp_dir_path = Path(temporary_folder)
    hut.get_bin_table(chrom_sizes_dict = get_chromosomes_sizes, bins = 2000,  output_dir = temp_dir_path)

    bin_table_path = temp_dir_path / FRAG_FILE

    assert bin_table_path.is_file()

def test_is_duplicated():
    pass

def test_is_poor_quality():
    pass

def test_is_unqualitative():
    pass

def test_is_unmapped():
    pass

def test_is_reverse():
    pass

def test_classify_reads():
    pass

def test_classify_reads_multi():
    pass

def test_classify_reads_multi():
    pass

def test_is_intra_chromosome():
    pass

def test_get_ordered_reads():
    pass

def test_is_weird():
    pass

def test_is_uncut():
    pass

def test_is_circle():
    pass

def test_get_cis_distance():
    pass

def test_get_event_stats():
    pass

def test_bam_iterator(bam_file = FORWARD_SORTED_BAM):

    bam_path = Path(bam_file)
    bam_iterator = hut.bam_iterator(bam_path)
    first_iteration = next(bam_iterator)

    # Check if the first iteration return a list of 1
    assert(len(first_iteration) == 1)
    # Check if first iteration has the right query name
    assert first_iteration[0].query_name  == "NS500150:487:HNLLNBGXC:1:11101:1047:14348"

def test_block_counter():
    pass

def test_chunk_bam():
    pass

def test_get_pair_cover():
    pass

def test_get_pair_ps():
    pass

def test_get_trans_ps():
    pass

def test_get_d1d2():
    pass

def test_get_d1d2_distance():
    pass

def test_compute_propentsity():
    pass

def test_decompose_propentsity():
    pass

def test_check_propensity():
    pass

def test_draw_read_couple():
    pass

def test_reattribute_reads():
    pass

def test_reattribute_reads_multiprocess():
    pass

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
    pass

def test_mad_smoothing():
    pass
