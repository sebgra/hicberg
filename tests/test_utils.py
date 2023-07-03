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




def test_attribute_xs():
    pass

def test_get_num_lines_alignment_file():
    pass

def test_get_restriction_table():
    pass

def test_write_frag_info():
    pass

def test_get_chromosomes_sizes():
    pass

def test_get_bin_table():
    pass

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

def test_sam_iterator():
    pass

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
