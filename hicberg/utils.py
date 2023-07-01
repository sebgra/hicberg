import re
import glob
import subprocess as sp
from os.path import join
from pathlib import Path, PurePath
import multiprocessing
from tools import partial
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




def attribute_xs():
    pass

def get_num_lines_alignment_file():
    pass

def get_restriction_table():
    pass

def write_frag_info():
    pass

def get_chromosomes_sizes():
    pass

def get_bin_table():
    pass

def is_duplicated():
    pass

def is_poor_quality():
    pass

def is_unqualitative():
    pass

def is_unmapped():
    pass

def is_reverse():
    pass

def classify_reads():
    pass

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

def sam_iterator():
    pass

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
