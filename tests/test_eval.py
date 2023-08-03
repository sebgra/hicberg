import os
import subprocess as sp
import glob
from pathlib import Path
import tempfile as tmpf
import shutil
from itertools import combinations

import pytest


import numpy as np
from scipy import spatial
from scipy.stats import median_abs_deviation, pearsonr
import bioframe as bf

import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

import cooler

from .conftest import random
import hicberg.eval as hev
import hicberg.io as hio

from .conftest import temporary_folder
from .test_align import test_hic_build_index, test_hic_align, test_hic_view, test_hic_sort

from .test_utils import test_get_chromosomes_sizes, test_classify_reads, test_hic_sort, test_get_bin_table
from .test_io import test_build_matrix, test_build_pairs

GENOME = "data_test/SC288_with_micron.fa"
CHROM_SIZES_DIC = "chromosome_sizes.npy"
MATRIX = "data_test/unrescued_map.cool"
FRAG_FILE = "fragments_fixed_sizes.txt"

FORWARD_SORTED_BAM = "data_test/alignments/1.sorted.bam"

FORWARD_UNMAPPED_BAM = "group0.1.bam"
REVERSE_UNMAPPED_BAM = "group0.2.bam"
FORWARD_UNIQUE_BAM = "group1.1.bam"
REVERSE_UNIQUE_BAM = "group1.2.bam"
FORWARD_MULTI_BAM = "group2.1.bam"
REVERSE_MULTI_BAM = "group2.2.bam"

NB_INTERVALS = 2
NB_INTERVALS_EMPTINESS = 10
POSITION = 5000
CHROMOSOME = "chr1"
BINS = 2000
FRIST_INTERVAL = (2 * BINS, 2* BINS + BINS)
LAST_INTERVAL = (230000 - BINS, 230000)
BOUNDARIES = (4000, 6000)
PROPORTIONS = {'chr2': 1, 'chr4': 1} # Seed set
INTERVALS_DICT = {'chr2': [(0, 1000), (1000, 2000)], 'chr4': [(0, 1000), (1000, 2000)], 
                'chr1': [(0, 1000), (1000, 2000), (4000, 6000)]} # Seed set

DRAWN_INTERVALS =   {'chr2': [(756000, 758000)], 'chr3': [(54000, 56000)]}



def test_get_interval_index(test_get_chromosomes_sizes):
    """
    Test if interval belonging is correctly computed.
    """

    chrom_size_dict = hio.load_dictionary(test_get_chromosomes_sizes)
    indexes = hev.get_interval_index(chromosome = CHROMOSOME, value = POSITION, intervals_dict = INTERVALS_DICT,  chrom_sizes_dict = chrom_size_dict)

    assert indexes[CHROMOSOME][1] == BOUNDARIES
    assert indexes[CHROMOSOME][0] == INTERVALS_DICT[CHROMOSOME][:2]

    assert indexes["chr8"][0] == [(None, None)]
    assert indexes["chr8"][1] == (None, None)


def test_select_reads(temporary_folder, test_classify_reads, test_build_matrix, test_get_chromosomes_sizes):
    hev.select_reads(bam_for = test_classify_reads[0], bam_rev = test_classify_reads[1], matrix_file = test_build_matrix, output_dir = temporary_folder, chromosome = CHROMOSOME)

    temp_dir_path = Path(temporary_folder)
    forward_in_path = temp_dir_path / "group1.1.in.bam"
    reverse_in_path = temp_dir_path / "group1.2.in.bam"
    forward_out_path = temp_dir_path / "group1.1.out.bam"
    reverse_out_path = temp_dir_path / "group1.2.out.bam"

    assert forward_in_path.is_file()
    assert reverse_in_path.is_file()
    assert forward_out_path.is_file()
    assert reverse_out_path.is_file()

def test_get_intervals_proportions(random, test_get_chromosomes_sizes):
    """
    Test if the dictionary of proportion is correctly computed.
    """
    proportions = hev.get_intervals_proportions(chrom_sizes_dict = test_get_chromosomes_sizes, nb_intervals = NB_INTERVALS)
    assert proportions == PROPORTIONS

def test_get_chromosomes_intervals(test_get_chromosomes_sizes):
    """
    Test if all the intervals possible are correctly computed.
    """
    intervals = hev.get_chromosomes_intervals(bins = BINS, chrom_sizes_dict = test_get_chromosomes_sizes, chromosome = CHROMOSOME)
    # print(f"Intervals: {intervals}")
    assert intervals[0] == FRIST_INTERVAL
    assert intervals[-1] == LAST_INTERVAL

def test_draw_intervals(test_get_chromosomes_sizes):
    """
    Test if interval drawing is correctly computed.
    """
    intervals = hev.draw_intervals(nb_intervals = NB_INTERVALS, bins = BINS, chrom_sizes_dict = test_get_chromosomes_sizes)
    assert intervals == DRAWN_INTERVALS


def test_draw_positions():
    pass

def test_is_overlapping():
    pass

def test_is_overlapping():
    pass

def test_get_boundaries(test_get_chromosomes_sizes):
    """
    Test if boundaries from genomic coordinates are correctly computed.
    """
    boundaries = hev.get_boundaries(position = POSITION, bins = BINS, chromosome = CHROMOSOME, chrom_sizes_dict = test_get_chromosomes_sizes)

    assert boundaries == BOUNDARIES

def test_emptyness(test_get_chromosomes_sizes):
    """
    Test if emptiness checking before considering genomic intervals is correctly computed.
    """

    clr = cooler.Cooler(MATRIX)
    dictionary_of_intervals = hev.draw_intervals(chrom_sizes_dict  = test_get_chromosomes_sizes, nb_intervals = NB_INTERVALS_EMPTINESS, bins = BINS)
    emptiness = hev.check_emptiness(intervals = dictionary_of_intervals, matrix = clr)

    assert not emptiness

def test_get_intervals_indexes():
    pass

def test_order_reads():
    pass

def test_get_bin_enrichment():
    pass

def test_check_empty_bins():
    pass


def test_get_intervals_gap():
    pass

def test_select_pattern_reads():
    pass

def test_get_interest_bins():
    pass

def test_get_interest_indexes():
    pass

def test_build_recovered_pairs():
    pass

def test_build_recovered_matrix():
    pass

def test_nan_sanitize():
    pass

def test_max_consecutive_nans():
    pass

def test_new_smoothing():
    pass

def test_compute_crop_coordinates():
    pass

def test_compute_crop_coordinates():
    pass

def test_compute_likelyhood():
    pass


def test_compute_comat_score():
    pass

def test_compute_cosine_score():
    pass

def test_compute_matrix_correlation():
    pass

def test_plot_1D_correlation():
    pass

def test_get_1D_correlation():
    pass

def test_sew_maps():
    pass

def test_plot_reconstructed_map():
    pass

def test_get_trans_bins():
    pass

def test_zoom_stripe():
    pass


def test_get_bin_indexes():
    pass

def test_check_reco_sum():
    pass

def test_grid_search():
    pass


