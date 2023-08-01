import os
import subprocess as sp
import glob
import pathlib
import tempfile as tmpf
import shutil
from itertools import combinations

from scipy import spatial
from scipy.stats import median_abs_deviation, pearsonr
import bioframe as bf

import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

import hicberg.eval as hev

from .conftest import temporary_folder
from .test_utils import test_get_chromosomes_sizes

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

POSITION = 5000
CHROMOSOME = "chr1"
BINS = 2000
BOUNDARIES = (4000, 6000)


def test_get_intervals_proportions():

    pass


def test_get_chromosomes_intervals():
    pass

def test_draw_intervals():
    pass

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


