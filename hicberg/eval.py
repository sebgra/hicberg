import os
import subprocess as sp
import glob
from pathlib import Path
import tempfile as tmpf
import shutil
from itertools import combinations

import numpy as np
from numpy.random import choice

from scipy import spatial
from scipy.stats import median_abs_deviation, pearsonr
import bioframe as bf

import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

import hicberg.io as hio


def get_intervals_proportions(chrom_sizes_dict : str = "chromosome_sizes.npy", nb_intervals : int = 1) -> dict[str, float]:
    """
    AI is creating summary for get_intervals_proportions

    Parameters
    ----------
    chrom_sizes_dict : str
        Path to a dictionary containing chromosome sizes as {chromosome : size} saved in .npy format. By default chromosome_sizes.npy
    nb_intervals : int, optional
        Number of intervals to draw, by default 1

    Returns
    -------
    dict[str, float]
        Dictionary containing proportion by intervals as {chromosome : proportion}.
    """ 
    
    chrom_sizes_path = Path(chrom_sizes_dict)
    chrom_sizes = hio.load_dictionary(chrom_sizes_path)

    # Get genome global length
    tot_genome_length = np.sum(list(chrom_sizes.values()))

    # Compute relative proportion of each chromosome through genome
    chr_proportions  = {k : (v / tot_genome_length)  for (k, v) in zip(chrom_sizes.keys(), chrom_sizes.values())}

    # Draw chromosomes considering their relative prpoportion in genome
    proportions_choice = choice(a = list(chr_proportions.keys()), size = nb_intervals, replace = True, p = list(chr_proportions.values()))

    # Get and return chromosomes picked and counts
    unique, counts = np.unique(proportions_choice, return_counts=True)

    return dict(zip(unique, counts))



def get_chromosomes_intervals():
    pass

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

    # TODO : implement get_intervals_proportion

    pass  

def draw_positions():
    pass

def is_overlapping():
    pass

def is_overlapping():
    pass

def get_boundaries(position : int = None, bins : int = 2000, chromosome : str = None, chrom_sizes_dict : str = "chromosome_sizes.npy") -> tuple[int, int]:
    """
    Return boundaries surrounding position considering regular binning of bins.

    Parameters
    ----------
    position : int, optional
        Position to get surounding boundaries from, by default None
    bins : int, optional
        Size of the desired bin, by default 2000., by default 2000
    chromosome : str, optional
        [description], by default None
    chrom_sizes_dict : str
        Path to a dictionary containing chromosome sizes as {chromosome : size} saved in .npy format. By default chromosome_sizes.npy

    Returns
    -------
    tuple[int, int]
        Tuple containing the boundaries of the position : (lower_bound, upper_bound).
    """

    chrom_sizes_path = Path(chrom_sizes_dict)
    chrom_sizes = hio.load_dictionary(chrom_sizes_path)

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


def get_intervals_indexes():
    pass

def order_reads():
    pass

def get_bin_enrichment():
    pass

def check_empty_bins():
    pass


def get_intervals_gap():
    pass

def select_pattern_reads():
    pass

def get_interest_bins():
    pass

def get_interest_indexes():
    pass

def build_recovered_pairs():
    pass

def build_recovered_matrix():
    pass

def nan_sanitize():
    pass

def max_consecutive_nans():
    pass

def new_smoothing():
    pass

def compute_crop_coordinates():
    pass

def compute_crop_coordinates():
    pass

def compute_likelyhood():
    pass


def compute_comat_score():
    pass

def compute_cosine_score():
    pass

def compute_matrix_correlation():
    pass

def plot_1D_correlation():
    pass

def get_1D_correlation():
    pass

def sew_maps():
    pass

def plot_reconstructed_map():
    pass

def get_trans_bins():
    pass

def zoom_stripe():
    pass


def get_bin_indexes():
    pass

def check_reco_sum():
    pass

def grid_search(mode : list[str] = [], chromosome : list[str] = [], trans_chromosome : list[str] = [], position : list[str] = [], trans_position : list[str] = [], stride : list[int] = []) -> None:
    pass



