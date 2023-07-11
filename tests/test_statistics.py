import pytest
from os.path import join
from pathlib import Path
import statistics 
import random
import itertools

import numpy as np
import scipy
from scipy.spatial.distance import pdist
import statsmodels.api as sm
import pandas as pd
import matplotlib.pyplot as plt

import pysam
from Bio import SeqIO, Restriction
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.Restriction import *

import cooler

import hicberg.statistics as hst
import hicberg.utils as hut


lowess = sm.nonparametric.lowess

GENOME = "data_test/SC288_with_micron.fa"
DEFAULT_ENZYME = ["DpnII"]
DICT_FIRST_KEY = "chr10"
DICT_FIRST_CHR_FIRST_POS = 0
DICT_FIRST_CHR_LAST_POS = 745751
LENGTH_DPNII = 2266
LENGTH_DPNII_HINFI = 4680
SUBSAMPLE_RATE = 1.0


@pytest.mark.parametrize("genome, enzyme, length", [(GENOME, DEFAULT_ENZYME, LENGTH_DPNII), (GENOME, ["DpnII", "HinfI"], LENGTH_DPNII_HINFI)])
def test_get_restriction_map(genome, enzyme, length):
    """
    Test if the function returns a dictionary with the right keys and values, i.e. the genome has been correctly digested.
    First restriction site for first chromosome should be 0 and last restriction site should be the length of the chromosome.
    """

    genome_path  = Path(genome)
    restriction_map = hst.get_restriction_map(genome = genome_path, enzyme = enzyme)


    # yield restriction_map

    # Check if the dictionary is not empty, if fisrt chromosome is right, and restrictions sites are in the right order.
    assert restriction_map[DICT_FIRST_KEY][0] == DICT_FIRST_CHR_FIRST_POS
    assert restriction_map[DICT_FIRST_KEY][-1] == DICT_FIRST_CHR_LAST_POS
    #Check if the length of the restriction map is right.
    assert len(restriction_map[DICT_FIRST_KEY]) == length

@pytest.mark.parametrize("genome, enzyme, length, rate", [(GENOME, DEFAULT_ENZYME, LENGTH_DPNII, 1.0), (GENOME, ["DpnII", "HinfI"], LENGTH_DPNII_HINFI, 1.0), (GENOME, DEFAULT_ENZYME, LENGTH_DPNII, 0.5), (GENOME, ["DpnII", "HinfI"], LENGTH_DPNII_HINFI, 0.5), (GENOME, DEFAULT_ENZYME, LENGTH_DPNII, 0.2), (GENOME, ["DpnII", "HinfI"], LENGTH_DPNII_HINFI, 0.2)])
def test_subsample_restriction_map(genome, enzyme, length, rate):

    genome_path  = Path(genome)
    restriction_map = hst.get_restriction_map(genome = genome_path, enzyme = enzyme)

    print(f"Rate : {rate}")
    subsampled_restriction_map = hut.subsample_restriction_map(restriction_map = restriction_map, rate = rate)


    assert subsampled_restriction_map[DICT_FIRST_KEY][0] == DICT_FIRST_CHR_FIRST_POS
    assert subsampled_restriction_map[DICT_FIRST_KEY][-1] == DICT_FIRST_CHR_LAST_POS
    #Check if the length of the restriction map is right. (+/-1 because of the random subsampling on even or odd numbers)
    assert float(len(subsampled_restriction_map[DICT_FIRST_KEY])) -1 <= length *  rate <= float(len(subsampled_restriction_map[DICT_FIRST_KEY])) + 1


def test_generate_xs():
    pass

def test_attribute_xs():
    pass

def test_get_dist_frags():
    pass

def test_generate_trans_ps():
    pass

def test_generate_probabilities():
    pass

def test_filter_poor_covered():
    pass

def test_generate_coverages():
    pass

def test_compute_d1d2():
    pass

def test_get_stats():
    pass
