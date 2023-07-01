from os.path import join
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


lowess = sm.nonparametric.lowess

def test_get_restriction_map():
    pass

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
