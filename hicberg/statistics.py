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


lowess = sm.nonparametric.lowess

def get_restriction_map(genome : str = None, enzyme : list[str] = ["DpnII"]) -> dict[str, np.ndarray[int]]:
    """
    Get ordered restriction map (including 0 and n) from a chromosome sequence.
    Return a dictionary where keys are chromosomes names and values are restrictions sites positions.

    Parameters
    ----------
    genome : str, optional
        Path to the genome to digest, by default None, by default None
    enzyme : list[str], optional
        Enzyme or list of enzyme to digest the genome with., by default None, by default ["DpnII"]

    Returns
    -------
    dict
        Dictionary of the product of digestion where keys are chromosomes names and values are restrictions sites positions.
    """

    genome_path = Path(genome)

    if not genome_path.is_file():

        raise FileNotFoundError(f"Genome file {genome} not found. Please provide a valid path to a genome file.")
    
    
    restriction_map_dictionary = dict()
    
    restriction_batch = Restriction.RestrictionBatch()
    for enz in enzyme:
        restriction_batch.add(enz)

    # parse sequence from fasta file
    for seq_record in SeqIO.parse(genome, "fasta"):

        # Get restriction map from the restriction batch.
        restriction_map = restriction_batch.search(seq_record.seq)

        # Convert dictionary values to numpy array
        restriction_map_array = np.sort(
            np.array([pos for sites in restriction_map.values() for pos in sites])
        )
        restriction_map_array = np.insert(
            restriction_map_array,
            [0, len(restriction_map_array)],
            [0, len(seq_record.seq)],
        )

        restriction_map_dictionary[seq_record.id] = restriction_map_array

    return restriction_map_dictionary


def generate_xs():
    pass

def attribute_xs():
    pass

def get_dist_frags(genome : str = None, restriction_map : dict = None, circular : str = None, subsample : float = 1.0) -> None:
    """
    Get the distribution of fragments' distance from a genome distribution .

    Parameters
    ----------
    genome : str, optional
        Path to the genome, by default None
    restriction_map : dict, optional
        Restriction map saved as a dictionary like chrom_name : list of restriction sites' position, by default None
    circular : str, optional
        Name of the chromosomes to consider as circular, by default None
    subsample : float, optional
        Set the proportion of restriction sites to condiser. Avoid memory overflow when restriction maps are very dense, by default 1.0

    Returns
    -------
    dict
        Dictionary of subsampled restriction map with keys as chromosome names and values as lists of restriction sites' position.
    """

    
    genome_path = Path(genome)

    if not genome_path.is_file():
            
        raise FileNotFoundError(f"Genome file {genome} not found. Please provide a valid path to a genome file.")
    
    if subsample == 1.0:

        restriction_map_dictionary = None

    return 0

def generate_trans_ps():
    pass

def generate_probabilities():
    pass

def filter_poor_covered():
    pass

def generate_coverages():
    pass

def compute_d1d2():
    pass

def get_stats():
    pass
