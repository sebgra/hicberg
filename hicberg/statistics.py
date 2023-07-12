from os import getcwd
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
import hicberg.utils as hut


lowess = sm.nonparametric.lowess

RESTRICTION_DICO = "dist.frag.npy"
XS = "xs.npy"
COVERAGE_DICO = "coverage.npy"

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


def generate_xs(chromosome_size : int, base : float = 1.1) -> dict[str, np.ndarray[int]]:
    """
    Generate xs array for computing and plotting P(s). Return xs array which is logspace.

    Parameters
    ----------
    chromosome_size : int
        Size of the chromosome to be binned in bp.
    base : float, optional
        Base of the logspace., by default 1.1

    Returns
    -------
    dict[str, np.ndarray[int]]
        Array of log bins related to the chromosome.
    """

    n_bins = np.divide(np.log1p(chromosome_size), np.log(base)).astype(int)
    xs = np.unique(
        np.logspace(0, n_bins + 1, base=base, num=n_bins + 2, endpoint=True, dtype=int)
    )


    return xs

def attribute_xs(xs : np.ndarray[int], distance : int) -> int:
    """
    Attibute genomic distance to the corresponding log bin of xs.

    Parameters
    ----------
    xs : np.ndarray[int]
        Array containing the log bins.
    distance : int
        Genomic distance in bp.

    Returns
    -------
    int
        Index of the corresponding bins where distance has to be attributed.
    """

    idx = np.searchsorted(xs, distance, side="right") - 1
    return idx

def get_dist_frags(genome : str = None, restriction_map : dict = None, circular : str = None, rate : float = 1.0, output_dir : str = None) -> None:
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
    rate : float, optional
        Set the proportion of restriction sites to condiser. Avoid memory overflow when restriction maps are very dense, by default 1.0
    output_dir : str, optional
        Path to the folder where to save the dictionary, by default None

    Returns
    -------
    dict
        Dictionary of subsampled restriction map with keys as chromosome names and values as lists of restriction sites' position.
    """

    if output_dir is None:

        folder_path = Path(getcwd())

    else:

        folder_path = Path(output_dir)

    if (0.0 > rate) or (rate > 1.0):
        raise ValueError("Subsampling rate must be between 0.0 and 1.0.")

    genome_path = Path(genome)

    if not genome_path.is_file():
            
        raise FileNotFoundError(f"Genome file {genome} not found. Please provide a valid path to a genome file.")
    
    dist_frag = dict()
    xs = dict()
    
    if rate != 1.0:

        restriction_map = hut.subsample_restriction_map(restriction_map = restriction_map, subsample = rate)

    for seq_record in SeqIO.parse(genome, "fasta"):

        if seq_record.id in restriction_map.keys():

            seq_name = seq_record.id
            map_size = restriction_map[seq_name].shape[0]

            forward_distances = pdist(
                np.reshape(restriction_map[seq_name], (map_size, 1))
            )
            max_size_vector = np.full(
                forward_distances.shape, np.max(restriction_map[seq_name])
            )
            backward_distances = max_size_vector - forward_distances
            pairwise_distances = np.minimum(forward_distances, backward_distances)
            pairwise_distances = np.delete(
                pairwise_distances, np.where(pairwise_distances == 0)
            )

            # freeing memory
            del forward_distances
            del backward_distances
            del max_size_vector

        else : 

            pairwise_distances = pdist(
                np.reshape(
                    restriction_map[seq_name],
                    (len(restriction_map[seq_name]), 1),
                )
            )
            pairwise_distances = np.delete(
                pairwise_distances, np.where(pairwise_distances == 0)
            )

        # Computing xs
        xs[seq_name] = generate_xs(len(seq_record.seq), base=1.1)
        dist_frag[seq_name] = np.zeros(xs[seq_name].shape)
        # Parse distances
        for distance in pairwise_distances:
            dist_frag[seq_name][attribute_xs(xs[seq_name], distance)] += 1



    # Save dictionaries
    np.save(
        folder_path / RESTRICTION_DICO,
        dist_frag,
    )

    print(f"Restriction map saved in {folder_path}")

def generate_trans_ps():
    pass

def generate_probabilities():
    pass

def filter_poor_covered():
    pass

def generate_coverages(genome : str = None, bins : int = 2000, forward_bam_file : str = "group1.1.bam", reverse_bam_file : str = "group1.2.bam", output_dir : str = None) -> None:
    """
    Take a genome and  both for and rev sam files for unambiguous group and return a dictionary containing the coverage in terms of reads overs chromosomes. .

    Parameters
    ----------
    genome : str, optional
        Path to the genome file to get coverage on ., by default None
    bins : int, optional
        Size of the desired bin., by default 2000
    forward_bam_file : str, optional
        Path to forward .bam alignment file, by default None, by default group1.1.bam
    reverse_bam_file : str, optional
        Path to reverse .bam alignment file, by default None, by default group1.2.bam
    output_dir : str, optional
        Path to the folder where to save the classified alignment files, by default None, by default None
    """        
    
    if output_dir is None:

        output_path = Path(getcwd())

    else:

        output_path = Path(output_dir)
    
    genome_path = Path(genome)

    if not genome_path.is_file():
                
        raise FileNotFoundError(f"Genome file {genome} not found. Please provide a valid path to a genome file.")   
    
    genome_parser = SeqIO.parse(genome, "fasta")
    genome_coverages = {seq_record.id : np.zeros(np.round(np.divide(len(seq_record.seq), bins) + 1).astype(int)) for seq_record in genome_parser}

    forward_bam_path = Path(forward_bam_file)
    reverse_bam_path = Path(reverse_bam_file)

    if not forward_bam_path.is_file():

        raise FileNotFoundError(f"Forward .bam file {forward_bam_file} not found. Please provide a valid path to a forward .bam file.")
    
    if not reverse_bam_path.is_file():

        raise FileNotFoundError(f"Reverse .bam file {reverse_bam_file} not found. Please provide a valid path to a reverse .bam file.")
    
    forward_bam_handler, reverse_bam_handler = pysam.AlignmentFile(forward_bam_path, "rb"), pysam.AlignmentFile(reverse_bam_path, "rb")

    for forward_read, reverse_read in zip(forward_bam_handler, reverse_bam_handler):

        genome_coverages[forward_read.reference_name][
            np.divide(forward_read.pos, bins).astype(int)
        ] += 1
        genome_coverages[reverse_read.reference_name][
            np.divide(reverse_read.pos, bins).astype(int)
        ] += 1

    # close files
    forward_bam_handler.close()
    reverse_bam_handler.close()

    # Smooth coverages
    smoothed_coverages = {seq_name : hut.mad_smoothing(coverage) for seq_name, coverage in genome_coverages.items()}

    np.save(output_path / COVERAGE_DICO, smoothed_coverages)



def compute_d1d2():
    pass

def get_stats():
    pass
