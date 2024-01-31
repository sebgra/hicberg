import time
import sys
from os import getcwd
from os.path import join
from pathlib import Path
import uuid
import multiprocessing as mp
from functools import partial

import itertools

import numpy as np
from numpy.random import choice
from scipy.spatial.distance import pdist
from scipy.stats import pearsonr
from scipy.ndimage import gaussian_filter
import statsmodels.api as sm
import pandas as pd
import matplotlib.pyplot as plt

import pysam
from Bio import SeqIO, Restriction
from Bio.Restriction import *

import cooler
import hicberg.utils as hut
import hicberg.io as hio

from hicberg import logger



lowess = sm.nonparametric.lowess

DIST_FRAG = "dist.frag.npy"
XS = "xs.npy"
COVERAGE_DICO = "coverage.npy"
D1D2 = "d1d2.npy"
UNCUTS = "uncuts.npy"
WEIRDS = "weirds.npy"
LOOPS = "loops.npy"
TRANS_PS = "trans_ps.npy"
RESTRICTION_MAP = "restriction_map.npy"
DENSITY_MAP = "density_map.npy"

# TODO : check if keeping is necessary
def generate_density_map_backup(matrix : str = "unrescued_map.cool", rounds : int = 1, magnitude : float = 1.0, output_dir : str = None) -> dict[str, np.ndarray[float]]:
    """
    Create density map from a Hi-C matrix. Return a dictionary where keys are chromosomes names and values are density maps.

    Parameters
    ----------
    matrix : str, optional
        Path to a cooler matrix, by default None
    rounds : int, optional
        Number of times the matrix has to be shaken, by default 1
    magnitude : float, optional
        Blending ratio between the native matrix and the shaken one, by default 1.0
    output_dir : str, optional
        Path to the folder where to save the density map, by default None
    Returns
    -------
    dict[str, np.ndarray[float]]
        Density maps as a dictionary where keys are chromosomes names couples as tuples and values are density maps.
    """

    logger.info("Start generating density map...")

    if output_dir is None:

        output_path = Path(getcwd())

    else : 

        output_path = Path(output_dir)

    matrix_path = output_path / matrix

    if not matrix_path.is_file():

        raise FileNotFoundError(f"Matrix file {matrix} not found. Please provide a valid path to a matrix file.")
    
    density_map = {}

    # Load cooler matrix

    matrix = hio.load_cooler(matrix_path)

    # Get chromosomes names

    chromosomes = matrix.chromnames

    chromosomes_combination = list(itertools.product(chromosomes, chromosomes))

    for combination in chromosomes_combination:

        # Get matrix for each chromosome pair

        matrix_chromosome = matrix.matrix(balance = False).fetch(combination[0], combination[1])

        # Get density map for each chromosome pair

        if combination[0] == combination[1]:

            density_map[combination] = hut.diffuse_matrix(matrix = matrix_chromosome, rounds = rounds, magnitude = magnitude, mode = "intra")

        elif combination[0] != combination[1]:

            density_map[combination] = hut.diffuse_matrix(matrix = matrix_chromosome, rounds = rounds, magnitude = magnitude, mode = "inter")

        blured_density  = dict()

    for key in density_map.keys():

        if key[0] == key[1]:

            blured_density[key] = gaussian_filter(density_map[key], sigma = 5)

        else :

            blured_density[key] = gaussian_filter(density_map[key], sigma = 10)

    np.save(output_path / DENSITY_MAP, blured_density)

    logger.info(f"Saved density map at : {output_path}")

    return blured_density


def generate_density_map(matrix : str = "unrescued_map.cool", size : int = 5, sigma : int = 2, n_mads : int = 2, nan_threshold : bool = False, output_dir : str = None) -> dict[str, np.ndarray[float]]:
    """
    Create density map from a Hi-C matrix. Return a dictionary where keys are chromosomes names and values are density maps.

    Parameters
    ----------
    matrix : str, optional
        Path to a cooler matrix, by default None
    rounds : int, optional
        Number of times the matrix has to be shaken, by default 1
    magnitude : float, optional
        Blending ratio between the native matrix and the shaken one, by default 1.0
    output_dir : str, optional
        Path to the folder where to save the density map, by default None
    Returns
    -------
    dict[str, np.ndarray[float]]
        Density maps as a dictionary where keys are chromosomes names couples as tuples and values are density maps.
    """

    logger.info("Start generating density map...")

    if output_dir is None:

        output_path = Path(getcwd())

    else : 

        output_path = Path(output_dir)

    matrix_path = output_path / matrix

    if not matrix_path.is_file():

        raise FileNotFoundError(f"Matrix file {matrix} not found. Please provide a valid path to a matrix file.")
    
    density_map = {}

    # Load cooler matrix

    matrix = hio.load_cooler(matrix_path)

    # Get chromosomes names

    chromosomes = matrix.chromnames

    chromosomes_combination = list(itertools.product(chromosomes, chromosomes))

    for combination in chromosomes_combination:

        # Get matrix for each chromosome pair

        matrix_chromosome = matrix.matrix(balance = True).fetch(combination[0], combination[1])

        # Get density map for each chromosome pair

        density_map[combination] = hut.get_local_density(matrix = matrix_chromosome, size  = size, sigma  = sigma, n_mads  = n_mads, nan_threshold  = nan_threshold)


    np.save(output_path / DENSITY_MAP, density_map)

    logger.info(f"Saved density map at : {output_path}")

    return density_map

def compute_density(cooler_file : str = None, threads : int = 2, kernel_size: int = 11, deviation : float = 0.5, output_dir : str = None) ->  None: 
    """
    Create density map from a Hi-C matrix. Return a dictionary where keys are chromosomes names and values are density maps.

    Parameters
    ----------
    cooler_file : str, optional
        [description], by default None
    threads : int, optional
        [description], by default 2
    output_dir : str, optional
        [description], by default None

    """
    logger.info("Start generating density map...")

    if output_dir is None:

        output_path = Path(getcwd())

    else : 

        output_path = Path(output_dir)

    matrix_path = output_path / cooler_file

    print(f"Matrix path : {matrix_path}")

    if not matrix_path.is_file():

        raise FileNotFoundError(f"Matrix file {matrix} not found. Please provide a valid path to a matrix file.")

    #Load cooler file
    matrix = hio.load_cooler(matrix = matrix_path)

    #Get chromosomes names
    chromosomes = matrix.chromnames

    #Get chromosomes couples
    chromosomes_couples = list(itertools.product(chromosomes, repeat = 2))

    # Get chromsomes maps
    chromosomes_maps = [matrix.matrix(balance = True).fetch(chrom1, chrom2) for chrom1, chrom2 in chromosomes_couples]

    pool = mp.Pool(processes=threads)

    results = pool.map(partial(hut.get_local_density, str(matrix_path), size  = kernel_size, sigma  = deviation, n_mads  = 2, nan_threshold  = False),
            chromosomes_couples)

    # Close the pool and wait for the work to finish
    pool.close()
    pool.join()

    results_dict =  {key : value for key, value in results}

    # TODO : to adjust to itertools meodification : so far itertools.product(chromosomes, repeat = 2) is used
    for chrom_pair in results_dict.copy().keys():
        if chrom_pair[0] == chrom_pair[1]:
            pass

        else :
            results_dict[(chrom_pair[1], chrom_pair[0])]  = results_dict[chrom_pair].T

    np.save(output_path / DENSITY_MAP, results_dict)

    logger.info(f"Saved density map at : {output_path}")



def get_restriction_map(genome : str = None, enzyme : list[str] = ["DpnII"], output_dir : str = None) -> dict[str, np.ndarray[int]]:
    """
    Get ordered restriction map (including 0 and n) from a chromosome sequence.
    Return a dictionary where keys are chromosomes names and values are restrictions sites positions.

    Parameters
    ----------
    genome : str, optional
        Path to the genome to digest, by default None, by default None
    enzyme : list[str], optional
        Enzyme or list of enzyme to digest the genome with. If integer passed, micro-C mode using Mnase is used, and the integer correspond to the size of nucleosomal fragment, by default None, by default ["DpnII"]
    output_dir : str, optional
        Path to the folder where to save the plot, by default None

    Returns
    -------
    dict
        Dictionary of the product of digestion where keys are chromosomes names and values are restrictions sites positions.
    """

    logger.info("Generating restriction map...")

    genome_path = Path(genome)

    if not genome_path.is_file():

        raise FileNotFoundError(f"Genome file {genome} not found. Please provide a valid path to a genome file.")

    if output_dir is None:

        output_path = Path(getcwd())

    else : 

        output_path = Path(output_dir)
    
    
    restriction_map_dictionary = dict()

    if len(enzyme) == 1 and enzyme[0].isnumeric():
        print(f"Enabled micro-c with enzyme : {enzyme}")
        enzyme = int(enzyme[0])

        for seq_record in SeqIO.parse(genome, "fasta"):

            # Get restriction map from the restriction batch.
            restriction_map = np.arange(0, len(seq_record.seq), enzyme)
            restriction_map = np.insert(
                restriction_map,
                [len(restriction_map)],
                [len(seq_record.seq)],
            )
            restriction_map_dictionary[seq_record.id] = restriction_map

    elif type(enzyme) == str or type(enzyme) == list or type(enzyme) == tuple:

        print(f"Non micro-C mode with enzyme : {enzyme}")

        restriction_batch = Restriction.RestrictionBatch()

        if len(enzyme) == 1:
            enzyme = enzyme[0].split(",")
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


    np.save(output_path / RESTRICTION_MAP, restriction_map_dictionary)

    logger.info(f"Saved restriction map at : {output_path}")
    

    return restriction_map_dictionary


def generate_xs(chromosome_size : int, base : float = 1.1) -> np.ndarray[int]:
    """
    Generate xs array for computing P(s). Return xs array which is log space.

    Parameters
    ----------
    chromosome_size : int
        Size of the chromosome to be binned in bp.
    base : float, optional
        Base of the log space., by default 1.1

    Returns
    -------
    np.ndarray[int]
        Array of log bins related to the chromosome.
    """

    
    n_bins = np.divide(np.log1p(chromosome_size), np.log(base)).astype(int)
    xs = np.unique(
            np.logspace(0, n_bins, num=n_bins + 1, base=base, dtype=int)
        )
    # xs[-1] = chromosome_size

    return xs

def log_bin_genome(genome :str, base : float = 1.1, output_dir : str = None) -> dict[str, np.ndarray[int]]:
    
    logger.info("Start log binning of genome...")

    genome_path = Path(genome)

    if not genome_path.is_file():

        raise FileNotFoundError(f"Genome file {genome} not found. Please provide a valid path to a genome file.")
    
    if output_dir is None:
            
        folder_path = Path(getcwd())

    else:
            
        folder_path = Path(output_dir)

    genome_parser = SeqIO.parse(genome, "fasta")
    xs_dict = {seq_record.id : generate_xs(chromosome_size = len(seq_record.seq), base = base) for seq_record in genome_parser}

    np.save(folder_path / XS, xs_dict)

    logger.info(f"Log binning of genome {genome} saved in {folder_path / XS}.")

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

def get_dist_frags(genome : str = None, restriction_map : dict = None, circular : str = "", rate : float = 1.0, output_dir : str = None) -> None:
    """
    Get the distribution of fragments' distance from a genome distribution .

    Parameters
    ----------
    genome : str, optional
        Path to the genome, by default None
    restriction_map : dict, optional
        Restriction map saved as a dictionary like chrom_name : list of restriction sites' position, by default None
    circular : str, optional
        Name of the chromosomes to consider as circular, by default ""
    rate : float, optional
        Set the proportion of restriction sites to consider. Avoid memory overflow when restriction maps are very dense, by default 1.0
    output_dir : str, optional
        Path to the folder where to save the dictionary, by default None

    Returns
    -------
    dict
        Dictionary of sub-sampled restriction map with keys as chromosome names and values as lists of restriction sites' position.
    """
    logger.info("Start generating distribution of fragments' distance...")

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

        restriction_map = hut.subsample_restriction_map(restriction_map = restriction_map, rate = rate)

    for seq_record in SeqIO.parse(genome, "fasta"):

        seq_name = seq_record.id

        
        if seq_record.id in circular:

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
        folder_path / DIST_FRAG,
        dist_frag,
    )

    logger.info(f"Saved restriction map at : {folder_path / DIST_FRAG}")

def generate_trans_ps(matrix : str = "unrescued_map.cool", chrom_sizes : str = "chromosome_sizes.npy", output_dir : str = None) -> None:
    

    logger.info("Start getting trans-P(s)")

    if output_dir is None:

        output_path = Path(getcwd())

    else:

        output_path = Path(output_dir)

    chrom_size_dict = hio.load_dictionary(output_path / chrom_sizes)

    matrix_path = Path(output_path, matrix)

    if not matrix_path.is_file():

        raise FileNotFoundError(f"Matrix file {matrix} not found. Please provide a valid path to a matrix file.")
    
    # logger.info(f"Loading matrix {matrix.name}...")
    
    matrix = cooler.Cooler(matrix_path.as_posix())

    chromosome_sets = itertools.product((chrom_size_dict.keys()), repeat=2)

    trans_ps = {}

    t_ps = np.zeros((len(chrom_size_dict.keys()) ** 2, 1))
    all_interaction_matrix = np.zeros((len(chrom_size_dict.keys()) ** 2, 1))
    n_frags_matrix = np.zeros((len(chrom_size_dict.keys()) ** 2, 1))


    for idx, s in enumerate(chromosome_sets):

        all_interactions = matrix.matrix(balance=False).fetch(s[0], s[1]).sum()
        n_frags = chrom_size_dict.get(s[0]) * chrom_size_dict.get(s[1])
        
        # trans_ps[s] = np.divide(all_interactions, np.multiply(n_frags, 4))
        trans_ps[s] = np.divide(all_interactions, np.multiply(n_frags, 4)) # Multiplied by 4 to balance 4 configurations of reads orientation (++/+-/-+/--)


        # t_ps[idx] = np.divide(all_interactions, np.multiply(n_frags, 4))
        
        all_interaction_matrix[idx] = all_interactions
        n_frags_matrix[idx] = n_frags

    # t_ps = t_ps.reshape(
    #     (len(restriction_map.keys()), (len(restriction_map.keys())))
    # )
    # np.fill_diagonal(t_ps, np.nan)

    # all_interaction_matrix = all_interaction_matrix.reshape(
    #     (len(restriction_map.keys()), (len(restriction_map.keys())))
    # )
    # np.fill_diagonal(all_interaction_matrix, np.nan)

    # n_frags_matrix = n_frags_matrix.reshape(
    #     (len(restriction_map.keys()), (len(restriction_map.keys())))
    # )
    # np.fill_diagonal(n_frags_matrix, np.nan)

    np.save(output_path / TRANS_PS, trans_ps)

    logger.info(f"Trans P(s) saved in {output_path}")


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
    
    logger.info("Start generating coverages...")

    if output_dir is None:

        output_path = Path(getcwd())

    else:

        output_path = Path(output_dir)
    
    genome_path = Path(genome)

    if not genome_path.is_file():
                
        raise FileNotFoundError(f"Genome file {genome} not found. Please provide a valid path to a genome file.")   
    
    genome_parser = SeqIO.parse(genome, "fasta")

    genome_coverages = {seq_record.id : np.zeros(np.round(np.divide(len(seq_record.seq), bins) + 1).astype(int)) for seq_record in genome_parser}

    forward_bam_path = Path(output_dir, forward_bam_file)
    reverse_bam_path = Path(output_dir, reverse_bam_file)

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

    logger.info(f"Coverage dictionary saved in {output_path}")



def generate_d1d2(forward_bam_file : str = "group1.1.bam", reverse_bam_file : str = "group1.2.bam", restriction_map : str = "restriction_map.npy", output_dir : str = None) -> None:
    """
    Compute d1d2 distance laws with the given alignments and restriction map.

    Parameters
    ----------
    forward_bam_file : str, optional
        Path to forward .bam alignment file, by default None, by default group1.1.bam, by default "group1.1.bam"
    reverse_bam_file : str, optional
        Path to reverse .bam alignment file, by default None, by default group1.1.bam, by default "group1.2.bam"
    restriction_map : str, optional
        Restriction map saved as a dictionary like chrom_name : list of restriction sites' position, by default "dist.frag.npy"
    output_dir : str, optional
        Path to the folder where to save the dictionary, by default None, by default None
    """    
    
    logger.info("Start generating d1d2 law...")
    
    if output_dir is None:

        output_path = Path(getcwd())

    else:

        output_path = Path(output_dir)

    forward_bam_path = Path(output_path, forward_bam_file)
    reverse_bam_path = Path(output_path, reverse_bam_file)

    if not forward_bam_path.is_file():
            
        raise FileNotFoundError(f"Forward .bam file {forward_bam_file} not found. Please provide a valid path to a forward .bam file.")
    
    if not reverse_bam_path.is_file():
            
        raise FileNotFoundError(f"Reverse .bam file {reverse_bam_file} not found. Please provide a valid path to a reverse .bam file.")
    
    forward_bam_handler, reverse_bam_handler = pysam.AlignmentFile(forward_bam_path, "rb"), pysam.AlignmentFile(reverse_bam_path, "rb")

    # Ensure that the restriction map is a dictionary to be loaded
    try:

        restriction_map = hio.load_dictionary(output_path / restriction_map)

    except : 

        print(f"restriction map not found")
        pass

    print(f"keys : {restriction_map.keys()}")

    # d1d2 = {seq_name : np.zeros((len(restriction_map[seq_name]), len(restriction_map[seq_name]))) for seq_name in restriction_map.keys()}
    list_d1d2 = []  # list containing the (d1+d2) i.e size of the fragment to sequence

    for forward_read, reverse_read in zip(forward_bam_handler, reverse_bam_handler):

        # print(f"forward_read.reference_name : {forward_read.reference_name}")

        r_sites_forward_read = restriction_map[forward_read.reference_name]
        r_sites_reverse_read = restriction_map[reverse_read.reference_name]

        if forward_read.flag == 0 or forward_read.flag == 256:

            index = np.searchsorted(r_sites_forward_read, forward_read.pos, side="right")

            # try : 
            distance_1 = np.subtract(r_sites_forward_read[index], forward_read.pos)

            # except : 

            #     distance_1 = np.subtract(r_sites_forward_read[index - 1], forward_read.pos)


        elif forward_read.flag == 16 or forward_read.flag == 272:

            index = np.searchsorted(r_sites_forward_read, forward_read.reference_end, side="left")

            # try : 
            distance_1 = np.abs(
                np.subtract(forward_read.reference_end, r_sites_forward_read[index])
            )
            # except : 

            #     distance_1 = np.abs(
            #         np.subtract(forward_read.reference_end, r_sites_forward_read[index - 1])
            #     )


        if reverse_read.flag == 0 or reverse_read.flag == 256:

            index = np.searchsorted(
                r_sites_reverse_read, reverse_read.reference_start, side="right"
            )  # right

            # try : 
            distance_2 = np.subtract(r_sites_reverse_read[index], reverse_read.reference_start)

            # except : 

            #     distance_2 = np.subtract(r_sites_reverse_read[index - 1], reverse_read.reference_start)


        elif reverse_read.flag == 16 or reverse_read.flag == 272:

            index = np.searchsorted(
                r_sites_reverse_read, reverse_read.reference_end, side="left"
            )  # left

            # try : 
            distance_2 = np.abs(
                np.subtract(reverse_read.reference_end, r_sites_reverse_read[index])
            )

            # except : 

            #     distance_2 = np.abs(
            #         np.subtract(reverse_read.reference_end, r_sites_reverse_read[index - 1])
            #     )

        # Correction for uncuts with no restriction sites inside
        if forward_read.reference_name == reverse_read.reference_name and np.add(
            distance_1, distance_2
        ) > np.abs(np.subtract(reverse_read.pos, forward_read.pos)):
            list_d1d2.append(np.abs(np.subtract(reverse_read.pos, forward_read.pos)))

        else:

            list_d1d2.append(np.add(distance_1, distance_2))

    # print(f"max(list_d1d2) : {int(max(list_d1d2))}")
    # print(f"list_d1d2 : {list_d1d2}")


    histo, bins = np.histogram(list_d1d2, int(max(list_d1d2)))

    np.save(output_path / D1D2, histo)

    logger.info(f"Saved d1d2 law at : {output_path / D1D2}")

def get_patterns(forward_bam_file : str = "group1.1.bam", reverse_bam_file : str = "group1.2.bam", xs : str = "xs.npy", chrom_sizes : str = "chromosome_sizes.npy", circular : str = "", output_dir : str = None) -> None:
    """
    Get the patterns distribution from read pairs alignment. .

    Parameters
    ----------
    forward_bam_file : str, optional
        Path to forward .bam alignment file, by default None, by default group1.1.bam, by default "group1.1.bam", by default "group1.1.bam"
    reverse_bam_file : str, optional
        Path to reverse .bam alignment file, by default None, by default group1.1.bam, by default "group1.1.bam", by default "group1.2.bam"
    xs : str, optional
        Path to the dictionary containing the xs values, by default "xs.npy"
    dist_frag : str, optional
        Path to the dictionary containing the inter fragment distances, by default "dist.frag.npy"
    circular : str, optional
        Name of the chromosomes to consider as circular, by default ""
    output_dir : str, optional
        Path to the folder where to save the dictionary, by default None, by default None, by default None
    """    

    logger.info("Start generating patterns distribution...")

    if output_dir is None:
            
        output_path = Path(getcwd())

    else:

        output_path = Path(output_dir)

    forward_bam_path = Path(output_path, forward_bam_file)
    reverse_bam_path = Path(output_path, reverse_bam_file)

    if not forward_bam_path.is_file():
        
        raise FileNotFoundError(f"Forward .bam file {forward_bam_file} not found. Please provide a valid path to a forward .bam file.")
    
    if not reverse_bam_path.is_file():

        raise FileNotFoundError(f"Reverse .bam file {reverse_bam_file} not found. Please provide a valid path to a reverse .bam file.")
    
    #Load xs

    xs = hio.load_dictionary(output_path / XS)
    # dist_frag = hio.load_dictionary(output_path / dist_frag)
    chrom_size_dict = hio.load_dictionary(output_path / chrom_sizes)

    # Create placeholders for the dictionaries

    weirds = {seq_name : np.zeros(xs.get(seq_name).shape) for seq_name in xs.keys()}
    uncuts = {seq_name : np.zeros(xs.get(seq_name).shape) for seq_name in xs.keys()}
    loops = {seq_name : np.zeros(xs.get(seq_name).shape) for seq_name in xs.keys()}

    # Create placeholder for area to divide logbins counts
    trapezoids_area = {seq_name : np.zeros(xs.get(seq_name).shape) for seq_name in xs.keys()}

    # Compute areas of trapezoids

    # print("Computinf trapezoids areas...")
    for chrom in xs.keys():
        xs_ = xs[chrom]
        chrom_size_ = chrom_size_dict[chrom]

        trapezoids_area[chrom] = [(2 * chrom_size_ - xs_[j+1] - xs_[j]) * (xs_[j+1] - xs_[j]) * 0.5 for j in range(len(xs_) - 1)]
        trapezoids_area[chrom].append( (((chrom_size_ - xs_[-1]) ** 2) / 2))

    forward_bam_handler, reverse_bam_handler = pysam.AlignmentFile(forward_bam_path, "rb"), pysam.AlignmentFile(reverse_bam_path, "rb")

    for forward_read, reverse_read in zip(forward_bam_handler, reverse_bam_handler):

        

        if hut.is_intra_chromosome(forward_read, reverse_read):
                

            if hut.is_weird(forward_read, reverse_read):

                weirds[forward_read.reference_name][
                    attribute_xs(
                        xs.get(forward_read.reference_name),
                        hut.get_cis_distance(forward_read, reverse_read, circular) + 1,
                    )
                ] += 1

            if hut.is_uncut(forward_read, reverse_read):

                uncuts[forward_read.reference_name][
                    attribute_xs(
                        xs.get(forward_read.reference_name),
                        hut.get_cis_distance(forward_read, reverse_read, circular) + 1,
                    )
                ] += 1

            if hut.is_circle(forward_read, reverse_read):
                loops[forward_read.reference_name][
                    attribute_xs(
                        xs.get(forward_read.reference_name),
                        hut.get_cis_distance(forward_read, reverse_read, circular) + 1,
                    )
                ] += 1
    
    forward_bam_handler.close()
    reverse_bam_handler.close()

    for chromosome in weirds.keys():

        weirds[chromosome] = np.divide(
            weirds[chromosome],
            np.multiply(trapezoids_area.get(chromosome), 2),
            out=np.zeros_like(weirds[chromosome]),
            where=trapezoids_area.get(chromosome) != 0,
        )

        uncuts[chromosome] = np.divide(
            uncuts[chromosome],
            trapezoids_area.get(chromosome),
            out=np.zeros_like(uncuts[chromosome]),
            where=trapezoids_area.get(chromosome) != 0,
        )

        loops[chromosome] = np.divide(
            loops[chromosome],
            trapezoids_area.get(chromosome),
            out=np.zeros_like(loops[chromosome]),
            where=trapezoids_area.get(chromosome) != 0,
        )

    np.save(output_path / WEIRDS, weirds)
    np.save(output_path / UNCUTS, uncuts)
    np.save(output_path / LOOPS, loops)

    logger.info(f"Saved {WEIRDS}, {UNCUTS} and {LOOPS} in {output_path}")

def get_pair_ps(read_forward : pysam.AlignedSegment, read_reverse : pysam.AlignedSegment, xs : dict, weirds :  dict, uncuts : dict, loops : dict, circular : str = "") -> float:
    """
    Take two reads and return the P(s) value depending on event type (intra-chromosomal case only).

    Parameters
    ----------
    read_forward : pysam.AlignedSegment
        Forward read to compare with the reverse read.
    read_reverse : pysam.AlignedSegment
        Reverse read to compare with the forward read.
    xs : dict
        Dictionary containing log binning values for each chromosome.
    weirds : dict
        Dictionary containing number of weird events considering distance for each chromosome.
    uncuts : dict
        Dictionary containing number of uncuts events considering distance for each chromosome.
    loops : dict
        Dictionary containing number of loops events considering distance for each chromosome.
    circular : str, optional
        Name of the chromosomes to consider as circular, by default None, by default "".

    Returns
    -------
    float
        P(s) of the pair considering the event type.
    """    
    
    if read_forward.query_name != read_reverse.query_name:
        raise ValueError("Reads are not coming from the same pair.")

    if not hut.is_intra_chromosome(read_forward, read_reverse):
        raise ValueError("Reads are not intra-chromosomal.")
    
    if hut.is_weird(read_forward, read_reverse):
        return weirds[read_forward.reference_name][
            attribute_xs(
                xs[read_forward.reference_name],
                hut.get_cis_distance(read_forward, read_reverse, circular),
            )
        ]

    elif hut.is_uncut(read_forward, read_reverse):
        return uncuts[read_forward.reference_name][
            attribute_xs(
                xs[read_forward.reference_name],
                hut.get_cis_distance(read_forward, read_reverse, circular),
            )
        ]

    elif hut.is_circle(read_forward, read_reverse):
        return loops[read_forward.reference_name][
            attribute_xs(
                xs[read_forward.reference_name],
                hut.get_cis_distance(read_forward, read_reverse, circular),
            )
        ]
    

def get_trans_ps(read_forward : pysam.AlignedSegment, read_reverse : pysam.AlignedSegment, trans_ps : dict) -> float:
    """
    Parameters
    ----------
    read_forward : pysam.AlignedSegment
        Forward read to compare with the reverse read.
    read_reverse : pysam.AlignedSegment
        Reverse read to compare with the forward read.
    trans_ps : dict
        Dictionary of trans-chromosomal P(s)

    Returns
    -------
    float
        Trans P(s) value

    """    
    if read_forward.query_name != read_reverse.query_name:
        raise ValueError("Reads are not coming from the same pair.")
    
    if  hut.is_intra_chromosome(read_forward, read_reverse):

        raise ValueError("Reads are not inter-chromosomal.")
    
    return trans_ps[
        tuple(sorted([read_forward.reference_name, read_reverse.reference_name]))
    ] 

def get_pair_cover(read_forward : pysam.AlignedSegment, read_reverse : pysam.AlignedSegment, coverage : dict, bins : int) -> float:
    """
    Get the coverage of a pair of reads.

    Parameters
    ----------
    read_forward : pysam.AlignedSegment
        Forward read to compare with the reverse read.
    read_reverse : pysam.AlignedSegment
        Reverse read to compare with the forward read.
    coverage : dict
        Dictionary containing the coverage of each chromosome.
    bins : int
        Size of the desired bin.

    Returns
    -------
    float
        Product of the coverage of the pair of reads.
    """    

    if read_forward.query_name != read_reverse.query_name:
        raise ValueError("Reads are not coming from the same pair.")
    

    if (read_forward.flag == 16 or read_forward.flag == 272) and( read_reverse.flag == 0 or read_reverse.flag == 256):

        return (
            coverage[read_forward.reference_name][
                int(read_forward.reference_end / bins)
            ]
            * coverage[read_reverse.reference_name][
                int(read_reverse.reference_start / bins)
            ]
        )

    elif( read_forward.flag == 0 or read_forward.flag == 256) and (read_reverse.flag == 16 or read_reverse.flag == 272):
        return (
            coverage[read_forward.reference_name][
                int(read_forward.reference_start / bins)
            ]
            * coverage[read_reverse.reference_name][
                int(read_reverse.reference_end / bins)
            ]
        )

    elif (read_forward.flag == 16 or read_forward.flag == 272) and (read_reverse.flag == 16 or read_reverse.flag == 272):
        return (
            coverage[read_forward.reference_name][
                int(read_forward.reference_end / bins)
            ]
            * coverage[read_reverse.reference_name][
                int(read_reverse.reference_end / bins)
            ]
        )

    else:
        return (
            coverage[read_forward.reference_name][
                int(read_forward.reference_start / bins)
            ]
            * coverage[read_reverse.reference_name][
                int(read_reverse.reference_start / bins)
            ]
        )


def get_d1d2(read_forward : pysam.AlignedSegment, read_reverse : pysam.AlignedSegment, restriction_map : dict = None, d1d2 : np.array = None) -> int:
    """
    Get the d1d2 value of a pair of reads.

    Parameters
    ----------
    read_forward : pysam.AlignedSegment
        Forward read to compare with the reverse read.
    read_reverse : pysam.AlignedSegment
        Reverse read to compare with the forward read.
    restriction_map : dict, optional
        Restriction map saved as a dictionary like chrom_name : list of restriction sites' position, by default None
    d1d2 : np.array, optional
        Distribution of d1d2 values, by default None

    Returns
    -------
    int
        d1d2 value
    """

    if read_forward.query_name != read_reverse.query_name:
        raise ValueError("Reads are not coming from the same pair.")
    
    # Get appropriate restriction sites vecctors in dicitionary
    r_sites_read_for = restriction_map[read_forward.reference_name]
    r_sites_read_rev = restriction_map[read_reverse.reference_name]

    if read_forward.flag == 0 or read_forward.flag == 256:

        index = np.searchsorted(r_sites_read_for, read_forward.pos, side="right")
        distance_1 = np.subtract(r_sites_read_for[index], read_forward.pos)

    if read_forward.flag == 16 or read_forward.flag == 272:

        index = np.searchsorted(r_sites_read_for, read_forward.reference_end, side="left")
        distance_1 = np.abs(
            np.subtract(read_forward.reference_end, r_sites_read_for[index])
        )

    if read_reverse.flag == 0 or read_reverse.flag == 256:

        index = np.searchsorted(
            r_sites_read_rev, read_reverse.reference_start, side="right"
        )  # right
        distance_2 = np.subtract(r_sites_read_rev[index], read_reverse.reference_start)

    if read_reverse.flag == 16 or read_reverse.flag == 272:

        index = np.searchsorted(
            r_sites_read_rev, read_reverse.reference_end, side="left"
        )  # left
        distance_2 = np.abs(
            np.subtract(read_reverse.reference_end, r_sites_read_rev[index])
        )

    # Correction for uncuts with no restriction sites inside

    if read_forward.reference_name == read_reverse.reference_name and np.add(
        distance_1, distance_2
    ) > np.abs(np.subtract(read_reverse.pos, read_forward.pos)):
        distance = np.abs(np.subtract(read_reverse.pos, read_forward.pos))

    else:
        distance = np.add(distance_1, distance_2)
        
    return d1d2[distance]

def get_density(read_forward : pysam.AlignedSegment, read_reverse : pysam.AlignedSegment, density_map : dict[(str, str) : np.array], bin_size : int = 2000) -> float:
    """
    Get density from density map dictionary.

    Parameters
    ----------
    read_forward : pysam.AlignedSegment
        Forward read to get density from.
    read_reverse : pysam.AlignedSegment
        Reverse read to get density from.
    density_map : dict
        Dictionary containing density maps for each chromosome couple.
    bin_size : int, optional
        Resolution of the matrix on which density map has been computed, by default 2000

    Returns
    -------
    float
        Density corresponding to the pair of reads.

    """    
    if read_forward.query_name != read_reverse.query_name:
        raise ValueError("Reads are not coming from the same pair.")

    if hut.is_reverse(read_forward):
        position_for = int(read_forward.reference_end // bin_size)

    elif not hut.is_reverse(read_forward):
        position_for = int(read_forward.reference_start // bin_size)

    if hut.is_reverse(read_reverse):
        position_rev = int(read_reverse.reference_end // bin_size)

    elif not hut.is_reverse(read_reverse):

        position_rev = int(read_reverse.reference_start // bin_size)

    # print(f"read_forward : {read_forward.reference_name}")
    # print(f"position_for : {position_for}")
    # print(f"read_reverse : {read_reverse.reference_name}")
    # print(f"position_rev : {position_rev}")
    # print(f"map : {density_map.get((read_forward.reference_name, read_reverse.reference_name))}")
    
    
    couple_density = density_map.get((read_forward.reference_name, read_reverse.reference_name))[position_for, position_rev]

    # print(f"couple_density : {couple_density}")

    return couple_density


def compute_propensity(read_forward : pysam.AlignedSegment, read_reverse : pysam.AlignedSegment, restriction_map : dict = None, xs : dict = None, weirds : dict = None, uncuts : dict = None, loops : dict = None, circular : str = "", trans_ps : dict = None,  coverage : dict = None, bins : int = 2000, d1d2 : dict = None, density_map : dict = None,  mode : str = "full") -> float:
    """
    Compute propensity for read pair to be selected among all plausible pairs related to multi-mapping reads.

    Parameters
    ----------
    read_forward : pysam.AlignedSegment
        Forward read to compare with the reverse read.
    read_reverse : pysam.AlignedSegment
        Reverse read to compare with the forward read.
    restriction_map : dict, optional
        Restriction map saved as a dictionary like chrom_name : list of restriction sites' position, by default None
    xs : dict
        Dictionary containing log binning values for each chromosome.
    weirds : dict
        Dictionary containing number of weird events considering distance for each chromosome.
    uncuts : dict
        Dictionary containing number of uncuts events considering distance for each chromosome.
    loops : dict
        Dictionary containing number of loops events considering distance for each chromosome.
    circular : str, optional
        Name of the chromosomes to consider as circular, by default None, by default "".
    trans_ps : dict
        Dictionary of trans-chromosomal P(s)
    coverage : dict
        Dictionary containing the coverage of each chromosome.
    bins : int
        Size of the desired bin, by default 2000
    d1d2 : np.array, optional
        Distribution of d1d2 values, by default None
    density : dict
        Dictionary of contact density by chromosome couple.
    mode : str, optional
        Mode to use to compute propensity among, by default "full"

    Returns
    -------
    float
        Propensity to use for read couple drawing
    """

    if read_forward.query_name != read_reverse.query_name:
        raise ValueError("Reads are not coming from the same pair.")
    
    if mode == "ps_only":

        if hut.is_intra_chromosome(read_forward, read_reverse):
    
            ps = get_pair_ps(
                read_forward,
                read_reverse,
                xs,
                weirds,
                uncuts,
                loops,
                circular,
            )

        else:
            
            ps = get_trans_ps(read_forward, read_reverse, trans_ps)

            # Avoid ps = 0 making the read unselectable. Value of 1 make the propensity unsensitive to P(s).
            if ps == 0:

                ps = 1

        return ps
    
    elif mode == "cover_only":

        cover = get_pair_cover(read_forward, read_reverse, coverage, bins=bins)

        # Avoid cover = 0 making the read unselectable. Value of 1 make the propensity unsensitive to coverage.
        if cover <= 0:
            cover = 1

        return cover
    
    elif mode == "d1d2_only":  # not functional so far

        try:
            d1d2 = get_d1d2(
            read_forward,
            read_reverse,
            restriction_map,
            d1d2,

        )

        except:

            d1d2 = 1

        return d1d2

    elif mode == "density_only":

        density = get_density(read_forward, read_reverse, density_map = density_map)

        return density


    elif mode == "no_ps":

        cover = get_pair_cover(read_forward, read_reverse, coverage, bins=bins)

        # Avoid cover = 0 making the read unselectable. Value of 1 make the propensity unsensitive to coverage.
        if cover <= 0:
            cover = 1

        try:
            d1d2 = get_d1d2(
            read_forward,
            read_reverse,
            restriction_map,
            d1d2,

        )

        except:

            d1d2 = 1

        density = get_density(read_forward, read_reverse, density_map = density_map)

        return cover * d1d2 * density
    
    elif mode == "no_cover":

        if hut.is_intra_chromosome(read_forward, read_reverse):
        
            ps = get_pair_ps(
                read_forward,
                read_reverse,
                xs,
                weirds,
                uncuts,
                loops,
                circular,
            )

        else:
            
            ps = get_trans_ps(read_forward, read_reverse, trans_ps)

            # Avoid ps = 0 making the read unselectable. Value of 1 make the propensity unsensitive to P(s).
            if ps == 0:

                ps = 1

        try:
            d1d2 = get_d1d2(
            read_forward,
            read_reverse,
            restriction_map,
            d1d2,
        )

        except:

            d1d2 = 1

        density = get_density(read_forward, read_reverse, density_map = density_map)

        return ps * d1d2 * density

    elif mode == "no_d1d2":

        if hut.is_intra_chromosome(read_forward, read_reverse):
        
            ps = get_pair_ps(
                read_forward,
                read_reverse,
                xs,
                weirds,
                uncuts,
                loops,
                circular,
            )

        else:
            
            ps = get_trans_ps(read_forward, read_reverse, trans_ps)

            # Avoid ps = 0 making the read unselectable. Value of 1 make the propensity unsensitive to P(s).
            if ps == 0:

                ps = 1

        cover = get_pair_cover(read_forward, read_reverse, coverage, bins = bins)

        # Avoid cover = 0 making the read unselectable. Value of 1 make the propensity unsensitive to coverage.
        if cover <= 0:
            cover = 1

        density = get_density(read_forward, read_reverse, density_map = density_map)

        return ps * cover * density
    
    elif mode == "no_density":

        if hut.is_intra_chromosome(read_forward, read_reverse):
        
            ps = get_pair_ps(
                read_forward,
                read_reverse,
                xs,
                weirds,
                uncuts,
                loops,
                circular,
            )

        else:
            
            ps = get_trans_ps(read_forward, read_reverse, trans_ps)

            # Avoid ps = 0 making the read unselectable. Value of 1 make the propensity unsensitive to P(s).
            if ps == 0:

                ps = 1

        try:
            d1d2 = get_d1d2(
            read_forward,
            read_reverse,
            restriction_map,
            d1d2,
        )

        except:

            d1d2 = 1


        cover = get_pair_cover(read_forward, read_reverse, coverage, bins=bins)

        # Avoid cover = 0 making the read unselectable. Value of 1 make the propensity unsensitive to coverage.
        if cover <= 0:
            cover = 1

        return ps * cover * d1d2
    
    elif mode == "no_cover_no_density":

        if hut.is_intra_chromosome(read_forward, read_reverse):
        
            ps = get_pair_ps(
                read_forward,
                read_reverse,
                xs,
                weirds,
                uncuts,
                loops,
                circular,
            )

        else:
            
            ps = get_trans_ps(read_forward, read_reverse, trans_ps)

            # Avoid ps = 0 making the read unselectable. Value of 1 make the propensity unsensitive to P(s).
            if ps == 0:

                ps = 1

        try:
            d1d2 = get_d1d2(
            read_forward,
            read_reverse,
            restriction_map,
            d1d2,
        )

        except:

            d1d2 = 1


        return ps * d1d2

        
    
    elif mode == "no_ps_no_d1d2":

        cover = get_pair_cover(read_forward, read_reverse, coverage, bins=bins)

        # Avoid cover = 0 making the read unselectable. Value of 1 make the propensity unsensitive to coverage.
        if cover <= 0:
            cover = 1

        density = get_density(read_forward, read_reverse, density_map = density_map)

        return cover * density

    elif mode == "no_ps_no_cover":

        try:
            d1d2 = get_d1d2(
            read_forward,
            read_reverse,
            restriction_map,
            d1d2,
        )

        except:

            d1d2 = 1

        density = get_density(read_forward, read_reverse, density_map = density_map)

        return d1d2 * density
    
    elif mode == "no_d1d2_no_cover":

        if hut.is_intra_chromosome(read_forward, read_reverse):
        
            ps = get_pair_ps(
                read_forward,
                read_reverse,
                xs,
                weirds,
                uncuts,
                loops,
                circular,
            )

        else:
            
            ps = get_trans_ps(read_forward, read_reverse, trans_ps)

            # Avoid ps = 0 making the read unselectable. Value of 1 make the propensity unsensitive to P(s).
            if ps == 0:

                ps = 1

        density = get_density(read_forward, read_reverse, density_map = density_map)

        return ps * density
    
    elif mode == "no_ps_no_density":

        cover = get_pair_cover(read_forward, read_reverse, coverage, bins=bins)

        try : 

            d1d2 = get_d1d2(
            read_forward,
            read_reverse,
            restriction_map,
            d1d2,
        )
            
        
        except :

            d1d2 = 1

        return cover * d1d2
    
    elif mode == "no_density_no_d1d2":

        if hut.is_intra_chromosome(read_forward, read_reverse):
        
            ps = get_pair_ps(
                read_forward,
                read_reverse,
                xs,
                weirds,
                uncuts,
                loops,
                circular,
            )

        else:
            
            ps = get_trans_ps(read_forward, read_reverse, trans_ps)

        cover = get_pair_cover(read_forward, read_reverse, coverage, bins=bins)

        return ps * cover

    elif mode == "random":

        return 1

    elif mode == "full":

        if hut.is_intra_chromosome(read_forward, read_reverse):
        
            ps = get_pair_ps(
                read_forward,
                read_reverse,
                xs,
                weirds,
                uncuts,
                loops,
                circular,
            )

        else:
            
            ps = get_trans_ps(read_forward, read_reverse, trans_ps)


            # Avoid ps = 0 making the read unselectable. Value of 1 make the propensity unsensitive to P(s).
            if ps == 0:

                ps = 1

        cover = get_pair_cover(read_forward, read_reverse, coverage, bins=bins)

        # Avoid cover = 0 making the read unselectable. Value of 1 make the propensity unsensitive to coverage.
        if cover <= 0:
            cover = 1

        try:
            d1d2 = get_d1d2(
            read_forward,
            read_reverse,
            restriction_map,
            d1d2,
        )

        except:

            d1d2 = 1

        density = get_density(read_forward, read_reverse, density_map = density_map)
    

        return ps * d1d2 * cover * density
    

def draw_read_couple(propensities : np.array) -> int:
    """
    Draw an index respecting distribution of propensities. This function is used to draw a couple of reads considering the propensity of each couple.

    Parameters
    ----------
    propensities : np.array
        Array containing all the propensities of each couple of reads.

    Returns
    -------
    int
        Index of the couple of reads drawn.
    """

    xk = np.arange(len(propensities))

    if  np.sum(propensities) > 0: 

        pk = np.divide(propensities, np.sum(propensities))

    elif np.sum(propensities) <= 0:

        pk = np.full(xk.shape, np.divide(1, len(propensities)))

    try:
        index = choice(xk, p=pk)

    # TODO : enventually remove
    except:
        
        
        print("-- draw_read_couple --")
        print(f"propensities : {propensities}")
        print(f"pk : {pk}")
        print(f"xk : {xk}")
        print(f"min pk : {np.min(pk)}")
        print(f"min propensity: {np.min(propensities)}")
        a
        

    return index

def reattribute_reads(reads_couple : tuple[str, str] = ("group2.1.bam", "group2.2.bam"), restriction_map : dict = None, xs : dict = "xs.npy", weirds : dict = "weirds.npy", uncuts : dict = "uncuts.npy", loops : dict = "loops.npy", circular : str = "", trans_ps : dict = "trans_ps.npy",  coverage : dict = "coverage.npy", bins : int = 2000, d1d2 : dict = "d1d2.npy", density_map : dict = "density_map.npy",  mode : str = "full", output_dir : str = None) -> None:
    """
    Re-attribute multi-mapping (ambiguous) reads considering sets of statistical laws.

    Parameters
    ----------
    reads_couple : tuple[str, str], optional
        Paths to ambiguous reads alignment files (.bam), by default ("group2.1.bam", "group2.2.bam")
    restriction_map : dict, optional
        Restriction map saved as a dictionary like chrom_name : list of restriction sites' position, by default None
    xs : dict
        Dictionary containing log binning values for each chromosome.
    weirds : dict
        Dictionary containing number of weird events considering distance for each chromosome.
    uncuts : dict
        Dictionary containing number of uncuts events considering distance for each chromosome.
    loops : dict
        Dictionary containing number of loops events considering distance for each chromosome.
    circular : str, optional
        Name of the chromosomes to consider as circular, by default None, by default "".
    trans_ps : dict
        Dictionary of trans-chromosomal P(s)
    coverage : dict
        Dictionary containing the coverage of each chromosome.
    bins : int
        Size of the desired bin, by default 2000
    d1d2 : np.array, optional
        Distribution of d1d2 values, by default None
    density : np.array, optional
        Dictionary containing density maps per chromosome couples as, by default None
    mode : str, optional
        Mode to use to compute propensity among, by default "full"
    output_dir : str, optional
        Path to the re-attributed ambiguous reads alignment files are saved, by default None
    """ 

    if output_dir is None:
            
        output_path = Path(getcwd())

    else:
            
        output_path = Path(output_dir)

    #TODO : Add selective loading of dictionaries depending on reconstruction mode
        
    

    #Reload dictionaries
    xs = hio.load_dictionary(output_path / xs)
    weirds  = hio.load_dictionary(output_path / weirds)
    uncuts = hio.load_dictionary(output_path / uncuts)
    loops = hio.load_dictionary(output_path / loops)
    trans_ps = hio.load_dictionary(output_path / trans_ps)
    coverage = hio.load_dictionary(output_path / coverage)
    d1d2 = None
    density = None

    if mode == "full" :

        d1d2 = hio.load_dictionary(output_path / "d1d2.npy")
        density = hio.load_dictionary(output_path / "density_map.npy")
    
    elif mode in ["d1d2_only", "no_ps", "no_cover", "no_density", "no_cover_no_density", "no_ps_no_cover", "no_ps_no_density"]: # TODO : to be completed

        d1d2 = hio.load_dictionary(output_path / "d1d2.npy")
        
    elif mode in ["no_d1d2_no_cover", "no_ps_no_cover",  "no_ps_no_d1d2", "no_d1d2", "no_cover", "no_ps", "density_only"]:
            
        density = hio.load_dictionary(output_path / "density_map.npy")

    forward_bam_path, reverse_bam_path = Path(reads_couple[0]), Path(reads_couple[1])
    file_id = time.time()
    id_for = uuid.uuid4()
    id_rev = uuid.uuid4()

    forward_bam_handler = pysam.AlignmentFile(forward_bam_path, "rb")
    reverse_bam_handler = pysam.AlignmentFile(reverse_bam_path, "rb")

    forward_out_bam_handler = pysam.AlignmentFile(output_path / f"forward_{id_for}_{file_id}_predicted.bam", "wb", template=forward_bam_handler)
    reverse_out_bam_handler = pysam.AlignmentFile(output_path / f"reverse_{id_rev}_{file_id}_predicted.bam", "wb", template=reverse_bam_handler)

    # Instanciate generators
    forward_generator = hut.bam_iterator(forward_bam_path)
    reverse_generator = hut.bam_iterator(reverse_bam_path)

    for forward_block, reverse_block in zip(forward_generator, reverse_generator):

        propensities = []

        combinations = list(itertools.product(tuple(forward_block), tuple(reverse_block)))

        for combination in combinations:

            # # TODO : enventually remove
            # if compute_propensity(read_forward = combination[0], read_reverse = combination[1], restriction_map = restriction_map, xs = xs, weirds = weirds, uncuts = uncuts, loops = loops, trans_ps = trans_ps, coverage = coverage, bins = bins, d1d2 = d1d2, density_map = density, mode = mode) < 0:
            #     pass

            propensities.append(compute_propensity(read_forward = combination[0], read_reverse = combination[1], restriction_map = restriction_map, xs = xs, weirds = weirds, uncuts = uncuts, loops = loops, trans_ps = trans_ps, coverage = coverage, bins = bins, d1d2 = d1d2, density_map = density, mode = mode))

        selected_couple_index = draw_read_couple(propensities)

        selected_read_forward, selected_read_reverse = combinations[selected_couple_index]
        selected_read_forward.set_tag("XL", len(forward_block))
        selected_read_reverse.set_tag("XL", len(reverse_block))

        forward_out_bam_handler.write(selected_read_forward)
        reverse_out_bam_handler.write(selected_read_reverse)
        
    
    forward_bam_handler.close()
    reverse_bam_handler.close()

    logger.info(f"Predictions written in {output_path}")


def pearson_score(original_matrix, rescued_matrix , markers):

    ori_matrix = original_matrix.matrix(balance=False)[:]
    reco_matrix = rescued_matrix.matrix(balance=False)[:]

    ori_vector = ori_matrix[markers]
    reco_vector = reco_matrix[markers]

    pearson_score = pearsonr(ori_vector.flatten(), reco_vector.flatten())

    return pearson_score[0]

# Benchamrk analysis functions
def get_top_pattern(file : str = None, top : int = 10, chromosome : str = None) -> pd.DataFrame:
    """
    Get top patterns from a dataframe

    Parameters
    ----------
    df : pd.DataFrame, optional
        Dataframe containing patterns given by Chromosight, by default None
    top : int, optional
        Percentage of top patterns to get, by default 10
    chromosome : str, optional
        Chromosome to consider, by default None

    Returns
    -------
    pd.DataFrame
        Dataframe containing top percentage patterns.
    """
    df = pd.read_csv(file, sep = "\t", header = 0)
    top_factor = (df.shape[0] * top) // 100

    if chromosome is not None:

        df = df.query(f"chrom1 == '{chromosome}' and chrom2 == '{chromosome}'")
    df_top = df.sort_values(by='score', ascending=False).head(top_factor).reset_index(drop=True)

    return df_top

