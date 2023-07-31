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

from .conftest import temporary_folder
from.test_utils import test_classify_reads, test_get_chromosomes_sizes, test_get_bin_table, test_chunk_bam
from .test_align import test_hic_build_index, test_hic_align, test_hic_view, test_hic_sort
from.test_io import test_build_matrix, test_build_pairs


import hicberg.statistics as hst
import hicberg.utils as hut
import hicberg.io as hio


lowess = sm.nonparametric.lowess

GENOME = "data_test/SC288_with_micron.fa"
XS = "xs.npy"
RESTRICTION_DICO = "dist.frag.npy"
COVERAGE_DICO = "coverage.npy"
D1D2 = "d1d2.npy"
UNCUTS = "uncuts.npy"
WEIRDS = "weirds.npy"
LOOPS = "loops.npy"
TRANS_PS = "trans_ps.npy"

DICT_FIRST_KEY = "chr10"
DICT_FIRST_CHR_LAST_POS = 745751
DICT_SECOND_KEY = "chr11"
DICT_SECOND_CHR_LAST_POS = 666816

HEADER = pysam.AlignmentHeader().from_dict({
            "HD": {"VN": "1.0", "SO": "unsorted"},
            "SQ": [
                {"SN": DICT_FIRST_KEY, "LN": DICT_FIRST_CHR_LAST_POS},
                {"SN": DICT_SECOND_KEY, "LN": DICT_SECOND_CHR_LAST_POS},
            ],
        })


BINS = 2000

DEFAULT_ENZYME = ["DpnII"]
DICT_FIRST_KEY = "chr10"
DICT_FIRST_CHR_FIRST_POS = 0
DICT_FIRST_CHR_LAST_POS = 745751
LENGTH_DPNII = 2266
LENGTH_DPNII_HINFI = 4680
SUBSAMPLE_RATE = 1.0

CIRCULAR = "chrM"

MODE = "full"

XS_BASE = 1.1
LENGTH_XS_CHR_FIRST = 127
XS_CHR_FIRST_VALUE = 1
XS_CHR_LAST_VALUE = 754677

DISTANCE = 100000
DISTANCE_XS_IDX = 104

PS_VALUE = 3.0
D1D2_VALUE = 18

@pytest.fixture(scope = "module")
@pytest.mark.parametrize("genome, enzyme, length", [(GENOME, DEFAULT_ENZYME, LENGTH_DPNII), (GENOME, ["DpnII", "HinfI"], LENGTH_DPNII_HINFI)])
def test_get_restriction_map(genome, enzyme, length):
    """
    Test if the function returns a dictionary with the right keys and values, i.e. the genome has been correctly digested.
    First restriction site for first chromosome should be 0 and last restriction site should be the length of the chromosome.
    """

    genome_path  = Path(genome)
    restriction_map = hst.get_restriction_map(genome = genome_path, enzyme = enzyme)

    # Check if the dictionary is not empty, if fisrt chromosome is right, and restrictions sites are in the right order.
    assert restriction_map[DICT_FIRST_KEY][0] == DICT_FIRST_CHR_FIRST_POS
    assert restriction_map[DICT_FIRST_KEY][-1] == DICT_FIRST_CHR_LAST_POS
    #Check if the length of the restriction map is right.
    assert len(restriction_map[DICT_FIRST_KEY]) == length
    


@pytest.fixture(scope = "module")
def test_get_restriction_map_mono():
    """
    Test if the function returns a dictionary with the right keys and values, i.e. the genome has been correctly digested.
    First restriction site for first chromosome should be 0 and last restriction site should be the length of the chromosome.
    """

    genome_path  = Path(GENOME)
    restriction_map = hst.get_restriction_map(genome = genome_path, enzyme = DEFAULT_ENZYME)

    yield restriction_map

    # Check if the dictionary is not empty, if fisrt chromosome is right, and restrictions sites are in the right order.
    assert restriction_map[DICT_FIRST_KEY][0] == DICT_FIRST_CHR_FIRST_POS
    assert restriction_map[DICT_FIRST_KEY][-1] == DICT_FIRST_CHR_LAST_POS
    #Check if the length of the restriction map is right.
    assert len(restriction_map[DICT_FIRST_KEY]) == LENGTH_DPNII


@pytest.fixture(scope = "module")
@pytest.mark.parametrize("genome, enzyme, length, rate", [(GENOME, DEFAULT_ENZYME, LENGTH_DPNII, 1.0), (GENOME, ["DpnII", "HinfI"], LENGTH_DPNII_HINFI, 1.0), (GENOME, DEFAULT_ENZYME, LENGTH_DPNII, 0.5), (GENOME, ["DpnII", "HinfI"], LENGTH_DPNII_HINFI, 0.5), (GENOME, DEFAULT_ENZYME, LENGTH_DPNII, 0.2), (GENOME, ["DpnII", "HinfI"], LENGTH_DPNII_HINFI, 0.2)])
def test_subsample_restriction_map(genome, enzyme, length, rate):
    """
    Test the subsample restriction map against the given enzyme and rate .
    """    

    genome_path  = Path(genome)
    restriction_map = hst.get_restriction_map(genome = genome_path, enzyme = enzyme)

    subsampled_restriction_map = hut.subsample_restriction_map(restriction_map = restriction_map, rate = rate)

    assert subsampled_restriction_map[DICT_FIRST_KEY][0] == DICT_FIRST_CHR_FIRST_POS
    assert subsampled_restriction_map[DICT_FIRST_KEY][-1] == DICT_FIRST_CHR_LAST_POS
    #Check if the length of the restriction map is right. (+/-1 because of the random subsampling on even or odd numbers)
    assert float(len(subsampled_restriction_map[DICT_FIRST_KEY])) -1 <= length *  rate <= float(len(subsampled_restriction_map[DICT_FIRST_KEY])) + 1

@pytest.fixture(scope = "module")
def test_generate_xs(temporary_folder):
    """
    Test if the log binning of the genome is correctly done .
    """

    temp_dir_path = Path(temporary_folder)
    genome_path = Path(GENOME)

    for seq_record in SeqIO.parse(genome_path, "fasta"):

        seq_name = seq_record.id

        if seq_name == DICT_FIRST_KEY:

            seq_len = len(seq_record.seq)
            break

    xs = hst.generate_xs(chromosome_size = seq_len, base = XS_BASE)


    yield xs

    # Check if the length of the list is right.
    assert len(xs) == LENGTH_XS_CHR_FIRST
    # Check if the first and last values are right.
    assert xs[0] == XS_CHR_FIRST_VALUE
    assert xs[-1] == XS_CHR_LAST_VALUE

@pytest.fixture(scope = "module")
def test_log_bin_genome(temporary_folder):
    """
    Test if the log binning of the genome is correctly done .
    """
    temp_dir_path = Path(temporary_folder)
    genome_path = Path(GENOME)
    xs = hst.log_bin_genome(genome = genome_path, base = XS_BASE, output_dir = temp_dir_path)

    xs_path = temp_dir_path / XS

    yield xs_path

    assert xs_path.is_file()


def test_attribute_xs(test_generate_xs):
    """
    Test if attribution of a distance to a corresponding log binned genomic segment is correctly done.
    """

    xs = test_generate_xs
    idx = hst.attribute_xs(xs = xs, distance = DISTANCE)

    assert idx == DISTANCE_XS_IDX 


@pytest.fixture(scope = "module")
def test_get_dist_frags(temporary_folder, test_get_restriction_map_mono):
    """
    Test if the distance between restriction sites is correctly computed.
    """
    temp_dir_path = Path(temporary_folder)
    print(f"temp_dir_path : {temp_dir_path}")
    
    hst.get_dist_frags(genome = GENOME, restriction_map = test_get_restriction_map_mono, circular = CIRCULAR, rate = SUBSAMPLE_RATE, output_dir = temp_dir_path)

    chrom_sizes_dictionary_path = temp_dir_path / RESTRICTION_DICO
    
    assert chrom_sizes_dictionary_path.is_file()

    yield chrom_sizes_dictionary_path

@pytest.fixture(scope = "module")
def test_generate_trans_ps(temporary_folder, test_get_restriction_map_mono, test_build_matrix):
    """
    Test if the transchromosomal P(s) are correctly generated.
    """
    temp_dir_path = Path(temporary_folder)

    print(f"Build matrix : {test_build_matrix}")
    print(f"File name : {test_build_matrix.name}")

    hst.generate_trans_ps(matrix = test_build_matrix, restriction_map = test_get_restriction_map_mono, output_dir = temp_dir_path)
    
    trans_ps_path = temp_dir_path / TRANS_PS

    yield trans_ps_path
    
    assert trans_ps_path.is_file()

def test_generate_probabilities():
    pass

def test_filter_poor_covered():
    pass

@pytest.fixture(scope = "module")
def test_generate_coverages(temporary_folder, test_classify_reads):
    """
    Test if the coverages are correctly generated.
    """
    temp_dir_path = Path(temporary_folder)
    hst.generate_coverages(genome = GENOME, forward_bam_file = test_classify_reads[0], reverse_bam_file = test_classify_reads[1], bins = BINS, output_dir = temp_dir_path)
    coverage_dictionary_path = temp_dir_path / COVERAGE_DICO

    yield coverage_dictionary_path
    
    assert coverage_dictionary_path.is_file()

@pytest.fixture(scope = "module")
def test_generate_d1d2(temporary_folder, test_classify_reads, test_get_restriction_map_mono):
    """
    Test if the d1d2 law is correctly generated.
    """
    temp_dir_path = Path(temporary_folder)
    hst.generate_d1d2(forward_bam_file = test_classify_reads[0], reverse_bam_file = test_classify_reads[1], restriction_map = test_get_restriction_map_mono, output_dir = temp_dir_path)
    d1d2_dictionary_path = temp_dir_path / D1D2

    yield d1d2_dictionary_path

    assert d1d2_dictionary_path.is_file()

# TODO : Add argument dist_frag

@pytest.fixture(scope = "module")
def test_get_patterns(temporary_folder, test_classify_reads, test_log_bin_genome, test_get_dist_frags):
    """
    Test if the patterns are correctly generated.
    """
    temp_dir_path = Path(temporary_folder)
    hst.get_patterns(forward_bam_file = test_classify_reads[0], reverse_bam_file = test_classify_reads[1], xs = test_log_bin_genome, output_dir = temp_dir_path)
    weirds_dictionary_path = temp_dir_path / WEIRDS
    uncuts_dictionary_path = temp_dir_path / UNCUTS
    loops_dictionary_path = temp_dir_path / LOOPS

    yield weirds_dictionary_path, uncuts_dictionary_path, loops_dictionary_path

    assert weirds_dictionary_path.is_file()
    assert uncuts_dictionary_path.is_file()
    assert loops_dictionary_path.is_file()

def test_get_pair_ps(temporary_folder, test_log_bin_genome, test_get_patterns):
    """
    Test if the get_pair_ps function is correctly returning the right P(s) regarding a pair of reads.
    """

    read_forward = pysam.AlignedSegment(header = HEADER)
    read_forward.query_name = "TEST"
    read_forward.reference_name = DICT_FIRST_KEY
    read_forward.reference_start = 100
    read_forward.cigarstring = "100M"
    read_forward.flag = 16

    read_reverse = pysam.AlignedSegment(header = HEADER)
    read_reverse.query_name = "TEST"
    read_reverse.reference_name = DICT_FIRST_KEY
    read_reverse.reference_start = 500
    read_reverse.cigarstring = "100M"
    read_reverse.flag = 0


    weirds_dictionary_path, uncuts_dictionary_path, loops_dictionary_path = test_get_patterns

    print(f"test_log_bin_genome : {test_log_bin_genome}")
    print(f"weirds_dictionary_path : {weirds_dictionary_path}")
    print(f"uncuts_dictionary_path : {uncuts_dictionary_path}")
    print(f"loops_dictionary_path : {loops_dictionary_path}")
    
    xs = hio.load_dictionary(test_log_bin_genome)
    weirds = hio.load_dictionary(weirds_dictionary_path)
    uncuts = hio.load_dictionary(uncuts_dictionary_path)
    loops = hio.load_dictionary(loops_dictionary_path)

    ps = hst.get_pair_ps(read_forward = read_forward, read_reverse = read_reverse, xs = xs, weirds = weirds, uncuts = uncuts, loops = loops)

    # assert ps == PS_VALUE

    # TODO : To be checked
    assert ps >= 0

def test_get_trans_ps(test_generate_trans_ps):
    """
    Test if the get_trans_ps function returns the right transchromosomal P(s) regarding a pair of reads.
    """

    read_forward = pysam.AlignedSegment(header = HEADER)
    read_forward.query_name = "TEST"
    read_forward.reference_name = DICT_FIRST_KEY
    read_forward.reference_start = 100
    read_forward.cigarstring = "100M"
    read_forward.flag = 16

    read_reverse = pysam.AlignedSegment(header = HEADER)
    read_reverse.query_name = "TEST"
    read_reverse.reference_name = DICT_SECOND_KEY
    read_reverse.reference_start = 10000
    read_reverse.cigarstring = "100M"
    read_reverse.flag = 0

    trans_ps = hio.load_dictionary(test_generate_trans_ps)

    tps = hst.get_trans_ps(read_forward = read_forward, read_reverse = read_reverse, trans_ps = trans_ps)



    assert tps <=  1.0

def test_get_pair_cover(test_generate_coverages):
    """
    Test if the get_pair_cover function returns the right coverage regarding a pair of reads.
    """

    read_forward = pysam.AlignedSegment(header = HEADER)
    read_forward.query_name = "TEST"
    read_forward.reference_name = DICT_FIRST_KEY
    read_forward.reference_start = 100
    read_forward.cigarstring = "100M"
    read_forward.flag = 16

    read_reverse = pysam.AlignedSegment(header = HEADER)
    read_reverse.query_name = "TEST"
    read_reverse.reference_name = DICT_SECOND_KEY
    read_reverse.reference_start = 10000
    read_reverse.cigarstring = "100M"
    read_reverse.flag = 0

    coverages = hio.load_dictionary(test_generate_coverages)

    cover = hst.get_pair_cover(read_forward = read_forward, read_reverse = read_reverse, coverage = coverages,  bins = BINS)

    assert cover > 0

def test_get_d1d2(test_get_restriction_map_mono, test_generate_d1d2):
    """
    Test if the get_d1d2 function returns the right d1d2 regarding a pair of reads.
    """

    read_forward = pysam.AlignedSegment(header = HEADER)
    read_forward.query_name = "TEST"
    read_forward.reference_name = DICT_FIRST_KEY
    read_forward.reference_start = 100
    read_forward.cigarstring = "100M"
    read_forward.flag = 16

    read_reverse = pysam.AlignedSegment(header = HEADER)
    read_reverse.query_name = "TEST"
    read_reverse.reference_name = DICT_SECOND_KEY
    read_reverse.reference_start = 10000
    read_reverse.cigarstring = "100M"
    read_reverse.flag = 0


    d1d2_law = hio.load_dictionary(test_generate_d1d2)

    print(f"TEST D1D2 LAW: {d1d2_law}")
    d1d2 = hst.get_d1d2(read_forward = read_forward, read_reverse = read_reverse, restriction_map = test_get_restriction_map_mono, d1d2 = d1d2_law)

    print(f"TEST D1D2: {d1d2}")


    assert d1d2 == D1D2_VALUE

def test_get_d1d2_distance():
    pass


def test_compute_propentsity(test_get_restriction_map_mono, test_log_bin_genome, test_get_patterns, test_generate_trans_ps, test_generate_coverages, test_generate_d1d2):
    """
    Test if the compute_propensity function returns the right propensity regarding a pair of reads.
    """
    read_forward = pysam.AlignedSegment(header = HEADER)
    read_forward.query_name = "TEST"
    read_forward.reference_name = DICT_FIRST_KEY
    read_forward.reference_start = 100
    read_forward.cigarstring = "100M"
    read_forward.flag = 16

    read_reverse = pysam.AlignedSegment(header = HEADER)
    read_reverse.query_name = "TEST"
    read_reverse.reference_name = DICT_SECOND_KEY
    read_reverse.reference_start = 10000
    read_reverse.cigarstring = "100M"
    read_reverse.flag = 0

    weirds_dictionary_path, uncuts_dictionary_path, loops_dictionary_path = test_get_patterns
    
    xs = hio.load_dictionary(test_log_bin_genome)
    weirds = hio.load_dictionary(weirds_dictionary_path)
    uncuts = hio.load_dictionary(uncuts_dictionary_path)
    loops = hio.load_dictionary(loops_dictionary_path)

    trans_ps = hio.load_dictionary(test_generate_trans_ps)
    coverage = hio.load_dictionary(test_generate_coverages)
    d1d2 = hio.load_dictionary(test_generate_d1d2)

    restriction_map = test_get_restriction_map_mono
    propensity = hst.compute_propentsity(read_forward = read_forward, read_reverse = read_reverse, restriction_map = restriction_map, xs = xs, weirds = weirds, uncuts = uncuts, loops = loops, trans_ps = trans_ps, coverage = coverage, bins = BINS, d1d2 = d1d2, mode = MODE)

    assert propensity >= 0

def test_decompose_propentsity():
    pass

def test_check_propensity():
    pass

def test_draw_read_couple(test_get_restriction_map_mono, test_log_bin_genome, test_get_patterns, test_generate_trans_ps, test_generate_coverages, test_generate_d1d2):
    """
    Test if the draw_read_couple function returns the right read couple index.
    """
    
    read_forward = pysam.AlignedSegment(header = HEADER)
    read_forward.query_name = "TEST"
    read_forward.reference_name = DICT_FIRST_KEY
    read_forward.reference_start = 100
    read_forward.cigarstring = "100M"
    read_forward.flag = 16

    read_reverse = pysam.AlignedSegment(header = HEADER)
    read_reverse.query_name = "TEST"
    read_reverse.reference_name = DICT_SECOND_KEY
    read_reverse.reference_start = 10000
    read_reverse.cigarstring = "100M"
    read_reverse.flag = 0

    read_reverse_2 = pysam.AlignedSegment(header = HEADER)
    read_reverse_2.query_name = "TEST"
    read_reverse_2.reference_name = DICT_SECOND_KEY
    read_reverse_2.reference_start = 50000
    read_reverse_2.cigarstring = "100M"
    read_reverse_2.flag = 0

    weirds_dictionary_path, uncuts_dictionary_path, loops_dictionary_path = test_get_patterns
    
    xs = hio.load_dictionary(test_log_bin_genome)
    weirds = hio.load_dictionary(weirds_dictionary_path)
    uncuts = hio.load_dictionary(uncuts_dictionary_path)
    loops = hio.load_dictionary(loops_dictionary_path)

    trans_ps = hio.load_dictionary(test_generate_trans_ps)
    coverage = hio.load_dictionary(test_generate_coverages)
    d1d2 = hio.load_dictionary(test_generate_d1d2)

    restriction_map = test_get_restriction_map_mono

    reads_forward = [read_forward]
    reads_reverse = [read_reverse, read_reverse_2]

    propensities = []

    all_couples = list(
            itertools.product(tuple(reads_forward), tuple(reads_reverse))
        )
    
    for couple in all_couples:
        
        propensity =  hst.compute_propentsity(read_forward = couple[0], read_reverse = couple[1], restriction_map = restriction_map, xs = xs, weirds = weirds, uncuts = uncuts, loops = loops, trans_ps = trans_ps, coverage = coverage, bins = BINS, d1d2 = d1d2, mode = MODE)
        propensities.append(propensity)

    read_couple_index = hst.draw_read_couple(propensities = propensities)

    assert read_couple_index == 1 or read_couple_index == 0




# @pytest.fixture(scope="session")
def test_reattribute_reads(temporary_folder, test_classify_reads, test_get_restriction_map_mono, test_log_bin_genome, test_get_patterns, test_generate_trans_ps, test_generate_coverages, test_generate_d1d2):
    """
    Test if the reattribute_reads function create prediction alignment files.
    """

    temp_dir_path = Path(temporary_folder)

    weirds_dictionary_path, uncuts_dictionary_path, loops_dictionary_path = test_get_patterns

    print(f"test_log_bin_genome: {test_log_bin_genome}")
    print(f"test_get_restriction_map_mono : {test_get_restriction_map_mono}")

    xs = test_log_bin_genome
    weirds = weirds_dictionary_path
    uncuts = uncuts_dictionary_path
    loops = loops_dictionary_path

    trans_ps = test_generate_trans_ps
    coverage = test_generate_coverages
    d1d2 = test_generate_d1d2

    restriction_map = test_get_restriction_map_mono

    forward_bam_file, reverse_bam_file = str(test_classify_reads[2]), str(test_classify_reads[3])

    hst.reattribute_reads(reads_couple = (forward_bam_file, reverse_bam_file), restriction_map = restriction_map, xs = xs, weirds = weirds, uncuts = uncuts, loops = loops, trans_ps = trans_ps, coverage = coverage, bins = BINS, d1d2 = d1d2, mode = MODE, output_dir = temp_dir_path)
    
    # yield temp_dir_path

    assert True