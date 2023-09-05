import subprocess as sp
from pathlib import Path
# from itertools import combinations
import pytest
import cooler

from .conftest import random
import hicberg.eval as hev
import hicberg.io as hio

from .conftest import temporary_folder
from .test_align import test_hic_build_index, test_hic_align, test_hic_view, test_hic_sort

from .test_utils import test_get_chromosomes_sizes, test_classify_reads, test_get_bin_table
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

BIN_INDEXES = [569, 1104, 1705, 1804, 1920, 4360, 4797, 5110, 5806, 5935]



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


def test_get_boundaries(test_get_chromosomes_sizes):
    """
    Test if boundaries from genomic coordinates are correctly computed.
    """
    boundaries = hev.get_boundaries(position = POSITION, bins = BINS, chromosome = CHROMOSOME, chrom_sizes_dict = test_get_chromosomes_sizes)

    assert boundaries == BOUNDARIES

def test_check_emptyness(test_get_chromosomes_sizes):
    """
    Test if emptiness checking before considering genomic intervals is correctly computed.
    """

    clr = cooler.Cooler(MATRIX)
    dictionary_of_intervals = hev.draw_intervals(chrom_sizes_dict  = test_get_chromosomes_sizes, nb_intervals = NB_INTERVALS_EMPTINESS, bins = BINS)
    emptiness = hev.check_emptiness(intervals = dictionary_of_intervals, matrix = clr)

    assert not emptiness

def test_get_bin_indexes(test_get_chromosomes_sizes):
    """
    Test if bin indexes are correctly computed.
    """
    clr = cooler.Cooler(MATRIX)
    dictionary_of_intervals = hev.draw_intervals(chrom_sizes_dict  = test_get_chromosomes_sizes, nb_intervals = NB_INTERVALS_EMPTINESS, bins = BINS)
    bin_indexes = hev.get_bin_indexes(dictionary = dictionary_of_intervals, matrix = clr)


    assert bin_indexes == BIN_INDEXES




