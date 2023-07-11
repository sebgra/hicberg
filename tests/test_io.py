import pytest
from os.path import join
from pathlib import Path

import glob
import shutil
import tempfile

import subprocess as sp

import numpy as np
import cooler

import hicberg.io as hio

from .conftest import temporary_folder
from .test_utils import test_get_chromosomes_sizes, test_classify_reads, test_get_bin_table
from .test_align import test_hic_build_index, test_hic_align, test_hic_view, test_hic_sort


FOLDER_TO_CREATE = "test_sample"
TEST_DICT = "data_test/chromosome_sizes.npy"
GROUP1_PAIRS = "group1.pairs"
UNRESCUED_MAP = "unrescued_map.cool"

DICT_FIRST_KEY = "chr10"
DICT_FIRST_SIZE = 745751


def test_create_folder(temporary_folder):
    """
    Test if the function creates a folder in the specified path.
    """

    temp_dir_path = Path(temporary_folder)
    hio.create_folder(sample_name = FOLDER_TO_CREATE, output_dir = temp_dir_path)

    assert (temp_dir_path / FOLDER_TO_CREATE).is_dir()

@pytest.fixture(scope = "module")
def test_build_pairs(temporary_folder, test_hic_build_index, test_hic_align, test_hic_view, test_hic_sort, test_classify_reads):
    """
    Test if the pairs file is correctly created.
    """    
    temp_dir_path = Path(temporary_folder)
    hio.build_pairs(output_dir = temp_dir_path, mode = False)
    assert (temp_dir_path / GROUP1_PAIRS).is_file() 

# @pytest.fixture(scope = "module")
def test_build_matrix(temporary_folder, test_get_bin_table, test_hic_build_index, test_hic_align, test_hic_view, test_hic_sort, test_classify_reads, test_build_pairs):
    """
    Test if the cool file is correctly created.
    """
    temp_dir_path = Path(temporary_folder)
    hio.build_matrix(output_dir = temp_dir_path, mode = False)
    
    assert (temp_dir_path / UNRESCUED_MAP).is_file()

def test_load_dictionary(test_get_chromosomes_sizes):
    """
    Test if the dictionary saved in .npy format is correctly loaded.
    """
    dictionary = hio.load_dictionary(dictionary = test_get_chromosomes_sizes)

    assert DICT_FIRST_KEY == list(dictionary.keys())[0]
    assert DICT_FIRST_SIZE == list(dictionary.values())[0]

def test_load_cooler():
    pass

def test_merge_predictions():
    pass