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
from .test_utils import test_get_chromosomes_sizes


FOLDER_TO_CREATE = "test_sample"
TEST_DICT = "data_test/chromosome_sizes.npy"

DICT_FIRST_KEY = "chr10"
DICT_FIRST_SIZE = 745751


def test_create_folder(temporary_folder):
    """
    Test if the function creates a folder in the specified path.
    """

    temp_dir_path = Path(temporary_folder)
    hio.create_folder(sample_name = FOLDER_TO_CREATE, output_dir = temp_dir_path)

    assert (temp_dir_path / FOLDER_TO_CREATE).is_dir()

def test_build_pairs():
    pass

def test_build_matrix():
    pass

def test_load_dictionary(test_get_chromosomes_sizes):
    
    dictionary = hio.load_dictionary(dictionary = test_get_chromosomes_sizes)

    assert DICT_FIRST_KEY == list(dictionary.keys())[0]
    assert DICT_FIRST_SIZE == list(dictionary.values())[0]

def test_load_cooler():
    pass

def test_merge_predictions():
    pass