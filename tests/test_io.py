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

FOLDER_TO_CREATE = "test_sample"

@pytest.fixture(scope="session")
def temporary_folder():
    fn = tempfile.TemporaryDirectory()
    # yiled temporary folder name to be used as Path object
    yield fn.name

@pytest.fixture(name = "create_folder")
def test_create_folder(temporary_folder):

    temp_dir_path = Path(temporary_folder)
    hio.create_folder(sample_name = FOLDER_TO_CREATE, output_dir = temp_dir_path)

    assert (temp_dir_path / FOLDER_TO_CREATE).is_dir()

def test_build_pairs():
    pass

def test_build_matrix():
    pass

def test_load_dictionary():
    pass

def test_load_cooler():
    pass

def test_merge_predictions():
    pass