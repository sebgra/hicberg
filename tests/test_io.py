import pytest
from pathlib import Path
import cooler

import hicberg.io as hio

from .conftest import temporary_folder
from .test_utils import test_get_chromosomes_sizes, test_classify_reads, test_get_bin_table, test_chunk_bam
from .test_align import test_hic_build_index, test_hic_align, test_hic_view, test_hic_sort


FOLDER_TO_CREATE = "test_sample"
TEST_DICT = "data_test/chromosome_sizes.npy"
GROUP1_PAIRS = "group1.pairs"
UNRESCUED_MAP = "unrescued_map.cool"

PREDICTED_BAM_FORWARD = "group2.1.rescued.bam"
PREDICTED_BAM_REVERSE = "group2.2.rescued.bam"


DICT_FIRST_KEY = "chr10"
DICT_FIRST_SIZE = 745751

#TODO : Improve with force mode
def test_create_folder(temporary_folder):
    """
    Test if the function creates a folder in the specified path.
    """

    temp_dir_path = Path(temporary_folder)
    hio.create_folder(sample_name = FOLDER_TO_CREATE, output_dir = temp_dir_path)

    assert (temp_dir_path / FOLDER_TO_CREATE).is_dir()

@pytest.fixture(scope = "session")
def test_build_pairs(temporary_folder, test_hic_build_index, test_hic_align, test_hic_view, test_hic_sort, test_classify_reads):
    """
    Test if the pairs file is correctly created.
    """    
    temp_dir_path = Path(temporary_folder)
    hio.build_pairs(output_dir = temp_dir_path, mode = False)
    assert (temp_dir_path / GROUP1_PAIRS).is_file() 

@pytest.fixture(scope = "session")
def test_build_matrix(temporary_folder, test_get_bin_table, test_hic_build_index, test_hic_align, test_hic_view, test_hic_sort, test_classify_reads, test_build_pairs):
    """
    Test if the cool file is correctly created.
    """
    temp_dir_path = Path(temporary_folder)
    hio.build_matrix(output_dir = temp_dir_path, mode = False)

    unrescued_map_path = temp_dir_path / UNRESCUED_MAP
    
    yield unrescued_map_path

    assert unrescued_map_path.is_file()


def test_load_dictionary(test_get_chromosomes_sizes):
    """
    Test if the dictionary saved in .npy format is correctly loaded.
    """
    dictionary = hio.load_dictionary(dictionary = test_get_chromosomes_sizes)

    assert DICT_FIRST_KEY == list(dictionary.keys())[0]
    assert DICT_FIRST_SIZE == list(dictionary.values())[0]

def test_load_cooler(test_build_matrix):

    matrix = hio.load_cooler(matrix = test_build_matrix)

    assert isinstance(matrix, cooler.Cooler)

def test_merge_predictions(temporary_folder, test_chunk_bam):

    temp_dir_path = Path(temporary_folder)

    folder_path = Path("data_test/alignments/chunks")

    hio.merge_predictions(output_dir = folder_path, clean = False)

    assert (folder_path / PREDICTED_BAM_FORWARD).is_file()
    assert (folder_path / PREDICTED_BAM_REVERSE).is_file()

    assert True

def test_tidy_folder():

    assert True 