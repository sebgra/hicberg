import pytest
from pathlib import Path
import tempfile
import hicberg.align as hal

def test_hic_align():


    pass

# @pytest.fixture()
def test_hic_build_index():

    # temp_dir = tempfile.TemporaryDirectory()
    temp_dir_path = Path("/home/sardine/Bureau/test")
    print(temp_dir_path)
    hal.hic_build_index(genome = "data_test/SC288_with_micron.fa", output = temp_dir_path, verbose = True)
    print(temp_dir_path.iterdir())
    assert  any(temp_dir_path.iterdir()) == True

def test_hic_view():

    pass


def test_hic_sort():

    pass

def test_hic_index():

    pass
