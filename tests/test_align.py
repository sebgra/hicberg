import pytest
from pathlib import Path
import tempfile
import hicberg.align as hal


@pytest.fixture()
def test_hic_build_index():

    # temp_dir = tempfile.TemporaryDirectory()
    temp_dir_path = Path("/home/sardine/Bureau/test")
    genome_path  = Path("data_test/SC288_with_micron.fa")
    # print(temp_dir_path)
    hal.hic_build_index(genome = genome_path, output = temp_dir_path, verbose = True)

    # print(f"genome_path.stem : {genome_path.stem}")

    # print(f"to be yielded : {temp_dir_path / genome_path.stem}")

    yield temp_dir_path / genome_path.stem

    # for f in temp_dir_path.iterdir():
    #     print(f"File : {f}")

    # print(temp_dir_path.iterdir())
    assert  any(temp_dir_path.iterdir()) == True

def test_hic_align(test_hic_build_index):

    temp_dir_path = Path("/home/sardine/Bureau/test")
    genome_path  = Path("data_test/SC288_with_micron.fa")
    fq_for_path = Path("data_test/forward_reads_test.fq.gz")
    fq_rev_path = Path("data_test/reverse_reads_test.fq.gz")
    
    hal.hic_align(genome = genome_path, index = test_hic_build_index, fq_for = fq_for_path, fq_rev = fq_rev_path, output = temp_dir_path, verbose = True)

    assert True

def test_hic_view():

    pass


def test_hic_sort():

    pass

def test_hic_index():

    pass
