import pytest
import os
from pathlib import Path
import tempfile
import hicberg.align as hal

@pytest.fixture(scope="session")
def temporary_folder():
    fn = tempfile.TemporaryDirectory()
    yield fn.name
    fn.cleanup()



@pytest.fixture(name = "bowtie2_index")
def test_hic_build_index(temporary_folder):

    temp_dir_path = Path(temporary_folder)
    genome_path  = Path("data_test/SC288_with_micron.fa")
    
    hal.hic_build_index(genome = genome_path, output = temp_dir_path, verbose = True)

    yield temp_dir_path / genome_path.stem

    assert  any(temp_dir_path.iterdir()) == True

def test_hic_align(bowtie2_index, temporary_folder):

    temp_dir_path = Path(temporary_folder)
    genome_path  = Path("data_test/SC288_with_micron.fa")
    fq_for_path = Path("data_test/forward_reads_test.fq.gz")
    fq_rev_path = Path("data_test/reverse_reads_test.fq.gz")

    hal.hic_align(genome = genome_path, index = bowtie2_index, fq_for = fq_for_path, fq_rev = fq_rev_path, output = temp_dir_path, verbose = True)

    assert True

def test_hic_view():

    pass


def test_hic_sort():

    pass

def test_hic_index():

    pass
