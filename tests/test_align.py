import pytest
import os
from pathlib import Path
import tempfile
import hicberg.align as hal

@pytest.fixture(scope="session")
def temporary_folder():
    fn = tempfile.TemporaryDirectory()
    # yiled temporary folder name to be used as Path object
    yield fn.name



@pytest.fixture(name = "bowtie2_index")
def test_hic_build_index(temporary_folder):

    temp_dir_path = Path(temporary_folder)
    genome_path  = Path("data_test/SC288_with_micron.fa")
    
    hal.hic_build_index(genome = genome_path, output = temp_dir_path, verbose = True)

    yield temp_dir_path / genome_path.stem

    assert  any(temp_dir_path.iterdir()) == True

@pytest.fixture(name = "bowtie2_alignment")
def test_hic_align(bowtie2_index, temporary_folder):

    temp_dir_path = Path(temporary_folder)
    genome_path  = Path("data_test/SC288_with_micron.fa")
    fq_for_path = Path("data_test/forward_reads_test.fq.gz")
    fq_rev_path = Path("data_test/reverse_reads_test.fq.gz")

    hal.hic_align(genome = genome_path, index = bowtie2_index, fq_for = fq_for_path, fq_rev = fq_rev_path, output = temp_dir_path, verbose = True)

    for_sam_path = temp_dir_path / '1.sam'
    rev_sam_path = temp_dir_path / '2.sam'

    # Check if the alignement files are created
    assert for_sam_path.is_file()
    assert rev_sam_path.is_file()

@pytest.fixture(name = "samtools_view")
def test_hic_view(temporary_folder, bowtie2_index, bowtie2_alignment):

    temp_dir_path = Path(temporary_folder)
    hal.hic_view(output = temp_dir_path, verbose = True)

    for_bam_path = temp_dir_path / '1.bam'
    rev_bam_path = temp_dir_path / '2.bam'

    # Check if the alignement files are created
    assert for_bam_path.is_file()
    assert rev_bam_path.is_file()

@pytest.fixture(name = "samtools_sort")
def test_hic_sort(temporary_folder, bowtie2_index, bowtie2_alignment, samtools_view):

    temp_dir_path = Path(temporary_folder)

    hal.hic_sort(output = temp_dir_path, verbose = True)

    for_sorted_bam_path = temp_dir_path / '1.sorted.bam'
    rev_sorted_bam_path = temp_dir_path / '2.sorted.bam'

    # Check if the alignement files are created
    assert for_sorted_bam_path.is_file()
    assert rev_sorted_bam_path.is_file()


def test_hic_index(temporary_folder, bowtie2_index, bowtie2_alignment, samtools_view, samtools_sort):

    temp_dir_path = Path(temporary_folder)

    hal.hic_sort(output = temp_dir_path, verbose = True)

    for_indexed_bam_path = temp_dir_path / '1.sorted.bam'
    rev_indexed_bam_path = temp_dir_path / '2.sorted.bam'

    # Check if the alignement files are created
    assert for_indexed_bam_path.is_file()
    assert rev_indexed_bam_path.is_file()
