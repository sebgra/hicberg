import pytest
from unittest import mock  # Import the mock module
import subprocess as sp
import logging
from pathlib import Path
import hicberg.align as hal

from .conftest import temporary_folder


GENOME = "data_test/sub_genome.fa"
FOR_FQ = "data_test/forward_reads_test.fq.gz"
REV_FQ = "data_test/reverse_reads_test.fq.gz"
FOR_SAM = "1.sam"
REV_SAM = "2.sam"
FOR_BAM = "1.bam"
REV_BAM = "2.bam"
FOR_SORTED_BAM = "1.sorted.bam"
REV_SORTED_BAM = "2.sorted.bam"

#####################################################
################### Fixtures ########################
#####################################################

@pytest.fixture(scope="function")
def capture_logs(caplog):
    """Fixture to capture log messages."""
    yield caplog


@pytest.fixture(scope = "session")
def hic_build_index_fixture(temporary_folder):
    """
    Fixture fot index building tests.
    """    

    temp_dir_path = Path(temporary_folder)
    genome_path  = Path(GENOME)
    
    hal.hic_build_index(genome = genome_path, output_dir = temp_dir_path, verbose = False)

    assert  any(temp_dir_path.iterdir()) == True
    yield temp_dir_path / genome_path.stem


def test_hic_build_index(temporary_folder):
    """
    Test if the bowtie2 index is correctly created
    """    

    temp_dir_path = Path(temporary_folder)
    genome_path  = Path(GENOME)
    
    hal.hic_build_index(genome = genome_path, output_dir = temp_dir_path, verbose = False)

    assert  any(temp_dir_path.iterdir()) == True

def test_hic_build_index_bowtie2_not_found(monkeypatch):
    """
    Test that hic_build_index raises a RuntimeError when bowtie2-build is not found.
    """
    # 
    def mock_check_output(cmd):
        """
        Mock function to simulate the behavior of subprocess.check_output.
        Raises OSError if the command is 'bowtie2-build' to simulate 
        the scenario where it's not found.
        """
        if cmd[0] == "bowtie2-build":
            raise OSError("Mock OSError: Command not found")
        else:
            return sp.check_output(cmd)  # Allow other commands to execute normally
    # Replace the actual subprocess.check_output with our mock function
    monkeypatch.setattr(sp, "check_output", mock_check_output)

    # Call the function and assert that it raises the expected RuntimeError
    with pytest.raises(RuntimeError) as excinfo:
        hal.hic_build_index(genome="genome.fasta", output_dir=".") 
    assert "bowtie2-build not found; check if it is installed and in $PATH\n install Bowtie2 with : conda install bowtie2" in str(excinfo.value)

def test_hic_build_index_genome_not_found(temporary_folder):
    """
    Test if the file not found ValueError is correctly raised
    """

    print("Test is effective")

    temp_dir_path = Path(temporary_folder)
    genome_path  = Path("wrong_path")

    with pytest.raises(ValueError) as excinfo:
        hal.hic_build_index(genome = genome_path, output_dir = temp_dir_path, verbose = True)

    assert str(excinfo.value) == f"Genome file {genome_path} not found"

def test_hic_build_index_logging(capture_logs, temporary_folder):
    """
    Test that hic_build_index logs the expected messages.
    """
    temp_dir_path = Path(temporary_folder)
    genome_path = Path(GENOME)

    # Mock sp.run
    with mock.patch("subprocess.run"):
        hal.hic_build_index(genome=genome_path, output_dir=temp_dir_path, verbose=True)

    # Get the logger used in hic_build_index
    hic_build_logger = logging.getLogger("hicberg")  # Replace your_module

    # Access logs from the specific logger
    captured_messages = [
        record.message for record in capture_logs.records if record.name == "hicberg"
    ]

    # Check for expected log messages
    assert "Start building index for alignment" in captured_messages
    assert f"Index built at {temp_dir_path / genome_path.stem}" in captured_messages
    assert f"bowtie2-build -q -f --threads 1 --large-index {genome_path} {temp_dir_path / genome_path.stem}" in captured_messages



@pytest.fixture(scope = "session")
def test_hic_align(hic_build_index_fixture, temporary_folder):
    """
    Test if the alignement is correctly performed.
    """

    temp_dir_path = Path(temporary_folder)
    genome_path  = Path(GENOME)
    fq_for_path = Path(FOR_FQ)
    fq_rev_path = Path(REV_FQ)

    hal.hic_align(index = hic_build_index_fixture, fq_for = fq_for_path, fq_rev = fq_rev_path, output_dir = temp_dir_path, verbose = True)

    for_sam_path = temp_dir_path / FOR_SAM
    rev_sam_path = temp_dir_path / REV_SAM

    # Check if the alignement files are created
    assert for_sam_path.is_file()
    assert rev_sam_path.is_file()

@pytest.fixture(scope = "session")
def test_hic_view(temporary_folder, hic_build_index_fixture, test_hic_align):
    """
    Test if the alignement compression is correctly performed.
    """
    temp_dir_path = Path(temporary_folder)
    hal.hic_view(output_dir = temp_dir_path, verbose = True)

    for_bam_path = temp_dir_path / FOR_BAM
    rev_bam_path = temp_dir_path / REV_BAM

    # Check if the alignement files are created
    assert for_bam_path.is_file()
    assert rev_bam_path.is_file()

@pytest.fixture(scope = "session")
def test_hic_sort(temporary_folder, hic_build_index_fixture, test_hic_align, test_hic_view):
    """
    Test if the alignement sorting by read ids is correctly performed.
    """

    temp_dir_path = Path(temporary_folder)

    hal.hic_sort(output_dir = temp_dir_path, verbose = True)

    for_sorted_bam_path = temp_dir_path / FOR_SORTED_BAM
    rev_sorted_bam_path = temp_dir_path / REV_SORTED_BAM


    yield for_sorted_bam_path, rev_sorted_bam_path

    # Check if the sorted alignement files are created
    # assert for_sorted_bam_path.is_file()
    # assert rev_sorted_bam_path.is_file()

@pytest.fixture(scope = "session")
def test_hic_index(temporary_folder, hic_build_index_fixture, test_hic_align, test_hic_view, test_hic_sort):
    """
    Test if the alignement indexing is correctly performed.
    """

    temp_dir_path = Path(temporary_folder)

    hal.hic_sort(output_dir = temp_dir_path, verbose = True)

    for_indexed_bam_path = temp_dir_path / FOR_SORTED_BAM
    rev_indexed_bam_path = temp_dir_path / REV_SORTED_BAM

    # Check if the sorted and indexed alignement files are created
    assert for_indexed_bam_path.is_file()
    assert rev_indexed_bam_path.is_file()
