from abc import abstractmethod
from pathlib import Path
from hicberg.core.command import BaseCommand

class Aligner(BaseCommand):
    @abstractmethod
    def index(self, genome: Path, output_index: Path, **kwargs):
        """
        Generate tool specific command for genome indexing if needed.

        Parameters
        ----------
        genome : Path
            Path to the reference genome file (e.g., FASTA format).
        output_index : Path
            Path to the output index file.
        """
        pass

    @abstractmethod
    def align(self, index: Path, reads: Path, output_align: Path, **kwargs):
        """
        Generate tool specific command for aligning reads to the reference genome.
        Parameters
        ----------
        index : Path
            Path to the index file.
        reads : Path
            Path to the reads file (e.g., FASTQ format).
        output_align : Path
            Path to the output alignment file.
        """
        pass
    

    