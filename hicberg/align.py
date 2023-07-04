from os import getcwd
from os.path import join
from pathlib import Path
import subprocess as sp
import click


def hic_build_index(genome : str, output : str = None , cpus : int = 1 , verbose : bool = False) -> None:
    """
    Building of bowtie2 index (.bt2l files) for read alignement.

    Parameters
    ----------
    genome : str
        Path to the genome file along which reads are going to be aligned.
    cpus : int, optional
        Number of threads allocated for the alignment, by default 1
    output : str, optional
        Path where the Bowtie2 index files should be stored, by default None
    verbose : bool, optional
        Set weither or not the shell command should be printed, by default False
    """

    try:

        sp.check_output(["bowtie2-build", "-h"])

    except OSError:

        raise RuntimeError(
            "bowtie2-build not found; check if it is installed and in $PATH\n install Bowtie2 with : conda install bowtie2"
        )
    
    genome_path = Path(genome)

    if not genome_path.is_file():
        
        raise ValueError(f"Genome file {genome} not found")

    if output is None:    
        output_path = Path(getcwd())

    # output_path  = Path(output)
    else : 
        output_path = Path(output)

        
    if not output_path.exists():

        raise ValueError(f"Output path {output_path} does not exist. Please provide existing ouput path.")
    
    sample = Path(genome).stem
    index_path = Path(output, sample)

    cmd_index = f"bowtie2-build -q -f  --large-index {genome} {index_path}"

    if verbose:

        print(cmd_index)

    sp.check_call([cmd_index], shell=True)

    print(f"Index built at {index_path}")


def hic_align(genome : str, index : str, fq_for : str, fq_rev : str, sensitivity : str = 'very-sensitive', max_alignment :  int = None, cpus : int = 1, output : str = None, verbose : bool = False) -> None:
    """
    AI is creating summary for hic_align

    Parameters
    ----------
    genome : str
        Path to the genome file along which reads are going to be aligned.
    index : str
        Path to the index of the genome along which reads are going to be aligned (path to .bt2l files). Default to None, index files are searched to sample_name/data/index/sample_name.
    fq_for : str
        Path to .fatsa containing set of reads to align (forward mate).
    fq_rev : str
        Path to .fatsa containing set of reads to align (forward mate).
    sensitivity : str, optional
        Sensitivity of the alignment., by default 'very_sensitive'
    max_alignment : int, optional
        Maximum number of alignments to be returned, by default None
    cpus : int, optional
        Number of threads allocated for the alignment, by default None
    output : str, optional
        Path where the alignement files (.sam) should be stored, by default None
    verbonse : bool, optional
        Set weither or not the shell command should be printed, by default False
    """    

    fq_for_path, fq_rev_path  = Path(fq_for), Path(fq_rev)

    if not fq_for_path.is_file() or not fq_rev_path.is_file():

        raise IOError(f"Wront path to fastq files : {fq_for_path} or {fq_rev_path} given. \
                    Pease provide existing files.")
    
    if output is None:    
        output_path = Path(getcwd())
    
    else : 
        output_path = Path(output)

    if not output_path.exists():

        raise ValueError(f"Output path {output_path} does not exist. Please provide existing ouput path.")

    index_path = Path(output_path / index)

    
    if max_alignment is None:
        
        cmd_alignment_rev = f"bowtie2 --{sensitivity} -p {cpus} -a -x {index_path} -S {output_path / '1.sam'} {fq_for}"
        cmd_alignment_for = f"bowtie2 --{sensitivity} -p {cpus} -a -x {index_path} -S {output_path / '2.sam'} {fq_rev}"

    else:
            
        cmd_alignment_rev = f"bowtie2 --{sensitivity} -p {cpus} -k {max_alignment} {sensitivity} -p {cpus}  -x {index_path} -S {output_path / '1.sam'} {fq_rev}"
        cmd_alignment_for = f"bowtie2 --{sensitivity} -p {cpus} -k {max_alignment} {sensitivity} -p {cpus}  -x {index_path} -S {output_path / '2.sam'} {fq_for}"

    if verbose :

        print(cmd_alignment_for)
        print(cmd_alignment_rev)

    
    sp.check_call([cmd_alignment_for], shell=True)
    sp.check_call([cmd_alignment_rev], shell=True)

    print(f"Alignement done at {output_path}")







def hic_view():

    pass


def hic_sort():

    pass

def hic_index():

    pass



