from os import getcwd
from os.path import join
from pathlib import Path
import subprocess as sp
import click
from hicberg import logger



def hic_build_index(genome : str, output_dir  : str = None , cpus : int = 1 , verbose : bool = False) -> None:
    """
    Building of bowtie2 index (.bt2l files) for read alignement.

    Parameters
    ----------
    genome : str
        Path to the genome file along which reads are going to be aligned.
    cpus : int, optional
        Number of threads allocated for the alignment, by default 1
    output_dir : str, optional
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

    if output_dir is None:    
        output_path = Path(getcwd())

    # output_path  = Path(output)
    else : 
        output_path = Path(output_dir)

        
    if not output_path.exists():

        raise ValueError(f"Output path {output_path} does not exist. Please provide existing ouput path.")
    
    sample = Path(genome).stem
    index_path = Path(output_dir, sample)

    cmd_index = f"bowtie2-build -q -f  --large-index {genome} {index_path}"

    if verbose:

        logger.info(cmd_index)

    sp.check_call([cmd_index], shell=True)

    logger.info(f"Index built at {index_path}")


    return index_path


def hic_align(genome : str, index : str, fq_for : str, fq_rev : str, sensitivity : str = 'very-sensitive', max_alignment :  int = None, cpus : int = 1, output_dir : str = None, verbose : bool = False) -> None:
    """
    lignement of reads from HiC experiments along an indexed genome.

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
        Number of threads allocated for the alignment, by default 1
    output_dir : str, optional
        Path where the alignement files (.sam) should be stored, by default None
    verbonse : bool, optional
        Set weither or not the shell command should be printed, by default False
    """    

    fq_for_path, fq_rev_path  = Path(fq_for), Path(fq_rev)

    if not fq_for_path.is_file() or not fq_rev_path.is_file():

        raise IOError(f"Wront path to fastq files : {fq_for_path} or {fq_rev_path} given. \
                    Pease provide existing files.")
    
    if output_dir is None:    
        output_path = Path(getcwd())
    
    else : 
        output_path = Path(output_dir)

    if not output_path.exists():

        raise ValueError(f"Output path {output_path} does not exist. Please provide existing ouput path.")

    index_path = Path(output_path / index)

    
    if max_alignment is None:
        
        cmd_alignment_for = f"bowtie2 --{sensitivity} -p {cpus} -a -x {index_path} -S {output_path / '1.sam'} {fq_rev}"
        cmd_alignment_rev = f"bowtie2 --{sensitivity} -p {cpus} -a -x {index_path} -S {output_path / '2.sam'} {fq_for}"

    elif max_alignment is not None:
            
        cmd_alignment_for = f"bowtie2 --{sensitivity} -p {cpus} -k {max_alignment}  -p {cpus}  -x {index_path} -S {output_path / '1.sam'} {fq_for}"
        cmd_alignment_rev = f"bowtie2 --{sensitivity} -p {cpus} -k {max_alignment}  -p {cpus}  -x {index_path} -S {output_path / '2.sam'} {fq_rev}"

    if verbose :

        logger.info(cmd_alignment_for)
        logger.info(cmd_alignment_rev)

    p_for = sp.Popen([cmd_alignment_for], shell=True, stdout = sp.PIPE, stderr = sp.PIPE)
    p_rev = sp.Popen([cmd_alignment_rev], shell=True, stdout = sp.PIPE,  stderr = sp.PIPE)

    stdout_for, stderr_for  = p_for.communicate()
    stdout_rev, stderr_rev = p_rev.communicate()


    if stdout_for : 
        logger.info(stdout_for.decode('ascii'))
    if stderr_for : 
        logger.info(stderr_for.decode('ascii'))

    if stdout_rev:
        logger.info(stdout_rev.decode('ascii'))
    if stderr_rev:
        logger.info(stderr_rev.decode('ascii'))

    logger.info(f"Alignement saved at {output_path}")




def hic_view(sam_for : str = "1.sam", sam_rev : str = "2.sam", cpus : int = 1, output_dir : str = None, verbose : bool = False) -> None:
    """
    Conversion of .sam alignement files to .bam alignement format (using samtools).

    Parameters
    ----------
    sam_for : str, optional
        Path to forward .sam alignment file, by default "1.sam"
    sam_rev : str, optional
        Path to reverse .sam alignment file, by default "2.sam"
    cpus : int, optional
        Number of threads allocated for the alignment, by default 1
    output_dir : str, optional
        Path where the alignement files (.bam) should be stored, by default None
    verbose : bool, optional
        Set weither or not the shell command should be printed, by default False
    """

    try:

        sp.check_output(["samtools", "--help"])

    except OSError:

        raise RuntimeError(
            "Samtools not found; check if it is installed and in $PATH\n install Samtools with : conda install samtools"
        )


    if output_dir is None:    
        output_path = Path(getcwd())
    
    else : 
        output_path = Path(output_dir)

    if not output_path.exists():

        raise ValueError(f"Output path {output_path} does not exist. Please provide existing ouput path.")
    

    cmd_view_for = f"samtools view -S -b {output_path / sam_for} -o {output_path / '1.bam'} --threads {cpus}"
    cmd_view_rev = f"samtools view -S -b {output_path / sam_rev} -o {output_path / '2.bam'} --threads {cpus}"

    if verbose:

        logger.info(cmd_view_for)
        logger.info(cmd_view_rev)

    sp.check_call([cmd_view_for], shell=True)
    sp.check_call([cmd_view_rev], shell=True)

    # Delete .sam files after .bam conversion
    (output_path / sam_for).unlink()
    (output_path / sam_rev).unlink()

    logger.info(f"Compressed alignemnet alignement done at {output_path}")




def hic_sort(bam_for : str = "1.bam", bam_rev : str = "2.bam", cpus : int = 1, output_dir : str = None, verbose : bool = False) -> None:
    """
    Sort .bam alignement files by read_name  (using samtools).

    Parameters
    ----------
    bam_for : str, optional
        Forward alignement file to be sorted, by default "1.bam"
    bam_rev : str, optional
        Reverse alignement file to be sorted, by default "2.bam"
    cpus : int, optional
        Number of threads allocated for the alignment, by default 1
    output_dir : str, optional
        Path where the alignement files (.bam) should be stored, by default None
    verbose : bool, optional
        Set weither or not the shell command should be printed, by default False

    """

    try:

        sp.check_output(["samtools", "--help"])

    except OSError:

        raise RuntimeError(
            "Samtools not found; check if it is installed and in $PATH\n install Samtools with : conda install samtools"
        )


    if output_dir is None:    
        output_path = Path(getcwd())
    
    else : 
        output_path = Path(output_dir)

    if not output_path.exists():

        raise ValueError(f"Output path {output_path} does not exist. Please provide existing ouput path.")
    

    cmd_sort_for = f"samtools sort -n {output_path / '1.bam'} -o {output_path / '1.sorted.bam'} --threads {cpus}"
    cmd_sort_rev = f"samtools sort -n {output_path / '2.bam'} -o {output_path / '2.sorted.bam'} --threads {cpus}"

    if verbose:

        logger.info(cmd_sort_for)
        logger.info(cmd_sort_rev)

    sp.check_call([cmd_sort_for], shell=True)
    sp.check_call([cmd_sort_rev], shell=True)

    (output_path / '1.bam').unlink()
    (output_path / '2.bam').unlink()

    logger.info(f"Sorted alignement done at {output_path}")



def hic_index(bam_for : str = "1.sorted.bam", bam_rev : str = "2.sorted.bam", cpus : int = 1, output_dir : str = None, verbose : bool = False) -> None:
    """
    Index a coordinate-sorted BGZIP-compressed SAM, BAM or CRAM file for fast random access.

    Parameters
    ----------
    bam_for : str, optional
        Forward alignement file to be indexed, by default "1.sorted.bam"
    bam_rev : str, optional
        Reverse alignement file to be indexed,, by default "2.sorted.bam"
    cpus : int, optional
        Number of threads allocated for the alignment, by default 1
    output_dir : str, optional
        Path where the alignement files (.bam) should be stored, by default None
    verbose : bool, optional
        Set weither or not the shell command should be printed, by default False

    Returns
    -------
    [type]
        [description]
    """

    try:

        sp.check_output(["samtools", "--help"])

    except OSError:

        raise RuntimeError(
            "Samtools not found; check if it is installed and in $PATH\n install Samtools with : conda install samtools"
        )

    if output_dir is None:    
        output_path = Path(getcwd())
    
    else : 
        output_path = Path(output_dir)

    if not output_path.exists():

        raise ValueError(f"Output path {output_path} does not exist. Please provide existing ouput path.")
    

    cmd_index_for = f"samtools index -b {bam_for} -@ {cpus}"
    cmd_index_rev = f"samtools index -b {bam_rev} -@ {cpus}"

    if verbose:

        logger.info(cmd_index_for)
        logger.info(cmd_index_rev)

    sp.check_call([cmd_index_for], shell=True)
    sp.check_call([cmd_index_rev], shell=True)

    logger.info(f"Indexed alignement done at {output_path}")



