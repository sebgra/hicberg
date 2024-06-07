import logging
import uuid
from os import getcwd, mkdir
from os.path import join
from pathlib import Path

from glob import glob
from shutil import rmtree

import subprocess as sp

import numpy as np
import cooler
import pysam as ps

from hicberg import logger

def create_folder(sample_name : str  = None, output_dir : str = None, force : bool = False) -> None:
    """
    Creates folder architecture to store results and intermediate files for the full HiC-BERG pipeline.

    Parameters
    ----------
    sample_name : str
        Name of the folder to be created.
    force : bool
        Set if existing folder has to be deleted before folder creation.
    output_dir : str
        Path where the folder will be created.

    Returns
    -------
    [str]
        Path of the folder created
    """

    logger.info(f"Creating folder {sample_name} in {output_dir}")
    
    if sample_name is None:

        sample_name = "sample"

    if output_dir is None:

        folder_path = Path(getcwd(), sample_name)

    else:

        folder_path = Path(output_dir, sample_name)

    if folder_path.exists() and force : 

        rmtree(folder_path)

    mkdir(folder_path)
    mkdir(folder_path / "index")
    mkdir(folder_path / "alignments")
    mkdir(folder_path / "statistics")
    mkdir(folder_path / "contacts")
    mkdir(folder_path / "contacts" / "matrices")
    mkdir(folder_path / "contacts" / "pairs")
    mkdir(folder_path / "plots")
    mkdir(folder_path / "plots" / "contact_maps")
    mkdir(folder_path / "plots" / "densities")
    mkdir(folder_path / "plots" / "ps")
    mkdir(folder_path / "plots" / "coverages")

    logger.info(f"Folder {sample_name} in {folder_path} created.")

    return folder_path.as_posix()

    

def build_pairs(bam_for : str = "group1.1.bam", bam_rev : str = "group1.2.bam", bam_for_rescued :str = "group2.1.rescued.bam", bam_rev_rescued : str = "group2.2.rescued.bam", mode : bool = False, output_dir : str = None) -> None:
    """
    Build pairs of reads from the aligned reads.

    Parameters
    ----------
    bam_for : str, optional
        Path to forward .bam file for the construction of the .pairs equivalent file (non rescued)., by default "group1.1.bam"
    bam_rev : str, optional
        Path to reverse .bam file for the construction of the .pairs equivalent file (non rescued)., by default None, by default "group1.2.bam"
    bam_for_rescued : str, optional
        Path to forward .bam file for the construction of the .pairs equivalent file (rescued)., by default "group2.1.rescued.bam"
    bam_rev_rescued : str, optional
        Path to reverse .bam file for the construction of the .pairs equivalent file (rescued)., by default "group2.2.rescued.bam"
    mode : bool, optional
        Choose wether the mode is rescued or unrescued to construct associated .pairs file, by default False
    output_dir : str, optional
        Path where the alignment files (.sam) should be stored, by default None
    """

    output_path = Path(output_dir)
    chromosome_sizes_path = Path(output_path / "chromosome_sizes.npy")
    chromosome_sizes = load_dictionary(chromosome_sizes_path)


    if not output_path.exists():

        raise ValueError(f"Output path {output_path} does not exist. Please provide existing output path.")

    if not mode:

        logger.info(f"Start building pairs file for unambiguously aligned reads")

        bam_for_path = Path(output_path / bam_for)
        bam_rev_path = Path(output_path / bam_rev)

        if not bam_for_path.exists():

            raise ValueError(f"Forward bam file {bam_for} not found")
        
        if not bam_rev_path.exists():
                
                raise ValueError(f"Reverse bam file {bam_rev} not found")
        
        bam_for_handler = ps.AlignmentFile(bam_for_path, "rb")
        bam_rev_handler = ps.AlignmentFile(bam_rev_path, "rb")

        with open(output_path / "group1.pairs", "w") as f_out:

            f_out.write("## pairs format v1.0\n")
            f_out.write("#columns: readID chr1 pos1 strand1 chr2 pos2 strand2\n")
            
            for chromosome, size in chromosome_sizes.items():

                f_out.write(f"#chromsize: {chromosome} {size}\n")

            for forward_read, reverse_read in zip(bam_for_handler, bam_rev_handler):

                if forward_read.query_name != reverse_read.query_name:

                    raise ValueError(f"Forward and reverse reads do not match. Please check the bam files.")
                
                f_out.write(f"{forward_read.query_name}\t{forward_read.reference_name}\t{forward_read.pos}\t{reverse_read.reference_name}\t{reverse_read.pos}\t{'+' if forward_read.flag == 0 or forward_read.flag == 256 else '-'}\t{'+' if reverse_read.flag == 0 or forward_read.flag == 256 else '-'}\n")

        f_out.close()
        bam_for_handler.close()
        bam_rev_handler.close()

    elif mode: 

        logger.info(f"Start building pairs file for ambiguously aligned reads")
        bam_for_path = Path(output_path / bam_for)
        bam_rev_path = Path(output_path / bam_rev)
        bam_for_path_rescued = Path(output_path / bam_for_rescued)
        bam_rev_path_rescued = Path(output_path / bam_rev_rescued)

        if not bam_for_path.exists():

            raise ValueError(f"Forward bam file {bam_for_rescued} not found")
        
        if not bam_rev_path.exists():
                
            raise ValueError(f"Reverse bam file {bam_rev_rescued} not found")
        
        if not bam_for_path_rescued.exists():
            
            raise ValueError(f"Forward rescued bam file {bam_for_rescued} not found")
        
        if not bam_rev_path_rescued.exists():
            
            raise ValueError(f"Reverse rescued bam file {bam_rev_rescued} not found")
        

        bam_for_handler = ps.AlignmentFile(bam_for_path, "rb")
        bam_rev_handler = ps.AlignmentFile(bam_rev_path, "rb")
        bam_for_handler_rescued = ps.AlignmentFile(bam_for_path_rescued, "rb")
        bam_rev_handler_rescued = ps.AlignmentFile(bam_rev_path_rescued, "rb")

        with open(output_path / "all_group.pairs", "w") as f_out:

            f_out.write("## pairs format v1.0\n")
            f_out.write("#columns: readID chr1 pos1 strand1 chr2 pos2 strand2\n")
            
            for chromosome, size in chromosome_sizes.items():

                f_out.write(f"#chromsize: {chromosome} {size}\n")

            for forward_read, reverse_read in zip(bam_for_handler, bam_rev_handler):

                if forward_read.query_name != reverse_read.query_name:

                    raise ValueError(f"Forward and reverse reads do not match. Please check the bam files.")

                f_out.write(f"{forward_read.query_name}\t{forward_read.reference_name}\t{forward_read.pos}\t{reverse_read.reference_name}\t{reverse_read.pos}\t{'+' if forward_read.flag == 0 or forward_read.flag == 256 else '-'}\t{'+' if reverse_read.flag == 0 or forward_read.flag == 256 else '-'}\n")

            for forward_read, reverse_read in zip(bam_for_handler_rescued, bam_rev_handler_rescued):

                if forward_read.query_name != reverse_read.query_name:

                    raise ValueError(f"Forward and reverse reads do not match. Please check the bam files.")

                f_out.write(f"{forward_read.query_name}\t{forward_read.reference_name}\t{forward_read.pos}\t{reverse_read.reference_name}\t{reverse_read.pos}\t{'+' if forward_read.flag == 0 or forward_read.flag == 256 else '-'}\t{'+' if reverse_read.flag == 0 or forward_read.flag == 256 else '-'}\n")

        f_out.close()
        bam_for_handler.close()
        bam_rev_handler.close()
        bam_for_handler_rescued.close()
        bam_rev_handler_rescued.close()

    logger.info(f"Pairs file successfully created in {output_path}")


def build_matrix(bins : str = "fragments_fixed_sizes.txt", pairs : str = "group1.pairs", mode : bool = False, balance : bool = False, cpus : int = 8, output_dir : str = None) -> None:
    """
    Take table of bins and .pairs file and build a matrix in .cool format.

    Parameters
    ----------
    bins : str
        Path to the bin file.
    pairs : str, optional
        Path to pairs file, by default "group1.pairs"
    mode : bool, optional
        Choose wether the mode is rescued or un-rescued to construct associated .cool file, by default False
    balance : bool, optional
        Set wether or not to balance the matrix, by default False
    output_dir : str, optional
        Path to the folder where to save the cooler matrix file, by default None, by default None
    """

    output_path = Path(output_dir)

    if not output_path.exists():
            
        raise ValueError(f"Output path {output_path} does not exist. Please provide existing ouput path.")
    
    pairs_path = output_path / Path(pairs)

    if not pairs_path.is_file():
            
        raise ValueError(f"Pairs file {pairs_path} not found. Please provide existing pairs file.")
    
    bins_path = output_path / Path(bins)

    if not bins_path.is_file():
                
        raise ValueError(f"Bins file {bins_path} not found. Please provide existing bins file.")
    
    if not mode:

        cool_path = output_path / "unrescued_map.cool"

        cooler_cmd = f"cooler cload pairs --zero-based -c1 2 -p1 3 -p2 4 -c2 5 {bins_path} {pairs_path} {cool_path}"

    elif mode:

        pairs_path = output_path /"all_group.pairs"

        cool_path = output_path / "rescued_map.cool"

        cooler_cmd = f"cooler cload pairs --zero-based -c1 2 -p1 3 -p2 4 -c2 5 {bins_path} {pairs_path} {cool_path}"

    balance_cmd = f"cooler balance --nproc {cpus} {cool_path}"

    sp.run(cooler_cmd, shell=True)

    if balance:
        sp.run(balance_cmd, shell=True)

    logger.info(f"Cooler matrix successfully created in {output_path}")

def load_dictionary(dictionary : str = None) -> dict:
    """
    Load dictionary save into numpy format (.npy).

    Parameters
    ----------
    dictionary : str, optional
        Path to a the dictionary saved in numpy (.npy) format, by default None

    Returns
    -------
    dict
        Python native dictionary
    """
    try:

        return np.load(dictionary, allow_pickle=True).item()
    
    except:

        return np.load(dictionary, allow_pickle=True)

def load_cooler(matrix : str = None) -> cooler.Cooler:
    """
    Load cooler matrix.

    Parameters
    ----------
    matrix : str, optional
        Path to a cooler matrix, by default None

    Returns
    -------
    cooler.Cooler
        Cooler matrix object
    """

    return cooler.Cooler(matrix.as_posix())

def merge_predictions(output_dir : str = None, clean : bool = True, stage = "prediction", cpus : int = 1) -> None:
    """
    Merge predictions of all chunks of ambiguous reads predictions.

    Parameters
    ----------
    output_dir : str, optional
        Path to the folder where to save the fused alignment file, by default None
    clean : bool, optional
        Set weither or not to remove temporary chunks, by default True
    stage : str, optional
        Set the stage of the merging, by default "prediction". Can be "prediction" or "classification"
    cpus : int, optional
        Number of cpus to use for the merging, by default 1
    """

    logger.info(f"Start merging predictions")
    if output_dir is None:
        output_path = Path(getcwd())

    else : 
        output_path = Path(output_dir)

    if stage == "prediction":

        forward_alignment_chunk_files = sorted(glob(str(output_path / "forward_*_predicted.bam")))
        reverse_alignment_chunk_files = sorted(glob(str(output_path / "reverse_*_predicted.bam")))

        forward_merge_cmd = f"samtools merge -f -n --threads {cpus} {output_path / 'group2.1.rescued.bam'} {' '.join(forward_alignment_chunk_files)}"
        reverse_merge_cmd = f"samtools merge -f -n --threads {cpus} {output_path / 'group2.2.rescued.bam'} {' '.join(reverse_alignment_chunk_files)}"

        # Launch merge
        sp.run(forward_merge_cmd, shell=True)

        # Launch merge
        sp.run(reverse_merge_cmd, shell=True)

        if clean:

            for forward_chunk, reverse_chunk in zip(forward_alignment_chunk_files, reverse_alignment_chunk_files):

                Path(forward_chunk).unlink()
                Path(reverse_chunk).unlink()

    elif stage == "classification":

        forward_unaligned_chunk_files = sorted(glob(str(output_path / "group_*_0.1.bam")))
        reverse_unaligned_chunk_files = sorted(glob(str(output_path / "group_*_0.2.bam")))
        forward_aligned_chunk_files = sorted(glob(str(output_path / "group_*_1.1.bam")))
        reverse_aligned_chunk_files = sorted(glob(str(output_path / "group_*_1.2.bam")))
        forward_multi_aligned_chunk_files = sorted(glob(str(output_path / "group_*_2.1.bam")))
        reverse_multi_aligned_chunk_files = sorted(glob(str(output_path / "group_*_2.2.bam")))


        forward_group0_merge_cmd = f"samtools merge -f -n --threads {cpus} {output_path / 'group0.1.bam'} {' '.join(forward_unaligned_chunk_files)}"
        reverse_group0_merge_cmd = f"samtools merge -f -n --threads {cpus} {output_path / 'group0.2.bam'} {' '.join(reverse_unaligned_chunk_files)}"
        forward_group1_merge_cmd = f"samtools merge -f -n --threads {cpus} {output_path / 'group1.1.bam'} {' '.join(forward_aligned_chunk_files)}"
        reverse_group1_merge_cmd = f"samtools merge -f -n --threads {cpus} {output_path / 'group1.2.bam'} {' '.join(reverse_aligned_chunk_files)}"
        forward_group2_merge_cmd = f"samtools merge -f -n --threads {cpus} {output_path / 'group2.1.bam'} {' '.join(forward_multi_aligned_chunk_files)}"
        reverse_group2_merge_cmd = f"samtools merge -f -n --threads {cpus} {output_path / 'group2.2.bam'} {' '.join(reverse_multi_aligned_chunk_files)}"

        # Unaligned reads
        sp.run(forward_group0_merge_cmd, shell=True)
        sp.run(reverse_group0_merge_cmd, shell=True)

        # Aligned once reads
        sp.run(forward_group1_merge_cmd, shell=True)
        sp.run(reverse_group1_merge_cmd, shell=True)

        # Multi-aligned reads
        sp.run(forward_group2_merge_cmd, shell=True)
        sp.run(reverse_group2_merge_cmd, shell=True)

        logger.info(f"Groups successfully merged in {output_path}")

        if clean:

            for forward_chunk, reverse_chunk in zip(forward_unaligned_chunk_files, reverse_unaligned_chunk_files):

                Path(forward_chunk).unlink()
                Path(reverse_chunk).unlink()

            for forward_chunk, reverse_chunk in zip(forward_aligned_chunk_files, reverse_aligned_chunk_files):

                Path(forward_chunk).unlink()
                Path(reverse_chunk).unlink()
    
            for forward_chunk, reverse_chunk in zip(forward_multi_aligned_chunk_files, reverse_multi_aligned_chunk_files):

                Path(forward_chunk).unlink()
                Path(reverse_chunk).unlink()

    elif stage == "benchmark":

        forward_out_chunk_files = sorted(glob(str(output_path / "chunk_for_*.out.bam")))
        reverse_out_chunk_files = sorted(glob(str(output_path / "chunk_rev_*.out.bam")))

        forward_out_merge_cmd = f"samtools merge -f -n --threads {cpus} {output_path / 'group1.1.out.bam'} {' '.join(forward_out_chunk_files)}"
        reverse_out_merge_cmd = f"samtools merge -f -n --threads {cpus} {output_path / 'group1.2.out.bam'} {' '.join(reverse_out_chunk_files)}"

        # Out reads
        sp.run(forward_out_merge_cmd, shell=True)
        sp.run(reverse_out_merge_cmd, shell=True)

        if clean:

            for forward_chunk, reverse_chunk in zip(forward_out_chunk_files, reverse_out_chunk_files):

                Path(forward_chunk).unlink()
                Path(reverse_chunk).unlink()

        logger.info(f"Out reads successfully merged in {output_path}")

def tidy_folder(output_dir : str = None) -> None:
    """
    Tidy all the files in the output folder.

    Parameters
    ----------
    output_dir : str, optional
        Path to the folder where to save the fused alignment file, by default None
    """ 

    if output_dir is None:
        output_path = Path(getcwd())

    else : 

        output_path = Path(output_dir)

    # Tidy folder
    files = [p for  p in output_path.glob("*")]

    for file in files :

        if Path(file).suffix == ".bt2l":

            Path(file).rename(output_path / "index" / Path(file).name)

        if Path(file).suffix == ".bam":

            Path(file).rename(output_path / "alignments" / Path(file).name)

        elif Path(file).suffix in [".npy", ".bed", ".bedgraph", ".bw"]:

            Path(file).rename(output_path / "statistics" / Path(file).name)

        elif Path(file).suffix == ".pairs":

            Path(file).rename(output_path / "contacts" / "pairs" / Path(file).name)

        elif Path(file).suffix == ".cool":

            Path(file).rename(output_path / "contacts" / "matrices" / Path(file).name)

        elif Path(file).suffix == ".pdf" or Path(file).suffix == ".svg":

            Path(file).rename(output_path / "plots" / Path(file).name)

    # Tidy plots
    plot_files = [p for  p in (output_path / "plots").glob("*.pdf")]

    for file in plot_files :

        if "density" in Path(file).name:

            Path(file).rename(output_path / "plots" / "densities" / Path(file).name)

        elif "coverage" in Path(file).name:

            Path(file).rename(output_path / "plots" / "coverages" / Path(file).name)
        
        elif "patterns" in Path(file).name:

            Path(file).rename(output_path / "plots" / "ps" / Path(file).name)

        elif Path(file).name.startswith("chr"):

            Path(file).rename(output_path / "plots" / "contact_maps" / Path(file).name)


