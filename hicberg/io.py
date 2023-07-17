from os import getcwd, mkdir
from os.path import join
from pathlib import Path

import glob
import shutil

import subprocess as sp

import numpy as np
import cooler
import pysam as ps


def create_folder(sample_name : str  = None, output_dir : str = None) -> None:
    """
    Creates folder architecture to store results and intermediate files for the full HiCBERG pipeline.

    Parameters
    ----------
    samlpe_name : str
        Name of the folder to be created.
    output_dir : str
        Path where the folder will be created.
    """

    if sample_name is None:

        sample_name = "sample"

    if output_dir is None:

        folder_path = Path(getcwd(), sample_name)

    else:

        folder_path = Path(output_dir, sample_name)

    mkdir(folder_path)

    

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
        Choose weither the mode is rescued or unrescued to construct associated .pairs file, by default False
    output_dir : str, optional
        Path where the alignement files (.sam) should be stored, by default None
    """

    output_path = Path(output_dir)

    if not output_path.exists():

        raise ValueError(f"Output path {output_path} does not exist. Please provide existing ouput path.")

    if not mode:

        bam_for_path = Path(output_path / bam_for)
        bam_rev_path = Path(output_path / bam_rev)

        if not bam_for_path.exists():

            raise ValueError(f"Forward bam file {bam_for} not found")
        
        if not bam_rev_path.exists():
                
                raise ValueError(f"Reverse bam file {bam_rev} not found")
        
        bam_for_handler = ps.AlignmentFile(bam_for_path, "rb")
        bam_rev_handler = ps.AlignmentFile(bam_rev_path, "rb")

        with open(output_path / "group1.pairs", "w") as f_out:

            for forward_read, reverse_read in zip(bam_for_handler, bam_rev_handler):

                if forward_read.query_name != reverse_read.query_name:

                    raise ValueError(f"Forward and reverse reads do not match. Please check the bam files.")

                f_out.write(f"{forward_read.query_name}\t{forward_read.reference_name}\t{forward_read.pos}\t{forward_read.flag}\t{reverse_read.reference_name}\t{reverse_read.pos}\t{reverse_read.flag}\n")

        f_out.close()


    elif mode: 
            
        bam_for_path = Path(bam_for_rescued)
        bam_rev_path = Path(bam_rev_rescued)
        bam_for_path_rescued = Path(bam_for_rescued)
        bam_rev_path_rescued = Path(bam_rev_rescued)

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

            for forward_read, reverse_read in zip(bam_for_handler, bam_rev_handler):

                if forward_read.query_name != reverse_read.query_name:

                    raise ValueError(f"Forward and reverse reads do not match. Please check the bam files.")

                f_out.write(f"{forward_read.query_name}\t{forward_read.reference_name}\t{forward_read.pos}\t{forward_read.flag}\t{reverse_read.reference_name}\t{reverse_read.pos}\t{reverse_read.flag}\n")

            for forward_read, reverse_read in zip(bam_for_handler_rescued, bam_rev_handler_rescued):

                if forward_read.query_name != reverse_read.query_name:

                    raise ValueError(f"Forward and reverse reads do not match. Please check the bam files.")

                f_out.write(f"{forward_read.query_name}\t{forward_read.reference_name}\t{forward_read.pos}\t{forward_read.flag}\t{reverse_read.reference_name}\t{reverse_read.pos}\t{reverse_read.flag}\n")

        f_out.close()

    print(f"Pairs file successfully created in {output_path}")


    

def build_matrix(bins : str = "fragments_fixed_sizes.txt", pairs : str = "group1.pairs", mode : bool = False, output_dir : str = None) -> None:
    """
    Take table of bins and .pairs file and build a matrix in .cool format.

    Parameters
    ----------
    bins : str
        Path to the bin file.
    pairs : str, optional
        Path to pairs file, by default "group1.pairs"
    mode : bool, optional
        Choose weither the mode is rescued or unrescued to construct associated .cool file, by default False
    output_dir : str, optional
        th to the folder where to save the cooler matrix file, by default None, by default None
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

        cooler_cmd = f"cooler cload pairs --zero-based -c1 2 -p1 3 -p2 5 -c2 6 {bins_path} {pairs_path} {cool_path}"

    elif mode:

        pairs_path = output_path /"all_group.pairs"

        cool_path = output_path / "rescued_map.cool"

        cooler_cmd = f"cooler cload pairs --zero-based -c1 2 -p1 3 -p2 5 -c2 6 {bins_path} {pairs_path} {cool_path}"

    balance_cmd = f"cooler balance {cool_path}"

    sp.check_call(cooler_cmd, shell=True)
    sp.check_call(balance_cmd, shell=True)

    print(f"Cooler matrix successfully created in {output_path}")

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

def merge_predictions():
    pass

