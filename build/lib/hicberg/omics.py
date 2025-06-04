import subprocess as sp
from pathlib import Path
import numpy as np

import hicberg.io as hio

from hicberg import logger


def preprocess_pairs(pairs_file : str = "all_group.pairs", threshold : int = 1000, output_dir : str = None) -> None:
    """
    Preprocess pairs file to remove pairs that are not in the same chromosome or are greater than a threshold.
    Retain columns are : chromosome, start, end, count.

    Parameters
    ----------
    pairs_file : str, optional
        Path to the pairs file, by default "all_group.pairs"
    threshold : int, optional
        Threshold distance beyond which pairs will not be kept, by default 1000
    output_dir : str, optional
        Path where the formatted pairs will be saved, by default None
    """    

    output_dir_path = Path(output_dir)
    if not output_dir_path.is_dir():
        raise IOError(f"Output directory {output_dir} not found. Please provide a valid path.")

    pairs_path = Path(output_dir, pairs_file)

    if not pairs_path.is_file():
            
        raise IOError(f"Pairs file {pairs_path} not found. Please provide a valid path.")
    
    pairs_handler = open(pairs_path, "r")

    processed_pairs_path = Path(output_dir_path , "preprocessed_pairs.pairs")

    with open(processed_pairs_path, "w") as f_out:

        for line in pairs_handler:

            if line.startswith("#"):
                continue

            read_id, chromosome_for, position_for, chromosome_rev, position_rev, strand_for, strand_rev = line.split("\t")

            if chromosome_for != chromosome_rev or np.abs(int(position_rev) - int(position_for)) < threshold:
                continue

            else: 

                if int(position_for) < int(position_rev):

                    f_out.write(f"{chromosome_for}\t{position_for}\t{position_rev}\t1\n")
                else :

                    f_out.write(f"{chromosome_for}\t{position_rev}\t{position_for}\t1\n")

    pairs_handler.close()

    logger.info(f"Formated paris saved at {processed_pairs_path}")


def format_chrom_sizes(chromosome_sizes : str = "chromosome_sizes.npy", output_dir : str = None) -> None:
    """
    Format chromosome sizes to bed and txt format.
    - bed format : chrom, start, end
    - txt format : chrom, size
    Parameters
    ----------
    chrom_sizes : str, optional
        Path to chromosome sizes file (.npy), by default "chromosome_sizes.npy"
    output_dir : str, optional
        Path where the formatted chromosome sizes will be saved, by default None

    """    
    output_dir_path = Path(output_dir)
    if not output_dir_path.is_dir():
        raise IOError(f"Output directory {output_dir} not found. Please provide a valid path.")

    chrom_size_path = Path(output_dir, chromosome_sizes)

    if not chrom_size_path.is_file():
            
        raise IOError(f"Pairs file {chrom_size_path.name} not found. Please provide a valid path.")
    
    chrom_size = hio.load_dictionary(chrom_size_path)

    chrom_size_bed_path = Path(output_dir_path / "chromosome_sizes.bed")
    chrom_size_txt_path = Path(output_dir_path / "chromosome_sizes.txt")
    

    with open(chrom_size_bed_path, 'w') as f_out:

        for k, v in chrom_size.items():
            f_out.write(f'{k}\t0\t{v}\n')

    f_out.close()

    with open(chrom_size_txt_path, 'w') as f_out:

        for k, v in chrom_size.items():
            f_out.write(f'{k}\t{v}\n')

    f_out.close()

    logger.info(f"Formated chromosome sizes saved at {chrom_size_bed_path} and {chrom_size_txt_path}")

def get_bed_coverage(chromosome_sizes : str = "chromosome_sizes.bed", pairs_file : str = "preprocessed_pairs.pairs", output_dir : str = None) -> None:
    """
    Get bed coverage from pairs file (using bedtools).

    Parameters
    ----------
    chromosome_sizes : str, optional
        Path to chromsomes sizes files (.bed format), by default "chromosome_sizes.bed"
    pairs_file : str, optional
        Path to processed pairs files (columns : chrom, start, end, count), by default "preprocessed_pairs.pairs"
    output_dir : str, optional
        Path where the coverage (.bed) will be saved, by default None
    """    

    output_dir_path = Path(output_dir)
    if not output_dir_path.is_dir():
        raise IOError(f"Output directory {output_dir} not found. Please provide a valid path.")

    chrom_size_path = Path(output_dir, chromosome_sizes)

    pairs_path = Path(output_dir, pairs_file)

    if not chrom_size_path.is_file():
            
        raise IOError(f"Pairs file {chrom_size_path} not found. Please provide a valid path.")
    
    if not pairs_path.is_file():
                
        raise IOError(f"Pairs file {pairs_path} not found. Please provide a valid path.")
    
    bed_coverage_path = Path(output_dir_path , "coverage.bed")
    
    bedtools_cmd = f"bedtools coverage -a {str(chrom_size_path)} -b {str(pairs_path)} -d"

    with open(bed_coverage_path, "w") as f_out:

        sp.run(bedtools_cmd, shell=True, stdout=f_out)

    f_out.close()

    logger.info(f"Saved data coverage at {bed_coverage_path}")



def get_bedgraph(bed_coverage : str = "coverage.bed", output_dir : str = None) -> None:
    """
    Convert bed coverage to bedgraph format.
    Format is : chrom, start, end, count.
    Start and end are different by 1bp (end  = start + 1).

    Parameters
    ----------
    bed_coverage : str, optional
        Path to coverage (.bed), by default "coverage.bed"
    output_dir : str, optional
        Path where the coverage (.bedgraph) will be saved, by default None
    """    
    output_dir_path = Path(output_dir)
    if not output_dir_path.is_dir():
        raise IOError(f"Output directory {output_dir} not found. Please provide a valid path.")

    bed_coverage_path = Path(output_dir, bed_coverage)

    if not bed_coverage_path.is_file():
            
        raise IOError(f"Pairs file {bed_coverage_path.name} not found. Please provide a valid path.")
    
    bed_handler = open(bed_coverage_path, "r")

    bedgraph_coverage_path = Path(output_dir_path, "coverage.bedgraph")

    with open(bedgraph_coverage_path, "w") as f_out:

        for line in bed_handler:

            chromosome, start, end, index, count = line.split("\t")

            if end == index:
                continue

            f_out.write(f"{chromosome}\t{int(index)}\t{int(index) + 1}\t{count}")

    f_out.close()
    bed_handler.close()
    

def bedgraph_to_bigwig(bedgraph_file : str = "coverage.bedgraph", chromosome_sizes : str = "chromosome_sizes.txt", output_dir : str = None) -> None:
    """
    Convert bedgraph to bigwig format.

    Parameters
    ----------
    bedgraph_file : str, optional
        Path to coverage (.bedgraph), by default "coverage.bedgraph"
    chromosome_sizes : str, optional
        Path to chromosome sizes file (chrom_id, size), by default "chromosome_sizes.txt"
    output_dir : str, optional
        [description], by default None

    Raises
    ------
    IOError
        [description]
    IOError
        [description]
    IOError
        [description]
    """    
    output_dir_path = Path(output_dir)
    if not output_dir_path.is_dir():
        raise IOError(f"Output directory {output_dir} not found. Please provide a valid path.")

    bedgraph_coverage_path = Path(output_dir, bedgraph_file)
    if not bedgraph_coverage_path.is_file():
        raise IOError(f"Pairs file {bedgraph_coverage_path.name} not found. Please provide a valid path.")
    
    chromosome_sizes_path = Path(output_dir, chromosome_sizes)
    if not bedgraph_coverage_path.is_file():
        raise IOError(f"Pairs file {chromosome_sizes_path.name} not found. Please provide a valid path.")
    
    output_bigwig_path = Path(output_dir, "signal.bw")
    
    bedgraphtobigwig_cmd = f"bedGraphToBigWig {bedgraph_coverage_path} {chromosome_sizes_path} {output_bigwig_path}"

    sp.run([bedgraphtobigwig_cmd], shell = True)

    logger.info(f"Saved data in BigWig format at {output_bigwig_path}")