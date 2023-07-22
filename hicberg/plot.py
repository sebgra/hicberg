from os import getcwd
from os.path import join
from pathlib import Path

from itertools import product

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

import cooler
import cooltools
import bioframe as bf
import pysam as ps

from hicberg.io import load_dictionary, load_cooler
from hicberg import logger


DIST_FRAG = "dist.frag.npy"
XS = "xs.npy"
COVERAGE = "coverage.npy"
D1D2 = "d1d2.npy"
UNCUTS = "uncuts.npy"
WEIRDS = "weirds.npy"
LOOPS = "loops.npy"
TRANS_PS = "trans_ps.npy"
CLR = "unrescued_map.cool"

def plot_laws(output_dir : str = None) -> None:
    """
    Plot P(s) patterns laws
    
    Parameters
    ----------
    output_dir : str, optional
        Path to the folder where to save the plot, by default None, by default None.
    """

    if output_dir is None:
        output_path = Path(getcwd())

    else : 

        output_path = Path(output_dir)

    # reload dictionaries

    xs = load_dictionary(output_path / XS)
    dist_frag = load_dictionary(output_path / DIST_FRAG)
    weirds = load_dictionary(output_path / WEIRDS)
    uncuts = load_dictionary(output_path / UNCUTS)
    loops = load_dictionary(output_path / LOOPS)

    for chromosome in xs.keys():

        plt.figure(figsize = (10, 10))

        plt.loglog(xs[chromosome], weirds[chromosome], "o", label="++/--")
        plt.loglog(xs[chromosome], uncuts[chromosome], "o", label="+-")
        plt.loglog(xs[chromosome], loops[chromosome], "o", label="-+")
        plt.title(f"Distribution of weirds, uncuts and loops events across chr11")
        plt.xlabel("Logarithmic binned genomic distances")
        plt.ylabel("Number of events")
        plt.grid()
        plt.legend()
        plt.savefig(output_path / f"patterns_distribution_{chromosome}.pdf", format = "pdf")
        plt.close()
     
    logger.info(f"Saved plots of patterns at : {output_path}")


def plot_trans_ps(output_dir : str = None) -> None:
    """
    Plot P(s) patterns laws
    
    Parameters
    ----------
    output_dir : str, optional
        Path to the folder where to save the plot, by default None, by default None.
    """

    if output_dir is None:
        output_path = Path(getcwd())

    else : 

        output_path = Path(output_dir)

    # reload dictionaries

    dist_frag = load_dictionary(output_path / DIST_FRAG)


    clr_unambiguous = load_cooler(output_path / CLR)
    map_table = clr_unambiguous.matrix(balance=False, as_pixels=True, join=True)[:]
    # chrm_sets = itertools.combinations_with_replacement(sorted(restriction_maps.keys()), 2)
    chrm_sets = product((dist_frag.keys()), repeat=2)

    t_ps = np.zeros((len(dist_frag.keys()) ** 2, 1))
    all_interaction_matrix = np.zeros((len(dist_frag.keys()) ** 2, 1))
    n_frags_matrix = np.zeros((len(dist_frag.keys()) ** 2, 1))

    trans_ps_dictionary = dict()

    # print(f"Size of t_ps : {t_ps.shape}")

    for idx, s in enumerate(chrm_sets):

        all_interactions = clr_unambiguous.matrix(balance=False).fetch(s[0], s[1]).sum()
        n_frags = len(dist_frag.get(str(s[0]))) * len(
            dist_frag.get(str(s[1]))
        )
        trans_ps_dictionary[s] = np.divide(all_interactions, np.multiply(n_frags, 4))

        # print(f"idx : {idx}, s : {s}, value : {np.divide(all_interactions, np.multiply(n_frags, 4))}, : all_interactions : {all_interactions}")

        t_ps[idx] = np.divide(all_interactions, np.multiply(n_frags, 4))
        all_interaction_matrix[idx] = all_interactions
        n_frags_matrix[idx] = n_frags

    t_ps = t_ps.reshape(
        (len(dist_frag.keys()), (len(dist_frag.keys())))
    )
    np.fill_diagonal(t_ps, np.nan)

    all_interaction_matrix = all_interaction_matrix.reshape(
        (len(dist_frag.keys()), (len(dist_frag.keys())))
    )
    np.fill_diagonal(all_interaction_matrix, np.nan)

    n_frags_matrix = n_frags_matrix.reshape(
        (len(dist_frag.keys()), (len(dist_frag.keys())))
    )
    np.fill_diagonal(n_frags_matrix, np.nan)

    plt.figure(figsize=(10, 10))

    plt.imshow(t_ps, cmap="Wistia", interpolation="None")
    plt.colorbar(fraction=0.046)
    plt.xticks(
        np.arange(len(list(dist_frag.keys()))),
        list(dist_frag.keys()),
        rotation="vertical",
    )
    plt.yticks(
        np.arange(len(list(dist_frag.keys()))),
        list(dist_frag.keys()),
    )
    plt.title("Pseudo P(s)")
    plt.savefig(output_path / f"pseudo_ps.pdf", format = "pdf")
    plt.close()

    logger.info(f"Saved pseudo P(s) of patterns at : {output_path}")

def plot_coverages(bins : int = 2000, output_dir : str = None ) -> None:
    """
    Plot coverages of chromosomes
    
    Parameters
    ----------
    bins : int, optional
        Size of the desired bin., by default 2000
    output_dir : str, optional
        Path to the folder where to save the plot, by default None, by default None.
    """
    
    if output_dir is None:
        output_path = Path(getcwd())

    else : 

        output_path = Path(output_dir)

    # reload dictionaries

    xs = load_dictionary(output_path / XS)
    coverage = load_dictionary(output_path / COVERAGE)


    for chromosome in xs.keys():

        
        plt.figure()
        plt.plot(coverage[chromosome], label="Covering smoothed")        
        plt.title(f"Covering across {chromosome} - bins of {bins} bp")
        plt.xlabel(f"Bin number")
        plt.ylabel("Number of reads")
        plt.legend()
        plt.grid()
        plt.savefig(output_path / f"coverage_{chromosome}.pdf", format = "pdf")
        plt.close()

    logger.info(f"Saved coverages at : {output_path}")



def plot_hic_matrix():
    pass

def plot_couple_repartition(forward_bam_file : str = "group2.1.rescued.bam", reverse_bam_file : str = "group2.1.rescued.bam",  output_dir : str = None ) -> None:
    """
    Plot read couples sizes distribution

    Parameters
    ----------
    forward_bam_file : str, optional
        Path to forward .bam alignment file, by default 1.sorted.bam
    reverse_bam_file : str, optional
        Path to reverse .bam alignment file, by default 2.sorted.bam
        Minimal read quality under which a Hi-C read pair will not be kept, by default 30
    output_dir : str, optional
        Path to the folder where to save the classified alignment files, by default None
    """

    if output_dir is None:
        output_path = Path(getcwd())

    else : 

        output_path = Path(output_dir)

    merged_forward_alignment_path = output_path / forward_bam_file
    merged_reverse_alignment_path = output_path / reverse_bam_file

    merged_forward_alignment_file_handler = ps.AlignmentFile(merged_forward_alignment_path, "rb")
    merged_reverse_alignment_file_handler = ps.AlignmentFile(merged_reverse_alignment_path, "rb")

     # Get the number of possible couples
    couple_lenght = list()

    for forward_read, reverse_read in zip(merged_forward_alignment_file_handler, merged_reverse_alignment_file_handler):

        couple_lenght.append(forward_read.get_tag("XL") * reverse_read.get_tag("XL"))

    _, bins_edges = np.histogram(couple_lenght, bins=max(couple_lenght))

    plt.figure()
    plt.vlines(
        x=np.mean(couple_lenght),
        ymin=0,
        ymax=max(_),
        color="red",
        label="mean",
        linestyles="dashed",
    )
    plt.vlines(
        x=np.median(couple_lenght),
        ymin=0,
        ymax=max(_),
        color="green",
        label="median",
        linestyles="dashed",
    )
    plt.vlines(
        x=np.percentile(couple_lenght, 99),
        ymin=0,
        ymax=max(_),
        color="purple",
        label="99 percentile",
        linestyles="dashed",
    )
    plt.loglog(_)
    plt.xlim(
        (2, (np.percentile(couple_lenght, 99) + np.std(couple_lenght)).astype(int))
    )
    plt.xlabel("Size of the set of possible couple")
    plt.ylabel("Number of couples")
    plt.title("Distribution of set of potential couple sizes")
    plt.legend()

    plt.savefig(
        join(output_path / f"Couple_sizes_distribution.pdf"),
        format="pdf",
    )
    plt.close()
    plt.close()

    logger.info(f"Saved couple size distribution at : {output_path}")

def plot_matrix_detailled():
    pass

def plot_matrix():
    pass

def plot_differences():
    pass

def plot_probability_maps():
    pass

def show_results():
    pass