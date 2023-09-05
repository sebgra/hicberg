from os import getcwd
from os.path import join
from pathlib import Path

from itertools import product, combinations

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
RESTRICTION_MAP = "restriction_map.npy"
DENSITY_MAP = "density_map.npy"


def plot_density(output_dir : str = None) -> None:
    """
    Plot density maps
    
    Parameters
    ----------
    output_dir : str, optional
        Path to the folder where to save the plots (one plot per chromosome couple), by default None.
    """

    if output_dir is None:
        output_path = Path(getcwd())

    else : 
        output_path = Path(output_dir)

    # reload dictionaries
    density_map = load_dictionary(output_path / DENSITY_MAP)

    for chromosome_couple in density_map.keys():

        matrix = density_map[chromosome_couple]

        plt.figure(figsize=(10, 10))
        plt.imshow(np.log10(matrix), cmap="afmhot_r", vmin = 0.0, vmax = 3.5)
        plt.title(f"Contatct density for  {chromosome_couple}")
        plt.savefig(output_path / f"density_{chromosome_couple}.pdf", format = "pdf")
        plt.close()

    logger.info(f"Saved plots of densities at : {output_path}")

def plot_benchmark(original_matrix : str = None, depleted_matrix : str = None, rescued_matrix : str = None, chromosomes : list[str] = None, output_dir : str = None) -> None:
    """
    Plot benchmark results (original, depleted and rescued matrices with associated log ratios). One plot per chromosome.

    Parameters
    ----------
    original_matrix : str, optional
        Path to the original matrix, by default None
    rescued_matrix : str, optional
        Path to the rescued matrix (reattributed reads), by default None
    chromosomes : list[str], optional
        List of chromosomes to plot, by default None
    output_dir : str, optional
        Path to where to save plots, by default None
    """

    if output_dir is None:  # if no output directory is provided, save in current directory     
        output_path = Path(getcwd())

    else : 

        output_path = Path(output_dir)    

    # define paths
    original_matrix_path = output_dir / original_matrix
    depleted_matrix_path = output_dir / depleted_matrix
    rescued_matrix_path = output_dir / rescued_matrix

    if not original_matrix_path.is_file():
        raise FileNotFoundError(f"Original matrix not found at {original_matrix_path}. Please provide a valid path.")

    if not depleted_matrix_path.is_file():
        raise FileNotFoundError(f"Depleted matrix not found at {depleted_matrix_path}. Please provide a valid path.")
    if not rescued_matrix_path.is_file():
        raise FileNotFoundError(f"Rescued matrix not found at {rescued_matrix_path}. Please provide a valid path.")
    
    
    # Relaod matricies

    original_matrix = load_cooler(original_matrix_path)
    depleted_matrix = load_cooler(depleted_matrix_path)
    rescued_matrix = load_cooler(rescued_matrix_path)

    print(f"Chromosomes : {chromosomes}")

    print(f"combinations :{list(combinations(chromosomes, 2))}")

    for chrm in chromosomes : 

        ori_matrix = original_matrix.matrix(balance=False).fetch(chrm)
        dep_matrix = depleted_matrix.matrix(balance=False).fetch(chrm)
        res_matrix = rescued_matrix.matrix(balance=False).fetch(chrm)
        ratio = np.divide(
            res_matrix,
            ori_matrix,
            out=np.ones(res_matrix.shape),
            where=ori_matrix != 0,
        )
        log_ratio = np.log10(ratio)

        # TODO : Adjust log non log and exponent
        plt.figure(figsize=(10, 10))
        plt.subplot(141)
        plt.imshow(ori_matrix ** 0.15, cmap = "afmhot_r", vmin = 0, vmax = np.max(ori_matrix ** 0.15))
        plt.title(f"Original map - {chrm}")
        plt.subplot(142)
        plt.imshow(dep_matrix ** 0.15, cmap = "afmhot_r", vmin = 0, vmax = np.max(ori_matrix ** 0.15))
        plt.title(f"Depleted map - {chrm}")
        plt.subplot(143)
        plt.imshow(res_matrix ** 0.15, cmap = "afmhot_r", vmin = 0, vmax = np.max(ori_matrix ** 0.15))
        plt.title(f"Rescued map - {chrm}")
        plt.subplot(144)
        plt.imshow(log_ratio, cmap = "bwr" , vmin = -1, vmax = 1) 
        plt.title(f"Log ratio - {chrm}")
        plt.colorbar(fraction=0.046)
        plt.savefig(output_path / f"benchmark_{chrm}.pdf", format = "pdf")
        plt.close()


def plot_d1d2(output_dir : str = None) -> None:
    """
    Plot d1d2 law
    
    Parameters
    ----------
    output_dir : str, optional
        Path to the folder where to save the plot, by default None, by default None.
    """

    if output_dir is None:
        output_path = Path(getcwd())

    else : 

        output_path = Path(output_dir)

    # reload dictionary
    d1d2 = load_dictionary(output_path / D1D2)

    histo, bins = np.histogram(d1d2, max(d1d2))

    plt.figure(figsize=(10, 10))
    plt.loglog(histo)
    plt.title("Log distribution of d1d2 distance")
    plt.xlabel("d1+d2")
    plt.ylabel("No. occurences")
    plt.savefig(output_path / f"d1d2.pdf", format = "pdf")
    plt.close()

    logger.info(f"Saved plots of d1d2 at : {output_path}")


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
        
        
        plt.figure(figsize=(10, 10))

        plt.loglog(xs[chromosome], weirds[chromosome], "o", label="++/--")
        plt.loglog(xs[chromosome], uncuts[chromosome], "o", label="+-")
        plt.loglog(xs[chromosome], loops[chromosome], "o", label="-+")
        plt.title(f"Distribution of weirds, uncuts and loops events across {chromosome}")
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
    # map_table = clr_unambiguous.matrix(balance=False, as_pixels=True, join=True)[:]
    # chrm_sets = itertools.combinations_with_replacement(sorted(restriction_maps.keys()), 2)
    chrm_sets = product((dist_frag.keys()), repeat=2)

    t_ps = np.zeros((len(dist_frag.keys()) ** 2, 1))
    all_interaction_matrix = np.zeros((len(dist_frag.keys()) ** 2, 1))
    n_frags_matrix = np.zeros((len(dist_frag.keys()) ** 2, 1))

    trans_ps_dictionary = dict()

    for idx, s in enumerate(chrm_sets):

        all_interactions = clr_unambiguous.matrix(balance=False).fetch(s[0], s[1]).sum()
        n_frags = len(dist_frag.get(str(s[0]))) * len(
            dist_frag.get(str(s[1]))
        )
        trans_ps_dictionary[s] = np.divide(all_interactions, np.multiply(n_frags, 4))


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

def plot_couple_repartition(forward_bam_file : str = "group2.1.rescued.bam", reverse_bam_file : str = "group2.2.rescued.bam",  output_dir : str = None ) -> None:
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
        Path to the folder where to save the plot, by default None
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

    plt.savefig(output_path / f"Couple_sizes_distribution.pdf",
        format="pdf",
    )
    plt.close()


    logger.info(f"Saved couple size distribution at : {output_path}")


def plot_matrix(unrescued_matrix : str = "unrescued_map.cool", rescued_matrix : str = "rescued_map.cool", restriction_map : str = "restriction_map.npy", genome : str = "", vmin : float = 0.0, vmax : float = 3.5, bins : int = 2000, output_dir : str = None) -> None:
    """
    Plot matrix with additional trackss

    Parameters
    ----------
    unrescued_matrix : str, optional
        Path to the unrescued map file, by default unrescued_map.cool
    rescued_matrix : str, optional
        Path to rescued map file, by default rescued_map.cool
    restriction_map : dict, optional
        Restriction map saved as a dictionary like chrom_name : list of restriction sites' position, by default dist.frag.npy
    genome : str, optional
        Path to the genome to digest, by default None, by default None
    vmin : float, optional
        Inferior limit for the colorscale, by default 0.0
    vmax : float, optional
        Superior limit for the colorscale, by default 3.5
    bins : int, optional
        Size of the desired bin., by default 2000
    output_dir : str, optional
        Path to the folder where to save the plot, by default None
    """

    if output_dir is None:

        output_path = Path(getcwd())

    else : 

        output_path = Path(output_dir)

    # Get the matrix
    unrescued_matrix = load_cooler(output_path / unrescued_matrix)
    rescued_matrix = load_cooler(output_path /rescued_matrix)

    genome_file = bf.load_fasta(genome, engine="pysam")

    restriction_map = load_dictionary(output_path / restriction_map)

    bins = unrescued_matrix.bins()[:]
    gc_cov = bf.frac_gc(bins[["chrom", "start", "end"]], genome_file)


    ### to make a list of chromosome start/ends in bins:
    chromstarts = []
    for i in rescued_matrix.chromnames:
        chromstarts.append(rescued_matrix.extent(i)[0])

    for i in rescued_matrix.chromnames:

        lower = rescued_matrix.extent(str(i))[0]
        upper = rescued_matrix.extent(str(i))[1]

        # Unrescued
        cis_coverage_unrescued, tot_coverage_unrescued = cooltools.coverage(
            unrescued_matrix, ignore_diags=False
        )

        # Rescued
        cis_coverage, tot_coverage = cooltools.coverage(
            rescued_matrix, ignore_diags=False
        )

        # Plot the matrix
        fig = plt.figure(figsize=(20, 20))
        gs = gridspec.GridSpec(2, 2, height_ratios=[10, 1], width_ratios=[1, 1])

        ax1 = plt.subplot(gs[0])
        divider1 = make_axes_locatable(ax1)
        cax1 = divider1.append_axes("right", size="5%", pad=0.1)
        im_unrescued = ax1.imshow(
            np.log10(unrescued_matrix.matrix(balance=False).fetch(i)), vmin = vmin, vmax = vmax,
            cmap = "afmhot_r",
        )
        fig.colorbar(im_unrescued, cax=cax1, label="corrected frequencies")
        ax1.set_title(
            f"Unrescued map of chromosome {i} \n binned at {int(rescued_matrix.binsize / 1000 )}kb",
            loc="center",
        )

        # Synchronize rescued and unrescued parts
        ax2 = plt.subplot(gs[1], sharex=ax1, sharey=ax1)

        # Rescued map
        divider2 = make_axes_locatable(ax2)
        cax2 = divider2.append_axes("right", size="5%", pad=0.1)
        im_rescued = ax2.imshow(
            np.log10(rescued_matrix.matrix(balance=False).fetch(i)), vmin = vmin, vmax = vmax,
            cmap = "afmhot_r",
        )
        fig.colorbar(im_rescued, cax=cax2, label="corrected frequencies")
        ax2.set_title(
            f"Rescued map of chromosome {i} \n binned at {int(unrescued_matrix.binsize / 1000 ) }kb",
            loc="center",
        )

        ax3 = divider1.append_axes("bottom", size="15%", pad=0.5, sharex=ax1)
        ax3.plot(tot_coverage_unrescued[lower:upper], label="total")
        ax3.plot(cis_coverage_unrescued[lower:upper], label="cis")
        ax3.set_ylabel("Coverage")
        ax3.legend(loc="lower left", bbox_to_anchor=(1, 0.5))
        ax3.set_xticks([])

        ax4 = divider1.append_axes("bottom", size="15%", pad=0.5, sharex=ax1)
        ax4.plot(list(gc_cov["GC"][lower:upper]), color="purple")
        ax4.set_ylabel("GC Content")

        ax5 = divider2.append_axes("bottom", size="15%", pad=0.5, sharex=ax2)
        ax5.plot(tot_coverage_unrescued[lower:upper], label="Before recovery")
        ax5.plot(tot_coverage[lower:upper], label="After recovery")
        ax5.set_xlim([0, len(unrescued_matrix.bins().fetch(str(i)))])
        ax5.set_ylabel("Coverage")
        ax5.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        ax5.set_xticks([])

        ax6 = divider2.append_axes("bottom", size="15%", pad=0.5, sharex=ax2)
        ax6.plot(list(gc_cov["GC"][lower:upper]), color="purple")
        ax6.set_ylabel("GC Content")

        plt.savefig(
            output_path / f"chr_{i}.pdf",
            format="pdf",
        )

        plt.close()

