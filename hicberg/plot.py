from os import getcwd
from os.path import join
from pathlib import Path

from itertools import product, combinations

import numpy as np
from scipy.ndimage import gaussian_filter1d
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as plc
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

import cooler
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

#Matplotlib parameters
plt.style.use('seaborn-v0_8')  
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


def plot_density(output_dir: str = None) -> None:
    """
    Plot density maps

    Parameters
    ----------
    output_dir : str, optional
        Path to the folder where to save the plots (one plot per chromosome couple), by default None.
    """

    if output_dir is None:
        output_path = Path(getcwd())

    else:
        output_path = Path(output_dir)

    # reload dictionaries
    density_map = load_dictionary(output_path / DENSITY_MAP)

    for chromosome_couple in density_map.keys():

        matrix = density_map[chromosome_couple]

        cmap = plt.get_cmap("seismic")
        cmap.set_bad(color="black")
        plt.figure(figsize=(10, 10))
        plt.imshow(np.log10(matrix), cmap=cmap, vmin=-1, vmax=1)
        plt.title(f"Contact density for  {chromosome_couple}")
        plt.colorbar(fraction=0.046)
        plt.savefig(
            output_path / f"density_{chromosome_couple[0]}-{chromosome_couple[1]}.pdf",
            format="pdf",
        )
        plt.close()

    logger.info(f"Saved plots of densities at : {output_path}")


def plot_benchmark(
    original_matrix: str = None,
    depleted_matrix: str = None,
    rescued_matrix: str = None,
    chromosomes: list[str] = None,
    output_dir: str = None,
) -> None:
    """
    Plot benchmark results (original, depleted and rescued matrices with associated log ratios). One plot per chromosome.

    Parameters
    ----------
    original_matrix : str, optional
        Path to the original matrix, by default None
    rescued_matrix : str, optional
        Path to the rescued matrix (re-attributed reads), by default None
    chromosomes : list[str], optional
        List of chromosomes to plot, by default None
    output_dir : str, optional
        Path to where to save plots, by default None
    """

    if (
        output_dir is None
    ):  # if no output directory is provided, save in current directory
        output_path = Path(getcwd())

    else:
        output_path = Path(output_dir)

    chromosomes = chromosomes if type(chromosomes) == list else chromosomes.split()

    # define paths
    original_matrix_path = output_dir / original_matrix
    depleted_matrix_path = output_dir / depleted_matrix
    rescued_matrix_path = output_dir / rescued_matrix

    if not original_matrix_path.is_file():
        raise FileNotFoundError(
            f"Original matrix not found at {original_matrix_path}. Please provide a valid path."
        )

    if not depleted_matrix_path.is_file():
        raise FileNotFoundError(
            f"Depleted matrix not found at {depleted_matrix_path}. Please provide a valid path."
        )
    if not rescued_matrix_path.is_file():
        raise FileNotFoundError(
            f"Rescued matrix not found at {rescued_matrix_path}. Please provide a valid path."
        )

    # Relaod matricies
    original_matrix = load_cooler(original_matrix_path)
    depleted_matrix = load_cooler(depleted_matrix_path)
    rescued_matrix = load_cooler(rescued_matrix_path)

    for chrm in chromosomes:

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
        plt.subplot(221)
        plt.imshow(
            ori_matrix**0.15, cmap="afmhot_r", vmin=0, vmax=np.max(ori_matrix**0.15)
        )
        plt.title(f"Original map - {chrm}")
        plt.subplot(222)
        plt.imshow(
            dep_matrix**0.15, cmap="afmhot_r", vmin=0, vmax=np.max(ori_matrix**0.15)
        )
        plt.title(f"Depleted map - {chrm}")
        plt.subplot(223)
        plt.imshow(
            res_matrix**0.15, cmap="afmhot_r", vmin=0, vmax=np.max(ori_matrix**0.15)
        )
        plt.title(f"Rescued map - {chrm}")
        plt.subplot(224)
        plt.imshow(log_ratio, cmap="bwr", vmin=-1, vmax=1)
        plt.title(f"Log ratio - {chrm}")
        plt.colorbar(fraction=0.046)
        plt.savefig(output_path / f"benchmark_{chrm}.pdf", format="pdf")
        plt.close()


def plot_d1d2(output_dir: str = None) -> None:
    """
    Plot d1d2 law

    Parameters
    ----------
    output_dir : str, optional
        Path to the folder where to save the plot, by default None, by default None.
    """

    if output_dir is None:
        output_path = Path(getcwd())

    else:

        output_path = Path(output_dir)

    # reload dictionary
    d1d2 = load_dictionary(output_path / D1D2)
# Ensure d1d2 is an array for consistent indexing and filtering
    d1d2 = np.array(d1d2)

    # Filter out non-positive values from d1d2 for histogram on log scale
    # max(d1d2) as bin count can be very large if d1d2 values are large.
    # It's usually better to pick a number of bins or a more sophisticated binning strategy
    # for a log-log histogram. Let's use `ax.hist` with `log=True` directly.

    # Ensure min_data_val is positive for log scale
    min_data_val = np.min(d1d2[d1d2 > 0]) if np.any(d1d2 > 0) else 1
    max_data_val = np.max(d1d2)

    if min_data_val <= 0 or max_data_val <= min_data_val:
        logger.error("d1d2 data range is not suitable for a log-log plot. "
                     "Ensure values are positive and span a range.")
        return

    # Create bins for the histogram in log space
    num_bins = 50 # You can adjust this number
    log_bins = np.logspace(np.log10(min_data_val), np.log10(max_data_val), num_bins)


    # --- Plotting with ax objects ---
    fig, ax = plt.subplots(figsize=(10, 10)) # Create figure and axes

    # Plot the histogram directly.
    # `log=True` here means the Y-axis will be logarithmic (counts are logged).
    # The X-axis scale will be set explicitly later.
    counts, bins, patches = ax.hist(
        d1d2,
        bins=log_bins, # Use log-spaced bins for the x-axis
        color='skyblue', # Choose a suitable color
        edgecolor='black',
        alpha=0.7
    )

    # Set both x and y axes to a logarithmic scale
    ax.set_xscale('log')
    ax.set_yscale('log') # This makes the y-axis (counts) logarithmic

    # Set x-axis limits based on the log bins
    ax.set_xlim(min_data_val, max_data_val)


    # Set labels and title directly on ax
    ax.set_title("Distribution of d1d2 distance", fontsize=16) # Removed "Log" from title since axes clarify
    ax.set_xlabel("d1+d2", fontsize=14)
    ax.set_ylabel("Number of occurrences", fontsize=14) # Sticking to "Number of occurrences"

    # Add grid lines for both major and minor ticks on both axes
    ax.grid(True, which="major", ls="-", alpha=0.7)
    ax.grid(True, which="minor", ls="--", alpha=0.3)

    # Set tick parameters
    ax.tick_params(axis='both', which='major', labelsize=13)
    ax.tick_params(axis='both', which='minor', labelsize=10) # Smaller labels for minor ticks

    # Ensure tight layout to prevent elements from overlapping
    fig.tight_layout()

    # Save and close the figure
    fig.savefig(output_path / f"d1d2.pdf", format="pdf")
    plt.close(fig)

    logger.info(f"Saved plots of d1d2 at : {output_path}")


def plot_laws(output_dir: str = None) -> None:
    """
    Plot P(s) patterns laws.

    Parameters
    ----------
    output_dir : str, optional
        Path to the folder where to save the plot, by default None, by default None.
    """

    if output_dir is None:
        output_path = Path(getcwd())

    else:
        output_path = Path(output_dir)

    # reload dictionaries

    xs = load_dictionary(output_path / XS)
    weirds = load_dictionary(output_path / WEIRDS)
    uncuts = load_dictionary(output_path / UNCUTS)
    loops = load_dictionary(output_path / LOOPS)
    
    for chromosome in xs.keys():
        
        # Interpolate curves
        interpolated_loops = gaussian_filter1d(loops[chromosome], 1.1)
        interpolated_uncuts = gaussian_filter1d(uncuts[chromosome], 1.1)
        interpolated_weirds = gaussian_filter1d(weirds[chromosome], 1.1)
        
        fig, ax = plt.subplots(figsize=(8, 6)) # Create figure and axes object

        ax.loglog(xs[chromosome], interpolated_weirds, label="++/--", color = colors[0])
        ax.loglog(xs[chromosome], weirds[chromosome], "o", label="+-",markersize=5, color = colors[0])

        
        ax.loglog(xs[chromosome], interpolated_uncuts, label="+-", color = colors[1])
        ax.loglog(xs[chromosome], uncuts[chromosome], "o", label="+-",markersize=5, color = colors[1])

        ax.loglog(xs[chromosome], interpolated_loops, label="-+", color = colors[2])
        ax.loglog(xs[chromosome], loops[chromosome], "o", label="+-",markersize=5, color = colors[2])


        ax.set_title(
            f"Distribution of weirds, uncuts and loops events across {chromosome}",
            fontsize=16  # Title font size set to 16
        )
        ax.set_xlabel("Logarithmic binned genomic distances", fontsize=14)
        ax.set_ylabel("Number of events", fontsize=14)

        ax.tick_params(axis='both', which='major', labelsize=13)

        ax.grid(True, linestyle='--', alpha=0.7)
        ax.legend(fontsize=15) 
        plt.savefig(
            output_path / f"patterns_distribution_{chromosome}.pdf", format="pdf"
        )
        plt.close()
    
    logger.info(f"Saved plots of patterns at : {output_path}")

def plot_trans_ps(output_dir: str = None) -> None:
    """
    Plot P(s) patterns laws

    Parameters
    ----------
    output_dir : str, optional
        Path to the folder where to save the plot, by default None, by default None.
    """

    if output_dir is None:
        output_path = Path(getcwd())

    else:
        output_path = Path(output_dir)

    # reload dictionaries
    dist_frag = load_dictionary(output_path / DIST_FRAG)
    clr_unambiguous = load_cooler(output_path / CLR)
    chrm_sets = product((dist_frag.keys()), repeat=2)

    t_ps = np.zeros((len(dist_frag.keys()) ** 2, 1))
    all_interaction_matrix = np.zeros((len(dist_frag.keys()) ** 2, 1))
    n_frags_matrix = np.zeros((len(dist_frag.keys()) ** 2, 1))

    trans_ps_dictionary = dict()

    for idx, s in enumerate(chrm_sets):

        all_interactions = clr_unambiguous.matrix(balance=False).fetch(s[0], s[1]).sum()
        n_frags = len(dist_frag.get(str(s[0]))) * len(dist_frag.get(str(s[1])))
        trans_ps_dictionary[s] = np.divide(all_interactions, np.multiply(n_frags, 4))

        t_ps[idx] = np.divide(all_interactions, np.multiply(n_frags, 4))
        all_interaction_matrix[idx] = all_interactions
        n_frags_matrix[idx] = n_frags

    t_ps = t_ps.reshape((len(dist_frag.keys()), (len(dist_frag.keys()))))
    np.fill_diagonal(t_ps, np.nan)

    all_interaction_matrix = all_interaction_matrix.reshape(
        (len(dist_frag.keys()), (len(dist_frag.keys())))
    )
    np.fill_diagonal(all_interaction_matrix, np.nan)

    n_frags_matrix = n_frags_matrix.reshape(
        (len(dist_frag.keys()), (len(dist_frag.keys())))
    )
    np.fill_diagonal(n_frags_matrix, np.nan)

    fig, ax = plt.subplots(figsize=(10, 10))

    im = ax.imshow(t_ps, cmap="Wistia", interpolation="None")
    plt.colorbar(im, ax=ax, fraction=0.046)

    # Set x-ticks and labels directly on the axes
    ax.set_xticks(
        np.arange(len(list(dist_frag.keys())))
    )
    ax.set_xticklabels(
        list(dist_frag.keys()),
        rotation="vertical",
    )
    ax.set_yticks(
        np.arange(len(list(dist_frag.keys())))
    )
    ax.set_yticklabels(
        list(dist_frag.keys()),
    )
    ax.tick_params(axis='both', which='major', labelsize=14)

    ax.set_title("Pseudo P(s)", fontsize=16)
    ax.grid(False)
    plt.savefig(output_path / f"pseudo_ps.pdf", format="pdf")
    plt.close()

    logger.info(f"Saved pseudo P(s) of patterns at : {output_path}")

    
def plot_coverages(bins: int = 2000, output_dir: str = None) -> None:
    """
    Plot coverages of chromosomes with x-axis converted to genomic coordinates,
    and the last tick being the maximum chromosome size rounded to the 10k ceil.

    Parameters
    ----------
    bins : int, optional
        Size of the desired bin in base pairs (bp), by default 2000.
    output_dir : str, optional
        Path to the folder where to save the plot, by default None.
    """

    if output_dir is None:
        output_path = Path(getcwd())
    else:
        output_path = Path(output_dir)

    xs = load_dictionary(output_path / XS)
    coverage = load_dictionary(output_path / COVERAGE)

    for chromosome in coverage.keys():
        fig, ax = plt.subplots(figsize=(12, 6))

        bin_indices = np.arange(len(coverage[chromosome]))
        genomic_coords = bin_indices * bins

        ax.plot(genomic_coords, coverage[chromosome], label="Covering smoothed")

        ax.set_title(f"{chromosome} coverage - bins of {bins / 1_000:.1f} kb", fontsize=16)
        ax.set_xlabel("Genomic position (kb)", fontsize=14)
        ax.set_ylabel("Number of reads", fontsize=14)
        ax.legend(fontsize=12)
        ax.grid(True, linestyle='--', alpha=0.7)

        # Calculate the actual maximum genomic coordinate (end of the last bin)
        # If there are no bins, max_genomic_coord is 0.
        if len(coverage[chromosome]) > 0:
            actual_max_genomic_coord = (len(coverage[chromosome])) * bins
        else:
            actual_max_genomic_coord = 0

        round_to = 10_000
        rounded_max_genomic_coord = np.ceil(actual_max_genomic_coord / round_to) * round_to

        ax.set_xlim(0, rounded_max_genomic_coord)

        tick_interval_bp = 100_000
        major_tick_locs_bp = np.arange(0, rounded_max_genomic_coord + 1, tick_interval_bp)

        if rounded_max_genomic_coord not in major_tick_locs_bp:
            major_tick_locs_bp = np.append(major_tick_locs_bp, rounded_max_genomic_coord)
            major_tick_locs_bp = np.sort(np.unique(major_tick_locs_bp))

        ax.set_xticks(major_tick_locs_bp)

        major_tick_labels = [f"{int(loc / 1_000)} kb" for loc in major_tick_locs_bp]
        ax.set_xticklabels(major_tick_labels, rotation=45, ha='right', fontsize=13)

        ax.tick_params(axis='y', which='major', labelsize=13)

        fig.tight_layout()

        fig.savefig(output_path / f"coverage_{chromosome}.pdf", format="pdf")
        plt.close(fig)

    logger.info(f"Saved coverages at: {output_path}")


def plot_couple_repartition(
    forward_bam_file: str = "group2.1.rescued.bam",
    reverse_bam_file: str = "group2.2.rescued.bam",
    output_dir: str = None,
) -> None:
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

    else:
        output_path = Path(output_dir)

    merged_forward_alignment_path = output_path / forward_bam_file
    merged_reverse_alignment_path = output_path / reverse_bam_file

    merged_forward_alignment_file_handler = ps.AlignmentFile(
        merged_forward_alignment_path, "rb"
    )
    merged_reverse_alignment_file_handler = ps.AlignmentFile(
        merged_reverse_alignment_path, "rb"
    )

    # Get the number of possible couples
    couple_length = list()

    for forward_read, reverse_read in zip(
        merged_forward_alignment_file_handler, merged_reverse_alignment_file_handler
    ):

        couple_length.append(forward_read.get_tag("XL") * reverse_read.get_tag("XL"))

    # Convert to numpy array for efficient calculations
    couple_length = np.array(couple_length)

    # Calculate statistics
    mean_val = np.mean(couple_length)
    median_val = np.median(couple_length)
    percentile_99_val = np.percentile(couple_length, 99)

    # --- Plotting with ax objects ---
    fig, ax = plt.subplots(figsize=(12, 8))

    # --- Changed: Determine histogram bins for a linear x-axis ---
    # You can either let ax.hist determine optimal bins, or specify a number:
    num_linear_bins = 100 # A reasonable number of bins for linear scale
    # Or calculate max data value to set a range for bins:
    min_data_val = np.min(couple_length)
    max_data_val = np.max(couple_length)
    linear_bins = np.linspace(min_data_val, max_data_val, num_linear_bins)

    # Plot the histogram using ax.hist()
    # No changes needed for y-axis scale (it's linear by default)
    counts, bins, patches = ax.hist(
        couple_length,
        bins=linear_bins, # Use linear-spaced bins
        color=colors[0],
        edgecolor='black',
        alpha=0.7,
        label="Distribution"
    )

    # --- Removed: ax.set_xscale('log') ---
    # The x-axis will now be linear by default.

    # Get the maximum y-value (count) from the histogram for positioning vlines
    max_y_val = np.max(counts) if counts.size > 0 else 0

    # Plot vertical lines for mean, median, and 99th percentile
    ax.vlines(
        x=mean_val,
        ymin=0,
        ymax=max_y_val,
        color=colors[1],
        label=f"Mean ({mean_val:.2f})",
        linestyles="dashed",
        linewidth=2
    )
    ax.vlines(
        x=median_val,
        ymin=0,
        ymax=max_y_val,
        color=colors[2],
        label=f"Median ({median_val:.2f})",
        linestyles="dashed",
        linewidth=2
    )
    ax.vlines(
        x=percentile_99_val,
        ymin=0,
        ymax=max_y_val,
        color=colors[3],
        label=f"99th Percentile ({percentile_99_val:.2f})",
        linestyles="dashed",
        linewidth=2
    )

    # --- Adjusted: Set x-axis limits for a linear scale ---
    # Start from 0 or slightly below the min data value, extend to upper calculated limit.
    lower_xlim = 0 # Start from 0 for counts
    upper_xlim_calc = np.ceil(percentile_99_val + np.std(couple_length))
    final_upper_xlim = max(upper_xlim_calc, max_data_val) # Ensure it covers all data
    ax.set_xlim(lower_xlim, final_upper_xlim)


    # Set labels and title
    ax.set_xlabel("Number of possible pairs", fontsize=14)
    ax.set_ylabel("Count", fontsize=14)
    ax.set_title("Distribution of set of potential couple number", fontsize=16)

    # Set legend and grid
    ax.legend(fontsize=12)
    # --- Adjusted: Grid for linear axes (no 'which="minor"' for x-axis if not needed) ---
    ax.grid(True, which="major", ls="-", alpha=0.7)
    # Minor grid on linear x-axis often not as useful as on log, so removed or can be adjusted.
    # ax.grid(True, which="minor", ls="--", alpha=0.3)

    # Set tick parameters
    ax.tick_params(axis='both', which='major', labelsize=13)
    # --- Removed: X-axis minor tick label size setting (not as relevant for linear) ---
    # ax.tick_params(axis='x', which='minor', labelsize=10)
    ax.tick_params(axis='y', which='minor', left=False)

    fig.tight_layout()

    fig.savefig(output_path / "Couple_number_distribution.pdf", format="pdf")
    plt.close(fig)

    logger.info(f"Saved couple number distribution at : {output_path}")


def plot_matrix(
    unrescued_matrix: str = "unrescued_map.cool",
    rescued_matrix: str = "rescued_map.cool",
    restriction_map: str = "restriction_map.npy",
    genome: str = "",
    vmin: float = 0.0,
    vmax: float = 3.5,
    bins: int = 2000,
    output_dir: str = None,
) -> None:
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

    else:
        output_path = Path(output_dir)

    # Get the matrix
    unrescued_matrix = load_cooler(output_path / unrescued_matrix)
    rescued_matrix = load_cooler(output_path / rescued_matrix)

    genome_file = bf.load_fasta(genome, engine="pysam")
    restriction_map = load_dictionary(output_path / restriction_map)
    bins = unrescued_matrix.bins()[:]
    gc_cov = bf.frac_gc(bins[["chrom", "start", "end"]], genome_file)

    ### to make a list of chromosome start/ends in bins:

    for i in rescued_matrix.chromnames:

        lower = rescued_matrix.extent(str(i))[0]
        upper = rescued_matrix.extent(str(i))[1]

        # Unrescued
        coverage_unrescued = np.sum(
            np.tril(unrescued_matrix.matrix(balance=False).fetch(i)), axis=1
        )
        median_coverage = np.repeat(
            np.median(coverage_unrescued), coverage_unrescued.shape[0]
        )
        # Rescued
        coverage_rescued = np.sum(
            np.tril(rescued_matrix.matrix(balance=False).fetch(i)), axis=1
        )

        # Plot the matrix
        fig = plt.figure(figsize=(20, 20))
        gs = gridspec.GridSpec(2, 2, height_ratios=[10, 1], width_ratios=[1, 1])

        ax1 = plt.subplot(gs[0])
        divider1 = make_axes_locatable(ax1)
        cax1 = divider1.append_axes("right", size="5%", pad=0.1)
        im_unrescued = ax1.imshow(
            np.log10(unrescued_matrix.matrix(balance=False).fetch(i)),
            vmin=vmin,
            vmax=vmax,
            cmap="afmhot_r",
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
            np.log10(rescued_matrix.matrix(balance=False).fetch(i)),
            vmin=vmin,
            vmax=vmax,
            cmap="afmhot_r",
        )
        fig.colorbar(im_rescued, cax=cax2, label="corrected frequencies")
        ax2.set_title(
            f"Rescued map of chromosome {i} \n binned at {int(unrescued_matrix.binsize / 1000 ) }kb",
            loc="center",
        )

        ax3 = divider1.append_axes("bottom", size="15%", pad=0.5, sharex=ax1)
        ax3.plot(coverage_unrescued)
        ax3.plot(median_coverage, linestyle="--", color="black")
        ax3.set_ylabel("Coverage")
        ax3.set_xticks([])
        ax3.set_title("Natural coverage")

        ax4 = divider1.append_axes("bottom", size="15%", pad=0.5, sharex=ax1)
        ax4.plot(list(gc_cov["GC"][lower:upper]), color="purple")
        ax4.set_ylabel("GC Content")

        ax5 = divider2.append_axes("bottom", size="15%", pad=0.5, sharex=ax2)
        ax5.plot(coverage_unrescued, label="Before HiC-BERG")
        ax5.plot(coverage_rescued, label="After HiC-BERG")
        ax5.plot(median_coverage, linestyle="--", color="black")
        ax5.set_title("Enhanced coverage")
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


def plot_pattern_reconstruction(
    table: pd.DataFrame = None,
    original_cool: str = None,
    rescued_cool: str = None,
    chromosome: str = None,
    threshold: float = 0.0,
    case: str = "",
    output_dir: str = None,
) -> None:
    """
    Create a plot of pattern reconstruction quality.

    Parameters
    ----------
    table : pd.DataFrame, optional
        Table containing either true positives, false positives or false negatives patterns, by default None
    original_cool : str, optional
        Path to the original matrix in .cool format, by default None
    rescued_cool : str, optional
        Path to the rescued matrix in .cool format, by default None
    chromosome : str, optional
        Selected chromosome, by default None
    threshold : float, optional
        Threshold for pattern score significance, by default 0.0
    case : str, optional
        Mode to consider, either true positives, false positives or false negatives, by default ""
    output_dir : str, optional
        Path to save plots, by default None
    """

    if output_dir is None:

        output_path = Path(getcwd())

    else:
        output_path = Path(output_dir)

    original_matrix = load_cooler(original_cool).matrix(balance=False)
    rescued_matrix = load_cooler(rescued_cool).matrix(balance=False)

    bin_size = load_cooler(original_cool).info["bin-size"]

    fig, ax = plt.subplots()
    plt.title(f"Reconstructed pattern {chromosome}\n {case}")
    # Use imshow to add the first set of data to the plot
    img1 = ax.imshow(
        original_matrix.fetch(chromosome) ** 0.15,
        cmap="afmhot_r",
        vmin=0,
        vmax=np.max(rescued_matrix.fetch(chromosome) ** 0.15),
    )

    if table is not None:
        colormap = plt.colormaps["Blues"]  # 'plasma' or 'viridis'
        colors = colormap(table["score"])
        norm = plc.Normalize(vmin=0.0, vmax=1.0)
        # Create a divider for the existing axes instance
        divider = make_axes_locatable(ax)

        # Append axes to the right of the main axes.
        cax1 = divider.append_axes("right", size="5%", pad=0.1)

        # Add the colorbar to the figure
        cbar1 = fig.colorbar(img1, cax=cax1)

        sc = ax.scatter(
            x=table["start1"] // bin_size,
            y=table["start2"] // bin_size,
            s=40,
            linewidth=2,
            color="none",
            edgecolors=colors,
        )
        sm = plt.cm.ScalarMappable(cmap=colormap)

        # Append axes to the bottom of the main axes.
        cax2 = divider.append_axes("bottom", size="5%", pad=0.4)

        # Add the second colorbar to the figure
        cbar2 = fig.colorbar(
            sm,
            cax=cax2,
            orientation="horizontal",
        )
        cbar2.set_label(f"Pattern score - threshold : {threshold}")

    fig.savefig(
        str(output_path / f"pattern_{case.replace(' ', '')}_{chromosome}.pdf"),
        format="pdf",
    )
