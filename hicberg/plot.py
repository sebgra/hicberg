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

from hicberg.io import load_dictionary, load_cooler
from hicberg import logger


DIST_FRAG = "dist.frag.npy"
XS = "xs.npy"
COVERAGE_DICO = "coverage.npy"
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


def plot_hic_matrix():
    pass

def plot_couple_repartition():
    pass

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