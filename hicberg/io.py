from os import getcwd, mkdir
from os.path import join
from pathlib import Path

import glob
import shutil

import subprocess as sp

import numpy as np
import cooler


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

    

def build_pairs():
    pass

def build_matrix():
    pass

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

    return np.load(dictionary, allow_pickle=True).item()

def load_cooler():
    pass

def merge_predictions():
    pass

