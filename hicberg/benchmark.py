from os import path
from pathlib import Path

from glob import glob
import tempfile as tmpf
import multiprocessing as mp
import subprocess as sp

from itertools import product

from shutil import rmtree
import click
import logging

import hicberg.io as hio
import hicberg.utils as hut
# import hicberg.align as hal
import hicberg.statistics as hst
import hicberg.eval as hev

def benchmark(output_dir : str = None, chromosome : str = "", position : int = 0, trans_chromosome : str = None, trans_position : int = None, strides : list[int] = [], mode : str = "full", auto : int = None, bins : int = None):
    """
    AI is creating summary for benchmark

    Parameters
    ----------
    output_dir : str, optional
        [description], by default None
    chromosome : str, optional
        [description], by default ""
    position : int, optional
        [description], by default 0
    trans_chromosome : str, optional
        [description], by default None
    trans_position : int, optional
        [description], by default None
    strides : list[int], optional
        [description], by default []
    mode : str, optional
        [description], by default "full"
    auto : int, optional
        [description], by default None
    bins : int, optional
        [description], by default None
    """ 

    #TODO : complete code and docstring   

    print(f"output_dir : {output_dir}")
    print(f"chromosome : {chromosome}")
    print(f"position : {position}")
    print(f"trans_chromosome : {trans_chromosome}")
    print(f"trans_position : {trans_position}")
    print(f"strides : {strides}")
    print(f"mode : {mode}")
    print(f"auto : {auto}")
    print(f"bins : {bins}")

    print("------ Formating inputs ------")

    OTHER = [[chromosome, position, trans_chromosome, trans_position, strides, auto, bins]]

    print(f"OTHER : {OTHER}")

    print('----- Creating grid search -----')


    print(f"----Mode : {mode}----")

    # CHROMOSOME = [c for c in chromosome.split(",")]
    CHROMOSOME = chromosome

    if trans_chromosome is not None:
        TRANS_CHROMOSOME = [str(t) for t in trans_chromosome.split(",")]
    else : 
        TRANS_CHROMOSOME = None

    if trans_position is not None:
        TRANS_POSITION = [int(p) for p in trans_position.split(",")]

    else :
        TRANS_POSITION = None


    NB_BINS = bins
    # POSITION = [int(p) for p in str(position).split(",")]
    POSITION = position
    BIN_SIZE = bins

    STRIDES = [int(s) for s in strides.split(",")]

    for m in mode.split(","):

        print(f"mode : {m}")
        print(f"output_dir : {output_dir}")
        print(f"chromosome : {CHROMOSOME}")
        print(f"position : {POSITION}")
        print(f"trans_chromosome : {TRANS_CHROMOSOME}")
        print(f"trans_position : {TRANS_POSITION}")
        print(f"strides : {STRIDES}")
        print(f"auto : {auto}")

        hev.select_reads(position = POSITION, chromosome = CHROMOSOME, strides = STRIDES, trans_chromosome = TRANS_CHROMOSOME, trans_position = TRANS_POSITION, auto = auto, nb_bins = NB_BINS, output_dir = output_dir)
        break


    return

