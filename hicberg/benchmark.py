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

    # position = [int(p) for p in position.split(",")]
    # trans_chromosome = [c for c in trans_chromosome.split(",")]
    # trans_position = [int(tp) for tp in trans_position.split(",")]
    # strides = [int(s) for s in strides.split(",")]

    MODE = [[m] for m in mode.split(",")]

    CHROMOSOME = [c for c in chromosome.split(",")]

    TRANS_CHROMOSOME = [[str(t)] for t in trans_chromosome.split(",")]

    # TRANS_POSITION = [[int(p) for p in P.split(",")] for P in trans_position]

    TRANS_POSITION = [[int(p)] for p in position.split(",")]


    NB_BINS = bins
    POSITION = [[int(p) for p in str(position).split(",")]]
    BIN_SIZE = bins

    STRIDES = [[int(s) for s in strides.split(",")]]

    print(f"output_dir : {output_dir}")
    print(f"chromosome : {CHROMOSOME}")
    print(f"position : {POSITION}")
    print(f"trans_chromosome : {TRANS_CHROMOSOME}")
    print(f"trans_position : {TRANS_POSITION}")
    print(f"strides : {STRIDES}")
    print(f"mode : {MODE}")
    print(f"auto : {auto}")
    print(f"bins : {bins}")

    print('----- Creating grid search -----')

    # stride = stride

    # trans_coordinates = list(zip(*trans_chromosome, trans_position))
    trans_coordinates = zip(TRANS_CHROMOSOME, TRANS_POSITION)

    GS = list(product(MODE, CHROMOSOME, POSITION, trans_coordinates, STRIDES))

    # print(list(product(MODE, CHROMOSOME, POSITION, trans_coordinates, STRIDES)))

    # if len(trans_chromosome) != 0:

    #     GS =  product(mode, chromosome, position, strides, trans_coordinates)

    # else:

    #     GS = list(product(mode, chromosome, position, strides))

    # print(f"GS : {GS}")

    for idx, g in enumerate(GS):

        # print(f"{idx} : {g}")

        print(f"GS[{idx}] : {g}")

        # if len(TRANS_CHROMOSOME) != 0:

    


        #     MODE, CHROMOSOME, POSITION, STRIDE = g[:-1]
        #     TRANS_CHROMOSOME, TRANS_POSITION = g[-1]

        # # conver list of position from CLI from str to int
        #     TRANS_POSITION = [int(x) for x in TRANS_POSITION.split(",")]

        # else : 

        #     MODE, CHROMOSOME, POSITION, STRIDE  =  g



    return

