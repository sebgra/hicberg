from os import path
from pathlib import Path

from glob import glob
import tempfile as tmpf
import multiprocessing as mp
import subprocess as sp

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

    position = [int(p) for p in position.split(",")]
    trans_chromosome = [c for c in trans_chromosome.split(",")]
    trans_position = [int(tp) for tp in trans_position.split(",")]
    strides = [int(s) for s in strides.split(",")]

    print(f"output_dir : {output_dir}")
    print(f"chromosome : {chromosome}")
    print(f"position : {position}")
    print(f"trans_chromosome : {trans_chromosome}")
    print(f"trans_position : {trans_position}")
    print(f"strides : {strides}")
    print(f"mode : {mode}")
    print(f"auto : {auto}")
    print(f"bins : {bins}")



    return

