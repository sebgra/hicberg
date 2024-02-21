from email.policy import default
import multiprocessing
from multiprocessing import Process
from functools import partial
import click
import hicberg.io as hio
import hicberg.utils as hut
import hicberg.align as hal
import hicberg.statistics as hst
import hicberg.pipeline as hpp
# import hicberg.benchmark as mbk
import hicberg.plot as hpl
import hicberg.benchmark as hbk

import ast

import numpy as np

from hicberg.version import __version__


class PythonLiteralOption(click.Option):
    def type_cast_value(self, ctx, value):
        try:
            return ast.literal_eval(value)
        except:
            raise click.BadParameter(value)
        


@click.group()
@click.version_option(version=__version__)
def cli(chain=True):
    """
    HiC-BERG: Rescuing of multi-reads in Hi-C maps.
    """
    ...


@click.command()
@click.argument('data', nargs = -1)
# @click.option("--genome", "-g", required = True, default = None, type = str, help = "Genome to perform analysis on.")
@click.option("--index", "-i", required = False, default = None, type = str, help = "Index of the genome.")
@click.option("--name", "-n", required = False, default = "sample", type = str, help = "Name of the analysis.")
# @click.option("--fq-for", required = True, default = None, type = str, help = "Forward fastq file to analyze.")
# @click.option("--fq-rev", required = True, default = None, type = str, help = "Reverse fastq file to analyze.")
@click.option("--rate", "-r", required = False, default = 1.0, type = float, help = "Rate to use for sub-sampling restriction map.")
@click.option("--cpus", "-t", required = False, default = 1, type = int, help = "Threads to use for analysis.")
@click.option("--kernel-size", "-K", required = False, default = 11, type = int, help = "Size of the gaussian kernel for contact density estimation.")
@click.option("--deviation", "-d", required = False, default = 0.5, type = float, help = "Standard deviation for contact density estimation.")
@click.option("--mode", "-m", required = False, default = "full", type = str, help = "Statistical model to use for ambiguous reads assignment.")
@click.option("--max-alignment", '-k', required = False, type = int, default = None, help = "Set the number of alignments to report in ambiguous reads case.")
@click.option("--sensitivity", "-s", required = False, type = click.Choice(["very-sensitive", "sensitive", "fast", "very-fast"]), default = "very-sensitive", help = "Set sensitivity level for Bowtie2")
@click.option("--bins", "-b", required = False, type = int, default = 2000, show_default = True, help = "Size of bins")
@click.option("--enzyme", "-e", required = False, type = str, multiple = True, help = "Enzymes to use for genome digestion.")
@click.option("--circular", "-c", required = False, type = str, default = "", help = "Name of the chromosome to consider as circular")
@click.option("--mapq", "-q", required = False, type = int, default = 35, help = "Minimum MAPQ to consider a read as valid")
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
@click.option("--start-stage", required = False, type = click.Choice(["fastq", "bam", "groups", "build", "stats", "rescue", "final"]), default = "fastq", help = "Stage to start the pipeline")
@click.option("--exit-stage", required = False, type = click.Choice(["None", "bam", "groups", "build", "stats", "rescue", "final"]), default = "None", help = "Stage to exit the pipeline")
@click.option("--force", "-f", is_flag = True, help = "Set if previous analysis files are deleted")
def pipeline_cmd(data, index, name, rate, mode, kernel_size, deviation, cpus, output, max_alignment, sensitivity, bins, enzyme, circular, mapq, start_stage, exit_stage, force):
    """
    Add documentation here


    """
    hpp.pipeline(genome = data[0], index = index, name = name, fq_for = data[1], fq_rev = data[2], output_dir = output, cpus = cpus, rate = rate, nb_chunks = 2 * cpus, mode = mode, kernel_size = kernel_size, deviation = deviation, max_alignment = max_alignment,  sensitivity = sensitivity, bins = bins, enzyme = enzyme, circular = circular, mapq = mapq, start_stage = start_stage, exit_stage = exit_stage, force = force)
    
@click.command()
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
@click.option("--force", "-f", is_flag = True, help = "Set if previous analysis files are deleted")
@click.option("--name", "-n", required = False, default = None, type = str, help = "Name of the output folder to create.")
def create_folder_cmd(output, force, name):
    hio.create_folder(sample_name=name, output_dir=output, force=force)

@click.command()
@click.argument('data', nargs = -1)
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
# @click.option("--genome", "-g", required = True, default = None, type = str, help = "Genome to perform analysis on.")
@click.option("--bins", "-b", required = False, type = int, default = 2000, help = "Size of bins")
def get_tables_cmd(data, bins, output):
    hut.get_chromosomes_sizes(genome = data[0], output_dir = output)
    hut.get_bin_table(bins = bins, output_dir = output)


@click.command()
@click.argument('data', nargs = -1)
@click.option("--index", "-i", required = False, default = None, type = str, help = "Index of the genome.")
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
# @click.option("--fq-for", required = True, default = None, type = str, help = "Forward fastq file to map.")
# @click.option("--fq-rev", required = True, default = None, type = str, help = "Reverse fastq file to map.")
@click.option("--sensitivity", "-s", required = False, type = click.Choice(["very-sensitive", "sensitive", "fast", "very-fast"]), default = "very-sensitive", help = "Set sensitivity level for Bowtie2")
@click.option("--max-alignment", '-k', required = False, type = int, default = None, help = "Set the number of alignments to report in ambiguous reads case.")
@click.option("--cpus", "-t", required = False, default = 1, type = int, help = "Threads to use for analysis.")
@click.option("--verbose", "-v", is_flag = True, help = "Set verbosity level.")
def alignment_cmd(data, index, max_alignment, sensitivity, output, cpus, verbose):

    if index is None:
        index = hal.hic_build_index(genome = data[0], output_dir = output, cpus = cpus, verbose = verbose)

    hal.hic_align(index = index, fq_for = data[1], fq_rev = data[2], sensitivity = sensitivity, max_alignment = max_alignment, output_dir = output, cpus = cpus, verbose = True)
    hal.hic_view(output_dir = output, cpus = cpus, verbose = verbose)
    hal.hic_sort(output_dir = output, cpus = cpus, verbose = verbose)


@click.command()
@click.option("--mapq", "-q", required = False, type = int, default = 35, help = "Minimum mapping quality to consider a read as valid")
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
def classify_cmd(mapq, output):
    hut.classify_reads(mapq = mapq, output_dir = output)


@click.command()
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
@click.option("--recover", "-r", required = False, default = False, is_flag = True, help = "Set if pairs are built after reads reassignment.")
def build_pairs_cmd(output, recover):
    hio.build_pairs(output_dir = output, mode = recover)

@click.command()
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
@click.option("--recover", "-r", required = False, default = False, is_flag = True, help = "Set if .cool matrix are built after reads reassignment.")
@click.option("--cpus", "-t", required = False, default = 1, type = int, help = "Threads to use for matrix building.")
def build_matrix_cmd(output, recover, cpus):
    hio.build_matrix(cpus = cpus, output_dir = output, mode = recover)

@click.command()
@click.argument('data', nargs = -1)
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
@click.option("--mode", "-m", required = False, default = "full", type = str, help = "Statistical model to use for ambiguous reads assignment.")
@click.option("--kernel-size", "-K", required = False, default = 11, type = int, help = "Size of the gaussian kernel for contact density estimation.")
@click.option("--deviation", "-d", required = False, default = 0.5, type = float, help = "Standard deviation for contact density estimation.")
# @click.option("--genome", "-g", required = True, default = None, type = str, help = "Genome to perform analysis on.")
@click.option("--rate", "-r", required = False, default = 1.0, type = float, help = "Rate to use for sub-sampling restriction map.")
@click.option("--enzyme", "-e", required = False, type = str, multiple = True, help = "Enzymes to use for genome digestion.")
@click.option("--circular", "-c", required = False, type = str, default = "", help = "Name of the chromosome to consider as circular")
@click.option("--bins", "-b", required = False, type = int, default = 2000, help = "Size of bins")
@click.option("--cpus", "-t", required = False, default = 1, type = int, help = "Threads to use for analysis.")
def statistics_cmd(data, mode, kernel_size, deviation,  rate, enzyme, circular, bins, output, cpus):
    restriction_map = hst.get_restriction_map(genome = data[0], enzyme = enzyme, output_dir = output)
    hst.get_dist_frags(genome = data[0], restriction_map = restriction_map, circular = circular, rate = rate, output_dir = output)
    hst.log_bin_genome(genome = data[0], output_dir = output)

    p1 = Process(target = hst.get_patterns(circular = circular, output_dir = output))
    p2 = Process(target = hst.generate_trans_ps(output_dir = output))
    p3 = Process(target = hst.generate_coverages(genome = data[0], bins = bins, output_dir = output))
    p4 = Process(target = hst.generate_d1d2(output_dir = output))

    # Launch processes
    for process in [p1, p2, p3, p4]:
        process.start()
        process.join()

    if mode in ["full", "density"]:

            hst.compute_density(kernel_size = kernel_size, deviation = deviation, threads  = cpus, output_dir  = output)
        

    

@click.command()
@click.argument('data', nargs = -1)
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
# @click.option("--genome", "-g", required = True, default = None, type = str, help = "Genome to perform analysis on.")
@click.option("--enzyme", "-e", required = False, type = str, multiple = True, help = "Enzymes to use for genome digestion.")
# @click.option("--nb-chunks", "-n", required = False, default = 1, type = int, help = "Number of chunks to split the alignment files into.")
@click.option("--mode", "-m", required = False, default = "full", type = str, help = "Statistical model to use for ambiguous reads assignment.")
@click.option("--cpus", "-t", required = False, default = 1, type = int, help = "Threads to use for analysis.")
def rescue_cmd(data, enzyme, mode, output, cpus):
    
    restriction_map = hst.get_restriction_map(genome = data[0], enzyme = enzyme, output_dir = output)
    hut.chunk_bam(nb_chunks = cpus, output_dir = output)
        
    # Get chunks as lists
    forward_chunks, reverse_chunks = hut.get_chunks(output)

    # Reattribute reads
    with multiprocessing.Pool(processes = cpus) as pool: # cpus

        results = pool.map(partial(hst.reattribute_reads, mode = mode, restriction_map = restriction_map, output_dir = output),
        zip(forward_chunks, reverse_chunks))
        pool.close()
        pool.join()

    hio.merge_predictions(output_dir = output, clean = True)


@click.command()
@click.argument('data', nargs = -1)
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
# @click.option("--genome", "-g", required = True, default = None, type = str, help = "Genome to perform analysis on.")
@click.option("--bins", "-b", required = False, default = 2000, type = int, help = "Size of bins")
def plot_cmd(data, bins, output):
    p1 = Process(target = hpl.plot_laws(output_dir = output))
    p2 = Process(target = hpl.plot_trans_ps(output_dir = output))
    p3 = Process(target = hpl.plot_coverages(bins = bins, output_dir = output))
    p4 = Process(target = hpl.plot_couple_repartition(output_dir = output))
    p5 = Process(target = hpl.plot_matrix(genome = data[0], output_dir = output))
    p6 = Process(target = hpl.plot_d1d2(output_dir = output))
    p7 = Process(target = hpl.plot_density, kwargs = dict(output_dir = output))

    # Launch processes
    for process in [p1, p2, p3, p4, p5, p6, p7]:
        process.start()
        process.join()

@click.command()
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
def tidy_cmd(output):
    hio.tidy_folder(output_dir = output)


@click.command()
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
@click.argument('data', nargs = -1)
# @click.option("--genome", "-g", required = True, default = None, type = str, help = "Genome to perform analysis on.")
@click.option("--chromosome", "-c", required = False, default = None, type = str, help = "Chromosome to get as source for duplication.")
@click.option("--position", "-p", required = False, default = None, type = int, help = "Position to get as source for duplication.")
@click.option("--trans-chromosome", "-C", required = False, default = None, type = str, help = "Chromosome to get as target for duplication.")
@click.option("--trans-position", "-P", required = False, default = None, type = str, help = "Position to get as target for duplication.")
@click.option("--bins", "-b", required = False, default = 1, type = int, help = "Number of bins to select from a genomic coordinates.")
@click.option("--strides", "-s", required = False, default = None, type = str, help = "Strides to apply from source genomic coordinates to define targets intervals.")
@click.option("--auto", "-a", required = False, default = None, type = int, help = "Automatically select auto intervals for duplication.")
@click.option("--kernel-size", "-K", required = False, default = 11, type = int, help = "Size of the gaussian kernel for contact density estimation.")
@click.option("--deviation", "-d", required = False, default = 0.5, type = float, help = "Standard deviation for contact density estimation.")
@click.option("--mode", "-m", required = False, default = "full", type = str, help = "Statistical model to use for ambiguous reads assignment.")
@click.option("--pattern", "-S", required = False, type = click.Choice(["loops", "borders", "hairpins", "-1"]), default = None, help = "Set pattern if benchmarking considering patterns")
@click.option("--threshold", "-t", required = False, type = float, default = 0.0, help = "Set pattern score threshold under which pattern are discarded")
@click.option("--jitter", "-j", required = False, type = int, default = 0, help = "Set jitter for pattern detection interval overlapping")
@click.option("--trend", "-T", is_flag = False, help = "Set if detrending of the contact map has to be performed")
@click.option("--top", "-k", required = False, type = int, default = 100, help = "Set the top k % of patterns to retain")
@click.option("--iterations", "-i", required = False, type = int, default = 3, help = "Set the number of iterations for benchmarking")
@click.option("--force", "-f", is_flag = True, help = "Set if previous analysis files have to be deleted")
@click.option("--cpus", required = False, default = 1, type = int, help = "Threads to use for analysis.")
def benchmark_cmd(data, chromosome, position, trans_chromosome, trans_position, bins, strides, auto, kernel_size, deviation, mode, pattern, threshold, jitter, trend, top, iterations, force, output, cpus):

    hbk.benchmark(output_dir = output, genome = data[0], chromosome = chromosome, position = position, trans_chromosome = trans_chromosome, trans_position = trans_position, strides = strides, mode = mode, force = force, bins = bins, auto = auto, kernel_size = kernel_size, deviation = deviation, pattern = pattern, threshold = threshold, jitter = jitter, trend = trend, top = top, iterations = iterations, cpus = cpus)

@click.command()
@click.argument('name', nargs = -1, metavar = "<name> <name>")
@click.option("--bins", "-b", required = False, type = int, default = 2000, show_default = True, help = "Size of bins")
def greet(bins, name):
    """Simple program that greets NAME."""
    click.echo(f'Hello, {name[0]} and {name[1]}!')
    print(f"bins : {bins}")




# Command group
cli.add_command(pipeline_cmd, name="pipeline")
cli.add_command(create_folder_cmd, name="create-folder")
cli.add_command(get_tables_cmd, name="get-tables")
cli.add_command(alignment_cmd, name="alignment")
cli.add_command(classify_cmd, name="classify")
cli.add_command(build_pairs_cmd, name="build-pairs")
cli.add_command(build_matrix_cmd, name="build-matrix")
cli.add_command(statistics_cmd, name="statistics")
cli.add_command(rescue_cmd, name="rescue")
cli.add_command(plot_cmd, name="plot")
cli.add_command(tidy_cmd, name="tidy")
cli.add_command(benchmark_cmd, name="benchmark")
cli.add_command(greet, name="greet")