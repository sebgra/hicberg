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

import ast

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
@click.option("--genome", "-g", required = True, default = None, type = str, help = "Genome to perform analysis on.")
@click.option("--fq-for", required = True, default = None, type = str, help = "Forward fastq file to analyze.")
@click.option("--fq-rev", required = True, default = None, type = str, help = "Reverse fastq file to analyze.")
@click.option("--rate", "-r", required = False, default = 1.0, type = float, help = "Rate to use for subsampling restriction map.")
@click.option("--cpus", "-t", required = False, default = 1, type = int, help = "Threads to use for analysis.")
@click.option("--mode", "-m", required = False, default = "full", type = str, help = "Statistical model to use for ambiguous reads assignment.")
@click.option("--max-alignment", '-k', required = False, type = int, default = None, help = "Set the number of alignments to report in ambiguous reads case.")
@click.option("--sensitivity", "-s", required = False, type = click.Choice(["very-sensitive", "sensitive", "fast", "very-fast"]), default = "very-sensitive", help = "Set sensitivity level for Bowtie2")
@click.option("--bins", "-b", required = False, type = int, default = 2000, help = "Size of bins")
@click.option("--enzyme", "-e", required = False, type = str, multiple = True, help = "Enzymes to use for genome digestion.")
@click.option("--circular", "-c", required = False, type = str, default = "", help = "Name of the chromosome to consider as circular")
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
@click.option("--start-stage", required = False, type = click.Choice(["fastq", "bam", "stats", "pairs", "rescue", "cool"]), default = "fastq", help = "Stage to start the pipeline")
@click.option("--exit-stage", required = False, type = click.Choice(["None", "bam", "stats", "pairs", "rescue", "cool"]), default = "None", help = "Stage to exit the pipeline")
@click.option("--force", "-f", is_flag = True, help = "Set if previous analysis files are deleted")
def pipeline_cmd(genome, fq_for, fq_rev, rate, mode, cpus, output, max_alignment, sensitivity, bins, enzyme, circular, start_stage, exit_stage, force):
    hpp.pipeline(genome = genome, fq_for = fq_for, fq_rev = fq_rev, output_dir = output, cpus = cpus, rate = rate, nb_chunks = cpus, mode = mode, max_alignment = max_alignment,  sensitivity = sensitivity, bins = bins, enzyme = enzyme, circular = circular, start_stage = start_stage, exit_stage = exit_stage, force = force)
    
@click.command()
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
@click.option("--force", "-f", is_flag = True, help = "Set if previous analysis files are deleted")
@click.option("--name", "-n", required = False, default = None, type = str, help = "Name of the output folder to create.")
def create_folder_cmd(output, force, name):
    hio.create_folder(sample_name=name, output_dir=output, force=force)

@click.command()
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
@click.option("--genome", "-g", required = True, default = None, type = str, help = "Genome to perform analysis on.")
@click.option("--bins", "-b", required = False, type = int, default = 2000, help = "Size of bins")
def get_tables_cmd(genome, bins, output):
    hut.get_chromosomes_sizes(genome = genome, output_dir = output)
    hut.get_bin_table(bins = bins, output_dir = output)


@click.command()
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
@click.option("--genome", "-g", required = True, default = None, type = str, help = "Genome to perform analysis on.")
@click.option("--fq-for", required = True, default = None, type = str, help = "Forward fastq file to map.")
@click.option("--fq-rev", required = True, default = None, type = str, help = "Reverse fastq file to map.")
@click.option("--sensitivity", "-s", required = False, type = click.Choice(["very-sensitive", "sensitive", "fast", "very-fast"]), default = "very-sensitive", help = "Set sensitivity level for Bowtie2")
@click.option("--max-alignment", '-k', required = False, type = int, default = None, help = "Set the number of alignments to report in ambiguous reads case.")
@click.option("--cpus", "-t", required = False, default = 1, type = int, help = "Threads to use for analysis.")
@click.option("--verbose", "-v", is_flag = True, help = "Set verbosity level.")
def alignment_cmd(genome, fq_for, fq_rev, max_alignment, sensitivity, output, cpus, verbose):
    index = hal.hic_build_index(genome = genome, output_dir = output, cpus = cpus, verbose = verbose)
    hal.hic_align(genome = genome, index = index, fq_for = fq_for, fq_rev = fq_rev, sensitivity = sensitivity, max_alignment = max_alignment, output_dir = output, cpus = cpus, verbose = True)
    hal.hic_view(output_dir = output, cpus = cpus, verbose = verbose)
    hal.hic_sort(output_dir = output, cpus = cpus, verbose = verbose)


@click.command()
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
def classify_cmd(output):
    hut.classify_reads(output_dir = output)


@click.command()
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
@click.option("--recover", "-r", required = False, default = False, is_flag = True, help = "Set if pairs are bulit after reads reassignment.")
def build_pairs_cmd(output, recover):
    hio.build_pairs(output_dir = output, mode = recover)

@click.command()
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
@click.option("--recover", "-r", required = False, default = False, is_flag = True, help = "Set if .cool matrix are bulit after reads reassignment.")
def build_matrix_cmd(output, recover):
    hio.build_matrix(output_dir = output, mode = recover)

@click.command()
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
@click.option("--genome", "-g", required = True, default = None, type = str, help = "Genome to perform analysis on.")
@click.option("--rate", "-r", required = False, default = 1.0, type = float, help = "Rate to use for subsampling restriction map.")
@click.option("--enzyme", "-e", required = False, type = str, multiple = True, help = "Enzymes to use for genome digestion.")
@click.option("--circular", "-c", required = False, type = str, default = "", help = "Name of the chromosome to consider as circular")
@click.option("--bins", "-b", required = False, type = int, default = 2000, help = "Size of bins")
def statistics_cmd(genome, rate, enzyme, circular, bins, output):
    restriction_map = hst.get_restriction_map(genome = genome, enzyme = enzyme, output_dir = output)
    hst.get_dist_frags(genome = genome, restriction_map = restriction_map, circular = circular, rate = rate, output_dir = output)
    hst.get_dist_frags(genome = genome, restriction_map = restriction_map, circular = circular, rate = rate, output_dir = output)
    hst.log_bin_genome(genome = genome, output_dir = output)

    p1 = Process(target = hst.get_patterns(circular = circular, output_dir = output))
    p2 = Process(target = hst.generate_trans_ps(restriction_map = restriction_map, output_dir = output))
    p3 = Process(target = hst.generate_coverages(genome = genome, bins = bins, output_dir = output))
    p4 = Process(target = hst.generate_d1d2(output_dir = output))

    # Launch processes
    for process in [p1, p2, p3, p4]:
        process.start()
        process.join()

@click.command()
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
@click.option("--genome", "-g", required = True, default = None, type = str, help = "Genome to perform analysis on.")
@click.option("--enzyme", "-e", required = False, type = str, multiple = True, help = "Enzymes to use for genome digestion.")
@click.option("--nb-chunks", "-n", required = False, default = 1, type = int, help = "Number of chunks to split the alignment files into.")
@click.option("--mode", "-m", required = False, default = "full", type = str, help = "Statistical model to use for ambiguous reads assignment.")
@click.option("--cpus", "-t", required = False, default = 1, type = int, help = "Threads to use for analysis.")
def rescue_cmd(genome, enzyme, nb_chunks, mode, output, cpus):
    
    restriction_map = hst.get_restriction_map(genome = genome, enzyme = enzyme, output_dir = output)
    hut.chunk_bam(nb_chunks = nb_chunks, output_dir = output)
        
    # Get chunks as lists
    forward_chunks, reverse_chunks = hut.get_chunks(output)

    # Reattribute reads
    with multiprocessing.Pool(processes = cpus) as pool:

        results = pool.map_async(partial(hst.reattribute_reads, mode = mode,  restriction_map = restriction_map, output_dir = output),
        zip(forward_chunks, reverse_chunks))
        pool.close()
        pool.join()

    hio.merge_predictions(output_dir = output, clean = True)


# @click.command()
# def plot_cmd():
#     pass

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
# cli.add_command(plot_cmd, name="plot")