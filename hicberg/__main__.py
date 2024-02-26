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
import hicberg.plot as hpl
import hicberg.benchmark as hbk
import ast

from hicberg.version import __version__

# CONTEXT_SETTINGS = dict(help_option_names = ['-h', '--help'])

CONTEXT_SETTINGS = {"help_option_names" :  ['-h', '--help'],
                    "max_content_width" : 120}

epilogs = {"general" : "For more information about hicberg, please visit :https://github.com/sebgra/hicberg", 
        "bowtie2" : "For more information about Bowtie2, please visit : http://bowtie-bio.sourceforge.net/bowtie2/index.shtml",
        "samtools" : "For more information about Samtools, please visit : http://www.htslib.org/doc/samtools.html",
        "cooler" : "For more information about cooler, please visit : https://cooler.readthedocs.io/en/latest/",
        "complete" : "For more information about hicberg and used components, please visit :\nhttps://github.com/sebgra/hicberg\n - Bowtie2 : http://bowtie-bio.sourceforge.net/bowtie2/index.shtml\n - Samtools : http://www.htslib.org/doc/samtools.html"}

class PythonLiteralOption(click.Option):
    def type_cast_value(self, ctx, value):
        try:
            return ast.literal_eval(value)
        except:
            raise click.BadParameter(value)
        

@click.group(context_settings = CONTEXT_SETTINGS, epilog = epilogs["complete"], options_metavar = "<options>")
@click.version_option(version=__version__)
def cli(chain=True):
    """
    HiC-BERG: Rescuing of multi-reads in Hi-C maps.
    """
    ...


@click.command(context_settings = CONTEXT_SETTINGS, epilog = epilogs["complete"], options_metavar = "<options>", )
@click.argument('data', nargs = -1, metavar = "<genome> <input1> <input2>")
@click.option("--index", "-i", required = False, default = None, type = str, show_default = True, metavar = "<str>", help = "Index of the genome.")
@click.option("--name", "-n", required = False, default = "sample", type = str, show_default = True, metavar = "<str>",  help = "Name of the analysis.")
@click.option("--rate", "-r", required = False, default = 1.0, type = float, show_default = True, metavar = "<float>", help = "Rate to use for sub-sampling restriction map.")
@click.option("--cpus", "-t", required = False, default = 1, type = int, show_default = True, metavar = "<int>", help = "Threads to use for analysis.")
@click.option("--kernel-size", "-K", required = False, default = 11, type = int, show_default = True, metavar = "<int>", help = "Size of the gaussian kernel for contact density estimation.")
@click.option("--deviation", "-d", required = False, default = 0.5, type = float, show_default = True, metavar = "<float>", help = "Standard deviation for contact density estimation.")
@click.option("--mode", "-m", required = False, default = "full", type = str, show_default = True, metavar = "<str>", help = "Statistical model to use for ambiguous reads assignment.")
@click.option("--max-alignment", '-k', required = False, type = int, default = None, metavar = "<int>", help = "Set the number of alignments to report in ambiguous reads case. Value of -1 reports all alignments.")
@click.option("--sensitivity", "-s", required = False, type = click.Choice(["very-sensitive", "sensitive", "fast", "very-fast"]), default = "very-sensitive", show_default = True, metavar = "<str>", help = "Set sensitivity level for Bowtie2.")
@click.option("--bins", "-b", required = False, type = int, default = 2000, show_default = True, metavar = "<int>", help = "Genomic resolution.")
@click.option("--enzyme", "-e", required = False, type = str, multiple = True, show_default = True, metavar = "<str>", help = "Enzymes to use for genome digestion.")
@click.option("--circular", "-c", required = False, type = str, default = "", show_default = True, metavar = "<str>", help = "Name of the chromosome to consider as circular.")
@click.option("--mapq", "-q", required = False, type = int, default = 35, show_default = True, metavar = "<int>", help = "Minimum mapping quality to consider a read as valid.")
@click.option("--output", "-o", required = False, default = None, type = str, metavar = "<str>", help = "Output folder to save results.")
@click.option("--start-stage", required = False, type = click.Choice(["fastq", "bam", "groups", "build", "stats", "rescue", "final"]), default = "fastq", show_default = True, metavar = "<str>", help = "Stage to start the pipeline.")
@click.option("--exit-stage", required = False, type = click.Choice(["None", "bam", "groups", "build", "stats", "rescue", "final"]), default = "None", show_default = True, metavar = "<str>", help = "Stage to exit the pipeline.")
@click.option("--force", "-f", is_flag = True, show_default = True, help = "Set if previous analysis files are deleted.")
def pipeline_cmd(data, index, name, rate, mode, kernel_size, deviation, cpus, output, max_alignment, sensitivity, bins, enzyme, circular, mapq, start_stage, exit_stage, force):
    """
    Hi-C pipeline to generate enhanced contact matrix from fastq files.
    """
    hpp.pipeline(genome = data[0], index = index, name = name, fq_for = data[1], fq_rev = data[2], output_dir = output, cpus = cpus, rate = rate, nb_chunks = 2 * cpus, mode = mode, kernel_size = kernel_size, deviation = deviation, max_alignment = max_alignment,  sensitivity = sensitivity, bins = bins, enzyme = enzyme, circular = circular, mapq = mapq, start_stage = start_stage, exit_stage = exit_stage, force = force)
    return


@click.command(context_settings = CONTEXT_SETTINGS, epilog = epilogs["general"], options_metavar = "<options>")
@click.option("--output", "-o", required = False, default = None, type = str, show_default = True, metavar = "<str>", help = "Output folder to save results. If not set, the current directory is used.")
@click.option("--name", "-n", required = False, default = None, type = str, show_default = True, metavar = "<str>",  help = "Name of the output folder to create. If not set, 'sample' is used.")
@click.option("--force", "-f", is_flag = True, help = "Set if previous analysis files are deleted.")
def create_folder_cmd(output, force, name):
    """
    Create a folder to save results. Folder will be set as <output>/<name>.
    """
    hio.create_folder(sample_name=name, output_dir=output, force=force)
    return

@click.command(context_settings = CONTEXT_SETTINGS, epilog = epilogs["general"], options_metavar = "<options>")
@click.argument('data', nargs = -1, metavar = "<genome>")
@click.option("--output", "-o", required = False, default = None, type = str, show_default = True, metavar = "<str>", help = "Output folder to save results. If not set, the current directory is used.")
@click.option("--bins", "-b", required = False, type = int, default = 2000, show_default = True, metavar = "<int>", help = "Genomic resolution.")
def get_tables_cmd(data, bins, output):
    """
    Create tables for the genome length detail and the bins.
    """
    hut.get_chromosomes_sizes(genome = data[0], output_dir = output)
    hut.get_bin_table(bins = bins, output_dir = output)
    return


@click.command(context_settings = CONTEXT_SETTINGS, epilog = epilogs["bowtie2"] + "\n"  + epilogs["samtools"], options_metavar = "<options>")
@click.argument('data', nargs = -1, metavar = "<genome> <input1> <input2>")
@click.option("--index", "-i", required = False, default = None, type = str, show_default = True, metavar = "<str>", help = "Index of the genome (path). If not set, the index is built.")
@click.option("--output", "-o", required = False, default = None, type = str, show_default = True, metavar = "<str>", help = "Output folder to save results. If not set, the current directory is used.")
@click.option("--sensitivity", "-s", required = False, type = click.Choice(["very-sensitive", "sensitive", "fast", "very-fast"]), default = "very-sensitive", show_default = True, metavar = "<str>", help = "Set sensitivity level for Bowtie2.")
@click.option("--max-alignment", '-k', required = False, type = int, default = None, show_default = True, metavar = "<int>", help = "Set the number of alignments to report in ambiguous reads case. If set to -1, all alignments are reported.")
@click.option("--cpus", "-t", required = False, default = 1, type = int, show_default = True, metavar = "<int>", help = "Threads to use for analysis.")
@click.option("--verbose", "-v", is_flag = True, help = "Set verbosity level.")
def alignment_cmd(data, index, max_alignment, sensitivity, output, cpus, verbose):
    """
    Perform alignment of Hi-C reads.
    """
    if index is None:
        index = hal.hic_build_index(genome = data[0], output_dir = output, cpus = cpus, verbose = verbose)

    hal.hic_align(index = index, fq_for = data[1], fq_rev = data[2], sensitivity = sensitivity, max_alignment = max_alignment, output_dir = output, cpus = cpus, verbose = True)
    hal.hic_view(output_dir = output, cpus = cpus, verbose = verbose)
    hal.hic_sort(output_dir = output, cpus = cpus, verbose = verbose)
    return


@click.command(context_settings = CONTEXT_SETTINGS, epilog = epilogs["general"], options_metavar = "<options>")
@click.option("--mapq", "-q", required = False, type = int, default = 35, help = "Minimum mapping quality to consider a read as non ambiguous.")
@click.option("--output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
def classify_cmd(mapq, output):
    """
    Perform classification of Hi-C reads (pairs).
    3 groups wil be defined and 2 alignment files (.bam) will be created per group:\n
    - Unmapped read pairs (group 0)\n
    - Read pairs with both reads mapping at only one position (group 1)\n
    - Read pairs with at least one read mapping at multiple positions (group 2)
    """
    hut.classify_reads(mapq = mapq, output_dir = output)
    return


@click.command(context_settings = CONTEXT_SETTINGS, epilog = epilogs["general"], options_metavar = "<options>")
@click.option("--output", "-o", required = False, default = None, type = str, show_default = True, metavar = "<str>", help = "Output folder to save results.")
@click.option("--recover", "-r", required = False, default = False, is_flag = True, help = "Set if pairs are built after reads reassignment. Therefore alignment files of group2 will be used.")
def build_pairs_cmd(output, recover):
    """
    Create pair files from a pair of alignment files.
    """
    hio.build_pairs(output_dir = output, mode = recover)
    return


@click.command(context_settings = CONTEXT_SETTINGS, epilog = epilogs["cooler"], options_metavar = "<options>")
@click.option("--output", "-o", required = False, default = None, type = str, show_default = True, metavar = "<str>", help = "Output folder to save results.")
@click.option("--recover", "-r", required = False, default = False, is_flag = True, help = "Set if .cool matrix are built after reads reassignment. Therefor pairs file from group2 will be used.")
@click.option("--cpus", "-t", required = False, default = 1, type = int, help = "Threads to use for matrix building.")
def build_matrix_cmd(output, recover, cpus):
    """
    Create matrix (.cool) from pairs files.
    """
    hio.build_matrix(cpus = cpus, output_dir = output, mode = recover)

    return

@click.command(context_settings = CONTEXT_SETTINGS, epilog = epilogs["general"], options_metavar = "<options>")
@click.argument('data', nargs = -1, metavar = "<genome>")
@click.option("--output", "-o", required = False, default = None, type = str, show_default = True, metavar = "<str>", help = "Output folder to save results.")
@click.option("--mode", "-m", required = False, default = "full", type = str, show_default = True, metavar = "<str>", help = "Statistical model to use for ambiguous reads assignment.")
@click.option("--kernel-size", "-K", required = False, default = 11, type = int, show_default = True, metavar = "<int>", help = "Size of the gaussian kernel for contact density estimation.")
@click.option("--deviation", "-d", required = False, default = 0.5, type = float, show_default = True, metavar = "<float>", help = "Standard deviation for contact density estimation.")
@click.option("--rate", "-r", required = False, default = 1.0, type = float, show_default = True, metavar = "<float>", help = "Rate to use for sub-sampling restriction map.")
@click.option("--enzyme", "-e", required = False, type = str, multiple = True, show_default = True, metavar = "<str>", help = "Enzymes to use for genome digestion.")
@click.option("--circular", "-c", required = False, type = str, default = "", show_default = True, metavar = "<str>", help = "Name of the chromosome to consider as circular.")
@click.option("--bins", "-b", required = False, type = int, default = 2000, show_default = True, metavar = "<int>", help = "Genomic resolution.")
@click.option("--cpus", "-t", required = False, default = 1, type = int, show_default = True, metavar = "<int>", help = "Threads to use for analysis.")
def statistics_cmd(data, mode, kernel_size, deviation,  rate, enzyme, circular, bins, output, cpus):
    """
    Extract statistics from non ambiguous Hi-C data.
    """
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
    return
        

@click.command(context_settings = CONTEXT_SETTINGS, epilog = epilogs["general"], options_metavar = "<options>")
@click.argument('data', nargs = -1, metavar = "<genome>")
@click.option("--output", "-o", required = False, default = None, type = str, show_default = True, metavar = "<str>", help = "Output folder to save results.")
@click.option("--enzyme", "-e", required = False, type = str, multiple = True, show_default = True, metavar = "<str>", help = "Enzymes to use for genome digestion.")
@click.option("--mode", "-m", required = False, default = "full", type = str, show_default = True, metavar = "<str>", help = "Statistical model to use for ambiguous reads assignment.")
@click.option("--cpus", "-t", required = False, default = 1, type = int, show_default = True, metavar = "<int>", help = "Threads to use for analysis.")
def rescue_cmd(data, enzyme, mode, output, cpus):
    """
    Reallocate ambiguous reads to the most plausible position according to model.
    """
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
    return 


@click.command(context_settings = CONTEXT_SETTINGS, epilog = epilogs["general"], options_metavar = "<options>")
@click.argument('data', nargs = -1, metavar = "<genome>")
@click.option("--output", "-o", required = False, default = None, type = str, show_default = True, metavar = "<str>", help = "Output folder to save results.")
@click.option("--bins", "-b", required = False, default = 2000, type = int, show_default = True, metavar = "<int>", help = "Genomic resolution.")
def plot_cmd(data, bins, output):
    """
    Plot results from analysis.
    """
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

    return

@click.command(context_settings = CONTEXT_SETTINGS, epilog = epilogs["general"], options_metavar = "<options>")
@click.option("--output", "-o", required = False, default = None, type = str, show_default = True, metavar = "<str>", help = "Output folder to save results.")
def tidy_cmd(output):
    """
    Tidy output folder.
    """
    hio.tidy_folder(output_dir = output)
    return


@click.command(context_settings = CONTEXT_SETTINGS, epilog = epilogs["general"], options_metavar = "<options>")
@click.argument('data', nargs = -1, metavar = "<genome>")
@click.option("--output", "-o", required = False, default = None, type = str, show_default = True, metavar = "<str>", help = "Output folder to save results.")
@click.option("--chromosome", "-c", required = False, default = None, type = str, show_default = True, metavar = "<str>", help = "Chromosome to get as source for duplication.")
@click.option("--position", "-p", required = False, default = None, type = int, show_default = True, metavar = "<int>", help = "Position to get as source for duplication.")
@click.option("--trans-chromosome", "-C", required = False, default = None, type = str, show_default = True, metavar = "<str>", help = "Chromosome to get as target for duplication.")
@click.option("--trans-position", "-P", required = False, default = None, type = str, show_default = True, metavar = "<int>", help = "Position to get as target for duplication.")
@click.option("--bins", "-b", required = False, default = 1, type = int, show_default = True, metavar = "<int>", help = "Number of bins to select from a genomic coordinates.")
@click.option("--strides", "-s", required = False, default = None, type = str, show_default = True, metavar = "<int>", help = "Strides to apply from source genomic coordinates to define targets intervals. Multiple strides must be coma separated.")
@click.option("--auto", "-a", required = False, default = None, type = int, show_default = True, metavar = "<int>", help = "Automatically select auto intervals for duplication.")
@click.option("--kernel-size", "-K", required = False, default = 11, type = int, show_default = True, metavar = "<int>", help = "Size of the gaussian kernel for contact density estimation.")
@click.option("--deviation", "-d", required = False, default = 0.5, type = float, show_default = True, metavar = "<float>", help = "Standard deviation for contact density estimation.")
@click.option("--mode", "-m", required = False, default = "full", type = str, show_default = True, metavar = "<str>", help = "Statistical model to use for ambiguous reads assignment. Multiple modes must be coma separated.")
@click.option("--pattern", "-S", required = False, type = click.Choice(["loops", "borders", "hairpins", "-1"]), default = None, show_default = True, metavar = "<str>", help = "Set pattern if benchmarking considering patterns.")
@click.option("--threshold", "-t", required = False, type = float, default = 0.0, show_default = True, metavar = "<float>", help = "Set pattern score threshold under which pattern are discarded.")
@click.option("--jitter", "-j", required = False, type = int, default = 0, show_default = True, metavar = "<int>", help = "Set jitter for pattern detection interval overlapping")
@click.option("--trend", "-T", is_flag = False, help = "Set if detrending of the contact map has to be performed.")
@click.option("--top", "-k", required = False, type = int, default = 100, show_default = True, metavar = "<int>", help = "Set the top k % of patterns to retain.")
@click.option("--iterations", "-i", required = False, type = int, default = 3, show_default = True, metavar = "<int>", help = "Set the number of iterations for benchmarking.")
@click.option("--force", "-f", is_flag = True, help = "Set if previous analysis files have to be deleted.")
@click.option("--cpus", required = False, default = 1, type = int, show_default = True, metavar = "<int>", help = "Threads to use for analysis.")
def benchmark_cmd(data, chromosome, position, trans_chromosome, trans_position, bins, strides, auto, kernel_size, deviation, mode, pattern, threshold, jitter, trend, top, iterations, force, output, cpus):
    """
    Perform benchmarking of the statistical model (this can be time consuming).
    """
    hbk.benchmark(output_dir = output, genome = data[0], chromosome = chromosome, position = position, trans_chromosome = trans_chromosome, trans_position = trans_position, strides = strides, mode = mode, force = force, bins = bins, auto = auto, kernel_size = kernel_size, deviation = deviation, pattern = pattern, threshold = threshold, jitter = jitter, trend = trend, top = top, iterations = iterations, cpus = cpus)
    return

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
