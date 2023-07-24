from email.policy import default
import click
import hicberg.io as hio
import hicberg.utils as hut
import hicberg.align as hal
import hicberg.statistics as hstst
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
@click.option("rate", "-r", required = False, default = 1.0, type = float, help = "Rate to use for subsampling restriction map.")
@click.option("cpus", "-t", required = False, default = 1, type = int, help = "Threads to use for analysis.")
@click.option("mode", "-m", required = False, default = "full", type = str, help = "Statistical model to use for ambiguous reads assignment.")
@click.option("output", "-o", required = False, default = None, type = str, help = "Output folder to save results.")
@click.option("--force", "-f", is_flag = True, help = "Set if previous analysis files are deleted")
def pipeline_cmd(genome, fq_for, fq_rev, rate, mode, cpus, output, force):

    hpp.pipeline(genome = genome, fq_for = fq_for, fq_rev = fq_rev, output_dir = output, cpus = cpus, rate = rate, nb_chunks = cpus, mode = mode, force = force)
    


cli.add_command(pipeline_cmd, name="pipeline")