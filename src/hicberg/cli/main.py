import click
from hicberg.version import __version__
# Use relative or absolute imports for your other modules
from hicberg.cli.commands import pipeline_cmd # etc.

@click.group()
@click.version_option(__version__)
def cli():
    """HiCberg: Rescuing multi-reads in Hi-C maps."""
    pass

# Add your commands
cli.add_command(pipeline_cmd, name="pipeline")