import click

@click.command()
@click.option('--input', '-i', required=True, help="Path to input reads.")
@click.option('--output', '-o', required=True, help="Output directory.")
def pipeline_cmd(input, output):
    """Run the full HiCberg pipeline."""
    click.echo(f"Starting pipeline for {input}. Results will be in {output}")
    # Later, this is where you will call:
    # engine = HicbergEngine(output_dir=output)
    # engine.execute_full_pipeline(input_data=input)