# Import the main click group function from the new dedicated module
from .cli import cli

if __name__ == "__main__":
    # Execute the click command group
    cli()