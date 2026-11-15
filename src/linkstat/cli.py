"""Command-line entry point. Subcommands are added as their modules are implemented (P1 steps 1.9-1.16)."""

import click

from . import __version__


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, prog_name="linkstat")
def main():
    """Chance-corrected co-localization of gene classes on long reads and contigs."""


if __name__ == "__main__":
    main()
