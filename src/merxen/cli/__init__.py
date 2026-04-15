"""Command-line entry points for MerXen."""

from __future__ import annotations

import click

from merxen.cli.run_comparison import compare_command
from merxen.cli.run_enrichment import enrich_command
from merxen.cli.run_qc import qc_command
from merxen.cli.run_segmentation import segment_command
from merxen.cli.run_visualization import visualize_command


@click.group(name="merxen")
def main() -> None:
    """MerXen spatial transcriptomics pipeline CLI."""


main.add_command(segment_command)
main.add_command(enrich_command)
main.add_command(qc_command)
main.add_command(compare_command)
main.add_command(visualize_command)
