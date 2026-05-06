"""CLI command for Scanpy/Squidpy clustering analysis."""

from __future__ import annotations

from pathlib import Path

import click

from merxen.analysis.clustering_squidpy import run_clustering_squidpy
from merxen.config import ClusteringSquidpyConfig, load_config_from_json


@click.command(name="clustering-squidpy")
@click.option(
    "--config",
    "config_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Path to JSON config validated against ClusteringSquidpyConfig.",
)
def clustering_squidpy_command(config_path: Path) -> None:
    """Run Scanpy/Squidpy QC, clustering, UMAP, and spatial plots."""
    cfg = load_config_from_json(config_path, ClusteringSquidpyConfig)
    assert isinstance(cfg, ClusteringSquidpyConfig)

    results = run_clustering_squidpy(cfg)

    click.echo("clustering_squidpy complete:")
    for sample_id, paths in results.items():
        click.echo(f"- {sample_id}")
        for key, value in paths.items():
            click.echo(f"  - {key}: {value}")
