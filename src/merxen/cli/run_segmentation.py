"""CLI command for segmentation (Cellpose + ProSeg)."""

from __future__ import annotations

import logging
from pathlib import Path

import click

from merxen.config import SegmentationConfig, load_config_from_json
from merxen.segmentation.pipeline import run_segmentation_pipeline

logger = logging.getLogger(__name__)


@click.command(name="segment")
@click.option(
    "--config",
    "config_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Path to JSON config validated against SegmentationConfig.",
)
@click.option(
    "--force-rerun",
    is_flag=True,
    default=False,
    help="Recompute outputs even when existing zarr outputs are present.",
)
def segment_command(config_path: Path, force_rerun: bool) -> None:
    """Run unified segmentation for one dataset."""
    cfg = load_config_from_json(config_path, SegmentationConfig)
    assert isinstance(cfg, SegmentationConfig)
    outputs = run_segmentation_pipeline(cfg, force_rerun=force_rerun)
    click.echo("Segmentation complete:")
    for key, value in outputs.items():
        click.echo(f"- {key}: {value}")
