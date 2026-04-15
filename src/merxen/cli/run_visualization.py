"""CLI command for generating visualization artifacts."""

from __future__ import annotations

from pathlib import Path

import click
import pandas as pd
import spatialdata as sd

from merxen.config import VisualizationConfig, load_config_from_json
from merxen.io.image_source import build_image_source, fetch_tile
from merxen.io.transcript_io import first_existing_col, to_pandas
from merxen.qc.gene_comparison import compute_gene_comparison_from_paths
from merxen.qc.metrics import compute_dataset_qc
from merxen.visualization.density_overview import plot_density_overview
from merxen.visualization.gene_scatter import plot_gene_scatter
from merxen.visualization.qc_plots import (
    plot_assignment_bar,
    plot_cell_metrics_violin,
    plot_geometry_histograms,
)
from merxen.visualization.sanity_plots import plot_sanity_overlay


@click.command(name="visualize")
@click.option(
    "--config",
    "config_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Path to JSON config validated against VisualizationConfig.",
)
def visualize_command(config_path: Path) -> None:
    """Generate visualization artifacts for a paired dataset comparison."""
    cfg = load_config_from_json(config_path, VisualizationConfig)
    assert isinstance(cfg, VisualizationConfig)
    cfg.output_dir.mkdir(parents=True, exist_ok=True)

    # Gene-level scatter plots from comparison tables.
    comparison = compute_gene_comparison_from_paths(
        xenium_zarr_path=cfg.xenium_zarr_path,
        merscope_zarr_path=cfg.merscope_zarr_path,
    )
    total_scatter = cfg.output_dir / f"{cfg.pair_id}_gene_scatter_total_normalized.png"
    assigned_scatter = (
        cfg.output_dir / f"{cfg.pair_id}_gene_scatter_assigned_normalized.png"
    )
    plot_gene_scatter(
        comparison["total_normalized_df"],
        total_scatter,
        title=f"{cfg.pair_id} Total Gene Counts (normalized)",
        x_label="Xenium (normalized)",
        y_label="MERSCOPE (normalized)",
        log_scale=True,
    )
    plot_gene_scatter(
        comparison["assigned_normalized_df"],
        assigned_scatter,
        title=f"{cfg.pair_id} Assigned Gene Counts (normalized)",
        x_label="Xenium (normalized)",
        y_label="MERSCOPE (normalized)",
        log_scale=True,
    )

    qc_records: list[dict[str, float | str]] = []
    dataset_items = [
        ("XENIUM", cfg.xenium_zarr_path),
        ("MERSCOPE", cfg.merscope_zarr_path),
    ]
    for dataset_name, zarr_path in dataset_items:
        qc = compute_dataset_qc(zarr_path, dataset_name=dataset_name)
        geom_plot = (
            cfg.output_dir / f"{cfg.pair_id}_{dataset_name.lower()}_geometry_hist.png"
        )
        cell_plot = (
            cfg.output_dir / f"{cfg.pair_id}_{dataset_name.lower()}_cell_violin.png"
        )
        plot_geometry_histograms(qc["geometry_metrics"], geom_plot)
        plot_cell_metrics_violin(qc["cell_metrics"], cell_plot)
        qc_records.append(
            {
                "dataset": dataset_name,
                "pct_assigned": float(qc["summary"]["pct_assigned"]),
            }
        )

        # Density + sanity overlays using the first points/image/shapes elements.
        sdata = sd.read_zarr(zarr_path)
        if len(sdata.points) > 0:
            points_key = list(sdata.points.keys())[0]
            points_pdf = to_pandas(sdata.points[points_key])
            x_col = first_existing_col(
                points_pdf,
                [
                    "x",
                    "x_micron",
                    "x_location",
                    "global_x",
                    "x_global_px",
                    "observed_x",
                ],
            )
            y_col = first_existing_col(
                points_pdf,
                [
                    "y",
                    "y_micron",
                    "y_location",
                    "global_y",
                    "y_global_px",
                    "observed_y",
                ],
            )
            if x_col is not None and y_col is not None:
                density_plot = (
                    cfg.output_dir
                    / f"{cfg.pair_id}_{dataset_name.lower()}_density_overview.png"
                )
                plot_density_overview(
                    points_pdf,
                    density_plot,
                    x_col=x_col,
                    y_col=y_col,
                    title=f"{dataset_name} Transcript Density",
                )

        if len(sdata.images) > 0:
            image_key = list(sdata.images.keys())[0]
            source = build_image_source(sdata.images[image_key], as_float32=False)
            height, width, _ = source["shape"]
            crop_size = int(min(1024, height, width))
            y0 = max(0, (height - crop_size) // 2)
            x0 = max(0, (width - crop_size) // 2)
            tile = fetch_tile(source, y0, y0 + crop_size, x0, x0 + crop_size)

            shapes = None
            if len(sdata.shapes) > 0:
                shapes_key = list(sdata.shapes.keys())[0]
                shapes = sdata.shapes[shapes_key]

            overlay_plot = (
                cfg.output_dir
                / f"{cfg.pair_id}_{dataset_name.lower()}_sanity_overlay.png"
            )
            plot_sanity_overlay(
                tile,
                overlay_plot,
                shapes=shapes,
                title=f"{dataset_name} Sanity Overlay",
            )

    assignment_df = pd.DataFrame(qc_records)
    assign_plot = cfg.output_dir / f"{cfg.pair_id}_assignment_rate_bar.png"
    plot_assignment_bar(assignment_df, assign_plot)

    click.echo("Visualization complete:")
    click.echo(f"- {total_scatter}")
    click.echo(f"- {assigned_scatter}")
    click.echo(f"- {assign_plot}")
