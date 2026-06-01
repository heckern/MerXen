"""Tests for gene scatter plotting."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from merxen.config import VisualizationConfig
from merxen.visualization.gene_scatter import plot_gene_abundance, plot_gene_scatter


def test_plot_gene_scatter_writes_file(tmp_path: Path) -> None:
    """Scatter plotting should write PNG and PDF outputs."""
    df = pd.DataFrame(
        {
            "gene": ["A", "B", "C"],
            "xenium": [1.0, 2.0, 4.0],
            "merscope": [1.2, 1.9, 4.1],
        }
    )
    out = tmp_path / "scatter.png"
    plot_gene_scatter(df, out, title="Test Plot")
    assert out.exists()
    assert out.with_suffix(".pdf").exists()


def test_plot_gene_abundance_writes_file(tmp_path: Path) -> None:
    """Single-dataset gene abundance plotting should write PNG and PDF outputs."""
    df = pd.DataFrame(
        {
            "gene": ["A", "B", "C"],
            "count": [10.0, 5.0, 1.0],
            "normalized": [0.625, 0.3125, 0.0625],
        }
    )
    out = tmp_path / "abundance.png"

    plot_gene_abundance(df, out, title="Single Dataset")

    assert out.exists()
    assert out.with_suffix(".pdf").exists()


def test_visualization_config_accepts_single_sample(tmp_path: Path) -> None:
    """Visualization config should support one-platform sample lists."""
    zarr_path = tmp_path / "sample.zarr"
    cfg = VisualizationConfig.model_validate(
        {
            "pair_id": "P1",
            "output_dir": str(tmp_path / "out"),
            "samples": [
                {
                    "sample_id": "P1_XENIUM",
                    "platform": "XENIUM",
                    "zarr_path": str(zarr_path),
                }
            ],
        }
    )

    assert len(cfg.samples) == 1
    assert cfg.samples[0].platform == "XENIUM"
