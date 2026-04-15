"""Transcript density overview plotting."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def density_hist2d(
    points_df: pd.DataFrame,
    *,
    x_col: str = "x",
    y_col: str = "y",
    bins: int = 400,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute 2D histogram density from point coordinates."""
    x = pd.to_numeric(points_df[x_col], errors="coerce").to_numpy(np.float64)
    y = pd.to_numeric(points_df[y_col], errors="coerce").to_numpy(np.float64)
    valid = np.isfinite(x) & np.isfinite(y)
    if not np.any(valid):
        raise ValueError("No valid point coordinates available for density histogram.")
    hist, x_edges, y_edges = np.histogram2d(x[valid], y[valid], bins=bins)
    return hist, x_edges, y_edges


def plot_density_overview(
    points_df: pd.DataFrame,
    output_path: Path | str,
    *,
    x_col: str = "x",
    y_col: str = "y",
    bins: int = 400,
    title: str = "Transcript Density",
) -> Path:
    """Render a full-field transcript density heatmap from point coordinates."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    hist, x_edges, y_edges = density_hist2d(
        points_df, x_col=x_col, y_col=y_col, bins=bins
    )

    fig, ax = plt.subplots(figsize=(7, 6))
    extent = (
        float(x_edges[0]),
        float(x_edges[-1]),
        float(y_edges[0]),
        float(y_edges[-1]),
    )
    im = ax.imshow(
        hist.T,
        origin="lower",
        extent=extent,
        aspect="auto",
        cmap="magma",
        interpolation="nearest",
    )
    fig.colorbar(im, ax=ax, label="Transcript count")
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=220)
    plt.close(fig)
    return output_path
