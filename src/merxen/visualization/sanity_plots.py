"""Sanity overlay plotting helpers."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_sanity_overlay(
    image: np.ndarray,
    output_path: Path | str,
    *,
    shapes: gpd.GeoDataFrame | None = None,
    points: pd.DataFrame | None = None,
    x_col: str = "x",
    y_col: str = "y",
    title: str = "Sanity Overlay",
    point_size: float = 1.5,
) -> Path:
    """Overlay shape boundaries and transcript points on an image crop."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6, 6))
    if image.ndim == 2:
        ax.imshow(image, cmap="gray")
    else:
        ax.imshow(image)

    if shapes is not None and len(shapes) > 0:
        shapes.boundary.plot(ax=ax, linewidth=0.5, color="#F97316", alpha=0.85)

    if points is not None and len(points) > 0:
        ax.scatter(
            points[x_col].to_numpy(float),
            points[y_col].to_numpy(float),
            s=point_size,
            c="#06B6D4",
            alpha=0.6,
            linewidths=0,
        )

    ax.set_title(title)
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(output_path, dpi=220)
    plt.close(fig)
    return output_path
