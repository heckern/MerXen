"""Sanity overlay plotting helpers."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _prepare_overlay_image(image: np.ndarray) -> tuple[np.ndarray, str | None]:
    """Convert microscopy image data into a display-safe grayscale or RGB array."""
    arr = np.asarray(image)

    if arr.ndim == 2:
        return _normalize_channel(arr), "gray"

    if arr.ndim != 3:
        raise ValueError(f"Expected 2D or 3D image data, got shape {arr.shape}")

    if arr.shape[-1] == 1:
        return _normalize_channel(arr[..., 0]), "gray"

    arr = arr.astype(np.float32, copy=False)
    if arr.shape[-1] == 2:
        arr = np.concatenate([arr, np.zeros_like(arr[..., :1])], axis=-1)
    elif arr.shape[-1] > 3:
        arr = arr[..., :3]

    return _normalize_channel(arr), None


def _normalize_channel(arr: np.ndarray) -> np.ndarray:
    """Percentile-normalize an array to the [0, 1] range for plotting."""
    arr = arr.astype(np.float32, copy=False)
    finite = np.isfinite(arr)
    if not np.any(finite):
        return np.zeros_like(arr, dtype=np.float32)

    values = arr[finite]
    lo, hi = np.percentile(values, (2, 98))
    if hi <= lo:
        hi = lo + 1.0
    scaled = np.asarray(np.clip((arr - lo) / (hi - lo), 0.0, 1.0), dtype=np.float32)
    scaled[~finite] = 0.0
    return scaled


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
    display_image, cmap = _prepare_overlay_image(image)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(display_image, cmap=cmap)

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
