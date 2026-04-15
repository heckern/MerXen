"""QC plotting utilities."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.patches import Rectangle


def plot_geometry_histograms(
    geometry_metrics: pd.DataFrame,
    output_path: Path | str,
    *,
    bins: int = 50,
) -> Path:
    """Plot geometry metric histograms for area/eccentricity/aspect ratio."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cols = [
        c for c in ["area", "eccentricity", "aspect_ratio"] if c in geometry_metrics
    ]
    ncols = max(1, len(cols))
    fig, axes = plt.subplots(1, ncols, figsize=(5 * ncols, 4))
    if ncols == 1:
        axes = [axes]

    for ax, col in zip(axes, cols, strict=False):
        vals = pd.to_numeric(geometry_metrics[col], errors="coerce").dropna()
        ax.hist(vals, bins=bins, color="#1f77b4", alpha=0.8)
        ax.set_title(col.replace("_", " ").title())
        ax.set_xlabel(col)
        ax.set_ylabel("Count")

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def plot_cell_metrics_violin(
    cell_metrics: pd.DataFrame,
    output_path: Path | str,
) -> Path:
    """Plot violin distributions of transcripts/cell and genes/cell."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cols = [c for c in ["transcripts_per_cell", "genes_per_cell"] if c in cell_metrics]
    if not cols:
        raise ValueError(
            "cell_metrics is missing required columns for violin plotting."
        )

    melted = cell_metrics[cols].melt(var_name="metric", value_name="value")
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.violinplot(data=melted, x="metric", y="value", fill=True, ax=ax)
    ax.set_xlabel("Metric")
    ax.set_ylabel("Value")
    ax.set_title("Per-cell QC Metrics")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def plot_assignment_bar(
    assignment_summaries: pd.DataFrame,
    output_path: Path | str,
    *,
    dataset_col: str = "dataset",
    pct_col: str = "pct_assigned",
) -> Path:
    """Plot per-dataset assignment rate bar chart."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if dataset_col not in assignment_summaries or pct_col not in assignment_summaries:
        raise KeyError(
            f"Expected columns '{dataset_col}' and '{pct_col}' in assignment_summaries."
        )

    fig, ax = plt.subplots(figsize=(6, 4))
    sns.barplot(data=assignment_summaries, x=dataset_col, y=pct_col, ax=ax)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Assigned Transcripts (%)")
    ax.set_xlabel("Dataset")
    ax.set_title("Transcript Assignment Rate")
    for patch in ax.patches:
        if not isinstance(patch, Rectangle):
            continue
        height = patch.get_height()
        ax.annotate(
            f"{height:.1f}%",
            (patch.get_x() + patch.get_width() / 2.0, height),
            ha="center",
            va="bottom",
            fontsize=9,
        )
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path
