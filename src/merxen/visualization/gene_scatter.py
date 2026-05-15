"""Per-gene scatter plotting with linear fits."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from merxen.plotting import prepare_plot_output, save_figure
from merxen.qc.gene_comparison import fit_linear


def plot_gene_scatter(
    df: pd.DataFrame,
    output_path: Path | str,
    *,
    title: str,
    x_label: str = "Xenium",
    y_label: str = "MERSCOPE",
    log_scale: bool = True,
) -> Path:
    """Plot gene-level scatter and fitted line for Xenium vs MERSCOPE counts."""
    output_path = prepare_plot_output(output_path)

    fig, ax = plt.subplots(figsize=(5, 5))
    if df.empty:
        ax.set_title(title + " (no overlapping genes)")
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        fig.tight_layout()
        save_figure(fig, output_path, dpi=220)
        plt.close(fig)
        return output_path

    x = df["xenium"].to_numpy(dtype=float)
    y = df["merscope"].to_numpy(dtype=float)

    if log_scale:
        eps = 1e-12
        x_plot = np.clip(x, eps, None)
        y_plot = np.clip(y, eps, None)
        ax.scatter(
            x_plot,
            y_plot,
            s=22,
            alpha=0.75,
            edgecolor="none",
            rasterized=True,
        )
        ax.set_xscale("log")
        ax.set_yscale("log")

        lo = float(min(x_plot.min(), y_plot.min()))
        hi = float(max(x_plot.max(), y_plot.max()))
        ax.plot([lo, hi], [lo, hi], "--", linewidth=1.2, label="y = x")

        lx = np.log10(x_plot)
        ly = np.log10(y_plot)
        slope, intercept, r2 = fit_linear(lx, ly)
        if np.isfinite(slope):
            x_fit_log = np.logspace(np.log10(lo), np.log10(hi), 200)
            y_fit_log = 10 ** (slope * np.log10(x_fit_log) + intercept)
            ax.plot(x_fit_log, y_fit_log, linewidth=1.5, label="best fit")
            ax.text(
                0.03,
                0.97,
                f"log10(y) = {slope:.3f} * log10(x) + {intercept:.3f}\nR² = {r2:.3f}",
                transform=ax.transAxes,
                va="top",
                ha="left",
                fontsize=9,
                bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.8},
            )
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
    else:
        ax.scatter(x, y, s=22, alpha=0.75, edgecolor="none", rasterized=True)
        lo = 0.0
        hi = float(max(x.max(), y.max()))
        ax.plot([lo, hi], [lo, hi], "--", linewidth=1.2, label="y = x")
        slope, intercept, r2 = fit_linear(x, y)
        if np.isfinite(slope):
            x_fit_linear = np.linspace(lo, hi, 200)
            y_fit_linear = slope * x_fit_linear + intercept
            ax.plot(x_fit_linear, y_fit_linear, linewidth=1.5, label="best fit")
            ax.text(
                0.03,
                0.97,
                f"y = {slope:.3f} * x + {intercept:.3e}\nR² = {r2:.3f}",
                transform=ax.transAxes,
                va="top",
                ha="left",
                fontsize=9,
                bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.8},
            )
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)

    ax.set_aspect("equal")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend(loc="lower right")
    fig.tight_layout()
    save_figure(fig, output_path, dpi=220)
    plt.close(fig)
    return output_path
