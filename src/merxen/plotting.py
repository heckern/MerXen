"""Shared Matplotlib output helpers."""

from __future__ import annotations

import os
from collections.abc import Sized
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection
from matplotlib.figure import Figure

BUNDLED_ARIAL_FONT_PATH = (
    Path(__file__).resolve().parent / "assets" / "fonts" / "arial.ttf"
)
ARIAL_FONT_ENV_VAR = "MERXEN_ARIAL_FONT_PATH"

_PLOT_OUTPUT_CONFIGURED = False


def configure_matplotlib_for_pdf(
    arial_font_path: Path | str | None = None,
) -> Path | None:
    """Configure Matplotlib so PDF outputs embed editable TrueType text."""
    global _PLOT_OUTPUT_CONFIGURED

    font_path = _resolve_arial_font_path(arial_font_path)
    if not _PLOT_OUTPUT_CONFIGURED:
        if font_path is not None:
            arial_prop = fm.FontProperties(fname=str(font_path))
            plt.rcParams["font.family"] = arial_prop.get_name()
            plt.rcParams.update({"mathtext.default": "regular"})
            fm.fontManager.addfont(str(font_path))
        matplotlib.rcParams["pdf.fonttype"] = 42
        _PLOT_OUTPUT_CONFIGURED = True
    return font_path


def prepare_plot_output(output_path: Path | str) -> Path:
    """Create the output directory and configure Matplotlib before plotting."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    configure_matplotlib_for_pdf()
    return path


def save_figure(
    fig: Figure,
    output_path: Path | str,
    *,
    dpi: int,
    bbox_inches: str | None = None,
    rasterize_min_points: int = 10,
    **savefig_kwargs: Any,
) -> Path:
    """Save a figure to PNG and to a sibling PDF path."""
    path = prepare_plot_output(output_path)
    rasterize_scatter_collections(fig, min_points=rasterize_min_points)

    save_args: dict[str, Any] = {"dpi": int(dpi), **savefig_kwargs}
    if bbox_inches is not None:
        save_args["bbox_inches"] = bbox_inches

    fig.savefig(path, **save_args)
    pdf_path = _pdf_output_path(path)
    if pdf_path != path:
        fig.savefig(pdf_path, **save_args)
    return path


def rasterize_scatter_collections(fig: Figure, *, min_points: int = 10) -> None:
    """Rasterize scatter point collections with more than a few points."""
    for ax in fig.axes:
        for collection in ax.collections:
            if not isinstance(collection, PathCollection):
                continue
            if _point_count(collection) >= min_points:
                collection.set_rasterized(True)


def _resolve_arial_font_path(arial_font_path: Path | str | None) -> Path | None:
    if arial_font_path is not None:
        path = Path(arial_font_path)
    elif os.environ.get(ARIAL_FONT_ENV_VAR):
        path = Path(os.environ[ARIAL_FONT_ENV_VAR])
    else:
        path = BUNDLED_ARIAL_FONT_PATH

    return path if path.exists() else None


def _pdf_output_path(output_path: Path) -> Path:
    if output_path.suffix.lower() == ".pdf":
        return output_path
    return output_path.with_suffix(".pdf")


def _point_count(collection: PathCollection) -> int:
    try:
        offsets = collection.get_offsets()
    except Exception:  # noqa: BLE001
        return 0
    if offsets is None or not isinstance(offsets, Sized):
        return 0
    return int(len(offsets))
