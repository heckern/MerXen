"""Cross-platform gene-level comparison utilities."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import spatialdata as sd

from merxen.io.transcript_io import assignment_mask, first_existing_col, to_pandas

logger = logging.getLogger(__name__)


def gene_totals_from_table(adata: Any, gene_label: str = "gene") -> pd.Series:
    """Aggregate per-gene totals from an AnnData table."""
    if gene_label in adata.var.columns:
        genes = adata.var[gene_label].astype(str)
    else:
        genes = adata.var_names.astype(str)

    totals = np.asarray(adata.X.sum(axis=0)).ravel()
    return pd.Series(totals, index=genes).groupby(level=0).sum().astype(float)


def gene_totals_from_points(
    sdata_obj: Any,
    assigned_only: bool = False,
) -> pd.Series:
    """Aggregate per-gene transcript counts from the first points table."""
    if len(sdata_obj.points) == 0:
        raise RuntimeError("No points found in spatialdata object.")

    points_key = list(sdata_obj.points.keys())[0]
    pts = sdata_obj.points[points_key]

    gene_col = first_existing_col(pts, ["gene", "feature_name", "target"])
    if gene_col is None:
        raise KeyError(
            f"Could not find gene column in points {points_key}: {list(pts.columns)}"
        )

    assign_col: str | None = None
    if assigned_only:
        assign_col = first_existing_col(pts, ["assignment", "cell", "cell_id"])
        if assign_col is None:
            raise KeyError(
                "Could not find assignment column in points "
                f"{points_key}: {list(pts.columns)}"
            )

    if hasattr(pts, "npartitions") and hasattr(pts, "partitions"):
        if assigned_only:
            work = pts[[gene_col, assign_col]]
            mask = work[assign_col].map_partitions(
                assignment_mask, meta=("assigned", "bool")
            )
            counts = work.loc[mask].groupby(gene_col).size().compute()
        else:
            counts = pts[[gene_col]].groupby(gene_col).size().compute()
    else:
        pdf = to_pandas(pts)
        if assigned_only:
            mask = assignment_mask(pdf[assign_col])
            counts = pdf.loc[mask].groupby(gene_col).size()
        else:
            counts = pdf.groupby(gene_col).size()

    counts.index = counts.index.astype(str)
    return counts.groupby(level=0).sum().astype(float)


def apply_dataset_filter(gene_counts: pd.Series, dataset_name: str) -> pd.Series:
    """Remove platform-specific control probes from gene counts."""
    idx = gene_counts.index.astype(str)
    if dataset_name.upper() == "XENIUM":
        keep = ~idx.str.contains("Blank", na=False)
    elif dataset_name.upper() == "MERSCOPE":
        keep = ~idx.str.contains(
            "UnassignedCodeword|NegControlCodeword|NegControlProbe",
            regex=True,
            na=False,
        )
    else:
        keep = pd.Series(True, index=gene_counts.index)
    return gene_counts.loc[keep].copy()


def normalize_counts(gene_counts: pd.Series) -> tuple[pd.Series, float]:
    """Normalize counts by total transcripts and return (normalized, total)."""
    total = float(gene_counts.sum())
    if total <= 0:
        raise ValueError("Cannot normalize: total count is <= 0.")
    return gene_counts / total, total


def compare_df(x_counts: pd.Series, m_counts: pd.Series) -> pd.DataFrame:
    """Construct an aligned comparison DataFrame across shared genes."""
    common = sorted(set(x_counts.index.astype(str)) & set(m_counts.index.astype(str)))
    return pd.DataFrame(
        {
            "gene": common,
            "xenium": x_counts.reindex(common).values,
            "merscope": m_counts.reindex(common).values,
        }
    )


def fit_linear(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Fit a 1D linear model and return slope, intercept, and R²."""
    if len(x) < 2:
        return np.nan, np.nan, np.nan
    slope, intercept = np.polyfit(x, y, 1)
    yhat = slope * x + intercept
    denom = np.sum((y - y.mean()) ** 2)
    r2 = np.nan if denom <= 0 else (1.0 - np.sum((y - yhat) ** 2) / denom)
    return float(slope), float(intercept), float(r2)


def compute_gene_comparison(
    xenium_sdata: Any,
    merscope_sdata: Any,
) -> dict[str, Any]:
    """Compute cross-platform gene totals, normalized counts, and fit metrics."""
    logger.info("[Gene Compare] Computing per-gene totals from points and tables")

    xenium_table = xenium_sdata.tables["table"]
    merscope_table = merscope_sdata.tables["table"]

    x_total_all = gene_totals_from_points(xenium_sdata, assigned_only=False)
    m_total_all = gene_totals_from_points(merscope_sdata, assigned_only=False)
    x_assigned_all = gene_totals_from_table(xenium_table)
    m_assigned_all = gene_totals_from_table(merscope_table)

    x_total = apply_dataset_filter(x_total_all, "XENIUM")
    m_total = apply_dataset_filter(m_total_all, "MERSCOPE")
    x_assigned = apply_dataset_filter(x_assigned_all, "XENIUM")
    m_assigned = apply_dataset_filter(m_assigned_all, "MERSCOPE")

    x_total_norm, x_total_sum = normalize_counts(x_total)
    m_total_norm, m_total_sum = normalize_counts(m_total)
    x_assigned_norm, x_assigned_sum = normalize_counts(x_assigned)
    m_assigned_norm, m_assigned_sum = normalize_counts(m_assigned)

    total_df = compare_df(x_total, m_total)
    assigned_df = compare_df(x_assigned, m_assigned)
    total_norm_df = compare_df(x_total_norm, m_total_norm)
    assigned_norm_df = compare_df(x_assigned_norm, m_assigned_norm)

    lx = np.log10(np.clip(total_norm_df["xenium"].to_numpy(float), 1e-12, None))
    ly = np.log10(np.clip(total_norm_df["merscope"].to_numpy(float), 1e-12, None))
    total_fit = fit_linear(lx, ly)

    lax = np.log10(np.clip(assigned_norm_df["xenium"].to_numpy(float), 1e-12, None))
    lay = np.log10(np.clip(assigned_norm_df["merscope"].to_numpy(float), 1e-12, None))
    assigned_fit = fit_linear(lax, lay)

    return {
        "total_counts_df": total_df,
        "assigned_counts_df": assigned_df,
        "total_normalized_df": total_norm_df,
        "assigned_normalized_df": assigned_norm_df,
        "totals": {
            "xenium_total_sum": x_total_sum,
            "merscope_total_sum": m_total_sum,
            "xenium_assigned_sum": x_assigned_sum,
            "merscope_assigned_sum": m_assigned_sum,
        },
        "fits": {
            "total_log10": {
                "slope": total_fit[0],
                "intercept": total_fit[1],
                "r2": total_fit[2],
            },
            "assigned_log10": {
                "slope": assigned_fit[0],
                "intercept": assigned_fit[1],
                "r2": assigned_fit[2],
            },
        },
    }


def compute_gene_comparison_from_paths(
    xenium_zarr_path: Path | str,
    merscope_zarr_path: Path | str,
) -> dict[str, Any]:
    """Load SpatialData zarrs and run cross-platform gene comparison."""
    xenium_sdata = sd.read_zarr(Path(xenium_zarr_path))
    merscope_sdata = sd.read_zarr(Path(merscope_zarr_path))
    return compute_gene_comparison(
        xenium_sdata=xenium_sdata, merscope_sdata=merscope_sdata
    )
