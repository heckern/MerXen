"""Downstream analysis modules for MerXen."""

from __future__ import annotations

from merxen.analysis.clustering_squidpy import (
    adata_from_spatialdata,
    add_qc_metrics,
    load_spatialdata_adata,
    plot_qc_histograms,
    plot_spatial_scatter,
    plot_umap,
    run_clustering_squidpy,
    run_scanpy_clustering,
    save_clustered_adata,
    save_qc_metrics,
)

__all__ = [
    "add_qc_metrics",
    "adata_from_spatialdata",
    "load_spatialdata_adata",
    "plot_qc_histograms",
    "plot_spatial_scatter",
    "plot_umap",
    "run_clustering_squidpy",
    "run_scanpy_clustering",
    "save_clustered_adata",
    "save_qc_metrics",
]
