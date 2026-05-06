"""Tests for the Scanpy/Squidpy clustering shim."""

from __future__ import annotations

from types import SimpleNamespace

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

from merxen.analysis.clustering_squidpy import (
    adata_from_spatialdata,
    run_scanpy_clustering,
)


def test_adata_from_spatialdata_adds_spatial_area_and_control_metrics() -> None:
    """SpatialData table extraction should add Squidpy-ready coordinates."""
    obs = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4"],
            "control_probe_counts": [1, 0, 2, 0],
        },
        index=["c1", "c2", "c3", "c4"],
    )
    var = pd.DataFrame(index=["GeneA", "Blank-1", "GeneB", "NegControlProbe-1"])
    adata = ad.AnnData(
        X=np.array(
            [
                [10, 1, 0, 2],
                [0, 0, 12, 1],
                [3, 4, 5, 0],
                [6, 0, 0, 0],
            ],
            dtype=np.int64,
        ),
        obs=obs,
        var=var,
    )
    adata.obsm["blank"] = pd.DataFrame(
        {"Blank-A": [5, 0, 1, 0]},
        index=adata.obs_names,
    )
    adata.uns["spatialdata_attrs"] = {"region": "MOSAIK_proseg"}

    gdf = gpd.GeoDataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4"],
            "geometry": [
                box(0, 0, 1, 1),
                box(2, 0, 3, 1),
                box(0, 2, 1, 3),
                box(2, 2, 3, 3),
            ],
        },
        geometry="geometry",
    )
    aligned_gdf = gdf.copy()
    aligned_gdf["geometry"] = aligned_gdf.geometry.translate(xoff=10.0)
    fake_sdata = SimpleNamespace(
        tables={"table": adata},
        shapes={
            "MOSAIK_proseg": gdf,
            "MOSAIK_proseg_aligned_nonrigid": aligned_gdf,
        },
    )

    out = adata_from_spatialdata(fake_sdata, platform="MERSCOPE")

    assert out.uns["merxen_clustering_squidpy"]["shape_key"] == (
        "MOSAIK_proseg_aligned_nonrigid"
    )
    np.testing.assert_allclose(out.obsm["spatial"][0], [10.5, 0.5])
    np.testing.assert_allclose(out.obs["cell_area"].to_numpy(float), 1.0)
    np.testing.assert_allclose(
        out.obs["control_counts"].to_numpy(float),
        [9.0, 1.0, 7.0, 0.0],
    )
    assert out.obs["nucleus_ratio"].isna().all()


def test_adata_from_spatialdata_adds_xenium_nucleus_ratio_from_shapes() -> None:
    """Xenium nucleus shapes should fill nucleus_area when tables lack it."""
    obs = pd.DataFrame(
        {"cell_id": ["x1", "x2"]},
        index=["x1", "x2"],
    )
    adata = ad.AnnData(
        X=np.array([[10, 1], [2, 8]], dtype=np.int64),
        obs=obs,
        var=pd.DataFrame(index=["GeneA", "GeneB"]),
    )
    adata.uns["spatialdata_attrs"] = {"region": "xenium_cell_boundaries"}
    cell_gdf = gpd.GeoDataFrame(
        {
            "cell_id": ["x1", "x2"],
            "geometry": [box(0, 0, 2, 2), box(4, 0, 6, 2)],
        },
        geometry="geometry",
    )
    nucleus_gdf = gpd.GeoDataFrame(
        {
            "cell_id": ["x1", "x2"],
            "geometry": [box(0, 0, 1, 1), box(4, 0, 5, 1)],
        },
        geometry="geometry",
    )
    fake_sdata = SimpleNamespace(
        tables={"table": adata},
        shapes={
            "xenium_cell_boundaries": cell_gdf,
            "xenium_nucleus": nucleus_gdf,
        },
    )

    out = adata_from_spatialdata(fake_sdata, platform="XENIUM")

    np.testing.assert_allclose(out.obs["cell_area"].to_numpy(float), [4.0, 4.0])
    np.testing.assert_allclose(out.obs["nucleus_area"].to_numpy(float), [1.0, 1.0])
    np.testing.assert_allclose(out.obs["nucleus_ratio"].to_numpy(float), [0.25, 0.25])


def test_run_scanpy_clustering_adds_umap_and_leiden() -> None:
    """The gentle Scanpy workflow should produce expected clustering fields."""
    rng = np.random.default_rng(1)
    adata = ad.AnnData(
        X=rng.poisson(lam=4, size=(12, 6)).astype(np.float32),
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(12)]),
        var=pd.DataFrame(index=[f"Gene{i}" for i in range(6)]),
    )
    adata.obsm["spatial"] = rng.normal(size=(12, 2))

    out = run_scanpy_clustering(
        adata,
        min_counts=1,
        min_cells=1,
        n_pcs=3,
        n_neighbors=3,
        random_seed=1,
    )

    assert "counts" in out.layers
    assert "X_umap" in out.obsm
    assert "leiden" in out.obs
