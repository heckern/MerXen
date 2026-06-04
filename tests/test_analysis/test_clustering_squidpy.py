"""Tests for the Scanpy/Squidpy clustering shim."""

from __future__ import annotations

import sys
import warnings
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from shapely.geometry import box

from merxen.analysis.clustering_squidpy import (
    AtlasMarkerSet,
    _add_spatial_scale_bar,
    _clean_spatial_axis,
    _make_neuron_split_marker_sets,
    _run_gpu_clustering,
    adata_from_spatialdata,
    collapse_atlas_label_to_broad_class,
    compute_group_gene_summary,
    load_atlas_marker_sets,
    plot_annotation_score_heatmap,
    plot_group_gene_dotplot,
    plot_spatial_cluster_grid,
    plot_spatial_scatter,
    remove_control_features,
    run_scanpy_clustering,
    score_clusters_by_atlas_markers,
)
from merxen.config import ClusteringSquidpyConfig


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


def test_adata_from_spatialdata_adds_ensembl_ids_from_original_table() -> None:
    """Gene IDs from one SpatialData table should annotate clustering tables."""
    obs = pd.DataFrame({"cell_id": ["x1", "x2"]}, index=["x1", "x2"])
    adata = ad.AnnData(
        X=np.array([[10, 1], [2, 8]], dtype=np.int64),
        obs=obs,
        var=pd.DataFrame({"gene": ["GeneA", "GeneB"]}, index=["GeneA", "GeneB"]),
    )
    adata.uns["spatialdata_attrs"] = {"region": "xenium_cell_boundaries"}
    original = ad.AnnData(
        X=np.ones((1, 2), dtype=np.float32),
        obs=pd.DataFrame(index=["cell0"]),
        var=pd.DataFrame(
            {
                "gene_ids": ["ENSG000001", "ENSG000002"],
                "feature_types": ["Gene Expression", "Gene Expression"],
            },
            index=["GeneA", "GeneB"],
        ),
    )
    cell_gdf = gpd.GeoDataFrame(
        {
            "cell_id": ["x1", "x2"],
            "geometry": [box(0, 0, 1, 1), box(2, 0, 3, 1)],
        },
        geometry="geometry",
    )
    fake_sdata = SimpleNamespace(
        tables={"table": adata, "table_original": original},
        shapes={"xenium_cell_boundaries": cell_gdf},
    )

    out = adata_from_spatialdata(fake_sdata, platform="XENIUM")

    assert list(out.var["ensembl_id"]) == ["ENSG000001", "ENSG000002"]
    assert out.uns["merxen_clustering_squidpy"]["ensembl_id_mapping"] == {
        "n_features": 2,
        "n_mapped": 2,
        "column": "ensembl_id",
    }


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
        normalize_exclude_highly_expressed=False,
        n_pcs=3,
        n_neighbors=3,
        umap_min_dist=0.2,
        umap_spread=1.5,
        random_seed=1,
        use_gpu=False,
    )

    assert "counts" in out.layers
    assert "X_umap" in out.obsm
    assert "leiden" in out.obs
    assert out.uns["merxen_clustering_params"]["umap_min_dist"] == 0.2
    assert out.uns["merxen_clustering_params"]["umap_spread"] == 1.5


def test_run_scanpy_clustering_can_start_from_counts_layer() -> None:
    """Branch reclustering should renormalize from raw counts, not parent log X."""
    rng = np.random.default_rng(11)
    counts = rng.poisson(lam=4, size=(12, 5)).astype(np.float32)
    adata = ad.AnnData(
        X=np.zeros_like(counts),
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(12)]),
        var=pd.DataFrame(index=[f"Gene{i}" for i in range(5)]),
        layers={"counts": counts.copy()},
    )
    adata.obsm["spatial"] = rng.normal(size=(12, 2))

    out = run_scanpy_clustering(
        adata,
        min_counts=1,
        min_cells=1,
        n_pcs=2,
        n_neighbors=3,
        random_seed=11,
        use_gpu=False,
        key_added="leiden_branch",
        input_layer="counts",
    )

    assert "leiden_branch" in out.obs
    assert "leiden" not in out.obs
    assert float(out.layers["counts"].sum()) > 0.0
    assert out.uns["merxen_clustering_params_leiden_branch"]["input_layer"] == (
        "counts"
    )


def test_run_scanpy_clustering_preserves_ensembl_ids() -> None:
    """Filtering and clustering should retain gene IDs needed by MapMyCells."""
    rng = np.random.default_rng(2)
    adata = ad.AnnData(
        X=rng.poisson(lam=4, size=(12, 4)).astype(np.float32),
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(12)]),
        var=pd.DataFrame(
            {
                "gene": ["GeneA", "GeneB", "Blank-1", "GeneC"],
                "ensembl_id": ["ENSG000001", "ENSG000002", "", "ENSG000003"],
            },
            index=["GeneA", "GeneB", "Blank-1", "GeneC"],
        ),
    )
    adata.obsm["spatial"] = rng.normal(size=(12, 2))

    out = run_scanpy_clustering(
        adata,
        min_counts=1,
        min_cells=1,
        n_pcs=2,
        n_neighbors=3,
        random_seed=2,
        use_gpu=False,
    )

    assert list(out.var["ensembl_id"]) == [
        "ENSG000001",
        "ENSG000002",
        "ENSG000003",
    ]


def test_remove_control_features_drops_blank_negative_and_unassigned() -> None:
    """Control-like features should be excluded from clustering inputs."""
    adata = ad.AnnData(
        X=np.ones((4, 5), dtype=np.float32),
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(4)]),
        var=pd.DataFrame(
            index=[
                "GeneA",
                "Blank-1",
                "NegControlProbe_00001",
                "UnassignedCodeword_0001",
                "GeneB",
            ]
        ),
    )

    filtered = remove_control_features(adata)

    assert list(filtered.var_names) == ["GeneA", "GeneB"]
    summary = filtered.uns["merxen_clustering_squidpy"]["control_feature_filter"]
    assert summary["n_features_before"] == 5
    assert summary["n_control_features_removed"] == 3
    assert summary["removed_control_features"] == [
        "Blank-1",
        "NegControlProbe_00001",
        "UnassignedCodeword_0001",
    ]


def test_load_atlas_marker_sets_joins_taxonomy_and_collapses_labels(
    tmp_path: Path,
) -> None:
    """MapMyCells marker keys should resolve through Allen taxonomy metadata."""
    marker_path = tmp_path / "markers.json"
    marker_path.write_text(
        """
{
  "CCN202210140_SUPC/CS1": ["ENSG1", "ENSG2"],
  "CCN202210140_SUPC/CS2": ["ENSG3"],
  "CCN202210140_SUPC/CS3": ["ENSG4"],
  "CCN202210140_SUBC/ignored": ["ENSG4"]
}
""".strip()
    )
    taxonomy_path = tmp_path / "cluster_annotation_term.csv"
    taxonomy_path.write_text(
        "\n".join(
            [
                "label,name,cluster_annotation_term_set_label",
                "CS1,Oligodendrocyte,CCN202210140_SUPC",
                "CS2,Upper-layer intratelencephalic,CCN202210140_SUPC",
                "CS3,MGE interneuron,CCN202210140_SUPC",
                "ignored,Ignored,CCN202210140_SUBC",
            ]
        )
        + "\n"
    )
    membership_path = tmp_path / "cluster_to_cluster_annotation_membership.csv"
    membership_path.write_text(
        "\n".join(
            [
                "cluster_annotation_term_label,cluster_annotation_term_set_label,"
                "cluster_alias,cluster_annotation_term_name",
                "CS2,CCN202210140_SUPC,1,Upper-layer intratelencephalic",
                "CS202210140_3820,CCN202210140_NEUR,1,VGLUT1",
                "CS3,CCN202210140_SUPC,2,MGE interneuron",
                "CS202210140_3810,CCN202210140_NEUR,2,GABA",
            ]
        )
        + "\n"
    )

    marker_sets = load_atlas_marker_sets(
        marker_path,
        taxonomy_path,
        cluster_membership_path=membership_path,
    )

    assert [marker_set.label_name for marker_set in marker_sets] == [
        "Oligodendrocyte",
        "Upper-layer intratelencephalic",
        "MGE interneuron",
    ]
    assert marker_sets[0].broad_class == "Oligodendrocytes"
    assert marker_sets[1].broad_class == "Neurons"
    assert marker_sets[1].neuron_split == "Excitatory"
    assert marker_sets[2].neuron_split == "Inhibitory"
    assert collapse_atlas_label_to_broad_class("Choroid plexus") == ("Choroid plexus")


def test_score_clusters_by_atlas_markers_resolves_ensembl_then_symbol() -> None:
    """Synthetic marker expression should recover known broad labels."""
    adata = ad.AnnData(
        X=np.array(
            [
                [9, 8, 1, 1],
                [8, 7, 1, 1],
                [1, 1, 8, 9],
                [1, 1, 7, 8],
            ],
            dtype=np.float32,
        ),
        obs=pd.DataFrame({"leiden_broad": ["0", "0", "1", "1"]}),
        var=pd.DataFrame(
            {
                "gene": ["GeneA", "GeneB", "GeneC", "GeneD"],
                "ensembl_id": ["ENSGA", "ENSGB", "ENSGC", "ENSGD"],
            },
            index=["GeneA", "GeneB", "GeneC", "GeneD"],
        ),
    )
    marker_sets = [
        AtlasMarkerSet(
            level="level",
            label_id="oligo",
            label_name="Oligodendrocyte",
            broad_class="Oligodendrocytes",
            marker_ids=("ENSGA", "ENSGB"),
        ),
        AtlasMarkerSet(
            level="level",
            label_id="astro",
            label_name="Astrocyte",
            broad_class="Astrocytes",
            marker_ids=("GeneC", "GeneD"),
        ),
    ]

    assignments, scores, markers = score_clusters_by_atlas_markers(
        adata,
        cluster_key="leiden_broad",
        marker_sets=marker_sets,
        min_marker_overlap=2,
    )

    label_by_cluster = dict(
        zip(assignments["cluster"], assignments["atlas_label"], strict=True)
    )
    assert label_by_cluster == {"0": "Oligodendrocyte", "1": "Astrocyte"}
    assert set(scores["atlas_label"]) == {"Oligodendrocyte", "Astrocyte"}
    assert set(markers["n_resolved_markers"]) == {2}


def test_score_clusters_by_atlas_markers_uses_marker_alias_lookup() -> None:
    """Reference gene metadata should bridge Ensembl markers to symbol panels."""
    adata = ad.AnnData(
        X=np.array(
            [
                [9, 8, 1, 1],
                [8, 7, 1, 1],
                [1, 1, 8, 9],
                [1, 1, 7, 8],
            ],
            dtype=np.float32,
        ),
        obs=pd.DataFrame({"leiden_broad": ["0", "0", "1", "1"]}),
        var=pd.DataFrame(index=["GeneA", "GeneB", "GeneC", "GeneD"]),
    )
    marker_sets = [
        AtlasMarkerSet(
            level="level",
            label_id="oligo",
            label_name="Oligodendrocyte",
            broad_class="Oligodendrocytes",
            marker_ids=("ENSGA", "ENSGB"),
        ),
        AtlasMarkerSet(
            level="level",
            label_id="astro",
            label_name="Astrocyte",
            broad_class="Astrocytes",
            marker_ids=("ENSGC", "ENSGD"),
        ),
    ]

    assignments, _, markers = score_clusters_by_atlas_markers(
        adata,
        cluster_key="leiden_broad",
        marker_sets=marker_sets,
        marker_alias_lookup={
            "ENSGA": "GeneA",
            "ENSGB": "GeneB",
            "ENSGC": "GeneC",
            "ENSGD": "GeneD",
        },
        min_marker_overlap=2,
    )

    label_by_cluster = dict(
        zip(assignments["cluster"], assignments["atlas_label"], strict=True)
    )
    assert label_by_cluster == {"0": "Oligodendrocyte", "1": "Astrocyte"}
    assert list(markers["n_resolved_markers"]) == [2, 2]


def test_score_clusters_by_atlas_markers_unknown_for_low_overlap() -> None:
    """Panels with too little marker overlap should still get stable outputs."""
    adata = ad.AnnData(
        X=np.ones((4, 2), dtype=np.float32),
        obs=pd.DataFrame({"leiden_broad": ["0", "0", "1", "1"]}),
        var=pd.DataFrame(index=["GeneA", "GeneB"]),
    )

    assignments, scores, markers = score_clusters_by_atlas_markers(
        adata,
        cluster_key="leiden_broad",
        marker_sets=[
            AtlasMarkerSet(
                level="level",
                label_id="missing",
                label_name="Microglia",
                broad_class="Microglia",
                marker_ids=("MissingGene",),
            )
        ],
        min_marker_overlap=2,
        unknown_label="Mixed/Unknown",
    )

    assert list(assignments["atlas_label"]) == ["Mixed/Unknown", "Mixed/Unknown"]
    assert scores.empty
    assert list(scores.columns) == [
        "cluster",
        "label_id",
        "atlas_label",
        "broad_class",
        "score",
        "n_markers",
        "resolved_markers",
    ]
    assert list(markers["n_resolved_markers"]) == [0]


def test_neuron_split_marker_sets_group_exc_inh_and_other() -> None:
    """Neuron supercluster marker sets should collapse to Exc/Inh/Other groups."""
    marker_sets = [
        AtlasMarkerSet(
            level="level",
            label_id="exc",
            label_name="Upper-layer intratelencephalic",
            broad_class="Neurons",
            marker_ids=("ENSG1",),
        ),
        AtlasMarkerSet(
            level="level",
            label_id="inh",
            label_name="MGE interneuron",
            broad_class="Neurons",
            marker_ids=("ENSG2",),
        ),
        AtlasMarkerSet(
            level="level",
            label_id="other",
            label_name="Unclassified neuron",
            broad_class="Neurons",
            marker_ids=("ENSG3",),
        ),
    ]

    split_sets = _make_neuron_split_marker_sets(marker_sets)

    markers_by_split = {
        marker_set.label_name: marker_set.marker_ids for marker_set in split_sets
    }
    assert markers_by_split == {
        "Excitatory": ("ENSG1",),
        "Inhibitory": ("ENSG2",),
        "Other": ("ENSG3",),
    }


def test_plot_annotation_score_heatmap_writes_png_and_pdf(tmp_path: Path) -> None:
    """Annotation heatmaps should be emitted as regular plot artifacts."""
    score_table = pd.DataFrame(
        {
            "cluster": ["0", "1"],
            "atlas_label": ["Astrocyte", "Microglia"],
            "score": [1.2, 0.9],
        }
    )

    output_path = plot_annotation_score_heatmap(
        score_table,
        tmp_path / "scores.png",
        title="Synthetic scores",
    )

    assert output_path.exists()
    assert output_path.with_suffix(".pdf").exists()


def test_clustering_squidpy_config_defaults_enable_hierarchical_mode() -> None:
    """Minimal stage configs should run broad annotation and subclustering."""
    cfg = ClusteringSquidpyConfig.model_validate(
        {
            "pair_id": "pair1",
            "output_dir": "/tmp/out",
            "samples": [
                {
                    "sample_id": "sample1",
                    "platform": "MERSCOPE",
                    "zarr_path": "/tmp/input.zarr",
                }
            ],
        }
    )

    assert cfg.hierarchical_enabled is True
    assert cfg.leiden_resolution == 0.5
    assert cfg.broad_round.leiden_resolution == 0.2
    assert cfg.subcluster_round.leiden_resolution == 0.5
    assert cfg.spatial_point_size == 0.5
    assert cfg.spatial_scatter_point_size == 2.0


def test_run_gpu_clustering_uses_chunked_pca_for_sparse_input(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Sparse GPU PCA should avoid rapids-singlecell's fragile sparse helper."""
    adata = ad.AnnData(
        X=sparse.csr_matrix(np.ones((20, 6), dtype=np.float32)),
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(20)]),
        var=pd.DataFrame(index=[f"Gene{i}" for i in range(6)]),
    )
    calls: dict[str, object] = {}

    fake_get = SimpleNamespace(
        anndata_to_GPU=lambda data: calls.setdefault("to_gpu", data),
        anndata_to_CPU=lambda data: calls.setdefault("to_cpu", data),
    )

    def fake_pca(data: ad.AnnData, **kwargs: object) -> None:
        calls["pca"] = kwargs
        data.obsm["X_pca"] = np.ones((data.n_obs, 3), dtype=np.float32)

    def fake_neighbors(data: ad.AnnData, **kwargs: object) -> None:
        calls["neighbors"] = kwargs

    fake_pp = SimpleNamespace(pca=fake_pca, neighbors=fake_neighbors)
    fake_tl = SimpleNamespace(
        umap=lambda data, **kwargs: calls.setdefault("umap", kwargs),
        leiden=lambda data, **kwargs: calls.setdefault("leiden", kwargs),
    )
    fake_rsc = SimpleNamespace(get=fake_get, pp=fake_pp, tl=fake_tl)
    monkeypatch.setitem(sys.modules, "rapids_singlecell", fake_rsc)

    gpu_used = _run_gpu_clustering(
        adata,
        max_pcs=3,
        n_pcs_for_neighbors=3,
        effective_neighbors=5,
        umap_min_dist=0.4,
        umap_spread=1.2,
        leiden_resolution=0.8,
        random_seed=7,
    )

    assert gpu_used is True
    assert calls["to_gpu"] is adata
    assert calls["to_cpu"] is adata
    assert calls["pca"] == {
        "n_comps": 3,
        "random_state": 7,
        "chunked": True,
        "chunk_size": adata.n_obs,
    }
    assert calls["neighbors"] == {
        "n_neighbors": 5,
        "n_pcs": 3,
        "use_rep": "X_pca",
        "random_state": 7,
    }


def test_plot_spatial_scatter_suppresses_squidpy_noise(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Image-less spatial scatter should not emit Squidpy library warnings."""
    adata = ad.AnnData(
        X=np.ones((5, 2), dtype=np.float32),
        obs=pd.DataFrame(
            {"leiden": pd.Categorical(["0", "1", "0", "1", "2"])},
            index=[f"cell{i}" for i in range(5)],
        ),
        var=pd.DataFrame(index=["Gene0", "Gene1"]),
    )
    adata.obsm["spatial"] = np.column_stack(
        [np.arange(adata.n_obs), np.arange(adata.n_obs)]
    )

    with warnings.catch_warnings(record=True) as recorded:
        warnings.simplefilter("always")
        output_path = plot_spatial_scatter(
            adata,
            tmp_path / "spatial.png",
            point_size=0.2,
        )

    captured = capsys.readouterr()
    warning_text = "\n".join(str(item.message) for item in recorded)
    output_text = f"{captured.out}\n{captured.err}"
    assert output_path.exists()
    assert output_path.with_suffix(".pdf").exists()
    assert "No data for colormapping provided via 'c'" not in warning_text
    assert "Please specify a valid `library_id`" not in output_text


def test_spatial_axis_cleanup_adds_scale_bar_without_coordinate_labels() -> None:
    """Spatial helper should hide raw coordinate axes and add a 200 um scale bar."""
    fig, ax = plt.subplots()
    ax.scatter([0, 100, 300], [0, 50, 100])
    ax.set_xlabel("spatial1")
    ax.set_ylabel("spatial2")

    _clean_spatial_axis(ax)
    _add_spatial_scale_bar(ax, length_um=200)

    assert ax.get_xlabel() == ""
    assert ax.get_ylabel() == ""
    assert not ax.get_xticks().size
    assert not ax.get_yticks().size
    assert any(text.get_text() == "200 um" for text in ax.texts)
    plt.close(fig)


def test_plot_spatial_cluster_grid_writes_png_and_pdf(tmp_path: Path) -> None:
    """Spatial cluster grid should highlight each Leiden cluster separately."""
    adata = ad.AnnData(
        X=np.ones((6, 2), dtype=np.float32),
        obs=pd.DataFrame(
            {"leiden": pd.Categorical(["0", "1", "0", "1", "2", "2"])},
            index=[f"cell{i}" for i in range(6)],
        ),
        var=pd.DataFrame(index=["Gene0", "Gene1"]),
    )
    adata.obsm["spatial"] = np.column_stack(
        [np.arange(adata.n_obs), np.arange(adata.n_obs)]
    )

    output_path = plot_spatial_cluster_grid(
        adata,
        tmp_path / "spatial_leiden_grid.png",
        point_size_highlight=0.4,
    )

    assert output_path.exists()
    assert output_path.with_suffix(".pdf").exists()


def test_group_gene_dotplot_writes_summary_plot(tmp_path: Path) -> None:
    """Branch dotplot helpers should summarize mean and fraction by group."""
    adata = ad.AnnData(
        X=np.array(
            [
                [3.0, 0.0, 0.0],
                [1.0, 2.0, 0.0],
                [0.0, 0.0, 5.0],
                [0.0, 1.0, 4.0],
            ],
            dtype=np.float32,
        ),
        obs=pd.DataFrame(
            {"leiden_subcluster": ["0", "0", "1", "1"]},
            index=[f"cell{i}" for i in range(4)],
        ),
        var=pd.DataFrame(index=["GeneA", "GeneB", "GeneC"]),
    )

    mean_expression, fraction_expression = compute_group_gene_summary(
        adata,
        ["GeneA", "GeneB", "GeneC"],
        groupby="leiden_subcluster",
    )
    output_path = plot_group_gene_dotplot(
        mean_expression,
        fraction_expression,
        tmp_path / "gene_dotplot.png",
    )

    assert mean_expression.loc["0", "GeneA"] == pytest.approx(2.0)
    assert mean_expression.loc["1", "GeneC"] == pytest.approx(4.5)
    assert fraction_expression.loc["0", "GeneB"] == pytest.approx(0.5)
    assert fraction_expression.loc["1", "GeneC"] == pytest.approx(1.0)
    assert output_path.exists()
    assert output_path.with_suffix(".pdf").exists()
