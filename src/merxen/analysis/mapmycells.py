"""Local MapMyCells annotation for clustered MerXen AnnData outputs."""

from __future__ import annotations

import json
import logging
import subprocess
import sys
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from merxen.config import MapMyCellsConfig
from merxen.memory import force_release, log_status

logger = logging.getLogger(__name__)


def prepare_mapmycells_query(
    input_h5ad: Path | str,
    output_h5ad: Path | str,
    *,
    query_layer: str | None = "counts",
    gene_id_column: str | None = None,
    obs_id_column: str | None = None,
) -> Path:
    """Write a MapMyCells-ready H5AD query file.

    The Squidpy clustering stage leaves normalized/log-transformed values in
    ``X`` and preserves raw counts in ``layers["counts"]``. MapMyCells expects
    the query matrix in ``X``, so this helper copies the selected layer into
    ``X`` before writing a local query file.

    Args:
        input_h5ad: Clustered AnnData from ``clustering_squidpy``.
        output_h5ad: Destination H5AD consumed by MapMyCells.
        query_layer: AnnData layer to copy into ``X``. Use ``None`` to keep the
            current ``X`` matrix.
        gene_id_column: Optional ``var`` column to use as gene identifiers.
        obs_id_column: Optional ``obs`` column to use as cell identifiers.

    Returns:
        Path to the written query H5AD.
    """
    input_h5ad = Path(input_h5ad)
    output_h5ad = Path(output_h5ad)
    output_h5ad.parent.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(input_h5ad)
    try:
        if query_layer is not None:
            if query_layer not in adata.layers:
                raise KeyError(
                    f"Requested query_layer={query_layer!r} not found in {input_h5ad}. "
                    f"Available layers: {list(adata.layers.keys())}"
                )
            adata.X = _copy_matrix(adata.layers[query_layer])

        if gene_id_column is not None:
            if gene_id_column not in adata.var.columns:
                raise KeyError(
                    f"Requested gene_id_column={gene_id_column!r} not found in "
                    f"{input_h5ad}. Available var columns: {list(adata.var.columns)}"
                )
            adata.var_names = pd.Index(
                adata.var[gene_id_column].astype(str),
                name=gene_id_column,
            )

        if obs_id_column is not None:
            if obs_id_column not in adata.obs.columns:
                raise KeyError(
                    f"Requested obs_id_column={obs_id_column!r} not found in "
                    f"{input_h5ad}. Available obs columns: {list(adata.obs.columns)}"
                )
            adata.obs_names = pd.Index(
                adata.obs[obs_id_column].astype(str),
                name=obs_id_column,
            )

        adata.var_names = pd.Index(
            adata.var_names.astype(str), name=adata.var_names.name
        )
        adata.obs_names = pd.Index(
            adata.obs_names.astype(str), name=adata.obs_names.name
        )
        adata.var_names_make_unique()
        adata.obs_names_make_unique()
        adata.write_h5ad(output_h5ad)
    finally:
        del adata
        force_release(note=f"after preparing MapMyCells query {input_h5ad.name}")

    return output_h5ad


def run_mapmycells(config: MapMyCellsConfig) -> dict[str, dict[str, Path]]:
    """Run local MapMyCells assignment for every sample in a pair.

    Args:
        config: Validated MapMyCells stage configuration.

    Returns:
        Mapping from sample ID to output artifact paths.
    """
    _require_existing_file(config.marker_lookup_path, "MapMyCells marker lookup")
    _require_existing_file(
        config.precomputed_stats_path, "MapMyCells precomputed stats"
    )
    if config.tmp_dir is not None:
        config.tmp_dir.mkdir(parents=True, exist_ok=True)

    config.output_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, dict[str, Path]] = {}

    for sample in config.samples:
        log_status(
            f"[{sample.sample_id}] Starting MapMyCells "
            f"(platform={sample.platform}, bootstrap_factor={config.bootstrap_factor})"
        )
        sample_dir = config.output_dir / sample.platform.lower()
        sample_dir.mkdir(parents=True, exist_ok=True)
        _require_existing_file(
            sample.anndata_path, f"clustered AnnData for {sample.sample_id}"
        )

        query_h5ad = prepare_mapmycells_query(
            sample.anndata_path,
            sample_dir / f"{sample.sample_id}_mapmycells_query.h5ad",
            query_layer=sample.query_layer,
            gene_id_column=sample.gene_id_column,
            obs_id_column=sample.obs_id_column,
        )
        extended_json = sample_dir / f"{sample.sample_id}_mapmycells_extended.json"
        csv_path = sample_dir / f"{sample.sample_id}_mapmycells.csv"
        log_path = sample_dir / f"{sample.sample_id}_mapmycells.log"
        command_path = sample_dir / f"{sample.sample_id}_mapmycells_command.json"
        annotated_h5ad = sample_dir / f"{sample.sample_id}_mapmycells_annotated.h5ad"

        command = build_mapmycells_command(
            config,
            query_h5ad=query_h5ad,
            extended_json=extended_json,
            csv_path=csv_path,
            log_path=log_path,
        )
        _write_command_manifest(command_path, command)
        _run_command(command)
        annotate_h5ad_with_mapmycells(
            sample.anndata_path,
            csv_path,
            annotated_h5ad,
        )

        results[sample.sample_id] = {
            "query_h5ad": query_h5ad,
            "extended_json": extended_json,
            "csv": csv_path,
            "log": log_path,
            "command_json": command_path,
            "annotated_h5ad": annotated_h5ad,
        }
        force_release(note=f"after MapMyCells {sample.sample_id}")

    manifest_path = config.output_dir / f"{config.pair_id}_mapmycells_manifest.json"
    _write_results_manifest(manifest_path, config, results)
    return results


def build_mapmycells_command(
    config: MapMyCellsConfig,
    *,
    query_h5ad: Path,
    extended_json: Path,
    csv_path: Path,
    log_path: Path,
) -> list[str]:
    """Build the ``cell_type_mapper`` command-line invocation."""
    command = [
        sys.executable,
        "-m",
        "cell_type_mapper.cli.from_specified_markers",
        "--query_path",
        str(query_h5ad),
        "--extended_result_path",
        str(extended_json),
        "--csv_result_path",
        str(csv_path),
        "--log_path",
        str(log_path),
        "--cloud_safe",
        _bool_arg(config.cloud_safe),
        "--query_markers.serialized_lookup",
        str(config.marker_lookup_path),
        "--precomputed_stats.path",
        str(config.precomputed_stats_path),
        "--type_assignment.normalization",
        config.normalization,
        "--type_assignment.bootstrap_iteration",
        str(config.bootstrap_iteration),
        "--type_assignment.bootstrap_factor",
        str(config.bootstrap_factor),
        "--type_assignment.n_processors",
        str(config.n_processors),
        "--flatten",
        _bool_arg(config.flatten),
    ]
    if config.drop_level is not None:
        command.extend(["--drop_level", config.drop_level])
    if config.chunk_size is not None:
        command.extend(["--type_assignment.chunk_size", str(config.chunk_size)])
    if config.rng_seed is not None:
        command.extend(["--type_assignment.rng_seed", str(config.rng_seed)])
    if config.max_gb is not None:
        command.extend(["--max_gb", str(config.max_gb)])
    if config.tmp_dir is not None:
        command.extend(["--tmp_dir", str(config.tmp_dir)])
    if config.verbose_csv:
        command.extend(["--verbose_csv", _bool_arg(config.verbose_csv)])
    command.extend(config.extra_args)
    return command


def annotate_h5ad_with_mapmycells(
    input_h5ad: Path | str,
    csv_path: Path | str,
    output_h5ad: Path | str,
) -> Path:
    """Attach MapMyCells CSV assignments to ``adata.obs`` and write H5AD."""
    input_h5ad = Path(input_h5ad)
    csv_path = Path(csv_path)
    output_h5ad = Path(output_h5ad)
    output_h5ad.parent.mkdir(parents=True, exist_ok=True)

    assignments = read_mapmycells_csv(csv_path)
    adata = ad.read_h5ad(input_h5ad)
    try:
        indexed = assignments.set_index(assignments.columns[0], drop=False)
        indexed.index = indexed.index.astype(str)
        indexed = indexed[~indexed.index.duplicated(keep="first")]
        aligned = indexed.reindex(adata.obs_names.astype(str))
        for column in aligned.columns:
            target = f"mapmycells_{column}"
            adata.obs[target] = aligned[column].to_numpy()

        matched = int(aligned.iloc[:, 0].notna().sum()) if not aligned.empty else 0
        adata.uns["merxen_mapmycells"] = {
            "csv_path": str(csv_path),
            "n_assignments": int(len(assignments)),
            "n_obs": int(adata.n_obs),
            "n_matched_obs": matched,
        }
        adata.write_h5ad(output_h5ad)
    finally:
        del adata
        force_release(note=f"after annotating MapMyCells output {input_h5ad.name}")

    return output_h5ad


def read_mapmycells_csv(csv_path: Path | str) -> pd.DataFrame:
    """Read the comment-prefixed MapMyCells CSV output."""
    csv_path = Path(csv_path)
    df = pd.read_csv(csv_path, comment="#", converters={0: str})
    if df.empty:
        raise ValueError(f"MapMyCells CSV is empty: {csv_path}")
    first_col = str(df.columns[0])
    df[first_col] = df[first_col].astype(str)
    return df


def _copy_matrix(matrix: Any) -> Any:
    if sparse.issparse(matrix):
        return matrix.copy()
    return np.array(matrix, copy=True)


def _run_command(command: list[str]) -> None:
    logger.info("Running MapMyCells command: %s", " ".join(command))
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            "Could not start MapMyCells. Install the Allen Institute "
            "cell_type_mapper package in the active environment."
        ) from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"MapMyCells failed with exit code {exc.returncode}"
        ) from exc


def _require_existing_file(path: Path, label: str) -> None:
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")


def _write_command_manifest(path: Path, command: list[str]) -> None:
    path.write_text(json.dumps({"command": command}, indent=2) + "\n")


def _write_results_manifest(
    path: Path,
    config: MapMyCellsConfig,
    results: dict[str, dict[str, Path]],
) -> None:
    payload = {
        "pair_id": config.pair_id,
        "marker_lookup_path": str(config.marker_lookup_path),
        "precomputed_stats_path": str(config.precomputed_stats_path),
        "bootstrap_factor": config.bootstrap_factor,
        "bootstrap_iteration": config.bootstrap_iteration,
        "n_processors": config.n_processors,
        "samples": {
            sample_id: {key: str(value) for key, value in paths.items()}
            for sample_id, paths in results.items()
        },
    }
    path.write_text(json.dumps(payload, indent=2) + "\n")


def _bool_arg(value: bool) -> str:
    return "True" if value else "False"
