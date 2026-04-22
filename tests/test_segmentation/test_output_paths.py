"""Tests for persistent segmentation output staging."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from merxen.config import SegmentationConfig
from merxen.segmentation.pipeline import run_segmentation_pipeline


def test_run_segmentation_pipeline_stages_persistent_outputs(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Persistent segmentation artifacts should be staged back into the work dir."""
    work_dir = tmp_path / "work" / "segment_out"
    persistent_root = tmp_path / "results"

    cfg = SegmentationConfig.model_validate(
        {
            "dataset": {
                "name": "P1_MERSCOPE",
                "platform": "MERSCOPE",
                "data_path": str(tmp_path / "input.zarr"),
                "channels": ["DAPI", "PolyT"],
                "output_dir": str(work_dir),
                "persistent_latest_zarr_path": str(
                    persistent_root / "latest" / "latest_spatialdata.zarr"
                ),
                "persistent_mask_path": str(
                    persistent_root / "segmentation" / "cellpose_masks_tiled.npy"
                ),
                "persistent_transcripts_path": str(
                    persistent_root / "segmentation" / "transcripts_for_proseg.csv"
                ),
            }
        }
    )

    points_df = pd.DataFrame({"x": [1.0], "y": [2.0], "gene": ["Gad1"]})

    monkeypatch.setattr(
        "merxen.segmentation.pipeline._load_dataset_sdata",
        lambda config: (object(), object(), 8, 8, np.eye(3), points_df),
    )
    monkeypatch.setattr(
        "merxen.segmentation.pipeline.run_tiled_cellpose",
        lambda *, output_mask_path, **kwargs: (
            output_mask_path.parent.mkdir(parents=True, exist_ok=True),
            np.save(output_mask_path, np.ones((4, 4), dtype=np.int32)),
            output_mask_path,
        )[-1],
    )
    monkeypatch.setattr(
        "merxen.segmentation.pipeline.build_cellpose_affine_to_microns",
        lambda *args, **kwargs: ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
    )
    monkeypatch.setattr(
        "merxen.segmentation.pipeline.write_proseg_csv_from_points",
        lambda *, csv_path, **kwargs: (
            csv_path.parent.mkdir(parents=True, exist_ok=True),
            csv_path.write_text("x_micron,y_micron,z_micron,feature_name,cell_id\n"),
            {"n_seeded": 1, "pct_seeded": 100.0},
        )[-1],
    )
    monkeypatch.setattr(
        "merxen.segmentation.pipeline.run_proseg_refinement",
        lambda *, output_path, **kwargs: (
            output_path.parent.mkdir(parents=True, exist_ok=True),
            output_path.mkdir(parents=True, exist_ok=True),
            output_path,
        )[-1],
    )
    monkeypatch.setattr(
        "merxen.segmentation.pipeline.convert_to_latest_zarr",
        lambda raw_path, latest_path: (
            latest_path.parent.mkdir(parents=True, exist_ok=True),
            latest_path.mkdir(parents=True, exist_ok=True),
            (latest_path / "marker.txt").write_text("latest"),
            latest_path,
        )[-1],
    )

    outputs = run_segmentation_pipeline(cfg, force_rerun=True)

    staged_latest = work_dir / "proseg_base_latest.zarr"
    staged_mask = work_dir / "cellpose_masks_tiled.npy"
    staged_transcripts = work_dir / "transcripts_for_proseg.csv"

    assert outputs["latest_output"] == staged_latest
    assert outputs["cellpose_mask_path"] == staged_mask
    assert outputs["transcripts_csv"] == staged_transcripts

    assert staged_latest.is_symlink()
    assert staged_mask.is_symlink()
    assert staged_transcripts.is_symlink()

    assert (
        staged_latest.resolve()
        == Path(cfg.dataset.persistent_latest_zarr_path).resolve()
    )
    assert staged_mask.resolve() == Path(cfg.dataset.persistent_mask_path).resolve()
    assert (
        staged_transcripts.resolve()
        == Path(cfg.dataset.persistent_transcripts_path).resolve()
    )

    assert not (work_dir / "proseg_base_raw.zarr").exists()
