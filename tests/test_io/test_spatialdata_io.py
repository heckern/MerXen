"""Tests for SpatialData write helpers."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

from merxen.io.spatialdata_io import write_spatialdata_zarr


def test_write_spatialdata_zarr_calls_write(tmp_path: Path) -> None:
    """write_spatialdata_zarr should delegate to sdata.write with the given path."""
    sdata = MagicMock()
    out = tmp_path / "out.zarr"

    write_spatialdata_zarr(sdata, out)

    sdata.write.assert_called_once_with(out)


def test_write_spatialdata_zarr_passes_overwrite_flag(tmp_path: Path) -> None:
    """write_spatialdata_zarr should forward the overwrite kwarg when supplied."""
    sdata = MagicMock()
    out = tmp_path / "out.zarr"

    write_spatialdata_zarr(sdata, out, overwrite=True)

    sdata.write.assert_called_once_with(out, overwrite=True)


def test_write_spatialdata_zarr_omits_overwrite_when_none(tmp_path: Path) -> None:
    """write_spatialdata_zarr should not pass overwrite when it is None."""
    sdata = MagicMock()
    out = tmp_path / "out.zarr"

    write_spatialdata_zarr(sdata, out, overwrite=None)

    sdata.write.assert_called_once_with(out)
