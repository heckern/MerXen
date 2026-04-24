"""Tests for image source normalization helpers."""

from __future__ import annotations

from types import SimpleNamespace

import dask.array as da
import numpy as np
import xarray as xr

from merxen.io.image_source import _get_image_dataarray, build_image_source


def test_get_image_dataarray_keeps_dask_backed_dataarray() -> None:
    """Plain DataArray inputs should not be probed as multiscale containers."""
    arr = xr.DataArray(
        da.from_array(np.arange(24).reshape(2, 3, 4), chunks=(1, 3, 4)),
        dims=("c", "y", "x"),
    )

    result = _get_image_dataarray(arr)

    assert result is arr


def test_get_image_dataarray_extracts_scale0_image() -> None:
    """Multiscale containers should still resolve to their scale0 image."""
    image = xr.DataArray(np.ones((1, 4, 5)), dims=("c", "y", "x"))
    multiscale = {"scale0": SimpleNamespace(ds={"image": image})}

    result = _get_image_dataarray(multiscale)

    assert result is image


def test_build_image_source_accepts_dask_backed_dataarray() -> None:
    """Image source construction should succeed for dask-backed DataArrays."""
    arr = xr.DataArray(
        da.from_array(np.arange(24).reshape(2, 3, 4), chunks=(1, 3, 4)),
        dims=("c", "y", "x"),
    )

    source = build_image_source(arr, as_float32=False)

    assert source["kind"] == "xarray"
    assert source["shape"] == (3, 4, 2)
