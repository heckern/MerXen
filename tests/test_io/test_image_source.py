"""Tests for image source normalization helpers."""

from __future__ import annotations

from types import SimpleNamespace

import dask.array as da
import numpy as np
import xarray as xr

from merxen.io.image_source import (
    MERSCOPE_ZPROJ_IMAGE_NAME,
    _get_image_dataarray,
    build_image_source,
    build_merscope_z_projection,
    image_to_cyx,
    max_project_image_elements,
)


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


def test_get_image_dataarray_extracts_s0_image() -> None:
    """SpatialData V2 multiscale containers should resolve from s0."""
    image = xr.DataArray(np.ones((1, 4, 5)), dims=("c", "y", "x"))
    multiscale = {"s0": SimpleNamespace(ds={"image": image})}

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


def test_image_to_cyx_preserves_channel_labels() -> None:
    """Channel-last image data should become channel-first without label loss."""
    arr = xr.DataArray(
        np.ones((3, 4, 2), dtype=np.uint16),
        dims=("y", "x", "c"),
        coords={"c": ["DAPI", "PolyT"]},
    )

    cyx = image_to_cyx(arr)

    assert tuple(cyx.dims) == ("c", "y", "x")
    assert list(cyx.coords["c"].values) == ["DAPI", "PolyT"]


def test_max_project_image_elements_preserves_channels() -> None:
    """Max projection should keep channel labels and max across planes."""
    z0 = xr.DataArray(
        np.ones((2, 3, 4), dtype=np.uint16),
        dims=("c", "y", "x"),
        coords={"c": ["DAPI", "PolyT"]},
    )
    z1 = xr.DataArray(
        np.full((2, 3, 4), 5, dtype=np.uint16),
        dims=("c", "y", "x"),
        coords={"c": ["DAPI", "PolyT"]},
    )

    projection = max_project_image_elements([z0, z1])

    assert tuple(projection.dims) == ("c", "y", "x")
    assert list(projection.coords["c"].values) == ["DAPI", "PolyT"]
    np.testing.assert_array_equal(
        np.asarray(projection),
        np.full((2, 3, 4), 5, dtype=np.uint16),
    )


def test_build_merscope_z_projection_prefers_existing_projection() -> None:
    """Projection-only sources should be accepted directly."""
    existing = xr.DataArray(
        np.ones((2, 3, 4), dtype=np.uint16),
        dims=("c", "y", "x"),
        coords={"c": ["DAPI", "PolyT"]},
    )

    projection = build_merscope_z_projection(
        {
            MERSCOPE_ZPROJ_IMAGE_NAME: existing,
            "run_z0": existing + 10,
        }
    )

    np.testing.assert_array_equal(np.asarray(projection), np.asarray(existing))
