"""Tests for MERSCOPE enrichment image handling."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import xarray as xr

from merxen.enrichment.enrich import MERSCOPE_ZPROJ_IMAGE_NAME, _copy_merscope_images


def _make_plane(fill_value: int) -> xr.DataArray:
    """Build a tiny plane image with MERSCOPE-style channel-last layout."""
    return xr.DataArray(
        np.full((4, 4, 2), fill_value, dtype=np.uint16),
        dims=("y", "x", "c"),
        coords={"c": ["DAPI", "PolyT"]},
    )


def test_copy_merscope_images_keeps_planes_and_projection() -> None:
    """Enrichment should keep the raw planes and add a projection."""
    src = SimpleNamespace(images={"run_z0": _make_plane(1), "run_z1": _make_plane(7)})
    dst = SimpleNamespace(images={})

    copied = _copy_merscope_images(dst, src, force=False)

    assert copied == 3
    assert set(dst.images) == {"run_z0", "run_z1", MERSCOPE_ZPROJ_IMAGE_NAME}
    np.testing.assert_array_equal(
        np.asarray(dst.images["run_z0"]),
        np.asarray(src.images["run_z0"]),
    )
    np.testing.assert_array_equal(
        np.asarray(dst.images["run_z1"]),
        np.asarray(src.images["run_z1"]),
    )
    projection = dst.images[MERSCOPE_ZPROJ_IMAGE_NAME]
    assert tuple(projection.dims) == ("c", "y", "x")
    np.testing.assert_array_equal(
        np.asarray(projection),
        np.full((2, 4, 4), 7, dtype=np.uint16),
    )
