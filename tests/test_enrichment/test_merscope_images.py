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


def test_copy_merscope_images_keeps_only_projection() -> None:
    """Enrichment should add only the projected MERSCOPE image."""
    src = SimpleNamespace(images={"run_z0": _make_plane(1), "run_z1": _make_plane(7)})
    dst = SimpleNamespace(images={})

    copied = _copy_merscope_images(dst, src, force=False)

    assert copied == 1
    assert set(dst.images) == {MERSCOPE_ZPROJ_IMAGE_NAME}
    projection = dst.images[MERSCOPE_ZPROJ_IMAGE_NAME]
    assert tuple(projection.dims) == ("c", "y", "x")
    assert list(projection.coords["c"].values) == ["DAPI", "PolyT"]
    np.testing.assert_array_equal(
        np.asarray(projection),
        np.full((2, 4, 4), 7, dtype=np.uint16),
    )


def test_copy_merscope_images_reuses_existing_projection() -> None:
    """Projection-only sources should be copied without adding plane layers."""
    projection = _make_plane(3).transpose("c", "y", "x")
    src = SimpleNamespace(images={MERSCOPE_ZPROJ_IMAGE_NAME: projection})
    dst = SimpleNamespace(images={})

    copied = _copy_merscope_images(dst, src, force=False)

    assert copied == 1
    assert set(dst.images) == {MERSCOPE_ZPROJ_IMAGE_NAME}
    np.testing.assert_array_equal(
        np.asarray(dst.images[MERSCOPE_ZPROJ_IMAGE_NAME]),
        np.asarray(projection),
    )
