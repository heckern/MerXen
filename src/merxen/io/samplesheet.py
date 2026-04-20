"""CSV samplesheet parsing and validation for multi-pair pipeline runs."""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class SamplePair:
    """A paired MERSCOPE + Xenium dataset for processing.

    Attributes:
        pair_id: Unique identifier for this sample pair.
        merscope_dir: Path to the raw MERSCOPE input folder.
        merscope_spatialdata_path: Optional path to a reusable built MERSCOPE
            SpatialData zarr.
        merscope_image_prefix: Prefix for matching z-plane image keys.
        merscope_z_range: Tuple of (z_start, z_end) inclusive.
        merscope_transform_path: Optional override path to the
            micron-to-mosaic transform CSV.
        merscope_channels: Channel names for MERSCOPE (e.g. ['DAPI', 'PolyT']).
        xenium_dir: Path to the raw Xenium output directory.
        xenium_spatialdata_path: Optional path to a reusable built Xenium
            SpatialData zarr.
        xenium_channels: Channel names for Xenium (e.g. ['DAPI', '18S']).
        xenium_min_qv: Minimum quality value for Xenium transcript filtering.
        merscope_voxel_layers: ProSeg voxel layers for MERSCOPE.
        xenium_voxel_layers: ProSeg voxel layers for Xenium.
    """

    pair_id: str
    merscope_z_range: tuple[int, int] = (0, 6)
    merscope_dir: Path | None = None
    merscope_spatialdata_path: Path | None = None
    merscope_image_prefix: str = ""
    merscope_transform_path: Path | None = None
    merscope_channels: list[str] = field(default_factory=lambda: ["DAPI", "PolyT"])
    xenium_dir: Path | None = None
    xenium_spatialdata_path: Path | None = None
    xenium_channels: list[str] = field(default_factory=lambda: ["DAPI", "18S"])
    xenium_min_qv: float = 20.0
    merscope_voxel_layers: int = 7
    xenium_voxel_layers: int = 2
    xenium_spec_path: Path | None = None


def parse_samplesheet(csv_path: Path) -> list[SamplePair]:
    """Parse a CSV samplesheet into a list of SamplePair objects.

    Expected columns:
        pair_id, merscope_dir, merscope_spatialdata_path, merscope_image_prefix,
        merscope_z_range, merscope_transform_path, merscope_channels,
        xenium_dir, xenium_spatialdata_path, xenium_channels, xenium_min_qv,
        merscope_voxel_layers, xenium_voxel_layers, xenium_spec_path

    Backward-compatible aliases:
        merscope_zarr_path -> merscope_spatialdata_path

    Args:
        csv_path: Path to the samplesheet CSV.

    Returns:
        List of validated SamplePair dataclass instances.

    Raises:
        FileNotFoundError: If csv_path does not exist.
        ValueError: If required columns are missing.
    """
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Samplesheet not found: {csv_path}")

    pairs = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        required_cols = {"pair_id"}
        if reader.fieldnames is None:
            raise ValueError(f"Empty samplesheet: {csv_path}")
        missing = required_cols - set(reader.fieldnames)
        if missing:
            raise ValueError(
                f"Samplesheet missing required columns: {missing}. "
                f"Found: {reader.fieldnames}"
            )

        for row_num, row in enumerate(reader, start=2):
            z_range_str = row.get("merscope_z_range", "0-6")
            z_parts = z_range_str.split("-")
            z_range = (int(z_parts[0]), int(z_parts[1]))
            merscope_spatialdata_raw = row.get("merscope_spatialdata_path") or row.get(
                "merscope_zarr_path"
            )

            pair = SamplePair(
                pair_id=row["pair_id"],
                merscope_dir=_optional_path(row.get("merscope_dir")),
                merscope_spatialdata_path=_optional_path(merscope_spatialdata_raw),
                merscope_image_prefix=row.get("merscope_image_prefix", ""),
                merscope_z_range=z_range,
                merscope_transform_path=_optional_path(
                    row.get("merscope_transform_path")
                ),
                merscope_channels=_parse_list(
                    row.get("merscope_channels", "DAPI,PolyT")
                ),
                xenium_dir=_optional_path(row.get("xenium_dir")),
                xenium_spatialdata_path=_optional_path(
                    row.get("xenium_spatialdata_path")
                ),
                xenium_channels=_parse_list(row.get("xenium_channels", "DAPI,18S")),
                xenium_min_qv=float(row.get("xenium_min_qv", "20.0")),
                merscope_voxel_layers=int(row.get("merscope_voxel_layers", "7")),
                xenium_voxel_layers=int(row.get("xenium_voxel_layers", "2")),
                xenium_spec_path=_optional_path(row.get("xenium_spec_path")),
            )
            pairs.append(pair)
            logger.info("Parsed sample pair %d: %s", row_num - 1, pair.pair_id)

    return pairs


def validate_samplesheet(pairs: list[SamplePair]) -> None:
    """Validate that all paths in a samplesheet exist.

    Args:
        pairs: List of SamplePair instances to validate.

    Raises:
        FileNotFoundError: If any required path does not exist.
    """
    errors = []
    for pair in pairs:
        if pair.merscope_dir is None and pair.merscope_spatialdata_path is None:
            errors.append(
                f"[{pair.pair_id}] Provide either merscope_dir or "
                "merscope_spatialdata_path."
            )
        if pair.xenium_dir is None and pair.xenium_spatialdata_path is None:
            errors.append(
                f"[{pair.pair_id}] Provide either xenium_dir or "
                "xenium_spatialdata_path."
            )
        if (
            pair.merscope_dir is not None
            and not pair.merscope_dir.exists()
            and pair.merscope_spatialdata_path is None
        ):
            errors.append(
                f"[{pair.pair_id}] MERSCOPE dir not found: {pair.merscope_dir}"
            )
        if (
            pair.merscope_spatialdata_path is not None
            and not pair.merscope_spatialdata_path.exists()
            and pair.merscope_dir is None
        ):
            errors.append(
                f"[{pair.pair_id}] MERSCOPE SpatialData zarr not found: "
                f"{pair.merscope_spatialdata_path}"
            )
        if (
            pair.xenium_dir is not None
            and not pair.xenium_dir.exists()
            and pair.xenium_spatialdata_path is None
        ):
            errors.append(f"[{pair.pair_id}] Xenium dir not found: {pair.xenium_dir}")
        if (
            pair.xenium_spatialdata_path is not None
            and not pair.xenium_spatialdata_path.exists()
            and pair.xenium_dir is None
        ):
            errors.append(
                f"[{pair.pair_id}] Xenium SpatialData zarr not found: "
                f"{pair.xenium_spatialdata_path}"
            )
        if (
            pair.merscope_transform_path is not None
            and not pair.merscope_transform_path.exists()
        ):
            errors.append(
                f"[{pair.pair_id}] MERSCOPE transform not found: "
                f"{pair.merscope_transform_path}"
            )
        if pair.xenium_spec_path is not None and not pair.xenium_spec_path.exists():
            errors.append(
                f"[{pair.pair_id}] Xenium spec not found: {pair.xenium_spec_path}"
            )
    if errors:
        raise FileNotFoundError("Samplesheet validation failed:\n" + "\n".join(errors))


def _parse_list(value: str) -> list[str]:
    """Parse a comma-separated string into a list of stripped strings."""
    return [v.strip() for v in value.split(",") if v.strip()]


def _optional_path(value: str | None) -> Path | None:
    """Convert a possibly-empty string field into an optional Path."""
    if value is None:
        return None
    cleaned = value.strip()
    if not cleaned:
        return None
    return Path(cleaned)
