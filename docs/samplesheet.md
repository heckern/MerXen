# Samplesheet format

The samplesheet is a CSV with one row per biological sample or adjacent-section
pair. In the default `--analysis_mode paired`, each row must contain one
MERSCOPE and one Xenium dataset. In `--analysis_mode merscope` or
`--analysis_mode xenium`, only the selected platform's source/cache columns are
required. A template lives at
[workflows/samplesheet.example.csv](../workflows/samplesheet.example.csv).

## Columns

| Column | Required | Description |
|--------|----------|-------------|
| `pair_id` | **yes** | Unique identifier for this row. Used as the top-level output directory name. |
| `merscope_dir` | required for MERSCOPE modes if no cache | Path to the raw MERSCOPE region export folder (contains `transcripts.parquet`, `cell_boundaries/`, `images/`, etc.). |
| `merscope_spatialdata_path` | required for MERSCOPE modes if no raw dir | Path to an existing (or desired) reusable MERSCOPE SpatialData zarr. If it exists, the build step is **skipped** unless `--force_spatialdata_build true` is passed to Nextflow. |
| `merscope_image_prefix` | no | Prefix used to match z-plane image keys when more than one run is present. |
| `merscope_z_range` | no | Inclusive z-layer range as `start-end`. Defaults to `0-6`. |
| `merscope_transform_path` | no | Override path to the `micron_to_mosaic_pixel_transform.csv`. If not set, MerXen looks inside `merscope_dir`. |
| `merscope_channels` | no | Comma-separated channel names for Cellpose. Defaults to `DAPI,PolyT`. |
| `merscope_voxel_layers` | no | ProSeg voxel layer count for MERSCOPE. Defaults to `7` (from `nextflow.config`). |
| `xenium_dir` | required for Xenium modes if no cache | Path to the raw Xenium export folder. |
| `xenium_spatialdata_path` | required for Xenium modes if no raw dir | Path to an existing (or desired) reusable Xenium SpatialData zarr. Cached the same way as the MERSCOPE path. |
| `xenium_channels` | no | Comma-separated channel names for Cellpose. Defaults to `DAPI,18S`. |
| `xenium_min_qv` | no | Minimum transcript quality value to retain. Defaults to `20`. |
| `xenium_voxel_layers` | no | ProSeg voxel layer count for Xenium. Defaults to `2`. |
| `xenium_spec_path` | no | Override path to `experiment.xenium` or `specs.json` used to derive the micron→pixel transform. |

### Aliases

`merscope_zarr_path` is accepted as an alias for `merscope_spatialdata_path`
for backwards compatibility.

## Validation rules

From [workflows/main.nf](../workflows/main.nf):

- Every row must have a non-empty `pair_id`.
- In `--analysis_mode paired`, every row must provide **either**
  `merscope_dir` **or** `merscope_spatialdata_path`, and **either**
  `xenium_dir` **or** `xenium_spatialdata_path`.
- In `--analysis_mode merscope`, only MERSCOPE source/cache columns are
  required.
- In `--analysis_mode xenium`, only Xenium source/cache columns are required.

The Python-side parser lives in
[src/merxen/io/samplesheet.py:52](../src/merxen/io/samplesheet.py#L52). It is
used by unit tests and can be invoked directly by scripts, but the Nextflow
workflow parses the CSV itself.

## Minimal example

Paired mode:

```csv
pair_id,merscope_dir,xenium_dir
EXAMPLE01,/path/to/merscope/EXAMPLE01/region_R1,/path/to/xenium/EXAMPLE01
```

This is enough to run the whole pipeline end-to-end with default parameters.
All other columns fall back to defaults.

Xenium-only mode:

```csv
pair_id,xenium_dir,xenium_spatialdata_path,xenium_channels,xenium_min_qv
EXAMPLE01,/path/to/xenium/EXAMPLE01,/path/to/cache/EXAMPLE01_xenium.zarr,"DAPI,18S",20
```

Run it with `--analysis_mode xenium`. MERSCOPE-only samplesheets are analogous:
provide `merscope_dir` or `merscope_spatialdata_path` and run with
`--analysis_mode merscope`.

## Full example

```csv
pair_id,merscope_dir,merscope_spatialdata_path,merscope_image_prefix,merscope_z_range,merscope_transform_path,merscope_channels,xenium_dir,xenium_spatialdata_path,xenium_channels,xenium_min_qv,merscope_voxel_layers,xenium_voxel_layers,xenium_spec_path
EXAMPLE01,/path/to/merscope/EXAMPLE01,/path/to/cache/EXAMPLE01_merscope.zarr,,0-6,,"DAPI,PolyT",/path/to/xenium/EXAMPLE01,/path/to/cache/EXAMPLE01_xenium.zarr,"DAPI,18S",20,7,2,
EXAMPLE02,/path/to/merscope/EXAMPLE02,,,,,"DAPI,PolyT",/path/to/xenium/EXAMPLE02,,"DAPI,18S",20,7,2,
```

Row 1 uses cached SpatialData zarrs if they already exist. Row 2 builds
everything fresh from the raw exports.

## Tips

- Quote any list-valued field that contains commas (`merscope_channels`,
  `xenium_channels`).
- Leave a field empty (`,,`) to fall back to the Nextflow default.
- Prefer absolute paths. Relative paths are resolved from the shell that ran
  `nextflow`, not from the work directory.
- Keep one samplesheet per batch. Nextflow runs all rows in parallel up to the
  executor limit.
