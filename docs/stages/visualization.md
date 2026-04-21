# Stage 6 — Visualization

Produces a fixed set of PNG plots for a pair: gene-level scatter plots,
per-dataset QC plots, transcript density overviews, and a sanity-check image
overlay. Runs once per pair, after comparison has completed.

## What it does

1. Rerun the gene comparison internally (re-opens the enriched zarrs) and
   plot log-log scatter plots of MERSCOPE vs Xenium normalized counts.
2. Recompute per-dataset QC metrics and plot geometry histograms and cell
   metric violins for both platforms.
3. For each dataset, plot a transcript density heatmap using the points
   table.
4. For each dataset, crop a central 1024×1024 region from the first image
   and overlay the first shape layer for a visual sanity check.
5. Plot an assignment-rate bar chart comparing the percentage of transcripts
   assigned across platforms.

## Nextflow process

[`VISUALIZE`](../../workflows/modules/visualization.nf) — one instance per
`pair_id`, downstream of `COMPARE`.

- **Input:** `tuple(pair_id, merscope_zarr, xenium_zarr)`.
- **CLI:** `merxen visualize --config visualize_config.json`.
- **Output:** `tuple(pair_id, visualize_out/)`.
- **publishDir:** `${outdir}/${pair_id}/visualization/` (copy mode).

## Python entry points

| Function | File |
|----------|------|
| CLI `visualize_command` | [cli/run_visualization.py:34](../../src/merxen/cli/run_visualization.py#L34) |
| `plot_gene_scatter` | [visualization/gene_scatter.py:14](../../src/merxen/visualization/gene_scatter.py#L14) |
| `plot_geometry_histograms` | [visualization/qc_plots.py:13](../../src/merxen/visualization/qc_plots.py#L13) |
| `plot_cell_metrics_violin` | [visualization/qc_plots.py:44](../../src/merxen/visualization/qc_plots.py#L44) |
| `plot_assignment_bar` | [visualization/qc_plots.py:70](../../src/merxen/visualization/qc_plots.py#L70) |
| `plot_density_overview` | [visualization/density_overview.py:29](../../src/merxen/visualization/density_overview.py#L29) |
| `plot_sanity_overlay` | [visualization/sanity_plots.py:13](../../src/merxen/visualization/sanity_plots.py#L13) |

## Config schema

`VisualizationConfig` — [config.py:186](../../src/merxen/config.py#L186).

| Field | Description |
|-------|-------------|
| `merscope_zarr_path` | Enriched MERSCOPE zarr. |
| `xenium_zarr_path` | Enriched Xenium zarr. |
| `output_dir` | Where `visualize_out/` is populated. |
| `pair_id` | Prefix for output filenames. |

## Outputs

Written under `visualize_out/`:

| Kind | File | Contents |
|------|------|----------|
| Gene scatter | `<pair_id>_gene_scatter_total_normalized.png` | MERSCOPE vs Xenium, all transcripts (normalized). |
| Gene scatter | `<pair_id>_gene_scatter_assigned_normalized.png` | MERSCOPE vs Xenium, transcripts assigned to cells. |
| Geometry | `<pair_id>_<platform>_geometry_hist.png` | Histograms of area, eccentricity, etc. |
| Cell metrics | `<pair_id>_<platform>_cell_violin.png` | Violin plots of transcripts/cell, genes/cell. |
| Density | `<pair_id>_<platform>_density_overview.png` | 2D histogram of transcript locations. |
| Sanity overlay | `<pair_id>_<platform>_sanity_overlay.png` | Central 1024×1024 image crop with polygons drawn on top. |
| Assignment rate | `<pair_id>_assignment_rate_bar.png` | Bar chart of `pct_assigned` per platform. |

## Notes

- The visualization stage does **not** read the CSVs produced by the
  comparison stage; it recomputes them. This keeps stages independent but
  means large zarrs are opened twice per run.
- The sanity overlay uses whichever image is keyed first in
  `sdata.images` and whichever shape layer is keyed first in
  `sdata.shapes`. For consistency it is worth running enrichment first so
  that both layers exist.
- Points coordinate columns are resolved with `first_existing_col` across
  `x`, `x_micron`, `x_location`, `global_x`, `x_global_px`, `observed_x`
  (and the corresponding `y_*`). If none are present the density plot is
  skipped silently.
