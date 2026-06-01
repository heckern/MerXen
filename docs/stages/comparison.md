# Stage 5 — Comparison

Cross-platform gene-level comparison. Runs once per pair in
`--analysis_mode paired`, after both MERSCOPE and Xenium have passed QC, and
quantifies how well the two platforms agree on the same 300-gene panel.
Single-platform runs skip this stage.

## What it does

For each dataset:

- **Total counts per gene** — from the transcript points, regardless of
  whether a transcript was assigned to a cell.
- **Assigned counts per gene** — from the primary counts table (only
  transcripts that landed in a segmented cell).

Then:

- Filter each platform's gene set to retain only panel genes (drops Xenium
  control codes and MERSCOPE blanks).
- CP10K-style normalization: scale to sum = 10 000.
- Align genes across platforms and emit four DataFrames
  (total, assigned, total-normalized, assigned-normalized).
- Fit a linear regression in log10 space on the normalized counts to
  produce slope / intercept / R².

## Nextflow process

[`COMPARE`](../../workflows/modules/comparison.nf) — one instance per
`pair_id`. Takes the enriched zarrs for both halves of the pair as input.

- **Input:** `tuple(pair_id, merscope_zarr, xenium_zarr)`.
- **CLI:** `merxen compare --config compare_config.json`.
- **Output:** `tuple(pair_id, compare_out/)`.
- **publishDir:** `${outdir}/${pair_id}/comparison/` (copy mode).

## Python entry points

| Function | File |
|----------|------|
| CLI `compare_command` | [cli/run_comparison.py:37](../../src/merxen/cli/run_comparison.py#L37) |
| `compute_gene_comparison_from_paths` | [qc/gene_comparison.py:187](../../src/merxen/qc/gene_comparison.py#L187) |
| `compute_gene_comparison` | [qc/gene_comparison.py:123](../../src/merxen/qc/gene_comparison.py#L123) |
| `gene_totals_from_points` | [qc/gene_comparison.py:29](../../src/merxen/qc/gene_comparison.py#L29) |
| `gene_totals_from_table` | [qc/gene_comparison.py:18](../../src/merxen/qc/gene_comparison.py#L18) |
| `apply_dataset_filter` | [qc/gene_comparison.py:76](../../src/merxen/qc/gene_comparison.py#L76) |
| `normalize_counts` | [qc/gene_comparison.py:92](../../src/merxen/qc/gene_comparison.py#L92) |
| `fit_linear` | [qc/gene_comparison.py:112](../../src/merxen/qc/gene_comparison.py#L112) |

## Config schema

`ComparisonConfig` — [config.py:177](../../src/merxen/config.py#L177).

| Field | Description |
|-------|-------------|
| `merscope_zarr_path` | Enriched MERSCOPE zarr. |
| `xenium_zarr_path` | Enriched Xenium zarr. |
| `output_dir` | Where `compare_out/` is populated. |
| `pair_id` | Prefix for output filenames. |

## Walkthrough

1. `compute_gene_comparison_from_paths` opens both enriched zarrs.
2. Pull four count vectors:
   - `x_total_all` / `m_total_all` — `gene_totals_from_points`, includes
     unassigned transcripts.
   - `x_assigned_all` / `m_assigned_all` — `gene_totals_from_table`, sums
     the `.X` matrix in each platform's primary table.
3. `apply_dataset_filter` strips platform-specific non-panel genes:
   - Xenium: drops codeword / control / unassigned patterns.
   - MERSCOPE: drops `Blank-*` probes.
4. `normalize_counts` rescales each vector so its entries sum to 10 000,
   returning the normalized vector and the original total.
5. `compare_df` joins the two platforms on gene name, keeping only genes
   present in both. Produces one DataFrame per count type.
6. `fit_linear` fits `log10(merscope) = slope × log10(xenium) + intercept`
   on the normalized total and normalized assigned counts, reporting R².

## Outputs

Written under `compare_out/`:

| File | Contents |
|------|----------|
| `<pair_id>_total_counts_compare.csv` | Gene × platform total counts. |
| `<pair_id>_assigned_counts_compare.csv` | Gene × platform counts from the table. |
| `<pair_id>_total_normalized_compare.csv` | CP10K-normalized total counts. |
| `<pair_id>_assigned_normalized_compare.csv` | CP10K-normalized assigned counts. |
| `<pair_id>_comparison_metrics.json` | Platform totals plus linear-fit metrics. |

The metrics JSON has the form:

```json
{
  "totals": {
    "xenium_total_sum": 1234567,
    "merscope_total_sum": 2345678,
    "xenium_assigned_sum": 987654,
    "merscope_assigned_sum": 1876543
  },
  "fits": {
    "total_log10":    {"slope": 0.98, "intercept": -0.03, "r2": 0.92},
    "assigned_log10": {"slope": 1.01, "intercept":  0.02, "r2": 0.95}
  },
  "n_genes": { "total_counts": 300, "assigned_counts": 300, ... }
}
```

## Interpreting the fits

- **R² near 1** — platforms agree on relative gene abundance.
- **Slope deviates from 1** — one platform compresses or stretches the
  dynamic range. Typical values fall between 0.8 and 1.2.
- **Assigned R² >> total R²** — platforms differ on transcript localization
  but agree on total panel loading. Usually a segmentation-quality story.
- **Total R² >> assigned R²** — segmentation is dropping different genes
  on each platform. Revisit ProSeg parameters.
