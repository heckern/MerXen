process CLUSTERING_SQUIDPY {
    tag "${pair_id}"

    publishDir { "${params.outdir}/${pair_id}/clustering_squidpy" }, mode: "copy", overwrite: true

    input:
    tuple val(pair_id),
        path(merscope_zarr, stageAs: "merscope_latest_input.zarr"),
        val(xenium_zarr)

    output:
    tuple val(pair_id), path("clustering_squidpy_out")

    script:
    """
    set -euo pipefail

    cat > clustering_squidpy_config.json <<JSON
{
  "pair_id": "${pair_id}",
  "output_dir": "clustering_squidpy_out",
  "samples": [
    {
      "sample_id": "${pair_id}_MERSCOPE",
      "platform": "MERSCOPE",
      "zarr_path": "${merscope_zarr}"
    },
    {
      "sample_id": "${pair_id}_XENIUM",
      "platform": "XENIUM",
      "zarr_path": "${xenium_zarr}"
    }
  ],
  "min_counts": ${params.clustering_squidpy_min_counts},
  "min_cells": ${params.clustering_squidpy_min_cells},
  "normalize_target_sum": ${params.clustering_squidpy_normalize_target_sum},
  "n_pcs": ${params.clustering_squidpy_n_pcs},
  "n_neighbors": ${params.clustering_squidpy_n_neighbors},
  "leiden_resolution": ${params.clustering_squidpy_leiden_resolution},
  "random_seed": ${params.clustering_squidpy_random_seed},
  "spatial_point_size": ${params.clustering_squidpy_spatial_point_size},
  "figure_dpi": ${params.clustering_squidpy_figure_dpi}
}
JSON

    merxen clustering-squidpy --config clustering_squidpy_config.json
    """
}
