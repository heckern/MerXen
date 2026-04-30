process ALIGN {
    tag "${pair_id}"

    publishDir { "${params.outdir}/${pair_id}/alignment" }, mode: "copy", overwrite: true

    input:
    tuple val(pair_id),
        path(merscope_zarr, stageAs: "merscope_latest_input.zarr"),
        path(xenium_zarr, stageAs: "xenium_latest_input.zarr")

    output:
    tuple val(pair_id),
        path("align_out/merscope_aligned.zarr"),
        path("align_out/xenium_aligned.zarr"),
        path("align_out/alignment_transform.json"),
        path("align_out/alignment_coords")

    script:
    """
    set -euo pipefail

    cat > align_config.json <<JSON
{
  "pair_id": "${pair_id}",
  "merscope_zarr_path": "${merscope_zarr}",
  "xenium_zarr_path": "${xenium_zarr}",
  "output_dir": "align_out",
  "spateo": {
    "mode": "${params.alignment_spateo_mode}",
    "device": "${params.alignment_device}",
    "dtype": "${params.alignment_dtype}",
    "selected_mode": "${params.alignment_selected_mode}",
    "max_iter": ${params.alignment_max_iter},
    "beta": ${params.alignment_beta},
    "lambda_vf": ${params.alignment_lambda_vf},
    "k": ${params.alignment_k},
    "partial_robust_level": ${params.alignment_partial_robust_level},
    "SVI_mode": ${params.alignment_svi_mode},
    "n_sampling": ${params.alignment_n_sampling},
    "sparse_calculation_mode": ${params.alignment_sparse_calculation_mode},
    "use_chunk": ${params.alignment_use_chunk},
    "chunk_capacity": ${params.alignment_chunk_capacity},
    "use_hvg": ${params.alignment_use_hvg},
    "n_top_genes": ${params.alignment_n_top_genes},
    "rbf_neighbors": ${params.alignment_rbf_neighbors},
    "rbf_smoothing": ${params.alignment_rbf_smoothing},
    "max_nonrigid_anchors": ${params.alignment_max_nonrigid_anchors}
  }
}
JSON

    export PYTORCH_CUDA_ALLOC_CONF="${params.alignment_pytorch_cuda_alloc_conf}"
    merxen align --config align_config.json
    """
}

process ALIGN_QC {
    tag "${pair_id}"

    publishDir { "${params.outdir}/${pair_id}/alignment_qc" }, mode: "copy", overwrite: true

    input:
    tuple val(pair_id),
        path(merscope_zarr, stageAs: "merscope_aligned_input.zarr"),
        path(xenium_zarr, stageAs: "xenium_aligned_input.zarr"),
        path(transform_json),
        path(coords_dir)

    output:
    tuple val(pair_id), path("alignment_qc_out")

    script:
    """
    set -euo pipefail

    cat > alignment_qc_config.json <<JSON
{
  "pair_id": "${pair_id}",
  "merscope_zarr_path": "${merscope_zarr}",
  "xenium_zarr_path": "${xenium_zarr}",
  "transform_json_path": "${transform_json}",
  "output_dir": "alignment_qc_out",
  "grid_rows": ${params.alignment_qc_grid_rows},
  "grid_cols": ${params.alignment_qc_grid_cols}
}
JSON

    merxen alignment-qc --config alignment_qc_config.json
    """
}
