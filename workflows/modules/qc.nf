process QC {
    tag "${pair_id}:${platform}"

    publishDir { "${params.outdir}/${pair_id}/${platform.toLowerCase()}/qc" }, mode: "symlink", overwrite: true

    input:
    tuple val(key), val(pair_id), val(platform), path(latest_zarr)

    output:
    tuple val(key), val(pair_id), val(platform), path(latest_zarr), path("qc_out")

    script:
    """
    set -euo pipefail

    cat > qc_config.json <<JSON
{
  "dataset_name": "${pair_id}_${platform}",
  "latest_zarr_path": "${latest_zarr}",
  "output_dir": "qc_out"
}
JSON

    merxen qc --config qc_config.json
    """
}
