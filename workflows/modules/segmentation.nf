process SEGMENT {
    tag "${pair_id}:${platform}"

    publishDir { "${params.outdir}/${pair_id}/${platform.toLowerCase()}/segmentation" }, mode: "symlink", overwrite: true

    input:
    tuple val(key), val(pair_id), val(platform), val(seg_config_json)

    output:
    tuple val(key), val(pair_id), val(platform), path("segment_out/proseg_base_latest.zarr"), path("segment_out/cellpose_masks_tiled.npy"), path("segment_out/transcripts_for_proseg.csv")

    script:
    """
    set -euo pipefail

    cat > segment_config.json <<'JSON'
${seg_config_json}
JSON

    merxen segment --config segment_config.json
    """
}
