process ENRICH {
    tag "${pair_id}:${platform}"

    publishDir { "${params.outdir}/${pair_id}/${platform.toLowerCase()}/enrichment" }, mode: "copy", overwrite: true

    input:
    tuple val(key), val(pair_id), val(platform), val(enrich_config_json), path(latest_zarr), path(mask_path)

    output:
    tuple val(key), val(pair_id), val(platform), path("latest_input.zarr"), path("enrich_out")

    script:
    """
    set -euo pipefail

    cp -r ${latest_zarr} latest_input.zarr
    cp ${mask_path} cellpose_masks_tiled.npy

    cat > enrich_config.json <<'JSON'
${enrich_config_json}
JSON

    merxen enrich --config enrich_config.json
    """
}
