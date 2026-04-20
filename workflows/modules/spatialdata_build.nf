process BUILD_SPATIALDATA {
    tag "${pair_id}:${platform}"

    publishDir { "${params.outdir}/${pair_id}/${platform.toLowerCase()}/spatialdata" }, mode: "copy", overwrite: true

    input:
    tuple val(key), val(pair_id), val(platform), val(build_config_json)

    output:
    tuple val(key), val(pair_id), val(platform), path("spatialdata_out/source_spatialdata.zarr")

    script:
    def forceFlag = params.force_spatialdata_build ? "--force-rerun" : ""
    """
    set -euo pipefail

    cat > build_config.json <<'JSON'
${build_config_json}
JSON

    merxen build-spatialdata --config build_config.json ${forceFlag}
    """
}
