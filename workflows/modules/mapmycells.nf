import groovy.json.JsonOutput

process MAPMYCELLS {
    tag "${pair_id}"

    publishDir { "${params.outdir}/${pair_id}/mapmycells" }, mode: "copy", overwrite: true

    input:
    tuple val(pair_id),
        path(clustering_out_dir, stageAs: "clustering_squidpy_input")

    output:
    tuple val(pair_id), path("mapmycells_out")

    script:
    def dropLevelJson = params.mapmycells_drop_level == null ? "null" : JsonOutput.toJson(params.mapmycells_drop_level.toString())
    def queryLayerJson = params.mapmycells_query_layer == null ? "null" : JsonOutput.toJson(params.mapmycells_query_layer.toString())
    def geneIdColumnJson = params.mapmycells_gene_id_column == null ? "null" : JsonOutput.toJson(params.mapmycells_gene_id_column.toString())
    def obsIdColumnJson = params.mapmycells_obs_id_column == null ? "null" : JsonOutput.toJson(params.mapmycells_obs_id_column.toString())
    def tmpDirJson = params.mapmycells_tmp_dir == null ? "null" : JsonOutput.toJson(params.mapmycells_tmp_dir.toString())
    """
    set -euo pipefail

    cat > mapmycells_config.json <<JSON
{
  "pair_id": "${pair_id}",
  "output_dir": "mapmycells_out",
  "samples": [
    {
      "sample_id": "${pair_id}_MERSCOPE",
      "platform": "MERSCOPE",
      "anndata_path": "${clustering_out_dir}/merscope/${pair_id}_MERSCOPE_clustered.h5ad",
      "query_layer": ${queryLayerJson},
      "gene_id_column": ${geneIdColumnJson},
      "obs_id_column": ${obsIdColumnJson}
    },
    {
      "sample_id": "${pair_id}_XENIUM",
      "platform": "XENIUM",
      "anndata_path": "${clustering_out_dir}/xenium/${pair_id}_XENIUM_clustered.h5ad",
      "query_layer": ${queryLayerJson},
      "gene_id_column": ${geneIdColumnJson},
      "obs_id_column": ${obsIdColumnJson}
    }
  ],
  "marker_lookup_path": "${params.mapmycells_marker_lookup_path}",
  "precomputed_stats_path": "${params.mapmycells_precomputed_stats_path}",
  "drop_level": ${dropLevelJson},
  "normalization": "${params.mapmycells_normalization}",
  "bootstrap_factor": ${params.mapmycells_bootstrap_factor},
  "bootstrap_iteration": ${params.mapmycells_bootstrap_iteration},
  "n_processors": ${params.mapmycells_n_processors},
  "chunk_size": ${params.mapmycells_chunk_size},
  "rng_seed": ${params.mapmycells_rng_seed},
  "max_gb": ${params.mapmycells_max_gb},
  "tmp_dir": ${tmpDirJson},
  "cloud_safe": ${params.mapmycells_cloud_safe},
  "flatten": ${params.mapmycells_flatten},
  "verbose_csv": ${params.mapmycells_verbose_csv}
}
JSON

    merxen mapmycells --config mapmycells_config.json
    """
}
