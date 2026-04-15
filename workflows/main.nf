nextflow.enable.dsl = 2

import groovy.json.JsonOutput

include { SEGMENT } from "./modules/segmentation"
include { ENRICH } from "./modules/enrichment"
include { QC } from "./modules/qc"
include { COMPARE } from "./modules/comparison"
include { VISUALIZE } from "./modules/visualization"

def parseChannels(rawValue, defaults) {
    if (rawValue == null) {
        return defaults
    }
    def values = rawValue
        .toString()
        .split(",")
        .collect { it.trim() }
        .findAll { it.length() > 0 }
    return values ? values : defaults
}

def parseRange(rawValue, fallbackStart = 0, fallbackEnd = 6) {
    if (rawValue == null || rawValue.toString().trim().isEmpty()) {
        return [fallbackStart as int, fallbackEnd as int]
    }
    def parts = rawValue.toString().split("-").collect { it.trim() }
    if (parts.size() != 2) {
        return [fallbackStart as int, fallbackEnd as int]
    }
    return [parts[0] as int, parts[1] as int]
}

def intOrDefault(rawValue, defaultValue) {
    if (rawValue == null || rawValue.toString().trim().isEmpty()) {
        return defaultValue as int
    }
    return rawValue as int
}

def floatOrDefault(rawValue, defaultValue) {
    if (rawValue == null || rawValue.toString().trim().isEmpty()) {
        return defaultValue as float
    }
    return rawValue as float
}

workflow {
    if (!params.samplesheet) {
        error "Missing required parameter: --samplesheet"
    }
    if (!params.proseg_binary) {
        error "Missing required parameter: --proseg_binary"
    }

    samplesheet_ch = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)

    segment_inputs_ch = samplesheet_ch.flatMap { row ->
        def pairId = row.pair_id?.toString()?.trim()
        if (!pairId) {
            error "Found samplesheet row with missing pair_id: ${row}"
        }

        def merscopeVoxelLayers = intOrDefault(
            row.merscope_voxel_layers,
            params.default_merscope_voxel_layers
        )
        def xeniumVoxelLayers = intOrDefault(
            row.xenium_voxel_layers,
            params.default_xenium_voxel_layers
        )
        def xeniumMinQv = floatOrDefault(row.xenium_min_qv, params.xenium_min_qv)
        def merscopeRange = parseRange(row.merscope_z_range, 0, 6)

        def baseConfig = [
            cellpose: [
                model_type: params.cellpose_model_type,
                gpu: params.cellpose_gpu,
                diameter: params.cellpose_diameter,
                flow_threshold: params.cellpose_flow_threshold,
                cellprob_threshold: params.cellpose_cellprob,
                tile_overlap: params.cellpose_tile_overlap,
                bsize: params.cellpose_bsize,
                factor_rescale: 1.0,
            ],
            proseg: [
                binary_path: params.proseg_binary,
                samples: params.proseg_samples,
                voxel_size: params.proseg_voxel_size,
                burnin_voxel_size: params.proseg_burnin_voxel_size,
                nuclear_reassignment_prob: params.proseg_nuclear_reassignment_prob,
                diffusion_probability: params.proseg_diffusion_probability,
                cell_compactness: params.proseg_cell_compactness,
                num_threads: params.proseg_num_threads,
                voxel_layers: 2,
            ],
            memory: [
                max_system_ram_gb: params.max_ram_gb,
                memory_warn_gb: params.warn_ram_gb,
                transcript_chunk_rows: params.transcript_chunk_rows,
            ],
        ]

        def merscopeConfig = baseConfig + [
            dataset: [
                name: "${pairId}_MERSCOPE",
                platform: "MERSCOPE",
                data_path: row.merscope_zarr_path,
                channels: parseChannels(row.merscope_channels, ["DAPI", "PolyT"]),
                output_dir: "segment_out",
                image_prefix: row.merscope_image_prefix ?: null,
                z_range: merscopeRange,
                transform_path: row.merscope_transform_path ?: null,
                proseg_overrides: [voxel_layers: merscopeVoxelLayers],
            ],
        ]

        def xeniumConfig = baseConfig + [
            dataset: [
                name: "${pairId}_XENIUM",
                platform: "XENIUM",
                data_path: row.xenium_dir,
                channels: parseChannels(row.xenium_channels, ["DAPI", "18S"]),
                output_dir: "segment_out",
                xenium_spec_path: row.xenium_spec_path ?: null,
                min_qv: xeniumMinQv,
                proseg_overrides: [voxel_layers: xeniumVoxelLayers],
            ],
        ]

        def merscopeKey = "${pairId}|MERSCOPE"
        def xeniumKey = "${pairId}|XENIUM"

        [
            tuple(
                merscopeKey,
                pairId,
                "MERSCOPE",
                JsonOutput.prettyPrint(JsonOutput.toJson(merscopeConfig)),
            ),
            tuple(
                xeniumKey,
                pairId,
                "XENIUM",
                JsonOutput.prettyPrint(JsonOutput.toJson(xeniumConfig)),
            ),
        ]
    }

    metadata_ch = samplesheet_ch.flatMap { row ->
        def pairId = row.pair_id?.toString()?.trim()
        def merscopeKey = "${pairId}|MERSCOPE"
        def xeniumKey = "${pairId}|XENIUM"

        [
            tuple(
                merscopeKey,
                pairId,
                "MERSCOPE",
                row.merscope_zarr_path.toString(),
                row.merscope_transform_path?.toString() ?: "",
            ),
            tuple(
                xeniumKey,
                pairId,
                "XENIUM",
                row.xenium_dir.toString(),
                row.xenium_spec_path?.toString() ?: "",
            ),
        ]
    }

    segment_results_ch = SEGMENT(segment_inputs_ch)

    enrich_inputs_ch = segment_results_ch
        .join(metadata_ch)
        .map { key, pairId, platform, latestZarr, maskPath, transcriptsCsv, pairMeta, platformMeta, originalDataPath, transformPath ->
            if (pairId != pairMeta || platform != platformMeta) {
                error "Internal channel mismatch for key=${key}: ${pairId}/${platform} vs ${pairMeta}/${platformMeta}"
            }

            def enrichConfig = [
                dataset_name: "${pairId}_${platform}",
                platform: platform,
                latest_zarr_path: "latest_input.zarr",
                mask_path: "cellpose_masks_tiled.npy",
                original_data_path: originalDataPath,
                output_dir: "enrich_out",
                transform_path: transformPath ? transformPath : null,
            ]

            tuple(
                key,
                pairId,
                platform,
                JsonOutput.prettyPrint(JsonOutput.toJson(enrichConfig)),
                latestZarr,
                maskPath,
            )
        }

    enrich_results_ch = ENRICH(enrich_inputs_ch)

    qc_inputs_ch = enrich_results_ch.map { key, pairId, platform, enrichedLatestZarr, enrichOutDir ->
        tuple(key, pairId, platform, enrichedLatestZarr)
    }

    qc_results_ch = QC(qc_inputs_ch)

    merscope_qc_ch = qc_results_ch
        .filter { key, pairId, platform, enrichedLatestZarr, qcOutDir -> platform == "MERSCOPE" }
        .map { key, pairId, platform, enrichedLatestZarr, qcOutDir -> tuple(pairId, enrichedLatestZarr) }

    xenium_qc_ch = qc_results_ch
        .filter { key, pairId, platform, enrichedLatestZarr, qcOutDir -> platform == "XENIUM" }
        .map { key, pairId, platform, enrichedLatestZarr, qcOutDir -> tuple(pairId, enrichedLatestZarr) }

    paired_zarrs_ch = merscope_qc_ch
        .join(xenium_qc_ch)
        .map { pairId, merscopeLatest, xeniumLatest ->
            tuple(pairId, file(merscopeLatest), file(xeniumLatest))
        }

    compare_results_ch = COMPARE(paired_zarrs_ch)

    compare_done_ch = compare_results_ch.map { pairId, compareOutDir ->
        tuple(pairId, true)
    }

    visualize_inputs_ch = paired_zarrs_ch
        .join(compare_done_ch)
        .map { pairId, merscopeLatest, xeniumLatest, doneFlag ->
            tuple(pairId, merscopeLatest, xeniumLatest)
        }

    VISUALIZE(visualize_inputs_ch)
}
