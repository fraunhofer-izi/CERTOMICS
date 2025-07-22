#!/usr/bin/env nextflow

include { RUN_NF_VALIDATION      } from './workflows/initiate_pipeline'
include { PARSE_PARAMETERS       } from './workflows/initiate_pipeline'
include { HANDLE_GEX_REFERENCE   } from './workflows/handle_references'
include { RUN_SECONDARY_ANALYSIS } from './workflows/secondary_analysis'
include { RUN_QUALITY_CONTROL    } from './workflows/quality_control'

workflow {
    // validation and Help message
    RUN_NF_VALIDATION (
        params.help,
        params.validate_params,
        projectDir.resolve('nextflow_schema.json')
    )

    // parse parameters
    PARSE_PARAMETERS ()
    safe_params = PARSE_PARAMETERS.out

    // build / manage references
    HANDLE_GEX_REFERENCE (
        safe_params.needs_gex,        // is reference needed?
        safe_params.gex_version,      // build version
        safe_params.gex_src_fa_url,   // (optional) url of fallback source_fa file
        safe_params.gex_src_gtf_url,  // (optional) url of fallback source_gtf file
        safe_params.gex_reference,    // (optional) path of prebuilt reference
        safe_params.gex_src_fa,       // (optional) path of source_fa
        safe_params.gex_src_gtf,      // (optional) path of source_gtf
        safe_params.gex_car_fa,       // (optional) path of car_fa
        safe_params.gex_car_gtf       // (optional) path of car_gtf
    )

    // Run analysis
    RUN_SECONDARY_ANALYSIS (
        safe_params.samples,
        HANDLE_GEX_REFERENCE.out,
        safe_params.vdj_reference,
        safe_params.feat_reference,
        safe_params.gex_car_fa,
        safe_params.gex_car_gtf
    )

    // Run QC
    RUN_QUALITY_CONTROL (
        safe_params.samples,
        RUN_SECONDARY_ANALYSIS.out.cellranger_web_summary,
        params.skip_qc,
        params.skip_fastqc,
        params.skip_fastq_screen,
        params.skip_multiqc,
        params.fastq_screen_config,
        params.multiqc_config,
    )
}
