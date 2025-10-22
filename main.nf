#!/usr/bin/env nextflow

include { RUN_NF_VALIDATION      } from './workflows/initiate_pipeline'
include { HANDLE_GEX_REFERENCE   } from './workflows/handle_references'
include { RUN_SECONDARY_ANALYSIS } from './workflows/secondary_analysis'
include { RUN_QUALITY_CONTROL    } from './workflows/quality_control'

workflow REFERENCE {
    take:
    ref_ver

    // paths
    src_fa
    src_fa_url
    src_gtf
    src_gtf_url
    car_fa
    car_gtf

    main:
    HANDLE_GEX_REFERENCE (
        ref_ver,
        (src_fa.getName() == 'NO_FILE') ? src_fa_url : src_fa,
        (src_gtf.getName() == 'NO_FILE') ? src_gtf_url : src_gtf,
        car_fa,
        car_gtf
    )

    emit:
    reference = HANDLE_GEX_REFERENCE.out
}

workflow ANALYSIS {
    take:
    samples
    gex_reference
    vdj_reference
    feat_reference
    car_fa
    car_gtf
    car_fa_multi
    scgate_model

    main:
    RUN_SECONDARY_ANALYSIS (
        samples,
        gex_reference,
        vdj_reference,
        feat_reference,
        car_fa,
        car_gtf,
        car_fa_multi,
        scgate_model
    )

    // Run QC
    RUN_QUALITY_CONTROL (
        samples,
        RUN_SECONDARY_ANALYSIS.out.cellranger_web_summary,
        params.skip_qc,
        params.skip_fastqc,
        params.skip_fastq_screen,
        params.skip_multiqc,
        params.fastq_screen_config,
        params.multiqc_config,
    )
}

workflow {
    // validation and Help message
    RUN_NF_VALIDATION (
        params.help,
        params.validate_params,
        projectDir.resolve('nextflow_schema.json')
    )

    // parse parameters
    /*
    PARSE_PARAMETERS ()
    safe_params = PARSE_PARAMETERS.out
    */

    if (params.pipeline_mode == null) {
        log.error('Parameter "pipeline_mode" cannot be null.')
        System.exit(1)
    } else if (params.pipeline_mode == 'reference') {
        // build / manage references
        def ref_ver = params.gene_expression_reference_version
        REFERENCE(
            ref_ver,
            file(params.gene_expression_source_fa),
            file(params.gene_expression_source_fa_url[ref_ver]),
            file(params.gene_expression_source_gtf),
            file(params.gene_expression_source_gtf_url[ref_ver]),
            file(params.gene_expression_car_fa),
            file(params.gene_expression_car_gtf)
        )
    } else if (params.pipeline_mode == 'analysis') {
        ANALYSIS(
            params.samples.collect { sampleMap -> Sample.create(sampleMap) },
            file(params.gene_expression_reference),
            file(params.vdj_reference),
            file(params.feature_reference),
            file(params.gene_expression_car_fa),
            file(params.gene_expression_car_gtf),
            file(params.multiple_car_fa),
            params.scGate_model
        )
    } else {
        log.error("Pipeline mode \"${params.pipeline_mode}\" is not supported")
        System.exit(1)
    }
}
