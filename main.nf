#!/usr/bin/env nextflow

include { NF_VALIDATION } from './workflows/initiate_pipeline'
include { BUILD_CUSTOM_REFERENCE } from './workflows/handle_references'
include { SECONDARY_ANALYSIS } from './workflows/secondary_analysis'
include { QUALITY_CONTROL } from './workflows/quality_control'
include { getNullFile; isNullFile; parseOptionalPath } from './modules/local/functions'

workflow REFERENCE {
    take:
    referenceVersion
    sourceFasta
    sourceGtf
    carFasta
    carGtf

    main:
    BUILD_CUSTOM_REFERENCE (
        referenceVersion,
        sourceFasta,
        sourceGtf,
        carFasta,
        carGtf
    )

    emit:
    reference = BUILD_CUSTOM_REFERENCE.out
}

workflow ANALYSIS {
    take:
    samples
    gexReference
    vdjReference
    featureReference
    carFasta
    carGtf
    multiCarFasta
    scGateModel

    main:
    SECONDARY_ANALYSIS (
        samples,
        gexReference,
        vdjReference,
        featureReference,
        carFasta,
        carGtf,
        multiCarFasta,
        scGateModel
    )

    // Run QC
    QUALITY_CONTROL (
        samples,
        params.skip_fastqc,
        params.skip_fastq_screen,
        params.skip_multiqc,
        params.fastq_screen_config,
        params.multiqc_config,
    )
}

workflow {
    // validation and Help message
    NF_VALIDATION (
        params.help,
        params.validate_params,
        projectDir.resolve('nextflow_schema.json')
    )

    // check / update parameters
    // samples
    samples = params.samples.collect { sampleMap -> Sample.create(sampleMap) }
    
    // prebuilt references
    gexReference = parseOptionalPath(params.gene_expression_reference)
    vdjReference = parseOptionalPath(params.vdj_reference)
    featureReference = parseOptionalPath(params.feature_reference)

    // custom reference data
    gexCarFasta = parseOptionalPath(params.gene_expression_car_fa)
    gexCarGtf = parseOptionalPath(params.gene_expression_car_gtf)
    referenceVersion = params.gene_expression_reference_version
    gexSourceFasta = (params.gene_expression_source_fa == null) ?
        file(params.gene_expression_source_fa_url[referenceVersion]) :
        file(params.gene_expression_source_fa, checkIfExists: true)
    gexSourceGtf = (params.gene_expression_source_gtf == null) ?
        file(params.gene_expression_source_gtf_url[referenceVersion]) :
        file(params.gene_expression_source_gtf, checkIfExists: true)
    
    // misc
    multiCarFasta = parseOptionalPath(params.multiple_car_fa)

    if (params.pipeline_mode == null) {
        error('Parameter "pipeline_mode" cannot be null.')
    } else if (params.pipeline_mode == 'reference') {
        REFERENCE (
            referenceVersion,
            gexSourceFasta,
            gexSourceGtf,
            gexCarFasta,
            gexCarGtf
        )
    } else if (params.pipeline_mode == 'analysis') {
        ANALYSIS(
            samples,
            gexReference,
            vdjReference,
            featureReference,
            gexCarFasta,
            gexCarGtf,
            multiCarFasta,
            params.scGate_model
        )
    } else if (params.pipeline_mode == 'full') {
        doBuildReference = samples.collect { sample -> sample.hasGeneExpressionLibrary() }.any() && isNullFile(gexReference)
        if (doBuildReference) {
            gexReference = REFERENCE (
                referenceVersion,
                gexSourceFasta,
                gexSourceGtf,
                gexCarFasta,
                gexCarGtf
            ).out
        }

        ANALYSIS (
            samples,
            gexReference,
            vdjReference,
            featureReference,
            gexCarFasta,
            gexCarGtf,
            multiCarFasta,
            params.scGate_model
        )
    } else {
        error("Unknown pipeline mode: '${params.pipeline_mode}'")
    }
}
