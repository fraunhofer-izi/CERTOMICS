#!/usr/bin/env nextflow
include { validateParameters ; paramsSummaryLog ; samplesheetToList } from 'plugin/nf-schema'
include { getNullFile ; isNullFile ; parseOptionalPath } from './modules/local/functions'

include { BUILD_CUSTOM_REFERENCE } from './workflows/handle_references'
include { SECONDARY_ANALYSIS } from './workflows/secondary_analysis'
include { QUALITY_CONTROL } from './workflows/quality_control'

def handleSourceUrls(sourceFile, sourceUrl, sourceUrlFallback) {
    if (sourceFile) {
        return channel.fromPath(sourceFile, checkIfExists: true)
    }
    if (sourceUrl) {
        return channel.fromPath(sourceUrl)
    }
    return channel.fromPath(sourceUrlFallback)
}

workflow REFERENCE {
    take:
    referenceVersion
    sourceFasta
    sourceGtf
    carFasta
    carGtf

    main:
    BUILD_CUSTOM_REFERENCE(
        referenceVersion,
        sourceFasta,
        sourceGtf,
        carFasta,
        carGtf,
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
    cellrangerClusterTemplate

    main:
    SECONDARY_ANALYSIS(
        samples,
        gexReference,
        vdjReference,
        featureReference,
        carFasta,
        carGtf,
        multiCarFasta,
        scGateModel,
        cellrangerClusterTemplate,
    )

    if (!params.skip_qc) {
        // Run QC
        QUALITY_CONTROL(
            samples,
            SECONDARY_ANALYSIS.out.cellranger_web_summary,
            params.skip_fastqc,
            params.skip_fastq_screen,
            params.skip_multiqc,
            params.multiqc_config,
        )
    }
}

workflow {
    // validation and Help message
    validateParameters()
    log.info(paramsSummaryLog(workflow))

    // read samples
    samples = channel.empty()
    if (params.samplesheet) {
        samples = channel.fromList(
                samplesheetToList(
                    params.samplesheet,
                    projectDir.resolve('assets/schemas/samples.json'),
                )
            )
            .map { sampleName, libraries ->
                Sample.create(
                    sampleName,
                    libraries.collect { id, path, type -> ['fastq_id': id, 'fastqs': path, 'feature_types': type] },
                )
            }
    }

    urlMap = samplesheetToList(
        projectDir.resolve('assets/reference_source_urls.yaml'),
        projectDir.resolve('assets/schemas/reference_source_urls.json'),
    ).collectEntries { version, fa, gtf -> [(version): ['fa': fa, 'gtf': gtf]] }

    // check / update parameters
    // prebuilt references
    gexReference = parseOptionalPath(params.gene_expression_reference)
    vdjReference = parseOptionalPath(params.vdj_reference)
    featureReference = parseOptionalPath(params.feature_reference)

    // custom reference data
    gexCarFasta = parseOptionalPath(params.gene_expression_car_fa)
    gexCarGtf = parseOptionalPath(params.gene_expression_car_gtf)
    referenceVersion = params.gene_expression_reference_version

    if (!urlMap.containsKey(referenceVersion)) {
        error("Missing source urls for reference version '${referenceVersion}'.")
    }

    gexSourceFasta = handleSourceUrls(
        params.gene_expression_source_fa,
        params.gene_expression_source_fa_url,
        urlMap[referenceVersion].fa,
    )

    gexSourceGtf = handleSourceUrls(
        params.gene_expression_source_gtf,
        params.gene_expression_source_gtf_url,
        urlMap[referenceVersion].gtf,
    )

    // misc
    multiCarFasta = parseOptionalPath(params.multiple_car_fa)
    cellrangerClusterTemplate = parseOptionalPath(params.cellranger_cluster_template)

    if (params.pipeline_mode == null) {
        error('Parameter "pipeline_mode" cannot be null.')
    }
    else if (params.pipeline_mode == 'reference') {
        REFERENCE(
            referenceVersion,
            gexSourceFasta,
            gexSourceGtf,
            gexCarFasta,
            gexCarGtf,
        )
    }
    else if (params.pipeline_mode == 'analysis') {
        ANALYSIS(
            samples,
            gexReference,
            vdjReference,
            featureReference,
            gexCarFasta,
            gexCarGtf,
            multiCarFasta,
            params.scGate_model,
            cellrangerClusterTemplate,
        )
    }
    else if (params.pipeline_mode == 'full') {
        doBuildReference = isNullFile(gexReference) && samples.collect { sample -> sample.hasGeneExpressionLibrary() }.any()
        if (doBuildReference) {
            REFERENCE(
                referenceVersion,
                gexSourceFasta,
                gexSourceGtf,
                gexCarFasta,
                gexCarGtf,
            )

            gexReference = REFERENCE.out.reference.first()
        }

        ANALYSIS(
            samples,
            gexReference,
            vdjReference,
            featureReference,
            gexCarFasta,
            gexCarGtf,
            multiCarFasta,
            params.scGate_model,
            cellrangerClusterTemplate,
        )
    }
    else {
        error("Unknown pipeline mode: '${params.pipeline_mode}'")
    }
}
