#!/usr/bin/env nextflow

include {
    getNullFile ;
    isNullFile ;
    getLibraryTypes ;
    getSampleLibraryPaths ;
    getSampleNames
} from '../modules/local/functions'

def getSampleListTypeBits(SampleList) {
    def libraryTypes = getLibraryTypes(SampleList.collect { sample -> sample.libraries }.flatten())

    def type_bits = 0
    type_bits += libraryTypes.gex ? 1 : 0
    type_bits += libraryTypes.vdj_b ? 10 : 0
    type_bits += libraryTypes.vdj_t ? 100 : 0
    type_bits += libraryTypes.feat ? 1000 : 0
    return type_bits
}

process CELLRANGER_MULTI {
    publishDir "${params.outdir}/cellranger_multi/${task.tag}"
    label 'module_cellranger'
    label 'big_task'
    fair true
    tag "${sampleName}"

    input:
    path clusterTemplate, stageAs: 'cluster/*', arity: '1'
    path gexReference, stageAs: 'references/gex/*', arity: '1'
    path vdjReference, stageAs: 'references/vdj/*', arity: '1'
    path featureReference, stageAs: 'references/feature/*', arity: '1'
    tuple path(libraryPaths, stageAs: 'libraries/lib*', arity: '1..*'), val(libraryList), val(sampleName)

    output:
    path 'output', emit: full
    path 'output/outs/per_sample_outs/*/web_summary.html', emit: webSummary
    path 'output/outs/per_sample_outs/*/count/sample_alignments.bam', emit: sampleAlignmentsBam
    path 'output/outs/per_sample_outs/*/count/sample_alignments.bam.bai', emit: sampleAlignmentsBai
    path 'output/outs/per_sample_outs/*/count/sample_filtered_feature_bc_matrix', emit: featureBcMatrix, optional: true
    path 'output/outs/multi/count/raw_feature_bc_matrix', emit: rawFeatureBcMatrix, optional: true
    path 'output/outs/per_sample_outs/*/vdj_t/filtered_contig_annotations.csv', emit: vdjTAnnotations, optional: true
    path 'output/outs/per_sample_outs/*/vdj_b/filtered_contig_annotations.csv', emit: vdjBAnnotations, optional: true

    script:
    libraryTypes = getLibraryTypes(libraryList)
    if (isNullFile(gexReference) && libraryTypes.gex) {
        error('Gene expression reference needed but not provided.')
    }
    if (isNullFile(vdjReference) && (libraryTypes.vdj_b || libraryTypes.vdj_t)) {
        error('VDJ reference needed but not provided')
    }
    if (isNullFile(featureReference) && libraryTypes.feat) {
        error('Feature reference needed but not provided.')
    }

    libs = ['[libraries]', 'fastq_id,fastqs,feature_types']
    [libraryList, libraryPaths]
        .transpose()
        .each { libraryObject, libraryPath ->
            libs.add(
                [
                    libraryObject.id,
                    "\$(realpath -s ${libraryPath})",
                    libraryObject.type,
                ].join(',')
            )
        }

    cr_args = params.cellranger_disable_ui ? ['--disable-ui'] : []
    if (params.cellranger_enable_cluster) {
        // cluster mode
        if (isNullFile(clusterTemplate)) {
            error('Trying to run CELLRANGER_MULTI process in cluster-mode without cluster template')
        }
        cr_args += [
            "--jobmode \"\$(realpath ${clusterTemplate})\"",
            "--mempercore ${params.cellranger_mem_per_core}",
            "--maxjobs ${params.cellranger_max_jobs}",
        ]
    }
    else {
        // local mode
        if (task.cpus) {
            cr_args += "--localcores ${task.cpus}"
        }
        if (task.memory) {
            cr_args += "--localmem ${task.memory.toGiga()}"
        }
    }

    """
    # Building multi.csv
    touch multi.csv
    CR_VERSION=\$(cellranger --version | sed -n 's/.*cellranger-\\([0-9]\\+\\).*/\\1/p')
    if [[ "\$CR_VERSION" == "8" ]]; then
        if ${libraryTypes.gex}; then
            echo "[gene-expression]\nreference,\$(realpath -s ${gexReference})\ncreate-bam,true" >> multi.csv
        fi
    elif [[ "\$CR_VERSION" == "7" ]]; then
        if ${libraryTypes.gex}; then
            echo "[gene-expression]\nreference,\$(realpath -s ${gexReference})" >> multi.csv
        fi
    else
        echo "Unsupported cellranger version: \$CR_VERSION"
        exit 1
    fi

    if ${libraryTypes.vdj_b || libraryTypes.vdj_t}; then
        echo "[vdj]\nreference,\$(realpath -s ${vdjReference})" >> multi.csv
    fi

    if ${libraryTypes.feat}; then
        echo "[feature]\nreference,\$(realpath -s ${featureReference})" >> multi.csv
    fi

    echo "${libs.join('\n')}" >> multi.csv

    # Running CR Multi
    cellranger multi --csv="multi.csv" --id=${sampleName} --output-dir=output ${cr_args.join(' ')}

    # Verify output
    pso="output/outs/per_sample_outs/"
    dirs=(\${pso}/*/)
    if [ \${#dirs[@]} -ne 1 ]; then
        echo "Unexpected output structure. Expected one directory in \${pso}" >&2
        exit 1
    fi
    """
}

process SEURAT_OBJECT {
    publishDir "${params.outdir}/seurat_object"
    label 'module_seurat'
    label 'big_task'

    input:
    path helperFunctionsScript, arity: '1'
    path featureBcMatrices, stageAs: 'feature_bc_matrix', arity: '1..*'
    path rawFeatureBcMatrices, stageAs: 'feature_raw_bc_matrix', arity: '1..*'
    val  vdjTAnnotations
    val  vdjBAnnotations
    path annotation, stageAs: 'annotation.gtf', arity: '1'
    val samples
    val scGateModel

    output:
    path 'seurat_object.rds'

    script:
    libraryTypes = getLibraryTypes(samples.collect { sample -> sample.libraries }.flatten())
    addMatrices = libraryTypes.gex || libraryTypes.feat
    addRawMatrices = libraryTypes.gex
    addAnnotationB = libraryTypes.vdj_b
    addAnnotationT = libraryTypes.vdj_t

    matrices = addMatrices ? featureBcMatrices.join(',') : 'none'
    matricesRaw = addRawMatrices ? rawFeatureBcMatrices.join(',') : 'none'
    annotationT = addAnnotationT ? vdjTAnnotations.join(',') : 'none'
    annotationB = addAnnotationB ? vdjBAnnotations.join(',') : 'none'

    """
    build_seurat_objects.R \
        ${helperFunctionsScript} \
        ${matrices} \
        ${matricesRaw} \
        ${annotationT} \
        ${annotationB} \
        ${samples.collect { sample -> "\"${sample.name}\"" }.join(',')} \
        "seurat_object.rds" \
        "${getSampleListTypeBits(samples)}" \
        "${annotation}" \
        "${scGateModel}"
    """
}

process CAR_METRICS {
    publishDir "${params.outdir}/car_metrics"
    label 'module_python3'

    input:
    path sampleAlignmentBams, stageAs: 'bams*/', arity: '1..*'
    path sampleAlignmentBais, stageAs: 'bams*/', arity: '1..*'
    path carFasta, stageAs: 'car.fa', arity: '1'
    path carGtf, stageAs: 'car.gtf', arity: '1'
    val samples

    output:
    path 'results_metrics_reads_CAR.csv', emit: metrics
    path 'results_coverage_against_CAR.csv', emit: coverage
    path 'results_coverage_against_CAR_unique.csv', emit: uniqueCoverage

    script:
    """
    CAR_quality.py \
        --sample_names ${samples.collect { sample -> sample.name }.join(' ')} \
        --bam_files ${sampleAlignmentBams.join(' ')} \
        --CAR_fasta_file ${carFasta} \
        --CAR_gtf_file ${carGtf}
    """
}

process QUARTO {
    publishDir "${params.outdir}/quarto"
    label 'module_quarto'

    input:
    path mapped_kmer_output, stageAs: 'mapped_kmer_outputs/*'
    path unmapped_kmer_output, stageAs: 'unmapped_kmer_outputs/*'
    path car_plot_qmd, arity: '1'
    path car_quality_plot_py, arity: '1'
    path helper_functions, arity: '1'
    path seurat_object, stageAs: 'seurat_object', arity: '1'
    path annotation, stageAs: 'annotation', arity: '1'
    path metrics_reads_car, arity: '1'
    path metrics_coverage_car, arity: '1'
    path metrics_coverage_car_unique, arity: '1'
    val samples

    output:
    path 'metrics_html'

    script:
    """
    export HOME=\$(realpath "quarto-cache")
    quarto render ${car_plot_qmd} \
        -P mapped_kmer_dir:"mapped_kmer_outputs" \
        -P unmapped_kmer_dir:"unmapped_kmer_outputs" \
        -P seurat_object:"${seurat_object}" \
        -P gtf:"${annotation}" \
        -P results_metrics_reads_CAR:"${metrics_reads_car}" \
        -P results_coverage_against_CAR:"${metrics_coverage_car}" \
        -P results_coverage_against_CAR_unique:"${metrics_coverage_car_unique}" \
        -P libraries:"${getSampleListTypeBits(samples)}" \
        --no-cache

    mkdir metrics_html
    mv *.html metrics_html/
    mv CAR_plot_files metrics_html/
    """
}
process MAKE_UNIQUE_KMERS {
    label 'module_car_identity'

    input:
    path carFasta, stageAs: 'car.fa'
    path kmerScript

    output:
    path "kmers_out"

    script:
    """
    python $kmerScript car.fa 31 kmers_out
    """
}


process CAR_IDENTITY_MAPPED {
    label 'module_car_identity'
    publishDir "${params.outdir}/car_identity/mapped", mode: 'copy'

    input:
    tuple path(bam), val(sample), path(kmers), val(ref_tg)

    output:
    path "${sample}.CAR_kmer_summary.txt"

    script:
    """
    samtools index -c $bam || true

    bash compare_all_CARs_against_reference.sh "$bam" "$sample" "$ref_tg" "$kmers"
    """
}
process CAR_IDENTITY_UNMAPPED {
    label 'module_car_identity'
    publishDir "${params.outdir}/car_identity/unmapped", mode: 'copy'

    cpus 8

    input:
    tuple path(bam), val(sample), path(kmers)

    output:
    path "${sample}.unmapped_kmc_intersect_summary.txt"

    script:
    """
    samtools index -c $bam || true

    bash compare_all_CARs_against_unmapped.sh "$bam" "$sample" "$kmers" 31 "$task.cpus" 
    """
}
process EXTRACT_REFERENCE_TG {
    label 'module_car_identity'

    input:
    path carFasta, stageAs: 'car.fa'

    output:
    stdout

    script:
    """
    grep '^>' car.fa | head -n 1 | sed 's/^>//' | cut -d' ' -f1 | tr -d '\\n'
    """
}


workflow SECONDARY_ANALYSIS {
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
    sample_libs = samples.map { sample ->
        tuple(
            sample.libraries.collect { library -> library.path },
            sample.libraries,
            sample.name,
        )
    }

    CELLRANGER_MULTI(
        cellrangerClusterTemplate,
        gexReference,
        vdjReference,
        featureReference,
        sample_libs,
    )

    doKmerWorkflow = !isNullFile(multiCarFasta)
    if (doKmerWorkflow) {
        secondaryBamsCh = CELLRANGER_MULTI.out.sampleAlignmentsBam
            .map { bam ->
                tuple(
                    bam,
                    bam.parent.parent.name
                )
            }
        // Stage CAR fasta
        carFastaCh = Channel.fromPath(carFasta)
        multiCarFastaCh = Channel.fromPath(multiCarFasta)
        // Extract reference CAR once
        ref_tg_ch = EXTRACT_REFERENCE_TG(carFastaCh)
        // Generate unique kmers once
        kmers_ch = MAKE_UNIQUE_KMERS(
            multiCarFastaCh,
            projectDir.resolve('bin/kmers_similarity.py')
        )
        car_identity_input_ch = secondaryBamsCh
            .combine(kmers_ch)
            .combine(ref_tg_ch)
        
        // Run per-sample CAR identity (mapped)
        CAR_IDENTITY_MAPPED(car_identity_input_ch)
        // Run per-sample CAR identity (unmapped)
        car_identity_unmapped_input_ch = secondaryBamsCh
            .combine(kmers_ch)

        CAR_IDENTITY_UNMAPPED(car_identity_unmapped_input_ch)
        }

    doCarWorkflow = !isNullFile(carFasta) && !isNullFile(carGtf)
    if (doCarWorkflow) {

        vdjT_ch = CELLRANGER_MULTI.out.vdjTAnnotations
            .map { it.toString() }
            .ifEmpty { Channel.value('none') }

        vdjB_ch = CELLRANGER_MULTI.out.vdjBAnnotations
                    .map { it.toString() }
                    .ifEmpty { Channel.value('none') }


        SEURAT_OBJECT(
            projectDir.resolve('bin/helper_functions.R'),
            CELLRANGER_MULTI.out.featureBcMatrix.collect(),
            CELLRANGER_MULTI.out.rawFeatureBcMatrix.collect(),
            vdjT_ch.collect(),
            vdjB_ch.collect(),
            carGtf,
            samples.collect(),
            scGateModel,
        )

        CAR_METRICS(
            CELLRANGER_MULTI.out.sampleAlignmentsBam.collect(),
            CELLRANGER_MULTI.out.sampleAlignmentsBai.collect(),
            carFasta,
            carGtf,
            samples.collect(),
        )

        QUARTO(
            doKmerWorkflow ? CAR_IDENTITY_MAPPED.out.collect() : getNullFile(),
            doKmerWorkflow ? CAR_IDENTITY_UNMAPPED.out.collect() : getNullFile(),
            projectDir.resolve('bin/CAR_plot.qmd'),
            projectDir.resolve('bin/CAR_quality_plot.py'),
            projectDir.resolve('bin/helper_functions.R'),
            SEURAT_OBJECT.out,
            carGtf,
            CAR_METRICS.out.metrics,
            CAR_METRICS.out.coverage,
            CAR_METRICS.out.uniqueCoverage,
            samples.collect(),
        )
    }

    emit:
    cellranger_web_summary = CELLRANGER_MULTI.out.webSummary
    cellranger_full = CELLRANGER_MULTI.out.full
    seurat_obj = doCarWorkflow ? SEURAT_OBJECT.out : getNullFile()
    quarto_out = doCarWorkflow ? QUARTO.out : getNullFile()
}
