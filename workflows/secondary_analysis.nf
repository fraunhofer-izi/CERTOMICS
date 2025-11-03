#!/usr/bin/env nextflow

include {
    getNullFile;
    isNullFile;
    getLibraryTypes;
    getSampleLibraryPaths;
    getSampleNames
} from '../modules/local/functions'

def getSampleListTypeBits (SampleList)  {
    def libraryTypes = getLibraryTypes(SampleList.collect { sample -> sample.libraries }.flatten())

    int type_bits = 0
    type_bits += libraryTypes.gex ? 1    : 0
    type_bits += libraryTypes.vdj_b ? 10   : 0
    type_bits += libraryTypes.vdj_t ? 100  : 0
    type_bits += libraryTypes.feat ? 1000 : 0
    return type_bits
}

process CELLRANGER_MULTI {
    publishDir "${params.outdir}/cellranger_multi/${task.tag}"
    label 'module_cellranger'
    label 'big_task'
    fair true
    tag "${sample.name}"

    input:
    path clusterTemplate, stageAs: 'cluster/*', arity: '1'
    path gexReference, stageAs: 'references/gex/*', arity: '1'
    path vdjReference, stageAs: 'references/vdj/*', arity: '1'
    path featureReference, stageAs: 'references/feature/*', arity: '1'
    tuple val(sample), path(libraries, stageAs: 'libraries/lib*', arity: '1..*')

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
    println("NAMES:" + gexReference.getSimpleName() + " " + isNullFile(gexReference))

    library_types = getLibraryTypes(sample.libraries)
    if (isNullFile(gexReference) && library_types.gex)
        error('Gene expression reference needed but not provided.')
    if (isNullFile(vdjReference) && (library_types.vdj_b || library_types.vdj_t))
        error('VDJ reference needed but not provided')
    if (isNullFile(featureReference) && library_types.feat)
        error('Feature reference needed but not provided.')
    
    // using ${libraries[index]} instead of ${library.path} to get staged path
    libs = ['[libraries]', 'fastq_id,fastqs,feature_types']
    sample.libraries.eachWithIndex { library, index ->
        libs.add([library.id, "\$(realpath -s ${libraries[index]})", library.type].join(','))
    }

    cr_args = params.cellranger_disable_ui ? ['--disable-ui'] : []
    if (params.cellranger_enable_cluster) {
        // cluster mode
        if (isNullFile(clusterTemplate))
            error('Trying to run CELLRANGER_MULTI process in cluster-mode without cluster template')
        cr_args += [
            "--jobmode \"\$(realpath ${clusterTemplate})\"",
            "--mempercore ${params.cellranger_mem_per_core}",
            "--maxjobs ${params.cellranger_max_jobs}"
        ]
    } else {
        // local mode
        if (task.cpus) cr_args += "--localcores ${task.cpus}"
        if (task.memory) cr_args += "--localmem ${task.memory.toGiga()}"
    }

    """
    # Building multi.csv
    touch multi.csv
    CR_VERSION=\$(cellranger --version | sed -n 's/.*cellranger-\\([0-9]\\+\\).*/\\1/p')
    if [[ "\$CR_VERSION" == "8" ]]; then
        if ${library_types.gex}; then
            echo "[gene-expression]\nreference,\$(realpath -s ${gexReference})\ncreate-bam,true" >> multi.csv
        fi
    elif [[ "\$CR_VERSION" == "7" ]]; then
        if ${library_types.gex}; then
            echo "[gene-expression]\nreference,\$(realpath -s ${gexReference})" >> multi.csv
        fi
    else
        echo "Unsupported cellranger version: \$CR_VERSION"
        exit 1
    fi

    if ${library_types.vdj_b || library_types.vdj_t}; then
        echo "[vdj]\nreference,\$(realpath -s ${vdjReference})" >> multi.csv
    fi

    if ${library_types.feat}; then
        echo "[feature]\nreference,\$(realpath -s ${featureReference})" >> multi.csv
    fi

    echo "${libs.join('\n')}" >> multi.csv

    # Running CR Multi
    cellranger multi --csv="multi.csv" --id=${sample.name} --output-dir=output ${cr_args.join(' ')}

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
    path vdjTAnnotations, stageAs: 'vdj_t_annotation', arity: '1..*'
    path vdjBAnnotations, stageAs: 'vdj_b_annotation', arity: '1..*'
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
    label 'module_python3'

    input:
    path sampleAlignmentBams, stageAs: 'bams*/', arity: '1..*'
    path sampleAlignmentBais, stageAs: 'bams*/', arity: '1..*'
    path carFasta, stageAs: 'car.fa', arity: '1'
    path carGtf, stageAs: 'car.gtf', arity: '1'
    path quantDirs, stageAs: 'quant_dirs/*'
    val samples
    
    output:
    path 'results_metrics_reads_CAR.csv', emit: metrics
    path 'results_coverage_against_CAR.csv', emit: coverage
    path 'results_coverage_against_CAR_unique.csv', emit: uniqueCoverage
    path "CAR_est_counts_matrix.csv", optional: true, emit: kallistoMatrix

    script:
    """
    CAR_quality.py \
        --sample_names ${samples.collect { sample -> sample.name }.join(' ')} \
        --bam_files ${sampleAlignmentBams.join(' ')} \
        --CAR_fasta_file ${carFasta} \
        --CAR_gtf_file ${carGtf}

    if ${!isNullFile(quantDirs)}; then
        Kallisto_Comparisons.py \
            --input-dir ${quantDirs.join(' ')} \
            --output CAR_est_counts_matrix.csv
    fi
    """
}

process QUARTO {
    publishDir "${params.outdir}/quarto"
    label 'module_quarto'

    input:
    path kallisto_matrix, arity: '0..*'
    path car_plot_qmd, arity: '1'
    path car_quality_plot_py, arity: '1'
    path helper_functions, arity: '1'
    path seurat_object, stageAs: 'seurat_object', arity: '1'
    path annotation, stageAs: 'annotation', arity: '1'
    path metrics_reads_car, arity: '1'
    path metrics_coverage_car, arity: '1'
    path metrics_coverage_car_unique, arity: '1'
    val  samples

    output:
    path 'metrics_html'

    script:
    """
    export HOME=\$(realpath "quarto-cache")
    quarto render ${car_plot_qmd} \
        -P kallisto_matrix:${kallisto_matrix} \
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


process KALLISTO_INDEX {
    label 'module_kallisto'
    input:
    path fasta

    output:
    path "*.idx"

    script:
    def fasta_basename = fasta.getBaseName()  // strips .fa/.fasta
    """
    kallisto index -i ${fasta_basename}.idx ${fasta}
    """
}

process KALLISTO_QUANT {
    label 'module_kallisto'
    fair true

    input:
    path index_file
    path libraries, stageAs: 'library', arity: '1..*'
    val sample_name

    output:
    path "quant_${sample_name}"

    script:
    """
    echo "Running kallisto for ${sample_name}"

    READS=\$(for lib in library*; do find -L "\$lib" -type f -name '*_R_*R2_*.fastq.gz' | head -n 1; done | sort | tr '\\n' ' ')

    echo "Using reads: \$READS"

    kallisto quant \\
        -i ${index_file} \\
        -o quant_${sample_name} \\
        --single -l 350 -s 50 \\
        --single-overhang \\
        \$READS
    """
}

workflow SECONDARY_ANALYSIS {
    take:
    // sample list
    samples

    // paths
    gexReference
    vdjReference
    featureReference

    //CAR references
    carFasta
    carGtf
    
     // for kallisto indexing & quantification
    multiCarFasta
    scGateModel

    // misc
    cellrangerClusterTemplate

    main:
    sample_libs = samples.map { sample ->
        def libraries = sample.libraries.collect { library -> file(library.path) }
        tuple(sample, libraries)
    }

    CELLRANGER_MULTI (
        cellrangerClusterTemplate,
        gexReference,
        vdjReference,
        featureReference,
        sample_libs
    )

    doKallistoWorkflow = !isNullFile(multiCarFasta)
    if (doKallistoWorkflow) {
        KALLISTO_INDEX(
            multiCarFasta
        )
        
        KALLISTO_QUANT(
            KALLISTO_INDEX.out,
            getSampleLibraryPaths(samples),
            getSampleNames(samples)
        )
    }

    doCarWorkflow = !isNullFile(carFasta) && !isNullFile(carGtf)
    if (doCarWorkflow) {
        SEURAT_OBJECT (
            projectDir.resolve('bin/helper_functions.R'),
            CELLRANGER_MULTI.out.featureBcMatrix.collect(),
            CELLRANGER_MULTI.out.rawFeatureBcMatrix.collect(),
            CELLRANGER_MULTI.out.vdjTAnnotations.collect(),
            CELLRANGER_MULTI.out.vdjBAnnotations.collect(),
            carGtf,
            samples.collect(),
            scGateModel
        )

        CAR_METRICS (
            CELLRANGER_MULTI.out.sampleAlignmentsBam.collect(),
            CELLRANGER_MULTI.out.sampleAlignmentsBai.collect(),
            carFasta,
            carGtf,
            samples.collect(),
            doKallistoWorkflow ? KALLISTO_QUANT.out.collect() : getNullFile()
        )

        QUARTO (
            CAR_METRICS.out.kallistoMatrix,
            projectDir.resolve('bin/CAR_plot.qmd'),
            projectDir.resolve('bin/CAR_quality_plot.py'),
            projectDir.resolve('bin/helper_functions.R'),
            SEURAT_OBJECT.out,
            carGtf,
            CAR_METRICS.out.metrics,
            CAR_METRICS.out.coverage,
            CAR_METRICS.out.uniqueCoverage,
            samples.collect()
        )
    }

    emit:
    cellranger_web_summary = CELLRANGER_MULTI.out.webSummary
    cellranger_full = CELLRANGER_MULTI.out.full
    seurat_obj = doCarWorkflow ? SEURAT_OBJECT.out : getNullFile()
    quarto_out = doCarWorkflow ? QUARTO.out : getNullFile()
}
