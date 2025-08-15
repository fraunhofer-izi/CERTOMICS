#!/usr/bin/env nextflow

def get_library_types (library_list) {
    def library_types = library_list.collect { library -> library.type }.unique()
    
    boolean gex = false
    boolean vdj_b = false
    boolean vdj_t = false
    boolean feat = false

    library_types.each() { type ->
        if ('Gene Expression'.equals(type)) {
            gex = true
        } else if ('VDJ-B'.equals(type)) {
            vdj_b = true
        } else if ('VDJ-T'.equals(type)) {
            vdj_t = true
        } else if ('Antibody Capture'.equals(type)) {
            feat = true
        } else {
            error "Invalid feature type: ${type}"
        }
    }

    return [gex, vdj_b, vdj_t, feat]
}



def get_sample_list_type_bits (sample_list)  {
    def library_types = get_library_types(sample_list.collect { sample -> sample.libraries }.flatten() )
    
    int type_bits = 0
    type_bits += library_types[0] ? 1    : 0
    type_bits += library_types[1] ? 10   : 0
    type_bits += library_types[2] ? 100  : 0
    type_bits += library_types[3] ? 1000 : 0
    return type_bits
}

process CELLRANGER_MULTI {
    publishDir "${params.outdir}/cellranger_multi/${task.tag}"
    label 'module_cellranger'
    label 'big_task'
    fair true
    tag "${sample.name}"

    input:
    path cluster_template, arity: '1'
    path libraries, stageAs: 'library', arity: '1..*'
    path gex_reference, stageAs: 'gex_reference'
    path vdj_reference, stageAs: 'vdj_reference'
    path feat_reference, stageAs: 'feat_reference'
    val  sample

    output:
    path 'output', emit: full
    path 'output/outs/per_sample_outs/*/web_summary.html', emit: web_summary
    path 'output/outs/per_sample_outs/*/count/sample_alignments.bam', emit: sample_alignments_bam
    path 'output/outs/per_sample_outs/*/count/sample_alignments.bam.bai', emit: sample_alignments_bai
    path 'output/outs/per_sample_outs/*/count/sample_filtered_feature_bc_matrix', emit: feature_bc_matrix, optional: true
    path 'output/outs/multi/count/raw_feature_bc_matrix', emit: feature_raw_bc_matrix, optional: true
    path 'output/outs/per_sample_outs/*/vdj_t/filtered_contig_annotations.csv', emit: vdj_t_annotations, optional: true
    path 'output/outs/per_sample_outs/*/vdj_b/filtered_contig_annotations.csv', emit: vdj_b_annotations, optional: true
    val  sample.name, emit: sample_name

    script:
    types = get_library_types(sample.libraries)
    libs = ['[libraries]', 'fastq_id,fastqs,feature_types']
    
    // using ${libraries[index]} instead of ${library.path} to get staged path
    sample.libraries.eachWithIndex { library, index ->
        libs.add([library.id, "\$(realpath -s ${libraries[index]})", library.type].join(','))
    }

    cr_args = params.cellranger_disable_ui ? ['--disable-ui'] : []
    if (params.cellranger_enable_cluster && cluster_template) {
        // cluster mode
        cr_args += [
            "--jobmode \"\$(realpath ${cluster_template})\"",
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
        if ${types[0]}; then
            echo "[gene-expression]\nreference,\$(realpath -s ${gex_reference})\ncreate-bam,true" >> multi.csv
        fi
    elif [[ "\$CR_VERSION" == "7" ]]; then
        if ${types[0]}; then
            echo "[gene-expression]\nreference,\$(realpath -s ${gex_reference})" >> multi.csv
        fi
    else
        echo "Unsupported cellranger version: \$CR_VERSION"
        exit 1
    fi

    if ${types[1] || types[2]}; then
        echo "[vdj]\nreference,\$(realpath -s ${vdj_reference})" >> multi.csv
    fi

    if ${types[3]}; then
        echo "[feature]\nreference,\$(realpath -s ${feat_reference})" >> multi.csv
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
    path helper_functions_script, arity: '1'
    path feature_bc_matrices, stageAs: 'feature_bc_matrix', arity: '1..*'
    path feature_raw_bc_matrix, stageAs: 'feature_raw_bc_matrix', arity: '1..*'
    path vdj_t_annotations, stageAs: 'vdj_t_annotation', arity: '1..*'
    path vdj_b_annotations, stageAs: 'vdj_b_annotation', arity: '1..*'
    path annotation, stageAs: 'annotation.gtf', arity: '1'
    val samples

    output:
    path 'seurat_object.rds'

    script:
    library_types = get_library_types(samples.collect { sample -> sample.libraries }.flatten() )
    add_matrices = library_types[0] || library_types[3]
    add_raw_matrices = library_types[0] 
    add_annotation_b = library_types[1]
    add_annotation_t = library_types[2]

    matrices = add_matrices ? feature_bc_matrices.join(',') : 'none'
    matrices_raw = add_raw_matrices ? feature_raw_bc_matrix.join(',') : 'none'
    annotation_t = add_annotation_t ? vdj_t_annotations.join(',') : 'none'
    annotation_b = add_annotation_b ? vdj_b_annotations.join(',') : 'none'

    """
    build_seurat_objects.R \
        ${helper_functions_script} \
        ${matrices} \
        ${matrices_raw} \
        ${annotation_t} \
        ${annotation_b} \
        ${samples.collect { sample -> "\"${sample.name}\"" }.join(',')} \
        "seurat_object.rds" \
        "${get_sample_list_type_bits(samples)}" \
        "${annotation}"
    """
}

process CAR_METRICS {
    label 'module_python3'

    input:
    path sample_alignments_bam, stageAs: 'dir*/sample_alignments.bam', arity: '1..*'
    path sample_alignments_bai, stageAs: 'dir*/sample_alignments.bam.bai', arity: '1..*'
    path car_fa, stageAs: 'car.fa', arity: '1'
    path car_gtf, stageAs: 'car.gtf', arity: '1'
    val samples
    path quant_dirs
    
    output:
    path 'results_metrics_reads_CAR.csv', emit: metrics
    path 'results_coverage_against_CAR.csv', emit: coverage
    path 'results_coverage_against_CAR_unique.csv', emit: coverage_unique
    path "CAR_est_counts_matrix.csv", emit: kallisto_matrix

    script:
    def kallisto_cmd = quant_dirs ? """
    Kallisto_Comparisons.py \\
        --input-dir ${quant_dirs.join(' ')} \\
        --output CAR_est_counts_matrix.csv
    """ : "touch CAR_est_counts_matrix.csv"

    """
    CAR_quality.py \
        --sample_names ${samples.collect { sample -> sample.name }.join(' ')} \
        --bam_files ${sample_alignments_bam.join(' ')} \
        --CAR_fasta_file ${car_fa} \
        --CAR_gtf_file ${car_gtf}

    ${kallisto_cmd}
    """
}

process QUARTO {
    publishDir "${params.outdir}/quarto"
    label 'module_quarto'

    input:
    path kallisto_matrix, arity: '1' 
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
        -P libraries:"${get_sample_list_type_bits(samples)}" \
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

workflow RUN_SECONDARY_ANALYSIS {
    take:
    // sample list
    samples

    // paths
    gex_reference
    vdj_reference
    feat_reference

    //CAR references
    car_fa
    car_gtf
    multiple_car_fa // for kallisto indexing & quantification

    main:
    CELLRANGER_MULTI (
        params.cellranger_cluster_template ?: [],
        samples.map { sample -> sample.libraries.collect { library -> library.path } },
        gex_reference.value ?: [],
        vdj_reference.value ?: [],
        feat_reference.value ?: [],
        samples
    )

    do_sub_workflow = (car_fa.value && car_gtf.value)
    if (multiple_car_fa.value) {
        KALLISTO_INDEX(multiple_car_fa)
        KALLISTO_QUANT(
            KALLISTO_INDEX.out, 
            samples.map { sample -> sample.libraries.collect { library -> library.path } },
            samples.map { sample -> sample.name })
    }

    if (do_sub_workflow) {
        SEURAT_OBJECT (
            projectDir.resolve('bin/helper_functions.R'),
            CELLRANGER_MULTI.out.feature_bc_matrix.collect(),
            CELLRANGER_MULTI.out.feature_raw_bc_matrix.collect(),
            CELLRANGER_MULTI.out.vdj_t_annotations.collect(),
            CELLRANGER_MULTI.out.vdj_b_annotations.collect(),
            car_gtf,
            samples.collect()
        )

        CAR_METRICS (
            CELLRANGER_MULTI.out.sample_alignments_bam.collect(),
            CELLRANGER_MULTI.out.sample_alignments_bai.collect(),
            car_fa,
            car_gtf,
            samples.collect(),
            multiple_car_fa.value ? KALLISTO_QUANT.out.collect() : []
        )

        QUARTO (
            CAR_METRICS.out.kallisto_matrix,
            projectDir.resolve('bin/CAR_plot.qmd'),
            projectDir.resolve('bin/CAR_quality_plot.py'),
            projectDir.resolve('bin/helper_functions.R'),
            SEURAT_OBJECT.out,
            car_gtf,
            CAR_METRICS.out.metrics,
            CAR_METRICS.out.coverage,
            CAR_METRICS.out.coverage_unique,
            samples.collect()
        )
    }

    emit:
    cellranger_web_summary = CELLRANGER_MULTI.out.web_summary
    cellranger_full = CELLRANGER_MULTI.out.full
    seurat_obj = do_sub_workflow ? SEURAT_OBJECT.out : []
    quarto_out = do_sub_workflow ? QUARTO.out : []
}
