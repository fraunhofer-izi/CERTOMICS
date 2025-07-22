#!/usr/bin/env nextflow

process FASTQC {
    publishDir "${params.outdir}/fastqc/${task.tag}"
    label 'module_fastqc'
    tag "${sample_name}"

    input:
    path fastq_paths, stageAs: 'fastq', arity: '1..*'
    val fastq_ids
    val sample_name

    output:
    path "fastqc"

    script:
    libraries = []
    fastq_paths.eachWithIndex { directory, index ->
        libraries.add(
            directory.resolve("${fastq_ids[index]}_*")
        )
    }

    """
    mkdir "fastqc"
    fastqc -o "fastqc" --noextract --threads ${task.cpus} ${libraries.join(' ')}
    """
}

process FASTQ_SCREEN {
    publishDir "${params.outdir}/fastq_screen/${task.tag}"
    label 'module_fastq_screen'
    tag "${sample_name}"

    input:
    path fastq_screen_config, arity: '1'
    path fastq_paths, stageAs: 'fastq', arity: '1..*'
    val fastq_ids
    val fastq_types
    val sample_name

    output:
    path "fastqs"

    script:
    libraries = []
    fastq_paths.eachWithIndex { directory, index ->
        if ('Gene Expression'.equals(fastq_types[index])) {
            libraries.add(
                directory.resolve("${fastq_ids[index]}_*")
            )
        }
    }

    """
    mkdir fastqs
    sed 's+{FQS_DIR}+${params.fastq_screen_database_dir}+g' ${fastq_screen_config} > fastq_screen_config_interpolated.conf
    fastq_screen ${libraries.join(' ')} \
        --conf fastq_screen_config_interpolated.conf \
        --threads ${task.cpus} \
        --outdir fastqs \
        --aligner bowtie2
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc"
    label 'module_multiqc'

    input:
    path config, stageAs: 'multiqc_config', arity: '1'
    path fastqc_out, stageAs: 'fastqc', arity: '0..*'
    path fastq_screen_out, stageAs: 'fastq_screen', arity: '0..*'
    path cellranger_out, stageAs: 'cellranger_multi', arity: '0..*'

    output:
    path 'multiqc'

    script:
    """
    multiqc ${fastqc_out} ${fastq_screen_out} ${cellranger_out} --config ${config} -o multiqc
    """
}

workflow RUN_QUALITY_CONTROL {
    take:
    // samples
    samples

    // cellranger multi
    cellranger_multi

    // booleans
    skip_qc
    skip_fastqc
    skip_fastq_screen
    skip_multiqc

    // paths
    fastqs_config
    multiqc_config

    main:

    if (!skip_qc) {
        path_ch = samples.map { sample -> sample.libraries.collect { library -> library.path } }
        type_ch = samples.map { sample -> sample.libraries.collect { library -> library.type } }
        id_ch   = samples.map { sample -> sample.libraries.collect { library -> library.id } }
        name_ch = samples.map { sample -> sample.name }

        if (!skip_fastqc) {
            FASTQC (path_ch, id_ch, name_ch)
        }

        if (!skip_fastq_screen) {
            FASTQ_SCREEN (fastqs_config, path_ch, id_ch, type_ch, name_ch)
        }

        if (!skip_multiqc) {
            MULTIQC (
                multiqc_config,
                skip_fastqc ? [] : FASTQC.out.collect(),
                skip_fastq_screen ? [] : FASTQ_SCREEN.out.collect(),
                cellranger_multi.collect()
            )
        }
    }

    emit:
    fastqc  = ( skip_qc || skip_fastqc       ) ? [] : FASTQC.out
    fastqs  = ( skip_qc || skip_fastq_screen ) ? [] : FASTQ_SCREEN.out
    multiqc = ( skip_qc || skip_multiqc      ) ? [] : MULTIQC.out
}
