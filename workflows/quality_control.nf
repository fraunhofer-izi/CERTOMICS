#!/usr/bin/env nextflow

process FASTQC {
    publishDir "${params.outdir}/fastqc/${task.tag}"
    label 'module_fastqc'
    tag "${sampleName}"

    input:
    path fastqPaths, stageAs: 'fastq', arity: '1..*'
    val fastqIds
    val sampleName

    output:
    path "fastqc"

    script:
    libraries = []
    fastqPaths.eachWithIndex { directory, index ->
        libraries.add(
            directory.resolve("${fastqIds[index]}_*")
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
    tag "${sampleName}"

    input:
    path fastqScreenConfig, arity: '1'
    path fastqPaths, stageAs: 'fastq', arity: '1..*'
    val fastqIds
    val fastqTypes
    val sampleName

    output:
    path "fastqs"

    script:
    libraries = []
    fastqPaths.eachWithIndex { directory, index ->
        if ('Gene Expression'.equals(fastqTypes[index])) {
            libraries.add(
                directory.resolve("${fastqIds[index]}_*")
            )
        }
    }

    """
    mkdir fastqs
    sed 's+{FQS_DIR}+${params.fastq_screen_database_dir}+g' ${fastqScreenConfig} > fastq_screen_config_interpolated.conf
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
    path multiQcConfig, stageAs: 'multiqc_config', arity: '1'
    path fastQcReports, stageAs: 'fastqc/*', arity: '0..*'
    path fastqScreenReports, stageAs: 'fastq_screen/*', arity: '0..*'
    
    output:
    path 'multiqc'

    script:
    """
    multiqc ${fastQcReports} ${fastqScreenReports} --config ${multiQcConfig} -o multiqc
    """
}

workflow QUALITY_CONTROL {
    take:
    // samples
    samples

    // booleans
    skipFastQc
    skipFastqScreen
    skipMultiQc

    // paths
    fastqScreenConfig
    multiQcConfig

    main:
    pathCh = samples.map { sample -> sample.libraries.collect { library -> library.path } }
    typeCh = samples.map { sample -> sample.libraries.collect { library -> library.type } }
    idCh   = samples.map { sample -> sample.libraries.collect { library -> library.id } }
    nameCh = samples.map { sample -> sample.name }
    
    if (!skipFastQc) {
        FASTQC (pathCh, idCh, nameCh)
    }

    if (!skipFastqScreen) {
        FASTQ_SCREEN (fastqScreenConfig, pathCh, idCh, typeCh, nameCh)
    }

    if (!skipMultiQc) {
        MULTIQC (
            multiQcConfig,
            skipFastQc ? Channel.empty() : FASTQC.out.collect(),
            skipFastqScreen ? Channel.empty() : FASTQ_SCREEN.out.collect()
        )
    }

    emit:
    fastqc  = skipFastQc ? Channel.empty() : FASTQC.out
    fastqs  = skipFastqScreen ? Channel.empty() : FASTQ_SCREEN.out
    multiqc = skipMultiQc ? Channel.empty() : MULTIQC.out
}
