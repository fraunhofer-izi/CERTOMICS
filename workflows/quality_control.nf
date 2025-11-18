#!/usr/bin/env nextflow

include { isNullFile; getNullFile } from '../modules/local/functions'

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

process BUILD_FASTQ_SCREEN_CONFIG {
    input:
    path dbDir, stageAs: 'databases/*', arity: '1..*'
    val dbName
    val dbFile
    
    output:
    path "config"

    script:
    configText = [dbDir, dbName, dbFile].transpose().collect { dir, name, db ->
        "DATABASE\t${name}\t${dir.resolve(db)}"
    }

    """
    mkdir config
    mv databases config
    echo "${configText.join('\n')}" > config/fastq_screen.conf
    """
}

process FASTQ_SCREEN {
    publishDir "${params.outdir}/fastq_screen/${task.tag}"
    label 'module_fastq_screen'
    tag "${sampleName}"

    input:
    path config, stageAs: 'configDir', arity: '1'
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
    fastq_screen ${libraries.join(' ')} \
        --conf ${config}/*.conf \
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
    path fastQcReports, stageAs: 'reports/fastqc?/*', arity: '0..*'
    path fastqScreenReports, stageAs: 'reports/fastq_screen?/*', arity: '0..*'
    path cellRangerReports, stageAs: 'reports/cellranger?/*', arity: '0..*'

    output:
    path 'multiqc'

    script:
    """
    multiqc ${cellRangerReports} ${fastQcReports} ${fastqScreenReports} --config ${multiQcConfig} -o multiqc
    """
}

workflow QUALITY_CONTROL {
    take:
    samples
    cellrangerReports
    skipFastQc
    skipFastqScreen
    skipMultiQc
    fastqScreenDatabases
    multiQcConfig

    main:
    pathCh = samples.map { sample -> sample.libraries.collect { library -> library.path } }
    typeCh = samples.map { sample -> sample.libraries.collect { library -> library.type } }
    idCh = samples.map { sample -> sample.libraries.collect { library -> library.id } }
    nameCh = samples.map { sample -> sample.name }

    if (!skipFastQc) {
        FASTQC(pathCh, idCh, nameCh)
    }

    if (!skipFastqScreen) {
        BUILD_FASTQ_SCREEN_CONFIG(
            fastqScreenDatabases.collect { name, dir, db -> dir },
            fastqScreenDatabases.collect { name, dir, db -> name },
            fastqScreenDatabases.collect { name, dir, db -> db },
        )

        FASTQ_SCREEN(
            BUILD_FASTQ_SCREEN_CONFIG.out,
            pathCh,
            idCh,
            typeCh,
            nameCh
        )
    }

    if (!skipMultiQc) {
        MULTIQC(
            multiQcConfig,
            skipFastQc ? getNullFile() : FASTQC.out.collect(),
            skipFastqScreen ? getNullFile() : FASTQ_SCREEN.out.collect(),
            cellrangerReports.collect(),
        )
    }

    emit:
    fastqc = skipFastQc ? channel.empty() : FASTQC.out
    fastqs = skipFastqScreen ? channel.empty() : FASTQ_SCREEN.out
    multiqc = skipMultiQc ? channel.empty() : MULTIQC.out
}
