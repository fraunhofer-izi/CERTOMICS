#!/usr/bin/env nextflow

include { isNullFile ; getNullFile } from '../modules/local/functions'

process FASTQC {
    publishDir "${params.outdir}/fastqc/${task.tag}"
    label 'module_fastqc'
    tag "${sampleName}"

    input:
    tuple path(libraries, stageAs: 'library*/', arity: '1..*'), val(libraryIds), val(sampleName)

    output:
    path "fastqc"

    script:
    """
    mkdir "fastqc"
    fastq_files=\$(find -L ${libraries} -type f -regextype posix-extended -regex '.*\\.(fastq|fq)(\\.gz)?\$')
    echo "Running fastqc for \$fastq_files"
    fastqc -o "fastqc" \
        --noextract \
        --threads ${task.cpus} \
        \$fastq_files
    """
}

process FASTQ_SCREEN {
    publishDir "${params.outdir}/fastq_screen/${task.tag}"
    label 'module_fastq_screen'
    tag "${sampleName}"

    input:
    tuple path(gexLibrary, stageAs: 'library'), val(fastqId), val(sampleName)

    output:
    path "fastqs"

    script:
    """
    fastq_screen ${gexLibrary}/${fastqId}* \
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
    multiQcConfig

    main:
    if (!skipFastQc) {
        fastqcInput = samples.map { sample ->
            tuple(
                sample.libraries.collect { library -> library.path },
                sample.libraries.collect { library -> library.id },
                sample.name,
            )
        }

        FASTQC(fastqcInput)
    }

    if (!skipFastqScreen) {
        fastqScreenInput = samples.map { sample ->
            tuple(
                sample.libraries.find { library -> library.type == 'Gene Expression' }.path,
                sample.libraries.find { library -> library.type == 'Gene Expression' }.id,
                sample.name,
            )
        }

        FASTQ_SCREEN(fastqScreenInput)
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
