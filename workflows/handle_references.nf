#!/usr/bin/env nextflow

include { isNullFile } from '../modules/local/functions'

process PREPARE_SOURCE_FILES {
    input:
    path fa, stageAs: 'src/fa/*'
    path gtf, stageAs: 'src/gtf/*'

    output:
    path 'source.fa', emit: fa
    path 'source.gtf', emit: gtf

    script:
    """
    if ${fa.getExtension() == 'gz'}; then
        echo "Unzipping FASTA"
        zcat ${fa} > 'source.fa'
    else
        mv ${fa} 'source.fa'
    fi

    if ${gtf.getExtension() == 'gz'}; then
        echo "Unzipping GTF"
        zcat ${gtf} > 'source.gtf'
    else
        mv ${gtf} 'source.gtf'
    fi
    """
}

process BUILD_REFERENCE {
    publishDir "${params.outdir}/gene_expression_reference"
    label 'module_cellranger'
    label 'big_task'
   
    input:
    path sourceFasta, stageAs: 'src/fa/*', arity: '1'
    path sourceGtf, stageAs: 'src/gtf/*', arity: '1'
    path carFasta, stageAs: 'car/fa/*', arity: '1'
    path carGtf, stageAs: 'car/gtf/*', arity: '1'
    val referenceVersion

    output:
    path "gex_reference"

    script:
    doCar = !isNullFile(carFasta) && !isNullFile(carGtf)
    if (referenceVersion == '2020') {
        """
        bash build_reference_2020.sh \
            ${sourceFasta} \
            ${sourceGtf} \
            ${doCar ? carFasta : 0} \
            ${doCar ? carGtf : 0} \
            ${task.cpus ? task.cpus : 0} \
            ${task.memory ? task.memory.toGiga() : 0}
        """
    } else if (referenceVersion == '2024') {
        """
        bash build_reference_2024.sh \
            ${sourceFasta} \
            ${sourceGtf} \
            ${doCar ? carFasta : 0} \
            ${doCar ? carGtf : 0} \
            ${task.cpus ? task.cpus : 0} \
            ${task.memory ? task.memory.toGiga() : 0}
        """
    } else {
        error("'invalid reference version: ${referenceVersion}")
    }
}

workflow BUILD_CUSTOM_REFERENCE {
    take:
    referenceVersion
    sourceFasta
    sourceGtf
    carFasta
    carGtf

    main:
    PREPARE_SOURCE_FILES (
        sourceFasta,
        sourceGtf
    )

    BUILD_REFERENCE (
        PREPARE_SOURCE_FILES.out.fa,
        PREPARE_SOURCE_FILES.out.gtf,
        carFasta,
        carGtf,
        referenceVersion
    )

    emit:
    reference = BUILD_REFERENCE.out
}
