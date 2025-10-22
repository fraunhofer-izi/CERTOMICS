#!/usr/bin/env nextflow

process PREPARE_SOURCE_FILES {
    input:
    path fa, stageAs: 'src/*'
    path gtf, stageAs: 'src/*'

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

process BUILD_GEX_REFERENCE {
    publishDir "${params.outdir}/gene_expression_reference"
    label 'module_cellranger'
    label 'big_task'
   
    input:
    path src_fa, stageAs: 'src/*', arity: '1'
    path src_gtf, stageAs: 'src/*', arity: '1'
    path car_fa, stageAs: 'car/fa/*'
    path car_gtf, stageAs: 'car/gtf/*'
    val ref_version

    output:
    path "gex_reference"

    script:
    def do_car = (
        car_fa.getSimpleName() != 'NO_FILE' &&
        car_gtf.getSimpleName() != 'NO_FILE'
    )

    if (ref_version == '2020') {
        """
        bash build_reference_2020.sh \
            ${src_fa} \
            ${src_gtf} \
            ${do_car ? car_fa : 0} \
            ${do_car ? car_gtf : 0} \
            ${task.cpus ? task.cpus : 0} \
            ${task.memory ? task.memory.toGiga() : 0}
        """
    } else if (ref_version == '2024') {
        """
        bash build_reference_2024.sh \
            ${src_fa} \
            ${src_gtf} \
            ${do_car ? car_fa : 0} \
            ${do_car ? car_gtf : 0} \
            ${task.cpus ? task.cpus : 0} \
            ${task.memory ? task.memory.toGiga() : 0}
        """
    } else {
        // alternative: allow the user to use custom templates with ref_version?
        error 'invalid reference version'
    }
}

workflow HANDLE_GEX_REFERENCE {
    take:
    // integer
    reference_version

    // paths
    src_fa
    src_gtf
    car_fa
    car_gtf

    main:
    PREPARE_SOURCE_FILES(
        src_fa,
        src_gtf
    )

    BUILD_GEX_REFERENCE(
        PREPARE_SOURCE_FILES.out.fa,
        PREPARE_SOURCE_FILES.out.gtf,
        car_fa,
        car_gtf,
        reference_version
    )

    emit:
    reference = BUILD_GEX_REFERENCE.out
}
