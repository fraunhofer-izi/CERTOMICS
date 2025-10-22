#!/usr/bin/env nextflow

process GET_GEX_SOURCE {
    input:
    path fa
    path gtf
    val fa_url
    val gtf_url

    output:
    path 'source.fa', emit: fa
    path 'source.gtf', emit: gtf

    script:
    """
    if [ -f "$fa" ]; then
        mv "$fa" "source.fa"
    else
        echo "Downloading FASTA"
        curl -sS "$fa_url" | zcat > "source.fa"
    fi

    if [ -f "$gtf" ]; then
        mv "$gtf" "source.gtf"
    else
        echo "Downloading GTF"
        curl -sS "$gtf_url" | zcat > "source.gtf"
    fi
    """
}

process BUILD_GEX_REFERENCE {
    publishDir "${params.outdir}/gene_expression_reference"
    label 'module_cellranger'
    label 'big_task'
   
    input:
    path src_fa, stageAs: 'src_fa', arity: '1'
    path src_gtf, stageAs: 'src_gtf', arity: '1'
    path car_fa, stageAs: 'car_fa'
    path car_gtf, stageAs: 'car_gtf'
    val ref_version

    output:
    path "gex_reference"

    shell:
    if (ref_version == '2020') {
        template 'build_reference_2020.sh'
    } else if (ref_version == '2024') {
        template 'build_reference_2024.sh'
    } else {
        // alternative: allow the user to use custom templates with ref_version?
        error 'invalid reference version'
    }
}

workflow HANDLE_GEX_REFERENCE {
    take:
    is_needed // boolean
    reference_version // integer

    // urls
    src_fa_url
    src_gtf_url

    // paths
    reference
    src_fa
    src_gtf
    car_fa
    car_gtf

    main:
    do_build = !reference.value
    if (is_needed.value && do_build) {
        do_get_fa  = !src_fa.value  || !file(src_fa.value).exists()
        do_get_gtf = !src_gtf.value || !file(src_gtf.value).exists()
        do_get_source = do_get_fa || do_get_gtf

        if (do_get_source) {
            GET_GEX_SOURCE (
                src_fa.value  ?: [],
                src_gtf.value ?: [],
                src_fa_url,
                src_gtf_url
            )
        }

        BUILD_GEX_REFERENCE (
            do_get_source ? GET_GEX_SOURCE.out.fa : src_fa,
            do_get_source ? GET_GEX_SOURCE.out.gtf : src_gtf,
            car_fa.value  ?: [],
            car_gtf.value ?: [],
            reference_version
        )
    }

    emit:
    is_needed ? (do_build ? BUILD_GEX_REFERENCE.out : reference) : []
}
