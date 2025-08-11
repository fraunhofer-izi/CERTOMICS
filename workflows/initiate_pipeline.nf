#!/usr/bin/env nextflow

include { paramsHelp         } from 'plugin/nf-validation'
include { paramsSummaryLog   } from 'plugin/nf-validation'
include { validateParameters } from 'plugin/nf-validation'

def verify_path (path, expectFile, expectDirectory) {
    if (!path) return 1

    def fileObject = file(path)
    if (!fileObject.exists()) return 2
    if (expectFile && !fileObject.isFile()) return 3
    if (expectDirectory && !fileObject.isDirectory()) return 4
    return 0
}

def verify_directory (path, parameter_name) {
    int verify_code = verify_path(path, false, true)
    if (verify_code == 0) {
        return path
    } else if (verify_code == 1) {
        error ("Falsy value ($path) path for $parameter_name")
    } else if (verify_code == 2) {
        error ("Path $path does not exist. Verify $parameter_name")
    } else if (verify_code == 3) {
        error ("Path $path is not a directory. Verify $parameter_name")
    } else {
        error ("Unexpected behaviour for path $path. Check $parameter_name")
    }
}

def verify_file (path, parameter_name) {
    int verify_code = verify_path(path, true, false)
    if (verify_code == 0) {
        return path
    } else if (verify_code == 1) {
        error ("Falsy value ($path) path for $parameter_name")
    } else if (verify_code == 2) {
        error ("Path $path does not exist. Verify $parameter_name")
    } else if (verify_code == 3) {
        error ("Path $path is not a file. Verify $parameter_name")
    } else {
        error ("Unexpected behaviour for path $path. Check $parameter_name")
    }
}

// helper method to verify gene_expression_source_* and gene_expression_source_*_url in combination
def verify_source_file_url (path, urlMap, urlKey, parameter_name) {
    if (path) {
        return [file: verify_file(path, parameter_name), url: null]
    } else if (urlMap) {
        log.info ("No file provided with $parameter_name. Downloading from URL instead.")
        // Do additional check for variable type?
        if (!urlMap.containsKey(urlKey))
            error ("parameter $parameter_name does not contain key $urlKey")
        return [file: null, url: urlMap[urlKey]]
    } else {
        error ("Both $parameter_name and ${parameter_name}_url are undefined.")
    }
}

workflow PARSE_PARAMETERS {
    main:
    samples = params.samples.collect { sampleMap -> Sample.create(sampleMap) }

    // decide which references are needed
    needs_gex  = samples.collect { sample -> sample.hasGeneExpressionLibrary() }.any()
    needs_vdj  = samples.collect { sample -> sample.hasVdjTLibrary() || sample.hasVdjBLibrary() }.any()
    needs_feat = samples.collect { sample -> sample.hasFeatureLibrary() }.any()

    // set default values for returned variables
    out_gene_expression_source_fa_url = null
    out_gene_expression_source_gtf_url = null
    out_gene_expression_reference = null
    out_gene_expression_source_fa = null
    out_gene_expression_source_gtf = null
    out_gene_expression_car_fa = null
    out_gene_expression_car_gtf = null
    out_vdj_reference = null
    out_feature_reference = null
    out_multiple_car_fa = null

    // overwrite defaults
    if (needs_gex) {
        // verify / check car files
        if (params.multiple_car_fa) { // if multiple_cars are given for kallisto process
            out_multiple_car_fa = verify_file(params.multiple_car_fa, 'multiple_car_fa')
        }
        if (params.gene_expression_car_fa && params.gene_expression_car_gtf) {
            // both car files provided
            out_gene_expression_car_fa = verify_file(
                params.gene_expression_car_fa,
                'gene_expression_car_fa')
            out_gene_expression_car_gtf = verify_file(
                params.gene_expression_car_gtf,
                'gene_expression_car_gtf')
        } else if (params.gene_expression_car_fa || params.gene_expression_car_gtf) {
            error ('Provided only one CAR file. Provide either both or none.')
        }

        if (params.gene_expression_reference) {
            // expecting to use prebuilt reference
            out_gene_expression_reference = verify_directory(
                params.gene_expression_reference,
                'gene_expression_reference')
        } else {
            // expecting to build reference at runtime
            // verify source files + urls
            source_fa_map = verify_source_file_url(
                params.gene_expression_source_fa,
                params.gene_expression_source_fa_url,
                params.gene_expression_reference_version,
                'gene_expression_source_fa')
            out_gene_expression_source_fa = source_fa_map.file
            out_gene_expression_source_fa_url = source_fa_map.url

            source_gtf_map = verify_source_file_url(
                params.gene_expression_source_gtf,
                params.gene_expression_source_gtf_url,
                params.gene_expression_reference_version,
                'gene_expression_source_gtf')
            out_gene_expression_source_gtf = source_gtf_map.file
            out_gene_expression_source_gtf_url = source_gtf_map.url
        }
    }

    if (needs_vdj) {
        if (params.vdj_reference) {
            out_vdj_reference = verify_directory(params.vdj_reference, 'vdj_reference')
        } else {
            error ('VDJ reference needed but not provided.')
        }
    }

    if (needs_feat) {
        if (params.feature_reference) {
            out_feature_reference = verify_file(params.feature_reference, 'feature_reference')
        } else {
            error ('Feature reference needed but not provided.')
        }
    }

    emit:
    samples = Channel.fromList(samples)
    needs_gex = needs_gex
    needs_vdj = needs_vdj
    needs_feat = needs_feat
    gex_version = params.gene_expression_reference_version
    gex_src_fa_url = out_gene_expression_source_fa_url
    gex_src_gtf_url = out_gene_expression_source_gtf_url
    gex_reference = out_gene_expression_reference
    gex_src_fa = out_gene_expression_source_fa
    gex_src_gtf = out_gene_expression_source_gtf
    gex_car_fa = out_gene_expression_car_fa
    gex_car_gtf = out_gene_expression_car_gtf
    vdj_reference = out_vdj_reference
    feat_reference = out_feature_reference
    multiple_car_fa = out_multiple_car_fa
}

// <copied with small modifications from  from https://github.com/nf-core/fetchngs/blob/1.12.0/subworkflows/nf-core/utils_nfvalidation_plugin/main.nf>
workflow RUN_NF_VALIDATION {
    take:
    print_help       // boolean: print help
    validate_params  // boolean: validate parameters
    schema_filename  // path: JSON schema file, null to use default value

    main:
    // Default values for strings
    workflow_command = 'nextflow run main.nf -profile <singularity/slurm/...> -params-file <parameter-file>'
    pre_help_text = 'LivingDrugOmics is an optimized single cell RNA sequencing omics pipeline for the purpose of high-resolution CAR T cell profiling. It supports the analysis of scRNA-seq data from common 10x Genomics single cell protocols, including gene expression, V(D)J repertoire and antibody/antigen recognition. Specific quality control metrics are incorporated for robust identification of CAR-positive cells.'
    post_help_text = "For more details, please visit our Wiki: (${workflow.manifest.docsUrl}) or contact the authors directly."

    // Print help message if needed
    if (print_help) {
        log.info([pre_help_text, paramsHelp(workflow_command, parameters_schema: schema_filename), post_help_text].join('\n\n'))
        System.exit(0)
    }

    if (validate_params) { validateParameters(parameters_schema: schema_filename) }
    log.info (paramsSummaryLog(workflow, parameters_schema: schema_filename))
}
// </>
