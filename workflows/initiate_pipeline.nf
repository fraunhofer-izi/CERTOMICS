#!/usr/bin/env nextflow

include { paramsHelp         } from 'plugin/nf-validation'
include { paramsSummaryLog   } from 'plugin/nf-validation'
include { validateParameters } from 'plugin/nf-validation'



// <copied with small modifications from  from https://github.com/nf-core/fetchngs/blob/1.12.0/subworkflows/nf-core/utils_nfvalidation_plugin/main.nf>
workflow NF_VALIDATION {
    take:
    print_help       // boolean: print help
    validate_params  // boolean: validate parameters
    schema_filename  // path: JSON schema file, null to use default value

    main:
    // Default values for strings
    workflow_command = 'nextflow run main.nf -profile <singularity/slurm/...> -params-file <parameter-file>'
    pre_help_text = 'CERTOMICS is an Nextfow-based pipeline offering enhanced CERTainty in immunophenotyping and data interpretation, tailored for single-cell multiOMICS profling of adoptive cellular immunotherapies. It supports the analysis of scRNA-seq data from common 10x Genomics single cell protocols, including gene expression, V(D)J repertoire and antibody/antigen recognition. Specific quality control metrics are incorporated for robust identification of CAR-positive cells.'
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
