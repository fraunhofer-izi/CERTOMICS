# LivingDrugOmics

LivingDrugOmics is an optimized single cell RNA sequencing omics pipeline for the purpose of high-resolution CAR T cell profiling. It supports the analysis of scRNA-seq data from common 10x Genomics single cell protocols, including gene expression, V(D)J repertoire and antibody/antigen recognition. Specific quality control metrics are incorporated for robust identification of CAR-positive cells. 

For further information and installation guidelines go to our [Website](TODO)

# Contacts

If you have any questions or require further assistance, please contact us. 

Emails: 
- Christina.kuhn@izi.fraunhofer.de
- David.Schmidt@izi.fraunhofer.de


# LivingDrugOmics/
├── main.nf                      # Main Nextflow entrypoint
├── nextflow.config              # Config file
|── nextflow_schema.json         # Nextflow schema for validation
├── lib/
│   ├── Sample.groovy            # Custom Sample class
│   └── Library.groovy           # Custom Library class
├── workflows/
│   ├── initiate_pipeline.nf     # Contains PARSE_PARAMETERS, RUN_NF_VALIDATION
│   ├── handle_references.nf     # Contains HANDLE_GEX_REFERENCE etc.
│   ├── secondary_analysis.nf    # Contains CELLRANGER, KALLISTO etc.
│   └── quality_control.nf       # (optional) QC-related workflows
├── assets/                      # config files, helper scripts etc.
├── bin/                         # python and R scripts, quarto script
└── params_file.yml              # Parameter input
