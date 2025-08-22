# CERTOMICS

CERTOMICS is a Nextflow-based pipeline offering enhanced **CERT**ainty in immunophenotyping and data interpretation, tailored for single-cell multi**OMICS** profiling of adoptive cellular immunotherapies. The pipeline standardises processing 10x Genomics single-cell multiomics data and integrates CAR-specific identification and quality control. Additionally, a curated repository of CAR construct sequences and annotation data is provided, serving as a resource to support the analysis and development of CAR T cell therapies.

For further information and **installation guidelines** go to our [Website](https://fraunhofer-izi.github.io/Living-Drugs-Wiki/Home/Pipelines/CERTOMICS-Pipeline/)

<img alt="Abb_CAR_T_cell_profiling_V2"
     src="https://github.com/user-attachments/assets/b49a424e-befc-47e0-a568-9c2523691c0c"
     width="50%" />


# Contacts

If you have any questions or require further assistance, please contact us. 

Emails: 
- Christina.kuhn@izi.fraunhofer.de
- David.Schmidt@izi.fraunhofer.de


# CERTOMICS/
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
