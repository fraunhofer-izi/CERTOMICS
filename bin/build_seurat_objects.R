#!/usr/bin/env Rscript

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries needed for script
.cran_packages = c("Seurat", "SeuratObject", "yaml", "dplyr", "doParallel", "parallel", "data.table", "Matrix", "stringr","scGate", "remotes","SoupX")
.bioc_packages = c("BiocParallel", "SingleCellExperiment", "scDblFinder", "scds","scRepertoire", "UCell", "clustifyr")
.github_packages = c("carmonalab/ProjecTILs") # Add ProjecTILs

# Installation of libraries (for local script testing)
.inst = .cran_packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}
.inst <- .bioc_packages %in% installed.packages()
if (any(!.inst)) {
  library(BiocManager)
  BiocManager::install(.bioc_packages[!.inst], ask = T)
}
.inst = sapply(.github_packages, function(pkg) {
  pkg_name <- strsplit(pkg, "/")[[1]][2] # Extract package name
  pkg_name %in% installed.packages()
})
if (any(!.inst)) {
  remotes::install_github(.github_packages[!.inst])
}
list.of.packages = c(.cran_packages, .bioc_packages, "ProjecTILs")

## Loading library
for (pack in list.of.packages) {
  suppressMessages(library(
    pack,
    quietly = TRUE,
    verbose = FALSE,
    character.only = TRUE
  ))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Logging Function for Debugging
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
log_file <- "build_seurat_object_log.txt"  # Define log file path
log_message <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_message <- paste0("[", timestamp, "] ", message, "\n")
  cat(full_message, file = log_file, append = TRUE)
  cat(full_message)  # Print to stdout for Nextflow logs
  flush.console()
}
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Validate Command-line Arguments
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
args = commandArgs(trailingOnly = TRUE)

validate_args <- function(args) {
  if (length(args) < 9) stop("Please provide at least eight arguments.")
  if (!file.exists(args[1])) stop("Helper script not found: ", args[1])

  filt <- unlist(strsplit(args[2], ","))
  raw  <- unlist(strsplit(args[3], ","))
  smps <- unlist(strsplit(args[6], ","))

  if (!all(file.exists(filt))) stop("One or more filtered 10x dirs do not exist.")
  if (!all(file.exists(raw)))  stop("One or more RAW 10x dirs do not exist.")
  
  if (length(filt) != length(smps))
    stop("The number of filtered 10x dirs does not match the number of samples.")
  if (length(raw)  != length(smps))
    stop("The number of RAW 10x dirs does not match the number of samples.")

  if (args[9] == "none") {
    message("No CAR are given")
    log_message("INFO: No CAR reference was given.")
  } else if (!file.exists(args[9])) {
    stop("GTF file not found: ", args[9])
  }
}
validate_args(args)
# log_message(paste("Arguments received:\n", paste(seq_along(args), args, sep = ": ", collapse = "\n")))

tryCatch({
  helper_script <- args[1]
  source(helper_script)
}, error = function(e) {
  stop(e)  # Stop execution after logging
})
# Parse arguments with error handling
tryCatch({
  cellranger_dirs_filt <- unlist(strsplit(args[2], ","))  # filtered per-sample MEX
  cellranger_dirs_raw  <- unlist(strsplit(args[3], ","))  # raw (unfiltered) MEX
  cellranger_vdjt      <- unlist(strsplit(args[4], ","))
  cellranger_vdjb      <- unlist(strsplit(args[5], ","))
  cellranger_samples   <- unlist(strsplit(args[6], ","))
  seurat_output        <- args[7]
  type_bits            <- args[8]
  gtf_path             <- args[9]
  scGate_model         <- args[10]

  flags <- parse_binary_string(type_bits)
  #Position 1: GEX (Gene Expression)
  #Position 2: VDJT (VDJ-T)
  #Position 3: VDJB (VDJ-B)
  #Position 4: AntibodyCapture = CITE-seq

  # Log the detected libraries
  log_message(paste("GEX:", flags$GEX))
  log_message(paste("VDJ-T:", flags$VDJT))
  log_message(paste("VDJ-B:", flags$VDJB))
  log_message(paste("AntibodyCapture:", flags$AntibodyCapture))

  # If neither VDJ-T nor VDJ-B are available, warn but do not stop
  if (!flags$VDJT && !flags$VDJB) {
    log_message("INFO: No VDJ-T or VDJ-B data provided. Proceeding with Gene Expression only.")
  } else if (flags$VDJT && flags$VDJB) {
    log_message("INFO: VDJ-T and VDJ-B data provided.")
  } else if  (flags$VDJT) {
    log_message("INFO: VDJ-T data provided.")
  } else if  (flags$VDJB) {
    log_message("INFO: VDJ-B data provided.")
  }
  log_message("Arguments successfully parsed.")

}, error = function(e) {
  stop(e)  # Stop execution after logging
})

# Check that the number of samples matches the directories
if (length(cellranger_dirs_filt) != length(cellranger_samples)) {
  stop("The number of filtered directories does not match the number of samples.")
}
if (length(cellranger_dirs_raw) != length(cellranger_samples)) {
  stop("The number of RAW directories does not match the number of samples.")
}
if (gtf_path == "none") {
  car_construct <- NULL
  message("No CAR are given")
} else {
  gtf <- read_gtf_file(gtf_path)
  if (is.null(gtf)) {
    stop("Failed to read the GTF file")
  }
  car_construct <- gtf[1, 1]
}

names(cellranger_dirs_filt) <- cellranger_samples
names(cellranger_dirs_raw)  <- cellranger_samples
if (any(cellranger_vdjt != "none")) names(cellranger_vdjt) <- cellranger_samples
if (any(cellranger_vdjb != "none")) names(cellranger_vdjb) <- cellranger_samples

# Validate scGate_model
allowed_models <- c("PBMC", "TME_HiRes")
if (!(scGate_model %in% allowed_models)) {
  stop(paste0(
    "Invalid scGate_model: ", scGate_model,
    ". Must be one of: ", paste(allowed_models, collapse = ", ")
  ))
} else {
  log_message(paste("Using scGate_model:", scGate_model))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load Rawcounts and create a merged Seurat object accoring to library used
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if (!any(unlist(flags))){
  stop("No libraries were provided. Check binary flag.")
}

# Load and Process Data
se.meta <- read_in_10x_parallel(cellranger_dirs_filt, flags$AntibodyCapture, cellranger_dirs_raw)
se.meta <- normalization_seurat(se.meta, flags$AntibodyCapture)
se.meta <- data_quality_mito_ribo_complexity(se.meta)
se.meta <- add_low_quality_info(se.meta)
se.meta <- cell_cycle_scores(se.meta)
se.meta <- annotation_exr_and_scgate(se.meta, car_construct, scGate_model)

if (flags$VDJT && flags$VDJB) {
  se.meta <- assign_both_vdj(obj = se.meta, batch_vdj = cellranger_vdjt, batch_vdb = cellranger_vdjb, TRUE)
} else if (flags$VDJT) {
  se.meta <- assign_vdj(se.meta, "vdj_t", cellranger_vdjt, TRUE)
} else if (flags$VDJB) {
  se.meta <- assign_vdj(se.meta, "vdj_b", cellranger_vdjb, TRUE)
}

DefaultAssay(se.meta) <- "RNA"
Idents(se.meta) <- "orig.ident"
saveRDS(se.meta, file = seurat_output)