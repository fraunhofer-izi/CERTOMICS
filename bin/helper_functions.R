# Seurat-based Single-Cell RNA-seq Processing Pipeline
################################################################################
# This script contains a set of functions to process and analyze single-cell RNA-seq data
# using the Seurat package and related tools. It covers data loading, quality control,
# cell filtering, normalization, cell type annotation, VDJ assignment, and cell cycle analysis.

################################################################################
# Function Descriptions:

### `label_cells_rm(obj)`
# - Labels cells for removal based on:
#   - `nFeature_RNA`, `nCount_RNA` (gene & UMI counts)
#   - Mitochondrial gene percentage (`Perc_of_mito_genes`)
#   - Complexity (`log10GenesPerUMI`)
#   - Doublet classification from scDblFinder
# - Outputs a summary table of cell retention per sample.

###  `Split_Object(object, split.by = "orig.ident", threads = 5)`
# - Splits a Seurat object by a metadata column (`orig.ident` by default).
# - Uses `BiocParallel::bplapply()` for parallel processing.
# - Returns a list of Seurat objects, one per sample or batch.

###  `count_cells_per_sample(obj, count.base = NULL, col.name = NULL)`
# - Counts the number of cells per sample (`orig.ident`).
# - If `count.base` is provided, it updates an existing tracking table.
# - Useful for tracking cell retention before and after filtering.

###  `assign_vdj(obj, vdj = "vdj_t", batch = NULL, present.bool = TRUE)`
# - Assigns VDJ sequences (TCR/BCR) to cells using `scRepertoire` and Cell Ranger outputs.
# - Combines and integrates clonotype data into the Seurat object.
# - Adds VDJ-T and VDJ-B availability to metadata.

### `parse_binary_string(binary_string)`
# - Converts a binary string into a boolean list.
# - Used to indicate the presence of different data modalities:
#   - `GEX`, `VDJ-T`, `VDJ-B`, `AntibodyCapture`.

### `read_in_10x_parallel(fltrd.dirs, ADT, raw.dirs)`
# - Reads 10x Genomics filtered feature barcode matrices into Seurat.
# - Supports Antibody-Derived Tags (ADT) if present.
# - Uses `scDblFinder()` to compute doublet scores.

### `normalization_seurat(se.meta, ADT)`
# - Normalizes RNA and ADT assays:
#   - RNA → LogNormalization
#   - ADT → Centered Log-Ratio (CLR) normalization
# - Ensures ADT data is stored as a sparse matrix.

### `annotation_exr_and_scgate(se.meta, car_construct)`
# - Annotates CAR-T cells based on raw expression (`car_by_expr()`).
# - Uses scGate to classify cell types.
# - Integrates CAR expression with CD4/CD8 annotations.

### `data_quality_mito_ribo_complexity(se.meta)`
# - Adds quality control metrics:
#   - Mitochondrial gene percentage (`MT-`)
#   - Ribosomal gene percentage (`RPL|RPS`)
#   - Log10 genes per UMI (complexity metric)

### `add_low_quality_info(se.meta)`
# - Removes low-quality cells based on QC thresholds.
# - Tracks cell counts before and after filtering.

### `cell_cycle_scores(se.meta)`
# - Scores cell cycle phases (`S`, `G2M`, `G1M`) using Seurat's `CellCycleScoring()`.
# - Uses clustering-based GSEA analysis for cell cycle annotation.

### `estimate_cc(se)`
# - Estimates cell cycle phase using gene set enrichment analysis (GSEA).
# - Uses the ProjecTILs package and clustifyr to assign phases.

### `car_by_expr(obj, car_construct)`
# - Annotates CAR-T cells based on gene expression.
# - Checks for CAR construct presence in RNA counts.

### `assign_both_vdj(obj, batch_vdj, batch_vdb, present.bool = TRUE)`
# - Assigns both VDJ-T and VDJ-B to cells.
# - Uses `scRepertoire::combineTCR()` and `combineBCR()`.
# - Integrates VDJ metadata into the Seurat object.

### `read_gtf_file(gtf_path)`
# - Reads a GTF annotation file and validates its structure.
# - Ensures correct column format (9 columns).
# - Handles file format errors gracefully.

### `estimate_cc(se)`
# - Runs cell cycle estimation using PCA-based clustering.
# - Uses clustifyr GSEA for phase assignment.

############### Setup #########################################################
# Filter cutoffs for low cell quality
nFeature_low_cutoff = 250
nFeature_high_cutoff = 8000
nCount_low_cutoff = 1000
nCount_high_cutoff = 100000
mt_cutoff = 15
complx_cutoff = 0.8

#Marker genes for cell cycle
# DOI: 10.1016/j.cell.2018.09.006 -> DOI: 10.1126/science.aad0501
s.genes = c("ATAD2", "BLM", "BRIP1", "CASP8AP2", "CCNE2", "CDC45", "CDC6", "CDCA7",
            "CHAF1B", "CLSPN", "DSCC1", "DTL", "E2F8", "EXO1", "FEN1", "GINS2", "GMNN", "HELLS",
            "MCM2", "MCM4", "MCM5", "MCM6", "MLF1IP", "MSH2", "NASP", "PCNA", "POLA1", "POLD3",
            "PRIM1", "RAD51", "RAD51AP1", "RFC2", "RPA2", "RRM1", "RRM2", "SLBP", "TIPIN",
            "TYMS", "UBR7", "UHRF1", "UNG", "USP1", "WDR76")

g2m.genes = c("ANLN", "ANP32E", "AURKA", "AURKB", "BIRC5", "BUB1", "CBX5", "CCNB2",
              "CDC20", "CDC25C", "CDCA2", "CDCA3", "CDCA8", "CDK1", "CENPA", "CENPE", "CENPF",
              "CKAP2", "CKAP2L", "CKAP5", "CKS1B", "CKS2", "CTCF", "DLGAP5", "ECT2", "FAM64A",
              "G2E3", "GAS2L3", "GTSE1", "HJURP", "HMGB2", "HMMR", "HN1", "KIF11", "KIF20B",
              "KIF23", "KIF2C", "LBR", "MKI67", "NCAPD2", "NDC80", "NEK2", "NUF2", "NUSAP1",
              "PSRC1", "RANGAP1", "SMC4", "TACC3", "TMPO", "TOP2A", "TPX2", "TTK", "TUBB4B", "UBE2C")

############### R Functions ####################################################

#label cells that will be removed by quality control cutoffs
label_cells_rm = function(obj) {
  obj@meta.data = obj@meta.data %>% mutate(
    KEEP_CELL = case_when(
      (nFeature_RNA < nFeature_low_cutoff) | (nFeature_RNA > nFeature_high_cutoff) |
      (nCount_RNA < nCount_low_cutoff) | (nCount_RNA > nCount_high_cutoff) |
      (Perc_of_mito_genes > mt_cutoff)  | (log10GenesPerUMI < complx_cutoff) |
      (scDblFinder_class == "doublet") ~ FALSE,  # Mark doublets as FALSE
      TRUE ~ TRUE
    )
  )
  table_output <- capture.output(table(obj$orig.ident, obj$KEEP_CELL))
  cat("Table of obj$orig.ident vs obj$KEEP_CELL:\n", paste(table_output, collapse = "\n"), "\n")
  flush.console()
  obj
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Read the GTF file with error handling
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
read_gtf_file <- function(gtf_path) {
  tryCatch({
    # Attempt to read the GTF file
    gtf <- read.table(
      gtf_path,
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE,
      comment.char = "#"
    )
    # Validate the structure of the file
    if (ncol(gtf) != 9) {
      stop("The GTF file does not have exactly 9 columns in all rows.")
    }
    return(gtf)
  }, error = function(e) {
    # Handle errors gracefully
    message("An error occurred while reading the GTF file: ", e$message)
    message("Please check the file format and ensure it has 9 tab-separated columns.")
    NULL  # Return NULL to indicate failure
  })
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Function to convert binary string to booleans for each dataset
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
parse_binary_string <- function(binary_string) {
  # Ensure the binary string is at least 4 characters long
  # Pad with "0"s to the right if it is shorter
  binary_string <- str_pad(binary_string, 4, pad = "0", side = "right")

  # Return a list of booleans for each position in the binary string
  # Each boolean indicates whether the corresponding dataset type is present
  list(
    GEX = substr(binary_string, 1, 1) == "1",
    VDJT = substr(binary_string, 2, 2) == "1",
    VDJB = substr(binary_string, 3, 3) == "1",
    AntibodyCapture = substr(binary_string, 4, 4) == "1"
  )
}
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Split Seurat object (BiocParallel)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Split_Object = function(object, split.by = "orig.ident", threads = 5) {

  library(parallel)
  library(BiocParallel)

  groupings <- FetchData(object = object, vars = split.by)[, 1]
  groupings <- unique(x = as.character(x = groupings))
  names(groupings) = groupings

  if (is.null(threads)) {
    bpparam = BiocParallel::MulticoreParam(workers = length(groupings))
  } else {
    bpparam = BiocParallel::MulticoreParam(workers = threads)
  }

  obj.list = BiocParallel::bplapply(groupings, function(grp) {
    cells <- which(x = object[[split.by, drop = TRUE]] == grp)
    cells <- colnames(x = object)[cells]
    se = subset(x = object, cells = cells)
    se@meta.data = droplevels(se@meta.data)
    se
  }, BPPARAM = bpparam)

  return(obj.list)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Track number of cell (pre-processing
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
count_cells_per_sample = function(obj = NULL, count.base = NULL, col.name = NULL){
  l = lapply(obj, function(x) {
  x@meta.data %>% dplyr::select(orig.ident)
  })
  df = do.call("rbind", l) %>%
    dplyr::count(orig.ident) %>%
    dplyr::rename(beforeFiltering = n)  # Rename count column to "beforeFiltering"
  if (is.null(count.base)) {
    return(df)  # Return df with renamed column
  } else {
    count.base[[col.name]] = df$beforeFiltering[match(count.base$orig.ident, df$orig.ident)]
    return(count.base)
  }
}
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Set CAR by expression
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
car_by_expr = function(obj, car_construct){
  if(car_construct %in% rownames(GetAssayData(obj, layer = c("counts"), assay = "RNA"))){
    car.ftr = FetchData(obj, c(car_construct), layer = "counts")
    obj$CAR_BY_EXPRS = as.factor(car.ftr[[car_construct]] > 0)
  } else {
    obj$CAR_BY_EXPRS = FALSE
    obj$CAR_BY_EXPRS = as.factor(obj$CAR_BY_EXPRS)
  }

  obj
}


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Assign VDJ
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
assign_both_vdj = function(obj, batch_vdj = NULL, batch_vdb = NULL, present.bool = TRUE){
  log_message("PROCESS: assign_vdj")
  # Read filtered contig annotation files and filter out empty or invalid files
  # Read filtered contig annotation files and filter out empty or invalid files
  contig_list_t <- lapply(batch_vdj, function(x) {
    tryCatch(read.csv(x), error=function(e) NULL)
  })
  contig_list_t = contig_list_t[lengths(contig_list_t) != 0]

  contig_list_b <- lapply(batch_vdb, function(x) {
    tryCatch(read.csv(x), error=function(e) NULL)
  })
  contig_list_b = contig_list_b[lengths(contig_list_b) != 0]

  # Initialize combined list
  combined <- list()

  # Check if contig_list_t is empty
  if (length(contig_list_t) > 0) {
    combined_t <- scRepertoire::combineTCR(
      contig_list_t,
      samples = paste0(names(contig_list_t))
    )
    combined <- c(combined, combined_t)
  } else {
    combined_t <- NULL
  }

  # Check if contig_list_b is empty
  if (length(contig_list_b) > 0) {
    combined_b <- scRepertoire::combineBCR(
      contig_list_b,
      samples = paste0(names(contig_list_b))
    )
    combined <- c(combined, combined_b)
  } else {
    combined_b <- NULL
  }

  # If both are empty, skip processing
  if (length(combined) == 0) {
    message("Both TCR and BCR lists are empty. Skipping processing.")
  } else {
    obj <- scRepertoire::combineExpression(
      combined, obj,
      cloneCall = "strict",
      group.by = "sample"
    )
    obj$cloneSize = droplevels(obj$cloneSize)
  }

  # Convert combined lists to data frames only if they are not empty
  if (!is.null(combined_t)) {
    combined_t = data.table::rbindlist(combined_t) %>% data.frame()
  }

  if (!is.null(combined_b)) {
    combined_b = data.table::rbindlist(combined_b) %>% data.frame()
  }

  # Annotate Seurat metadata if present.bool is TRUE
  if (present.bool == TRUE) {
    if (!is.null(combined_t)) {
      obj$VDJ_T_AVAIL = combined_t$CTnt[match(rownames(obj@meta.data), combined_t$barcode)]
      obj$VDJ_T_AVAIL = ifelse(is.na(obj$VDJ_T_AVAIL), FALSE, TRUE)
    } else {
      obj$VDJ_T_AVAIL <- FALSE
    }

    if (!is.null(combined_b)) {
      obj$VDJ_B_AVAIL = combined_b$CTnt[match(rownames(obj@meta.data), combined_b$barcode)]
      obj$VDJ_B_AVAIL = ifelse(is.na(obj$VDJ_B_AVAIL), FALSE, TRUE)
    } else {
      obj$VDJ_B_AVAIL <- FALSE
    }
  }
  log_message("DONE: assign_vdj")
  obj

}

assign_vdj = function(obj, vdj = "vdj_t", batch = NULL, present.bool = TRUE){
  log_message("PROCESS: assign_vdj")

  # Ensure the provided 'vdj' argument is valid ("vdj_t" or "vdj_b")
  stopifnot(vdj %in% c("vdj_t", "vdj_b"))

  fltrd.vdj = batch

  # Read filtered contig annotation files and filter out empty or invalid files
  contig_list <- lapply(fltrd.vdj, function(x) {
    tryCatch(read.csv(x), error=function(e) NULL)
  })
  contig_list = contig_list[lengths(contig_list) != 0]

  # Combine TCR or BCR data depending on the 'vdj' type
  if(vdj == "vdj_t"){
    cat("vdj-t library given")
    combined <- scRepertoire::combineTCR(
      contig_list,
      samples = paste0(names(contig_list))
    )
  } else {
    cat("vdj-b library given")
    combined <- scRepertoire::combineBCR(
      contig_list,
      samples = paste0(names(contig_list))
    )
  }
  combined = data.table::rbindlist(combined) %>% data.frame()

  # Annotate Seurat metadata if present.bool is TRUE
  if(present.bool == TRUE){
    if(vdj == "vdj_t"){
      # Mark availability of VDJ-T data in the Seurat metadata
      obj$VDJ_T_AVAIL = combined$CTnt[match(rownames(obj@meta.data), combined$barcode)]
      obj$VDJ_T_AVAIL = ifelse(is.na(obj$VDJ_T_AVAIL), FALSE, TRUE)
    } else {
      # Mark availability of VDJ-B data in the Seurat metadata
      obj$VDJ_B_AVAIL = combined$CTnt[match(rownames(obj@meta.data), combined$barcode)]
      obj$VDJ_B_AVAIL = ifelse(is.na(obj$VDJ_B_AVAIL), FALSE, TRUE)
    }
  }

  combined$sample = obj$orig.ident[match(combined$barcode, rownames(obj@meta.data))]
  combined = combined[!is.na(combined$sample), ]
  combined = split(combined, combined$sample)

  obj <- scRepertoire::combineExpression(
    combined, obj,
    cloneCall = "strict",
    group.by = "sample"
  )
  obj$cloneSize = droplevels(obj$cloneSize)
  log_message("DONE: assign_vdj")
  obj

}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Read in 10x data GEX (and ADT)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
read_in_10x_parallel <- function(fltrd.dirs, ADT, raw.dirs){
  log_message("PROCESS: Read_in_10x_parallel.")
  bpparam = BiocParallel::MulticoreParam(workers = 16)
  seurat.l = BiocParallel::bplapply(names(fltrd.dirs), function(i) {
    id = i
    # Read 10X data from the specified directory for the current dataset (RNA and ADT counts)
    log_message(paste("Processing 10X data for:", id))
    data = Seurat::Read10X(data.dir = fltrd.dirs[names(fltrd.dirs) == id], gene.column = 2)
    # Create a Seurat object using the gene expression (RNA) counts, setting the project name as the dataset ID
    if (!ADT) {
      seu.obj = Seurat::CreateSeuratObject(counts = data, project = id) # if no adt
    } else {
      seu.obj = Seurat::CreateSeuratObject(counts = data[[1]], project = id) # if adt
      # Add ADT (Antibody-Derived Tag) counts as a new assay to the Seurat object
      seu.obj[['ADT']] = Seurat::CreateAssayObject(counts = data[[2]]) #get adt counts
    }
    # Rename cells by removing the "multi_" prefix from their names
    seu.obj = Seurat::RenameCells(seu.obj, new.names = gsub("multi_", "", colnames(seu.obj)))
    seu.obj@meta.data$orig.ident = id

    # Add doublet removel scores based on scDblFinder
    log_message("PROCESS: scDblFinder")
    sce = scDblFinder(Seurat::GetAssayData(seu.obj, assay="RNA", layer="counts"), verbose=TRUE)
    df = data.frame(sce@colData) %>% dplyr::select(scDblFinder.score, scDblFinder.class)
    colnames(df) = c("scDblFinder_score", "scDblFinder_class")
    seu.obj = Seurat::AddMetaData(seu.obj, df)

    # Perform soupX
    log_message("PROCESS: SoupX")
    seu.obj$soup_group <- get_soup_groups(seu.obj)
    raw_dir <- raw.dirs[[id]]
    if (is.null(raw_dir) || !dir.exists(raw_dir)) {
      stop(sprintf("Raw directory missing for %s: %s", id, as.character(raw_dir)))
    }
    seu.obj <- make_soup(seu.obj, raw_dir = raw_dir, ADT)

    seu.obj
  }, BPPARAM = bpparam)
  
  log_message("PROCESS: Combine Seurat.l")
  # If more than one Seurat object was created, merge them into a single Seurat object
  if(length(seurat.l) > 1) {
    se.meta = merge(
      seurat.l[[1]], y = seurat.l[2:length(seurat.l)],
      add.cell.ids = names(seurat.l))
  } else {
    se.meta = seurat.l[[1]]
  }

  # Convert the original identity metadata to a factor for consistency
  se.meta@meta.data$orig.ident = factor(se.meta@meta.data$orig.ident)
  # Join multiple layers (if any) in the RNA assay for the Seurat object
  se.meta[["RNA"]] <- SeuratObject::JoinLayers(se.meta[["RNA"]])

  DefaultAssay(se.meta) = "RNA"
  Idents(se.meta) = "orig.ident"

  if (ADT) {  ####### Process ADT
    # Switch the default assay to ADT (Antibody-Derived Tags) for further processing
    DefaultAssay(se.meta) = "ADT"
    # Clean up the feature names in the ADT assay by removing the "totalseqC-" prefix and converting to uppercase
    adt.ftrs = gsub("totalseqC-", "", rownames(se.meta))
    rownames(se.meta@assays$ADT@counts) = toupper(adt.ftrs)
    # Set the raw counts as the data for the ADT assay (i.e., no transformations yet)
    se.meta@assays$ADT@data = se.meta@assays$ADT@counts
    rownames(se.meta@assays$ADT@meta.features) = rownames(se.meta@assays$ADT@counts)
    # Switch the default assay back to RNA for further processing
    DefaultAssay(se.meta) = "RNA"
  }
  log_message("DONE: read_in_10x_parallel.")
  return(se.meta)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Normalization
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
normalization_seurat <- function(se.meta, ADT){
  log_message("PROCESS: normalization_seurat")
  # Normalize the RNA assay data using log-normalization (standard in Seurat)
  se.meta = NormalizeData(se.meta, assay = 'RNA', normalization.method = "LogNormalize")
  if (ADT){
    # Normalize the ADT assay data using CLR (Centered Log-Ratio) normalization
    se.meta = NormalizeData(se.meta, assay = "ADT", normalization.method = 'CLR', margin = 2)
    # Ensure the ADT data is stored as a sparse matrix (dgCMatrix), which is memory-efficient for sparse data
    slot(object = se.meta[["ADT"]], name = 'data') = as(se.meta@assays$ADT@data, "dgCMatrix")
  }
  log_message("DONE: normalization_seurat")
  return(se.meta)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Annotation Cell types
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
annotation_exr_and_scgate <- function(se.meta, car_construct){
  log_message("PROCESS: Annotation_exr_and_scgate")
  # Annotate CAR-T cells based on raw counts by calling a custom function
  if (!is.null(car_construct)){
    se.meta <- car_by_expr(obj = se.meta, car_construct)
  }
  # Annotate with scGATE (based on normalized counts)
  scGate_models_DB <- get_scGateDB()
  current_model <-  scGate_models_DB$human$PBMC
  se.meta_annotated <- scGate(se.meta, model = current_model, save.levels = FALSE, multi.asNA = FALSE)
  #annotate T cells with and without
  if (!is.null(car_construct)){
    se.meta_annotated$T_CAR <- ifelse(
      (se.meta_annotated$is.pure_CD8T == "Pure" & se.meta_annotated$CAR_BY_EXPRS == "TRUE") |
      (se.meta_annotated$is.pure_CD4T == "Pure" & se.meta_annotated$CAR_BY_EXPRS == "TRUE") |
      (se.meta_annotated$CAR_BY_EXPRS == "TRUE" & se.meta_annotated$is.pure_gdT == "Pure"),
      TRUE,  # Set T_CAR to TRUE if any of the conditions are met
      FALSE  # Set T_CAR to FALSE otherwise
    )
  }
  se.meta <- se.meta_annotated
  log_message("DONE: annotation_exr_and_scgate")
  return(se.meta)
}
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Add data quality of cells regarding mito and ribo genes and complexity
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
data_quality_mito_ribo_complexity <- function(se.meta){
  log_message("PROCESS: data_quality_mito_ribo_complexity")
  # Calculate the percentage of mitochondrial genes for each cell based on the "^MT-" pattern in the gene names
  se.meta[["Perc_of_mito_genes"]] = Seurat::PercentageFeatureSet(se.meta, pattern = "^MT-")
  # Calculate the percentage of ribosomal genes (RPL and RPS) for each cell based on the pattern "^RPL|^RPS"
  se.meta[["Perc_of_ribosomal_genes"]] = Seurat::PercentageFeatureSet(se.meta, pattern = "^RPL|^RPS")
  # Calculate the log10 ratio of the number of unique genes per UMI for quality control
  se.meta@meta.data$log10GenesPerUMI = log10(se.meta$nFeature_RNA) / log10(se.meta$nCount_RNA)
  log_message("DONE: data_quality_mito_ribo_complexity")
  return(se.meta)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Remove low quality cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
add_low_quality_info <- function(se.meta){
  log_message("PROCESS: add_low_quality_info")
  # Track the number of cells per sample and store the results in the cell.track variable
  cell.track = count_cells_per_sample(c(se.meta))
  # Label cells based on some criteria using the custom function `label_cells_rm`, which possibly removes unwanted cells
  se.meta = label_cells_rm(se.meta)
  se.meta.subset <- tryCatch({
    subset(se.meta, subset = KEEP_CELL == TRUE)  # Attempt subsetting
    cell.track = count_cells_per_sample(c(se.meta.subset), cell.track, "afterFiltering")
  }, error = function(e) {
    log_message("WARNING: No cells passing filtering step. Skipping subsetting.")
    return(NULL)  # Return NULL to handle gracefully in downstream code
  })
  log_message("DONE: add_low_quality_info")
  return(se.meta)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Cell cyle scoring
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cell_cycle_scores <- function(se.meta){
  log_message("PROCESS: cell_cycle_scores")
  ###### Add Cell cylce information seurat scores
  se.meta = CellCycleScoring(
    se.meta, s.features = s.genes, g2m.features = g2m.genes,
    assay = 'RNA', search = TRUE
  )
  ##### Cell cycle based on clustering
  # Load the cell cycle object, typically from the ProjecTILs package, containing relevant cell cycle gene sets
  data(cell.cycle.obj)
  obj.l = Split_Object(se.meta, split.by = "orig.ident", threads = 4)
  # Run the cell cycle estimation in parallel using the custom 'estimate_cc' function on each split object, utilizing 4 cores
  cc.res = parallel::mclapply(obj.l, function(se){
    estimate_cc(se)
  }, mc.cores = 4)
  # Combine the results from each parallel execution into a single data frame using 'rbind'
  cc.res = do.call("rbind", cc.res)
  rownames(cc.res) = cc.res$barcode
  se.meta = AddMetaData(se.meta, cc.res)
  se.meta$barcode = NULL

  se.meta@meta.data = se.meta@meta.data %>%
    dplyr::mutate_if(is.character, as.factor)

  # Annotate cells that were missed by the cluster approach
  cutoff = 0.20
  # Annotate cells with their respective cell cycle phases based on the S and G2M scores
  se.meta$CellCycle_Phase = dplyr::case_when(
    se.meta$S.Score > cutoff & se.meta$G2M.Score <= cutoff ~ "S",
    se.meta$G2M.Score > cutoff & se.meta$S.Score <= cutoff ~ "G2M",
    (se.meta$S.Score > cutoff & se.meta$G2M.Score > cutoff) & (se.meta$S.Score > se.meta$G2M.Score) ~ "S",
    (se.meta$S.Score > cutoff & se.meta$G2M.Score > cutoff) & (se.meta$G2M.Score > se.meta$S.Score) ~ "G2M",
    TRUE ~ se.meta$CellCycle_Phase
  )
  se.meta$CellCycle_Phase = factor(se.meta$CellCycle_Phase, levels = c("G1M", "S", "G2M"))
  se.meta@meta.data$CellCycle = factor(ifelse(se.meta$CellCycle_Phase == "G1M", "FALSE", "TRUE"))

  log_message("PROCESS: cell_cycle_scores")
  return(se.meta)
}


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Estimate cell cylce phase
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#estimate the cell cycle phase of cells in single-cell RNA sequencing (scRNA-seq) data using the Seurat and ProjecTILs R packages,
#alongside the clustifyr package for gene set enrichment analysis (GSEA)
estimate_cc = function(se){

  suppressWarnings({
    suppressMessages({

      data(cell.cycle.obj) # ProjecTILs package
      DefaultAssay(se) = "RNA"

      se = se %>%
        FindVariableFeatures(verbose = F) %>%
        ScaleData(verbose = F) %>%
        RunPCA(npcs = 20, verbose = F) %>%
        FindNeighbors(reduction = "pca", dims = 1:20, verbose = F)

      tmp =  tryCatch(
        FindClusters(se, resolution = 2, verbose = F),
        error=function(e) "error"
      )
      if(class(tmp) != "Seurat") {
        se = FindClusters(se, resolution = 1, verbose = F)
      } else {
        se = tmp
      }

      quiet <- function(x) {
        sink(tempfile())
        on.exit(sink())
        invisible(force(x))
      }

      cc.phase = quiet(
        clustifyr::run_gsea(
          GetAssayData(object = se, assay = "RNA", slot = "data"),
          query_genes = cell.cycle.obj$human$cycling,
          cluster_ids =  se@meta.data[["seurat_clusters"]], n_perm = 1000
        )
      )

      cc.phase$pval_adj = p.adjust(cc.phase$pval, method = "BH")
      cc.cl = rownames(cc.phase[cc.phase$pval < 0.05, ])

      cc.cells = (se@meta.data[["seurat_clusters"]] %in% cc.cl)
      se@meta.data$CellCycle = factor(cc.cells)

      pd.cc = se@meta.data[cc.cells, ]
      pd.cc = pd.cc %>% dplyr::mutate(
        CellCycle_Phase = dplyr::case_when(
          G2M.Score > S.Score ~ "G2M",
          TRUE ~ "S"
        )
      )
      se$CellCycle_Phase = pd.cc$CellCycle_Phase[match(rownames(se@meta.data), rownames(pd.cc))]
      se$CellCycle_Phase[is.na(se$CellCycle_Phase)] = "G1M"
      se$CellCycle_Phase = factor(se$CellCycle_Phase, levels = c("G1M", "S", "G2M"))

      res = se@meta.data %>% dplyr::select(CellCycle, CellCycle_Phase)
      res$barcode = rownames(res)
      res
    })
  })
}

get_soup_groups <- function(sobj){
  sobj <- NormalizeData(sobj, verbose = FALSE)
  sobj <- FindVariableFeatures(
    object = sobj, nfeatures = 2000, verbose = FALSE, selection.method = 'vst'
  )
  sobj <- ScaleData(sobj, verbose = FALSE)
  sobj <- RunPCA(sobj, npcs = 20, verbose = FALSE)
  sobj <- FindNeighbors(sobj, dims = 1:20, verbose = FALSE)
  sobj <- FindClusters(sobj, resolution = 0.5, verbose = FALSE)
  return(sobj@meta.data[['seurat_clusters']])
}

make_soup <- function(sobj,raw_dir, ADT){
  stopifnot(dir.exists(raw_dir))
  sample_id = as.character(sobj$orig.ident[1])
  print(sample_id)
  print(raw_dir)
  data = Seurat::Read10X(data.dir = raw_dir, gene.column = 2)
  # Create a Seurat object using the gene expression (RNA) counts, setting the project name as the dataset ID
  if (!ADT) {
    raw = Seurat::CreateSeuratObject(counts = data, project = sample_id) # if no adt
  } else {
    raw = Seurat::CreateSeuratObject(counts = data[[1]], project = sample_id) # if adt
    # Add ADT (Antibody-Derived Tag) counts as a new assay to the Seurat object
    raw[['ADT']] = Seurat::CreateAssayObject(counts = data[[2]]) #get adt counts
  }
  raw_counts <- GetAssayData(raw, layer = "counts")
  fltrd_counts <- GetAssayData(sobj, layer = "counts")
  sc <- SoupChannel(raw_counts, fltrd_counts)
  sc <- setClusters(sc, sobj$soup_group)
  sc  <- autoEstCont(sc)
  out = adjustCounts(sc, roundToInt = TRUE)

  #optional keep original
  sobj[["original.counts"]] <- CreateAssayObject(
    counts = GetAssayData(sobj, layer = "counts")
  )
  sobj@assays$RNA$counts <- out

  return(sobj)
}
