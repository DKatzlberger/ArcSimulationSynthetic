## max. code width ============================================================

## Title: Arc Institute in vivo pertubation screens
## Author: DKatzlberger
## Date: 21.01.2026

# "Preparing a Seurat object for synthetic benchmarking of",
# "differential gene expression methods in Perturb-seq data.",
# "",
# "First, the Seurat object is subset to the specified celltype.", 
# "Second guides are assigned to cells based on maximum CRISPR counts.",
# "Third, the Seurat object is subset to the analysed knockout.",
# "Finally, the Seurat object is converted to a SingleCellExperiment object.\n",

## Reproducibility ============================================================
seed <- 42
set.seed(seed)


## Settings ===================================================================
if (!interactive()) {
  # Commands
  full_args  <- commandArgs(trailingOnly = FALSE)
  trail_args <- commandArgs(trailingOnly = TRUE)

  # Script
  script_path <- sub("^--file=", "", grep("^--file=", full_args, value = TRUE)[1])
  script_name <- tools::file_path_sans_ext(basename(script_path))
  script_id   <- sub(pattern = "-.*", replacement = "", script_name)

  # Paths
  project_root <- trail_args[1]
  input_path   <- gsub(pattern = "\\\\ ", replacement = " ", trail_args[2])
  output_path  <- "."

  # Params
  celltype <- as.character(trail_args[3])
  knockout <- as.character(trail_args[4])

  knockout      <- if (identical(knockout, "NULL")) NULL else knockout
  knockout_name <- if (is.null(knockout)) "NULL" else knockout

  # Output
  output_name <- paste0(script_id, "-")
  output_name <- paste0(output_name, "celltype_", celltype, "_knockout_", knockout_name)

} else {
  # Script
  script_path <- rstudioapi::getActiveDocumentContext()$path
  script_name <- tools::file_path_sans_ext(basename(script_path))
  script_id   <- sub(pattern = "-.*", replacement = "", script_name)

  # Paths
  project_root <- "/home/scc_e_928652/ArcSimulationSynthetic"
  input_path   <- "data_generated/GSE213511/rna_integrated.rds"
  output_path  <- file.path("results/manual/ArcSimSyn_offtarget_none", script_name)

  if (!dir.exists(output_path)) dir.create(path = output_path, recursive = TRUE)

  # Params
  celltype <- "Gran. P"
  knockout <- "NULL"

  knockout      <- if (identical(knockout, "NULL")) NULL else knockout
  knockout_name <- if (is.null(knockout)) "NULL" else knockout

  # Output
  output_name <- paste0(script_id, "-")
  output_name <- paste0(output_name, "celltype_", celltype, "_knockout_", knockout_name)
}

## r-utils ====================================================================
source(file.path(project_root, "R", "utils.R"))


## Environment ================================================================
scc_info <- scc_status()

script_startup(
  root = project_root,
  renv = TRUE,
  bioc = TRUE,
  packages = c(
    "SummarizedExperiment",
    "SingleCellExperiment",
    "Seurat"
  )
)

## Script info ================================================================
showf(fmt = "\n%s\n", line80(paste0(script_name, ".R")), type = "cat")


## Dataset ====================================================================
.read_rds  <- fun_with_ram(base::readRDS)
.read_h5ad <- fun_with_ram(anndataR::read_h5ad)
.check_obj <- fun_with_ram(check_obj)

ext <- tools::file_ext(input_path)
txt <- sprintf("'%s' as 'Seurat'", input_path)

seu_obj <- switch(EXPR = ext,
  rds  = .read_rds(input_path, info = txt),
  h5ad = .read_h5ad(input_path, as = "Seurat", info = txt),
  showf("Unsupported file extension: %s", ext, type = "stop")
)

# Check input
txt <- "checking input obj."
.check_obj(seu_obj, info = txt)


## Prepare ====================================================================
### Celltype ------------------------------------------------------------------
.subset <- fun_with_ram(base::subset)

txt     <- sprintf("subsetting to celltype '%s'", celltype)
keep    <- colnames(seu_obj)[seu_obj$celltype %in% celltype]
seu_obj <- .subset(seu_obj, cells = keep, info = txt)

### Genotype, guide & on_guide ------------------------------------------------
.assign_guide_by_max_count <- fun_with_ram(assign_guide_by_max_count)

seu_obj <- .assign_guide_by_max_count(seu_obj)

# Check
has_KO  <- !is.null(knockout) && any(knockout %in% seu_obj$genotype)
has_NTC <- "NTC" %in% seu_obj$genotype

if (!has_KO && !has_NTC) showf("Neither '%s' nor 'NTC' assigned to any cells.", knockout, type = "stop")
if (!is.null(knockout) && !has_KO) showf("'%s' not assigned to any cells.", knockout, type = "stop")
if (!has_NTC) showf("'%s' not assigned to any cells.", "NTC", type = "stop")

### Cells: subset by knockout, NTC & feature counts ---------------------------
if (is.null(knockout)) {
  txt  <- "subsetting to 'NTC'"
  keep <- colnames(seu_obj)[seu_obj$genotype %in% "NTC"]
} else {
  txt  <- sprintf("subsetting to genotype '%s' and 'NTC'", paste(knockout, collapse = ", "))
  keep <- colnames(seu_obj)[seu_obj$genotype %in% c(knockout, "NTC")]
}

seu_obj <- .subset(seu_obj, cells = keep, info = txt)
seu_obj <- subset(seu_obj, subset = nFeature_RNA >= 1000 & nFeature_RNA <= 5000)

### Genes: subset by detection ------------------------------------------------
matr <- SeuratObject::LayerData(seu_obj, assay = "RNA", layer = "counts")
keep <- Matrix::rowSums(matr > 0) >= 10

seu_obj <- subset(seu_obj, features = rownames(matr)[keep])

### Convert to SingleCellExperiment -------------------------------------------
matr <- SeuratObject::LayerData(seu_obj, assay = "RNA", layer = "counts")
meta <- seu_obj@meta.data

sce_obj <- SingleCellExperiment::SingleCellExperiment(
  assays  = list(counts = matr),
  colData = meta
)

# Create rowData
sce_rd <- do.call(
  what = S4Vectors::DataFrame,
  args = c(
    list(
      sim_mode = "synthetic",
      celltype = as.character(celltype),
      knockout = if (is.null(knockout)) rep("NULL", nrow(sce_obj)) else knockout,
      n_genes = nrow(sce_obj),
      n_cells = ncol(sce_obj),
      gene = rownames(sce_obj)
    )
  )
)

# Add the rowData
SummarizedExperiment::rowData(sce_obj) <- NULL
SummarizedExperiment::rowData(sce_obj) <- sce_rd


## Save =======================================================================
.write_h5ad <- fun_with_ram(anndataR::write_h5ad)

# Integer check
stopifnot(all(SingleCellExperiment::counts(sce_obj) %% 1 == 0))

# Check output
txt <- "checking output obj."
.check_obj(sce_obj, info = txt)

txt <- sprintf("saved as %s", paste0(output_name, ".h5ad"))
.write_h5ad(
  object = sce_obj, 
  path = file.path(
    output_path,
    paste0(output_name, ".h5ad")
  ),
  mode = "w",
  info = txt
)


## Script info ================================================================
showf("\n◷ %s ====================================================\n", date(), type = "cat")

