## max. code width ============================================================

## Title: Arc Institute in vivo pertubation screens
## Author: DKatzlberger
## Date: 21.01.2026

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
  groups <- strsplit(trail_args[3], split = ",", fixed = TRUE)[[1]]

  # Output
  input_name  <- tools::file_path_sans_ext(basename(input_path))
  output_name <- paste0(script_id, "-", input_name)
  output_name <- paste0(output_name)

} else {
  # Commands
  script_path <- rstudioapi::getActiveDocumentContext()$path
  script_name <- tools::file_path_sans_ext(basename(script_path))
  script_id   <- sub(pattern = "-.*", replacement = "", script_name)
  
  # Paths
  project_root <- "/home/scc_e_928652/ArcSimulationSynthetic"
  input_path   <- "results/manual/ArcSimSyn_singlecell_baseline/02-splatter/02-01-celltype_Gran. P_knockout_NULL_facLoc_2.h5ad"
  output_path  <- file.path("results/manual/ArcSimSyn_singlecell_baseline", script_name)

  stopifnot(dirname(dirname(input_path)) == dirname(output_path))
  if (!dir.exists(output_path)) dir.create(path = output_path, recursive = TRUE)

  # Params
  groups <- c("sample", "guide")

  # Output
  input_name  <- tools::file_path_sans_ext(basename(input_path))
  output_name <- paste0(script_id, "-", input_name)
  output_name <- paste0(output_name)
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
    "SingleCellExperiment"
  )
)


## Script info ================================================================
showf(fmt = "\n%s\n", line80(paste0(script_name, ".R")), type = "cat")


## Dataset ====================================================================
.read_rds  <- fun_with_ram(base::readRDS)
.read_h5ad <- fun_with_ram(anndataR::read_h5ad)
.check_obj <- fun_with_ram(check_obj)

ext <- tools::file_ext(input_path)
txt <- sprintf(fmt = "'%s' as 'SingleCellExperiment'", input_path)

sce_obj <- switch(EXPR = ext,
  rds  = .read_rds(input_path, info = txt),
  h5ad = .read_h5ad(input_path, as = "SingleCellExperiment", info = txt),
  showf("Unsupported file extension: %s", ext, type = "stop")
)

# Check input
txt <- "checking input obj."
.check_obj(sce_obj, info = txt)

# Gene truth (pipeline)
sce_rd <- SingleCellExperiment::rowData(sce_obj)
sce_rd <- as.data.frame(sce_rd)


## Pseudobulk =================================================================
### Aggregate -----------------------------------------------------------------
.aggregate_expression <- fun_with_ram(scuttle::aggregateAcrossCells)

txt     <- sprintf("pseudobulk by %s", paste(groups, collapse = ", "))
sce_obj <- .aggregate_expression(
  sce_obj, 
  ids = colData(sce_obj)[, groups], 
  statistics = "sum",
  use.assay.type = "counts",
  info = txt
)

### Ground truth --------------------------------------------------------------
sce_rd_new  <- sce_rd[rownames(sce_obj), , drop = FALSE]

# Log counts
sce_obj <- scater::logNormCounts(sce_obj)

cells_ntc <- which(sce_obj$group == "NTC")
cells_tar <- which(sce_obj$group == "KO")

log_counts <- SummarizedExperiment::assay(sce_obj, "logcounts")
obs_logFC  <- Matrix::rowMeans(log_counts[, cells_tar]) - Matrix::rowMeans(log_counts[, cells_ntc])

sce_rd_new$obs_logFC   <- obs_logFC
sce_rd_new$total_reads <- Matrix::rowSums(SingleCellExperiment::counts(sce_obj))

# Assign back rowData
SummarizedExperiment::rowData(sce_obj) <- NULL
SummarizedExperiment::rowData(sce_obj) <- sce_rd_new
sce_obj$sizeFactor <- NULL
sce_obj$ExpLibSize <- NULL


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
