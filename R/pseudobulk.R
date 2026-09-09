## max. code width ============================================================

## Title: Arc Institute in vivo pertubation screens
## Author: DKatzlberger
## Date: 21.01.2026

## Reproducibility ============================================================
seed <- 42
set.seed(seed)

## Settings ===================================================================
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  project_root <- args[1]
  input_path   <- args[2]
  output_path  <- "."

  # Meta
  facLoc <- args[3]
  nOn    <- args[4]
  koEff  <- args[5]
  nGuide <- args[6]

  file_name <- paste0("pseudobulk_facLoc_", facLoc, "_nOn_", nOn, "_koEff_", koEff, "_nGuide_", nGuide)

} else {
  project_root <- "/home/scc_e_928652/ArcSimulationSynthetic"
  input_path   <- "results/manual/offtarget_mixed/subset/subset_facLoc_2_nOn_10_koEff_25_nGuide_50.h5ad"
  output_path  <- "results/manual/offtarget_mixed/pseudobulk"

  stopifnot(dirname(dirname(input_path)) == dirname(output_path))

  # Meta
  facLoc <- 2
  nOn    <- 10
  koEff  <- 25
  nGuide <- 50

  file_name <- paste0("pseudobulk_facLoc_", facLoc, "_nOn_", nOn, "_koEff_", koEff, "_nGuide_", nGuide)
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

## Dataset ====================================================================
.read_h5ad <- fun_with_ram(anndataR::read_h5ad)

txt     <- sprintf("%s as 'SingleCellExperiment'", input_path)
sce_obj <- .read_h5ad(input_path, as = "SingleCellExperiment", info = txt)

# Gene truths (simulations)
sce_rd <- SingleCellExperiment::rowData(sce_obj)
sce_rd <- as.data.frame(sce_rd)

## Pseudobulk =================================================================
.aggregate_expression <- fun_with_ram(scuttle::aggregateAcrossCells)

groups  <- c("sample", "guide")
txt     <- sprintf("pseudobulk by %s", paste(groups, collapse = ", "))
sce_obj <- .aggregate_expression(
  sce_obj, 
  ids = colData(sce_obj)[, groups], 
  statistics = "sum",
  use.assay.type = "counts",
  info = txt
)

## Observed annotations =======================================================
sce_rd_new  <- sce_rd[rownames(sce_obj), , drop = FALSE]

# Log counts
sce_obj <- scater::logNormCounts(sce_obj)

cells_ntc <- which(sce_obj$group == "non-targeting")
cells_tar <- which(sce_obj$group == "targeting")

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
txt         <- sprintf("saved as %s", paste0(file_name, ".h5ad"))

# Integer check
stopifnot(all(SingleCellExperiment::counts(sce_obj) %% 1 == 0))

.write_h5ad(
  object = sce_obj, 
  path = file.path(
    output_path,
    paste0(file_name, ".h5ad")
  ),
  mode = "w",
  info = txt
)
