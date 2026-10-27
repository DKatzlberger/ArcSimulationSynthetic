## max. code width ============================================================

## Title: Arc Institute in vivo pertubation screens
## Author: DKatzlberger
## Date: 21.01.2026

## Reproducibility ============================================================
seed <- 42
set.seed(seed)

## Settings ===================================================================
if (!interactive()) {
  full_args  <- commandArgs(trailingOnly = FALSE)
  trail_args <- commandArgs(trailingOnly = TRUE)

  # Commands
  script_path <- sub("^--file=", "", grep("^--file=", full_args, value = TRUE)[1])
  script_name <- tools::file_path_sans_ext(basename(script_path))

  project_root <- trail_args[1]
  input_path   <- trail_args[2]
  output_path  <- "."

  # Meta
  # celltype <- as.character(trail_args[3])
  # facLoc   <- as.numeric(trail_args[4])
  # nOn      <- as.numeric(trail_args[5])
  # koEff    <- as.numeric(trail_args[6])

  # Params
  nGuide <- as.numeric(args[7])

  # Output file
  input_name <- tools::file_path_sans_ext(basename(input_path))
  file_name  <- paste0(input_name, "_nGuide_", nGuide)

} else {
  # Commands
  script_path  <- rstudioapi::getActiveDocumentContext()$path
  script_name  <- tools::file_path_sans_ext(basename(script_path))

  project_root <- "/home/scc_e_928652/ArcSimulationSynthetic"
  input_path   <- "results/manual/ArcSimSyn_offtarget_mixed/scdesign3/celltype_Gran. P_facLoc_2_nOn_10_koEff_100.h5ad"
  output_path  <- file.path("results/manual/ArcSimSyn_offtarget_mixed", script_name)
  stopifnot(dirname(dirname(input_path)) == dirname(output_path))

  # Meta
  # celltype <- "Gran. P"
  # facLoc   <- 2
  # nOn      <- 10
  # koEff    <- 25

  # Params
  nGuide <- 300

  # Output file
  input_name <- tools::file_path_sans_ext(basename(input_path))
  file_name  <- paste0(input_name, "_nGuide_", nGuide)
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
.read_rds  <- fun_with_ram(base::readRDS)
.read_h5ad <- fun_with_ram(anndataR::read_h5ad)

ext <- tools::file_ext(input_path)
txt <- sprintf("%s as 'SingleCellExperiment'", input_path)

sce_obj <- switch(EXPR = ext,
  rds  = .read_rds(input_path, info = txt),
  h5ad = .read_h5ad(input_path, as = "SingleCellExperiment", info = txt),
  showf("Unsupported file extension: %s", ext, type = "stop")
)

# Gene truth (pipeline)
sce_rd <- SingleCellExperiment::rowData(sce_obj)
sce_rd <- as.data.frame(sce_rd)

# Meta
celltype <- unique(sce_rd$celltype)
facLoc   <- unique(sce_rd$facLoc)
nOn      <- unique(sce_rd$nOn)
koEff    <- unique(sce_rd$koEff)

## Subset =====================================================================
cells_ntc <- which(sce_obj$group == "non-targeting")
cells_tar <- which(sce_obj$group == "targeting")

# Check n guides
guide_idx    <- split(cells_tar, sce_obj$guide[cells_tar])
guide_counts <- lengths(guide_idx)

stopifnot(!any(guide_counts < nGuide))

# Stratified sampling
keep_tar <- unlist(
  lapply(guide_idx, function(idx) {
    
    # Cells for this guide, split by sample
    sample_idx <- split(idx, sce_obj$sample[idx])
    
    # Allocate nGuide cells across samples
    sample_counts <- lengths(sample_idx)
    
    # Proportional allocation
    n_sample <- floor(
      sample_counts / sum(sample_counts) * nGuide
    )
    
    # Make sure every represented sample gets at least 1 cell
    n_sample[n_sample == 0] <- 1
    
    # Adjust rounding so total = exactly nGuide
    while (sum(n_sample) > nGuide) {
      eligible <- which(n_sample > 1)
      j <- eligible[which.max(n_sample[eligible])]
      n_sample[j] <- n_sample[j] - 1
    }
    
    while (sum(n_sample) < nGuide) {
      eligible <- which(n_sample < sample_counts)
      j <- eligible[which.max(sample_counts[eligible] - n_sample[eligible])]
      n_sample[j] <- n_sample[j] + 1
    }
    
    # Sample from each sample
    unlist(
      Map(
        function(idx, n) sample(idx, n),
        sample_idx,
        n_sample
      )
    )
  })
)

# Subset
sce_obj <- sce_obj[, c(cells_ntc, keep_tar)]

## Observed annotations =======================================================
sce_rd_new  <- sce_rd[rownames(sce_obj), , drop = FALSE]

sce_rd_new$nGuide <- nGuide

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

## Save sim. ==================================================================
.write_h5ad <- fun_with_ram(anndataR::write_h5ad)
txt         <- sprintf("saved as %s", paste0(file_name, ".h5ad"))

# Integer check
counts  <- SingleCellExperiment::counts(sce_obj)
non_int <- (counts %% 1 != 0)
non_int_genes <- rownames(counts)[rowSums(non_int) > 0]

stopifnot(length(non_int_genes) == 0)

.write_h5ad(
  object = sce_obj, 
  path = file.path(
    output_path,
    paste0(file_name, ".h5ad")
  ),
  mode = "w",
  info = txt
)