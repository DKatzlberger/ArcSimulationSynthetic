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
  nGuide <- as.numeric(args[6])

  file_name <- paste0("subset_facLoc_", facLoc, "_nOn_", nOn, "_koEff_", koEff, "_nGuide_", nGuide)

} else {
  project_root <- "/home/scc_e_928652/ArcSimulationSynthetic"
  input_path   <- "results/manual/offtarget_mixed/scdesign3/scdesign3_facLoc_2_nOn_10_koEff_25.h5ad"
  output_path  <- "results/manual/offtarget_mixed/subset"

  stopifnot(dirname(dirname(input_path)) == dirname(output_path))

  # Meta
  facLoc <- 2
  nOn    <- 10
  koEff  <- 25
  nGuide <- 50

  file_name <- paste0("subset_facLoc_", facLoc, "_nOn_", nOn, "_koEff_", koEff, "_nGuide_", nGuide)
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