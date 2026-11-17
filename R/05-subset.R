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
  nGuide <- as.numeric(trail_args[3])

  # Output
  input_name  <- tools::file_path_sans_ext(basename(input_path))
  output_name <- paste0(script_id, "-", input_name)
  output_name <- paste0(output_name, "_nGuide_", nGuide)

} else {
  # Script
  script_path <- rstudioapi::getActiveDocumentContext()$path
  script_name <- tools::file_path_sans_ext(basename(script_path))
  script_id   <- sub(pattern = "-.*", replacement = "", script_name)

  # Paths
  project_root <- "/home/scc_e_928652/ArcSimulationSynthetic"
  input_path   <- "results/manual/ArcSimSyn_offtarget_mixed/scdesign3_offtarget_mixed/celltype_Gran. P_knockout_NULL_facLoc_2_nOn_10_koEff_100.h5ad"
  output_path  <- file.path("results/manual/ArcSimSyn_offtarget_mixed", script_name)

  stopifnot(dirname(dirname(input_path)) == dirname(output_path))
  if (!dir.exists(output_path)) dir.create(path = output_path, recursive = TRUE)

  # Params
  nGuide <- 50

  # Output
  input_name  <- tools::file_path_sans_ext(basename(input_path))
  output_name <- paste0(script_id, "-", input_name)
  output_name <- paste0(output_name, "_nGuide_", nGuide)
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


## Subset =====================================================================
### Sampling ------------------------------------------------------------------
cells_ntc <- which(sce_obj$group == "NTC")
cells_tar <- which(sce_obj$group == "KO")

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
stopifnot(length(which(sce_obj$group == "KO")) / length(names(guide_counts)) == nGuide)

### Ground truth --------------------------------------------------------------
sce_rd_new <- sce_rd[rownames(sce_obj), , drop = FALSE]

sce_rd_new$nGuide <- nGuide

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
    paste0(
      output_name, 
      ".h5ad"
    )
  ),
  mode = "w",
  info = txt
)


## Script info ================================================================
showf("\n◷ %s ====================================================\n", date(), type = "cat")
