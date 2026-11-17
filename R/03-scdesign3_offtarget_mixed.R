## max. code width ============================================================

## Title: Arc Institute in vivo pertubation screens
## Author: DKatzlberger
## Date: 21.01.2026

## This script learns marginal distributions and gene correlations from 
## single-cell RNA-seq data via 'scDesign3' (ZINB linear regression models).

## Simulates a Perturb-seq dataset with targeting and non-targeting (NTC)
## cells. Knockouts have 3 unique guide RNAs (2 on-target, 1 off-target).
## The off-target guide shares on-target genes but impacts extra genes.

## Genes are split into: downstream, off-target, and background sets.
## Mean parameters for all sets are conditioned on sample info 'I'.

## Downstream genes depend on knockout expression 'C'. Off-target genes
## depend on 'C' only for off-target guides. For on-target guides, they
## act as background. Background genes depend only on sample info 'I'.
## The marginal distributions of genes are grouped for on-target guides
## and off target guides. 

## Next the gene-gene correlations are learned using gaussian copulas. This 
## learns network-level dependencies for within samples 'I' and knockout 
## indication 'K' across genes.

## The mean of targeting cells are modified by modifying 'C' by sampling form
## expression values from NTC cells (currently outdated).

## Then expression counts for each are sampled using the fitted paramters and 
## the gene-gene correlation matrix to generate a Perturb-Seq dataset.

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
  fracOff  <- as.numeric(trail_args[3])
  koEff    <- as.numeric(trail_args[4])
  nDown    <- Inf
  nBack    <- Inf

  # nDown <- 50
  # nBack <- 50

  # Output
  input_name  <- tools::file_path_sans_ext(basename(input_path))
  output_name <- paste0(script_id, "-", input_name)
  output_name <- paste0(output_name, "_fracOff_", fracOff, "_koEff_", koEff)

} else {
  # Script
  script_path <- rstudioapi::getActiveDocumentContext()$path
  script_name <- tools::file_path_sans_ext(basename(script_path))
  script_id   <- sub(pattern = "-.*", replacement = "", script_name)

  # Paths
  project_root <- "/home/scc_e_928652/ArcSimulationSynthetic"
  input_path   <- "results/manual/ArcSimSyn_offtarget_mixed/02-splatter/02-01-celltype_Gran. P_knockout_NULL_facLoc_2.h5ad"
  output_path  <- file.path("results/manual/ArcSimSyn_offtarget_mixed", script_name)

  stopifnot(dirname(dirname(input_path)) == dirname(output_path))
  if (!dir.exists(output_path)) dir.create(path = output_path, recursive = TRUE)

  # Params
  fracOff <- 50
  koEff   <- 100
  nDown   <- 20
  nBack   <- 20

  # Output
  input_name  <- tools::file_path_sans_ext(basename(input_path))
  output_name <- paste0(script_id, "-", input_name)
  output_name <- paste0(output_name, "_fracOff_", fracOff, "_koEff_", koEff)
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
    "scDesign3"
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


## scDesign3 ==================================================================
### Prepare -------------------------------------------------------------------
knockout <- sce_rd$gene[sce_rd$status == "Knockout"]
stopifnot(length(knockout) == 1L, !is.na(knockout))

counts   <- SingleCellExperiment::counts(sce_obj)
if (anyNA(counts)) showf("Counts contain NA values.", type = "stop")
if (any(!is.finite(counts))) showf("Counts contain non-finite values.", type = "stop")

# Library size (pipeline)
lib_size <- sce_obj$lib_size_orig
max_lib  <- quantile(lib_size, 0.999)
stopifnot(all(is.finite(lib_size)), all(lib_size > 0))

Y <- counts[rownames(counts) != knockout, , drop = FALSE]
C <- as.numeric(counts[knockout, ])
stopifnot(all(is.finite(C)), all(C >= 0))

K <- factor(ifelse(sce_obj$genotype == knockout, "1", "0"), levels = c("0", "1"))
I <- factor(sce_obj$sample)

# Interaction between knockout and sample
K_I <- factor(paste0("K", K, "_I", I))

sce_obj <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = Y),
  colData = data.frame(
    C = C,
    K = K,
    I = I,
    K_I = K_I,
    library = lib_size,
    celltype = "macrophages",
    row.names = colnames(Y)
  )
)

# Downstream genes
genes_downstream  <- sce_rd$gene[sce_rd$status == "Downstream"]
downstream_reads  <- sce_rd$total_reads[sce_rd$status == "Downstream"]
genes_downstream  <- genes_downstream[order(downstream_reads, decreasing = TRUE)]
genes_downstream  <- genes_downstream[seq_len(min(nDown, length(genes_downstream)))]

# On- & off-target
nOff <- round(length(genes_downstream) * (fracOff / 100))
genes_off_target <- if (nOff > 0) genes_downstream[seq_len(nOff)] else character(0)
genes_downstream <- if (nOff < length(genes_downstream)) genes_downstream[seq.int(nOff + 1, length(genes_downstream))] else character(0)
stopifnot(!any(genes_downstream %in% genes_off_target))

# Background
genes_background <- sce_rd$gene[sce_rd$status == "Background"]
background_reads <- sce_rd$total_reads[sce_rd$status == "Background"]
genes_background <- genes_background[order(background_reads, decreasing = TRUE)]
genes_background <- genes_background[seq_len(min(nBack, length(genes_background)))]

### Marginals -----------------------------------------------------------------
.estimate_marginals <- fun_with_ram(estimate_marginals)

txt <- sprintf("on-target downstream: conditioned by %s expression 'C' and sample info 'I'", knockout)
marg_on_downstream <- .estimate_marginals(
  sce_obj = sce_obj[genes_downstream, ],
  assay_use = "counts",
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I", "library"),
  mu_formula = "K + I + offset(log(library))",
  sigma_formula = "1",
  family_use = "zinb",
  corr_by = "K_I",
  n_cores = 1,
  verbose = FALSE,
  info = txt
)

txt <- sprintf("off-target downstream: conditioned by %s expression 'C' and sample info 'I'", knockout)
marg_off_downstream <- .estimate_marginals(
  sce_obj = sce_obj[genes_off_target, ],
  assay_use = "counts",
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I", "library"),
  mu_formula = "K + I + offset(log(library))",
  sigma_formula = "1",
  family_use = "zinb",
  corr_by = "K_I",
  n_cores = 1,
  verbose = FALSE,
  info = txt
)

txt <- sprintf("off-target background: conditioned by sample info 'I'")
marg_off_background <- .estimate_marginals(
  sce_obj = sce_obj[genes_off_target, ],
  assay_use = "counts",
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I", "library"),
  mu_formula = "I + offset(log(library))",
  sigma_formula = "1",
  family_use = "zinb",
  corr_by = "K_I",
  n_cores = 1,
  verbose = FALSE,
  info = txt
)

txt <- sprintf("actual background: conditioned by sample info 'I'")
marg_background <- .estimate_marginals(
  sce_obj = sce_obj[genes_background, ],
  assay_use = "counts",
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I", "library"),
  mu_formula = "I + offset(log(library))",
  sigma_formula = "1",
  family_use = "zinb",
  corr_by = "K_I",
  n_cores = 1,
  verbose = FALSE,
  info = txt
)

# Validate marginals
marg_on_downstream  <- marg_on_downstream[sapply(marg_on_downstream, function(x) is_valid_fit(x$fit))]
marg_off_downstream <- marg_off_downstream[sapply(marg_off_downstream, function(x) is_valid_fit(x$fit))]
marg_off_background <- marg_off_background[sapply(marg_off_background, function(x) is_valid_fit(x$fit))]
marg_background     <- marg_background[sapply(marg_background, function(x) is_valid_fit(x$fit))]

# Build complete marginal sets 
marg_on  <- c(marg_on_downstream, marg_off_background, marg_background)
marg_off <- c(marg_on_downstream, marg_off_downstream, marg_background)

### Correlation ---------------------------------------------------------------
.estimate_correlation <- fun_with_ram(estimate_correlation)

# On-target copula
copula_on_downstream <- .estimate_correlation(
  sce_obj = sce_obj[names(marg_on), ],
  marginals = marg_on,
  covariates = c("C", "K", "I", "K_I", "library"),
  corr_by = "K_I",
  family_use = "zinb",
  copula = "gaussian",
  corr_features = rownames(sce_obj[names(marg_on), ]) %in% genes_downstream,
  n_cores = 1,
  verbose = FALSE
)

copula_on_background <- .estimate_correlation(
  sce_obj = sce_obj[names(marg_on), ],
  marginals = marg_on,
  covariates = c("C", "K", "I", "K_I", "library"),
  corr_by = "K_I",
  family_use = "zinb",
  copula = "gaussian",
  corr_features = rownames(sce_obj[names(marg_on), ]) %in% c(genes_off_target, genes_background),
  n_cores = 1,
  verbose = FALSE
)

stopifnot(identical(names(copula_on_downstream$copula_list), names(copula_on_background$copula_list)))

# Combine copula
copula_on <- copula_on_downstream
copula_on$copula_list <- Map(
  function(R_down, R_bg) {

    stopifnot(identical(dimnames(R_down), dimnames(R_bg)))
    # Combine the two correlation blocks
    R <- R_down + R_bg - diag(nrow(R_down))
    dimnames(R) <- dimnames(R_down)
    R
  },
  copula_on_downstream$copula_list,
  copula_on_background$copula_list
)

# Update important features: downstream OR background
copula_on$important_feature <- copula_on_downstream$important_feature | copula_on_background$important_feature

# Off-target copula
copula_off_downstream <- .estimate_correlation(
  sce_obj = sce_obj[names(marg_off), ],
  marginals = marg_off,
  covariates = c("C", "K", "I", "K_I", "library"),
  corr_by = "K_I",
  family_use = "zinb",
  copula = "gaussian",
  corr_features = rownames(sce_obj[names(marg_on), ]) %in% c(genes_downstream, genes_off_target),
  n_cores = 1,
  verbose = FALSE
)

copula_off_background <- .estimate_correlation(
  sce_obj = sce_obj[names(marg_off), ],
  marginals = marg_off,
  covariates = c("C", "K", "I", "K_I", "library"),
  corr_by = "K_I",
  family_use = "zinb",
  copula = "gaussian",
  corr_features = rownames(sce_obj[names(marg_on), ]) %in% genes_background,
  n_cores = 1,
  verbose = FALSE
)

stopifnot(identical(names(copula_off_downstream$copula_list), names(copula_off_background$copula_list)))

# Combine copula
copula_off <- copula_off_downstream
copula_off$copula_list <- Map(
  function(R_down, R_bg) {

    stopifnot(identical(dimnames(R_down), dimnames(R_bg)))
    # Combine the two correlation blocks
    R <- R_down + R_bg - diag(nrow(R_down))
    dimnames(R) <- dimnames(R_down)
    R
  },
  copula_off_downstream$copula_list,
  copula_off_background$copula_list
)

# Update important features: downstream OR background
copula_off$important_feature <- copula_off_downstream$important_feature | copula_off_background$important_feature

### Modify C (intervention) ---------------------------------------------------
.modify_covariate <- fun_with_ram(modify_covariate)

txt <- "modifying cells for on-target effect"
new_cov_on <- .modify_covariate(
  sce_obj = sce_obj[names(marg_on), ],
  covariates = c("C", "K", "I", "K_I", "library"),
  corr_by = "K_I",       
  ko_eff = koEff,
  mode = "partial",
  info = txt
)

txt <- "modifying cells for off-target effect"
new_cov_off <- .modify_covariate(
  sce_obj = sce_obj[names(marg_off), ],
  covariates = c("C", "K", "I", "K_I", "library"),
  corr_by = "K_I",      
  ko_eff = koEff,
  mode = "partial",
  info = txt
)

txt <- "modifying cells for non-targeting control"
new_cov_ntc <- .modify_covariate(
  sce_obj = sce_obj[names(marg_on), ],
  covariates = c("C", "K", "I", "K_I", "library"),
  corr_by = "K_I",      
  ko_eff = NULL,
  mode = "partial",
  info = txt
)

### Params --------------------------------------------------------------------
.extract_params <- fun_with_ram(extract_params)

params_on <- .extract_params(
  sce_obj = sce_obj[names(marg_on), ],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I", "library"),
  new_covariate = new_cov_on[new_cov_on$K == "1", , drop = FALSE],
  marginals = marg_on,
  family_use = "zinb",
  corr_by = "K_I"
)

params_off <- .extract_params(
  sce_obj = sce_obj[names(marg_off), ],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I", "library"),
  new_covariate = new_cov_off[new_cov_off$K == "1", , drop = FALSE],
  marginals = marg_off,
  family_use = "zinb",
  corr_by = "K_I"
)

params_ntc <- .extract_params(
  sce_obj = sce_obj[names(marg_on), ],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I", "library"),
  new_covariate = new_cov_ntc[new_cov_ntc$K == "0", , drop = FALSE],
  marginals = marg_on,
  family_use = "zinb",
  corr_by = "K_I"
)

### Counts --------------------------------------------------------------------
#### Targeting guides ---------------------------------------------------------
.sim_pipe <- fun_with_ram(sim_pipe)

new_cov_on   <- new_cov_on[new_cov_on$K == "1", , drop = FALSE]
new_cov_on$G <- NA_integer_

for (s in unique(new_cov_on$I)) {
  idx <- which(new_cov_on$I == s)
  new_cov_on$G[idx] <- sample(rep(1:3, length.out = length(idx)))
}

stopifnot(!anyNA(new_cov_on$G))

# Allocate cells per guide
cells_tar_1 <- rownames(new_cov_on[new_cov_on$G == 1, , drop = FALSE])
cells_tar_2 <- rownames(new_cov_on[new_cov_on$G == 2, , drop = FALSE])
cells_tar_3 <- rownames(new_cov_on[new_cov_on$G == 3, , drop = FALSE])

new_cov_on$G <- NULL

# On-target guides
txt <- "simulating cells for on-target guide 1"
sim_tar_1 <- .sim_pipe(
  sce_obj = sce_obj[names(marg_on), cells_tar_1],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I", "library"),
  params = lapply(params_on, function(x) x[cells_tar_1, , drop = FALSE]),
  copulas = copula_on,
  corr_by = "K_I",
  family_use = "zinb",
  new_covariate = new_cov_on[cells_tar_1, , drop = FALSE],
  knockout = knockout,
  seed = 1000,
  info = txt
)

sim_tar_1$group     <- "KO"
sim_tar_1$genotype  <- knockout
sim_tar_1$guide     <- paste0(knockout, "-1")
sim_tar_1$on_guide  <- paste0(knockout, "-1")
sim_tar_1$on_target <- TRUE

txt <- "simulating cells for on-target guide 2"
sim_tar_2 <- .sim_pipe(
  sce_obj = sce_obj[names(marg_on), cells_tar_2],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I", "library"),
  params = lapply(params_on, function(x) x[cells_tar_2, , drop = FALSE]),
  copulas = copula_on,
  corr_by = "K_I",
  family_use = "zinb",
  new_covariate = new_cov_on[cells_tar_2, , drop = FALSE],
  knockout = knockout,
  seed = 2000,
  info = txt
)

sim_tar_2$group     <- "KO"
sim_tar_2$genotype  <- knockout
sim_tar_2$guide     <- paste0(knockout, "-2")
sim_tar_2$on_guide  <- paste0(knockout, "-2")
sim_tar_2$on_target <- TRUE

# Off-target guides
txt <- "simulating cells for off-target guide 1"
sim_tar_3 <- .sim_pipe(
  sce_obj = sce_obj[names(marg_off), cells_tar_3],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I", "library"),
  params = lapply(params_off, function(x) x[cells_tar_3, , drop = FALSE]),
  copulas = copula_off,
  corr_by = "K_I",
  family_use = "zinb",
  new_covariate = new_cov_off[cells_tar_3, , drop = FALSE],
  knockout = knockout,
  seed = 3000,
  info = txt
)

sim_tar_3$group     <- "KO"
sim_tar_3$genotype  <- knockout
sim_tar_3$guide     <- paste0(knockout, "-3")
sim_tar_3$on_guide  <- paste0(knockout, "-3")
sim_tar_3$on_target <- FALSE

#### Non-targeting guides -----------------------------------------------------
new_cov_ntc   <- new_cov_ntc[new_cov_ntc$K == "0", , drop = FALSE]
new_cov_ntc$G <- NA_integer_

for (s in unique(new_cov_ntc$I)) {
    idx <- which(new_cov_ntc$I == s)
    new_cov_ntc$G[idx] <- sample(rep(1:3, length.out = length(idx)))
}

stopifnot(!anyNA(new_cov_ntc$G))

# Allocate cells per guide
cells_ntc_1 <- rownames(new_cov_ntc[new_cov_ntc$G == 1, , drop = FALSE])
cells_ntc_2 <- rownames(new_cov_ntc[new_cov_ntc$G == 2, , drop = FALSE])
cells_ntc_3 <- rownames(new_cov_ntc[new_cov_ntc$G == 3, , drop = FALSE])

new_cov_ntc$G <- NULL

txt <- "simulating cells for NTC guide 1"
sim_ntc_1 <- .sim_pipe(
  sce_obj = sce_obj[names(marg_on), cells_ntc_1],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I", "library"),
  params = lapply(params_ntc, function(x) x[cells_ntc_1, , drop = FALSE]),
  copulas = copula_on,
  corr_by = "K_I",
  family_use = "zinb",
  new_covariate = new_cov_ntc[cells_ntc_1, , drop = FALSE],
  knockout = knockout,
  seed = 4000,
  info = txt
)

sim_ntc_1$group     <- "NTC"
sim_ntc_1$genotype  <- "NTC"
sim_ntc_1$guide     <- "mouseNTC-1"
sim_ntc_1$on_guide  <- "NTC"
sim_ntc_1$on_target <- TRUE

txt <- "simulating cells for NTC guide 2"
sim_ntc_2 <- .sim_pipe(
  sce_obj = sce_obj[names(marg_on), cells_ntc_2],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I", "library"),
  params = lapply(params_ntc, function(x) x[cells_ntc_2, , drop = FALSE]),
  copulas = copula_on,
  corr_by = "K_I",
  family_use = "zinb",
  new_covariate = new_cov_ntc[cells_ntc_2, , drop = FALSE],
  knockout = knockout,
  seed = 5000,
  info = txt
)

sim_ntc_2$group     <- "NTC"
sim_ntc_2$genotype  <- "NTC"
sim_ntc_2$guide     <- "mouseNTC-2"
sim_ntc_2$on_guide  <- "NTC"
sim_ntc_2$on_target <- TRUE

txt <- "simulating cells for NTC guide 3"
sim_ntc_3 <- .sim_pipe(
  sce_obj = sce_obj[names(marg_on), cells_ntc_3],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I", "library"),
  params = lapply(params_ntc, function(x) x[cells_ntc_3, , drop = FALSE]),
  copulas = copula_on,
  corr_by = "K_I",
  family_use = "zinb",
  new_covariate = new_cov_ntc[cells_ntc_3, , drop = FALSE],
  knockout = knockout,
  seed = 6000,
  info = txt
)

sim_ntc_3$group     <- "NTC"
sim_ntc_3$genotype  <- "NTC"
sim_ntc_3$guide     <- "mouseNTC-3"
sim_ntc_3$on_guide  <- "NTC"
sim_ntc_3$on_target <- TRUE

### Ground truth --------------------------------------------------------------
sce_obj <- BiocGenerics::cbind(sim_tar_1, sim_tar_2, sim_tar_3, sim_ntc_1, sim_ntc_2, sim_ntc_3)

sce_rd_new            <- sce_rd[rownames(sce_obj), , drop = FALSE]
sce_rd_new$guide_mode <- "mixed off-target effect"
sce_rd_new$fracOff    <- fracOff
sce_rd_new$koEff      <- koEff

# Log counts 
sce_obj   <- scater::logNormCounts(sce_obj)
cells_ntc <- which(sce_obj$group == "NTC")
cells_tar <- which(sce_obj$group == "KO")

log_counts <- SummarizedExperiment::assay(sce_obj, "logcounts")
obs_logFC  <- Matrix::rowMeans(log_counts[, cells_tar]) - Matrix::rowMeans(log_counts[, cells_ntc])

sce_rd_new$obs_logFC   <- obs_logFC
sce_rd_new$total_reads <- Matrix::rowSums(SingleCellExperiment::counts(sce_obj))

# Gene response
sce_rd_new$status                                              <- "Background"
sce_rd_new$status[rownames(sce_rd_new) == knockout]            <- "Knockout"
sce_rd_new$status[rownames(sce_rd_new) %in% genes_downstream]  <- "On-target"
sce_rd_new$status[rownames(sce_rd_new) %in% genes_off_target]  <- "Off-target"

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


## Visualize ==================================================================
### Panel A -------------------------------------------------------------------
#### Plot data ----------------------------------------------------------------
rd <- SummarizedExperiment::rowData(sce_obj)

# Col order
guide_levels <- c(
  sort(grep(knockout, unique(sce_obj$guide), value = TRUE)),
  sort(grep("NTC",  unique(sce_obj$guide), value = TRUE))
)
sce_obj$guide <- factor(sce_obj$guide, levels = guide_levels)

col_ord <- order(sce_obj$genotype, sce_obj$guide)

# Row order
get_sorted_indices <- function(status_name, row_data) {
  idx <- which(row_data$status == status_name)
  if (length(idx) == 0) return(integer(0))
  idx[order(row_data$sim_logFC[idx], decreasing = TRUE, na.last = TRUE)]
}

ko_idx   <- get_sorted_indices("Knockout", rd)
on_idx   <- get_sorted_indices("On-target", rd)
off_idx  <- get_sorted_indices("Off-target", rd)
back_idx <- get_sorted_indices("Background", rd)

top_ontarget   <- on_idx[order(rd$total_reads[on_idx], decreasing = TRUE)][1:pmin(length(on_idx), 24)]
top_offtarget  <- off_idx[order(rd$total_reads[off_idx], decreasing = TRUE)][1:pmin(length(off_idx), 15)]
top_background <- back_idx[order(rd$total_reads[back_idx], decreasing = TRUE)][1:pmin(length(back_idx), 15)]

top_ontarget   <- top_ontarget[order(match(top_ontarget, on_idx))]
top_offtarget  <- top_offtarget[order(match(top_offtarget, off_idx))]
top_background <- top_background[order(match(top_background, back_idx))]

row_ord <- c(ko_idx, top_ontarget, top_offtarget, top_background)

# Matrix
mat <- as.matrix(SummarizedExperiment::assay(sce_obj, "logcounts"))
mat <- mat[row_ord, col_ord, drop = FALSE]

#### Panel A01: Heatmap -------------------------------------------------------
col_mat <- circlize::colorRamp2(
  c(0, max(mat, na.rm = TRUE) * 0.4, max(mat, na.rm = TRUE)), 
  c("#000000", "#A52A2A", "#FFD700") 
)

# Top column annotations
col_guides <- setNames(
  wesanderson::wes_palette(
    "Darjeeling2", 
    n = length(levels(sce_obj$guide)), 
    type = "continuous"
  ),
  levels(sce_obj$guide)
)

top_col_anno <- cmplx_anno(
  `On-target` = sce_obj$on_target[col_ord],
  `Guide` = sce_obj$guide[col_ord],
  col = list(
    `On-target` = c(
      "TRUE"  = wesanderson::wes_palette("Moonrise3")[1],
      "FALSE" = wesanderson::wes_palette("Moonrise3")[2]
    ),
    `Guide` = col_guides
  ),
  which = "column",
  annotation_name_side = "left",
  show_annotation_name = TRUE
)

# Bottom column annotations
lib_sizes       <- Matrix::colSums(SummarizedExperiment::assay(sce_obj, "counts"))
bottom_col_anno <- cmplx_anno(
  `Sample` = sce_obj$sample[col_ord],
  `Group` = sce_obj$group[col_ord],
  `Lib. size` = cmplx_bar(
    lib_sizes[col_ord], 
    which = "column", 
    axis_side = "right",
    flip = TRUE
  ),
  col = list(
    `Sample` = c(
      "mouse1" = wesanderson::wes_palette("Rushmore1")[3], 
      "mouse2" = wesanderson::wes_palette("Rushmore1")[4], 
      "mouse3" = wesanderson::wes_palette("Rushmore1")[5]
    ),
    `Group` = c(
      "KO"  = wesanderson::wes_palette("Moonrise2")[3],
      "NTC" = wesanderson::wes_palette("Moonrise2")[1]
    )
  ),
  which = "column",
  annotation_name_side = "right",
  show_annotation_name = TRUE,
  annotation_name_rot = c(
    `Sample` = 0,     
    `Group`  = 0,    
    `Lib. size` = 270   
  )
)

# Left Row Annotations
status     <- rd$status[row_ord]
sim_logfc  <- rd$sim_logFC[row_ord]
obs_logfc  <- rd$obs_logFC[row_ord]

finite_logfc <- c(sim_logfc, obs_logfc)
finite_logfc <- finite_logfc[is.finite(finite_logfc)]
max_logfc    <- max(abs(finite_logfc), na.rm = TRUE)

col_logfc <- circlize::colorRamp2(c(-max_logfc, 0, max_logfc), c("#313695", "#F7F7F7", "#B2182B"))

left_row_anno <- cmplx_anno(
  `Gene status`  = status,
  `Sim. log2FC`  = sim_logfc,
  `Obs. log2FC`  = obs_logfc,
  col = list(
    `Gene status` = colors_gene_status,
    `Sim. log2FC` = col_logfc,
    `Obs. log2FC` = col_logfc
  ),
  which = "row",
  show_annotation_name = TRUE,
  annotation_name_side = "bottom"
)

# Right row annotations
total_reads    <- rd$total_reads[row_ord]
right_row_anno <- cmplx_anno(
  `Total reads` = cmplx_bar(
    total_reads, 
    which = "row", 
    axis_side = "top"
  ), 
  which = "row",
  annotation_name_side = "top", 
  annotation_name_rot = 0,
  show_annotation_name = TRUE
)

# Titles
col_title <- unique(sapply(sce_obj$genotype[col_ord], function(g) {
  sprintf("%s\n(n = %s)", g, table(sce_obj$genotype)[g])
}))

row_title <- sprintf("%s%% off-target genes with a knockout efficiencey of %s%%", fracOff, koEff)

# Heatmap 
p_panel_A01 <- cmplx_ht(
  mat,
  name = "Log-norm.\ncounts",
  col = col_mat,
  
  # Row settings
  left_annotation = left_row_anno,
  right_annotation = right_row_anno,
  cluster_rows = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_title = row_title,
  
  # Col settings
  top_annotation = top_col_anno,
  bottom_annotation = bottom_col_anno,
  column_split = sce_obj$genotype[col_ord],
  column_title = col_title, 
  cluster_columns = FALSE,          
  show_column_names = FALSE,
  use_raster = TRUE
)

#### Patchwork & save ---------------------------------------------------------
p_panel_A <- 
  p_panel_A01

# Save
.cmplx_save <- fun_with_ram(cmplx_save)
txt         <- sprintf("saved as %s", paste0(output_name, "_A.svg"))
.cmplx_save(
  p_panel_A,
  file = file.path(
    output_path,
    paste0(
      output_name, 
      "_A.svg"
    )
  ),
  width = 15,
  height = 10,
  units = "cm",
  info = txt
)

### Panel B -------------------------------------------------------------------
#### Plot data ----------------------------------------------------------------
rd <- SingleCellExperiment::rowData(sce_obj)

get_sorted_indices <- function(status_name, row_data) {
  idx <- which(row_data$status == status_name)
  if (length(idx) == 0) return(integer(0))
  idx[order(row_data$sim_logFC[idx], decreasing = TRUE, na.last = TRUE)]
}

ko_idx   <- get_sorted_indices("Knockout", rd)
on_idx   <- get_sorted_indices("On-target", rd)
off_idx  <- get_sorted_indices("Off-target", rd)
back_idx <- get_sorted_indices("Background", rd)

top_ontarget   <- on_idx[order(rd$total_reads[on_idx], decreasing = TRUE)][1:pmin(length(on_idx), 24)]
top_offtarget  <- off_idx[order(rd$total_reads[off_idx], decreasing = TRUE)][1:pmin(length(off_idx), 15)]
top_background <- back_idx[order(rd$total_reads[back_idx], decreasing = TRUE)][1:pmin(length(back_idx), 15)]

top_ontarget   <- top_ontarget[order(match(top_ontarget, on_idx))]
top_offtarget  <- top_offtarget[order(match(top_offtarget, off_idx))]
top_background <- top_background[order(match(top_background, back_idx))]

selected_genes <- c(ko_idx, top_ontarget, top_offtarget, top_background)
selected_genes <- rownames(sce_obj)[selected_genes]

# Norm. counts
mat <- as.matrix(SingleCellExperiment::logcounts(sce_obj))
mat <- mat[selected_genes, , drop = FALSE]

# Plot data
p_data <- as.data.frame(t(mat))
p_data$guide <- sce_obj$guide

p_data <- reshape(
  p_data, 
  varying = selected_genes, 
  v.names = "expression", 
  timevar = "gene", 
  times = selected_genes, 
  direction = "long",
  new.row.names = 1:(nrow(p_data) * length(selected_genes))
)

ko_names   <- rownames(sce_obj)[ko_idx]
on_names   <- rownames(sce_obj)[top_ontarget]
off_names  <- rownames(sce_obj)[top_offtarget]

p_data$status <- "Background"
p_data$status[p_data$gene %in% ko_names]  <- "Knockout"
p_data$status[p_data$gene %in% on_names]  <- "On-target"
p_data$status[p_data$gene %in% off_names] <- "Off-target"

# Standardize levels to perfectly map palette keys across backgrounds & geoms
p_data$status <- factor(p_data$status, levels = c("Knockout", "On-target", "Off-target", "Background"))
p_data$gene   <- factor(p_data$gene, levels = selected_genes)

guide_levels <- c(
  sort(grep(knockout,  unique(p_data$guide), value = TRUE)),
  sort(grep("NTC",  unique(p_data$guide), value = TRUE))
)
p_data$guide <- factor(p_data$guide, levels = guide_levels)

# Background cols
legend_cols    <- colors_gene_status
present_status <- unique(as.character(p_data$status))
legend_cols    <- legend_cols[names(legend_cols) %in% present_status]

bg_cols_with_alpha <- paste0(legend_cols, "40")
names(bg_cols_with_alpha) <- names(legend_cols)

#### Panel B01: Boxplot -------------------------------------------------------
p_panel_B01 <- ggplot2::ggplot(
) +
ggplot2::geom_rect(
  data = unique(p_data[, c("gene", "status")]),
  mapping = ggplot2::aes(
    fill = status
  ), 
  xmin = -Inf, 
  xmax = Inf, 
  ymin = -Inf,
  ymax = Inf, 
  show.legend = TRUE 
) +
ggplot2::scale_fill_manual(
  name = "Gene response",
  values = bg_cols_with_alpha
) +
ggplot2::guides(
  fill = ggplot2::guide_legend(
    override.aes = list(
      fill = legend_cols
    ) 
  )
) +
ggnewscale::new_scale_fill() + 
ggplot2::geom_violin(
  data = p_data,
  mapping = ggplot2::aes(
    x = guide,
    y = expression,
    fill = guide
  ),
  alpha = 0.7,
  scale = "width", 
  color = "transparent",
  show.legend = TRUE
) +
ggplot2::geom_boxplot(
  data = p_data,
  mapping = ggplot2::aes(
    x = guide,
    y = expression
  ),
  width = 0.15, 
  linewidth = 0.15,
  fill = "white", 
  color = "grey30", 
  outlier.shape = NA
) +
ggplot2::facet_wrap(
  ~ gene, 
  scales = "free_y", 
  ncol = 5
) +
ggplot2::scale_fill_manual(
  name = "Guide",
  values = col_guides
) +
ggplot2::scale_y_continuous(
  n.breaks = 4,
  labels = scales::label_number(accuracy = 1)
) +
ggplot2::labs(
  x = "Guide",
  y = "Log counts"
) +
theme_DK(
  show_grid = FALSE,
  rotate = 45
)

#### Patchwork & save ---------------------------------------------------------
.ggsave <- fun_with_ram(ggplot2::ggsave)

# Plot
p_panel_B <-
  p_panel_B01 +
  patchwork::plot_layout(
    guides = "collect"
  ) 

# Save
txt <- sprintf("saved as %s", paste0(output_name, "_B.svg"))
.ggsave(
  p_panel_B,
  file = file.path(
    output_path,
    paste0(
      output_name,
      "_B.svg"
    )
  ),
  width = 15,
  height = 10,
  units = "cm",
  info = txt
)

### Panel C -------------------------------------------------------------------
#### Plot data ----------------------------------------------------------------
SummarizedExperiment::assay(sce_obj, "logcpm") <- log2(scuttle::calculateCPM(sce_obj) + 1)

counts <- SummarizedExperiment::assay(sce_obj, "counts")
logcpm <- SummarizedExperiment::assay(sce_obj, "logcpm")

# Gene response
gene_status <- setNames(rep("On-target", nrow(counts)), rownames(counts))

gene_status[intersect(genes_downstream, rownames(counts))] <- "On-target"
gene_status[intersect(genes_off_target, rownames(counts))] <- "Off-target"
gene_status[intersect(genes_background, rownames(counts))] <- "Background"

gene_status_levels <- c("On-target", "Off-target", "Background")
gene_status <- factor(gene_status, levels = gene_status_levels)

# Groups
g <- split(seq_along(sce_obj$group), interaction(sce_obj$group, sce_obj$guide, drop = TRUE, sep = "_"))

# Gene summary
gene_summary <- do.call(rbind, lapply(names(g), function(grp) {
  cells <- g[[grp]]

  data.frame(
    gene = rownames(logcpm),
    gene_status = gene_status[rownames(logcpm)],
    mean_counts = rowMeans(counts[, cells, drop = FALSE]),
    mean_logCPM = rowMeans(logcpm[, cells, drop = FALSE]),
    var_logCPM = apply(logcpm[, cells, drop = FALSE], 1, var),
    group = sce_obj$group[cells[1]],
    guide = sce_obj$guide[cells[1]]
  )
}))

gene_zero <- do.call(rbind, lapply(names(g), function(grp) {
  cells <- g[[grp]]

  data.frame(
    gene = rownames(counts),
    gene_status = gene_status[rownames(counts)],
    pct_zero = Matrix::rowMeans(counts[, cells, drop = FALSE] == 0) * 100,
    group = sce_obj$group[cells[1]],
    guide = sce_obj$guide[cells[1]]
  )
}))

# Cell summary
lib_sizes <- do.call(rbind, lapply(
  gene_status_levels,
  function(status) {

    genes <- names(gene_status)[gene_status == status]
    if (length(genes) == 0) return(NULL)

    lib <- Matrix::colSums(counts[genes, , drop = FALSE])

    data.frame(
      cell = colnames(counts),
      gene_status = status,
      total_counts = as.numeric(lib),
      group = sce_obj$group,
      guide = sce_obj$guide
    )
  }
))
lib_sizes$gene_status <- factor(lib_sizes$gene_status, levels = gene_status_levels)

cell_zero <- do.call(rbind, lapply(
  gene_status_levels,
  function(status) {

    genes <- names(gene_status)[gene_status == status]
    if (length(genes) == 0) return(NULL)

    data.frame(
      cell = colnames(counts),
      gene_status = status,
      pct_zero = Matrix::colMeans(counts[genes, , drop = FALSE] == 0) * 100,
      group = sce_obj$group,
      guide = sce_obj$guide
    )
  }
))
cell_zero$gene_status <- factor(cell_zero$gene_status, levels = gene_status_levels)

#### Panel C01: Distribution of mean ------------------------------------------
p_panel_C01 <- ggplot2::ggplot(
  data = gene_summary,
  mapping = ggplot2::aes(
    x = group,
    y = mean_logCPM,
    color = guide
  )
) +
ggplot2::geom_boxplot(
  linewidth = 0.3,
  coef = Inf
) +
ggplot2::scale_color_manual(
  name = "Guide",
  values = col_guides
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(gene_status)
) +
ggplot2::labs(
  x = "Group",
  y = expression("Mean log"[2] * "(CPM + 1)"),
  color = "Guide"
) +
theme_DK(
  rotate = 45,
  show_grid = FALSE
)

#### Panel C02: Distribution of variance --------------------------------------
p_panel_C02 <- ggplot2::ggplot(
  data = gene_summary,
  mapping = ggplot2::aes(
    x = group,
    y = var_logCPM,
    color = guide
  )
) +
ggplot2::geom_boxplot(
  linewidth = 0.3,
  coef = Inf
) +
ggplot2::scale_color_manual(
  name = "Guide",
  values = col_guides
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(gene_status)
) +
ggplot2::labs(
  x = "Group",
  y = expression("Variance log"[2] * "(CPM + 1)"),
  color = "Guide"
) +
theme_DK(
  rotate = 45,
  show_grid = FALSE
)

#### Panel C03: Mean variance relationship ------------------------------------
p_panel_C03 <- ggplot2::ggplot(
  data = gene_summary,
  mapping = ggplot2::aes(
    x = mean_logCPM,
    y = var_logCPM,
    color = guide
  )
) +
ggplot2::scale_color_manual(
  name = "Guide",
  values = col_guides
) +
ggplot2::geom_smooth(
  formula = "y ~ x",
  method = "loess",
  se = FALSE,
  linewidth = 0.2,
  show.legend = FALSE
) +
ggplot2::labs(
  x = expression("Mean log"[2] * "(CPM + 1)"),
  y = expression("Variance log"[2] * "(CPM + 1)"),
  color = "Guide"
) +
theme_DK(
  rotate = 45,
  show_grid = FALSE
)

#### Panel C04: Library size --------------------------------------------------
p_panel_C04 <- ggplot2::ggplot(
  data = lib_sizes,
  mapping = ggplot2::aes(
    x = group,
    y = total_counts,
    color = guide
  )
) +
ggplot2::geom_boxplot(
  linewidth = 0.3,
  coef = Inf
) +
ggplot2::scale_color_manual(
  name = "Guide",
  values = col_guides
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(gene_status)
) +
ggplot2::labs(
  x = "Group",
  y = "Total counts per cell",
  color = "Guide"
) +
theme_DK(
  rotate = 45,
  show_grid = FALSE
)

#### Panel C05: Percentage zero per cell --------------------------------------
p_panel_C05 <- ggplot2::ggplot(
  data = cell_zero,
  mapping = ggplot2::aes(
    x = group,
    y = pct_zero,
    color = guide
  )
) +
ggplot2::geom_boxplot(
  linewidth = 0.3,
  coef = Inf
) +
ggplot2::scale_color_manual(
  name = "Guide",
  values = col_guides
) +
ggplot2::scale_y_continuous(
  limits = c(0, 100),
  breaks = seq(0, 100, 25)
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(gene_status)
) +
ggplot2::labs(
  x = "Group",
  y = "Percentage zeros per cell",
  color = "Guide"
) +
theme_DK(
  rotate = 45,
  show_grid = FALSE
)

#### Panel C06: Percentage zero per gene --------------------------------------
p_panel_C06 <- ggplot2::ggplot(
  data = gene_zero,
  mapping = ggplot2::aes(
    x = group,
    y = pct_zero,
    color = guide
  )
) +
ggplot2::geom_boxplot(
  linewidth = 0.3,
  coef = Inf
) +
ggplot2::scale_color_manual(
  name = "Guide",
  values = col_guides
) +
ggplot2::scale_y_continuous(
  limits = c(0, 100),
  breaks = seq(0, 100, 25)
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(gene_status)
) +
ggplot2::labs(
  x = "Group",
  y = "Percentage zero per gene",
  color = "Guide"
) +
theme_DK(
  rotate = 45,
  show_grid = FALSE
)

#### Panel C07: Zero mean relationship ----------------------------------------
p_panel_C07 <- ggplot2::ggplot(
  data = cbind(
    gene_zero,
    mean_logCPM = gene_summary$mean_counts
  ),
  mapping = ggplot2::aes(
    x = mean_logCPM,
    y = pct_zero,
    color = guide
  )
) +
ggplot2::geom_smooth(
  formula = "y ~ x",
  method = "loess",
  se = FALSE,
  linewidth = 0.2,
  show.legend = FALSE
) +
ggplot2::scale_color_manual(
  name = "Guide",
  values = col_guides
) +
ggplot2::scale_y_continuous(
  limits = c(0, 100),
  breaks = seq(0, 100, 25)
) +
ggplot2::labs(
  x = "Mean count",
  y = "Percentage zero",
  color = "Guide"
) +
theme_DK(
  rotate = 45,
  show_grid = FALSE
)

#### Patchwork & save ---------------------------------------------------------
.ggsave <- fun_with_ram(ggplot2::ggsave)

# Plot
p_panel_C <-
  p_panel_C01 +
  p_panel_C02 +
  p_panel_C03 +
  p_panel_C04 +
  p_panel_C05 +
  p_panel_C06 +
  p_panel_C07 +
  patchwork::guide_area() +
  patchwork::plot_layout(
    ncol = 4,
    nrow = 2,
    guides = "collect"
  ) 

# Save
txt <- sprintf("saved as %s", paste0(output_name, "_C.svg"))
.ggsave(
  p_panel_C,
  file = file.path(
    output_path,
    paste0(
      output_name,
      "_C.svg"
    )
  ),
  width = 18,
  height = 8,
  units = "cm",
  info = txt
)


## Script info ================================================================
showf("\n◷ %s ====================================================\n", date(), type = "cat")
