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

  formula <- args[7]
  ref     <- args[8]

  # Meta
  facLoc <- args[3]
  nOn    <- args[4]
  koEff  <- args[5]
  nGuide <- args[6]

  file_name <- paste0("sceptre_facLoc_", facLoc, "_nOn_", nOn, "_koEff_", koEff, "_nGuide_", nGuide, "_", formula)

} else {
  project_root <- "/home/scc_e_928652/ArcSimulationSynthetic"
  input_path   <- "results/manual/offtarget_mixed/subset/subset_facLoc_2_nOn_10_koEff_25_nGuide_50.h5ad"
  output_path  <- "results/manual/offtarget_mixed/sceptre"

  stopifnot(dirname(dirname(input_path)) == dirname(output_path))

  formula <- "~genotype"
  ref     <- "NTC"

  # Meta
  facLoc <- 2
  nOn    <- 10
  koEff  <- 25
  nGuide <- 50

  file_name <- paste0("sceptre_facLoc_", facLoc, "_nOn_", nOn, "_koEff_", koEff, "_nGuide_", nGuide, "_", formula)
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
    "sceptre"
  )
)

## Dataset ====================================================================
.read_h5ad <- fun_with_ram(anndataR::read_h5ad)

txt     <- sprintf("%s as 'SingleCellExperiment'", input_path)
sce_obj <- .read_h5ad(input_path, as = "SingleCellExperiment", info = txt)

# Gene truths (simulations)
sce_rd <- SingleCellExperiment::rowData(sce_obj)
sce_rd <- as.data.frame(sce_rd)

## Formula ====================================================================
.parse_model_formula <- fun_with_ram(parse_model_formula)

covs <- .parse_model_formula(formula)

target_covariate <- covs$target_covariate
extra_covariates <- covs$extra_covariates

## DGE ========================================================================
matr <- SingleCellExperiment::counts(sce_obj)
meta <- SingleCellExperiment::colData(sce_obj)

# Target
guide_ids <- unique(meta$guide)

guide_tars <- sapply(guide_ids, function(g) {
  unique(meta[[target_covariate]][meta$guide == g])
})

guide_tars <- ifelse(
  guide_tars == ref, 
  "non-targeting", 
  sprintf("%s%s", target_covariate, guide_tars)
)

guide_tars <- data.frame(
  grna_id     = guide_ids,
  grna_target = guide_tars,
  stringsAsFactors = FALSE
)

guide_matr <- sapply(guide_ids, function(g) {
  as.numeric(meta$guide == g)
})
colnames(guide_matr) <- guide_ids
rownames(guide_matr) <- rownames(meta)

guide_matr <- t(guide_matr)
guide_matr <- guide_matr * 100

# Extra covariates
if (length(extra_covariates) > 0) {
  extra_covariates <- as.data.frame(
    meta[, extra_covariates, drop = FALSE]
  ) 
} else {
  extra_covariates <- data.frame()
}

# Fit
sceptre_obj <- sceptre::import_data(
  response_matrix = matr,
  grna_matrix = guide_matr,
  grna_target_data_frame = guide_tars,
  moi = "low",
  extra_covariates = data.frame(),
  response_names = NA_character_,
  use_ondisc = FALSE,
  directory_to_write = NULL
)

pos_ctrl_pair  <- sceptre::construct_positive_control_pairs(
  sceptre_object = sceptre_obj
)

discovery_pair <- sceptre::construct_trans_pairs(
  sceptre_object         = sceptre_obj,
  positive_control_pairs = pos_ctrl_pair
)

sceptre_obj <- sceptre::set_analysis_parameters(
  sceptre_object         = sceptre_obj,
  discovery_pairs        = discovery_pair,
  positive_control_pairs = pos_ctrl_pair,
)

# Fit
fit <- sceptre::assign_grnas(sceptre_object = sceptre_obj)
fit <- sceptre::run_qc(sceptre_object = fit)
fit <- sceptre::run_discovery_analysis(sceptre_object = fit)

# Save 
.saveRDS  <- fun_with_ram(saveRDS)
txt       <- "saving 'sceptre' fit obj."

.saveRDS(
  fit, 
  file = file.path(
    output_path,
    paste0(file_name, "_fit.rds")
  ),
  info = txt
)

## Statistics =================================================================
.extract_sceptre <- fun_with_ram(extract_sceptre)

# Keep only target covariate
txt <- "getting 'sceptre' stats"
dge <- .extract_sceptre(fit, info = txt)
dge <- dge[grepl(paste0("^", target_covariate), dge$coef), ]

# Add genotype
dge$genotype <- sub(paste0("^", target_covariate), "", dge$coef)

# Adding gene truth
dge <- merge(
  x = dge, 
  y = sce_rd, 
  by.x = "gene", 
  by.y = "gene", 
  all.x = TRUE
)

# Add information
dge$formula     <- formula
dge$method      <- "sceptre"
dge$workflow    <- "singlecell"
dge$method_type <- "semi-parametric"

# Save
txt <- "saving 'sceptre' dge obj. (incl. gene truth)"
.saveRDS(
  dge, 
  file = file.path(
    output_path,
    paste0(file_name, "_dge.rds")
  )
)
