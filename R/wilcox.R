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

  file_name <- paste0("wilcox_facLoc_", facLoc, "_nOn_", nOn, "_koEff_", koEff, "_nGuide_", nGuide, "_", formula)

} else {
  project_root <- "/home/scc_e_928652/ArcSimulationSynthetic"
  input_path   <- "results/manual/offtarget_mixed/subset/subset_facLoc_2_nOn_10_koEff_25_nGuide_50.h5ad"
  output_path  <- "results/manual/offtarget_mixed/wilcox"

  stopifnot(dirname(dirname(input_path)) == dirname(output_path))

  formula <- "~genotype"
  ref     <- "NTC"

  # Meta
  facLoc <- 2
  nOn    <- 10
  koEff  <- 25
  nGuide <- 50

  file_name <- paste0("wilcox_facLoc_", facLoc, "_nOn_", nOn, "_koEff_", koEff, "_nGuide_", nGuide, "_", formula)
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
random_covariate <- covs$random_covariates

## DGE ========================================================================
seu_obj <- Seurat::as.Seurat(sce_obj, counts = "counts")
seu_obj <- Seurat::NormalizeData(seu_obj, verbose = FALSE)

Seurat::Idents(seu_obj) <- target_covariate
all_groups <- unique(Seurat::Idents(seu_obj))
test_group <- setdiff(all_groups, ref)

fit <- Seurat::FindMarkers(
  object = seu_obj,
  ident.1 = test_group,  
  ident.2 = ref,
  test.use = "wilcox",
  verbose = FALSE
)

fit$coef <- paste0("genotype", test_group)

# Save
.saveRDS  <- fun_with_ram(saveRDS)
txt       <- "saving 'wilcox' fit obj."

.saveRDS(
  fit,
  file = file.path(
    output_path,
    paste0(file_name, "_fit.rds")
  ),
  info = txt
)

## Statistics =================================================================
.extract_seurat <- fun_with_ram(extract_seurat)

# Keep only target covariate
txt <- "getting 'seurat wilcox' stats"
dge <- .extract_seurat(fit, info = txt)
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
dge$method      <- "wilcox"
dge$workflow    <- "singlecell"
dge$method_type <- "non-parametric"

# Save
txt <- "saving 'wilcox' dge obj. (incl. gene truth)"
.saveRDS(
  dge, 
  file = file.path(
    output_path,
    paste0(file_name, "_dge.rds")
  ),
  info = txt
)
