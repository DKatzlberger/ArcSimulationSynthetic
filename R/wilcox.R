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
  # facLoc <- args[3]
  # nOn    <- args[4]
  # koEff  <- args[5]
  # nGuide <- args[6]

  # Params
  formula <- trail_args[7]
  ref     <- trail_args[8]
  gene_qc <- TRUE

  # Output file
  input_name <- tools::file_path_sans_ext(basename(input_path))
  file_name  <- paste0(input_name, "_", formula)

} else {
  # Commands
  script_path  <- rstudioapi::getActiveDocumentContext()$path
  script_name  <- tools::file_path_sans_ext(basename(script_path))

  project_root <- "/home/scc_e_928652/ArcSimulationSynthetic"
  input_path   <- "results/manual/ArcSimSyn_offtarget_mixed/subset/celltype_Gran. P_facLoc_2_nOn_10_koEff_100_nGuide_300.h5ad"
  output_path  <- file.path("results/manual/ArcSimSyn_offtarget_mixed", script_name)
  stopifnot(dirname(dirname(input_path)) == dirname(output_path))

  # Meta
  # facLoc <- 2
  # nOn    <- 10
  # koEff  <- 25
  # nGuide <- 50

  # Params
  formula <- "~genotype + sample"
  ref     <- "NTC"
  gene_qc <- TRUE

  # Output file
  input_name <- tools::file_path_sans_ext(basename(input_path))
  file_name  <- paste0(input_name, "_", formula)
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
.read_rds  <- fun_with_ram(base::readRDS)
.read_h5ad <- fun_with_ram(anndataR::read_h5ad)

ext <- tools::file_ext(input_path)
txt <- sprintf("%s as 'SingleCellExperiment'", input_path)

sce_obj <- switch(ext,
  rds  = .read_rds(input_path, info = txt),
  h5ad = .read_h5ad(input_path, as = "SingleCellExperiment", info = txt),
  showf("Unsupported file extension: %s", ext, type = "stop")
)

# Gene truth (pipeline)
sce_rd <- SingleCellExperiment::rowData(sce_obj)
sce_rd <- as.data.frame(sce_rd)

## Formula ====================================================================
.parse_model_formula <- fun_with_ram(parse_model_formula)

covs <- .parse_model_formula(formula)

target_covariate <- covs$target_covariate
extra_covariates <- covs$extra_covariates
random_covariate <- covs$random_covariates

# Fixed formula
fixed_formula <- paste("~", paste(c(target_covariate, extra_covariates), collapse = " + "))
fixed_formula <- as.formula(fixed_formula)

## Model matrix ===============================================================
matr <- SingleCellExperiment::counts(sce_obj)
meta <- SingleCellExperiment::colData(sce_obj)

meta[[target_covariate]] <- factor(meta[[target_covariate]])
meta[[target_covariate]] <- relevel(meta[[target_covariate]], ref = ref)

mm <- model.matrix(fixed_formula, data = meta)

# N cells per group
groups <- meta[[target_covariate]]
target <- setdiff(levels(groups), ref)
stopifnot(length(target) == 1)

sce_rd$n_cells_target <- sum(groups == target, na.rm = TRUE)
sce_rd$n_cells_ref    <- sum(groups == ref, na.rm = TRUE)

## Wilcox prep. ===============================================================
seu_obj <- Seurat::as.Seurat(sce_obj, counts = "counts", data = NULL)
seu_obj <- Seurat::NormalizeData(seu_obj, verbose = FALSE)

Seurat::Idents(seu_obj) <- target_covariate
all_groups <- unique(Seurat::Idents(seu_obj))
test_group <- setdiff(all_groups, ref)

## Fit ========================================================================
t0 <- Sys.time()
fit <- Seurat::FindMarkers(
  object = seu_obj,
  ident.1 = test_group,  
  ident.2 = ref,
  test.use = "wilcox",
  min.pct = if (gene_qc) 0.01 else 0,
  logfc.threshold = if (gene_qc) 0.1 else 0,
  verbose = FALSE
)
fit$coef <- paste0("genotype", test_group)

t1 <- Sys.time()

## Statistics =================================================================
.extract_seurat <- fun_with_ram(extract_seurat)

# Keep only target covariate
txt <- "getting 'seurat wilcox' stats"
dge <- .extract_seurat(fit, info = txt)
dge <- dge[grepl(paste0("^", target_covariate), dge$coef), ]

# Add target covariate
dge[[target_covariate]] <- sub(pattern = paste0("^", target_covariate), replacement = "", dge$coef)

# Adding gene truth
dge <- merge(
  x = dge, 
  y = sce_rd, 
  by = "gene", 
  all.x = TRUE
)

sce_cols  <- setdiff(names(sce_rd), "gene")
na_sce_rd <- anyNA(dge[, sce_cols, drop = FALSE])
if (na_sce_rd) showf("Gene meta could not be assigned.", type = "stop")

# Add information
dge$formula     <- formula
dge$method      <- "wilcox"
dge$workflow    <- "singlecell"
dge$method_type <- "non-parametric"
dge$gene_qc     <- gene_qc
dge$time_min    <- as.numeric(difftime(t1, t0, units = "mins"))

## Save =======================================================================
.save_rds <- fun_with_ram(base::saveRDS)
txt       <- "saving 'wilcox' dge obj. (incl. gene truth)"

.save_rds(
  dge, 
  file = file.path(
    output_path,
    paste0(file_name, "_dge.rds")
  ),
  info = txt
)