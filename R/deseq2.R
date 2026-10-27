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
  input_path   <- "results/manual/ArcSimSyn_offtarget_mixed/pseudobulk/celltype_Gran. P_facLoc_2_nOn_10_koEff_100_nGuide_300.h5ad"
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
    "DESeq2"
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

# Model matrix
mm <- model.matrix(fixed_formula, data = meta)

# N cells per group
groups <- meta[[target_covariate]]
target <- setdiff(levels(groups), ref)
stopifnot(length(target) == 1)

sce_rd$n_cells_target <- sum(groups == target, na.rm = TRUE)
sce_rd$n_cells_ref    <- sum(groups == ref, na.rm = TRUE)

## QC: filter genes ===========================================================
.filterByExpr <- fun_with_ram(edgeR::filterByExpr)

if (gene_qc) {
  dge  <- edgeR::DGEList(matr)
  keep <- .filterByExpr(dge, design = mm)
  matr <- matr[keep, ]
}

## Fit ========================================================================
.DESeqDataSetFromMatrix <- fun_with_ram(DESeq2::DESeqDataSetFromMatrix)
.estimateSizeFactors    <- fun_with_ram(DESeq2::estimateSizeFactors)
.DESeq                  <- fun_with_ram(DESeq2::DESeq)

t0 <- Sys.time()
fit <- DESeq2::DESeqDataSetFromMatrix(countData = matr, colData = meta, design = as.formula(formula))
fit <- DESeq2::estimateSizeFactors(fit, type = "poscounts")
fit <- DESeq2::DESeq(fit)
t1 <- Sys.time()

## Statistics =================================================================
.extract_deseq2 <- fun_with_ram(extract_deseq2)

# Keep only target covariate
txt <- "getting 'deseq2' stats"
dge <- .extract_deseq2(fit, gene_qc = gene_qc, info = txt)
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
dge$method      <- "deseq2"
dge$workflow    <- "pseudobulk"
dge$method_type <- "parametric"
dge$gene_qc     <- gene_qc
dge$time_min    <- as.numeric(difftime(t1, t0, units = "mins"))

## Save =======================================================================
.save_rds <- fun_with_ram(base::saveRDS)
txt      <- "saving 'deseq2' dge obj. (incl. gene truth)"

.save_rds(
  dge, 
  file = file.path(
    output_path,
    paste0(file_name, "_dge.rds")
  ),
  info = txt
)
