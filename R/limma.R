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
    "limma"
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
mm  <- model.matrix(fixed_formula, data = meta)

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
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
} else {
  dge <- edgeR::DGEList(matr)
}

## Fit ========================================================================
t0 <- Sys.time()
dge  <- edgeR::calcNormFactors(dge)
vobj <- limma::voom(dge, design = mm, plot = FALSE)
fit  <- limma::lmFit(vobj, design = mm)
fit  <- limma::eBayes(fit)
t1 <- Sys.time()

## Statistics =================================================================
.extract_limma <- fun_with_ram(extract_limma)

# Keep only target covariate
txt <- "getting 'limma' stats"
dge <- .extract_limma(fit, info = txt)
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
dge$method      <- "limma"
dge$workflow    <- "pseudobulk"
dge$method_type <- "parametric"
dge$gene_qc     <- gene_qc
dge$time_min    <- as.numeric(difftime(t1, t0, units = "mins"))

## Save =======================================================================
.save_rds <- fun_with_ram(base::saveRDS)
txt      <- "saving 'limma' dge obj. (incl. gene truth)"

.save_rds(
  dge, 
  file = file.path(
    output_path,
    paste0(file_name, "_dge.rds")
  ),
  info = txt
)


## Aggregation ================================================================
# if (target_covariate == "guide") {

#   # Filter for target covariate
#   pattern  <- paste0("^", target_covariate)
#   agg_data <- dge[grepl(pattern, dge$coef), ]

#   # This will be a fun!
#   if (length(unique(agg_data$coef)) > 1) {

#     # Transcriptome wide impact
#     guides   <- unique(agg_data$coef)
#     n_guides <- length(guides)

#     # Correlation (sampling, biology)
#     cor_sampling <- matrix(0, nrow = n_guides, ncol = n_guides, dimnames = list(guides, guides))
#     diag(cor_sampling) <- 1.0

#     cor_biology <- matrix(NA, nrow = n_guides, ncol = n_guides, dimnames = list(guides, guides))
#     diag(cor_biology) <- 1.0

#     for (i in 1:(length(guides) - 1)) {
#       for (j in (i + 1):length(guides)) {

#         guide_a <- guides[i]
#         guide_b <- guides[j]

#         results1 <- prepare_trade_input(agg_data, guide_a)
#         results2 <- prepare_trade_input(agg_data, guide_b)

#         trade_bivariate <- TRADEtools::TRADE(
#           mode = "bivariate",
#           results1 = results1,
#           results2 = results2,
#           genes_exclude = NULL,
#           estimate_sampling_covariance = TRUE,
#           covariance_matrix_set = "combined",
#           component_varexplained_threshold = 0,
#           weight_nocorr = 1,
#           n_sample = NULL
#         )

#         # Sampling correlation
#         trade_corr_samp <- trade_bivariate$V
#         pairwise_samp   <- trade_corr_samp[1, 2]

#         cor_sampling[guide_a, guide_b] <- pairwise_samp
#         cor_sampling[guide_b, guide_a] <- pairwise_samp

#         # Biological correlation
#         trade_corr_bio <- trade_bivariate$correlation_matrix
#         pairwise_bio   <- trade_corr_bio[1, 2]

#         cor_biology[guide_a, guide_b] <- pairwise_bio
#         cor_biology[guide_b, guide_a] <- pairwise_bio
#       }
#     }

#     # Sort the data to fit to cor matrix
#     agg_data <- agg_data[order(agg_data$gene, agg_data$coef), ]

#     # Inverse variance meta
#     gene_list <- split(agg_data, agg_data$gene)
#     agg_list  <- lapply(gene_list, function(x) {

#       n_estimates <- nrow(x)
#       guide_ids   <- x$coef

#       stopifnot(length(guide_ids) == n_estimates)
      
#       # Stats
#       R_bio  <- cor_biology[guide_ids, guide_ids, drop = FALSE]
#       R_samp <- cor_sampling[guide_ids, guide_ids, drop = FALSE]

#       logfc_vec <- x$logFC
#       se_vec    <- x$se

#       # Biological reliability
#       if (n_estimates == 1) {

#           reliability <- 1
#       } else {

#         diag(R_bio) <- NA
#         guide_concordance <- apply(R_bio, 1, median, na.rm = TRUE)
#         reliability       <- pmax(guide_concordance, 0.05)
#       }

#       # Sampling quality
#       D <- diag(se_vec, nrow = n_estimates)
#       Sigma <- D %*% R_samp %*% D

#       scale <- diag(1 / sqrt(reliability))
#       Sigma_weighted <- scale %*% Sigma %*% scale

#       Sigma_inv <- chol2inv(chol(Sigma_weighted))
#       X         <- rep(1, n_estimates)

#       sum_weights <- as.numeric(crossprod(X, Sigma_inv %*% X))
#       gls_weights <- as.numeric(Sigma_inv %*% X) / sum_weights

#       # Meta stats
#       meta_logFC <- sum(gls_weights * logfc_vec)
#       meta_se    <- sqrt(as.numeric(crossprod(gls_weights, Sigma %*% gls_weights)))

#       meta_z  <- meta_logFC / meta_se
#       meta_p  <- 2 * pnorm(-abs(meta_z))

#       # Results
#       res <- x[1, , drop = FALSE]

#       coef <- gsub(paste0("^", target_covariate), "", res$coef[1])
#       coef <- gsub("-.*$", "", coef)
      
#       res$coef <- paste0("genotype", coef)
#       res$logFC  <- meta_logFC
#       res$se     <- meta_se
#       res$pvalue <- meta_p

#       res
#     })

#     # Combine across genes
#     dge_agg <- do.call(rbind, agg_list)
#     rownames(dge_agg) <- NULL

#     # Multiple test correction
#     dge_agg$FDR <- p.adjust(dge_agg$pvalue, method = "BH")
#   }

#   dge <- dge_agg
# }
# head(dge)

