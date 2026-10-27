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
  full_args  <- commandArgs(trailingOnly = FALSE)
  trail_args <- commandArgs(trailingOnly = TRUE)

  # Commands
  script_path <- sub("^--file=", "", grep("^--file=", full_args, value = TRUE)[1])
  script_name <- tools::file_path_sans_ext(basename(script_path))

  project_root <- trail_args[1]
  input_path   <- trail_args[2]
  output_path  <- "."

  # Meta
  # celltype <- as.character(args[3])
  # facLoc   <- as.numeric(args[4])

  # Params
  nOn    <- as.numeric(trail_args[5])
  koEff  <- as.numeric(trail_args[6])

  family_use <- "zinb"
  nBack      <- Inf

  # Output file
  input_name <- tools::file_path_sans_ext(basename(input_path))
  file_name  <- paste0(input_name , "_nOn_", nOn, "_koEff_", koEff)

} else {
  # Commands
  script_path  <- rstudioapi::getActiveDocumentContext()$path
  script_name  <- tools::file_path_sans_ext(basename(script_path))

  project_root <- "/home/scc_e_928652/ArcSimulationSynthetic"
  input_path   <- "results/manual/ArcSimSyn_offtarget_mixed/splatter/celltype_Gran. P_facLoc_2.h5ad"
  output_path  <- file.path("results/manual/ArcSimSyn_offtarget_mixed", script_name)
  stopifnot(dirname(dirname(input_path)) == dirname(output_path))

  # Meta 
  # celltype <- "Gran. P"
  # facLoc   <- 2

  # Params
  nOn    <- 10
  koEff  <- 100

  family_use <- "zinb"
  nBack      <- 10

  # Output file
  input_name <- tools::file_path_sans_ext(basename(input_path))
  file_name  <- paste0(input_name , "_nOn_", nOn, "_koEff_", koEff)
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

## scDesign3 ==================================================================
knockout <- sce_rd$gene[sce_rd$status == "Knockout"]

counts <- SingleCellExperiment::counts(sce_obj)
Y <- counts[rownames(counts) != knockout, , drop = FALSE]
C <- counts[knockout, ]
K <- factor(ifelse(sce_obj$genotype == knockout, "1", "0"), levels = c("0", "1"))
I <- sce_obj$sample

# Interaction between knockout and sample
K_I <- factor(paste0("K", K, "_I", I))

sce_obj <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = Y),
  colData = data.frame(
    C = C,
    K = K,
    I = I,
    K_I = K_I,
    celltype = "macrophages",
    row.names = colnames(Y)
  )
)

# Downstream genes
genes_downstream  <- sce_rd$gene[sce_rd$status == "Downstream"]
downstream_reads  <- sce_rd$total_reads[sce_rd$status == "Downstream"]
sorted_downstream <- genes_downstream[order(downstream_reads, decreasing = TRUE)]

# On- & off-target
genes_downstream <- sorted_downstream[1:nOn]
genes_off_target <- sorted_downstream[(nOn + 1):(nOn + nOn)]
stopifnot(!any(genes_downstream %in% genes_off_target))

# Background
genes_background <- sce_rd$gene[sce_rd$status == "Background"]
background_reads <- sce_rd$total_reads[sce_rd$status == "Background"]
genes_background <- genes_background[order(background_reads, decreasing = TRUE)]

genes_background <- genes_background[1:min(nBack, length(genes_background))]

## Marginals ==================================================================
.estimate_marginals <- fun_with_ram(estimate_marginals)

txt <- sprintf("on-target downstream: conditioned by %s expression 'C' and sample info 'I'", knockout)
marg_on_downstream <- .estimate_marginals(
  sce_obj = sce_obj[genes_downstream, ],
  assay_use = "counts",
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I"),
  mu_formula = "C + I",
  sigma_formula = "1",
  family_use = family_use,
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
  covariates = c("C", "K", "I", "K_I"),
  mu_formula = "C + I",
  sigma_formula = "1",
  family_use = family_use,
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
  covariates = c("C", "K", "I", "K_I"),
  mu_formula = "I",
  sigma_formula = "1",
  family_use = family_use,
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
  covariates = c("C", "K", "I", "K_I"),
  mu_formula = "I",
  sigma_formula = "1",
  family_use = family_use,
  corr_by = "K_I",
  n_cores = 1,
  verbose = FALSE,
  info = txt
)

# Validate marginals
marg_on_downstream  <- marg_on_downstream[sapply(marg_on_downstream, function(x) .is_valid_fit(x$fit))]
marg_off_downstream <- marg_off_downstream[sapply(marg_off_downstream, function(x) .is_valid_fit(x$fit))]
marg_off_background <- marg_off_background[sapply(marg_off_background, function(x) .is_valid_fit(x$fit))]
marg_background     <- marg_background[sapply(marg_background, function(x) .is_valid_fit(x$fit))]

## Build complete marginal sets ===============================================
marg_on  <- c(marg_on_downstream, marg_off_background, marg_background)
marg_off <- c(marg_on_downstream, marg_off_downstream, marg_background)

## Correlation ================================================================
.estimate_correlation <- fun_with_ram(estimate_correlation)

# Estimate from the real downstream effect
copula <- .estimate_correlation(
  sce_obj = sce_obj[names(marg_on), ],
  marginals = marg_on,
  covariates = c("C", "K", "I", "K_I"),
  corr_by = "K_I",
  family_use = family_use,
  copula = "gaussian",
  n_cores = 1,
  verbose = FALSE
)

## Modify C (intervention) ====================================================
.modify_covariate <- fun_with_ram(modify_covariate)

txt <- "modifying cells for on-target effect"
new_cov_on <- .modify_covariate(
  sce_obj = sce_obj[names(marg_on), ],
  covariates = c("C", "K", "I", "K_I"),
  corr_by = "K_I",       
  ko_eff = koEff,
  mode = "partial",
  info = txt
)

txt <- "modifying cells for off-target effect"
new_cov_off <- .modify_covariate(
  sce_obj = sce_obj[names(marg_off), ],
  covariates = c("C", "K", "I", "K_I"),
  corr_by = "K_I",      
  ko_eff = koEff,
  mode = "partial",
  info = txt
)

txt <- "modifying cells for non-targeting control"
new_cov_ntc <- .modify_covariate(
  sce_obj = sce_obj[names(marg_on), ],
  covariates = c("C", "K", "I", "K_I"),
  corr_by = "K_I",      
  ko_eff = NULL,
  mode = "partial",
  info = txt
)

## Parameters ==================================================================
.extract_params <- fun_with_ram(extract_params)

params_tar_on <- .extract_params(
  sce_obj = sce_obj[names(marg_on), ],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I"),
  new_covariate = new_cov_on[new_cov_on$K == "1", , drop = FALSE],
  marginals = marg_on,
  family_use = family_use,
  corr_by = "K_I"
)

params_tar_off <- .extract_params(
  sce_obj = sce_obj[names(marg_off), ],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I"),
  new_covariate = new_cov_off[new_cov_off$K == "1", , drop = FALSE],
  marginals = marg_off,
  family_use = family_use,
  corr_by = "K_I"
)

params_ntc <- .extract_params(
  sce_obj = sce_obj[names(marg_on), ],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I"),
  new_covariate = new_cov_ntc[new_cov_ntc$K == "0", , drop = FALSE],
  marginals = marg_on,
  family_use = family_use,
  corr_by = "K_I"
)

## Simulate counts ============================================================
.sim_pipe <- fun_with_ram(sim_pipe)

# On-target guide 1
txt <- "simulating cells for on-target guide 1"
sim_tar_on_1 <- .sim_pipe(
  sce_obj = sce_obj[names(marg_on), ],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I"),
  params = params_tar_on,
  copulas = copula,
  corr_by = "K_I",
  family_use = family_use,
  new_covariate = new_cov_on[new_cov_on$K == "1", , drop = FALSE],
  knockout = knockout,
  info = txt
)

sim_tar_on_1$group     <- "targeting"
sim_tar_on_1$genotype  <- knockout
sim_tar_on_1$guide     <- paste0(knockout, "-On1")
sim_tar_on_1$on_target <- TRUE

# On-target guide 2
txt <- "simulating cells for on-target guide 2"
sim_tar_on_2 <- .sim_pipe(
  sce_obj = sce_obj[names(marg_on), ],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I"),
  params = params_tar_on,
  copulas = copula,
  corr_by = "K_I",
  family_use = family_use,
  new_covariate = new_cov_on[new_cov_on$K == "1", , drop = FALSE],
  knockout = knockout,
  info = txt
)

sim_tar_on_2$group     <- "targeting"
sim_tar_on_2$genotype  <- knockout
sim_tar_on_2$guide     <- paste0(knockout, "-On2")
sim_tar_on_2$on_target <- TRUE

# Off-target guide
txt <- "simulating cells for off-target guide 1"
sim_tar_off <- .sim_pipe(
  sce_obj = sce_obj[names(marg_off), ],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I"),
  params = params_tar_off,
  copulas = copula,
  corr_by = "K_I",
  family_use = family_use,
  new_covariate = new_cov_off[new_cov_off$K == "1", , drop = FALSE],
  knockout = knockout,
  info = txt
)

sim_tar_off$group     <- "targeting"
sim_tar_off$genotype  <- knockout
sim_tar_off$guide     <- paste0(knockout, "-Off1")
sim_tar_off$on_target <- FALSE

# Non-targeting guide
txt <- "simulating cells for NTC guide 1"
sim_ntc <- .sim_pipe(
  sce_obj = sce_obj[names(marg_on), ],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I"),
  params = params_ntc,
  copulas = copula,
  corr_by = "K_I",
  family_use = family_use,
  new_covariate = new_cov_ntc[new_cov_ntc$K == "0", , drop = FALSE],
  knockout = knockout,
  info = txt
)

sim_ntc$group     <- "non-targeting"
sim_ntc$genotype  <- "NTC"
sim_ntc$guide     <- "NTC-1"
sim_ntc$on_target <- TRUE

## Observed annotations =======================================================
sce_obj <- BiocGenerics::cbind(sim_ntc, sim_tar_on_1, sim_tar_on_2, sim_tar_off)

sce_rd_new <- sce_rd[rownames(sce_obj), , drop = FALSE]

sce_rd_new$guide_mode <- "mixed guides"
sce_rd_new$nOn        <- nOn
sce_rd_new$koEff      <- koEff

# Log counts
sce_obj <- scater::logNormCounts(sce_obj)

cells_ntc <- which(sce_obj$group == "non-targeting")
cells_tar <- which(sce_obj$group == "targeting")

log_counts <- SummarizedExperiment::assay(sce_obj, "logcounts")
obs_logFC  <- Matrix::rowMeans(log_counts[, cells_tar]) - Matrix::rowMeans(log_counts[, cells_ntc])

sce_rd_new$obs_logFC   <- obs_logFC
sce_rd_new$total_reads <- Matrix::rowSums(SingleCellExperiment::counts(sce_obj))

## DEG status =================================================================
sce_rd_new$status                                              <- "Background"
sce_rd_new$status[rownames(sce_rd_new) == knockout]            <- "Knockout"
sce_rd_new$status[rownames(sce_rd_new) %in% genes_downstream]  <- "On-target"
sce_rd_new$status[rownames(sce_rd_new) %in% genes_off_target]  <- "Off-target"

# Assign back rowData
SummarizedExperiment::rowData(sce_obj) <- NULL
SummarizedExperiment::rowData(sce_obj) <- sce_rd_new
sce_obj$sizeFactor <- NULL
sce_obj$ExpLibSize <- NULL

## Save sim. ==================================================================
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

## Heatmap ====================================================================
rd <- SummarizedExperiment::rowData(sce_obj)

# Col order
guide_levels <- c(
  sort(grep("-Off", unique(sce_obj$guide), value = TRUE)),
  sort(grep("-On",  unique(sce_obj$guide), value = TRUE)),
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

top_ontarget   <- on_idx[order(rd$total_reads[on_idx], decreasing = TRUE)][1:pmin(length(on_idx), 15)]
top_offtarget  <- off_idx[order(rd$total_reads[off_idx], decreasing = TRUE)][1:pmin(length(off_idx), 15)]
top_background <- back_idx[order(rd$total_reads[back_idx], decreasing = TRUE)][1:pmin(length(back_idx), 9)]

top_ontarget   <- top_ontarget[order(match(top_ontarget, on_idx))]
top_offtarget  <- top_offtarget[order(match(top_offtarget, off_idx))]
top_background <- top_background[order(match(top_background, back_idx))]

row_ord <- c(ko_idx, top_ontarget, top_offtarget, top_background)

# Matrix
mat <- as.matrix(SummarizedExperiment::assay(sce_obj, "logcounts"))
mat <- mat[row_ord, col_ord, drop = FALSE]

col_mat <- circlize::colorRamp2(
  c(0, max(mat, na.rm = TRUE) * 0.4, max(mat, na.rm = TRUE)), 
  c("#000000", "#A52A2A", "#FFD700") 
)

# Top column annotations
col_guides <- setNames(
  wesanderson::wes_palette(
    "Cavalcanti1", 
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
      "targeting"     = wesanderson::wes_palette("Moonrise2")[3],
      "non-targeting" = wesanderson::wes_palette("Moonrise2")[1]
    )
  ),
  which = "column",
  annotation_name_side = "right",
  show_annotation_name = TRUE,
  annotation_name_rot = c(
    `Sample`   = 0,     
    `Group`    = 0,    
    `Lib size` = 270   
  )
)

# Left Row Annotations
status     <- rd$status[row_ord]
sim_logfc  <- rd$sim_logFC[row_ord]
obs_logfc  <- rd$obs_logFC[row_ord]

max_logfc <- pmax(abs(sim_logfc), abs(obs_logfc), na.rm = TRUE)
max_logfc <- max(max_logfc, na.rm = TRUE)
col_logfc <- circlize::colorRamp2(
  c(-max_logfc, 0, max_logfc), 
  c("#313695", "#F7F7F7", "#B2182B")
)

left_row_anno <- cmplx_anno(
  `Gene status`  = status,
  `Sim. log2FC`  = sim_logfc,
  `Obs. log2FC`  = obs_logfc,
  col = list(
    `Gene status` = c(
      "Knockout"   = wesanderson::wes_palette("Chevalier1")[2],  
      "On-target"  = wesanderson::wes_palette("Chevalier1")[4],
      "Off-target" = wesanderson::wes_palette("Chevalier1")[1],
      "Background" = wesanderson::wes_palette("Chevalier1")[3]
    ),
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

row_title <- sprintf("Genes with a facLoc %s and\na knockout efficiencey of %s%%", facLoc, koEff)

# Heatmap 
ht <- cmplx_ht(
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

# Save
.cmplx_save <- fun_with_ram(cmplx_save)
txt         <- sprintf("saved as %s", paste0(file_name, "_ht.svg"))

.cmplx_save(
  ht,
  file = file.path(
    output_path,
    paste0(file_name, "_ht.svg")
  ),
  width = 15,
  height = 10,
  units = "cm",
  info = txt
)

## Boxplots ===================================================================
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

top_ontarget   <- on_idx[order(rd$total_reads[on_idx], decreasing = TRUE)][1:pmin(length(on_idx), 15)]
top_offtarget  <- off_idx[order(rd$total_reads[off_idx], decreasing = TRUE)][1:pmin(length(off_idx), 15)]
top_background <- back_idx[order(rd$total_reads[back_idx], decreasing = TRUE)][1:pmin(length(back_idx), 9)]

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
  sort(grep("-Off", unique(p_data$guide), value = TRUE)),
  sort(grep("-On",  unique(p_data$guide), value = TRUE)),
  sort(grep("NTC",  unique(p_data$guide), value = TRUE))
)
p_data$guide <- factor(p_data$guide, levels = guide_levels)

# Background cols
legend_cols <- c(
  "Knockout"   = wesanderson::wes_palette("Chevalier1")[2],  
  "On-target"  = wesanderson::wes_palette("Chevalier1")[4],  
  "Off-target" = wesanderson::wes_palette("Chevalier1")[1],
  "Background" = wesanderson::wes_palette("Chevalier1")[3]
)

bg_cols_with_alpha <- paste0(legend_cols, "40")
names(bg_cols_with_alpha) <- names(legend_cols)

# Plot
p <- ggplot2::ggplot() +
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
  name = "Gene status",
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

# Save
.ggsave <- fun_with_ram(ggplot2::ggsave)
txt     <- sprintf("saved as %s", paste0(file_name, "_bx.svg"))

.ggsave(
  p,
  file = file.path(
    output_path,
    paste0(file_name, "_bx.svg")
  ),
  width = 10,
  height = 10,
  units = "cm",
  info = txt
)
