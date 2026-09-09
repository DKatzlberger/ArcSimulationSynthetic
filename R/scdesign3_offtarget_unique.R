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
  family_use   <- "zinb"

  # Meta
  facLoc <- as.numeric(args[3])
  nOn    <- as.numeric(args[4])
  koEff  <- as.numeric(args[5])

  nBack  <- Inf

  file_name <- paste0("scdesign3_facLoc_", facLoc, "_nOn_", nOn, "_koEff_", koEff)

} else {
  project_root <- "/home/scc_e_928652/ArcSimulationSynthetic"
  input_path   <- "results/manual/splatter/splatter_facLoc_2.h5ad"
  output_path  <- "results/manual/ArcSimSyn_offtarget_unique/scdesign3"
  family_use   <- "zinb"

  # Meta
  facLoc <- 2
  nOn    <- 10
  koEff  <- 100

  nBack  <- 10

  file_name <- paste0("scdesign3_facLoc_", facLoc, "_nOn_", nOn, "_koEff_", koEff)
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
.read_h5ad <- fun_with_ram(anndataR::read_h5ad)

txt     <- sprintf("'%s' as 'SingleCellExperiment'", input_path)
sce_obj <- .read_h5ad(input_path, as = "SingleCellExperiment", info = txt)

# Create C, Y, K
sce_rd <- SingleCellExperiment::rowData(sce_obj)
sce_rd <- as.data.frame(sce_rd)

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
genes_downstream <- sce_rd$gene[sce_rd$status == "Downstream"]
downstream_reads <- sce_rd$total_reads[sce_rd$status == "Downstream"]
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

## scDesign3 ==================================================================
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


txt <- sprintf("on-target background: conditioned by sample info 'I'")
marg_on_background <- .estimate_marginals(
  sce_obj = sce_obj[genes_downstream, ],
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
marg_on_background  <- marg_on_background[sapply(marg_on_background, function(x) .is_valid_fit(x$fit))]
marg_off_downstream <- marg_off_downstream[sapply(marg_off_downstream, function(x) .is_valid_fit(x$fit))]
marg_off_background <- marg_off_background[sapply(marg_off_background, function(x) .is_valid_fit(x$fit))]
marg_background <- marg_background[sapply(marg_background, function(x) .is_valid_fit(x$fit))]

## Build complete marginal sets ===============================================
# On target
marg_on <- c(
  marg_on_downstream,
  marg_off_background,
  marg_background
)

# Off target
marg_off <- c(
  marg_on_background,
  marg_off_downstream,
  marg_background
)

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

new_cov_on <- modify_covariate(
  sce_obj = sce_obj[names(marg_on), ],
  covariates = c("C", "K", "I", "K_I"),
  corr_by = "K_I",       
  ko_eff = koEff
)

new_cov_off <- modify_covariate(
  sce_obj = sce_obj[names(marg_off), ],
  covariates = c("C", "K", "I", "K_I"),
  corr_by = "K_I",      
  ko_eff = koEff
)

new_cov_ntc <- modify_covariate(
  sce_obj = sce_obj[names(marg_on), ],
  covariates = c("C", "K", "I", "K_I"),
  corr_by = "K_I",      
  ko_eff = NULL
)

## Parameters ==================================================================
.extract_params <- fun_with_ram(extract_params)

params_tar_on <- extract_params(
  sce_obj = sce_obj[names(marg_on), ],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I"),
  new_covariate = new_cov_on[new_cov_on$K == "1", , drop = FALSE],
  marginals = marg_on,
  family_use = family_use,
  corr_by = "K_I"
)

params_tar_off <- extract_params(
  sce_obj = sce_obj[names(marg_off), ],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I"),
  new_covariate = new_cov_off[new_cov_off$K == "1", , drop = FALSE],
  marginals = marg_off,
  family_use = family_use,
  corr_by = "K_I"
)

params_ntc <- extract_params(
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
sim_tar_on_1 <- sim_pipe(
  sce_obj = sce_obj[names(marg_on), ],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I"),
  params = params_tar_on,
  copulas = copula,
  corr_by = "K_I",
  family_use = family_use,
  new_covariate = new_cov_on[new_cov_on$K == "1", , drop = FALSE],
  knockout = knockout
)

sim_tar_on_1$group     <- "targeting"
sim_tar_on_1$genotype  <- knockout
sim_tar_on_1$guide     <- paste0(knockout, "-On1")
sim_tar_on_1$on_target <- TRUE

# On-target guide 2
sim_tar_on_2 <- sim_pipe(
  sce_obj = sce_obj[names(marg_on), ],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I"),
  params = params_tar_on,
  copulas = copula,
  corr_by = "K_I",
  family_use = family_use,
  new_covariate = new_cov_on[new_cov_on$K == "1", , drop = FALSE],
  knockout = knockout
)

sim_tar_on_2$group     <- "targeting"
sim_tar_on_2$genotype  <- knockout
sim_tar_on_2$guide     <- paste0(knockout, "-On2")
sim_tar_on_2$on_target <- TRUE

# Off-target guide
sim_tar_off <- sim_pipe(
  sce_obj = sce_obj[names(marg_off), ],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I"),
  params = params_tar_off,
  copulas = copula,
  corr_by = "K_I",
  family_use = family_use,
  new_covariate = new_cov_off[new_cov_off$K == "1", , drop = FALSE],
  knockout = knockout
)

sim_tar_off$group     <- "targeting"
sim_tar_off$genotype  <- knockout
sim_tar_off$guide     <- paste0(knockout, "-Off1")
sim_tar_off$on_target <- FALSE

sim_ntc <- sim_pipe(
  sce_obj = sce_obj[names(marg_on), ],
  celltype = "celltype",
  covariates = c("C", "K", "I", "K_I"),
  params = params_ntc,
  copulas = copula,
  corr_by = "K_I",
  family_use = family_use,
  new_covariate = new_cov_ntc[new_cov_ntc$K == "0", , drop = FALSE],
  knockout = knockout
)

sim_ntc$group     <- "non-targeting"
sim_ntc$genotype  <- "NTC"
sim_ntc$guide     <- "NTC-1"
sim_ntc$on_target <- TRUE

## Observed annotations =======================================================
sce_obj <- BiocGenerics::cbind(sim_ntc, sim_tar_on_1, sim_tar_on_2, sim_tar_off)

sce_rd_new  <- sce_rd[rownames(sce_obj), , drop = FALSE]

sce_rd_new$sim_mode <- "unique guides"
sce_rd_new$nOn      <- nOn
sce_rd_new$koEff    <- koEff

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
status <- rd$status[row_ord]

true_logfc <- round(unique(rd$true_logFC), 2)
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
  sprintf("%s (n = %s)", g, table(sce_obj$genotype)[g])
}))

row_title <- sprintf(
  "Genes with a center log2FC of ±%s (facLoc %s) and\na knockout efficiencey of %s%%", 
  true_logfc, facLoc, koEff
)

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
ggplot2::ggsave(
  p,
  file = file.path(
    output_path,
    paste0(file_name, "_bx.svg")
  ),
  width = 10,
  height = 10,
  units = "cm"
)

## Demographics ===============================================================
p_ncells <- ggplot2::ggplot(
  data = as.data.frame(
    SummarizedExperiment::colData(sce_obj)
  ), 
  mapping = ggplot2::aes(
    x = sample, 
    fill = guide
  )
) +
ggplot2::geom_bar(
  position = "dodge", 
  width = 0.6
) +
ggplot2::scale_fill_manual(
  name = NULL,
  values = col_guides
) +
ggplot2::geom_text(
  mapping = ggplot2::aes(
    label = ggplot2::after_stat(count)
  ),
  position = ggplot2::position_dodge(
    width = 0.6
  ),
  stat = "count",
  angle = 90,    
  hjust = -0.2,  
  vjust = 0.5,  
  size = 1.5
) +
ggplot2::scale_y_continuous(
  expand = ggplot2::expansion(mult = c(0, 0.55)) 
) +
ggplot2::labs(
  x = "Sample",
  y = "Nr. of cells"
) +
theme_DK(
  lgd_position = "inside.top.right",
  show_grid = FALSE,
  show_border = FALSE
)

p_ngenes <- ggplot2::ggplot(
  data = as.data.frame(
    SummarizedExperiment::rowData(sce_obj)
  ), 
  mapping = ggplot2::aes(
    x = paste("sim. true logFC =", round(unique(rd$true_logFC))), 
    fill = factor(status, levels = c("Background", "On-target", "Off-target", "Knockout"))
  )
) +
ggplot2::geom_bar(
  mapping = ggplot2::aes(
    y = ggplot2::after_stat(count) / sum(ggplot2::after_stat(count))
  ),
  position = "dodge",
  width = 0.6
) +
ggplot2::scale_fill_manual(
  name = NULL,
  values = c(
    "Knockout"   = wesanderson::wes_palette("Chevalier1")[2],  
    "On-target"  = wesanderson::wes_palette("Chevalier1")[4],
    "Off-target" = wesanderson::wes_palette("Chevalier1")[1],
    "Background" = wesanderson::wes_palette("Chevalier1")[3]
  ),
) +
ggplot2::geom_text(
  mapping = ggplot2::aes(
    y = ggplot2::after_stat(count) / sum(ggplot2::after_stat(count)), 
    label = ggplot2::after_stat(count) 
  ),
  stat = "count",
  position = ggplot2::position_dodge(
    width = 0.6
  ),
  vjust = -0.5, 
  size = 1.5
) +
ggplot2::scale_y_continuous(
  expand = ggplot2::expansion(mult = c(0, 0.55)),
  labels = scales::percent
) +
ggplot2::labs(
  x = "Gene status",
  y = "Pct. of genes"

) + 
theme_DK(
  lgd_position = "inside.top.right",
  show_grid = FALSE,
  show_border = FALSE,
  show_ticks_x = FALSE
)

# Density
rd  <- SummarizedExperiment::rowData(sce_obj)

p_data <- data.frame(
  status = rep(rd$status, 2),
  logFC  = c(rd$sim_logFC, rd$obs_logFC),
  source = rep(c("Sim.", "Obs."), each = nrow(rd))
)
p_data <- p_data[p_data$status != "Knockout", ]

true_val       <- round(unique(rd$true_logFC), 2)
true_val_label <- paste(round(true_val))
vline_data <- data.frame(
  status = "On-target",
  x_pos  = c(true_val, -true_val)
)

p_density <- ggplot2::ggplot(
  data = p_data, 
  mapping = ggplot2::aes(
    x = logFC,
    color = source,
    fill = source
  )
) +
ggplot2::geom_density(
  alpha = 0.5
) + 
ggplot2::scale_fill_manual(
  name = "log2FC",
  values = c(
    "Sim." = "firebrick3", 
    "Obs." = "grey50"
  )
) +
ggplot2::scale_color_manual(
  name = "log2FC",
  values = c(
    "Sim." = "firebrick3", 
    "Obs." = "grey50"
  )
) +
ggplot2::geom_vline(
  data = vline_data,
  mapping = ggplot2::aes(
    xintercept = x_pos, 
    linetype = true_val_label
  ),
  color = "black",        
  linewidth = 0.3,
  inherit.aes = FALSE    
) +
ggplot2::scale_linetype_manual(
  name = "Center",
  values = setNames("dashed", true_val_label) 
) +
ggplot2::facet_wrap(
  ~ status, 
  scales = "free_y"
) +
ggplot2::labs(
  x = "log2FC",
  y = "Density",
  fill = "Source"
) +
theme_DK(
  lgd_position = "right",
  show_grid = FALSE
)

# UMAP
num_neighbors <- min(15, floor(ncol(sce_obj) / 2))
sce_obj <- scater::runUMAP(sce_obj, n_neighbors = num_neighbors)

p_data <- SingleCellExperiment::reducedDim(sce_obj, "UMAP")
colnames(p_data) <- c("UMAP1", "UMAP2")
p_data <- cbind(p_data, SingleCellExperiment::colData(sce_obj))
p_data <- as.data.frame(p_data)

p_sample <- ggplot2::ggplot(
  data = p_data, 
  mapping = ggplot2::aes(
    x = UMAP1, 
    y = UMAP2, 
    color = sample
  )
) +
ggplot2::geom_point(
  size = 0.1
) +
ggplot2::scale_color_manual(
  values = c(
    "mouse1" = wesanderson::wes_palette("Rushmore1")[3],
    "mouse2" = wesanderson::wes_palette("Rushmore1")[4],
    "mouse3" = wesanderson::wes_palette("Rushmore1")[5]
  )
) +
ggplot2::scale_x_continuous(
  labels = scales::label_number(accuracy = 1)
) +
ggplot2::scale_y_continuous(
  labels = scales::label_number(accuracy = 1)
) +
ggplot2::facet_wrap(
  ~ "Sample", 
  scales = "free_y"
) +
ggplot2::labs(
  x = "UMAP1",
  y = "UMAP2",
  color = "Sample"
) + 
theme_DK(
  lgd_position = "inside.top.right",
  show_grid = FALSE,
  aspect.ratio = 1
)

p_genotype <- ggplot2::ggplot(
  data = p_data, 
  mapping = ggplot2::aes(
    x = UMAP1, 
    y = UMAP2, 
    color = genotype
  )
) +
ggplot2::geom_point(
  size = 0.1
) +
ggplot2::scale_x_continuous(
  labels = scales::label_number(accuracy = 1)
) +
ggplot2::scale_y_continuous(
  labels = scales::label_number(accuracy = 1)
) +
ggplot2::facet_wrap(
  ~ "Genotype", 
  scales = "free_y"
) +
ggplot2::labs(
  x = "UMAP1",
  y = "UMAP2",
  color = "Genotype"
) + 
theme_DK(
  lgd_position = "inside.top.right",
  show_grid = FALSE,
  aspect.ratio = 1
)

p_guide <- ggplot2::ggplot(
  data = p_data, 
  mapping = ggplot2::aes(
    x = UMAP1, 
    y = UMAP2, 
    color = guide
  )
) +
ggplot2::geom_point(
  size = 0.1
) +
ggplot2::scale_color_manual(
  values = col_guides
) + 
ggplot2::scale_x_continuous(
  labels = scales::label_number(accuracy = 1)
) +
ggplot2::scale_y_continuous(
  labels = scales::label_number(accuracy = 1)
) +
ggplot2::facet_wrap(
  ~ "Guides", 
  scales = "free_y"
) +
ggplot2::labs(
  x = "UMAP1",
  y = "UMAP2",
  color = "Guides"
) + 
theme_DK(
  lgd_position = "inside.top.right",
  show_grid = FALSE,
  aspect.ratio = 1
)

p_on_target <- ggplot2::ggplot(
  data = p_data, 
  mapping = ggplot2::aes(
    x = UMAP1, 
    y = UMAP2, 
    color = on_target
  )
) +
ggplot2::geom_point(
  size = 0.1
) +
ggplot2::scale_x_continuous(
  labels = scales::label_number(accuracy = 1)
) +
ggplot2::scale_y_continuous(
  labels = scales::label_number(accuracy = 1)
) +
ggplot2::scale_color_manual(
  values = c(
    "TRUE"  = wesanderson::wes_palette("Moonrise3")[1],
    "FALSE" = wesanderson::wes_palette("Moonrise3")[2]
  )
) +
ggplot2::facet_wrap(
  ~ "On-target", 
  scales = "free_y"
) +
ggplot2::labs(
  x = "UMAP1",
  y = "UMAP2",
  color = "On-target"
) + 
theme_DK(
  lgd_position = "inside.top.right",
  show_grid = FALSE,
  aspect.ratio = 1
)

p_umap <- patchwork::wrap_plots(
  p_sample,
  p_genotype,
  p_on_target,
  p_guide
) + patchwork::plot_layout(
  guides = "collect"
) & 
ggplot2::guides(
  color = ggplot2::guide_legend(
    override.aes = list(size = 0.5)
  ) 
)

p <- p_ncells + p_umap + p_ngenes + p_density +
  patchwork::plot_layout(
    design = "
      AB
      CD
    ",
    widths = c(0.5, 1) 
  )

# Save
ggplot2::ggsave(
  p,
  file = file.path(
    output_path,
    paste0(file_name, "_um.svg")
  ),
  width = 12,
  height = 10,
  units = "cm"
)