#!/usr/bin/env Rscript

## max. code width ============================================================

## Title: Arc Institute in vivo pertubation screens
## Author: DKatzlberger
## Date: 21.01.2026

## This script generates a synthetic Perturb-seq dataset using real
## non-targeting control (NTC) cells as the basis for the simulation.
## The library 'splatter' is used to learn the distribution of the data.

## The input Seurat object is filtered for NTC cells with 1,000-5,000
## detected genes, followed by normalization and selection of the 3,000
## most highly variable genes.

## The NTC expression data are converted to a SingleCellExperiment and
## used to estimate Splatter parameters. Splatter then simulates 80%
## non-targeting and 20% targeting cells across three batches.

## The parameter 'facLoc' controls the magnitude of simulated
## differential expression, with 20% of genes designated as DEGs.
## In addition, one highly expressed gene is selected as an explicit
## knockout. Its expression in targeting cells is reduced to 1% of
## the mean expression observed in NTC cells.

## Ground-truth, simulated, and observed log fold changes are stored
## for each gene. Genes are classified as Background, Downstream, or
## Knockout according to their simulated status.

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
  facLoc <-  as.numeric(args[3])

  file_name <- paste0("splatter_facLoc_", facLoc)

} else {
  project_root <- "/home/scc_e_928652/ArcSimulationSynthetic"
  input_path   <- "data/JB22_PM/rna_integrated.rds"
  output_path  <- "results/manual/splatter"

  # Meta
  facLoc <- 2

  file_name <- paste0("splatter_facLoc_", facLoc)
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
    "splatter"
  )
)

## Dataset ====================================================================
.readRDS <- fun_with_ram(base::readRDS)

txt      <- sprintf("'%s' as 'Seurat'", input_path)
seu_obj  <- .readRDS(input_path, info = txt)

# Genotype
.assign_guide_by_max_count <- fun_with_ram(assign_guide_by_max_count)
seu_obj <- .assign_guide_by_max_count(seu_obj)

# Select cells
seu_obj <- subset(seu_obj, subset = genotype == "NTC")
seu_obj <- subset(seu_obj, subset = nFeature_RNA >= 1000 & nFeature_RNA <= 5000)

# Select 3000 hvgs
seu_obj <- Seurat::NormalizeData(seu_obj, verbose = FALSE)
seu_obj <- Seurat::FindVariableFeatures(seu_obj, nfeatures = 3000, verbose = FALSE)

seu_hvg <- Seurat::VariableFeatures(seu_obj)
seu_obj <- seu_obj[seu_hvg, ]

# Single cell obj.
sce_obj <- Seurat::as.SingleCellExperiment(seu_obj, assay = "RNA")
SingleCellExperiment::counts(sce_obj) <- as.matrix(SingleCellExperiment::counts(sce_obj))

## Splatter ===================================================================
.splatEstimate <- fun_with_ram(splatter::splatEstimate)

params <- .splatEstimate(sce_obj)
ncells <- splatter::getParams(params, "nCells")$nCells

# Set sim. parameters
params <- splatter::setParams(
  params,
  seed = seed,
  # DEGs
  group.prob = c(0.8, 0.2),
  de.prob = 0.2,
  de.facLoc = facLoc,
  de.facScale = 1,

  # Samples
  batchCells = round(ncells * c(0.325, 0.35, 0.325)),
  batch.facLoc = 0.1,
  batch.facScale = 0.1
)

# Simulate cells
.splatSimulate <- fun_with_ram(splatter::splatSimulate)

sce_obj <- .splatSimulate(
  params,
  method = "groups",
  verbose = FALSE
)

# Remove hidden layers
S4Vectors::metadata(sce_obj)$Params             <- NULL
SummarizedExperiment::assay(sce_obj, "Dropout") <- NULL

# Rename columns
sce_obj$group <- ifelse(sce_obj$Group == "Group1", "non-targeting", "targeting")
sce_obj$group <- factor(sce_obj$group, levels = c("targeting", "non-targeting"))
sce_obj$Group <- NULL

sce_obj$sample <- factor(gsub("Batch", "mouse", sce_obj$Batch), levels = c("mouse1", "mouse2", "mouse3"))
sce_obj$Batch  <- NULL

## Knockout ===================================================================
cells_ntc <- which(sce_obj$group == "non-targeting")
cells_tar <- which(sce_obj$group == "targeting")

reads_ntc <- Matrix::rowSums(SingleCellExperiment::counts(sce_obj)[, cells_ntc, drop = FALSE])
knockout  <- which.max(reads_ntc)

# Mean (knockout)
mu_ntc <- mean(SingleCellExperiment::counts(sce_obj)[knockout, cells_ntc])
mu_tar <- mu_ntc * 0.01
disp   <- splatter::getParam(params, "bcv.common")

SingleCellExperiment::counts(sce_obj)[knockout, cells_tar] <- rnbinom(length(cells_tar), mu = mu_tar, size = 1/disp^2)

# Genotype
sce_obj$genotype <- ifelse(sce_obj$group == "non-targeting", "NTC", rownames(sce_obj)[knockout])

## Ground-truth annotations ===================================================
sce_rd_old <- SummarizedExperiment::rowData(sce_obj)
sce_rd_new <- S4Vectors::DataFrame(row.names = rownames(sce_obj))

sce_rd_new$gene       <- sce_rd_old$Gene
sce_rd_new$facLoc     <- facLoc
sce_rd_new$true_logFC <- facLoc * log2(exp(1))
sce_rd_new$sim_logFC  <- log2(sce_rd_old$DEFacGroup2 / sce_rd_old$DEFacGroup1)

sce_rd_new$sim_logFC[knockout] <- log2(mu_tar / (sce_rd_old$GeneMean[knockout] * sce_rd_old$DEFacGroup1[knockout]))

## Observed annotations =======================================================
sce_obj <- scater::logNormCounts(sce_obj)

log_counts <- SummarizedExperiment::assay(sce_obj, "logcounts")
obs_logFC  <- Matrix::rowMeans(log_counts[, cells_tar]) - Matrix::rowMeans(log_counts[, cells_ntc])

sce_rd_new$obs_logFC   <- obs_logFC
sce_rd_new$total_reads <- Matrix::rowSums(SingleCellExperiment::counts(sce_obj))

## DEG status =================================================================
sce_rd_new$status                             <- "Background"
sce_rd_new$status[sce_rd_new$sim_logFC != 0]  <- "Downstream"
sce_rd_new$status[knockout]                   <- "Knockout"

# Assign back rowData
SummarizedExperiment::rowData(sce_obj) <- NULL
SummarizedExperiment::rowData(sce_obj) <- sce_rd_new
sce_obj$sizeFactor <- NULL
sce_obj$ExpLibSize <- NULL

## Save =======================================================================
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
col_ord <- order(sce_obj$genotype, sce_obj$sample)

# Row order
get_sorted_indices <- function(status_name, row_data) {
  idx <- which(row_data$status == status_name)
  if (length(idx) == 0) return(integer(0))
  idx[order(row_data$sim_logFC[idx], decreasing = TRUE, na.last = TRUE)]
}

ko_idx   <- get_sorted_indices("Knockout", rd)
down_idx <- get_sorted_indices("Downstream", rd)
back_idx <- get_sorted_indices("Background", rd)

top_downstream <- down_idx[order(rd$total_reads[down_idx], decreasing = TRUE)][1:20]
top_background <- back_idx[order(rd$total_reads[back_idx], decreasing = TRUE)][1:19]

top_downstream <- top_downstream[order(match(top_downstream, down_idx))]
top_background <- top_background[order(match(top_background, back_idx))]

row_ord  <- c(ko_idx, top_downstream, top_background)

# Matrix
mat <- as.matrix(SummarizedExperiment::assay(sce_obj, "logcounts"))
mat <- mat[row_ord, col_ord, drop = FALSE]

col_mat <- circlize::colorRamp2(
  c(0, max(mat, na.rm = TRUE) * 0.4, max(mat, na.rm = TRUE)), 
  c("#000000", "#A52A2A", "#FFD700") 
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
  `Gene status` = status,
  `Sim. log2FC` = sim_logfc,
  `Obs. log2FC` = obs_logfc,
  col = list(
    `Gene status` = c(
      "Knockout"   = wesanderson::wes_palette("Chevalier1")[2],
      "Downstream" = wesanderson::wes_palette("Chevalier1")[4],  
      "Background" = wesanderson::wes_palette("Chevalier1")[3]
    ),
    `Sim. log2FC` = col_logfc,
    `Obs. log2FC` = col_logfc
  ),
  which = "row",
  annotation_name_side = "bottom",
  show_annotation_name = TRUE
)

# Right Row Annotations
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
  "Genes with a center log2FC of ±%s (facLoc %s)", 
  true_logfc, facLoc
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
  top_annotation = NULL,
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
down_idx <- get_sorted_indices("Downstream", rd)
back_idx <- get_sorted_indices("Background", rd)

top_downstream <- down_idx[order(rd$total_reads[down_idx], decreasing = TRUE)[1:20]]
top_background <- back_idx[order(rd$total_reads[back_idx], decreasing = TRUE)[1:19]]

selected_genes <- c(ko_idx, top_downstream, top_background)
selected_genes <- rownames(sce_obj)[selected_genes]

# Norm. counts
mat <- as.matrix(SingleCellExperiment::logcounts(sce_obj))
mat <- mat[selected_genes, , drop = FALSE]

# Plot data
p_data <- as.data.frame(t(mat))
p_data$group <- sce_obj$group
p_data$genotype <- sce_obj$genotype

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
down_names <- rownames(sce_obj)[top_downstream]

p_data$status                              <- "Background"
p_data$status[p_data$gene %in% ko_names]   <- "Knockout"
p_data$status[p_data$gene %in% down_names] <- "Downstream"

p_data$gene   <- factor(p_data$gene, levels = selected_genes)
p_data$status <- factor(p_data$status, levels = c("Knockout", "Downstream", "Background"))
p_data$group  <- factor(p_data$group, levels = c("targeting", "non-targeting"))

# Background cols
legend_cols <- c(
  "Knockout"   = wesanderson::wes_palette("Chevalier1")[2],  
  "Downstream" = wesanderson::wes_palette("Chevalier1")[4],  
  "Background" = wesanderson::wes_palette("Chevalier1")[3]
)

bg_cols_with_alpha <- paste0(legend_cols, "40")
names(bg_cols_with_alpha) <- names(legend_cols)

# Plot
p <- ggplot2::ggplot(
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
  alpha = 0.4,
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
    x = genotype,
    y = expression,
    fill = group
  ),
  alpha = 0.6,
  scale = "width", 
  color = "transparent",
  show.legend = TRUE
) +
ggplot2::scale_fill_manual(
  name = "Group",
  values = c(
    "targeting"     = wesanderson::wes_palette("Moonrise2")[3],
    "non-targeting" = wesanderson::wes_palette("Moonrise2")[1]
  )
) +
ggplot2::geom_boxplot(
  data = p_data,
  mapping = ggplot2::aes(
    x = genotype,
    y = expression
  ),
  width = 0.1, 
  linewidth = 0.1,
  fill = "white", 
  color = "grey30", 
  outlier.shape = NA
) +
ggplot2::facet_wrap(
  ~ gene, 
  scales = "free_y", 
  ncol = 5
) +
ggplot2::scale_y_continuous(
  n.breaks = 4,
  labels = scales::label_number(accuracy = 1)
) +
ggplot2::labs(
  x = "Genotype",
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
    fill = group
  )
) +
ggplot2::geom_bar(
  position = "dodge", 
  width = 0.6
) +
ggplot2::scale_fill_manual(
  name = NULL,
  values = c(
    "targeting"     = wesanderson::wes_palette("Moonrise2")[3],
    "non-targeting" = wesanderson::wes_palette("Moonrise2")[1]
  )
) +
ggplot2::geom_text(
  mapping = ggplot2::aes(
    label = ggplot2::after_stat(count)
  ),
  position = ggplot2::position_dodge(
    width = 0.6
  ),
  stat = "count",
  vjust = -0.5,
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
    fill = status
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
    "Downstream" = wesanderson::wes_palette("Chevalier1")[4],  
    "Background" = wesanderson::wes_palette("Chevalier1")[3]
  )
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

true_val   <- unique(rd$true_logFC)
vline_data <- data.frame(
  status = "Downstream",
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

# Patchwork
p_umap <- patchwork::wrap_plots(
  p_sample,
  p_genotype
)+ patchwork::plot_layout(
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
