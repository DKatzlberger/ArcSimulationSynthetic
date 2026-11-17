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
  real_data <- gsub(pattern = "\\\\ ", replacement = " ", trail_args[3])
  celltype  <- as.character(trail_args[4])
  knockout  <- c("Brd9", "Kmt2d")

  # Output
  input_name  <- tools::file_path_sans_ext(basename(input_path))
  output_name <- paste0(script_id, "-", input_name)
  output_name <- paste0(output_name)

} else {
  # Commands
  script_path <- rstudioapi::getActiveDocumentContext()$path
  script_name <- tools::file_path_sans_ext(basename(script_path))
  script_id   <- sub(pattern = "-.*", replacement = "", script_name)

  # Paths
  project_root <- "/home/scc_e_928652/ArcSimulationSynthetic"
  input_path   <- "results/manual/ArcSimSyn_offtarget_mixed/03-scdesign3_offtarget_mixed/03-02-01-celltype_Gran. P_knockout_NULL_facLoc_2_fracOff_50_koEff_100.h5ad"
  output_path  <- file.path("results/manual/ArcSimSyn_offtarget_mixed", script_name)

  stopifnot(dirname(dirname(input_path)) == dirname(output_path))
  if (!dir.exists(output_path)) dir.create(path = output_path, recursive = TRUE)

  # Params
  real_data <- "data/GSE213511/rna_integrated.rds"
  celltype  <- "Gran. P"
  knockout  <- c("Brd9", "Kmt2d")

  # Output
  input_name  <- tools::file_path_sans_ext(basename(input_path))
  output_name <- paste0(script_id, "-", input_name)
  output_name <- paste0(output_name)
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


## Compare ====================================================================
### Sim. data -----------------------------------------------------------------
sim_obj <- sce_obj

### Real data -----------------------------------------------------------------
ext <- tools::file_ext(real_data)
txt <- sprintf("%s as 'Seurat'", real_data)

real_obj <- switch(EXPR = ext,
  rds  = .read_rds(real_data, info = txt),
  h5ad = .read_h5ad(real_data, as = "Seurat", info = txt),
  showf("Unsupported file extension: %s", ext, type = "stop")
)

### Celltype ------------------------------------------------------------------
.subset <- fun_with_ram(base::subset)

txt      <- sprintf("subsetting to celltype '%s'", celltype)
keep     <- colnames(real_obj)[real_obj$celltype %in% celltype]
real_obj <- .subset(real_obj, cells = keep, info = txt)

### Genotype, guide & on_guide ------------------------------------------------
.assign_guide_by_max_count <- fun_with_ram(assign_guide_by_max_count)

real_obj <- .assign_guide_by_max_count(real_obj)

# Check
has_KO  <- !is.null(knockout) && any(knockout %in% real_obj$genotype)
has_NTC <- "NTC" %in% real_obj$genotype

if (!has_KO && !has_NTC) showf("Neither '%s' nor 'NTC' assigned to any cells.", knockout, type = "stop")
if (!is.null(knockout) && !has_KO) showf("'%s' not assigned to any cells.", knockout, type = "stop")
if (!has_NTC) showf("'%s' not assigned to any cells.", "NTC", type = "stop")

### Cells: subset by knockout, NTC & feature counts ---------------------------
if (is.null(knockout)) {
  txt  <- "subsetting to 'NTC'"
  keep <- colnames(real_obj)[real_obj$genotype %in% "NTC"]
} else {
  txt  <- sprintf("subsetting to genotype '%s' and 'NTC'", paste(knockout, collapse = ", "))
  keep <- colnames(real_obj)[real_obj$genotype %in% c(knockout, "NTC")]
}

real_obj <- .subset(real_obj, cells = keep, info = txt)
real_obj <- subset(real_obj, subset = nFeature_RNA >= 1000 & nFeature_RNA <= 5000)

### Genes: subset by detection ------------------------------------------------
matr <- SeuratObject::LayerData(real_obj, assay = "RNA", layer = "counts")
keep <- Matrix::rowSums(matr > 0) >= 10

real_obj <- subset(real_obj, features = rownames(matr)[keep])

### Convert to SingleCellExperiment -------------------------------------------
matr <- SeuratObject::LayerData(real_obj, assay = "RNA", layer = "counts")
meta <- real_obj@meta.data

real_obj <- SingleCellExperiment::SingleCellExperiment(
  assays  = list(counts = matr),
  colData = meta
)

# Check input
txt <- "checking real obj."
.check_obj(real_obj, info = txt)


## Visualize ==================================================================
### Panel A -------------------------------------------------------------------
#### Plot data ----------------------------------------------------------------
SummarizedExperiment::assay(sim_obj, "logCPM")  <- log2(scuttle::calculateCPM(sim_obj) + 1)
SummarizedExperiment::assay(real_obj, "logCPM") <- log2(scuttle::calculateCPM(real_obj) + 1)

real_counts <- SummarizedExperiment::assay(real_obj, "counts")
sim_counts  <- SummarizedExperiment::assay(sim_obj, "counts")

real_logCPM <- SummarizedExperiment::assay(real_obj, "logCPM")
sim_logCPM  <- SummarizedExperiment::assay(sim_obj, "logCPM")

# Groups
real_g <- split(seq_along(real_obj$genotype), interaction(real_obj$genotype, real_obj$guide, drop = TRUE, sep = "_"))
sim_g  <- split(seq_along(sim_obj$genotype), interaction(sim_obj$genotype, sim_obj$guide, drop = TRUE, sep = "_"))

# Gene summary
real_gene_summary <- do.call(rbind, lapply(names(real_g), function(grp) {
  cells <- real_g[[grp]]

  data.frame(
    gene = rownames(real_logCPM),
    mean_counts = rowMeans(real_counts[, cells, drop = FALSE]),
    mean_logCPM = rowMeans(real_logCPM[, cells, drop = FALSE]),
    var_logCPM = apply(real_logCPM[, cells, drop = FALSE], 1, var),
    genotype = real_obj$genotype[cells[1]],
    guide = real_obj$guide[cells[1]],
    dataset = "Real"
  )
}))

sim_gene_summary <- do.call(rbind, lapply(names(sim_g), function(grp) {
  cells <- sim_g[[grp]]

  data.frame(
    gene = rownames(sim_logCPM),
    mean_counts = rowMeans(sim_counts[, cells, drop = FALSE]),
    mean_logCPM = rowMeans(sim_logCPM[, cells, drop = FALSE]),
    var_logCPM = apply(sim_logCPM[, cells, drop = FALSE], 1, var),
    genotype = sim_obj$genotype[cells[1]],
    guide = sim_obj$guide[cells[1]],
    dataset = "Simulated"
  )
}))

real_gene_zero <- do.call(rbind, lapply(names(real_g), function(grp) {
  cells <- real_g[[grp]]

  data.frame(
    gene = rownames(real_counts),
    pct_zero = Matrix::rowMeans(
      real_counts[, cells, drop = FALSE] == 0
    ) * 100,
    genotype = real_obj$genotype[cells[1]],
    guide = real_obj$guide[cells[1]],
    dataset = "Real"
  )
}))

# Gene zero
real_gene_zero <- do.call(rbind, lapply(names(real_g), function(grp) {
  cells <- real_g[[grp]]

  data.frame(
    gene = rownames(real_counts),
    pct_zero = Matrix::rowMeans(real_counts[, cells, drop = FALSE] == 0) * 100,
    genotype = real_obj$genotype[cells[1]],
    guide = real_obj$guide[cells[1]],
    dataset = "Real"
  )
}))

sim_gene_zero <- do.call(rbind, lapply(names(sim_g), function(grp) {
  cells <- sim_g[[grp]]

  data.frame(
    gene = rownames(sim_counts),
    pct_zero = Matrix::rowMeans(sim_counts[, cells, drop = FALSE] == 0) * 100,
    genotype = sim_obj$genotype[cells[1]],
    guide = sim_obj$guide[cells[1]],
    dataset = "Simulated"
  )
}))

# Library size
real_lib_sizes <- data.frame(
  total_counts = Matrix::colSums(real_counts),
  genotype = real_obj$genotype,
  guide = real_obj$guide,
  dataset = "Real"
)

sim_lib_sizes <- data.frame(
  total_counts = Matrix::colSums(sim_counts),
  genotype = sim_obj$genotype,
  guide = sim_obj$guide,
  dataset = "Simulated"
)

# Cell zero
real_cell_zero <- data.frame(
  pct_zero = Matrix::colMeans(real_counts == 0) * 100,
  genotype = real_obj$genotype,
  guide = real_obj$guide,
  dataset = "Real"
)

sim_cell_zero <- data.frame(
  pct_zero = Matrix::colMeans(sim_counts == 0) * 100,
  genotype = sim_obj$genotype,
  guide = sim_obj$guide,
  dataset = "Simulated"
)

# Combine 
gene_summary <- rbind(real_gene_summary, sim_gene_summary)
lib_sizes <- rbind(real_lib_sizes, sim_lib_sizes)
gene_zero <- rbind(real_gene_zero, sim_gene_zero)
cell_zero <- rbind(real_cell_zero, sim_cell_zero)

#### Panel A01: Distribution of mean ------------------------------------------
p_panel_A01 <- ggplot2::ggplot(
  data = gene_summary,
  mapping = ggplot2::aes(
    x = genotype,
    y = mean_logCPM,
    color = guide
  )
) +
ggplot2::geom_boxplot(
  linewidth = 0.3,
  coef = Inf,
  show.legend = FALSE
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(dataset),
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Genotype",
  y = expression("Mean log"[2] * "(CPM + 1)"),
  color = "Guide"
) +
theme_DK(
  rotate = 45,
  show_grid = FALSE
)

#### Panel A02: Distribution of variance --------------------------------------
p_panel_A02 <- ggplot2::ggplot(
  data = gene_summary,
  mapping = ggplot2::aes(
    x = genotype,
    y = var_logCPM,
    color = guide
  )
) +
ggplot2::geom_boxplot(
  linewidth = 0.3,
  coef = Inf,
  show.legend = FALSE
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(dataset),
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Genotype",
  y = expression("Variance log"[2] * "(CPM + 1)"),
  color = "Guide"
) +
theme_DK(
  rotate = 45,
  show_grid = FALSE
)

#### Panel A03: Mean variance relationship ------------------------------------
p_panel_A03 <- ggplot2::ggplot(
  data = gene_summary,
  mapping = ggplot2::aes(
    x = mean_logCPM,
    y = var_logCPM,
    color = guide
  )
) +
# ggplot2::geom_point(
#   shape = 16,
#   size = 0.2,
#   stroke = 0,
#   show.legend = FALSE
# ) +
ggplot2::geom_smooth(
  formula = "y ~ x",
  method = "loess",
  se = FALSE,
  linewidth = 0.2,
  show.legend = FALSE
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(dataset)
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

#### Panel A04: Library size --------------------------------------------------
p_panel_A04 <- ggplot2::ggplot(
  data = lib_sizes,
  mapping = ggplot2::aes(
    x = genotype,
    y = total_counts,
    color = guide
  )
) +
ggplot2::geom_boxplot(
  linewidth = 0.3,
  coef = Inf,
  show.legend = FALSE
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(dataset),
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Genotype",
  y = "Total counts per cell"
) +
theme_DK(
  rotate = 45,
  show_grid = FALSE
)

#### Panel A05: Percentage zero per cell --------------------------------------
p_panel_A05 <- ggplot2::ggplot(
  data = cell_zero,
  mapping = ggplot2::aes(
    x = genotype,
    y = pct_zero,
    color = guide
  )
) +
ggplot2::geom_boxplot(
  linewidth = 0.3,
  coef = Inf,
  show.legend = FALSE
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(dataset),
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Genotype",
  y = "Percentage zeros per cell"
) +
ggplot2::scale_y_continuous(
  limits = c(0, 100),
  breaks = seq(0, 100, 25)
) +
theme_DK(
  rotate = 45,
  show_grid = FALSE
)

#### Panel A06: Percentage zero per gene --------------------------------------
p_panel_A06 <- ggplot2::ggplot(
  data = gene_zero,
  mapping = ggplot2::aes(
    x = genotype,
    y = pct_zero,
    color = guide
  )
) +
ggplot2::geom_boxplot(
  linewidth = 0.3,
  coef = Inf,
  show.legend = FALSE
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(dataset),
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Genotype",
  y = "Percentage zero per gene",
  color = "Guide"
) +
ggplot2::scale_y_continuous(
  limits = c(0, 100),
  breaks = seq(0, 100, 25)
) +
theme_DK(
  rotate = 45,
  show_grid = FALSE
)

#### Panel A07: Percentage of zeros versus mean expression --------------------
p_panel_A07 <- ggplot2::ggplot(
  data = cbind(
    gene_zero,
    mean_counts = gene_summary$mean_counts
  ),
  mapping = ggplot2::aes(
    x = mean_counts,
    y = pct_zero,
    color = guide
  )
) +
# ggplot2::geom_point(
#   shape = 16,
#   size = 0.2,
#   stroke = 0,
#   show.legend = FALSE
# ) +
ggplot2::geom_smooth(
  formula = "y ~ x",
  method = "loess",
  se = FALSE,
  linewidth = 0.2,
  show.legend = FALSE
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(dataset)
) +
ggplot2::scale_y_continuous(
  limits = c(0, 100),
  breaks = seq(0, 100, 25)
) +
ggplot2::labs(
  x = "Mean count",
  y = "Percentage zero"
) +
theme_DK(
  rotate = 45,
  show_grid = FALSE
)

#### Panel A08: Distribution of percentage of zeros by gene -------------------
p_panel_A08 <- ggplot2::ggplot(
  data = gene_zero,
  mapping = ggplot2::aes(
    x = pct_zero,
    fill = guide
  )
) +
ggplot2::geom_density(
  alpha = 0.3,
  linewidth = 0.1,
  show.legend = FALSE
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(dataset)
) +
ggplot2::scale_y_continuous(
  limits = c(0, 100),
  breaks = seq(0, 100, 25)
) +
ggplot2::labs(
  x = "Percentage zero per gene",
  y = "Density",
  fill = "Guide"
) +
theme_DK(
  rotate = 45,
  show_grid = FALSE
)

#### Patchwork & save ---------------------------------------------------------
.ggsave <- fun_with_ram(ggplot2::ggsave)

# Plot
p_panel_A <-
  p_panel_A01 +
  p_panel_A02 +
  p_panel_A03 +
  p_panel_A04 +
  p_panel_A05 +
  p_panel_A06 +
  p_panel_A07 +
  patchwork::guide_area() +
  patchwork::plot_layout(
    ncol = 4,
    nrow = 2,
    guides = "collect"
  ) 

# Save
txt <- sprintf("saved as %s", paste0(output_name, "_A.svg"))
.ggsave(
  p_panel_A,
  file = file.path(
    output_path,
    paste0(
      output_name,
      "_A.svg"
    )
  ),
  width = 18,
  height = 8,
  units = "cm",
  info = txt
)


## Script info ================================================================
showf("\n◷ %s ====================================================\n", date(), type = "cat")
