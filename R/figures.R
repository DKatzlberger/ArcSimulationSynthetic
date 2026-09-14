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

} else {
  project_root <- "/home/scc_e_928652/ArcSimulationSynthetic"
  input_path   <- "results/nextflow/ArcSimSyn_offtarget_mixed"
  output_path  <- "results/nextflow/ArcSimSyn_offtarget_mixed/figures"

  alpha <- 0.05

  file_name <- paste0("_alpha_", alpha)
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
    "patchwork"
  )
)

## Dataset ====================================================================
readDGE <- function(
  path,
  pattern
){

  # Expected files
  dat_comb <- 700

  # Formula: ~genotype
  n_geno_limma    <- dat_comb * 2   # 1400
  n_geno_edger    <- dat_comb * 2   # 1400
  n_geno_deseq2   <- dat_comb * 2   # 1400
  n_geno_dream    <- dat_comb * 5   # 3500
  n_geno_nebula   <- dat_comb * 4   # 2800
  n_geno_sceptre  <- dat_comb * 2   # 1400
  n_geno_wilcox   <- dat_comb * 1   # 700
  n_geno_scmageck <- dat_comb * 1   # 700

  exp_genotype <- 
    n_geno_limma + n_geno_edger + 
    n_geno_deseq2 + n_geno_dream + 
    n_geno_nebula + n_geno_sceptre + 
    n_geno_wilcox + n_geno_scmageck

  # Formula: ~guide
  n_guide_limma   <- dat_comb * 2   # 1400
  n_guide_edger   <- dat_comb * 2   # 1400
  n_guide_deseq2  <- dat_comb * 2   # 1400
  n_guide_dream   <- dat_comb * 3   # 2100
  n_guide_nebula  <- dat_comb * 3   # 2100

  exp_guide <- 
    n_guide_limma + n_guide_edger + 
    n_guide_deseq2 + n_guide_dream + 
    n_guide_nebula

  # Total expected
  exp_total <- exp_genotype + exp_guide   # 21700

  # Found files
  exp_found <- list.files(
    path = path, 
    pattern = "_dge\\.rds$", 
    full.names = TRUE, 
    recursive = TRUE
  )

  # Check for duplication
  any_duplicated_paths <- any(duplicated(exp_found))
  if (any_duplicated_paths) {
    n_dup_paths <- sum(duplicated(exp_found))
    showf("Note: There are %s duplicated DGE files.", n_dup_paths)
  }

  showf("Expected %s DGE files found %s (pattern: total).", exp_total, length(exp_found), type = "message")
  showf(" ~genotype: %s", exp_genotype)
  showf(" ~guide:    %s", exp_guide)

  # Find DGE statistics as RDS
  file_paths <- list.files(
    path = path, 
    pattern = pattern, 
    full.names = TRUE, 
    recursive = TRUE
  )

  # Check for duplication
  any_duplicated_paths <- any(duplicated(file_paths))
  if (any_duplicated_paths) {
    n_dup_paths <- sum(duplicated(file_paths))
    showf("Note: There are %s duplicated DGE files.", n_dup_paths)
  }

  # Pre-allocate a list to store the data frames
  loaded_data <- vector("list", length(file_paths))
  bad_files <- c()

  # Load the files using a loop to track progress
  for (i in seq_along(file_paths)) {
    current_path <- file_paths[i]
    
    loaded_data[[i]] <- tryCatch({
      as.data.frame(readRDS(current_path))
    },
    error = function(e) {
      bad_files <<- c(bad_files, current_path)
      NULL
    })
    
    # Print progress every 500 files
    if (i %% 1000 == 0) {
      txt <- sprintf("%s/%s", i, length(file_paths))
      showf("Loaded %s", txt, type = "message")
    }
  }

  if (length(bad_files) > 0) {
    showf(
      "NOTE: Not all files could be loaded. %s failed files.",
      length(bad_files),
      type = "message"
    )
  }

  # Column standardization
  showf("Standardizing columns.", type = "message")
  
  valid_data <- loaded_data[!vapply(loaded_data, is.null, logical(1))]
  if (length(valid_data) > 0) {

    all_col_lists  <- lapply(valid_data, colnames)
    union_cols     <- unique(unlist(all_col_lists))
    intersect_cols <- Reduce(intersect, all_col_lists)
    
    non_overlapping_cols <- setdiff(union_cols, intersect_cols)
    
    if (length(non_overlapping_cols) > 0) {
      showf(
        "Removing '%s' non-overlapping column(s).", 
        paste(non_overlapping_cols, collapse = ", "), 
        type = "message"
      )
      
      loaded_data <- lapply(loaded_data, function(df) {
        if (is.null(df)) return(NULL)
        return(df[, intersect_cols, drop = FALSE])
      })
    } else {
      showf("All files have matching columns.", type = "message")
    }
  }

  # Combine into a dataframe
  do.call(rbind, loaded_data)
}

dge <- readDGE(
  path = input_path,
  pattern = "~genotype.*_dge\\.rds$"
)

dge$sim_mode <- tolower(dge$sim_mode)
dge$facLoc   <- factor(dge$facLoc, levels = c(0.1, 0.5, 1, 1.5, 2))
dge$nOn      <- factor(dge$nOn)
dge$koEff    <- factor(dge$koEff)
dge$nGuide   <- factor(dge$nGuide)

dge$formula  <- factor(
  dge$formula, 
  levels = c(
    "~genotype", 
    "~genotype + sample", 
    "~genotype + (1|sample)",
    "~genotype + sample + (1|guide)", 
    "~genotype + (1|sample) + (1|guide)"
  )
)

## Stats. =====================================================================
group_vars <- c("sim_mode", "workflow", "method_type", "method", "formula", "facLoc", "nOn", "koEff", "nGuide")

# Confusion matrix
stats <- aggregate(
  cbind(
    # TN: FDR bigger than alpha
    # FN: FDR bigger than alpha
    # TP: FDR smaller than alpha
    # FP: FDR smaller than alpha
    TN     = dge$status %in% c("Background", "Off-target") & dge$FDR >= alpha,
    FN     = dge$status %in% c("Knockout", "On-target")    & dge$FDR >= alpha,
    TP     = dge$status %in% c("Knockout", "On-target")    & dge$FDR < alpha,
    FP     = dge$status %in% c("Background", "Off-target") & dge$FDR < alpha,

    TN_off = dge$status == "Off-target"                    & dge$FDR >= alpha,
    TN_bg  = dge$status == "Background"                    & dge$FDR >= alpha,
    TP_ko  = dge$status == "Knockout"                      & dge$FDR < alpha,
    TP_on  = dge$status == "On-target"                     & dge$FDR < alpha,
    FP_off = dge$status == "Off-target"                    & dge$FDR < alpha,
    FP_bg  = dge$status == "Background"                    & dge$FDR < alpha
  ),
  by  = dge[group_vars],
  FUN = function(x) sum(x, na.rm = TRUE)
)

# Pct. DEGs per status
stats$pct_ko  <- stats$TP_ko  / (stats$TP + stats$FP)
stats$pct_on  <- stats$TP_on  / (stats$TP + stats$FP)
stats$pct_off <- stats$FP_off / (stats$TP + stats$FP)
stats$pct_bg  <- stats$FP_bg  / (stats$TP + stats$FP)

status_long <- reshape(
  stats,
  varying = c("pct_ko", "pct_on", "pct_off", "pct_bg"),
  v.names = "pct",
  timevar = "status",
  times = c("Knockout", "On-target", "Off-target", "Background"),
  direction = "long"
)

# Ratio between FPs (Background vs. Off-target)
stats$fpr_off_target <- ifelse((stats$TN_off + stats$FP_off) == 0, 0, stats$FP_off / (stats$TN_off + stats$FP_off))
stats$fpr_background <- ifelse((stats$TN_bg + stats$FP_bg) == 0, 0, stats$FP_bg / (stats$TN_bg + stats$FP_bg))

epsilon <- 0.001
stats$fpr_ratio <- log2((stats$fpr_off_target + epsilon) / (stats$fpr_background + epsilon))

# Classic metric
stats$sensitivity <- ifelse((stats$TP + stats$FN) == 0, 0, stats$TP / (stats$TP + stats$FN))
stats$specificity <- ifelse((stats$FP + stats$TN) == 0, 0, stats$TN / (stats$FP + stats$TN))
stats$precision   <- ifelse((stats$TP + stats$FP) == 0, 0, stats$TP / (stats$TP + stats$FP))
stats$fdr_total   <- ifelse((stats$TP + stats$FP) == 0, 0, stats$FP / (stats$TP + stats$FP))

# Higher metric
f1_numerator   <- stats$precision * stats$sensitivity
f1_denominator <- stats$precision + stats$sensitivity
stats$f1       <- ifelse(f1_denominator == 0, 0, 2 * f1_numerator / f1_denominator)

beta              <- 0.5
fbeta_numerator   <- (1 + beta^2) * (stats$precision * stats$sensitivity)
fbeta_denominator <- ((beta^2 * stats$precision) + stats$sensitivity)
stats$fbeta       <- ifelse(fbeta_denominator == 0, 0, fbeta_numerator / fbeta_denominator)

mcc_numerator   <- (as.numeric(stats$TP) * as.numeric(stats$TN)) - (as.numeric(stats$FP) * as.numeric(stats$FN))
mcc_denominator <- as.numeric(stats$TP + stats$FP) * as.numeric(stats$TP + stats$FN) * as.numeric(stats$TN + stats$FP) * as.numeric(stats$TN + stats$FN)
stats$mcc       <- ifelse(mcc_denominator == 0, 0, mcc_numerator / sqrt(mcc_denominator))

higher_metric <- aggregate(
  cbind(f1, fbeta, mcc) ~ method + formula + workflow + method_type + sim_mode,
  data = stats,
  FUN = mean,
  na.rm = TRUE
)

higher_metric_long <- reshape(
  higher_metric,
  varying = c("f1", "fbeta", "mcc"),
  v.names = "value",
  timevar = "metric",
  times = c("F1", "F-beta", "MCC"),
  direction = "long"
)

higher_metric_long$is_best <- as.logical(
  ave(
    higher_metric_long$value, 
    higher_metric_long$metric, 
    FUN = function(x) x == max(x, na.rm = TRUE)
  )
)


## Panel A ====================================================================
## Panel A0: Annotations ======================================================

# Sim. mode
p_sim_mode <- ggplot2::ggplot(
  data = stats,
  mapping = ggplot2::aes(
    x = method,
    y = "Sim. mode",
    fill = sim_mode
  )
) +
ggplot2::geom_tile(
  color = NA,
  linewidth = 0.1,
  show.legend = TRUE
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::scale_fill_manual(
  name = "Sim. mode", 
  values = c(
    "mixed guides"  = RColorBrewer::brewer.pal(8, "Dark2")[3], 
    "unique guides" = RColorBrewer::brewer.pal(8, "Dark2")[4]
  ),
  limits = c(
    "mixed guides", 
    "unique guides"
  ),
  drop = FALSE
) +
ggplot2::labs(
  x = NULL,
  y = NULL
) +
theme_DK(
  show_facet_x = TRUE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_borders = FALSE,
  show_grid = FALSE,
  show_axis = FALSE
) +
ggplot2::theme(
  strip.text.x = ggplot2::element_text(
    angle = 90,
    hjust = 0,
    vjust = 0.5
  ),
  axis.title.y = ggplot2::element_text(
    angle = 0,
    vjust = 0.5,
    hjust = 0.5
  ),
  plot.margin = ggplot2::margin(0, 0, 0, 0, "pt")
)

# Workflows
p_workflow <- ggplot2::ggplot(
  data = stats,
  mapping = ggplot2::aes(
    x = method,
    y = "Workflow",
    fill = workflow
  )
) +
ggplot2::geom_tile(
  color = "black",
  linewidth = 0.1
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::scale_fill_brewer(
  name = "Workflow",
  palette = "Set2"
) +
ggplot2::labs(
  x = "Method",
  y = NULL
) +
theme_DK(
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_borders = FALSE,
  show_grid = FALSE,
  show_axis = FALSE
) +
ggplot2::theme(
  axis.title.y = ggplot2::element_text(
    angle = 0,
    vjust = 0.5,
    hjust = 0.5
  ),
  plot.margin = ggplot2::margin(0, 0, 0, 0, "pt")
)

# Method type
p_method_type <- ggplot2::ggplot(
  data = stats,
  mapping = ggplot2::aes(
    x = method,
    y = "Method type",
    fill = method_type
  )
) +
ggplot2::geom_tile(
  color = "black",
  linewidth = 0.1
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::scale_fill_brewer(
  name = "Method type",
  palette = "Set1"
) +
ggplot2::labs(
  x = NULL,
  y = NULL
) +
theme_DK(
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_borders = FALSE,
  show_grid = FALSE,
  show_axis = FALSE
) +
ggplot2::theme(
  axis.title.y = ggplot2::element_text(
    angle = 0,
    vjust = 0.5,
    hjust = 0.5
  ),
  plot.margin = ggplot2::margin(0, 0, 0, 0, "pt")
)

# Patchwork
p_fig_1_panel_A0 <- wrap_plots(
  p_sim_mode,
  p_workflow,
  p_method_type,
  ncol = 1
)

## Panel A1: Gene status ======================================================
p_fig_1_panel_A1 <- ggplot2::ggplot(
  data = status_long, 
  mapping = ggplot2::aes(
    x = method,
  )
) +
ggplot2::geom_bar(
  mapping = ggplot2::aes(
    fill = status,
    weight = pct
  ),
  stat = "count",
  position = "fill",
  width = 0.8
) +
ggplot2::scale_y_continuous(
  breaks = seq(0, 1, by = 0.25)
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula), 
  scales = "free_x",
  space = "free"
) +
ggplot2::scale_fill_manual(
  name = "Gene status",
  values = c(
    "Knockout"   = wesanderson::wes_palette("Chevalier1")[2],  
    "On-target"  = wesanderson::wes_palette("Chevalier1")[4],  
    "Off-target" = wesanderson::wes_palette("Chevalier1")[1],
    "Background" = wesanderson::wes_palette("Chevalier1")[3]
  )
) +
ggplot2::labs(
  x = "Method", 
  y = "Pct. of DEGs"
) +
theme_DK(
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel A2: Power ============================================================
mean_metric <- aggregate(
  sensitivity ~ method + formula, 
  data = stats, 
  FUN = mean
)
colnames(mean_metric)[colnames(mean_metric) == "sensitivity"] <- "value"

mean_metric$x_pos <- as.numeric(as.factor(mean_metric$method))

winner_val  <- max(mean_metric$value, na.rm = TRUE)
global_best <- mean_metric[mean_metric$value == winner_val, ]

global_best$xmin <- global_best$x_pos - 0.4
global_best$xmax <- global_best$x_pos + 0.4
global_best$ymin <- -Inf
global_best$ymax <- Inf

# Plot
p_fig_1_panel_A2 <- ggplot2::ggplot(
  data = stats, 
  mapping = ggplot2::aes(
    x = method, 
    y = sensitivity
  )
) +
ggplot2::geom_rect(
  data = global_best,
  mapping = ggplot2::aes(
    xmin = xmin, xmax = xmax, 
    ymin = ymin, ymax = ymax
  ),
  fill = "pink",
  inherit.aes = FALSE
) +
ggplot2::stat_summary(
  fun = "mean",
  geom = "col",
  position = ggplot2::position_dodge(
    width = 0.6
  ),
  width = 0.8,
  color = "black",
  fill = "grey80",
  linewidth = 0.1
) +
ggplot2::stat_summary(
  fun.data = "mean_sdl",
  fun.args = list(mult = 1),
  geom = "errorbar",
  position = ggplot2::position_dodge(
    width = 0.6
  ),
  color = "black",
  width = 0.2,
  linewidth = 0.2
) +
ggplot2::scale_y_continuous(
  limits = function(x) c(min(0, x), max(1, x)),
  breaks = seq(0, 1, by = 0.25)
  # labels = scales::percent_format(
  #   accuracy = 1
  # )
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula), 
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Method", 
  y = "Power"
) +
theme_DK(
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel A3: Realized FDR =====================================================
mean_metric <- aggregate(
  fdr_total ~ method + formula, 
  data = stats, 
  FUN = mean
)
colnames(mean_metric)[colnames(mean_metric) == "fdr_total"] <- "value"

mean_metric$x_pos <- as.numeric(as.factor(mean_metric$method))

winner_val  <- min(mean_metric$value, na.rm = TRUE)
global_best <- mean_metric[mean_metric$value == winner_val, ]

global_best$xmin <- global_best$x_pos - 0.4
global_best$xmax <- global_best$x_pos + 0.4
global_best$ymin <- -Inf
global_best$ymax <- Inf

# Plot
p_fig_1_panel_A3 <- ggplot2::ggplot(
  data = stats, 
  mapping = ggplot2::aes(
    x = method, 
    y = fdr_total
  )
) +
ggplot2::geom_rect(
  data = global_best,
  mapping = ggplot2::aes(
    xmin = xmin, xmax = xmax, 
    ymin = ymin, ymax = ymax
  ),
  fill = "pink",
  inherit.aes = FALSE
) +
ggplot2::stat_summary(
  fun = "mean",
  geom = "col",
  position = ggplot2::position_dodge(
    width = 0.6
  ),
  width = 0.8,
  color = "black",
  fill = "grey80",
  linewidth = 0.1
) +
ggplot2::stat_summary(
  fun.data = "mean_sdl",
  fun.args = list(mult = 1),
  geom = "errorbar",
  position = ggplot2::position_dodge(
    width = 0.6
  ),
  color = "black",
  width = 0.2,
  linewidth = 0.2
) +
ggplot2::geom_hline(
  ggplot2::aes(
    yintercept = alpha, 
    linetype = paste0("FDR ", alpha) 
  ),
  color = "black",     
  linewidth = 0.3
) +
ggplot2::scale_linetype_manual(
  name = "Threshold",
  values = setNames("solid", paste0("FDR ", alpha))
) +
ggplot2::scale_y_continuous(
  limits = function(x) c(min(0, x), max(1, x)),
  breaks = seq(0, 1, by = 0.25)
  # labels = scales::percent_format(
  #   accuracy = 1
  # )
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula), 
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Method", 
  y = "Realized FDR"
) +
theme_DK(
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel A4: Higher metric ====================================================
p_fig_1_panel_A4 <- ggplot2::ggplot(
  data = higher_metric_long,
  mapping = ggplot2::aes(
    x = method,
    y = metric,
    fill = is_best
  )
) +
ggplot2::geom_tile(
  color = "black",  
  linewidth = 0.1
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::scale_fill_manual(
  name = NULL,
  values = c(
    "TRUE"  = "pink",
    "FALSE" = "white"
  ),
  breaks = c("TRUE"),
  labels = c(
    "TRUE"  = "Best method"
  )
) +
ggplot2::labs(
  x = NULL,
  y = NULL
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_label_x = TRUE,
  show_ticks_x = TRUE,
  show_borders = FALSE,
  show_grid = FALSE,
  show_axis = FALSE
) +
ggplot2::theme(
  axis.title.y = ggplot2::element_text(
    angle = 0,
    vjust = 0.5,
    hjust = 0.5
  ),
  plot.margin = ggplot2::margin(0, 0, 0, 0, "pt")
)

## Patchwork & save ===========================================================
p_fig_1_panel_A <-
  p_fig_1_panel_A0 + 
  p_fig_1_panel_A1 + 
  p_fig_1_panel_A2 +
  p_fig_1_panel_A3 +
  p_fig_1_panel_A4 +
  patchwork::plot_layout(
    ncol = 1,
    heights = c(0.05, 0.05, 0.05, 0.3, 0.3, 0.3, 0.15),
    guides = "collect"
  ) 

# Save
ggplot2::ggsave(
  p_fig_1_panel_A,
  file = file.path(
    output_path,
    paste0(
      "panel_A",
      file_name,
      ".svg"
    )
  ),
  width = 8,
  height = 9.5,
  units = "cm"
)

saveRDS(
  p_fig_1_panel_A,
  file = file.path(
    output_path,
    paste0(
      "panel_A",
      file_name,
      ".rds"
    )
  ),
)


## Panel B ====================================================================
## Panel B0: Annotations ======================================================

# Sim. mode
p_sim_mode <- ggplot2::ggplot(
  data = stats,
  mapping = ggplot2::aes(
    x = nGuide,
    y = "Sim. mode",
    fill = sim_mode
  )
) +
ggplot2::geom_tile(
  color = NA,
  linewidth = 0.1,
  show.legend = TRUE
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::scale_fill_manual(
  name = "Sim. mode", 
  values = c(
    "mixed guides"  = RColorBrewer::brewer.pal(8, "Dark2")[3], 
    "unique guides" = RColorBrewer::brewer.pal(8, "Dark2")[4]
  ),
  limits = c(
    "mixed guides", 
    "unique guides"
  ),
  drop = FALSE
) +
ggplot2::labs(
  x = NULL,
  y = NULL
) +
theme_DK(
  show_facet_x = TRUE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_borders = FALSE,
  show_grid = FALSE,
  show_axis = FALSE
) +
ggplot2::theme(
  strip.text.x = ggplot2::element_text(
    angle = 90,
    hjust = 0,
    vjust = 0.5
  ),
  axis.title.y = ggplot2::element_text(
    angle = 0,
    vjust = 0.5,
    hjust = 0.5
  ),
  plot.margin = ggplot2::margin(0, 0, 0, 0, "pt")
)

# Patchwork
p_fig_1_panel_B0 <- wrap_plots(
  p_sim_mode,
  ncol = 1
)

## Panel B1: Power ============================================================
p_fig_1_panel_B1 <- ggplot2::ggplot(
  data = stats, 
  mapping = ggplot2::aes(
    x = nGuide, 
    y = sensitivity,
    color = method,
    group = method
  )
) +
ggplot2::stat_summary(
  fun = mean, 
  geom = "point", 
  size = 0.1
) +
ggplot2::stat_summary(
  fun = mean, 
  geom = "line", 
  linewidth = 0.3
) +
ggplot2::scale_y_continuous(
  limits = c(0, 1)
  # labels = scales::percent_format(
  #   accuracy = 1
  # )
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula), 
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Cells per guide", 
  y = "Power",
  color = "Method"
) +
theme_DK(
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel B2: Realized FDR =====================================================
p_fig_1_panel_B2 <- ggplot2::ggplot(
  data = stats, 
  mapping = ggplot2::aes(
    x = nGuide, 
    y = fdr_total,
    color = method,
    group = method
  )
) +
ggplot2::stat_summary(
  fun = mean, 
  geom = "point", 
  size = 0.1
) +
ggplot2::stat_summary(
  fun = mean, 
  geom = "line", 
  linewidth = 0.3
) +
ggplot2::geom_hline(
  ggplot2::aes(
    yintercept = alpha, 
    linetype = paste0("FDR ", alpha) 
  ),
  color = "black",     
  linewidth = 0.3
) +
ggplot2::scale_linetype_manual(
  name = "Threshold",
  values = setNames("solid", paste0("FDR ", alpha))
) +
ggplot2::scale_y_continuous(
  limits = c(0, 1)
  # labels = scales::percent_format(
  #   accuracy = 1
  # )
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula), 
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Cells per targeting guide", 
  y = "Realized FDR",
  color = "Method"
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_grid = FALSE
)

## Patchwork & save ===========================================================
p_fig_1_panel_B <-
  p_fig_1_panel_B0 + 
  p_fig_1_panel_B1 + 
  p_fig_1_panel_B2 +
  patchwork::plot_layout(
    ncol = 1,
    heights = c(0.05, 0.3, 0.3),
    guides = "collect"
  ) 

# Save
ggplot2::ggsave(
  p_fig_1_panel_B,
  file = file.path(
    output_path,
    paste0(
      "panel_B",
      file_name,
      ".svg"
    )
  ),
  width = 9,
  height = 7,
  units = "cm"
)

saveRDS(
  p_fig_1_panel_B,
  file = file.path(
    output_path,
    paste0(
      "panel_B",
      file_name,
      ".rds"
    )
  ),
)


## Panel C ====================================================================
## Panel C0: Annotations ======================================================

# Sim. mode
p_sim_mode <- ggplot2::ggplot(
  data = stats,
  mapping = ggplot2::aes(
    x = method,
    y = "Sim. mode",
    fill = sim_mode
  )
) +
ggplot2::geom_tile(
  color = NA,
  linewidth = 0.1,
  show.legend = TRUE
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::scale_fill_manual(
  name = "Sim. mode", 
  values = c(
    "mixed guides"  = RColorBrewer::brewer.pal(8, "Dark2")[3], 
    "unique guides" = RColorBrewer::brewer.pal(8, "Dark2")[4]
  ),
  limits = c(
    "mixed guides", 
    "unique guides"
  ),
  drop = FALSE
) +
ggplot2::labs(
  x = NULL,
  y = NULL
) +
theme_DK(
  show_facet_x = TRUE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_borders = FALSE,
  show_grid = FALSE,
  show_axis = FALSE
) +
ggplot2::theme(
  strip.text.x = ggplot2::element_text(
    angle = 90,
    hjust = 0,
    vjust = 0.5
  ),
  axis.title.y = ggplot2::element_text(
    angle = 0,
    vjust = 0.5,
    hjust = 0.5
  ),
  plot.margin = ggplot2::margin(0, 0, 0, 0, "pt")
)

# Workflows
p_workflow <- ggplot2::ggplot(
  data = stats,
  mapping = ggplot2::aes(
    x = method,
    y = "Workflow",
    fill = workflow
  )
) +
ggplot2::geom_tile(
  color = "black",
  linewidth = 0.1
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::scale_fill_brewer(
  name = "Workflow",
  palette = "Set2"
) +
ggplot2::labs(
  x = "Method",
  y = NULL
) +
theme_DK(
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_borders = FALSE,
  show_grid = FALSE,
  show_axis = FALSE
) +
ggplot2::theme(
  axis.title.y = ggplot2::element_text(
    angle = 0,
    vjust = 0.5,
    hjust = 0.5
  ),
  plot.margin = ggplot2::margin(0, 0, 0, 0, "pt")
)

# Method type
p_method_type <- ggplot2::ggplot(
  data = stats,
  mapping = ggplot2::aes(
    x = method,
    y = "Method type",
    fill = method_type
  )
) +
ggplot2::geom_tile(
  color = "black",
  linewidth = 0.1
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::scale_fill_brewer(
  name = "Method type",
  palette = "Set1"
) +
ggplot2::labs(
  x = NULL,
  y = NULL
) +
theme_DK(
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_borders = FALSE,
  show_grid = FALSE,
  show_axis = FALSE
) +
ggplot2::theme(
  axis.title.y = ggplot2::element_text(
    angle = 0,
    vjust = 0.5,
    hjust = 0.5
  ),
  plot.margin = ggplot2::margin(0, 0, 0, 0, "pt")
)

# Patchwork
p_fig_1_panel_C0 <- wrap_plots(
  p_sim_mode,
  p_workflow,
  p_method_type,
  ncol = 1
)

## Panel C1: Total FP =========================================================
p_fig_1_panel_C1 <- ggplot2::ggplot(
  data = stats, 
  mapping = ggplot2::aes(
    x = method, 
    y = FP
  )
) +
ggplot2::stat_summary(
  fun = "mean",
  geom = "col",
  position = ggplot2::position_dodge(
    width = 0.6
  ),
  width = 0.8,
  color = "black",
  fill = "grey80",
  linewidth = 0.1
) +
ggplot2::stat_summary(
  fun.data = "mean_sdl",
  fun.args = list(mult = 1),
  geom = "errorbar",
  position = ggplot2::position_dodge(
    width = 0.6
  ),
  color = "black",
  width = 0.2,
  linewidth = 0.2
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula), 
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Method", 
  y = "Total FPs"
) +
theme_DK(
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel C2: FP ratio =========================================================
p_fig_1_panel_C2 <- ggplot2::ggplot(
  data = stats, 
  mapping = ggplot2::aes(
    x = method, 
    y = fpr_ratio,
  )
) +
ggplot2::stat_summary(
  mapping = ggplot2::aes(
    fill = ggplot2::stage(
      start = NULL, 
      after_stat = ggplot2::after_stat(y >= 0)
    )
  ),
  fun = "mean",
  geom = "col",
  position = ggplot2::position_dodge(
    width = 0.6
  ),
  width = 0.8,
  linewidth = 0.1
) +
ggplot2::scale_fill_manual(
  name = "Dominant FPR",
  values = c(
    "TRUE" = wesanderson::wes_palette("Chevalier1")[1],
    "FALSE" = wesanderson::wes_palette("Chevalier1")[3]
  ),
  labels = c(
    "TRUE" = "Off-target", 
    "FALSE" = "Background"
  )
) +
ggplot2::geom_hline(
  yintercept = 0, 
  color = "black", 
  linewidth = 0.3
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula), 
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Method", 
  y = expression(log[2] ~ FPR[Off] ~ "/" ~ FPR[Bg])
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_grid = FALSE
)

## Patchwork & save ===========================================================
p_fig_1_panel_C <-
  p_fig_1_panel_C0 + 
  p_fig_1_panel_C1 + 
  p_fig_1_panel_C2 +
  patchwork::plot_layout(
    ncol = 1,
    heights = c(0.05, 0.05, 0.05, 0.3, 0.3),
    guides = "collect"
  ) 

# Save
ggplot2::ggsave(
  p_fig_1_panel_C,
  file = file.path(
    output_path,
    paste0(
      "panel_C",
      file_name,
      ".svg"
    )
  ),
  width = 8,
  height = 10,
  units = "cm"
)

saveRDS(
  p_fig_1_panel_C,
  file = file.path(
    output_path,
    paste0(
      "panel_C",
      file_name,
      ".rds"
    )
  ),
)


## Panel D ====================================================================
## Panel D0: Annotations ======================================================

# Sim. mode
p_sim_mode <- ggplot2::ggplot(
  data = stats,
  mapping = ggplot2::aes(
    x = nGuide,
    y = "Sim. mode",
    fill = sim_mode
  )
) +
ggplot2::geom_tile(
  color = NA,
  linewidth = 0.1,
  show.legend = TRUE
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::scale_fill_manual(
  name = "Sim. mode", 
  values = c(
    "mixed guides"  = RColorBrewer::brewer.pal(8, "Dark2")[3], 
    "unique guides" = RColorBrewer::brewer.pal(8, "Dark2")[4]
  ),
  limits = c(
    "mixed guides", 
    "unique guides"
  ),
  drop = FALSE
) +
ggplot2::labs(
  x = NULL,
  y = NULL
) +
theme_DK(
  show_facet_x = TRUE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_borders = FALSE,
  show_grid = FALSE,
  show_axis = FALSE
) +
ggplot2::theme(
  strip.text.x = ggplot2::element_text(
    angle = 90,
    hjust = 0,
    vjust = 0.5
  ),
  axis.title.y = ggplot2::element_text(
    angle = 0,
    vjust = 0.5,
    hjust = 0.5
  ),
  plot.margin = ggplot2::margin(0, 0, 0, 0, "pt")
)

# Patchwork
p_fig_1_panel_D0 <- wrap_plots(
  p_sim_mode,
  ncol = 1
)
## Panel D1: Off-target FPR ===================================================
p_fig_1_panel_D1 <- ggplot2::ggplot(
  data = stats, 
  mapping = ggplot2::aes(
    x = nGuide, 
    y = fpr_off_target,
    color = method,
    group = method
  )
) +
ggplot2::stat_summary(
  fun = mean, 
  geom = "point", 
  size = 0.1
) +
ggplot2::stat_summary(
  fun = mean, 
  geom = "line", 
  linewidth = 0.3
) +
ggplot2::scale_y_continuous(
  limits = c(0, 1)
  # labels = scales::percent_format(
  #   accuracy = 1
  # )
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula), 
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Cells per guide", 
  y = expression(FPR[Off]),
  color = "Method"
) +
theme_DK(
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel D2: Background FPR ===================================================
p_fig_1_panel_D2 <- ggplot2::ggplot(
  data = stats, 
  mapping = ggplot2::aes(
    x = nGuide, 
    y = fpr_background,
    color = method,
    group = method
  )
) +
ggplot2::stat_summary(
  fun = mean, 
  geom = "point", 
  size = 0.1
) +
ggplot2::stat_summary(
  fun = mean, 
  geom = "line", 
  linewidth = 0.3
) +
ggplot2::scale_y_continuous(
  limits = c(0, 1)
  # labels = scales::percent_format(
  #   accuracy = 1
  # )
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula), 
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Cells per targeting guide", 
  y = expression(FPR[Bg]),
  color = "Method"
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_grid = FALSE
)

## Patchwork & save ===========================================================
p_fig_1_panel_D <-
  p_fig_1_panel_D0 + 
  p_fig_1_panel_D1 + 
  p_fig_1_panel_D2 +
  patchwork::plot_layout(
    ncol = 1,
    heights = c(0.05, 0.3, 0.3),
    guides = "collect"
  ) 

# Save
ggplot2::ggsave(
  p_fig_1_panel_D,
  file = file.path(
    output_path,
    paste0(
      "panel_D",
      file_name,
      ".svg"
    )
  ),
  width = 9,
  height = 7,
  units = "cm"
)

saveRDS(
  p_fig_1_panel_D,
  file = file.path(
    output_path,
    paste0(
      "panel_D",
      file_name,
      ".rds"
    )
  ),
)


