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
  input_path   <- gsub(pattern = "\\\\ ", replacement = " ", trail_args[-1])
  output_path  <- "."

  # Params
  target_alpha <- 0.05

  # Output
  output_name <- paste0("_alpha_", target_alpha)

} else {
  # Script
  script_path <- rstudioapi::getActiveDocumentContext()$path
  script_name <- tools::file_path_sans_ext(basename(script_path))
  script_id   <- sub(pattern = "-.*", replacement = "", script_name)

  # Paths
  project_root <- "/home/scc_e_928652/ArcSimulationSynthetic"
  input_path   <- list.files(path = "results/nextflow/ArcSimSyn_offtarget_mixed/LabMeeting_2026-11-18", pattern = "_dge\\.rds$", recursive = TRUE, full.names = TRUE)
  output_path  <- "results/nextflow/ArcSimSyn_offtarget_mixed/LabMeeting_2026-11-18/15-visualize"

  stopifnot(dirname(dirname(input_path)) == dirname(output_path))
  if (!dir.exists(output_path)) dir.create(path = output_path, recursive = TRUE)

  # Params
  target_alpha <- 0.05

  # Output
  output_name <- paste0("_alpha_", target_alpha)
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


## Script info ================================================================
showf(fmt = "\n%s\n", line80(paste0(script_name, ".R")), type = "cat")


## Dataset ====================================================================
readDGE <- function(file_paths) {

  # Validate input files
  file_paths <- as.character(file_paths)
  file_paths <- file_paths[file.exists(file_paths)]
  if (length(file_paths) == 0) showf("No DGE files found in the input.", type = "stop")
  
  # Load files
  showf("Loading %s DGE files.", length(file_paths), type = "message")
  loaded_data <- vector("list", length(file_paths))
  bad_files   <- character()

  for (i in seq_along(file_paths)) {
    current_path <- file_paths[i]

    loaded_data[[i]] <- tryCatch(
      expr = as.data.frame(readRDS(current_path)),
      error = function(e) {
        bad_files <<- c(bad_files, current_path)
        showf(
          "Failed to load %s: %s",
          current_path,
          conditionMessage(e),
          type = "message"
        )
        NULL
      }
    )

    if (i %% 1000 == 0) {
      showf(
        "Loaded %s/%s files.",
        i,
        length(file_paths),
        type = "message"
      )
    }
  }

  if (length(bad_files) > 0) {
    showf(
      "%s DGE files could not be loaded.",
      length(bad_files),
      type = "message"
    )
  }

  # Keep successfully loaded data frames
  valid_data <- Filter(Negate(is.null), loaded_data)

  if (length(valid_data) == 0) {
    stop("None of the DGE files could be loaded successfully.")
  }

  # Standardize columns
  showf("Standardizing columns.", type = "message")
  all_col_lists  <- lapply(valid_data, colnames)
  union_cols     <- unique(unlist(all_col_lists))
  intersect_cols <- Reduce(intersect, all_col_lists)

  non_overlapping_cols <- setdiff(union_cols, intersect_cols)

  if (length(non_overlapping_cols) > 0) {
    showf(
      "Removing non-overlapping columns: %s",
      paste(non_overlapping_cols, collapse = ", "),
      type = "message"
    )
  }

  valid_data <- lapply(
    valid_data,
    function(df) df[, intersect_cols, drop = FALSE]
  )

  # Combine data frames
  dge <- do.call(rbind, valid_data)
  rownames(dge) <- NULL
  dge
}

.readDGE <- fun_with_ram(readDGE)
dge <- readDGE(input_path)


## Levels & order ============================================================
method_levels  <- c("limma", "edger", "deseq2", "dream", "nebula", "wilcox", "sceptre", "scmageck")
formula_levels <- c("~genotype", "~genotype + sample", "~genotype + (1|sample)", "~genotype + sample + (1|guide)", "~genotype + (1|sample) + (1|guide)")

dge$facLoc  <- factor(if ("facLoc"  %in% names(dge)) dge$facLoc  else "1")
dge$facOff  <- factor(if ("facOff"  %in% names(dge)) dge$facOff  else "1")
dge$koEff   <- factor(if ("koEff"   %in% names(dge)) dge$koEff   else "1")
dge$nGuide  <- factor(if ("nGuide"  %in% names(dge)) dge$nGuide  else "1")

dge$method  <- factor(if ("method"  %in% names(dge)) dge$method  else "1", levels = method_levels)
dge$formula <- factor(if ("formula" %in% names(dge)) dge$formula else "1", levels = formula_levels)


## Data Panel A, B, C & D =====================================================
vars_1     <- c("sim_mode", "workflow", "method_type", "method", "gene_qc", "formula", "time_min")
vars_2     <- c("facLoc", "facOff", "koEff", "nGuide")
group_vars <- c(vars_1, vars_2)

# Confusion matrix
stats <- aggregate(
  cbind(
    # TN: FDR bigger than alpha
    # FN: FDR bigger than alpha
    # TP: FDR smaller than alpha
    # FP: FDR smaller than alpha
    TN     = dge$status %in% c("Background", "Off-target")               & !is.na(dge$FDR) & dge$FDR >= target_alpha,
    FN     = dge$status %in% c("Knockout", "On-target", "Downstream")    & !is.na(dge$FDR) & dge$FDR >= target_alpha,
    TP     = dge$status %in% c("Knockout", "On-target", "Downstream")    & !is.na(dge$FDR) & dge$FDR < target_alpha,
    FP     = dge$status %in% c("Background", "Off-target")               & !is.na(dge$FDR) & dge$FDR < target_alpha,

    TN_off = dge$status == "Off-target"                                  & !is.na(dge$FDR) & dge$FDR >= target_alpha,
    TN_bg  = dge$status == "Background"                                  & !is.na(dge$FDR) & dge$FDR >= target_alpha,

    TP_ko  = dge$status == "Knockout"                                    & !is.na(dge$FDR) & dge$FDR < target_alpha,
    TP_on  = dge$status == "On-target"                                   & !is.na(dge$FDR) & dge$FDR < target_alpha,
    TP_dwn = dge$status == "Downstream"                                  & !is.na(dge$FDR) & dge$FDR < target_alpha,

    FP_off = dge$status == "Off-target"                                  & !is.na(dge$FDR) & dge$FDR < target_alpha,
    FP_bg  = dge$status == "Background"                                  & !is.na(dge$FDR) & dge$FDR < target_alpha,

    n_tested = as.integer(!is.na(dge$FDR))
  ),
  by  = dge[group_vars],
  FUN = sum
)

# Pct. DEGs per status
stats$pct_ko  <- stats$TP_ko  / (stats$TP + stats$FP)
stats$pct_on  <- stats$TP_on  / (stats$TP + stats$FP)
stats$pct_dwn <- stats$TP_dwn / (stats$TP + stats$FP)
stats$pct_off <- stats$FP_off / (stats$TP + stats$FP)
stats$pct_bg  <- stats$FP_bg  / (stats$TP + stats$FP)

status_long <- reshape(
  data = stats,
  varying = c("pct_ko", "pct_on", "pct_dwn", "pct_off", "pct_bg"),
  v.names = "pct",
  timevar = "status",
  times = c("Knockout", "On-target", "Downstream", "Off-target", "Background"),
  direction = "long"
)

# FPR
stats$fpr_total     <- ifelse(test = (stats$TN + stats$FP) == 0, yes = NA_real_, no = stats$FP / (stats$TN + stats$FP))
stats$fpr_off_share <- ifelse(test = (stats$TN + stats$FP) == 0, yes = NA_real_, no = stats$FP_off / (stats$TN + stats$FP))
stats$fpr_bg_share  <- ifelse(test = (stats$TN + stats$FP) == 0, yes = NA_real_, no = stats$FP_bg / (stats$TN + stats$FP))

# Ratio
stats$fpr_off   <- ifelse(test = (stats$TN_off + stats$FP_off) == 0, yes = NA_real_, no = stats$FP_off / (stats$TN_off + stats$FP_off))
stats$fpr_bg    <- ifelse(test = (stats$TN_bg  + stats$FP_bg)  == 0, yes = NA_real_, no = stats$FP_bg / (stats$TN_bg + stats$FP_bg))

stats$fpr_off_share <- ifelse(test = is.na(stats$fpr_off + stats$fpr_bg) | stats$fpr_off + stats$fpr_bg == 0, yes = NA_real_, no = stats$fpr_off / stats$fpr_off + stats$fpr_bg)
stats$fpr_bg_share  <- ifelse(test = is.na(stats$fpr_off + stats$fpr_bg) | stats$fpr_off + stats$fpr_bg == 0, yes = NA_real_, no = stats$fpr_bg / stats$fpr_off + stats$fpr_bg)

stats$fpr_ratio <- ifelse(test = !is.na(stats$fpr_off) & !is.na(stats$fpr_bg) & stats$fpr_off > 0 & stats$fpr_bg > 0, yes = log2(stats$fpr_off / stats$fpr_bg), no = NA_real_)

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
## Panel A00: Annotations =====================================================

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
  color = "black",
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
    "synthetic"   = RColorBrewer::brewer.pal(8, "Dark2")[3], 
    "permutation" = RColorBrewer::brewer.pal(8, "Dark2")[4]
  )
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
ggplot2::scale_fill_manual(
  name = "Workflow",
  values = c(
    "singlecell" = "#66C2A5",
    "pseudobulk" = "#FC8D62"
  )
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

# QC: gene filtering
p_gene_qc <- ggplot2::ggplot(
  data = stats,
  mapping = ggplot2::aes(
    x = method,
    y = "Gene filter",
    fill = gene_qc
  )
) +
ggplot2::geom_tile(
  color = "black",
  linewidth = 0.1
) +
ggplot2::scale_fill_manual(
  name = "Gene filter", 
  values = c(
    "TRUE"  = "#D4A017",
    "FALSE" = "darkred"
  ),
  labels = c(
    "TRUE"  = "method-specific",
    "FALSE" = "intersection"
  )
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
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
  plot.margin = ggplot2::margin(0, 0, 0, 0)
)

# Patchwork
p_panel_A00 <- wrap_plots(
  p_sim_mode,
  p_workflow,
  p_method_type,
  p_gene_qc,
  ncol = 1
)

## Panel A01: Tested ==========================================================
p_panel_A01 <- ggplot2::ggplot(
  data = stats,
  mapping = ggplot2::aes(
    x = method,
    y = n_tested
  )
) +
ggplot2::stat_summary(
  fun = "mean",
  geom = "col",
  color = "black",
  fill = "grey80", 
  linewidth = 0.1,
  width = 0.8,
  position = ggplot2::position_dodge(
    width = 0.6
  )
) +
ggplot2::stat_summary(
  fun.data = "mean_sdl",
  fun.args = list(mult = 1),
  geom = "errorbar",
  color = "black",
  width = 0.2,
  linewidth = 0.2,
  position = ggplot2::position_dodge(
    width = 0.6
  )
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Method",
  y = "Nr. tests"
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel A02: Gene status =====================================================
p_panel_A02 <- ggplot2::ggplot(
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
  name = "Gene response",
  values = colors_gene_status,
  breaks = c(
    "Knockout",
    "On-target",
    "Off-target",
    "Downstream",
    "Background"
  )[c(
    "Knockout",
    "On-target",
    "Off-target",
    "Downstream",
    "Background"
  ) %in% status_long$status[
    !is.na(status_long$pct) & status_long$pct > 0
  ]]
) +
ggplot2::labs(
  x = "Method", 
  y = "Pct. of pos."
) +
theme_DK(
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel A03: Power ===========================================================
p_panel_A03 <- ggplot2::ggplot(
  data = stats,
  mapping = ggplot2::aes(
    x = method, 
    y = sensitivity
  )
) +
ggplot2::stat_summary(
  data = local({
    winners <- aggregate(
      sensitivity ~ formula + method,
      data = stats,
      FUN = mean,
      na.rm = TRUE
    )

    winners$best <- factor(
      as.numeric(
        winners$sensitivity == max(winners$sensitivity, na.rm = TRUE)
      ),
      levels = c(0, 1),
      labels = c("FALSE", "TRUE")
    )

    merge(
      stats,
      winners[c("formula", "method", "best")],
      by = c("formula", "method")
    )
  }),
  mapping = ggplot2::aes(
    fill = best
  ),
  fun = mean,
  geom = "col",
  position = ggplot2::position_dodge(
    width = 0.6
  ),
  width = 0.8,
  color = "black",
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
ggplot2::scale_fill_manual(
  values = c(
    "FALSE" = "grey80", 
    "TRUE" = "pink"
  ),
  guide = "none"
) +
ggplot2::scale_y_continuous(
  limits = function(x) c(min(0, x), max(1, x)),
  breaks = seq(0, 1, by = 0.25)
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

## Panel A04: FDR =============================================================
p_panel_A04 <- ggplot2::ggplot(
  data = stats, 
  mapping = ggplot2::aes(
    x = method, 
    y = fdr_total
  )
) +
ggplot2::stat_summary(
  data = local({
    winners <- aggregate(
      fdr_total ~ formula + method,
      data = stats,
      FUN = mean,
      na.rm = TRUE
    )

    winners$best <- factor(
      as.numeric(
        winners$fdr_total == min(winners$fdr_total, na.rm = TRUE)
      ),
      levels = c(0, 1),
      labels = c("FALSE", "TRUE")
    )

    merge(
      stats,
      winners[c("formula", "method", "best")],
      by = c("formula", "method")
    )
  }),
  mapping = ggplot2::aes(
    fill = best
  ),
  fun = "mean",
  geom = "col",
  position = ggplot2::position_dodge(
    width = 0.6
  ),
  width = 0.8,
  color = "black",
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
ggplot2::scale_fill_manual(
  values = c(
    "FALSE" = "grey80", 
    "TRUE" = "pink"
  ),
  guide = "none"
) +
ggplot2::scale_y_continuous(
  limits = function(x) c(min(0, x), max(1, x)),
  breaks = seq(0, 1, by = 0.25)
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula), 
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Method", 
  y = "FDR"
) +
theme_DK(
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel A05: Higher metric ===================================================
p_panel_A05 <- ggplot2::ggplot(
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
  x = "Method",
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
theme_DK(
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel A06: Time (min) ======================================================
p_panel_A06 <- ggplot2::ggplot(
  data = aggregate(time_min ~ formula + method, data = stats, FUN = mean),
  mapping = ggplot2::aes(
    x = method,
    y = "Time (min)",
    fill = time_min
  )
) +
ggplot2::geom_tile(
  color = "black",
  linewidth = 0.1
) +
ggplot2::scale_fill_gradientn(
  name = "Time (min)",
  colours = c("#14213D", "#007F86", "#4CC9A6", "#F9C74F"),
  limits = c(0, NA),
  breaks = function(x) c(0, max(x)),
  labels = function(x) round(x, 1)
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Method",
  y = NULL
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_label_x = TRUE,
  show_ticks_x = TRUE,
  show_borders = FALSE,
  show_axis = FALSE,
  show_grid = FALSE
) +
ggplot2::theme(
  plot.margin = ggplot2::margin(0, 0, 0, 0)
)

## Patchwork & save ===========================================================
panel_dir <- file.path(output_path, "panel_A")
if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)

p_panel_A <-
  p_panel_A00 + 
  p_panel_A01 + 
  p_panel_A02 +
  p_panel_A03 +
  p_panel_A04 +
  p_panel_A05 +
  p_panel_A06 +
  patchwork::plot_layout(
    ncol = 1,
    heights = c(0.05, 0.05, 0.05, 0.05, 0.25, 0.3, 0.3, 0.3, 0.15, 0.05),
    guides = "collect"
  ) 

# Save
purrr::walk(
  c("svg", "pdf", "png"),
  \(fmt) ggplot2::ggsave(
    plot = p_panel_A,
    file = file.path(
      panel_dir,
      paste0(
        "panel_A", 
        output_name, 
        ".", fmt
      )
    ),
    width = 10,
    height = 12,
    units = "cm",
    dpi = 300
  )
)

saveRDS(
  object = p_panel_A,
  file = file.path(
    panel_dir,
    paste0(
      "panel_A",
      output_name,
      ".rds"
    )
  ),
)


## Panel B ====================================================================
## Panel B00: Annotations =====================================================

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
  color = "black",
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
    "synthetic"   = RColorBrewer::brewer.pal(8, "Dark2")[3], 
    "permutation" = RColorBrewer::brewer.pal(8, "Dark2")[4]
  )
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
p_panel_B00 <- wrap_plots(
  p_sim_mode,
  ncol = 1
)

## Panel B01: Power ===========================================================
p_panel_B01 <- ggplot2::ggplot(
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

## Panel B02: FDR =============================================================
p_panel_B02 <- ggplot2::ggplot(
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
ggplot2::scale_y_continuous(
  limits = c(0, 1)
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula), 
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Cells per targeting guide", 
  y = "FDR",
  color = "Method"
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_grid = FALSE
)

## Patchwork & save ===========================================================
panel_dir <- file.path(output_path, "panel_B")
if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)

p_panel_B <-
  p_panel_B00 + 
  p_panel_B01 + 
  p_panel_B02 +
  patchwork::plot_layout(
    ncol = 1,
    heights = c(0.05, 0.3, 0.3),
    guides = "collect"
  ) 

# Save
purrr::walk(
  c("svg", "pdf", "png"),
  \(fmt) ggplot2::ggsave(
    plot = p_panel_B,
    file = file.path(
      panel_dir,
      paste0(
        "panel_B", 
        output_name, 
        ".", fmt
      )
    ),
    width = 9,
    height = 7,
    units = "cm",
    dpi = 300
  )
)

saveRDS(
  object = p_panel_B,
  file = file.path(
    panel_dir,
    paste0(
      "panel_B",
      output_name,
      ".rds"
    )
  )
)


if (!all(is.na(stats$fpr_ratio))) {

## Panel C ====================================================================
## Panel C00: Annotations =====================================================

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
  color = "black",
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
    "synthetic"   = RColorBrewer::brewer.pal(8, "Dark2")[3], 
    "permutation" = RColorBrewer::brewer.pal(8, "Dark2")[4]
  )
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
ggplot2::scale_fill_manual(
  name = "Workflow",
  values = c(
    "singlecell" = "#66C2A5",
    "pseudobulk" = "#FC8D62"
  )
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
p_panel_C00 <- wrap_plots(
  p_sim_mode,
  p_workflow,
  p_method_type,
  ncol = 1
)

## Panel C01: Tested ==========================================================
p_panel_C01 <- ggplot2::ggplot(
  data = stats,
  mapping = ggplot2::aes(
    x = method,
    y = n_tested
  )
) +
ggplot2::stat_summary(
  fun = "mean",
  geom = "col",
  color = "black",
  fill = "grey80", 
  linewidth = 0.1,
  width = 0.8,
  position = ggplot2::position_dodge(
    width = 0.6
  )
) +
ggplot2::stat_summary(
  fun.data = "mean_sdl",
  fun.args = list(mult = 1),
  geom = "errorbar",
  color = "black",
  width = 0.2,
  linewidth = 0.2,
  position = ggplot2::position_dodge(
    width = 0.6
  )
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Method",
  y = "Nr. tests"
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel C02: FPR =============================================================
p_panel_C02 <- ggplot2::ggplot(
  data = stats, 
  mapping = ggplot2::aes(
    x = method, 
    y = fpr_total
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
  y = "FPR"
) +
theme_DK(
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel C03: FPR contribution ================================================
p_panel_C03 <- ggplot2::ggplot(
  data = reshape(
    data = stats,
    varying = c("fpr_bg_share", "fpr_off_share"),
    v.names = "value",
    timevar = "status",
    times = c("Background", "Off-target"),
    direction = "long"
  ),
  mapping = ggplot2::aes(
    x = method,
    y = value,
    fill = status
  )
) +
ggplot2::geom_col(
  position = "stack",
  width = 0.8,
  linewidth = 0.1
) +
ggplot2::scale_y_continuous(
  breaks = seq(0, 1, by = 0.25),
  limits = c(0, 1)
) +
ggplot2::scale_fill_manual(
  name = "Gene response",
  values = c(
    "Background" = wesanderson::wes_palette("Chevalier1")[3],
    "Off-target" = wesanderson::wes_palette("Chevalier1")[1]
  )
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Method",
  y = "Relative\ncontribution"
) +
theme_DK(
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel C04: FPR ratio =======================================================
p_panel_C04 <- ggplot2::ggplot(
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
panel_dir <- file.path(output_path, "panel_C")
if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)

p_panel_C <-
  p_panel_C00 + 
  p_panel_C01 + 
  p_panel_C02 +
  p_panel_C03 +
  # p_panel_C04 +
  patchwork::plot_layout(
    ncol = 1,
    heights = c(0.05, 0.05, 0.05, 0.25, 0.3, 0.3),
    guides = "collect"
  ) 

# Save
purrr::walk(
  c("svg", "pdf", "png"),
  \(fmt) ggplot2::ggsave(
    plot = p_panel_C,
    file = file.path(
      panel_dir,
      paste0(
        "panel_C", 
        output_name, 
        ".", fmt
      )
    ),
    width = 8,
    height = 10,
    units = "cm",
    dpi = 300
  )
)

saveRDS(
  object = p_panel_C,
  file = file.path(
    panel_dir,
    paste0(
      "panel_C",
      output_name,
      ".rds"
    )
  )
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
    "synthetic"   = RColorBrewer::brewer.pal(8, "Dark2")[3], 
    "permutation" = RColorBrewer::brewer.pal(8, "Dark2")[4]
  )
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
p_panel_D0 <- wrap_plots(
  p_sim_mode,
  ncol = 1
)
## Panel D1: Off-target FPR ===================================================
p_panel_D1 <- ggplot2::ggplot(
  data = stats, 
  mapping = ggplot2::aes(
    x = nGuide, 
    y = fpr_off,
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
p_panel_D2 <- ggplot2::ggplot(
  data = stats, 
  mapping = ggplot2::aes(
    x = nGuide, 
    y = fpr_bg,
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
panel_dir <- file.path(output_path, "panel_D")
if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)

p_panel_D <-
  p_panel_D0 + 
  p_panel_D1 + 
  p_panel_D2 +
  patchwork::plot_layout(
    ncol = 1,
    heights = c(0.05, 0.3, 0.3),
    guides = "collect"
  ) 

# Save
purrr::walk(
  c("svg", "pdf", "png"),
  \(fmt) ggplot2::ggsave(
    plot = p_panel_D,
    file = file.path(
      panel_dir,
      paste0(
        "panel_D", 
        output_name, 
        ".", fmt
      )
    ),
    width = 9,
    height = 7,
    units = "cm",
    dpi = 300
  )
)

saveRDS(
  object = p_panel_D,
  file = file.path(
    panel_dir,
    paste0(
      "panel_D",
      output_name,
      ".rds"
    )
  ),
)


}

## Script info ================================================================
showf("\n◷ %s ====================================================\n", date(), type = "cat")
