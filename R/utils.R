## max. code width ============================================================

## Title: Arc Institute in vivo pertubation screens
## Author: DKatzlberger
## Date: 21.01.2026

## Log/Printing ===============================================================

showf <- function(fmt, ..., type = c(
  "message", "print", "cat", "warning", "stop"
)) {
  type   <- match.arg(type)
  txt    <- sprintf(fmt, ...)
  caller <- sys.call(-1L)

  switch(
    type,
    message = message(txt),
    print   = print(txt),
    cat     = cat(txt),
    warning = warning(simpleWarning(txt, call = caller)),
    stop    = stop(simpleError(txt, call = caller))
  )

  invisible(txt)
}

fun_with_ram <- function(fun) {

  # Helper: extend '-' to line 80
  .line80 <- function(x) {
    paste(x, strrep("-", max(1, 78 - nchar(x))))
  }

  # Helper: job id  
  .scc_jobid <- function() {
    cg  <- readLines("/proc/self/cgroup", warn = FALSE)
    hit <- grep("job_[0-9]+", cg, value = TRUE)
    if (length(hit) == 0) return(NA_character_)
    sub(".*job_([0-9]+).*", "\\1", hit[1])
  }

  # Helper: Get system cgroup memory
  .get_system_mem <- function() {
    job_id  <- .scc_jobid()
    limit_gb <- NA_real_
    used_gb  <- NA_real_
    
    if (!is.na(job_id) && job_id != "") {
      x <- tryCatch(
        system(sprintf("scontrol show job %s", job_id), intern = TRUE),
        error = function(e) NULL
      )
      if (!is.null(x)) {
        line <- grep("AllocTRES=", x, value = TRUE)
        if (length(line) > 0) {
          limit_str <- sub(".*mem=([^,\\s]+).*", "\\1", line, perl = TRUE)
          
          num_val <- as.numeric(gsub("[^0-9.]", "", limit_str))
          unit    <- toupper(gsub("[0-9.]", "", limit_str))
          
          if (grepl("G", unit)) {
            limit_gb <- num_val          
          } else if (grepl("M", unit)) {
            limit_gb <- num_val / 1024   
          } else if (grepl("T", unit)) {
            limit_gb <- num_val * 1024  
          } else {
            limit_gb <- num_val         
          }
        }
      }
    }
    
    if (file.exists("/proc/self/status")) {
      status <- readLines("/proc/self/status", warn = FALSE)
      vsize  <- grep("^VmRSS:", status, value = TRUE)
      if (length(vsize) > 0) {
        kb <- as.numeric(gsub("[^0-9]", "", vsize))
        used_gb <- kb / (1024^2) # KB -> MB -> GB
      }
    }
    
    list(limit = limit_gb, usage = used_gb)
  }

  # Extract fun. name
  fun_name <- deparse(substitute(fun))
  fun_name <- paste("fun.:", fun_name)

  force(fun)
  force(fun_name)

  function(...) {

    args_list <- list(...)
    
    if (!is.null(names(args_list)) && "info" %in% names(args_list)) {
      info <- args_list[["info"]]
      args_list[["info"]] <- NULL
    } else {
      info <- NULL
    }

    cnd <- rlang::catch_cnd({
      t0 <- Sys.time()
      messages <- character()
      warnings <- character()
      std_logs <- character()  

      mem <- peakRAM::peakRAM({
        std_logs <- capture.output({
          visible_output <- withCallingHandlers(
            withVisible(do.call(fun, args_list, envir = parent.frame())),

            message = function(m) {
              messages <<- c(messages, conditionMessage(m))
              invokeRestart("muffleMessage")
            },

            warning = function(w) {
              warnings <<- c(warnings, conditionMessage(w))
              invokeRestart("muffleWarning")
            }
          )
        }, type = "output")
      })
      t1 <- Sys.time()

      output_obj <- visible_output$value

      # Stats
      time <- as.numeric(difftime(t1, t0, units = "secs"))
      size <- as.numeric(object.size(output_obj)) / 1024^2
      peak <- as.numeric(mem[[grep("^Peak", names(mem), value = TRUE)[1]]][1])

      # Read allocated system memory
      sys_mem <- .get_system_mem()

      # Report
      showf("\n\n", type = "cat")
      showf("%s\n", .line80(fun_name), type = "cat")
      if (!is.null(info)) showf("info: %s\n", info, type = "cat")
      showf("- time: %.2f seconds\n", time, type = "cat")
      showf("- size: %.2f MB\n", size, type = "cat")
      showf("- peak: %.2f MB\n", peak, type = "cat")
      showf("- used: %.2f GB (%.2f GB)\n", sys_mem$usage, sys_mem$limit, type = "cat")
      
      # Std. out
      std_logs <- std_logs[nchar(trimws(std_logs)) > 0]
      if (length(std_logs)) {
        showf("\noutput:\n", type = "cat")
        showf("%s\n", paste(std_logs, collapse = "\n"), type = "cat")
      }

      # Messages
      if (length(messages)) {
        showf("\nmessages:\n", type = "cat")
        showf("%s", paste(messages, collapse = ""), type = "cat")
      }

      # Warnings
      if (length(warnings)) {
        showf("\nwarnings:\n", type = "cat")
        showf("%s\n", paste(warnings, collapse = "\n"), type = "cat")
      }

      showf("%s\n", .line80("✓ executed"), type = "cat")

      # Return
      if (visible_output$visible) {
        return(output_obj)
      } else {
        return(invisible(output_obj))
      }
    }, classes = "error")

    if (!is.null(cnd)) {
      showf("\n\n", type = "cat")
      showf("%s\n", .line80(fun_name), type = "cat")
      if (!is.null(info)) showf("info: %s\n", info, type = "cat")
      
      tb <- rlang::trace_back(cnd)
      showf("\nerror:\n", type = "cat")
      showf("%s\n", conditionMessage(cnd), type = "cat")

      showf("\ntraceback:\n", type = "cat")
      tb     <- rlang::trace_back(cnd)
      std_tb <- capture.output(print(tb))
      showf("%s\n", paste(std_tb, collapse = "\n"), type = "cat")
      
      showf("%s\n", .line80("✗ failed"), type = "cat")

      # Script stop
      old_ops <- options(show.error.messages = FALSE)
      on.exit(options(old_ops))
      stop(structure(list(message = ""), class = c("simpleError", "error", "condition")))
    }
  }
}

## Environment ================================================================
scc_status <- function(verbose = TRUE) {

  .scc_jobid <- function() {
    cg <- readLines("/proc/self/cgroup", warn = FALSE)
    hit <- grep("job_[0-9]+", cg, value = TRUE)
    if (length(hit) == 0) return(NA_character_)
    sub(".*job_([0-9]+).*", "\\1", hit[1])
  }

  .slurm_mem <- function(jobid) {
    if (is.na(jobid)) return(NA_character_)
    x <- system(
      sprintf("scontrol show job %s", jobid),
      intern = TRUE
    )
    line <- grep("AllocTRES=", x, value = TRUE)
    sub(".*mem=([^,]+).*", "\\1", line)
  }

  # User & node
  job_id   <- .scc_jobid()
  username <- Sys.getenv("USER")
  node     <- Sys.info()[["nodename"]]
  n_cores  <- as.integer(system("nproc", intern = TRUE))
  memory   <- .slurm_mem(job_id)

  if (verbose){
    cat(sprintf("%s@%s\n", username, node))
    cat(sprintf("- scc job id: %s\n", job_id))
    cat(sprintf("- cores available: %d\n", n_cores))
    cat(sprintf("- mem available: %s\n", memory))
  }

  invisible(
    list(
      job_id = job_id,
      n_cores = n_cores,
      mem = .slurm_mem(job_id)
    )
  )
}

script_startup <- function(
  root,
  renv = FALSE,
  bioc = FALSE,
  packages = NULL
){

  # renv
  t0 <- Sys.time()
  if (renv) {

    # Environment variables
    csbg_mnt <- "/mnt/ceph/groups_hdd/SCCGroup/computational_systems_biology_group/dkatzlberger"
    Sys.setenv(RENV_PATHS_RENV = ".renv")
    Sys.setenv(RENV_PATHS_ROOT = file.path(csbg_mnt, ".renv"))
    Sys.setenv(RENV_PATHS_CACHE = file.path(csbg_mnt, ".renv/cache"))
    Sys.setenv(RENV_PATHS_SANDBOX = file.path(csbg_mnt, ".renv/sandbox"))

    Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
    Sys.setenv(RENV_CONFIG_CACHE_SYMLINKS = "TRUE")
    Sys.setenv(RENV_CONFIG_INSTALL_TRANSACTIONAL = "FALSE")
    Sys.setenv(RENV_CONFIG_AUTO_SNAPSHOT = "FALSE")
    Sys.setenv(RENV_CONFIG_PKG_TYPES = "source")

    # Activate renv
    renv::load(project = root)

    project <- tryCatch(renv::project(), error = function(e) NULL)
    version <- as.character(packageVersion("renv"))

    # BiocManager 
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        renv::install("BiocManager", prompt = FALSE)
        # BiocManager repositories 
        suppressMessages(options(repos = BiocManager::repositories()))
        renv::install("BiocVersion", prompt = FALSE)
      }
    }
  }
  t1 <- Sys.time()

  # Libraries (only if provided)
  if (!is.null(packages)) {

    if (!is.character(packages)) showf("`packages` must be a character vector\n", type = "stop")

    suppressWarnings(
      suppressPackageStartupMessages({
        for (pkg in unique(packages)) {
          library(pkg, character.only = TRUE, quietly = TRUE)
        }
      })
    )
  }

  t2 <- Sys.time()
  
  # Report
  renv_time <- as.numeric(difftime(t1, t0, units = "secs"))
  lib_time  <- as.numeric(difftime(t2, t1, units = "secs"))
  showf("- renv activation: %.2f seconds\n", renv_time, type = "cat")
  showf("- lib. attachment: %.2f seconds\n", lib_time, type = "cat")

  # Info
  showf("NOTE: R runs through pixi. Make sure all variables are set!\n", type = "cat")
  showf("NOTE: renv sandbox is deactivated!\n", type = "cat")
  
  
  # Return
  invisible(TRUE)
}

## CRISPR =====================================================================
assign_guide_by_max_count <- function(
  obj,
  assay = "CRISPR",
  ntc_pattern = "^NTCms|^NonTargeting|^BRDN|^NTC",
  verbose = TRUE
) {

  # Extract crispr counts
  counts <- obj@assays[[assay]]$counts

  # Assign guides
  obj$guide <- apply(counts, 2, function(x) {
    if (max(x) == 0) return(NA)
    rownames(counts)[which.max(x)]
  })

  # Find cells with NTC pattern
  is_ntc <- grepl(ntc_pattern, obj$guide, ignore.case = TRUE)

  # Assign new columns
  obj$genotype <- ifelse(is_ntc, "NTC", sub("-.*$", "", obj$guide))
  obj$on_guide <- ifelse(is_ntc, "NTC", obj$guide)

  # Validate guide assignment 
  valid_cells <- !is.na(obj$guide)
  genotypes_per_guide <- tapply(obj$genotype[valid_cells], obj$guide[valid_cells], function(x) length(unique(x)))
  problematic_guides  <- names(genotypes_per_guide[genotypes_per_guide > 1])

  if (length(problematic_guides) > 0) {
    showf("Identical guides were assigned to different genotypes", type = "stop")
  }
  
  # Verbose
  n_valid_cells <- sum(!is.na(obj$guide))
  n_na_cells    <- sum(is.na(obj$guide))
  n_genotypes   <- length(unique(obj$genotype)) 

  if (verbose) {
    showf("Assigned %s unique genotypes to %s cells.", n_genotypes, n_valid_cells, type = "message")
    if (n_na_cells > 0) showf("NOTE: %s cells couldn't be assigned.", n_na_cells, type = "message")
  }

  # Return
  obj
}

make_guide_target_df <- function(
  data,
  guide_col,
  target_col, 
  ref,
  verbose = TRUE
) {

  data <- as.data.frame(data)
  
  # Unique grna targets
  obj <- unique(data[, c(guide_col, target_col)])
  colnames(obj) <- c("guide_id", "guide_target")
  
  # Targeting vs. non-targeting
  obj$guide_target <- ifelse(
    obj$guide_target == ref, 
    "non-targeting", 
    sprintf("%s%s", target_col, obj$guide_target)
  )

  # Verbose
  n_guides  <- length(unique(obj$guide_id))
  n_targets <- length(unique(obj$guide_target))

  if (verbose) {
    showf("Found %d unique guides mapped to %d unique genotypes.", n_guides, n_targets, type = "message")
    showf("Recoding reference value '%s' to 'non-targeting'.", ref, type = "message")
  }

  # Return
  rownames(obj) <- NULL
  obj
}

make_pseudo_guide_matr <- function(
  data,
  guide_col,
  pseudocount = 1,
  verbose = TRUE
){

  data <- as.data.frame(data)

  cell_ids  <- rownames(data)
  guide_ids <- unique(data[[guide_col]])

  # Create matrix
  obj <- sapply(guide_ids, function(g) {
    as.numeric(data[[guide_col]] == g)
  })

  colnames(obj) <- guide_ids
  rownames(obj) <- cell_ids

  # Add pseudocount
  obj <- obj * pseudocount
  obj <- t(obj)

  # Verbose
  n_cells  <- length(cell_ids)
  n_guides <- length(guide_ids)

  if (verbose) {
    showf("Created pseudocounts for %d guides and %d cells.", n_guides, n_cells, type = "message")
    showf("Added a pseudocount of %d to assigned cells.", pseudocount, type = "message")
  }

  # Return
  obj
}

make_extra_cov_df <- function(
  data,
  target_col
){

  data <- as.data.frame(data)

  # Create extra covariate df
  obj <- data[, target_col, drop = FALSE]
  
  # Return
  obj
}

## scDesign3 ==================================================================
.is_valid_fit <- function(fit) {
  !is.null(fit) && inherits(fit, c("gamlss", "gam", "glm", "lm"))
}

.model_family <- function(
  fit
){

  # Extract model type
  fam <- tryCatch({
    if (inherits(fit, "gamlss")) {
      if (startsWith(toupper(fit$family[1]), "ZINB")) {
        "zinb"
      } else if (startsWith(toupper(fit$family[1]), "ZIP")) {
        "zip"
      } else {
        "nb"
      }
    } else if (!is.null(fit$family$family) && grepl("Negative Binomial", fit$family$family)) {
      "nb"
    } else if (inherits(fit, "lm")) {
      "lm"
    } else {
      class(fit)[1]
    }
  }, error = function(e) "unknown")

  # Return
  fam
}

.get_predictors <- function(fit) {

  coefs <- tryCatch(coef(fit), error = function(e) NULL)
  coefs <- coefs[setdiff(names(coefs), "(Intercept)")]
  if (!length(coefs)) return("'1'")

  quoted_names <- paste0("'", names(coefs), "'")
  padded_names <- sprintf("%-10s", quoted_names)

  txt <- paste(sprintf("%s %-6.3f", padded_names, coefs), collapse = "  ")

  # Return
  txt
}

estimate_marginals <- function(
  sce_obj,
  assay_use,
  celltype = "celltype",
  covariates,
  mu_formula,
  sigma_formula,
  family_use,
  corr_by,
  n_cores,
  verbose = TRUE
){

  # Input
  input <- scDesign3::construct_data(
    sce = sce_obj,
    assay_use = assay_use,
    celltype = celltype,
    pseudotime = NULL,
    spatial = NULL,
    other_covariates = covariates,
    corr_by = corr_by
  )
  
  # Marginals
  marginals <- scDesign3::fit_marginal(
    data = input,
    predictor = "gene",
    mu_formula = mu_formula,
    sigma_formula = sigma_formula,
    family_use = family_use,
    n_cores = n_cores,
    usebam = TRUE
  )
   
  # Verbose
  if (verbose) { 
    for (g in names(marginals)) {

      if (!.is_valid_fit(marginals[[g]]$fit)) {
        message(sprintf("%-15s -> failed", g)) 
        next
      }

      fam <- .model_family(marginals[[g]]$fit)
      pre <- .get_predictors(marginals[[g]]$fit)
      message(sprintf("%-15s -> fam: %-10s | covariate: %s", g, fam, pre))
    }
  }

  # Return
  marginals
}

estimate_correlation <- function(
  sce_obj,
  celltype = "celltype",
  marginals,
  covariates,
  corr_by,
  copula,
  family_use,
  n_cores = 1,
  verbose = TRUE
){

  # Input
  input <- scDesign3::construct_data(
    sce = sce_obj,
    assay_use = "counts",
    celltype = celltype,
    pseudotime = NULL,
    spatial = NULL,
    other_covariates = covariates,
    corr_by = corr_by
  )

  # Copula
  copula <- scDesign3::fit_copula(
    sce = sce_obj,
    assay_use = "counts",
    input_data = input$dat,
    marginal_list = marginals,
    family_use = family_use,
    copula = copula,
    DT = FALSE,
    n_cores = n_cores
  )

  # Verbose
  if (verbose) { 
    for (g in names(marginals)) {

      if (!.is_valid_fit(marginals[[g]]$fit)) {
        message(sprintf("%-15s -> failed", g)) 
        next
      }

      fam <- .model_family(marginals[[g]]$fit)
      pre <- .get_predictors(marginals[[g]]$fit)
      message(sprintf("%-15s -> fam: %-10s | covariate: %s", g, fam, pre))
    }
  }

  # Return
  copula
}

modify_covariate <- function(
  sce_obj,
  covariates,
  celltype = "celltype",
  corr_by,
  ko_eff,
  mode = c("all", "partial"),
  seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  mode <- match.arg(mode)

  # Construct input data
  input <- scDesign3::construct_data(
    sce = sce_obj,
    assay_use = "counts",
    celltype = celltype,
    other_covariates = covariates,
    pseudotime = NULL,
    spatial = NULL,
    corr_by = corr_by
  )

  new_covariate <- input$dat

  # Identify KO and NTC cells
  ko_idx  <- which(new_covariate$K == 1)
  ntc_idx <- which(new_covariate$K == 0)

  # No modification
  if (is.null(ko_eff) || ko_eff == 100) {
    showf("No modification applied (ko_eff = %s).", if (is.null(ko_eff)) "NULL" else "100", type = "message")
    return(new_covariate)
  }

  # Validate efficiency
  if (!is.numeric(ko_eff) ||
      length(ko_eff) != 1 ||
      is.na(ko_eff) ||
      ko_eff < 0 ||
      ko_eff > 100) {
    showf("ko_eff must be NULL or a single number between 0 and 100.", type = "stop")
  }

  # Ensure both groups are present
  if (length(ko_idx) == 0 || length(ntc_idx) == 0) {
    showf("Both KO (K = 1) and NTC (K = 0) cells must be present.", type = "stop")
  }

  if (mode == "all") {
    showf("'ALL' mode: Modifying all %d KO cells using NTC counts scaled by %.0f.", length(ko_idx), 1 - ko_eff / 100, type = "message")

    # PAPER BEHAVIOR:
    # Modify all KO cells using randomly sampled NTC counts.
    # Scale counts according to ko_eff.
    sampled_ntc <- sample(
      new_covariate$C[ntc_idx],
      size = length(ko_idx),
      replace = TRUE
    )

    new_covariate$C[ko_idx] <- round(
      sampled_ntc * (1 - ko_eff / 100)
    )

  } else if (mode == "partial") {
    # NEW BEHAVIOR:
    # Calculate number of KO cells to modify
    n_modify <- round(length(ko_idx) * ko_eff / 100)

    showf("'PARTIAL' mode: Modifying %d of %d KO cells (%.1f%%).", n_modify, length(ko_idx), ko_eff, type = "message")

    # Randomly select KO cells to modify
    selected_ko <- sample(
      ko_idx,
      size = n_modify,
      replace = FALSE
    )

    # Assign randomly sampled NTC counts directly
    sampled_ntc <- sample(
      new_covariate$C[ntc_idx],
      size = n_modify,
      replace = TRUE
    )

    new_covariate$C[selected_ko] <- sampled_ntc
  }

  new_covariate
}

extract_params <- function(
  sce_obj,
  covariates,
  celltype = "celltype",
  marginals,
  family_use,
  corr_by,
  new_covariate,
  n_cores = 1
){

  # Input
  input <- scDesign3::construct_data(
    sce = sce_obj,
    assay_use = "counts",
    celltype = celltype,
    other_covariates = covariates,
    pseudotime = NULL,
    spatial = NULL,
    corr_by = corr_by
  )

  # Model params
  params <- scDesign3::extract_para(
    sce = sce_obj,
    marginal_list = marginals,
    family_use = family_use,
    new_covariate = new_covariate,
    data = input$dat,
    n_cores = n_cores
  )

  # Return
  params
}

.simulate_counts <- function(
  sce_obj,
  covariates,
  celltype = "celltype",
  params,
  copulas,
  corr_by,
  family_use,
  new_covariate,
  seed = NULL,
  n_cores = 1
){

  if (!is.null(seed)) set.seed(seed)

  # Input
  input <- scDesign3::construct_data(
    sce = sce_obj,
    assay_use = "counts",
    celltype = celltype,
    other_covariates = covariates,
    pseudotime = NULL,
    spatial = NULL,
    corr_by = corr_by
  )

  # Counts
  counts <- scDesign3::simu_new(
    sce = sce_obj,
    mean_mat = params$mean_mat,
    sigma_mat = params$sigma_mat,
    zero_mat = params$zero_mat,
    copula_list = copulas$copula_list,
    family_use = family_use,
    new_covariate = new_covariate,
    input_data = input$dat,
    filtered_gene = NULL,
    n_cores = n_cores
  )

  # Return
  counts
}

sim_pipe <- function(
  sce_obj,
  covariates,
  celltype = "celltype",
  params,
  copulas,
  corr_by,
  family_use,
  new_covariate,
  knockout,
  n_cores = 1,
  seed = NULL
){

  # Input
  if (is.null(new_covariate)) {

    input <- scDesign3::construct_data(
      sce = sce_obj,
      assay_use = "counts",
      celltype = celltype,
      other_covariates = covariates,
      pseudotime = NULL,
      spatial = NULL,
      corr_by = corr_by
    )

    meta <- input$dat
    C    <- sce_obj$C
  } else {

    meta <- new_covariate
    C    <- new_covariate$C
  }

  # Simulate count matrix
  counts <- .simulate_counts(
    sce_obj = sce_obj,
    celltype = celltype,
    covariates = covariates,
    params = params,
    copulas = copula,
    corr_by = corr_by,
    family_use = family_use,
    new_covariate = new_covariate,
    seed = seed,
    n_cores = n_cores
  )

  # Add knockout expr
  C <- Matrix::Matrix(
    matrix(C, nrow = 1, dimnames = list(knockout, colnames(counts))),
    sparse = TRUE
  )

  counts <- rbind(counts, C)

  # Unique cell barcodes
  if (!is.null(seed)) set.seed(seed + sample(1:10000, 1))

  n_cells     <- ncol(counts)
  nucleotides <- c("A", "C", "G", "T")

  barcodes <- replicate(n_cells, {
    paste0(sample(nucleotides, 18, replace = TRUE), collapse = "")
  })
  barcodes <- paste0(barcodes, "-1")

  colnames(counts) <- barcodes

  # colData
  col_data <- data.frame(
    celltype = meta$celltype,
    sample = meta$I,
    stringsAsFactors = FALSE
  )

  rownames(col_data) <- barcodes

  # Single cell obj.
  sce_obj <- SingleCellExperiment::SingleCellExperiment(
    assays  = list(counts = counts),
    colData = col_data
  )

  # Return
  sce_obj
}

## DGE analysis ===============================================================
parse_model_formula <- function(
  formula,
  verbose = TRUE
) {

  # From string to formula
  f <- as.formula(formula)

  # Random effects
  bars <- reformulas::findbars(f)
  if (is.null(bars)) {
    rand_vars <- NULL
  } else {
    rand_vars <- unique(unname(unlist(lapply(bars, function(b) all.vars(b[[3]])))))
  }
  
  # Fixed effects
  fixed_vars <- labels(terms(reformulas::nobars(f)))
  if (length(fixed_vars) > 0) {
    target_covariate <- fixed_vars[1]
    extra_covariates <- if (length(fixed_vars) > 1) fixed_vars[-1] else NULL
  } else {
    target_covariate <- NULL
    extra_covariates <- NULL
  }

  # Verbose
  if (verbose) {
    showf("Formula: %s", formula, type = "message")
    showf("Target covariate: %s", target_covariate, type = "message")
    showf("Extra covariates: %s", paste(extra_covariates, collapse = ", "), type = "message")
    showf("Random covariate: %s", paste(rand_vars, collapse = ", "), type = "message")
  }
  
  # Return
  list(
    target_covariate  = target_covariate,
    extra_covariates  = extra_covariates,
    random_covariates = rand_vars
  )
}

validate_dge <- function(
  x
){
  
  # Check
  if (nrow(x) == 0) {
    showf("Statistic has zero rows. Maybe wrong extraction fun.", type = "stop")
  }

  # Universal req. columns
  required_cols <- c("coef", "gene", "logFC", "se", "pvalue")
  missing_cols  <- setdiff(required_cols, colnames(x))

  if (length(missing_cols) > 0) {
    showf(
      paste0("Statistic is missing req. column(s): ", 
      paste(missing_cols, collapse = ", ")), 
      type = "stop"
    )
  }

  # Check for NA columns
  required_cols <- c("pvalue")
  na_cols       <- colSums(is.na(x[required_cols]))

  if (any(na_cols > 0)) {
    showf(
      paste0("NOTE: Statistic has NA in column(s): ", 
      paste(names(na_cols)[na_cols > 0], collapse = ", ")), 
      type = "message"
    )
  }

  # Return
  invisible(TRUE)
}

extract_nebula <- function(
  fit,
  intercept = FALSE,
  adjust_p = TRUE,
  adjust_m = "BH"
){

  stats <- fit$summary

  # Identify coefficients
  coef_names <- unique(sub(
      "^(logFC|se|p)_",
      "",
      grep("^(logFC|se|p)_", colnames(stats), value = TRUE)
  ))

  # Build long-format table
  out <- do.call(
    rbind,
    lapply(coef_names, function(coef){

      data.frame(
        coef = coef,
        gene = stats$gene,
        logFC = stats[[paste0("logFC_", coef)]],
        se = stats[[paste0("se_", coef)]],
        pvalue = stats[[paste0("p_", coef)]],
        stringsAsFactors = FALSE
      )
    })
  )

  # Optionally remove intercept
  if (!intercept) out <- out[out$coef != "(Intercept)", ]

  # Optionally adjust p-values
  if (adjust_p) {
    out$FDR <- ave(
      out$pvalue,
      out$coef,
      FUN = function(z){ p.adjust(z, method = adjust_m)}
    )
  }

  # Validate
  rownames(out) <- NULL
  validate_dge(out)

  # Return
  out
}

extract_sceptre <- function(
  fit,
  intercept = FALSE,
  adjust_p = TRUE,
  adjust_m = "BH"
) {

  # Extract the results table
  stats <- if (inherits(fit, "sceptre_object")) {
    sceptre::get_result(fit, analysis = "run_discovery_analysis")
  } else {
    showf("fit must be a 'sceptre_object'", type = "stop")
  }
  
  out <- data.frame(
    coef = stats$grna_target,
    gene = stats$response_id,
    logFC = stats$log_2_fold_change,
    se = stats$se_fold_change,
    pvalue = stats$p_value,
    stringsAsFactors = FALSE
  )
  
  if (adjust_p) {
    out$FDR <- ave(
      out$pvalue,
      out$coef,
      FUN = function(z) { p.adjust(z, method = adjust_m) }
    )
  }
  
  # Validate
  rownames(out) <- NULL
  validate_dge(out)

  # Return
  out
}

extract_limma <- function(
  fit,
  intercept = FALSE,
  adjust_p = TRUE,
  adjust_m = "BH"
){

  coefs <- colnames(fit$coefficients)

  if (!intercept) {
    coefs <- setdiff(coefs, "(Intercept)")
  }

  out <- do.call(
    rbind,
    lapply(coefs, function(coef){

      idx <- match(coef, colnames(fit$coefficients))

      res <- data.frame(
        coef = coef,
        gene = rownames(fit$coefficients),
        logFC = fit$coefficients[, idx],
        se = fit$stdev.unscaled[, idx] * sqrt(fit$s2.post),
        pvalue = fit$p.value[, idx],
        stringsAsFactors = FALSE
      )

      res
    })
  )

  # Optionally adjust p-values
  if (adjust_p) {
    out$FDR <- ave(
      out$pvalue,
      out$coef,
      FUN = function(z){ p.adjust(z, method = adjust_m)}
    )
  }

  # Validate
  rownames(out) <- NULL
  validate_dge(out)

  # Return
  out
}

extract_edger <- function(
  fit,
  intercept = FALSE,
  adjust_p = TRUE,
  adjust_m = "BH"
){

  coefs <- colnames(fit$design)

  if (!intercept) {
    coefs <- setdiff(coefs, "(Intercept)")
  }

  out <- do.call(
    rbind,
    lapply(coefs, function(coef){

      qlf <- edgeR::glmQLFTest(fit, coef = coef)
      
      res_table <- qlf$table
      coef_idx  <- match(coef, colnames(fit$design))

      res <- data.frame(
        coef = coef,
        gene = rownames(qlf),
        logFC = res_table$logFC,
        se = NA_real_,
        pvalue = res_table$PValue,
        stringsAsFactors = FALSE
      )

      res
    })
  )

  # Optionally adjust p-values
  if (adjust_p) {
    out$FDR <- ave(
      out$pvalue,
      out$coef,
      FUN = function(z){ p.adjust(z, method = adjust_m)}
    )
  }

  # Validate
  rownames(out) <- NULL
  validate_dge(out)

  # Return
  out
}

extract_deseq2 <- function(
  fit,
  gene_qc = FALSE,
  intercept = FALSE,
  adjust_p = TRUE,
  adjust_m = "BH"
){
  coefs <- DESeq2::resultsNames(fit)

  if (!intercept) {
    coefs <- setdiff(coefs, "Intercept")
  }

  out <- do.call(
    rbind,
    lapply(coefs, function(coef){
      
      res_obj <- DESeq2::results(
        fit, 
        name = coef,
        cooksCutoff = gene_qc, 
        independentFiltering = gene_qc
      )

      # Rename to limma format
      clean_coef <- coef
      if (clean_coef == "Intercept") {
        clean_coef <- "(Intercept)"
      } else {
        clean_coef <- gsub("_vs_.*$", "", clean_coef)
        clean_coef <- gsub("_", "", clean_coef)
      }
      
      res <- data.frame(
        coef = clean_coef,
        gene = rownames(res_obj),
        logFC = res_obj$log2FoldChange,
        se = res_obj$lfcSE,
        pvalue = res_obj$pvalue,
        stringsAsFactors = FALSE
      )

      res
    })
  )

  # Optionally adjust p-values
  if (adjust_p) {
    out$FDR <- ave(
      out$pvalue,
      out$coef,
      FUN = function(z){ p.adjust(z, method = adjust_m)}
    )
  }

  # Validate
  rownames(out) <- NULL
  validate_dge(out)

  # Return
  out
}

extract_seurat <- function(
  fit,
  adjust_p = TRUE,
  adjust_m = "BH"
){

  out <- data.frame(
    coef = fit$coef,
    gene = rownames(fit),
    logFC = fit$avg_log2FC,
    se = NA_real_,
    pvalue = fit$p_val,
    stringsAsFactors = FALSE
  )

  # Optionally adjust p-values
  if (adjust_p) {
    out$FDR <- ave(
      out$pvalue,
      out$coef,
      FUN = function(z){ p.adjust(z, method = adjust_m)}
    )
  }

  # Validate
  rownames(out) <- NULL
  validate_dge(out)

  # Return
  out
}

extract_scmageck <- function(
  fit,
  intercept = FALSE,
  adjust_p = TRUE,
  adjust_m = "BH"
){

  score <- fit$score
  pval  <- fit$pvalue

  # coefficients (= perturbations)
  coefs <- score$Perturbedgene

  # Remove NegCtrl unless requested (analogous to removing Intercept)
  if (!intercept) {
    coefs <- setdiff(coefs, "NegCtrl")
  }

  # Response genes
  genes <- setdiff(colnames(score), "Perturbedgene")

  out <- do.call(
    rbind,
    lapply(coefs, function(coef){

      idx <- match(coef, score$Perturbedgene)

      res <- data.frame(
        coef = coef,
        gene = genes,
        logFC = as.numeric(score[idx, genes]),
        se = NA_real_,
        pvalue = as.numeric(pval[idx, genes]),
        stringsAsFactors = FALSE
      )

      res
    })
  )

  # Optionally adjust p-values
  if (adjust_p) {
    out$FDR <- ave(
      out$pvalue,
      out$coef,
      FUN = function(z){ p.adjust(z, method = adjust_m)}
    )
  }

  # Validate
  rownames(out) <- NULL
  validate_dge(out)

  # Return
  out
}

## Meta analysis ==============================================================
prepare_trade_input <- function(stats, guide) {

  # Rename columns
  stats <- stats[stats$coef == guide, c("gene", "logFC", "se", "pvalue")]
  stats <- stats[order(stats$gene), ]
  colnames(stats) <- c("gene", "log2FoldChange", "lfcSE", "pvalue")
  
  # Genes as rownames
  rownames(stats) <- stats$gene
  stats$gene <- NULL
  
  stats
}


.robust_corr_weights <- function(
    x,
    top_genes = 500
){
  stats <- x$stats

  # Create matrices for logFC and SE
  fc_mat <- xtabs(logFC ~ coef + gene, data = stats)
  se_mat <- xtabs(se ~ coef + gene, data = stats)
  
  # Calculate the weights matrix (inverse-variance)
  w_mat <- 1 / (se_mat^2)
  w_mat[is.infinite(w_mat)] <- 0 # Handle division by zero edge cases

  guides <- rownames(fc_mat)
  guide_info <- unique(stats[, c("target", "coef")])

  out <- vector("list", length(guides))
  
  for(i in seq_along(guides)){
    guide <- guides[i]
    target <- guide_info$target[match(guide, guide_info$coef)]
    
    # Correctly look for other guides belonging to the SAME target
    siblings <- guide_info$coef[guide_info$target == target & guide_info$coef != guide]

    weight_with_se    <- NA_real_
    weight_without_se <- NA_real_

    if(length(siblings) > 0){
      
      # With se (Inverse-Variance Weighted Consensus)
      sib_weights <- w_mat[siblings, , drop = FALSE]
      sib_logfc   <- fc_mat[siblings, , drop = FALSE]
      
      consensus_with_se <- colSums(sib_logfc * sib_weights, na.rm = TRUE) / 
                           colSums(sib_weights, na.rm = TRUE)
      consensus_with_se[is.nan(consensus_with_se)] <- 0
      
      # Top genes determined by error-aware consensus magnitude
      keep_with_se <- order(abs(consensus_with_se), decreasing = TRUE)
      keep_with_se <- keep_with_se[seq_len(min(top_genes, length(keep_with_se)))]

      # Weighted Spearman Ranks
      guide_ranks_w <- rank(fc_mat[guide, keep_with_se])
      cons_ranks_w  <- rank(consensus_with_se[keep_with_se])
      
      # Joint analytical weights (Guide SE precision * Mean Sibling SE precision)
      joint_weights <- w_mat[guide, keep_with_se] * 
                       colMeans(sib_weights[, keep_with_se, drop = FALSE], na.rm = TRUE)
      
      weight_with_se <- tryCatch({
        cov.wt(
          cbind(guide_ranks_w, cons_ranks_w), 
          wt = joint_weights, 
          cor = TRUE
        )$cor[1, 2] 
      }, error = function(e) NA_real_)
    
      # Without se (Standard Unweighted Consensus)
      consensus_without_se <- colMeans(sib_logfc, na.rm = TRUE)
      consensus_without_se[is.nan(consensus_without_se)] <- 0
      
      # Top genes determined by unweighted raw logFC consensus magnitude
      keep_without_se <- order(abs(consensus_without_se), decreasing = TRUE)
      keep_without_se <- keep_without_se[seq_len(min(top_genes, length(keep_without_se)))]

      # Standard Spearman Correlation (using standard cor())
      weight_without_se <- suppressWarnings(
        cor(
          fc_mat[guide, keep_without_se],
          consensus_without_se[keep_without_se],
          method = "spearman",
          use = "pairwise.complete.obs"
        )
      )
    }

    out[[i]] <- data.frame(
      target            = target,
      coef              = guide,
      weight_with_se    = weight_with_se,
      weight_without_se = weight_without_se
    )
  }

  x$weights <- do.call(rbind, out)
  
  return(x)
}



.agg_stouffer <- function(
    x,
    weights = NULL,
    signed = TRUE
){

  if (is.null(weights)) {
    weights <- rep(1, nrow(x))
  }

  if(length(weights) != nrow(x)){
    stop("weights must have length nrow(x)")
  }

  p <- pmin(
    pmax(
      x$pvalue,
      .Machine$double.xmin
    ),
    1 - 1e-16
  )

  if(signed){

    zscore <- qnorm(
      1 - p / 2
    ) * sign(x$logFC)

  } else {

    zscore <- qnorm(
      1 - p
    )
  }

  x$..zscore <- zscore
  x$..weight <- weights

  res <- lapply(
    split(
      x,
      list(
        x$target,
        x$gene
      ),
      drop = TRUE
    ),
    function(d){

      z <- d$..zscore
      w <- d$..weight

      keep <- is.finite(z) & is.finite(w)

      z <- z[keep]
      w <- w[keep]

      if(length(z) == 0 || sum(w^2) == 0){

        return(
          data.frame(
            target = d$target[1],
            gene = d$gene[1],
            stat = NA_real_,
            pvalue = NA_real_
          )
        )
      }

      Z <- sum(w * z) / sqrt(sum(w^2))

      pval <- if(signed){
        2 * pnorm(-abs(Z))
      } else {
        pnorm(-Z)
      }

      data.frame(
        target = d$target[1],
        gene = d$gene[1],
        stat = Z,
        pvalue = pval
      )
    }
  )

  do.call(rbind, res)
}

.aggregate_pvals <- function(
    x,
    method = c(
      "stouffer",
      "fisher",
      "cauchy",
      "hmp",
      "minp"
    ),
    weights = NULL,
    signed = TRUE
){

  # Method selection
  method <- match.arg(method)

  # Extract input obj.
  stats <- x$stats

  # Required columns
  req_cols <- c("target", "gene", "pvalue")
  if (signed) { req_cols <- c(req_cols, "logFC") }

  miss_cols <- setdiff(req_cols, names(stats))
  if(length(miss_cols) > 0){
    stop("Missing required columns: ", paste(miss_cols, collapse = ", "))
  }

  # Aggregation
  out <- switch(
    method,

    stouffer = .agg_stouffer(
      stats,
      weights = weights,
      signed = signed
    ),

    fisher = .agg_fisher(stats),
    cauchy = .agg_cauchy(stats),
    hmp = .agg_hmp(stats),
    minp = .agg_minp(stats)
  )

  out$FDR <- p.adjust(
    out$pvalue,
    method = "BH"
  )

  out$method <- method
  rownames(out) <- NULL
  x$agg <- out

  x
}

## ggplot2 ====================================================================
.plot_pval_dist <- function(
  data,
  x_var,
  color_var = NULL,
  fill_var = NULL,
  facet_col = NULL,
  facet_levels = NULL,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  color_name = NULL,
  fill_name = NULL,
  caption = NULL,
  bins = 50,
  fill_bins = 15,
  fill_limits = NULL   
) {

  # Convert to dataframe
  df <- as.data.frame(data)

  # Validate x_var
  if (!(x_var %in% names(df))) {
    stop("Data must contain column: ", x_var)
  }

  # Faceting
  if (!is.null(facet_col)) {
    if (!(facet_col %in% names(df))){
      stop("Data must contain column: ", facet_col)
    }
    if (is.null(facet_levels)) {
      facet_levels <- head(unique(df[[facet_col]]), 9)
    }
    df <- df[df[[facet_col]] %in% facet_levels, , drop = FALSE]
    df[[facet_col]] <- factor(df[[facet_col]], levels = facet_levels)
  }

  # Fill
  fill_colors <- NULL
  if (!is.null(fill_var)) {
    if (!(fill_var %in% names(df))) {
      stop("Data must contain column: ", fill_var)
    }

    # Determine range for fill scaling
    t_vals <- df[[fill_var]]

    if (is.null(fill_limits)) {
      max_abs <- max(abs(t_vals), na.rm = TRUE)
      if (max_abs < 1e-6) max_abs <- 1
      neg_limit <- -max_abs
      pos_limit <- max_abs
    } else {
      neg_limit <- fill_limits[1]
      pos_limit <- fill_limits[2]
    }

    # Ensure fill_bins is odd
    if (fill_bins %% 2 == 0) {
      stop("fill_bins must be odd to ensure a center bin at 0.")
    }

    # Define bin centers & breaks
    bin_centers <- seq(neg_limit, pos_limit, length.out = fill_bins)
    bin_width  <- (pos_limit - neg_limit) / (fill_bins - 1)
    bin_breaks <- c(-Inf, head(bin_centers, -1) + bin_width / 2, Inf)
    bin_labels <- format(round(bin_centers, 1), nsmall = 1)
    bin_labels[1] <- paste0("< ", round(neg_limit, 1))
    bin_labels[length(bin_labels)] <- paste0("> ", round(pos_limit, 1))

    df$fill_bin <- cut(
      t_vals,
      breaks = bin_breaks,
      labels = bin_labels,
      include.lowest = TRUE
    )

    # Number of bins below/above zero for color gradients
    zero_idx <- which.min(abs(bin_centers))
    n_below <- zero_idx - 1
    n_above <- length(bin_centers) - zero_idx

    fill_colors <- c(
      colorRampPalette(c("#002082", "white"))(n_below + 1),
      colorRampPalette(c("white", "#a50026"))(n_above + 1)[-1]
    )
    names(fill_colors) <- bin_labels
  } 

  # Color
  if (!is.null(color_var)) {
    if (!(color_var %in% names(df))) {
      stop("Data must contain column: ", color_var)
    }

    df[[color_var]] <- as.factor(df[[color_var]])
  }

  # Plot
  p <- ggplot2::ggplot(
    data = df, 
    mapping = ggplot2::aes(
      x = .data[[x_var]]
    )
  )

 # Histogram
 # Conditional Geoms & Aesthetics Injection
  if (!is.null(fill_var) && !is.null(color_var)) {
    p <- p + ggplot2::geom_histogram(
      mapping = ggplot2::aes(
        fill = fill_bin, 
        color = .data[[color_var]]
      ),
      bins = bins,
      position = "stack",
      linewidth = 0.2
    )

  } else if (!is.null(fill_var)) {
    p <- p + ggplot2::geom_histogram(
      mapping = ggplot2::aes(
        fill = fill_bin
      ),
      bins = bins,
      position = "stack",
      color = "black",
      linewidth = 0.2
    )

  } else if (!is.null(color_var)) {
    p <- p + ggplot2::geom_histogram(
      mapping = ggplot2::aes(
        color = .data[[color_var]]
      ),
      bins = bins,
      position = "stack",
      fill = "grey80",
      linewidth = 0.2
    )

  } else {
    p <- p + ggplot2::geom_histogram(
      bins = bins,
      position = "stack",
      fill = "grey80",
      color = "black",
      linewidth = 0.2
    )
  }

  # Fill scale
  if (!is.null(fill_var)) {
    p <- p + ggplot2::scale_fill_manual(
      name = fill_name,
      values = fill_colors,
      guide = ggplot2::guide_legend(
        override.aes = list(
          color = "black", 
          linewidth = 0.1
        )
      )
    )
  }

  # Color scale
  if (!is.null(color_var)) {
    p <- p + ggplot2::scale_color_brewer(
      name = color_name,
      palette = "Set1",
      guide = ggplot2::guide_legend(
        override.aes = list(
          fill = "white", 
          linewidth = 0.1
        )
      )
    )
  }

  # Facets
  if (do_facet) {
    p <- p + ggplot2::facet_wrap(
        as.formula(paste("~", facet_col)
      )
    )
  }

  # Scales
  p <- p + ggplot2::scale_y_continuous(
    breaks = function(x) seq(
      floor(min(x)),
      ceiling(max(x)),
      by = 1
    )
  )

  # Labels
  p <- p + ggplot2::labs(
    title = title,
    x = ifelse(is.null(x_label), x_var, x_label),
    y = ifelse(is.null(y_label), "count", y_label),
    caption = caption
  )

  # Return
  p
}

theme_DK <- function(
  base_size = 5,
  lgd_key_prop = 1,
  lgd_position = "right",
  make_void = FALSE,
  show_axis = TRUE,
  show_grid = TRUE,
  show_facet_x = TRUE,
  show_facet_y = TRUE,
  show_label_x = TRUE,
  show_label_y = TRUE,
  show_ticks_x = TRUE,
  show_ticks_y = TRUE,
  show_borders = TRUE,
  zero_margins = TRUE,
  show_identity = FALSE,
  xintercept = NULL,
  yintercept = NULL,
  rotate = NULL,
  base_family = "sans",
  ...
) {

  # Unit conversion helpers 
  mm_to_pt <- function(mm) mm / 0.352778
  pt_to_mm <- function(pt) pt * 0.352778

  # Set internal params. for void
  if (make_void) {
    show_axis <- FALSE
    show_grid <- FALSE
    show_label_x <- FALSE
    show_label_y <- FALSE
    show_ticks_x <- FALSE
    show_ticks_y <- FALSE
    show_borders <- FALSE
  }

  # Base theme
  p <- ggplot2::theme(
    plot.background  = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    strip.background = ggplot2::element_rect(fill = "white", color = NA),

    # Facets
    strip.placement = "inside",
    strip.text.x = ggplot2::element_text(size = base_size, family = base_family, margin = ggplot2::margin(b = 1, t = 0, unit = "mm")),
    strip.text.y = ggplot2::element_text(size = base_size, family = base_family, margin = ggplot2::margin(b = 1, t = 0, unit = "mm")),
    panel.spacing.x = grid::unit(1, "mm"),
    panel.spacing.y = grid::unit(1, "mm"),

    axis.line  = ggplot2::element_line(color = "black", linewidth = 0.3),
    axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = ggplot2::unit(1, "mm"),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = if (show_grid) ggplot2::element_line(color = "grey80", linewidth = 0.3) else ggplot2::element_blank(),

    # Text
    axis.text  = ggplot2::element_text(size = base_size, family = base_family),
    axis.title = ggplot2::element_text(size = base_size, family = base_family),
    plot.title = ggplot2::element_text(size = base_size, hjust = 0.5, family = base_family),
    plot.subtitle = ggplot2::element_text(size = base_size, hjust = 0.5, family = base_family),
    plot.caption = ggplot2::element_text(size = base_size, hjust = 0, family = base_family),
    ...
  )

  # Legend sizing
  lgd_key_size  <- base_size * lgd_key_prop

  p <- p + ggplot2::theme(
    legend.title = ggplot2::element_text(size = base_size, family = base_family, margin = ggplot2::margin(b = 1, unit = "mm")),
    legend.text = ggplot2::element_text(size = base_size, family = base_family),

    legend.key.height = ggplot2::unit(pt_to_mm(lgd_key_size), "mm"),
    legend.key.width = ggplot2::unit(pt_to_mm(base_size), "mm"),
    legend.key.spacing = ggplot2::unit(0.5, "mm"),
    legend.box.spacing = ggplot2::unit(1, "mm"),
    legend.background = ggplot2::element_rect(fill = NA, color = NA)
  )

  # Legend position
  theme_lgd_pos    <- lgd_position
  theme_lgd_inside <- NULL
  theme_lgd_just   <- NULL

  if (lgd_position %in% c("inside.top.left", "inside.top.right", "inside.bottom.left", "inside.bottom.right")) {
    theme_lgd_pos <- "inside"
    
    if (lgd_position == "inside.top.left") {
      theme_lgd_inside <- c(0.02, 0.98)
      theme_lgd_just   <- c(0, 1)
    } else if (lgd_position == "inside.top.right") {
      theme_lgd_inside <- c(0.98, 0.98)
      theme_lgd_just   <- c(1, 1)
    } else if (lgd_position == "inside.bottom.left") {
      theme_lgd_inside <- c(0.02, 0.02)
      theme_lgd_just   <- c(0, 0)
    } else if (lgd_position == "inside.bottom.right") {
      theme_lgd_inside <- c(0.98, 0.02)
      theme_lgd_just   <- c(1, 0)
    }
  }

  p <- p + ggplot2::theme(
    legend.position = theme_lgd_pos,
    legend.position.inside = theme_lgd_inside,
    legend.justification = theme_lgd_just
  )

  # Adaptive margins (zero margins)
  if (zero_margins) {
    p <- p + ggplot2::theme(
      legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
      plot.margin = ggplot2::margin(1, 1, 1, 1, unit = "mm")
    )
  }

  # Adaptive facets (x and y facet)
  if (!show_facet_x) p <- p + ggplot2::theme(strip.text.x = ggplot2::element_blank())
  if (!show_facet_y) p <- p + ggplot2::theme(strip.text.y = ggplot2::element_blank())

  # Adaptive boarders
  if (show_borders) { p <- p + ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.line = ggplot2::element_blank()
    )
  }

  # Adaptive axis
  if (!show_axis) {
    p <- p + ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(), 
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 0.5, b = 0, unit = "mm")), # pull x labels closer
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 0.5, l = 0, unit = "mm"))  # pull y labels closer
    )
  }

  # Adaptive rotation
  if (!is.null(rotate)) {
    angle <- rotate %% 360
    # Determine justification heuristically
    if (angle == 0) {
      hjust <- 0.5; vjust <- 0.5
    } else if (angle == 90) {
      hjust <- 1; vjust <- 0.5
    } else if (angle == 270) {
      hjust <- 0; vjust <- 0.5
    } else if (angle > 0 && angle < 90) {
      hjust <- 1; vjust <- 1
    } else if (angle > 90 && angle < 180) {
      hjust <- 1; vjust <- 0
    } else if (angle > 180 && angle < 270) {
      hjust <- 0; vjust <- 0
    } else { # e.g., 315°
      hjust <- 0; vjust <- 1
    }

    p <- p + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = angle,
        hjust = hjust,
        vjust = vjust,
        family = base_family
      )
    )
  }

  if (!show_label_x) {
    p <- p + ggplot2::theme(
      axis.title.x = ggplot2::element_blank()
    )
  }

  if (!show_label_y) {
    p <- p + ggplot2::theme(
      axis.title.y = ggplot2::element_blank()
    )
  }

  if (!show_ticks_x) {
    p <- p + ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank()
    )
  }

  if (!show_ticks_y) {
    p <- p + ggplot2::theme(
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    )
  }

  # Extra layers
  layers <- list()

  # Identity line y = x
  if (show_identity) {
    layers <- append(
      layers,
      list(
        ggplot2::geom_abline(
          slope = 1,
          intercept = 0,
          linetype = "solid",
          linewidth = 0.3,
          color = "black"
        )
      )
    )
  }

   # Vertical intercept(s)
  if (!is.null(xintercept)) {
    layers <- append(
      layers,
      list(
        ggplot2::geom_vline(
          xintercept = xintercept,
          linetype = "solid",
          linewidth = 0.3,
          color = "black"
        )
      )
    )
  }

  # Horizontal intercept(s)
  if (!is.null(yintercept)) {
    layers <- append(
      layers,
      list(
        ggplot2::geom_hline(
          yintercept = yintercept,
          linetype = "solid",
          linewidth = 0.3,
          color = "black"
        )
      )
    )
  }

  # Return theme + layers
  return(c(list(p), layers))
}

add_corner_axis <- function(p, base_size = 5) {

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  xlab <- p$labels$x %||% "x"
  ylab <- p$labels$y %||% "y"

  # Apply coord_fixed
  build <- ggplot2::ggplot_build(p)

  xr <- build$layout$panel_params[[1]]$x.range
  yr <- build$layout$panel_params[[1]]$y.range

  dx <- diff(xr)
  dy <- diff(yr)

  x0 <- xr[1]
  y0 <- yr[1]

  # proportional arrow lengths
  arrow_x <- 0.1 * dx
  arrow_y <- 0.1 * dy

  # proportional label offsets
  offset_x <- 0.03 * dy
  offset_y <- 0.03 * dx

  p <- p + ggplot2::annotate(
    "segment",
    x = x0, 
    xend = x0 + arrow_x,
    y = y0, 
    yend = y0,
    arrow = grid::arrow(
      length = grid::unit(1, "mm"),
      type = "closed"
    ),
    linewidth = 0.3
  ) +
  ggplot2::annotate(
    "segment",
    x = x0, 
    xend = x0,
    y = y0, 
    yend = y0 + arrow_y,
    arrow = grid::arrow(
      length = grid::unit(1, "mm"),
      type = "closed"
    ),
    linewidth = 0.3
  ) +
  ggplot2::annotate(
    "text",
    x = x0 + arrow_x / 2,
    y = y0 - offset_x,
    label = xlab,
    hjust = 0.5,
    vjust = 1,
    size = base_size / 3
  ) +
  ggplot2::annotate(
    "text",
    x = x0 - offset_y,
    y = y0 + arrow_y / 2,
    label = ylab,
    angle = 90,
    hjust = 0.5,
    vjust = 0,
    size = base_size / 3
  )

  return(p)
}

## ComplexHeatmap =============================================================
cmplx_gp <- function(
  base_size = 5, 
  base_family = "sans", 
  ...
) {
  grid::gpar(
    fontsize = base_size,
    fontfamily = base_family,
    fontface = "plain",
    ...
  )
}

cmplx_bar <- function(
  x,
  base_size = 5,
  lgd_key_prop = 1,
  base_family = "sans",
  which = c("row", "column"),
  axis_side = NULL,
  flip = FALSE,
  ...
) {

  which <- match.arg(which)
  
  if (is.null(axis_side)) {
    axis_side <- if (which == "row") "top" else "left"
  } else {
    if (which == "row") {
      axis_side <- match.arg(axis_side, c("top", "bottom"))
    } else {
      axis_side <- match.arg(axis_side, c("left", "right"))
    }
  }

  gp_base <- cmplx_gp(base_size, base_family, lwd = 0.5)
  anno_height <- grid::unit(base_size * 2, "mm")

  max_val <- max(x, na.rm = TRUE)
  tick_positions <- seq(from = 0, to = max_val, length.out = 3) 

  tick_labels <- scales::label_number(
    scale_cut = scales::cut_short_scale(),
    accuracy = 1
  )(tick_positions)

  fun <- ComplexHeatmap::anno_barplot(
    x,
    which = which,
    border = FALSE,
    axis = TRUE,
    gp = grid::gpar(
      fill = "grey40", 
      col = NA
    ),
    axis_param = list(
      gp = gp_base,
      labels_rot = 90,
      side = axis_side,
      at = tick_positions,
      labels = tick_labels,
      direction = if (flip) "reverse" else "normal" 
    ),
    width = if (which == "row") anno_height else NULL,
    height = if (which == "column") anno_height else NULL
  )

  # Return
  fun
}

cmplx_anno <- function(
  ...,
  base_size = 5,
  lgd_key_prop = 1,
  base_family = "sans",
  which = c("column", "row")
) {

  which <- match.arg(which)

  dots <- list(...)
  name_side <- dots$annotation_name_side

  # Only validate if provided
  if (!is.null(name_side)) {
    
    if (which == "column" && !name_side %in% c("left", "right")) {
      stop("For column annotations, annotation_name_side must be 'left' or 'right'")
    }

    if (which == "row" && !name_side %in% c("top", "bottom")) {
      stop("For row annotations, annotation_name_side must be 'top' or 'bottom'")
    }
  }

  # Unit conversion helpers 
  mm_to_pt <- function(mm) mm / 0.352778
  pt_to_mm <- function(pt) pt * 0.352778

  # Annotation scaling
  gp_base <- cmplx_gp(base_size, base_family)
  anno_height <- grid::unit(pt_to_mm(base_size), "mm")

  # Legend scaling
  lgd_key_size <- base_size * lgd_key_prop
  lgd_unit <- grid::unit(pt_to_mm(lgd_key_size), "mm")

  anno <- ComplexHeatmap::HeatmapAnnotation(
    ...,

    which = which,
    gap = grid::unit(0.7, "mm"),
    simple_anno_size = anno_height,
    annotation_name_gp = gp_base,

    # Legend styling
    annotation_legend_param = list(
      title_position = "topleft",
      title_gp = gp_base,
      labels_gp = gp_base,
      gap = grid::unit(0.5, "mm"),
      grid_width = lgd_unit,
      grid_height = lgd_unit,
      legend_height = lgd_unit
    )
  )

  # Return
  anno
}

cmplx_ht <- function(
  mat,
  ...,
  base_size = 5,
  lgd_key_prop = 1,
  base_family = "sans",
  cell_size = NULL
) {

  stopifnot(is.matrix(mat) || is.data.frame(mat))
  mat <- as.matrix(mat)

  # Dimensions
  n_row <- nrow(mat)
  n_col <- ncol(mat)

  # Unit conversion helpers 
  mm_to_pt <- function(mm) mm / 0.352778
  pt_to_mm <- function(pt) pt * 0.352778

  # Annotation scaling
  gp_base <- cmplx_gp(base_size, base_family)
  anno_height <- grid::unit(pt_to_mm(base_size), "mm")

  # Legend scaling
  lgd_key_size <- base_size * lgd_key_prop
  lgd_unit <- grid::unit(pt_to_mm(lgd_key_size), "mm")

  # Conditionally enforce square tiles
  size_args <- list()
  if (!is.null(cell_size)) {
    if (is.numeric(cell_size)) {
      cell_size <- grid::unit(cell_size, "mm")
    }
    size_args$width  <- cell_size * n_col
    size_args$height <- cell_size * n_row
  }

  ht <- do.call(ComplexHeatmap::Heatmap, c(list(
    mat,
    ...,

    # Heatmap text
    row_names_gp = gp_base,
    column_names_gp = gp_base,
    column_title_gp = gp_base,
    row_title_gp = gp_base,

    # Heatmap legend
    heatmap_legend_param = list(
      title_gp = gp_base,
      labels_gp = gp_base,
      grid_width = lgd_unit,
      legend_height = lgd_unit
    ),

    # Dendrogram styling
    row_dend_gp = grid::gpar(lwd = 0.3),
    column_dend_gp = grid::gpar(lwd = 0.3)
  ), size_args))

  # Return
  return(ht)
}

cmplx_save <- function(
  plot,
  file,
  width,
  height,
  units = "in",
  device = NULL,
  res = 300,
  ...
) {

  to_inches <- function(x, units) {
    switch(
      units,
      "in" = x,
      "cm" = x / 2.54,
      "mm" = x / 25.4,
      "px" = x / res,
      stop("Unsupported units: ", units)
    )
  }

  width_in  <- to_inches(width, units)
  height_in <- to_inches(height, units)

  if (is.null(device)) {
    ext <- tools::file_ext(file)
    device <- switch(
      tolower(ext),
      pdf = grDevices::pdf,
      svg = grDevices::svg,
      png = grDevices::png,
      jpeg = grDevices::jpeg,
      jpg = grDevices::jpeg,
      tiff = grDevices::tiff,
      stop("Unsupported file type: ", ext)
    )
  }

  # raster vs vector handling
  if (identical(device, grDevices::png) ||
      identical(device, grDevices::jpeg) ||
      identical(device, grDevices::tiff)) {

    device(
      file,
      width = width_in,
      height = height_in,
      units = "in",
      res  = res,
      ...
    )

  } else {

    device(
      file,
      width = width_in,
      height = height_in,
      ...
    )
  }

  ComplexHeatmap::ht_opt(
    legend_gap = grid::unit(1, "mm")
  )

  grid::grid.newpage()
  ComplexHeatmap::draw(
    plot,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    legend_grouping = "original",
    align_heatmap_legend = "heatmap_top",
    align_annotation_legend = "heatmap_top",
    merge_legend = TRUE
  )

  invisible(grDevices::dev.off())
}
