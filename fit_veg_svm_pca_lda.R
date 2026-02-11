# fit_veg_svm_pca_lda.R
# Lightweight copy of the PCA-LDA weight training + SVM evaluation pipeline
# - DOES NOT modify the original `fit_veg_mixture_mesma.R`
# - Does NOT use spline-based weighting or MESMA library search
# - Uses an OOB tuning/validation split to evaluate PCA-LDA thresholding and final SVM
# - Saves `svm_pca_lda_model.rds` and `svm_threshold_results.csv` on success

# Usage: Rscript fit_veg_svm_pca_lda.R [input_csv] [oob_fraction] [cost]

library(dplyr)
library(magrittr)
library(e1071)  # SVM
library(MASS)    # lda
library(lubridate)
library(ggplot2) # plotting (prototype plots)
# Respect MESMA_SEED for reproducibility (default 42)
MESMA_SEED <- as.integer(Sys.getenv("MESMA_SEED", unset = "42"))
if (!is.finite(MESMA_SEED)) MESMA_SEED <- 42L
set.seed(MESMA_SEED)
# --- CONFIG ---
INPUT_CSV <- "C:/Users/yolan/Downloads/LS_S2_Harmonized_Timeseries_training.csv"
INFERENCE_CSV <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/landsat_lower_inference.csv"
OUTPUT_DIR <- "C:/MAP/svm_results"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
OOB_TUNING_FRACTION <- ifelse(length(commandArgs(trailingOnly = TRUE)) >= 2, as.numeric(commandArgs(trailingOnly = TRUE)[2]), 0.10)
SVM_COST <- ifelse(length(commandArgs(trailingOnly = TRUE)) >= 3, as.numeric(commandArgs(trailingOnly = TRUE)[3]), 1)
PCA_LDA_THRESHOLD_CANDIDATES <- c(0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)
TEMPORAL_AGGREGATION_DAYS <- 10L
TEMPORAL_BUDGET <- ceiling(365 / TEMPORAL_AGGREGATION_DAYS)
PRUNE_ZERO_MIN_FEATURES <- 1
ENABLE_LDA_L2_NORMALIZATION <- TRUE
RAW_BANDS <- c("blue","green","red","nir","swir1","swir2")
SOIL_LINE_SLOPE <- 1.0

cat(sprintf("[SVM PIPELINE] Input CSV: %s\n", INPUT_CSV))
cat(sprintf("[SVM PIPELINE] OOB fraction: %.2f, SVM cost: %s\n", OOB_TUNING_FRACTION, SVM_COST))

# --- Minimal helper functions copied/derived from main fit script ---
pheno_doy <- function(d) {
  d <- tryCatch(as.Date(d), error = function(e) NA)
  month <- lubridate::month(d)
  ifelse(is.na(d), NA_integer_,
         ifelse(month >= 3,
                as.integer(d - as.Date(paste0(lubridate::year(d), "-03-01"))) + 1L,
                as.integer(d - as.Date(paste0(lubridate::year(d) - 1, "-03-01"))) + 1L))
}

# Assign phenological year (March 1 = start of phenological year)
assign_pheno_year <- function(d) {
  d <- tryCatch(as.Date(d), error = function(e) NA)
  ifelse(is.na(d), NA_integer_, ifelse(lubridate::month(d) >= 3, lubridate::year(d), lubridate::year(d) - 1))
}

doy_to_pentad <- function(doy) {
  pmin(ceiling(doy / TEMPORAL_AGGREGATION_DAYS), TEMPORAL_BUDGET)
}

build_pentad_matrix <- function(dly_year, avail_idx) {
  if (is.null(dly_year) || nrow(dly_year) == 0) return(NULL)
  if (!"doy" %in% names(dly_year) || any(is.na(dly_year$doy))) {
    dly_year$doy <- pheno_doy(dly_year$date)
  }
  dly_year$pentad <- doy_to_pentad(dly_year$doy)
  K <- length(avail_idx)
  pentad_mat <- matrix(NA_real_, nrow = TEMPORAL_BUDGET, ncol = K)
  colnames(pentad_mat) <- avail_idx
  for (p in 1:TEMPORAL_BUDGET) {
    subset_p <- dly_year[dly_year$pentad == p, ]
    if (nrow(subset_p) == 0) next
    t_start <- (p - 1) * TEMPORAL_AGGREGATION_DAYS + 1
    t_end <- min(p * TEMPORAL_AGGREGATION_DAYS, TEMPORAL_BUDGET * TEMPORAL_AGGREGATION_DAYS)
    t_center <- (t_start + t_end) / 2
    for (j in seq_along(avail_idx)) {
      idx <- avail_idx[j]
      if (!idx %in% names(subset_p)) next
      v <- subset_p[[idx]]
      v <- v[is.finite(v)]
      if (length(v) == 0) next
      doys <- subset_p$doy
      doys <- doys[is.finite(doys)]
      if (length(doys) != length(v)) {
        n_use <- min(length(doys), length(v))
        doys <- doys[seq_len(n_use)]
        v <- v[seq_len(n_use)]
      }
      if (length(v) == 1 || length(unique(doys)) < 2) {
        pentad_mat[p, j] <- mean(v, na.rm = TRUE)
      } else {
        x <- doys - t_center
        model <- tryCatch(stats::lm(v ~ x), error = function(e) NULL)
        if (!is.null(model)) {
          b0 <- tryCatch(stats::coef(model)[[1]], error = function(e) NA_real_)
          if (is.finite(b0)) pentad_mat[p, j] <- as.numeric(b0) else pentad_mat[p, j] <- mean(v, na.rm = TRUE)
        } else {
          pentad_mat[p, j] <- mean(v, na.rm = TRUE)
        }
      }
    }
  }
  for (j in 1:K) {
    vals <- pentad_mat[, j]
    if (any(is.na(vals))) {
      if (all(is.na(vals))) {
        pentad_mat[, j] <- 0
      } else {
        idx_present <- which(!is.na(vals))
        if (length(idx_present) >= 2) pentad_mat[, j] <- approx(idx_present, vals[idx_present], xout = 1:TEMPORAL_BUDGET, rule = 2)$y else pentad_mat[, j] <- vals[idx_present[1]]
      }
    }
  }
  pentad_mat
}

# Minimal normalization function (z-score) for indices/bands
normalize_mesma_data <- function(df, cols) {
  INDEX_SCALES <- list()
  present <- intersect(cols, names(df))
  for (col in present) {
    vals <- df[[col]]
    vals <- vals[is.finite(vals)]
    if (length(vals) > 0) {
      mu <- mean(vals)
      sigma <- sd(vals)
      if (!is.finite(sigma) || sigma < 1e-10) sigma <- 1.0
      INDEX_SCALES[[col]] <- list(mean = mu, sd = sigma)
    }
  }
  for (col in names(INDEX_SCALES)) {
    if (col %in% names(df)) {
      params <- INDEX_SCALES[[col]]
      mu <- params$mean; sigma <- params$sd
      if (is.finite(sigma) && sigma > 0) df[[col]] <- (df[[col]] - mu) / sigma
    }
  }
  list(df = df, INDEX_SCALES = INDEX_SCALES)
}

safe_lda_call <- function(X_pca, y, min_n_pcs = 2) {
  if (is.null(X_pca) || ncol(X_pca) < min_n_pcs) {
    cat(sprintf("safe_lda_call: Not enough PCs (have=%d, min=%d).\n", ncol(X_pca), min_n_pcs))
    return(NULL)
  }
  curr_n_pcs <- ncol(X_pca)
  for (attempt in seq_len(max(1, curr_n_pcs - min_n_pcs + 1))) {
    tryCatch({
      res <- MASS::lda(X_pca[, 1:curr_n_pcs, drop = FALSE], grouping = y)
      return(res)
    }, warning = function(w) {
      warn_msg <- conditionMessage(w)
      if (grepl("collinear", warn_msg, ignore.case = TRUE)) {
        curr_n_pcs <- curr_n_pcs - 1
      } else return(res)
    }, error = function(e) return(NULL))
  }
  cat("safe_lda_call: Exhausted retries; LDA could not be computed without collinearity.\n")
  NULL
}

train_feature_pipeline <- function(df, class_col, feature_cols, use_lda = TRUE) {
  cat(sprintf("\n=== Training Feature Pipeline for Class: %s ===\n", class_col))
  X_raw <- list(); y_labels <- c()
  class_values <- unique(na.omit(df[[class_col]]))
  traces_by_class <- lapply(class_values, function(cv) {
    df_class <- df[df[[class_col]] == cv, , drop = FALSE]
    if (nrow(df_class) == 0) return(list())
    split(df_class, list(df_class$location_id, df_class$pheno_year), drop = TRUE)
  })
  traces <- unlist(traces_by_class, recursive = FALSE)
  cat("  Building trace matrix...\n")
  for(sub in traces) {
    if(nrow(sub) < 5) next
    mat <- build_pentad_matrix(sub, feature_cols)
    if(is.null(mat)) next
    vec <- as.numeric(mat); vec[!is.finite(vec)] <- NA
    X_raw[[length(X_raw)+1]] <- vec
    lbl <- names(sort(table(sub[[class_col]]), decreasing=TRUE))[1]
    y_labels <- c(y_labels, lbl)
  }
  if (length(X_raw) < 10) return(NULL)
  X_mat <- do.call(rbind, X_raw)
  if (isTRUE(ENABLE_LDA_L2_NORMALIZATION)) {
    cat("  L2-normalizing training samples (per-observation) for LDA...\n")
    X_mat <- t(apply(X_mat, 1, function(r) { r_clean <- r; r_clean[is.na(r_clean)] <- 0; nrm <- sqrt(sum(r_clean^2)); if (!is.finite(nrm) || nrm < 1e-9) return(r); r / nrm }))
  } else { cat("  L2-normalization for LDA skipped (using raw values)...\n") }
  n_bins <- TEMPORAL_BUDGET
  n_idx <- length(feature_cols)
  global_means <- numeric(n_idx); global_sds <- numeric(n_idx); names(global_means) <- feature_cols; names(global_sds) <- feature_cols
  X_z <- X_mat
  cat("  Computing Z-score parameters...\n")
  for(k in seq_along(feature_cols)) {
    col_idx_start <- (k-1)*n_bins + 1; col_idx_end <- k*n_bins
    vals <- X_mat[, col_idx_start:col_idx_end]
    mu <- mean(vals, na.rm=TRUE); sigma <- sd(vals, na.rm=TRUE)
    if(sigma == 0 || is.na(sigma)) sigma <- 1
    global_means[k] <- mu; global_sds[k] <- sigma
    X_z[, col_idx_start:col_idx_end] <- (vals - mu) / sigma
  }
  X_z[!is.finite(X_z)] <- 0
  cat("  Computing PCA (optionally LDA) weights...\n")
  vars <- apply(X_z, 2, var); keep_cols <- vars > 1e-9
  X_pca_in <- X_z[, keep_cols, drop=FALSE]
  pca_res <- prcomp(X_pca_in, center = FALSE, scale. = FALSE)
  cum_var <- cumsum(pca_res$sdev^2) / sum(pca_res$sdev^2); n_pcs <- which(cum_var > 0.95)[1]
  if(is.na(n_pcs)) n_pcs <- ncol(pca_res$x)

  # PCA-only weighting (no LDA)
  if (exists("use_lda") && identical(use_lda, FALSE)) {
    # Importance per PC (variance explained)
    prop <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
    # Use first n_pcs
    prop_n <- prop[1:n_pcs]
    R <- pca_res$rotation[, 1:n_pcs, drop = FALSE]
    # Feature importance = abs(R) * PC importance summed across PCs
    weights_clean <- as.numeric(abs(R) %*% prop_n)
    final_weights <- numeric(ncol(X_z)); final_weights[keep_cols] <- weights_clean
    return(list(means = global_means, sds = global_sds, weights = final_weights, indices = feature_cols))
  }

  # Default: PCA -> LDA pipeline (legacy behaviour)
  class_counts <- table(y_labels); n_min <- if (length(class_counts) > 0) min(class_counts) else 0; n_classes <- length(unique(y_labels))
  min_pcs_for_lda <- max(1, n_classes - 1)
  max_pcs_for_lda <- min(20, max(min_pcs_for_lda, n_min - 2))
  if (n_pcs > max_pcs_for_lda) { n_pcs <- max_pcs_for_lda }
  if (n_pcs < min_pcs_for_lda) { warning("Not enough degrees of freedom for LDA"); return(NULL) }

  lda_res <- safe_lda_call(pca_res$x[, 1:n_pcs, drop=FALSE], as.factor(y_labels), min_n_pcs = min_pcs_for_lda)
  if (is.null(lda_res)) { cat("[FALLBACK] LDA failed; using uniform weights.\n"); final_weights <- rep(1, ncol(X_z)); return(list(means = global_means, sds = global_sds, weights = final_weights, indices = feature_cols)) }

  W_pc <- lda_res$scaling; R <- pca_res$rotation[, 1:n_pcs, drop=FALSE]; W_std <- R %*% W_pc
  svd_vals <- lda_res$svd; prop <- svd_vals / sum(svd_vals)
  if (ncol(W_std) > 1) { n_dim <- min(length(prop), ncol(W_std)); weights_clean <- rowSums(abs(W_std[, 1:n_dim, drop=FALSE]) %*% diag(prop[1:n_dim], nrow=n_dim)) } else { weights_clean <- abs(W_std[, 1]) }
  final_weights <- numeric(ncol(X_z)); final_weights[keep_cols] <- weights_clean
  return(list(means = global_means, sds = global_sds, weights = final_weights, indices = feature_cols))
}

# Helper: build matrix (rows = traces) given weights mask - used for SVM training/eval
build_matrix_from_df <- function(df_in, indices, params, weights) {
  traces <- unique(df_in[, c("location_id", "pheno_year", "Veg")])
  X_list <- list(); y_vec <- c()
  for (j in seq_len(nrow(traces))) {
    lid <- traces$location_id[j]; pyr <- traces$pheno_year[j]
    sub <- df_in[df_in$location_id == lid & df_in$pheno_year == pyr, ]
    mat <- build_pentad_matrix(sub, indices)
    if (is.null(mat)) next
    vec <- as.numeric(mat)
    for (k in seq_along(indices)) {
      idx_start <- (k-1)*TEMPORAL_BUDGET + 1; idx_end <- k*TEMPORAL_BUDGET
      vec[idx_start:idx_end] <- (vec[idx_start:idx_end] - params$means[k]) / params$sds[k]
    }
    vec[!is.finite(vec)] <- 0
    vecw <- vec * weights
    X_list[[length(X_list) + 1]] <- vecw
    y_vec <- c(y_vec, as.character(traces$Veg[j]))
  }
  if (length(X_list) == 0) return(list(X = NULL, y = NULL))
  X <- do.call(rbind, X_list)
  return(list(X = X, y = y_vec))
}

# Plotting helper: generate prototype plots (similar to original pipeline)
plot_vegetation_prototypes <- function(lib, indices = NULL, out_dir = "prototype_plots", prefix = "veg_prototypes", save_png = TRUE, dpi = 150) {
  if (is.null(lib) || length(lib) == 0) return(NULL)
  if (is.null(indices)) indices <- avail
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required for prototype plotting")
  if (!dir.exists(file.path(OUTPUT_DIR, out_dir)) && save_png) dir.create(file.path(OUTPUT_DIR, out_dir), recursive = TRUE, showWarnings = FALSE)

  rows <- list()
  for (v in names(lib)) {
    variants <- lib[[v]]
    for (var in variants) {
      vec <- as.numeric(var$vec)
      vid <- if (!is.null(var$id)) var$id else if (!is.null(var$variant_id)) var$variant_id else paste0(v, "_unknown")
      n_idx <- length(indices)
      expected_len <- n_idx * TEMPORAL_BUDGET
      if (length(vec) < expected_len) next
      for (k in seq_len(n_idx)) {
        idx_name <- indices[k]
        start <- (k-1) * TEMPORAL_BUDGET + 1
        end <- k * TEMPORAL_BUDGET
        vals <- vec[start:end]
        df_tmp <- data.frame(pentad = seq_len(TEMPORAL_BUDGET), value = vals, Veg = v, variant_id = vid, index = idx_name, stringsAsFactors = FALSE)
        rows[[length(rows) + 1]] <- df_tmp
      }
    }
  }
  if (length(rows) == 0) return(NULL)
  proto_df <- do.call(rbind, rows)
  plots <- list()
  for (idx in unique(proto_df$index)) {
    df_idx <- proto_df[proto_df$index == idx, , drop = FALSE]

    # Base plot: individual variant traces
    p <- ggplot2::ggplot(df_idx, ggplot2::aes(x = pentad, y = value, color = Veg, group = variant_id)) +
      ggplot2::geom_line(alpha = 0.6, size = 0.7) +
      ggplot2::labs(title = sprintf("Vegetation prototypes: %s", idx), x = "Pentad", y = sprintf("%s (endmember center)", idx)) +
      ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_color_brewer(palette = "Set1")

    # Compute mean prototype per vegetation (across variants)
    mean_df <- stats::aggregate(value ~ Veg + pentad, data = df_idx, FUN = mean)

    # Fit candidate models per vegetation and pick the best by AIC
    fit_lines <- list()
    ann <- data.frame(Veg = character(0), model = character(0), x = numeric(0), y = numeric(0), stringsAsFactors = FALSE)
    for (v in unique(mean_df$Veg)) {
      dsub <- mean_df[mean_df$Veg == v, , drop = FALSE]
      if (nrow(dsub) < 4) next  # not enough points to fit robust models
      x <- as.numeric(dsub$pentad)
      y <- as.numeric(dsub$value)

      models <- list()

      # Linear
      lin <- tryCatch(stats::lm(y ~ x), error = function(e) NULL)
      if (!is.null(lin)) models$`Linear` <- lin

      # Exponential: y = a * exp(b * x) + c
      exp_fit <- tryCatch({
        start <- list(a = max(y) - min(y), b = -0.01, c = min(y))
        stats::nls(y ~ a * exp(b * x) + c, start = start, control = stats::nls.control(maxiter = 200))
      }, error = function(e) NULL)
      if (!is.null(exp_fit)) models$`Exponential` <- exp_fit

      # Michaelis-Menten: y = Vmax * x / (K + x) + c
      mm_fit <- tryCatch({
        start <- list(Vmax = max(y) - min(y), K = max(1, median(x)), c = min(y))
        stats::nls(y ~ Vmax * x / (K + x) + c, start = start, control = stats::nls.control(maxiter = 200))
      }, error = function(e) NULL)
      if (!is.null(mm_fit)) models$`Michaelis-Menten` <- mm_fit

      if (length(models) == 0) next

      # Compute AIC and select best model
      aics <- sapply(models, function(m) { tryCatch(AIC(m), error = function(e) Inf) })
      best_name <- names(which.min(aics))[1]
      best_model <- models[[best_name]]

      # Predict over pentad range
      xp <- seq(min(x), max(x), by = 1)
      preds <- tryCatch({
        if (inherits(best_model, "lm")) {
          stats::predict(best_model, newdata = data.frame(x = xp))
        } else {
          stats::predict(best_model, newdata = data.frame(x = xp))
        }
      }, error = function(e) rep(NA_real_, length(xp)))

      fit_lines[[length(fit_lines) + 1]] <- data.frame(pentad = xp, pred = preds, Veg = v, model = best_name, stringsAsFactors = FALSE)

      # Store annotation near the last pentad
      ann <- rbind(ann, data.frame(Veg = v, model = best_name, x = max(xp), y = tail(preds, 1), stringsAsFactors = FALSE))
    }

    # Overlay best-fit lines and annotate model names
    if (length(fit_lines) > 0) {
      fit_df_all <- do.call(rbind, fit_lines)
      p <- p + ggplot2::geom_line(data = fit_df_all, ggplot2::aes(x = pentad, y = pred, color = Veg), size = 1.1)
      # Add model name annotations (placed slightly to the right of the last pentad)
      if (nrow(ann) > 0) {
        for (ri in seq_len(nrow(ann))) {
          # offset label slightly to the right
          p <- p + ggplot2::annotate("text", x = ann$x[ri] + 0.5, y = ann$y[ri], label = ann$model[ri], color = "black", hjust = 0, size = 3)
        }
      }
    }

    plots[[idx]] <- p
    if (save_png) {
      fn <- file.path(OUTPUT_DIR, out_dir, sprintf("%s_%s.png", prefix, idx))
      ggplot2::ggsave(filename = fn, plot = p, width = 8, height = 4, dpi = dpi)
    }
  }
  invisible(plots)
}

# Apply stored normalization params (training) to a dataframe
apply_stored_normalization <- function(df_in, norm_params) {
  if (is.null(norm_params$INDEX_SCALES) || length(norm_params$INDEX_SCALES) == 0) return(df_in)
  df <- df_in
  for (col in names(norm_params$INDEX_SCALES)) {
    if (col %in% names(df)) {
      params <- norm_params$INDEX_SCALES[[col]]
      if (is.list(params) && all(c("mean","sd") %in% names(params))) {
        mu <- params$mean; sigma <- params$sd
        if (!is.finite(sigma) || sigma <= 0) sigma <- 1.0
        df[[col]] <- (df[[col]] - mu) / sigma
      }
    }
  }
  if ("PPI" %in% names(df) && !"PPI_raw" %in% names(df)) df$PPI_raw <- df$PPI
  df
}

# SVM-based threshold evaluator used on OOB set
evaluate_threshold_svm <- function(threshold_quantile, params, train_df, oob_df, indices) {
  weights <- params$weights
  positive_weights <- weights[weights > 0]
  if (length(positive_weights) == 0) stop("No positive weights found in params")
  if (threshold_quantile > 0) { thr <- quantile(positive_weights, threshold_quantile, na.rm = TRUE); test_weights <- weights; test_weights[test_weights < thr] <- 0 } else { thr <- 0; test_weights <- weights }
  n_zeroed <- sum(test_weights == 0); n_kept <- sum(test_weights > 0)
  train_res <- build_matrix_from_df(train_df, indices, params, test_weights)
  oob_res <- build_matrix_from_df(oob_df, indices, params, test_weights)
  if (is.null(train_res$X) || is.null(oob_res$X)) stop("Failed to build train or oob matrices")
  keep_cols <- which(apply(abs(train_res$X), 2, function(x) any(x > 0)))
  if (length(keep_cols) == 0) stop("No features retained after masking/pruning")
  X_train <- train_res$X[, keep_cols, drop = FALSE]; y_train <- as.factor(train_res$y)
  X_oob <- oob_res$X[, keep_cols, drop = FALSE]; y_oob <- as.factor(oob_res$y)
  mu <- colMeans(X_train); sdv <- apply(X_train, 2, sd); sdv[sdv == 0] <- 1
  X_train_s <- sweep(sweep(X_train, 2, mu, "-"), 2, sdv, "/")
  X_oob_s <- sweep(sweep(X_oob, 2, mu, "-"), 2, sdv, "/")
  df_train_svm <- as.data.frame(X_train_s); df_train_svm$target <- y_train
  svm_model <- tryCatch(e1071::svm(target ~ ., data = df_train_svm, kernel = "linear", cost = SVM_COST, scale = FALSE), error = function(e) stop(sprintf("SVM training failed: %s", e$message)))
  preds <- predict(svm_model, as.data.frame(X_oob_s))
  accuracy <- mean(preds == y_oob)
  diag_fracs <- sapply(unique(y_oob), function(cls) { idx <- which(y_oob == cls); if (length(idx) == 0) return(NA_real_); mean(preds[idx] == cls) })
  names(diag_fracs) <- unique(y_oob)
  list(score = mean(diag_fracs, na.rm = TRUE), accuracy = accuracy, n_zeroed = n_zeroed, n_kept = n_kept, threshold_value = thr)
}

# --- MAIN ---
if (!file.exists(INPUT_CSV)) stop(sprintf("Input CSV not found: %s", INPUT_CSV))
raw_df <- read.csv(INPUT_CSV, stringsAsFactors = FALSE)

# Minimal normalization and filtering similar to main script
# Ensure date column is Date
if ("date" %in% names(raw_df)) raw_df$date <- as.Date(raw_df$date)

# Map common vegetation column names and validate key columns before proceeding
if (!"Veg" %in% names(raw_df) && "vegetation" %in% names(raw_df)) {
  raw_df$Veg <- raw_df$vegetation
  cat("[NOTICE] Mapped 'vegetation' -> 'Veg'\n")
}
# If still missing, stop with clear message
if (!"Veg" %in% names(raw_df)) stop(sprintf("Input CSV '%s' missing required column: 'Veg' (or 'vegetation')", INPUT_CSV))
if (!"location_id" %in% names(raw_df)) stop(sprintf("Input CSV '%s' missing required column: 'location_id'", INPUT_CSV))

# Keep selected Vegs
raw_df$Veg <- tolower(as.character(raw_df$Veg))
# Treat certain classes as herbs (match behavior in fit_veg_mixture_mesma.R)
raw_df$Veg <- ifelse(trimws(raw_df$Veg) %in% c("herbs", "alhagi", "salicornia", "halocnemum", "phragmites"), "herbs", raw_df$Veg)
selected_vegs <- c("herbs", "populus", "tamarix", "barren")
df <- raw_df[raw_df$Veg %in% selected_vegs, ]
if (nrow(df) == 0) stop("No training rows after Veg filtering")

# Compute indices from raw bands (same logic as original fit script)
compute_indices_from_bands <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  eps <- 1e-9
  has_bands <- intersect(RAW_BANDS, names(df))
  if (length(has_bands) == 0) return(df)

  if (all(c('nir','red') %in% names(df))) df$DVI <- as.numeric(df$nir) - as.numeric(df$red)
  if (all(c('nir','red') %in% names(df))) df$WDVI <- as.numeric(df$nir) - SOIL_LINE_SLOPE * as.numeric(df$red)
  if (all(c('nir','red') %in% names(df))) df$OSAVI <- (as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + as.numeric(df$red) + 0.16)
  if (all(c('red','green','blue') %in% names(df))) df$MCARI <- ((as.numeric(df$red) - as.numeric(df$green)) - 0.2*(as.numeric(df$red) - as.numeric(df$blue))) * (as.numeric(df$red) / (as.numeric(df$green) + eps))
  if (all(c('green','red') %in% names(df))) df$PRI <- (as.numeric(df$green) - as.numeric(df$red)) / (as.numeric(df$green) + as.numeric(df$red) + eps)
  if (all(c('nir','red') %in% names(df))) df$NIRv <- as.numeric(df$nir) * ((as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + as.numeric(df$red) + eps))
  if (all(c('red','blue','nir') %in% names(df))) df$PSRI <- (as.numeric(df$red) - as.numeric(df$blue)) / (as.numeric(df$nir) + eps)
  if (all(c('nir','swir2') %in% names(df))) df$NBR <- (as.numeric(df$nir) - as.numeric(df$swir2)) / (as.numeric(df$nir) + as.numeric(df$swir2) + eps)

  if (all(c('green','red','nir','swir1','swir2','blue') %in% names(df))) {
    df$TCB <- 0.3029 * as.numeric(df$blue) + 0.2786 * as.numeric(df$green) + 0.4733 * as.numeric(df$red) + 0.5599 * as.numeric(df$nir) + 0.508 * as.numeric(df$swir1) + 0.1872 * as.numeric(df$swir2)
    df$TCG <- -0.2941 * as.numeric(df$blue) - 0.243 * as.numeric(df$green) - 0.5424 * as.numeric(df$red) + 0.7276 * as.numeric(df$nir) + 0.0713 * as.numeric(df$swir1) - 0.1608 * as.numeric(df$swir2)
    df$TCW <- 0.1511 * as.numeric(df$blue) + 0.1973 * as.numeric(df$green) + 0.3283 * as.numeric(df$red) + 0.3407 * as.numeric(df$nir) - 0.7117 * as.numeric(df$swir1) - 0.4559 * as.numeric(df$swir2)
    df$GVI <- df$TCG
  }

  if (all(c('nir','red') %in% names(df))) df$MSAVI2 <- (2 * as.numeric(df$nir) + 1 - sqrt(pmax(0, (2 * as.numeric(df$nir) + 1)^2 - 8 * (as.numeric(df$nir) - as.numeric(df$red))))) / 2
  if (all(c('nir','swir1') %in% names(df))) df$NDMI <- (as.numeric(df$nir) - as.numeric(df$swir1)) / (as.numeric(df$nir) + as.numeric(df$swir1) + eps)
  if (all(c('nir','red','blue') %in% names(df))) df$EVI <- 2.5 * ((as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + 6 * as.numeric(df$red) - 7.5 * as.numeric(df$blue) + 1 + eps))
  if (all(c('swir1','swir2') %in% names(df))) df$NDTI <- (as.numeric(df$swir1) - as.numeric(df$swir2)) / (as.numeric(df$swir1) + as.numeric(df$swir2) + eps)
  # snow-index (NDSI) removed; pipeline uses NDDI (dust) only
  if (all(c('red','nir') %in% names(df))) df$NDDI <- (as.numeric(df$red) - as.numeric(df$nir)) / (as.numeric(df$red) + as.numeric(df$nir) + eps)
  if (all(c('swir1','red','swir2') %in% names(df))) df$SATVI <- ((as.numeric(df$swir1) - as.numeric(df$red)) / (as.numeric(df$swir1) + as.numeric(df$red) + 0.5 + eps)) * 1.5 - (as.numeric(df$swir2) / 2)
  if (all(c('nir','green') %in% names(df))) df$CIG <- (as.numeric(df$nir) / (as.numeric(df$green) + eps)) - 1
  if (all(c('swir1','red','nir','blue') %in% names(df))) {
    term1 <- as.numeric(df$swir1) + as.numeric(df$red)
    term2 <- as.numeric(df$nir) + as.numeric(df$blue)
    df$BSI <- (term1 - term2) / (term1 + term2 + eps)
  }
  if (all(c('green','nir') %in% names(df))) df$NDWI <- (as.numeric(df$green) - as.numeric(df$nir)) / (as.numeric(df$green) + as.numeric(df$nir) + eps)
  if (all(c('swir1','nir') %in% names(df))) df$NDBI <- (as.numeric(df$swir1) - as.numeric(df$nir)) / (as.numeric(df$swir1) + as.numeric(df$nir) + eps)
  if (all(c('swir1','nir') %in% names(df))) df$MSI <- as.numeric(df$swir1) / (as.numeric(df$nir) + eps)
  if (all(c('green','red','blue') %in% names(df))) df$VARI <- (as.numeric(df$green) - as.numeric(df$red)) / (as.numeric(df$green) + as.numeric(df$red) - as.numeric(df$blue) + eps)
  if (all(c('nir','blue','red') %in% names(df))) df$SIPI <- (as.numeric(df$nir) - as.numeric(df$blue)) / (as.numeric(df$nir) - as.numeric(df$red) + eps)
  if (all(c('nir','red','blue') %in% names(df))) {
    rb <- (2 * as.numeric(df$red) - as.numeric(df$blue))
    df$ARVI <- (as.numeric(df$nir) - rb) / (as.numeric(df$nir) + rb + eps)
  }
  if (all(c('nir','green') %in% names(df))) df$GNDVI <- (as.numeric(df$nir) - as.numeric(df$green)) / (as.numeric(df$nir) + as.numeric(df$green) + eps)
  if ('NIRv' %in% names(df)) df$NIRv <- df$NIRv * 1.3
  df
}

df <- compute_indices_from_bands(df)
# Ensure NDVI is not used (keep backup) to be consistent with MESMA fitter
if ("NDVI" %in% names(df)) {
  df$NDVI_orig <- df$NDVI
  df$NDVI <- NULL
  cat("[NOTICE] NDVI column removed from dataframe to prevent usage (backup stored in NDVI_orig)\n")
}

# Simple index normalization (use column set from data)
cols_to_norm <- intersect(names(df), c("WDVI","GVI","NIRv","PSRI","MSAVI2","NDMI","EVI","NDTI","SATVI","CIG","BSI","TCW","TCB","TCG","NDWI","NDBI","MSI","SIPI","ARVI","GNDVI"))
norm_res <- normalize_mesma_data(df, cols_to_norm)
df <- norm_res$df

# Ensure pheno_year and doy
if (!"pheno_year" %in% names(df)) df$pheno_year <- assign_pheno_year(df$date)
if (!"doy" %in% names(df)) df$doy <- pheno_doy(df$date)

# Create df_train (filter to TRAIN_YEARS if available in environment else use all)
TRAIN_YEARS <- if (exists("TRAIN_YEARS")) get("TRAIN_YEARS") else NULL
if (!is.null(TRAIN_YEARS) && length(TRAIN_YEARS) > 0) {
  df_train <- df[df$pheno_year %in% TRAIN_YEARS, , drop = FALSE]
} else df_train <- df

# Stratified OOB split by location/Veg
set.seed(MESMA_SEED)
loc_veg_summary <- df_train %>% group_by(location_id, Veg) %>% summarise(n = n(), .groups = "drop") %>% arrange(location_id, desc(n)) %>% group_by(location_id) %>% slice(1) %>% ungroup()
unique_vegs <- unique(loc_veg_summary$Veg)
oob_locs_list <- vector("list", length(unique_vegs))
for (i in seq_along(unique_vegs)) {
  v <- unique_vegs[i]
  v_locs <- loc_veg_summary$location_id[loc_veg_summary$Veg == v]
  n_v <- length(v_locs)
  n_oob <- ceiling(n_v * OOB_TUNING_FRACTION)
  if (n_oob > 0 && n_v > 1) {
    n_oob <- min(n_oob, n_v - 1)
    selected <- sample(v_locs, n_oob)
    oob_locs_list[[i]] <- data.frame(location_id = selected, Veg = v, stringsAsFactors = FALSE)
  } else oob_locs_list[[i]] <- data.frame()
}
oob_locs_df <- do.call(rbind, oob_locs_list)
if (is.null(oob_locs_df) || nrow(oob_locs_df) == 0) stop("No OOB locations could be selected; increase dataset or OOB fraction")
oob_location_ids <- as.character(oob_locs_df$location_id)
df_train_oob <- df_train[df_train$location_id %in% oob_location_ids, , drop = FALSE]
df_train_model <- df_train[!df_train$location_id %in% oob_location_ids, , drop = FALSE]
cat(sprintf("[OOB SPLIT] OOB rows: %d (locations=%d); Model rows: %d (locations=%d)\n", nrow(df_train_oob), length(unique(df_train_oob$location_id)), nrow(df_train_model), length(unique(df_train_model$location_id))))

# Build multi-class training table for PCA-LDA (exclude validation if any)
train_for_pipeline <- df_train_model
multi_class_data <- mutate(train_for_pipeline, target_class = tolower(as.character(Veg))) %>% filter(!is.na(target_class) & target_class != "")
avail <- intersect(names(df), cols_to_norm)  # feature set to use
cat(sprintf("[PIPELINE] Using %d features for PCA-LDA: %s\n", length(avail), paste(head(avail, 10), collapse=", ")))

# Train PCA-LDA feature pipeline
MESMA_PARAMS_INITIAL <- train_feature_pipeline(multi_class_data, "target_class", avail, use_lda = FALSE)
if (is.null(MESMA_PARAMS_INITIAL)) stop("PCA-LDA training failed or returned NULL")
MESMA_PARAMS_INITIAL$weights[is.na(MESMA_PARAMS_INITIAL$weights)] <- 0
MESMA_PARAMS_INITIAL$weights[!is.finite(MESMA_PARAMS_INITIAL$weights)] <- 0

# Evaluate threshold candidates using OOB SVM evaluation
threshold_results <- data.frame(threshold_quantile = numeric(), threshold_value = numeric(), score = numeric(), accuracy = numeric(), n_zeroed = numeric(), n_kept = numeric(), stringsAsFactors = FALSE)
for (thresh_q in PCA_LDA_THRESHOLD_CANDIDATES) {
  res <- tryCatch({ evaluate_threshold_svm(thresh_q, MESMA_PARAMS_INITIAL, df_train_model, df_train_oob, avail) }, error = function(e) { cat(sprintf("[WARN] Threshold %.2f evaluation failed: %s\n", thresh_q, e$message)); NULL })
  if (is.null(res)) next
  threshold_results <- rbind(threshold_results, data.frame(threshold_quantile = thresh_q, threshold_value = res$threshold_value, score = res$score, accuracy = res$accuracy, n_zeroed = res$n_zeroed, n_kept = res$n_kept, stringsAsFactors = FALSE))
  cat(sprintf("  Threshold q=%.2f (val=%.4f): score=%.3f, accuracy=%.3f, zeroed=%d, kept=%d\n", thresh_q, res$threshold_value, res$score, res$accuracy, res$n_zeroed, res$n_kept))
}
if (nrow(threshold_results) == 0) stop("No successful threshold evaluations")
best_idx <- which.max(threshold_results$score + 0.001 * threshold_results$n_zeroed / pmax(1, max(threshold_results$n_zeroed)))
optimal_threshold <- threshold_results$threshold_quantile[best_idx]
optimal_threshold_value <- threshold_results$threshold_value[best_idx]
cat(sprintf("[RESULT] Optimal threshold quantile=%.2f (value=%.4f), score=%.3f, accuracy=%.3f\n", optimal_threshold, optimal_threshold_value, threshold_results$score[best_idx], threshold_results$accuracy[best_idx]))

# Train final SVM using optimal threshold and evaluate on OOB
pos_w <- MESMA_PARAMS_INITIAL$weights[MESMA_PARAMS_INITIAL$weights > 0]
final_weights <- MESMA_PARAMS_INITIAL$weights
if (optimal_threshold > 0) {
  thr <- quantile(pos_w, optimal_threshold, na.rm = TRUE)
  final_weights[final_weights < thr] <- 0
}
train_res <- build_matrix_from_df(df_train_model, avail, MESMA_PARAMS_INITIAL, final_weights)
oob_res <- build_matrix_from_df(df_train_oob, avail, MESMA_PARAMS_INITIAL, final_weights)
keep_cols <- which(apply(abs(train_res$X), 2, function(x) any(x > 0)))
X_train <- train_res$X[, keep_cols, drop = FALSE]; y_train <- as.factor(train_res$y)
X_oob <- oob_res$X[, keep_cols, drop = FALSE]; y_oob <- as.factor(oob_res$y)
mu <- colMeans(X_train); sdv <- apply(X_train, 2, sd); sdv[sdv == 0] <- 1
X_train_s <- sweep(sweep(X_train, 2, mu, "-"), 2, sdv, "/")
X_oob_s <- sweep(sweep(X_oob, 2, mu, "-"), 2, sdv, "/")

cat("[SVM] Training final linear SVM (probability estimates enabled)...\n")
df_train_svm <- as.data.frame(X_train_s); df_train_svm$target <- y_train
svm_model <- e1071::svm(target ~ ., data = df_train_svm, kernel = "linear", cost = SVM_COST, scale = FALSE, probability = TRUE)
# Evaluate on OOB with probabilities
preds <- predict(svm_model, as.data.frame(X_oob_s), probability = TRUE)
prob_mat <- attr(preds, "probabilities")
if (is.null(prob_mat)) {
  acc <- mean(preds == y_oob)
} else {
  acc <- mean(preds == y_oob)
}
cat(sprintf("[SVM] Final OOB accuracy = %.4f (n=%d)\n", acc, length(y_oob)))

# --- Confusion matrix and metrics for OOB ---
tryCatch({
  cm <- table(True = y_oob, Predicted = preds)
  # Save raw confusion matrix
  cm_file <- file.path(OUTPUT_DIR, "svm_oob_confusion_matrix.csv")
  write.csv(as.data.frame.matrix(cm), file = cm_file, row.names = TRUE)
  cat(sprintf("[SAVE] Wrote OOB confusion matrix to %s\n", cm_file))

  # Save row-normalized confusion matrix (recall per class)
  cm_row_file <- file.path(OUTPUT_DIR, "svm_oob_confusion_matrix_row_normalized.csv")
  write.csv(as.data.frame.matrix(prop.table(cm, 1)), file = cm_row_file, row.names = TRUE)
  cat(sprintf("[SAVE] Wrote row-normalized OOB confusion matrix to %s\n", cm_row_file))

  # Per-class precision and recall
  recall <- ifelse(rowSums(cm) > 0, diag(cm) / rowSums(cm), NA_real_)
  precision <- ifelse(colSums(cm) > 0, diag(cm) / colSums(cm), NA_real_)
  metrics <- data.frame(class = names(recall), recall = as.numeric(recall), precision = as.numeric(precision))
  write.csv(metrics, file = file.path(OUTPUT_DIR, "svm_oob_confusion_metrics.csv"), row.names = FALSE)
  cat("[SAVE] Wrote OOB confusion metrics to svm_oob_confusion_metrics.csv\n")

  # Plot heatmap using ggplot2 if available
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    cm_df <- as.data.frame(cm)
    names(cm_df) <- c("True", "Predicted", "Count")
    p <- ggplot2::ggplot(cm_df, ggplot2::aes(x = Predicted, y = True, fill = Count)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = Count), color = "black") +
      ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "SVM OOB Confusion Matrix", x = "Predicted", y = "True")
    ggsave_file <- file.path(OUTPUT_DIR, "svm_oob_confusion_matrix.png")
    ggplot2::ggsave(ggsave_file, plot = p, width = 6, height = 5, dpi = 150)
    cat(sprintf("[SAVE] Wrote OOB confusion matrix plot to %s\n", ggsave_file))
  } else {
    cat("[WARN] ggplot2 not available; skipping confusion matrix plot\n")
  }
}, error = function(e) {
  cat(sprintf("[WARN] Failed to compute or save confusion matrix: %s\n", e$message))
})

# --- Location-level bootstrap percentages (OOB) ---
# Computes the percent of unique locations predicted as each veg class with bootstrap CIs
tryCatch({
  # Reconstruct trace-level location ids corresponding to OOB predictions
  traces_oob <- unique(df_train_oob[, c("location_id", "pheno_year", "Veg")])
  if (length(preds) != nrow(traces_oob)) {
    warning("OOB predictions length does not match reconstructed OOB traces; skipping location-level bootstrap")
  } else {
    # Majority vote per location (across years if multiple)
    df_pred <- data.frame(location_id = traces_oob$location_id, pred = as.character(preds), stringsAsFactors = FALSE)
    library(dplyr)
    loc_pred <- df_pred %>% dplyr::group_by(location_id) %>% dplyr::summarise(pred = names(sort(table(pred), decreasing = TRUE))[1], .groups = 'drop')

    all_classes <- sort(unique(c(as.character(y_oob), as.character(preds))))
    nloc <- nrow(loc_pred)
    if (nloc == 0) stop("No locations found in OOB traces for bootstrap")

    BOOTSTRAP_B <- ifelse(exists("BOOTSTRAP_B"), BOOTSTRAP_B, 1000)
    set.seed(MESMA_SEED)
    bs_mat <- matrix(0, nrow = length(all_classes), ncol = BOOTSTRAP_B)
    rownames(bs_mat) <- all_classes
    for (b in seq_len(BOOTSTRAP_B)) {
      samp <- sample(loc_pred$location_id, size = nloc, replace = TRUE)
      samp_preds <- loc_pred$pred[match(samp, loc_pred$location_id)]
      tab <- table(factor(samp_preds, levels = all_classes))
      bs_mat[, b] <- as.numeric(tab) / nloc
    }
    mean_prop <- rowMeans(bs_mat, na.rm = TRUE)
    lower <- apply(bs_mat, 1, quantile, probs = 0.025, na.rm = TRUE)
    upper <- apply(bs_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
    out_df <- data.frame(class = rownames(bs_mat), prop = mean_prop, lower = lower, upper = upper)
    write.csv(out_df, file = file.path(OUTPUT_DIR, "svm_oob_location_fraction_bootstrap.csv"), row.names = FALSE)
    cat("[SAVE] Wrote location-level bootstrap fractions to svm_oob_location_fraction_bootstrap.csv\n")

    # Plot
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p2 <- ggplot2::ggplot(out_df, ggplot2::aes(x = class, y = prop)) +
        ggplot2::geom_col(fill = "steelblue") +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.2) +
        ggplot2::labs(title = "OOB: % of locations predicted per class (bootstrap 95% CI)", x = "Class", y = "Proportion of locations") +
        ggplot2::theme_minimal()
      ggplot2::ggsave(file.path(OUTPUT_DIR, "svm_oob_location_fraction.png"), plot = p2, width = 6, height = 4, dpi = 150)
      cat("[SAVE] Wrote OOB location fraction plot to svm_oob_location_fraction.png\n")
    } else {
      cat("[WARN] ggplot2 not available; skipping location-level fraction plot\n")
    }
  }
}, error = function(e) {
  cat(sprintf("[WARN] Failed to compute location-level bootstrap fractions: %s\n", e$message))
})

# Save training normalization parameters for inference use
TRAINING_NORM_PARAMS <- list(INDEX_SCALES = norm_res$INDEX_SCALES)
saveRDS(TRAINING_NORM_PARAMS, file = file.path(OUTPUT_DIR, "training_norm_params.rds"))
cat("[SAVE] Saved training normalization params to training_norm_params.rds\n")

# Build simple prototype library (medoid per class) and generate prototype plots to match original pipeline
res_lib <- list()
unique_classes <- unique(multi_class_data$target_class)
for (cls in unique_classes) {
  cls_data <- multi_class_data[multi_class_data$target_class == cls, ]
  if (nrow(cls_data) == 0) next
  cls_traces <- unique(cls_data[, c("location_id", "pheno_year")])
  cls_vecs <- list()
  for (j in seq_len(nrow(cls_traces))) {
    lid <- cls_traces$location_id[j]
    pyr <- cls_traces$pheno_year[j]
    sub <- cls_data[cls_data$location_id == lid & cls_data$pheno_year == pyr, ]
    mat <- build_pentad_matrix(sub, avail)
    if (is.null(mat)) next
    vec <- as.numeric(mat)
    # z-score using training moments
    for (k in seq_along(avail)) {
      idx_start <- (k-1)*TEMPORAL_BUDGET + 1; idx_end <- k*TEMPORAL_BUDGET
      vec[idx_start:idx_end] <- (vec[idx_start:idx_end] - MESMA_PARAMS_INITIAL$means[k]) / MESMA_PARAMS_INITIAL$sds[k]
    }
    vec[!is.finite(vec)] <- 0
    cls_vecs[[length(cls_vecs) + 1]] <- vec
  }
  if (length(cls_vecs) == 0) next
  cls_mat <- do.call(rbind, cls_vecs)
  if (nrow(cls_mat) > 1) {
    dist_mat <- as.matrix(dist(cls_mat)); medoid_idx <- which.min(rowSums(dist_mat)); medoid_vec <- cls_mat[medoid_idx, ]
  } else medoid_vec <- cls_mat[1, ]
  res_lib[[cls]] <- list(list(vec = medoid_vec, variant_id = paste0(cls, "_medoid")))
}

# Generate prototype plots (matches original behavior)
tryCatch({
  plot_vegetation_prototypes(res_lib, indices = avail, out_dir = file.path(OUTPUT_DIR, "prototype_plots"))
  cat(sprintf("[PLOTS] Generated prototype plots (dir=prototype_plots)\n"))
}, error = function(e) cat(sprintf("[WARN] Prototype plotting failed: %s\n", e$message)))

# If inference CSV exists, load and apply stored normalization so SVM uses same inference data handling
if (file.exists(INFERENCE_CSV)) {
  cat(sprintf("[INFERENCE] Found INFERENCE_CSV: %s - loading and applying stored normalization\n", INFERENCE_CSV))
  df_inf <- read.csv(INFERENCE_CSV, stringsAsFactors = FALSE)
  if ("date" %in% names(df_inf)) df_inf$date <- as.Date(df_inf$date)
  if (!"pheno_year" %in% names(df_inf)) df_inf$pheno_year <- assign_pheno_year(df_inf$date)
  if (!"doy" %in% names(df_inf)) df_inf$doy <- pheno_doy(df_inf$date)
  df_inf <- apply_stored_normalization(df_inf, TRAINING_NORM_PARAMS)
  saveRDS(df_inf, file = file.path(OUTPUT_DIR, "df_inference_processed.rds"))
  cat("[INFERENCE] Saved processed inference data to df_inference_processed.rds\n")

  # Build per-trace pentad vectors for inference and predict class probabilities using final SVM
  traces_inf <- unique(df_inf[, c("location_id", "pheno_year")])
  pred_rows <- list()
  if (nrow(traces_inf) == 0) {
    cat("[INFERENCE] No traces found in inference CSV after grouping; skipping prediction\n")
  } else {
    for (j in seq_len(nrow(traces_inf))) {
      lid <- traces_inf$location_id[j]
      pyr <- traces_inf$pheno_year[j]
      sub <- df_inf[df_inf$location_id == lid & df_inf$pheno_year == pyr, , drop = FALSE]
      mat <- build_pentad_matrix(sub, avail)
      if (is.null(mat)) next
      vec <- as.numeric(mat)
      # Z-score using training moments
      for (k in seq_along(avail)) {
        idx_start <- (k-1)*TEMPORAL_BUDGET + 1; idx_end <- k*TEMPORAL_BUDGET
        vec[idx_start:idx_end] <- (vec[idx_start:idx_end] - MESMA_PARAMS_INITIAL$means[k]) / MESMA_PARAMS_INITIAL$sds[k]
      }
      vec[!is.finite(vec)] <- 0
      vecw <- vec * final_weights
      # Ensure same columns as training keep_cols
      if (length(keep_cols) == 0) {
        cat("[INFERENCE] No keep_cols available; cannot make predictions\n")
        next
      }
      xrow <- vecw[keep_cols]
      # Scale using training mu/sdv
      xrow_s <- (xrow - mu) / sdv
      df_x <- as.data.frame(t(xrow_s))
      preds_inf <- tryCatch(predict(svm_model, df_x, probability = TRUE), error = function(e) { cat(sprintf("[INFERENCE] Prediction failed for loc=%s year=%s: %s\n", lid, pyr, e$message)); NA })
      if (is.na(preds_inf)) next
      probs <- attr(preds_inf, "probabilities")
      if (is.null(probs)) {
        # Fall back to hard prediction only
        pred_class <- as.character(preds_inf)
        # Create zeroed probability vector and set predicted class = 1
        lvl <- levels(as.factor(y_train))
        pv <- setNames(rep(0, length(lvl)), lvl); pv[pred_class] <- 1
      } else {
        pv <- as.numeric(probs[1, ])
        names(pv) <- colnames(probs)
        pred_class <- colnames(probs)[which.max(pv)]
      }
      row_out <- c(location_id = lid, pheno_year = pyr, predicted = pred_class, as.list(pv))
      pred_rows[[length(pred_rows) + 1]] <- row_out
    }
    if (length(pred_rows) > 0) {
      pred_df <- do.call(rbind.data.frame, lapply(pred_rows, function(x) { as.data.frame(as.list(x), stringsAsFactors = FALSE) }))
      # Ensure probability columns are numeric
      prob_cols <- setdiff(names(pred_df), c("location_id", "pheno_year", "predicted"))
      for (cname in prob_cols) pred_df[[cname]] <- as.numeric(as.character(pred_df[[cname]]))
      write.csv(pred_df, file = file.path(OUTPUT_DIR, "svm_inference_cover.csv"), row.names = FALSE)
      cat(sprintf("[INFERENCE] Saved SVM inference cover predictions to svm_inference_cover.csv (rows=%d)\n", nrow(pred_df)))
    } else {
      cat("[INFERENCE] No valid inference traces produced predictions\n")
    }
  }

} else {
  cat(sprintf("[INFERENCE] INFERENCE_CSV not found at %s - skipping inference processing\n", INFERENCE_CSV))
}

# Save final artifacts (SVM model, thresholds, and training splits)
saveRDS(svm_model, file = file.path(OUTPUT_DIR, "svm_pca_lda_model.rds"))
write.csv(threshold_results, file = file.path(OUTPUT_DIR, "svm_threshold_results.csv"), row.names = FALSE)
saveRDS(df_train_model, file = file.path(OUTPUT_DIR, "df_train_model.rds"))
saveRDS(df_train_oob, file = file.path(OUTPUT_DIR, "df_train_oob.rds"))
cat("[SAVE] Saved svm model, threshold results, and train/oob datasets\n")

# Exit normally
invisible(NULL)
