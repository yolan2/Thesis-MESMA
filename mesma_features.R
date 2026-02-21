# meia_features.R — small set of canonical MESMA helpers copied from
# `fit_veg_mixture_mesma.R` so lightweight pipelines (SVM, tests) can reuse
# identical L2 / PCA / LDA representations without sourcing the full MESMA
# training script (which executes heavy top-level work).

# L2-normalize a feature vector per observation (whole-vector)
l2_normalize_perindex <- function(vec, n_indices = NULL, n_bins = NULL) {
  v <- as.numeric(vec)
  v_clean <- v
  v_clean[!is.finite(v_clean)] <- 0
  nrm <- sqrt(sum(v_clean^2))
  if (!is.finite(nrm) || nrm < 1e-9) return(v)
  v / nrm
}

mesma_apply_representation_vec <- function(vec_raw, n_base_idx, n_bins, l2_normalize) {
  if (isTRUE(l2_normalize)) {
    return(l2_normalize_perindex(vec_raw, n_base_idx, n_bins))
  }
  vec_raw
}

mesma_apply_representation_mat <- function(mat_raw, n_base_idx, n_bins, l2_normalize) {
  if (!isTRUE(l2_normalize)) return(mat_raw)
  t(apply(mat_raw, 1, function(r) {
    l2_normalize_perindex(r, n_base_idx, n_bins)
  }))
}

mesma_zscore_vec_by_index <- function(vec, indices, means, sds, n_bins, eps_sigma = NULL) {
  if (is.null(eps_sigma)) {
    eps_sigma <- if (exists("EPS_SIGMA", inherits = TRUE)) get("EPS_SIGMA", inherits = TRUE) else 1e-8
  }

  if (length(indices) == 0) return(vec)

  idx_pos <- match(indices, names(means))
  out <- vec
  for (k in seq_along(indices)) {
    pos <- idx_pos[k]
    if (is.na(pos)) next
    idx_start <- (k - 1) * n_bins + 1
    idx_end <- k * n_bins
    mu <- means[pos]
    sigma <- sds[pos]
    if (!is.finite(sigma) || sigma < eps_sigma) sigma <- eps_sigma
    out[idx_start:idx_end] <- (out[idx_start:idx_end] - mu) / sigma
  }
  out[!is.finite(out)] <- 0
  out
}

mesma_zscore_mat_by_index <- function(mat, indices, means, sds, n_bins, eps_sigma = NULL) {
  if (is.null(eps_sigma)) {
    eps_sigma <- if (exists("EPS_SIGMA", inherits = TRUE)) get("EPS_SIGMA", inherits = TRUE) else 1e-8
  }

  if (is.null(mat) || nrow(mat) == 0 || length(indices) == 0) return(mat)

  idx_pos <- match(indices, names(means))
  out <- mat
  for (k in seq_along(indices)) {
    pos <- idx_pos[k]
    if (is.na(pos)) next
    idx_start <- (k - 1) * n_bins + 1
    idx_end <- k * n_bins
    mu <- means[pos]
    sigma <- sds[pos]
    if (!is.finite(sigma) || sigma < eps_sigma) sigma <- eps_sigma
    out[, idx_start:idx_end] <- (out[, idx_start:idx_end] - mu) / sigma
  }
  out[!is.finite(out)] <- 0
  out
}

# Canonical pentad builder (copies MESMA behavior)
build_pentad_matrix <- function(dly_year, avail_idx) {
  if (is.null(dly_year) || nrow(dly_year) == 0) return(NULL)

  if (!"doy" %in% names(dly_year) || any(is.na(dly_year$doy))) {
    # Use the centralized `pheno_doy()` implementation (May 1 = day 1)
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

  pentad_mat
}

# Safe LDA wrapper (same behavior as canonical MESMA)
safe_lda_call <- function(X_pca, y, min_n_pcs = 2) {
  if (is.null(X_pca) || ncol(X_pca) < min_n_pcs) return(NULL)
  curr_n_pcs <- ncol(X_pca)
  for (attempt in seq_len(max(1, curr_n_pcs - min_n_pcs + 1))) {
    tryCatch({
      res <- MASS::lda(X_pca[, 1:curr_n_pcs, drop = FALSE], grouping = y)
      return(res)
    }, warning = function(w) {
      warn_msg <- conditionMessage(w)
      if (grepl("collinear", warn_msg, ignore.case = TRUE)) curr_n_pcs <- curr_n_pcs - 1 else return(res)
    }, error = function(e) return(NULL))
  }
  NULL
}

# Canonical PCA->LDA trainer (copied from fit_veg_mixture_mesma.R)
train_feature_pipeline <- function(df, class_col, feature_cols, use_lda = TRUE) {
  # This is a near-identical copy of the MESMA implementation used for training
  # PCA/LDA weights. It intentionally reproduces L2-normalization, z-scoring,
  # PCA computation, PC pruning, and LDA weighting.
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
  X_mat_raw <- do.call(rbind, X_raw)
  n_bins_local <- TEMPORAL_BUDGET
  n_idx_local <- length(feature_cols)
  l2_only_mode <- isTRUE(ENABLE_LDA_L2_NORMALIZATION)
  if (l2_only_mode) {
    cat(sprintf("  L2-normalizing training samples (per-observation) for %d indices (%d pentads each)...\n", n_idx_local, n_bins_local))
    X_mat <- t(apply(X_mat_raw, 1, function(r) { l2_normalize_perindex(r, n_idx_local, n_bins_local) }))
    all_feature_cols <- feature_cols
    cat(sprintf("  L2 normalization ENABLED: using %d indices.\n", n_idx_local))
  } else {
    X_mat <- X_mat_raw
    all_feature_cols <- feature_cols
    cat(sprintf("  L2 normalization DISABLED: using %d raw indices.\n", n_idx_local))
  }
  apply_zscore <- if (exists("ENABLE_ZSCORE_AFTER_L2")) isTRUE(ENABLE_ZSCORE_AFTER_L2) else TRUE
  n_total_indices <- length(all_feature_cols)
  global_means <- numeric(n_total_indices); global_sds <- numeric(n_total_indices); names(global_means) <- all_feature_cols; names(global_sds) <- all_feature_cols
  X_z <- X_mat
  if (apply_zscore) {
    cat(sprintf("  Computing Z-score parameters for %d indices...\n", n_total_indices))
    for(k in seq_along(all_feature_cols)) {
      col_idx_start <- (k-1)*n_bins_local + 1; col_idx_end <- k*n_bins_local
      vals <- X_mat[, col_idx_start:col_idx_end]
      mu <- mean(vals, na.rm=TRUE); sigma <- sd(vals, na.rm=TRUE)
      if(sigma == 0 || is.na(sigma)) sigma <- 1
      global_means[k] <- mu; global_sds[k] <- sigma
      X_z[, col_idx_start:col_idx_end] <- (vals - mu) / sigma
    }
    X_z[!is.finite(X_z)] <- 0
  } else {
    cat("  Z-scoring DISABLED: using L2-normalized features directly\n")
    global_means[] <- 0; global_sds[] <- 1; X_z <- X_mat; X_z[!is.finite(X_z)] <- 0
  }
  cat("  Computing PCA-LDA weights...\n")
  vars <- apply(X_z, 2, var); keep_cols <- vars > 1e-9; X_pca_in <- X_z[, keep_cols, drop=FALSE]
  pca_res <- prcomp(X_pca_in, center = FALSE, scale. = FALSE)
  cum_var <- cumsum(pca_res$sdev^2) / sum(pca_res$sdev^2); n_pcs <- which(cum_var > 0.95)[1]
  if(is.na(n_pcs)) n_pcs <- ncol(pca_res$x)
  class_counts <- table(y_labels); n_min <- if (length(class_counts) > 0) min(class_counts) else 0; n_classes <- length(unique(y_labels))
  min_pcs_for_lda <- max(1, n_classes - 1)
  max_pcs_for_lda <- min(20, max(min_pcs_for_lda, n_min - 2))
  if (n_pcs > max_pcs_for_lda) n_pcs <- max_pcs_for_lda
  if (n_pcs < min_pcs_for_lda) stop("Not enough degrees of freedom for LDA")
  lda_res <- safe_lda_call(pca_res$x[, 1:n_pcs, drop = FALSE], as.factor(y_labels), min_n_pcs = min_pcs_for_lda)
  if (is.null(lda_res)) stop("LDA could not be computed")
  W_pc <- lda_res$scaling; R <- pca_res$rotation[, 1:n_pcs, drop = FALSE]; W_std <- R %*% W_pc
  svd_vals <- lda_res$svd; prop <- svd_vals / sum(svd_vals)
  if (ncol(W_std) > 1) { n_dim <- min(length(prop), ncol(W_std)); weights_clean <- rowSums(abs(W_std[, 1:n_dim, drop=FALSE]) %*% diag(prop[1:n_dim], nrow=n_dim)) } else { weights_clean <- abs(W_std[, 1]) }
  final_weights <- numeric(ncol(X_z)); final_weights[keep_cols] <- weights_clean
  n_full_features <- ncol(X_z); n_lda_dims <- ncol(W_pc)
  std_to_lda <- matrix(0, nrow = n_full_features, ncol = n_lda_dims); std_to_lda[keep_cols, ] <- W_std
  std_to_lda_abs_denom <- colSums(abs(std_to_lda)); std_to_lda_abs_denom[!is.finite(std_to_lda_abs_denom) | std_to_lda_abs_denom < 1e-12] <- 1
  return(list(means = global_means, sds = global_sds, weights = final_weights, pca_lda = list(feature_space = "pca_lda", keep_cols = keep_cols, std_to_lda = std_to_lda, std_to_lda_abs_denom = std_to_lda_abs_denom), indices = all_feature_cols, base_indices = feature_cols, l2_normalize = l2_only_mode, zscore_applied = apply_zscore))
}
