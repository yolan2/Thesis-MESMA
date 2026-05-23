# preprocess_data.R
# This script preprocesses the raw data for MESMA analysis.
# It loads the raw CSV, applies all preprocessing steps, and saves the preprocessed data.

# --- Source shared config and helpers ---
# All tunable parameters are defined in mesma_config.R (single source of truth).
# Override any parameter AFTER the source() call below if needed for this run.
stopifnot(file.exists("mesma_config.R"))
source("mesma_config.R")

stopifnot(file.exists("ppi_helpers.R"))
source("ppi_helpers.R")

stopifnot(file.exists("mesma_helpers.R"))
source("mesma_helpers.R")
set_mesma_seed()

# --- required packages -----------------------------------------------------
# load in bulk using require/ lapply
pkgs <- c(
  "dplyr", "purrr", "readr", "ggplot2",
  "magrittr", # pipes
  "future", "future.apply" # parallel backend
)
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed", p))
  }
  library(p, character.only = TRUE)
}


# --- Shared plot constants ---

SPECIES_COLORS <- c("Herbs" = "#9ACD32", "Populus" = "#006400", "Tamarix" = "#D95F02")

# --- Shared utility: normalize vegetation names (lowercase + trim) ---
normalize_veg_name <- function(x) {
  x <- as.character(x)
  x[!nzchar(trimws(x))] <- NA_character_
  tolower(trimws(x))
}

# Compute_haversine_distance_matrix is defined in mesma_helpers.R to keep behaviour consistent across scripts.


# --- Inform user of chosen temporal fill/interpolation method ---
# Logging must occur after the helper above is defined.
interp_method <- resolve_interpolation_method()
cat(sprintf("[NOTICE] Temporal fill option at startup: %s\n", interp_method))

# --- Simple Whittaker smoothing for a 1‑D vector ---
# NA values are treated as missing (weight zero) and will be filled by the
# smoothing operation.  The smoothing strength is controlled by `lambda`,
# which defaults to the global `WHITTAKER_LAMBDA` (if defined) or 1600.
whittaker_smooth <- function(y,
                              lambda = ifelse(exists("WHITTAKER_LAMBDA"), WHITTAKER_LAMBDA, 1600),
                              differences = 2) {
  n <- length(y)
  # weight 1 for finite observations, 0 for NA
  w <- as.numeric(is.finite(y))
  z <- y
  z[!is.finite(z)] <- 0
  D <- diff(diag(n), differences = differences)
  C <- diag(w) + lambda * crossprod(D)
  z_hat <- tryCatch(
    solve(C, w * z),
    error = function(e) {
      warning("whittaker_smooth: linear solve failed, returning original data")
      return(z)
    })
  as.numeric(z_hat)
}

# Shared exponential variogram fitter now lives in mesma_helpers.R
# to ensure consistent behaviour across scripts.

# --- Shared utility: get dominant veg type per location ---
get_dominant_veg_per_location <- function(df) {
  df %>%
    dplyr::group_by(location_id, Veg) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(location_id, dplyr::desc(n)) %>%
    dplyr::group_by(location_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
}

# --- Shared utility: deduplicate coordinates preferring finite lat/lon ---
deduplicate_coords <- function(df) {
  ok <- is.finite(df$lat) & is.finite(df$lon)
  df <- df[order(df$location_id, -as.integer(ok)), , drop = FALSE]
  df[!duplicated(df$location_id), , drop = FALSE]
}

# --- Shared utility: left_join and fill missing Veg from geo source ---
join_and_fill_veg <- function(df, geo_map, join_by = "location_id") {
  joined <- dplyr::left_join(df, geo_map, by = join_by, suffix = c("", ".geo"))
  if ("Veg.geo" %in% names(joined)) {
    joined$Veg <- ifelse(is.na(joined$Veg) | joined$Veg == "", joined$Veg.geo, joined$Veg)
    joined$Veg.geo <- NULL
  }
  joined
}

# --- Shared utility: Olofsson-style area-proportion bias correction ----------
# Applies bias correction per Olofsson et al. 2014.
# Returns data frame with adj_coef, adj_se, etc., or NULL if not applicable.
bias_correct_timeseries <- function(full_data, bias_correction) {
  if (is.null(bias_correction) ||
      is.null(bias_correction$norm_mat) ||
      is.null(bias_correction$class_counts)) {
    return(NULL)
  }
  norm_mat     <- bias_correction$norm_mat      # rows = true class, cols = predicted
  class_counts <- bias_correction$class_counts  # named integer n_{true,.}

  # Stored matrix is P(predicted | true), but Olofsson correction needs
  # P(true | map class). Reconstruct approximate counts and normalize by column.
  count_mat <- sweep(norm_mat, 1, as.numeric(class_counts[rownames(norm_mat)]), `*`)
  count_mat[!is.finite(count_mat)] <- 0
  pred_counts <- colSums(count_mat, na.rm = TRUE)
  post_mat <- count_mat
  for (cj in seq_len(ncol(post_mat))) {
    if (is.finite(pred_counts[cj]) && pred_counts[cj] > 0) {
      post_mat[, cj] <- post_mat[, cj] / pred_counts[cj]
    } else {
      post_mat[, cj] <- 0
    }
  }
  map_classes  <- colnames(post_mat)
  true_classes <- rownames(post_mat)

  years <- sort(unique(full_data$year))
  result_rows <- vector("list", length(years) * nrow(post_mat))
  row_idx <- 0L

  for (yr in years) {
    yr_data <- full_data[full_data$year == yr, , drop = FALSE]

    # W_i: raw MAP proportions keyed by lower-cased Veg name.
    # Keep the original absolute scale instead of renormalizing to vegetation-only.
    W_raw <- setNames(as.numeric(yr_data$global_coef),
                      tolower(as.character(yr_data$Veg)))
    W_raw[!is.finite(W_raw)] <- 0

    shared  <- intersect(names(W_raw), map_classes)
    if (length(shared) == 0) next
    W_sub   <- W_raw[shared]

    for (true_j in true_classes) {
      adj_p <- 0
      se_sq <- 0
      for (map_i in shared) {
        wi  <- W_sub[[map_i]]
        qij <- post_mat[true_j, map_i]
        ni  <- pred_counts[[map_i]]
        if (!is.finite(wi) || !is.finite(qij)) next
        adj_p <- adj_p + wi * qij
        if (!is.na(ni) && ni > 1L) {
          se_sq <- se_sq + wi^2 * qij * (1 - qij) / (ni - 1L)
        }
      }
      se <- sqrt(max(0, se_sq))

      # Preserve original Veg capitalisation
      orig_rows <- yr_data[tolower(as.character(yr_data$Veg)) == true_j, , drop = FALSE]
      orig_veg  <- if (nrow(orig_rows) > 0) as.character(orig_rows$Veg[1]) else true_j

      row_idx <- row_idx + 1L
      result_rows[[row_idx]] <- data.frame(
        year     = yr,
        Veg      = orig_veg,
        adj_coef = adj_p,
        adj_se   = se,
        adj_025  = pmax(0, adj_p - 1.96 * se),
        adj_975  = pmin(1, adj_p + 1.96 * se),
        stringsAsFactors = FALSE
      )
    }

    # Pass-through uncorrected values for classes absent from confusion matrix
    missing_classes <- setdiff(tolower(names(W_raw)), map_classes)
    for (cls in missing_classes) {
      # Find original row to preserve Veg capitalisation
      orig_rows <- yr_data[tolower(as.character(yr_data$Veg)) == cls, , drop = FALSE]
      if (nrow(orig_rows) == 0) next
      row_idx <- row_idx + 1L
      result_rows[[row_idx]] <- data.frame(
        year     = yr,
        Veg      = as.character(orig_rows$Veg[1]),
        adj_coef = W_raw[[cls]],
        adj_se   = NA_real_,
        adj_025  = NA_real_,
        adj_975  = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }

  result_rows <- result_rows[!vapply(result_rows, is.null, logical(1))]
  if (length(result_rows) == 0) return(NULL)
  do.call(rbind, result_rows)
}

# --- Shared utility: plot all inference results for a given method ---
# Replaces 4 nearly identical PPI/MSAVI/NDVI/NoIndex inference plot blocks (~400 lines).
plot_inference_method_results <- function(full_data, method, file_prefix,
                                          use_excluded_years_shade = TRUE,
                                          include_species_plots = TRUE,
                                          bias_correction = NULL) {
  if (is.null(full_data) || nrow(full_data) == 0) return(invisible(NULL))

  veg_norm <- normalize_veg_name(full_data$Veg)
  inf_veg <- full_data[!veg_norm %in% c("barren"), ]
  inf_barren <- full_data[veg_norm %in% c("barren"), ]

  # Helper: base layers for time series plots
  .ts_layers <- function() {
    layers <- list()
    if (use_excluded_years_shade) layers <- c(layers, list(add_excluded_years_shade(is_date = FALSE)))
    layers <- c(layers, list(add_year_lines(is_date = FALSE)))
    layers
  }

  # 1. Vegetation time series
  if (!is.null(inf_veg) && nrow(inf_veg) > 0) {
    p_ts <- ggplot(inf_veg, aes(x = year, y = global_coef, color = Veg, group = Veg)) +
      .ts_layers() +
      geom_line(linewidth = 1) +
      geom_point(show.legend = FALSE) +
      geom_ribbon(aes(ymin = coef_025, ymax = coef_975, fill = Veg), alpha = 0.15, color = NA) +
       labs(title = paste0(method, ": Vegetation Fractions"),
           x = "Year", y = "Total Normalized Fraction", color = "Veg", fill = "Veg") +
      scale_x_continuous(limits = c(1984, NA)) +
      theme_minimal()

    # --- Olofsson-style bias correction overlay (dashed lines + SE band) ---
    adj_veg <- bias_correct_timeseries(full_data, bias_correction)
    if (!is.null(adj_veg) && nrow(adj_veg) > 0) {
      adj_veg <- adj_veg[adj_veg$Veg %in% unique(inf_veg$Veg), , drop = FALSE]
      if (nrow(adj_veg) > 0) {
        p_ts <- p_ts +
          geom_line(data = adj_veg,
                    aes(x = year, y = adj_coef, color = Veg, group = Veg),
                    linewidth = 1, linetype = "dashed", inherit.aes = FALSE) +
          geom_ribbon(data = adj_veg,
                      aes(x = year, ymin = adj_025, ymax = adj_975,
                          fill = Veg, group = Veg),
                      alpha = 0.10, color = NA, inherit.aes = FALSE) +
          labs(subtitle = "Solid = raw MESMA; dashed = bias-corrected (Olofsson 2014); band = \u00b11.96 SE")
        readr::write_csv(adj_veg,
                         file.path(OUT_DIR, paste0("inference_", file_prefix, "_bias_corrected.csv")))
        cat(sprintf("[BIAS CORRECTION] Saved bias-corrected time series: %s\n",
                    file.path(OUT_DIR, paste0("inference_", file_prefix, "_bias_corrected.csv"))))
      }
    }

    ggsave(file.path(OUT_DIR, paste0("inference_", file_prefix, "_normalized_timeseries.png")), p_ts, width = 8, height = 6)
    readr::write_csv(inf_veg, file.path(OUT_DIR, paste0("inference_", file_prefix, "_normalized_timeseries.csv")))
    cat(sprintf("Saved inference %s-normalized time series plot to: %s\n", method, file.path(OUT_DIR, paste0("inference_", file_prefix, "_normalized_timeseries.png"))))
  }

  if (is.null(inf_barren) || nrow(inf_barren) == 0) {
    cat(sprintf("[INFERENCE] No barren rows found in %s inference aggregation results.\n", method))
    return(invisible(NULL))
  }

  # 2. Barren cover
  p_barren <- ggplot(inf_barren, aes(x = year, y = global_coef)) +
    .ts_layers() +
    geom_line(color = "saddlebrown", linewidth = 1) +
    geom_point(color = "saddlebrown") +
    geom_ribbon(aes(ymin = coef_025, ymax = coef_975), alpha = 0.15, fill = "saddlebrown", color = NA) +
    labs(title = paste0(method, ": Barren Fraction"), x = "Year", y = "Barren Fraction") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
    scale_x_continuous(limits = c(1984, NA)) +
    theme_minimal()
  ggsave(file.path(OUT_DIR, paste0("inference_", file_prefix, "_barren_cover.png")), p_barren, width = 8, height = 6)
  readr::write_csv(inf_barren, file.path(OUT_DIR, paste0("inference_", file_prefix, "_barren_cover.csv")))
  cat(sprintf("Saved inference %s-based barren cover plot to: %s\n", method, file.path(OUT_DIR, paste0("inference_", file_prefix, "_barren_cover.png"))))

  # 4. Species-level plots
  if (include_species_plots) {
    species_data <- full_data[veg_norm %in% c("herbs","populus","tamarix"), ]
    if (nrow(species_data) > 0) {
      species_data$Veg <- dplyr::case_when(
        tolower(species_data$Veg) == "herbs" ~ "Herbs",
        tolower(species_data$Veg) == "populus" ~ "Populus",
        tolower(species_data$Veg) == "tamarix" ~ "Tamarix",
        TRUE ~ species_data$Veg
      )

      p_sp_ts <- ggplot(species_data, aes(x = year, y = global_coef, color = Veg, group = Veg)) +
        .ts_layers() +
        geom_line(linewidth = 1) +
        geom_point(show.legend = FALSE) +
        geom_ribbon(aes(ymin = coef_025, ymax = coef_975, fill = Veg), alpha = 0.12, color = NA) +
        scale_color_manual(values = SPECIES_COLORS) +
        scale_fill_manual(values = SPECIES_COLORS) +
           labs(title = paste0(method, " Species"),
             x = "Year", y = "Total Normalized Fraction", color = "Veg", fill = "Veg") +
        scale_x_continuous(limits = c(1984, NA)) +
        theme_minimal()
      ggsave(file.path(OUT_DIR, paste0("inference_", file_prefix, "_species_separate.png")), p_sp_ts, width = 8, height = 6)
      readr::write_csv(species_data, file.path(OUT_DIR, paste0("inference_", file_prefix, "_species_separate.csv")))
      cat(sprintf("Saved inference %s species plot to: %s\n", method, file.path(OUT_DIR, paste0("inference_", file_prefix, "_species_separate.png"))))

      # Species stacked proportion
      df_wide_sp <- tryCatch({ tidyr::pivot_wider(species_data |> dplyr::select(year, Veg, global_coef), names_from = Veg, values_from = global_coef) }, error = function(e) NULL)
      if (!is.null(df_wide_sp)) {
        sp_cols <- intersect(c("Herbs","Populus","Tamarix"), names(df_wide_sp))
        if (length(sp_cols) >= 2) {
          for (col in sp_cols) {
            df_wide_sp[[col]] <- as.numeric(df_wide_sp[[col]])
            df_wide_sp[[col]][is.na(df_wide_sp[[col]])] <- 0
          }
          df_wide_sp$total <- rowSums(df_wide_sp[, sp_cols, drop = FALSE], na.rm = TRUE)
          df_prop_sp <- df_wide_sp |> dplyr::filter(is.finite(total) & total > 0) |>
            dplyr::mutate(across(all_of(sp_cols), ~ . / total)) |>
            tidyr::pivot_longer(cols = all_of(sp_cols), names_to = "Veg", values_to = "prop")

          p_sp_stacked <- ggplot(df_prop_sp, aes(x = year, y = prop, fill = Veg)) +
            geom_area() +
            scale_fill_manual(values = SPECIES_COLORS) +
              labs(title = paste0(method, " Species Share"),
                 x = "Year", y = "Proportion", fill = "Veg") +
            scale_x_continuous(limits = c(1984, NA)) +
            theme_minimal()
          ggsave(file.path(OUT_DIR, paste0("inference_", file_prefix, "_species_stacked.png")), p_sp_stacked, width = 8, height = 6)
          readr::write_csv(df_prop_sp, file.path(OUT_DIR, paste0("inference_", file_prefix, "_species_stacked.csv")))
          cat(sprintf("Saved inference %s species stacked plot to: %s\n", method, file.path(OUT_DIR, paste0("inference_", file_prefix, "_species_stacked.png"))))
        }
      }
    }
  }

  invisible(NULL)
}

# RAW_BANDS defined in mesma_config.R
# NOTE: Keep default args literal (not RAW_BANDS) so this remains worker-safe in parallel futures.
# Blue band is always retained.
normalize_band_names <- function(df, bands = c("blue", "green", "red", "nir", "swir1", "swir2")) {
  if (is.null(df)) stop("normalize_band_names: df is NULL")
  if (nrow(df) == 0) return(df)
  current_names <- names(df)
  for (b in bands) {
    candidates <- c(b, toupper(b), tools::toTitleCase(b), paste0('band_', b), toupper(paste0('band_', b)), paste0('Band_', b))
    for (cand in candidates) {
      if (cand %in% current_names && !(b %in% current_names)) {
        names(df)[names(df) == cand] <- b
        current_names <- names(df)
        break
      }
    }
  }
  df
}

# Compute all supported indices from the current raw-band columns.
# This intentionally overwrites any pre-existing band-derived columns so they
# stay aligned with the latest band values (e.g. after sensor bias correction).
# IMPORTANT: This must be worker-safe (parallel futures).
# NOTE: Indices that depend on a soil-line slope (e.g. WDVI) will be computed only
# if a finite slope is available (either passed in or present as SOIL_LINE_SLOPE).
# L2-normalize a feature vector per observation (whole-vector)
#
# Input: vec with n_indices * n_bins values.
# Output: vec / ||vec||_2 (with NA treated as 0 for the norm).
#
# The extra parameters are ignored but retained for API stability.
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

# optional matrix version (unused here but kept for API parity)
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

# (rest of script continues...)
# with small sample sizes typical in endmember bundles (n=3-7)
compute_bundle_covariance <- function(Mv, verbose = FALSE) {
  # Mv: n x p matrix of endmember variant vectors (rows = variants, cols = features)
  # Returns: list(C = diagonal covariance, A = transformation matrix for sampling, mu = mean)

  Mv <- as.matrix(Mv)
  n <- nrow(Mv)
  p <- ncol(Mv)

  if (n < 2 || p < 1) {
    return(NULL)
  }

  # Compute mean
  mu <- colMeans(Mv, na.rm = TRUE)

  # Compute per-band variance (diagonal covariance)
  vars <- apply(Mv, 2, var, na.rm = TRUE)
  vars[!is.finite(vars) | vars < 1e-12] <- 1e-12

  # Diagonal covariance matrix
  C <- diag(vars, nrow = p)

  # Transformation matrix for sampling: A = diag(sqrt(vars))
  # For sampling: x = mu + A * z where z ~ N(0, I)
  A <- diag(sqrt(vars), nrow = p)

  if (verbose) {
    cat(sprintf("    [BUNDLE] n=%d variants, p=%d bands, var range=[%.4g, %.4g]\n",
                n, p, min(vars), max(vars)))
  }

  return(list(
    C = C,
    A = A,
    mu = mu,
    method = "diagonal"
  ))
}

safe_mul_vec <- function(a, b, allow_recycle = TRUE, caller = NULL) {
  la <- length(a); lb <- length(b)
  if (la == 0 || lb == 0) return(numeric(0))
  if (la == lb) return(a * b)
  if (la == 1) return(rep(a, lb) * b)
  if (lb == 1) return(a * rep(b, la))
  if (allow_recycle && (la %% lb == 0 || lb %% la == 0)) {
    if (la < lb) a <- rep(a, length.out = lb) else b <- rep(b, length.out = la)
    return(a * b)
  }
  caller_text <- if (is.null(caller)) "safe_mul_vec" else paste0(caller, ": ")
  stop(sprintf("%sIncompatible lengths for multiplication: %d vs %d", caller_text, la, lb))
}

safe_dot <- function(a, b, na.rm = TRUE) {
  if (length(a) == 0 || length(b) == 0) return(0)
  prod <- safe_mul_vec(a, b, allow_recycle = TRUE, caller = "safe_dot")
  if (na.rm) sum(prod, na.rm = TRUE) else sum(prod)
}

safe_col_weighted_avg <- function(mat, wts) {
  if (is.null(mat) || nrow(mat) == 0) return(rep(0, ifelse(is.null(mat), 1, ncol(mat))))
  n <- nrow(mat)
  if (length(wts) == 0) wts <- rep(1, n)
  if (length(wts) != n) {
    if (length(wts) == 1 || n %% length(wts) == 0 || length(wts) %% n == 0) {
      wts <- rep(wts, length.out = n)
      warning(sprintf("safe_col_weighted_avg: adjusted weights vector to length %d", n))
    } else {
      stop(sprintf("safe_col_weighted_avg: weights length (%d) incompatible with rows in mat (%d)", length(wts), n))
    }
  }
  if (sum(wts, na.rm = TRUE) == 0) {
    return(as.numeric(colMeans(mat, na.rm = TRUE)))
  }
  wts <- as.numeric(wts) / sum(wts, na.rm = TRUE)
  as.numeric(colSums(mat * wts, na.rm = TRUE))
}

# OUTPUT_DIR and OUT_DIR defined in mesma_config.R

if (!file.exists(INPUT_CSV)) {
  stop(paste0("Required input CSV not found: ", INPUT_CSV))
}
raw_df <- tryCatch(
  {
    readr::read_csv(INPUT_CSV, show_col_types = FALSE)
  },
  error = function(e) stop(paste0("Failed to read INPUT_CSV: ", e$message))
)
if (nrow(raw_df) < 2) stop("INPUT_CSV contains fewer than 2 rows")
if (!"location_id" %in% names(raw_df) || (!"prediction_date" %in% names(raw_df) && !"date" %in% names(raw_df))) {
  stop("INPUT_CSV must contain 'location_id' and 'prediction_date' (or 'date') columns")
}

date_col <- if ("prediction_date" %in% names(raw_df)) "prediction_date" else "date"

# Canonicalize vegetation labels (typos, herbs grouping, agriculture aliases)
raw_df <- canonicalize_veg_labels(raw_df)

# === GEE INPUT MAPPING BLOCK ===
# Map common GEE-exported column names to the names expected by this script
# 1) Map 'vegetation' -> 'Veg'
if ("vegetation" %in% names(raw_df) && !"Veg" %in% names(raw_df)) {
  raw_df$Veg <- raw_df$vegetation
}
# 3) Ensure 'location_id' is character (GEE sometimes exports numeric)
if ("location_id" %in% names(raw_df) && !is.character(raw_df$location_id)) {
  raw_df$location_id <- as.character(raw_df$location_id)
}
# Ensure band columns like Blue/Green -> lower-case handled by normalize_band_names later
# === END GEE INPUT MAPPING BLOCK ===

df <- raw_df

# --- drop observations with phenology year before 1984 ----------------------
# applying the cutoff here ensures downstream scripts (january_averages,
# fit_mesma, etc.) never see the unwanted years.  We compute a temporary
# pheno_year if it is not already present so that filtering works on raw
# input CSVs as well.
cutoff <- 1984
if ("date" %in% names(df) || "pheno_year" %in% names(df)) {
  ph <- NULL
  if ("pheno_year" %in% names(df)) {
    ph <- df$pheno_year
  } else {
    ph <- ifelse(lubridate::month(df$date) >= 3,
                 lubridate::year(df$date),
                 lubridate::year(df$date) - 1)
  }
  n_before <- nrow(df)
  keep <- is.na(ph) | ph >= cutoff
  df <- df[keep, , drop = FALSE]
  if (n_before != nrow(df)) {
    cat(sprintf("[FILTER] dropped %d rows with pheno_year < %d\n", n_before - nrow(df), cutoff))
  }
}

# Preserve zenith.angle if present (e.g. from metadata), otherwise allow recalculation
if ("zenith.angle" %in% names(df)) {
  cat("[NOTICE] Preserving existing 'zenith.angle' in input data.\n")
} else {
  df$zenith.angle <- NA_real_
}

# Remove large outliers robustly per (location_id, pheno_year) where possible, otherwise per-location.
# Uses spline-based outlier detection for groups with sufficient data, otherwise falls back to MAD.
remove_large_outliers <- function(df, candidates = NULL, mad_thresh = OUTLIER_MAD_THRESHOLD) {
  if (!isTRUE(ENABLE_OUTLIER_REMOVAL)) return(df)

  # determine global interpolation method (same logic as build_pentad_matrix)
  interp_method <- NULL
  if (exists("MESMA_PARAMS") && !is.null(MESMA_PARAMS$interpolate_inference)) {
    interp_method <- MESMA_PARAMS$interpolate_inference
  } else if (exists("INTERPOLATE_INFERENCE")) {
    interp_method <- INTERPOLATE_INFERENCE
  }
  interp_method <- get_interpolate_method(interp_method)

  if (is.null(candidates)) {
    # Use indices we know are meaningful: OPTIMAL_INDICES + RAW_BANDS, if present
    candidates <- intersect(unique(c(OPTIMAL_INDICES, RAW_BANDS)), names(df))
  } else {
    candidates <- intersect(candidates, names(df))
  }
  if (length(candidates) == 0) {
    cat("[OUTLIER] No candidate indices found for outlier detection; skipping\n")
    return(df)
  }
  if (!"location_id" %in% names(df)) {
    cat("[OUTLIER] 'location_id' missing from data; skipping outlier removal\n")
    return(df)
  }

  # Assign pheno_year if possible for better grouping
  if (!"pheno_year" %in% names(df) && "date" %in% names(df)) {
    df$pheno_year <- assign_pheno_year(df$date)
  }

  grp <- interaction(df$location_id, ifelse(is.na(df$pheno_year), "NA", as.character(df$pheno_year)), drop = TRUE)
  groups <- split(seq_len(nrow(df)), grp)
  removed_idx <- logical(nrow(df))
  n_groups <- length(groups)

  # if Whittaker smoothing is active, replace spline detection with a simple MAD
  # threshold per group/column (no date requirement)
  if (identical(interp_method, "whittaker")) {
    for (g in seq_along(groups)) {
      rows <- groups[[g]]
      sub <- df[rows, , drop = FALSE]
      if (length(rows) < 5) {
        removed_idx[rows] <- TRUE
        next
      }
      for (col in candidates) {
        if (!is.numeric(sub[[col]])) next
        v <- sub[[col]]
        finite_idx <- is.finite(v)
        if (sum(finite_idx) < 3) next
        medv <- stats::median(v[finite_idx], na.rm = TRUE)
        madv <- stats::mad(v[finite_idx], na.rm = TRUE)
        if (!is.finite(madv) || madv <= 0) next
        mask <- rep(FALSE, length(v))
        mask[finite_idx] <- abs(v[finite_idx] - medv) > mad_thresh * madv
        removed_idx[rows[mask]] <- TRUE
      }
    }

    if (any(removed_idx, na.rm = TRUE)) {
      n_removed <- sum(removed_idx, na.rm = TRUE)
      cat(sprintf("[OUTLIER] Removed %d observations across %d groups (MAD-only due to whittaker)\n",
                  n_removed, n_groups))
      df <- df[!removed_idx, , drop = FALSE]
    }
    return(df)
  }

  # --- original spline-based algorithm (for linear/none) --------------------
  for (g in seq_along(groups)) {
    rows <- groups[[g]]
    sub <- df[rows, , drop = FALSE]
    # Remove location-years with fewer than 5 observations entirely
    if (length(rows) < 5) {
      removed_idx[rows] <- TRUE
      next
    }
    out_mask <- rep(FALSE, nrow(sub))

    # Check if we have date for spline - require at least 10 observations
    has_date <- "date" %in% names(sub) && any(!is.na(sub$date))
    if (!has_date || length(rows) < 10) next  # skip outlier removal if spline fitting is not possible

    # Compute DOY
    sub$doy <- as.numeric(format(sub$date, "%j"))
    for (col in candidates) {
      if (!is.numeric(sub[[col]])) next
      colv <- sub[[col]]
      finite_idx <- is.finite(colv) & is.finite(sub$doy)
      if (sum(finite_idx) < 5) next  # not enough for spline
      tryCatch({
        # --- Iterative Spline Fitting for Robustness ---
        # Pass 1: Initial fit to identify gross outliers
        x <- sub$doy[finite_idx]
        y <- colv[finite_idx]

        n_unique <- length(unique(x))
        fit1 <- stats::smooth.spline(x, y, df = min(OUTLIER_SPLINE_MAX_DF, length(x)/2, n_unique - 1))
        pred1 <- predict(fit1, x)$y
        res1 <- y - pred1
        mad1 <- stats::mad(res1, na.rm = TRUE)

        if (!is.finite(mad1) || mad1 <= 1e-6) next

        # Temporarily exclude gross outliers for the second pass (1.5x threshold)
        keep_mask <- abs(res1 - stats::median(res1, na.rm = TRUE)) <= (mad_thresh * 1.5 * mad1)

        # Pass 2: Refit on cleaner data if we have enough points left
        if (sum(keep_mask) >= 5) {
          n_unique2 <- length(unique(x[keep_mask]))
          fit2 <- stats::smooth.spline(x[keep_mask], y[keep_mask], df = min(OUTLIER_SPLINE_MAX_DF, sum(keep_mask)/2, n_unique2 - 1))
          # Predict against the refined trend for ALL points
          pred_final <- predict(fit2, x)$y
        } else {
          # Fallback to first pass if too many points removed
          pred_final <- pred1
        }

        # Final outlier detection against the robust trend
        residuals <- y - pred_final
        med_res <- stats::median(residuals, na.rm = TRUE)
        mad_res <- stats::mad(residuals, na.rm = TRUE)

        if (!is.finite(mad_res) || mad_res <= 0) next

        # Flag outliers
        this_mask <- rep(FALSE, length(colv))
        this_mask[finite_idx] <- abs(residuals - med_res) > mad_thresh * mad_res
        out_mask <- out_mask | this_mask
      }, error = function(e) {
        # Skip this column if spline fitting fails
      })
    }

    # Mark for removal
    if (any(out_mask, na.rm = TRUE)) {
      removed_idx[rows[which(out_mask)]] <- TRUE
    }
  }

  if (any(removed_idx, na.rm = TRUE)) {
    n_removed <- sum(removed_idx, na.rm = TRUE)
    cat(sprintf("[OUTLIER] Removed %d observations across %d groups\n", n_removed, n_groups))
    df <- df[!removed_idx, , drop = FALSE]
  }

  df
}

df <- normalize_band_names(df)

if (!"date" %in% names(df) && "prediction_date" %in% names(df)) df$date <- as.Date(df$prediction_date)
if ("date" %in% names(df)) {
  df$date <- as.Date(df$date)
}

# =============================================================================
# SENSOR BIAS CORRECTION — applied in-place before any index computation,
# outlier removal, or PPI baseline derivation.
# Harmonises LANDSAT_89 (OLI) raw bands to the ETM+ radiometric scale using
# per-band affine coefficients (slope + intercept) from satellite_bias_check.R:
#   ETM+ ≈ slope * OLI + intercept
# =============================================================================
{
  df <- apply_oli_etm_bias_correction(df, dataset_label = "training", log_prefix = "[BIAS CORR]")
}
# =============================================================================

# Recompute all band-derived indices from the bias-corrected raw bands so early
# contamination filtering and downstream normalization use the corrected values.
before_cols <- names(df)
df <- compute_indices_from_bands(df)

new_cols <- setdiff(names(df), before_cols)
if (length(new_cols) > 0) {
  cat(sprintf("[NOTICE] Computed indices from raw bands before early filtering: %s\n", paste(new_cols, collapse=", ")))
}

# === CRITICAL: Filter out dust contamination BEFORE year filtering and PPI baseline calculation ===
# This ensures the PPI baseline is computed from clean observations only
if ("NDDI" %in% names(df)) {
  dust_count <- sum(df$NDDI > NDDI_DUST_THRESHOLD, na.rm = TRUE)
  total_before <- nrow(df)
  df <- df[!(df$NDDI > NDDI_DUST_THRESHOLD), , drop = FALSE]
  total_after <- nrow(df)
  filtered <- total_before - total_after
  cat(sprintf("[EARLY FILTERING] Filtered out %d observations with dust (NDDI > %s) contamination\n", filtered, .nddi_thresh_fmt()))
  cat(sprintf("[EARLY FILTERING] Dataset after contamination filtering: %d rows from %d locations\n",
              total_after, length(unique(df$location_id))))

  # Additionally remove extreme outliers (robust MAD-based), operating per location-year where possible
  df <- remove_large_outliers(df)

  # Defer soil line calculation until after all filtering is complete
  cat("[SOIL LINE] Soil line computation deferred until after final filtering and normalization; it will be calculated later before index computation.\n")

  # Print average images per location for each year, excluding years with < 10 observations
  if ("date" %in% names(df)) {
    year_stats <- df %>%
      dplyr::mutate(.year = lubridate::year(date)) %>%
      dplyr::group_by(.year) %>%
      dplyr::summarise(total_images = dplyr::n(),
                       n_locations = dplyr::n_distinct(location_id),
                       avg_images_per_location = total_images / n_locations,
                       .groups = "drop") %>%
      dplyr::filter(total_images >= 10)
    if (nrow(year_stats) == 0) {
      cat("[DATA STATS] No years with >= 10 observations to summarize.\n")
    } else {
      # only create a diagnostic plot (no table printed)
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        p_years <- ggplot2::ggplot(year_stats, ggplot2::aes(x = .year, y = avg_images_per_location)) +
          ggplot2::geom_line() + ggplot2::geom_point() +
          ggplot2::theme_minimal() +
          ggplot2::labs(title = "Images per Location by Year",
                        x = "Year", y = "Avg images / location") +
          ggplot2::scale_x_continuous(limits = c(1984, NA))
        print(p_years)
        if (exists("OUTPUT_DIR") && !is.null(OUTPUT_DIR)) {
          try(
            ggplot2::ggsave(file.path(OUTPUT_DIR, "avg_images_per_location.png"), p_years, width = 6, height = 4),
            silent = TRUE)
        }
      }
    }
  } else if ("year" %in% names(df)) {
    year_stats <- df %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(total_images = dplyr::n(),
                       n_locations = dplyr::n_distinct(location_id),
                       avg_images_per_location = total_images / n_locations,
                       .groups = "drop") %>%
      dplyr::filter(total_images >= 10)
    if (nrow(year_stats) == 0) {
      cat("[DATA STATS] No years with >= 10 observations to summarize (using 'year' column).\n")
    } else {
      # only create a diagnostic plot (no table printed)
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        p_years <- ggplot2::ggplot(year_stats, ggplot2::aes(x = year, y = avg_images_per_location)) +
          ggplot2::geom_line() + ggplot2::geom_point() +
          ggplot2::theme_minimal() +
          ggplot2::labs(title = "Images per Location by Year",
                        x = "Year", y = "Avg images / location") +
          ggplot2::scale_x_continuous(limits = c(1984, NA))
        print(p_years)
        if (exists("OUTPUT_DIR") && !is.null(OUTPUT_DIR)) {
          try(
            ggplot2::ggsave(file.path(OUTPUT_DIR, "avg_images_per_location.png"), p_years, width = 6, height = 4),
            silent = TRUE)
        }
      }
    }
  }

} else {
  cat("[WARNING] NDDI not found in data; skipping early contamination filtering\n")
}
# ========================================================================================================

# Robustly prune collinear features based on correlation matrix
prune_collinear_features <- function(df, features, threshold = 0.95) {
  if (is.null(features) || length(features) < 2) return(features)
  
  common_feats <- intersect(names(df), features)
  if (length(common_feats) < 2) return(common_feats)
  
  cat(sprintf("[FEATURE PRUNE] Checking %d features for collinearity > %.3f...\n", length(common_feats), threshold))
  
  # Compute correlation matrix
  # Use spearman for robustness to non-linearity? or pearson. Pearson is standard for strict linear redundancy.
  cm <- cor(df[, common_feats], use = "pairwise.complete.obs", method = "pearson")
  
  if (any(is.na(cm))) {
    cm[is.na(cm)] <- 0
  }
  
  diag(cm) <- 0 # Ignore self-correlation
  
  dropped <- character(0)
  
  # Greedy removal: prioritize keeping features earlier in the list (assuming user preference order)
  for (i in seq_along(common_feats)) {
    f1 <- common_feats[i]
    if (f1 %in% dropped) next
    
    for (j in (i + 1):length(common_feats)) {
      if (j > length(common_feats)) break
      f2 <- common_feats[j]
      if (f2 %in% dropped) next
      
      score <- abs(cm[f1, f2])
      if (score > threshold) {
        cat(sprintf("  Dropping '%s' (corr %.3f with '%s')\n", f2, score, f1))
        dropped <- c(dropped, f2)
      }
    }
  }
  
  kept <- setdiff(common_feats, dropped)
  cat(sprintf("[FEATURE PRUNE] Kept %d features, dropped %d.\n", length(kept), length(dropped)))
  return(kept)
}

normalize_mesma_data <- function(df, cols = unique(c(OPTIMAL_INDICES, RAW_BANDS)), lat_default = 40.2) {
  cat("Applying comprehensive MESMA data normalization...\n")
  
  if (!"zenith.angle" %in% names(df)) df$zenith.angle <- NA_real_
  if (!"PPI" %in% names(df)) stop("[PPI] Missing required column 'PPI' in input data")
  if (all(!is.finite(df$PPI))) stop("[PPI] All PPI values are non-finite; refusing to continue")
  
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
      if (is.list(params) && all(c("mean", "sd") %in% names(params))) {
        mu <- params$mean
        sigma <- params$sd
        if (is.finite(sigma) && sigma > 0) {
          df[[col]] <- (df[[col]] - mu) / sigma
        }
      }
    }
  }
  
  list(
    df = df,
    INDEX_SCALES = INDEX_SCALES
  )
}

apply_stored_normalization <- function(df, norm_params, cols = unique(c(OPTIMAL_INDICES, RAW_BANDS)), lat_default = 40.2) {
  cat("Applying stored normalization parameters to data...\n")
  
  if (!"zenith.angle" %in% names(df)) df$zenith.angle <- NA_real_
  if (is.null(norm_params$INDEX_SCALES) || length(norm_params$INDEX_SCALES) == 0) {
    stop("[INDEX_SCALES] Missing stored normalization parameters (INDEX_SCALES)")
  }

  # Ensure PPI is present with at least some finite values
  if (!"PPI" %in% names(df) || all(!is.finite(df$PPI))) {
    if (all(c("nir", "red") %in% names(df)) && !"DVI" %in% names(df)) {
      df$DVI <- as.numeric(df$nir) - as.numeric(df$red)
    }
    if (exists("add_ppi_columns")) {
      dvi_soil_vec <- compute_dvi_soil_per_location(df)
      df <- add_ppi_columns(df, dvi_soil = dvi_soil_vec)
      cat("[PPI] Added PPI to data before applying stored normalization (per-location dvi_soil + per-location M).\n")
    }
  }
  if (!"PPI" %in% names(df) || all(!is.finite(df$PPI))) {
    stop("[PPI] Missing or all PPI values non-finite after attempted auto-add; refusing to continue")
  }
  if (!"PPI_raw" %in% names(df)) df$PPI_raw <- df$PPI
  
  # Apply stored INDEX_SCALES (mean/sd) from training normalization to df
  cat(sprintf("[apply_stored_normalization] Applying INDEX_SCALES to %d indices\n", length(norm_params$INDEX_SCALES)))
  for (col in names(norm_params$INDEX_SCALES)) {
    if (col %in% names(df)) {
      params <- norm_params$INDEX_SCALES[[col]]
      if (!is.list(params) || !all(c("mean", "sd") %in% names(params))) {
        stop(sprintf("[INDEX_SCALES] Invalid params for index '%s' (expected list(mean, sd))", col))
      }
      mu <- params$mean
      sigma <- params$sd
      if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0) {
        stop(sprintf("[INDEX_SCALES] Non-finite or non-positive mean/sd for index '%s'", col))
      }
      df[[col]] <- (df[[col]] - mu) / sigma
    }
  }
  
  df
}

apply_secondary_normalization <- function(vec, idx_names, norm_params) {
  if (is.null(norm_params$INDEX_SCALES) || length(norm_params$INDEX_SCALES) == 0) {
    stop("[INDEX_SCALES] Missing INDEX_SCALES in normalization params")
  }
  
  if (length(vec) != length(idx_names)) {
    stop("[INDEX_SCALES] apply_secondary_normalization: vec/idx_names length mismatch")
  }
  
  for (i in seq_along(idx_names)) {
    idx_name <- idx_names[i]
    if (!idx_name %in% names(norm_params$INDEX_SCALES)) {
      stop(sprintf("[INDEX_SCALES] Missing normalization params for index '%s'", idx_name))
    }
    params <- norm_params$INDEX_SCALES[[idx_name]]
    if (!is.list(params) || !all(c("mean", "sd") %in% names(params))) {
      stop(sprintf("[INDEX_SCALES] Invalid params format for index '%s'", idx_name))
    }
    mu <- params$mean
    sigma <- params$sd
    if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0) {
      stop(sprintf("[INDEX_SCALES] Non-finite or non-positive mean/sd for index '%s'", idx_name))
    }
    vec[i] <- (vec[i] - mu) / sigma
  }
  
  vec
}

## No external GeoJSON: construct location map directly from CSV lat/lon
if (all(c("lon", "lat") %in% names(df))) {
  df$location_id <- make_location_id(df$lon, df$lat)
  # Ensure mapping columns exist (fill with NA if missing)
  if (!"Veg" %in% names(df)) df$Veg <- NA_character_
  # Build a minimal gpts_map from unique lat/lon combos in the CSV
  # Use per-location aggregation: pick the first non-missing lat/lon and the
  # first non-missing Veg for that location (avoids losing Veg when
  # the first row happens to have NA)
  gpts_map <- df |>
    dplyr::group_by(location_id) |>
    dplyr::summarise(
      lat = if ("lat" %in% names(df)) { v <- na.omit(lat); if (length(v)>0) v[1] else NA_real_ } else NA_real_,
      lon = if ("lon" %in% names(df)) { v <- na.omit(lon); if (length(v)>0) v[1] else NA_real_ } else NA_real_,
      Veg = if ("Veg" %in% names(df)) { v <- na.omit(as.character(Veg)); if (length(v)>0) tolower(v[1]) else NA_character_ } else NA_character_,
      .groups = "drop"
    )
  gpts_map$location_row <- as.character(seq_len(nrow(gpts_map)))
  gpts_map$location_id_seq <- gpts_map$location_row
  cat(sprintf("[NOTICE] Constructed gpts_map from %d unique lat/lon combos in CSV\n", nrow(gpts_map)))
} else {
  cat("[NOTICE] No 'lat'/'lon' columns found in CSV; GeoJSON mapping disabled.\n")
  gpts_map <- data.frame(location_id = character(0), location_row = character(0), Veg = character(0), stringsAsFactors = FALSE)
}

if (nrow(gpts_map) == 0) {
  cat("[NOTICE] No location mapping points found; continuing without derived soil baselines.\n")
}

if ("location_id" %in% names(df) && "location_id" %in% names(gpts_map)) {
  if (!is.character(df$location_id)) df$location_id <- as.character(df$location_id)
  if (!is.character(gpts_map$location_id)) gpts_map$location_id <- as.character(gpts_map$location_id)

  # Normalize IDs to row-number format so training CSV `location_id` (which uses
  # row numbers) will match GeoJSON row IDs (`location_id_seq`).
  df$location_id <- trimws(df$location_id)
  df$location_id[df$location_id == ""] <- NA_character_
  gpts_map$location_id <- trimws(gpts_map$location_id)
  gpts_map$location_id[gpts_map$location_id == ""] <- NA_character_

  # If CSV has 'L_123' style IDs, strip the 'L_' and coerce numeric IDs to
  # canonical integer string form (e.g., '156'). This ensures '156' matches
  # the GeoJSON row id '156'.
  df$location_id <- ifelse(grepl('^L_?[0-9]+$', toupper(df$location_id)), sub('^L_', '', toupper(df$location_id)), df$location_id)
  if (any(grepl('^[0-9]+$', df$location_id, perl=TRUE), na.rm = TRUE)) {
    df$location_id[grepl('^[0-9]+$', df$location_id)] <- as.character(as.integer(df$location_id[grepl('^[0-9]+$', df$location_id)]))
  }
  # Ensure GeoJSON ids are canonical integer strings as well
  if (any(grepl('^[0-9]+$', gpts_map$location_id, perl=TRUE), na.rm = TRUE)) {
    gpts_map$location_id[grepl('^[0-9]+$', gpts_map$location_id)] <- as.character(as.integer(gpts_map$location_id[grepl('^[0-9]+$', gpts_map$location_id)]))
  }
  cat("[NOTICE] Normalized CSV 'location_id' to row-number format where possible.\n")

  if (!"Veg" %in% names(df)) df$Veg <- NA_character_
  pre_non_na <- sum(!is.na(df$Veg) & df$Veg != "")

  joined <- join_and_fill_veg(df, gpts_map)

  post_non_na <- sum(!is.na(joined$Veg) & joined$Veg != "")

  # The row-number mapping fallback is not used here. We
  # only ever join on the canonical "location_id" field and do not attempt any
  # secondary join strategies.  This simplifies behaviour and avoids confusing
  # implicit matches when IDs are missing or misaligned.

  if ("lat" %in% names(gpts_map)) {
    if ("lat.geo" %in% names(joined)) {
      joined$lat <- joined$lat.geo
      joined$lat.geo <- NULL
      cat("[NOTICE] Replaced CSV latitude with GeoJSON latitude where available\n")
    }
  }
  df <- joined

  matched_locs <- length(intersect(na.omit(unique(as.character(df$location_id))), na.omit(unique(as.character(gpts_map$location_id)))))
  cat(sprintf("[NOTICE] GeoJSON join results - Veg before=%d after=%d; matched location_id strings=%d\n", pre_non_na, post_non_na, matched_locs))

  # === TRAINING DATA DIAGNOSTIC ===
  cat("\n=== TRAINING DATA DIAGNOSTIC ===\n")
  cat(sprintf("Total rows in joined data: %d\n", nrow(df)))
  cat(sprintf("Rows with Veg='barren': %d\n", sum(tolower(df$Veg) == "barren", na.rm = TRUE)))
  cat(sprintf("Rows with non-missing Veg: %d\n", sum(!is.na(df$Veg) & df$Veg != "")))

  cat("\nVegetation counts by type:\n")
  veg_types <- sort(unique(tolower(na.omit(df$Veg))))
  for (vt in veg_types) {
    cat(sprintf("  %s: %d rows\n", vt, sum(tolower(df$Veg) == vt, na.rm = TRUE)))
  }
  cat("=========================================\n\n")

}

# Dust contamination (NDDI > NDDI_DUST_THRESHOLD; default 0.18) already filtered in [EARLY FILTERING] above.
# Proceed with outlier removal, soil line estimation, and index recomputation.
df <- remove_large_outliers(df)
if (!exists("compute_soil_line_slope", mode = "function")) {
  stop("[SOIL LINE] Required helper 'compute_soil_line_slope' not available. Ensure mesma_helpers.R is loaded.")
}
compute_soil_line_slope(df)
df <- compute_indices_from_bands(df)

# STEP 2: Recalculate PPI from the bias-corrected raw bands, now that Veg is available
cat("[NOTICE] Retaining all years for PPI baseline calculation and trend analysis. Training subset will be filtered later.\n")
if (exists("add_ppi_columns")) {
  if ("PPI_raw" %in% names(df)) df$PPI_raw <- NULL
  if ("ppi_norm" %in% names(df)) df$ppi_norm <- NULL
  dvi_soil_vec <- compute_dvi_soil_per_location(df)
  df <- add_ppi_columns(df, dvi_soil = dvi_soil_vec)
  cat("[PPI] Recomputed PPI from bias-corrected raw bands before normalization (per-location dvi_soil + per-location M).\n")
} else {
  stop("[PPI] add_ppi_columns not available; cannot compute PPI")
}

if (!"PPI" %in% names(df) || all(!is.finite(df$PPI))) {
  stop("[PPI] Training data is missing PPI or all PPI values are non-finite after PPI computation")
}
n_ppi_na <- sum(!is.finite(df$PPI))
if (n_ppi_na > 0) {
  cat(sprintf("[PPI] %d of %d rows have non-finite PPI (missing DVI or edge-case inputs); these will be excluded from PPI-dependent modelling\n", n_ppi_na, nrow(df)))
}

# STEP 3: Apply normalization to all indices, including the new PPI
cat("\n=== APPLYING TRAINING DATA NORMALIZATION ===\n")

# --- FEATURE PRUNING (Optional) ---
if (exists("ENABLE_FEATURE_PRUNING") && isTRUE(ENABLE_FEATURE_PRUNING)) {
   thresh <- if(exists("FEATURE_PRUNING_THRESHOLD")) FEATURE_PRUNING_THRESHOLD else 0.95
   OPTIMAL_INDICES <- prune_collinear_features(df, OPTIMAL_INDICES, threshold = thresh)
}

# Backup raw values before normalization (z-scoring)
if ("PPI" %in% names(df) && !"PPI_raw" %in% names(df)) df$PPI_raw <- df$PPI
if ("MSAVI" %in% names(df)) { df$MSAVI_raw <- df$MSAVI; cat("[NOTICE] Backed up raw MSAVI values to 'MSAVI_raw' before normalization.\n") }
if ("nir" %in% names(df)) df$nir_raw <- df$nir
if ("red" %in% names(df)) df$red_raw <- df$red

norm_result <- normalize_mesma_data(df, cols = unique(c(OPTIMAL_INDICES, RAW_BANDS, "PPI")), lat_default = 40.2)
df <- norm_result$df
INDEX_SCALES <- norm_result$INDEX_SCALES

df <- backup_and_normalize_ppi(df)

# Filter to only include selected vegetation types
selected_vegs <- c("herbs", "populus", "tamarix", "barren")
df <- df[tolower(df$Veg) %in% selected_vegs, ]
cat(sprintf("Filtered training data to selected vegetation types: %s\n", paste(selected_vegs, collapse = ", ")))

# Pre-filter: remove vegetation observations whose feature signature is too similar to barren
# Uses the same feature space (indices + raw bands) that the endmember library will use,
# so barren-like observations are removed before endmember building.
barren_rows <- tolower(df$Veg) == "barren"
prefilter_features <- intersect(c(OPTIMAL_INDICES, RAW_BANDS), names(df))
if (any(barren_rows) && length(prefilter_features) > 0) {
  barren_mat <- as.matrix(df[barren_rows, prefilter_features, drop = FALSE])
  barren_mat[!is.finite(barren_mat)] <- NA
  barren_mean <- colMeans(barren_mat, na.rm = TRUE)
  barren_mean_norm <- barren_mean / sqrt(sum(barren_mean^2))

  if (all(is.finite(barren_mean_norm))) {
    veg_rows_idx <- which(!barren_rows)
    veg_mat <- as.matrix(df[veg_rows_idx, prefilter_features, drop = FALSE])
    veg_mat[!is.finite(veg_mat)] <- 0

    veg_norms <- sqrt(rowSums(veg_mat^2))
    veg_mat_normed <- veg_mat / ifelse(veg_norms > 0, veg_norms, 1)

    sims <- as.numeric(veg_mat_normed %*% barren_mean_norm)
    too_similar <- sims > BARREN_SIM_THRESHOLD

    if (any(too_similar)) {
      remove_idx <- veg_rows_idx[too_similar]
      for (v in unique(df$Veg[remove_idx])) {
        n_rem <- sum(df$Veg[remove_idx] == v)
        cat(sprintf("[BARREN PRE-FILTER] Removed %d/%d %s observations (cosine similarity to barren > %.2f)\n",
                    n_rem, sum(df$Veg == v), v, BARREN_SIM_THRESHOLD))
      }
      df <- df[-remove_idx, , drop = FALSE]
    } else {
      cat("[BARREN PRE-FILTER] No vegetation observations exceeded barren similarity threshold\n")
    }
  }
}

cat(sprintf("Remaining samples: %d\n", nrow(df)))

TRAINING_NORM_PARAMS <- list(
  INDEX_SCALES = INDEX_SCALES,
  INDEX_SCALES_SECONDARY = list()  # Will be populated after location mapping
)
cat(sprintf("Stored normalization params: INDEX_SCALES for %d indices\n",
            length(INDEX_SCALES)))
cat("(Secondary normalization will be computed after location mapping)\n")
cat("===========================================\n\n")

# Save preprocessed data
saveRDS(df, file = "preprocessed_data.rds")
saveRDS(gpts_map, file = "gpts_map.rds")
saveRDS(TRAINING_NORM_PARAMS, file = "training_norm_params.rds")

cat("Preprocessed data saved to preprocessed_data.rds, gpts_map.rds, and training_norm_params.rds\n")