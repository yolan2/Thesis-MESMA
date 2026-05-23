
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

# Compute a per-location DVI soil baseline for PPI.
# Baseline is defined as the median of the lowest `quantile_p` fraction of DVI
# observations within each location. This is deterministic, per-location,
# and avoids any constant/default soil baseline.
compute_dvi_soil_per_location <- function(df, quantile_p = 0.10, min_samples = 5L) {
  stopifnot(!is.null(df), nrow(df) > 0L, "location_id" %in% names(df))
  if (!"DVI" %in% names(df) && all(c("nir","red") %in% names(df))) {
    df <- df %>% mutate(DVI = as.numeric(nir) - as.numeric(red))
  }
  stopifnot("DVI" %in% names(df))

  soil_df <- df %>%
    filter(is.finite(DVI)) %>%
    group_by(location_id) %>%
    filter(n() >= min_samples) %>%
    summarise(
      dvi_soil = {
        vals <- DVI
        q <- suppressWarnings(quantile(vals, probs = quantile_p,
                                       na.rm = TRUE, names = FALSE, type = 7))
        median(vals[vals <= q], na.rm = TRUE)
      },
      .groups = "drop"
    ) %>%
    filter(is.finite(dvi_soil))

  out <- df %>% left_join(soil_df, by = "location_id") %>% pull(dvi_soil)
  if (any(is.na(out) & is.finite(df$DVI))) {
    bad_locs <- unique(as.character(df$location_id[is.na(out) & is.finite(df$DVI)]))
    stop(sprintf("[PPI] Cannot compute per-location dvi_soil for %d rows across %d locations (example locs: %s)",
                 sum(is.na(out) & is.finite(df$DVI)), length(bad_locs),
                 paste(head(bad_locs, 10), collapse = ", ")))
  }
  out
}

mesma_prepare_feature_columns <- function(df, required_cols, optional_cols = NULL, context = "data") {
  if (is.null(df)) stop(sprintf("[%s] df is NULL", context))
  if (!is.data.frame(df)) stop(sprintf("[%s] df is not a data.frame", context))

  required_cols <- unique(as.character(required_cols))
  required_cols <- required_cols[!is.na(required_cols) & nzchar(required_cols)]

  optional_cols <- if (is.null(optional_cols)) character(0) else unique(as.character(optional_cols))
  optional_cols <- optional_cols[!is.na(optional_cols) & nzchar(optional_cols)]

  missing_required <- setdiff(required_cols, names(df))
  if (length(missing_required) > 0) {
    cat(sprintf("\n[%s FEATURE SPACE CHECK]\n", toupper(context)))
    cat(sprintf("  Required columns (%d): %s\n", length(required_cols), paste(required_cols, collapse = ", ")))
    cat(sprintf("  Missing required (%d): %s\n", length(missing_required), paste(missing_required, collapse = ", ")))
    stop(sprintf("%s is missing %d required column(s).", context, length(missing_required)))
  }

  # Ensure optional columns exist (filled with NA_real_)
  if (length(optional_cols) > 0) {
    for (col in optional_cols) {
      if (!col %in% names(df)) df[[col]] <- NA_real_
    }
  }

  df
}

# --- Source shared config and helpers ---
# All tunable parameters are defined in mesma_config.R (single source of truth).
# Override any parameter AFTER the source() call below if needed for this run.
stopifnot(file.exists("mesma_config.R"))
source("mesma_config.R")
if (exists("INTERPOLATE_INFERENCE")) {
  if (!tolower(as.character(INTERPOLATE_INFERENCE)) %in%
      c("linear","whittaker","none","true","false","")) {
    stop(sprintf("[CONFIG] INTERPOLATE_INFERENCE has invalid value '%s'", INTERPOLATE_INFERENCE))
  }
}

stopifnot(file.exists("ppi_helpers.R"))
source("ppi_helpers.R")

stopifnot(file.exists("mesma_helpers.R"))
source("mesma_helpers.R")
set_mesma_seed()


# --- Shared plot constants ---

SPECIES_COLORS <- c("Herbs" = "#9ACD32", "Populus" = "#006400", "Tamarix" = "#D95F02")

# --- Shared utility: normalize vegetation names (lowercase + trim) ---

# --- Shared utility: compute pairwise Haversine distance matrix (km) ---
# vectorised version avoids explicit loops for speed
compute_haversine_distance_matrix <- function(coords) {
  rad <- pi / 180
  lat <- coords[,2] * rad
  lon <- coords[,1] * rad
  dlat <- outer(lat, lat, "-")
  dlon <- outer(lon, lon, "-")
  a <- sin(dlat/2)^2 + outer(cos(lat), cos(lat)) * sin(dlon/2)^2
  2 * 6371 * asin(pmin(1, sqrt(a)))
}


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

# fit_exponential_variogram has been relocated to mesma_helpers.R to
# provide a single implementation shared with january_averages.R and
# preprocess_data.R.  That file is sourced at the top of this script.

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

# --- Shared utility: plot all inference results for a given method ---
# Replaces 4 nearly identical PPI/MSAVI/NDVI/NoIndex inference plot blocks (~400 lines).
plot_inference_method_results <- function(full_data, method, file_prefix,
                                          use_excluded_years_shade = TRUE,
                                          include_species_plots = TRUE) {
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
# The default list is kept literal to avoid evaluation surprises in parallel workers.
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

# Remove large per-location-year outliers using spline residuals (or MAD when
# Whittaker smoothing is active).  Mirrors the definition in preprocess_data.R.
remove_large_outliers <- function(df, candidates = NULL, mad_thresh = OUTLIER_MAD_THRESHOLD) {
  if (!isTRUE(ENABLE_OUTLIER_REMOVAL)) return(df)

  interp_method <- resolve_interpolation_method()

  if (is.null(candidates)) {
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
  if (!"pheno_year" %in% names(df) && "date" %in% names(df)) {
    df$pheno_year <- assign_pheno_year(df$date)
  }

  grp <- interaction(df$location_id,
                     ifelse(is.na(df$pheno_year), "NA", as.character(df$pheno_year)),
                     drop = TRUE)
  groups <- split(seq_len(nrow(df)), grp)
  removed_idx <- logical(nrow(df))
  n_groups <- length(groups)

  if (identical(interp_method, "whittaker")) {
    for (g in seq_along(groups)) {
      rows <- groups[[g]]
      sub <- df[rows, , drop = FALSE]
      if (length(rows) < 5) { removed_idx[rows] <- TRUE; next }
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
      cat(sprintf("[OUTLIER] Removed %d observations across %d groups (MAD-only due to whittaker)\n",
                  sum(removed_idx, na.rm = TRUE), n_groups))
      df <- df[!removed_idx, , drop = FALSE]
    }
    return(df)
  }

  for (g in seq_along(groups)) {
    rows <- groups[[g]]
    sub <- df[rows, , drop = FALSE]
    if (length(rows) < 5) { removed_idx[rows] <- TRUE; next }
    out_mask <- rep(FALSE, nrow(sub))
    has_date <- "date" %in% names(sub) && any(!is.na(sub$date))
    if (!has_date || length(rows) < 10) next
    sub$doy <- as.numeric(format(sub$date, "%j"))
    for (col in candidates) {
      if (!is.numeric(sub[[col]])) next
      colv <- sub[[col]]
      finite_idx <- is.finite(colv) & is.finite(sub$doy)
      if (sum(finite_idx) < 5) next
      tryCatch({
        x <- sub$doy[finite_idx]; y <- colv[finite_idx]
        n_unique <- length(unique(x))
        fit1 <- stats::smooth.spline(x, y, df = min(OUTLIER_SPLINE_MAX_DF, length(x)/2, n_unique - 1))
        pred1 <- predict(fit1, x)$y
        res1 <- y - pred1
        mad1 <- stats::mad(res1, na.rm = TRUE)
        if (!is.finite(mad1) || mad1 <= 1e-6) return(NULL)
        keep_mask <- abs(res1 - stats::median(res1, na.rm = TRUE)) <= (mad_thresh * 1.5 * mad1)
        if (sum(keep_mask) >= 5) {
          n_unique2 <- length(unique(x[keep_mask]))
          fit2 <- stats::smooth.spline(x[keep_mask], y[keep_mask],
                                       df = min(OUTLIER_SPLINE_MAX_DF, sum(keep_mask)/2, n_unique2 - 1))
          pred_final <- predict(fit2, x)$y
        } else {
          pred_final <- pred1
        }
        residuals <- y - pred_final
        med_res <- stats::median(residuals, na.rm = TRUE)
        mad_res <- stats::mad(residuals, na.rm = TRUE)
        if (!is.finite(mad_res) || mad_res <= 0) return(NULL)
        this_mask <- rep(FALSE, length(colv))
        this_mask[finite_idx] <- abs(residuals - med_res) > mad_thresh * mad_res
        out_mask <<- out_mask | this_mask
      }, error = function(e) {})
    }
    if (any(out_mask, na.rm = TRUE)) removed_idx[rows[which(out_mask)]] <- TRUE
  }

  if (any(removed_idx, na.rm = TRUE)) {
    cat(sprintf("[OUTLIER] Removed %d observations across %d groups\n",
                sum(removed_idx, na.rm = TRUE), n_groups))
    df <- df[!removed_idx, , drop = FALSE]
  }
  df
}

# Compute all supported indices from raw bands.
# IMPORTANT: This must be worker-safe (parallel futures).
# NOTE: Indices that depend on a soil-line slope (e.g. WDVI) will be computed only
# if a finite slope is available (either passed in or present as SOIL_LINE_SLOPE).


# L2-normalize a feature vector per observation (whole-vector).
# Input: vec with n_indices * n_bins values.
# Output: vec / ||vec||_2 (with NA treated as 0 for the norm).
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
  valid <- !is.na(idx_pos)
  if (!any(valid)) return(vec)

  # Vectorize the z-scoring
  expanded_mu <- rep(means[idx_pos[valid]], each = n_bins)
  expanded_sigma <- rep(sds[idx_pos[valid]], each = n_bins)
  expanded_sigma[!is.finite(expanded_sigma) | expanded_sigma < eps_sigma] <- eps_sigma

  # Column indices to update
  k_valid <- which(valid)
  idx_starts <- (k_valid - 1) * n_bins + 1
  idx_ends <- k_valid * n_bins
  idx_all <- unlist(mapply(seq, idx_starts, idx_ends))

  out <- vec
  out[idx_all] <- (out[idx_all] - expanded_mu) / expanded_sigma
  out
}

mesma_zscore_mat_by_index <- function(mat, indices, means, sds, n_bins, eps_sigma = NULL) {
  if (is.null(eps_sigma)) {
    eps_sigma <- if (exists("EPS_SIGMA", inherits = TRUE)) get("EPS_SIGMA", inherits = TRUE) else 1e-8
  }

  if (is.null(mat) || nrow(mat) == 0 || length(indices) == 0) return(mat)

  idx_pos <- match(indices, names(means))
  valid <- !is.na(idx_pos)
  if (!any(valid)) return(mat)

  # Vectorize the z-scoring
  expanded_mu <- rep(means[idx_pos[valid]], each = n_bins)
  expanded_sigma <- rep(sds[idx_pos[valid]], each = n_bins)
  expanded_sigma[!is.finite(expanded_sigma) | expanded_sigma < eps_sigma] <- eps_sigma

  # Column indices to update
  k_valid <- which(valid)
  idx_starts <- (k_valid - 1) * n_bins + 1
  idx_ends <- k_valid * n_bins
  idx_all <- unlist(mapply(seq, idx_starts, idx_ends))

  out <- mat
  out[, idx_all] <- (out[, idx_all] - expanded_mu) / expanded_sigma
  out
}

# Compute diagonal covariance for endmember bundle sampling
# Uses per-band variance only (no cross-band covariance) which is more reliable
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
# Load preprocessed data from preprocess_data.R

cat("[NOTICE] Loading preprocessed data from preprocess_data.R outputs...\n")

df <- readRDS("preprocessed_data.rds")
# safety net: ensure no rows with phenology year before 1984 remain
cutoff <- 1984
if ("date" %in% names(df) || "pheno_year" %in% names(df)) {
  if (!"pheno_year" %in% names(df) && "date" %in% names(df)) {
    df$pheno_year <- ifelse(lubridate::month(df$date) >= 3,
                             lubridate::year(df$date),
                             lubridate::year(df$date) - 1)
  }
  n_before <- nrow(df)
  df <- df[is.na(df$pheno_year) | df$pheno_year >= cutoff, , drop = FALSE]
  if (n_before != nrow(df)) {
    cat(sprintf("[FILTER] fit_mesma removed %d rows with pheno_year < %d from preprocessed data\n", n_before - nrow(df), cutoff))
  }
}
gpts_map <- readRDS("gpts_map.rds")
TRAINING_NORM_PARAMS <- readRDS("training_norm_params.rds")

INDEX_SCALES <- TRAINING_NORM_PARAMS$INDEX_SCALES

# Reconstruct selected feature list (formerly built in preprocessing block)
if (!is.null(INDEX_SCALES) && length(INDEX_SCALES) > 0) {
  avail <- intersect(names(INDEX_SCALES), names(df))
} else {
  avail <- intersect(unique(c(OPTIMAL_INDICES, RAW_BANDS, "PPI")), names(df))
}
if (length(avail) == 0) {
  stop("No feature indices available after loading preprocessed data")
}

cat(sprintf("Selected %d indices from preprocessed data: %s\n",
            length(avail), paste(avail, collapse = ", ")))

cat(sprintf("[NOTICE] Loaded preprocessed data: %d rows, %d columns\n", nrow(df), ncol(df)))
cat(sprintf("[NOTICE] Loaded gpts_map with %d locations\n", nrow(gpts_map)))
cat(sprintf("[NOTICE] Loaded normalization params for %d indices\n", length(INDEX_SCALES)))

timing_info <- list()
timing_info$start_time <- Sys.time()

cat("Starting vegetation mixture analysis with MESMA approach...\n")

cat("\n")
cat(sprintf("Dataset size: %d rows, %d columns\n", nrow(df), ncol(df)))
cat(sprintf("Number of locations: %d\n", length(unique(df$location_id))))
cat(sprintf("Date range: %s to %s\n", min(df$date, na.rm = TRUE), max(df$date, na.rm = TRUE)))
cat(sprintf("Temporal aggregation: %d days (%d bins)\n", TEMPORAL_AGGREGATION_DAYS, TEMPORAL_BUDGET))

cat("\n=== TRAIN/INFERENCE DATA CONFIGURATION ===\n")
cat(sprintf("Training years (config): %s\n", paste(TRAIN_YEARS, collapse = ", ")))
if (isTRUE(ENABLE_LDA_L2_NORMALIZATION)) {
  cat("Feature representation: L2-normalized only\n")
} else {
  cat("Feature representation: Raw only\n")
}

if (exists("TRAIN_YEARS") && !is.null(TRAIN_YEARS) && length(TRAIN_YEARS) > 0) {
  cat(sprintf("Filtering training data to phenological years (March-February): %s\n", paste(TRAIN_YEARS, collapse = ", ")))
  df_train <- df[df$pheno_year %in% TRAIN_YEARS, , drop = FALSE]
} else {
  df_train <- df
}
# Apply balancing/downsampling to df_train at the LOCATION level for balanced training library.
# Downsample non-barren classes to the minimum number of unique locations among non-barren classes.
# Barren is kept as-is.
set.seed(get_mesma_seed(4))
loc_class <- df_train %>%
  dplyr::filter(.data$Veg != "barren") %>%
  dplyr::distinct(.data$location_id, .data$Veg)
loc_counts_per_class <- tapply(loc_class$location_id, loc_class$Veg, function(x) length(unique(x)))
if (length(loc_counts_per_class) > 0) {
  min_locs <- min(loc_counts_per_class)
  if (is.finite(min_locs) && min_locs > 0) {
    kept_locs <- loc_class %>%
      dplyr::group_by(.data$Veg) %>%
      dplyr::slice_sample(n = min_locs) %>%
      dplyr::ungroup() %>%
      dplyr::pull(.data$location_id)

    df_train_barren <- df_train %>% dplyr::filter(.data$Veg == "barren")
    df_train_non_barren <- df_train %>% dplyr::filter(.data$Veg != "barren" & .data$location_id %in% kept_locs)
    df_train <- dplyr::bind_rows(df_train_non_barren, df_train_barren)

    cat(sprintf(
      "[BALANCE] Downsampled NON-barren classes to %d locations/class; kept barren as-is (total=%d rows, %d locations)\n",
      min_locs, nrow(df_train), length(unique(df_train$location_id))
    ))
  }
}
cat(sprintf(
  "Training dataset (Initial): %d rows from %d locations\n",
  nrow(df_train), length(unique(df_train$location_id))
))

# --- STRATIFIED TRAIN/VALIDATION SPLIT ---
if (nrow(df_train) > 0) {
  cat(sprintf("[SPLIT] Performing stratified %.0f/%.0f train/validation split based on location_id and Veg...\n",
              100 * (1 - VALIDATION_FRACTION), 100 * VALIDATION_FRACTION))
  set.seed(get_mesma_seed(123))

  # Get dominant veg type per location for stratified splitting
  loc_veg_summary <- get_dominant_veg_per_location(df_train)

  # Stratified validation split: hold out VALIDATION_FRACTION of locations per Veg class
  val_locs_list <- vector("list", length(unique(loc_veg_summary$Veg)))
  unique_vegs_split <- unique(loc_veg_summary$Veg)
  for (i in seq_along(unique_vegs_split)) {
    v <- unique_vegs_split[i]
    v_locs <- loc_veg_summary$location_id[loc_veg_summary$Veg == v]
    n_v <- length(v_locs)
    n_val <- ceiling(n_v * VALIDATION_FRACTION)
    if (n_val > 0 && n_v > 1) {
      n_val <- min(n_val, n_v - 1)  # keep at least 1 location for training
      if (n_val > 0) {
        selected <- sample(v_locs, n_val)
        val_locs_list[[i]] <- data.frame(location_id = selected, Veg = v, stringsAsFactors = FALSE)
      }
    }
  }
  val_locs_df <- do.call(rbind, val_locs_list)
  validation_location_ids <- if (!is.null(val_locs_df) && nrow(val_locs_df) > 0) {
    as.character(val_locs_df$location_id)
  } else {
    character(0)
  }
  df_validation <- df_train[df_train$location_id %in% validation_location_ids, , drop = FALSE]
  df_train <- df_train[!df_train$location_id %in% validation_location_ids, , drop = FALSE]
  df_inference <- NULL

  cat(sprintf("[SPLIT] Validation set: %d rows from %d locations\n",
              nrow(df_validation), length(unique(df_validation$location_id))))
  cat(sprintf("[SPLIT] Training set: %d rows from %d locations\n",
              nrow(df_train), length(unique(df_train$location_id))))

  # ==========================================================================
  # CREATE OOB HOLDOUT FOR THRESHOLD TUNING (from training data)
  # ==========================================================================
  cat("\n[OOB SPLIT] Creating OOB holdout from training data for threshold/cluster tuning...\n")
  set.seed(get_mesma_seed(43))

  train_loc_veg <- get_dominant_veg_per_location(df_train)

  oob_locs_list <- vector("list", length(unique(train_loc_veg$Veg)))
  unique_train_vegs <- unique(train_loc_veg$Veg)

  for (i in seq_along(unique_train_vegs)) {
    v <- unique_train_vegs[i]
    v_locs <- train_loc_veg$location_id[train_loc_veg$Veg == v]
    n_v <- length(v_locs)
    n_oob <- ceiling(n_v * OOB_TUNING_FRACTION)

    if (n_oob > 0 && n_v > 1) {
      n_oob <- min(n_oob, n_v - 1)
      if (n_oob > 0) {
        selected <- sample(v_locs, n_oob)
        oob_locs_list[[i]] <- data.frame(location_id = selected, Veg = v, stringsAsFactors = FALSE)
      }
    }
  }

  oob_locs_df <- do.call(rbind, oob_locs_list)
  if (is.null(oob_locs_df) || nrow(oob_locs_df) == 0) {
    stop("[OOB SPLIT] No OOB locations could be selected")
  }

  oob_location_ids <- as.character(oob_locs_df$location_id)
  df_train_oob <- df_train[df_train$location_id %in% oob_location_ids, , drop = FALSE]
  df_train_model <- df_train[!df_train$location_id %in% oob_location_ids, , drop = FALSE]

  cat(sprintf("[OOB SPLIT] OOB tuning set: %d rows from %d locations\n",
              nrow(df_train_oob), length(unique(df_train_oob$location_id))))
  cat(sprintf("[OOB SPLIT] Model training set: %d rows from %d locations\n",
              nrow(df_train_model), length(unique(df_train_model$location_id))))
}

# Note: Snow/dust contamination filtering already applied early in the pipeline (before PPI baseline calculation)
# See [EARLY FILTERING] section above for details

df_test <- df
cat(sprintf(
  "Testing dataset: %d rows from %d locations\n",
  nrow(df_test), length(unique(df_test$location_id))
))

# Look for missing vegetation in original df (before train/test split) for augmentation
if ("Veg" %in% names(df) && length(ALLOWED_VEG) > 0) {
  # Check which classes are missing from df_train and add from df (original full dataset)
  original_df <- df  # df still has all data before it was assigned to df_train
  missing_vegs <- sapply(ALLOWED_VEG, function(v) {
    sum(tolower(df$Veg) == v, na.rm = TRUE)
  })
  missing_names <- names(missing_vegs)[missing_vegs == 0]
  if (length(missing_names) > 0) {
    safe_veg_subset <- function(d, veg_name) {
      if (!is.data.frame(d)) return(data.frame())
      if (!("Veg" %in% names(d))) return(data.frame())
      if (!isTRUE(nrow(d) > 0)) return(data.frame())
      idx <- !is.na(d$Veg) & tolower(d$Veg) == veg_name
      d[idx, , drop = FALSE]
    }
    for (mv in missing_names) {
      # Look in df_validation and df_inference combined for missing class samples
      cand <- rbind(
        if (exists("df_validation", inherits = FALSE)) safe_veg_subset(df_validation, mv) else data.frame(),
        if (exists("df_inference", inherits = FALSE)) safe_veg_subset(df_inference, mv) else data.frame()
      )
      if (nrow(cand) > 0) {
        add_n <- min(nrow(cand), max(5L, as.integer(floor(nrow(df) / 10))))
        add_n <- max(1L, add_n)
        to_add <- cand[seq_len(add_n), , drop = FALSE]
        df <- rbind(df, to_add)
        cat(sprintf("[NOTICE] Added %d samples from non-training years for Veg='%s' to ensure representation in training set\n", nrow(to_add), mv))
      } else {
        cat(sprintf("[WARNING] No samples found anywhere for Veg='%s'; cannot add examples to training set\n", mv))
      }
    }
    cat(sprintf("Training dataset after augmentation: %d rows from %d locations\n", nrow(df), length(unique(df$location_id))))
    df_train <- df
  }
}


cat("Using training data for vegetation library construction\n")
cat("=====================================\n\n")

if (isTRUE(PARALLEL_ENABLE)) {
  cleanup_parallel <- setup_parallel_backend(workers = PARALLEL_WORKERS)
} else {
  cleanup_parallel <- function() {}
}
on.exit(cleanup_parallel(), add = TRUE)


safe_lda_call <- function(X_pca, y, min_n_pcs = 2) {
  if (is.null(X_pca) || ncol(X_pca) < min_n_pcs) {
    cat(sprintf("safe_lda_call: Not enough PCs (have=%d, min=%d).\n", ncol(X_pca), min_n_pcs))
    return(NULL)
  }
  curr_n_pcs <- ncol(X_pca)
  while (curr_n_pcs >= min_n_pcs) {
    lda_res <- NULL
    warn_msg <- NULL
    withCallingHandlers({
      lda_res <- tryCatch({
        MASS::lda(X_pca[, 1:curr_n_pcs, drop = FALSE], grouping = y)
      }, error = function(e) {
        e
      })
    }, warning = function(w) {
      warn_msg <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    })

    if (inherits(lda_res, "error")) {
      cat(sprintf("safe_lda_call: LDA error: %s\n", lda_res$message))
      return(NULL)
    }

    if (is.null(warn_msg)) {
      return(lda_res)
    }

    if (grepl("collinear", warn_msg, ignore.case = TRUE)) {
      cat(sprintf("safe_lda_call: Received LDA warning '%s' -> reducing PCs from %d to %d and retrying.\n",
                  warn_msg, curr_n_pcs, curr_n_pcs - 1))
      curr_n_pcs <- curr_n_pcs - 1
      next
    }
    cat(sprintf("safe_lda_call: LDA warning (non-collinearity): %s\n", warn_msg))
    return(lda_res)
  }

  cat("safe_lda_call: Exhausted retries; LDA could not be computed without collinearity.\n")
  NULL
}



# ---------------------------------------------------------------------------
# Helper: apply endmember-level LDA grouping to y_labels.
# Returns updated y_labels vector (compound "class__location" when enabled,
# with small-group fallback to coarse class label).
# ---------------------------------------------------------------------------
.apply_endmember_grouping <- function(y_class, y_location) {
  use_endmember_lda <- exists("LDA_ENDMEMBER_GROUPING") && isTRUE(LDA_ENDMEMBER_GROUPING)
  if (!use_endmember_lda) return(y_class)

  min_per_group <- if (exists("LDA_ENDMEMBER_MIN_SAMPLES_PER_GROUP"))
                     as.integer(LDA_ENDMEMBER_MIN_SAMPLES_PER_GROUP) else 2L
  compound     <- paste(y_class, y_location, sep = "__")
  group_counts <- table(compound)
  small_groups <- names(group_counts)[group_counts < min_per_group]
  if (length(small_groups) > 0) {
    cat(sprintf("  [ENDMEMBER-LDA] %d location-class groups below min_samples=%d; merging back to class label\n",
                length(small_groups), min_per_group))
  }
  y_labels <- ifelse(compound %in% small_groups, y_class, compound)
  n_groups <- length(unique(y_labels))
  cat(sprintf("  [ENDMEMBER-LDA] Using %d compound endmember groups for LDA (vs %d coarse classes)\n",
              n_groups, length(unique(y_class))))
  y_labels
}

# ---------------------------------------------------------------------------
# Helper: compute z-score parameters from X_mat and return z-scored matrix
# plus the per-index means/sds vectors.
# ---------------------------------------------------------------------------
.train_z_score <- function(X_mat, feature_cols, n_bins) {
  n_total_indices <- length(feature_cols)
  global_means <- setNames(numeric(n_total_indices), feature_cols)
  global_sds   <- setNames(numeric(n_total_indices), feature_cols)
  X_z <- X_mat

  apply_zscore <- if (exists("ENABLE_ZSCORE_AFTER_L2")) isTRUE(ENABLE_ZSCORE_AFTER_L2) else TRUE

  if (apply_zscore) {
    cat(sprintf("  Computing Z-score parameters for %d indices...\n", n_total_indices))
    for (k in seq_along(feature_cols)) {
      col_start <- (k - 1) * n_bins + 1
      col_end   <- k * n_bins
      vals  <- X_mat[, col_start:col_end]
      mu    <- mean(vals, na.rm = TRUE)
      sigma <- sd(vals,   na.rm = TRUE)
      if (is.na(sigma) || sigma == 0) sigma <- 1
      global_means[k] <- mu
      global_sds[k]   <- sigma
      X_z[, col_start:col_end] <- (vals - mu) / sigma
    }
    X_z[!is.finite(X_z)] <- 0
  } else {
    cat("  Z-scoring DISABLED: using L2-normalized features directly\n")
    # Identity transform: un-zscoring is a no-op downstream
    global_means[] <- 0
    global_sds[]   <- 1
    X_z[!is.finite(X_z)] <- 0
  }

  list(X_z = X_z, means = global_means, sds = global_sds, zscore_applied = apply_zscore)
}

# ---------------------------------------------------------------------------
# Helper: run PCA then LDA on z-scored matrix; return LDA weights projected
# back into original feature space.
# ---------------------------------------------------------------------------
.train_pca_lda <- function(X_z, y_labels) {
  if (!exists("USE_LDA_SPACE_SOLVER")) {
    stop("[CONFIG] USE_LDA_SPACE_SOLVER is not set. Please define it in mesma_config.R.")
  }
  use_lda_space <- isTRUE(USE_LDA_SPACE_SOLVER)

  cat("  Computing PCA-LDA weights...\n")
  vars      <- apply(X_z, 2, var)
  keep_cols <- vars > 1e-9
  X_pca_in  <- X_z[, keep_cols, drop = FALSE]

  pca_res <- prcomp(X_pca_in, center = FALSE, scale. = FALSE)
  cum_var <- cumsum(pca_res$sdev^2) / sum(pca_res$sdev^2)
  n_pcs   <- which(cum_var > 0.95)[1]
  if (is.na(n_pcs)) n_pcs <- ncol(pca_res$x)

  class_counts     <- table(y_labels)
  n_min            <- if (length(class_counts) > 0) min(class_counts) else 0
  n_classes        <- length(unique(y_labels))
  min_pcs_for_lda  <- max(1, n_classes - 1)
  max_pcs_for_lda  <- min(20, max(min_pcs_for_lda, n_min - 2))

  if (n_pcs > max_pcs_for_lda) {
    warning(sprintf(
      "PCA->LDA: Reducing n_pcs from %d to %d to satisfy p < n_min constraint (smallest class has %d samples, max_pcs=20)",
      n_pcs, max_pcs_for_lda, n_min
    ))
    n_pcs <- max_pcs_for_lda
  }
  if (n_pcs < min_pcs_for_lda) {
    stop(sprintf(
      "[LDA] Not enough degrees of freedom for LDA (n_pcs=%d < min_required=%d, n_min=%d). Need more samples in smallest class.",
      n_pcs, min_pcs_for_lda, n_min
    ))
  }

  lda_res <- safe_lda_call(pca_res$x[, 1:n_pcs, drop = FALSE], as.factor(y_labels),
                            min_n_pcs = min_pcs_for_lda)
  if (is.null(lda_res)) {
    stop("[LDA] LDA could not be computed (collinearity / too few samples / invalid feature space)")
  }

  R      <- pca_res$rotation[, 1:n_pcs, drop = FALSE]
  W_pc   <- lda_res$scaling
  W_std  <- R %*% W_pc
  svd_v  <- lda_res$svd
  prop   <- svd_v^2 / sum(svd_v^2)
  n_dim  <- min(length(prop), ncol(W_std))

  if (ncol(W_std) > 1) {
    weights_clean <- rowSums(abs(W_std[, 1:n_dim, drop = FALSE]) %*% diag(prop[1:n_dim], nrow = n_dim))
  } else {
    weights_clean <- abs(W_std[, 1])
  }

  final_weights <- numeric(ncol(X_z))
  final_weights[keep_cols] <- weights_clean
  cat(sprintf("LDA weights (no normalization): min=%.4f, max=%.4f, mean=%.4f\n",
      min(final_weights[final_weights > 0], na.rm = TRUE),
      max(final_weights, na.rm = TRUE),
      mean(final_weights, na.rm = TRUE)))

  res <- list(final_weights = final_weights, use_lda_space = use_lda_space)

  if (use_lda_space) {
    lda_comp_weights <- prop[1:n_dim]
    # Build full-space basis: map zero-variance columns to zero rows
    W_full <- matrix(0, nrow = ncol(X_z), ncol = ncol(W_std))
    W_full[keep_cols, ] <- W_std
    res$lda_basis             <- W_full
    res$lda_component_weights <- lda_comp_weights
    # Store PCA loadings (rotation matrix for non-zero-variance columns only)
    # The plot script will compute X_pca_in the same way and use these loadings directly
    res$pca_rotation <- pca_res$rotation[, 1:n_pcs, drop = FALSE]
    res$pca_sdev <- pca_res$sdev[1:n_pcs]
    cat(sprintf("  LDA-space solver enabled: %d discriminant components, weights sum=%.3f\n",
                n_dim, sum(lda_comp_weights)))
  }

  res
}

train_feature_pipeline <- function(df, class_col, feature_cols) {
  cat(sprintf("\n=== Training Feature Pipeline for Class: %s ===\n", class_col))

  if (!exists(".fit_landsat_filter_plot_done", envir = .GlobalEnv) ||
      !isTRUE(get(".fit_landsat_filter_plot_done", envir = .GlobalEnv))) {
    obs_by_loc_year <- df |>
      dplyr::group_by(location_id, pheno_year) |>
      dplyr::summarise(n_obs = dplyr::n(), .groups = "drop")

    if (nrow(obs_by_loc_year) > 0) {
      min_required <- as.integer(MIN_PENTADS_PER_TRAIN_SAMPLE)
      keep_loc_year <- obs_by_loc_year |>
        dplyr::filter(n_obs >= min_required)

      n_locations_before <- dplyr::n_distinct(obs_by_loc_year$location_id)
      n_locations_after <- if (nrow(keep_loc_year) > 0) dplyr::n_distinct(keep_loc_year$location_id) else 0L

      cat(sprintf("[FILTER] Fit minimum Landsat images per location-year: %d\n", min_required))
      cat(sprintf("[FILTER] Locations surviving fit image filter: %d / %d\n",
                  n_locations_after, n_locations_before))

      if (requireNamespace("ggplot2", quietly = TRUE)) {
        survival_df <- data.frame(
          stage = factor(c("Before filter", sprintf(">= %d images", min_required)),
                         levels = c("Before filter", sprintf(">= %d images", min_required))),
          n_locations = c(n_locations_before, n_locations_after)
        )

        p_survival <- ggplot2::ggplot(survival_df, ggplot2::aes(x = stage, y = n_locations, fill = stage)) +
          ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
          ggplot2::geom_text(ggplot2::aes(label = n_locations), vjust = -0.3) +
          ggplot2::theme_minimal() +
          ggplot2::labs(
            title = "Fit Filter Survival",
            x = NULL,
            y = "Number of locations"
          )

        print(p_survival)
        # Plot saving disabled; OUTPUT_DIR ggsave call is skipped
      }
    }

    assign(".fit_landsat_filter_plot_done", TRUE, envir = .GlobalEnv)
  }
  
  X_raw      <- list()
  y_labels   <- c()   # fine-grained labels used for LDA
  y_class    <- c()   # coarse class labels (always vegetation type)
  y_location <- c()   # location id per sample

  class_values <- unique(na.omit(df[[class_col]]))
  traces_by_class <- lapply(class_values, function(cv) {
    df_class <- df[df[[class_col]] == cv, , drop = FALSE]
    if (nrow(df_class) == 0) return(list())
    split(df_class, list(df_class$location_id, df_class$pheno_year), drop = TRUE)
  })
  traces <- unlist(traces_by_class, recursive = FALSE)

  cat("  Building trace matrix...\n")
  for(sub in traces) {
    if(nrow(sub) < MIN_PENTADS_PER_TRAIN_SAMPLE) next

    mat <- build_pentad_matrix(sub, feature_cols)
    if(is.null(mat)) next

    vec <- as.numeric(mat)
    vec[!is.finite(vec)] <- NA

    X_raw[[length(X_raw)+1]] <- vec
    cls  <- names(sort(table(sub[[class_col]]), decreasing=TRUE))[1]
    loc  <- as.character(sub$location_id[1])
    y_class    <- c(y_class,    cls)
    y_location <- c(y_location, loc)
    y_labels   <- c(y_labels,   cls)   # default: coarse class label
  }

  # --- Endmember-level grouping ---
  y_labels <- .apply_endmember_grouping(y_class, y_location)

  if (length(X_raw) < 10) return(NULL)
  X_mat_raw <- do.call(rbind, X_raw)

  n_bins       <- TEMPORAL_BUDGET
  n_idx_local  <- length(feature_cols)
  l2_only_mode <- isTRUE(ENABLE_LDA_L2_NORMALIZATION)

  if (l2_only_mode) {
    cat(sprintf("  L2-normalizing training samples (per-observation) for %d indices (%d pentads each)...\n",
                n_idx_local, n_bins))
    X_mat <- t(apply(X_mat_raw, 1, function(r) l2_normalize_perindex(r, n_idx_local, n_bins)))
    cat(sprintf("  L2 normalization ENABLED: using %d indices.\n", n_idx_local))
  } else {
    X_mat <- X_mat_raw
    cat(sprintf("  L2 normalization DISABLED: using %d raw indices.\n", n_idx_local))
  }

  # --- Z-score ---
  zscore_res   <- .train_z_score(X_mat, feature_cols, n_bins)
  X_z          <- zscore_res$X_z
  global_means <- zscore_res$means
  global_sds   <- zscore_res$sds
  apply_zscore <- zscore_res$zscore_applied

  # --- PCA-LDA weights ---
  lda_res_out <- .train_pca_lda(X_z, y_labels)

  # Build return list
  res <- list(
    means                 = global_means,
    sds                   = global_sds,
    weights               = lda_res_out$final_weights,
    indices               = feature_cols,
    base_indices          = feature_cols,
    l2_normalize          = l2_only_mode,
    zscore_applied        = apply_zscore,
    lda_basis             = lda_res_out$lda_basis,
    lda_component_weights = lda_res_out$lda_component_weights
  )
  return(res)
}

doy_to_pentad <- function(doy) {
  pmin(ceiling(doy / TEMPORAL_AGGREGATION_DAYS), TEMPORAL_BUDGET)
}

# ---------------------------------------------------------------------------
# Helper: aggregate daily observations into a pentad matrix (no interpolation).
# Returns a [TEMPORAL_BUDGET x length(avail_idx)] matrix with NAs for missing
# pentads. The final pentad row is always NA (intentionally excluded).
# Returns NULL if the input is empty or all observations are in the last pentad.
# ---------------------------------------------------------------------------
.build_pentad_from_observations <- function(dly_year, avail_idx) {
  if (is.null(dly_year) || nrow(dly_year) == 0) return(NULL)

  # CRITICAL: Use phenological DOY (March 1 = day 1), not calendar DOY
  # This ensures temporal alignment when data spans phenological years (March-February)
  if (!"doy" %in% names(dly_year) || any(is.na(dly_year$doy))) {
    # Local implementation (keeps this function self-contained for future workers)
    local_pheno_doy <- function(d) {
      d <- as.Date(d)
      yr <- as.integer(format(d, "%Y"))
      march1 <- as.Date(paste0(yr, "-03-01"))
      pheno_start <- ifelse(d >= march1, march1, as.Date(paste0(yr - 1L, "-03-01")))
      as.integer(d - as.Date(pheno_start) + 1L)
    }
    dly_year$doy <- local_pheno_doy(dly_year$date)
  }

  dly_year$pentad <- doy_to_pentad(dly_year$doy)

  # ---------------------------------------------------------------------------
  # **IMPORTANT**: do not include the last pentad in any analysis.  The final
  # pentad (equal to `TEMPORAL_BUDGET`) is often incomplete/reliable and must
  # be treated as missing.  To guarantee it is never used we drop any
  # observations from that pentad and leave the corresponding row of the
  # output matrix as NA.  Later steps implicitly convert NA -> 0 and constant
  # columns are pruned, so the last pentad will have no influence.
  # ---------------------------------------------------------------------------
  if (any(dly_year$pentad >= TEMPORAL_BUDGET, na.rm = TRUE)) {
    n_dropped <- sum(dly_year$pentad >= TEMPORAL_BUDGET, na.rm = TRUE)
    loc_id <- if ("location_id" %in% names(dly_year)) paste(unique(dly_year$location_id), collapse = ",") else "unknown"
    dly_year <- dly_year[dly_year$pentad < TEMPORAL_BUDGET, , drop = FALSE]
    if (nrow(dly_year) == 0) {
      # All observations fell in the final (incomplete) pentad — nothing usable.
      warning(sprintf(
        "[build_pentad_matrix] All %d observations for location '%s' fell in the final pentad (>= TEMPORAL_BUDGET=%d) and were dropped; returning NULL.",
        n_dropped, loc_id, TEMPORAL_BUDGET
      ))
      return(NULL)
    }
  }

  K <- length(avail_idx)
  pentad_mat <- matrix(NA_real_, nrow = TEMPORAL_BUDGET, ncol = K)
  colnames(pentad_mat) <- avail_idx

  # iterate only up to the penultimate pentad; the final row is intentionally
  # left as NA and will be zeroed out after interpolation
  for (p in seq_len(TEMPORAL_BUDGET - 1)) {
    subset_p <- dly_year[dly_year$pentad == p, ]
    if (nrow(subset_p) == 0) next

    # Simple mean of the bin's valid values
    for (j in seq_along(avail_idx)) {
      idx <- avail_idx[j]
      if (!idx %in% names(subset_p)) next

      v <- subset_p[[idx]]
      v <- v[is.finite(v)]
      if (length(v) == 0) next

      pentad_mat[p, j] <- mean(v, na.rm = TRUE)
    }
  }

  pentad_mat
}

# ---------------------------------------------------------------------------
# Helper: apply temporal interpolation / smoothing to a pentad matrix.
# `interp_method` must be one of "linear", "whittaker", or "none"/NULL.
# Always clears the final pentad row (TEMPORAL_BUDGET) after filling.
# ---------------------------------------------------------------------------
.apply_temporal_interpolation <- function(pentad_mat, interp_method) {
  K <- ncol(pentad_mat)
  if (!is.null(interp_method) && interp_method != "none") {
    if (interp_method == "linear") {
      for (j in seq_len(K)) {
        vals <- pentad_mat[, j]
        if (any(is.na(vals))) {
          if (all(is.na(vals))) {
            pentad_mat[, j] <- 0
          } else {
            idx_present <- which(!is.na(vals))
            if (length(idx_present) >= 2) {
              pentad_mat[, j] <- approx(idx_present, vals[idx_present], xout = seq_len(TEMPORAL_BUDGET), rule = 2)$y
            } else {
              pentad_mat[, j] <- vals[idx_present[1]]
            }
          }
        }
      }
    } else if (interp_method == "whittaker") {
      for (j in seq_len(K)) {
        pentad_mat[, j] <- whittaker_smooth(pentad_mat[, j])
      }
    } else {
      stop(sprintf("Unsupported interpolation method '%s'", interp_method))
    }
  }
  # Always clear the final pentad row — it is intentionally excluded from analysis
  pentad_mat[TEMPORAL_BUDGET, ] <- NA_real_
  pentad_mat
}

build_pentad_matrix <- function(dly_year, avail_idx, interpolate = NULL) {
  interp_method <- if (is.null(interpolate)) {
    resolve_interpolation_method(fallback = "linear")
  } else {
    get_interpolate_method(interpolate)
  }
  pentad_mat <- .build_pentad_from_observations(dly_year, avail_idx)
  if (is.null(pentad_mat)) return(NULL)
  .apply_temporal_interpolation(pentad_mat, interp_method)
}

 





.run_map <- function(X, FUN, show_pb = TRUE) {
  f_FUN <- FUN
  
  if (!PARALLEL_ENABLE) {
    lapply(X, function(x) { f_FUN(x) })
  } else {
    # Parallel execution uses the global future plan configured earlier (setup_parallel_backend).
    # Do NOT create/tear down clusters per call: that adds major overhead and can look like a hang.
    if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
      stop("Packages 'future' and 'future.apply' are required for parallel processing")
    }

    options(future.globals.maxSize = max(getOption("future.globals.maxSize", 2e9), 12e9))

    current_plan <- future::plan()
    plan_is_sequential <- inherits(current_plan, "sequential") || inherits(current_plan, "uniprocess")
    if (isTRUE(plan_is_sequential)) {
      if (!exists(".RUN_MAP_WARNED_SEQUENTIAL", envir = .GlobalEnv)) {
        assign(".RUN_MAP_WARNED_SEQUENTIAL", TRUE, envir = .GlobalEnv)
        cat("[PARALLEL] NOTE: PARALLEL_ENABLE=TRUE but future plan is sequential; running sequentially.\n")
        cat("[PARALLEL]       If you want parallel, ensure setup_parallel_backend() ran successfully.\n")
      }
      return(lapply(X, f_FUN))
    }

    future.apply::future_lapply(X, f_FUN, future.seed = TRUE)
  }
}



cat("===============================================\n\n")

lib_df <- df
# No artificial per-vegetation sampling applied; use full lib_df for variant construction

vegs <- unique(na.omit(lib_df$Veg))
vegs <- vegs[vegs != ""]
vegs <- vegs[tolower(vegs) %in% c("herbs", "populus", "tamarix", "barren")]  # FIXED: case-insensitive matching

# Filter lib_df to only include the selected vegetation types
lib_df <- lib_df[tolower(lib_df$Veg) %in% tolower(vegs), ]


lib <- list()
for (v in vegs) {
  lib[[v]] <- list(n_samples = 0)
}

timing_info$lib_construction_done <- Sys.time()

cat("=== Building raw index library ===\n")

feature_cols <- avail
cat(sprintf("Raw index features: %s\n", paste(feature_cols, collapse=", ")))

X_all <- as.matrix(lib_df[, feature_cols, drop = FALSE])
for (j in seq_len(ncol(X_all))) {
  col_vals <- X_all[, j]
  if (any(!is.finite(col_vals))) {
    mu_j <- mean(col_vals[is.finite(col_vals)], na.rm = TRUE)
    if (!is.finite(mu_j)) mu_j <- 0
    X_all[!is.finite(col_vals), j] <- mu_j
  }
}
mu_all <- colMeans(X_all)
feature_sds <- apply(X_all, 2, sd)
feature_sds[feature_sds <= 1e-10] <- 1.0

n_features <- length(feature_cols)
cat(sprintf("Raw index library: %d features from %d training samples\n", n_features, nrow(X_all)))

cat("Computing raw index templates per vegetation type...\n")
raw_lib_templates <- list()

for (vname in vegs) {
  dveg <- lib_df[lib_df$Veg == vname, , drop = FALSE]

  if (!"date" %in% names(dveg)) {
    cat(sprintf("[ERROR] %s: 'date' column missing from dveg!\n", vname))
    next
  }

  X_v <- as.matrix(dveg[, feature_cols, drop = FALSE])
  for (j in seq_len(ncol(X_v))) {
    col_vals <- X_v[, j]
    if (any(!is.finite(col_vals))) {
      X_v[!is.finite(col_vals), j] <- mu_all[j]
    }
  }
  
  X_v_c <- sweep(X_v, 2, mu_all, "-")
  X_v_std <- sweep(X_v_c, 2, feature_sds, "/")
  
  doy_vec <- pheno_doy(dveg$date)  # Use phenological DOY
  pentad_vec <- doy_to_pentad(doy_vec)

  if (all(is.na(pentad_vec))) {
    cat(sprintf("  [ERROR] %s: All pentad values are NA\n", vname))
    next
  }

  T_medoid <- matrix(NA_real_, nrow = TEMPORAL_BUDGET, ncol = n_features)
  n_filled <- 0

  for (p in seq_len(TEMPORAL_BUDGET)) {
    rows_p <- which(pentad_vec == p)
    if (length(rows_p) > 0) {
      sub <- X_v_std[rows_p, , drop = FALSE]
      if (nrow(sub) == 1) {
        T_medoid[p, ] <- sub[1, ]
      } else {
        # Use median center for medoid selection within pentad
        center_med <- apply(sub, 2, median, na.rm = TRUE)
        dists <- rowSums(sweep(sub, 2, center_med, "-")^2)
        T_medoid[p, ] <- sub[which.min(dists), ]
      }
      n_filled <- n_filled + 1
    }
  }

  if (n_filled == 0) {
    cat(sprintf("  [ERROR] %s: No pentads filled!\n", vname))
    next
  }
  
  
  raw_lib_templates[[vname]] <- list(T = T_medoid, n_samples = nrow(dveg))
  cat(sprintf("  %s: T_medoid range [%.4f, %.4f], mean=%.4f (filled %d/%d pentads)\n",
              vname, min(T_medoid, na.rm=TRUE), max(T_medoid, na.rm=TRUE), mean(T_medoid, na.rm=TRUE),
              n_filled, TEMPORAL_BUDGET))
}

cat("Raw index library templates computed.\n")

# --- GLOBAL HELPER: Bootstrap Median Vector ---
# Defined at global scope so all bootstrap functions can use it
boot_median_vec <- function(x, n_boot) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(rep(NA_real_, n_boot))
  if (length(x) < 3) return(rep(median(x), n_boot)) # Too few to bootstrap meaningfully

  # Fast vectorized bootstrap of median
  replicate(n_boot, median(sample(x, length(x), replace = TRUE), na.rm = TRUE))
}

# --- GLOBAL HELPER: Spatial Block Bootstrap over Locations ---
# Builds coarse spatial blocks from per-location lat/lon and resamples blocks.
# This preserves spatial dependence better than i.i.d. location resampling.
collect_location_coords <- function(locations, df_tasks = NULL, all_coefs = NULL) {
  locations <- trimws(as.character(locations))
  locations <- locations[!is.na(locations) & locations != ""]
  if (length(locations) == 0) return(data.frame(location_id = character(0), lat = numeric(0), lon = numeric(0)))

  cand <- NULL

  # Prefer df_tasks if it has usable coordinates
  if (!is.null(df_tasks) && all(c("location_id", "lat", "lon") %in% names(df_tasks))) {
    cand <- df_tasks[, c("location_id", "lat", "lon"), drop = FALSE]
  } else if (!is.null(all_coefs) && all(c("location_id", "lat", "lon") %in% names(all_coefs))) {
    cand <- all_coefs[, c("location_id", "lat", "lon"), drop = FALSE]
  }

  if (!is.null(cand)) {
    cand$location_id <- trimws(as.character(cand$location_id))
    cand$lat <- suppressWarnings(as.numeric(cand$lat))
    cand$lon <- suppressWarnings(as.numeric(cand$lon))
    cand <- cand[!is.na(cand$location_id) & cand$location_id != "", , drop = FALSE]

    # Choose first non-missing coord per location (prefer finite lat/lon)
    cand <- deduplicate_coords(cand)
  }

  # Fall back to global gpts_map if available
  if ((is.null(cand) || nrow(cand) == 0 || all(is.na(cand$lat)) || all(is.na(cand$lon))) && exists("gpts_map", envir = globalenv())) {
    gm <- get("gpts_map", envir = globalenv())
    if (!is.null(gm) && all(c("location_id", "lat", "lon") %in% names(gm))) {
      gm2 <- gm[, c("location_id", "lat", "lon"), drop = FALSE]
      gm2$location_id <- trimws(as.character(gm2$location_id))
      gm2$lat <- suppressWarnings(as.numeric(gm2$lat))
      gm2$lon <- suppressWarnings(as.numeric(gm2$lon))
      gm2 <- gm2[!is.na(gm2$location_id) & gm2$location_id != "", , drop = FALSE]
      gm2 <- deduplicate_coords(gm2)
      cand <- gm2
    }
  }

  if (is.null(cand) || nrow(cand) == 0) {
    return(data.frame(location_id = locations, lat = NA_real_, lon = NA_real_))
  }

  # Align to the requested locations list
  cand <- cand[match(locations, cand$location_id), , drop = FALSE]
  if (nrow(cand) == 0) {
    return(data.frame(location_id = locations, lat = NA_real_, lon = NA_real_))
  }
  cand
}

# estimate_autocorrelation_range now lives in mesma_helpers.R and is
# sourced at startup.  That centralisation ensures january_averages.R and the
# main MESMA script behave identically.


# Spatial autocorrelation and block bootstrap helpers are defined in
# mesma_helpers.R and sourced earlier in this script.  The definitions have
# See the helper file for the shared implementation.
# common implementation.

# PPI bootstrap: thin wrapper around the generic location_bootstrap_index
location_bootstrap_ppi <- function(all_coefs, df_tasks, B = BOOTSTRAP_B, seed = 123) {
  ppi_max <- if (exists("PPI_FULL_VEG_COVER")) get("PPI_FULL_VEG_COVER") else 0.4
  if (ppi_max <= 0) ppi_max <- 0.4
  result <- location_bootstrap_index(all_coefs, df_tasks, index_name = "PPI",
                                     B = B, seed = seed, index_max = ppi_max)
  if (!is.null(result) && "Veg" %in% names(result)) {
    result$Veg <- ifelse(normalize_veg_name(result$Veg) == "phragmites", "herbs", result$Veg)
  }
  result
}

# Shared helper for index-based location bootstrap (summer detrended medians).
# When noindex=TRUE, uses raw coefficients as absolute cover (no index scaling).
location_bootstrap_index <- function(all_coefs, df_tasks = NULL, index_name = "NOINDEX",
                                     B = BOOTSTRAP_B, seed = 123,
                                     index_max = 0.6,
                                     summer_months = 6:9,
                                     detrend_poly_degree = 3,
                                     min_detrend_n = 50,
                                     noindex = FALSE) {
  set.seed(seed)
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    warning(sprintf("dplyr required for %s bootstrap", index_name))
    return(NULL)
  }

  # --- NOINDEX mode: uses raw coefficients as absolute cover ---
  if (isTRUE(noindex)) {
    merged_veg <- filter_valid_vegetation(all_coefs, exclude_barren = TRUE)
    merged_veg$rel_coef <- merged_veg$coef

    veg_types <- unique(merged_veg$Veg[!tolower(merged_veg$Veg) %in% c("barren")])
    years <- sort(unique(merged_veg$pheno_year[!is.na(merged_veg$pheno_year)]))
    locations <- unique(merged_veg$location_id)
    n_locs <- length(locations)
    if (length(years) == 0 || n_locs == 0 || nrow(merged_veg) == 0) { warning("[NOINDEX BOOTSTRAP] No valid data after filtering - skipping bootstrap"); return(all_coefs) }
    if (!is.numeric(B) || length(B) != 1 || B < 1) { warning("[NOINDEX BOOTSTRAP] Invalid B parameter - skipping bootstrap"); return(all_coefs) }
    B <- as.integer(B)

    .loc_coords <- collect_location_coords(locations = locations, df_tasks = df_tasks, all_coefs = all_coefs)
    loc_means <- tapply(merged_veg$coef, merged_veg$location_id, mean, na.rm = TRUE)
    coords_for_range <- .loc_coords[match(names(loc_means), .loc_coords$location_id), , drop = FALSE]
    .block_km <- block_km_if_significant(coords_for_range, as.numeric(loc_means))
    cat(sprintf("[NOINDEX BOOTSTRAP] Block size: %.1f km (%s)\n", .block_km, if (.block_km > 0) "spatial" else "i.i.d."))

    unique_loc_years <- merged_veg |> dplyr::distinct(location_id, pheno_year)
    veg_boot_res <- list()
    for (v in c(veg_types, "barren")) {
      mat <- matrix(NA_real_, nrow = B, ncol = length(years))
      if (ncol(mat) > 0) colnames(mat) <- as.character(years)
      veg_boot_res[[v]] <- mat
    }

    for (b in seq_len(B)) {
      boot_locs <- spatial_block_sample_locations(locations, .loc_coords, n_draw = n_locs, block_km = .block_km)
      loc_counts <- table(boot_locs)
      sub_veg <- merged_veg[merged_veg$location_id %in% names(loc_counts), ]
      if (nrow(sub_veg) > 0) {
        sub_veg$weight <- as.integer(loc_counts[sub_veg$location_id])
        agg_veg <- sub_veg |> dplyr::group_by(pheno_year, Veg) |> dplyr::summarize(total_cover = sum(rel_coef * weight, na.rm = TRUE), .groups = "drop")
        for (i in 1:nrow(agg_veg)) {
          v <- agg_veg$Veg[i]; y <- as.character(agg_veg$pheno_year[i])
          if (v %in% veg_types && y %in% colnames(veg_boot_res[[v]])) veg_boot_res[[v]][b, y] <- agg_veg$total_cover[i] / n_locs
        }
        loc_year_sums <- sub_veg |> dplyr::group_by(location_id, pheno_year) |> dplyr::summarize(total_veg = sum(rel_coef, na.rm = TRUE), .groups = "drop")
        loc_year_sums$weight <- as.integer(loc_counts[as.character(loc_year_sums$location_id)])
        loc_year_sums$barren <- pmax(0, pmin(1, 1.0 - loc_year_sums$total_veg))
        barren_by_year <- loc_year_sums |> dplyr::group_by(pheno_year) |> dplyr::summarize(mean_barren = sum(barren * weight, na.rm = TRUE) / sum(weight, na.rm = TRUE), .groups = "drop")
        for (i in 1:nrow(barren_by_year)) {
          y <- as.character(barren_by_year$pheno_year[i])
          if (y %in% colnames(veg_boot_res[["barren"]])) veg_boot_res[["barren"]][b, y] <- barren_by_year$mean_barren[i]
        }
      }
    }
    return(compile_bootstrap_results(veg_boot_res, years, unique_loc_years, method_name = "location_bootstrap_noindex"))
  }

  index_name <- as.character(index_name)
  raw_col <- paste0(index_name, "_raw")
  idx_col <- NULL
  if (raw_col %in% names(df_tasks) && any(!is.na(df_tasks[[raw_col]]))) {
    idx_col <- raw_col
    cat(sprintf("[%s BOOTSTRAP] Using raw %s values (%s) for weighting.\n", index_name, index_name, raw_col))
  } else if (index_name %in% names(df_tasks)) {
    idx_col <- index_name
    cat(sprintf("[%s BOOTSTRAP] Using standard %s column.\n", index_name, index_name))
  } else {
    warning(sprintf("%s column missing from df_tasks", index_name))
    return(NULL)
  }

  # Ensure temporal columns exist
  if (!"month" %in% names(df_tasks)) {
    if ("date" %in% names(df_tasks)) df_tasks$month <- as.integer(format(as.Date(df_tasks$date), "%m")) else df_tasks$month <- 1L
  }
  if (!"doy" %in% names(df_tasks)) {
    if ("date" %in% names(df_tasks)) df_tasks$doy <- as.numeric(format(as.Date(df_tasks$date), "%j")) else df_tasks$doy <- 150
  }

  summer_df <- df_tasks |> dplyr::filter(month %in% summer_months)
  detrended_col <- paste0(tolower(index_name), "_detrended")
  if (nrow(summer_df) > min_detrend_n) {
    tryCatch({
      f <- as.formula(paste(idx_col, "~ poly(doy,", detrend_poly_degree, ")"))
      seasonal_model <- lm(f, data = summer_df)
      summer_df$seasonal_trend <- predict(seasonal_model, newdata = summer_df)
      global_seasonal_mean <- mean(summer_df$seasonal_trend, na.rm = TRUE)
      summer_df[[detrended_col]] <- as.numeric(summer_df[[idx_col]] - (summer_df$seasonal_trend - global_seasonal_mean))
      cat(sprintf("[%s BOOTSTRAP] Computed detrended summer %s (N=%d).\n", index_name, index_name, nrow(summer_df)))
    }, error = function(e) {
      cat(sprintf("[%s BOOTSTRAP] Detrending failed: %s. Using raw %s.\n", index_name, e$message, index_name))
      summer_df[[detrended_col]] <- as.numeric(summer_df[[idx_col]])
    })
  } else {
    summer_df[[detrended_col]] <- as.numeric(summer_df[[idx_col]])
    cat(sprintf("[%s BOOTSTRAP] Insufficient summer data for detrending; using raw values.\n", index_name))
  }

  # Build per-location-year bootstrap matrix (split/lapply avoids dplyr scoping pitfalls)
  if (nrow(summer_df) == 0) {
    boot_combined <- data.frame(location_id = character(0), pheno_year = integer(0), boot_list = I(list()))
  } else {
    summer_df$grp_key <- paste(summer_df$location_id, summer_df$pheno_year, sep = "|||")
    grp <- split(summer_df, summer_df$grp_key, drop = TRUE)
    grp_names <- names(grp)
    boot_lists <- vector("list", length(grp))
    for (gi in seq_along(grp)) {
      g <- grp[[gi]]
      vec <- as.numeric(g[[detrended_col]])
      if (sum(is.finite(vec)) > 0) boot_lists[[gi]] <- boot_median_vec(vec, as.integer(B)) else boot_lists[[gi]] <- rep(NA_real_, as.integer(B))
    }
    parsed <- do.call(rbind, strsplit(grp_names, "\\|\\|\\|"))
    boot_combined <- data.frame(location_id = parsed[, 1], pheno_year = as.integer(parsed[, 2]), boot_list = I(boot_lists), stringsAsFactors = FALSE)
  }
  boot_combined$key <- paste(boot_combined$location_id, boot_combined$pheno_year, sep = "_")
  n_rows <- nrow(boot_combined)
  boot_mat <- matrix(NA_real_, nrow = n_rows, ncol = as.integer(B))
  for (i in 1:n_rows) {
    s_vec <- boot_combined$boot_list[[i]]
    if (!is.null(s_vec) && !all(is.na(s_vec))) boot_mat[i, ] <- s_vec
  }
  key_map <- setNames(seq_len(n_rows), boot_combined$key)

  # Prepare coefficient data (exclude barren rows)
  merged_veg <- filter_valid_vegetation(all_coefs, exclude_barren = TRUE)
  merged_veg$key <- paste(merged_veg$location_id, merged_veg$pheno_year, sep = "_")
  merged_veg$row_idx <- key_map[merged_veg$key]
  merged_veg <- merged_veg[!is.na(merged_veg$row_idx), ]

  # Relative coef normalization (same logic as PPI bootstrap)
  merged_veg <- merged_veg |>
    dplyr::group_by(location_id, pheno_year) |>
    dplyr::mutate(
      sum_veg_coef = sum(coef[tolower(Veg) != "barren"], na.rm = TRUE),
      rel_coef = ifelse(sum_veg_coef > 1e-6, coef / sum_veg_coef, coef)
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(location_id, pheno_year) |>
    dplyr::mutate(
      total_rel = sum(rel_coef, na.rm = TRUE),
      rel_coef = ifelse(total_rel > 1e-6 & abs(total_rel - 1) > 1e-9, rel_coef / total_rel, rel_coef)
    ) |>
    dplyr::select(-total_rel) |>
    dplyr::ungroup()

  veg_types <- unique(merged_veg$Veg[!tolower(merged_veg$Veg) %in% c("barren")])
  years <- sort(unique(merged_veg$pheno_year[!is.na(merged_veg$pheno_year)]))
  locations <- unique(merged_veg$location_id)
  n_locs <- length(locations)

  if (length(years) == 0 || n_locs == 0 || nrow(merged_veg) == 0) {
    warning(sprintf("[%s BOOTSTRAP] No valid data after filtering - skipping bootstrap", index_name))
    return(NULL)
  }
  if (!is.numeric(B) || length(B) != 1 || B < 1) {
    warning(sprintf("[%s BOOTSTRAP] Invalid B parameter - skipping bootstrap", index_name))
    return(NULL)
  }
  B <- as.integer(B)

  unique_loc_years <- merged_veg |> dplyr::distinct(location_id, pheno_year, row_idx)

  veg_boot_res <- list()
  all_veg_plus_barren <- c(veg_types, "Barren")
  for (v in all_veg_plus_barren) {
    mat <- matrix(NA_real_, nrow = as.integer(B), ncol = as.integer(length(years)))
    if (is.matrix(mat) && ncol(mat) > 0) colnames(mat) <- as.character(years)
    veg_boot_res[[v]] <- mat
  }

  # Precompute per-location coordinates and data-driven block size
  .loc_coords <- collect_location_coords(locations = locations, df_tasks = df_tasks, all_coefs = all_coefs)
  loc_means <- tapply(merged_veg$coef, merged_veg$location_id, mean, na.rm = TRUE)
  loc_ids_for_range <- names(loc_means)
  coords_for_range <- .loc_coords[match(loc_ids_for_range, .loc_coords$location_id), , drop = FALSE]
  .block_km <- block_km_if_significant(coords_for_range, as.numeric(loc_means))
  cat(sprintf("[%s BOOTSTRAP] Block size: %.1f km (%s)\n", index_name, .block_km, if (.block_km > 0) "spatial" else "i.i.d."))

  for (b in seq_len(B)) {
    boot_locs <- spatial_block_sample_locations(locations, .loc_coords, n_draw = n_locs, block_km = .block_km)
    loc_counts <- table(boot_locs)
    sampled_locs <- names(loc_counts)

    raw_vals <- boot_mat[, b]
    norm_vals <- pmin(pmax(raw_vals / index_max, 0), 1)

    sub_veg <- merged_veg[merged_veg$location_id %in% sampled_locs, ]
    if (nrow(sub_veg) > 0) {
      sub_veg$current_norm <- norm_vals[sub_veg$row_idx]
      sub_veg$abs_cover <- sub_veg$rel_coef * sub_veg$current_norm
      sub_veg$weight <- as.integer(loc_counts[sub_veg$location_id])

      agg_veg <- sub_veg |>
        dplyr::group_by(pheno_year, Veg) |>
        dplyr::summarize(total_cover = sum(abs_cover * weight, na.rm = TRUE), .groups = "drop")

      for (i in 1:nrow(agg_veg)) {
        v <- agg_veg$Veg[i]
        y <- as.character(agg_veg$pheno_year[i])
        if (v %in% veg_types && y %in% colnames(veg_boot_res[[v]])) {
          veg_boot_res[[v]][b, y] <- agg_veg$total_cover[i] / n_locs
        }
      }
    }

    sub_barren <- unique_loc_years[unique_loc_years$location_id %in% sampled_locs, ]
    if (nrow(sub_barren) > 0) {
      sub_barren$current_norm <- norm_vals[sub_barren$row_idx]
      sub_barren$barren_frac <- 1.0 - sub_barren$current_norm
      sub_barren$weight <- as.integer(loc_counts[sub_barren$location_id])

      agg_barren <- sub_barren |>
        dplyr::group_by(pheno_year) |>
        dplyr::summarize(
          total_barren = sum(barren_frac * weight, na.rm = TRUE),
          n_obs_year = sum(weight, na.rm = TRUE),
          .groups = "drop"
        )

      for (i in 1:nrow(agg_barren)) {
        y <- as.character(agg_barren$pheno_year[i])
        if (y %in% colnames(veg_boot_res[["Barren"]])) {
          veg_boot_res[["Barren"]][b, y] <- agg_barren$total_barren[i] / agg_barren$n_obs_year[i]
        }
      }
    }
  }

  compile_bootstrap_results(veg_boot_res, years, unique_loc_years,
                            method_name = paste0("location_bootstrap_", tolower(index_name)))
}


# Location-based bootstrap for global aggregation
location_bootstrap_aggregate <- function(all_coefs, B = BOOTSTRAP_B, seed = 123) {
  set.seed(seed)

  veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
  years <- sort(unique(all_coefs$pheno_year[!is.na(all_coefs$pheno_year)]))
  results_list <- list()

  # Guard against empty data
  if (length(years) == 0 || length(veg_types) == 0) {
    warning("[BOOTSTRAP] No valid years or veg types - returning empty results")
    return(data.frame())
  }

  for (veg in veg_types) {
    # Filter data for this veg
    veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), ]
    if (nrow(veg_data) == 0) next

    # --- POOLED VARIANCE CALCULATION ---
    # Calculate robust SD of the coefficients across the entire time series (or per year)
    # to use as a fallback proxy for spatial variability when N is small.

    # Calculate SD per year where we have at least 3 points
    sd_by_year <- tapply(veg_data$coef, veg_data$pheno_year, function(x) {
      if (length(x) >= 3) sd(x, na.rm = TRUE) else NA_real_
    })

    # Use median of yearly SDs as the "typical spatial variability" for this class
    pooled_sd <- median(sd_by_year, na.rm = TRUE)

    cat(sprintf("[BOOTSTRAP] Veg '%s': Pooled Spatial SD = %.4f (used for N < 8)\n", veg, pooled_sd))

    boot_means <- matrix(NA_real_, nrow = B, ncol = length(years))
    if (length(years) > 0) {
      colnames(boot_means) <- as.character(years)
    }
    
    n_eff_vec <- numeric(length(years))

    for (i in seq_along(years)) {
      yr <- years[i]
      yr_data <- veg_data[veg_data$pheno_year == yr, ]
      n_obs <- nrow(yr_data)

      # Estimate effective sample size using pairwise spatial correlation
      # Default: n_eff = n_obs (assume independence if no spatial info)
      n_eff_est <- n_obs
      if (n_obs >= 3) {
        # Check for pre-existing coordinates
        if (all(c("lat", "lon") %in% names(yr_data))) {
          lat <- yr_data$lat
          lon <- yr_data$lon
        } else if (all(c("latitude", "longitude") %in% names(yr_data))) {
          lat <- yr_data$latitude
          lon <- yr_data$longitude
        } else {
          # Fallback to parsing location_id
          parts <- strsplit(as.character(yr_data$location_id), "_")
          lon_str <- sapply(parts, function(x) if (length(x) >= 2) x[length(x)-1] else NA_character_)
          lat_str <- sapply(parts, function(x) if (length(x) >= 1) x[length(x)] else NA_character_)
          lon <- suppressWarnings(as.numeric(lon_str))
          lat <- suppressWarnings(as.numeric(lat_str))
        }

        # If we have NAs, try to fill them from gpts_map global if available
        if ((any(is.na(lon)) || any(is.na(lat))) && exists("gpts_map", envir = globalenv())) {
           gpts <- get("gpts_map", envir = globalenv())
           if (all(c("location_id", "lat", "lon") %in% names(gpts))) {
             missing_idx <- which(is.na(lon) | is.na(lat))
             m_ids <- as.character(yr_data$location_id[missing_idx])
             match_idx <- match(m_ids, as.character(gpts$location_id))
             lon[missing_idx] <- ifelse(is.na(match_idx), lon[missing_idx], gpts$lon[match_idx])
             lat[missing_idx] <- ifelse(is.na(match_idx), lat[missing_idx], gpts$lat[match_idx])
             # Warn about locations that still could not be matched
             still_missing <- missing_idx[is.na(match_idx)]
             if (length(still_missing) > 0) {
               warning(sprintf(
                 "[collect_location_coords] %d location(s) not found in gpts_map and remain without coordinates: %s",
                 length(still_missing),
                 paste(m_ids[is.na(match_idx)], collapse = ", ")
               ))
             }
           }
        }

        if (all(is.na(lon)) || all(is.na(lat))) {
          n_eff_est <- n_obs  # No spatial info: assume independence
        } else {
          # Compute pairwise great-circle distances (km) using Haversine
          valid <- which(!is.na(lon) & !is.na(lat))
          if (length(valid) >= 3) {
            coords <- cbind(lon[valid], lat[valid])
            coefs_valid <- yr_data$coef[valid]

            # Distance matrix in km (Haversine approximation)
            dist_mat <- compute_haversine_distance_matrix(coords)

            # Estimate spatial correlation range via empirical variogram
            # Fit exponential model: C(d) = exp(-d / range)
            # Use Moran-style approach: mean pairwise correlation as function of distance
            dists <- dist_mat[upper.tri(dist_mat)]
            coef_diffs_sq <- outer(coefs_valid, coefs_valid, function(a, b) (a - b)^2)
            gamma_vals <- coef_diffs_sq[upper.tri(coef_diffs_sq)] / 2  # semivariance

            if (length(dists) > 0 && var(coefs_valid, na.rm = TRUE) > 0) {
              total_var <- var(coefs_valid, na.rm = TRUE)
              # Estimate range by finding distance at which semivariance reaches ~63% of sill
              # Simple robust estimate: fit exponential variogram gamma(d) = sill * (1 - exp(-d/range))
              # Use median-based bins for robustness
              n_bins <- min(10, max(3, length(dists) %/% 5))
              bin_breaks <- unique(quantile(dists, probs = seq(0, 1, length.out = n_bins + 1)))
              if (length(bin_breaks) >= 3) {
                bin_mid <- (bin_breaks[-length(bin_breaks)] + bin_breaks[-1]) / 2
                bin_gamma <- numeric(length(bin_mid))
                for (bb in seq_along(bin_mid)) {
                  in_bin <- dists >= bin_breaks[bb] & dists < bin_breaks[bb + 1]
                  bin_gamma[bb] <- if (sum(in_bin) > 0) median(gamma_vals[in_bin]) else NA
                }

                range_est <- fit_exponential_variogram(bin_mid, bin_gamma, total_var, dists)
                if (!is.null(range_est)) {
                  # Compute mean pairwise correlation: C(d) = exp(-d / range)
                  # Effective n = n / (1 + (n-1) * mean_corr)  [Kish formula]
                  mean_corr <- mean(exp(-dists / range_est))
                  n_eff_est <- max(1, n_obs / (1 + (n_obs - 1) * mean_corr))
                } else {
                  n_eff_est <- n_obs
                }
              } else {
                n_eff_est <- n_obs
              }
            } else {
              n_eff_est <- n_obs
            }
          } else {
            n_eff_est <- n_obs
          }
        }
      }
      n_eff_vec[i] <- n_eff_est

      if (n_obs > 0) {
        # DECISION: Small Sample Size vs Large Sample Size
        if (n_obs < 15) {
             # --- SMALL SAMPLE: USE POOLED VARIANCE ---
             mu <- mean(yr_data$coef, na.rm = TRUE)
             n_eff_i <- n_eff_vec[i]

             if (n_obs >= 2) {
               sample_sd <- sd(yr_data$coef, na.rm = TRUE)
               if (is.na(sample_sd) || sample_sd == 0) sample_sd <- pooled_sd
               se_mean <- sample_sd / sqrt(n_eff_i)
             } else {
               se_mean <- pooled_sd / sqrt(n_eff_i)
             }

             boot_means[, i] <- rnorm(B, mean = mu, sd = se_mean)

        } else {
             # --- SUFFICIENT SAMPLE: USE RESAMPLING ---
             for (b in 1:B) {
               boot_sample <- sample(yr_data$coef, n_obs, replace = TRUE)
               boot_means[b, i] <- mean(boot_sample, na.rm = TRUE)
             }
        }
      }
    }
    
    boot_result <- data.frame(
      year = years,
      Veg = veg,
      n_locations = sapply(years, function(y) sum(veg_data$pheno_year == y & !is.na(veg_data$coef))),
      global_coef = apply(boot_means, 2, mean, na.rm = TRUE),
      se = apply(boot_means, 2, sd, na.rm = TRUE),
      coef_025 = apply(boot_means, 2, quantile, 0.025, na.rm = TRUE),
      coef_975 = apply(boot_means, 2, quantile, 0.975, na.rm = TRUE),
      n_eff = n_eff_vec,
      method = "robust_location_bootstrap"
    )

    # Print spatial effective sample size diagnostics
    cat(sprintf("[BOOTSTRAP] Veg '%s': effective sample size (n_eff) across years — mean=%.1f, median=%.1f, min=%.1f, max=%.1f (of n_obs mean=%.1f)\n",
        veg,
        mean(n_eff_vec, na.rm = TRUE), median(n_eff_vec, na.rm = TRUE),
        min(n_eff_vec, na.rm = TRUE), max(n_eff_vec, na.rm = TRUE),
        mean(sapply(years, function(y) sum(veg_data$pheno_year == y & !is.na(veg_data$coef)))))
    )

    # NOTE: No post-hoc autocorrelation inflation applied.
    # For the large-N path (n>=15), the bootstrap already captures between-location variance.
    # For the small-N path (n<15), spatial autocorrelation is accounted for via n_eff
    # in the SE calculation above (dividing by sqrt(n_eff) instead of sqrt(n_obs)).

    boot_result$coef_025 <- pmax(0, boot_result$coef_025)
    boot_result$coef_975 <- pmin(1, boot_result$coef_975)
    
    # Clamp the median to prevent negative coverages
    boot_result$global_coef <- pmax(0, boot_result$global_coef)
    
    results_list[[veg]] <- boot_result
  }
  
  dplyr::bind_rows(results_list)
}









# Helper: ensure the variant similarity heatmap is generated once and early


# Helper: load inference dataset and apply filtering only AFTER the variant similarity heatmap exists ✅
load_and_prepare_inference_data <- function() {
  # Ensure df_inf always exists (prevents 'object df_inf not found' when inference is disabled/missing)
  df_inf <- NULL

  # Clear any cached inference data; this function always reloads from CSV
  if (exists("df_tasks_inference", inherits = TRUE)) {
    df_tasks_inference <<- NULL
  }

  # Always reload inference data from CSV (do not use cached version)
  if (exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0) {
    cat("[INFO] Clearing cached inference data to force reload from CSV.\n")
  }

  # Load and filter (copied logic)
  if (exists("INFERENCE_CSV")) {
    if (isTRUE(is.na(INFERENCE_CSV)) || !nzchar(as.character(INFERENCE_CSV))) {
      cat("[INFERENCE] INFERENCE_CSV is NA/empty; skipping inference loading.\n")
      return(invisible(NULL))
    }
    if (file.exists(INFERENCE_CSV)) {
      cat(sprintf("Loading inference data from %s...\n", INFERENCE_CSV))
      df_inf <- tryCatch(read.csv(INFERENCE_CSV), error = function(e) { cat(sprintf("[WARNING] Error reading inference CSV: %s\n", e$message)); NULL })
      if (!is.null(df_inf)) {
        cat(sprintf("Loaded %d rows from inference file.\n", nrow(df_inf)))
        
        df_inf <- canonicalize_veg_labels(df_inf)
        
        # IMPORTANT: Reconstruct location_id from lat/lon coordinates
        # Do NOT use the existing location_id column (if present) as it may be incorrect
        if ("location_id" %in% names(df_inf)) {
          df_inf$location_id <- NULL  # Remove to force reconstruction
          cat("[NOTICE] Removed existing 'location_id' column - will reconstruct from lat/lon\n")
        }
        
        # Normalize coordinate column names
        lon_candidates <- names(df_inf)[grepl("(^|_)(lon|longitude|x)(_|$)", names(df_inf), ignore.case = TRUE)]
        lat_candidates <- names(df_inf)[grepl("(^|_)(lat|latitude|y)(_|$)", names(df_inf), ignore.case = TRUE)]
        
        if (length(lon_candidates) > 0 && length(lat_candidates) > 0) {
          lon_col <- lon_candidates[1]
          lat_col <- lat_candidates[1]
          
          # Ensure numeric
          df_inf[[lon_col]] <- as.numeric(df_inf[[lon_col]])
          df_inf[[lat_col]] <- as.numeric(df_inf[[lat_col]])
          
          # Standardize to 'lon' and 'lat' column names
          if (lon_col != "lon") {
            df_inf$lon <- df_inf[[lon_col]]
            cat(sprintf("[NOTICE] Using '%s' as longitude\n", lon_col))
          }
          if (lat_col != "lat") {
            df_inf$lat <- df_inf[[lat_col]]
            cat(sprintf("[NOTICE] Using '%s' as latitude\n", lat_col))
          }
          
          # Reconstruct location_id from coordinates
          df_inf$location_id <- make_location_id(df_inf$lon, df_inf$lat)
          cat(sprintf("[NOTICE] Reconstructed location_id from lat/lon: %d unique locations\n",
                      length(unique(df_inf$location_id[!is.na(df_inf$location_id)]))))

          # Apply MAX_INFERENCE_LOCATIONS limit once per inference CSV file
          if (exists("MAX_INFERENCE_LOCATIONS") && is.numeric(MAX_INFERENCE_LOCATIONS)) {
            unique_locations <- unique(df_inf$location_id[!is.na(df_inf$location_id)])
            n_unique <- length(unique_locations)
            inf_file_label <- if (exists("INFERENCE_CSV") && nzchar(as.character(INFERENCE_CSV))) basename(as.character(INFERENCE_CSV)) else "<unknown file>"

            if (n_unique > MAX_INFERENCE_LOCATIONS) {
              set.seed(get_mesma_seed(123))  # Deterministic sampling for reproducibility
              sampled_locations <- sample(unique_locations, MAX_INFERENCE_LOCATIONS, replace = FALSE)

              # Filter the dataframe to only include sampled locations
              df_inf <- df_inf[df_inf$location_id %in% sampled_locations, ]

              cat(sprintf("[INFERENCE LOADING] %s: reduced from %d to %d locations (MAX_INFERENCE_LOCATIONS=%d per inference file)\n",
                          inf_file_label, n_unique, MAX_INFERENCE_LOCATIONS, MAX_INFERENCE_LOCATIONS))
              cat(sprintf("[INFERENCE LOADING] Filtered dataset: %d rows from %d locations\n",
                          nrow(df_inf), length(unique(df_inf$location_id[!is.na(df_inf$location_id)]))))
            } else {
              cat(sprintf("[INFERENCE LOADING] %s: using all %d locations (<= MAX_INFERENCE_LOCATIONS=%d per inference file)\n",
                          inf_file_label, n_unique, MAX_INFERENCE_LOCATIONS))
            }
          }

          # Add Veg from gpts_map if available
          if (exists("gpts_map") && "Veg" %in% names(gpts_map)) {
            gpts_map$location_id <- as.character(gpts_map$location_id)
            df_inf$location_id <- as.character(df_inf$location_id)
            df_inf <- df_inf |>
              dplyr::left_join(
                gpts_map |> dplyr::select(location_id, Veg),
                by = "location_id"
              )
            cat(sprintf("[NOTICE] Added Veg column to inference data from gpts_map.\n"))
          }
        } else {
          stop("ERROR: Inference CSV must contain latitude and longitude columns for location_id reconstruction")
        }
      }
    } else {
      stop(sprintf("[INFERENCE] Inference file not found at: %s", INFERENCE_CSV))
    }
  } else {
    stop("[INFERENCE] INFERENCE_CSV variable not defined.")
  }

  if (!is.null(df_inf) && nrow(df_inf) > 0) {
    df_inf <- normalize_band_names(df_inf)
    if ("...1" %in% names(df_inf) && !"location_id" %in% names(df_inf)) {
      if (is.character(df_inf$...1) || is.numeric(df_inf$...1)) names(df_inf)[names(df_inf) == "...1"] <- "location_id"
    }

    # derive date column
    if ("prediction_date" %in% names(df_inf)) df_inf$date <- as.Date(df_inf$prediction_date) else if ("date" %in% names(df_inf)) df_inf$date <- as.Date(df_inf$date) else {
      for (col in names(df_inf)) {
        if (inherits(df_inf[[col]], "Date")) { df_inf$date <- df_inf[[col]]; break }
        if (is.character(df_inf[[col]]) && all(grepl("^\\d{4}-\\d{2}-\\d{2}", na.omit(df_inf[[col]][1:min(10, nrow(df_inf))])))) { df_inf$date <- as.Date(df_inf[[col]]); break }
      }
    }

    if ("location_id" %in% names(df_inf) && "date" %in% names(df_inf)) {
      df_inf$location_id <- as.character(df_inf$location_id)

      # =======================================================================
      # SENSOR BIAS CORRECTION (inference) — mirrors preprocess_data.R block.
      # Harmonises LANDSAT_89 (OLI) raw bands to the ETM+ radiometric scale
      # using per-band affine coefficients (slope + intercept) so inference
      # features are on the same scale as the bias-corrected training prototypes.
      # =======================================================================
      {
        df_inf <- apply_oli_etm_bias_correction(df_inf, dataset_label = "inference", log_prefix = "[BIAS CORR INF]")
      }
      # =======================================================================

      before_cols <- names(df_inf)
      df_inf <- compute_indices_from_bands(df_inf)
      new_cols <- setdiff(names(df_inf), before_cols)
      if (length(new_cols) > 0) cat(sprintf("[NOTICE] Computed indices from raw bands for inference: %s\n", paste(new_cols, collapse=", ")))

      if ("NDDI" %in% names(df_inf)) {
        dust_count <- sum(df_inf$NDDI > NDDI_DUST_THRESHOLD, na.rm = TRUE)
        total_before <- nrow(df_inf)
        df_inf <- df_inf[!(df_inf$NDDI > NDDI_DUST_THRESHOLD), , drop = FALSE]
        total_after <- nrow(df_inf)
        filtered <- total_before - total_after
        cat(sprintf("[INFERENCE FILTERING] Filtered out %d observations with dust (NDDI > %s) contamination\n", filtered, .nddi_thresh_fmt()))
        cat(sprintf("[INFERENCE FILTERING] Inference dataset after contamination filtering: %d rows from %d locations\n", total_after, length(unique(df_inf$location_id))))
        df_inf <- remove_large_outliers(df_inf)
      } else {
        cat("[WARNING] NDDI not found in inference data; skipping contamination filtering\n")
        # Indices were computed above; continue.
      }

      if (!"Veg" %in% names(df_inf)) df_inf$Veg <- NA_character_
      if (!"pheno_year" %in% names(df_inf)) df_inf$pheno_year <- assign_pheno_year(df_inf$date)
      if (!"doy" %in% names(df_inf)) df_inf$doy <- pheno_doy(df_inf$date)
      if (!"zenith.angle" %in% names(df_inf)) df_inf$zenith.angle <- NA_real_
      df_inf$DVI_max <- 0.7

      # IMPORTANT: Use ALL years for inference - inference should cover the full temporal range
      if (exists("TRAIN_YEARS") && !is.null(TRAIN_YEARS)) {
        cat(sprintf("[NOTICE] TRAIN_YEARS is set to %s for training, but inference uses ALL available years.\n", paste(TRAIN_YEARS, collapse=", ")))
        cat(sprintf("[NOTICE] Inference dataset has %d rows from %d locations (all years)\n", nrow(df_inf), length(unique(df_inf$location_id))))
        cat("[NOTICE] Trend computations will use this inference dataset (INFERENCE_CSV) and will NOT use training data.\n")
      }

      # -----------------------------------------------------------------------------
      # TREND (INFERENCE DATA ONLY): compute per-pheno_year means for key indices
      # Do not use training data for trend computation — use the INFERENCE_CSV file provided
      # -----------------------------------------------------------------------------
      tryCatch({
        # only track PPI trends
        trend_indices <- c("PPI")
        trend_indices <- intersect(trend_indices, names(df_inf))
        # Detrend PPI using saved seasonal model if available
        if ("PPI" %in% trend_indices && exists("INDEX_SEASONAL_MODELS") && exists("INDEX_SEASONAL_MEANS")) {
          if ("PPI" %in% names(INDEX_SEASONAL_MODELS)) {
            seasonal_model <- INDEX_SEASONAL_MODELS[["PPI"]]
            global_seasonal_mean <- INDEX_SEASONAL_MEANS[["PPI"]]
            summer_inf <- df_inf |> dplyr::filter(lubridate::month(date) %in% c(6,7,8,9))
            if (nrow(summer_inf) > 0) {
              summer_inf$.orig_row <- seq_len(nrow(df_inf))[which(lubridate::month(df_inf$date) %in% c(6,7,8,9))]
              summer_inf$doy <- lubridate::yday(summer_inf$date)
              summer_inf$seasonal_trend <- predict(seasonal_model, newdata = summer_inf)
              summer_inf$ppi_detrended <- summer_inf$PPI - (summer_inf$seasonal_trend - global_seasonal_mean)
              # assign back
              df_inf$ppi_detrended <- NA_real_
              df_inf$ppi_detrended[summer_inf$.orig_row] <- summer_inf$ppi_detrended
              if (!"ppi_detrended" %in% trend_indices) trend_indices <- c(trend_indices, "ppi_detrended")
            } else {
              df_inf$ppi_detrended <- df_inf$PPI
              if (!"ppi_detrended" %in% trend_indices) trend_indices <- c(trend_indices, "ppi_detrended")
            }
          } else {
            df_inf$ppi_detrended <- df_inf$PPI
            if (!"ppi_detrended" %in% trend_indices) trend_indices <- c(trend_indices, "ppi_detrended")
          }
        }

        # Aggregate summer (Jun-Sep) per pheno_year
        summer_inf_all <- df_inf |> dplyr::filter(lubridate::month(date) %in% c(6,7,8,9))
        if (nrow(summer_inf_all) > 0 && length(trend_indices) > 0) {
          sum_by_year <- summer_inf_all |> dplyr::group_by(pheno_year) |>
            dplyr::summarise(across(dplyr::all_of(trend_indices), ~ mean(.x, na.rm = TRUE)), n = dplyr::n(), .groups = "drop")

          if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
          # Use fixed filenames (no date/time) so outputs are deterministic and do not embed system timestamps
          csvfile <- file.path(OUT_DIR, "inference_trend_summary.csv")
          tryCatch({ utils::write.csv(sum_by_year, file = csvfile, row.names = FALSE); cat(sprintf("[TREND] Saved inference trend summary CSV: %s\n", csvfile)) }, error = function(e) cat(sprintf("[TREND] Failed to save inference trend CSV: %s\n", e$message)))

          # Plot trends
          library(ggplot2)
          
          plot_dt <- tidyr::pivot_longer(sum_by_year, cols = dplyr::starts_with("ppi"), names_to = "index", values_to = "mean_val")
          if (nrow(plot_dt) > 0) {
            p_trend <- ggplot(plot_dt, aes(x = pheno_year, y = mean_val)) +
              add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
              geom_line() + geom_point() + facet_wrap(~index, scales = "free_y") + theme_minimal() +
              scale_x_continuous(limits = c(1984, NA)) +
              coord_cartesian(ylim = c(0, NA))
            # Use fixed filename (no date/time)
            pfile <- file.path(OUT_DIR, "inference_trends.png")
            tryCatch({ ggplot2::ggsave(pfile, plot = p_trend, width = 10, height = 6); cat(sprintf("[TREND] Saved inference trend plot: %s\n", pfile)) }, error = function(e) cat(sprintf("[TREND] Failed to save trend plot: %s\n", e$message)))
          }
        } else {
          cat("[TREND] No summer observations in inference data or no indices available; skipping trend summary\n")
        }
      }, error = function(e) {
        cat(sprintf("[TREND] Error during inference trend computation: %s\n", e$message))
      })

      # Determine available indices from training parameters
      if (exists("TRAINING_NORM_PARAMS") && !is.null(TRAINING_NORM_PARAMS) && !is.null(TRAINING_NORM_PARAMS$INDEX_SCALES)) {
        avail <- names(TRAINING_NORM_PARAMS$INDEX_SCALES)
      } else {
        avail <- OPTIMAL_INDICES
      }

      # Ensure avail includes all indices required by the trained model (MESMA_PARAMS)
      if (exists("MESMA_PARAMS") && !is.null(MESMA_PARAMS) && !is.null(MESMA_PARAMS$indices)) {
        avail <- unique(c(avail, MESMA_PARAMS$indices))
      }

      dvi_soil_arg <- NA_real_
      
      # Defensive: ensure `ppi_inf_res` exists so downstream checks do not error when
      # the auto-add branch is skipped (was causing "object 'ppi_inf_res' not found").
      ppi_inf_res <- NULL
      
      # If PPI not in df_inf, try to add it
      if (exists("add_ppi_columns") && "PPI" %in% avail && !"PPI" %in% names(df_inf)) {
        dvi_soil_vec_inf <- compute_dvi_soil_per_location(df_inf, quantile_p = 0.10)
        df_inf <- add_ppi_columns(df_inf, dvi_soil = dvi_soil_vec_inf)
        cat("[PPI] Added PPI to inference data (per-location dvi_soil + per-location M).\n")
      }

      if (!"PPI" %in% names(df_inf) || all(!is.finite(df_inf$PPI))) {
        stop("[PPI] Missing or all PPI values non-finite in inference data after preprocessing")
      }

      # One-time feature-space check + fill optional cols with NA
      required_indices <- if (exists("MESMA_PARAMS") && !is.null(MESMA_PARAMS) && !is.null(MESMA_PARAMS$indices)) MESMA_PARAMS$indices else avail
      df_inf <- mesma_prepare_feature_columns(df_inf, required_cols = required_indices, optional_cols = avail, context = "inference")

      if ("prediction_date" %in% names(df_inf)) df_inf$prediction_date <- as.Date(df_inf$prediction_date)
      if ("reference_date" %in% names(df_inf)) df_inf$reference_date <- as.Date(df_inf$reference_date)

      # Backup raw PPI before normalization (inference)
      if ("PPI" %in% names(df_inf) && !"PPI_raw" %in% names(df_inf)) df_inf$PPI_raw <- df_inf$PPI


      # (nir/red backups are left in place in case other routines require them)
      if ("nir" %in% names(df_inf) && !"nir_raw" %in% names(df_inf)) df_inf$nir_raw <- df_inf$nir
      if ("red" %in% names(df_inf) && !"red_raw" %in% names(df_inf)) df_inf$red_raw <- df_inf$red

      if (exists("TRAINING_NORM_PARAMS") && !is.null(TRAINING_NORM_PARAMS)) { 
        cat("\n=== APPLYING STORED NORMALIZATION TO INFERENCE DATA ===\n")
        df_inf <- apply_stored_normalization(df_inf, TRAINING_NORM_PARAMS, cols = avail, lat_default = 40.2)
        cat("=======================================================\n") 
      } else { 
        warning("TRAINING_NORM_PARAMS not found, applying fresh normalization to inference data (scale factors may differ from training!)")
        inf_norm_result <- normalize_mesma_data(df_inf, cols = avail, lat_default = 40.2)
        df_inf <- inf_norm_result$df 
      }

      df_inf <- backup_and_normalize_ppi(df_inf, label = "Inference")

      df_tasks_inference <<- df_inf

      # Compute average value of first pentad of first index for inference data
      if (exists("MESMA_PARAMS") && !is.null(MESMA_PARAMS$indices) && length(MESMA_PARAMS$indices) > 0) {
        first_index <- MESMA_PARAMS$indices[1]
        col_name <- paste0(first_index, "_1")
        if (col_name %in% names(df_inf)) {
          avg_inference <- mean(df_inf[[col_name]], na.rm = TRUE)
          cat(sprintf("[INFO] Average value of first pentad of first index (inference, after normalization): %.6f\n", avg_inference))
        } else {
          cat("[WARNING] Column for first pentad of first index not found in inference data\n")
        }
      }

      n_infer_loc_years <- nrow(unique(df_inf[c("location_id", "pheno_year")] ))
      if (exists("df_train") && !is.null(df_train) && nrow(df_train) > 0) {
          if (!"pheno_year" %in% names(df_train) && "date" %in% names(df_train))
              df_train$pheno_year <- assign_pheno_year(df_train$date)
          n_train_loc_years <- nrow(unique(df_train[c("location_id", "pheno_year")]))
      } else {
          n_train_loc_years <- 0
      }
      cat("(NOTICE) Inference dataset location-years:", n_infer_loc_years, "\n")
      cat("(NOTICE) Training dataset location-years:", n_train_loc_years, "\n")
      if (n_train_loc_years > 0 && n_infer_loc_years == n_train_loc_years)
          cat("(WARNING) Training and inference datasets have the same number of location-years ",
              sprintf("(%d). This may be expected if IDs are independent; no automatic filtering will be applied.\n", n_train_loc_years))

      cat("Keeping df_tasks as training data for separate processing.\n")
    } else {
      cat("[WARNING] Inference data missing 'location_id' or 'date' column. Skipping.\n")
      cat(sprintf("Columns found: %s\n", paste(names(df_inf), collapse=", ")))
    }
  }
  df_inf_state <- if (exists("df_tasks_inference", envir = globalenv())) df_tasks_inference else NULL
  inf_ids <- if (!is.null(df_inf_state) && is.data.frame(df_inf_state) && nrow(df_inf_state) > 0 && "location_id" %in% names(df_inf_state)) unique(df_inf_state$location_id) else character(0)
  inference_location_ids <<- inf_ids
  INFERENCE_LOAD_DEFERRED <<- FALSE
}



  # =============================================================================
  # FCLS: Fully Constrained Least Squares Solver
  # Enforces BOTH non-negativity AND sum-to-one constraints simultaneously
  # Uses iterative active-set method (Heinz & Chang, 2001)
  # =============================================================================
  
  solve_weights_fcls <- function(E, y, feature_weights = NULL, max_iter = 500, tol = 1e-8) {
    if (is.null(E) || ncol(E) < 1) return(NULL)

    # Ensure numeric
    E_fit <- as.matrix(E)
    y_fit <- as.numeric(y)
    n_endmembers <- ncol(E_fit)
    n_bands <- nrow(E_fit)

    # Handle length mismatch
    if (length(y_fit) != n_bands) {
      if (length(y_fit) > n_bands) y_fit <- y_fit[1:n_bands] else y_fit <- c(y_fit, rep(0, n_bands - length(y_fit)))
    }

    # Mask missing values instead of just zero-padding
    valid <- is.finite(y_fit)
    if (!all(valid)) {
      y_fit[!valid] <- 0
      E_fit[!valid, ] <- 0
    }
    E_fit[!is.finite(E_fit)] <- 0

    # Base feature weights (from PCA-LDA or uniform)
    if (!is.null(feature_weights) && length(feature_weights) == n_bands) {
      feature_weights <- as.numeric(feature_weights)
      feature_weights[!is.finite(feature_weights)] <- 0
      base_weights <- pmax(feature_weights, 0)
    } else {
      base_weights <- rep(1, n_bands)
    }

    # === AUGMENTED NNLS (FCLS) — Heinz & Chang 2001 ===
    # Problem: min ||Ex - y||^2  s.t. x >= 0, sum(x) = 1
    # Enforce sum-to-one by augmenting the system with a scaled row:
    #   E_aug = [E_w; delta * 1^T],  y_aug = [y_w; delta]
    # then solve: min ||E_aug x - y_aug||^2  s.t. x >= 0
    # Normalise afterward to guarantee sum(x) = 1 exactly.
      E_w <- E_fit * base_weights
      y_w <- y_fit * base_weights

      # delta: scale to match spectral magnitude so the constraint row is
      # influential but not overwhelming (sqrt of mean squared E value * 100).
      delta <- sqrt(mean(E_w^2)) * 100
      if (!is.finite(delta) || delta < 1e-8) delta <- 1.0

      E_aug <- rbind(E_w, delta * rep(1, n_endmembers))
      y_aug <- c(y_w, delta)

      res_nnls <- nnls::nnls(E_aug, y_aug)
      w <- res_nnls$x
      w[!is.finite(w)] <- 0
      w_sum <- sum(w)
      if (w_sum > 0) w <- w / w_sum else w <- rep(1 / n_endmembers, n_endmembers)

      # Calculate RMSE
      pred  <- as.numeric(E_fit %*% w)
      resid <- y_fit - pred
      rmse  <- sqrt(mean(resid^2))

      return(list(w = w, rmse = rmse, residuals = resid, loss_type = "rmse"))
    }
  

  # Batch FCLS Solver using Augmented NNLS (Optimization for GA/Grid Search)
  # E: Features x Endmembers matrix
  # Y: Samples x Features matrix
  # Returns: Samples x Endmembers weight matrix
  solve_batch_fcls <- function(E, Y,
                                  feature_weights = NULL,
                                  lda_basis = NULL,
                                  lda_component_weights = NULL,
                                  solver = "fcls",
                                  iwlmm_iters = 3) {
    # Ensure matrices
    E <- as.matrix(E)
    Y <- as.matrix(Y)
    if (!is.null(lda_basis)) lda_basis <- as.matrix(lda_basis)
    if (!is.null(lda_component_weights)) lda_component_weights <- as.numeric(lda_component_weights)
    # possible solvers: "fcls" (default), "sparse", "iwlmm"
    # lda options as before
    if (!is.null(lda_basis)) {
      # save original for missing-data reliability
      Y_orig <- Y

      # Ensure LDA basis rows align with feature dimension used by E and Y.
      # Mismatches can occur after feature pruning/masking.
      n_feat_E <- nrow(E)
      n_feat_Y <- ncol(Y)
      n_feat_LDA <- nrow(lda_basis)
      use_lda_projection <- TRUE

      if (!identical(n_feat_E, n_feat_Y)) {
        warning(sprintf("solve_batch_fcls: E/Y feature mismatch before LDA projection (E rows=%d, Y cols=%d); disabling LDA projection for this call", n_feat_E, n_feat_Y))
        use_lda_projection <- FALSE
      } else if (!identical(n_feat_LDA, n_feat_E)) {
        n_common <- min(n_feat_LDA, n_feat_E)
        if (n_common < 1) {
          warning("solve_batch_fcls: No common feature dimensions between lda_basis and E/Y; disabling LDA projection for this call")
          use_lda_projection <- FALSE
        } else {
          warning(sprintf("solve_batch_fcls: lda_basis/E/Y dimension mismatch (lda=%d, E=%d, Y=%d); trimming to first %d features", n_feat_LDA, n_feat_E, n_feat_Y, n_common))
          lda_basis <- lda_basis[seq_len(n_common), , drop = FALSE]
          E <- E[seq_len(n_common), , drop = FALSE]
          Y <- Y[, seq_len(n_common), drop = FALSE]
          Y_orig <- Y_orig[, seq_len(n_common), drop = FALSE]
        }
      }

      if (use_lda_projection) {
        E_safe <- E
        E_safe[!is.finite(E_safe)] <- 0
        lda_basis[!is.finite(lda_basis)] <- 0
        E <- t(lda_basis) %*% E_safe
        Y_t <- t(Y)
        Y_t_safe <- Y_t
        Y_t_safe[!is.finite(Y_t_safe)] <- 0
        Y_t <- t(lda_basis) %*% Y_t_safe
        Y <- t(Y_t)
        if (!is.null(lda_component_weights)) {
          feature_weights <- lda_component_weights
        } else {
          feature_weights <- NULL
        }
        # compute reliability weights whenever any input features are missing
        # (this is mandatory when projecting to LDA space)
        if (any(!is.finite(Y_orig))) {
          absA <- abs(lda_basis)
          denom <- colSums(absA)
          # avoid zero denom
          denom[denom == 0] <- 1
          # construct matrix of reliabilities
          rel_mat <- matrix(0, nrow = nrow(Y_orig), ncol = ncol(absA))
          for (k in seq_len(nrow(Y_orig))) {
            valid <- is.finite(Y_orig[k, ])
            if (any(valid)) {
              rel_mat[k, ] <- colSums(absA[valid, , drop = FALSE]) / denom
            } else {
              rel_mat[k, ] <- 0
            }
          }
          # apply to projected data
          Y <- Y * rel_mat
        }
      }
    }

    if (is.null(E) || ncol(E) < 1) return(NULL)
    n_endmembers <- ncol(E)
    n_samples <- nrow(Y)
    if (n_samples == 0) return(matrix(0, 0, n_endmembers))

    # Convert Y to Features x Samples for matrix math
    Y_t <- t(Y)

    # helper: run augmented NNLS (FCLS) and return weight matrix
    run_fcls <- function(E_fit, Y_fit, wts) {
      # delta: match spectral magnitude for sum-to-one enforcement
      delta <- sqrt(mean(E_fit^2, na.rm=TRUE)) * 100
      if (!is.finite(delta) || delta < 1e-8) delta <- 1.0
      aug_row <- delta * rep(1, n_endmembers)
      
      w_out <- matrix(0, nrow = n_samples, ncol = n_endmembers)
      
      for (i in seq_len(n_samples)) {
        y_i <- Y_fit[, i]
        
        # Determine valid features (non-NA in both y and E)
        valid_y <- is.finite(y_i)
        valid_E <- apply(E_fit, 1, function(x) all(is.finite(x)))
        valid <- valid_y & valid_E
        
        E_i <- E_fit
        E_i[!is.finite(E_i)] <- 0 # Impute any residual before masking so nnls doesn't crash on NAs in unmasked parts
        
        if (!all(valid)) {
          y_i[!valid] <- 0
          E_i[!valid, ] <- 0
        }
        
        E_aug <- rbind(E_i, aug_row)
        y_aug <- c(y_i, delta)
        
        res   <- nnls::nnls(E_aug, y_aug)
        w     <- res$x
        w[!is.finite(w)] <- 0
        s <- sum(w)
        if (s > 0) w <- w / s else w <- rep(1 / n_endmembers, n_endmembers)
        w_out[i, ] <- w
      }
      w_out
    }

    # determine E_fit and Y_fit based on feature_weights
    if (!is.null(feature_weights) && length(feature_weights) == nrow(E)) {
      w <- pmax(feature_weights, 0); w[!is.finite(w)] <- 0
      E_fit_base <- E * w
      Y_fit_base <- Y_t * w
    } else {
      E_fit_base <- E
      Y_fit_base <- Y_t
    }

    run_iwlmm <- function(E_in, Y_in_t, fw) {
      Wmat <- run_fcls(E_in, Y_in_t, fw)
      for (iter in 1:iwlmm_iters) {
        wgt <- 1 / (abs(Wmat) + 1e-6)
        feat_w <- colMeans(wgt)
        if (!is.null(fw)) feat_w <- feat_w * fw
        if (exists("IWLMM_FEAT_W_MIN") && !is.na(IWLMM_FEAT_W_MIN)) feat_w <- pmax(feat_w, IWLMM_FEAT_W_MIN)
        if (exists("IWLMM_FEAT_W_MAX") && !is.na(IWLMM_FEAT_W_MAX)) feat_w <- pmin(feat_w, IWLMM_FEAT_W_MAX)
        E_iter <- E * feat_w
        Y_iter <- Y_t * feat_w
        Wmat <- run_fcls(E_iter, Y_iter, feat_w)
      }
      Wmat
    }

    sparsify <- function(Wmat) {
      k <- if (exists("SPARSE_UNMIX_K", inherits = TRUE)) SPARSE_UNMIX_K else 2L
      em_names <- colnames(Wmat)
      barren_col <- if (!is.null(em_names)) which(tolower(em_names) == "barren") else integer(0)
      veg_cols <- setdiff(seq_len(ncol(Wmat)), barren_col)
      for (i in 1:nrow(Wmat)) {
        wi <- Wmat[i, ]
        if (length(veg_cols) > k) {
          veg_order <- order(wi[veg_cols], decreasing = TRUE)
          drop <- veg_cols[veg_order[(k + 1):length(veg_cols)]]
          wi[drop] <- 0
        }
        s <- sum(wi)
        if (s > 0) wi <- wi / s
        Wmat[i, ] <- wi
      }
      Wmat
    }

    if (solver == "fcls") {
      return(run_fcls(E_fit_base, Y_fit_base, feature_weights))
    } else if (solver == "sparse") {
      Wmat <- run_fcls(E_fit_base, Y_fit_base, feature_weights)
      return(sparsify(Wmat))
    } else if (solver == "iwlmm") {
      return(run_iwlmm(E_fit_base, Y_fit_base, feature_weights))
    } else if (solver == "iwlmm_sparse") {
      # IWLMM iterations first, then sparsify to top-k endmembers
      Wmat <- run_iwlmm(E_fit_base, Y_fit_base, feature_weights)
      return(sparsify(Wmat))
    } else {
      stop(sprintf("Unknown solver '%s' requested", solver))
    }
  }






  compute_diagnostics <- function(y, E, w, mesma_result = NULL) {
    if (!is.matrix(E) || ncol(E) == 0 || length(w) != ncol(E)) {
      return(data.frame(
        condition_number = NA_real_,
        residual_sum_of_squares = NA_real_,
        r_squared = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    
    pred <- as.numeric(E %*% w)
    residuals <- y - pred
    rss <- sum(residuals^2)
    tss <- sum((y - mean(y))^2)
    r_squared <- if (tss > 0) 1 - rss/tss else 0
    
    cond_num <- tryCatch({
      s <- svd(E)$d
      if (length(s) > 0 && min(s) > 0) max(s) / min(s) else stop("compute_diagnostics: singular values invalid for condition number calculation")
    }, error = function(e) stop(sprintf("compute_diagnostics: failed computing SVD for condition number: %s", e$message)))
    
    data.frame(
      condition_number = cond_num,
      residual_sum_of_squares = rss,
      r_squared = r_squared,
      stringsAsFactors = FALSE
    )
  }

  cat("Building MESMA endmember library...\n")


  

  # ==========================================================================
  # Validation and inference datasets already created during stratified split
  # Validation pipeline reconstructed in the dedicated section below

  cat("\n=== DATA DISTRIBUTION ANALYSIS ===\n")
  if (exists("df_tasks") && nrow(df_tasks) > 0) {
    cat(sprintf("Total locations in df_tasks: %d\n", length(unique(df_tasks$location_id))))

    sample_sizes <- df_tasks |> 
      dplyr::group_by(location_id, pheno_year) |> 
      dplyr::summarize(n_obs = dplyr::n(), .groups = "drop")

    # Also compute distribution excluding barren (so min sample size is not driven by barren rows)
    sample_sizes_no_barren <- NULL
    if ("Veg" %in% names(df_tasks)) {
      sample_sizes_no_barren <- df_tasks |>
        dplyr::filter(!is.na(.data$Veg) & normalize_veg_name(.data$Veg) != "barren") |>
        dplyr::group_by(location_id, pheno_year) |>
        dplyr::summarize(n_obs = dplyr::n(), .groups = "drop")
    }

    if (nrow(sample_sizes) > 0) {
      sample_sizes$n_obs <- as.numeric(sample_sizes$n_obs)
      cat("\nObservations per location-year distribution:\n")
      cat(sprintf("  Min:    %d\n", as.integer(min(sample_sizes$n_obs, na.rm = TRUE))))
      cat(sprintf("  Q1:     %d\n", as.integer(stats::quantile(sample_sizes$n_obs, 0.25, na.rm = TRUE))))
      cat(sprintf("  Median: %d\n", as.integer(median(sample_sizes$n_obs, na.rm = TRUE))))
      cat(sprintf("  Q3:     %d\n", as.integer(stats::quantile(sample_sizes$n_obs, 0.75, na.rm = TRUE))))
      cat(sprintf("  Max:    %d\n", as.integer(max(sample_sizes$n_obs, na.rm = TRUE))))
      cat(sprintf("  Mean:   %.1f\n", mean(sample_sizes$n_obs, na.rm = TRUE)))

      if (!is.null(sample_sizes_no_barren) && nrow(sample_sizes_no_barren) > 0) {
        sample_sizes_no_barren$n_obs <- as.numeric(sample_sizes_no_barren$n_obs)
        cat("\nObservations per location-year distribution (excluding barren):\n")
        cat(sprintf("  Min (no barren): %d\n", as.integer(min(sample_sizes_no_barren$n_obs, na.rm = TRUE))))
      }

      problem_threshold <- MIN_OBS_PER_LOC_YEAR  # Minimum observations for meaningful MESMA
      n_problem <- sum(sample_sizes$n_obs < problem_threshold, na.rm = TRUE)
      cat(sprintf("\nLocation-years with < %d observations: %d (%.1f%%)\n",
                  problem_threshold, n_problem, 100 * n_problem / nrow(sample_sizes)))
      if (n_problem > 0) {
        cat("\nSample of problematic location-years:\n")
        problem_cases <- sample_sizes[sample_sizes$n_obs < problem_threshold, , drop = FALSE]
        print(utils::head(problem_cases, 10))
      }

      if ("doy" %in% names(df_tasks)) {
        doy_coverage <- df_tasks |> 
          dplyr::group_by(location_id, pheno_year) |> 
          dplyr::summarize(
            n_unique_doys = length(unique(doy[!is.na(doy)])),
            doy_span = if (any(!is.na(doy))) as.integer(max(doy, na.rm = TRUE) - min(doy, na.rm = TRUE)) else 0,
            .groups = "drop"
          )
        cat("\nDOY coverage per location-year:\n")
        cat(sprintf("  Median unique DOYs: %d\n", as.integer(median(doy_coverage$n_unique_doys, na.rm = TRUE))))
        cat(sprintf("  Median DOY span:    %d days\n", as.integer(median(doy_coverage$doy_span, na.rm = TRUE))))
        n_poor_coverage <- sum(doy_coverage$n_unique_doys < 30, na.rm = TRUE)
        cat(sprintf("  Location-years with < 30 unique DOYs: %d (%.1f%%)\n",
                    n_poor_coverage, 100 * n_poor_coverage / nrow(doy_coverage)))
      }
    } else {
      cat("No location-year sample sizes available to summarize.\n")
    }
  }

  # NOTE: location_list check moved to after validation split creates it (around line 5705)





  evaluate_all_combinations <- function(
    y,
    top_variants,
    lambda = 0,
    early_stop_rmse = 0,
    feature_weights = NULL,
    solver = "fcls"
  ) {
    # lambda: additional penalty per active component (used for sparse unmixing).
    # The penalty will be applied as rmse + lambda * n_active_components, where
    # active components are coefficients > SPARSE_COEF_THRESHOLD.  Barren is
    # explicitly excluded from this count so that sparsity incentives act only
    # on vegetation classes.  A lambda of zero disables the penalty.
    if (length(top_variants) == 0) return(NULL)

    full_veg_names <- names(top_variants)
    n_veg_full <- length(full_veg_names)
    if (n_veg_full == 0) return(NULL)

    y_target <- y

    score_from_solution <- function(res) {
      if (is.null(res) || is.null(res$residuals)) return(Inf)
      base_rmse <- as.numeric(res$rmse)
      if (!is.null(lambda) && lambda > 0) {
        coefs <- res$w
        n_active <- sum(coefs > 0, na.rm = TRUE)
        # barren + first veg endmember are free; penalty starts from the second veg endmember
        penalized <- max(0L, n_active - 2L)
        base_rmse <- base_rmse + lambda * penalized
      }
      base_rmse
    }

    # Solve a combination of variant indices for a SUBSET of vegetation types
    # veg_subset: character vector of vegetation type names to include
    # variant_indices: integer vector of variant indices (one per veg type in subset)
    solve_combo_subset <- function(veg_subset, variant_indices) {
      cols <- list()
      ids <- character(length(veg_subset))
      names(ids) <- veg_subset

      for (v_idx in seq_along(veg_subset)) {
        v <- veg_subset[v_idx]
        idx <- variant_indices[v_idx]
        cand <- top_variants[[v]][[idx]]
        if (is.null(cand) || is.null(cand$vec)) return(NULL)
        cols[[length(cols) + 1]] <- as.numeric(cand$vec)
        ids[v] <- cand$id
      }
      if (length(cols) == 0) return(NULL)

      E <- do.call(cbind, cols)
      if (is.null(E) || ncol(E) < 1) return(NULL)

      fw_combo <- if (!is.null(feature_weights) && length(feature_weights) == nrow(E)) feature_weights else NULL
      res <- solve_weights_fcls(E, y_target, feature_weights = fw_combo)
      if (is.null(res)) return(NULL)

      # Apply IWLMM iterations if requested
      if (solver %in% c("iwlmm", "iwlmm_sparse")) {
        n_iter <- if (exists("IWLMM_MAX_ITER", inherits = TRUE)) IWLMM_MAX_ITER else 5
        w_cur <- res$w
        for (iter_i in seq_len(n_iter)) {
          w_prev <- w_cur
          feat_w <- 1 / (abs(w_cur) + 1e-6)
          if (!is.null(fw_combo)) feat_w <- feat_w * fw_combo[seq_along(feat_w)]
          if (exists("IWLMM_FEAT_W_MIN") && !is.na(IWLMM_FEAT_W_MIN)) feat_w <- pmax(feat_w, IWLMM_FEAT_W_MIN)
          if (exists("IWLMM_FEAT_W_MAX") && !is.na(IWLMM_FEAT_W_MAX)) feat_w <- pmin(feat_w, IWLMM_FEAT_W_MAX)
          E_iter <- E * feat_w
          res_iter <- solve_weights_fcls(E_iter, y_target, feature_weights = NULL)
          if (is.null(res_iter)) break
          w_cur <- res_iter$w
          tol_val <- if (exists("IWLMM_TOL", inherits = TRUE)) IWLMM_TOL else 1e-4
          if (max(abs(w_cur - w_prev)) < tol_val) break
        }
        # recompute residuals with original (unweighted) E
        pred <- as.numeric(E %*% w_cur)
        resid <- y_target - pred
        res$w <- w_cur
        res$residuals <- resid
        res$rmse <- sqrt(mean(resid^2))
      }


      # Expand weights to full veg list (zeros for excluded types)
      w_full <- rep(0, n_veg_full)
      names(w_full) <- full_veg_names
      w_full[veg_subset] <- res$w

      ids_full <- rep(NA_character_, n_veg_full)
      names(ids_full) <- full_veg_names
      ids_full[veg_subset] <- ids

      res$w <- w_full
      res$ids <- ids_full
      res$E <- E
      res$veg_subset <- veg_subset
      res$score <- score_from_solution(res)
      return(res)
    }

    # Build list of veg subsets to search over.
    # When sparse unmixing is active, also try smaller subsets (1 and 2 veg classes)
    # so the lambda penalty can select the sparsest model that fits well.
    # Barren is always included inside solve_combo_subset; subsets here refer to veg only.
    veg_only_names <- full_veg_names[tolower(full_veg_names) != "barren"]
    barren_names   <- full_veg_names[tolower(full_veg_names) == "barren"]

    sparse_search <- solver %in% c("sparse", "iwlmm_sparse") && lambda > 0
    max_veg_k <- if (sparse_search && exists("SPARSE_UNMIX_K", inherits = TRUE)) SPARSE_UNMIX_K else length(veg_only_names)

    subsets_to_try <- list()
    for (k in seq_len(max_veg_k)) {
      for (combo in combn(veg_only_names, k, simplify = FALSE)) {
        subsets_to_try[[length(subsets_to_try) + 1]] <- c(barren_names, combo)
      }
    }

    n_variants_per_veg <- sapply(top_variants, length)
    global_best_res <- NULL
    global_best_score <- Inf

    search_subset <- function(veg_subset) {
      n_variants_subset <- n_variants_per_veg[veg_subset]
      total_combos <- as.numeric(prod(n_variants_subset))

      if (total_combos <= 500) {
        combos <- expand.grid(lapply(n_variants_subset, seq_len), KEEP.OUT.ATTRS = FALSE)
        for (i in seq_len(nrow(combos))) {
          r <- solve_combo_subset(veg_subset, as.integer(combos[i, ]))
          if (!is.null(r) && r$score < global_best_score) {
            global_best_score <<- r$score
            global_best_res <<- r
          }
        }
      } else {
        TOPK_COARSE <- 3
        coarse_idx <- lapply(n_variants_subset, function(n) seq_len(min(TOPK_COARSE, n)))
        combos <- expand.grid(coarse_idx, KEEP.OUT.ATTRS = FALSE)
        best_idx <- rep(1, length(veg_subset))

        for (i in seq_len(nrow(combos))) {
          r <- solve_combo_subset(veg_subset, as.integer(combos[i, ]))
          if (!is.null(r) && r$score < global_best_score) {
            global_best_score <<- r$score
            global_best_res <<- r
            best_idx <- as.integer(combos[i, ])
          }
        }

        curr_idx <- best_idx
        improved <- TRUE
        while (improved) {
          improved <- FALSE
          for (k in seq_along(veg_subset)) {
            for (v_opt in seq_len(n_variants_subset[k])) {
              if (v_opt == curr_idx[k]) next
              t_idx <- curr_idx; t_idx[k] <- v_opt
              r <- solve_combo_subset(veg_subset, t_idx)
              if (!is.null(r) && r$score < global_best_score) {
                global_best_score <<- r$score
                global_best_res <<- r
                curr_idx <- t_idx
                improved <- TRUE
              }
            }
          }
        }
      }
    }

    for (veg_subset in subsets_to_try) {
      search_subset(veg_subset)
    }
    
    if (is.null(global_best_res)) return(NULL)

    # [MODIFIED] Return only the single best model as top_models per request
    # Previous logic sorted accumulated candidates and kept near-optimal ones.
    top_models <- list(global_best_res)

    return(list(
      w = global_best_res$w,
      rmse = global_best_res$rmse,
      score = global_best_score,
      ids = global_best_res$ids,
      residuals = global_best_res$residuals,
      E_best = global_best_res$E,
      top_models = top_models
    ))
  }

  
  build_mesma_library_weighted <- function(df_train, indices, params, allowed_veg, precomputed_clusters = NULL, generate_proto_plots = NULL) {
    # MESMA library: treats barren and all vegetation types as equal endmembers
    # Returns a library with all endmember types (barren + veg types) in one unified structure
    # If precomputed_clusters is provided (list with class -> k mapping), skip cluster optimization

    if (!is.null(precomputed_clusters)) {
      cat("\n[LIBRARY BUILD] Using precomputed cluster counts (skipping optimization)...\n")
      cat(sprintf("[LIBRARY BUILD] Precomputed clusters: %s\n", paste(names(precomputed_clusters), precomputed_clusters, sep="=", collapse=", ")))
    } else {
      cat("\n[LIBRARY BUILD] Starting Global Combinatorial Optimization for Cluster Counts (including Barren)...\n")
    }

    # Validate params structure
    if (is.null(params) || is.null(params$means) || is.null(params$sds)) {
      stop("[ERROR] build_mesma_library_weighted: params$means or params$sds is NULL!")
    }
    if (length(params$means) != length(indices)) {
      cat(sprintf("[ERROR] params$means has length %d but indices has length %d\n", length(params$means), length(indices)))
      cat(sprintf("  indices: %s\n", paste(indices, collapse = ", ")))
      cat(sprintf("  params$means names: %s\n", paste(names(params$means), collapse = ", ")))
      stop("[ERROR] Length mismatch between params$means and indices")
    }

    # -------------------------------------------------------------------------
    # PASS 1: Load, Normalize, and Store Data (with Metadata)
    # Separate storage for model training and OOB evaluation
    # -------------------------------------------------------------------------
    expected_cols <- length(indices) * TEMPORAL_BUDGET
    storage <- list()       # Store unweighted Z-score matrices (for model training)
    storage_meta <- list()  # Store metadata (location_id) for spatial bootstrapping
    storage_oob <- list()   # Store OOB data for cluster optimization evaluation
    storage_oob_meta <- list()

    # Verify OOB data is available (required, no fallback)
    df_oob_state <- if (exists("df_train_oob", envir = globalenv())) df_train_oob else NULL
    if (is.null(df_oob_state) || !is.data.frame(df_oob_state) || nrow(df_oob_state) == 0) {
      stop("[ERROR] OOB holdout data (df_train_oob) is required but not found. Ensure OOB split was performed during training data preparation.")
    }

    df_oob <- df_oob_state
    oob_locs <- if (exists("oob_location_ids", envir = globalenv())) oob_location_ids else character(0)
    has_oob_data <- TRUE  # OOB data verified above
    cat(sprintf("[LIBRARY BUILD] OOB data available: %d rows from %d locations\n",
                nrow(df_oob), length(oob_locs)))
    cat("[LIBRARY BUILD] Will use OOB data for cluster optimization evaluation\n")

    # Unified class list
    target_classes <- unique(c("barren", allowed_veg))
    valid_classes <- c()

    # Helper function to build storage from a dataframe
    build_storage_from_df <- function(df_source, storage_list, meta_list, source_name = "train") {
      l2_normalize <- isTRUE(params$l2_normalize)
      base_indices <- if (!is.null(params$base_indices)) params$base_indices else indices
      n_base_idx <- length(base_indices)
      expected_base_cols <- n_base_idx * TEMPORAL_BUDGET

      for(v in target_classes) {
        veg_data <- dplyr::filter(df_source, .data$Veg == v)
        if(nrow(veg_data) == 0) next

        veg_list <- list()
        loc_list <- character(0)

        traces <- unique(veg_data[, c("location_id", "pheno_year")])

        for(i in seq_len(nrow(traces))) {
          lid <- traces$location_id[i]
          pyr <- traces$pheno_year[i]
          sub <- veg_data[veg_data$location_id == lid & veg_data$pheno_year == pyr, ]
          # Always build from base indices
          mat <- build_pentad_matrix(sub, base_indices)
          if(!is.null(mat)) {
            veg_list[[length(veg_list) + 1]] <- as.numeric(mat)
            loc_list <- c(loc_list, as.character(lid))
          }
        }

        if(length(veg_list) > 0) {
          veg_mat_raw <- do.call(rbind, veg_list)

          if(ncol(veg_mat_raw) == expected_base_cols) {
            veg_mat <- mesma_apply_representation_mat(veg_mat_raw, n_base_idx, TEMPORAL_BUDGET, l2_normalize)
            veg_mat <- mesma_zscore_mat_by_index(veg_mat, indices, params$means, params$sds, TEMPORAL_BUDGET)

            storage_list[[v]] <- veg_mat
            meta_list[[v]] <- data.frame(location_id = loc_list, stringsAsFactors = FALSE)
          }
        }
      }
      return(list(storage = storage_list, meta = meta_list))
    }

    # Build storage from training data (excluding OOB if available)
    if (has_oob_data) {
      # Use df_train_model for building endmembers
      df_train_model_local <- if (exists("df_train_model", envir = globalenv())) df_train_model else NULL
      if (is.null(df_train_model_local) || !is.data.frame(df_train_model_local) || nrow(df_train_model_local) == 0) {
        stop("[LIBRARY BUILD] df_train_model is missing/empty; expected after OOB split")
      }
      cat(sprintf("[LIBRARY BUILD] Building endmember storage from model data: %d rows\n", nrow(df_train_model_local)))
      result <- build_storage_from_df(df_train_model_local, storage, storage_meta, "model")
      storage <- result$storage
      storage_meta <- result$meta

      # Build OOB storage for evaluation
      cat(sprintf("[LIBRARY BUILD] Building OOB storage for evaluation: %d rows\n", nrow(df_oob)))
      result_oob <- build_storage_from_df(df_oob, storage_oob, storage_oob_meta, "oob")
      storage_oob <- result_oob$storage
      storage_oob_meta <- result_oob$meta

      # Log OOB storage sizes
      for (v in names(storage_oob)) {
        cat(sprintf("[LIBRARY BUILD] OOB storage for %s: %d samples\n", v, nrow(storage_oob[[v]])))
      }
    }

    valid_classes <- names(storage)

    if (length(valid_classes) == 0) {
      warning("No valid data found for library building.")
      return(list())
    }

    # Compute average value of first pentad of first index for training data
    if (length(valid_classes) > 0) {
      all_vals <- c()
      for (v in valid_classes) {
        if (!is.null(storage[[v]]) && nrow(storage[[v]]) > 0) {
          all_vals <- c(all_vals, storage[[v]][, 1])
        }
      }
      if (length(all_vals) > 0) {
        avg_train <- mean(all_vals, na.rm = TRUE)
        cat(sprintf("[INFO] Average value of first pentad of first index (training, after normalization): %.6f\n", avg_train))
      }
    }

    # Weights for clustering (PCA-LDA weights for distance calculations)
    # Storage matrices are built to match `indices` exactly, so they always have `expected_cols` columns.
    n_storage_cols <- expected_cols
    w_vec <- if(!is.null(params$weights) && length(params$weights) == n_storage_cols) {
      pmax(params$weights, 0)
    } else {
      rep(1, n_storage_cols)
    }
    w_vec[w_vec < 1e-9] <- 1e-9 # Prevent div by zero
    
    min_cluster_size <- MIN_CLUSTER_SIZE

    # Greedy EAR-based endmember extraction (minimizes reconstruction error)
    ear_extract_all_levels_greedy <- function(data_mat, max_k, w_vec = NULL) {
      n <- nrow(data_mat)
      if (n == 0) return(NULL)

      
      # Apply weighting
      if (!is.null(w_vec) && length(w_vec) == ncol(data_mat)) {
        data_w <- sweep(data_mat, 2, w_vec, "*")
      } else {
        data_w <- data_mat
      }

      # Impute NAs in data_w for distance computation only.
      # NA features propagate through distance sums and can make which.min() return
      # integer(0), producing 0-row prototype matrices and downstream NULL crashes.
      # We replace NAs with the column median so that samples with missing pentads
      # still receive sensible (if imperfect) distances; the stored prototypes
      # are taken from data_mat (never from data_w) so imputation is invisible
      # to the returned endmembers.
      col_medians_w <- apply(data_w, 2, function(x) median(x, na.rm = TRUE))
      col_medians_w[!is.finite(col_medians_w)] <- 0
      for (ci in seq_len(ncol(data_w))) {
        na_rows <- !is.finite(data_w[, ci])
        if (any(na_rows)) data_w[na_rows, ci] <- col_medians_w[ci]
      }

      # Determine max_k effectively
      max_k <- min(max_k, n)
      
      # Precompute Distance Matrix (O(N^2))
      dist_mat <- as.matrix(dist(data_w))^2
      candidates <- 1:n

      selected <- integer(0)
      min_sq_dists <- rep(Inf, n)

      results <- list()

      for (i in 1:max_k) {
        if (i == 1) {
          # Step 1: Find Medoid (sample closest to median)
          centroid <- apply(data_w, 2, median, na.rm = TRUE)
          dists_to_centroid <- rowSums(sweep(data_w, 2, centroid, "-")^2, na.rm = TRUE)
          best_local_idx <- which.min(dists_to_centroid)
          if (length(best_local_idx) == 0) best_local_idx <- 1L
          best_idx <- candidates[best_local_idx]
          
        } else {
          remaining_local_idx <- seq_along(candidates) 
          remaining_local_idx <- remaining_local_idx[!candidates[remaining_local_idx] %in% selected]
          
          # We want to maximize total gain = sum(pmax(0, current_err - new_err))
          # new_err = dist(pt, cand)
          # gain = pmax(0, min_sq_dists - D[pt, cand])
          
          # Vectorized over all candidates
          # Calculate gains for all remaining candidates in one matrix op if poss,
          # or loop over candidates but use the precomputed dist_mat
          
          best_gain <- -1
          best_local_idx <- -1
          
          # Fast scan using precomputed distance matrix
          # Total Error with current set = sum(min_sq_dists)
          # If we add candidate C: new_error_sum = sum(min(min_sq_dists[j], dist_mat[j, C]))
          
          # We can compute this for all C in 'remaining' quickly
          current_err_sum <- sum(min_sq_dists)
          
          # This loop is effectively vectorized compared to re-calculating distances
          for(cand_loc in remaining_local_idx) {
             d_col <- dist_mat[, cand_loc]
             # Total error if we add this candidate
             new_total <- sum(pmin(min_sq_dists, d_col))
             gain <- current_err_sum - new_total
             
             if(gain > best_gain) {
               best_gain <- gain
               best_local_idx <- cand_loc
             }
          }
          best_idx <- candidates[best_local_idx]
        }
        
        selected <- c(selected, best_idx)
        
        # Update distances
        # dist_mat is [local_idx, local_idx]
        # best_local_idx corresponds to the selected candidate
        best_idx_num <- as.integer(best_idx)
        if (i == 1) {
           # Initial distances logic
           min_sq_dists <- dist_mat[, best_idx_num]
        } else {
           min_sq_dists <- pmin(min_sq_dists, dist_mat[, best_idx_num])
        }
        
        # Save snapshot
        results[[as.character(i)]] <- data_mat[selected, , drop = FALSE]
      }
      return(results)
    }

    # Endmember extraction wrapper: returns endmembers for a specific k
    ear_extract_endmembers <- function(data_mat, k, w_vec = NULL) {
      if (k < 1) return(NULL)
      res <- ear_extract_all_levels_greedy(data_mat, k, w_vec)
      return(res[[as.character(k)]])
    }
    
    optimize_library <- function(n_boot = 5) {
      cat(sprintf("\n  --- Running EAR-Based Endmember Selection with Cross-Validation (%d folds) ---\n", n_boot))
      
      # Pre-compute EAR-based endmembers on full model data
      boot_endmember_cache <- list()
      k_ranges <- list()

      cat("    Pre-computing EAR-based endmembers on full model data...\n")

      for(v in valid_classes) {
        n_total <- nrow(storage[[v]])
        max_k <- floor(n_total / min_cluster_size)
        if (v == "barren") {
           # Allow multiple barren endmembers controlled by RAW_BARREN_N_PROTOTYPES
           max_k_barren <- if (exists("RAW_BARREN_N_PROTOTYPES") && is.finite(RAW_BARREN_N_PROTOTYPES) && RAW_BARREN_N_PROTOTYPES >= 1) as.integer(RAW_BARREN_N_PROTOTYPES) else MAX_K_EAR
           k_candidates <- 1:max_k_barren
           k_candidates <- k_candidates[k_candidates <= max_k]
           if (length(k_candidates) == 0) k_candidates <- 1
        } else {
           k_candidates <- 1:MAX_K_EAR
           k_candidates <- k_candidates[k_candidates <= max_k]
           if (length(k_candidates) == 0) k_candidates <- 1
        }
        k_ranges[[v]] <- k_candidates
      }

      # Build endmembers from full storage (model data)
      boot_endmember_cache[[1]] <- list()

      for(v in valid_classes) {
         boot_endmember_cache[[1]][[v]] <- list()

         veg_mat_train <- storage[[v]]

         # Determine max K needed for this class
         max_k_needed <- max(k_ranges[[v]])

         # Compute ALL levels at once
         all_levels <- ear_extract_all_levels_greedy(veg_mat_train, max_k_needed, w_vec)

         for(k in k_ranges[[v]]) {
           k_char <- as.character(k)
           if(!is.null(all_levels[[k_char]])) {
             boot_endmember_cache[[1]][[v]][[k_char]] <- all_levels[[k_char]]
           }
         }
      }

      # --- Random Search with OOB + Training Evaluation ---
      cat(sprintf("    Optimizing library structure using OOB holdout + training data for evaluation...\n"))
      cat(sprintf("    (OOB classes available: %s)\n", paste(names(storage_oob), collapse = ", ")))

      # 3a. Prepare Test Matrices using OOB data

      # Build OOB evaluation set from storage_oob
      oob_samples <- list(); oob_lbls <- c()
      for (v in valid_classes) {
        if (v %in% names(storage_oob) && !is.null(storage_oob[[v]]) && nrow(storage_oob[[v]]) > 0) {
          oob_samples[[length(oob_samples)+1]] <- storage_oob[[v]]
          oob_lbls <- c(oob_lbls, rep(v, nrow(storage_oob[[v]])))
        }
      }

      if (length(oob_samples) == 0) {
        stop("[ERROR] No OOB samples found in storage_oob. Cannot proceed with cluster optimization.")
      }

      oob_set <- list(Y = do.call(rbind, oob_samples), labels = oob_lbls)
      cat(sprintf("    OOB evaluation set: %d samples across %d classes\n",
                  nrow(oob_set$Y), length(unique(oob_lbls))))

      # Build training evaluation set from storage (model training data, excluding OOB)
      train_samples <- list(); train_lbls <- c()
      for (v in valid_classes) {
        if (v %in% names(storage) && !is.null(storage[[v]]) && nrow(storage[[v]]) > 0) {
          train_samples[[length(train_samples)+1]] <- storage[[v]]
          train_lbls <- c(train_lbls, rep(v, nrow(storage[[v]])))
        }
      }
      if (length(train_samples) > 0) {
        train_set <- list(Y = do.call(rbind, train_samples), labels = train_lbls)
        cat(sprintf("    Training evaluation set: %d samples across %d classes\n",
                    nrow(train_set$Y), length(unique(train_lbls))))
      } else {
        train_set <- NULL
        cat("    [WARNING] No training samples available for combined evaluation; using OOB only\n")
      }

      # 3b. Fitness Function - Evaluate on OOB holdout AND training data
      # Optimizes VEGETATION ACCURACY (populus, tamarix, herbs) - excludes barren
      # Returns min of OOB and training accuracy — config must perform well on both
      score_on_dataset <- function(target_set, M, col_names, veg_classes) {
          Y_test <- target_set$Y
          labels_test <- target_set$labels
          # Drop rows where ALL features are non-finite (completely missing samples)
          finite_rows <- which(rowSums(is.finite(Y_test)) > 0)
          if (length(finite_rows) == 0) return(0)
          Y_test <- Y_test[finite_rows, , drop = FALSE]
          labels_test <- labels_test[finite_rows]
          # Replace remaining NAs with 0 so NNLS doesn't receive non-finite input
          Y_test[!is.finite(Y_test)] <- 0
          n_test <- nrow(Y_test)

          # Batch process using solve_batch_fcls for speed
          solver_mode <- if (!is.null(params$solver)) params$solver else "fcls"
          all_coefs <- solve_batch_fcls(M, Y_test,
                                        params$weights,
                                        lda_basis             = params$lda_basis,
                                        lda_component_weights = params$lda_component_weights,
                                        solver = solver_mode)

          # Compute average correctly predicted fraction for vegetation classes
          # Row-normalized after barren subtraction: veg_frac_true / sum(veg_fracs)
          veg_norm_frac_sums <- setNames(rep(0, length(veg_classes)), veg_classes)
          veg_counts <- setNames(rep(0L, length(veg_classes)), veg_classes)

          for(j in 1:n_test) {
             true_label <- labels_test[j]
             coefs <- all_coefs[j, ]
             sums <- tapply(coefs, col_names, sum)

             for (vc in valid_classes) if (!(vc %in% names(sums))) sums[[vc]] <- 0

             # Only count vegetation samples (exclude barren)
             if (true_label %in% veg_classes) {
               # Sum of vegetation fractions (excluding barren)
               veg_total <- sum(sapply(veg_classes, function(vc) sums[[vc]]), na.rm = TRUE)

               # Row-normalized fraction for true class (after barren subtraction)
               if (veg_total > 1e-10) {
                 norm_frac <- sums[[true_label]] / veg_total
               } else {
                 norm_frac <- 0
               }

               veg_norm_frac_sums[true_label] <- veg_norm_frac_sums[true_label] + norm_frac
               veg_counts[true_label] <- veg_counts[true_label] + 1L
             }
          }

          veg_diag_fracs <- ifelse(veg_counts > 0, veg_norm_frac_sums / veg_counts, NA)
          mean_veg_frac <- mean(veg_diag_fracs, na.rm = TRUE)
          if (is.finite(mean_veg_frac)) return(mean_veg_frac) else return(0)
      }

      evaluate_config <- function(combo_list) {
          if(is.null(oob_set)) return(-1.0)

          # Build library matrix using ALL k prototypes per class
          cols <- list()
          col_names <- c()
          for(v in valid_classes) {
             k_val <- as.character(combo_list[[v]])
             endmembers <- boot_endmember_cache[[1]][[v]][[k_val]]
             if(is.null(endmembers)) return(-1.0)
             if (is.null(dim(endmembers))) endmembers <- matrix(endmembers, nrow = 1)
             # Each row is a prototype; transpose so columns = endmembers
             for (ri in seq_len(nrow(endmembers))) {
               cols[[length(cols) + 1]] <- endmembers[ri, ]
               col_names <- c(col_names, v)
             }
          }

          veg_classes <- setdiff(valid_classes, "barren")

          # Complexity penalty: lambda * sum(k) across all classes
          total_k <- sum(sapply(combo_list, as.integer))
          penalty <- CLUSTER_COMPLEXITY_LAMBDA * total_k

          if (length(cols) == 0) return(-1.0)
          M <- do.call(cbind, cols)
          if (is.null(M) || !is.matrix(M) || ncol(M) == 0) return(-1.0)

          oob_score <- score_on_dataset(oob_set, M, col_names, veg_classes)

          combo_score <- if (!is.null(train_set)) {
            train_score <- score_on_dataset(train_set, M, col_names, veg_classes)
            min(oob_score, train_score) - penalty
          } else {
            oob_score - penalty
          }

          return(combo_score)
      }

      # 3c. Exhaustive Grid Search over all cluster combinations
      # Skip if precomputed clusters are provided
      if (!is.null(precomputed_clusters)) {
        cat("[CLUSTER] Using precomputed cluster counts - skipping optimization\n")
        best_combo <- precomputed_clusters
        best_mean_score <- NA  # Not computed when using precomputed
        cat(sprintf("    Using precomputed combo: %s\n", paste(names(best_combo), best_combo, sep="=", collapse=", ")))
      } else {
        # Build full grid of all combinations
        grid <- expand.grid(k_ranges, KEEP.OUT.ATTRS = FALSE)
        colnames(grid) <- valid_classes
        n_combos <- nrow(grid)

        cat(sprintf("      [Grid Search] Evaluating all %d cluster combinations...\n", n_combos))

        best_mean_score <- -1
        best_combo <- NULL

        for(i in 1:n_combos) {
           combo <- as.list(grid[i, ])

           score <- evaluate_config(combo)

           if(score > best_mean_score) {
               best_mean_score <- score
               best_combo <- combo
           }

           # Progress report
           if(i %% 50 == 0 || i == n_combos) {
              cat(sprintf("      [Combo %d/%d] Current Best Score: %.4f\n", i, n_combos, best_mean_score))
           }
        }

        # If no valid combination was found, return NULL to trigger hard failure
        if (is.null(best_combo)) {
          cat("[ERROR] No valid combinations found during grid search.\n")
          return(NULL)
        }

        cat(sprintf("    Best OOB+Train Score: %.4f using combo: %s\n", best_mean_score, paste(names(best_combo), best_combo, sep="=", collapse=", ")))

        # Store optimal cluster counts globally for reuse (avoids re-optimization after threshold tuning)
        assign("OPTIMAL_CLUSTER_COUNTS", best_combo, envir = globalenv())
        cat(sprintf("[CLUSTER] Stored optimal cluster counts: %s\n", paste(names(best_combo), best_combo, sep="=", collapse=", ")))
      }

      # NOTE: Confusion matrix is computed AFTER threshold optimization (Step 3)
      # to use the final optimized weights. Skipping here.
      cat("\n    [CLUSTER] Confusion matrix will be computed after threshold optimization (Step 3)\n")

      # Return best combination (or NULL if none found)
      return(best_combo)
    }
    
    # -------------------------------------------------------------------------
    # MAIN EXECUTION
    # -------------------------------------------------------------------------
    
    # 1. Run Optimization to get Best K-Combo
    best_combo <- optimize_library(n_boot = 5)
    
    # Fail hard if optimization failed to find a valid combo
    if (is.null(best_combo) || length(best_combo) == 0) {
      stop("[ERROR] optimize_library() failed to find any valid cluster combination. Aborting library build.")
    } else {
      cat(sprintf("[CLUSTER OPT] Best combination found: %s\n", 
                  paste(names(best_combo), "=", sapply(best_combo, as.integer), collapse = ", ")))
    }
    
    # 2. Extract Final Endmembers on FULL Dataset
    # IMPORTANT: Now that we've found optimal cluster sizes using OOB evaluation,
    # we merge OOB data back into the training storage for final endmember extraction
    cat("\n[LIBRARY BUILD] Extracting final endmembers on full dataset...\n")

    # Merge OOB data back into storage for final clustering
    storage_final <- storage
    if (length(storage_oob) > 0 && any(sapply(storage_oob, function(x) !is.null(x) && nrow(x) > 0))) {
      cat("[LIBRARY BUILD] Merging OOB data back into storage for final endmember extraction...\n")
      for (v in names(storage_oob)) {
        if (!is.null(storage_oob[[v]]) && nrow(storage_oob[[v]]) > 0) {
          if (v %in% names(storage_final) && !is.null(storage_final[[v]])) {
            n_before <- nrow(storage_final[[v]])
            storage_final[[v]] <- rbind(storage_final[[v]], storage_oob[[v]])
            n_after <- nrow(storage_final[[v]])
            cat(sprintf("  %s: merged %d OOB samples (total: %d -> %d)\n",
                        v, nrow(storage_oob[[v]]), n_before, n_after))
          } else {
            storage_final[[v]] <- storage_oob[[v]]
            cat(sprintf("  %s: added %d OOB samples (was empty)\n", v, nrow(storage_oob[[v]])))
          }
        }
      }
    } else {
      cat("[LIBRARY BUILD] No OOB data to merge, using storage as-is\n")
    }

    final_lib_cache <- list()

    for(v in valid_classes) {
       k_opt <- as.numeric(best_combo[[v]])
       veg_mat <- storage_final[[v]]

       # Check availability
       if(nrow(veg_mat) < k_opt) {
          k_opt <- max(1, nrow(veg_mat))
          cat(sprintf("  [WARNING] Reducing k for '%s' to %d (insufficient samples)\n", v, k_opt))
       }

       # EAR-based endmember extraction (brightness-invariant)
       final_endmembers <- ear_extract_endmembers(veg_mat, k_opt, w_vec)

       # Note: Barren-similar observations removed during early data pre-filtering

       if (!is.null(final_endmembers)) {
          final_lib_cache[[v]] <- final_endmembers
       }
    }
    
    # 3. Construct Final Library Structure (Standard Format)
    res_lib <- list()
    for(v in names(final_lib_cache)) {
       mat <- final_lib_cache[[v]]
       res_lib[[v]] <- list()
       for(i in 1:nrow(mat)) {
          # Construct dummy variant object
          # Re-normalization to T_medoid structure if needed?
          # The rest of the pipeline expects 'vec' (unweighted, z-scored) or 'T' (medoid).
          # Here we provide 'vec'.
          # Use storage_final for n_samples to reflect merged OOB data
          n_total <- if (v %in% names(storage_final)) nrow(storage_final[[v]]) else nrow(storage[[v]])
          res_lib[[v]][[length(res_lib[[v]]) + 1]] <- list(
             vec = mat[i, ],
             id = paste0(v, "_opt_", i),
             n_samples = floor(n_total / nrow(mat)) # Approximate n_samples
          )
       }
    }

    # -------------------------------------------------------------------------
    # (Barren-similar observations removed during early data pre-filtering)

    # Optional: Generate prototype plots (one plot per index/band) showing endmember centers across pentads
    # NOTE: Plots use the exact feature vectors stored in the model/library (no plot-only transforms).
    #
    # If `feature_weights` is supplied (a numeric vector of length
    # length(indices)*TEMPORAL_BUDGET), the background of each pentad is shaded
    # proportionally to the corresponding weight (higher weight = darker),
    # providing a visual cue for temporally-important bins.  If weights are not
    # provided, a binary 'significant' mask may still be used for light shading.
    plot_vegetation_prototypes <- function(lib, indices = NULL, out_dir = if (exists("OUT_DIR")) OUT_DIR else ".", prefix = "veg_prototypes", save_png = TRUE, dpi = 150, feature_weights = NULL, lda_weights = NULL, show_medians = TRUE, variant_alpha = NULL) {
      if (is.null(lib) || length(lib) == 0) return(NULL)
      if (is.null(indices)) {
        if (exists("MESMA_PARAMS") && !is.null(MESMA_PARAMS) && !is.null(MESMA_PARAMS$indices)) indices <- MESMA_PARAMS$indices else indices <- OPTIMAL_INDICES
      }
      if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required for prototype plotting")
      if (!dir.exists(out_dir) && save_png) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      # Precompute pentad-to-month mapping for x-axis (pheno year: Mar=1)
      # Month start days in phenological DOY (March 1 = day 1)
      pheno_month_starts <- c(1, 32, 62, 93, 123, 154, 184, 215, 245, 276, 306, 337)
      pheno_month_labels <- c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan","Feb")
      # Convert month-start pheno DOY to pentad position: pentad = ceil(doy / agg_days)
      month_breaks_pentad <- ceiling(pheno_month_starts / TEMPORAL_AGGREGATION_DAYS)
      # Exclude the last (partial/reserved) pentad from plots — it is always
      # sparse and creates a visual artefact at the right edge of the plot.
      n_plot_bins <- TEMPORAL_BUDGET - 1L

      # Keep only unique breaks within the plotted temporal budget
      keep <- month_breaks_pentad >= 1 & month_breaks_pentad <= n_plot_bins & !duplicated(month_breaks_pentad)
      month_breaks_pentad <- month_breaks_pentad[keep]
      month_labels <- pheno_month_labels[keep]

      # Build per-index weight mask from feature_weights (0 = pruned / non-significant)
      # weights may be arbitrary positive numbers; we'll normalize later per-index when plotting.
      # `orig_indices` preserves the full ordered index list so that slices into `vec`
      # (which was built from the original full set) remain correctly aligned even after
      # excluded indices are dropped from `indices`.
      orig_indices <- indices
      weight_mask <- NULL      # pruning mask (0 = pruned): used for fill colour
      lda_weight_mask <- NULL  # raw LDA weights: used for bar heights (always positive)
      if (!is.null(feature_weights) && length(feature_weights) == length(indices) * TEMPORAL_BUDGET) {
        weight_mask <- list()
        for (k in seq_along(indices)) {
          w_start <- (k - 1) * TEMPORAL_BUDGET + 1
          weight_mask[[indices[k]]] <- feature_weights[w_start:(w_start + n_plot_bins - 1L)]
        }
        # Drop fully-excluded indices (all pentad weights zero) — do NOT plot these
        excluded <- names(weight_mask)[vapply(weight_mask, function(m) all(m == 0), logical(1))]
        if (length(excluded) > 0) {
          cat(sprintf("[PROTO PLOT] Skipping %d fully-excluded indices: %s\n", length(excluded), paste(excluded, collapse = ", ")))
          indices <- setdiff(indices, excluded)
          weight_mask <- weight_mask[indices]
        }
      }
      if (!is.null(lda_weights) && length(lda_weights) == length(orig_indices) * TEMPORAL_BUDGET) {
        lda_weight_mask <- list()
        for (k in seq_along(orig_indices)) {
          w_start <- (k - 1) * TEMPORAL_BUDGET + 1
          lda_weight_mask[[orig_indices[k]]] <- lda_weights[w_start:(w_start + n_plot_bins - 1L)]
        }
      }

      # Build long dataframe for plotting
      rows <- list()
      for (v in names(lib)) {
        variants <- lib[[v]]
        for (var in variants) {
          vec <- as.numeric(var$vec)
          vid <- if (!is.null(var$id)) var$id else if (!is.null(var$variant_id)) var$variant_id else paste0(v, "_unknown")
          # vec was constructed from orig_indices — use its length as the ground truth
          expected_len <- length(orig_indices) * TEMPORAL_BUDGET
          if (length(vec) < expected_len) {
            warning(sprintf("Skipping variant %s for veg %s: length(vec)=%d != expected=%d", vid, v, length(vec), expected_len))
            next
          }
          # Iterate over orig_indices so slice positions are always correct;
          # skip any index that was removed (fully excluded or not in the active set).
          for (k in seq_along(orig_indices)) {
            idx_name <- orig_indices[k]
            if (!idx_name %in% indices) next   # excluded index — don't plot
            start <- (k - 1) * TEMPORAL_BUDGET + 1
            end   <- (k - 1) * TEMPORAL_BUDGET + n_plot_bins   # drop last pentad
            vals <- vec[start:end]
            sig <- if (!is.null(weight_mask) && idx_name %in% names(weight_mask)) weight_mask[[idx_name]][1:n_plot_bins] > 0 else rep(TRUE, n_plot_bins)
            wt <- if (!is.null(weight_mask) && idx_name %in% names(weight_mask)) weight_mask[[idx_name]][1:n_plot_bins] else rep(NA_real_, n_plot_bins)

            df_tmp <- data.frame(pentad = seq_len(n_plot_bins), value = vals, Veg = v, variant_id = vid, index = idx_name, significant = sig, weight = wt, stringsAsFactors = FALSE)
            rows[[length(rows) + 1]] <- df_tmp
          }
        }
      }

      if (length(rows) == 0) return(NULL)
      proto_df <- do.call(rbind, rows)

      plots <- list()
      for (idx in unique(proto_df$index)) {
        df_idx <- proto_df[proto_df$index == idx, , drop = FALSE]

        # Prepare a veg -> color mapping: prefer user-provided VEG_CALIBRATION_COLORS (case-insensitive),
        # otherwise fall back to RColorBrewer Set1 palette.
        veg_levels <- unique(df_idx$Veg)
        veg_palette <- NULL
        if (exists("VEG_CALIBRATION_COLORS", envir = globalenv())) {
          supplied <- get("VEG_CALIBRATION_COLORS", envir = globalenv())
          # Match by case-insensitive names
          matched <- sapply(veg_levels, function(v) {
            nm <- names(supplied)
            im <- which(tolower(nm) == tolower(v))
            if (length(im) > 0) supplied[im[1]] else NA_character_
          }, USE.NAMES = FALSE)
          if (!all(is.na(matched))) {
            veg_palette <- setNames(matched, veg_levels)
          }
        }

        # Fallback: a simple Set1 palette sized to available veg levels
        if (is.null(veg_palette) || any(is.na(veg_palette))) {
          if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
            stop("RColorBrewer required for default veg palette")
          }
          nveg <- max(3, length(veg_levels))
          brewer_cols <- RColorBrewer::brewer.pal(n = nveg, name = "Set1")
          brewer_cols <- brewer_cols[seq_len(length(veg_levels))]
          names(brewer_cols) <- veg_levels
          if (is.null(veg_palette)) veg_palette <- brewer_cols else {
            na_idx <- which(is.na(veg_palette))
            if (length(na_idx) > 0) veg_palette[na_idx] <- brewer_cols[na_idx]
          }
        }

        # Compute y-axis limits to include all data (including negative values)
        y_min <- min(df_idx$value, na.rm = TRUE)
        y_max <- max(df_idx$value, na.rm = TRUE)
        # Add 5% padding on each side
        y_range <- y_max - y_min
        y_pad <- if (y_range > 0) y_range * 0.05 else 0.1

        # Determine if we have per-bin weight information
        has_weight <- "weight" %in% names(df_idx) && any(is.finite(df_idx$weight))

        p <- ggplot2::ggplot(df_idx, ggplot2::aes(x = pentad, y = value, group = variant_id))
        p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.4)

        # Shade pruned (zero-weight) pentads with a light grey rectangle so they
        # remain distinguishable but don't dominate the view.  The weight panel
        # below carries the continuous weight information.
        if (has_weight) {
          wts <- if (!is.null(weight_mask) && idx %in% names(weight_mask)) {
                    weight_mask[[idx]]
                 } else {
                    stats::aggregate(weight ~ pentad, df_idx, FUN = function(x) x[1])$weight
                 }
          if (length(wts) == TEMPORAL_BUDGET) {
            pruned <- wts == 0
            if (any(pruned)) {
              pr_df <- data.frame(xmin = which(pruned) - 0.5,
                                  xmax = which(pruned) + 0.5,
                                  ymin = y_min - y_pad, ymax = y_max + y_pad)
              p <- p + ggplot2::geom_rect(data = pr_df, inherit.aes = FALSE,
                                          ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                                          fill = "grey85", alpha = 0.65)
            }
          }
        } else {
          # fallback: grey shading for non-significant pentads
          has_nonsig <- "significant" %in% names(df_idx) && any(!df_idx$significant)
          if (has_nonsig) {
            nonsig_pentads <- unique(df_idx$pentad[!df_idx$significant])
            shade_df <- data.frame(xmin = nonsig_pentads - 0.5, xmax = nonsig_pentads + 0.5,
                                   ymin = y_min - y_pad, ymax = y_max + y_pad)
            p <- p + ggplot2::geom_rect(data = shade_df, inherit.aes = FALSE,
                                        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                                        fill = "grey90", alpha = 0.5)
          }
        }

        if (is.null(variant_alpha)) {
          variant_alpha <- if (isTRUE(show_medians)) 0.15 else 0.6
        }
        variant_size <- if (isTRUE(show_medians)) 0.45 else 0.65

        # Always draw all variant curves; optionally overlay a per-class median curve.
        p <- p + ggplot2::geom_line(ggplot2::aes(color = Veg), alpha = variant_alpha, size = variant_size)

        if (isTRUE(show_medians)) {
          # Compute median explicitly so we can add markers and keep full control of styling.
          med_df <- dplyr::summarise(
            dplyr::group_by(df_idx, .data$Veg, .data$pentad),
            value = stats::median(.data$value, na.rm = TRUE),
            .groups = "drop"
          )
          med_pts <- med_df[med_df$pentad %in% month_breaks_pentad, , drop = FALSE]

          p <- p +
            ggplot2::geom_line(data = med_df, inherit.aes = FALSE,
                               ggplot2::aes(x = pentad, y = value, color = Veg, group = Veg),
                               alpha = 0.95, size = 1.25) +
            ggplot2::geom_point(data = med_pts, inherit.aes = FALSE,
                                ggplot2::aes(x = pentad, y = value, color = Veg),
                                alpha = 0.95, size = 1.3)
        }

        y_lab <- sprintf("%s (L2-normalized feature space)", idx)

        subtitle_txt <- if (isTRUE(show_medians)) {
          "Thin lines: variants; thick lines: class median"
        } else {
          "Variants only (no median overlay)"
        }

           p <- p +
             ggplot2::labs(title = sprintf("Prototypes: %s", idx),
                           subtitle = subtitle_txt,
                           x = NULL,
                           y = y_lab) +
             ggplot2::scale_x_continuous(breaks = month_breaks_pentad, labels = month_labels,
                                         limits = c(0.5, n_plot_bins + 0.5)) +
             ggplot2::theme_minimal() +
             ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                            plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9),
                            axis.text.x = ggplot2::element_blank(),
                            axis.ticks.x = ggplot2::element_blank()) +
             ggplot2::scale_color_manual(values = veg_palette) +
             ggplot2::coord_cartesian(ylim = c(y_min - y_pad, y_max + y_pad))

        # ── Weight strip panel ─────────────────────────────────────────────────
        # A bar chart showing the per-pentad permutation weight below the main
        # prototype panel.  Pruned bins (weight = 0) are coloured grey so they
        # are visually de-emphasized; surviving bins use blue.
        w_panel <- NULL
        if (has_weight && !is.null(weight_mask) && idx %in% names(weight_mask)) {
          wts_plot <- weight_mask[[idx]]   # pruning mask (0/1 or 0/positive)
          # Use raw LDA weights for bar heights if available, else fall back to mask
          raw_wts <- if (!is.null(lda_weight_mask) && idx %in% names(lda_weight_mask)) {
            lda_weight_mask[[idx]]
          } else {
            wts_plot
          }
          if (length(wts_plot) == n_plot_bins) {
            wmax <- max(raw_wts, na.rm = TRUE)
            wnorm_plot <- if (wmax > 0) raw_wts / wmax else rep(0, n_plot_bins)
            w_df <- data.frame(
              pentad = seq_len(n_plot_bins),
              weight = wnorm_plot,
              pruned = wts_plot == 0
            )
            w_panel <- ggplot2::ggplot(w_df, ggplot2::aes(x = pentad)) +
              ggplot2::geom_col(ggplot2::aes(y = weight, fill = pruned), width = 0.85, show.legend = FALSE) +
              ggplot2::scale_fill_manual(values = c("FALSE" = "#2166ac", "TRUE" = "grey70")) +
              ggplot2::scale_x_continuous(breaks = month_breaks_pentad, labels = month_labels,
                                          limits = c(0.5, n_plot_bins + 0.5)) +
              ggplot2::scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"),
                                          limits = c(0, 1.05)) +
              ggplot2::labs(x = "Month", y = "Rel. weight") +
              ggplot2::theme_minimal() +
              ggplot2::theme(
                axis.title.y     = ggplot2::element_text(size = 7),
                axis.text.y      = ggplot2::element_text(size = 6),
                panel.grid.minor = ggplot2::element_blank()
              )
          }
        }

        # Compose with patchwork if available; otherwise attach x-axis to main plot
        final_plot <- tryCatch({
          if (!is.null(w_panel) && requireNamespace("patchwork", quietly = TRUE)) {
            patchwork::wrap_plots(p, w_panel, ncol = 1, heights = c(3, 1))
          } else {
            # Fallback: restore x-axis labels directly on the prototype panel
            p + ggplot2::theme(axis.text.x = ggplot2::element_text(),
                               axis.ticks.x = ggplot2::element_line()) +
              ggplot2::labs(x = "Month")
          }
        }, error = function(e) {
          p + ggplot2::theme(axis.text.x = ggplot2::element_text(),
                             axis.ticks.x = ggplot2::element_line()) +
            ggplot2::labs(x = "Month")
        })

        plots[[idx]] <- final_plot
        if (save_png) {
          fn <- file.path(out_dir, sprintf("%s_%s.png", prefix, idx))
          ggplot2::ggsave(filename = fn, plot = final_plot, width = 8,
                          height = if (!is.null(w_panel)) 5 else 4, dpi = dpi)
        }
      }
      invisible(plots)
    }

    if (is.null(generate_proto_plots)) {
      generate_proto_plots <- exists("GENERATE_PROTO_PLOTS") && isTRUE(GENERATE_PROTO_PLOTS)
    }

    if (isTRUE(generate_proto_plots)) {
      tryCatch({
        # Use params$indices (passed to this function) instead of MESMA_PARAMS which may not exist yet
        plot_indices <- if (!is.null(params) && !is.null(params$indices)) params$indices else indices
        plot_weights     <- if (!is.null(params) && !is.null(params$weights)) params$weights else NULL
        plot_lda_weights <- plot_weights
        out_base <- file.path(if (exists("OUT_DIR")) OUT_DIR else ".", "prototype_plots")
        plot_vegetation_prototypes(res_lib, indices = plot_indices, out_dir = out_base, feature_weights = plot_weights, lda_weights = plot_lda_weights, show_medians = TRUE)
        cat(sprintf("[NOTICE] Generated prototype plots (with medians) to %s\n", out_base))

        if (exists("GENERATE_PROTO_PLOTS_VARIANTS_ONLY") && isTRUE(GENERATE_PROTO_PLOTS_VARIANTS_ONLY)) {
          out_variants <- file.path(if (exists("OUT_DIR")) OUT_DIR else ".", "prototype_plots_variants_only")
          plot_vegetation_prototypes(res_lib, indices = plot_indices, out_dir = out_variants, feature_weights = plot_weights, lda_weights = plot_lda_weights, show_medians = FALSE)
          cat(sprintf("[NOTICE] Generated prototype plots (variants-only) to %s\n", out_variants))
        }
      }, error = function(e) {
        cat(sprintf("[WARN] Failed to generate prototype plots: %s\n", e$message))
      })
    }

    return(res_lib)
  }


  precompute_optimized_library_weighted <- function(mesma_lib, grid_type = "full", feature_weights = NULL) {
    opt_lib <- list()

    # Determine pruning indices if requested
    keep_idx_global_w <- NULL
    if (!is.null(feature_weights) && exists("PRUNE_ZERO_WEIGHT_FEATURES") && isTRUE(PRUNE_ZERO_WEIGHT_FEATURES)) {
      zero_mask_glob <- feature_weights == 0
      n_zero_glob <- sum(zero_mask_glob, na.rm = TRUE)
      frac_zero_glob <- n_zero_glob / length(feature_weights)
      if (n_zero_glob > 0) {
        cat(sprintf("[INFO] MESMA pruning requested: %d/%d features zeroed (%.1f%%)\n", n_zero_glob, length(feature_weights), 100*frac_zero_glob))

        # Helpful diagnostics: summarize which indices were effectively ignored (all pentads zero-weighted)
        if (exists("MESMA_PARAMS") && !is.null(MESMA_PARAMS) &&
            !is.null(MESMA_PARAMS$indices) &&
            length(feature_weights) == (TEMPORAL_BUDGET * length(MESMA_PARAMS$indices))) {
          idxs <- MESMA_PARAMS$indices
          bins <- TEMPORAL_BUDGET
          zero_by_idx <- vapply(seq_along(idxs), function(k) {
            rng <- ((k - 1) * bins + 1):(k * bins)
            sum(feature_weights[rng] == 0, na.rm = TRUE)
          }, numeric(1))
          names(zero_by_idx) <- idxs

          all_zero <- names(zero_by_idx)[zero_by_idx >= bins]
          if (length(all_zero) > 0) {
            cat(sprintf("[INFO] Indices with ALL %d pentads zero-weighted: %s\n", bins, paste(all_zero, collapse = ", ")))
          }
        }

        if ((length(feature_weights) - n_zero_glob) > 0) {
          keep_idx_global_w <- which(!zero_mask_glob)
          cat(sprintf("[INFO] Will prune %d features for MESMA, keeping %d features\n", n_zero_glob, length(keep_idx_global_w)))
        } else {
          cat("[WARN] Skipping MESMA pruning: all features are zero-weighted.\n")
        }
      }
    }

    for(v in names(mesma_lib)) {
        if(is.null(mesma_lib[[v]])) next
        variants_raw <- mesma_lib[[v]]
        if (!is.list(variants_raw) || (!is.null(variants_raw$vec) && !is.list(variants_raw$vec))) {
          variants <- list(variants_raw)
        } else {
          variants <- variants_raw
        }

        vecs <- lapply(variants, function(x) as.numeric(if (is.list(x) && !is.null(x$vec)) x$vec else x))
        ids <- sapply(variants, function(x) {
          if (is.list(x) && !is.null(x$id)) return(x$id)
          if (is.list(x) && !is.null(x$variant_id)) return(x$variant_id)
          return(paste0(v, "_1"))
        })

        M <- do.call(rbind, vecs)

        pruned_info <- NULL
        if (!is.null(keep_idx_global_w) && length(keep_idx_global_w) > 0) {
          if (max(keep_idx_global_w) <= ncol(M)) {
            kept_names <- NULL
            if (!is.null(colnames(M))) kept_names <- colnames(M)[keep_idx_global_w]
            M <- M[, keep_idx_global_w, drop = FALSE]
            pruned_info <- list(kept_idx = keep_idx_global_w, kept_names = kept_names, n_kept = ncol(M))
            cat(sprintf("[INFO] Pruned MESMA features for veg '%s': kept %d columns\n", v, ncol(M)))
          } else {
            cat(sprintf("[WARN] MESMA prune requested but index mismatch for veg '%s'; skipping prune\n", v))
          }
        }

        M_norm <- t(apply(M, 1, function(r) {
          nrm <- sqrt(sum(r^2))
          if (!is.finite(nrm) || nrm == 0) return(rep(0, length(r))) else return(r / nrm)
        }))

        # Note: Barren-similar observations removed during early data pre-filtering
        # This ensures accuracy scores reflect performance of filtered endmembers

        opt_lib[[v]] <- list(
          M = M,
          M_norm = M_norm,
          ids = ids,
          pruned_info = pruned_info
        )
    }
    
    opt_lib
  }





  precompute_compressed_templates <- function(mesma_lib, grid_name = "full") {
    template_db <- list()
    for (veg in names(mesma_lib)) {
      template_db[[veg]] <- list()
      variants <- mesma_lib[[veg]]
      for (variant in variants) {
        vid <- if (!is.null(variant$variant_id)) variant$variant_id else variant$id
        if (is.null(vid)) next
        template_db[[veg]][[vid]] <- list()
        
        raw_src <- if (!is.null(variant$raw_mat)) variant$raw_mat else if (!is.null(variant$vec)) variant$vec else NULL
        if (is.null(raw_src)) next
        
        compressed_vec <- as.numeric(raw_src)
        compressed_vec[!is.finite(compressed_vec)] <- NA_real_
        template_db[[veg]][[vid]][[grid_name]] <- compressed_vec
      }
    }
    template_db
  }

print_weights_summary <- function(stage_name, params) {
  if (is.null(params) || is.null(params$weights) || is.null(params$indices)) {
    cat(sprintf("[WEIGHTS %s] No weights available\n", stage_name))
    return(invisible(NULL))
  }

  n_indices <- length(params$indices)
  wlen <- length(params$weights)
  expected <- n_indices * TEMPORAL_BUDGET

  if (wlen == expected) {
    w_mat <- matrix(params$weights, nrow = TEMPORAL_BUDGET, ncol = n_indices)
    per_index_mean <- colMeans(w_mat, na.rm = TRUE)
    per_index_max <- apply(w_mat, 2, max, na.rm = TRUE)
    ord <- order(per_index_mean, decreasing = TRUE)

    cat(sprintf("[WEIGHTS %s] Per-index mean weights (top %d):\n", stage_name, min(8, n_indices)))
    topn <- min(8, n_indices)
    for (i in seq_len(topn)) {
      ii <- ord[i]
      cat(sprintf("  %s: mean=%.4f, max=%.4f\n", params$indices[ii], per_index_mean[ii], per_index_max[ii]))
    }

    cat(sprintf("[WEIGHTS %s] Full per-index means: %s\n", stage_name, paste(sprintf("%s=%.4f", params$indices, per_index_mean), collapse=", ")))
  } else {
    cat(sprintf("[WEIGHTS %s] weights length (%d) != expected (%d) -> printing sample values\n", stage_name, wlen, expected))
    cat(sprintf("[WEIGHTS %s] sample weights (first 20): %s\n", stage_name, paste(sprintf("%.4f", head(params$weights, 20)), collapse=", ")))
  }
  invisible(NULL)
}

  cat("Building MESMA library with Z-score/PCA/LDA weighting...\n")


  # ==========================================================================
  # OOB-OPTIMIZED PCA-LDA THRESHOLD SELECTION
  # Step 1: Train initial PCA-LDA on df_train_model (excluding OOB tuning set)
  # Step 2: Find optimal threshold for zeroing weights using OOB data
  # Step 3: Remove variables with zero weights and re-run PCA-LDA
  # ==========================================================================

  # Create multi-class classification data: Barren vs Veg A vs Veg B ...
  # This maximizes separation between ALL classes in the feature space.
  # Use df_train_model for initial fitting if OOB split was performed
  train_for_pipeline <- NULL
  if (exists("df_train_model") && !is.null(df_train_model) && nrow(df_train_model) > 0) {
    train_for_pipeline <- df_train_model
    cat("[PCA/LDA] Using df_train_model (OOB holdout excluded) for initial PCA-LDA fitting\n")
  } else if (exists("df_train") && !is.null(df_train) && nrow(df_train) > 0) {
    train_for_pipeline <- df_train
    cat("[PCA/LDA] Using df_train for PCA-LDA fitting (no OOB split)\n")
  } else if (exists("df") && !is.null(df) && nrow(df) > 0) {
    train_for_pipeline <- df
  } else {
    stop("No training data available for PCA/LDA training pipeline")
  }

  # Exclude validation locations
  if (exists("df_validation") && !is.null(df_validation) && nrow(df_validation) > 0) {
    val_locs <- unique(df_validation$location_id)
    before_n <- nrow(train_for_pipeline)
    train_for_pipeline <- train_for_pipeline[!(train_for_pipeline$location_id %in% val_locs), , drop = FALSE]
    after_n <- nrow(train_for_pipeline)
    cat(sprintf("[PCA/LDA] Excluded %d rows belonging to %d validation locations from PCA/LDA training (rows: %d -> %d)\n", before_n - after_n, length(val_locs), before_n, after_n))
  } else {
    cat(sprintf("[PCA/LDA] Training rows used for PCA/LDA: %d\n", nrow(train_for_pipeline)))
  }

  # One-time training feature-space check (ensures all columns in `avail` exist)
  train_for_pipeline <- mesma_prepare_feature_columns(train_for_pipeline, required_cols = avail, optional_cols = NULL, context = "training")

  multi_class_data <- dplyr::mutate(train_for_pipeline, target_class = tolower(as.character(Veg)))
  multi_class_data <- dplyr::filter(multi_class_data, !is.na(target_class) & target_class != "")
  multi_class_data <- dplyr::select(multi_class_data, dplyr::all_of(c("location_id", "pheno_year", "date", "doy", "Veg", "target_class", avail)))

  cat("Training feature pipeline (Multi-Class: Barren vs Veg A vs Veg B...)\n")
  cat(sprintf("[NOTICE] Using feature set for training: %s\n", paste(avail, collapse=", ")))

  class_counts <- table(multi_class_data$target_class)
  cat("[MESMA] Training data class distribution:\n")
  print(class_counts)

  # Step 1: Initial PCA-LDA training on model data (excluding OOB)
  cat("\n=== STEP 1: Initial PCA-LDA Training (on model data, excluding OOB) ===\n")
  MESMA_PARAMS_INITIAL <- train_feature_pipeline(multi_class_data, "target_class", avail)
  # set solver mode from config default if not already defined
  if (!is.null(MESMA_PARAMS_INITIAL) && is.null(MESMA_PARAMS_INITIAL$solver)) {
    MESMA_PARAMS_INITIAL$solver <- if (exists("DEFAULT_UNMIX_SOLVER")) DEFAULT_UNMIX_SOLVER else "fcls"
  }
  # override with global toggles if requested
  if (isTRUE(exists("USE_IWLMM") && USE_IWLMM) && isTRUE(exists("USE_SPARSE_UNMIXING") && USE_SPARSE_UNMIXING)) {
    MESMA_PARAMS_INITIAL$solver <- "iwlmm_sparse"
  } else if (isTRUE(exists("USE_IWLMM") && USE_IWLMM)) {
    MESMA_PARAMS_INITIAL$solver <- "iwlmm"
  } else if (isTRUE(exists("USE_SPARSE_UNMIXING") && USE_SPARSE_UNMIXING)) {
    MESMA_PARAMS_INITIAL$solver <- "sparse"
  }
  cat(sprintf("[SOLVER] Active solver: %s (USE_IWLMM=%s, USE_SPARSE_UNMIXING=%s)\n",
              MESMA_PARAMS_INITIAL$solver,
              isTRUE(exists("USE_IWLMM") && USE_IWLMM),
              isTRUE(exists("USE_SPARSE_UNMIXING") && USE_SPARSE_UNMIXING)))

  if (!is.null(MESMA_PARAMS_INITIAL) && !is.null(MESMA_PARAMS_INITIAL$weights)) {
    MESMA_PARAMS_INITIAL$weights[is.na(MESMA_PARAMS_INITIAL$weights)] <- 0
    MESMA_PARAMS_INITIAL$weights[!is.finite(MESMA_PARAMS_INITIAL$weights)] <- 0
    print_weights_summary("INITIAL_PCA_LDA", MESMA_PARAMS_INITIAL)
    
    # Save PCA-LDA training parameters to training_norm_params.rds for use in plot_pca_lda_scatter.R
    # This ensures the plot uses the exact same PCA-LDA space as the fitting
    if (file.exists("training_norm_params.rds")) {
      TRAINING_NORM_PARAMS <- readRDS("training_norm_params.rds")
    } else {
      TRAINING_NORM_PARAMS <- list(INDEX_SCALES = NULL)
    }
    
    # Augment with PCA-LDA training parameters
    TRAINING_NORM_PARAMS$l2_normalize <- MESMA_PARAMS_INITIAL$l2_normalize
    TRAINING_NORM_PARAMS$zscore_applied <- MESMA_PARAMS_INITIAL$zscore_applied
    TRAINING_NORM_PARAMS$means <- MESMA_PARAMS_INITIAL$means
    TRAINING_NORM_PARAMS$sds <- MESMA_PARAMS_INITIAL$sds
    TRAINING_NORM_PARAMS$lda_basis <- MESMA_PARAMS_INITIAL$lda_basis
    TRAINING_NORM_PARAMS$lda_component_weights <- MESMA_PARAMS_INITIAL$lda_component_weights
    # Store PCA parameters for projection in plot scripts
    if (!is.null(MESMA_PARAMS_INITIAL$pca_rotation) && !is.null(MESMA_PARAMS_INITIAL$pca_sdev)) {
      TRAINING_NORM_PARAMS$pca_loadings <- list(
        rotation = MESMA_PARAMS_INITIAL$pca_rotation,
        sdev = MESMA_PARAMS_INITIAL$pca_sdev
      )
      cat("[INFO] Saved PCA loadings for plot_pca_lda_scatter.R\n")
    }
    # Store the trained indices to ensure plot_pca_lda_scatter.R uses the same features
    TRAINING_NORM_PARAMS$indices <- MESMA_PARAMS_INITIAL$indices
    
    saveRDS(TRAINING_NORM_PARAMS, "training_norm_params.rds")
    cat("[INFO] Updated training_norm_params.rds with PCA-LDA basis for plot_pca_lda_scatter.R\n")
  }

  # Step 2: Optimal cluster sizing FIRST (before threshold optimization)
  # This ensures cluster optimization uses the full feature set
  cat("\n=== STEP 2: Optimal Cluster Sizing (before threshold optimization) ===\n")
  cat("Building initial MESMA library with full feature weights for cluster optimization...\n")

  # Use indices from params (includes L2norm if enabled) instead of original avail
  indices_for_library <- if (!is.null(MESMA_PARAMS_INITIAL$indices)) MESMA_PARAMS_INITIAL$indices else avail

  # Build library with initial (unthresholded) weights to determine optimal cluster counts
  mesma_lib_initial <- build_mesma_library_weighted(df_train, indices_for_library, MESMA_PARAMS_INITIAL, ALLOWED_VEG, generate_proto_plots = FALSE)

  # Store the optimal cluster counts discovered during library building
  if (exists("OPTIMAL_CLUSTER_COUNTS", envir = globalenv())) {
    cat("[CLUSTER] Optimal cluster counts determined from full-feature optimization\n")
  }

  pruned_indices <- avail
  MESMA_PARAMS   <- MESMA_PARAMS_INITIAL

  if (!is.null(MESMA_PARAMS) && !is.null(MESMA_PARAMS$weights)) {
    print_weights_summary("MESMA_FINAL", MESMA_PARAMS)

    # ── Per-index summary table (post-Phase-2 LDA rebuild) ──────────────────
    tryCatch({
      idx_final  <- MESMA_PARAMS$indices
      w_final    <- MESMA_PARAMS$weights

      if (!is.null(idx_final) && length(w_final) == length(idx_final) * TEMPORAL_BUDGET) {
        tbl <- data.frame(
          Index        = idx_final,
          Mean_LDA_wt  = NA_real_,
          stringsAsFactors = FALSE
        )
        for (k in seq_along(idx_final)) {
          rng   <- (k - 1L) * TEMPORAL_BUDGET + seq_len(TEMPORAL_BUDGET)
          wvals <- w_final[rng]
          tbl$Mean_LDA_wt[k] <- mean(abs(wvals), na.rm = TRUE)
        }
        tbl <- tbl[order(tbl$Mean_LDA_wt, decreasing = TRUE), ]
        cat("\n=== Post-Phase-2 LDA: per-index summary ===\n")
        cat(sprintf("%-12s  %12s\n", "Index", "Mean_LDA_wt"))
        cat(strrep("-", 26), "\n")
        for (i in seq_len(nrow(tbl))) {
          cat(sprintf("%-12s  %11.4f\n", tbl$Index[i], tbl$Mean_LDA_wt[i]))
        }
        cat(strrep("-", 26), "\n\n")
      }
    }, error = function(e) {
      cat(sprintf("[WARN] Could not print per-index summary table: %s\n", e$message))
    })
  }

  # Store pruning results
  assign("PRUNED_INDICES", pruned_indices, envir = globalenv())

  # Step 3c: Re-fit PCA-LDA on full dataset (training + OOB) using the pruned feature set.
  # The permutation test determined which features to keep using training-only data.
  # Now we refit on all available data so the final library and scoring use consistent params.
  cat("\n=== STEP 3c: Re-fitting PCA-LDA on full dataset (training + OOB) ===\n")
  final_indices <- if (!is.null(MESMA_PARAMS$indices)) MESMA_PARAMS$indices else pruned_indices
  if (exists("df_train") && !is.null(df_train) && nrow(df_train) > 0) {
    multi_class_data_full <- dplyr::mutate(df_train, target_class = tolower(as.character(Veg)))
    multi_class_data_full <- dplyr::filter(multi_class_data_full, !is.na(target_class) & target_class != "")
    multi_class_data_full <- dplyr::select(multi_class_data_full,
                                            dplyr::all_of(c("location_id", "pheno_year", "date", "doy", "Veg", "target_class", final_indices)))
    MESMA_PARAMS_FULL <- tryCatch(
      train_feature_pipeline(multi_class_data_full, "target_class", final_indices),
      error = function(e) {
        cat(sprintf("[STEP 3c] Re-fit on full dataset failed (%s); keeping training-only params.\n", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(MESMA_PARAMS_FULL)) {
      # Propagate solver/behaviour fields from the pruned params
      for (.field in c("solver", "l2_normalize", "interpolate_inference")) {
        if (!is.null(MESMA_PARAMS[[.field]]) && is.null(MESMA_PARAMS_FULL[[.field]]))
          MESMA_PARAMS_FULL[[.field]] <- MESMA_PARAMS[[.field]]
      }
      MESMA_PARAMS <- MESMA_PARAMS_FULL
      cat("[STEP 3c] PCA-LDA successfully re-fitted on full dataset.\n")
    }
  } else {
    cat("[STEP 3c] df_train not available; keeping training-only PCA-LDA params.\n")
  }

  # NOTE: STEP 4 (confusion matrix / validation matrix) is computed later,
  # after the OOB cluster optimization that happens during final library build.

  # MESMA UNMIXING: All endmembers (barren + veg types) treated equally
  cat("[MODE] MESMA unmixing ENABLED (barren and vegetation types treated as equals)\n")

  # Log loss function choice (now always RMSE)
  cat("[MODE] Using standard RMSE loss for FCLS\n")



  cat("Building final MESMA library (barren + all vegetation types)...\n")

  # Use precomputed cluster counts from Step 2 (if available) to skip re-optimization
  precomputed_clusters <- NULL
  if (exists("OPTIMAL_CLUSTER_COUNTS", envir = globalenv())) {
    precomputed_clusters <- get("OPTIMAL_CLUSTER_COUNTS", envir = globalenv())
  }

  # Use indices from MESMA_PARAMS (includes L2norm if enabled)
  indices_for_final_library <- if (!is.null(MESMA_PARAMS$indices)) MESMA_PARAMS$indices else avail

  mesma_lib <- build_mesma_library_weighted(df_train, indices_for_final_library, MESMA_PARAMS, ALLOWED_VEG, precomputed_clusters = precomputed_clusters)

  cat("Pre-computing optimized library for MESMA...\n")
  OPTIMIZED_LIBRARY <- precompute_optimized_library_weighted(mesma_lib, grid_type = "full", feature_weights = (if (!is.null(MESMA_PARAMS) && !is.null(MESMA_PARAMS$weights)) MESMA_PARAMS$weights else NULL) )

  # Summarize pruning if any
  if (exists("PRUNE_ZERO_WEIGHT_FEATURES") && isTRUE(PRUNE_ZERO_WEIGHT_FEATURES)) {
    total_pruned <- 0
    total_kept <- 0
    vegs_pruned <- 0
    for (v in names(OPTIMIZED_LIBRARY)) {
      libv <- OPTIMIZED_LIBRARY[[v]]
      if (!is.null(libv$pruned_info) && !is.null(libv$pruned_info$n_kept)) {
        total_kept <- total_kept + libv$pruned_info$n_kept
        vegs_pruned <- vegs_pruned + 1
      }
    }
    if (vegs_pruned > 0) {
      cat(sprintf("[INFO] Pruning summary: %d vegetation types pruned; avg kept features per pruned veg: %.1f\n", vegs_pruned, total_kept / vegs_pruned))
    }
  }

  assign("MESMA_PARAMS", MESMA_PARAMS, envir = globalenv())
  assign("mesma_lib", mesma_lib, envir = globalenv())
  assign("OPTIMIZED_LIBRARY", OPTIMIZED_LIBRARY, envir = globalenv())

  cat(sprintf("[NOTICE] MESMA feature count: avail=%d, params_indices=%d\n",
              length(avail), length(MESMA_PARAMS$indices)))
  cat(sprintf("[NOTICE] MESMA feature list: %s\n", paste(MESMA_PARAMS$indices, collapse=", ")))

  compressed_templates_accessor <- precompute_compressed_templates(mesma_lib, "full")
  assign("compressed_templates_accessor", compressed_templates_accessor, envir = globalenv())
  assign(".COMPRESSED_TEMPLATES_ACCESSOR", compressed_templates_accessor, envir = globalenv())

  # === Shared confusion matrix function (used for both training and validation) ===
  # Pipeline: OPTIMIZED_LIBRARY endmembers -> energy normalization -> batch FCLS
  # with feature weights. Observations: z-score -> prune to match library.
  compute_confusion_matrix <- function(df_data, label = "Training", interpolate = NULL) {
    # Always echo an attempt to compute, even if we exit early.
    if (is.null(interpolate)) {
        if (!is.null(MESMA_PARAMS$interpolate_inference)) {
            interpolate <- MESMA_PARAMS$interpolate_inference
        } else {
            interpolate <- INTERPOLATE_INFERENCE
        }
    }
    interpolate <- get_interpolate_method(interpolate)
    cat(sprintf("[CONFUSION MATRIX] Starting '%s' (interpolate=%s)\n", label, interpolate))

    if (is.null(df_data) || nrow(df_data) == 0) {
      cat(sprintf("[CONFUSION MATRIX] %s skipped: no data rows\n", label))
      return(invisible(NULL))
    }
    if (!exists("OPTIMIZED_LIBRARY") || is.null(OPTIMIZED_LIBRARY)) {
      cat(sprintf("[CONFUSION MATRIX] %s skipped: library not available\n", label))
      return(invisible(NULL))
    }
    if (is.null(MESMA_PARAMS) || is.null(MESMA_PARAMS$weights)) {
      cat(sprintf("[CONFUSION MATRIX] %s skipped: MESMA_PARAMS or weights missing\n", label))
      return(invisible(NULL))
    }

    # always log the method so that behaviour is clear even in quiet runs
    cat(sprintf("[NOTICE] %s confusion matrix: temporal fill method = %s\n",
                label, interpolate))

    tryCatch({
      df_cm <- dplyr::mutate(df_data, target_class = tolower(as.character(Veg)))
      df_cm <- dplyr::filter(df_cm, !is.na(target_class) & target_class != "")

      unique_classes <- names(OPTIMIZED_LIBRARY)
      veg_classes <- setdiff(unique_classes, "barren")

      cat(sprintf("[CONFUSION MATRIX] Using %d rows from %d locations\n",
                  nrow(df_cm), length(unique(df_cm$location_id))))

      # 1. Build endmember matrix E from OPTIMIZED_LIBRARY (multi-prototype per class)
      E_cols <- list()
      col_class_labels <- c()
      pruned_info <- NULL
      for (cls in unique_classes) {
        lib_cls <- OPTIMIZED_LIBRARY[[cls]]
        if (is.null(lib_cls) || is.null(lib_cls$M)) next
        if (is.null(pruned_info)) pruned_info <- lib_cls$pruned_info
        M_cls <- lib_cls$M
        if (is.null(dim(M_cls))) M_cls <- matrix(M_cls, nrow = 1)
        for (ri in seq_len(nrow(M_cls))) {
          E_cols[[length(E_cols) + 1]] <- M_cls[ri, ]
          col_class_labels <- c(col_class_labels, cls)
        }
      }
      if (length(E_cols) < 2) {
        cat(sprintf("[CONFUSION MATRIX] %s skipped: not enough endmember columns (found %d)\n",
                    label, length(E_cols)))
        return(invisible(NULL))
      }

      E <- do.call(cbind, E_cols)

      # 2. Build observation vectors (z-score + prune)
      base_indices_cm <- if (!is.null(MESMA_PARAMS$base_indices)) MESMA_PARAMS$base_indices else MESMA_PARAMS$indices
      indices_cm <- MESMA_PARAMS$indices
      n_base_idx_cm <- length(base_indices_cm)
      l2_normalize_cm <- isTRUE(MESMA_PARAMS$l2_normalize)


      traces <- unique(df_cm[, c("location_id", "pheno_year", "target_class")])
      obs_vecs <- list()
      obs_labels <- c()

      for (j in seq_len(nrow(traces))) {
        lid <- traces$location_id[j]
        pyr <- traces$pheno_year[j]
        true_cls <- traces$target_class[j]

        # Filter by target_class too, so the pentad matrix is built from
        # single-class rows only (matches build_storage_from_df behaviour).
        sub <- df_cm[df_cm$location_id == lid & df_cm$pheno_year == pyr & df_cm$target_class == true_cls, ]
        mat <- build_pentad_matrix(sub, base_indices_cm, interpolate = interpolate)
        if (!is.null(mat)) {
          # Ensure matrix columns align exactly with indices_cm (order and selection)
          # This fixes potential mismatches where base_indices differs from indices_cm.
          if (all(indices_cm %in% colnames(mat))) {
             mat <- mat[, indices_cm, drop=FALSE]
             n_proc_idx <- length(indices_cm)
          } else {
             n_proc_idx <- n_base_idx_cm
          }

          vec_raw <- as.numeric(mat)
          vec_raw <- mesma_apply_representation_vec(vec_raw, n_proc_idx, TEMPORAL_BUDGET, l2_normalize_cm)
          vec <- mesma_zscore_vec_by_index(vec_raw, indices_cm, MESMA_PARAMS$means, MESMA_PARAMS$sds, TEMPORAL_BUDGET)
          if (j <= 2 && grepl("Training", label)) {
             idx1_name <- if (length(indices_cm) >= 1) indices_cm[1] else "?"
             mu1 <- if (idx1_name %in% names(MESMA_PARAMS$means)) MESMA_PARAMS$means[[idx1_name]] else NA
             sd1 <- if (idx1_name %in% names(MESMA_PARAMS$sds))   MESMA_PARAMS$sds[[idx1_name]]   else NA
             cat(sprintf("\n[DEBUG %s] Sample %d (True=%s)\n", label, j, true_cls))
             cat(sprintf("  Index[1]='%s'  mean=%.4f  sd=%.4f\n", idx1_name, mu1, sd1))
             cat(sprintf("  RAW Vec[1:5] (pre-zscore): [%s]\n",
                         paste(round(vec_raw[1:min(5,length(vec_raw))],4), collapse=",")))
             cat(sprintf("  Z-scored Vec[1:5]:          [%s]\n",
                         paste(round(vec[1:min(5,length(vec))],3), collapse=",")))
             # Reverse-engineer what the library E[,1] was before z-scoring
             E_raw_approx <- E[1:min(5,nrow(E)),1] * sd1 + mu1
             cat(sprintf("  E[,1] class='%s'  z-scored[1:5]=[%s]  back-calc raw=[%s]\n",
                         col_class_labels[1],
                         paste(round(E[1:min(5,nrow(E)),1],3), collapse=","),
                         paste(round(E_raw_approx,4), collapse=",")))
          }
          if (!is.null(pruned_info) && !is.null(pruned_info$kept_idx)) {
            vec <- vec[pruned_info$kept_idx]
          }
          obs_vecs[[length(obs_vecs) + 1]] <- vec
          obs_labels <- c(obs_labels, true_cls)
        }
      }
      if (length(obs_vecs) == 0) {
        cat(sprintf("[CONFUSION MATRIX] %s skipped: no observation vectors constructed (interpolate=%s)\n",
                    label, interpolate))
        return(invisible(NULL))
      }

      Y <- do.call(rbind, obs_vecs)

      # 3. FCLS with pruned weights
      # MESMA_PARAMS is already trimmed to the kept-index feature space.
      # pruned_info$kept_idx (if set) reflects per-bin pruning *within* that space.
      w <- MESMA_PARAMS$weights
      lda_basis_cm <- MESMA_PARAMS$lda_basis
      lda_cw_cm    <- MESMA_PARAMS$lda_component_weights
      if (!is.null(pruned_info) && !is.null(pruned_info$kept_idx)) {
        w <- w[pruned_info$kept_idx]
        if (!is.null(lda_basis_cm)) lda_basis_cm <- lda_basis_cm[pruned_info$kept_idx, , drop = FALSE]
      }
      solver_mode <- if (!is.null(MESMA_PARAMS$solver)) MESMA_PARAMS$solver else "fcls"
      all_coefs <- solve_batch_fcls(E, Y,
                                    w,
                                    lda_basis             = lda_basis_cm,
                                    lda_component_weights = lda_cw_cm,
                                    solver = solver_mode)
      if (is.null(all_coefs)) return(invisible(NULL))

      # 4. Aggregate coefficients by class and build confusion matrix
      cm_labels <- unique_classes
      frac_mat <- matrix(0, nrow = length(cm_labels), ncol = length(cm_labels))
      rownames(frac_mat) <- cm_labels; colnames(frac_mat) <- cm_labels
      class_counts <- setNames(rep(0, length(cm_labels)), cm_labels)
      correct_veg <- 0; total_veg <- 0

      for (j in seq_len(nrow(Y))) {
        true_cls <- obs_labels[j]
        coefs <- tapply(all_coefs[j, ], col_class_labels, sum)
        for (cc in cm_labels) if (!(cc %in% names(coefs))) coefs[[cc]] <- 0
        
        # Standard abundance accumulation
        if (true_cls %in% cm_labels) {
          for (cc in cm_labels) frac_mat[true_cls, cc] <- frac_mat[true_cls, cc] + coefs[[cc]]
          class_counts[true_cls] <- class_counts[true_cls] + 1
        }
        
        # Hard classification stats
        pred_cls <- names(which.max(coefs))
        if (true_cls %in% veg_classes) {
          total_veg <- total_veg + 1
          if (pred_cls == true_cls) correct_veg <- correct_veg + 1
        }
      }

      # 5. Normalize and display
      avg_pred_frac <- frac_mat
      for (r in seq_len(nrow(frac_mat))) {
        if (class_counts[rownames(frac_mat)[r]] > 0)
          avg_pred_frac[r, ] <- frac_mat[r, ] / class_counts[rownames(frac_mat)[r]]
      }
      row_sums <- rowSums(avg_pred_frac, na.rm = TRUE)
      for (ri in seq_len(nrow(avg_pred_frac))) {
        if (row_sums[ri] > 0) avg_pred_frac[ri, ] <- avg_pred_frac[ri, ] / row_sums[ri]
      }

      cat(sprintf("\n[CONFUSION MATRIX] %s Predicted Fractions - Row Normalized (All Classes):\n", label))
      print(round(avg_pred_frac, 3))

      veg_rows <- intersect(veg_classes, rownames(avg_pred_frac))
      veg_cols <- intersect(veg_classes, colnames(avg_pred_frac))
      avg_cor_pred_veg_frac <- NA
      if (length(veg_rows) > 0 && length(veg_cols) > 0) {
        # Subset veg classes and re-normalize locally
        veg_matrix <- avg_pred_frac[veg_rows, veg_cols, drop = FALSE]
        row_sums_veg <- rowSums(veg_matrix, na.rm = TRUE)
        for (ri in seq_len(nrow(veg_matrix))) {
          if (row_sums_veg[ri] > 0) veg_matrix[ri, ] <- veg_matrix[ri, ] / row_sums_veg[ri]
        }

        cat(sprintf("\n[CONFUSION MATRIX] %s Predicted Fractions - Vegetation Only (Row Normalized):\n", label))
        print(round(veg_matrix, 3))

        veg_diag_fracs <- sapply(veg_rows, function(cls) {
          if (cls %in% colnames(veg_matrix)) veg_matrix[cls, cls] else NA
        })
        avg_cor_pred_veg_frac <- mean(veg_diag_fracs, na.rm = TRUE)
      }

      veg_class_accuracy <- if (total_veg > 0) correct_veg / total_veg else NA
      cat(sprintf("\n[CONFUSION MATRIX] %s Avg Row-Normalized Veg Fraction (diagonal): %.1f%%\n", label, 100 * avg_cor_pred_veg_frac))
      cat(sprintf("[CONFUSION MATRIX] %s Veg Classification Accuracy: %.1f%% (%d/%d correct)\n",
                  label, 100 * veg_class_accuracy, correct_veg, total_veg))

      # --- Save numeric matrices and make PNG visualizations if OUTPUT_DIR exists or can be created ---
      if (exists("OUTPUT_DIR") && !is.null(OUTPUT_DIR)) {
        if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
        lab_clean <- tolower(gsub("[^A-Za-z0-9]+", "_", label))
        all_csv <- file.path(OUTPUT_DIR, sprintf("%s_confusion_matrix_all.csv", lab_clean))
        write.csv(avg_pred_frac, file = all_csv, row.names = TRUE)
        cat(sprintf("[SAVE] Wrote %s\n", all_csv))
        if (length(veg_rows) > 0 && length(veg_cols) > 0) {
          veg_csv <- file.path(OUTPUT_DIR, sprintf("%s_confusion_matrix_veg.csv", lab_clean))
          write.csv(veg_matrix, file = veg_csv, row.names = TRUE)
          cat(sprintf("[SAVE] Wrote %s\n", veg_csv))
        }

        if (requireNamespace("ggplot2", quietly = TRUE)) {
          plot_mat <- function(mat, title_text, fname) {
            df_plot <- as.data.frame(as.table(mat))
            names(df_plot) <- c("True", "Predicted", "Frac")
            p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Predicted, y = True, fill = Frac)) +
                 ggplot2::geom_tile() +
                 ggplot2::geom_text(ggplot2::aes(label = round(Frac,3)), color = "black", size = 3) +
                 ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
                 ggplot2::theme_minimal() +
                 ggplot2::labs(title = title_text, x = "Predicted", y = "True")
            ggplot2::ggsave(fname, plot = p, width = 6, height = 5, dpi = 150)
            cat(sprintf("[SAVE] Wrote confusion matrix plot to %s\n", fname))
          }
          png_all <- file.path(OUTPUT_DIR, sprintf("%s_confusion_matrix_all.png", lab_clean))
          plot_mat(avg_pred_frac, sprintf("%s Confusion Matrix", label), png_all)
          if (length(veg_rows) > 0 && length(veg_cols) > 0) {
            png_veg <- file.path(OUTPUT_DIR, sprintf("%s_confusion_matrix_veg.png", lab_clean))
            plot_mat(veg_matrix, sprintf("%s Veg Confusion Matrix", label), png_veg)
          }
        } else {
          cat("[WARN] ggplot2 not available; skipping confusion matrix plots\n")
        }
        # Save class counts (n_{i.}) so callers can compute SE for bias correction
        counts_csv <- file.path(OUTPUT_DIR, sprintf("%s_class_counts.csv", lab_clean))
        write.csv(data.frame(class_name = names(class_counts), n = as.integer(class_counts)),
                  file = counts_csv, row.names = FALSE)
        cat(sprintf("[SAVE] Wrote %s\n", counts_csv))
      }
      list(norm_mat = avg_pred_frac, class_counts = class_counts)
    }, error = function(e) {
      cat(sprintf("[CONFUSION MATRIX] Error computing %s confusion matrix: %s\n", label, e$message))
      NULL
    })
  }

  # === STEP 4: Compute Training Confusion Matrix ===
  # Use df_train_model (the data used to build the library) for internal consistency.
  # Fall back to df_train if df_train_model is unavailable.
  .cm_train_data <- if (exists("df_train_model", envir = globalenv()) &&
                         !is.null(df_train_model) && nrow(df_train_model) > 0) {
    cat("[CONFUSION MATRIX] Using df_train_model for training confusion matrix (library-consistent)\n")
    df_train_model
  } else if (exists("df_train") && !is.null(df_train) && nrow(df_train) > 0) {
    cat("[CONFUSION MATRIX] df_train_model not found, falling back to df_train\n")
    df_train
  } else {
    NULL
  }

  if (!is.null(.cm_train_data) && nrow(.cm_train_data) > 0) {
    cat("\n=== STEP 4: Computing Confusion Matrix with Final Weights (All Training Data) ===\n")
    cat("[CONFUSION MATRIX] computing training interpolated matrix\n")
    compute_confusion_matrix(.cm_train_data, "Training (interpolated)", interpolate = "linear")

    cat("[CONFUSION MATRIX] computing training non-interpolated matrix\n")
    compute_confusion_matrix(.cm_train_data, "Training (non-interpolated)", interpolate = "none")
  }

  # --- early validation confusion matrix for user visibility -------------
  validation_confusion_done <- FALSE
  if (exists("df_validation") && !is.null(df_validation) && nrow(df_validation) > 0) {
    cat("\n=== VALIDATION CONFUSION MATRIX (pre-processing) ===\n")
    # always print both versions so user can compare directly
    cat("[CONFUSION MATRIX] computing validation interpolated matrix\n")
    compute_confusion_matrix(df_validation, "Validation (interpolated)", interpolate = "linear")
    cat("[CONFUSION MATRIX] computing validation non-interpolated matrix\n")
    compute_confusion_matrix(df_validation, "Validation (non-interpolated)", interpolate = FALSE)
    validation_confusion_done <- TRUE
  } else {
    cat("\n[WARNING] No validation data available to compute confusion matrix early\n")
  }


  # Do not exit yet: we want to compute validation confusion too when in
  # VALIDATE_ONLY mode.  The actual early termination will happen after the
  # validation confusion block below (now an additional sanity check).

  

# ==========================================================================
# Training is performed externally — this script focuses on inference & validation.
# validation_location_ids already set during stratified train/validation split
# NOTE: do NOT create or populate `df_tasks` here — inference code will build task tables from the inference input.
# ==========================================================================

# ============================================================================
# OOB FRACTION PREDICTION UNCERTAINTY
# Uses OOB validation residuals to add realistic prediction error to MC draws.
# For each true class, we store the distribution of prediction residuals
# (predicted_fraction - true_fraction) observed in OOB validation.
# During MC, we sample from these residuals to perturb the predicted fractions.
# ============================================================================

# Global storage for OOB fraction residuals per true class
.OOB_FRACTION_RESIDUALS <- NULL  # List: true_class -> matrix of residuals (n_samples x n_classes)


# Sample a residual vector for a given dominant predicted class
# Returns a named vector of residuals to add to predicted fractions
sample_oob_residual <- function(dominant_class) {
  if (!exists(".OOB_FRACTION_RESIDUALS", envir = globalenv())) {
    return(NULL)
  }

  residuals_by_class <- get(".OOB_FRACTION_RESIDUALS", envir = globalenv())
  if (is.null(residuals_by_class)) return(NULL)

  # Match dominant class (case-insensitive)
  cls_lower <- tolower(dominant_class)
  available_classes <- tolower(names(residuals_by_class))
  match_idx <- which(available_classes == cls_lower)

  if (length(match_idx) == 0) {
    return(NULL)
  }

  resid_mat <- residuals_by_class[[match_idx]]
  if (is.null(resid_mat) || nrow(resid_mat) == 0) {
    return(NULL)
  }

  # Sample one row uniformly at random
  row_idx <- sample.int(nrow(resid_mat), 1)
  residual <- as.numeric(resid_mat[row_idx, ])
  names(residual) <- colnames(resid_mat)

  return(residual)
}


  fit_one_task <- function(task_data) {
    if (is.null(task_data) || nrow(task_data) == 0) return(NULL)

    loc <- as.character(task_data$location_id[1])
    yr <- as.integer(task_data$pheno_year[1])
    

    # Extract lat/lon if available
    lat_val <- NA_real_
    lon_val <- NA_real_
    if ("lat" %in% names(task_data)) lat_val <- as.numeric(task_data$lat[1]) else if ("latitude" %in% names(task_data)) lat_val <- as.numeric(task_data$latitude[1])
    if ("lon" %in% names(task_data)) lon_val <- as.numeric(task_data$lon[1]) else if ("longitude" %in% names(task_data)) lon_val <- as.numeric(task_data$longitude[1])

    # Calculate summer observations
    if ("month" %in% names(task_data)) {
      summer_months <- 6:9
      n_summer_valid <- sum(task_data$month %in% summer_months, na.rm = TRUE)
    } else if ("date" %in% names(task_data)) {
      task_data$month <- as.integer(format(task_data$date, "%m"))
      summer_months <- 6:9
      n_summer_valid <- sum(task_data$month %in% summer_months, na.rm = TRUE)
    } else {
      n_summer_valid <- NA_integer_
    }
    low_data_flag <- isTRUE(n_summer_valid < 10)

    # Use MESMA parameters
    PARAMS <- MESMA_PARAMS
    
    if (is.null(PARAMS) || is.null(PARAMS$indices) || is.null(PARAMS$means) || is.null(PARAMS$sds)) {
      cat(sprintf("[ERROR] loc=%s yr=%d: PARAMS structure is incomplete (PARAMS=%s, indices=%s, means=%s, sds=%s)\n",
                  loc, yr,
                  !is.null(PARAMS), 
                  !is.null(PARAMS$indices), 
                  !is.null(PARAMS$means), 
                  !is.null(PARAMS$sds)))
      return(NULL)
    }
    
    # For inference/validation data: temporal filling is disabled by default to avoid
    # filling gaps in held-out data.  The option can be set to
    # "linear"/TRUE, "whittaker" or "none" via the global `INTERPOLATE_INFERENCE`
    # flag or `PARAMS$interpolate_inference`.
    base_indices_for_build <- if (!is.null(PARAMS$base_indices)) PARAMS$base_indices else PARAMS$indices
    interp_method_for_inference <- if (!is.null(PARAMS$interpolate_inference)) {
      get_interpolate_method(PARAMS$interpolate_inference)
    } else {
      get_interpolate_method(INTERPOLATE_INFERENCE)
    }
    cat(sprintf("[NOTICE] temporal fill method for inference/validation: %s\n", interp_method_for_inference))
    raw_mat <- build_pentad_matrix(task_data, base_indices_for_build, interpolate = interp_method_for_inference)
    if (is.null(raw_mat)) {
      return(NULL)
    }

    y_raw <- as.numeric(raw_mat)

    # Detect if library features were pruned (PRUNE_ZERO_WEIGHT_FEATURES)
    # so we can apply the same pruning after z-scoring
    .pruned_kept_idx <- NULL
    if (exists("OPTIMIZED_LIBRARY") && !is.null(OPTIMIZED_LIBRARY)) {
      for (vv in names(OPTIMIZED_LIBRARY)) {
        pi <- OPTIMIZED_LIBRARY[[vv]]$pruned_info
        if (!is.null(pi) && !is.null(pi$kept_idx)) { .pruned_kept_idx <- pi$kept_idx; break }
      }
    }

    # === CREATE VALIDITY MASK ===
    # Mask for valid (non-NA) observations - DO NOT replace NA with 0
    valid_mask <- is.finite(y_raw)
    n_valid <- sum(valid_mask)

    # Require at least one valid observation
    if (n_valid == 0) {
      return(NULL)
    }

    # ===== MESMA UNMIXING =====
      # Get indices from params
      n_bins <- TEMPORAL_BUDGET
      indices <- PARAMS$indices
      base_indices <- if (!is.null(PARAMS$base_indices)) PARAMS$base_indices else indices
      n_base_idx <- length(base_indices)
      l2_normalize <- isTRUE(PARAMS$l2_normalize)

      # Build feature vector based on representation mode
      y_work <- mesma_apply_representation_vec(y_raw, n_base_idx, n_bins, l2_normalize)
      y_norm <- mesma_zscore_vec_by_index(y_work, indices, PARAMS$means, PARAMS$sds, n_bins)

      # --- LDA-space projection ---
      # When USE_LDA_SPACE_SOLVER is TRUE the observation and every library
      # template are projected into the LDA discriminant subspace *before*
      # masking.  After projection the feature dimension equals the number of
      # LDA components; valid_mask must therefore be recomputed in that space.
      lda_basis_inf <- PARAMS$lda_basis
      lda_cw_inf    <- PARAMS$lda_component_weights

      if (!is.null(lda_basis_inf) && ncol(lda_basis_inf) >= 1) {
        # lda_basis_inf: (n_features x n_components)
        # Apply pruning to basis rows first if the library was pruned
        lda_basis_use <- lda_basis_inf
        if (!is.null(.pruned_kept_idx) && nrow(lda_basis_use) > length(.pruned_kept_idx)) {
          lda_basis_use <- lda_basis_use[.pruned_kept_idx, , drop = FALSE]
          y_norm_for_proj <- y_norm[.pruned_kept_idx]
          valid_mask      <- valid_mask[.pruned_kept_idx]
        } else {
          y_norm_for_proj <- y_norm
        }
        # Zero-out missing pentads in obs before projection so they don't contribute
        y_norm_for_proj[!valid_mask] <- 0
        lda_basis_use[!is.finite(lda_basis_use)] <- 0
        y_norm_full_pca_lda <- as.numeric(t(lda_basis_use) %*% y_norm_for_proj)
        # In LDA space every component is always "observed"
        valid_mask    <- rep(TRUE, length(y_norm_full_pca_lda))
        weights_for_mask <- if (!is.null(lda_cw_inf)) lda_cw_inf else rep(1, length(y_norm_full_pca_lda))
        # Store the (pruned) basis for template projection below
        .lda_basis_for_templates <- lda_basis_use
      } else {
        # No LDA projection: stay in original feature space
        if (!is.null(.pruned_kept_idx) && length(.pruned_kept_idx) < length(y_norm)) {
          y_norm_full_pca_lda <- y_norm[.pruned_kept_idx]
          valid_mask          <- valid_mask[.pruned_kept_idx]
        } else {
          y_norm_full_pca_lda <- y_norm
        }
        weights_for_mask <- PARAMS$weights
        if (!is.null(weights_for_mask) && !is.null(.pruned_kept_idx) &&
            length(weights_for_mask) > length(.pruned_kept_idx)) {
          weights_for_mask <- weights_for_mask[.pruned_kept_idx]
        }
        .lda_basis_for_templates <- NULL
      }

      # Mask the observation
      y_norm_masked <- y_norm_full_pca_lda[valid_mask]

      # Weights to pass to solver
      if (!is.null(weights_for_mask) && length(weights_for_mask) == length(y_norm_full_pca_lda)) {
        weights_masked <- weights_for_mask[valid_mask]
      } else {
        weights_masked <- rep(1, length(y_norm_masked))
      }

      # Check if we have sufficient signal (using unweighted norm for validation)
      y_norm_val <- sqrt(sum(y_norm_masked^2, na.rm = TRUE))

      if (is.na(y_norm_val) || y_norm_val < 1e-9) {
        return(NULL)
      }

      # Use the masked, z-score normalized observation (NOT pre-weighted)
      # Weighting will be applied inside solve_weights_fcls via feature_weights parameter
      y_for_unmixing <- y_norm_masked

      # Perform MESMA using all endmembers (barren + all veg types)
      # Barren is kept as from MESMA - no replacement
      veg_kept <- names(mesma_lib)
      # SAFEGUARD: Ensure a 'barren' endmember exists
      if (!"barren" %in% veg_kept) {
        warning("No 'barren' endmember found in MESMA library - ensure barren endmembers are present.")
      }
      
      if (is.null(OPTIMIZED_LIBRARY) || length(OPTIMIZED_LIBRARY) == 0) {
        cat(sprintf("[ERROR] loc=%s yr=%d: OPTIMIZED_LIBRARY is NULL or empty\n", loc, yr))
        return(NULL)
      }
      
      top_variants <- list()

      for(v in veg_kept) {
        lib <- OPTIMIZED_LIBRARY[[v]]
        if(is.null(lib)) {
          next
        }

        # Build masked templates from raw library matrix, applying the SAME PCA-LDA transform
        # used in training (if available) so the observation and templates lie in identical
        # feature space before similarity ranking and unmixing.
        # NOTE: lib$M is already pruned (if PRUNE_ZERO_WEIGHT_FEATURES is enabled)
        # during precompute_optimized_library_weighted, so no need to prune again here
        M_full <- lib$M

        # Project templates into LDA space if active, otherwise keep as-is
        if (!is.null(.lda_basis_for_templates)) {
          # M_full: (n_variants x n_features); basis: (n_features x n_components)
          # Result: (n_variants x n_components)
          M_full_trans <- M_full %*% .lda_basis_for_templates
        } else {
          M_full_trans <- M_full
        }

        # Mask to the valid observation columns (all TRUE in LDA space)
        lib_M_masked <- M_full_trans[, valid_mask, drop = FALSE]

        # Guard against NA/Inf in templates after masking
        lib_M_masked[!is.finite(lib_M_masked)] <- 0

        lib_ids_kept <- lib$ids

        if (nrow(lib_M_masked) == 0) {
          next
        }

        # L2-normalize each row (used for both similarity ranking and unmixing)
        if (nrow(lib_M_masked) == 1) {
          row_norm <- sqrt(sum(lib_M_masked^2, na.rm = TRUE))
          lib_M_norm_masked <- if (is.finite(row_norm) && row_norm >= 1e-9) lib_M_masked / row_norm else lib_M_masked
        } else {
          lib_M_norm_masked <- t(apply(lib_M_masked, 1, function(row) {
            row_norm <- sqrt(sum(row^2, na.rm = TRUE))
            if (!is.finite(row_norm) || row_norm < 1e-9) row else row / row_norm
          }))
        }

        # Compute similarities using weighted vectors (for ranking only)
        w_masked <- weights_masked
        if (is.null(w_masked) || length(w_masked) != length(y_for_unmixing)) w_masked <- rep(1, length(y_for_unmixing))
        y_for_sim <- y_for_unmixing * sqrt(pmax(w_masked, 0))
        denom_sim <- sqrt(sum(y_for_sim^2, na.rm = TRUE)); if (denom_sim < 1e-12) denom_sim <- 1
        y_sim_norm <- y_for_sim / denom_sim
        sims <- as.numeric(lib_M_norm_masked %*% y_sim_norm)
        # Include all variants for consideration (remove similarity-based top-K filtering)
        best_idx <- order(sims, decreasing=TRUE)
        
        top_variants[[v]] <- lapply(best_idx, function(i) {
          # Use L2-normalized masked template for unmixing
          masked_vec <- lib_M_norm_masked[i, ]
          list(vec = masked_vec, id = lib_ids_kept[i], similarity = sims[i])
        })
      }

      # Remove empty vegetation types
      empty_vegs <- names(top_variants)[sapply(top_variants, function(x) is.null(x) || length(x) == 0)]
      if (length(empty_vegs) > 0) {
        for (ev in empty_vegs) top_variants[[ev]] <- NULL
      }
      

      if (length(top_variants) == 0) {
        return(NULL)
      }
      
      # Evaluate all combinations of endmembers
      # Pass unweighted observation; weighting happens inside solve_weights_fcls
      
      # decide penalty strength (may be zero)
      sparse_lambda <- 0
      if (exists("SPARSE_MIXING_LAMBDA", inherits = TRUE) && is.finite(get("SPARSE_MIXING_LAMBDA", inherits = TRUE))) {
        sparse_lambda <- get("SPARSE_MIXING_LAMBDA", inherits = TRUE)
      }
      if (sparse_lambda != 0) {
        cat(sprintf("[INFO] Applying sparse unmixing penalty lambda=%.6f (barren included)\n", sparse_lambda))
      }
      active_solver <- if (!is.null(PARAMS$solver)) PARAMS$solver else "fcls"
      if (!isTRUE(getOption(".fit_one_task_solver_logged"))) {
        cat(sprintf("[SOLVER] fit_one_task using solver: %s\n", active_solver))
        options(.fit_one_task_solver_logged = TRUE)
      }
      best_result <- tryCatch({
        evaluate_all_combinations(
          y_for_unmixing,
          top_variants,
          lambda = sparse_lambda,
          feature_weights = weights_masked,
          solver = active_solver
        )
      }, error = function(e) {
        cat(sprintf("[ERROR fit_one_task] evaluate_all_combinations failed for loc=%s year=%d: %s\n", as.character(loc), as.integer(yr), e$message))
        NULL
      })
      
      if (is.null(best_result)) {
        return(NULL)
      }
      

      # Extract coefficients and create output
      chosen_ids <- best_result$ids
      coefs <- best_result$w
      rmse <- best_result$rmse
      residuals <- best_result$residuals
      E_best_masked <- if (!is.null(best_result$E_best)) best_result$E_best else NULL

      
      # Calculate model selection uncertainty from top models
      coef_sd_vec <- rep(NA, length(coefs))
      names(coef_sd_vec) <- names(coefs)
      
      if (!is.null(best_result$top_models) && length(best_result$top_models) > 1) {
         # Stack weights from all top models: n_veg x n_models
         w_matrix <- do.call(cbind, lapply(best_result$top_models, function(m) m$w))
         
         if (!is.null(w_matrix) && ncol(w_matrix) > 1) {
            # Compute SD for each vegetation type across the top models
            coef_sd_vec <- apply(w_matrix, 1, sd, na.rm=TRUE)
         }
      }

      # --- Monte Carlo error propagation through unmixing (optional) ---
      mc_coef_sd <- NULL
      mc_coef_q025 <- NULL
      mc_coef_q975 <- NULL

      run_monte_carlo <- isTRUE(ENABLE_UNCERTAINTY) && exists("ENABLE_MONTE_CARLO") && isTRUE(ENABLE_MONTE_CARLO)
      n_mc <- if (exists("MC_N_DRAWS") && is.finite(MC_N_DRAWS)) as.integer(MC_N_DRAWS) else 0L

      # E_best_masked may have fewer columns than length(coefs) when a veg subset
      # was selected (excluded types get coef=0 but no endmember column).
      # Identify active subset so MC draws match the endmember matrix dimensions.
      mc_veg_subset <- if (!is.null(best_result$veg_subset)) best_result$veg_subset else names(coefs)
      n_mc_subset <- length(mc_veg_subset)

      if (run_monte_carlo && n_mc >= 10 && !is.null(E_best_masked) && is.matrix(E_best_masked) && ncol(E_best_masked) == n_mc_subset) {
        # Noise model: y is z-scored; use residual scale by default.
        mc_sigma <- rmse
        if (exists("MC_NOISE_SD") && is.finite(MC_NOISE_SD) && MC_NOISE_SD > 0) mc_sigma <- as.numeric(MC_NOISE_SD)
        if (exists("MC_NOISE_SCALE") && is.finite(MC_NOISE_SCALE) && MC_NOISE_SCALE > 0) mc_sigma <- mc_sigma * as.numeric(MC_NOISE_SCALE)
        if (!is.finite(mc_sigma) || mc_sigma < 0) mc_sigma <- 0

        # Optional: endmember bundles (sample endmembers from per-band variance estimated from the candidate pool)
        enable_bundles <- exists("ENABLE_ENDMEMBER_BUNDLES") && isTRUE(ENABLE_ENDMEMBER_BUNDLES)
        bundle_params <- NULL
        if (enable_bundles) {
          bundle_params <- list()
          verbose_bundle <- FALSE

          for (v in mc_veg_subset) {
            cand_list <- top_variants[[v]]
            if (is.null(cand_list) || length(cand_list) < 3) next
            Mv <- do.call(rbind, lapply(cand_list, function(z) as.numeric(z$vec)))
            if (is.null(Mv) || nrow(Mv) < 3) next

            bundle_result <- compute_bundle_covariance(Mv, verbose = verbose_bundle)

            if (!is.null(bundle_result) && !is.null(bundle_result$A)) {
              bundle_params[[v]] <- list(
                mu = bundle_result$mu,
                A = bundle_result$A
              )
            }
          }
        }

        sample_from_bundle <- function(mu, A) {
          z <- stats::rnorm(length(mu))
          x <- as.numeric(mu + A %*% z)
          x[!is.finite(x)] <- 0
          nrm <- sqrt(sum(x^2))
          if (!is.finite(nrm) || nrm < 1e-9) return(mu)
          x / nrm
        }

        # Check if OOB fraction uncertainty is enabled and residuals are available
        enable_oob_frac_uncertainty <- exists("ENABLE_OOB_FRACTION_UNCERTAINTY") &&
                                        isTRUE(ENABLE_OOB_FRACTION_UNCERTAINTY) &&
                                        exists(".OOB_FRACTION_RESIDUALS", envir = globalenv())

        # Determine dominant class for OOB residual sampling (class with highest coefficient)
        dominant_class <- NULL
        if (enable_oob_frac_uncertainty && length(coefs) > 0) {
          dominant_class <- names(which.max(coefs))
        }

        # w_draws sized to full coef vector; inactive types stay 0
        w_draws <- matrix(NA_real_, nrow = length(coefs), ncol = n_mc)
        rownames(w_draws) <- names(coefs)

        for (b in seq_len(n_mc)) {
          y_mc <- y_for_unmixing
          if (length(mc_sigma) == 1 && is.finite(mc_sigma) && mc_sigma > 0) {
            y_mc <- y_mc + stats::rnorm(length(y_mc), mean = 0, sd = mc_sigma)
          }

          E_mc <- E_best_masked
          if (enable_bundles && !is.null(bundle_params) && length(bundle_params) > 0) {
            for (j in seq_along(mc_veg_subset)) {
              v <- mc_veg_subset[j]
              bp <- bundle_params[[v]]
              if (!is.null(bp) && !is.null(bp$mu) && !is.null(bp$A)) {
                E_mc[, j] <- sample_from_bundle(bp$mu, bp$A)
              }
            }
          }

          res_mc <- solve_weights_fcls(E_mc, y_mc, feature_weights = weights_masked)
          if (!is.null(res_mc) && !is.null(res_mc$w) && length(res_mc$w) == n_mc_subset) {
            # Expand subset weights back to full veg vector
            w_mc_full <- rep(0, length(coefs))
            names(w_mc_full) <- names(coefs)
            w_mc_full[mc_veg_subset] <- as.numeric(res_mc$w)

            # Apply OOB fraction residual perturbation if enabled
            if (enable_oob_frac_uncertainty && !is.null(dominant_class)) {
              oob_resid <- sample_oob_residual(dominant_class)
              if (!is.null(oob_resid)) {
                # Match residual names to coefficient names and add perturbation
                for (vname in names(w_mc_full)) {
                  vname_lower <- tolower(vname)
                  resid_names_lower <- tolower(names(oob_resid))
                  match_idx <- which(resid_names_lower == vname_lower)
                  if (length(match_idx) > 0) {
                    w_mc_full[vname] <- w_mc_full[vname] + oob_resid[match_idx[1]]
                  }
                }
                # Re-enforce sum-to-one and non-negativity constraints
                w_mc_full[w_mc_full < 0] <- 0
                w_sum <- sum(w_mc_full)
                if (w_sum > 1e-9) {
                  w_mc_full <- w_mc_full / w_sum
                }
              }
            }

            w_draws[, b] <- w_mc_full
          }
        }

        if (ncol(w_draws) > 1) {
          mc_coef_sd <- apply(w_draws, 1, stats::sd, na.rm = TRUE)
          mc_coef_q025 <- apply(w_draws, 1, stats::quantile, probs = 0.025, na.rm = TRUE)
          mc_coef_q975 <- apply(w_draws, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
        }
      }

      # Build coefficient dataframe with variant-level detail
      # The veg names are stored as names of chosen_ids (which is a named vector)
      # Values of chosen_ids are variant IDs like "populus_opt_1", but can be NA if veg type not used
      early_zero_applied <- FALSE

      veg_names <- names(chosen_ids)
      if (is.null(veg_names) || length(veg_names) == 0) {
        stop(sprintf("ERROR loc=%s yr=%d: Vegetation names missing from chosen_ids; refusing to continue", loc, yr))
      }
      if (any(is.na(veg_names) | trimws(as.character(veg_names)) == "")) {
        stop(sprintf("ERROR loc=%s yr=%d: Vegetation names contain NA/empty values; refusing to continue", loc, yr))
      }
      coef_df <- data.frame(
        location_id = loc,
        pheno_year = yr,
        lat = lat_val,
        lon = lon_val,
        Veg = veg_names,
        variant_id = as.character(chosen_ids),
        coef = coefs,
        rmse = rmse,
        coef_median = NA_real_,
        coef_025 = if (!is.null(mc_coef_q025)) as.numeric(mc_coef_q025) else NA,
        coef_975 = if (!is.null(mc_coef_q975)) as.numeric(mc_coef_q975) else NA,
        coef_sd = {
          base_sd <- as.numeric(coef_sd_vec)
          if (!is.null(mc_coef_sd)) {
            base_sd <- sqrt(pmax(base_sd, 0)^2 + pmax(as.numeric(mc_coef_sd), 0)^2)
          }
          base_sd
        },
        interval = if (!is.null(mc_coef_q975) && !is.null(mc_coef_q025)) as.numeric(mc_coef_q975) - as.numeric(mc_coef_q025) else NA,
        n_obs = nrow(task_data),
        n_summer_obs = n_summer_valid,
        low_data_flag = low_data_flag,
        inseparable_variant_flag = FALSE,
        inseparable_variant_details = NA_character_,
        stringsAsFactors = FALSE
      )

      # DIAGNOSTIC: Check coef_df immediately after creation

      # --- Mark inseparable variants (drop instead of assign Veg = 'unknown') if detected in similarity tables ---
      # Wrapped in tryCatch to prevent non-critical metadata lookups from crashing the task
      tryCatch({
          sim_tbl <- NULL
          if (exists("INSEPARABLE_VARIANT_INFO") && !is.null(INSEPARABLE_VARIANT_INFO$similarity_table)) {
            sim_tbl <- INSEPARABLE_VARIANT_INFO$similarity_table
          } else if (exists("VARIANT_SIMILARITY_TABLE") && !is.null(VARIANT_SIMILARITY_TABLE)) {
            sim_tbl <- VARIANT_SIMILARITY_TABLE
          }
          if (!is.null(sim_tbl) && nrow(sim_tbl) > 0) {
            cols <- names(sim_tbl)
            variant_id_cols <- intersect(cols, c("variant_id", "variant", "id", "ids", "var_id", "variant1", "variant_a"))
            other_variant_cols <- intersect(cols, c("other_variant_id", "other_variant", "variant2", "other_id", "variant_b"))
            ids_in_tbl <- character(0)
            for (cname in variant_id_cols) ids_in_tbl <- c(ids_in_tbl, as.character(sim_tbl[[cname]]))
            for (cname in other_variant_cols) ids_in_tbl <- c(ids_in_tbl, as.character(sim_tbl[[cname]]))
            ids_in_tbl <- unique(na.omit(ids_in_tbl))
            if (length(ids_in_tbl) > 0) {
              mask <- coef_df$variant_id %in% ids_in_tbl
              mask[is.na(mask)] <- FALSE
              if (any(mask, na.rm = TRUE)) {
                coef_df$inseparable_variant_flag[mask] <- TRUE
                # Build detail strings: collect related vegs and other variant ids
                details <- vapply(coef_df$variant_id[mask], function(vid) {
                  rows <- apply(sim_tbl, 1, function(r) any(vid == as.character(r), na.rm = TRUE))
                  if (any(rows)) {
                    row_tbl <- sim_tbl[rows, , drop = FALSE]
                    other_vegs <- character(0)
                    if ("other_veg" %in% names(row_tbl)) other_vegs <- c(other_vegs, as.character(row_tbl$other_veg))
                    if ("veg" %in% names(row_tbl)) other_vegs <- c(other_vegs, as.character(row_tbl$veg))
                    other_variant_ids <- character(0)
                    for (cname in other_variant_cols) if (cname %in% names(row_tbl)) other_variant_ids <- c(other_variant_ids, as.character(row_tbl[[cname]]))
                    other_variant_ids <- unique(na.omit(other_variant_ids))
                    info_parts <- character(0)
                    if (length(other_vegs) > 0) info_parts <- c(info_parts, paste(unique(na.omit(other_vegs)), collapse = ";"))
                    if (length(other_variant_ids) > 0) info_parts <- c(info_parts, paste(other_variant_ids, collapse = ";"))
                    paste(info_parts, collapse = "|")
                  } else {
                    NA_character_
                  }
                }, FUN.VALUE = "")
                coef_df$inseparable_variant_details[mask] <- details

              }
            }
          }
      }, error = function(e) {
      })

      # Barren is kept as from MESMA - no replacement
      # Aggregate by vegetation type (sum coefficients for same veg type)
      # Must happen AFTER all filtering but BEFORE creating diagnostics/coef_agg usage
      if (nrow(coef_df) > 0) {
        # DIAGNOSTIC: Check coef values before aggregation
        coef_agg <- aggregate(coef ~ Veg, data = coef_df, FUN = sum)
        # DIAGNOSTIC: Check coef values after aggregation
      } else {
        # CRITICAL FIX: When coef_df is empty (e.g., MESMA picked 100% barren),
        # create coef_agg with vegetation types from library, initialized to 0
        # This allows PPI scaling logic to distribute vegetation cover among available types
        veg_types_nonbarren <- setdiff(veg_kept, "barren")
        if (length(veg_types_nonbarren) > 0) {
          coef_agg <- data.frame(
            Veg = veg_types_nonbarren,
            coef = rep(0, length(veg_types_nonbarren)),
            stringsAsFactors = FALSE
          )
        } else {
          coef_agg <- data.frame(Veg = character(0), coef = numeric(0), stringsAsFactors = FALSE)
        }
      }

      veg_coefs_mask <- tolower(coef_agg$Veg) != "barren"
      veg_coefs_mask[is.na(veg_coefs_mask)] <- FALSE
      sum_original_veg_coefs <- sum(coef_agg$coef[veg_coefs_mask], na.rm = TRUE)

      scale_factor <- 1.0

      # compute MESMA-derived barren fraction for diagnostics/output
      mesma_barren <- 1 - sum(coef_agg$coef[veg_coefs_mask], na.rm = TRUE)
      mesma_barren <- pmax(0, pmin(1, mesma_barren))
      barren_fraction <- mesma_barren

      # --- PPI-based barren cap ---
      # Clamp barren fraction to what PPI allows: if peak summer PPI indicates
      # substantial vegetation cover, MESMA should not assign more barren than
      # (1 - ppi_veg_cover).  Freed fraction is redistributed to vegetation coefs.
      if (isTRUE(exists("ENABLE_PPI_BARREN_CAP") && ENABLE_PPI_BARREN_CAP) &&
          "PPI" %in% names(task_data)) {
        ppi_full <- if (exists("PPI_FULL_VEG_COVER") && is.finite(PPI_FULL_VEG_COVER) && PPI_FULL_VEG_COVER > 0) PPI_FULL_VEG_COVER else 0.584568
        summer_months_ppi <- if (exists("SUMMER_DETREND_MONTHS")) SUMMER_DETREND_MONTHS else 6:9
        ppi_vals <- task_data$PPI
        if ("month" %in% names(task_data)) {
          ppi_vals <- task_data$PPI[task_data$month %in% summer_months_ppi]
        } else if ("date" %in% names(task_data)) {
          mo <- as.integer(format(task_data$date, "%m"))
          ppi_vals <- task_data$PPI[mo %in% summer_months_ppi]
        }
        peak_ppi <- suppressWarnings(max(ppi_vals[is.finite(ppi_vals)], na.rm = TRUE))
        if (is.finite(peak_ppi)) {
          ppi_veg_cover <- min(peak_ppi / ppi_full, 1.0)
          ppi_barren_min <- max(0.0, 1.0 - ppi_veg_cover)  # minimum barren allowed by PPI
          if (barren_fraction > ppi_barren_min + 1e-6) {
            freed <- barren_fraction - ppi_barren_min
            barren_fraction <- ppi_barren_min
            # Compute scale factor and apply to both coef_agg and coef_df (via scale_factor)
            if (sum_original_veg_coefs > 1e-9) {
              scale_factor <- (sum_original_veg_coefs + freed) / sum_original_veg_coefs
              coef_agg$coef[veg_coefs_mask] <- coef_agg$coef[veg_coefs_mask] * scale_factor
            } else if (sum(veg_coefs_mask) > 0) {
              # No veg coefs yet — distribute freed equally among veg types in coef_agg only
              n_veg <- sum(veg_coefs_mask)
              coef_agg$coef[veg_coefs_mask] <- freed / n_veg
              # scale_factor stays 1.0; coef_df has no veg rows to scale
            }
            # Update barren row in coef_agg to match
            barren_idx_agg <- which(tolower(coef_agg$Veg) == "barren")
            if (length(barren_idx_agg) > 0) coef_agg$coef[barren_idx_agg] <- barren_fraction
          }
        }
      }

      # ensure barren row exists in coef_agg (value from MESMA)
      barren_idx <- which(tolower(coef_agg$Veg) == "barren")
      if (length(barren_idx) == 0) {
        barren_row <- data.frame(Veg = "barren", coef = barren_fraction)
        coef_agg <- rbind(coef_agg, barren_row)
      }
      # if barren already existed, leave it untouched (already updated by cap block if needed)

      # --- Propagate scaling to coef_df (no equal distribution fallback) ---

      # Wrap scaling in a defensive tryCatch to diagnose errors
      tryCatch({
        # Apply the same scale_factor to coef_df vegetation fractions
        # NO rebuilding coef_df from coef_agg - if MESMA returned empty, keep it empty
        if (nrow(coef_df) > 0 && is.finite(sum_original_veg_coefs) && sum_original_veg_coefs > 1e-9) {
          v_rows <- tolower(coef_df$Veg) != "barren"
          v_rows[is.na(v_rows)] <- FALSE
          coef_df$coef[v_rows] <- coef_df$coef[v_rows] * scale_factor
          if ("coef_sd" %in% names(coef_df)) coef_df$coef_sd[v_rows] <- coef_df$coef_sd[v_rows] * abs(scale_factor)
          if ("coef_025" %in% names(coef_df)) coef_df$coef_025[v_rows] <- coef_df$coef_025[v_rows] * scale_factor
          if ("coef_975" %in% names(coef_df)) coef_df$coef_975[v_rows] <- coef_df$coef_975[v_rows] * scale_factor
        }
        # If MESMA returned 0 vegetation, coef_df veg fractions stay at 0 - no fallback
      }, error = function(e) {
        stop(e)
      })
      
      if (nrow(coef_df) > 0) {
        barren_row_df <- coef_df[1, , drop=FALSE] 
        barren_row_df[] <- NA 
        barren_row_df$location_id <- loc
        barren_row_df$pheno_year <- yr
        barren_row_df$lat <- lat_val
        barren_row_df$lon <- lon_val
        barren_row_df$Veg <- "barren"
        barren_row_df$variant_id <- "barren_ppi"
        barren_row_df$coef <- barren_fraction
        if ("coef_sd" %in% names(barren_row_df)) barren_row_df$coef_sd <- 0
        if ("coef_025" %in% names(barren_row_df)) barren_row_df$coef_025 <- NA
        if ("coef_975" %in% names(barren_row_df)) barren_row_df$coef_975 <- NA
        coef_df <- rbind(coef_df, barren_row_df)
      } else {
        # No vegetation detected, use PPI barren fraction
        if (isTRUE(is.finite(barren_fraction) && barren_fraction > 0)) {
          # Return barren only result
          coef_df <- data.frame(
            location_id = loc,
            pheno_year = yr,
            lat = lat_val,
            lon = lon_val,
            Veg = "barren",
            variant_id = "barren_ppi",
            coef = barren_fraction,
            rmse = if(exists("rmse") && is.numeric(rmse)) rmse else 0,
            coef_025 = NA,
            coef_975 = NA,
            coef_sd = 0,
            interval = NA,
            n_obs = nrow(task_data),
            n_summer_obs = n_summer_valid,
            low_data_flag = low_data_flag,
            inseparable_variant_flag = FALSE,
            inseparable_variant_details = NA_character_,
            stringsAsFactors = FALSE
          )
        } else {
          # No vegetation and no barren - fail
          return(NULL)
        }
      }

      # Create diagnostics dataframe
      diag_df <- data.frame(
        location_id = loc,
        pheno_year = yr,
        stringsAsFactors = FALSE
      )
      # Add each vegetation type's total fraction
      for (v in coef_agg$Veg) {
        if (!is.na(v)) {
          v_mask <- coef_agg$Veg == v
          v_mask[is.na(v_mask)] <- FALSE
          diag_df[[paste0(v, "_fraction")]] <- coef_agg$coef[v_mask][1]
        }
      }
      diag_df$barren_fraction_ppi_based <- barren_fraction # Add barren fraction

      # Reconstruct E_best for returning
      E_best_list <- list()
      for(v in names(chosen_ids)){
        vid <- chosen_ids[[v]]
        # Defensive: ensure we have candidate variants for this veg
        if (is.null(top_variants[[v]]) || length(top_variants[[v]]) == 0) {
          next
        }
        found <- FALSE
        for(variant in top_variants[[v]]){
          # Defensive equality check to avoid NA logicals causing crashes
          if (!is.na(variant$id) && !is.na(vid) && isTRUE(variant$id == vid)) {
            E_best_list[[v]] <- variant$vec
            found <- TRUE
            break
          }
        }
        # Fallback: if chosen variant not found (e.g., NA mismatch), use first candidate and log
        if (!found) {
          E_best_list[[v]] <- top_variants[[v]][[1]]$vec
        }
      }
      E_best <- do.call(cbind, E_best_list)

      return(tryCatch({
        list(
          coef_df = coef_df,
          diagnostics = diag_df,
          uncertainty = NULL,
          residuals = residuals,
          y_hat = y_for_unmixing - residuals,
          y_obs = y_for_unmixing,
          E_best = E_best,
          top_variants = top_variants,
          weights_masked = weights_masked,
          valid_mask = valid_mask
        )
      }, error = function(e) {
        stop(e)
      }))

      # ===== END MESMA UNMIXING =====
  # These variables are already set in the earlier conditional blocks above.
  # No need to reassign them here
  # Single-stage MESMA is used, so mesma_lib and OPTIMIZED_LIBRARY are set accordingly

  }

# New function: Process all years for a single location with multi-year bootstrap
  fit_one_location <- function(location_data) {
    if (is.null(location_data) || nrow(location_data) == 0) {
      return(NULL)
    }

    loc <- as.character(location_data$location_id[1])

    # Get all phenological years for this location
    years <- sort(unique(location_data$pheno_year))
    years <- years[!is.na(years)]
    
    if (length(years) == 0) {
      return(NULL)
    }

    # Process each year individually first (to get point estimates and chosen variants)
    year_results <- list()
    y_vecs_by_year <- list()
    chosen_ids_by_year <- list()
    w_hat_by_year <- list()

    # Wrap per-location processing in tryCatch so a single error does not discard valid year results
    out <- tryCatch({
      for (yr in years) {
        year_data <- location_data[location_data$pheno_year == yr, , drop = FALSE]
        
        res_yr <- tryCatch({
            fit_one_task(year_data)
          }, error = function(e) {
            msg <- sprintf("loc=%s yr=%d: ERROR in fit_one_task: %s", loc, yr, e$message)
            message(paste("[ERROR fit_one_task CRASH]", msg))
            # Capture a compact call stack to pinpoint the failing if() condition
            calls <- tryCatch(sys.calls(), error = function(e2) NULL)
            if (!is.null(calls)) {
              tail_calls <- tail(calls, 12)
              call_str <- paste(vapply(tail_calls, function(x) paste(deparse(x), collapse = ""), character(1)), collapse = " | ")
            }
            NULL
          })

        if (!is.null(res_yr)) {
          tryCatch({
            year_results[[as.character(yr)]] <- res_yr
          }, error = function(e) {
            stop(e)  # Re-throw to be caught by outer handler
          })

          # Store data needed for multi-year bootstrap - wrapped in tryCatch to isolate errors
          tryCatch({
            # For inference/validation data, the temporal fill method is controlled
            # by the global `INTERPOLATE_INFERENCE` flag or
            # `MESMA_PARAMS$interpolate_inference` (can be "linear", "whittaker" or "none").
            # Build from the same raw feature columns (base_indices) as fit_one_task(),
            # then apply representation + z-scoring to match the library feature space.
            n_bins <- TEMPORAL_BUDGET
            base_indices <- if (!is.null(MESMA_PARAMS$base_indices)) MESMA_PARAMS$base_indices else MESMA_PARAMS$indices
            n_base_idx <- length(base_indices)
            l2_normalize <- isTRUE(MESMA_PARAMS$l2_normalize)

            interp_flag <- if (!is.null(MESMA_PARAMS$interpolate_inference)) {
                get_interpolate_method(MESMA_PARAMS$interpolate_inference)
            } else {
                get_interpolate_method(INTERPOLATE_INFERENCE)
            }
            raw_mat_yr <- build_pentad_matrix(year_data, base_indices, interpolate = interp_flag)
            if (!is.null(raw_mat_yr)) {
              vec_raw <- as.numeric(raw_mat_yr)
              y_work <- mesma_apply_representation_vec(vec_raw, n_base_idx, n_bins, l2_normalize)
              y_norm <- mesma_zscore_vec_by_index(y_work, MESMA_PARAMS$indices, MESMA_PARAMS$means, MESMA_PARAMS$sds, n_bins)

              expected_length <- length(MESMA_PARAMS$indices) * n_bins
              if (length(y_norm) != expected_length) {
                # skip; length mismatch
              } else {
                y_vecs_by_year[[length(y_vecs_by_year) + 1]] <- y_norm
              }

              # Get chosen variants and weights from this year's result
              if (!is.null(res_yr$chosen_variants)) {
                chosen_ids_by_year[[yr]] <- res_yr$chosen_variants
              }
              if (!is.null(res_yr$vegetation_proportions)) {
                w_hat_by_year[[yr]] <- res_yr$vegetation_proportions
              }
            }
          }, error = function(e) {
            invisible(NULL)
          })
        } else {
          # fit_one_task returned NULL for this year
        }
      }

      # If we have results for multiple years and ENABLE_UNCERTAINTY, do multi-year bootstrap
      if (length(year_results) == 0) {
        return(NULL)
      }

      # Return results for all years
      return(year_results)

    }, error = function(e) {
      cat(sprintf("[ERROR fit_one_location] loc=%s: %s\n", loc, e$message))
      if (length(year_results) > 0) {
        return(year_results)
      }
      return(NULL)
    })

    out
  }

  # Training is disabled; do not short-circuit validation here. Initialize
  # validation containers and let the validation routines run normally below.
  validation_coefs <- data.frame()
  validation_results_list <- list()

  cat("Starting main processing loop...\n")
  cat("[NOTICE] Validation confusion matrix was computed earlier (or will be computed before inference)\n")
  

if (dir.exists(TEMP_RESULTS_DIR)) {
  cat(sprintf("Removing existing temporary results directory: %s\n", TEMP_RESULTS_DIR))
  unlink(TEMP_RESULTS_DIR, recursive = TRUE)
}
cat(sprintf("Creating temporary results directory: %s\n", TEMP_RESULTS_DIR))
# create directory and parent directories if needed; suppress warnings for existing paths
dir.create(TEMP_RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(TEMP_RESULTS_DIR)) {
  stop(sprintf("Failed to create temporary results directory: %s", TEMP_RESULTS_DIR))
}

# Main processing executes at top-level; helper functions remain in scope.

    # Assign df_tasks for the main processing loop (inference data, not training data)
    if (exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0) {
      df_tasks <- df_tasks_inference
    } else {
      df_tasks <- data.frame()  # Empty if no inference data
    }

    n_train_loc_years <- 0L
    n_infer_loc_years <- 0L
  if (exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0) {
    if (!"pheno_year" %in% names(df_tasks_inference) && "date" %in% names(df_tasks_inference)) df_tasks_inference$pheno_year <- assign_pheno_year(df_tasks_inference$date)
    n_infer_loc_years <- nrow(unique(df_tasks_inference[c("location_id", "pheno_year")]))
  } else if (exists("df_tasks") && !is.null(df_tasks) && nrow(df_tasks) > 0) {
    if (!"pheno_year" %in% names(df_tasks) && "date" %in% names(df_tasks)) df_tasks$pheno_year <- assign_pheno_year(df_tasks$date)
    n_infer_loc_years <- nrow(unique(df_tasks[c("location_id", "pheno_year")]))
  }

  cat(sprintf("Dataset location-years: %d\n", n_infer_loc_years))
  
  cat("Preparing locations for batched processing (multi-year bootstrap)...\n")

  # Group by LOCATION (not location-year)
  target_locations <- if (!is.null(df_tasks) && is.data.frame(df_tasks) &&
                          nrow(df_tasks) > 0 && "location_id" %in% names(df_tasks)) {
    unique(as.character(df_tasks$location_id))
  } else {
    character(0)
  }
  target_locations <- target_locations[!is.na(target_locations) & nzchar(target_locations)]

  n_locs_to_process <- length(target_locations)
  loc_batches <- list()
  n_batches <- 0L

  # Skip processing if no locations to process
  if (n_locs_to_process == 0) {
    cat("[INFO] No locations to process - skipping main processing loop\n")
    results_list <- list()
  } else {
  loc_batches <- split(target_locations, ceiling(seq_along(target_locations) / BATCH_SIZE))
  n_batches <- length(loc_batches)
  # Collect full per-task results for downstream reporting (variant trajectories, diagnostics, uncertainty)
  results_list <- list()

  start_time <- Sys.time()

  for (i in seq_along(loc_batches)) {
    batch_locs <- loc_batches[[i]]
    batch_df <- df_tasks[df_tasks$location_id %in% batch_locs, ]
    batch_location_list <- split(batch_df, batch_df$location_id)
    # Suppress any verbose output from per-location processing
    batch_results <- suppress_output_safely(
      .run_map(batch_location_list, fit_one_location, show_pb = FALSE)
    )

    results_list <- aggregate_batch_results(batch_results, results_list, save_csv_dir = TEMP_RESULTS_DIR)

  cat(sprintf("\r  [Batch %d/%d complete]  ", i, n_batches))
  }
  cat("\n")
  
  end_time <- Sys.time()
  processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(
    "\nProcessing finished in %.2f seconds (%.2f minutes)\n",
    processing_time, processing_time / 60
  ))
  } # End of n_locs_to_process > 0 conditional

  # ==========================================================================
  # VALIDATION PROCESSING (minimal unmix for confusion matrix)
  # ========================================================================== 
  cat("\n=== STARTING VALIDATION PROCESSING ===\n")
  # Re-initialise the parallel backend with a conservative worker count to avoid
  # SOCK worker crashes (Windows OOM / socket desync) during validation unmixing.
  if (exists("cleanup_parallel") && is.function(cleanup_parallel)) {
    tryCatch(cleanup_parallel(), error = function(e) NULL)
  }
  val_workers <- min(PARALLEL_WORKERS, max(1L, parallel::detectCores() - 2L), 15L)
  cleanup_parallel <- setup_parallel_backend(workers = val_workers)
  cat(sprintf("[VALIDATION] Parallel backend reset to %d workers for validation unmixing.\n", val_workers))

  validation_coefs <- data.frame()
  if (exists("df_validation") && !is.null(df_validation) && nrow(df_validation) > 0) {
    # We perform a lightweight unmix of the held-out validation rows so that
    # a proper confusion matrix (predicted vs true) can be computed.  This is
    # the "correct" validation matrix the user expects.  No additional
    # diagnostics (accuracy stats, PPI bias, RMSE, etc.) are calculated.
    
    # Build brief task list matching inference logic
    df_val_proc <- df_validation
    if (!"pheno_year" %in% names(df_val_proc) && "date" %in% names(df_val_proc)) {
      df_val_proc$pheno_year <- assign_pheno_year(df_val_proc$date)
    }
    if (!"task_key" %in% names(df_val_proc)) {
      df_val_proc$task_key <- paste(df_val_proc$location_id, df_val_proc$pheno_year, sep = "_")
    }
    val_task_list <- split(df_val_proc, df_val_proc$task_key)

    cat(sprintf("[VALIDATION] Performing MESMA unmix on %d validation rows (%d locations)\n",
                nrow(df_val_proc), length(unique(df_val_proc$location_id))))

    # Execute in batches to mirror inference; reuse BATCH_SIZE
    val_locs <- unique(df_val_proc$location_id)
    val_loc_batches <- split(val_locs, ceiling(seq_along(val_locs) / BATCH_SIZE))
    validation_results_list <- list()
    for (batch in val_loc_batches) {
      batch_df <- df_val_proc[df_val_proc$location_id %in% batch, ]
      batch_list <- split(batch_df, batch_df$location_id)
      results <- suppress_output_safely(.run_map(batch_list, fit_one_location, show_pb = FALSE))
      validation_results_list <- aggregate_batch_results(results, validation_results_list)
    }
    # combine into data.frame
    validation_coefs <- do.call(rbind, lapply(validation_results_list, function(r) {
      if (!is.null(r$coef_df)) r$coef_df else NULL
    }))
    if (is.null(validation_coefs)) validation_coefs <- data.frame()
    cat(sprintf("[VALIDATION] Obtained %d predicted coefficient rows from %d locations\n",
                if(!is.null(validation_coefs)) nrow(validation_coefs) else 0,
                length(val_locs)))

    # -----------------------------------------------------------------------
    # CORRELATION: avg correctly predicted fraction vs. n filled pentads
    # Use INTERPOLATED-ONLY location-years for accuracy, but count observations
    # BEFORE interpolation on the x-axis.
    # -----------------------------------------------------------------------
    tryCatch({
      if (!is.null(validation_coefs) && nrow(validation_coefs) > 0 &&
          exists("df_validation") && !is.null(df_validation) && nrow(df_validation) > 0) {

        cat("\n[VALIDATION] Computing correlation for interpolated-only samples: correct fraction ~ n_filled_pentads_before\n")

        # Count filled pentads BEFORE interpolation and detect whether temporal
        # interpolation actually added information for that location-year.
        base_idx_corr <- if (!is.null(MESMA_PARAMS$base_indices)) MESMA_PARAMS$base_indices else MESMA_PARAMS$indices
        interp_method_corr <- resolve_interpolation_method(fallback = "linear")

        val_tasks_for_count <- df_validation %>%
          dplyr::mutate(pheno_year = if ("pheno_year" %in% names(.)) pheno_year else assign_pheno_year(as.Date(date))) %>%
          dplyr::group_by(location_id, pheno_year) %>%
          dplyr::group_split()

        val_pentad_counts <- do.call(rbind, lapply(val_tasks_for_count, function(sub) {
          lid <- sub$location_id[1]
          pyr <- sub$pheno_year[1]

          # Always compute pre-interpolation count for x-axis.
          mat_before <- build_pentad_matrix(sub, base_idx_corr, interpolate = "none")
          n_filled_before <- if (!is.null(mat_before)) {
            valid_rows_before <- which(apply(mat_before[-nrow(mat_before), , drop = FALSE], 1, function(r) any(is.finite(r))))
            length(valid_rows_before)
          } else {
            NA_integer_
          }

          # Interpolated-only filter: retain groups where interpolation adds at
          # least one finite pentad value that was missing before.
          mat_after <- build_pentad_matrix(sub, base_idx_corr, interpolate = interp_method_corr)
          was_interpolated <- FALSE
          n_filled_after <- NA_integer_
          if (!is.null(mat_before) && !is.null(mat_after)) {
            before_core <- mat_before[-nrow(mat_before), , drop = FALSE]
            after_core  <- mat_after[-nrow(mat_after), , drop = FALSE]
            n_filled_after <- length(which(apply(after_core, 1, function(r) any(is.finite(r)))))
            newly_filled <- (!is.finite(before_core)) & is.finite(after_core)
            was_interpolated <- any(newly_filled, na.rm = TRUE)
          }

          data.frame(location_id = lid,
                     pheno_year = pyr,
                     n_filled_pentads_before = n_filled_before,
                     n_filled_pentads_after = n_filled_after,
                     was_interpolated = was_interpolated,
                     stringsAsFactors = FALSE)
        }))

        # true_veg per location-year comes from df_validation (not from coef_df)
        val_true_veg <- df_validation %>%
          dplyr::group_by(location_id, pheno_year) %>%
          dplyr::summarise(true_veg = {
            tv <- tolower(unique(as.character(Veg)))
            tv <- tv[!is.na(tv) & tv != ""]
            if (length(tv) == 1) tv else NA_character_
          }, .groups = "drop")

        # Per location-year: correct fraction = coef of true_veg / sum(all veg coefs)
        correct_frac_df2 <- validation_coefs %>%
          dplyr::filter(tolower(as.character(Veg)) != "barren") %>%
          dplyr::left_join(val_true_veg, by = c("location_id", "pheno_year")) %>%
          dplyr::group_by(location_id, pheno_year, true_veg) %>%
          dplyr::summarise(
            veg_total = sum(coef, na.rm = TRUE),
            correct_coef = sum(coef[tolower(as.character(Veg)) == true_veg[1]], na.rm = TRUE),
            .groups = "drop"
          ) %>%
          dplyr::mutate(
            correct_frac = dplyr::if_else(veg_total > 1e-10, correct_coef / veg_total, NA_real_)
          ) %>%
          dplyr::filter(!is.na(true_veg) & !is.na(correct_frac))

        # Merge with pentad counts and keep interpolated-only location-years.
        val_corr_df <- correct_frac_df2 %>%
          dplyr::inner_join(val_pentad_counts, by = c("location_id", "pheno_year"))

        val_corr_df <- val_corr_df %>%
          dplyr::filter(was_interpolated, true_veg != "barren")

        if (nrow(val_corr_df) >= 3) {
          r_val <- cor(val_corr_df$n_filled_pentads_before, val_corr_df$correct_frac,
                       use = "complete.obs", method = "pearson")
          cat(sprintf("[VALIDATION] Pearson r (interpolated-only; correct_frac ~ n_filled_pentads_before): %.3f  (n=%d)\n",
                      r_val, nrow(val_corr_df)))

          p_val <- tryCatch(
            cor.test(val_corr_df$n_filled_pentads_before, val_corr_df$correct_frac,
                     method = "pearson")$p.value,
            error = function(e) NA_real_
          )
          if (!is.na(p_val)) cat(sprintf("[VALIDATION] p-value: %.4g\n", p_val))

          # Scatter plot
          p_corr <- ggplot2::ggplot(val_corr_df,
                                    ggplot2::aes(x = n_filled_pentads_before, y = correct_frac,
                                                 colour = true_veg)) +
            ggplot2::geom_jitter(alpha = 0.45, size = 1.8,
                                 width = 0.25, height = 0.015) +
            ggplot2::geom_smooth(method = "lm", se = TRUE, colour = "black",
                                 linewidth = 0.8, linetype = "dashed") +
            ggplot2::annotate("text", x = Inf, y = -Inf,
                              label = sprintf("r = %.3f\np = %.3g\nn = %d",
                                             r_val,
                                             if (!is.na(p_val)) p_val else NaN,
                                             nrow(val_corr_df)),
                              hjust = 1.05, vjust = -0.3, size = 3.5) +
            ggplot2::labs(
              title = sprintf("Validation (interpolated-only): avg correctly predicted fraction vs. filled pentads before interpolation [%s]", interp_method_corr),
              x = "N filled pentads (before interpolation)",
              y = "Correctly predicted fraction (normalised)",
              colour = "True veg"
            ) +
            ggplot2::theme_bw()

          corr_plot_path <- file.path(OUTPUT_DIR, "validation_correct_frac_vs_preinterp_pentads_interpolated_only.png")
          tryCatch({
            ggplot2::ggsave(corr_plot_path, p_corr, width = 7, height = 5, dpi = 150)
            cat(sprintf("[VALIDATION] Saved correlation plot: %s\n", corr_plot_path))
          }, error = function(e) cat(sprintf("[VALIDATION] Could not save plot: %s\n", e$message)))
        } else {
          cat(sprintf("[VALIDATION] Too few interpolated rows (%d) to compute correlation\n", nrow(val_corr_df)))
        }
      }
    }, error = function(e) {
      cat(sprintf("[VALIDATION] Correlation analysis failed: %s\n", e$message))
    })
    # -----------------------------------------------------------------------

  } else {
    cat("[VALIDATION] No df_validation available; skipping unmixing\n")
  }
  
  # ==========================================================================
  # INFERENCE PROCESSING (separate from validation)
  # ==========================================================================
  cat("\n=== STARTING INFERENCE PROCESSING ===\n")
  
  # Generate variant similarity heatmap before any inference processing
  cat("[INFERENCE] Generating variant similarity heatmap...\n")
  ensure_variant_similarity_heatmap(force = TRUE)
  
  # Load inference data from INFERENCE_CSV if not already loaded
  cat("[INFERENCE] Loading inference data from separate CSV file...\n")
  load_and_prepare_inference_data()

  if (exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0) {
    # Prepare inference task list - use df_tasks_inference loaded from INFERENCE_CSV
    df_tasks_inference_proc <- df_tasks_inference
    
    if ("Veg" %in% names(df_tasks_inference_proc)) {
      df_tasks_inference_proc$Veg <- tolower(as.character(df_tasks_inference_proc$Veg))
    }
    
    if ("date" %in% names(df_tasks_inference_proc)) {
      df_tasks_inference_proc$date <- as.Date(df_tasks_inference_proc$date)
      if (!"pheno_year" %in% names(df_tasks_inference_proc)) {
        df_tasks_inference_proc$pheno_year <- assign_pheno_year(df_tasks_inference_proc$date)
      }
      df_tasks_inference_proc$doy <- pheno_doy(df_tasks_inference_proc$date)
      df_tasks_inference_proc$doy[df_tasks_inference_proc$doy < 1 | df_tasks_inference_proc$doy > 366] <- NA_integer_
      if (any(is.na(df_tasks_inference_proc$doy))) stop("[DOY] Missing/invalid DOY values in inference data")
    }
    
    # Get unique inference locations
    inference_locations <- df_tasks_inference_proc |>
      dplyr::filter(!is.na(.data$location_id) & trimws(.data$location_id) != "") |>
      dplyr::distinct(.data$location_id)
    inference_locations$location_id <- trimws(as.character(inference_locations$location_id))
    
    # Do not re-sample here: MAX_INFERENCE_LOCATIONS is already applied once per inference file during loading.

    # Build inference_loc_years with defensive checks
    required_cols <- c("location_id", "pheno_year")
    if (!all(required_cols %in% names(df_tasks_inference_proc))) {
      missing_cols <- setdiff(required_cols, names(df_tasks_inference_proc))
      stop(sprintf("[INFERENCE ERROR] Missing required columns in df_tasks_inference_proc: %s", paste(missing_cols, collapse = ", ")))
    }

    tmp_inf <- df_tasks_inference_proc
    tmp_inf$location_id <- as.character(tmp_inf$location_id)
    tmp_inf$pheno_year <- as.integer(tmp_inf$pheno_year)
    tmp_inf$location_id <- trimws(tmp_inf$location_id)

    valid_mask <- !is.na(tmp_inf$location_id) & tmp_inf$location_id != "" & !is.na(tmp_inf$pheno_year) & tmp_inf$pheno_year > 0
    if (any(valid_mask)) {
      inference_loc_years <- unique(tmp_inf[valid_mask, c("location_id", "pheno_year")])
      inference_loc_years <- data.frame(location_id = inference_loc_years[,1], pheno_year = inference_loc_years[,2], stringsAsFactors = FALSE)
    } else {
      inference_loc_years <- data.frame(location_id = character(0), pheno_year = integer(0), stringsAsFactors = FALSE)
    }

    cat(sprintf("[INFERENCE] Processing %d locations (%d location-year pairs)\n",
                nrow(inference_locations), nrow(inference_loc_years)))
    
    # Process inference in batches
    inference_location_list <- as.character(inference_locations$location_id)
    inference_loc_batches <- split(inference_location_list, 
                                    ceiling(seq_along(inference_location_list) / BATCH_SIZE))
    
    cat(sprintf("[INFERENCE] Processing %d locations in %d batches (approx %d locations/batch)...\n",
                length(inference_location_list), length(inference_loc_batches), BATCH_SIZE))

    n_inference_batches <- length(inference_loc_batches)
    progress_targets <- seq.int(10L, 100L, by = 10L)
    next_progress_idx <- 1L
    
    inference_results_list <- list()
    inference_start_time <- Sys.time()
    
    for (i in seq_along(inference_loc_batches)) {
      batch_locs <- inference_loc_batches[[i]]
      batch_df <- df_tasks_inference_proc[df_tasks_inference_proc$location_id %in% batch_locs, ]
      batch_location_list <- split(batch_df, batch_df$location_id)

      # Always emit a lightweight batch heartbeat so long batches don't look hung
      batch_t0 <- Sys.time()
      
      batch_results <- .run_map(batch_location_list, fit_one_location, show_pb = FALSE)

      # Store results using shared helper
      n_before <- length(inference_results_list)
      inference_results_list <- aggregate_batch_results(batch_results, inference_results_list)
      n_stored <- length(inference_results_list) - n_before
      
      # Emit progress at each 10% completion mark.
      if (n_inference_batches > 0) {
        pct_done <- as.integer(floor(100 * i / n_inference_batches))
        while (next_progress_idx <= length(progress_targets) && pct_done >= progress_targets[[next_progress_idx]]) {
          cat(sprintf("[INFERENCE] %d%% complete (%d/%d batches)\n",
                      progress_targets[[next_progress_idx]], i, n_inference_batches))
          flush.console()
          next_progress_idx <- next_progress_idx + 1L
        }
      }
    }

    # Ensure 100% is always reported once if not already printed
    if (n_inference_batches > 0 && next_progress_idx <= length(progress_targets)) {
      cat(sprintf("[INFERENCE] 100%% complete (%d/%d batches)\n", n_inference_batches, n_inference_batches))
    }

    cat("\n")
    
    inference_end_time <- Sys.time()
    inference_processing_time <- as.numeric(difftime(inference_end_time, inference_start_time, units = "secs"))
    cat(sprintf("[INFERENCE] Processing finished in %.2f seconds (%.2f minutes)\n",
                inference_processing_time, inference_processing_time / 60))
    
    # Combine inference results
    inference_coefs <- do.call(rbind, lapply(inference_results_list, function(r) {
      if (!is.null(r$coef_df)) r$coef_df else NULL
    }))

    # Remove rows with NA or empty Veg values
    if (!is.null(inference_coefs) && nrow(inference_coefs) > 0) {
      n_before <- nrow(inference_coefs)
      inference_coefs <- filter_valid_vegetation(inference_coefs, exclude_barren = FALSE)
      na_veg_count <- n_before - nrow(inference_coefs)
      if (na_veg_count > 0) {
        cat(sprintf("[INFERENCE] Removing %d rows with NA/empty Veg values\n", na_veg_count))
      }
    }

    cat(sprintf("[INFERENCE] Processed %d inference coefficient rows\n",
                if(!is.null(inference_coefs)) nrow(inference_coefs) else 0))

    # Calculate exact point-wise PPI normals for CSV output to ensure absolute cover fraction trends are preserved downstream
    inference_coefs_csv <- inference_coefs
    if (!is.null(inference_coefs_csv) && nrow(inference_coefs_csv) > 0 && exists("df_tasks_inference_proc") && "ppi_norm" %in% names(df_tasks_inference_proc)) {
      cat("[INFERENCE] Adjusting pixel-level absolute PPI-normalized coefficients (saving as proper 'coef') for CSV output\n")
      tryCatch({
        summer_ppi <- df_tasks_inference_proc |> 
          dplyr::filter(!is.na(pheno_year)) |>
          dplyr::group_by(location_id, pheno_year) |>
          dplyr::summarise(ppi_scaler = median(ppi_norm, na.rm = TRUE), .groups = "drop") |>
          dplyr::mutate(ppi_scaler = pmin(pmax(ifelse(is.na(ppi_scaler), 0, ppi_scaler), 0), 1))
        
        inference_coefs_csv <- inference_coefs_csv |>
          dplyr::left_join(summer_ppi, by = c("location_id", "pheno_year")) |>
          dplyr::group_by(location_id, pheno_year) |>
          dplyr::mutate(
            sum_veg_coef = sum(coef[tolower(Veg) != "barren"], na.rm = TRUE),
            rel_coef = ifelse(sum_veg_coef > 1e-6, coef / sum_veg_coef, coef),
            total_rel = sum(rel_coef[tolower(Veg) != "barren"], na.rm = TRUE),
            rel_coef = ifelse(total_rel > 1e-6 & abs(total_rel - 1) > 1e-9 & tolower(Veg) != "barren", 
                              rel_coef / total_rel, rel_coef),
            # Overwrite the actual main coefficient directly with the absolute normalized value
            coef = as.numeric(ifelse(tolower(Veg) == "barren", 
                                   1.0 - ifelse(is.na(ppi_scaler), 0, ppi_scaler), 
                                   rel_coef * ifelse(is.na(ppi_scaler), 0, ppi_scaler)))
          ) |>
          dplyr::select(-sum_veg_coef, -rel_coef, -total_rel, -ppi_scaler) |>
          dplyr::ungroup()
      }, error = function(e) {
        cat(sprintf("[INFERENCE ERROR] Could not substitute 'coef' column: %s\n", e$message))
      })
    }

    # Always export canonical long-format inference fractions for downstream scripts
    tryCatch({
      save_inference_results_csv(
        inference_coefs_csv,
        out_csv = "inference_results/inference_results.csv",
        source_csv = if (exists("INFERENCE_CSV")) INFERENCE_CSV else NA_character_
      )
    }, error = function(e) {
      cat(sprintf("[INFERENCE] Failed to write canonical inference CSV: %s\n", e$message))
    })

    # --- INFERENCE: PPI-based barren estimation (always enabled) ---
    if (!is.null(inference_coefs) && nrow(inference_coefs) > 0 && exists("df_tasks_inference_proc") && ("PPI" %in% names(df_tasks_inference_proc) || "PPI_raw" %in% names(df_tasks_inference_proc))) {
      cat("[INFERENCE] Running PPI-based barren estimation (inference mode)\n")
      ppi_inf_full <- location_bootstrap_ppi(inference_coefs, df_tasks_inference_proc, B = BOOTSTRAP_B, seed = get_mesma_seed(123))
      if (is.null(ppi_inf_full) || nrow(ppi_inf_full) == 0) {
        cat("[INFERENCE] PPI inference aggregation returned no results (no matching loc-year PPI values).\n")
      } else {
        plot_inference_method_results(ppi_inf_full, "PPI", "ppi",
                                      use_excluded_years_shade = TRUE,
                                      include_species_plots = TRUE)
      }
    } else {
      cat("[INFERENCE] PPI data not available for inference; skipping PPI barren estimation.\n")
    }

    # --- INFERENCE: only PPI aggregation runs ---


    # --- INFERENCE: Aggregate bootstrap ---
    if (!is.null(inference_coefs) && nrow(inference_coefs) > 0) {
      cat("[INFERENCE] Running aggregate bootstrap\n")
      aggregate_inf_full <- tryCatch({
        location_bootstrap_aggregate(inference_coefs, B = BOOTSTRAP_B, seed = get_mesma_seed(123))
      }, error = function(e) { cat(sprintf("[AGGREGATE BOOTSTRAP] failed: %s\n", e$message)); NULL })
      if (!is.null(aggregate_inf_full) && nrow(aggregate_inf_full) > 0) {
        plot_inference_method_results(aggregate_inf_full, "Aggregate", "aggregate",
                                      use_excluded_years_shade = FALSE,
                                      include_species_plots = TRUE)
      } else {
        cat("[INFERENCE] Aggregate bootstrap returned no results.\n")
      }
    } else {
      cat("[INFERENCE] No inference coefficients available; skipping aggregate bootstrap.\n")
    }

  } else {
    cat("[INFERENCE] No inference data to process\n")
    inference_coefs <- NULL
  }
  
  # Combine validation and inference results for downstream processing
  all_coefs <- rbind(validation_coefs, inference_coefs)
  cat(sprintf("[COMBINED] Total coefficient rows: %d (validation: %d, inference: %d)\n",
              if(!is.null(all_coefs)) nrow(all_coefs) else 0,
              if(!is.null(validation_coefs)) nrow(validation_coefs) else 0,
              if(!is.null(inference_coefs)) nrow(inference_coefs) else 0))

  # Defensive: ensure certain columns are numeric (avoid coercion warnings later)
  if (!is.null(all_coefs) && nrow(all_coefs) > 0) {
    num_cols <- c("coef", "coef_sd", "coef_025", "coef_975", "pheno_year")
    for (nc in num_cols) {
      if (nc %in% names(all_coefs)) {
        if (!is.numeric(all_coefs[[nc]])) {
          before_na <- sum(is.na(all_coefs[[nc]]))
          all_coefs[[nc]] <- suppressWarnings(as.numeric(as.character(all_coefs[[nc]])))
          after_na <- sum(is.na(all_coefs[[nc]]))
          if (after_na > before_na) {
            warning(sprintf("Coerced column '%s' to numeric; NAs increased by %d", nc, after_na - before_na))
          }
        }
      }
    }
  }

  # Training results CSV saving is disabled by default (disabled by config)
  if (!is.null(all_coefs) && nrow(all_coefs) > 0) {
    cat("[INFO] Training results CSV saving skipped (training outputs removed by config)\n")
  }

  # Calculate total number of year results processed
  n_year_results <- 0
  if (!is.null(all_coefs) && nrow(all_coefs) > 0) {
    n_year_results <- length(unique(paste(all_coefs$location_id, all_coefs$pheno_year, sep = "_")))
  }

  if (n_year_results > 0 && exists("processing_time")) {
    cat(sprintf("Average time per year-result: %.2f seconds\n", processing_time / n_year_results))
  } else {
    cat("Average time per year-result: N/A (0 results)\n")
  }

  # Ensure required library/templates exist for visualization
  ensure_library_and_templates()

  df <- do.call(rbind, lapply(results, as.data.frame))
  list(best_params=best_params,
       best_score=best_score,
       best_W=best_W,
       all_results=df)
# End of script