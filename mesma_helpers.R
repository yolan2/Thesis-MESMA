
apply_oli_etm_bias_correction <- function(df,
                                          dataset_label = "training",
                                          satellite_col = "satellite",
                                          bias_csv = NULL,
                                          raw_bands = c("blue", "green", "red", "nir", "swir1", "swir2"),
                                          enable = NULL,
                                          log_prefix = "[BIAS CORR]") {
  if (is.null(bias_csv)) {
    # Prefer the one in OUT_DIR (where satellite_bias_check.R actually writes if mesma_config.R is sourced)
    if (exists("OUT_DIR", inherits = TRUE)) {
      bias_csv <- file.path(get("OUT_DIR", inherits = TRUE), "bias_stats_features.csv")
    }
    # Fallback to local folder if OUT_DIR version doesn't exist
    if (is.null(bias_csv) || !file.exists(bias_csv)) {
      bias_csv <- "satellite_bias_check/bias_stats_features.csv"
    }
  }

  if (is.null(enable)) {
    enable <- if (exists("ENABLE_BAND_BIAS_CORRECTION", inherits = TRUE)) {
      isTRUE(get("ENABLE_BAND_BIAS_CORRECTION", inherits = TRUE))
    } else {
      TRUE
    }
  }

  if (!isTRUE(enable)) {
    cat(sprintf("%s Band-level bias correction disabled (ENABLE_BAND_BIAS_CORRECTION=FALSE).\n", log_prefix))
    return(df)
  }

  if (!(satellite_col %in% names(df))) {
    cat(sprintf("%s No '%s' column in %s data; bias correction skipped.\n",
                log_prefix, satellite_col, dataset_label))
    return(df)
  }

  if (!file.exists(bias_csv)) {
    stop(sprintf("%s %s not found — run satellite_bias_check.R to generate it.", log_prefix, bias_csv))
  }

  bias_tbl <- tryCatch(data.table::fread(bias_csv), error = function(e) NULL)
  if (is.null(bias_tbl) || !all(c("band", "ols_slope", "ols_intercept") %in% names(bias_tbl))) {
    stop(sprintf("%s Bias CSV is missing required columns (band / ols_slope / ols_intercept). Re-run satellite_bias_check.R to regenerate it.",
                 log_prefix))
  }
  if (!("ols_significant" %in% names(bias_tbl))) {
    bias_tbl[, ols_significant := TRUE]
    cat(sprintf("%s 'ols_significant' not found in bias CSV; defaulting to TRUE for all bands.\n", log_prefix))
  }

  bias_tbl[, band_lower := tolower(trimws(as.character(band)))]
  bias_bands <- bias_tbl[band_lower %in% raw_bands]
  oli_mask <- tolower(trimws(as.character(df[[satellite_col]]))) == "landsat_89"

  cat(sprintf("%s %d/%d %s rows are LANDSAT_89; applying affine (slope+intercept) correction\n",
              log_prefix, sum(oli_mask, na.rm = TRUE), nrow(df), dataset_label))

  for (i in seq_len(nrow(bias_bands))) {
    b <- bias_bands$band_lower[i]
    if (!(b %in% names(df))) next

    if (!isTRUE(bias_bands$ols_significant[i])) {
      p_val <- if ("ols_p_value" %in% names(bias_bands)) suppressWarnings(as.numeric(bias_bands$ols_p_value[i])) else NA_real_
      p_msg <- if (is.finite(p_val)) sprintf("p=%.6f", p_val) else "p=NA"
      cat(sprintf("%s   %-6s  skipped (no significant relation, %s)\n",
                  log_prefix, b, p_msg))
      next
    }

    slope <- suppressWarnings(as.numeric(bias_bands$ols_slope[i]))
    intcp <- suppressWarnings(as.numeric(bias_bands$ols_intercept[i]))
    if (!is.finite(slope) || !is.finite(intcp)) {
      stop(sprintf("%s Non-finite affine terms for band '%s'; check satellite_bias_check.R output.", log_prefix, b))
    }

    df[[b]] <- as.numeric(df[[b]])
    df[[b]][oli_mask] <- slope * df[[b]][oli_mask] + intcp
    cat(sprintf("%s   %-6s  slope=%.6f  intercept=%+.6f  (%d OLI rows)\n",
                log_prefix, b, slope, intcp, sum(oli_mask & is.finite(df[[b]]))))
  }

  cat(sprintf("%s Band correction complete; indices will be recomputed from corrected bands.\n", log_prefix))
  df
}

compute_soil_line_slope <- function(input_df, min_samples = NULL, assign_global_dvi = TRUE) {
  if (is.null(min_samples)) {
    if (exists("MIN_ENDMEMBER_SAMPLES", inherits = TRUE)) {
      min_samples <- get("MIN_ENDMEMBER_SAMPLES", inherits = TRUE)
    } else {
      min_samples <- 5L
    }
  }

  if (is.null(input_df) || nrow(input_df) == 0) {
    stop("[SOIL LINE] No input data provided to compute_soil_line_slope")
  }
  if (!all(c("nir", "red", "Veg") %in% names(input_df))) {
    stop("[SOIL LINE] Cannot compute soil line slope: required columns 'nir', 'red', and 'Veg' are missing")
  }

  veg_norm <- tolower(trimws(as.character(input_df$Veg)))
  bare_soil_df <- input_df[veg_norm == "barren" & is.finite(input_df$nir) & is.finite(input_df$red), , drop = FALSE]

  if (nrow(bare_soil_df) <= as.integer(min_samples)) {
    stop(sprintf(
      "[SOIL LINE] Not enough bare soil pixels to estimate SOIL_LINE_SLOPE (need > %d, have %d)",
      as.integer(min_samples),
      as.integer(nrow(bare_soil_df))
    ))
  }

  soil_line_model <- tryCatch(lm(nir ~ red, data = bare_soil_df), error = function(e) e)
  if (inherits(soil_line_model, "error")) {
    stop("[SOIL LINE] Linear fit failed; cannot estimate SOIL_LINE_SLOPE")
  }

  slope <- as.numeric(coef(soil_line_model)[2])
  if (!is.finite(slope)) {
    stop("[SOIL LINE] Estimated slope is non-finite")
  }

  assign("SOIL_LINE_SLOPE", slope, envir = globalenv())
  cat(sprintf("[SOIL LINE] Calculated SOIL_LINE_SLOPE=%.4f from %d bare soil pixels\n", slope, nrow(bare_soil_df)))

  if (isTRUE(assign_global_dvi)) {
    dvi_soil_calc <- mean(bare_soil_df$nir - bare_soil_df$red, na.rm = TRUE)
    if (is.finite(dvi_soil_calc)) {
      cat(sprintf("[SOIL LINE] Computed training DVI soil baseline (local only): dvi_soil = %.6f\n", dvi_soil_calc))
    }
  }

  invisible(slope)
}

# --- Spatial autocorrelation / block-bootstrap utilities ---
# Shared helper functions for estimating spatial dependence and performing
# spatially-aware "moving block" resampling of locations.  Used by both the
# MESMA fitting script and the january_averages workflow so that behaviour is
# identical regardless of the entry point.

# compute pairwise Haversine distance matrix (km) for an Nx2 lon/lat input.
compute_haversine_distance_matrix <- function(coords) {
  rad <- pi / 180
  lat <- coords[,2] * rad
  lon <- coords[,1] * rad
  dlat <- outer(lat, lat, "-")
  dlon <- outer(lon, lon, "-")
  a <- sin(dlat/2)^2 + outer(cos(lat), cos(lat)) * sin(dlon/2)^2
  2 * 6371 * asin(pmin(1, sqrt(a)))
}

# Shared utility: fit exponential variogram via NLS with fallback
# (previously duplicated in multiple scripts).  Returns estimated range (r)
# or NULL if fit failed; used by estimate_autocorrelation_range.
fit_exponential_variogram <- function(bin_mid, bin_gamma, total_var, dists) {
  valid_bins <- !is.na(bin_gamma)
  if (sum(valid_bins) < 2) return(NULL)
  tryCatch({
    nls_fit <- nls(g ~ s * (1 - exp(-d / r)),
                   data = data.frame(d = bin_mid[valid_bins], g = bin_gamma[valid_bins]),
                   start = list(s = total_var, r = median(dists)),
                   lower = list(s = total_var * 0.1, r = max(dists) * 0.01),
                   upper = list(s = total_var * 3, r = max(dists) * 2),
                   algorithm = "port",
                   control = list(maxiter = 50, warnOnly = TRUE))
    as.numeric(coef(nls_fit)["r"])
  }, error = function(e) {
    thresh_idx <- which(bin_gamma[valid_bins] >= 0.5 * total_var)
    if (length(thresh_idx) > 0) bin_mid[valid_bins][thresh_idx[1]]
    else max(dists)
  })
}

# Collect unique coordinates for a list of locations.  Prefers the provided
# df_tasks/all_coefs frame but will fall back to the global gpts_map if
# necessary.  Coordinates are deduplicated (keeping finite values first).
collect_location_coords <- function(locations, df_tasks = NULL, all_coefs = NULL) {
  locations <- trimws(as.character(locations))
  locations <- locations[!is.na(locations) & locations != ""]
  if (length(locations) == 0) return(data.frame(location_id = character(0), lat = numeric(0), lon = numeric(0)))

  cand <- NULL
  if (!is.null(df_tasks) && all(c("location_id","lat","lon") %in% names(df_tasks))) {
    cand <- df_tasks[, c("location_id","lat","lon"), drop = FALSE]
  } else if (!is.null(all_coefs) && all(c("location_id","lat","lon") %in% names(all_coefs))) {
    cand <- all_coefs[, c("location_id","lat","lon"), drop = FALSE]
  }

  if (!is.null(cand)) {
    cand$location_id <- trimws(as.character(cand$location_id))
    cand$lat <- suppressWarnings(as.numeric(cand$lat))
    cand$lon <- suppressWarnings(as.numeric(cand$lon))
    cand <- cand[!is.na(cand$location_id) & cand$location_id != "", , drop = FALSE]
    cand <- deduplicate_coords(cand)
  }

  # fall back to global gpts_map if coords still inadequate
  if ((is.null(cand) || nrow(cand) == 0 || all(is.na(cand$lat)) || all(is.na(cand$lon))) &&
      exists("gpts_map", envir = globalenv())) {
    gm <- get("gpts_map", envir = globalenv())
    if (!is.null(gm) && all(c("location_id","lat","lon") %in% names(gm))) {
      cand <- gm[, c("location_id","lat","lon"), drop = FALSE]
      cand$location_id <- trimws(as.character(cand$location_id))
      cand$lat <- suppressWarnings(as.numeric(cand$lat))
      cand$lon <- suppressWarnings(as.numeric(cand$lon))
      cand <- cand[!is.na(cand$location_id) & cand$location_id != "", , drop = FALSE]
      cand <- deduplicate_coords(cand)
    }
  }

  if (is.null(cand) || nrow(cand) == 0) {
    return(data.frame(location_id = locations, lat = NA_real_, lon = NA_real_))
  }
  cand <- cand[match(locations, cand$location_id), , drop = FALSE]
  if (nrow(cand) == 0) {
    return(data.frame(location_id = locations, lat = NA_real_, lon = NA_real_))
  }
  cand
}

# Estimate spatial autocorrelation range using an exponential variogram.
estimate_autocorrelation_range <- function(coords_df, values, fallback_km = 30.0) {
  if (is.null(coords_df) || nrow(coords_df) == 0) return(fallback_km)
  if (!all(c("lat", "lon") %in% names(coords_df))) return(fallback_km)

  lat <- suppressWarnings(as.numeric(coords_df$lat))
  lon <- suppressWarnings(as.numeric(coords_df$lon))
  values <- suppressWarnings(as.numeric(values))

  valid <- which(is.finite(lat) & is.finite(lon) & is.finite(values))
  if (length(valid) < 5) return(fallback_km)  # need enough pairs

  coords <- cbind(lon[valid], lat[valid])
  vals <- values[valid]
  if (var(vals, na.rm = TRUE) == 0) return(fallback_km)

  dist_mat <- compute_haversine_distance_matrix(coords)
  dists <- dist_mat[upper.tri(dist_mat)]
  coef_diffs_sq <- outer(vals, vals, function(a, b) (a - b)^2)
  gamma_vals <- coef_diffs_sq[upper.tri(coef_diffs_sq)] / 2
  if (length(dists) == 0) return(fallback_km)

  total_var <- var(vals, na.rm = TRUE)
  n_bins <- min(10, max(3, length(dists) %/% 5))
  bin_breaks <- unique(quantile(dists, probs = seq(0, 1, length.out = n_bins + 1)))
  if (length(bin_breaks) < 3) return(fallback_km)

  bin_mid <- (bin_breaks[-length(bin_breaks)] + bin_breaks[-1]) / 2
  bin_gamma <- numeric(length(bin_mid))
  for (bb in seq_along(bin_mid)) {
    in_bin <- dists >= bin_breaks[bb] & dists < bin_breaks[bb + 1]
    bin_gamma[bb] <- if (sum(in_bin) > 0) median(gamma_vals[in_bin]) else NA
  }
  valid_bins <- !is.na(bin_gamma)
  if (sum(valid_bins) < 2) return(fallback_km)

  range_est <- fit_exponential_variogram(bin_mid, bin_gamma, total_var, dists)
  if (is.null(range_est)) return(fallback_km)
  range_est <- as.numeric(range_est)
  if (!is.finite(range_est) || range_est <= 0) return(fallback_km)

  range_est <- max(1, min(range_est, 500))
  cat(sprintf("[SPATIAL] Estimated autocorrelation range: %.1f km (used as block size)\n", range_est))
  range_est
}

# Moran's I test for significant spatial autocorrelation; returns TRUE if
# block bootstrap should be used.
test_spatial_autocorrelation <- function(coords_df, values, alpha = 0.05, n_perm = 199) {
  if (is.null(coords_df) || nrow(coords_df) == 0) return(FALSE)
  if (!all(c("lat", "lon") %in% names(coords_df))) return(FALSE)
  lat  <- suppressWarnings(as.numeric(coords_df$lat))
  lon  <- suppressWarnings(as.numeric(coords_df$lon))
  vals <- suppressWarnings(as.numeric(values))
  ok   <- is.finite(lat) & is.finite(lon) & is.finite(vals)
  if (sum(ok) < 5) return(FALSE)
  lat  <- lat[ok]; lon <- lon[ok]; vals <- vals[ok]
  if (var(vals) == 0) return(FALSE)

  dist_mat <- as.matrix(compute_haversine_distance_matrix(cbind(lon, lat)))
  if (!is.matrix(dist_mat) || nrow(dist_mat) < 2 || ncol(dist_mat) < 2) return(FALSE)
  w <- 1 / dist_mat
  diag(w) <- 0
  w[!is.finite(w)] <- 0
  rs <- rowSums(w, na.rm = TRUE); rs[rs == 0] <- 1
  w  <- w / rs

  moran_stat <- function(x) {
    n  <- length(x); xc <- x - mean(x)
    s0 <- sum(w, na.rm = TRUE); if (s0 == 0) return(0)
    (n / s0) * sum(w * outer(xc, xc), na.rm = TRUE) / sum(xc^2)
  }

  obs_i  <- moran_stat(vals)
  perm_i <- replicate(n_perm, moran_stat(sample(vals)))
  p_val  <- (sum(perm_i >= obs_i) + 1L) / (n_perm + 1L)
  cat(sprintf("[SPATIAL] Moran's I = %.4f, p = %.3f (%d perms) -> block bootstrap: %s\n",
              obs_i, p_val, n_perm, if (p_val < alpha) "YES" else "NO"))
  p_val < alpha
}

# Return estimated block size (km) if spatial autocorrelation is significant,
# else zero so the bootstrap falls back to i.i.d. sampling.
block_km_if_significant <- function(coords_df, values, alpha = 0.05, n_perm = 199,
                                    fallback_km = BOOTSTRAP_BLOCK_KM) {
  sig <- tryCatch(
    test_spatial_autocorrelation(coords_df, values, alpha = alpha, n_perm = n_perm),
    error = function(e) { warning("[SPATIAL] Moran test error: ", e$message); FALSE }
  )
  if (!isTRUE(sig)) {
    cat("[SPATIAL] No significant spatial autocorrelation -> using regular (i.i.d.) bootstrap\n")
    return(0)
  }
  estimate_autocorrelation_range(coords_df, values, fallback_km = fallback_km)
}

# Spatial block bootstrap for a vector of location IDs.  Returns a resampled
# subset of the same length as n_draw, either i.i.d. or respecting spatial
# blocks of size block_km.  Includes guard against too few blocks causing
# artificially narrow confidence intervals.
spatial_block_sample_locations <- function(locations, coords_df, n_draw,
                                           block_km = NULL,
                                           max_missing_frac = BOOTSTRAP_BLOCK_MAX_MISSING_FRAC) {
  locations <- trimws(as.character(locations))
  locations <- locations[!is.na(locations) & locations != ""]
  n_draw <- as.integer(n_draw)
  if (length(locations) == 0 || n_draw < 1) return(character(0))

  if (!isTRUE(exists("ENABLE_SPATIAL_BLOCK_BOOTSTRAP")) || !isTRUE(ENABLE_SPATIAL_BLOCK_BOOTSTRAP)) {
    return(sample(locations, n_draw, replace = TRUE))
  }

  if (is.null(coords_df) || nrow(coords_df) == 0 || !all(c("location_id","lat","lon") %in% names(coords_df))) {
    stop("[SPATIAL BOOTSTRAP] Coordinate data (location_id, lat, lon) is required but missing or malformed.")
  }

  coords_df$location_id <- trimws(as.character(coords_df$location_id))
  coords_df$lat <- suppressWarnings(as.numeric(coords_df$lat))
  coords_df$lon <- suppressWarnings(as.numeric(coords_df$lon))

  coords_df <- coords_df[match(locations, coords_df$location_id), , drop = FALSE]
  if (nrow(coords_df) == 0) {
    stop("[SPATIAL BOOTSTRAP] No coordinate records matched the supplied location IDs.")
  }

  ok <- is.finite(coords_df$lat) & is.finite(coords_df$lon)
  missing_frac <- mean(!ok)
  if (!is.finite(missing_frac) || missing_frac > max_missing_frac || sum(ok) < 3) {
    stop(sprintf(
      "[SPATIAL BOOTSTRAP] Too many locations lack valid coordinates (%.0f%% missing, need at least 3 valid). Check your coordinate data.",
      missing_frac * 100
    ))
  }

  if (is.null(block_km)) block_km <- BOOTSTRAP_BLOCK_KM
  block_km <- as.numeric(block_km)
  if (!is.finite(block_km) || block_km <= 0) return(sample(locations, n_draw, replace = TRUE))

  mean_lat <- mean(coords_df$lat[ok])
  km_per_deg_lat <- 111.32
  km_per_deg_lon <- 111.32 * cos(mean_lat * pi / 180)
  km_per_deg_lon <- max(1e-6, km_per_deg_lon)

  cell_lat_deg <- block_km / km_per_deg_lat
  cell_lon_deg <- block_km / km_per_deg_lon
  cell_lat_deg <- max(1e-8, cell_lat_deg)
  cell_lon_deg <- max(1e-8, cell_lon_deg)

  cell_id <- rep(NA_character_, length(locations))
  cell_id[ok] <- paste0(
    floor(coords_df$lat[ok] / cell_lat_deg), "_", floor(coords_df$lon[ok] / cell_lon_deg)
  )

  blocks <- split(locations[ok], cell_id[ok], drop = TRUE)
  block_ids <- names(blocks)

  min_blocks_needed <- max(3L, as.integer(ceiling(n_draw / 3)))
  if (length(block_ids) < 2) {
    stop(sprintf(
      "[SPATIAL BOOTSTRAP] All %d locations fall into a single spatial block (block_km = %.1f km). Increase block size or check coordinate spread.",
      n_draw, block_km
    ))
  }
  if (length(block_ids) < min_blocks_needed) {
    stop(sprintf(
      "[SPATIAL BOOTSTRAP] Only %d spatial blocks for %d locations (need at least %d). Block size %.1f km is too large relative to the spatial extent — reduce BOOTSTRAP_BLOCK_KM or check coordinates.",
      length(block_ids), n_draw, min_blocks_needed, block_km
    ))
  }

  pool <- character(0)
  avg_block_size <- mean(lengths(blocks))
  n_blocks_draw <- max(1L, ceiling(n_draw / max(1, avg_block_size)))
  sampled_blocks <- sample(block_ids, size = n_blocks_draw, replace = TRUE)
  pool <- unlist(blocks[sampled_blocks], use.names = FALSE)

  while (length(pool) < n_draw) {
    sb <- sample(block_ids, size = 1L, replace = TRUE)
    pool <- c(pool, blocks[[sb]])
  }

  sample(pool, n_draw, replace = FALSE)
}

analyze_library_similarity <- function(mesma_lib, compressed_templates_accessor, grid_type = "full") {
  cat("\n=== INTER-CLASS VARIANT SIMILARITY ANALYSIS ===\n")

  # DIAGNOSTICS
  cat(sprintf("[DIAGNOSTIC] names(mesma_lib): %s\n", paste(names(mesma_lib), collapse=", ")))
  cat(sprintf("[DIAGNOSTIC] names(compressed_templates_accessor): %s\n", paste(names(compressed_templates_accessor), collapse=", ")))
  
  # Collect compressed 'full' templates in a safe manner
  all_ids <- character()
  all_vegs <- character()
  all_vecs <- list()
  all_nsamples <- numeric()

  for (veg in names(mesma_lib)) {
    if (is.null(mesma_lib[[veg]])) next
    # Try to access raw_lib_templates from global environment for fallback
    rt <- if (exists("raw_lib_templates", envir = globalenv()) && !is.null(get("raw_lib_templates", envir = globalenv())[[veg]])) get("raw_lib_templates", envir = globalenv())[[veg]] else NULL
    for (variant in filter_variants_by_min_samples(mesma_lib[[veg]], min_samples = MIN_ENDMEMBER_SAMPLES, veg = veg, raw_template = rt)) {
      vid <- if (!is.null(variant$variant_id)) variant$variant_id else if (!is.null(variant$id)) variant$id else NA_character_
      
      # If vid is still NA, try to construct one
      if (is.na(vid)) {
         # heuristic to generate a unique ID
         vid <- paste0(veg, "_", length(all_ids) + 1)
      }

      vec <- NULL
      
      # Robust lookup strategy
      # 1. Try exact match
      if (!is.null(compressed_templates_accessor[[veg]]) && !is.null(compressed_templates_accessor[[veg]][[vid]]) && !is.null(compressed_templates_accessor[[veg]][[vid]][["full"]])) {
         vec <- compressed_templates_accessor[[veg]][[vid]][["full"]]
      }
      
      # 2. Try case-insensitive veg match if failed
      if (is.null(vec)) {
         c_vegs <- names(compressed_templates_accessor)
         # Find veg key that matches case-insensitively
         v_match <- c_vegs[tolower(c_vegs) == tolower(veg)]
         if (length(v_match) > 0) {
            used_veg_key <- v_match[1]
            acc_veg <- compressed_templates_accessor[[used_veg_key]]
            
            # Try exact vid match in the resolved veg
            if (!is.null(acc_veg[[vid]]) && !is.null(acc_veg[[vid]][["full"]])) {
               vec <- acc_veg[[vid]][["full"]]
            } else {
               # 3. Try case-insensitive vid match
               c_vids <- names(acc_veg)
               vid_match <- c_vids[tolower(c_vids) == tolower(vid)]
               if (length(vid_match) > 0) {
                  vec <- acc_veg[[vid_match[1]]][["full"]]
               }
            }
         }
      }

      # 3. Fallback to using raw_mat from variant if available
      if (is.null(vec) && !is.null(variant$raw_mat)) {
         vec <- as.numeric(variant$raw_mat)
         vec[!is.finite(vec)] <- NA_real_
      }
      # 4. Fallback to using vec from variant if available
      if (is.null(vec) && !is.null(variant$vec)) {
         vec <- as.numeric(variant$vec)
         vec[!is.finite(vec)] <- NA_real_
      }


      if (!is.null(vec) && length(vec) > 0 && !is.na(vid)) {
        all_vecs[[length(all_vecs) + 1]] <- as.numeric(vec)
        all_ids <- c(all_ids, vid)
        all_vegs <- c(all_vegs, veg)
        all_nsamples <- c(all_nsamples, if (!is.null(variant$n_samples)) variant$n_samples else 0)
      }
    }
  }

  cat(sprintf("[DIAGNOSTIC] Collected %d variant vectors for similarity analysis.\n", length(all_vecs)))

  # Barren similarity filtering is applied at the observation level (before endmember building)
  # so it is not repeated here. See [BARREN PRE-FILTER] in fit_veg_mixture_mesma.R.

  # Print variant counts per vegetation type
  if (length(all_vecs) > 0) {
    veg_counts <- table(all_vegs)
    cat("[INFO] Variants per vegetation type (after filtering):\n")
    print(veg_counts)
    
    if (min(veg_counts) < 2 && length(veg_counts) > 1) {
       cat("[INFO] Note: Some vegetation types have few variants. All variants will be retained for visualization.\n")
    }
  }

  if (length(all_vecs) < 2) {
    cat("[WARN] Not enough variant templates for similarity heatmap.\n")
    write.csv(data.frame(info = "insufficient_variants", n_variants = length(all_vecs)), file.path(OUT_DIR, "variant_similarity_matrix_info.csv"), row.names = FALSE)
    cat("==============================================\n\n")
    return(invisible(NULL))
  }

  # Normalize lengths by truncating to shortest vector if needed
  lens <- sapply(all_vecs, length)
  if (length(unique(lens)) > 1) {
    min_len <- min(lens)
    cat(sprintf("[WARN] Variant vectors have differing lengths; truncating to min length %d\n", min_len))
    all_vecs <- lapply(all_vecs, function(v) v[1:min_len])
  }

  n_v <- length(all_vecs)
  
  # Vectorized Similarity Calculation
  vec_mat <- do.call(rbind, all_vecs)
  vec_mat[!is.finite(vec_mat)] <- 0
  row_norms <- sqrt(rowSums(vec_mat^2))
  row_norms[row_norms == 0] <- 1
  vec_mat_norm <- vec_mat / row_norms
  sim_mat <- tcrossprod(vec_mat_norm)
  sim_mat[sim_mat > 1.01] <- 1.01
  sim_mat[sim_mat < -0.01] <- -0.01
  
  rownames(sim_mat) <- all_ids
  colnames(sim_mat) <- all_ids

  sim_df <- as.data.frame(as.table(sim_mat))
  colnames(sim_df) <- c("Var1", "Var2", "Similarity")
  
  # AGENT FIX: Ensure ordered_ids is defined
  ord_idx <- order(all_vegs, all_ids)
  ordered_ids <- all_ids[ord_idx]
  ordered_vegs <- all_vegs[ord_idx]
  
  # AGENT CHANGE: Add sample counts to weight tile thickness and allocate contiguous positions so there are no gaps between variant bars
  count_map <- setNames(all_nsamples, all_ids)

  # Determine ordered unique ids and counts
  counts_ordered <- as.numeric(count_map[ordered_ids])
  counts_ordered[is.na(counts_ordered)] <- 0
  total_counts <- sum(counts_ordered)
  if (!is.finite(total_counts) || total_counts <= 0) {
    n_ord <- length(ordered_ids)
    widths_prop <- rep(1 / n_ord, n_ord)
  } else {
    widths_prop <- counts_ordered / total_counts
    # Avoid zero widths (give tiny epsilon to zero-count variants so they remain visible)
    eps <- 1e-6
    widths_prop[widths_prop == 0] <- eps
    widths_prop <- widths_prop / sum(widths_prop)
  }

  # Compute continuous positions in [0,1]
  cumw <- cumsum(widths_prop)
  lefts <- c(0, cumw[-length(cumw)])
  rights <- cumw
  names(lefts) <- ordered_ids
  names(rights) <- ordered_ids

  # For Y axis we use reversed ordering (top-to-bottom)
  lefts_rev <- rev(lefts)
  rights_rev <- rev(rights)
  names(lefts_rev) <- rev(ordered_ids)
  names(rights_rev) <- rev(ordered_ids)

  # Map rectangle coordinates into sim_df
  sim_df$xmin <- lefts[as.character(sim_df$Var1)]
  sim_df$xmax <- rights[as.character(sim_df$Var1)]
  sim_df$ymin <- lefts_rev[as.character(sim_df$Var2)]
  sim_df$ymax <- rights_rev[as.character(sim_df$Var2)]

  # Compute centers for axis ticks/labels
  x_centers <- (lefts + rights) / 2
  y_centers <- (lefts_rev + rights_rev) / 2

  # Boundaries at vegetation type changes (in cumulative units)
  if (length(ordered_vegs) > 1) {
    change_idx <- which(ordered_vegs[-1] != ordered_vegs[-length(ordered_vegs)])
    vlines_x <- cumsum(widths_prop)[change_idx]
    hlines_y <- cumsum(widths_prop)[change_idx]
  } else {
    vlines_x <- numeric(0)
    hlines_y <- numeric(0)
  }

  # Build heatmap using fixed rectangles (no gaps between adjacent variants)
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  p_heat <- ggplot(sim_df) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Similarity), color = NA) +
    # Add thick black lines at vegetation boundaries
    geom_vline(xintercept = vlines_x, color = "black", size = 1.2) +
    geom_hline(yintercept = hlines_y, color = "black", size = 1.2) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.9, limits = c(-0.01, 1.01), na.value = "white") +
    scale_x_continuous(expand = c(0,0), breaks = x_centers, labels = ordered_ids) +
    scale_y_continuous(expand = c(0,0), breaks = y_centers, labels = rev(ordered_ids)) +
    coord_equal() +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6), axis.text.y = element_text(size = 6)) +
    labs(title = "Variant Cosine Similarity", x = NULL, y = NULL)

  tryCatch({
    ggsave(file.path(OUT_DIR, "variant_similarity_heatmap.png"), p_heat, width = 12, height = 10)
    write.csv(sim_mat, file.path(OUT_DIR, "variant_similarity_matrix.csv"))
    cat(sprintf("Saved similarity heatmap to: %s\n", file.path(OUT_DIR, "variant_similarity_heatmap.png")))
  }, error = function(e) {
    cat(sprintf("[ERROR] Failed to generate similarity heatmap: %s\n", e$message))
  })

  # == AGENT FIX: Print #loc-years per variant ==
  cat("\n=== VARIANT SAMPLE SIZES (Loc-Years) ===\n")
  sample_stats <- data.frame(
    Vegetation = all_vegs,
    VariantID = all_ids,
    LocYears = all_nsamples
  )
  # Sort by vegetation then ID for display
  sample_stats <- sample_stats[order(sample_stats$Vegetation, sample_stats$VariantID), ]
  print(sample_stats, row.names = FALSE)
  write.csv(sample_stats, file.path(OUT_DIR, "variant_sample_sizes.csv"), row.names = FALSE)
  cat("=========================================\n\n")

  cat("==============================================\n\n")
  invisible(TRUE)
}



# global flag controlling whether years 1992-1999 (and their shading)
# should be treated as 'excluded'.  By default we now leave this FALSE so that
# plots and analyses include all years; set EXCLUDE_PRE2000 <- TRUE in your
# environment to restore the original exclusions.
if (!exists("EXCLUDE_PRE2000", inherits = TRUE)) EXCLUDE_PRE2000 <- FALSE

# Helper: add a shaded rectangle covering excluded years (1992-1999) for ggplot2 time-series plots.
# Usage: + add_excluded_years_shade(is_date = TRUE)  # x axis is Date
#        + add_excluded_years_shade(is_date = FALSE) # x axis is numeric (year)
add_excluded_years_shade <- function(start_year = 1992, end_year = 1999, is_date = FALSE, fill = "grey70", alpha = 0.35) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(list())
  if (!isTRUE(EXCLUDE_PRE2000)) {
    # shading disabled by global flag; return empty list so callers can concatenate safely
    return(list())
  }
  # Two shaded regions: 1992-1999 and 2007-2009
  if (is_date) {
    xmin1 <- as.Date(paste0(start_year, "-01-01"))
    xmax1 <- as.Date(paste0(end_year, "-12-31"))
    xmin2 <- as.Date("2007-01-01")
    xmax2 <- as.Date("2009-12-31")
    return(list(
      ggplot2::annotate("rect", xmin = xmin1, xmax = xmax1, ymin = -Inf, ymax = Inf, fill = fill, alpha = alpha),
      ggplot2::annotate("rect", xmin = xmin2, xmax = xmax2, ymin = -Inf, ymax = Inf, fill = fill, alpha = alpha)
    ))
  } else {
    return(list(
      ggplot2::annotate("rect", xmin = start_year, xmax = end_year, ymin = -Inf, ymax = Inf, fill = fill, alpha = alpha),
      ggplot2::annotate("rect", xmin = 2007, xmax = 2009, ymin = -Inf, ymax = Inf, fill = fill, alpha = alpha)
    ))
  }
}

# Helper: add vertical lines at year boundaries for ggplot2 time-series plots.
# Usage: + add_year_lines(is_date = TRUE)  # x axis is Date
#        + add_year_lines(is_date = FALSE) # x axis is numeric (year)
add_year_lines <- function(start_year = 1985, end_year = 2025, is_date = FALSE, color = "grey50", linetype = "dashed", alpha = 0.5) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  years <- seq(start_year, end_year, by = 1)
  if (is_date) {
    xintercepts <- as.Date(paste0(years, "-01-01"))
    return(ggplot2::geom_vline(xintercept = as.numeric(xintercepts), color = color, linetype = linetype, alpha = alpha))
  } else {
    return(ggplot2::geom_vline(xintercept = years, color = color, linetype = linetype, alpha = alpha))
  }
}


ensure_library_and_templates <- function(force = FALSE) {
  # Attempt to construct or expose `mesma_lib` and `compressed_templates_accessor`
  # so visualizations can run reliably. If these objects already exist, this function is a no-op.
  # If both objects already exist in globalenv, nothing to do
  if (!isTRUE(force) && exists("mesma_lib", envir = globalenv()) && (exists("compressed_templates_accessor", envir = globalenv()) || exists(".COMPRESSED_TEMPLATES_ACCESSOR", envir = globalenv()))) {
    return(invisible(TRUE))
  }

  cat("[INFO] Ensuring MESMA library and compressed templates are available for visualization...\n")

  # Attempt to use local variables if available in parent frame
  ms <- NULL; ct <- NULL
  if (exists("mesma_lib", envir = globalenv())) ms <- get("mesma_lib", envir = globalenv())
  if (exists("compressed_templates_accessor", envir = globalenv())) ct <- get("compressed_templates_accessor", envir = globalenv())
  if (is.null(ct) && exists(".COMPRESSED_TEMPLATES_ACCESSOR", envir = globalenv())) ct <- get(".COMPRESSED_TEMPLATES_ACCESSOR", envir = globalenv())
  if (is.null(ms) && exists("mesma_lib", envir = parent.frame())) ms <- get("mesma_lib", envir = parent.frame())
  if (is.null(ct) && exists("compressed_templates_accessor", envir = parent.frame())) ct <- get("compressed_templates_accessor", envir = parent.frame())
  if (is.null(ct) && exists(".COMPRESSED_TEMPLATES_ACCESSOR", envir = parent.frame())) ct <- get(".COMPRESSED_TEMPLATES_ACCESSOR", envir = parent.frame())

  # If still missing, attempt to build them from available inputs


  if (is.null(ms) || is.null(ct)) {
    cat("[WARN] After attempts, 'mesma_lib' or 'compressed_templates_accessor' still not available. Some visualizations may be skipped.\n")
    return(invisible(FALSE))
  }

  invisible(TRUE)
}

ensure_variant_similarity_heatmap <- function(force = FALSE) {
  if (exists("VARIANT_SIMILARITY_HEATMAP_DONE", envir = globalenv()) && isTRUE(get("VARIANT_SIMILARITY_HEATMAP_DONE", envir = globalenv())) && !isTRUE(force)) {
    return(invisible(TRUE))
  }

  # Ensure that core objects are present (try to create them if missing)
  ensure_library_and_templates(force = force)

  ms <- NULL; ct <- NULL
  if (exists("mesma_lib", envir = globalenv())) ms <- get("mesma_lib", envir = globalenv())
  if (exists("compressed_templates_accessor", envir = globalenv())) ct <- get("compressed_templates_accessor", envir = globalenv())
  if (is.null(ct) && exists(".COMPRESSED_TEMPLATES_ACCESSOR", envir = globalenv())) ct <- get(".COMPRESSED_TEMPLATES_ACCESSOR", envir = globalenv())

  if (is.null(ms) || is.null(ct)) {
    cat("[WARN] Cannot generate variant similarity heatmap: 'mesma_lib' or 'compressed_templates_accessor' not available yet\n")
    return(invisible(FALSE))
  }

  tryCatch({
    analyze_library_similarity(ms, ct, grid_type = "full")
    assign("VARIANT_SIMILARITY_HEATMAP_DONE", TRUE, envir = globalenv())
    cat("[INFO] Variant similarity heatmap generated (ensure_variant_similarity_heatmap).\n")
  }, error = function(e) {
    cat(sprintf("[ERROR] Failed to generate variant similarity heatmap: %s\n", e$message))
  })
  invisible(TRUE)
}


# --- Reproducible RNG helpers -------------------------------------------------
# Provide a single, overridable MESMA seed so all scripts and functions can be
# made comparable by setting the MESMA_SEED environment variable (or calling
# set_mesma_seed()). Defaults to 42 for backward-compatibility.
#
# Usage:
#  - call set_mesma_seed() early in a top-level script
#  - use get_mesma_seed(offset) when a deterministic distinct seed is needed
#  - existing code that uses `seed` or option("mesma.seed") will continue to work
set_mesma_seed <- function(base = NULL, announce = TRUE, set_env_vars = TRUE) {
  # Resolve base seed: explicit arg -> MESMA_SEED env var -> default 42
  if (is.null(base) || !is.finite(base)) {
    ev <- Sys.getenv("MESMA_SEED", unset = NA_character_)
    if (!is.na(ev) && nzchar(ev)) {
      base <- suppressWarnings(as.integer(ev))
      if (is.na(base)) base <- 42L
    } else {
      base <- 42L
    }
  }
  base <- as.integer(base)

  # Use L'Ecuyer for better parallel reproducibility and future::future.seed support
  try({ RNGkind("L'Ecuyer-CMRG") }, silent = TRUE)
  set.seed(base)

  # record seed in options and a global variable for compatibility
  options(mesma.seed = base)
  assign("seed", base, envir = globalenv())

  # Ensure common env vars used by scripts are seeded (if not explicitly provided)
  if (isTRUE(set_env_vars)) {
    if (nzchar(Sys.getenv("FVC_SAMPLING_SEED", "")) == FALSE) Sys.setenv(FVC_SAMPLING_SEED = as.character(base))
  }

  if (isTRUE(announce)) message(sprintf("[MESMA] RNG initialized with MESMA_SEED=%d (override with MESMA_SEED env var)" , base))
  invisible(base)
}

# Return a deterministic seed derived from the master seed (useful for offsets)
get_mesma_seed <- function(offset = 0L) {
  base <- as.integer(getOption("mesma.seed", Sys.getenv("MESMA_SEED", unset = "42")))
  if (!is.finite(base) || length(base) == 0) base <- 42L
  offset <- as.integer(offset)
  # keep result in integer range
  as.integer((base + offset) %% .Machine$integer.max)
}


# --- NDDI (dust) threshold configuration ------------------------------------
# Make the NDDI dust-contamination threshold tunable via environment variable
# `MESMA_NDDI_THRESHOLD` or the global `NDDI_DUST_THRESHOLD` variable.
# Default is 0.28 (updated per user request).
get_nddi_threshold <- function() {
  v <- as.numeric(Sys.getenv("MESMA_NDDI_THRESHOLD", unset = "0.28"))
  if (!is.finite(v)) v <- 0.28
  v
}
# Cached default value so scripts can reference `NDDI_DUST_THRESHOLD` directly
NDDI_DUST_THRESHOLD <- get_nddi_threshold()

# Helper to format threshold for messages
.nddi_thresh_fmt <- function() sprintf("%.3f", as.numeric(get_nddi_threshold()))

# --- Shared vegetation label canonicalization helpers -------------------------

# Canonical list of herb-group labels (map these to "herbs")
HERBS_GROUP <- c("herbs", "alhagi", "salicornia", "halocnemum", "phragmites")

# Agriculture alias mapping (map these to "agriculture")
AGRICULTURE_ALIASES <- c("agri", "agric", "agriculture", "agricultural")

# Canonicalize vegetation labels in a data.frame column.
# Fixes typos, maps herb-group species -> "herbs", maps agriculture aliases -> "agriculture".
# `veg_col` is determined automatically if NULL.
canonicalize_veg_labels <- function(df, veg_col = NULL) {
  if (is.null(veg_col)) {
    veg_col <- intersect(c("Veg", "vegetation"), names(df))[1]
  }
  if (is.na(veg_col) || !veg_col %in% names(df)) return(df)

  v <- tolower(trimws(as.character(df[[veg_col]])))
  # Fix typos
  v[v == "tamairx"] <- "tamarix"
  # Map herb-group
  v[v %in% HERBS_GROUP] <- "herbs"
  # Map agriculture aliases
  v[v %in% AGRICULTURE_ALIASES] <- "agriculture"
  df[[veg_col]] <- v
  df
}

# --- Bootstrap results compilation helper ------------------------------------
# Compiles a list of [B x n_years] bootstrap matrices into a single data.frame
# with columns: year, Veg, global_coef, se, coef_025, coef_975, method, n_locations.
compile_bootstrap_results <- function(veg_boot_res, years, unique_loc_years, method_name) {
  final_results <- list()
  for (v in names(veg_boot_res)) {
    mat <- veg_boot_res[[v]]
    if (is.null(mat) || ncol(mat) == 0 || is.null(colnames(mat)) || length(colnames(mat)) == 0) next
    df_res <- data.frame(
      year = as.integer(colnames(mat)),
      Veg = v,
      global_coef = apply(mat, 2, mean, na.rm = TRUE),
      se = apply(mat, 2, sd, na.rm = TRUE),
      coef_025 = apply(mat, 2, quantile, 0.025, na.rm = TRUE),
      coef_975 = apply(mat, 2, quantile, 0.975, na.rm = TRUE),
      method = method_name
    )
    n_locs_per_year <- sapply(years, function(y) sum(unique_loc_years$pheno_year == y))
    df_res$n_locations <- n_locs_per_year[match(df_res$year, years)]
    final_results[[v]] <- df_res
  }
  dplyr::bind_rows(final_results)
}

# --- PPI norm/backup helper --------------------------------------------------
# Backs up raw PPI to PPI_raw, then creates ppi_norm clamped to [0,1].
backup_and_normalize_ppi <- function(df, label = "") {
  prefix <- if (nzchar(label)) paste0(label, ": ") else ""

  # Backup raw PPI
  if (!"PPI_raw" %in% names(df)) {
    df$PPI_raw <- df$PPI
    cat(sprintf("[NOTICE] %sBacked up raw PPI values to 'PPI_raw' before normalization.\n", prefix))
  }

  # Compute ppi_norm
  if (!"ppi_norm" %in% names(df)) df$ppi_norm <- NA_real_
  if ("PPI_raw" %in% names(df) && any(is.finite(df$PPI_raw))) {
    df$ppi_norm <- pmin(pmax(df$PPI_raw, 0), 1)
    cat(sprintf("[PPI NORM] %sCreated 'ppi_norm' from 'PPI_raw' and clamped to [0,1]\n", prefix))
  } else if ("PPI" %in% names(df) && any(is.finite(df$PPI))) {
    df$ppi_norm <- pmin(pmax(df$PPI, 0), 1)
    warning(sprintf("%s'PPI_raw' not found - computed 'ppi_norm' from 'PPI' (may be z-scored); values were clamped to [0,1].", prefix))
  } else {
    df$ppi_norm <- NA_real_
    cat(sprintf("[PPI NORM] %sNo PPI or PPI_raw available to compute 'ppi_norm' (all NA)\n", prefix))
  }
  df
}


# Compute a per-location DVI soil baseline for PPI.
# Baseline is defined as the median of the lowest quantile_p fraction of DVI
# observations within each location. This is deterministic, per-location,
# and avoids any constant/default soil baseline.
compute_dvi_soil_per_location <- function(df, quantile_p = 0.10, min_samples = 5L) {
  stopifnot(!is.null(df), nrow(df) > 0L, "location_id" %in% names(df))
  if (!"DVI" %in% names(df) && all(c("nir", "red") %in% names(df))) {
    df <- df %>% dplyr::mutate(DVI = as.numeric(nir) - as.numeric(red))
  }
  stopifnot("DVI" %in% names(df))

  soil_df <- df %>%
    dplyr::filter(is.finite(DVI)) %>%
    dplyr::group_by(location_id) %>%
    dplyr::filter(dplyr::n() >= min_samples) %>%
    dplyr::summarise(
      dvi_soil = {
        vals <- DVI
        q <- suppressWarnings(quantile(vals, probs = quantile_p,
                                       na.rm = TRUE, names = FALSE, type = 7))
        median(vals[vals <= q], na.rm = TRUE)
      },
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(dvi_soil))

  out <- df %>% dplyr::left_join(soil_df, by = "location_id") %>% dplyr::pull(dvi_soil)
  if (any(is.na(out) & is.finite(df$DVI))) {
    bad_locs <- unique(as.character(df$location_id[is.na(out) & is.finite(df$DVI)]))
    stop(sprintf("[PPI] Cannot compute per-location dvi_soil for %d rows across %d locations (example locs: %s)",
                 sum(is.na(out) & is.finite(df$DVI)), length(bad_locs),
                 paste(head(bad_locs, 10), collapse = ", ")))
  }
  out
}

apply_stored_normalization <- function(df, norm_params, cols = unique(c(OPTIMAL_INDICES, RAW_BANDS)), lat_default = 40.2) {
  cat("Applying stored normalization parameters to data...\n")

  if (!"zenith.angle" %in% names(df)) df[["zenith.angle"]] <- NA_real_
  if (is.null(norm_params[["INDEX_SCALES"]]) || length(norm_params[["INDEX_SCALES"]]) == 0) {
    stop("[INDEX_SCALES] Missing stored normalization parameters (INDEX_SCALES)")
  }

  if (!"PPI" %in% names(df) || all(!is.finite(df[["PPI"]]))) {
    if (all(c("nir", "red") %in% names(df)) && !"DVI" %in% names(df)) {
      df[["DVI"]] <- as.numeric(df[["nir"]]) - as.numeric(df[["red"]])
    }
    if (exists("add_ppi_columns")) {
      dvi_soil_vec <- compute_dvi_soil_per_location(df)
      df <- add_ppi_columns(df, dvi_soil = dvi_soil_vec)
      cat("[PPI] Added PPI to data before applying stored normalization (per-location dvi_soil + per-location M).\n")
    }
  }
  if (!"PPI" %in% names(df) || all(!is.finite(df[["PPI"]]))) {
    stop("[PPI] Missing or all PPI values non-finite after attempted auto-add; refusing to continue")
  }
  if (!"PPI_raw" %in% names(df)) df[["PPI_raw"]] <- df[["PPI"]]

  cat(sprintf("[apply_stored_normalization] Applying INDEX_SCALES to %d indices\n", length(norm_params[["INDEX_SCALES"]])))
  n_scaled <- 0L
  for (col in names(norm_params[["INDEX_SCALES"]])) {
    if (col %in% names(df)) {
      params <- norm_params[["INDEX_SCALES"]][[col]]
      if (!is.list(params) || !all(c("mean", "sd") %in% names(params))) {
        stop(sprintf("[INDEX_SCALES] Invalid params for index '%s' (expected list(mean, sd))", col))
      }
      mu <- params[["mean"]]
      sigma <- params[["sd"]]
      if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0) {
        stop(sprintf("[INDEX_SCALES] Non-finite or non-positive mean/sd for index '%s'", col))
      }
      df[[col]] <- (df[[col]] - mu) / sigma
      n_scaled <- n_scaled + 1L
    }
  }
  cat(sprintf("[apply_stored_normalization] Z-scored %d feature columns\n", n_scaled))

  df
}

normalize_veg_name <- function(x) {
  x <- as.character(x)
  x[!nzchar(trimws(x))] <- NA_character_
  tolower(trimws(x))
}
