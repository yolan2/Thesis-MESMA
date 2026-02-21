
# --- Consolidated shared helpers (moved here) ---------------------------------
# General small utilities used across multiple scripts. These were de-duplicated
# from inline copies in other files and are kept here as the canonical source.

# Additional helpers migrated from fit_veg_mixture_mesma.R to reduce redundancy.
# These functions were previously defined in the main script but are generic
# enough to be reused by other tools; moving them here centralizes maintenance.

safe_as_numeric <- function(x) {
  as.numeric(as.character(x))
}

make_location_id <- function(lon, lat) {
  lon <- as.numeric(lon)
  lat <- as.numeric(lat)
  invalid_mask <- !is.finite(lon) | !is.finite(lat)
  res <- sprintf("L_%0.6f_%0.6f", round(lat, 6), round(lon, 6))
  res[invalid_mask] <- NA_character_
  res
}

assign_pheno_year <- function(d) {
  d <- as.Date(d)
  # phenological year begins on May 1 (boundary = 30 April).
  # Dates in Jan-Apr are assigned to the previous calendar year.
  ifelse(is.na(d), NA_integer_, ifelse(lubridate::month(d) >= 5, lubridate::year(d), lubridate::year(d) - 1))
}

# Perform Whittaker smoothing on a 1‑D series.
#
# The smoothed vector `z` is obtained by minimising
#   \sum_i (y_i - z_i)^2 + \lambda \sum_i (\Delta^2 z_i)^2
# where \Delta^2 z_i = z_i - 2 z_{i-1} + z_{i-2} is the second
# difference.  The first term enforces fidelity to the original data
# and the second term penalises roughness; because we penalise the
# *second* difference the result is a piecewise‑linear (once
# differentiable) trend rather than merely encouraging consecutive
# points to be close.  A larger `lambda` yields a smoother (flatter)
# curve.  This is the standard second‑order Whittaker smoother; if you
# prefer the first‑order formulation
# \sum (\Delta z_i)^2 replace `diff(..., differences = 2)` with a
# first‑difference matrix.
#
# Arguments:
#   y: numeric vector of observations (can contain NAs).
#   lambda: smoothing penalty (higher -> smoother).
#
whittaker_smooth <- function(y, lambda = 500) {
  n <- length(y)
  if (n < 3) return(rep(NA_real_, n))

  yy <- as.numeric(y)
  w <- as.numeric(is.finite(yy))
  yy[!is.finite(yy)] <- 0

  W <- diag(w, nrow = n, ncol = n)
  D <- diff(diag(n), differences = 2)
  A <- W + lambda * crossprod(D)
  b <- W %*% yy

  tryCatch(
    as.numeric(base::solve(A, b)),
    error = function(e) {
      tryCatch(as.numeric(base::qr.solve(A, b)), error = function(e2) rep(NA_real_, n))
    }
  )
}

# Iteratively reweighted Whittaker smoother (IRW)
# Uses Huber-style weights to reduce influence of outliers.
# Returns a list with: `z` (smoothed values), `weights` (final weights), `residuals`.
whittaker_smooth_irw <- function(y, lambda = 500, max_iter = 10L, k = 1.345, tol = 1e-6) {
  n <- length(y)
  if (n < 3) return(list(z = rep(NA_real_, n), weights = rep(0, n), residuals = rep(NA_real_, n)))

  yy <- as.numeric(y)
  finite_mask <- is.finite(yy)
  yy[!finite_mask] <- 0

  # initial weights: 1 for finite, 0 for NA
  w <- as.numeric(finite_mask)
  z_prev <- rep(0, n)

  D <- diff(diag(n), differences = 2)
  DtD <- crossprod(D)

  for (it in seq_len(as.integer(max_iter))) {
    W <- diag(w, nrow = n, ncol = n)
    A <- W + lambda * DtD
    b <- W %*% yy
    z <- tryCatch(as.numeric(base::solve(A, b)), error = function(e) as.numeric(base::qr.solve(A, b)))

    resid <- yy - z
    # robust scale estimate (MAD) on finite residuals; fallback to sd
    rfin <- resid[finite_mask]
    s <- stats::mad(rfin, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) s <- stats::sd(rfin, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) s <- 1

    u <- resid / (k * s)
    # Huber weight: psi(u)/u where psi(u)=u for |u|<=1, =sign(u) otherwise
    w_new <- rep(0, n)
    small <- abs(u) <= 1
    w_new[small] <- 1
    w_new[!small & finite_mask] <- 1 / abs(u[!small & finite_mask])
    w_new[!finite_mask] <- 0

    # check convergence (max change in z)
    if (max(abs(z - z_prev), na.rm = TRUE) < tol) {
      return(list(z = z, weights = w_new, residuals = resid))
    }

    w <- as.numeric(w_new)
    z_prev <- z
  }

  # final solve (in case loop ended by iter count)
  W <- diag(w, nrow = n, ncol = n)
  A <- W + lambda * DtD
  b <- W %*% yy
  z <- tryCatch(as.numeric(base::solve(A, b)), error = function(e) as.numeric(base::qr.solve(A, b)))
  resid <- yy - z
  list(z = z, weights = w, residuals = resid)
}

remove_large_outliers_whittaker <- function(df, candidates, lambda = 500, min_group_rows = 3L, min_rows_with_date = 5L) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if (is.null(candidates) || length(candidates) == 0) {
    cat("[WHITTAKER] No candidate indices found for smoothing; skipping\n")
    return(df)
  }
  candidates <- intersect(candidates, names(df))
  if (length(candidates) == 0) {
    cat("[WHITTAKER] No candidate indices found for smoothing; skipping\n")
    return(df)
  }
  if (!"location_id" %in% names(df)) {
    cat("[WHITTAKER] 'location_id' missing from data; skipping smoothing\n")
    return(df)
  }

  if (!"pheno_year" %in% names(df) && "date" %in% names(df)) {
    df$pheno_year <- assign_pheno_year(df$date)
  }

  grp <- interaction(df$location_id, ifelse(is.na(df$pheno_year), "NA", as.character(df$pheno_year)), drop = TRUE)
  groups <- split(seq_len(nrow(df)), grp)
  total_smoothed <- 0L
  n_groups <- length(groups)

  for (g in seq_along(groups)) {
    rows <- groups[[g]]
    sub <- df[rows, , drop = FALSE]
    if (length(rows) < min_group_rows) next

    has_date <- "date" %in% names(sub) && any(!is.na(sub$date))
    if (!has_date || length(rows) < min_rows_with_date) next

    sub$doy <- as.numeric(format(sub$date, "%j"))

    for (col in candidates) {
      if (!is.numeric(sub[[col]])) next
      colv <- sub[[col]]

      idx_all <- which(is.finite(sub$doy))
      if (length(idx_all) < min_rows_with_date) next

      ord <- order(sub$doy[idx_all], na.last = NA)
      idx_ord <- idx_all[ord]
      y_ord <- as.numeric(colv[idx_ord])
      if (sum(is.finite(y_ord)) < min_rows_with_date) next

      tryCatch({
        fit_res <- whittaker_smooth_irw(y_ord, lambda = lambda)
        fit_ord <- fit_res$z
        if (is.null(fit_ord) || all(!is.finite(fit_ord))) next

        smooth_ord <- is.finite(fit_ord)
        if (any(smooth_ord, na.rm = TRUE)) {
          colv[idx_ord[smooth_ord]] <- fit_ord[smooth_ord]
          total_smoothed <- total_smoothed + sum(smooth_ord, na.rm = TRUE)
          sub[[col]] <- colv
        }
      }, error = function(e) {
        # If smoothing fails for this column/group, keep original values.
      })
    }

    df[rows, names(sub)] <- sub
  }

  if (total_smoothed > 0) {
    cat(sprintf("[WHITTAKER] Smoothed %d values across %d groups\n", total_smoothed, n_groups))
  }

  df
}

pheno_doy <- function(d) {
  d <- tryCatch(as.Date(d), error = function(e) NA)
  month <- lubridate::month(d)
  ifelse(is.na(d), NA_integer_,
         # use May 1 as Day 1 of the phenological year
         ifelse(month >= 5,
                as.integer(d - as.Date(paste0(lubridate::year(d), "-05-01"))) + 1L,
                as.integer(d - as.Date(paste0(lubridate::year(d) - 1, "-05-01"))) + 1L
         )
  )
}

# Compute a per-location DVI soil baseline for PPI (deterministic, per-location)
compute_dvi_soil_per_location <- function(df, quantile_p = 0.10, min_samples = 5L) {
  if (is.null(df) || nrow(df) == 0) stop("[PPI] compute_dvi_soil_per_location: empty df")
  if (!"location_id" %in% names(df)) stop("[PPI] compute_dvi_soil_per_location: missing location_id")
  if (!"DVI" %in% names(df) && all(c("nir", "red") %in% names(df))) df$DVI <- as.numeric(df$nir) - as.numeric(df$red)
  if (!"DVI" %in% names(df)) stop("[PPI] compute_dvi_soil_per_location: missing DVI (and nir/red not available)")

  locs <- unique(as.character(df$location_id))
  dvi_soil_vec <- rep(NA_real_, nrow(df))
  for (loc in locs) {
    idx <- which(as.character(df$location_id) == loc)
    vals <- df$DVI[idx]
    vals <- vals[is.finite(vals)]
    if (length(vals) < as.integer(min_samples)) next
    q <- suppressWarnings(as.numeric(stats::quantile(vals, probs = quantile_p, na.rm = TRUE, names = FALSE, type = 7)))
    if (!is.finite(q)) next
    low_vals <- vals[vals <= q]
    if (length(low_vals) < 2L) next
    soil <- suppressWarnings(as.numeric(stats::median(low_vals, na.rm = TRUE)))
    if (!is.finite(soil)) next
    dvi_soil_vec[idx] <- soil
  }

  need_idx <- is.finite(df$DVI) & !is.finite(dvi_soil_vec)
  if (any(need_idx)) {
    bad_locs <- unique(as.character(df$location_id[need_idx]))
    stop(sprintf("[PPI] Cannot compute per-location dvi_soil for %d rows across %d locations (example locs: %s)",
                 sum(need_idx), length(bad_locs), paste(head(bad_locs, 10), collapse = ", ")))
  }
  dvi_soil_vec
}

# Calculate common spectral indices (data.table-style implementation; requires data.table when df is a data.table)
calculate_indices <- function(df) {
  eps <- 1e-9

  # keep implementation identical to existing scripts (data.table idiom)
  df[, `:=`(
    DVI   = nir - red,
    OSAVI = (nir - red) / (nir + red + 0.16),
    MCARI = ((red - green) - 0.2*(red - blue)) * (red / (green + eps)),
    NIRv  = (nir * ((nir - red) / (nir + red + eps))) * 1.3,
    PSRI  = (red - blue) / (nir + eps),
    NBR   = (nir - swir2) / (nir + swir2 + eps),
    TCW   = (swir1 - swir2) / (swir1 + swir2 + eps),
    NDVI   = (nir - red) / (nir + red + eps),
    MSAVI2 = (2 * nir + 1 - sqrt(pmax(0, (2 * nir + 1)^2 - 8 * (nir - red)))) / 2,
    MSAVI  = (2 * nir + 1 - sqrt(pmax(0, (2 * nir + 1)^2 - 8 * (nir - red)))) / 2,
    NDMI   = (nir - swir1) / (nir + swir1 + eps),
    TCB    = 0.3029 * blue + 0.2786 * green + 0.4733 * red + 0.5599 * nir + 0.508 * swir1 + 0.1872 * swir2,
    GVI    = -0.2941 * blue - 0.243 * green - 0.5424 * red + 0.7276 * nir + 0.0713 * swir1 - 0.1608 * swir2,
    SATVI  = (swir1 - red) / (swir1 + red + 0.5) * (1 + 0.5),
    EVI    = 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)
  )]

  # Add PPI if available (caller scripts may provide ppi/calculate_solar_zenith)
  if (exists("ppi") && exists("calculate_solar_zenith") && "fraction_veg" %in% names(df)) {
    dvi_soil_val <- df$DVI[df$fraction_veg == 0][1]
    if (is.finite(dvi_soil_val)) {
      zenith_rad <- calculate_solar_zenith(lat = 40, doy = 180, hour = 10.5)
      M_val <- suppressWarnings(max(df$DVI, na.rm = TRUE))
      if (is.finite(M_val)) df$PPI <- ppi(dvi = df$DVI, zenith.angle = zenith_rad, M = M_val, dvi.soil = dvi_soil_val)
    }
  }

  return(df)
}

# ------------------------------------------------------------------------------
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

  diag_count <- 0
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

      # Diagnostic for the first few attempts
      if (diag_count < 5 && exists("TESTING_MODE") && isTRUE(TESTING_MODE)) {
         cat(sprintf("[DIAGNOSTIC] Checking veg='%s', vid='%s' -> Found? %s\n", veg, vid, !is.null(vec)))
         if (is.na(vid) || vid == "NA") {
             cat(sprintf("   [DEBUG] variant class: %s, names: %s\n", class(variant), paste(names(variant), collapse=", ")))
         }
         diag_count <- diag_count + 1
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
    labs(title = "Pairwise Cosine Similarity of All Variants", x = NULL, y = NULL)

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



# Helper: add a shaded rectangle covering excluded years (1992-1999) for ggplot2 time-series plots.
# Usage: + add_excluded_years_shade(is_date = TRUE)  # x axis is Date
#        + add_excluded_years_shade(is_date = FALSE) # x axis is numeric (year)
add_excluded_years_shade <- function(start_year = 1992, end_year = 1999, is_date = FALSE, fill = "grey70", alpha = 0.35) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
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
# Draws major vertical grid lines (default: every 10 yrs). By default we DO NOT draw
# the minor (5-year) vertical lines anymore, but we still provide x-axis labels every
# 5 years so plots show ticks/labels at 5-year intervals. Set `minor_by` to a value
# to re-enable minor vlines if needed.
# Usage: + add_year_lines(is_date = TRUE)  # x axis is Date
#        + add_year_lines(is_date = FALSE) # x axis is numeric (year)
# Optional arguments:
#  minor_by, major_by - spacing for minor/major lines (integers). Set to NULL to disable.
#  minor_size, major_size - line widths for minor/major lines.
#  minor_color, major_color, minor_linetype, major_linetype, minor_alpha, major_alpha - styles.
add_year_lines <- function(start_year = 1985,
                           end_year = 2025,
                           is_date = FALSE,
                           color = "grey50",
                           linetype = "dashed",
                           alpha = 0.5,
                           minor_by = NULL,
                           major_by = 10L,
                           minor_size = 0.4,
                           major_size = 0.9,
                           minor_color = color,
                           major_color = color,
                           minor_linetype = linetype,
                           major_linetype = "solid",
                           minor_alpha = alpha,
                           major_alpha = alpha) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)

  # Generate sequences for minor/major lines; allow disabling by setting to NULL or <=0
  minor_seq <- if (!is.null(minor_by) && minor_by > 0) seq(start_year, end_year, by = as.integer(minor_by)) else integer(0)
  major_seq <- if (!is.null(major_by) && major_by > 0) seq(start_year, end_year, by = as.integer(major_by)) else integer(0)

  # Avoid drawing duplicate lines where major years overlap minor years
  minor_seq <- setdiff(minor_seq, major_seq)

  layers <- list()
  if (length(minor_seq) > 0) {
    if (is_date) {
      xints <- as.Date(paste0(minor_seq, "-01-01"))
      layers <- c(layers, list(ggplot2::geom_vline(xintercept = as.numeric(xints), color = minor_color, linetype = minor_linetype, alpha = minor_alpha, size = minor_size)))
    } else {
      layers <- c(layers, list(ggplot2::geom_vline(xintercept = minor_seq, color = minor_color, linetype = minor_linetype, alpha = minor_alpha, size = minor_size)))
    }
  }

  if (length(major_seq) > 0) {
    if (is_date) {
      xints <- as.Date(paste0(major_seq, "-01-01"))
      layers <- c(layers, list(ggplot2::geom_vline(xintercept = as.numeric(xints), color = major_color, linetype = major_linetype, alpha = major_alpha, size = major_size)))
    } else {
      layers <- c(layers, list(ggplot2::geom_vline(xintercept = major_seq, color = major_color, linetype = major_linetype, alpha = major_alpha, size = major_size)))
    }
  }

  # Add axis breaks/labels at 5-year intervals (labels only; no minor vlines).
  # Using the same visual styling as the existing major-year labels (i.e. default axis text).
  if (!is.null(start_year) && !is.null(end_year)) {
    if (is_date) {
      date_breaks <- seq(as.Date(paste0(start_year, "-01-01")), as.Date(paste0(end_year, "-01-01")), by = "5 years")
      layers <- c(layers, list(ggplot2::scale_x_date(breaks = date_breaks, date_labels = "%Y")))
    } else {
      layers <- c(layers, list(ggplot2::scale_x_continuous(breaks = seq(start_year, end_year, by = 5))))
    }
  }

  # Return list of layers so callers can add them directly to ggplot
  return(layers)
}


# Colorblind-friendly discrete palette (Okabe-Ito) — use for categorical veg colors
MESMA_OKABE_ITO <- c(
  "black" = "#000000",
  "orange" = "#E69F00",
  "skyblue" = "#56B4E9",
  "bluishgreen" = "#009E73",
  "yellow" = "#F0E442",
  "blue" = "#0072B2",
  "vermillion" = "#D55E00",
  "purple" = "#CC79A7"
)

# A consistent theme for MESMA time-series plots (subtle grid, panel border)
theme_mesma <- function(base_size = 11, base_family = "") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(ggplot2::theme_minimal(base_size = base_size))
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(fill = NA, color = "#d9d9d9", size = 0.5),
      panel.grid.major = ggplot2::element_line(color = "#efefef", size = 0.45),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = base_size * 0.95),
      axis.text = ggplot2::element_text(size = base_size * 0.85),
      plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.05),
      legend.position = "bottom",
      legend.background = ggplot2::element_rect(fill = NA, colour = NA),
      legend.key = ggplot2::element_rect(fill = NA, color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
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
# set_mesma_seed()). Defaults to 123 (was 42) for reproducibility.
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
      if (is.na(base)) base <- 123L
    } else {
      base <- 123L
    }
  }
  base <- as.integer(base)

  # Use L'Ecuyer for better parallel reproducibility and future::future.seed support
  try({ RNGkind("L'Ecuyer-CMRG") }, silent = TRUE)
  set.seed(base)

  options(mesma.seed = base)
  # keep legacy code working: many places check for a global 'seed'
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
  base <- as.integer(getOption("mesma.seed", Sys.getenv("MESMA_SEED", unset = "123")))
  if (!is.finite(base) || length(base) == 0) base <- 123L
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
HERBS_GROUP <- c("herbs", "salicornia", "halocnemum", "phragmites")

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
  # Fix common typos and canonicalize typos to "tamarix"
  v[v %in% c("tamairx", "tamarix")] <- "tamarix"
  # Map herb-group
  v[v %in% HERBS_GROUP] <- "herbs"
  # Map agriculture aliases
  v[v %in% AGRICULTURE_ALIASES] <- "agriculture"
  df[[veg_col]] <- v
  df
}

# --- Woody aggregation support removed ----------------------------------
# The old `aggregate_woody_bootstrap()` helper has been removed because the
# concept of "woody aggregation" is no longer used.  Previous callers should
# now operate on the raw bootstrap results directly.

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
# `ppi_max` is the normalizing constant (PPI_FULL_VEG_COVER).
backup_and_normalize_ppi <- function(df, ppi_max, label = "") {
  prefix <- if (nzchar(label)) paste0(label, ": ") else ""

  # Backup raw PPI
  if (!"PPI_raw" %in% names(df)) {
    df$PPI_raw <- df$PPI
    cat(sprintf("[NOTICE] %sBacked up raw PPI values to 'PPI_raw' before normalization.\n", prefix))
  }

  # Compute ppi_norm
  if (!"ppi_norm" %in% names(df)) df$ppi_norm <- NA_real_
  if ("PPI_raw" %in% names(df) && any(is.finite(df$PPI_raw))) {
    df$ppi_norm <- pmin(pmax(df$PPI_raw / ppi_max, 0), 1)
    cat(sprintf("[PPI NORM] %sCreated 'ppi_norm' from 'PPI_raw' and clamped to [0,1] using PPI_FULL_VEG_COVER=%.3f\n", prefix, ppi_max))
  } else if ("PPI" %in% names(df) && any(is.finite(df$PPI))) {
    df$ppi_norm <- pmin(pmax(df$PPI / ppi_max, 0), 1)
    warning(sprintf("%s'PPI_raw' not found - computed 'ppi_norm' from 'PPI' (may be z-scored); values were clamped to [0,1].", prefix))
  } else {
    df$ppi_norm <- NA_real_
    cat(sprintf("[PPI NORM] %sNo PPI or PPI_raw available to compute 'ppi_norm' (all NA)\n", prefix))
  }
  df
}

