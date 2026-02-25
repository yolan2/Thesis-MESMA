
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

  # --- Detect high-similarity populus/tamarix variants (cross-class similarity > 0.95) ---
  # These variants are spectrally indistinguishable and should be labeled as "woody_unknown"
  WOODY_CROSS_SIM_THRESHOLD <- 0.95
  populus_indices <- which(tolower(all_vegs) == "populus")
  tamarix_indices <- which(tolower(all_vegs) == "tamarix")

  woody_unknown_variants <- character(0)

  if (length(populus_indices) > 0 && length(tamarix_indices) > 0) {
    cat(sprintf("\n[CROSS-CLASS SIMILARITY] Checking %d populus vs %d tamarix variants (threshold: %.2f)\n",
                length(populus_indices), length(tamarix_indices), WOODY_CROSS_SIM_THRESHOLD))

    # Extract cross-class similarities from the similarity matrix
    for (p_idx in populus_indices) {
      for (t_idx in tamarix_indices) {
        sim_val <- sim_mat[p_idx, t_idx]
        if (!is.na(sim_val) && sim_val > WOODY_CROSS_SIM_THRESHOLD) {
          # Both variants are indistinguishable - add both to woody_unknown list
          p_id <- all_ids[p_idx]
          t_id <- all_ids[t_idx]
          if (!(p_id %in% woody_unknown_variants)) {
            woody_unknown_variants <- c(woody_unknown_variants, p_id)
            cat(sprintf("  - %s (populus) has similarity %.3f with %s (tamarix)\n", p_id, sim_val, t_id))
          }
          if (!(t_id %in% woody_unknown_variants)) {
            woody_unknown_variants <- c(woody_unknown_variants, t_id)
            cat(sprintf("  - %s (tamarix) has similarity %.3f with %s (populus)\n", t_id, sim_val, p_id))
          }
        }
      }
    }

    if (length(woody_unknown_variants) > 0) {
      cat(sprintf("[CROSS-CLASS SIMILARITY] Found %d variants with cross-class similarity > %.2f\n",
                  length(woody_unknown_variants), WOODY_CROSS_SIM_THRESHOLD))
      cat(sprintf("  Variants: %s\n", paste(woody_unknown_variants, collapse = ", ")))
    } else {
      cat(sprintf("[CROSS-CLASS SIMILARITY] No populus/tamarix variants exceed similarity threshold of %.2f\n",
                  WOODY_CROSS_SIM_THRESHOLD))
    }
  } else {
    cat("[CROSS-CLASS SIMILARITY] Skipping: need both populus and tamarix variants present\n")
  }

  # Store in global environment for use during MESMA fitting
  assign("WOODY_UNKNOWN_VARIANTS", woody_unknown_variants, envir = globalenv())
  cat(sprintf("[CROSS-CLASS SIMILARITY] Assigned WOODY_UNKNOWN_VARIANTS to global env (%d variants)\n",
              length(woody_unknown_variants)))

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

# --- Woody aggregation helper for bootstrap results --------------------------
# Sums bootstrap matrices for populus + tamarix + woody_unknown into a "woody" category.
# `veg_boot_res` is a named list of [B x n_years] matrices. Modifies in place and returns.
aggregate_woody_bootstrap <- function(veg_boot_res, label = "BOOTSTRAP") {
  woody_types <- c("populus", "tamarix", "woody_unknown")
  woody_mats <- veg_boot_res[tolower(names(veg_boot_res)) %in% woody_types]
  if (length(woody_mats) >= 1) {
    woody_mat <- Reduce(`+`, lapply(woody_mats, function(m) { m[is.na(m)] <- 0; m }))
    all_na_mask <- Reduce(`&`, lapply(woody_mats, is.na))
    woody_mat[all_na_mask] <- NA_real_
    colnames(woody_mat) <- colnames(woody_mats[[1]])
    veg_boot_res[["woody"]] <- woody_mat
    cat(sprintf("[%s] Added 'woody' category (populus + tamarix + woody_unknown) with combined bootstrap CIs\n", label))
  }
  veg_boot_res
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

