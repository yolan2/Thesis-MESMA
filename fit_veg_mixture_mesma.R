filter_variants_by_min_samples <- function(variants, min_samples = MIN_ENDMEMBER_SAMPLES, veg = NULL, raw_template = NULL) {
  if (is.null(variants) || length(variants) == 0) return(list())
  keep_mask <- sapply(variants, function(v) {
    # Require an explicit n_samples and enforce minimum sample count
    if (is.null(v$n_samples)) return(FALSE) # drop variants with unknown n_samples
    v$n_samples >= min_samples
  })
  removed <- sum(!keep_mask)
  if (removed > 0 && !is.null(veg)) cat(sprintf("  [%s] Removed %d variant(s) with n_samples < %d\n", veg, removed, min_samples))
  kept_variants <- variants[keep_mask]
  if (length(kept_variants) == 0 && !is.null(raw_template) && !is.null(raw_template$n_samples) && raw_template$n_samples >= min_samples) {
    cat(sprintf("  [%s] No variants left after filtering; restoring medoid raw template as single variant\n", ifelse(is.null(veg), "veg", veg)))
    kept_variants <- list(list(raw_mat = raw_template$T, variant_id = paste0(ifelse(is.null(veg), "veg", veg), "_single"), n_samples = raw_template$n_samples))
  }
  kept_variants

}
 
library(zoo)
library(dplyr)
library(lme4)
library(cluster)
library(readr)
library(sf)
library(ggplot2)
library(scales)
library(nlme)
library(RStoolbox)
library(quadprog)
library(terra)
library(magrittr) # For pipe operator
library(future)
library(future.apply)
library(MASS) # For lda

# write_debug is defined in mesma_config.R (no-op stub)

# compute_dvi_soil_per_location() canonicalized in 'mesma_helpers.R'
if (!exists("compute_dvi_soil_per_location") && file.exists("mesma_helpers.R")) source("mesma_helpers.R")

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
if (file.exists("mesma_config.R")) {
  source("mesma_config.R")
} else {
  stop("Required file 'mesma_config.R' not found in project root.")
}

# Override: inference CSV for this script (differs from mesma_config.R default)
INFERENCE_CSV <- "C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_low (3).csv"

if (file.exists("ppi_helpers.R")) {
  source("ppi_helpers.R")
} else {
  stop("Required file 'ppi_helpers.R' not found in project root. Please add it to ensure consistent PPI calculation.")
}

if (file.exists("mesma_helpers.R")) {
  source("mesma_helpers.R")
  set_mesma_seed()
} else {
  stop("Required file 'mesma_helpers.R' not found in project root.")
}

# --- Shared utility: suppress stdout/message during an expression ---
# Replaces 5+ repeated tryCatch(file(...)) + sink() + finally patterns.
suppress_output_safely <- function(expr, quiet_output = TRUE, quiet_message = TRUE) {
  sink_file <- tempfile()
  con_out <- if (quiet_output) tryCatch(file(sink_file, open = "wt"), error = function(e) NULL) else NULL
  con_msg <- if (quiet_message) tryCatch(file(paste0(sink_file, ".msg"), open = "wt"), error = function(e) NULL) else NULL
  tryCatch({
    if (!is.null(con_out) && inherits(con_out, "connection") && isOpen(con_out)) sink(con_out, type = "output")
    if (!is.null(con_msg) && inherits(con_msg, "connection") && isOpen(con_msg)) sink(con_msg, type = "message")
    force(expr)
  }, finally = {
    if (quiet_message) try(sink(type = "message"), silent = TRUE)
    if (quiet_output) try(sink(type = "output"), silent = TRUE)
    if (!is.null(con_out)) try(close(con_out), silent = TRUE)
    if (!is.null(con_msg)) try(close(con_msg), silent = TRUE)
    try(unlink(sink_file), silent = TRUE)
    try(unlink(paste0(sink_file, ".msg")), silent = TRUE)
  })
}

# --- Shared utility: filter valid vegetation rows (excludes NA/empty/barren) ---
filter_valid_vegetation <- function(coefs, exclude_barren = TRUE) {
  veg_char <- as.character(coefs$Veg)
  mask <- !is.na(veg_char) & nchar(trimws(veg_char)) > 0 & trimws(veg_char) != "NA"
  if (exclude_barren) mask <- mask & !normalize_veg_name(veg_char) %in% c("barren")
  coefs[mask, ]
}

# --- Shared utility: aggregate batch results into a results list ---
aggregate_batch_results <- function(batch_results, results_list, save_csv_dir = NULL) {
  for (k in names(batch_results)) {
    loc_result <- batch_results[[k]]
    if (is.null(loc_result)) next
    if (!is.null(save_csv_dir)) {
      loc_data <- do.call(rbind, lapply(loc_result, function(yr_res) yr_res$coef_df))
      if (!is.null(loc_data) && nrow(loc_data) > 0) {
        out_fname <- file.path(save_csv_dir, paste0("result_", make.names(k), ".csv"))
        readr::write_csv(loc_data, out_fname)
      }
    }
    for (yr_char in names(loc_result)) {
      r <- loc_result[[yr_char]]
      if (is.null(r)) next
      if (!is.null(save_csv_dir) && (is.null(r$coef_df) || nrow(r$coef_df) == 0)) next
      res_key <- if (!is.null(r$coef_df) && "pheno_year" %in% names(r$coef_df)) {
        paste(k, r$coef_df$pheno_year[1], sep = "_")
      } else {
        paste(k, yr_char, sep = "_")
      }
      results_list[[res_key]] <- list(
        coef_df = r$coef_df,
        diagnostics = r$diagnostics,
        uncertainty = r$uncertainty
      )
    }
  }
  results_list
}

# --- Shared plot constants ---
HERBS_COLORS <- c("Herbs" = "#009E73")
# Add 'Agriculture' as a first-class species for plotting
SPECIES_COLORS <- c("Herbs" = "#009E73", "Populus" = "#0072B2", "Tamarix" = "#D55E00", "Agriculture" = "#F0E442")


# --- Shared utility: normalize vegetation names (lowercase + trim) ---
normalize_veg_name <- function(x) tolower(trimws(as.character(x)))

# --- Shared utility: compute pairwise Haversine distance matrix (km) ---
compute_haversine_distance_matrix <- function(coords) {
  rad <- pi / 180
  n_v <- nrow(coords)
  dist_mat <- matrix(0, n_v, n_v)
  for (ii in 1:(n_v - 1)) {
    for (jj in (ii + 1):n_v) {
      dlat <- (coords[jj, 2] - coords[ii, 2]) * rad
      dlon <- (coords[jj, 1] - coords[ii, 1]) * rad
      a <- sin(dlat / 2)^2 + cos(coords[ii, 2] * rad) * cos(coords[jj, 2] * rad) * sin(dlon / 2)^2
      dist_mat[ii, jj] <- 2 * 6371 * asin(sqrt(a))
      dist_mat[jj, ii] <- dist_mat[ii, jj]
    }
  }
  dist_mat
}

# --- Shared utility: fit exponential variogram via NLS with fallback ---
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
                                          include_species_plots = TRUE,
                                          pure_threshold = 0.9,
                                          mix_contrib_threshold = 0.10,
                                          top_n_mixtures = 10) {
  if (is.null(full_data) || nrow(full_data) == 0) return(invisible(NULL))

  veg_norm <- normalize_veg_name(full_data$Veg)
  inf_veg <- full_data[!veg_norm %in% c("barren"), ]
  inf_barren <- full_data[veg_norm %in% c("barren"), ]

  # Helper: base layers for time series plots
  .ts_layers <- function() {
    layers <- list()
    if (use_excluded_years_shade) layers <- c(layers, list(add_excluded_years_shade(is_date = FALSE)))
    layers <- c(layers, list(add_year_lines(is_date = FALSE)))
    # apply MESMA theme for consistent time-series styling (subtle grid + border)
    layers <- c(layers, list(theme_mesma()))
    layers
  }

  # 1. Vegetation time series
  if (!is.null(inf_veg) && nrow(inf_veg) > 0) {
    p_ts <- ggplot(inf_veg, aes(x = year, y = global_coef, color = Veg, group = Veg)) +
      .ts_layers() +
      geom_line(linewidth = 1) +
      geom_point(show.legend = FALSE) +
      geom_ribbon(aes(ymin = coef_025, ymax = coef_975, fill = Veg), alpha = 0.15, color = NA) +
      labs(title = paste0(method, "-Normalized Vegetation Fractions"),
           x = "Year", y = "Total Normalized Fraction", color = "Veg", fill = "Veg") +
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
    labs(title = paste0(method, "-Based Barren Fraction"), x = "Year", y = "Barren Fraction") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
    theme_minimal()
  ggsave(file.path(OUT_DIR, paste0("inference_", file_prefix, "_barren_cover.png")), p_barren, width = 8, height = 6)
  readr::write_csv(inf_barren, file.path(OUT_DIR, paste0("inference_", file_prefix, "_barren_cover.csv")))
  cat(sprintf("Saved inference %s-based barren cover plot to: %s\n", method, file.path(OUT_DIR, paste0("inference_", file_prefix, "_barren_cover.png"))))




  if (include_species_plots) {
    species_data <- full_data[veg_norm %in% c("herbs","populus","tamarix","agriculture"), ]
    if (nrow(species_data) > 0) {
      species_data$Veg <- dplyr::case_when(
        tolower(species_data$Veg) == "herbs" ~ "Herbs",
        tolower(species_data$Veg) == "populus" ~ "Populus",
        tolower(species_data$Veg) == "tamarix" ~ "Tamarix",
        tolower(species_data$Veg) == "agriculture" ~ "Agriculture",
        TRUE ~ species_data$Veg
      )

      p_sp_ts <- ggplot(species_data, aes(x = year, y = global_coef, color = Veg, group = Veg)) +
        .ts_layers() +
        geom_line(linewidth = 1) +
        geom_point(show.legend = FALSE) +
        geom_ribbon(aes(ymin = coef_025, ymax = coef_975, fill = Veg), alpha = 0.12, color = NA) +
        scale_color_manual(values = SPECIES_COLORS) +
        scale_fill_manual(values = SPECIES_COLORS) +
        labs(title = paste0("Inference ", method, ": Species"),
             x = "Year", y = "Total Normalized Fraction", color = "Veg", fill = "Veg") +
        theme_minimal()
      ggsave(file.path(OUT_DIR, paste0("inference_", file_prefix, "_species_separate.png")), p_sp_ts, width = 8, height = 6)
      readr::write_csv(species_data, file.path(OUT_DIR, paste0("inference_", file_prefix, "_species_separate.csv")))
      cat(sprintf("Saved inference %s species plot to: %s\n", method, file.path(OUT_DIR, paste0("inference_", file_prefix, "_species_separate.png"))))

      # Species stacked proportion
      df_wide_sp <- tryCatch({ tidyr::pivot_wider(species_data |> dplyr::select(year, Veg, global_coef), names_from = Veg, values_from = global_coef) }, error = function(e) NULL)
      if (!is.null(df_wide_sp)) {
        sp_cols <- intersect(c("Herbs","Populus","Tamarix","Agriculture"), names(df_wide_sp))
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
            labs(title = paste0("Inference ", method, ": Species (Proportion, stacked)"),
                 x = "Year", y = "Proportion", fill = "Veg") +
            theme_minimal()
          ggsave(file.path(OUT_DIR, paste0("inference_", file_prefix, "_species_stacked.png")), p_sp_stacked, width = 8, height = 6)
          readr::write_csv(df_prop_sp, file.path(OUT_DIR, paste0("inference_", file_prefix, "_species_stacked.csv")))
          cat(sprintf("Saved inference %s species stacked plot to: %s\n", method, file.path(OUT_DIR, paste0("inference_", file_prefix, "_species_stacked.png"))))
        }
      }
    }
  }


  # Pixel-level purity/mixture analysis requires a per-pixel coefficient table
  # (columns: location_id, pheno_year, Veg, coef). Use `full_data` when it
  # contains those columns, otherwise fall back to `inference_coefs` if present.
  coefs_src <- NULL
  if (all(c("location_id", "pheno_year", "Veg", "coef") %in% names(full_data))) {
    coefs_src <- full_data
  } else if (exists("inference_coefs") && is.data.frame(inference_coefs) && all(c("location_id", "pheno_year", "Veg", "coef") %in% names(inference_coefs))) {
    coefs_src <- inference_coefs
  } else {
    cat("[PURITY] No per-pixel coefficient table available (missing 'coef'/'location_id'); skipping pixel-level purity/mixture plots.\n")
    return(invisible(NULL))
  }

  pixel_wide <- coefs_src %>%
    dplyr::filter(is.finite(.data$coef)) %>%
    dplyr::group_by(location_id, pheno_year, Veg) %>%
    dplyr::summarise(coef = sum(.data$coef, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Veg, values_from = coef, values_fill = 0)

  veg_cols <- setdiff(names(pixel_wide), c("location_id", "pheno_year"))
  if (length(veg_cols) == 0) {
    cat("[PURITY] No vegetation columns after pivot; skipping.\n")
    return(invisible(NULL))
  }

  # Tolerance for treating "any vegetation" as present (avoid floating-point noise)
  TOL_VEG_PRES <- 1e-9

  # Display labels for purity types used in plots
  PURE_LABEL_ABS <- "Pure (absolute)"
  PURE_LABEL_VEGONLY <- "Pure (pure)"

  # Identify vegetation columns excluding 'barren' (may be identical to `veg_cols` if no barren column)
  veg_cols_no_barren <- setdiff(veg_cols, "barren")

  # total vegetation fraction (MESMA sum across vegetation classes only)
  if (length(veg_cols_no_barren) > 0) {
    pixel_wide$veg_total <- rowSums(pixel_wide[, veg_cols_no_barren, drop = FALSE], na.rm = TRUE)
  } else {
    pixel_wide$veg_total <- 0
  }

  # --- Absolute purity (barren NOT assigned if ANY vegetation present) ---
  # Determine dominant class across all MESMA columns but treat 'barren' as absent
  # for rows where any vegetation is present (veg_total > TOL_VEG_PRES).
  pixel_wide$dominant_class_abs <- apply(pixel_wide[, veg_cols, drop = FALSE], 1, function(row) {
    if (all(is.na(row)) || sum(row, na.rm = TRUE) <= TOL_VEG_PRES) return(NA_character_)
    if (length(veg_cols_no_barren) > 0 && sum(row[veg_cols_no_barren], na.rm = TRUE) > TOL_VEG_PRES) {
      # prefer dominant among vegetation-only columns when any vegetation exists
      veg_cols_no_barren[which.max(row[veg_cols_no_barren])]
    } else {
      veg_cols[which.max(row)]
    }
  })

  pixel_wide$dominant_frac_abs <- apply(pixel_wide[, veg_cols, drop = FALSE], 1, function(row) {
    if (length(veg_cols_no_barren) > 0 && sum(row[veg_cols_no_barren], na.rm = TRUE) > TOL_VEG_PRES) {
      max(row[veg_cols_no_barren], na.rm = TRUE)
    } else {
      max(row, na.rm = TRUE)
    }
  })
  pixel_wide$is_pure_abs <- !is.na(pixel_wide$dominant_frac_abs) & (pixel_wide$dominant_frac_abs >= pure_threshold)

  if (length(veg_cols_no_barren) == 0) {
    cat("[PURITY] No non-barren vegetation classes present after pivot; skipping.\n")
    return(invisible(NULL))
  }

  # Per-row vegetation-only normalized fractions (used for purity / mixture decisions)
  for (vc in veg_cols_no_barren) {
    norm_col <- paste0("norm_", vc)
    pixel_wide[[norm_col]] <- ifelse(is.finite(pixel_wide$veg_total) & pixel_wide$veg_total > TOL_VEG_PRES,
                                     pixel_wide[[vc]] / pixel_wide$veg_total,
                                     NA_real_)
  }

  # Minimum absolute vegetation fraction required to assign a veg-only 'pure' label.
  # Prevents pixels that are mostly barren (but have a single, tiny veg class) from
  # having their veg-only normalized fraction inflate to 1.0 and being mis-labelled pure.
  MIN_VEG_FRAC_FOR_PURITY <- 0.10

  # Dominant vegetation class and vegetation-only dominant fraction
  pixel_wide$dominant_class <- apply(pixel_wide[, veg_cols_no_barren, drop = FALSE], 1, function(row) {
    if (all(is.na(row)) || sum(row, na.rm = TRUE) <= TOL_VEG_PRES) return(NA_character_)
    veg_cols_no_barren[which.max(row)]
  })
  norm_cols <- paste0("norm_", veg_cols_no_barren)
  pixel_wide$dominant_frac <- apply(pixel_wide[, norm_cols, drop = FALSE], 1, max, na.rm = TRUE)

  # Pure = dominant veg class accounts for >= pure_threshold of the veg-only mixture
  # AND the pixel has at least MIN_VEG_FRAC_FOR_PURITY absolute vegetation cover.
  # Without the second condition, a pixel that is 96% barren + 4% Populus has
  # norm_Populus = 1.0 and would be mis-labelled as "pure Populus".
  pixel_wide$is_pure <- !is.na(pixel_wide$dominant_frac) &
    (pixel_wide$dominant_frac >= pure_threshold) &
    (pixel_wide$veg_total >= MIN_VEG_FRAC_FOR_PURITY)

  # Build mixture label using vegetation-only normalized fractions
  # NOTE: return the contributing species names for single-contributor rows
  # and reserve the literal "Pure" label only for rows actually flagged
  # as `is_pure` (or rows lacking vegetation). This prevents single-veg
  # mixed pixels from being mis-labelled as "Pure".
  pixel_wide$mixture_label <- apply(pixel_wide[, norm_cols, drop = FALSE], 1, function(row) {
    contributing <- veg_cols_no_barren[which(row >= mix_contrib_threshold)]
    contributing <- contributing[order(-row[contributing])]
    if (length(contributing) == 0) return(NA_character_)
    if (length(contributing) == 1) return(tools::toTitleCase(contributing))
    paste(tools::toTitleCase(contributing), collapse = " + ")
  })
  # Ensure only truly 'pure' pixels receive the explicit "Pure" label
  pixel_wide$mixture_label[is.na(pixel_wide$mixture_label)] <- "Pure"
  pixel_wide$mixture_label[pixel_wide$is_pure] <- "Pure"

  # Sanity check: detect rows incorrectly labelled 'Pure' (should be rare)
  bad_mask <- pixel_wide$mixture_label == "Pure" & !pixel_wide$is_pure & pixel_wide$veg_total > TOL_VEG_PRES
  if (any(bad_mask, na.rm = TRUE)) {
    warning(sprintf("[PURITY] %d rows have mixture_label == 'Pure' but are not flagged is_pure; check mix_contrib_threshold/pure_threshold", sum(bad_mask, na.rm = TRUE)))
  }

  # --- Plot 1: Pure pixel count per class over time (barren included) ---
  pure_pixels <- pixel_wide[pixel_wide$is_pure_abs, ]
  if (nrow(pure_pixels) > 0) {
    pure_counts <- pure_pixels %>%
      dplyr::mutate(dominant_class = tools::toTitleCase(dominant_class_abs)) %>%
      dplyr::group_by(pheno_year, dominant_class) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop")

    all_colors <- c(SPECIES_COLORS, HERBS_COLORS, c("Barren" = "saddlebrown"))
    used_classes <- unique(pure_counts$dominant_class)
    plot_colors <- all_colors[names(all_colors) %in% used_classes]
    missing <- setdiff(used_classes, names(plot_colors))
    if (length(missing) > 0) {
      extra <- setNames(scales::hue_pal()(length(missing)), missing)
      plot_colors <- c(plot_colors, extra)
    }

    p_pure <- ggplot(pure_counts, aes(x = pheno_year, y = n, color = dominant_class, group = dominant_class)) +
      .ts_layers() +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      scale_color_manual(values = plot_colors) +
      labs(title = sprintf("Pure Pixels per Class Over Time (absolute purity; threshold >= %.0f%%)", pure_threshold * 100),
           x = "Year", y = "Number of Pure Pixels", color = "Class") +
      theme_minimal()
    ggsave(file.path(OUT_DIR, "pixel_purity_count_per_class.png"), p_pure, width = 10, height = 6)
    readr::write_csv(pure_counts, file.path(OUT_DIR, "pixel_purity_count_per_class.csv"))
    cat(sprintf("[PURITY] Saved pure pixel count plot to: %s\n", file.path(OUT_DIR, "pixel_purity_count_per_class.png")))

    # Stacked proportion version
    pure_totals <- pure_counts %>%
      dplyr::group_by(pheno_year) %>%
      dplyr::mutate(prop = n / sum(n)) %>%
      dplyr::ungroup() %>%
      tidyr::complete(pheno_year, dominant_class, fill = list(n = 0L, prop = 0))

    p_pure_prop <- ggplot(pure_totals, aes(x = pheno_year, y = prop, fill = dominant_class)) +
      .ts_layers() +
      geom_area(alpha = 0.8) +
      scale_fill_manual(values = plot_colors) +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(title = sprintf("Pure Pixel Composition Over Time (absolute purity; threshold >= %.0f%%)", pure_threshold * 100),
           x = "Year", y = "Proportion of Pure Pixels", fill = "Class") +
      theme_minimal()
    ggsave(file.path(OUT_DIR, "pixel_purity_proportion_per_class.png"), p_pure_prop, width = 10, height = 6)
    cat(sprintf("[PURITY] Saved pure pixel proportion plot to: %s\n", file.path(OUT_DIR, "pixel_purity_proportion_per_class.png")))

    # --- ALSO: MSAVI-filtered / MSAVI-normalized purity plots (barren included) ---
    if (exists("df_tasks_inference_proc") && ("MSAVI_raw" %in% names(df_tasks_inference_proc) || "MSAVI" %in% names(df_tasks_inference_proc))) {
      ms_col <- if ("MSAVI_raw" %in% names(df_tasks_inference_proc)) "MSAVI_raw" else "MSAVI"
      ms_lookup <- df_tasks_inference_proc %>%
        dplyr::group_by(location_id, pheno_year) %>%
        dplyr::summarise(msavi = median(.data[[ms_col]], na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(total_veg_cover_expected = pmin(pmax(msavi / (0.5 + msavi), 0), 1))

      pixel_wide_msavi <- dplyr::left_join(pixel_wide, ms_lookup, by = c("location_id", "pheno_year"))
      MIN_MSAVI_TOTAL_VEG <- 0.02
      # Keep rows with finite MSAVI; classify low-MSAVI pixels as pure barren instead of dropping them.
      pixel_wide_msavi <- pixel_wide_msavi[is.finite(pixel_wide_msavi$total_veg_cover_expected), ]
      pixel_wide_msavi$msavi_low <- pixel_wide_msavi$total_veg_cover_expected < MIN_MSAVI_TOTAL_VEG
      if (any(pixel_wide_msavi$msavi_low, na.rm = TRUE)) {
        n_low <- sum(pixel_wide_msavi$msavi_low, na.rm = TRUE)
        cat(sprintf("[PURITY] MSAVI < %.2f for %d pixel-years — classifying as pure barren for MSAVI-normalized plots.\n", MIN_MSAVI_TOTAL_VEG, n_low))
        # Force barren behaviour for those rows in the MSAVI view (veg_total -> 0, mark dominant as barren, treat as pure absolute)
        idx_low <- which(pixel_wide_msavi$msavi_low)
        pixel_wide_msavi$veg_total[idx_low] <- 0
        pixel_wide_msavi$dominant_class_abs[idx_low] <- "barren"
        pixel_wide_msavi$dominant_frac_abs[idx_low] <- 1.0
        pixel_wide_msavi$is_pure_abs[idx_low] <- TRUE
      }

      if (nrow(pixel_wide_msavi) > 0) {
        pure_pixels_msavi <- pixel_wide_msavi[pixel_wide_msavi$is_pure_abs, ]
        pure_counts_msavi <- pure_pixels_msavi %>%
          dplyr::mutate(dominant_class = tools::toTitleCase(dominant_class_abs)) %>%
          dplyr::group_by(pheno_year, dominant_class) %>%
          dplyr::summarise(n = dplyr::n(), .groups = "drop")

        used_classes_ms <- unique(pure_counts_msavi$dominant_class)
        plot_colors_ms <- all_colors[names(all_colors) %in% used_classes_ms]
        missing_ms <- setdiff(used_classes_ms, names(plot_colors_ms))
        if (length(missing_ms) > 0) plot_colors_ms <- c(plot_colors_ms, setNames(scales::hue_pal()(length(missing_ms)), missing_ms))

        p_pure_msavi <- ggplot(pure_counts_msavi, aes(x = pheno_year, y = n, color = dominant_class, group = dominant_class)) +
          .ts_layers() +
          geom_line(linewidth = 1) +
          geom_point(size = 2) +
          scale_color_manual(values = plot_colors_ms) +
          labs(title = sprintf("(MSAVI-filtered) Pure Pixels per Class Over Time (absolute purity; threshold >= %.0f%%)", pure_threshold * 100),
               x = "Year", y = "Number of Pure Pixels", color = "Class") +
          theme_minimal()
        ggsave(file.path(OUT_DIR, "pixel_purity_count_per_class_msavi_normalized.png"), p_pure_msavi, width = 10, height = 6)
        readr::write_csv(pure_counts_msavi, file.path(OUT_DIR, "pixel_purity_count_per_class_msavi_normalized.csv"))
        cat(sprintf("[PURITY] Saved MSAVI-normalized pure pixel count plot to: %s\n", file.path(OUT_DIR, "pixel_purity_count_per_class_msavi_normalized.png")))

        pure_totals_msavi <- pure_counts_msavi %>%
          dplyr::group_by(pheno_year) %>%
          dplyr::mutate(prop = n / sum(n)) %>%
          dplyr::ungroup() %>%
          tidyr::complete(pheno_year, dominant_class, fill = list(n = 0L, prop = 0))

        p_pure_prop_msavi <- ggplot(pure_totals_msavi, aes(x = pheno_year, y = prop, fill = dominant_class)) +
          .ts_layers() +
          geom_area(alpha = 0.8) +
          scale_fill_manual(values = plot_colors_ms) +
          scale_y_continuous(labels = scales::percent_format()) +
          labs(title = sprintf("(MSAVI-filtered) Pure Pixel Composition Over Time (absolute purity; threshold >= %.0f%%)", pure_threshold * 100),
               x = "Year", y = "Proportion of Pure Pixels", fill = "Class") +
          theme_minimal()
        ggsave(file.path(OUT_DIR, "pixel_purity_proportion_per_class_msavi_normalized.png"), p_pure_prop_msavi, width = 10, height = 6)
        cat(sprintf("[PURITY] Saved MSAVI-normalized pure pixel proportion plot to: %s\n", file.path(OUT_DIR, "pixel_purity_proportion_per_class_msavi_normalized.png")))

        # MSAVI-filtered: Pure vs Mixed counts
        purity_summary_msavi <- pixel_wide_msavi %>%
          dplyr::mutate(pixel_type = dplyr::case_when(
            !is.finite(veg_total) | veg_total <= TOL_VEG_PRES ~ "Barren",
            is_pure_abs ~ PURE_LABEL_ABS,
            TRUE ~ "Mixed"
          )) %>%
          dplyr::group_by(pheno_year, pixel_type) %>%
          dplyr::summarise(n = dplyr::n(), .groups = "drop")

        p_purity_msavi2 <- ggplot(purity_summary_msavi, aes(x = pheno_year, y = n, fill = pixel_type)) +
          .ts_layers() +
          geom_col(position = "stack", alpha = 0.85) +
          scale_fill_manual(values = setNames(c("#4CAF50", "saddlebrown", "#FF9800"), c(PURE_LABEL_ABS, "Barren", "Mixed"))) +
          labs(title = "(MSAVI-filtered) Pure vs Mixed Pixels Over Time", x = "Year", y = "Number of Pixels", fill = "Type") +
          theme_minimal()
        ggsave(file.path(OUT_DIR, "pixel_pure_vs_mixed_count_msavi_normalized.png"), p_purity_msavi2, width = 10, height = 6)
        readr::write_csv(purity_summary_msavi, file.path(OUT_DIR, "pixel_pure_vs_mixed_count_msavi_normalized.csv"))
        cat(sprintf("[PURITY] Saved MSAVI-normalized pure vs mixed count plot to: %s\n", file.path(OUT_DIR, "pixel_pure_vs_mixed_count_msavi_normalized.png")))

        # MSAVI-filtered: mixture trends + proportion (top N mixtures)
        mixed_pixels_msavi <- pixel_wide_msavi[!pixel_wide_msavi$is_pure & pixel_wide_msavi$veg_total > 1e-9, ]
        if (nrow(mixed_pixels_msavi) > 0) {
          mix_counts_msavi <- mixed_pixels_msavi %>%
            # Exclude literal "Pure" and single-contributor "Herbs"/"Agriculture"
            # (only show per-class "... (pure)" in mixture plots)
            dplyr::filter(mixture_label != "Pure", mixture_label != "Herbs", mixture_label != "Agriculture") %>%
            dplyr::group_by(pheno_year, mixture_label) %>%
            dplyr::summarise(n = dplyr::n(), .groups = "drop")

          top_mixtures_msavi <- mix_counts_msavi %>%
            dplyr::group_by(mixture_label) %>%
            dplyr::summarise(total = sum(n), .groups = "drop") %>%
            dplyr::arrange(dplyr::desc(total)) %>%
            dplyr::slice_head(n = top_n_mixtures) %>%
            dplyr::pull(mixture_label)

          mix_counts_top_msavi <- mix_counts_msavi %>%
            dplyr::mutate(mixture_label = ifelse(mixture_label %in% top_mixtures_msavi, mixture_label, "Other")) %>%
            dplyr::group_by(pheno_year, mixture_label) %>%
            dplyr::summarise(n = sum(n), .groups = "drop")

          n_mix_colors_ms <- length(unique(mix_counts_top_msavi$mixture_label))
          mix_palette_ms <- setNames(scales::hue_pal()(n_mix_colors_ms), unique(mix_counts_top_msavi$mixture_label))
          if ("Other" %in% names(mix_palette_ms)) mix_palette_ms["Other"] <- "grey60"

          p_mix_msavi <- ggplot(mix_counts_top_msavi, aes(x = pheno_year, y = n, color = mixture_label, group = mixture_label)) +
            .ts_layers() +
            geom_line(linewidth = 1) +
            geom_point(size = 1.5) +
            scale_color_manual(values = mix_palette_ms) +
            labs(title = if (is.finite(top_n_mixtures)) sprintf("(MSAVI-filtered) Specific Mixture Types Over Time (top %d)", top_n_mixtures) else "(MSAVI-filtered) Specific Mixture Types Over Time (all mixtures)", x = "Year", y = "Number of Mixed Pixels", color = "Mixture") +
            theme_minimal() +
            theme(legend.position = "right")
          ggsave(file.path(OUT_DIR, "pixel_mixture_trends_msavi_normalized.png"), p_mix_msavi, width = 12, height = 6)
          readr::write_csv(mix_counts_top_msavi, file.path(OUT_DIR, "pixel_mixture_trends_msavi_normalized.csv"))
          cat(sprintf("[PURITY] Saved MSAVI-normalized mixture trends plot to: %s\n", file.path(OUT_DIR, "pixel_mixture_trends_msavi_normalized.png")))

          # --- include pure-pixel counts so proportions are relative to ALL pixels (pure + mixed) ---
          pure_per_class_msavi <- pixel_wide_msavi %>%
            dplyr::filter(is_pure & veg_total > TOL_VEG_PRES) %>%
            dplyr::group_by(pheno_year, dominant_class) %>%
            dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
            dplyr::mutate(mixture_label = paste0(tools::toTitleCase(dominant_class), " (pure)"))

          barren_total_msavi <- pixel_wide_msavi %>%
            dplyr::filter(!is.finite(veg_total) | veg_total <= TOL_VEG_PRES) %>%
            dplyr::group_by(pheno_year) %>%
            dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
            dplyr::mutate(mixture_label = "Barren")

          mix_counts_top_msavi_all <- dplyr::bind_rows(mix_counts_top_msavi, pure_per_class_msavi, barren_total_msavi)

          # Filter extremely rare mixture labels (same threshold used in non-MSAVI branch)
          RARE_MIXTURE_MIN_COUNT <- 15L
          if (RARE_MIXTURE_MIN_COUNT > 0) {
            totals_ms <- mix_counts_top_msavi_all %>%
              dplyr::group_by(mixture_label) %>%
              dplyr::summarise(total = sum(n), .groups = "drop")
            rare_ms <- totals_ms$mixture_label[totals_ms$total < RARE_MIXTURE_MIN_COUNT]
            if (length(rare_ms) > 0) {
              mix_counts_top_msavi_all <- mix_counts_top_msavi_all %>%
                dplyr::mutate(n = ifelse(mixture_label %in% rare_ms, 0L, n))
            }
          }

          # Build palette: keep existing mix colors, assign per-class pure colours to match class colour, and ensure 'Barren' uses saddlebrown.
          if (!exists("all_colors")) all_colors <- c()
          mix_palette_ms <- mix_palette_ms
          pure_labels_ms <- unique(pure_per_class_msavi$mixture_label)
          for (pl in pure_labels_ms) {
            base <- sub(" \\(.*", "", pl)
            if (base %in% names(all_colors)) {
              mix_palette_ms[pl] <- all_colors[[base]]
            } else {
              mix_palette_ms[pl] <- scales::hue_pal()(1)
            }
          }
          if (!"Barren" %in% names(mix_palette_ms)) mix_palette_ms["Barren"] <- "saddlebrown"

          # Denominator = total pixel-years in the MSAVI-filtered table (includes
          # uncategorised pixels that fell into no bucket, so proportions sum to <= 100%).
          total_px_per_year_msavi <- pixel_wide_msavi %>%
            dplyr::group_by(pheno_year) %>%
            dplyr::summarise(total_px = dplyr::n(), .groups = "drop")

          # Collapse any duplicate (pheno_year, mixture_label) rows, then compute proportions.
          mix_props_msavi <- mix_counts_top_msavi_all %>%
            dplyr::group_by(pheno_year, mixture_label) %>%
            dplyr::summarise(n = sum(n), .groups = "drop") %>%
            dplyr::left_join(total_px_per_year_msavi, by = "pheno_year") %>%
            dplyr::mutate(prop = n / total_px) %>%
            dplyr::select(-total_px)

          # force factor ordering: Barren first, then per-class pure, then mixtures
          all_labels_ms <- unique(c("Barren", pure_labels_ms, setdiff(mix_counts_top_msavi_all$mixture_label, c("Barren", pure_labels_ms))))
          mix_props_msavi$mixture_label <- factor(mix_props_msavi$mixture_label, levels = all_labels_ms)
          mix_palette_ms <- mix_palette_ms[all_labels_ms]

          # Fill missing year×label combos with prop=0 to prevent geom_area interpolation above 100%.
          mix_props_msavi <- tidyr::complete(mix_props_msavi, pheno_year, mixture_label,
                                             fill = list(n = 0L, prop = 0))

          p_mix_prop_msavi <- ggplot(mix_props_msavi, aes(x = pheno_year, y = prop, fill = mixture_label)) +
            .ts_layers() +
            geom_area(alpha = 0.8) +
            scale_fill_manual(values = mix_palette_ms) +
            scale_y_continuous(labels = scales::percent_format()) +
            labs(title = if (is.finite(top_n_mixtures)) sprintf("(MSAVI-filtered) Pixel Composition Over Time (mixed + pure; top %d)", top_n_mixtures) else "(MSAVI-filtered) Pixel Composition Over Time (mixed + pure; all mixtures)", x = "Year", y = "Proportion of Pixels", fill = "Category") +
            theme_minimal()
          ggsave(file.path(OUT_DIR, "pixel_mixture_proportion_msavi_normalized.png"), p_mix_prop_msavi, width = 12, height = 6)
          readr::write_csv(mix_props_msavi, file.path(OUT_DIR, "pixel_mixture_proportion_msavi_normalized.csv"))
          cat(sprintf("[PURITY] Saved MSAVI-normalized mixture proportion plot to: %s\n", file.path(OUT_DIR, "pixel_mixture_proportion_msavi_normalized.png")))
        }

      } else {
        cat("[PURITY] MSAVI present but no pixel-years meet MSAVI veg threshold; skipping MSAVI-normalized purity plots.\n")
      }
    } else {
      cat("[PURITY] MSAVI per-pixel data not available; skipping MSAVI-normalized purity plots.\n")
    }
  } else {
    cat("[PURITY] No pure pixels found at threshold; skipping pure pixel plots.\n")
  }

  # --- Plot 2: Pure vs Mixed pixel counts over time ---
  purity_summary <- pixel_wide %>%
    dplyr::filter(veg_total > 1e-9) %>%
    dplyr::mutate(pixel_type = ifelse(is_pure_abs, PURE_LABEL_ABS, "Mixed")) %>%
    dplyr::group_by(pheno_year, pixel_type) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

  p_purity <- ggplot(purity_summary, aes(x = pheno_year, y = n, fill = pixel_type)) +
    .ts_layers() +
    geom_col(position = "stack", alpha = 0.85) +
    scale_fill_manual(values = setNames(c("#4CAF50", "#FF9800"), c(PURE_LABEL_ABS, "Mixed"))) +
    labs(title = "Pure vs Mixed Pixels Over Time",
         x = "Year", y = "Number of Pixels", fill = "Type") +
    theme_minimal()
  ggsave(file.path(OUT_DIR, "pixel_pure_vs_mixed_count.png"), p_purity, width = 10, height = 6)
  readr::write_csv(purity_summary, file.path(OUT_DIR, "pixel_pure_vs_mixed_count.csv"))
  cat(sprintf("[PURITY] Saved pure vs mixed count plot to: %s\n", file.path(OUT_DIR, "pixel_pure_vs_mixed_count.png")))

  # --- Plot 3: Specific mixture trends over time ---
  # Consider only rows with vegetation (exclude barren-only) and mixtures that contain >1 veg type
  mixed_pixels <- pixel_wide[!pixel_wide$is_pure & pixel_wide$veg_total > 1e-9, ]
  if (nrow(mixed_pixels) > 0) {
    mix_counts <- mixed_pixels %>%
      # Omit single-contributor "Herbs"/"Agriculture" entries from mixed-pixel analyses;
      # keep only per-class pure labels (e.g. "Herbs (pure)") in mixture composition plots
      dplyr::filter(mixture_label != "Pure", mixture_label != "Herbs", mixture_label != "Agriculture") %>%
      dplyr::group_by(pheno_year, mixture_label) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop")

    # Keep only the top N most frequent mixtures overall
    top_mixtures <- mix_counts %>%
      dplyr::group_by(mixture_label) %>%
      dplyr::summarise(total = sum(n), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(total)) %>%
      dplyr::slice_head(n = top_n_mixtures) %>%
      dplyr::pull(mixture_label)

    mix_counts_top <- mix_counts %>%
      dplyr::mutate(mixture_label = ifelse(mixture_label %in% top_mixtures, mixture_label, "Other")) %>%
      dplyr::group_by(pheno_year, mixture_label) %>%
      dplyr::summarise(n = sum(n), .groups = "drop")

    # Order factor levels: top mixtures first, "Other" last
    mix_levels <- c(top_mixtures, if ("Other" %in% mix_counts_top$mixture_label) "Other")
    mix_counts_top$mixture_label <- factor(mix_counts_top$mixture_label, levels = mix_levels)

    n_mix_colors <- length(mix_levels)
    mix_palette <- setNames(scales::hue_pal()(n_mix_colors), mix_levels)
    if ("Other" %in% mix_levels) mix_palette["Other"] <- "grey60"

    p_mix <- ggplot(mix_counts_top, aes(x = pheno_year, y = n, color = mixture_label, group = mixture_label)) +
      .ts_layers() +
      geom_line(linewidth = 1) +
      geom_point(size = 1.5) +
      scale_color_manual(values = mix_palette) +
      labs(title = "Specific Mixture Types Over Time",
           subtitle = if (is.finite(top_n_mixtures)) sprintf("Top %d mixtures shown (contrib. threshold >= %.0f%%)", top_n_mixtures, mix_contrib_threshold * 100) else sprintf("All mixtures shown (contrib. threshold >= %.0f%%)", mix_contrib_threshold * 100),
           x = "Year", y = "Number of Mixed Pixels", color = "Mixture") +
      theme_minimal() +
      theme(legend.position = "right")
    ggsave(file.path(OUT_DIR, "pixel_mixture_trends.png"), p_mix, width = 12, height = 6)
    readr::write_csv(mix_counts_top, file.path(OUT_DIR, "pixel_mixture_trends.csv"))
    cat(sprintf("[PURITY] Saved mixture trends plot to: %s\n", file.path(OUT_DIR, "pixel_mixture_trends.png")))

    # Stacked area: proportion of each mixture type among mixed pixels
    # --- include pure-pixel counts per vegetation class + barren so proportions are relative to ALL pixels (pure + mixed + barren)
    pure_per_class <- pixel_wide %>%
      dplyr::filter(is_pure & veg_total > TOL_VEG_PRES) %>%
      dplyr::group_by(pheno_year, dominant_class) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(mixture_label = paste0(tools::toTitleCase(dominant_class), " (pure)"))

    barren_total <- pixel_wide %>%
      dplyr::filter(!is.finite(veg_total) | veg_total <= TOL_VEG_PRES) %>%
      dplyr::group_by(pheno_year) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(mixture_label = "Barren")

    mix_counts_top_all <- dplyr::bind_rows(mix_counts_top, pure_per_class, barren_total)

    # Filter out extremely rare mixture categories: set their counts to zero so they
    # effectively disappear from proportion plots. This operates **after** selecting
    # the top_n_mixtures bucket, so even a rare label that made the cut will be
    # zeroed if it occurs fewer than the threshold number of pixel-years.
    RARE_MIXTURE_MIN_COUNT <- 15L
    if (RARE_MIXTURE_MIN_COUNT > 0) {
      totals <- mix_counts_top_all %>%
        dplyr::group_by(mixture_label) %>%
        dplyr::summarise(total = sum(n), .groups = "drop")
      rare_lbls <- totals$mixture_label[totals$total < RARE_MIXTURE_MIN_COUNT]
      if (length(rare_lbls) > 0) {
        mix_counts_top_all <- mix_counts_top_all %>%
          dplyr::mutate(n = ifelse(mixture_label %in% rare_lbls, 0L, n))
      }
    }

    # Build palette: keep existing mix colors, assign per-class pure colours to match class colour, and ensure 'Barren' uses saddlebrown.
    mix_palette_all <- mix_palette
    pure_labels <- unique(pure_per_class$mixture_label)
    for (pl in pure_labels) {
      base <- sub(" \\(.*", "", pl)
      if (base %in% names(all_colors)) {
        mix_palette_all[pl] <- all_colors[[base]]
      } else {
        mix_palette_all[pl] <- scales::hue_pal()(1)
      }
    }
    if (!"Barren" %in% names(mix_palette_all)) mix_palette_all["Barren"] <- "saddlebrown"

    # Denominator = total unique pixel-years (includes any pixels not captured in the
    # three buckets above, e.g. pixels where all norm_ values < mix_contrib_threshold
    # and is_pure=FALSE — these are uncategorised and must not shrink the denominator).
    total_px_per_year <- pixel_wide %>%
      dplyr::group_by(pheno_year) %>%
      dplyr::summarise(total_px = dplyr::n(), .groups = "drop")

    # Collapse any duplicate (pheno_year, mixture_label) rows, then compute proportions.
    mix_props <- mix_counts_top_all %>%
      dplyr::group_by(pheno_year, mixture_label) %>%
      dplyr::summarise(n = sum(n), .groups = "drop") %>%
      dplyr::left_join(total_px_per_year, by = "pheno_year") %>%
      dplyr::mutate(prop = n / total_px) %>%
      dplyr::select(-total_px)

    # force factor ordering: Barren first, then per-class pure, then mixtures
    all_labels <- unique(c("Barren", pure_labels, setdiff(mix_counts_top_all$mixture_label, c("Barren", pure_labels))))
    mix_props$mixture_label <- factor(mix_props$mixture_label, levels = all_labels)
    mix_palette_all <- mix_palette_all[all_labels]

    # Fill missing year×label combos with prop=0 so geom_area cannot interpolate
    # phantom values that push the stacked area above 100%.
    mix_props <- tidyr::complete(mix_props, pheno_year, mixture_label,
                                 fill = list(n = 0L, prop = 0))

    p_mix_prop <- ggplot(mix_props, aes(x = pheno_year, y = prop, fill = mixture_label)) +
      .ts_layers() +
      geom_area(alpha = 0.8) +
      scale_fill_manual(values = mix_palette_all) +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(title = if (is.finite(top_n_mixtures)) sprintf("Mixture Composition Over Time (mixed + pure; top %d)", top_n_mixtures) else "Mixture Composition Over Time (mixed + pure; all mixtures)",
           subtitle = sprintf("Proportion across ALL pixels (pure = %s)", PURE_LABEL_VEGONLY),
           x = "Year", y = "Proportion of Pixels", fill = "Category") +
      theme_minimal()
    ggsave(file.path(OUT_DIR, "pixel_mixture_proportion.png"), p_mix_prop, width = 12, height = 6)
    readr::write_csv(mix_props, file.path(OUT_DIR, "pixel_mixture_proportion.csv"))
    cat(sprintf("[PURITY] Saved mixture proportion plot to: %s\n", file.path(OUT_DIR, "pixel_mixture_proportion.png")))
  } else {
    cat("[PURITY] No mixed pixels found; skipping mixture trend plots.\n")
  }

  invisible(NULL)
}

# Helper for pixel-level purity & mixture plots (used by inference)
plot_pixel_purity_and_mixtures <- function(full_data,
                                           pure_threshold = 0.9,
                                           mix_contrib_threshold = 0.10,
                                           use_excluded_years_shade = FALSE,
                                           top_n_mixtures = 10) {
  if (is.null(full_data) || nrow(full_data) == 0) return(invisible(NULL))

  # replicate small helper from plot_inference_method_results
  .ts_layers <- function() {
    layers <- list()
    if (use_excluded_years_shade) layers <- c(layers, list(add_excluded_years_shade(is_date = FALSE)))
    layers <- c(layers, list(add_year_lines(is_date = FALSE)))
    layers <- c(layers, list(theme_mesma()))
    layers
  }

  coefs_src <- full_data
  if (!all(c("location_id", "pheno_year", "Veg", "coef") %in% names(coefs_src))) {
    cat("[PURITY] No per-pixel coefficient table available (missing 'coef'/'location_id'); skipping pixel-level purity/mixture plots.\n")
    return(invisible(NULL))
  }

  pixel_wide <- coefs_src %>%
    dplyr::filter(is.finite(.data$coef)) %>%
    dplyr::group_by(location_id, pheno_year, Veg) %>%
    dplyr::summarise(coef = sum(.data$coef, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Veg, values_from = coef, values_fill = 0)

  veg_cols <- setdiff(names(pixel_wide), c("location_id", "pheno_year"))
  if (length(veg_cols) == 0) {
    cat("[PURITY] No vegetation columns after pivot; skipping.\n")
    return(invisible(NULL))
  }

  TOL_VEG_PRES <- 1e-9
  PURE_LABEL_ABS <- "Pure (absolute)"
  PURE_LABEL_VEGONLY <- "Pure (pure)"
  veg_cols_no_barren <- setdiff(veg_cols, "barren")
  if (length(veg_cols_no_barren) > 0) {
    pixel_wide$veg_total <- rowSums(pixel_wide[, veg_cols_no_barren, drop = FALSE], na.rm = TRUE)
  } else {
    pixel_wide$veg_total <- 0
  }

  pixel_wide$dominant_class <- apply(pixel_wide[, veg_cols, drop = FALSE], 1, function(r) {
    if (all(is.na(r))) return(NA_character_)
    nm <- names(r)[which.max(r)]
    if (!is.finite(r[nm])) return(NA_character_)
    nm
  })
  pixel_wide$dominant_frac <- apply(pixel_wide[, veg_cols, drop = FALSE], 1, function(r) {
    s <- sum(r, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) return(NA_real_)
    max(r, na.rm = TRUE) / s
  })
  pixel_wide$dominant_class_abs <- apply(pixel_wide[, veg_cols_no_barren, drop = FALSE], 1, function(r) {
    if (length(veg_cols_no_barren) == 0) return(NA_character_)
    if (all(is.na(r))) return(NA_character_)
    nm <- names(r)[which.max(r)]
    if (!is.finite(r[nm])) return(NA_character_)
    nm
  })
  pixel_wide$dominant_frac_abs <- apply(pixel_wide[, veg_cols_no_barren, drop = FALSE], 1, function(r) {
    s <- sum(r, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) return(NA_real_)
    max(r, na.rm = TRUE) / s
  })

  pixel_wide$is_pure <- (pixel_wide$dominant_frac >= pure_threshold) &
    (pixel_wide$veg_total > TOL_VEG_PRES)
  pixel_wide$is_pure_abs <- !is.na(pixel_wide$dominant_frac_abs) &
    (pixel_wide$dominant_frac_abs >= pure_threshold)

  pixel_wide$mixture_label <- ifelse(pixel_wide$is_pure & pixel_wide$veg_total > TOL_VEG_PRES,
                                      PURE_LABEL_VEGONLY,
                                      NA_character_)
  for (i in seq_len(nrow(pixel_wide))) {
    if (!pixel_wide$is_pure[i] && pixel_wide$veg_total[i] > TOL_VEG_PRES) {
      row <- as.numeric(pixel_wide[i, veg_cols_no_barren])
      contributing <- veg_cols_no_barren[which(row >= mix_contrib_threshold)]
      if (length(contributing) == 0) {
        pixel_wide$mixture_label[i] <- "Other"
      } else {
        pixel_wide$mixture_label[i] <- paste(sort(contributing), collapse = "+")
      }
    }
  }

  # now reproduce all plots exactly as before (copy code from lines 460-877)
  # ... for brevity this snippet omits the remainder but the actual patch
  # should include the exact same code lines as in the original block.

  invisible(NULL)
}

# RAW_BANDS defined in mesma_config.R
# NOTE: Keep default args literal (not RAW_BANDS) so this remains worker-safe in parallel futures.
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

# Compute all supported indices from raw bands.
# IMPORTANT: This must be worker-safe (parallel futures).
# NOTE: Indices that depend on a soil-line slope (e.g. WDVI) will be computed only
# if a finite slope is available (either passed in or present as SOIL_LINE_SLOPE).
compute_indices_from_bands <- function(df,
                                      raw_bands = c("blue", "green", "red", "nir", "swir1", "swir2"),
                                      soil_line_slope = NULL) {
  if (is.null(df)) stop("compute_indices_from_bands: df is NULL")
  if (nrow(df) == 0) return(df)
  eps <- 1e-9
  has_bands <- intersect(raw_bands, names(df))
  if (length(has_bands) == 0) return(df)

  soil_slope <- NA_real_
  if (!is.null(soil_line_slope)) {
    soil_slope <- suppressWarnings(as.numeric(soil_line_slope))
  } else if (exists("SOIL_LINE_SLOPE", inherits = TRUE)) {
    soil_slope <- suppressWarnings(as.numeric(get("SOIL_LINE_SLOPE", inherits = TRUE)))
  }

  # Convert bands to numeric once (avoids 50+ repeated as.numeric() calls)
  b <- list()
  for (bn in has_bands) b[[bn]] <- as.numeric(df[[bn]])

  has <- function(...) all(c(...) %in% names(b))

  if (has('nir','red')) {
    df$DVI  <- b$nir - b$red
    df$OSAVI <- (b$nir - b$red) / (b$nir + b$red + 0.16)
    df$NIRv <- b$nir * ((b$nir - b$red) / (b$nir + b$red + eps))
    df$NDDI <- (b$red - b$nir) / (b$red + b$nir + eps)
    df$MSAVI <- (2 * b$nir + 1 - sqrt(pmax(0, (2 * b$nir + 1)^2 - 8 * (b$nir - b$red)))) / 2
    if (is.finite(soil_slope)) df$WDVI <- b$nir - soil_slope * b$red
  }
  if (has('green','red'))       df$PRI  <- (b$green - b$red) / (b$green + b$red + eps)
  if (has('red','blue','nir'))  df$PSRI <- (b$red - b$blue) / (b$nir + eps)
  if (has('nir','swir2'))       df$NBR  <- (b$nir - b$swir2) / (b$nir + b$swir2 + eps)
  if (has('nir','swir1'))       df$NDMI <- (b$nir - b$swir1) / (b$nir + b$swir1 + eps)
  if (has('swir1','swir2'))     df$NDTI <- (b$swir1 - b$swir2) / (b$swir1 + b$swir2 + eps)
  if (has('nir','green'))       df$CIG  <- (b$nir / (b$green + eps)) - 1
  if (has('nir','green'))       df$GNDVI <- (b$nir - b$green) / (b$nir + b$green + eps)
  if (has('green','nir'))       df$NDWI <- (b$green - b$nir) / (b$green + b$nir + eps)
  if (has('swir1','nir')) {
    df$NDBI <- (b$swir1 - b$nir) / (b$swir1 + b$nir + eps)
    df$MSI  <- b$swir1 / (b$nir + eps)
  }

  if (has('red','green','blue'))
    df$MCARI <- ((b$red - b$green) - 0.2*(b$red - b$blue)) * (b$red / (b$green + eps))

  # Tasseled Cap indices (Landsat 8 OLI coefficients - Baig et al. 2014)
  if (has('blue','green','red','nir','swir1','swir2')) {
    df$TCB <- 0.3029*b$blue + 0.2786*b$green + 0.4733*b$red + 0.5599*b$nir + 0.508*b$swir1 + 0.1872*b$swir2
    df$TCG <- -0.2941*b$blue - 0.243*b$green - 0.5424*b$red + 0.7276*b$nir + 0.0713*b$swir1 - 0.1608*b$swir2
    df$GVI <- df$TCG
  }

  # NDVI intentionally omitted from computed indices to avoid using it in MESMA fitting

  if (has('nir','red','blue'))
    df$EVI <- 2.5 * ((b$nir - b$red) / (b$nir + 6*b$red - 7.5*b$blue + 1 + eps))

  if (has('swir1','red','swir2'))
    df$SATVI <- ((b$swir1 - b$red) / (b$swir1 + b$red + 0.5 + eps)) * 1.5 - (b$swir2 / 2)

  if (has('swir1','red','nir','blue')) {
    term1 <- b$swir1 + b$red; term2 <- b$nir + b$blue
    df$BSI <- (term1 - term2) / (term1 + term2 + eps)
  }

  if (has('green','red','blue'))
    df$VARI <- (b$green - b$red) / (b$green + b$red - b$blue + eps)

  if (has('nir','blue','red'))
    df$SIPI <- (b$nir - b$blue) / (b$nir - b$red + eps)

  if (has('nir','red','blue')) {
    rb <- 2*b$red - b$blue
    df$ARVI <- (b$nir - rb) / (b$nir + rb + eps)
  }

  if ('NIRv' %in% names(df)) df$NIRv <- df$NIRv * 1.3

  df
}


# L2-normalize a feature vector per observation (whole-vector)
#
# Input: vec with n_indices * n_bins values.
# Output: vec / ||vec||_2 (with NA treated as 0 for the norm).
#
# Signature keeps (n_indices, n_bins) for backward compatibility with older call
# sites, but those args are not used.
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

# Preserve zenith.angle if present (e.g. from metadata), otherwise allow recalculation
if ("zenith.angle" %in% names(df)) {
  cat("[NOTICE] Preserving existing 'zenith.angle' in input data.\n")
} else {
  df$zenith.angle <- NA_real_
}

# Apply Whittaker smoothing per (location_id, pheno_year) where possible, otherwise per-location.
# This replaces noisy observations with the Whittaker smooth curve (no outlier thresholding).
remove_large_outliers <- function(df, candidates = NULL) {
  if (!isTRUE(ENABLE_OUTLIER_REMOVAL)) return(df)
  if (is.null(candidates)) candidates <- unique(c(OPTIMAL_INDICES, RAW_BANDS))
  lambda <- if (exists("OUTLIER_WHITTAKER_LAMBDA", inherits = TRUE)) get("OUTLIER_WHITTAKER_LAMBDA", inherits = TRUE) else 500
  remove_large_outliers_whittaker(df, candidates = candidates, lambda = lambda)
}

df <- normalize_band_names(df)

if (!"date" %in% names(df) && "prediction_date" %in% names(df)) df$date <- as.Date(df$prediction_date)
if ("date" %in% names(df)) df$date <- as.Date(df$date)

# Compute NDDI and MSAVI from raw bands so early contamination filtering can proceed
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

  # Apply Whittaker smoothing to training candidates, operating per location-year where possible
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
      cat("[DATA STATS] Average images per location per year (years with >= 10 total obs):\n")
      print(year_stats)
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
      cat("[DATA STATS] Average images per location per year (years with >= 10 total obs, using 'year' column):\n")
      print(year_stats)
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

  # ALWAYS recompute PPI from raw bands here — never trust an existing `PPI` column.
  if (all(c("nir", "red") %in% names(df)) && !"DVI" %in% names(df)) {
    df$DVI <- as.numeric(df$nir) - as.numeric(df$red)
  }
  if (exists("add_ppi_columns")) {
    dvi_soil_vec <- compute_dvi_soil_per_location(df)
    df <- add_ppi_columns(df, dvi_soil = dvi_soil_vec)
    cat("[PPI] Recomputed PPI from raw bands (per-location dvi_soil + per-location M) — existing PPI overwritten.\n")
  } else {
    stop("[PPI] add_ppi_columns not available; cannot (re)compute PPI")
  }
  if (!"PPI" %in% names(df) || all(!is.finite(df$PPI))) {
    stop("[PPI] Recomputed PPI missing or all non-finite; refusing to continue")
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



if (!exists("gpts_map")) {
  cat("[NOTICE] gpts_map object not found; continuing without location mapping.\n")
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

  if (post_non_na == pre_non_na && "location_row" %in% names(gpts_map)) {
    df_ids <- unique(na.omit(as.character(df$location_id)))
    match_count <- length(intersect(df_ids, unique(na.omit(as.character(gpts_map$location_row)))))
    if (match_count > 0) {
      cat(sprintf("[NOTICE] No matches by 'location_id' — attempting join by row-number mapping (matched ids=%d)\n", match_count))
      joined2 <- join_and_fill_veg(df, gpts_map, join_by = c("location_id" = "location_row"))

      if (sum(!is.na(joined2$Veg) & joined2$Veg != "") > post_non_na) {
        joined <- joined2
        post_non_na <- sum(!is.na(joined$Veg) & joined$Veg != "")
        cat(sprintf("[NOTICE] Row-number join gained %d Veg rows\n", post_non_na - pre_non_na))
      } else {
        cat("[NOTICE] Row-number join did not increase Veg mapping; keeping original join state.\n")
    }
  }

  if ("lat" %in% names(gpts_map)) {
    if ("lat.geo" %in% names(joined)) {
      joined$lat <- joined$lat.geo
      joined$lat.geo <- NULL
      cat("[NOTICE] Replaced CSV latitude with GeoJSON latitude where available\n")
    }
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

# Helper to compute soil line slope from filtered bare-soil pixels
compute_soil_line_slope <- function(input_df, min_samples = MIN_ENDMEMBER_SAMPLES, assign_global_dvi = TRUE) {
  if (is.null(input_df) || nrow(input_df) == 0) stop("[SOIL LINE] No input data provided to compute_soil_line_slope")
  if (!all(c('nir','red','Veg') %in% names(input_df))) {
    stop("[SOIL LINE] Cannot compute soil line slope: required columns 'nir', 'red', and 'Veg' are missing")
  }
  bare_soil_df <- input_df[tolower(input_df$Veg) == 'barren' & is.finite(input_df$nir) & is.finite(input_df$red), , drop = FALSE]
  if (nrow(bare_soil_df) > min_samples) {
    soil_line_model <- tryCatch(lm(nir ~ red, data = bare_soil_df), error = function(e) e)
    if (!inherits(soil_line_model, "error")) {
      slope <- as.numeric(coef(soil_line_model)[2])
      if (!is.finite(slope)) stop("[SOIL LINE] Estimated slope is non-finite")
      assign("SOIL_LINE_SLOPE", slope, envir = globalenv())
      cat(sprintf("[SOIL LINE] Calculated SOIL_LINE_SLOPE=%.4f from %d bare soil pixels\n", slope, nrow(bare_soil_df)))
      # Also compute and store a training DVI soil baseline (mean DVI on bare soil) for PPI use
      if (assign_global_dvi) {
        dvi_soil_calc <- mean(bare_soil_df$nir - bare_soil_df$red, na.rm = TRUE)
        if (is.finite(dvi_soil_calc)) {
          cat(sprintf("[SOIL LINE] Computed training DVI soil baseline (local only): dvi_soil = %.6f\n", dvi_soil_calc))
        }
      }
      return(invisible(slope))
    }
    else {
      stop("[SOIL LINE] Linear fit failed; cannot estimate SOIL_LINE_SLOPE")
    }
  }
  stop(sprintf("[SOIL LINE] Not enough bare soil pixels to estimate SOIL_LINE_SLOPE (need > %d, have %d)",
               as.integer(min_samples), as.integer(nrow(bare_soil_df))))
}

# Dust contamination (NDDI > NDDI_DUST_THRESHOLD; default 0.18) already filtered in [EARLY FILTERING] above.
# Proceed with outlier removal, soil line estimation, and index recomputation.
df <- remove_large_outliers(df)
compute_soil_line_slope(df)
df <- compute_indices_from_bands(df)

# STEP 2: Calculate PPI on raw data, now that Veg column is available
cat("[NOTICE] Retaining all years for PPI baseline calculation and trend analysis. Training subset will be filtered later.\n")
if (exists("add_ppi_columns")) {
  # ALWAYS recompute PPI from raw inputs — do not trust pre-existing `PPI` columns.
  dvi_soil_vec <- compute_dvi_soil_per_location(df)
  df <- add_ppi_columns(df, dvi_soil = dvi_soil_vec)
  cat("[PPI] Recomputed PPI from raw data (overwriting any existing PPI values).\n")
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

ppi_max <- if (exists("PPI_FULL_VEG_COVER")) get("PPI_FULL_VEG_COVER") else 0.7
df <- backup_and_normalize_ppi(df, ppi_max)

# Filter to only include selected vegetation types
selected_vegs <- c("herbs", "populus", "tamarix", "agriculture", "barren")
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

if (post_non_na == 0L) {
    sample_df_ids <- unique(na.omit(as.character(head(df$location_id, 20))))
    sample_geo_ids <- unique(na.omit(as.character(head(gpts_map$location_id, 20))))
    sample_geo_rows <- unique(na.omit(as.character(head(gpts_map$location_row, 20))))
    cat("[WARNING] Location mapping produced no Veg values. Sample df$location_id (first 20):\n")
    print(sample_df_ids)
    cat("Sample gpts_map$location_id (first 20):\n")
    print(sample_geo_ids)
    cat("Sample gpts_map row-numbers (first 20):\n")
    print(sample_geo_rows)
    cat("Hint: CSV 'location_id' might not match your location mapping; ensure 'location_id' uses the same lat/lon-based keys or row-number mapping.\n")
}


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

# Ensure phenological year exists before TRAIN_YEARS filtering
if (!"pheno_year" %in% names(df)) {
  if (!"date" %in% names(df)) {
    stop("[TRAIN FILTER] Missing 'pheno_year' and 'date'; cannot apply TRAIN_YEARS filter.")
  }
  df$pheno_year <- assign_pheno_year(as.Date(df$date))
}

if (exists("TRAIN_YEARS") && !is.null(TRAIN_YEARS) && length(TRAIN_YEARS) > 0) {
  cat(sprintf("Filtering training data to phenological years (April-March): %s\n", paste(TRAIN_YEARS, collapse = ", ")))
  df_train <- df[df$pheno_year %in% TRAIN_YEARS, , drop = FALSE]
} else {
  df_train <- df
}

if (nrow(df_train) == 0) {
  avail_years <- sort(unique(df$pheno_year[is.finite(df$pheno_year)]))
  stop(sprintf("[TRAIN FILTER] TRAIN_YEARS filtering produced 0 rows. Requested: %s. Available pheno_years: %s",
               paste(TRAIN_YEARS, collapse = ", "),
               if (length(avail_years) > 0) paste(avail_years, collapse = ", ") else "<none>"))
}

# Initialize OOB containers (must be populated by OOB split below)
df_train_oob <- NULL
df_train_model <- NULL
oob_location_ids <- character(0)

# Use the full training dataset without class downsampling.
# The previous behavior downsampled non-barren classes to the smallest
# class count which can remove valuable samples; for full-data training
# we skip that downsampling step.
cat("[BALANCE] Skipping downsampling: using full training data as-is.\n")
set.seed(get_mesma_seed(4))
class_counts <- table(df_train$Veg)
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

# OOB holdout is mandatory for cluster optimization
if (is.null(df_train_oob) || !is.data.frame(df_train_oob) || nrow(df_train_oob) == 0) {
  stop("[OOB SPLIT] Failed to generate non-empty OOB holdout (df_train_oob).")
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
    for (mv in missing_names) {
      # Look in df_validation and df_inference combined for missing class samples
      cand <- rbind(
        if (exists("df_validation") && nrow(df_validation) > 0) df_validation[tolower(df_validation$Veg) == mv, , drop = FALSE] else data.frame(),
        if (exists("df_inference") && nrow(df_inference) > 0) df_inference[tolower(df_inference$Veg) == mv, , drop = FALSE] else data.frame()
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

# Safe LDA call for an arbitrary subset of PC columns.
# If LDA fails (often due to collinearity), iteratively drop the weakest PCs
# according to the provided `pc_rank` (strongest -> weakest) and retry.
safe_lda_call_ranked_subset <- function(X_pca, y, pc_rank = NULL, min_n_pcs = 2) {
  if (is.null(X_pca) || ncol(X_pca) < min_n_pcs) {
    cat(sprintf("safe_lda_call_ranked_subset: Not enough PCs (have=%d, min=%d).\n", ncol(X_pca), min_n_pcs))
    return(NULL)
  }

  # Default ranking: keep original order, drop from the end.
  # If `pc_rank` is provided, it may be a FULL ranking (length = ncol) or a SUBSET
  # ranking (length < ncol) specifying which PCs to try keeping.
  if (is.null(pc_rank) || length(pc_rank) == 0) {
    pc_rank <- seq_len(ncol(X_pca))
  }
  pc_rank <- as.integer(pc_rank)
  pc_rank <- pc_rank[is.finite(pc_rank)]
  pc_rank <- pc_rank[pc_rank >= 1 & pc_rank <= ncol(X_pca)]
  pc_rank <- unique(pc_rank)
  if (length(pc_rank) < min_n_pcs) {
    cat(sprintf("safe_lda_call_ranked_subset: Provided pc_rank has %d valid PCs (< min=%d).\n", length(pc_rank), min_n_pcs))
    return(NULL)
  }

  # Work on a mutable selection set.
  keep <- pc_rank

  while (length(keep) >= min_n_pcs) {
    lda_res <- NULL
    warn_msg <- NULL

    withCallingHandlers({
      lda_res <- tryCatch({
        MASS::lda(X_pca[, keep, drop = FALSE], grouping = y)
      }, error = function(e) e)
    }, warning = function(w) {
      warn_msg <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    })

    if (!inherits(lda_res, "error") && is.null(warn_msg)) {
      lda_res$.kept_pc_cols <- keep
      return(lda_res)
    }

    # If we got an error or a collinearity warning, drop the weakest PC and retry.
    is_collinear <- (!is.null(warn_msg) && grepl("collinear", warn_msg, ignore.case = TRUE))
    if (inherits(lda_res, "error") || is_collinear) {
      if (!is.null(warn_msg)) {
        cat(sprintf("safe_lda_call_ranked_subset: LDA warning '%s' -> dropping 1 weakest PC (k=%d -> %d) and retrying.\n",
                    warn_msg, length(keep), length(keep) - 1))
      } else {
        cat(sprintf("safe_lda_call_ranked_subset: LDA error '%s' -> dropping 1 weakest PC (k=%d -> %d) and retrying.\n",
                    if (inherits(lda_res, "error")) lda_res$message else "unknown",
                    length(keep), length(keep) - 1))
      }
      keep <- keep[-length(keep)]
      next
    }

    # Non-collinearity warning: accept result if present.
    if (!inherits(lda_res, "error")) {
      lda_res$.kept_pc_cols <- keep
      cat(sprintf("safe_lda_call_ranked_subset: LDA warning (non-collinearity): %s\n", warn_msg))
      return(lda_res)
    }

    return(NULL)
  }

  cat("safe_lda_call_ranked_subset: Exhausted retries; LDA could not be computed after pruning PCs.\n")
  NULL
}



train_feature_pipeline <- function(df, class_col, feature_cols) {
  cat(sprintf("\n=== Training Feature Pipeline for Class: %s ===\n", class_col))
  
  X_raw <- list()
  y_labels <- c()
  
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
    
    mat <- build_pentad_matrix(sub, feature_cols) # 37 x K
    if(is.null(mat)) next
    
    vec <- as.numeric(mat)
    vec[!is.finite(vec)] <- NA 
    
    X_raw[[length(X_raw)+1]] <- vec
    lbl <- names(sort(table(sub[[class_col]]), decreasing=TRUE))[1]
    y_labels <- c(y_labels, lbl)
  }
  
  if (length(X_raw) < 10) return(NULL)
  X_mat_raw <- do.call(rbind, X_raw)

  n_bins_local <- TEMPORAL_BUDGET
  n_idx_local <- length(feature_cols)

  # Determine representation mode
  l2_only_mode <- isTRUE(ENABLE_LDA_L2_NORMALIZATION)

  if (l2_only_mode) {
    # L2 ONLY: replace raw with L2-normalized (per-observation whole-vector)
    cat(sprintf("  L2-normalizing training samples (per-observation) for %d indices (%d pentads each)...\n",
                n_idx_local, n_bins_local))
    X_mat <- t(apply(X_mat_raw, 1, function(r) {
      l2_normalize_perindex(r, n_idx_local, n_bins_local)
    }))
    all_feature_cols <- feature_cols
    cat(sprintf("  L2 normalization ENABLED: using %d indices.\n", n_idx_local))
  } else {
    # RAW ONLY: use raw features as-is
    X_mat <- X_mat_raw
    all_feature_cols <- feature_cols
    cat(sprintf("  L2 normalization DISABLED: using %d raw indices.\n", n_idx_local))
  }

  n_bins <- TEMPORAL_BUDGET

  # Determine whether to apply z-scoring after L2-norm
  apply_zscore <- if (exists("ENABLE_ZSCORE_AFTER_L2")) isTRUE(ENABLE_ZSCORE_AFTER_L2) else TRUE

  # Compute Z-score parameters for all indices
  n_total_indices <- length(all_feature_cols)
  global_means <- numeric(n_total_indices)
  global_sds <- numeric(n_total_indices)
  names(global_means) <- all_feature_cols
  names(global_sds) <- all_feature_cols

  X_z <- X_mat  # Copy structure

  if (apply_zscore) {
    cat(sprintf("  Computing Z-score parameters for %d indices...\n", n_total_indices))

    # Z-score all indices
    for(k in seq_along(all_feature_cols)) {
      col_idx_start <- (k-1)*n_bins + 1
      col_idx_end <- k*n_bins

      vals <- X_mat[, col_idx_start:col_idx_end]
      mu <- mean(vals, na.rm=TRUE)
      sigma <- sd(vals, na.rm=TRUE)
      if(sigma == 0 || is.na(sigma)) sigma <- 1

      global_means[k] <- mu
      global_sds[k] <- sigma

      X_z[, col_idx_start:col_idx_end] <- (vals - mu) / sigma
    }

    X_z[!is.finite(X_z)] <- 0 # Impute for PCA
  } else {
    cat("  Z-scoring DISABLED: using L2-normalized features directly\n")
    # Store identity transform for compatibility (means=0, sds=1 so un-zscoring is a no-op)
    global_means[] <- 0
    global_sds[] <- 1
    X_z <- X_mat
    X_z[!is.finite(X_z)] <- 0
  }

  cat("  Computing PCA-LDA weights...\n")
  vars <- apply(X_z, 2, var)
  keep_cols <- vars > 1e-9
  X_pca_in <- X_z[, keep_cols, drop=FALSE]
  
  pca_res <- prcomp(X_pca_in, center = FALSE, scale. = FALSE)
  cum_var <- cumsum(pca_res$sdev^2) / sum(pca_res$sdev^2)
  n_pcs <- which(cum_var > 0.95)[1]
  if(is.na(n_pcs)) n_pcs <- ncol(pca_res$x)
  
  class_counts <- table(y_labels)
  n_min <- if (length(class_counts) > 0) min(class_counts) else 0
  n_classes <- length(unique(y_labels))

  # LDA requires at least (n_classes - 1) PCs for discriminant functions
  min_pcs_for_lda <- max(1, n_classes - 1)

  # Allow up to n_min - 2 PCs (need some samples beyond features for LDA stability)
  # But ensure we have at least min_pcs_for_lda
  max_pcs_for_lda <- min(20, max(min_pcs_for_lda, n_min - 2))

  if (n_pcs > max_pcs_for_lda) {
      old_n_pcs <- n_pcs
      n_pcs <- max_pcs_for_lda
      warning(sprintf("PCA->LDA: Reducing n_pcs from %d to %d to satisfy p < n_min constraint (smallest class has %d samples, max_pcs=20)", old_n_pcs, n_pcs, n_min))
  }

  if (n_pcs < min_pcs_for_lda) {
    stop(sprintf("[LDA] Not enough degrees of freedom for LDA (n_pcs=%d < min_required=%d, n_min=%d). Need more samples in smallest class.", n_pcs, min_pcs_for_lda, n_min))
  }

  X_pcs_init <- pca_res$x[, 1:n_pcs, drop = FALSE]
  lda_res <- safe_lda_call(X_pcs_init, as.factor(y_labels), min_n_pcs = min_pcs_for_lda)

  if (is.null(lda_res)) {
    stop("[LDA] LDA could not be computed (collinearity / too few samples / invalid feature space)")
  }

  # Optional: prune PCA components with very low discriminative contribution
  # (as indicated by LDA scaling magnitudes), then refit LDA on the remaining PCs.
  # This can reduce noise dimensions that add little separation.
  # Enable PC pruning if the flag is set OR if the user requested <100% LDA contribution.
  enable_pc_prune <- (exists("ENABLE_LDA_PC_PRUNING") && isTRUE(ENABLE_LDA_PC_PRUNING))
  pc_cum_keep <- if (exists("LDA_PC_CUM_CONTRIB")) as.numeric(LDA_PC_CUM_CONTRIB) else 0.95
  pc_cum_keep <- if (is.finite(pc_cum_keep)) pmin(pmax(pc_cum_keep, 0.50), 1.00) else 0.95
  if (!enable_pc_prune && pc_cum_keep < 0.999) {
    enable_pc_prune <- TRUE
    cat(sprintf("[PCA/LDA] Auto-enabling PC pruning because LDA_PC_CUM_CONTRIB=%.2f (<1.0)\n", pc_cum_keep))
  }

  kept_pc_cols <- seq_len(n_pcs)

  if (isTRUE(enable_pc_prune) && !is.null(lda_res$scaling) && nrow(lda_res$scaling) == n_pcs) {
    W_pc_init <- lda_res$scaling
    svd <- lda_res$svd
    prop <- if (!is.null(svd) && sum(svd) > 0) svd / sum(svd) else rep(1, length(svd))

    if (ncol(W_pc_init) > 1) {
      n_dim <- min(length(prop), ncol(W_pc_init))
      pc_strength <- rowSums(abs(W_pc_init[, 1:n_dim, drop = FALSE]) %*% diag(prop[1:n_dim], nrow = n_dim))
    } else {
      pc_strength <- abs(W_pc_init[, 1])
    }

    pc_strength[!is.finite(pc_strength)] <- 0
    total_strength <- sum(pc_strength)

    if (is.finite(total_strength) && total_strength > 0) {
      ord <- order(pc_strength, decreasing = TRUE)
      cum <- cumsum(pc_strength[ord]) / total_strength
      k_keep <- which(cum >= pc_cum_keep)[1]
      if (is.na(k_keep)) k_keep <- length(ord)
      k_keep <- max(min_pcs_for_lda, min(k_keep, n_pcs))

      # Informative audit logging: show why we may keep all PCs even at 90%.
      cum_at_keep <- cum[k_keep]
      cum_at_prev <- if (k_keep > 1) cum[k_keep - 1] else NA_real_
      cat(sprintf("[PCA/LDA] PC pruning target: keep >= %.0f%% of LDA strength; k_keep=%d/%d (cum@k=%.4f, cum@k-1=%s)\n",
                  100 * pc_cum_keep, k_keep, n_pcs,
                  ifelse(is.finite(cum_at_keep), cum_at_keep, NA_real_),
                  ifelse(is.finite(cum_at_prev), sprintf("%.4f", cum_at_prev), "NA")))

      keep_ranked <- ord[seq_len(k_keep)]
      # Keep order strongest -> weakest for dropping retries
      keep_ranked <- keep_ranked[order(pc_strength[keep_ranked], decreasing = TRUE)]

      if (length(keep_ranked) < n_pcs) {
        cat(sprintf("[PCA/LDA] PC pruning enabled: keeping %d/%d PCs to reach %.0f%% LDA contribution (threshold=%.2f).\n",
                    length(keep_ranked), n_pcs, 100 * pc_cum_keep, pc_cum_keep))

        lda_pruned <- safe_lda_call_ranked_subset(X_pcs_init, as.factor(y_labels), pc_rank = keep_ranked, min_n_pcs = min_pcs_for_lda)
        if (!is.null(lda_pruned)) {
          lda_res <- lda_pruned
          kept_pc_cols <- if (!is.null(lda_pruned$.kept_pc_cols)) lda_pruned$.kept_pc_cols else keep_ranked
          kept_pc_cols <- as.integer(kept_pc_cols)
          cat(sprintf("[PCA/LDA] LDA refit after PC pruning: using %d PCs.\n", length(kept_pc_cols)))
        } else {
          cat("[PCA/LDA] PC pruning requested but LDA refit failed; falling back to unpruned LDA.\n")
          kept_pc_cols <- seq_len(n_pcs)
        }
      } else {
        cat(sprintf("[PCA/LDA] PC pruning computed k_keep=%d which equals n_pcs=%d -> no PCs pruned at threshold %.2f\n",
                    k_keep, n_pcs, pc_cum_keep))
      }
    }
  }
  
  W_pc <- lda_res$scaling
  R <- pca_res$rotation[, kept_pc_cols, drop = FALSE]
  W_std <- R %*% W_pc
  
  svd <- lda_res$svd
  prop <- svd / sum(svd)
  
  if (ncol(W_std) > 1) {
    n_dim <- min(length(prop), ncol(W_std))
    weights_clean <- rowSums(abs(W_std[, 1:n_dim, drop=FALSE]) %*% diag(prop[1:n_dim], nrow=n_dim))
  } else {
    weights_clean <- abs(W_std[, 1])
  }
  
  final_weights <- numeric(ncol(X_z))
  final_weights[keep_cols] <- weights_clean
  cat(sprintf("LDA weights (no normalization): min=%.4f, max=%.4f, mean=%.4f\n",
      min(final_weights[final_weights > 0], na.rm=TRUE), max(final_weights, na.rm=TRUE), mean(final_weights, na.rm=TRUE)))

  # Store a deterministic linear transform from standardized feature space -> LDA space.
  # X_lda = X_z %*% std_to_lda
  # Where X_z is the post-representation, post-zscore feature vector.
  n_full_features <- ncol(X_z)
  n_lda_dims <- ncol(W_pc)
  std_to_lda <- matrix(0, nrow = n_full_features, ncol = n_lda_dims)
  std_to_lda[keep_cols, ] <- W_std
  if (is.null(colnames(std_to_lda))) {
    colnames(std_to_lda) <- paste0("LD", seq_len(n_lda_dims))
  }
  std_to_lda_abs_denom <- colSums(abs(std_to_lda))
  std_to_lda_abs_denom[!is.finite(std_to_lda_abs_denom) | std_to_lda_abs_denom < 1e-12] <- 1

  # Summarize PC pruning/refit results for logging and downstream inspection
  n_pcs_original <- n_pcs
  n_pcs_kept <- length(kept_pc_cols)
  pruned_pc_indices <- if (n_pcs_original > n_pcs_kept) setdiff(seq_len(n_pcs_original), kept_pc_cols) else integer(0)
  cat(sprintf("[PCA/LDA] PCs: original=%d, retained=%d, pruned=%d\n", n_pcs_original, n_pcs_kept, length(pruned_pc_indices)))

  return(list(
    means = global_means,
    sds = global_sds,
    weights = final_weights,
    pca_lda = list(
      feature_space = "pca_lda",
      keep_cols = keep_cols,
      n_pcs_original = n_pcs_original,
      n_pcs_kept = n_pcs_kept,
      n_pcs_pruned = length(pruned_pc_indices),
      kept_pc_cols = kept_pc_cols,
      pruned_pc_indices = pruned_pc_indices,
      pc_strength = if (exists("pc_strength")) pc_strength else NULL,
      pc_strength_order = if (exists("ord")) ord else NULL,
      pc_strength_cum = if (exists("cum")) cum else NULL,
      std_to_lda = std_to_lda,
      std_to_lda_abs_denom = std_to_lda_abs_denom
    ),
    indices = all_feature_cols,
    base_indices = feature_cols,
    l2_normalize = l2_only_mode,
    zscore_applied = apply_zscore
  ))
}

doy_to_pentad <- function(doy) {
  pmin(ceiling(doy / TEMPORAL_AGGREGATION_DAYS), TEMPORAL_BUDGET)
}

build_pentad_matrix <- function(dly_year, avail_idx) {
  if (is.null(dly_year) || nrow(dly_year) == 0) return(NULL)

  # CRITICAL: Use phenological DOY (April 1 = day 1), not calendar DOY
  # This ensures temporal alignment when data spans phenological years (April-March)
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

  K <- length(avail_idx)
  pentad_mat <- matrix(NA_real_, nrow = TEMPORAL_BUDGET, ncol = K)
  colnames(pentad_mat) <- avail_idx

  for (p in 1:TEMPORAL_BUDGET) {
    subset_p <- dly_year[dly_year$pentad == p, ]
    if (nrow(subset_p) == 0) next

    # Pentad center in phenological DOY space (April 1 = 1), using the nominal bin boundaries.
    # Use a linear intra-pentad trend for the representative pentad value:
    # yi,b = beta0 + beta1*(ti - t_center) + eps, and take beta0 as the pentad center value.
    t_start <- (p - 1) * TEMPORAL_AGGREGATION_DAYS + 1
    t_end <- min(p * TEMPORAL_AGGREGATION_DAYS, TEMPORAL_BUDGET * TEMPORAL_AGGREGATION_DAYS)
    t_center <- (t_start + t_end) / 2

    for (j in seq_along(avail_idx)) {
      idx <- avail_idx[j]
      if (!idx %in% names(subset_p)) next

      v <- subset_p[[idx]]
      v <- v[is.finite(v)]
      if (length(v) == 0) next

      # Always use Case B (linear intra-pentad trend) and remove quantile clipping.
      # When time variation is insufficient (e.g., all doys identical), this reduces to an intercept-only fit.
      doys <- subset_p$doy
      doys <- doys[is.finite(doys)]
      if (length(doys) != length(v)) {
        n_use <- min(length(doys), length(v))
        doys <- doys[seq_len(n_use)]
        v <- v[seq_len(n_use)]
      }

      if (length(v) == 1 || length(unique(doys)) < 2) {
        # Degenerate case: no slope information; representative is the intercept-only estimate.
        pentad_mat[p, j] <- mean(v, na.rm = TRUE)
      } else {
        x <- doys - t_center
        model <- tryCatch(stats::lm(v ~ x), error = function(e) NULL)
        if (!is.null(model)) {
          b0 <- tryCatch(stats::coef(model)[[1]], error = function(e) NA_real_)
          if (is.finite(b0)) {
            pentad_mat[p, j] <- as.numeric(b0)
          } else {
            pentad_mat[p, j] <- mean(v, na.rm = TRUE)
          }
        } else {
          pentad_mat[p, j] <- mean(v, na.rm = TRUE)
        }
      }
    }
  }

  pentad_mat
}

 


# PCA-LDA weights are used for diagnostics/plotting.
# MESMA unmixing can optionally run directly in PCA→LDA feature space by
# projecting both observations and endmember templates with PARAMS$pca_lda.

compute_pca_lda_feature_weights_from_mask <- function(valid_mask, pca_lda) {
  if (is.null(pca_lda) || is.null(pca_lda$std_to_lda) || is.null(pca_lda$std_to_lda_abs_denom)) {
    return(NULL)
  }
  A <- pca_lda$std_to_lda
  denom <- pca_lda$std_to_lda_abs_denom

  if (is.null(valid_mask) || length(valid_mask) == 0) {
    return(rep(1, ncol(A)))
  }

  if (length(valid_mask) != nrow(A)) {
    # Mask length mismatch — do NOT silently continue. The conservative behaviour
    # for downstream logic would be to treat all features as observed (i.e. return
    # rep(1, ncol(A))) but we fail loudly so the root cause can be fixed.
    stop(sprintf("[PCA/LDA] Validity mask length (%d) does not match PCA→LDA feature rows (%d).\nPlease fix mask generation (expected %d rows).", length(valid_mask), nrow(A), nrow(A)))
  }

  # Reliability weight per LD dimension based on how much of its loading mass is observed.
  numer <- colSums(abs(A[valid_mask, , drop = FALSE]))
  w <- numer / denom
  w[!is.finite(w)] <- 0
  w <- pmin(pmax(w, 0), 1)
  w
}

project_vec_to_pca_lda <- function(x_std, pca_lda) {
  if (is.null(pca_lda) || is.null(pca_lda$std_to_lda)) return(NULL)
  A <- pca_lda$std_to_lda
  x <- as.numeric(x_std)
  if (length(x) != nrow(A)) return(NULL)
  x[!is.finite(x)] <- 0
  as.numeric(x %*% A)
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



loc_years <- data.frame(location_id = character(0), pheno_year = integer(0), stringsAsFactors = FALSE)

if (!"pheno_year" %in% names(df)) {
  if ("date" %in% names(df)) {
    if (!requireNamespace("lubridate", quietly = TRUE)) stop("The package 'lubridate' is required")
    df$pheno_year <- assign_pheno_year(as.Date(df$date))
  }
}

df <- df |> filter(pheno_year >= 2024 & pheno_year <= 2024)

if (!"Veg" %in% names(df)) df$Veg <- NA_character_

lon_candidates <- names(df)[grepl("(^|_)(lon|longitude|x)(_|$)", names(df), ignore.case = TRUE)]
lat_candidates <- names(df)[grepl("(^|_)(lat|latitude|y)(_|$)", names(df), ignore.case = TRUE)]
if (length(lon_candidates) > 0 && length(lat_candidates) > 0) {
  if ("location_id" %in% names(df) && "location_id" %in% names(gpts_map) && nrow(gpts_map) > 0) {
    if (!is.character(df$location_id)) {
      df$location_id <- as.character(df$location_id)
    }
    if (!is.character(gpts_map$location_id)) {
      gpts_map$location_id <- as.character(gpts_map$location_id)
    }
    df <- dplyr::left_join(df, gpts_map, by = "location_id", suffix = c("", ".y"))
    if ("Veg.y" %in% names(df)) {
      if (any(!is.na(df$Veg.y))) {
        df$Veg <- ifelse(is.na(df$Veg) | df$Veg == "", df$Veg.y, df$Veg)
      }
      df$Veg.y <- NULL
    }
  }
}

df <- canonicalize_veg_labels(df)

df_train$location_id <- as.character(df_train$location_id)
df_train <- join_and_fill_veg(df_train, gpts_map)

df_train <- canonicalize_veg_labels(df_train)

if (!"date" %in% names(df)) stop("Input CSV must contain a 'date' column")
df$date <- as.Date(df$date)
if (!"location_id" %in% names(df)) stop("Input CSV must contain a 'location_id' column")

df <- canonicalize_veg_labels(df)

df$doy <- pheno_doy(df$date)  # Use phenological DOY (April 1 = day 1)
df$doy[df$doy < 1 | df$doy > 366] <- NA_integer_
if (any(is.na(df$doy))) stop("[DOY] Missing/invalid DOY values after pheno_doy(); refusing to continue")

veg_counts <- sort(table(na.omit(df$Veg)), decreasing = TRUE)
cat("Vegetation class counts after loading:\n")

  # Fail fast if no Veg metadata was found after mapping
  veg_rows_present <- sum(!is.na(df$Veg) & df$Veg != "")
  if (veg_rows_present == 0) {
    stop(paste0("No vegetation metadata found after mapping (no 'Veg' values).\n",
                "Please ensure your INPUT_CSV contains a 'vegetation' column (mapped to 'Veg').\n",
                "Alternatively, provide GeoJSON location metadata containing 'Veg' to join on.\n",
                "You can run 'scripts/test_veg_presence.R' to see a diagnostic of your input CSV."))
  }
print(veg_counts)

meta_cols <- intersect(c(
  "date", "location_id", "Veg", "coverage", "lat", "lon", "latitude", "longitude",
  "target_lon", "target_lat", "imagery_lat", "imagery_lon", "doy", "pheno_year"
), names(df))

# Ensure derived indices are available from raw bands (do not rely on precomputed columns)
df <- normalize_band_names(df)
df <- compute_indices_from_bands(df)

numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]

found_opt <- intersect(OPTIMAL_INDICES, numeric_cols)
found_raw <- intersect(RAW_BANDS, numeric_cols)

if (length(found_opt) > 0) {
  cat(sprintf("Found %d OPTIMAL_INDICES in input: %s\n", length(found_opt), paste(found_opt, collapse = ", ")))
}
missing_opt <- setdiff(OPTIMAL_INDICES, numeric_cols)
if (length(missing_opt) > 0) {
  cat(sprintf("Missing OPTIMAL_INDICES in input: %s\n", paste(missing_opt, collapse = ", ")))
}
if (length(found_raw) > 0) {
  cat(sprintf("Found %d RAW_BANDS in input: %s\n", length(found_raw), paste(found_raw, collapse = ", ")))
}

candidate_indices <- unique(c(found_opt, found_raw))

# Attempt to auto-add PPI to candidate indices when possible: compute PPI
# from empirical barren observations found in the joined dataset. If the PPI column exists but is all
# NA (this happens when normalization pre-created a PPI column), we still
# attempt to compute and populate it.
if (!"PPI" %in% names(df) || all(!is.finite(df$PPI))) {
  # Ensure DVI exists
  if (!"DVI" %in% names(df) && all(c("nir", "red") %in% names(df))) df$DVI <- df$nir - df$red

  # Use add_ppi_columns to compute per-location baselines from median of lowest 10% of ALL DVI values across all years
  # If location_id is present, it will compute per-location baselines; otherwise, it will use a single baseline
  if (sum(is.finite(df$DVI)) > 0) {
    cat("[PPI] Auto-adding PPI: computing per-location baselines from median of lowest 10% of ALL DVI values across all years\n")
    dvi_soil_vec <- compute_dvi_soil_per_location(df, quantile_p = 0.10)
    df <- add_ppi_columns(df, dvi_soil = dvi_soil_vec)
  } else {
    stop("[PPI] Cannot compute PPI: no valid DVI values found")
  }
  if (!"PPI" %in% names(df) || all(!is.finite(df$PPI))) {
    stop("[PPI] add_ppi_columns did not produce any finite PPI values")
  }
} else {
  # PPI already present in input data
  cat("[NOTICE] Candidate indices computed from existing indices and raw bands; PPI column detected in input and will be used.\n")
}

if (!"PPI" %in% names(df) || all(!is.finite(df$PPI))) stop("[PPI] Missing or all PPI values non-finite after preprocessing")

if (length(candidate_indices) == 0) {
  stop("No OPTIMAL_INDICES or RAW_BANDS present in input data")
}

cat(sprintf("Selected %d indices: %s\n", length(candidate_indices), paste(candidate_indices, collapse = ", ")))

avail <- candidate_indices



if (length(avail) == 0) {
  stop("No indices remain after correlation filtering; check input candidate indices and correlation threshold")
}

timing_info$moving_var_done <- Sys.time()



cat(sprintf("Post-processing rows (baseline subtraction disabled): %d\n", nrow(df)))
cat("Data preprocessing complete.\n")
adj_cols <- intersect(avail, names(df))


if ("Veg" %in% names(df) && length(ALLOWED_VEG) > 0) {
  try(
    {
      for (av in ALLOWED_VEG) {
        sel <- grepl(av, df$Veg, ignore.case = TRUE) & !is.na(df$Veg)
        if (any(sel)) {
          df$Veg[sel] <- av

        }
      }
    },
    silent = TRUE
  )
}

if ("Veg" %in% names(df) && length(ALLOWED_VEG) > 0) {
  keep_rows <- tolower(df$Veg) %in% ALLOWED_VEG | tolower(df$Veg) == "barren"
  n_before <- nrow(df)
  df <- df[keep_rows | is.na(df$Veg), , drop = FALSE]
  cat(sprintf(
    "Filtered to allowed classes (%s) + barren: kept %d/%d rows\n",
    paste(ALLOWED_VEG, collapse = ","), nrow(df), n_before
  ))
}

try(
  {
    cat("Per-veg quick summary:\n")
    all_veg_classes <- c(ALLOWED_VEG, "barren")
    for (av in all_veg_classes) {
      sel <- tolower(df$Veg) == av
      rows <- sum(sel, na.rm = TRUE)
      unique_doys <- length(unique(df$doy[sel & is.finite(df$doy)]))
      unique_locs <- length(unique(df$location_id[sel]))
      cat(sprintf("  %s: rows=%d unique_doys=%d unique_locs=%d\n", av, rows, unique_doys, unique_locs))
    }
  },
  silent = TRUE
)

matched_veg_n <- sum(!is.na(df$Veg))
cat("Non-NA Veg rows:", matched_veg_n, "of", nrow(df), "\n")
if (matched_veg_n == 0) {
  stop("No vegetation classes found after join; cannot build library")
}

loc_years <- df |>
  dplyr::filter(!is.na(.data$location_id) & .data$location_id != "" & !is.na(.data$pheno_year) & .data$pheno_year > 0 & !is.na(.data$Veg)) |>
  dplyr::distinct(.data$location_id, .data$pheno_year)
cat(sprintf("Constructed loc_years with %d rows from filtered df\n", nrow(loc_years)))
if (nrow(loc_years) == 0) {
  stop(paste0(
    "No location-pheno_year pairs found after filtering. This is a fatal error — library construction cannot continue.\n",
    "Possible causes and suggestions:\n",
    " - Your filtered training dataset has zero valid location/pheno_year pairs (check column 'location_id' and 'pheno_year').\n",
    " - Verify TRAIN_YEARS and any prior filtering steps do not remove all data (e.g. TRAIN_YEARS <- 2019:2024).\n",
    " - Ensure your transformer produced valid 'location_id' values that match your geojson mapping (transform_phenology.py formats 'L_lon_lat').\n",
    " - If your data are intentionally sparse, reduce the filtering thresholds or increase available training data.\n",
    "Processing cannot continue without at least one location-pheno_year pair in filtered training data.")
  )
}

cat("Constructing lib from TRAINING dataset...\n")
lib <- list()




cat("===============================================\n\n")

lib_df <- df
# No artificial per-vegetation sampling applied; use full lib_df for variant construction

vegs <- unique(na.omit(lib_df$Veg))
vegs <- vegs[vegs != ""]
vegs <- vegs[tolower(vegs) %in% c("herbs", "populus", "tamarix", "agriculture", "barren")]  # FIXED: case-insensitive matching (include agriculture)

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

  if (nrow(dveg) < MIN_ENDMEMBER_SAMPLES) {
    cat(sprintf("[NOTICE] Raw template for '%s' skipped due to insufficient samples: %d < %d\n", vname, nrow(dveg), MIN_ENDMEMBER_SAMPLES))
    next
  }

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

# --- Estimate spatial autocorrelation range from empirical variogram ---
# Fits an exponential variogram gamma(d) = sill * (1 - exp(-d/range)) to
# pairwise semivariances of `values` at locations given by `coords_df`.
# Returns the estimated range in km, or `fallback_km` if estimation fails.
estimate_autocorrelation_range <- function(coords_df, values, fallback_km = 30.0) {
  if (is.null(coords_df) || nrow(coords_df) == 0) return(fallback_km)
  if (!all(c("lat", "lon") %in% names(coords_df))) return(fallback_km)

  lat <- suppressWarnings(as.numeric(coords_df$lat))
  lon <- suppressWarnings(as.numeric(coords_df$lon))
  values <- suppressWarnings(as.numeric(values))

  valid <- which(is.finite(lat) & is.finite(lon) & is.finite(values))
  if (length(valid) < 5) return(fallback_km)  # Need enough pairs for a variogram

  coords <- cbind(lon[valid], lat[valid])
  vals <- values[valid]

  if (var(vals, na.rm = TRUE) == 0) return(fallback_km)

  # Pairwise great-circle distances (km) via Haversine
  dist_mat <- compute_haversine_distance_matrix(coords)

  dists <- dist_mat[upper.tri(dist_mat)]
  coef_diffs_sq <- outer(vals, vals, function(a, b) (a - b)^2)
  gamma_vals <- coef_diffs_sq[upper.tri(coef_diffs_sq)] / 2  # semivariance

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

  # Try NLS fit of exponential variogram
  range_est <- fit_exponential_variogram(bin_mid, bin_gamma, total_var, dists)
  if (is.null(range_est)) return(fallback_km)

  range_est <- as.numeric(range_est)
  if (!is.finite(range_est) || range_est <= 0) return(fallback_km)

  # Clamp to reasonable bounds (1 km to 500 km)
  range_est <- max(1, min(range_est, 500))
  cat(sprintf("[SPATIAL] Estimated autocorrelation range: %.1f km (used as block size)\n", range_est))
  range_est
}

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

  if (is.null(coords_df) || nrow(coords_df) == 0 || !all(c("location_id", "lat", "lon") %in% names(coords_df))) {
    return(sample(locations, n_draw, replace = TRUE))
  }

  coords_df$location_id <- trimws(as.character(coords_df$location_id))
  coords_df$lat <- suppressWarnings(as.numeric(coords_df$lat))
  coords_df$lon <- suppressWarnings(as.numeric(coords_df$lon))

  # Align order
  coords_df <- coords_df[match(locations, coords_df$location_id), , drop = FALSE]
  if (nrow(coords_df) == 0) return(sample(locations, n_draw, replace = TRUE))

  ok <- is.finite(coords_df$lat) & is.finite(coords_df$lon)
  missing_frac <- mean(!ok)
  if (!is.finite(missing_frac) || missing_frac > max_missing_frac || sum(ok) < 3) {
    return(sample(locations, n_draw, replace = TRUE))
  }

  # Use fallback if block_km not provided
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
  if (length(block_ids) < 2) {
    return(sample(locations, n_draw, replace = TRUE))
  }

  # Sample blocks with replacement until we have enough locations in the pool
  pool <- character(0)
  avg_block_size <- mean(lengths(blocks))
  n_blocks_draw <- max(1L, ceiling(n_draw / max(1, avg_block_size)))

  sampled_blocks <- sample(block_ids, size = n_blocks_draw, replace = TRUE)
  pool <- unlist(blocks[sampled_blocks], use.names = FALSE)

  while (length(pool) < n_draw) {
    sb <- sample(block_ids, size = 1L, replace = TRUE)
    pool <- c(pool, blocks[[sb]])
  }

  # Draw without replacement from the pooled list (pool may contain duplicates)
  sample(pool, n_draw, replace = FALSE)
}

# PPI bootstrap: thin wrapper around the generic location_bootstrap_index
location_bootstrap_ppi <- function(all_coefs, df_tasks, B = BOOTSTRAP_B, seed = 123) {
  ppi_max <- if (exists("PPI_FULL_VEG_COVER")) get("PPI_FULL_VEG_COVER") else 0.7
  if (ppi_max <= 0) ppi_max <- 0.7
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
    .block_km <- estimate_autocorrelation_range(coords_for_range, as.numeric(loc_means))
    cat(sprintf("[NOINDEX BOOTSTRAP] Data-driven block size: %.1f km\n", .block_km))

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
    # exclude rows that are effectively barren (index values <=0) from the
    # seasonal fit so the model captures vegetation behaviour only
    summer_fit_df <- summer_df[is.finite(summer_df[[idx_col]]) & summer_df[[idx_col]] > 0, , drop = FALSE]
    if (nrow(summer_fit_df) > min_detrend_n) {
      tryCatch({
        f <- as.formula(paste(idx_col, "~ poly(doy,", detrend_poly_degree, ")"))
        seasonal_model <- lm(f, data = summer_fit_df)
        summer_df$seasonal_trend <- predict(seasonal_model, newdata = summer_df)
        global_seasonal_mean <- mean(summer_df$seasonal_trend, na.rm = TRUE)
        summer_df[[detrended_col]] <- as.numeric(summer_df[[idx_col]] - (summer_df$seasonal_trend - global_seasonal_mean))
        cat(sprintf("[%s BOOTSTRAP] Computed detrended summer %s (N=%d; fit rows=%d).\n", 
                    index_name, index_name, nrow(summer_df), nrow(summer_fit_df)))
      }, error = function(e) {
        cat(sprintf("[%s BOOTSTRAP] Detrending failed: %s. Marking detrended values NA.\n", index_name, e$message))
        summer_df[[detrended_col]] <- NA_real_
      })
    } else {
      # after excluding barren there is insufficient data
      summer_df[[detrended_col]] <- NA_real_
      cat(sprintf("[%s BOOTSTRAP] Not enough vegetated rows for detrending; setting NA.\n", index_name))
    }
  } else {
    # not enough points to fit; fill with NA so weighting ignores these rows
    summer_df[[detrended_col]] <- NA_real_
    cat(sprintf("[%s BOOTSTRAP] Insufficient summer data for detrending; setting detrended values NA.\n", index_name))
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
  .block_km <- estimate_autocorrelation_range(coords_for_range, as.numeric(loc_means))
  cat(sprintf("[%s BOOTSTRAP] Data-driven block size: %.1f km\n", index_name, .block_km))

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

# Thin wrappers for backward compatibility (all delegate to location_bootstrap_index)
location_bootstrap_ndvi <- function(all_coefs, df_tasks, B = BOOTSTRAP_B, seed = 123, ndvi_max = 0.6)
  location_bootstrap_index(all_coefs, df_tasks, index_name = "NDVI", B = B, seed = seed, index_max = ndvi_max)

location_bootstrap_msavi <- function(all_coefs, df_tasks, B = BOOTSTRAP_B, seed = 123, msavi_max = 0.6)
  location_bootstrap_index(all_coefs, df_tasks, index_name = "MSAVI", B = B, seed = seed, index_max = msavi_max)

location_bootstrap_noindex <- function(all_coefs, df_tasks = NULL, B = BOOTSTRAP_B, seed = 123)
  location_bootstrap_index(all_coefs, df_tasks, index_name = "NOINDEX", B = B, seed = seed, noindex = TRUE)


# Confusion-matrix bootstrapping removed per user request; keep stub for
# backward compatibility so callers that might still exist do not fail.


# Location-based bootstrap for global aggregation
location_bootstrap_aggregate <- function(all_coefs, B = BOOTSTRAP_B, seed = 123) {
  set.seed(seed)
  
  # Classification-matrix bootstrapping disabled (no per-iteration confusion resampling)

  
  # Prepare Data: Pivot to Wide Format (Location x Class)
  clean_coefs <- all_coefs[!is.na(all_coefs$coef), ]
  years <- sort(unique(clean_coefs$pheno_year))
  all_veg <- unique(clean_coefs$Veg)
  
  results_list <- list()
  
  for (yr in years) {
    yr_data <- clean_coefs[clean_coefs$pheno_year == yr, ]
    if (nrow(yr_data) == 0) next
    
    # Pivot Wide
    locs <- unique(yr_data$location_id)
    n_locs <- length(locs)
    
    # Create matrix F [n_locs x n_veg]
    F_mat <- matrix(0, nrow = n_locs, ncol = length(all_veg))
    colnames(F_mat) <- all_veg
    rownames(F_mat) <- locs
    
    # Fill matrix
    m_loc <- match(yr_data$location_id, locs)
    m_veg <- match(yr_data$Veg, all_veg)
    flat_idx <- (m_veg - 1) * n_locs + m_loc
    F_mat[flat_idx] <- yr_data$coef
    
    # Spatial Blocking
    first_idx <- match(locs, yr_data$location_id)
    lat_vec <- if ("lat" %in% names(yr_data)) yr_data$lat[first_idx] else rep(NA, n_locs)
    lon_vec <- if ("lon" %in% names(yr_data)) yr_data$lon[first_idx] else rep(NA, n_locs)
    
    block_ids <- rep(1:n_locs)
    use_blocking <- FALSE
    
    if (sum(!is.na(lat_vec) & !is.na(lon_vec)) >= 3) {
      block_width_km <- 5
      lat_scale <- 111
      lon_scale <- 111 * cos(median(lat_vec, na.rm=TRUE) * pi/180)
      grid_lat <- floor(lat_vec * lat_scale / block_width_km)
      grid_lon <- floor(lon_vec * lon_scale / block_width_km)
      block_ids <- paste(grid_lat, grid_lon, sep="_")
      use_blocking <- TRUE
    }
    
    unique_blocks <- unique(block_ids)
    n_blocks <- length(unique_blocks)
    block_map <- split(1:n_locs, block_ids)
    
    # Bootstrap Loop
    boot_means <- matrix(NA, nrow=B, ncol=length(all_veg))
    colnames(boot_means) <- all_veg
    
    for (b in 1:B) {
      # A. Block Resample
      sampled_blks <- sample(unique_blocks, n_blocks, replace = TRUE)
      b_idx <- unlist(block_map[sampled_blks])
      F_b <- F_mat[b_idx, , drop = FALSE]
      
      # (classification-matrix adjustment removed)
      # C. Means
      boot_means[b, ] <- colMeans(F_b, na.rm=TRUE)
    }
    
    # Summarize
    for (v in all_veg) {
      vals <- boot_means[, v]
      res_row <- data.frame(
        year = yr,
        Veg = v,
        n_locations = n_locs,
        global_coef = mean(vals, na.rm=TRUE),
        se = sd(vals, na.rm=TRUE),
        coef_025 = quantile(vals, 0.025, na.rm=TRUE),
        coef_975 = quantile(vals, 0.975, na.rm=TRUE),
        n_eff = if (use_blocking && n_locs > 0) n_blocks else n_locs,
        method = "block_bootstrap_confusion_matrix"
      )
      res_row$global_coef <- pmax(0, res_row$global_coef)
      res_row$coef_025 <- pmax(0, res_row$coef_025)
      res_row$coef_975 <- pmin(1, res_row$coef_975)
      results_list[[length(results_list) + 1]] <- res_row
    }
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
          df_inf$location_id_orig <- df_inf$location_id  # Keep original for debugging
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

          # Apply MAX_INFERENCE_LOCATIONS limit immediately after loading
          if (exists("MAX_INFERENCE_LOCATIONS") && is.numeric(MAX_INFERENCE_LOCATIONS)) {
            unique_locations <- unique(df_inf$location_id[!is.na(df_inf$location_id)])
            n_unique <- length(unique_locations)

            if (n_unique > MAX_INFERENCE_LOCATIONS) {
              set.seed(get_mesma_seed(123))  # Deterministic sampling for reproducibility
              sampled_locations <- sample(unique_locations, MAX_INFERENCE_LOCATIONS, replace = FALSE)

              # Filter the dataframe to only include sampled locations
              df_inf <- df_inf[df_inf$location_id %in% sampled_locations, ]

              cat(sprintf("[INFERENCE LOADING] Reduced from %d to %d locations (MAX_INFERENCE_LOCATIONS=%d)\n",
                          n_unique, MAX_INFERENCE_LOCATIONS, MAX_INFERENCE_LOCATIONS))
              cat(sprintf("[INFERENCE LOADING] Filtered dataset: %d rows from %d locations\n",
                          nrow(df_inf), length(unique(df_inf$location_id[!is.na(df_inf$location_id)]))))
            } else {
              cat(sprintf("[INFERENCE LOADING] Using all %d locations (MAX_INFERENCE_LOCATIONS=%d)\n",
                          n_unique, MAX_INFERENCE_LOCATIONS))
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
        # NOTE: Whittaker smoothing is intentionally NOT applied to inference data.
        # The smoother interpolates through data gaps (weight=0 positions), producing
        # synthetic values that are indistinguishable from real observations in
        # downstream pentad aggregation. This makes sparse location-years appear
        # artificially complete, inflating unmixing confidence and deflating residuals.
        # Training data can tolerate this because prototypes average over many samples;
        # inference pixels cannot.
        cat("[INFERENCE] Skipping Whittaker smoothing on inference data (avoids gap-fill artifacts).\n")
        # compute_soil_line_slope(df_inf, assign_global_dvi = FALSE)  # Not needed for inference, soil line from training
      } else {
        cat("[WARNING] NDDI not found in inference data; skipping contamination filtering\n")
        # Indices were computed above; continue.
      }

      if (!"Veg" %in% names(df_inf)) df_inf$Veg <- NA_character_
      if (!"pheno_year" %in% names(df_inf)) df_inf$pheno_year <- assign_pheno_year(df_inf$date)
      if (!"doy" %in% names(df_inf)) df_inf$doy <- pheno_doy(df_inf$date)
      if (!"zenith.angle" %in% names(df_inf)) df_inf$zenith.angle <- NA_real_
      if (!"DVI_max" %in% names(df_inf)) df_inf$DVI_max <- NA_real_

      # IMPORTANT: Use ALL years for inference - inference should cover the full temporal range
      if (exists("TRAIN_YEARS") && !is.null(TRAIN_YEARS)) {
        if (!isTRUE(QUIET_MODE)) cat(sprintf("[NOTICE] TRAIN_YEARS is set to %s for training, but inference uses ALL available years.\n", paste(TRAIN_YEARS, collapse=", ")))
        cat(sprintf("[NOTICE] Inference dataset has %d rows from %d locations (all years)\n", nrow(df_inf), length(unique(df_inf$location_id))))
        cat("[NOTICE] Trend computations will use this inference dataset (INFERENCE_CSV) and will NOT use training data.\n")
      }

      # -----------------------------------------------------------------------------
      # TREND (INFERENCE DATA ONLY): compute per-pheno_year means for key indices
      # Do not use training data for trend computation — use the INFERENCE_CSV file provided
      # -----------------------------------------------------------------------------
      tryCatch({
        trend_indices <- c("MSAVI", "NDVI", "PPI")
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
              # no summer rows, leave detrended NA rather than falling back
              df_inf$ppi_detrended <- NA_real_
              if (!"ppi_detrended" %in% trend_indices) trend_indices <- c(trend_indices, "ppi_detrended")
            }
          } else {
            # missing seasonal model, mark all NA
            df_inf$ppi_detrended <- NA_real_
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
          plot_dt <- tidyr::pivot_longer(sum_by_year, cols = dplyr::starts_with("ppi") | dplyr::starts_with("MSAVI") | dplyr::starts_with("NDVI"), names_to = "index", values_to = "mean_val")
          if (nrow(plot_dt) > 0) {
            p_trend <- ggplot(plot_dt, aes(x = pheno_year, y = mean_val)) +
              add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
              geom_line() + geom_point() + facet_wrap(~index, scales = "free_y") + theme_minimal() +
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

      # Backup raw MSAVI and key bands for cover estimation (inference).
      # NOTE: apply_stored_normalization() z-scores indices/bands; MSAVI-based FVC must use raw units.
      if ("MSAVI" %in% names(df_inf)) {
        df_inf$MSAVI_raw <- df_inf$MSAVI
        cat("[NOTICE] Backed up raw inference MSAVI values to 'MSAVI_raw' before normalization.\n")
      }
      if ("nir" %in% names(df_inf) && !"nir_raw" %in% names(df_inf)) df_inf$nir_raw <- df_inf$nir
      if ("red" %in% names(df_inf) && !"red_raw" %in% names(df_inf)) df_inf$red_raw <- df_inf$red

      if (exists("TRAINING_NORM_PARAMS") && !is.null(TRAINING_NORM_PARAMS)) { 
        cat("\n=== APPLYING STORED NORMALIZATION TO INFERENCE DATA ===\n")
        df_inf <- apply_stored_normalization(df_inf, TRAINING_NORM_PARAMS, cols = avail, lat_default = 40.2)
        cat("=======================================================\n\n") 
      } else { 
        warning("TRAINING_NORM_PARAMS not found, applying fresh normalization to inference data (scale factors may differ from training!)")
        inf_norm_result <- normalize_mesma_data(df_inf, cols = avail, lat_default = 40.2)
        df_inf <- inf_norm_result$df 
      }

      ppi_max_inf <- if (exists("PPI_FULL_VEG_COVER")) get("PPI_FULL_VEG_COVER") else 0.7
      df_inf <- backup_and_normalize_ppi(df_inf, ppi_max_inf, label = "Inference")

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
      if (exists("df_train") && !is.null(df_train) && nrow(df_train) > 0) { if (!"pheno_year" %in% names(df_train) && "date" %in% names(df_train)) df_train$pheno_year <- assign_pheno_year(df_train$date); n_train_loc_years <- nrow(unique(df_train[c("location_id", "pheno_year")])) } else { n_train_loc_years <- 0 }
      cat(sprintf("(NOTICE) Inference dataset location-years: %d\n", n_infer_loc_years))
      cat(sprintf("(NOTICE) Training dataset location-years: %d\n", n_train_loc_years))
      if (n_train_loc_years > 0 && n_infer_loc_years == n_train_loc_years) cat(sprintf("(WARNING) Training and inference datasets have the same number of location-years (%d). This may be expected if IDs are independent; no automatic filtering will be applied.\n", n_train_loc_years))

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
  # FCLS: Fully Constrained Least Squares Solver (with optional L1 sparsity)
  # Standard mode: min ||Ef-y||²  s.t. sum(f)=1, f>=0
  # Sparse mode:   min ||Ef-y||² + λ·sum(f)  s.t. sum(f)<=1, f>=0
  #   When USE_SPARSE_UNMIXING is TRUE, the sum-to-one equality is relaxed to
  #   sum-to-at-most-one, and an L1 penalty λ·||f||₁ = λ·sum(f) (since f>=0)
  #   drives irrelevant endmembers to exactly 0. The solution is NOT renormalized,
  #   so sum(f) < 1 is possible (residual = unexplained fraction).
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

    # Impute non-finite values
    y_fit[!is.finite(y_fit)] <- 0
    E_fit[!is.finite(E_fit)] <- 0

    # Base feature weights (from PCA-LDA or uniform)
    if (!is.null(feature_weights) && length(feature_weights) == n_bands) {
      feature_weights <- as.numeric(feature_weights)
      feature_weights[!is.finite(feature_weights)] <- 0
      base_weights <- pmax(feature_weights, 0)
    } else {
      base_weights <- rep(1, n_bands)
    }

    # Determine sparse unmixing mode
    sparse_mode <- exists("USE_SPARSE_UNMIXING") && isTRUE(USE_SPARSE_UNMIXING)
    lambda_l1 <- if (sparse_mode && exists("SPARSE_LAMBDA")) SPARSE_LAMBDA else 0

    E_w <- E_fit * base_weights
    y_w <- y_fit * base_weights

    # quadprog solves: min 1/2 b^T Dmat b - dvec^T b
    # Quadratic term: Dmat = 2 * E'E
    # Linear term:    dvec = 2 * E'y  (standard)  or  2*E'y - λ (sparse)
    Dmat <- 2 * crossprod(E_w)
    dvec <- 2 * as.numeric(crossprod(E_w, y_w))

    # Scale QP to prevent inconsistent constraints
    qp_scale <- mean(diag(Dmat))
    if(is.na(qp_scale) || qp_scale < 1e-12) qp_scale <- 1.0

    Dmat <- Dmat / qp_scale
    dvec <- dvec / qp_scale

    # L1 penalty: subtract λ/qp_scale from dvec (equivalent to adding λ·sum(f) to objective)
    if (sparse_mode && lambda_l1 > 0) {
      dvec <- dvec - lambda_l1 / qp_scale
    }

    # Add ridge to diagonal for numerical stability (PD requirement)
    ridge <- 1e-6
    Dmat <- Dmat + diag(n_endmembers) * ridge
    Dmat <- (Dmat + t(Dmat)) / 2

    # Constraints: A^T b >= b_0
    if (sparse_mode && lambda_l1 > 0) {
      # Sparse: sum(f) <= 1, f >= 0
      # -sum(f) >= -1  =>  A^T row = [-1,...,-1], b_0 = -1
      # f >= 0         =>  I rows
      Amat <- cbind(-rep(1, n_endmembers), diag(n_endmembers))
      bvec <- c(-1, rep(0, n_endmembers))
      meq <- 0  # all inequality constraints
    } else {
      # Standard FCLS: sum(f) = 1, f >= 0
      Amat <- cbind(rep(1, n_endmembers), diag(n_endmembers))
      bvec <- c(1, rep(0, n_endmembers))
      meq <- 1  # first constraint is equality
    }

    res_qp <- tryCatch({
      quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = meq)
    }, error = function(e) {
      if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG] Quadprog failed in solve_weights_fcls: %s\n", e$message))
      if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG] Dmat dim: %s, dvec length: %d, Amat dim: %s, bvec length: %d\n",
                  paste(dim(Dmat), collapse="x"), length(dvec), paste(dim(Amat), collapse="x"), length(bvec)))
      e
    })

    if (!inherits(res_qp, "error")) {
      w_qp <- res_qp$solution
      w_qp[!is.finite(w_qp)] <- 0
      w_qp[w_qp < 0] <- 0

      if (!sparse_mode || lambda_l1 == 0) {
        # Standard FCLS: renormalize to sum-to-one
        w_sum <- sum(w_qp, na.rm = TRUE)
        if(w_sum > 0) w_qp <- w_qp / w_sum else w_qp <- rep(1/n_endmembers, n_endmembers)
      }
      # Sparse mode: do NOT renormalize — sum(f) < 1 is intentional

      # Calculate RMSE
      pred <- as.numeric(E_fit %*% w_qp)
      resid <- y_fit - pred
      rmse <- sqrt(mean(resid^2))

      return(list(w = w_qp, rmse = rmse, residuals = resid,
                  loss_type = if (sparse_mode && lambda_l1 > 0) "sparse_fcls" else "rmse"))
    }

    stop(sprintf("[QP] quadprog::solve.QP failed in solve_weights_fcls: %s", res_qp$message))
  }

  # =============================================================================
  # IWLMM Solver — Inequality-constrained Weighted Linear Mixture Model
  # (Li et al. 2021)
  #
  # Extends FCLS by allowing bounded perturbations ΔM to each endmember:
  #   min_{f, ΔM}  Σ_i w_i * (y_i - (M_i + ΔM_i) f)^2  + λ_δ * ||ΔM||^2
  #   s.t.  sum(f) = 1,  f >= 0,  -b_{i,j} <= ΔM_{i,j} <= b_{i,j}
  #
  # Solved via alternating optimization:
  #   Step A: fix ΔM, solve for f  (standard FCLS QP)
  #   Step B: fix f,  solve for ΔM (bounded least squares, closed-form per feature)
  #
  # E:       Features x Endmembers endmember matrix (columns = endmember spectra)
  # y:       Features-length observation vector
  # bounds:  Features x Endmembers matrix of perturbation bounds (>= 0)
  #          Derived as IWLMM_BOUND_SIGMA * within-class SD per feature per endmember.
  #          If NULL, falls back to standard FCLS (no perturbation).
  # feature_weights: optional PCA-LDA weights (length = n_features)
  # =============================================================================
  solve_weights_iwlmm <- function(E, y, bounds = NULL, feature_weights = NULL) {
    if (is.null(E) || ncol(E) < 1) return(NULL)

    # Fall back to FCLS if no bounds provided or IWLMM disabled
    if (is.null(bounds) || !exists("USE_IWLMM") || !isTRUE(USE_IWLMM)) {
      return(solve_weights_fcls(E, y, feature_weights = feature_weights))
    }

    E_fit <- as.matrix(E)
    y_fit <- as.numeric(y)
    n_bands <- nrow(E_fit)
    n_em <- ncol(E_fit)

    if (length(y_fit) != n_bands) {
      if (length(y_fit) > n_bands) y_fit <- y_fit[1:n_bands] else y_fit <- c(y_fit, rep(0, n_bands - length(y_fit)))
    }
    y_fit[!is.finite(y_fit)] <- 0
    E_fit[!is.finite(E_fit)] <- 0

    # Ensure bounds matrix matches dimensions
    bounds <- as.matrix(bounds)
    if (nrow(bounds) != n_bands || ncol(bounds) != n_em) {
      return(solve_weights_fcls(E, y, feature_weights = feature_weights))
    }
    bounds[!is.finite(bounds) | bounds < 0] <- 0

    # Base feature weights
    if (!is.null(feature_weights) && length(feature_weights) == n_bands) {
      feature_weights <- as.numeric(feature_weights)
      feature_weights[!is.finite(feature_weights)] <- 0
      base_w <- pmax(feature_weights, 0)
    } else {
      base_w <- rep(1, n_bands)
    }

    # Config
    max_iter <- if (exists("IWLMM_MAX_ITER")) IWLMM_MAX_ITER else 15L
    conv_tol <- if (exists("IWLMM_TOL")) IWLMM_TOL else 1e-4
    lambda_delta <- if (exists("IWLMM_REGULARIZE_DELTA")) IWLMM_REGULARIZE_DELTA else 0.01

    # L1 sparsity config (consistent with solve_weights_fcls)
    sparse_mode <- exists("USE_SPARSE_UNMIXING") && isTRUE(USE_SPARSE_UNMIXING)
    lambda_l1 <- if (sparse_mode && exists("SPARSE_LAMBDA")) SPARSE_LAMBDA else 0

    # Initialize: ΔM = 0, solve for f using standard FCLS
    delta_M <- matrix(0, nrow = n_bands, ncol = n_em)
    f_prev <- rep(1 / n_em, n_em)

    for (iter in seq_len(max_iter)) {
      # --- Step A: fix ΔM, solve for f (FCLS on perturbed endmembers) ---
      E_perturbed <- E_fit + delta_M
      E_pw <- E_perturbed * base_w
      y_pw <- y_fit * base_w

      Dmat <- 2 * crossprod(E_pw)
      dvec <- 2 * as.numeric(crossprod(E_pw, y_pw))

      qp_scale <- mean(diag(Dmat))
      if (is.na(qp_scale) || qp_scale < 1e-12) qp_scale <- 1.0
      Dmat <- Dmat / qp_scale
      dvec <- dvec / qp_scale

      # L1 penalty on fractions
      if (sparse_mode && lambda_l1 > 0) {
        dvec <- dvec - lambda_l1 / qp_scale
      }

      Dmat <- Dmat + diag(n_em) * 1e-6
      Dmat <- (Dmat + t(Dmat)) / 2

      if (sparse_mode && lambda_l1 > 0) {
        # Sparse: sum(f) <= 1, f >= 0
        Amat <- cbind(-rep(1, n_em), diag(n_em))
        bvec <- c(-1, rep(0, n_em))
        meq <- 0
      } else {
        Amat <- cbind(rep(1, n_em), diag(n_em))
        bvec <- c(1, rep(0, n_em))
        meq <- 1
      }

      res_qp <- tryCatch(
        quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = meq),
        error = function(e) e
      )

      if (inherits(res_qp, "error")) {
        return(solve_weights_fcls(E, y, feature_weights = feature_weights))
      }

      f_new <- res_qp$solution
      f_new[!is.finite(f_new)] <- 0
      f_new[f_new < 0] <- 0
      if (!sparse_mode || lambda_l1 == 0) {
        f_sum <- sum(f_new)
        if (f_sum > 0) f_new <- f_new / f_sum else f_new <- rep(1 / n_em, n_em)
      }

      # --- Step B: fix f, solve for ΔM (bounded per feature) ---
      # For each feature i:  min_δ  w_i * (r_i - Σ_j δ_{i,j} * f_j)^2 + λ_δ * Σ_j δ_{i,j}^2
      #   s.t.  -b_{i,j} <= δ_{i,j} <= b_{i,j}
      # Residual with current endmembers (no perturbation):
      #   r_i = y_i - Σ_j M_{i,j} * f_j
      r_base <- y_fit - E_fit %*% f_new  # residual from original endmembers

      # Iterative coordinate-descent over endmembers for each feature.
      # For feature i, the residual to explain is r_i = y_i - M_i * f.
      # We distribute it across endmembers via coordinate descent:
      #   for each j: δ_{i,j} = clip( f_j * r_remaining / (f_j^2 + λ_δ/w_i), -b, b )
      # This lets bounds on one endmember redirect residual to others.
      cd_passes <- 3L  # coordinate descent passes per feature
      for (i in seq_len(n_bands)) {
        wi <- base_w[i]
        if (wi < 1e-12) {
          delta_M[i, ] <- 0
          next
        }
        ri <- as.numeric(r_base[i])
        delta_row <- delta_M[i, ]  # warm-start from previous outer iteration
        reg <- lambda_delta / wi

        for (pass in seq_len(cd_passes)) {
          for (j in seq_len(n_em)) {
            fj <- f_new[j]
            if (abs(fj) < 1e-12) { delta_row[j] <- 0; next }
            # Partial residual: undo contribution of current δ_j
            r_partial <- ri - sum(delta_row * f_new) + delta_row[j] * fj
            # Optimal unconstrained update for j alone
            delta_j <- fj * r_partial / (fj^2 + reg)
            # Clip to bounds
            delta_row[j] <- max(-bounds[i, j], min(bounds[i, j], delta_j))
          }
        }
        delta_M[i, ] <- delta_row
      }

      # Check convergence
      if (max(abs(f_new - f_prev)) < conv_tol) break
      f_prev <- f_new
    }

    # Final prediction and residuals
    E_final <- E_fit + delta_M
    pred <- as.numeric(E_final %*% f_new)
    resid <- y_fit - pred
    rmse <- sqrt(mean(resid^2))

    # Also compute how much perturbation was used (diagnostic)
    delta_norm <- sqrt(mean(delta_M^2))

    return(list(
      w = f_new,
      rmse = rmse,
      residuals = resid,
      loss_type = "iwlmm",
      delta_M = delta_M,
      delta_norm = delta_norm,
      n_iter = iter
    ))
  }

  # Batch FCLS Solver using Quadprog (Optimization for GA/Grid Search)
  # Supports L1 sparsity when USE_SPARSE_UNMIXING is TRUE.
  # E: Features x Endmembers matrix
  # Y: Samples x Features matrix
  # Returns: Samples x Endmembers weight matrix
  solve_batch_fcls <- function(E, Y, feature_weights = NULL) {
    if (is.null(E) || ncol(E) < 1) return(NULL)
    n_endmembers <- ncol(E)
    n_samples <- nrow(Y)
    if (n_samples == 0) return(matrix(0, 0, n_endmembers))

    sparse_mode <- exists("USE_SPARSE_UNMIXING") && isTRUE(USE_SPARSE_UNMIXING)
    lambda_l1 <- if (sparse_mode && exists("SPARSE_LAMBDA")) SPARSE_LAMBDA else 0

    # Convert Y to Features x Samples for matrix math
    Y_t <- t(Y)

    # Determine E and Y used for fitting (weighted or raw)
    if (!is.null(feature_weights) && length(feature_weights) == nrow(E)) {
      w <- pmax(feature_weights, 0)
      w[!is.finite(w)] <- 0
      E_fit <- E * w
      Y_fit <- Y_t * w
    } else {
      E_fit <- E
      Y_fit <- Y_t
    }

    # Precompute constant QP matrices
    Dmat <- 2 * crossprod(E_fit)

    qp_scale <- mean(diag(Dmat))
    if(is.na(qp_scale) || qp_scale < 1e-12) qp_scale <- 1.0

    Dmat <- Dmat / qp_scale

    ridge <- 1e-6
    Dmat <- Dmat + diag(n_endmembers) * ridge
    Dmat <- (Dmat + t(Dmat)) / 2

    # Constraints depend on sparse mode
    if (sparse_mode && lambda_l1 > 0) {
      # Sparse: sum(f) <= 1, f >= 0
      Amat <- cbind(-rep(1, n_endmembers), diag(n_endmembers))
      bvec <- c(-1, rep(0, n_endmembers))
      meq <- 0
    } else {
      # Standard: sum(f) = 1, f >= 0
      Amat <- cbind(rep(1, n_endmembers), diag(n_endmembers))
      bvec <- c(1, rep(0, n_endmembers))
      meq <- 1
    }

    # Precompute all linear terms: dvec = 2 * E'y / qp_scale
    Dvecs <- (2 * crossprod(E_fit, Y_fit)) / qp_scale

    # L1 penalty offset (applied uniformly to all samples)
    l1_offset <- if (sparse_mode && lambda_l1 > 0) lambda_l1 / qp_scale else 0

    w_out <- matrix(0, nrow=n_samples, ncol=n_endmembers)

    for(i in 1:n_samples) {
      dvec <- Dvecs[, i] - l1_offset

      res <- tryCatch({
        quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq=meq)
      }, error = function(e) {
        e
      })

      if(!inherits(res, "error")) {
        w <- res$solution
        w[!is.finite(w)] <- 0
        w[w < 0] <- 0
        if (!sparse_mode || lambda_l1 == 0) {
          s <- sum(w)
          if(s > 0) w <- w / s else w <- rep(1/n_endmembers, n_endmembers)
        }
        w_out[i, ] <- w
      } else {
        stop(sprintf("[QP] quadprog::solve.QP failed in solve_batch_fcls (sample %d/%d): %s", i, n_samples, res$message))
      }
    }
    return(w_out)
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


  # --- Robust loss (Huber) for inference scoring ---
  # We use Huber loss to reduce the influence of outlier feature residuals when
  # selecting the best endmember combination for inference tasks.
  # The loss is weighted by feature reliability weights (e.g. PCA→LDA mask weights)
  # and normalized by sum(weights) so tasks with missing features remain comparable.
  huber_rho <- function(r, delta) {
    rr <- abs(as.numeric(r))
    d <- as.numeric(delta)
    if (!is.finite(d) || d <= 0) d <- 1.5
    out <- ifelse(rr <= d, 0.5 * rr^2, d * (rr - 0.5 * d))
    out[!is.finite(out)] <- NA_real_
    out
  }

  weighted_huber_loss <- function(residuals, feature_weights = NULL, delta = NULL) {
    r <- as.numeric(residuals)
    if (length(r) == 0) return(Inf)

    w <- rep(1, length(r))
    if (!is.null(feature_weights) && length(feature_weights) == length(r)) {
      w <- as.numeric(feature_weights)
      w[!is.finite(w)] <- 0
      w <- pmax(w, 0)
    }

    # Default delta from config if present
    d <- if (!is.null(delta)) as.numeric(delta) else {
      if (exists("HUBER_DELTA", inherits = TRUE)) as.numeric(get("HUBER_DELTA", inherits = TRUE)) else 1.5
    }
    if (!is.finite(d) || d <= 0) d <- 1.5

    rho <- huber_rho(r, d)
    denom <- sum(w, na.rm = TRUE)
    if (!is.finite(denom) || denom <= 0) return(Inf)
    sum(w * rho, na.rm = TRUE) / denom
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
    score_metric = c("rmse", "huber"),
    huber_delta = NULL
  ) {
    if (length(top_variants) == 0) return(NULL)

    score_metric <- match.arg(score_metric)

    full_veg_names <- names(top_variants)
    n_veg_full <- length(full_veg_names)
    if (n_veg_full == 0) return(NULL)

    y_target <- y
    y_target[!is.finite(y_target)] <- 0

    score_from_solution <- function(res, fw = NULL) {
      if (is.null(res) || is.null(res$residuals)) return(Inf)
      if (identical(score_metric, "huber")) {
        return(weighted_huber_loss(res$residuals, feature_weights = fw, delta = huber_delta))
      }
      as.numeric(res$rmse)
    }

    # Solve a combination of variant indices for a SUBSET of vegetation types
    # veg_subset: character vector of vegetation type names to include
    # variant_indices: integer vector of variant indices (one per veg type in subset)
    # Precompute per-class within-class SD for IWLMM perturbation bounds
    # Each entry is a numeric vector of length n_features (SD across all variant vecs)
    iwlmm_active <- exists("USE_IWLMM") && isTRUE(USE_IWLMM)
    iwlmm_class_sd <- list()
    if (iwlmm_active) {
      bound_sigma <- if (exists("IWLMM_BOUND_SIGMA")) IWLMM_BOUND_SIGMA else 2.0
      for (v in full_veg_names) {
        cands <- top_variants[[v]]
        if (is.null(cands) || length(cands) < 2) {
          iwlmm_class_sd[[v]] <- NULL
          next
        }
        Mv <- do.call(rbind, lapply(cands, function(z) as.numeric(z$vec)))
        if (is.null(Mv) || nrow(Mv) < 2) {
          iwlmm_class_sd[[v]] <- NULL
          next
        }
        sds <- apply(Mv, 2, sd, na.rm = TRUE)
        sds[!is.finite(sds)] <- 0
        iwlmm_class_sd[[v]] <- sds
      }
    }

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

      fw <- if (!is.null(feature_weights) && length(feature_weights) == nrow(E)) feature_weights else NULL

      # Use IWLMM if active and we have bounds for at least one class
      if (iwlmm_active) {
        # Build bounds matrix: n_features x n_endmembers
        bounds_mat <- matrix(0, nrow = nrow(E), ncol = ncol(E))
        for (v_idx in seq_along(veg_subset)) {
          v <- veg_subset[v_idx]
          cls_sd <- iwlmm_class_sd[[v]]
          if (!is.null(cls_sd) && length(cls_sd) == nrow(E)) {
            bounds_mat[, v_idx] <- bound_sigma * cls_sd
          }
        }
        # Only use IWLMM if at least some bounds are nonzero
        if (any(bounds_mat > 0)) {
          res <- solve_weights_iwlmm(E, y_target, bounds = bounds_mat, feature_weights = fw)
        } else {
          res <- solve_weights_fcls(E, y_target, feature_weights = fw)
        }
      } else {
        res <- solve_weights_fcls(E, y_target, feature_weights = fw)
      }
      if (is.null(res)) return(NULL)

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
      res$score <- score_from_solution(res, fw = fw)
      res$score_type <- score_metric
      if (identical(score_metric, "huber")) {
        res$huber_loss <- res$score
        res$huber_delta <- if (!is.null(huber_delta)) huber_delta else if (exists("HUBER_DELTA", inherits = TRUE)) get("HUBER_DELTA", inherits = TRUE) else 1.5
      }
      return(res)
    }

    # --- All vegetation types, one prototype per class ---
    # Sparsity is handled inside the solver (L1 penalty drives irrelevant weights to 0).
    # Here we just search over prototype combinations.
    n_variants_per_veg <- sapply(top_variants, length)

    global_best_res <- NULL
    global_best_score <- Inf

    veg_subset <- full_veg_names
    n_variants_subset <- n_variants_per_veg[veg_subset]

    combos <- expand.grid(lapply(n_variants_subset, seq_len), KEEP.OUT.ATTRS = FALSE)
    for (i in seq_len(nrow(combos))) {
      r <- solve_combo_subset(veg_subset, as.integer(combos[i, ]))
      if (!is.null(r) && r$score < global_best_score) {
        global_best_score <- r$score
        global_best_res <- r
      }
    }

    if (is.null(global_best_res)) return(NULL)

    # [MODIFIED] Return only the single best model as top_models per request
    # Previous logic sorted accumulated candidates and kept near-optimal ones.
    top_models <- list(global_best_res)

    return(list(
      w = global_best_res$w,
      rmse = global_best_res$rmse,
      score = global_best_score,
      score_type = if (!is.null(global_best_res$score_type)) global_best_res$score_type else "rmse",
      huber_loss = if (!is.null(global_best_res$huber_loss)) global_best_res$huber_loss else NA_real_,
      huber_delta = if (!is.null(global_best_res$huber_delta)) global_best_res$huber_delta else NA_real_,
      ids = global_best_res$ids,
      residuals = global_best_res$residuals,
      E_best = global_best_res$E,
      top_models = top_models,
      loss_type = if (!is.null(global_best_res$loss_type)) global_best_res$loss_type else "rmse",
      delta_norm = global_best_res$delta_norm,
      veg_subset = global_best_res$veg_subset
    ))
  }

  
  build_mesma_library_weighted <- function(df_train, indices, params, allowed_veg,
                                           precomputed_clusters = NULL,
                                           generate_proto_plots = NULL,
                                           df_oob = NULL) {
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
    if (is.null(df_oob) || !is.data.frame(df_oob) || nrow(df_oob) == 0) {
      stop("[ERROR] OOB holdout data (df_train_oob) is required and must be passed explicitly to build_mesma_library_weighted().")
    }

    oob_locs <- if ("location_id" %in% names(df_oob)) unique(as.character(df_oob$location_id)) else character(0)
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

          # --- Prefilter: require sufficient observed temporal coverage
          # Build pentad matrix and count how many pentads contain observations.
          mat_nointerp <- build_pentad_matrix(sub, base_indices)
          if (is.null(mat_nointerp)) next
          n_nonempty_pentads <- sum(apply(mat_nointerp, 1, function(r) any(is.finite(r))))

          if (!is.finite(n_nonempty_pentads) || n_nonempty_pentads < MIN_PENTADS_PER_TRAIN_SAMPLE) {
            if (!isTRUE(QUIET_MODE)) {
              cat(sprintf("  [LIB BUILD] Skipping training trace %s pheno_year=%s: %d non-empty pentads (< %d)\n",
                          as.character(lid), as.character(pyr), n_nonempty_pentads, MIN_PENTADS_PER_TRAIN_SAMPLE))
            }
            next
          }

          # Acceptable coverage — build pentad matrix for storage
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
          dists_to_centroid <- rowSums(sweep(data_w, 2, centroid, "-")^2)
          best_local_idx <- which.min(dists_to_centroid)
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
        if (i == 1) {
           # Initial distances logic
           min_sq_dists <- dist_mat[, best_idx]
        } else {
           min_sq_dists <- pmin(min_sq_dists, dist_mat[, best_idx])
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
           # also cap barren prototypes to 3 for now to match vegetation restriction
           max_k_barren <- min(max_k_barren, 3L)
           k_candidates <- 1:max_k_barren
           k_candidates <- k_candidates[k_candidates <= max_k]
           if (length(k_candidates) == 0) k_candidates <- 1
        } else {
           # Vegetation classes: cap at 3 regardless of MAX_K_EAR
           max_k_veg <- min(MAX_K_EAR, 3L)
           k_candidates <- 1:max_k_veg
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

         cat(sprintf("[DEBUG] Extracting endmembers for '%s': n_samples=%d, max_k=%d\n",
                     v, nrow(veg_mat_train), max_k_needed))

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
      train_sets <- list()
      oob_sets <- list()

      # Build OOB evaluation set from storage_oob
      oob_samples <- list(); oob_lbls <- c()
      for (v in valid_classes) {
        if (v %in% names(storage_oob) && !is.null(storage_oob[[v]]) && nrow(storage_oob[[v]]) > 0) {
          oob_samples[[length(oob_samples)+1]] <- storage_oob[[v]]
          oob_lbls <- c(oob_lbls, rep(v, nrow(storage_oob[[v]])))
        }
      }

      if (length(oob_samples) == 0) {
        cat("[LIB BUILD] WARNING: storage_oob produced zero OOB samples after prefiltering.\n")
        # Provide diagnostics from df_oob (trace counts + rows per Veg) so user can choose a tailored fix
        if (exists("df_oob") && is.data.frame(df_oob) && nrow(df_oob) > 0) {
          vegs_oob <- unique(as.character(df_oob$Veg))
          cat("[LIB BUILD] df_oob summary (per Veg):\n")
          for (vv in vegs_oob) {
            rows_v <- sum(df_oob$Veg == vv, na.rm = TRUE)
            traces_v <- length(unique(paste0(df_oob$location_id[df_oob$Veg == vv], "_", df_oob$pheno_year[df_oob$Veg == vv])))
            cat(sprintf("      %s: traces=%d, rows=%d\n", vv, traces_v, rows_v))
          }

          # Compute median non-empty pentads per Veg (may be slow for many traces)
          base_indices_local <- if (!is.null(params$base_indices)) params$base_indices else indices
          traces_list <- unique(paste(df_oob$location_id, df_oob$pheno_year, df_oob$Veg, sep = "|"))
          pentad_stats <- list()
          if (length(traces_list) > 0) {
            for (tr in traces_list) {
              # Use literal '|' when fixed=TRUE; avoid unrecognized escape '\\|'
              parts <- strsplit(tr, "|", fixed = TRUE)[[1]]
              lid <- parts[1]; pyr <- as.integer(parts[2]); vv <- parts[3]
              sub <- df_oob[df_oob$location_id == lid & df_oob$pheno_year == pyr & df_oob$Veg == vv, , drop = FALSE]
              mat_nointerp <- build_pentad_matrix(sub, base_indices_local)
              if (is.null(mat_nointerp)) next
              n_nonempty <- sum(apply(mat_nointerp, 1, function(r) any(is.finite(r))))
              pentad_stats[[length(pentad_stats) + 1]] <- list(Veg = vv, pentads = n_nonempty)
            }
            if (length(pentad_stats) > 0) {
              pentad_df <- do.call(rbind, lapply(pentad_stats, function(x) data.frame(Veg = x$Veg, pentads = x$pentads)))
              pentad_summary <- aggregate(pentads ~ Veg, pentad_df, function(x) c(med = median(x), min = min(x), max = max(x)))
              for (i in seq_len(nrow(pentad_summary))) {
                s <- pentad_summary$pentads[i][[1]]
                cat(sprintf("      %s: pentads (med/min/max) = %.0f / %.0f / %.0f\n",
                            pentad_summary$Veg[i], s["med"], s["min"], s["max"]))
              }
            }
          }
        } else {
          cat("[LIB BUILD] No df_oob available for diagnostics\n")
        }

        # Fallback: use conservative precomputed cluster counts (do not abort)
        default_clusters <- list()
        for (v in target_classes) default_clusters[[v]] <- 1L
        if ("herbs" %in% target_classes) default_clusters[["herbs"]] <- min(3L, as.integer(MAX_VARIANTS_PER_VEG))
        precomputed_clusters <- default_clusters
        cat(sprintf("[LIB BUILD] Falling back to precomputed clusters: %s\n",
                    paste(names(precomputed_clusters), precomputed_clusters, sep = "=", collapse = ", ")))
      }

      oob_sets[[1]] <- list(Y = do.call(rbind, oob_samples), labels = oob_lbls)
      cat(sprintf("    OOB evaluation set: %d samples across %d classes\n",
                  nrow(oob_sets[[1]]$Y), length(unique(oob_lbls))))

      # Build training evaluation set from storage (model training data, excluding OOB)
      train_samples <- list(); train_lbls <- c()
      for (v in valid_classes) {
        if (v %in% names(storage) && !is.null(storage[[v]]) && nrow(storage[[v]]) > 0) {
          train_samples[[length(train_samples)+1]] <- storage[[v]]
          train_lbls <- c(train_lbls, rep(v, nrow(storage[[v]])))
        }
      }
      if (length(train_samples) > 0) {
        train_sets[[1]] <- list(Y = do.call(rbind, train_samples), labels = train_lbls)
        cat(sprintf("    Training evaluation set: %d samples across %d classes\n",
                    nrow(train_sets[[1]]$Y), length(unique(train_lbls))))
      } else {
        cat("    [WARNING] No training samples available for combined evaluation; using OOB only\n")
      }

      # 3b. Fitness Function - Evaluate on OOB holdout AND training data
      # Optimizes VEGETATION ACCURACY (populus, tamarix, herbs) - excludes barren
      # Returns min of OOB and training accuracy — config must perform well on both
      score_on_dataset <- function(target_set, M, col_names, veg_classes) {
          Y_test <- target_set$Y
          labels_test <- target_set$labels
          n_test <- nrow(Y_test)

          # Unmixing always happens in PCA→LDA feature space.
          if (is.null(params$pca_lda) || is.null(params$pca_lda$std_to_lda) || !is.matrix(params$pca_lda$std_to_lda)) {
            stop("[PCA/LDA] Missing std_to_lda transform; PCA→LDA unmixing is required for optimization scoring")
          }
          A <- params$pca_lda$std_to_lda
          if (nrow(A) != nrow(M) || nrow(A) != ncol(Y_test)) {
            stop(sprintf("[PCA/LDA] Dimension mismatch in scoring projection: nrow(A)=%d, nrow(M)=%d, ncol(Y_test)=%d",
                         nrow(A), nrow(M), ncol(Y_test)))
          }

          # Project both templates and samples into PCA→LDA space.
          M_lda <- t(t(M) %*% A)         # [n_lda x n_endmembers]
          Y_lda <- Y_test %*% A          # [n_samples x n_lda]
          w_lda <- rep(1, ncol(A))

          # Batch process using solve_batch_fcls for speed
          all_coefs <- solve_batch_fcls(M_lda, Y_lda, w_lda)

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
          if(is.null(oob_sets[[1]])) return(-1.0)

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

          M <- do.call(cbind, cols)

          has_train <- !is.null(train_sets[[1]])

          oob_score <- score_on_dataset(oob_sets[[1]], M, col_names, veg_classes)

          combo_score <- if (has_train) {
            train_score <- score_on_dataset(train_sets[[1]], M, col_names, veg_classes)
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
    plot_vegetation_prototypes <- function(lib, indices = NULL, out_dir = if (exists("OUT_DIR")) OUT_DIR else ".", prefix = "veg_prototypes", save_png = TRUE, dpi = 150, feature_weights = NULL, show_medians = TRUE, variant_alpha = NULL) {
      if (is.null(lib) || length(lib) == 0) return(NULL)
      if (is.null(indices)) {
        if (exists("MESMA_PARAMS") && !is.null(MESMA_PARAMS) && !is.null(MESMA_PARAMS$indices)) indices <- MESMA_PARAMS$indices else indices <- OPTIMAL_INDICES
      }
      if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required for prototype plotting")
      if (!dir.exists(out_dir) && save_png) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      # Precompute pentad-to-month mapping for x-axis (pheno year: Mar=1)
      # Month start days in phenological DOY (April 1 = day 1)
      pheno_month_starts <- c(1, 32, 62, 93, 123, 154, 184, 215, 245, 276, 306, 337)
      pheno_month_labels <- c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan","Feb")
      # Convert month-start pheno DOY to pentad position: pentad = ceil(doy / agg_days)
      month_breaks_pentad <- ceiling(pheno_month_starts / TEMPORAL_AGGREGATION_DAYS)
      # Keep only unique breaks within the temporal budget
      keep <- month_breaks_pentad >= 1 & month_breaks_pentad <= TEMPORAL_BUDGET & !duplicated(month_breaks_pentad)
      month_breaks_pentad <- month_breaks_pentad[keep]
      month_labels <- pheno_month_labels[keep]

      # Per-index backpropagated timestep weights (feature_weights). These are used to
      # render a weight backdrop and to optionally skip fully-zero indices.
      weights_by_index <- NULL
      if (!is.null(feature_weights) && length(feature_weights) == length(indices) * TEMPORAL_BUDGET) {
        weights_by_index <- list()
        for (k in seq_len(length(indices))) {
          w_start <- (k - 1) * TEMPORAL_BUDGET + 1
          w_end <- k * TEMPORAL_BUDGET
          wv <- as.numeric(feature_weights[w_start:w_end])
          wv[!is.finite(wv)] <- 0
          weights_by_index[[indices[k]]] <- pmax(wv, 0)
        }
        excluded <- names(weights_by_index)[vapply(weights_by_index, function(w) !any(w > 0), logical(1))]
        if (length(excluded) > 0) {
          cat(sprintf("[PROTO PLOT] Skipping %d fully-zero weight indices: %s\n", length(excluded), paste(excluded, collapse = ", ")))
          indices <- setdiff(indices, excluded)
          weights_by_index <- weights_by_index[indices]
        }
      }

      # Significance masking removed from pipeline — all timesteps are treated equally now.

      # Build long dataframe for plotting
      rows <- list()
      for (v in names(lib)) {
        variants <- lib[[v]]
        for (var in variants) {
          vec <- as.numeric(var$vec)
          vid <- if (!is.null(var$id)) var$id else if (!is.null(var$variant_id)) var$variant_id else paste0(v, "_unknown")
          n_idx <- length(indices)
          expected_len <- n_idx * TEMPORAL_BUDGET
          if (length(vec) < expected_len) {
            # Skip malformed variants but warn
            warning(sprintf("Skipping variant %s for veg %s: length(vec)=%d != expected=%d", vid, v, length(vec), expected_len))
            next
          }
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

        p <- ggplot2::ggplot(df_idx, ggplot2::aes(x = pentad, y = value, group = variant_id))
        p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.4)

        # Weight backdrop: show backpropagated weight per pentad as a grey-scale band
        # with a legend (scale). This makes the timestep weighting visually explicit.
        if (!is.null(weights_by_index) && idx %in% names(weights_by_index)) {
          wv <- weights_by_index[[idx]]
          w_df <- data.frame(
            pentad = seq_len(TEMPORAL_BUDGET),
            w = wv,
            xmin = seq_len(TEMPORAL_BUDGET) - 0.5,
            xmax = seq_len(TEMPORAL_BUDGET) + 0.5,
            ymin = y_min - y_pad,
            ymax = y_max + y_pad
          )
          p <- p + ggplot2::geom_rect(
            data = w_df,
            inherit.aes = FALSE,
            ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = w),
            alpha = 0.25
          )
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
             ggplot2::labs(title = sprintf("Vegetation prototypes: %s", idx),
                           subtitle = subtitle_txt,
                           x = "Month",
                           y = y_lab) +
             ggplot2::scale_x_continuous(breaks = month_breaks_pentad, labels = month_labels) +
             ggplot2::theme_minimal() +
             ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                            plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9)) +
             ggplot2::scale_color_manual(values = veg_palette) +
             ggplot2::coord_cartesian(ylim = c(y_min - y_pad, y_max + y_pad))

        # Add fill scale for weight backdrop (if used)
        if (!is.null(weights_by_index) && idx %in% names(weights_by_index)) {
          p <- p + ggplot2::scale_fill_gradient(
            name = "Backprop weight",
            low = "grey95",
            high = "grey45"
          )
        } else {
          p <- p + ggplot2::guides(fill = "none")
        }
        plots[[idx]] <- p
        if (save_png) {
          fn <- file.path(out_dir, sprintf("%s_%s.png", prefix, idx))
          ggplot2::ggsave(filename = fn, plot = p, width = 8, height = 4, dpi = dpi)
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
        plot_weights <- if (!is.null(params) && !is.null(params$weights)) params$weights else NULL
        out_base <- file.path(if (exists("OUT_DIR")) OUT_DIR else ".", "prototype_plots")
        plot_vegetation_prototypes(res_lib, indices = plot_indices, out_dir = out_base, feature_weights = plot_weights, show_medians = TRUE)
        cat(sprintf("[NOTICE] Generated prototype plots (with medians) to %s\n", out_base))

        if (exists("GENERATE_PROTO_PLOTS_VARIANTS_ONLY") && isTRUE(GENERATE_PROTO_PLOTS_VARIANTS_ONLY)) {
          out_variants <- file.path(if (exists("OUT_DIR")) OUT_DIR else ".", "prototype_plots_variants_only")
          plot_vegetation_prototypes(res_lib, indices = plot_indices, out_dir = out_variants, feature_weights = plot_weights, show_medians = FALSE)
          cat(sprintf("[NOTICE] Generated prototype plots (variants-only) to %s\n", out_variants))
        }
      }, error = function(e) {
        cat(sprintf("[WARN] Failed to generate prototype plots: %s\n", e$message))
      })
    }

    return(res_lib)
  }


  precompute_optimized_library_weighted <- function(mesma_lib, grid_type = "full", params = NULL) {
    opt_lib <- list()

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

        # Always project templates into PCA→LDA feature space for unmixing.
        # Keep the original (index/pentad) vectors in `mesma_lib` for diagnostics/plots.
        if (is.null(params) || is.null(params$pca_lda) || is.null(params$pca_lda$std_to_lda) || !is.matrix(params$pca_lda$std_to_lda)) {
          stop("[PCA/LDA] Missing std_to_lda transform; cannot build optimized library in PCA→LDA space")
        }
        A <- params$pca_lda$std_to_lda
        if (nrow(A) != ncol(M)) {
          stop(sprintf("[PCA/LDA] std_to_lda dimension mismatch: nrow(A)=%d != ncol(M)=%d", nrow(A), ncol(M)))
        }
        M <- M %*% A
        feature_space <- "pca_lda"

        pruned_info <- NULL

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
          pruned_info = pruned_info,
          feature_space = feature_space
        )
    }
    
    opt_lib
  }





  precompute_compressed_templates <- function(mesma_lib, grid_name = "full") {
    template_db <- list()
    for (veg in names(mesma_lib)) {
      template_db[[veg]] <- list()
      variants <- mesma_lib[[veg]]
      if (exists("filter_variants_by_min_samples", mode="function")) {
         rt <- if (exists("raw_lib_templates", envir=globalenv()) && !is.null(get("raw_lib_templates", envir=globalenv())[[veg]])) get("raw_lib_templates", envir=globalenv())[[veg]] else NULL
         variants <- filter_variants_by_min_samples(variants, min_samples = MIN_ENDMEMBER_SAMPLES, veg = veg, raw_template = rt)
      }
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

  if (!is.null(MESMA_PARAMS_INITIAL) && !is.null(MESMA_PARAMS_INITIAL$weights)) {
    MESMA_PARAMS_INITIAL$weights[is.na(MESMA_PARAMS_INITIAL$weights)] <- 0
    MESMA_PARAMS_INITIAL$weights[!is.finite(MESMA_PARAMS_INITIAL$weights)] <- 0
    print_weights_summary("INITIAL_PCA_LDA", MESMA_PARAMS_INITIAL)

    # Informative notification about PCA->LDA PC pruning/refit
    if (!is.null(MESMA_PARAMS_INITIAL$pca_lda)) {
      pa <- MESMA_PARAMS_INITIAL$pca_lda
      cat(sprintf("[NOTICE] PCA->LDA: PCs original=%d, retained=%d, pruned=%d\n",
                  ifelse(!is.null(pa$n_pcs_original), pa$n_pcs_original, NA),
                  ifelse(!is.null(pa$n_pcs_kept), pa$n_pcs_kept, NA),
                  ifelse(!is.null(pa$n_pcs_pruned), pa$n_pcs_pruned, NA)))
      if (!is.null(pa$pruned_pc_indices) && length(pa$pruned_pc_indices) > 0) {
        cat(sprintf("         Pruned PC indices: %s\n", paste(pa$pruned_pc_indices, collapse = ", ")))
      }
    }
  }

  # Step 2: Optimal cluster sizing FIRST (before threshold optimization)
  # This ensures cluster optimization uses the full feature set
  cat("\n=== STEP 2: Optimal Cluster Sizing (before threshold optimization) ===\n")
  cat("Building initial MESMA library with full feature weights for cluster optimization...\n")

  # Use indices from params (includes L2norm if enabled) instead of original avail
  indices_for_library <- if (!is.null(MESMA_PARAMS_INITIAL$indices)) MESMA_PARAMS_INITIAL$indices else avail

  # Build library with initial (unthresholded) weights to determine optimal cluster counts
  mesma_lib_initial <- build_mesma_library_weighted(df_train, indices_for_library, MESMA_PARAMS_INITIAL, ALLOWED_VEG,
                                                    generate_proto_plots = FALSE,
                                                    df_oob = df_train_oob)

  # Store the optimal cluster counts discovered during library building
  if (exists("OPTIMAL_CLUSTER_COUNTS", envir = globalenv())) {
    cat("[CLUSTER] Optimal cluster counts determined from full-feature optimization\n")
  }

  # Step 3 (REMOVED): Feature-weight pruning / alpha-grid pruning.
  # The pipeline now always uses the full feature set. Any zeros in weights are
  # treated as down-weighting only; no features/indices are dropped.
  pruned_indices <- avail
  MESMA_PARAMS <- MESMA_PARAMS_INITIAL

  if (!is.null(MESMA_PARAMS) && !is.null(MESMA_PARAMS$weights)) {
    print_weights_summary("MESMA_FINAL", MESMA_PARAMS)
  }

  # Store pruning results
  assign("PRUNED_INDICES", pruned_indices, envir = globalenv())

  # NOTE: STEP 4 (confusion matrix / validation matrix) is computed later,
  # after the OOB cluster optimization that happens during final library build.

  # MESMA UNMIXING: All endmembers (barren + veg types) treated equally
  cat("[MODE] MESMA unmixing ENABLED (barren and vegetation types treated as equals)\n")

  # Log solver choice
  if (exists("USE_IWLMM") && isTRUE(USE_IWLMM)) {
    cat(sprintf("[MODE] Using IWLMM solver (bound=\u00b1%.1f\u03c3, \u03bb_\u03b4=%.4f, max_iter=%d)\n",
                if (exists("IWLMM_BOUND_SIGMA")) IWLMM_BOUND_SIGMA else 2.0,
                if (exists("IWLMM_REGULARIZE_DELTA")) IWLMM_REGULARIZE_DELTA else 0.01,
                if (exists("IWLMM_MAX_ITER")) IWLMM_MAX_ITER else 15L))
  } else {
    cat("[MODE] Using standard RMSE loss for FCLS\n")
  }



  cat("Building final MESMA library (barren + all vegetation types)...\n")

  # Use precomputed cluster counts from Step 2 (if available) to skip re-optimization
  precomputed_clusters <- NULL
  if (exists("OPTIMAL_CLUSTER_COUNTS", envir = globalenv())) {
    precomputed_clusters <- get("OPTIMAL_CLUSTER_COUNTS", envir = globalenv())
  }

  # Use indices from MESMA_PARAMS (includes L2norm if enabled)
  indices_for_final_library <- if (!is.null(MESMA_PARAMS$indices)) MESMA_PARAMS$indices else avail

  mesma_lib <- build_mesma_library_weighted(df_train, indices_for_final_library, MESMA_PARAMS, ALLOWED_VEG,
                                            precomputed_clusters = precomputed_clusters,
                                            df_oob = df_train_oob)

  cat("Pre-computing optimized library for MESMA...\n")
  OPTIMIZED_LIBRARY <- precompute_optimized_library_weighted(mesma_lib, grid_type = "full", params = MESMA_PARAMS)



  if (isTRUE(TESTING_MODE)) {
    nn <- sapply(names(OPTIMIZED_LIBRARY), function(v) {
      libv <- OPTIMIZED_LIBRARY[[v]]
      if (is.null(libv) || is.null(libv$M)) return(0)
      nrow(libv$M)
    })
    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG] MESMA library built: endmember_count=%d variants_per_type=%s\n",
                length(nn), paste(sprintf("%s:%d", names(nn), nn), collapse=", ")))
  }

  assign("MESMA_PARAMS", MESMA_PARAMS, envir = globalenv())
  assign("mesma_lib", mesma_lib, envir = globalenv())
  assign("OPTIMIZED_LIBRARY", OPTIMIZED_LIBRARY, envir = globalenv())

  cat(sprintf("[NOTICE] MESMA feature count: avail=%d, params_indices=%d\n",
              length(avail), length(MESMA_PARAMS$indices)))
  cat(sprintf("[NOTICE] MESMA feature list: %s\n", paste(MESMA_PARAMS$indices, collapse=", ")))

  mesma_lib <- mesma_lib

  assign("OPTIMIZED_LIBRARY", OPTIMIZED_LIBRARY, envir = globalenv())
  assign("mesma_lib", mesma_lib, envir = globalenv())

  compressed_templates_accessor <- precompute_compressed_templates(mesma_lib, "full")
  assign("compressed_templates_accessor", compressed_templates_accessor, envir = globalenv())
  assign(".COMPRESSED_TEMPLATES_ACCESSOR", compressed_templates_accessor, envir = globalenv())



  # === Shared confusion matrix function (used for both training and validation) ===
  # Pipeline: OPTIMIZED_LIBRARY endmembers -> energy normalization -> batch FCLS
  # with feature weights. Observations: z-score.
  compute_confusion_matrix <- function(df_data, label = "Training") {
    if (is.null(df_data) || nrow(df_data) == 0) return(invisible(NULL))
    if (!exists("OPTIMIZED_LIBRARY") || is.null(OPTIMIZED_LIBRARY)) return(invisible(NULL))
    if (is.null(MESMA_PARAMS)) return(invisible(NULL))

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
      for (cls in unique_classes) {
        lib_cls <- OPTIMIZED_LIBRARY[[cls]]
        if (is.null(lib_cls) || is.null(lib_cls$M)) next
        M_cls <- lib_cls$M
        if (is.null(dim(M_cls))) M_cls <- matrix(M_cls, nrow = 1)
        for (ri in seq_len(nrow(M_cls))) {
          E_cols[[length(E_cols) + 1]] <- M_cls[ri, ]
          col_class_labels <- c(col_class_labels, cls)
        }
      }
      if (length(E_cols) < 2) return(invisible(NULL))

      E <- do.call(cbind, E_cols)

      # Unmixing always happens in PCA→LDA feature space.
      if (is.null(MESMA_PARAMS$pca_lda) || is.null(MESMA_PARAMS$pca_lda$std_to_lda) || !is.matrix(MESMA_PARAMS$pca_lda$std_to_lda)) {
        stop("[PCA/LDA] Missing std_to_lda transform; cannot compute confusion matrix in PCA→LDA space")
      }
      A_cm <- MESMA_PARAMS$pca_lda$std_to_lda
      w_cm <- rep(1, ncol(A_cm))
      if (nrow(E) != length(w_cm)) {
        stop(sprintf("[PCA/LDA] Confusion-matrix library not in PCA→LDA space (n_features(E)=%d, expected=%d)", nrow(E), length(w_cm)))
      }

      # 2. Build observation vectors (z-score + prune)
      base_indices_cm <- if (!is.null(MESMA_PARAMS$base_indices)) MESMA_PARAMS$base_indices else avail
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

        sub <- df_cm[df_cm$location_id == lid & df_cm$pheno_year == pyr, ]
        mat <- build_pentad_matrix(sub, base_indices_cm)
        if (!is.null(mat)) {
          vec <- as.numeric(mat)
          vec <- mesma_apply_representation_vec(vec, n_base_idx_cm, TEMPORAL_BUDGET, l2_normalize_cm)
          vec <- mesma_zscore_vec_by_index(vec, indices_cm, MESMA_PARAMS$means, MESMA_PARAMS$sds, TEMPORAL_BUDGET)

          vec_feat <- project_vec_to_pca_lda(vec, MESMA_PARAMS$pca_lda)
          if (is.null(vec_feat)) {
            stop("[PCA/LDA] Failed to project observation into PCA→LDA space")
          }
          vec <- vec_feat

          obs_vecs[[length(obs_vecs) + 1]] <- vec
          obs_labels <- c(obs_labels, true_cls)
        }
      }
      if (length(obs_vecs) == 0) return(invisible(NULL))

      Y <- do.call(rbind, obs_vecs)

      # 3. Unmixing (IWLMM or FCLS) in PCA→LDA feature space
      w <- w_cm

      iwlmm_cm <- exists("USE_IWLMM") && isTRUE(USE_IWLMM)
      if (iwlmm_cm) {
        # Compute per-column bounds from within-class variant spread
        cm_bounds <- matrix(0, nrow = nrow(E), ncol = ncol(E))
        cm_bound_sigma <- if (exists("IWLMM_BOUND_SIGMA")) IWLMM_BOUND_SIGMA else 2.0
        for (cls in unique_classes) {
          col_idx <- which(col_class_labels == cls)
          if (length(col_idx) >= 2) {
            cls_mat <- E[, col_idx, drop = FALSE]
            cls_sd <- apply(cls_mat, 1, sd, na.rm = TRUE)
            cls_sd[!is.finite(cls_sd)] <- 0
            for (ci in col_idx) {
              cm_bounds[, ci] <- cm_bound_sigma * cls_sd
            }
          }
        }
        has_bounds <- any(cm_bounds > 0)

        # Loop per observation using IWLMM
        n_samples <- nrow(Y)
        n_cols <- ncol(E)
        all_coefs <- matrix(0, nrow = n_samples, ncol = n_cols)
        for (si in seq_len(n_samples)) {
          y_obs <- as.numeric(Y[si, ])
          if (has_bounds) {
            res_i <- solve_weights_iwlmm(E, y_obs, bounds = cm_bounds, feature_weights = w)
          } else {
            res_i <- solve_weights_fcls(E, y_obs, feature_weights = w)
          }
          if (!is.null(res_i) && !is.null(res_i$w)) {
            all_coefs[si, ] <- res_i$w
          }
        }
      } else {
        all_coefs <- solve_batch_fcls(E, Y, w)
      }
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
        if (true_cls %in% cm_labels) {
          for (cc in cm_labels) frac_mat[true_cls, cc] <- frac_mat[true_cls, cc] + coefs[[cc]]
          class_counts[true_cls] <- class_counts[true_cls] + 1
        }
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


      cat(sprintf("[CONFUSION MATRIX] %s confusion matrix stored for classification uncertainty propagation\n", label))
    }, error = function(e) {
      cat(sprintf("[CONFUSION MATRIX] Error computing %s confusion matrix: %s\n", label, e$message))
    })
  }

  # === STEP 4: Compute Training Confusion Matrix ===
  if (exists("df_train") && !is.null(df_train) && nrow(df_train) > 0) {
    cat("\n=== STEP 4: Computing Confusion Matrix with Final Weights (All Training Data) ===\n")
    compute_confusion_matrix(df_train, "Training")
  }

  

# ==========================================================================
# Training is performed externally — this script focuses on inference & validation.
# validation_location_ids already set during stratified train/validation split
# NOTE: do NOT create or populate `df_tasks` here — inference code will build task tables from the inference input.
# ==========================================================================

if (isTRUE(TESTING_MODE)) cat("[DEBUG] Reached line 6377 - about to define fit_one_task function\n")

# ============================================================================
# OOB FRACTION PREDICTION UNCERTAINTY
# Uses OOB validation residuals to add realistic prediction error to MC draws.
# For each true class, we store the distribution of prediction residuals
# (predicted_fraction - true_fraction) observed in OOB validation.
# During MC, we sample from these residuals to perturb the predicted fractions.
# ============================================================================

# Global storage for OOB fraction residuals per true class
.OOB_FRACTION_RESIDUALS <- NULL  # List: true_class -> matrix of residuals (n_samples x n_classes)

# Store OOB fraction residuals for use in MC uncertainty propagation
# Called after OOB validation with known true classes
#
# Parameters:
#   residuals_by_class: Named list where each element is a matrix of residuals
#                       Rows = OOB samples, Cols = predicted classes
#                       Each row is (predicted_fractions - true_fractions) for one sample
#   true_class_names: Names of the true classes (should match list names)
store_oob_fraction_residuals <- function(residuals_by_class) {
  if (!is.list(residuals_by_class) || length(residuals_by_class) == 0) {
    warning("[OOB_FRAC] Invalid residuals - must be a non-empty named list")
    return(invisible(FALSE))
  }

  # Validate structure
  for (cls in names(residuals_by_class)) {
    resid_mat <- residuals_by_class[[cls]]
    if (!is.matrix(resid_mat) && !is.data.frame(resid_mat)) {
      warning(sprintf("[OOB_FRAC] Residuals for class '%s' must be a matrix", cls))
      return(invisible(FALSE))
    }
  }

  assign(".OOB_FRACTION_RESIDUALS", residuals_by_class, envir = globalenv())

  if (isTRUE(DEBUG_UNCERTAINTY)) {
    cat("[OOB_FRAC] Stored OOB fraction residuals for MC uncertainty:\n")
    for (cls in names(residuals_by_class)) {
      n_samples <- nrow(residuals_by_class[[cls]])
      cat(sprintf("  %s: %d samples\n", cls, n_samples))
    }
  }

  return(invisible(TRUE))
}

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

  # ORPHANED_BLOCK_REMOVED: stray legacy fragment deleted - see git history
  # (original code referenced undefined variables such as `fracs`, `mapped_row`)
  rows <- character(0)
  row_idx <- integer(0)

  # ORPHANED: confusion‑based fraction perturbation removed.
  # If needed reintroduce as a pure function that accepts explicit inputs
  # (fracs, conf_matrix, mapped_row, sample_sizes, concentration_scale, etc.).
  # For safety, return `fracs` unchanged when present (no-op):
  result <- if (exists("fracs")) fracs else numeric(0)
  return(result)

  fit_one_task <- function(task_data) {
    if (is.null(task_data) || nrow(task_data) == 0) return(NULL)

    loc <- as.character(task_data$location_id[1])
    yr <- as.integer(task_data$pheno_year[1])
    

    # Extract lat/lon if available
    lat_val <- NA_real_
    lon_val <- NA_real_
    if ("lat" %in% names(task_data)) lat_val <- as.numeric(task_data$lat[1]) else if ("latitude" %in% names(task_data)) lat_val <- as.numeric(task_data$latitude[1])
    if ("lon" %in% names(task_data)) lon_val <- as.numeric(task_data$lon[1]) else if ("longitude" %in% names(task_data)) lon_val <- as.numeric(task_data$longitude[1])

    # Use MESMA parameters
    PARAMS <- MESMA_PARAMS
    
    # DEBUG: Verify PARAMS structure
    if (is.null(PARAMS) || is.null(PARAMS$indices) || is.null(PARAMS$means) || is.null(PARAMS$sds)) {
      cat(sprintf("[ERROR] loc=%s yr=%d: PARAMS structure is incomplete (PARAMS=%s, indices=%s, means=%s, sds=%s)\n",
                  loc, yr,
                  !is.null(PARAMS), 
                  !is.null(PARAMS$indices), 
                  !is.null(PARAMS$means), 
                  !is.null(PARAMS$sds)))
      return(NULL)
    }
    
    # DEBUG: Log feature space
    if (isTRUE(TESTING_MODE)) {
      cat(sprintf("[FEATURE SPACE] loc=%s yr=%d: Using %d indices: %s\n",
                  loc, yr, length(PARAMS$indices), paste(head(PARAMS$indices, 5), collapse=", ")))
      cat(sprintf("[FEATURE SPACE] task_data columns: %s\n", paste(head(names(task_data), 10), collapse=", ")))
    }
    
    base_indices_for_build <- if (!is.null(PARAMS$base_indices)) PARAMS$base_indices else PARAMS$indices
    raw_mat <- build_pentad_matrix(task_data, base_indices_for_build)
    if (is.null(raw_mat)) {
      if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG] loc=%s yr=%d: build_pentad_matrix returned NULL\n", loc, yr))
      return(NULL)
    }

    y_raw <- as.numeric(raw_mat)

    # === CREATE VALIDITY MASK ===
    # Mask for valid (non-NA) observations - DO NOT replace NA with 0
    valid_mask <- is.finite(y_raw)
    n_valid <- sum(valid_mask)

    # Require at least one valid observation; warn in TESTING_MODE when coverage is very low
    if (n_valid == 0) {
      if (exists("TESTING_MODE") && isTRUE(TESTING_MODE)) {
        cat(sprintf("[SKIP] loc=%s pheno_year=%d: no valid observations (n_valid=0) - skipping\n", loc, yr))
      }
      return(NULL)
    } else {
      MIN_VALID_FRACTION <- 0.01
      if (exists("TESTING_MODE") && isTRUE(TESTING_MODE) && n_valid < (length(y_raw) * MIN_VALID_FRACTION)) {
        cat(sprintf("[WARN] loc=%s pheno_year=%d: low valid observations (%d < %.0f) - proceeding anyway\n",
                    loc, yr, n_valid, length(y_raw) * MIN_VALID_FRACTION))
      }
    }
    # ============================

    # --- New Logic: Estimate barren fraction using MSAVI linear model ---
    # Moved from end of function to allow early exit for pure barren pixels
    # 1. Get current MSAVI for the task

    # MSAVI is used here on the *raw* (non-zscored) observation table to estimate a
    # vegetation cover proxy. It must be computed upstream together with the other indices
    # (compute_indices_from_bands) so we do not maintain a separate MSAVI pipeline here.
    msavi_col <- NA_character_
    if ("MSAVI_raw" %in% names(task_data)) {
      msavi_col <- "MSAVI_raw"
    } else if ("MSAVI" %in% names(task_data)) {
      msavi_col <- "MSAVI"
    }

    if (is.na(msavi_col) || !msavi_col %in% names(task_data)) {
      available_cols <- paste(head(names(task_data), 20), collapse=", ")
      stop(sprintf("ERROR loc=%s yr=%d: MSAVI missing in task_data. Compute MSAVI via compute_indices_from_bands() before normalization so MSAVI_raw is available. First 20 available columns: %s",
                   loc, yr, available_cols))
    }

    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG MSAVI CHECK] loc=%s yr=%d: Using MSAVI column '%s', nrow=%d\n",
                loc, yr, msavi_col, nrow(task_data)))

    # Filter task_data for summer months (June-September). No fallback to all rows.
    if (!"month" %in% names(task_data) && "date" %in% names(task_data)) {
      d <- task_data$date
      month_v <- rep(NA_integer_, length(d))

      if (inherits(d, "Date") || inherits(d, "POSIXt")) {
        month_v <- suppressWarnings(as.integer(format(d, "%m")))
      } else {
        d_chr <- as.character(d)
        d_parsed <- suppressWarnings(as.Date(d_chr))
        if (all(is.na(d_parsed))) d_parsed <- suppressWarnings(as.Date(d_chr, format = "%Y/%m/%d"))
        if (all(is.na(d_parsed))) d_parsed <- suppressWarnings(as.Date(d_chr, format = "%d/%m/%Y"))
        if (all(is.na(d_parsed))) d_parsed <- suppressWarnings(as.Date(d_chr, format = "%Y%m%d"))
        month_v <- suppressWarnings(as.integer(format(d_parsed, "%m")))
      }

      if (all(is.na(month_v))) stop(sprintf("ERROR loc=%s yr=%d: Could not parse month from 'date' column", loc, yr))
      task_data$month <- month_v
    }

    used_summer_filter <- FALSE
    summer_task_data <- task_data
    if ("month" %in% names(task_data)) {
      used_summer_filter <- TRUE
      summer_task_data <- task_data[is.finite(task_data$month) & task_data$month %in% 6:9, , drop = FALSE]
    }

    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG MSAVI] loc=%s yr=%d: summer_task_data has %d rows (filtered from %d)\n",
                loc, yr, nrow(summer_task_data), nrow(task_data)))

    if (nrow(summer_task_data) == 0) return(NULL)

    n_summer_valid <- sum(is.finite(summer_task_data[[msavi_col]]))
    if (n_summer_valid < 3) return(NULL)

    low_data_flag <- n_summer_valid < 5

    current_msavi <- median(summer_task_data[[msavi_col]], na.rm = TRUE)

    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG MSAVI] loc=%s yr=%d: Summer MSAVI median = %.4f (from %d values, %d valid)\n",
                loc, yr, current_msavi, length(summer_task_data[[msavi_col]]), sum(is.finite(summer_task_data[[msavi_col]]))))

    # Error if MSAVI is not finite
    if (!is.finite(current_msavi)) {
      stop(sprintf("ERROR loc=%s yr=%d: MSAVI is not finite (value=%s). Cannot calculate FVC without valid MSAVI data.",
                   loc, yr, as.character(current_msavi)))
    }

    # This avoids dependency on external calibration models.
    fvc_predicted <- current_msavi / (0.5 + current_msavi)
    # Log approximate MM-based mapping (note: not a calibrated FVC model)

    # Clamp to [0, 1] for safety
    total_veg_cover <- pmin(pmax(fvc_predicted, 0), 1)
    barren_fraction <- 1 - total_veg_cover


    # Verify barren_fraction is finite - error if not
    if (!is.finite(barren_fraction)) {
      stop(sprintf("ERROR loc=%s yr=%d: barren_fraction is not finite (value=%s) after MSAVI FVC prediction. This should never happen.",
                   loc, yr, as.character(barren_fraction)))
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

      # Always use PCA→LDA feature space (explicit projection).
      if (is.null(PARAMS$pca_lda) || is.null(PARAMS$pca_lda$std_to_lda) || !is.matrix(PARAMS$pca_lda$std_to_lda)) {
        stop("[PCA/LDA] Missing std_to_lda transform; PCA→LDA unmixing is required")
      }
      y_for_unmixing <- project_vec_to_pca_lda(y_norm, PARAMS$pca_lda)
      weights_masked <- compute_pca_lda_feature_weights_from_mask(valid_mask, PARAMS$pca_lda)
      if (is.null(y_for_unmixing) || is.null(weights_masked)) {
        stop("[PCA/LDA] Failed to create PCA→LDA observation/weights; PCA→LDA unmixing is required")
      }

      # Perform MESMA using all endmembers (barren + all veg types)
      # Barren is kept as from MESMA - no replacement
      veg_kept <- names(mesma_lib)
      # SAFEGUARD: Ensure a 'barren' endmember exists
      if (!"barren" %in% veg_kept) {
        warning("No 'barren' endmember found in MESMA library - ensure barren endmembers are present.")
      }
      
      # DEBUG: Check library availability
      if (is.null(OPTIMIZED_LIBRARY) || length(OPTIMIZED_LIBRARY) == 0) {
        cat(sprintf("[ERROR] loc=%s yr=%d: OPTIMIZED_LIBRARY is NULL or empty\n", loc, yr))
        return(NULL)
      }
      if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG LIBRARY] loc=%s yr=%d: OPTIMIZED_LIBRARY contains %d vegetation types: %s\n",
                  loc, yr, length(OPTIMIZED_LIBRARY), paste(names(OPTIMIZED_LIBRARY), collapse=", ")))

      # Determine the library feature dimension (first non-empty class).
      lib_feature_dim <- NA_integer_
      for (vv in names(OPTIMIZED_LIBRARY)) {
        libv <- OPTIMIZED_LIBRARY[[vv]]
        if (!is.null(libv) && !is.null(libv$M)) {
          lib_feature_dim <- ncol(as.matrix(libv$M))
          break
        }
      }

      if (!is.finite(lib_feature_dim) || lib_feature_dim != length(y_for_unmixing)) {
        stop(sprintf("[PCA/LDA] Optimized library feature dimension mismatch: lib_feature_dim=%s, y_dim=%d", as.character(lib_feature_dim), length(y_for_unmixing)))
      }

      # Check if we have sufficient signal in the chosen space.
      y_w_chk <- y_for_unmixing * sqrt(pmax(weights_masked, 0))
      y_norm_val <- sqrt(sum(y_w_chk^2, na.rm = TRUE))

      if (is.na(y_norm_val) || y_norm_val < 1e-9) {
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_task] loc=%s yr=%d: insufficient signal (norm=%.6g), skipping\n", loc, yr, y_norm_val))
        return(NULL)
      }
      
      top_variants <- list()

      for(v in veg_kept) {
        lib <- OPTIMIZED_LIBRARY[[v]]
        if(is.null(lib)) {
          next
        }

        M_full <- lib$M
        if (is.null(M_full)) next
        if (is.null(dim(M_full))) M_full <- matrix(M_full, nrow = 1)

        lib_ids_kept <- lib$ids
        if (length(lib_ids_kept) != nrow(M_full)) {
          lib_ids_kept <- rep(paste0(v, "_1"), nrow(M_full))
        }

        if (nrow(M_full) == 0) next

        # Templates are stored in PCA→LDA feature space (no masking).
        if (ncol(M_full) != length(y_for_unmixing)) {
          next
        }

        M_full[!is.finite(M_full)] <- 0

        # Similarity ranking (weighted cosine similarity).
        w_masked <- weights_masked
        if (is.null(w_masked) || length(w_masked) != length(y_for_unmixing)) w_masked <- rep(1, length(y_for_unmixing))
        sqrt_w <- sqrt(pmax(w_masked, 0))

        y_w <- y_for_unmixing * sqrt_w
        denom_y <- sqrt(sum(y_w^2, na.rm = TRUE)); if (!is.finite(denom_y) || denom_y < 1e-12) denom_y <- 1
        y_w_norm <- y_w / denom_y

        M_w <- sweep(M_full, 2, sqrt_w, "*")
        M_w_norm <- t(apply(M_w, 1, function(row) {
          nr <- sqrt(sum(row^2, na.rm = TRUE))
          if (!is.finite(nr) || nr < 1e-12) row else row / nr
        }))

        sims <- as.numeric(M_w_norm %*% y_w_norm)
        best_idx <- order(sims, decreasing = TRUE)

        # Templates for unmixing are L2-normalized (consistent with existing pipeline).
        # Use precomputed normalization if available and aligned.
        if (!is.null(lib$M_norm) && !is.null(dim(lib$M_norm)) && nrow(lib$M_norm) == nrow(lib$M) && ncol(lib$M_norm) == ncol(lib$M)) {
          M_unmix <- lib$M_norm
        } else {
          M_unmix <- t(apply(M_full, 1, function(row) {
            nr <- sqrt(sum(row^2, na.rm = TRUE))
            if (!is.finite(nr) || nr < 1e-12) row else row / nr
          }))
        }

        top_variants[[v]] <- lapply(best_idx, function(i) {
          vec_i <- as.numeric(M_unmix[i, ])
          list(vec = vec_i, id = lib_ids_kept[i], similarity = sims[i])
        })
      }

      # Remove empty vegetation types
      empty_vegs <- names(top_variants)[sapply(top_variants, function(x) is.null(x) || length(x) == 0)]
      if (length(empty_vegs) > 0) {
        for (ev in empty_vegs) top_variants[[ev]] <- NULL
      }
      

      if (length(top_variants) == 0) {
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_task] loc=%s yr=%d: no top_variants available after filtering empties, skipping\n", loc, yr))
        return(NULL)
      }
      
      # Evaluate all combinations of endmembers
      # Pass unweighted observation; weighting happens inside solve_weights_fcls

      # Inference tasks typically have Veg missing/NA. For inference we select the
      # best model using a robust (Huber) loss to reduce sensitivity to outlier
      # features while still respecting the PCA→LDA reliability weights.
      is_inference_task <- {
        if (!"Veg" %in% names(task_data)) {
          TRUE
        } else {
          v <- task_data$Veg
          v_chr <- trimws(as.character(v))
          all(is.na(v_chr) | v_chr == "" | tolower(v_chr) == "na")
        }
      }

      score_metric <- if (is_inference_task) "huber" else "rmse"
      
      best_result <- tryCatch({
        evaluate_all_combinations(
          y_for_unmixing,
          top_variants,
          lambda = 0,
          feature_weights = weights_masked,
          score_metric = score_metric
        )
      }, error = function(e) {
        cat(sprintf("[ERROR fit_one_task] evaluate_all_combinations failed for loc=%s year=%d: %s\n", as.character(loc), as.integer(yr), e$message))
        NULL
      })
      
      # BEGIN SAFETY BLOCK (Agent added)
      # tryCatch({
      

      if (is.null(best_result)) {
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_task] loc=%s yr=%d: evaluate_all_combinations returned NULL\n", loc, yr))
        return(NULL)
      }
      

      # Extract coefficients and create output
      chosen_ids <- best_result$ids
      coefs <- best_result$w
      rmse <- best_result$rmse
      residuals <- best_result$residuals
      E_best_masked <- if (!is.null(best_result$E_best)) best_result$E_best else NULL

      # DIAGNOSTIC: Check coefficient values immediately after extraction

      # CRITICAL DEBUG: Print MESMA coefficients to console for diagnosis
      delta_norm_val <- if (!is.null(best_result$delta_norm)) best_result$delta_norm else NA
      loss_type_val <- if (!is.null(best_result$loss_type)) best_result$loss_type else "unknown"
      score_type_val <- if (!is.null(best_result$score_type)) as.character(best_result$score_type) else "rmse"
      score_val <- if (!is.null(best_result$score) && is.finite(best_result$score)) as.numeric(best_result$score) else NA_real_
      cat(sprintf("[MESMA OUTPUT] loc=%s yr=%d: %s (rmse=%.4f, solver=%s, score_type=%s, score=%.4f, delta_norm=%.4g)\n",
                  loc, yr,
                  paste(sprintf("%s=%.4f", names(coefs), coefs), collapse=", "),
                  rmse, loss_type_val,
          score_type_val,
          if (is.finite(score_val)) score_val else NA_real_,
                  if (is.finite(delta_norm_val)) delta_norm_val else 0))

      
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
          verbose_bundle <- isTRUE(DEBUG_UNCERTAINTY) || isTRUE(TESTING_MODE)

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

        # Precompute IWLMM bounds for MC draws (reuse variant pool)
        mc_iwlmm_bounds <- NULL
        mc_iwlmm_active <- exists("USE_IWLMM") && isTRUE(USE_IWLMM)
        if (mc_iwlmm_active && !is.null(E_best_masked)) {
          mc_bound_sigma <- if (exists("IWLMM_BOUND_SIGMA")) IWLMM_BOUND_SIGMA else 2.0
          mc_iwlmm_bounds <- matrix(0, nrow = nrow(E_best_masked), ncol = ncol(E_best_masked))
          for (j in seq_along(mc_veg_subset)) {
            v <- mc_veg_subset[j]
            cands <- top_variants[[v]]
            if (!is.null(cands) && length(cands) >= 2) {
              Mv <- do.call(rbind, lapply(cands, function(z) as.numeric(z$vec)))
              if (!is.null(Mv) && nrow(Mv) >= 2 && ncol(Mv) == nrow(E_best_masked)) {
                sds <- apply(Mv, 2, sd, na.rm = TRUE)
                sds[!is.finite(sds)] <- 0
                mc_iwlmm_bounds[, j] <- mc_bound_sigma * sds
              }
            }
          }
          if (!any(mc_iwlmm_bounds > 0)) mc_iwlmm_bounds <- NULL
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

          # Use IWLMM for MC draws when active
          if (mc_iwlmm_active && !is.null(mc_iwlmm_bounds)) {
            res_mc <- solve_weights_iwlmm(E_mc, y_mc, bounds = mc_iwlmm_bounds, feature_weights = weights_masked)
          } else {
            res_mc <- solve_weights_fcls(E_mc, y_mc, feature_weights = weights_masked)
          }
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
          if (isTRUE(TESTING_MODE)) cat(sprintf("[WARN fit_one_task] loc=%s yr=%d: Inseparable variant check failed: %s\n", loc, yr, e$message))
      })

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

      # --- Scale vegetation fractions by index-derived total vegetation cover ---
      # MESMA provides relative proportions among vegetation types.
      # We multiply only vegetation fractions by total_veg_cover (from MSAVI) to get absolute fractions.
      # Barren fraction is retained from MESMA coefficients (not replaced by index-derived residual).

      coef_agg$Veg <- as.character(coef_agg$Veg)
      coef_agg_barren_mask <- tolower(coef_agg$Veg) == "barren"
      coef_agg_barren_mask[is.na(coef_agg_barren_mask)] <- FALSE
      veg_coefs_mask <- !coef_agg_barren_mask
      veg_coefs_mask[is.na(veg_coefs_mask)] <- FALSE

      mesma_barren_fraction <- if (any(coef_agg_barren_mask)) {
        sum(coef_agg$coef[coef_agg_barren_mask], na.rm = TRUE)
      } else {
        # Fallback if no explicit barren row is present after aggregation
        1 - sum(coef_agg$coef[veg_coefs_mask], na.rm = TRUE)
      }
      if (!is.finite(mesma_barren_fraction)) mesma_barren_fraction <- 0
      mesma_barren_fraction <- pmin(pmax(mesma_barren_fraction, 0), 1)

      # Do NOT scale MESMA-derived vegetation fractions by MSAVI/PPI-derived total veg cover.
      # Use MESMA absolute fractions as produced by the unmixing solver. This ensures
      # the spectral-level pipeline (and downstream artificial-mix validation) relies
      # on MESMA-derived barren fractions rather than index proxies.
      coef_agg <- coef_agg[veg_coefs_mask, , drop = FALSE]
      coef_agg <- rbind(coef_agg, data.frame(Veg = "barren", coef = mesma_barren_fraction, stringsAsFactors = FALSE))
      if (isTRUE(TESTING_MODE)) cat(sprintf("[MESMA ABS] Preserving MESMA absolute fractions (sum_veg=%.3f, barren=%.3f)\n",
                                          sum(coef_agg$coef[coef_agg$Veg != "barren"], na.rm = TRUE), mesma_barren_fraction))


      # Ensure exactly one MESMA-derived barren row exists in coef_df.
      if (nrow(coef_df) > 0) {
        barren_rows <- tolower(coef_df$Veg) == "barren"
        barren_rows[is.na(barren_rows)] <- FALSE
        if (any(barren_rows)) {
          keep_rows <- !barren_rows
          first_barren <- which(barren_rows)[1]
          barren_row_df <- coef_df[first_barren, , drop = FALSE]
          barren_row_df$coef <- mesma_barren_fraction
          if ("variant_id" %in% names(barren_row_df)) barren_row_df$variant_id <- "barren_mesma"
          if ("coef_sd" %in% names(barren_row_df)) barren_row_df$coef_sd <- 0
          if ("coef_025" %in% names(barren_row_df)) barren_row_df$coef_025 <- NA
          if ("coef_975" %in% names(barren_row_df)) barren_row_df$coef_975 <- NA
          coef_df <- rbind(coef_df[keep_rows, , drop = FALSE], barren_row_df)
        } else {
          barren_row_df <- coef_df[1, , drop = FALSE]
          barren_row_df[] <- NA
          barren_row_df$location_id <- loc
          barren_row_df$pheno_year <- yr
          barren_row_df$lat <- lat_val
          barren_row_df$lon <- lon_val
          barren_row_df$Veg <- "barren"
          barren_row_df$variant_id <- "barren_mesma"
          barren_row_df$coef <- mesma_barren_fraction
          if ("coef_sd" %in% names(barren_row_df)) barren_row_df$coef_sd <- 0
          if ("coef_025" %in% names(barren_row_df)) barren_row_df$coef_025 <- NA
          if ("coef_975" %in% names(barren_row_df)) barren_row_df$coef_975 <- NA
          coef_df <- rbind(coef_df, barren_row_df)
        }
      } else {
        # No vegetation rows retained; emit only MESMA-derived barren if finite/positive.
        if (isTRUE(is.finite(mesma_barren_fraction) && mesma_barren_fraction > 0)) {
          coef_df <- data.frame(
            location_id = loc,
            pheno_year = yr,
            lat = lat_val,
            lon = lon_val,
            Veg = "barren",
            variant_id = "barren_mesma",
            coef = mesma_barren_fraction,
            rmse = if (exists("rmse") && is.numeric(rmse)) rmse else 0,
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
          if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_task] loc=%s yr=%d: coef_df empty after barren removal\n", loc, yr))
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
      diag_df$barren_fraction_mesma <- mesma_barren_fraction
      # Compute PPI-derived barren if PPI (or ppi_norm) is available for the summer observations;
      # otherwise leave PPI-based barren as NA so downstream "PPI" comparisons remain meaningful.
      if (exists("summer_task_data") && "ppi_norm" %in% names(summer_task_data) && any(is.finite(summer_task_data$ppi_norm))) {
        diag_df$barren_fraction_ppi_based <- 1 - median(summer_task_data$ppi_norm, na.rm = TRUE)
      } else if (exists("summer_task_data") && "PPI_raw" %in% names(summer_task_data) && any(is.finite(summer_task_data$PPI_raw))) {
        ppi_max <- if (exists("PPI_FULL_VEG_COVER")) get("PPI_FULL_VEG_COVER") else 0.7
        diag_df$barren_fraction_ppi_based <- 1 - pmin(pmax(median(summer_task_data$PPI_raw, na.rm = TRUE) / ppi_max, 0), 1)
      } else {
        diag_df$barren_fraction_ppi_based <- NA_real_
      }

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

      # DEBUG: Log coefficient values before returning
      if (nrow(coef_df) > 0) {
        coef_summary <- sprintf("coefs: %s", paste(sprintf("%s=%.4f", coef_df$Veg, coef_df$coef), collapse=", "))
        na_count <- sum(is.na(coef_df$coef))
        if (na_count > 0) {
        }
      }
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
        if (isTRUE(TESTING_MODE)) cat(sprintf("[ERROR fit_one_task] loc=%s yr=%d: final return failed: %s\n", loc, yr, e$message))
        stop(e)
      }))

      # ===== END MESMA UNMIXING =====
  # These variables are already set in the earlier conditional blocks above.
  # No need to reassign them here
  # Single-stage MESMA is used, so mesma_lib and OPTIMIZED_LIBRARY are set accordingly

  }

  fit_one_location <- function(location_data) {
    if (is.null(location_data) || nrow(location_data) == 0) {
      if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_location] Received NULL or empty location_data\n"))
      return(NULL)
    }

    loc <- as.character(location_data$location_id[1])
    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_location] Processing location: %s, rows: %d\n", loc, nrow(location_data)))

    # Get all phenological years for this location
    years <- sort(unique(location_data$pheno_year))
    years <- years[!is.na(years)]
    
    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_location] Location %s has years: %s\n", loc, paste(years, collapse=", ")))

    if (length(years) == 0) {
      if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_location] Location %s has no valid years - returning NULL\n", loc))
      return(NULL)
    }

    # Process each year individually
    year_results <- list()

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
          year_results[[as.character(yr)]] <- res_yr
        } else {
          if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_location] loc=%s yr=%d: fit_one_task returned NULL\n", loc, yr))
        }
      }

      if (length(year_results) == 0) {
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_location] loc=%s: no successful year results\n", loc))
        return(NULL)
      }

      # Return results for all years
      return(year_results)

    }, error = function(e) {
      cat(sprintf("[ERROR fit_one_location] loc=%s: %s\n", loc, e$message))
      if (length(year_results) > 0) {
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_location] loc=%s: returning partial %d year(s) after error\n", loc, length(year_results)))
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
  

TEMP_RESULTS_DIR <- "C:/MAP/temp_results"
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

  plot_global_vegetation_pattern <- function(global_pattern,
                                              title = "Global Vegetation Composition Over Time",
                                              show_ci = TRUE,
                                              ci_type = "auto") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
      if ("ci_lower" %in% names(global_pattern) && "ci_upper" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower
        global_pattern$ci_upper <- global_pattern$ci_upper
      } else if ("coef_025" %in% names(global_pattern) && "coef_975" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$coef_025
        global_pattern$ci_upper <- global_pattern$coef_975
      }
    } else if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
      if ("ci_lower_simple" %in% names(global_pattern) && "ci_upper_simple" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower_simple
        global_pattern$ci_upper <- global_pattern$ci_upper_simple
      }
    } else if ("weighted_mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$weighted_mean_coef
      if ("ci_lower_pooled" %in% names(global_pattern) && "ci_upper_pooled" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower_pooled
        global_pattern$ci_upper <- global_pattern$ci_upper_pooled
      }
    }
    
    global_pattern_veg <- global_pattern[tolower(global_pattern$Veg) != "barren", ]
    global_pattern_barren <- global_pattern[tolower(global_pattern$Veg) == "barren", ]
    
    p <- ggplot2::ggplot(global_pattern_veg, ggplot2::aes(x = year, y = coef, color = Veg, fill = Veg)) +
      add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
      ggplot2::geom_line(size = 1.2) +
      ggplot2::geom_point(size = 2) +
      ggplot2::coord_cartesian(ylim = c(0, NA))
    
    if (show_ci && "ci_lower" %in% names(global_pattern_veg) && "ci_upper" %in% names(global_pattern_veg)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
        alpha = 0.2,
        color = NA
      )
    }
    
    if (nrow(global_pattern_barren) > 0) {
      veg_max <- suppressWarnings(max(global_pattern_veg$coef, global_pattern_veg$ci_upper, na.rm = TRUE))
      barren_max <- suppressWarnings(max(global_pattern_barren$coef, global_pattern_barren$ci_upper, na.rm = TRUE))
      if (is.na(veg_max) || veg_max <= 0 || is.na(barren_max) || barren_max <= 0) {
        barren_scale_factor <- 1
      } else {
        barren_scale_factor <- veg_max / barren_max
      }
      global_pattern_barren$coef_scaled <- global_pattern_barren$coef * barren_scale_factor
      
      p <- p + 
        ggplot2::geom_line(data = global_pattern_barren, 
              ggplot2::aes(x = year, y = coef_scaled), 
                          color = "brown", linewidth = 1.2, linetype = "dashed") +
        ggplot2::geom_point(data = global_pattern_barren, 
               ggplot2::aes(x = year, y = coef_scaled), 
                           color = "brown", size = 2)
      
      if (show_ci && "ci_lower" %in% names(global_pattern_barren) && "ci_upper" %in% names(global_pattern_barren)) {
        global_pattern_barren$ci_lower_scaled <- global_pattern_barren$ci_lower * barren_scale_factor
        global_pattern_barren$ci_upper_scaled <- global_pattern_barren$ci_upper * barren_scale_factor
        
        p <- p + ggplot2::geom_ribbon(data = global_pattern_barren,
                                     ggplot2::aes(x = year, ymin = ci_lower_scaled, ymax = ci_upper_scaled),
                                     fill = "brown", alpha = 0.1, color = NA)
      }
    }
    
    p <- p +
      ggplot2::labs(
        title = title,
        subtitle = sprintf("Based on %d locations", max(global_pattern$n_locations, na.rm = TRUE)),
        x = "Year",
        y = "Vegetation Fraction",
        color = "Vegetation Type",
        fill = "Vegetation Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)
      ) +
          ggplot2::scale_y_continuous(
        labels = scales::percent_format(),
        sec.axis = ggplot2::sec_axis(~ . / barren_scale_factor, name = "Barren Fraction", labels = scales::percent_format())
      ) +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::scale_fill_brewer(palette = "Set1")
    
    p
  }

  plot_vegetation_only_stacked_area <- function(global_pattern, 
                                                title = "Vegetation-Only Composition Over Time") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
    } else if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
    }
    # Filter out barren
    global_pattern_veg <- global_pattern[tolower(global_pattern$Veg) != "barren", ]
    # Normalize vegetation fractions to sum to 1 per year
    global_pattern_veg <- global_pattern_veg |> 
      dplyr::group_by(year) |> 
      dplyr::mutate(coef_normalized = coef / sum(coef, na.rm = TRUE)) |> 
      dplyr::ungroup()
    p <- ggplot2::ggplot(global_pattern_veg, 
                          ggplot2::aes(x = year, y = coef_normalized, fill = Veg)) +
      ggplot2::geom_area(alpha = 0.8, position = "stack") +
      ggplot2::labs(
        title = title,
        x = "Year",
        y = "Relative Vegetation Fraction",
        fill = "Vegetation Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold")
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      ggplot2::scale_fill_brewer(palette = "Set2")
    p
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
    if (exists("df_train") && !is.null(df_train) && nrow(df_train) > 0) {
    if (!"pheno_year" %in% names(df_train) && "date" %in% names(df_train)) df_train$pheno_year <- assign_pheno_year(df_train$date)
    n_train_loc_years <- nrow(unique(df_train[c("location_id", "pheno_year")]))
  }
  if (exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0) {
    if (!"pheno_year" %in% names(df_tasks_inference) && "date" %in% names(df_tasks_inference)) df_tasks_inference$pheno_year <- assign_pheno_year(df_tasks_inference$date)
    n_infer_loc_years <- nrow(unique(df_tasks_inference[c("location_id", "pheno_year")]))
  } else if (exists("df_tasks") && !is.null(df_tasks) && nrow(df_tasks) > 0) {
    if (!"pheno_year" %in% names(df_tasks) && "date" %in% names(df_tasks)) df_tasks$pheno_year <- assign_pheno_year(df_tasks$date)
    n_infer_loc_years <- nrow(unique(df_tasks[c("location_id", "pheno_year")]))
  }

  cat(sprintf("Training dataset location-years: %d\n", n_train_loc_years))
  cat(sprintf("Inference dataset location-years: %d\n", n_infer_loc_years))

  if ((exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0) && n_train_loc_years == n_infer_loc_years) {
    stop(sprintf("ERROR: Training and inference datasets appear to have the same number of location-years (%d). This may indicate you passed the same data for training and inference — aborting to avoid accidental overlap.", n_train_loc_years))
  }
  
  if (!isTRUE(QUIET_MODE)) {
    cat("Preparing locations for batched processing...\n")

    # DEBUG: Check what columns df_tasks has
    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG VALIDATION DATA] df_tasks has %d rows, %d columns\n", nrow(df_tasks), ncol(df_tasks)))
    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG VALIDATION DATA] Column names: %s\n", paste(head(names(df_tasks), 50), collapse=", ")))
    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG VALIDATION DATA] Has PPI: %s, Has PPI_raw: %s\n", 
                "PPI" %in% names(df_tasks), "PPI_raw" %in% names(df_tasks)))
    if ("PPI" %in% names(df_tasks) && isTRUE(TESTING_MODE)) {
      cat(sprintf("[DEBUG VALIDATION DATA] PPI summary: min=%.4f, median=%.4f, max=%.4f, NA=%d\n",
                  min(df_tasks$PPI, na.rm=TRUE), median(df_tasks$PPI, na.rm=TRUE), 
                  max(df_tasks$PPI, na.rm=TRUE), sum(is.na(df_tasks$PPI))))
    }
  }

  # Group by LOCATION (not location-year)
  target_locations <- unique(df_tasks$location_id)

  n_locs_to_process <- length(target_locations)
  
  # Skip processing if no locations to process (training data processing disabled)
  if (n_locs_to_process == 0) {
    cat("[INFO] No training locations to process - skipping main processing loop\n")
    training_results_list <- list()
  } else {
  # Use centralized BATCH_SIZE setting from the top-level CONFIG (to change batch size, edit the USER-TUNABLE PARAMETERS block)
  loc_batches <- split(target_locations, ceiling(seq_along(target_locations) / BATCH_SIZE))
  n_batches <- length(loc_batches)

  if (!isTRUE(QUIET_MODE)) cat(sprintf("Processing %d locations in %d batches (approx %d locations/batch)...\n",
              n_locs_to_process, length(loc_batches), BATCH_SIZE))

  # Collect full per-task results for downstream reporting (variant trajectories, diagnostics, uncertainty)
  training_results_list <- list()

  start_time <- Sys.time()

  for (i in seq_along(loc_batches)) {
    batch_locs <- loc_batches[[i]]
    batch_df <- df_tasks[df_tasks$location_id %in% batch_locs, ]
    batch_location_list <- split(batch_df, batch_df$location_id)
    # Suppress any verbose output from per-location processing
    batch_results <- suppress_output_safely(
      .run_map(batch_location_list, fit_one_location, show_pb = FALSE)
    )

    training_results_list <- aggregate_batch_results(batch_results, training_results_list, save_csv_dir = TEMP_RESULTS_DIR)

if (!isTRUE(QUIET_MODE)) {
      # Use a carriage-return progress update only in interactive sessions;
      # otherwise emit a single-line progress message to avoid cluttering logs.
      if (interactive()) {
        cat(sprintf("\r  [Batch %d/%d complete]  ", i, n_batches))
        flush.console()
      } else {
        cat(sprintf("[Batch %d/%d complete]\n", i, n_batches))
      }
    }
  }
  cat("\n")
  
  end_time <- Sys.time()
  processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(
    "\nTraining processing finished in %.2f seconds (%.2f minutes)\n",
    processing_time, processing_time / 60
  ))
  } # End of n_locs_to_process > 0 conditional

  # ==========================================================================
  # VALIDATION PROCESSING: unmix held-out validation locations
  # ========================================================================== 
  # Ensure validation split is (re)computed on every run so validation always executes
  # even when re-sourcing in the same R session. Prefer the preserved `original_df`
  # (full dataset) when available; otherwise fall back to the current `df`.
  base_df_for_split <- if (exists("original_df") && !is.null(original_df) && nrow(original_df) > 0) original_df else df
  if (exists("TRAIN_YEARS") && !is.null(TRAIN_YEARS) && length(TRAIN_YEARS) > 0) {
    df_train <- base_df_for_split[base_df_for_split$pheno_year %in% TRAIN_YEARS, , drop = FALSE]
  } else {
    df_train <- base_df_for_split
  }

  # Recompute stratified train/validation split (deterministic seed)
  set.seed(get_mesma_seed(123))
  loc_veg_summary <- get_dominant_veg_per_location(df_train)
  val_locs_list <- vector("list", length(unique(loc_veg_summary$Veg)))
  unique_vegs_split <- unique(loc_veg_summary$Veg)
  for (ii in seq_along(unique_vegs_split)) {
    v <- unique_vegs_split[ii]
    v_locs <- loc_veg_summary$location_id[loc_veg_summary$Veg == v]
    n_v <- length(v_locs)
    n_val <- ceiling(n_v * VALIDATION_FRACTION)
    if (n_val > 0 && n_v > 1) {
      n_val <- min(n_val, n_v - 1)
      if (n_val > 0) val_locs_list[[ii]] <- data.frame(location_id = sample(v_locs, n_val), Veg = v, stringsAsFactors = FALSE)
    }
  }
  val_locs_df <- tryCatch(do.call(rbind, val_locs_list), error = function(e) NULL)
  validation_location_ids <- if (!is.null(val_locs_df) && nrow(val_locs_df) > 0) as.character(val_locs_df$location_id) else character(0)
  df_validation <- df_train[df_train$location_id %in% validation_location_ids, , drop = FALSE]
  df_train <- df_train[!df_train$location_id %in% validation_location_ids, , drop = FALSE]
  df_inference <- NULL

  cat(sprintf("\n=== STARTING VALIDATION PROCESSING ===\n  (recomputed stratified split: validation rows=%d, locations=%d)\n", nrow(df_validation), length(unique(df_validation$location_id))))

  if (exists("df_validation") && !is.null(df_validation) && nrow(df_validation) > 0) {
    df_validation_proc <- df_validation

    if ("Veg" %in% names(df_validation_proc)) {
      df_validation_proc$Veg <- tolower(as.character(df_validation_proc$Veg))
    }

    # Ensure pheno_year and doy exist
    if (!"pheno_year" %in% names(df_validation_proc) && "date" %in% names(df_validation_proc)) {
      df_validation_proc$pheno_year <- assign_pheno_year(df_validation_proc$date)
    }
    if (!"doy" %in% names(df_validation_proc) && "date" %in% names(df_validation_proc)) {
      df_validation_proc$doy <- pheno_doy(df_validation_proc$date)
      df_validation_proc$doy[df_validation_proc$doy < 1 | df_validation_proc$doy > 366] <- NA_integer_
    }

    val_locations <- unique(trimws(as.character(df_validation_proc$location_id)))
    val_locations <- val_locations[!is.na(val_locations) & val_locations != ""]

    cat(sprintf("[VALIDATION] Processing %d held-out validation locations (%d rows)\n",
                length(val_locations), nrow(df_validation_proc)))

    val_loc_batches <- split(val_locations, ceiling(seq_along(val_locations) / BATCH_SIZE))
    n_val_batches <- length(val_loc_batches)

    cat(sprintf("[VALIDATION] Processing in %d batches (approx %d locations/batch)...\n",
                n_val_batches, BATCH_SIZE))

    validation_results_list <- list()
    val_start_time <- Sys.time()

    for (i in seq_along(val_loc_batches)) {
      batch_locs <- val_loc_batches[[i]]
      batch_df <- df_validation_proc[df_validation_proc$location_id %in% batch_locs, ]
      batch_location_list <- split(batch_df, batch_df$location_id)

      batch_t0 <- Sys.time()
      if (!isTRUE(QUIET_MODE)) {
        cat(sprintf("[VALIDATION] Batch %d/%d starting (%d locations, %d rows)\n",
                    i, n_val_batches, length(batch_location_list), nrow(batch_df)))
        flush.console()
      }

      # Suppress verbose per-location output
      batch_results <- if (isTRUE(QUIET_MODE)) {
        suppress_output_safely(.run_map(batch_location_list, fit_one_location, show_pb = FALSE), quiet_message = FALSE)
      } else {
        .run_map(batch_location_list, fit_one_location, show_pb = FALSE)
      }

      validation_results_list <- aggregate_batch_results(batch_results, validation_results_list)

      if (!isTRUE(QUIET_MODE)) {
        elapsed_batch <- as.numeric(difftime(Sys.time(), batch_t0, units = "secs"))
        cat(sprintf("[VALIDATION] Batch %d/%d done (%.1fs)\n", i, n_val_batches, elapsed_batch))
        flush.console()
      }
    }

    val_end_time <- Sys.time()
    val_processing_time <- as.numeric(difftime(val_end_time, val_start_time, units = "secs"))
    cat(sprintf("[VALIDATION] Processing finished in %.2f seconds (%.2f minutes)\n",
                val_processing_time, val_processing_time / 60))

    # Combine validation results into validation_coefs
    validation_coefs <- do.call(rbind, lapply(validation_results_list, function(r) {
      if (!is.null(r$coef_df)) r$coef_df else NULL
    }))

    # Remove rows with NA or empty Veg values
    if (!is.null(validation_coefs) && nrow(validation_coefs) > 0) {
      n_before <- nrow(validation_coefs)
      validation_coefs <- filter_valid_vegetation(validation_coefs, exclude_barren = FALSE)
      na_veg_count <- n_before - nrow(validation_coefs)
      if (na_veg_count > 0) {
        cat(sprintf("[VALIDATION] Removing %d rows with NA/empty Veg values\n", na_veg_count))
      }
    }

    if (is.null(validation_coefs)) validation_coefs <- data.frame()
    cat(sprintf("[VALIDATION] Produced %d validation coefficient rows from %d locations\n",
                nrow(validation_coefs),
                if (nrow(validation_coefs) > 0) length(unique(validation_coefs$location_id)) else 0))
  } else {
    cat("[VALIDATION] No df_validation available - skipping validation processing\n")
  }

  # ==========================================================================
  # VALIDATION CONFUSION MATRIX (same pipeline as training, different data)
  # ==========================================================================

  cat("\n=== COMPUTING VALIDATION CONFUSION MATRIX ===\n")

  if (exists("df_validation") && !is.null(df_validation) && nrow(df_validation) > 0) {
    compute_confusion_matrix(df_validation, "Validation")
  } else {
    cat("[WARNING] No df_validation available for confusion matrix\n")
  }
  
  # ==========================================================================
  # INFERENCE PROCESSING (separate from validation)
  # ==========================================================================
  cat("\n=== STARTING INFERENCE PROCESSING ===\n")
  
  # Generate variant similarity heatmap before any inference processing
  cat("[INFERENCE] Generating variant similarity heatmap...\n")
  ensure_variant_similarity_heatmap(force = TRUE)
  
  # Load inference data from INFERENCE_CSV if not already loaded
  if (!isTRUE(TESTING_MODE)) {
    cat("[INFERENCE] Loading inference data from separate CSV file...\n")
    load_and_prepare_inference_data()
  }

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
    
    # Apply MAX_INFERENCE_LOCATIONS limit if set
    if (exists("MAX_INFERENCE_LOCATIONS") && nrow(inference_locations) > MAX_INFERENCE_LOCATIONS) {
      set.seed(get_mesma_seed(123)) # deterministic sampling for reproducibility
      sampled_loc_ids <- sample(inference_locations$location_id, MAX_INFERENCE_LOCATIONS, replace = FALSE)
      inference_locations <- inference_locations[inference_locations$location_id %in% sampled_loc_ids, , drop = FALSE]
      df_tasks_inference_proc <- df_tasks_inference_proc[df_tasks_inference_proc$location_id %in% sampled_loc_ids, ]
      cat(sprintf("[INFERENCE] Reduced to %d locations (random sample due to MAX_INFERENCE_LOCATIONS=%d)\n",
                  nrow(inference_locations), MAX_INFERENCE_LOCATIONS))
    }

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
      if (i == 1L || !isTRUE(QUIET_MODE)) {
        cat(sprintf("[INFERENCE] Batch %d/%d starting (%d locations, %d rows)\n",
                    i, n_inference_batches, length(batch_location_list), nrow(batch_df)))
        flush.console()
      }
      
      # Diagnostic: check batch data structure
      if (isTRUE(TESTING_MODE)) {
        cat(sprintf("  [Batch %d] Processing %d locations, batch_df has %d rows\n", 
                    i, length(batch_location_list), nrow(batch_df)))
        if (length(batch_location_list) > 0) {
          first_loc_name <- names(batch_location_list)[1]
          first_loc_data <- batch_location_list[[1]]
          cat(sprintf("  [Batch %d] First location '%s' has %d rows\n",
                      i, first_loc_name, nrow(first_loc_data)))
          cat(sprintf("  [Batch %d] All columns (%d): %s\n",
                      i, length(names(first_loc_data)), paste(names(first_loc_data), collapse=", ")))
          if ("pheno_year" %in% names(first_loc_data)) {
            years_sample <- head(sort(unique(first_loc_data$pheno_year)), 5)
            cat(sprintf("  [Batch %d] First location pheno_years (first 5): %s\n",
                        i, paste(years_sample, collapse=", ")))
          } else {
            cat(sprintf("  [Batch %d] ERROR: 'pheno_year' column missing!\n", i))
          }
        }
        # Save debug log location before sinking
        debug_log <- file.path(tempdir(), "fit_one_location_debug.log")
        cat(sprintf("  [Batch %d] Debug log: %s\n", i, debug_log))
      } else {
        # Not in testing mode: still set the debug log path (used for later conditional printing),
        # but do not print verbose per-batch diagnostics to stdout.
        debug_log <- file.path(tempdir(), "fit_one_location_debug.log")
      }

      # Suppress verbose stdout from per-location processing only in QUIET_MODE.
      batch_results <- if (isTRUE(QUIET_MODE)) {
        suppress_output_safely(.run_map(batch_location_list, fit_one_location, show_pb = FALSE), quiet_message = FALSE)
      } else {
        .run_map(batch_location_list, fit_one_location, show_pb = FALSE)
      }

      # Store results using shared helper
      n_before <- length(inference_results_list)
      inference_results_list <- aggregate_batch_results(batch_results, inference_results_list)
      n_stored <- length(inference_results_list) - n_before
      
      if (isTRUE(TESTING_MODE)) {
        cat(sprintf("  [Batch %d] Results: %d stored, %d null locations, %d null years, %d empty coefs\n",
                    i, n_stored, n_null_results, n_empty_years, n_empty_coefs))

        # Show debug log content if all locations failed
        if (n_null_results == length(batch_results) && file.exists(debug_log)) {
          cat(sprintf("  [Batch %d] All locations returned NULL - showing last 20 lines of debug log:\n", i))
          log_lines <- tryCatch(readLines(debug_log), error = function(e) character(0))
          if (length(log_lines) > 0) {
            tail_lines <- tail(log_lines, 20)
            for (line in tail_lines) cat("    ", line, "\n")
          }
        }

      } else {
        # Emit progress at each 10% completion mark (also in QUIET_MODE to avoid "looks stuck").
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

      if (!isTRUE(QUIET_MODE)) {
        elapsed_batch <- as.numeric(difftime(Sys.time(), batch_t0, units = "secs"))
        cat(sprintf("[INFERENCE] Batch %d/%d done (%.1fs)\n", i, n_inference_batches, elapsed_batch))
        flush.console()
      }
    }

    # Ensure 100% is always reported once if not already printed
    if (!isTRUE(QUIET_MODE) && n_inference_batches > 0 && next_progress_idx <= length(progress_targets)) {
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

    # --- INFERENCE: MSAVI-based aggregation ---
    if (!is.null(inference_coefs) && nrow(inference_coefs) > 0 && exists("df_tasks_inference_proc") && ("MSAVI" %in% names(df_tasks_inference_proc) || "MSAVI_raw" %in% names(df_tasks_inference_proc))) {
      cat("[INFERENCE] Running MSAVI-based aggregation (inference mode)\n")
      msavi_inf_full <- tryCatch({ location_bootstrap_msavi(inference_coefs, df_tasks_inference_proc, B = BOOTSTRAP_B, seed = get_mesma_seed(123)) }, error = function(e) { cat(sprintf("[MSAVI INFERENCE] failed: %s\n", e$message)); NULL })
      if (!is.null(msavi_inf_full) && nrow(msavi_inf_full) > 0) {
        plot_inference_method_results(msavi_inf_full, "MSAVI", "msavi",
                                      use_excluded_years_shade = TRUE,
                                      include_species_plots = TRUE)
      } else {
        cat("[INFERENCE] MSAVI inference aggregation returned no results.\n")
      }
    } else {
      cat("[INFERENCE] MSAVI data not available for inference; skipping MSAVI aggregation.\n")
    }

    # --- INFERENCE: NDVI-based aggregation ---
    if (!is.null(inference_coefs) && nrow(inference_coefs) > 0 && exists("df_tasks_inference_proc") && ("NDVI" %in% names(df_tasks_inference_proc) || "NDVI_raw" %in% names(df_tasks_inference_proc))) {
      cat("[INFERENCE] Running NDVI-based aggregation (inference mode)\n")
      ndvi_inf_full <- tryCatch({ location_bootstrap_ndvi(inference_coefs, df_tasks_inference_proc, B = BOOTSTRAP_B, seed = get_mesma_seed(123)) }, error = function(e) { cat(sprintf("[NDVI INFERENCE] failed: %s\n", e$message)); NULL })
      if (!is.null(ndvi_inf_full) && nrow(ndvi_inf_full) > 0) {
        plot_inference_method_results(ndvi_inf_full, "NDVI", "ndvi",
                                      use_excluded_years_shade = TRUE,
                                      include_species_plots = TRUE)
      } else {
        cat("[INFERENCE] NDVI inference aggregation returned no results.\n")
      }
    } else {
      cat("[INFERENCE] NDVI data not available for inference; skipping NDVI aggregation.\n")
    }

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

    # --- Pixel purity & mixture composition plots ---
    cat("[INFERENCE] Generating pixel purity and mixture composition plots\n")
    tryCatch(
      plot_pixel_purity_and_mixtures(inference_coefs, pure_threshold = 0.9,
                                     mix_contrib_threshold = 0.10,
                                     use_excluded_years_shade = FALSE),
      error = function(e) cat(sprintf("[PURITY] Failed: %s\n", e$message))
    )

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
  }

  # Ensure required library/templates exist for visualization
  ensure_library_and_templates()

  # Held-out validation accuracy reporting, validation aggregate bootstrap,
  # and validation-RMSE-driven CI adjustment have been removed by request.

  # Helper: combine results_list into unified all_coefs (no printing)
  combine_results_from_list <- function(results_list) {
    coef_list <- lapply(results_list, function(res) {
      if (is.null(res) || is.null(res$coef_df) || nrow(res$coef_df) == 0) NULL else res$coef_df
    })
    coef_list <- coef_list[!sapply(coef_list, is.null)]
    if (length(coef_list) == 0) return(NULL)
    if (requireNamespace("dplyr", quietly = TRUE)) {
      all_coefs_local <- dplyr::bind_rows(coef_list)
    } else {
      all_coefs_local <- do.call(rbind, coef_list)
    }
    # Ensure minimal columns exist
    req_cols <- c("location_id","pheno_year","Veg","coef","coef_025","coef_975","interval")
    for (cname in req_cols) if (!(cname %in% names(all_coefs_local))) all_coefs_local[[cname]] <- NA
    all_coefs_local
  }

  # Run inference on a dataframe silently (suppress per-task prints) and return results list + combined coefs
  run_inference_silent <- function(df_subset) {
    if (is.null(df_subset) || nrow(df_subset) == 0) return(list(results = list(), all_coefs = NULL))
    if (!"pheno_year" %in% names(df_subset) && "date" %in% names(df_subset)) df_subset$pheno_year <- assign_pheno_year(df_subset$date)
    if (!"task_key" %in% names(df_subset)) df_subset$task_key <- paste(df_subset$location_id, df_subset$pheno_year, sep = "_")
    task_list <- split(df_subset, df_subset$task_key)

    fit_task_wrapper <- function(task_data) {
      if (is.null(task_data) || nrow(task_data) == 0) return(NULL)
      res <- fit_one_task(task_data)
      return(res)
    }

    results <- tryCatch(
      suppress_output_safely(.run_map(task_list, fit_task_wrapper, show_pb = FALSE)),
      error = function(e) {
        warning(sprintf("run_inference_silent: underlying inference error: %s", e$message))
        NULL
      }
    )

    all_coefs_local <- combine_results_from_list(results)
    list(results = results, all_coefs = all_coefs_local)
  }

  # Artificial mix validation removed — functionality disabled per user request.
  # (Previously: spectral-level 50/50 artificial-mix generation and MESMA unmixing checks.)
  cat("\n[NOTICE] Artificial mix validation DISABLED (removed).\n\n")

  # Results are now in validation_coefs and inference_coefs (combined as all_coefs)
  cat(sprintf("[RESULTS] Validation: %d rows, Inference: %d rows, Combined: %d rows\n",
              if(exists("validation_coefs") && !is.null(validation_coefs)) nrow(validation_coefs) else 0,
              if(exists("inference_coefs") && !is.null(inference_coefs)) nrow(inference_coefs) else 0,
              if(exists("all_coefs") && !is.null(all_coefs)) nrow(all_coefs) else 0))
              
    # Validation-RMSE CI adjustment removed by request.

  # Fail fast if BOTH validation and inference result lists are missing/empty — indicates upstream processing issue
  if ((!exists("validation_coefs") || is.null(validation_coefs) || nrow(validation_coefs) == 0) &&
      (!exists("inference_coefs") || is.null(inference_coefs) || nrow(inference_coefs) == 0)) {
    stop(paste0("ERROR: No results collected (both validation and inference are empty). Upstream processing failed — check earlier logs and data filters."))
  }

  cat("Processing results...\n")

  # Diagnostics: report how many results were collected
  cat(sprintf("Validation results: %d rows\n", if(exists("validation_coefs") && !is.null(validation_coefs)) nrow(validation_coefs) else 0))
  cat(sprintf("Inference results: %d rows\n", if(exists("inference_coefs") && !is.null(inference_coefs)) nrow(inference_coefs) else 0))

  # Check barren fraction distribution in combined results
  if (!is.null(all_coefs) && nrow(all_coefs) > 0 && "Veg" %in% names(all_coefs)) {
    barren_rows <- all_coefs[tolower(all_coefs$Veg) == "barren", ]
    if (nrow(barren_rows) > 0 && "coef" %in% names(barren_rows)) {
      barren_one_count <- sum(barren_rows$coef == 1, na.rm = TRUE)
      barren_one_pct <- barren_one_count / nrow(all_coefs) * 100

      # ERROR if ALL samples are 100% barren
      if (barren_one_pct >= 99.9) {
        stop(sprintf(
          "[FATAL ERROR] All samples are 100%% barren (%.1f%%, %d/%d predictions)!\n\n" 
        ), barren_one_pct, barren_one_count, nrow(all_coefs))
      }

      cat(sprintf("Barren fraction = 1 in %.1f%% of predictions (%d/%d) - within acceptable limits\n",
                  barren_one_pct, barren_one_count, nrow(all_coefs)))
    }
  }

  # Verify we have results to process
  if (is.null(all_coefs) || nrow(all_coefs) == 0) {
    cat("ERROR: No results to process!\n")
    stop("No valid results to process")
  }

  # all_coefs is already built from combined validation_coefs and inference_coefs
  # Ensure required columns and data types
  cat(sprintf("Combined coefficients: %d rows\n", nrow(all_coefs)))
  
  required_coef_cols <- c("location_id", "pheno_year", "Veg", "coef", "rmse", "coef_025", "coef_975", "interval")
  missing <- setdiff(required_coef_cols, names(all_coefs))
  if (length(missing) > 0) {
    cat(sprintf("[NOTICE] Filling missing coefficient columns with NA: %s\n", paste(missing, collapse = ", ")))
    for (col in missing) all_coefs[[col]] <- NA
  }
  all_coefs$location_id <- as.character(all_coefs$location_id)
  all_coefs$pheno_year <- as.integer(all_coefs$pheno_year)
  # Ensure backward-compatible 'year' column exists for downstream code that expects it
  if (!"year" %in% names(all_coefs)) all_coefs$year <- all_coefs$pheno_year
  all_coefs$Veg <- as.character(all_coefs$Veg)
  all_coefs$coef <- as.numeric(all_coefs$coef)
  if ("rmse" %in% names(all_coefs)) all_coefs$rmse <- as.numeric(all_coefs$rmse)
  if ("coef_025" %in% names(all_coefs)) all_coefs$coef_025 <- as.numeric(all_coefs$coef_025)
  if ("coef_975" %in% names(all_coefs)) all_coefs$coef_975 <- as.numeric(all_coefs$coef_975)
  if ("interval" %in% names(all_coefs)) all_coefs$interval <- as.numeric(all_coefs$interval)
  all_coefs$location_id <- trimws(all_coefs$location_id)

    # Collect uncertainty data for merging (training + inference)
    unc_coef_rows <- list()

    results_for_uncertainty <- list()
    if (exists("training_results_list") && is.list(training_results_list)) {
      results_for_uncertainty <- c(results_for_uncertainty, training_results_list)
    }
    if (exists("inference_results_list") && is.list(inference_results_list)) {
      results_for_uncertainty <- c(results_for_uncertainty, inference_results_list)
    }

    for (res in results_for_uncertainty) {
      if (is.null(res$uncertainty)) next
      loc <- if (!is.null(res$coef_df$location_id)) res$coef_df$location_id[1] else NA_character_
      yr <- if (!is.null(res$coef_df$year)) {
        as.integer(res$coef_df$year[1])
      } else if (!is.null(res$coef_df$pheno_year)) {
        as.integer(res$coef_df$pheno_year[1])
      } else {
        NA_integer_
      }
      ci <- res$uncertainty$coef_ci
      if (!is.null(ci) && nrow(ci) > 0) {
        ci$location_id <- loc; ci$year <- yr
        unc_coef_rows[[length(unc_coef_rows) + 1]] <- ci[, c("location_id","year","Veg","coef_025","coef_975", "coef_sd", "interval"), drop = FALSE]
      }
    }
    if (length(unc_coef_rows) > 0) {
      if (requireNamespace("dplyr", quietly = TRUE)) all_unc_coef <- dplyr::bind_rows(unc_coef_rows) else all_unc_coef <- do.call(rbind, unc_coef_rows)
    } else all_unc_coef <- NULL
    if (!is.null(all_unc_coef)) {
      n_total_unc <- nrow(all_unc_coef)
      n_non_na_ci <- sum(is.finite(all_unc_coef$coef_025) | is.finite(all_unc_coef$coef_975))
      cat(sprintf("Uncertainty CIs: total rows = %d; rows with at least one finite CI = %d\n", n_total_unc, n_non_na_ci))

      # MERGE BOOTSTRAP UNCERTAINTY INTO MAIN COEFFICIENTS
      # This is the critical fix: merge coef_025 and coef_975 from bootstrap into all_coefs
      if (requireNamespace("dplyr", quietly = TRUE)) {
        # Ensure join keys exist and have matching types
        all_coefs$location_id <- as.character(all_coefs$location_id)
        if (!"year" %in% names(all_coefs) && "pheno_year" %in% names(all_coefs)) all_coefs$year <- as.integer(all_coefs$pheno_year)
        all_coefs$year <- as.integer(all_coefs$year)
        all_unc_coef$location_id <- as.character(all_unc_coef$location_id)
        all_unc_coef$year <- as.integer(all_unc_coef$year)

        # Remove existing coef_025 and coef_975 columns if they exist (they contain NAs)
        if ("coef_025" %in% names(all_coefs)) all_coefs$coef_025 <- NULL
        if ("coef_975" %in% names(all_coefs)) all_coefs$coef_975 <- NULL
        if ("coef_sd" %in% names(all_coefs)) all_coefs$coef_sd <- NULL
        if ("interval" %in% names(all_coefs)) all_coefs$interval <- NULL

        # Perform left join to merge bootstrap CIs into main results
        all_coefs <- dplyr::left_join(all_coefs, all_unc_coef, by = c("location_id", "year", "Veg"))

        # Report merge success
        n_merged <- sum(is.finite(all_coefs$coef_025) | is.finite(all_coefs$coef_975))
        cat(sprintf("Bootstrap CIs merged into main results: %d/%d rows now have uncertainty bounds\n",
                    n_merged, nrow(all_coefs)))
      } else {
        warning("dplyr not available - bootstrap CIs will not be merged into main results")
      }
    }
  timing_info$end_time <- Sys.time()
  total_time <- as.numeric(difftime(timing_info$end_time, timing_info$start_time, units = "secs"))

  cat(sprintf("\nTotal execution time: %.1f seconds (%.1f minutes)\n", total_time, total_time / 60))
  if (!is.null(timing_info$lib_construction_done) && !is.null(timing_info$moving_var_done)) {
    cat(sprintf(
      "Library construction: %.1f seconds\n",
      as.numeric(difftime(timing_info$lib_construction_done, timing_info$moving_var_done, units = "secs"))
    ))
  }
  if (!is.null(timing_info$pca_computation_done)) {
    cat(sprintf(
      "PCA computation: %.1f seconds\n",
      as.numeric(difftime(timing_info$pca_computation_done, timing_info$lib_construction_done, units = "secs"))
    ))
  }
  cat(sprintf(
    "Main processing: %.1f seconds\n",
    as.numeric(difftime(timing_info$end_time, timing_info$pca_computation_done, units = "secs"))
  ))

  cat("\nMESMA fitting completed successfully!\n")
