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

# Compute a per-location DVI soil baseline for PPI.
# Baseline is defined as the median of the lowest `quantile_p` fraction of DVI
# observations within each location. This is deterministic, per-location, and
# avoids any constant/default soil baseline.
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
INFERENCE_CSV <- "C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_kon (1).csv"

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
HERBS_WOODY_COLORS <- c("Herbs" = "#2E8B57", "Woody" = "#8B4513")
SPECIES_COLORS <- c("Herbs" = "#2E8B57", "Populus" = "#228B22", "Tamarix" = "#CD853F", "Woody Unknown" = "#808080")
WOODY_TYPES_COLORS <- c("Tamarix" = "#CD853F", "Populus" = "#228B22", "Woody Unknown" = "#808080")

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
                                          include_woody_types_plot = FALSE) {
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

  # 3. Herbs vs Woody
  herbs_woody <- full_data[veg_norm %in% c("herbs", "woody"), ]
  if (nrow(herbs_woody) == 0) {
    cat(sprintf("No herbs/woody data available for inference %s herbs vs woody plot.\n", method))
    return(invisible(NULL))
  }
  herbs_woody$Veg <- ifelse(tolower(herbs_woody$Veg) == "herbs", "Herbs", "Woody")
  p_hw <- ggplot(herbs_woody, aes(x = year, y = global_coef, color = Veg, group = Veg)) +
    .ts_layers() +
    geom_line(linewidth = 1) +
    geom_point(show.legend = FALSE) +
    geom_ribbon(aes(ymin = coef_025, ymax = coef_975, fill = Veg), alpha = 0.15, color = NA) +
    scale_color_manual(values = HERBS_WOODY_COLORS) +
    scale_fill_manual(values = HERBS_WOODY_COLORS) +
    labs(title = paste0("Inference ", method, ": Herbs vs Woody Vegetation"),
         x = "Year", y = "Total Normalized Fraction", color = "Type", fill = "Type") +
    theme_minimal()
  ggsave(file.path(OUT_DIR, paste0("inference_", file_prefix, "_herbs_vs_woody.png")), p_hw, width = 8, height = 6)
  readr::write_csv(herbs_woody, file.path(OUT_DIR, paste0("inference_", file_prefix, "_herbs_vs_woody.csv")))
  cat(sprintf("Saved inference %s herbs vs woody plot to: %s\n", method, file.path(OUT_DIR, paste0("inference_", file_prefix, "_herbs_vs_woody.png"))))

  # 4. Species-level plots (no "Woody" aggregation)
  if (include_species_plots) {
    species_data <- full_data[veg_norm %in% c("herbs","populus","tamarix","woody_unknown"), ]
    if (nrow(species_data) > 0) {
      species_data$Veg <- dplyr::case_when(
        tolower(species_data$Veg) == "herbs" ~ "Herbs",
        tolower(species_data$Veg) == "populus" ~ "Populus",
        tolower(species_data$Veg) == "tamarix" ~ "Tamarix",
        tolower(species_data$Veg) == "woody_unknown" ~ "Woody Unknown",
        TRUE ~ species_data$Veg
      )

      p_sp_ts <- ggplot(species_data, aes(x = year, y = global_coef, color = Veg, group = Veg)) +
        .ts_layers() +
        geom_line(linewidth = 1) +
        geom_point(show.legend = FALSE) +
        geom_ribbon(aes(ymin = coef_025, ymax = coef_975, fill = Veg), alpha = 0.12, color = NA) +
        scale_color_manual(values = SPECIES_COLORS) +
        scale_fill_manual(values = SPECIES_COLORS) +
        labs(title = paste0("Inference ", method, ": Species (no Woody aggregation)"),
             x = "Year", y = "Total Normalized Fraction", color = "Veg", fill = "Veg") +
        theme_minimal()
      ggsave(file.path(OUT_DIR, paste0("inference_", file_prefix, "_species_separate.png")), p_sp_ts, width = 8, height = 6)
      readr::write_csv(species_data, file.path(OUT_DIR, paste0("inference_", file_prefix, "_species_separate.csv")))
      cat(sprintf("Saved inference %s species (no woody aggregation) plot to: %s\n", method, file.path(OUT_DIR, paste0("inference_", file_prefix, "_species_separate.png"))))

      # Species stacked proportion
      df_wide_sp <- tryCatch({ tidyr::pivot_wider(species_data |> dplyr::select(year, Veg, global_coef), names_from = Veg, values_from = global_coef) }, error = function(e) NULL)
      if (!is.null(df_wide_sp)) {
        sp_cols <- intersect(c("Herbs","Populus","Tamarix","Woody Unknown"), names(df_wide_sp))
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

  # 5. Stacked area (Proportion) and Woody/Herbs ratio
  df_wide <- tryCatch({
    tidyr::pivot_wider(herbs_woody |> dplyr::select(year, Veg, global_coef), names_from = Veg, values_from = global_coef)
  }, error = function(e) NULL)
  if (!is.null(df_wide) && all(c("Herbs","Woody") %in% names(df_wide))) {
    df_wide$Herbs <- as.numeric(df_wide$Herbs)
    df_wide$Woody <- as.numeric(df_wide$Woody)
    df_wide$total <- rowSums(df_wide[, c("Herbs","Woody")], na.rm = TRUE)

    df_prop <- df_wide |> dplyr::filter(is.finite(total) & total > 0) |>
      dplyr::mutate(Herbs = Herbs/total, Woody = Woody/total) |>
      tidyr::pivot_longer(cols = c("Herbs","Woody"), names_to = "Veg", values_to = "prop")
    p_stacked <- ggplot(df_prop, aes(x = year, y = prop, fill = Veg)) +
      geom_area() +
      scale_fill_manual(values = HERBS_WOODY_COLORS) +
      labs(title = paste0("Inference ", method, ": Herbs vs Woody (Proportion, stacked)"),
           x = "Year", y = "Proportion", fill = "Type") +
      theme_minimal()
    ggsave(file.path(OUT_DIR, paste0("inference_", file_prefix, "_herbs_vs_woody_stacked.png")), p_stacked, width = 8, height = 6)
    readr::write_csv(df_prop, file.path(OUT_DIR, paste0("inference_", file_prefix, "_herbs_vs_woody_stacked.csv")))
    cat(sprintf("Saved inference %s herbs vs woody stacked plot to: %s\n", method, file.path(OUT_DIR, paste0("inference_", file_prefix, "_herbs_vs_woody_stacked.png"))))

    df_ratio <- df_wide |> dplyr::mutate(ratio = ifelse(is.finite(Herbs) & Herbs > 0, Woody / Herbs, NA_real_))
    p_ratio <- ggplot(df_ratio, aes(x = year, y = ratio)) +
      geom_line(color = "#8B4513", linewidth = 1) +
      geom_point() +
      labs(title = paste0("Inference ", method, ": Woody / Herbs Ratio"),
           x = "Year", y = "Woody / Herbs") +
      theme_minimal()
    ggsave(file.path(OUT_DIR, paste0("inference_", file_prefix, "_woody_over_herbs.png")), p_ratio, width = 8, height = 6)
    readr::write_csv(df_ratio, file.path(OUT_DIR, paste0("inference_", file_prefix, "_woody_over_herbs.csv")))
    cat(sprintf("Saved inference %s woody/herbs ratio plot to: %s\n", method, file.path(OUT_DIR, paste0("inference_", file_prefix, "_woody_over_herbs.png"))))
  } else {
    cat(sprintf("Cannot create inference %s stacked/ratio plot: missing Herbs/Woody rows.\n", method))
  }

  # 6. Tamarix vs Populus vs Woody_unknown stacked plot (PPI only)
  if (include_woody_types_plot) {
    woody_types <- full_data[veg_norm %in% c("tamarix", "populus", "woody_unknown"), ]
    if (nrow(woody_types) > 0) {
      woody_types$Veg <- dplyr::case_when(
        tolower(woody_types$Veg) == "tamarix" ~ "Tamarix",
        tolower(woody_types$Veg) == "populus" ~ "Populus",
        tolower(woody_types$Veg) == "woody_unknown" ~ "Woody Unknown",
        TRUE ~ woody_types$Veg
      )
      df_wide_woody <- tryCatch({
        tidyr::pivot_wider(woody_types |> dplyr::select(year, Veg, global_coef), names_from = Veg, values_from = global_coef)
      }, error = function(e) NULL)
      if (!is.null(df_wide_woody)) {
        woody_cols <- intersect(c("Tamarix", "Populus", "Woody Unknown"), names(df_wide_woody))
        if (length(woody_cols) >= 2) {
          for (col in woody_cols) {
            df_wide_woody[[col]] <- as.numeric(df_wide_woody[[col]])
            df_wide_woody[[col]][is.na(df_wide_woody[[col]])] <- 0
          }
          df_wide_woody$total <- rowSums(df_wide_woody[, woody_cols, drop = FALSE], na.rm = TRUE)
          df_prop_woody <- df_wide_woody |>
            dplyr::filter(is.finite(total) & total > 0) |>
            dplyr::mutate(across(all_of(woody_cols), ~ . / total)) |>
            tidyr::pivot_longer(cols = all_of(woody_cols), names_to = "Veg", values_to = "prop")
          p_woody_stacked <- ggplot(df_prop_woody, aes(x = year, y = prop, fill = Veg)) +
            geom_area() +
            scale_fill_manual(values = WOODY_TYPES_COLORS) +
            labs(title = paste0("Inference ", method, ": Tamarix vs Populus vs Woody Unknown (Proportion, stacked)"),
                 x = "Year", y = "Proportion", fill = "Type") +
            theme_minimal()
          ggsave(file.path(OUT_DIR, paste0("inference_", file_prefix, "_tamarix_populus_woody_unknown_stacked.png")), p_woody_stacked, width = 8, height = 6)
          readr::write_csv(df_prop_woody, file.path(OUT_DIR, paste0("inference_", file_prefix, "_tamarix_populus_woody_unknown_stacked.csv")))
          cat(sprintf("Saved inference %s Tamarix vs Populus vs Woody Unknown stacked plot to: %s\n", method, file.path(OUT_DIR, paste0("inference_", file_prefix, "_tamarix_populus_woody_unknown_stacked.png"))))
        } else {
          cat(sprintf("Cannot create inference %s Tamarix/Populus/Woody_unknown stacked plot: need at least 2 of these veg types.\n", method))
        }
      } else {
        cat(sprintf("Cannot create inference %s Tamarix/Populus/Woody_unknown stacked plot: pivot failed.\n", method))
      }
    } else {
      cat(sprintf("No Tamarix/Populus/Woody_unknown data available for inference %s woody types stacked plot.\n", method))
    }
  }

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
    df$TCW <- 0.1511*b$blue + 0.1973*b$green + 0.3283*b$red + 0.3407*b$nir - 0.7117*b$swir1 - 0.4559*b$swir2
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
# NOTE: We previously L2-normalized per-index (blockwise). That was incorrect for
# this project: we want per-observation brightness/scale invariance.
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

# Remove large outliers robustly per (location_id, pheno_year) where possible, otherwise per-location.
# Uses spline-based outlier detection for groups with sufficient data, otherwise falls back to MAD.
remove_large_outliers <- function(df, candidates = NULL, mad_thresh = OUTLIER_MAD_THRESHOLD) {
  if (!isTRUE(ENABLE_OUTLIER_REMOVAL)) return(df)
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

  joined <- dplyr::left_join(df, gpts_map, by = "location_id", suffix = c("", ".geo"))
  if ("Veg.geo" %in% names(joined)) {
    joined$Veg <- ifelse(is.na(joined$Veg) | joined$Veg == "", joined$Veg.geo, joined$Veg)
    joined$Veg.geo <- NULL
  }


  post_non_na <- sum(!is.na(joined$Veg) & joined$Veg != "")

  if (post_non_na == pre_non_na && "location_row" %in% names(gpts_map)) {
    df_ids <- unique(na.omit(as.character(df$location_id)))
    match_count <- length(intersect(df_ids, unique(na.omit(as.character(gpts_map$location_row)))))
    if (match_count > 0) {
      cat(sprintf("[NOTICE] No matches by 'location_id' — attempting join by row-number mapping (matched ids=%d)\n", match_count))
      joined2 <- dplyr::left_join(df, gpts_map, by = c("location_id" = "location_row"), suffix = c("", ".geo"))
      if ("Veg.geo" %in% names(joined2)) {
        joined2$Veg <- ifelse(is.na(joined2$Veg) | joined2$Veg == "", joined2$Veg.geo, joined2$Veg)
        joined2$Veg.geo <- NULL
      }

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
  if (!"PPI" %in% names(df) || all(!is.finite(df$PPI))) {
    dvi_soil_vec <- compute_dvi_soil_per_location(df)
    df <- add_ppi_columns(df, dvi_soil = dvi_soil_vec)
    cat("[PPI] Added PPI to dataset before normalization (per-location dvi_soil + per-location M).\n")
  }
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

ppi_max <- if (exists("PPI_FULL_VEG_COVER")) get("PPI_FULL_VEG_COVER") else 0.4
df <- backup_and_normalize_ppi(df, ppi_max)

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

if (exists("TRAIN_YEARS") && !is.null(TRAIN_YEARS) && length(TRAIN_YEARS) > 0) {
  cat(sprintf("Filtering training data to phenological years (March-February): %s\n", paste(TRAIN_YEARS, collapse = ", ")))
  df_train <- df[df$pheno_year %in% TRAIN_YEARS, , drop = FALSE]
} else {
  df_train <- df
}
# Apply balancing/downsampling to df_train for balanced training library

# AGENT: Apply balancing/downsampling to df_train ONLY (ensure balanced training library)
set.seed(get_mesma_seed(4))
class_counts <- table(df_train$Veg)
if (length(class_counts) > 0) {
  # Downsample vegetation classes to the minimum count among NON-barren classes.
  # This prevents barren from driving the downsampling target.
  non_barren_counts <- class_counts[names(class_counts) != "barren"]
  if (length(non_barren_counts) > 0) {
    min_count_non_barren <- min(non_barren_counts)
    if (is.finite(min_count_non_barren) && min_count_non_barren > 0) {
      df_train_non_barren <- df_train %>%
        dplyr::filter(.data$Veg != "barren") %>%
        dplyr::group_by(.data$Veg) %>%
        dplyr::slice_sample(n = min_count_non_barren) %>%
        dplyr::ungroup()

      df_train_barren <- df_train %>% dplyr::filter(.data$Veg == "barren")
      df_train <- dplyr::bind_rows(df_train_non_barren, df_train_barren)

      cat(sprintf(
        "[BALANCE] Downsampled NON-barren classes to %d samples/class; kept barren as-is (total=%d)\n",
        min_count_non_barren, nrow(df_train)
      ))
    }
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

  # Get unique location-Veg pairs.
  # Note: A location typically has one dominant veg type, but might vary.
  # We'll assign each location to its most frequent Veg type for stratification.
  loc_veg_summary <- df_train %>%
    dplyr::group_by(location_id, Veg) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(location_id, dplyr::desc(n)) %>%
    dplyr::group_by(location_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

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

  train_loc_veg <- df_train %>%
    dplyr::group_by(location_id, Veg) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(location_id, dplyr::desc(n)) %>%
    dplyr::group_by(location_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

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

  lda_res <- safe_lda_call(pca_res$x[, 1:n_pcs, drop=FALSE], as.factor(y_labels), min_n_pcs = min_pcs_for_lda)

  if (is.null(lda_res)) {
    stop("[LDA] LDA could not be computed (collinearity / too few samples / invalid feature space)")
  }
  
  W_pc <- lda_res$scaling
  R <- pca_res$rotation[, 1:n_pcs, drop=FALSE]
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

  return(list(
    means = global_means,
    sds = global_sds,
    weights = final_weights,
    indices = all_feature_cols,
    base_indices = feature_cols,
    l2_normalize = l2_only_mode,
    zscore_applied = apply_zscore
  ))
}

doy_to_pentad <- function(doy) {
  pmin(ceiling(doy / TEMPORAL_AGGREGATION_DAYS), TEMPORAL_BUDGET)
}

build_pentad_matrix <- function(dly_year, avail_idx, interpolate = TRUE) {
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

  K <- length(avail_idx)
  pentad_mat <- matrix(NA_real_, nrow = TEMPORAL_BUDGET, ncol = K)
  colnames(pentad_mat) <- avail_idx

  for (p in 1:TEMPORAL_BUDGET) {
    subset_p <- dly_year[dly_year$pentad == p, ]
    if (nrow(subset_p) == 0) next

    # Pentad center in phenological DOY space (March 1 = 1), using the nominal bin boundaries.
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

  # Only interpolate missing values if interpolate=TRUE (for training endmembers)
  # For inference/validation data, keep NAs to use only actual observations
  if (interpolate) {
    for (j in 1:K) {
      vals <- pentad_mat[, j]
      if (any(is.na(vals))) {
        if (all(is.na(vals))) {
          pentad_mat[, j] <- 0
        } else {
          idx_present <- which(!is.na(vals))
          if (length(idx_present) >= 2) {
            pentad_mat[, j] <- approx(idx_present, vals[idx_present], xout = 1:TEMPORAL_BUDGET, rule = 2)$y
          } else {
            pentad_mat[, j] <- vals[idx_present[1]]
          }
        }
      }
    }
  }

  pentad_mat
}

 

# apply_pca_lda_transform removed (was a NO-OP returning y unchanged).
# PCA-LDA weights are applied only in the solver via feature_weights.


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
df_train <- dplyr::left_join(df_train, gpts_map, by = "location_id", suffix = c("", ".geo"))
if ("Veg.geo" %in% names(df_train)) {
  df_train$Veg <- ifelse(is.na(df_train$Veg) | df_train$Veg == "", df_train$Veg.geo, df_train$Veg)
  df_train$Veg.geo <- NULL
}

df_train <- canonicalize_veg_labels(df_train)

if (!"date" %in% names(df)) stop("Input CSV must contain a 'date' column")
df$date <- as.Date(df$date)
if (!"location_id" %in% names(df)) stop("Input CSV must contain a 'location_id' column")

df <- canonicalize_veg_labels(df)

df$doy <- pheno_doy(df$date)  # Use phenological DOY (March 1 = day 1)
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
    okc <- is.finite(cand$lat) & is.finite(cand$lon)
    cand <- cand[order(cand$location_id, -as.integer(okc)), , drop = FALSE]
    cand <- cand[!duplicated(cand$location_id), , drop = FALSE]
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
      okg <- is.finite(gm2$lat) & is.finite(gm2$lon)
      gm2 <- gm2[order(gm2$location_id, -as.integer(okg)), , drop = FALSE]
      gm2 <- gm2[!duplicated(gm2$location_id), , drop = FALSE]
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
  range_est <- tryCatch({
    nls_fit <- nls(g ~ s * (1 - exp(-d / r)),
                   data = data.frame(d = bin_mid[valid_bins], g = bin_gamma[valid_bins]),
                   start = list(s = total_var, r = median(dists)),
                   lower = list(s = total_var * 0.1, r = max(dists) * 0.01),
                   upper = list(s = total_var * 3, r = max(dists) * 2),
                   algorithm = "port",
                   control = list(maxiter = 50, warnOnly = TRUE))
    coef(nls_fit)["r"]
  }, error = function(e) {
    # Fallback: distance at which semivariance first exceeds 0.5 * sill
    thresh_idx <- which(bin_gamma[valid_bins] >= 0.5 * total_var)
    if (length(thresh_idx) > 0) bin_mid[valid_bins][thresh_idx[1]]
    else max(dists)  # All correlated -> conservative
  })

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
    veg_boot_res <- aggregate_woody_bootstrap(veg_boot_res, label = "NOINDEX BOOTSTRAP")
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

  veg_boot_res <- aggregate_woody_bootstrap(veg_boot_res, label = paste0(index_name, " BOOTSTRAP"))
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

# Precompute B Dirichlet draws for every location-year so that each bootstrap
# replicate gets a fresh classification perturbation.
# Returns list(draws = named list of [B x K] matrices, alpha_info = list of alpha metadata).
precompute_dirichlet_draws <- function(all_coefs, B, conf_matrix, sample_sizes = NULL,
                                       concentration_scale = DIRICHLET_CONCENTRATION_SCALE) {
  loc_yr_keys <- unique(paste(all_coefs$location_id, all_coefs$pheno_year, sep = "___"))
  row_names <- tolower(rownames(conf_matrix))
  col_names <- colnames(conf_matrix)
  K <- ncol(conf_matrix)

  draws <- list()

  for (key in loc_yr_keys) {
    parts <- strsplit(key, "___")[[1]]
    loc_id <- parts[1]
    p_year <- as.numeric(parts[2])

    mask <- all_coefs$location_id == loc_id & all_coefs$pheno_year == p_year
    loc_yr_data <- all_coefs[mask, ]
    if (nrow(loc_yr_data) == 0) next

    # Build fraction vector
    fracs <- setNames(pmax(as.numeric(loc_yr_data$coef), 0),
                      tolower(as.character(loc_yr_data$Veg)))
    fracs[!is.finite(fracs)] <- 0
    frac_sum <- sum(fracs)
    if (frac_sum <= 0) next
    fracs <- fracs / frac_sum

    # Map fraction names to confusion-matrix rows
    mapped_row <- rep(NA_character_, length(fracs))
    for (i in seq_along(fracs)) {
      nm <- names(fracs)[i]
      nm2 <- sub("^frac_", "", nm)
      if (nm %in% row_names) mapped_row[i] <- nm
      else if (nm2 %in% row_names) mapped_row[i] <- nm2
    }

    keep <- is.finite(fracs) & fracs > 0 & !is.na(mapped_row) & fracs >= 1e-6
    if (!any(keep)) next

    w <- fracs[keep]; w <- w / sum(w)
    rows <- mapped_row[keep]
    row_idx <- match(rows, row_names)

    # Weighted mean confusion row
    mu <- rep(0, K)
    for (k in seq_along(row_idx)) {
      mu <- mu + w[k] * as.numeric(conf_matrix[row_idx[k], ])
    }
    mu[!is.finite(mu)] <- 0
    mu[mu < 1e-6] <- 1e-6
    mu <- mu / sum(mu)

    # Concentration from validation support
    n_eff_val <- 20
    if (!is.null(sample_sizes) && length(sample_sizes) > 0) {
      ss_names <- tolower(names(sample_sizes))
      ss <- as.numeric(sample_sizes)
      ss[!is.finite(ss) | ss < 0] <- 0
      ss_match <- ss[match(rows, ss_names)]
      ss_match[is.na(ss_match)] <- 0
      if (sum(ss_match) > 0) n_eff_val <- sum(w * ss_match)
    }

    alpha <- max(concentration_scale * n_eff_val, 5)
    alpha_vec <- alpha * mu

    # Vectorized B Dirichlet draws via rgamma
    gamma_mat <- matrix(rgamma(B * K, shape = rep(alpha_vec, each = B), rate = 1),
                        nrow = B, ncol = K)
    row_sums <- rowSums(gamma_mat)
    zero_rows <- row_sums < 1e-10
    if (any(zero_rows)) {
      gamma_mat[zero_rows, ] <- matrix(rep(mu, sum(zero_rows)),
                                       nrow = sum(zero_rows), byrow = TRUE)
      row_sums[zero_rows] <- 1
    }
    dirichlet_mat <- gamma_mat / row_sums
    colnames(dirichlet_mat) <- tolower(col_names)

    draws[[key]] <- dirichlet_mat
  }

  draws
}

# Location-based bootstrap for global aggregation
location_bootstrap_aggregate <- function(all_coefs, B = BOOTSTRAP_B, seed = 123) {
  set.seed(seed)

  # =========================================================================
  # CLASSIFICATION UNCERTAINTY: Precompute B Dirichlet draws per location-year
  # so each bootstrap replicate gets a fresh classification perturbation.
  # =========================================================================
  dirichlet_draws <- NULL
  if (isTRUE(ENABLE_CLASSIFICATION_UNCERTAINTY) &&
      exists(".CONFUSION_MATRIX", envir = globalenv())) {
    conf_matrix <- get(".CONFUSION_MATRIX", envir = globalenv())
    sample_sizes <- if (exists(".VALIDATION_SAMPLE_SIZES", envir = globalenv())) {
      get(".VALIDATION_SAMPLE_SIZES", envir = globalenv())
    } else NULL

    message("[BOOTSTRAP] Classification uncertainty enabled; precomputing B Dirichlet draws per location-year.")
    dirichlet_draws <- precompute_dirichlet_draws(all_coefs, B, conf_matrix, sample_sizes)
    message(sprintf("[BOOTSTRAP] Precomputed Dirichlet draws for %d location-years", length(dirichlet_draws)))
  } else if (isTRUE(ENABLE_CLASSIFICATION_UNCERTAINTY)) {
    message("[BOOTSTRAP] Classification uncertainty enabled but no confusion matrix available; skipping.")
  }

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
              bin_breaks <- quantile(dists, probs = seq(0, 1, length.out = n_bins + 1))
              bin_breaks <- unique(bin_breaks)
              if (length(bin_breaks) >= 3) {
                bin_mid <- (bin_breaks[-length(bin_breaks)] + bin_breaks[-1]) / 2
                bin_gamma <- numeric(length(bin_mid))
                for (bb in seq_along(bin_mid)) {
                  in_bin <- dists >= bin_breaks[bb] & dists < bin_breaks[bb + 1]
                  if (sum(in_bin) > 0) bin_gamma[bb] <- median(gamma_vals[in_bin])
                  else bin_gamma[bb] <- NA
                }
                valid_bins <- !is.na(bin_gamma)
                if (sum(valid_bins) >= 2) {
                  # Fit exponential variogram: gamma(d) = sill * (1 - exp(-d/range))
                  # Using NLS with reasonable starting values
                  fit_ok <- tryCatch({
                    nls_fit <- nls(g ~ s * (1 - exp(-d / r)),
                                   data = data.frame(d = bin_mid[valid_bins], g = bin_gamma[valid_bins]),
                                   start = list(s = total_var, r = median(dists)),
                                   lower = list(s = total_var * 0.1, r = max(dists) * 0.01),
                                   upper = list(s = total_var * 3, r = max(dists) * 2),
                                   algorithm = "port",
                                   control = list(maxiter = 50, warnOnly = TRUE))
                    range_est <- coef(nls_fit)["r"]
                    TRUE
                  }, error = function(e) FALSE)

                  if (!fit_ok) {
                    # Fallback: use distance at which semivariance first exceeds 0.5 * total_var
                    thresh_idx <- which(bin_gamma[valid_bins] >= 0.5 * total_var)
                    if (length(thresh_idx) > 0) {
                      range_est <- bin_mid[valid_bins][thresh_idx[1]]
                    } else {
                      range_est <- max(dists)  # All correlated -> conservative
                    }
                  }

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
        # Build per-replicate Dirichlet-perturbed coefficients if available.
        # perturbed_coef_mat is [B x n_obs]: each row b gives Dirichlet draw b
        # for each location's coefficient of this veg type.
        perturbed_coef_mat <- NULL
        if (!is.null(dirichlet_draws) && n_obs > 0) {
          loc_keys <- paste(yr_data$location_id, yr_data$pheno_year, sep = "___")
          veg_lower <- tolower(as.character(veg))
          perturbed_coef_mat <- matrix(NA_real_, nrow = B, ncol = n_obs)
          for (j in seq_len(n_obs)) {
            draw_mat <- dirichlet_draws[[loc_keys[j]]]
            if (!is.null(draw_mat)) {
              col_idx <- match(veg_lower, colnames(draw_mat))
              if (!is.na(col_idx)) {
                perturbed_coef_mat[, j] <- draw_mat[, col_idx]
                next
              }
            }
            perturbed_coef_mat[, j] <- yr_data$coef[j]
          }
        }

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

             if (!is.null(perturbed_coef_mat)) {
               # Each replicate gets a Dirichlet-perturbed mean + parametric noise
               for (b in 1:B) {
                 perturbed_mu <- mean(perturbed_coef_mat[b, ], na.rm = TRUE)
                 boot_means[b, i] <- rnorm(1, mean = perturbed_mu, sd = se_mean)
               }
             } else {
               boot_means[, i] <- rnorm(B, mean = mu, sd = se_mean)
             }

        } else {
             # --- SUFFICIENT SAMPLE: USE RESAMPLING ---
             if (!is.null(perturbed_coef_mat)) {
               for (b in 1:B) {
                 idx <- sample(seq_len(n_obs), n_obs, replace = TRUE)
                 boot_means[b, i] <- mean(perturbed_coef_mat[b, idx], na.rm = TRUE)
               }
             } else {
               for (b in 1:B) {
                 boot_sample <- sample(yr_data$coef, n_obs, replace = TRUE)
                 boot_means[b, i] <- mean(boot_sample, na.rm = TRUE)
               }
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
        df_inf <- remove_large_outliers(df_inf)
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

      ppi_max_inf <- if (exists("PPI_FULL_VEG_COVER")) get("PPI_FULL_VEG_COVER") else 0.4
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

    # Check if Huber loss is enabled
    use_huber <- exists("USE_HUBER_LOSS") && isTRUE(USE_HUBER_LOSS)

    if (use_huber) {
      # === HUBER LOSS via IRLS (Iteratively Reweighted Least Squares) ===
      # Huber loss: L(r) = 0.5*r^2 if |r| <= delta, else delta*(|r| - 0.5*delta)
      # IRLS weights: w_i = 1 if |r_i| <= delta, else delta/|r_i|

      huber_delta <- if (exists("HUBER_DELTA")) HUBER_DELTA else 1.345
      huber_max_iter <- if (exists("HUBER_MAX_ITER")) HUBER_MAX_ITER else 20
      huber_tol <- if (exists("HUBER_TOL")) HUBER_TOL else 1e-4

      # Initialize with uniform weights for IRLS
      irls_weights <- rep(1, n_bands)
      w_qp_prev <- rep(1/n_endmembers, n_endmembers)

      for (irls_iter in 1:huber_max_iter) {
        # Combine base weights with IRLS weights
        combined_weights <- base_weights * irls_weights
        E_w <- E_fit * combined_weights
        y_w <- y_fit * combined_weights

        # Solve weighted QP
        Dmat <- 2 * crossprod(E_w)
        dvec <- 2 * crossprod(E_w, y_w)

        qp_scale <- mean(diag(Dmat))
        if(is.na(qp_scale) || qp_scale < 1e-12) qp_scale <- 1.0

        Dmat <- Dmat / qp_scale
        dvec <- dvec / qp_scale

        ridge <- 1e-6
        Dmat <- Dmat + diag(n_endmembers) * ridge
        Dmat <- (Dmat + t(Dmat)) / 2

        Amat <- cbind(rep(1, n_endmembers), diag(n_endmembers))
        bvec <- c(1, rep(0, n_endmembers))

        res_qp <- tryCatch({
          quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
        }, error = function(e) e)

        if (inherits(res_qp, "error")) {
          stop(sprintf("[QP] quadprog::solve.QP failed in solve_weights_fcls (Huber IRLS): %s", res_qp$message))
        }

        w_qp <- res_qp$solution
        w_qp[!is.finite(w_qp)] <- 0
        w_qp[w_qp < 0] <- 0
        w_sum <- sum(w_qp, na.rm = TRUE)
        if(w_sum > 0) w_qp <- w_qp / w_sum else w_qp <- rep(1/n_endmembers, n_endmembers)

        # Calculate residuals
        pred <- as.numeric(E_fit %*% w_qp)
        resid <- y_fit - pred

        # Check convergence
        if (max(abs(w_qp - w_qp_prev)) < huber_tol) {
          break
        }
        w_qp_prev <- w_qp

        # Update IRLS weights based on Huber loss
        # Scale residuals by MAD for robust delta scaling
        mad_resid <- mad(resid, na.rm = TRUE)
        if (!is.finite(mad_resid) || mad_resid < 1e-10) mad_resid <- 1
        scaled_resid <- abs(resid) / mad_resid

        # Huber weights: 1 for small residuals, delta/|r| for large ones
        irls_weights <- ifelse(scaled_resid <= huber_delta,
                               1,
                               huber_delta / pmax(scaled_resid, 1e-10))
        irls_weights[!is.finite(irls_weights)] <- 1
      }

      # Final solution
      pred <- as.numeric(E_fit %*% w_qp)
      resid <- y_fit - pred

      # Compute Huber loss instead of RMSE
      mad_resid <- mad(resid, na.rm = TRUE)
      if (!is.finite(mad_resid) || mad_resid < 1e-10) mad_resid <- 1
      scaled_resid <- abs(resid) / mad_resid
      huber_losses <- ifelse(scaled_resid <= huber_delta,
                             0.5 * (resid/mad_resid)^2,
                             huber_delta * (scaled_resid - 0.5 * huber_delta))
      loss <- sqrt(mean(huber_losses))  # Return sqrt for consistency with RMSE scale

      return(list(w = w_qp, rmse = loss, residuals = resid, loss_type = "huber"))

    } else {
      # === STANDARD RMSE (Original behavior) ===
      E_w <- E_fit * base_weights
      y_w <- y_fit * base_weights

      # --- DEDICATED QP SOLVER (ALWAYS used) ---
      # Problem: min ||Ex - y||^2  s.t. sum(x)=1, x>=0
      # ||Ex - y||^2 = x'E'Ex - 2y'Ex + y'y
      # quadprog solves: min 1/2 b^T Dmat b - dvec^T b
      # Dmat = 2 * E'E, dvec = 2 * E'y

      Dmat <- 2 * crossprod(E_w)
      dvec <- 2 * crossprod(E_w, y_w)

      # Scale QP to prevent inconsistent constraints
      qp_scale <- mean(diag(Dmat))
      if(is.na(qp_scale) || qp_scale < 1e-12) qp_scale <- 1.0

      Dmat <- Dmat / qp_scale
      dvec <- dvec / qp_scale

      # Add stronger ridge to diagonal for numerical stability (PD requirement)
      ridge <- 1e-6
      Dmat <- Dmat + diag(n_endmembers) * ridge

      # Ensure Dmat is symmetric
      Dmat <- (Dmat + t(Dmat)) / 2

      # Constraints: A^T b >= b_0
      # 1. sum(x) = 1  => use meq=1. Row 1 of A^T is [1, 1, ..., 1].
      # 2. x >= 0      => I * x >= 0. Rows 2..N+1 of A^T are I.

      Amat <- cbind(rep(1, n_endmembers), diag(n_endmembers))
      bvec <- c(1, rep(0, n_endmembers))

      res_qp <- tryCatch({
        quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
      }, error = function(e) {
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG] Quadprog failed in solve_weights_fcls: %s\n", e$message))
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG] Dmat dim: %s, dvec length: %d, Amat dim: %s, bvec length: %d\n",
                    paste(dim(Dmat), collapse="x"), length(dvec), paste(dim(Amat), collapse="x"), length(bvec)))
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG] Dmat eigenvalues: %s\n", paste(eigen(Dmat)$values[1:5], collapse=", ")))
        e
      })

      if (!inherits(res_qp, "error")) {
        w_qp <- res_qp$solution
        w_qp[!is.finite(w_qp)] <- 0
        w_qp[w_qp < 0] <- 0
        w_sum <- sum(w_qp, na.rm = TRUE)
        if(w_sum > 0) w_qp <- w_qp / w_sum else w_qp <- rep(1/n_endmembers, n_endmembers)

        # Calculate RMSE for consistency
        pred <- as.numeric(E_fit %*% w_qp)
        resid <- y_fit - pred
        rmse <- sqrt(mean(resid^2))

        return(list(w = w_qp, rmse = rmse, residuals = resid, loss_type = "rmse"))
      }

      stop(sprintf("[QP] quadprog::solve.QP failed in solve_weights_fcls: %s", res_qp$message))
    }
  }

  # Batch FCLS Solver using Quadprog (Optimization for GA/Grid Search)
  # E: Features x Endmembers matrix
  # Y: Samples x Features matrix
  # Returns: Samples x Endmembers weight matrix
  solve_batch_fcls <- function(E, Y, feature_weights = NULL) {
    if (is.null(E) || ncol(E) < 1) return(NULL)
    n_endmembers <- ncol(E)
    n_samples <- nrow(Y)
    if (n_samples == 0) return(matrix(0, 0, n_endmembers))

    # Convert Y to Features x Samples for matrix math
    Y_t <- t(Y)

    # Determine E and Y used for fitting (weighted or raw)
    if (!is.null(feature_weights) && length(feature_weights) == nrow(E)) {
      w <- pmax(feature_weights, 0)
      w[!is.finite(w)] <- 0

      # Apply weights directly: E_w = W * E, Y_w = W * Y
      # Since W is diagonal, just multiply rows
      E_fit <- E * w
      Y_fit <- Y_t * w
    } else {
      E_fit <- E
      Y_fit <- Y_t
    }

    # Precompute constant QP matrices
    # Dmat = 2 * E'E
    Dmat <- 2 * crossprod(E_fit)

    # Scale QP problem to improve numerical stability for quadprog
    qp_scale <- mean(diag(Dmat))
    if(is.na(qp_scale) || qp_scale < 1e-12) qp_scale <- 1.0

    Dmat <- Dmat / qp_scale

    # Add stronger ridge
    ridge <- 1e-6 # Relative to scaled matrix
    Dmat <- Dmat + diag(n_endmembers) * ridge
    # Ensure symmetric
    Dmat <- (Dmat + t(Dmat)) / 2

    # Amat: [1s; I]
    Amat <- cbind(rep(1, n_endmembers), diag(n_endmembers))
    bvec <- c(1, rep(0, n_endmembers))

    # Precompute all linear terms: dvec = 2 * E'y
    # Dvecs: Endmembers x Samples
    # Must also be scaled by same factor
    Dvecs <- (2 * crossprod(E_fit, Y_fit)) / qp_scale

    w_out <- matrix(0, nrow=n_samples, ncol=n_endmembers)

    # Loop over samples - solve.QP is fast when Dmat is precomputed
    for(i in 1:n_samples) {
      dvec <- Dvecs[, i]

      res <- tryCatch({
        quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq=1)
      }, error = function(e) {
        e
      })

      if(!inherits(res, "error")) {
        w <- res$solution
        w[!is.finite(w)] <- 0
        w[w < 0] <- 0
        s <- sum(w)
        if(s > 0) w <- w / s else w <- rep(1/n_endmembers, n_endmembers)
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
    feature_weights = NULL
  ) {
    if (length(top_variants) == 0) return(NULL)

    full_veg_names <- names(top_variants)
    n_veg_full <- length(full_veg_names)
    if (n_veg_full == 0) return(NULL)

    y_target <- y
    y_target[!is.finite(y_target)] <- 0

    score_from_solution <- function(res) {
      if (is.null(res) || is.null(res$residuals)) return(Inf)
      as.numeric(res$rmse)
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

      res <- solve_weights_fcls(E, y_target, feature_weights = if (!is.null(feature_weights) && length(feature_weights) == nrow(E)) feature_weights else NULL)
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
      res$score <- score_from_solution(res)
      return(res)
    }

    # --- Always use ALL vegetation types, one prototype per class ---
    # Try all prototype combinations across classes to find the best fit.
    n_variants_per_veg <- sapply(top_variants, length)

    global_best_res <- NULL
    global_best_score <- Inf

    # Single fixed subset: all classes must participate
    veg_subset <- full_veg_names
    n_in_subset <- length(veg_subset)
    n_variants_subset <- n_variants_per_veg[veg_subset]
    total_combos <- as.numeric(prod(n_variants_subset))

    if (total_combos <= 500) {
      # Full exhaustive search (typical: max 4^4 = 256 with MAX_K_EAR=4)
      combos <- expand.grid(lapply(n_variants_subset, seq_len), KEEP.OUT.ATTRS = FALSE)
      for (i in seq_len(nrow(combos))) {
        r <- solve_combo_subset(veg_subset, as.integer(combos[i, ]))
        if (!is.null(r) && r$score < global_best_score) {
          global_best_score <- r$score
          global_best_res <- r
        }
      }
    } else {
      # Coarse search + coordinate descent for large search spaces
      TOPK_COARSE <- 3
      coarse_idx <- lapply(n_variants_subset, function(n) seq_len(min(TOPK_COARSE, n)))
      combos <- expand.grid(coarse_idx, KEEP.OUT.ATTRS = FALSE)
      best_idx <- rep(1, n_in_subset)

      for (i in seq_len(nrow(combos))) {
        r <- solve_combo_subset(veg_subset, as.integer(combos[i, ]))
        if (!is.null(r) && r$score < global_best_score) {
          global_best_score <- r$score
          global_best_res <- r
          best_idx <- as.integer(combos[i, ])
        }
      }

      # Coordinate descent refinement
      curr_idx <- best_idx
      improved <- TRUE
      while (improved) {
        improved <- FALSE
        for (k in seq_len(n_in_subset)) {
          for (v_opt in seq_len(n_variants_subset[k])) {
            if (v_opt == curr_idx[k]) next
            t_idx <- curr_idx
            t_idx[k] <- v_opt
            r <- solve_combo_subset(veg_subset, t_idx)
            if (!is.null(r) && r$score < global_best_score) {
              global_best_score <- r$score
              global_best_res <- r
              curr_idx <- t_idx
              improved <- TRUE
            }
          }
        }
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
        stop("[ERROR] No OOB samples found in storage_oob. Cannot proceed with cluster optimization.")
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

          # Batch process using solve_batch_fcls for speed
          all_coefs <- solve_batch_fcls(M, Y_test, params$weights)

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
    # Post-construction filtering REMOVED - all filtering now happens pre-construction
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
      # Month start days in phenological DOY (March 1 = day 1)
      pheno_month_starts <- c(1, 32, 62, 93, 123, 154, 184, 215, 245, 276, 306, 337)
      pheno_month_labels <- c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan","Feb")
      # Convert month-start pheno DOY to pentad position: pentad = ceil(doy / agg_days)
      month_breaks_pentad <- ceiling(pheno_month_starts / TEMPORAL_AGGREGATION_DAYS)
      # Keep only unique breaks within the temporal budget
      keep <- month_breaks_pentad >= 1 & month_breaks_pentad <= TEMPORAL_BUDGET & !duplicated(month_breaks_pentad)
      month_breaks_pentad <- month_breaks_pentad[keep]
      month_labels <- pheno_month_labels[keep]

      # Build per-index significance mask from feature_weights (0 = pruned / non-significant)
      sig_mask <- NULL
      if (!is.null(feature_weights) && length(feature_weights) == length(indices) * TEMPORAL_BUDGET) {
        sig_mask <- list()
        for (k in seq_len(length(indices))) {
          w_start <- (k - 1) * TEMPORAL_BUDGET + 1
          w_end <- k * TEMPORAL_BUDGET
          sig_mask[[indices[k]]] <- feature_weights[w_start:w_end] > 0
        }
        # Drop fully-excluded indices (all pentad weights zero)
        excluded <- names(sig_mask)[vapply(sig_mask, function(m) !any(m), logical(1))]
        if (length(excluded) > 0) {
          cat(sprintf("[PROTO PLOT] Skipping %d fully-excluded indices: %s\n", length(excluded), paste(excluded, collapse = ", ")))
          indices <- setdiff(indices, excluded)
          sig_mask <- sig_mask[indices]
        }
      }

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
            sig <- if (!is.null(sig_mask) && idx_name %in% names(sig_mask)) sig_mask[[idx_name]] else rep(TRUE, TEMPORAL_BUDGET)
            df_tmp <- data.frame(pentad = seq_len(TEMPORAL_BUDGET), value = vals, Veg = v, variant_id = vid, index = idx_name, significant = sig, stringsAsFactors = FALSE)
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

        # Determine which pentads are non-significant for this index
        has_nonsig <- "significant" %in% names(df_idx) && any(!df_idx$significant)

        p <- ggplot2::ggplot(df_idx, ggplot2::aes(x = pentad, y = value, group = variant_id))
        p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.4)

        if (has_nonsig) {
          # Add grey shading behind non-significant pentads
          nonsig_pentads <- unique(df_idx$pentad[!df_idx$significant])
          shade_df <- data.frame(xmin = nonsig_pentads - 0.5, xmax = nonsig_pentads + 0.5,
                                 ymin = y_min - y_pad, ymax = y_max + y_pad)
          p <- p + ggplot2::geom_rect(data = shade_df, inherit.aes = FALSE,
                                      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                                      fill = "grey90", alpha = 0.5)
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

        if (frac_zero_glob <= PRUNE_ZERO_WEIGHT_MAX_FRAC && (length(feature_weights) - n_zero_glob) >= PRUNE_ZERO_MIN_FEATURES) {
          keep_idx_global_w <- which(!zero_mask_glob)
          cat(sprintf("[INFO] Will prune %d features for MESMA, keeping %d features\n", n_zero_glob, length(keep_idx_global_w)))
        } else {
          cat(sprintf("[WARN] Skipping MESMA pruning: zero fraction %.2f exceeds max allowed %.2f or resulting features < %d\n", frac_zero_glob, PRUNE_ZERO_WEIGHT_MAX_FRAC, PRUNE_ZERO_MIN_FEATURES))
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

  # Step 3: Alpha-grid pruning using precomputed permutation p-values.
  # IMPORTANT: This step does NOT run permutation tests anymore.
  # It only loads p-values from permutation_importance_results.csv and optimizes
  # (index_alpha, pentad_alpha) via grid search on OOB + training evaluation.
  pruned_indices <- avail
  MESMA_PARAMS <- MESMA_PARAMS_INITIAL

  if (exists("df_train_oob") && !is.null(df_train_oob) && nrow(df_train_oob) > 0 &&
      !is.null(MESMA_PARAMS_INITIAL) && !is.null(MESMA_PARAMS_INITIAL$weights)) {
    perm_path <- "permutation_importance_results.csv"
    if (!file.exists(perm_path)) {
      stop(sprintf("Missing '%s'. Alpha-grid pruning requires precomputed permutation p-values.", perm_path))
    }

    perm_results <- tryCatch({
      read.csv(perm_path, stringsAsFactors = FALSE)
    }, error = function(e) {
      stop(sprintf("Failed to read '%s': %s", perm_path, e$message))
    })

    if (!"p_value" %in% names(perm_results) || !"index" %in% names(perm_results)) {
      stop(sprintf("'%s' must contain columns: index, p_value (and ideally stage, pentad)", perm_path))
    }

    # Use feature list from params for pruning
    perm_indices <- if (!is.null(MESMA_PARAMS_INITIAL$indices)) MESMA_PARAMS_INITIAL$indices else avail

    # Split permutation results into index-level and pentad-level tables
    if (!"stage" %in% names(perm_results)) {
      perm_results$stage <- if ("pentad" %in% names(perm_results) && any(is.finite(suppressWarnings(as.numeric(perm_results$pentad))))) "pentad" else "index"
    }

    index_results <- perm_results[tolower(perm_results$stage) == "index", , drop = FALSE]
    if (nrow(index_results) == 0) {
      # Fallback: treat rows with missing pentad as index-level
      if ("pentad" %in% names(perm_results)) {
        idx_mask <- is.na(suppressWarnings(as.numeric(perm_results$pentad)))
        index_results <- perm_results[idx_mask, , drop = FALSE]
      }
    }
    if (nrow(index_results) == 0) {
      stop(sprintf("'%s' contains no index-level rows (stage=='index' or pentad missing).", perm_path))
    }

    # Keep only indices present in current run
    index_results$index <- as.character(index_results$index)
    index_results <- index_results[index_results$index %in% perm_indices & is.finite(as.numeric(index_results$p_value)), , drop = FALSE]

    # Pentad-level results
    pentad_results <- perm_results[tolower(perm_results$stage) == "pentad", , drop = FALSE]
    if (nrow(pentad_results) == 0 && "pentad" %in% names(perm_results)) {
      # Fallback: any finite pentad is treated as pentad-level
      pj <- suppressWarnings(as.integer(perm_results$pentad))
      pentad_results <- perm_results[is.finite(pj), , drop = FALSE]
      pentad_results$pentad <- pj[is.finite(pj)]
    }
    if (nrow(pentad_results) > 0) {
      pentad_results$index <- as.character(pentad_results$index)
      pentad_results$pentad <- suppressWarnings(as.integer(pentad_results$pentad))
      pentad_results$p_value <- suppressWarnings(as.numeric(pentad_results$p_value))
      pentad_results <- pentad_results[pentad_results$index %in% perm_indices & is.finite(pentad_results$pentad) & is.finite(pentad_results$p_value), , drop = FALSE]
    }

    # Prepare OOB data
    oob_for_eval <- dplyr::mutate(df_train_oob, target_class = tolower(as.character(Veg)))
    oob_for_eval <- dplyr::filter(oob_for_eval, !is.na(target_class) & target_class != "")

    n_oob_samples <- nrow(unique(oob_for_eval[, c("location_id", "pheno_year")]))
    if (n_oob_samples < PERMUTATION_MIN_SAMPLES) {
      stop(sprintf("Only %d OOB samples (< min=%d). Need more OOB for alpha-grid scoring.", n_oob_samples, PERMUTATION_MIN_SAMPLES))
    }

    # Helper: build endmember matrix from initial library
    build_endmember_matrix <- function(mesma_library, unique_classes) {
      cols <- list(); col_names <- character(0)
      for (cls in names(mesma_library)) {
        if (!(cls %in% unique_classes)) next
        cls_variants <- mesma_library[[cls]]
        if (is.null(cls_variants) || length(cls_variants) == 0) next
        for (variant in cls_variants) {
          vec <- variant$vec
          if (!is.null(vec) && length(vec) > 0) {
            cols[[length(cols) + 1]] <- vec
            col_names <- c(col_names, cls)
          }
        }
      }
      if (length(cols) < 2) stop("[ALPHA OPT] Insufficient endmembers in initial library")
      list(E = do.call(cbind, cols), col_names = col_names)
    }

    build_Y_from_df <- function(df_in, params, n_bins) {
      traces <- unique(df_in[, c("location_id", "pheno_year", "target_class")])
      if (nrow(traces) == 0) stop("[ALPHA OPT] No traces found")

      indices <- params$indices
      base_indices <- if (!is.null(params$base_indices)) params$base_indices else indices
      n_base_idx <- length(base_indices)
      l2_normalize <- isTRUE(params$l2_normalize)

      vecs <- vector("list", nrow(traces))
      labels <- as.character(traces$target_class)
      for (j in seq_len(nrow(traces))) {
        lid <- traces$location_id[j]
        pyr <- traces$pheno_year[j]
        sub <- df_in[df_in$location_id == lid & df_in$pheno_year == pyr, ]
        mat <- build_pentad_matrix(sub, base_indices)
        if (is.null(mat)) { vecs[[j]] <- NULL; next }
        vec_raw <- as.numeric(mat)
        vec <- mesma_apply_representation_vec(vec_raw, n_base_idx, n_bins, l2_normalize)
        vec <- mesma_zscore_vec_by_index(vec, indices, params$means, params$sds, n_bins)
        vecs[[j]] <- vec
      }
      keep <- vapply(vecs, function(x) !is.null(x) && length(x) > 0, logical(1))
      vecs <- vecs[keep]; labels <- labels[keep]
      if (length(vecs) == 0) stop("[ALPHA OPT] No vectors built")
      list(Y = do.call(rbind, vecs), labels = labels)
    }

    compute_score_from_Y <- function(E, col_names, Y, labels, unique_classes, weights) {
      all_coefs <- solve_batch_fcls(E, Y, weights)
      if (is.null(all_coefs)) stop("[ALPHA OPT] FCLS solver returned NULL")

      veg_classes <- setdiff(unique_classes, "barren")
      veg_norm_frac_sums <- setNames(rep(0, length(veg_classes)), veg_classes)
      veg_counts <- setNames(rep(0L, length(veg_classes)), veg_classes)

      total <- nrow(Y)
      for (j in seq_len(total)) {
        true_cls <- labels[j]
        if (!(true_cls %in% veg_classes)) next
        coefs <- all_coefs[j, ]
        class_sums <- tapply(coefs, col_names, sum)
        for (uc in unique_classes) if (!(uc %in% names(class_sums))) class_sums[[uc]] <- 0
        veg_total <- sum(sapply(veg_classes, function(vc) class_sums[[vc]]), na.rm = TRUE)
        norm_frac <- if (veg_total > 1e-10) class_sums[[true_cls]] / veg_total else 0
        veg_norm_frac_sums[true_cls] <- veg_norm_frac_sums[true_cls] + norm_frac
        veg_counts[true_cls] <- veg_counts[true_cls] + 1L
      }

      veg_diag_fracs <- ifelse(veg_counts > 0, veg_norm_frac_sums / veg_counts, NA)
      mean(veg_diag_fracs, na.rm = TRUE)
    }

    apply_pruning_rules <- function(weights, index_names, n_bins, idx_prune, pent_prune_df) {
      if (is.null(weights) || length(weights) == 0) return(weights)

      # Stage 1: prune entire indices
      if (!is.null(idx_prune) && length(idx_prune) > 0) {
        for (idx_name in idx_prune) {
          idx_i <- which(index_names == idx_name)
          if (length(idx_i) == 1) {
            s <- (idx_i - 1) * n_bins + 1
            e <- idx_i * n_bins
            weights[s:e] <- 0
          }
        }
      }

      # Stage 2: prune specific pentads (only those explicitly listed)
      if (!is.null(pent_prune_df) && nrow(pent_prune_df) > 0) {
        for (i in seq_len(nrow(pent_prune_df))) {
          idx_name <- as.character(pent_prune_df$index[i])
          pentad_j <- as.integer(pent_prune_df$pentad[i])
          idx_i <- which(index_names == idx_name)
          if (length(idx_i) == 1 && is.finite(pentad_j)) {
            pos <- (idx_i - 1) * n_bins + pentad_j
            if (pos >= 1 && pos <= length(weights)) weights[pos] <- 0
          }
        }
      }

      weights
    }

    unique_classes <- unique(oob_for_eval$target_class)
    if (length(unique_classes) < 2) stop("[ALPHA OPT] Need >=2 classes in OOB")

    endm <- build_endmember_matrix(mesma_lib_initial, unique_classes)
    oob_built <- build_Y_from_df(oob_for_eval, MESMA_PARAMS_INITIAL, TEMPORAL_BUDGET)
    Y_base <- oob_built$Y
    oob_labels <- oob_built$labels

    # Training evaluation set (exclude OOB locations)
    train_eval_df <- dplyr::mutate(
      dplyr::filter(multi_class_data, !location_id %in% unique(oob_for_eval$location_id)),
      target_class = tolower(as.character(Veg))
    )
    train_eval_df <- dplyr::filter(train_eval_df, !is.na(target_class) & target_class != "")
    train_built <- build_Y_from_df(train_eval_df, MESMA_PARAMS_INITIAL, TEMPORAL_BUDGET)
    Y_train_eval <- train_built$Y
    train_eval_labels <- train_built$labels

    # Alpha grid search
    index_alphas <- seq(0.01, 0.5, by = 0.01)
    pentad_alphas <- seq(0.01, 0.5, by = 0.01)

    alpha_grid_results <- data.frame(
      idx_alpha = numeric(), pent_alpha = numeric(),
      oob_score = numeric(), train_score = numeric(),
      combined_score = numeric(), n_features = integer(),
      n_indices = integer(),
      stringsAsFactors = FALSE
    )

    # Make sure missing p-values never cause accidental pruning.
    # Indices missing from index_results are treated as p=0 (always eligible to keep).
    p_by_index <- setNames(rep(0, length(perm_indices)), perm_indices)
    if (nrow(index_results) > 0) {
      tmp <- tapply(as.numeric(index_results$p_value), index_results$index, function(x) suppressWarnings(min(x, na.rm = TRUE)))
      tmp <- tmp[names(tmp) %in% names(p_by_index)]
      p_by_index[names(tmp)] <- tmp
    }

    for (ia in index_alphas) {
      idx_keep <- names(p_by_index)[is.finite(p_by_index) & p_by_index < ia]
      idx_prune <- setdiff(perm_indices, idx_keep)

      for (pa in pentad_alphas) {
        pent_prune <- if (nrow(pentad_results) > 0) {
          pentad_results[pentad_results$p_value >= pa, c("index", "pentad"), drop = FALSE]
        } else {
          data.frame(index = character(), pentad = integer(), stringsAsFactors = FALSE)
        }

        w <- apply_pruning_rules(MESMA_PARAMS_INITIAL$weights, perm_indices, TEMPORAL_BUDGET, idx_prune, pent_prune)
        n_feat <- sum(w > 0)
        n_idx <- sum(sapply(seq_along(perm_indices), function(k) {
          any(w[((k - 1) * TEMPORAL_BUDGET + 1):(k * TEMPORAL_BUDGET)] > 0)
        }))

        if (n_idx < 2) next

        oob_score <- compute_score_from_Y(endm$E, endm$col_names, Y_base, oob_labels, unique_classes, w)
        train_score <- compute_score_from_Y(endm$E, endm$col_names, Y_train_eval, train_eval_labels, unique_classes, w)
        combined <- min(oob_score, train_score)

        alpha_grid_results <- rbind(alpha_grid_results, data.frame(
          idx_alpha = ia, pent_alpha = pa,
          oob_score = oob_score, train_score = train_score,
          combined_score = combined, n_features = n_feat,
          n_indices = n_idx,
          stringsAsFactors = FALSE
        ))
      }
    }

    if (nrow(alpha_grid_results) == 0) stop("[ALPHA OPT] No alpha combinations evaluated")
    best_row <- which.max(alpha_grid_results$combined_score)
    best_idx_alpha <- alpha_grid_results$idx_alpha[best_row]
    best_pent_alpha <- alpha_grid_results$pent_alpha[best_row]

    # Apply optimized alphas
    PERMUTATION_ALPHA <- best_idx_alpha
    PERMUTATION_PENTAD_ALPHA <- best_pent_alpha

    idx_keep_final <- names(p_by_index)[is.finite(p_by_index) & p_by_index < PERMUTATION_ALPHA]
    idx_prune_final <- setdiff(perm_indices, idx_keep_final)
    pent_prune_final <- if (nrow(pentad_results) > 0) {
      pentad_results[pentad_results$p_value >= PERMUTATION_PENTAD_ALPHA, c("index", "pentad"), drop = FALSE]
    } else {
      data.frame(index = character(), pentad = integer(), stringsAsFactors = FALSE)
    }

    pruned_weights <- apply_pruning_rules(MESMA_PARAMS_INITIAL$weights, perm_indices, TEMPORAL_BUDGET, idx_prune_final, pent_prune_final)

    # Determine kept indices (any non-zero weight)
    indices_to_keep <- character(0)
    for (k in seq_along(perm_indices)) {
      s <- (k - 1) * TEMPORAL_BUDGET + 1
      e <- k * TEMPORAL_BUDGET
      if (any(pruned_weights[s:e] > 0)) indices_to_keep <- c(indices_to_keep, perm_indices[k])
    }
    indices_to_keep <- unique(indices_to_keep)

    if (length(indices_to_keep) >= PRUNE_ZERO_MIN_FEATURES) {
      pruned_indices <- indices_to_keep
      MESMA_PARAMS <- MESMA_PARAMS_INITIAL
      MESMA_PARAMS$weights <- pruned_weights
      avail <- pruned_indices
    } else {
      # Keep full feature set if pruning is too aggressive
      pruned_indices <- avail
      MESMA_PARAMS <- MESMA_PARAMS_INITIAL
    }

    write.csv(alpha_grid_results, file = "alpha_grid_results.csv", row.names = FALSE)
  }

  if (!is.null(MESMA_PARAMS) && !is.null(MESMA_PARAMS$weights)) {
    print_weights_summary("MESMA_FINAL", MESMA_PARAMS)
  }

  # Store pruning results
  assign("PRUNED_INDICES", pruned_indices, envir = globalenv())

  # NOTE: STEP 4 (confusion matrix / validation matrix) is computed later,
  # after the OOB cluster optimization that happens during final library build.

  # MESMA UNMIXING: All endmembers (barren + veg types) treated equally
  cat("[MODE] MESMA unmixing ENABLED (barren and vegetation types treated as equals)\n")

  # Log loss function choice
  if (exists("USE_HUBER_LOSS") && isTRUE(USE_HUBER_LOSS)) {
    cat(sprintf("[MODE] Using HUBER LOSS for FCLS (delta=%.3f, robust to outliers)\n",
                if (exists("HUBER_DELTA")) HUBER_DELTA else 1.345))
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
  # with feature weights. Observations: z-score -> prune to match library.
  compute_confusion_matrix <- function(df_data, label = "Training") {
    if (is.null(df_data) || nrow(df_data) == 0) return(invisible(NULL))
    if (!exists("OPTIMIZED_LIBRARY") || is.null(OPTIMIZED_LIBRARY)) return(invisible(NULL))
    if (is.null(MESMA_PARAMS) || is.null(MESMA_PARAMS$weights)) return(invisible(NULL))

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
      if (length(E_cols) < 2) return(invisible(NULL))

      E <- do.call(cbind, E_cols)

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
          if (!is.null(pruned_info) && !is.null(pruned_info$kept_idx)) {
            vec <- vec[pruned_info$kept_idx]
          }
          obs_vecs[[length(obs_vecs) + 1]] <- vec
          obs_labels <- c(obs_labels, true_cls)
        }
      }
      if (length(obs_vecs) == 0) return(invisible(NULL))

      Y <- do.call(rbind, obs_vecs)

      # 3. FCLS with pruned weights
      w <- MESMA_PARAMS$weights
      if (!is.null(pruned_info) && !is.null(pruned_info$kept_idx)) {
        w <- w[pruned_info$kept_idx]
      }
      all_coefs <- solve_batch_fcls(E, Y, w)
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
# TRAINING DISABLED
# This repository/script no longer performs model *training*. Training is expected
# to be performed offline and any required artifacts (normalization params,
# libraries, models) should be provided on disk. The script will continue to
# perform inference and validation where possible.
# ==========================================================================

# validation_location_ids already set during stratified train/validation split

# NOTE: do NOT create or populate `df_tasks` here — inference code will build
# task tables from the inference input. Training/task-generation code was
# removed; keep this script focused on inference + validation only.

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

# ============================================================================
# CLASSIFICATION UNCERTAINTY: DIRICHLET PERTURBATION
# Propagates classification/misclassification uncertainty through bootstrap by
# resampling fractional covers from a Dirichlet distribution whose mean is
# derived from the confusion matrix.
#
# IMPORTANT: This implementation is mixture-aware: it does NOT condition on a
# single dominant/true class. Instead, it uses the current fraction vector as
# weights over confusion-matrix rows.
# ============================================================================

# Global storage for confusion matrix and validation sample sizes
.CONFUSION_MATRIX <- NULL
.VALIDATION_SAMPLE_SIZES <- NULL

# Store the confusion matrix for use in Dirichlet perturbation
# Called after validation accuracy is computed
store_confusion_matrix <- function(conf_matrix, sample_sizes = NULL) {
  if (!is.matrix(conf_matrix) || nrow(conf_matrix) != ncol(conf_matrix)) {
    stop("[DIRICHLET] Invalid confusion matrix - must be square matrix")
  }

  # Ensure row-normalized (each row sums to 1)
  row_sums <- rowSums(conf_matrix, na.rm = TRUE)
  conf_matrix_norm <- conf_matrix
  for (i in seq_len(nrow(conf_matrix_norm))) {
    if (row_sums[i] > 1e-9) {
      conf_matrix_norm[i, ] <- conf_matrix_norm[i, ] / row_sums[i]
    }
  }

  assign(".CONFUSION_MATRIX", conf_matrix_norm, envir = globalenv())

  if (!is.null(sample_sizes)) {
    assign(".VALIDATION_SAMPLE_SIZES", sample_sizes, envir = globalenv())
  }

  if (isTRUE(DEBUG_UNCERTAINTY)) {
    cat(sprintf("[DIRICHLET] Stored %dx%d confusion matrix for classification uncertainty\n",
                nrow(conf_matrix_norm), ncol(conf_matrix_norm)))
    cat(sprintf("[DIRICHLET] Classes: %s\n", paste(rownames(conf_matrix_norm), collapse = ", ")))
    if (!is.null(sample_sizes)) {
      cat(sprintf("[DIRICHLET] Validation sample sizes: %s\n",
                  paste(sprintf("%s=%d", names(sample_sizes), sample_sizes), collapse = ", ")))
    }
  }

  return(invisible(TRUE))
}

# Apply Dirichlet perturbation to fractional covers based on true class
# This simulates classification uncertainty by drawing from a Dirichlet
# distribution centered on the confusion matrix row for the true class
#
# Parameters:
#   fractions: Named numeric vector of fractional covers (must sum to ~1)
#   true_class: The true vegetation class for this location
#   conf_matrix: Row-normalized confusion matrix (optional, uses global if NULL)
#   sample_sizes: Named vector of validation sample sizes per class (optional)
#   concentration_scale: Multiplier for sample size to get Dirichlet alpha
#
# Returns:
#   Perturbed fractional covers (same names, still sum to 1)
#
apply_dirichlet_perturbation <- function(fractions, true_class,
                                          conf_matrix = NULL,
                                          sample_sizes = NULL,
                                          concentration_scale = DIRICHLET_CONCENTRATION_SCALE) {
  # Use global confusion matrix if not provided
 if (is.null(conf_matrix)) {
    if (exists(".CONFUSION_MATRIX", envir = globalenv())) {
      conf_matrix <- get(".CONFUSION_MATRIX", envir = globalenv())
    } else {
      stop("[DIRICHLET] No confusion matrix available")
    }
  }

  if (is.null(sample_sizes) && exists(".VALIDATION_SAMPLE_SIZES", envir = globalenv())) {
    sample_sizes <- get(".VALIDATION_SAMPLE_SIZES", envir = globalenv())
  }

  # Validate inputs
  if (is.null(conf_matrix) || !is.matrix(conf_matrix)) {
    stop("[DIRICHLET] Invalid confusion matrix")
  }

  # Match true_class to confusion matrix row (case-insensitive)
  row_names <- tolower(rownames(conf_matrix))
  true_class_lower <- tolower(true_class)
  row_idx <- which(row_names == true_class_lower)

  if (length(row_idx) == 0) {
    stop(sprintf("[DIRICHLET] True class '%s' not found in confusion matrix", as.character(true_class)))
  }

  # Get the confusion row for the true class (this is our Dirichlet mean)
  mu <- as.numeric(conf_matrix[row_idx, ])
  class_names <- colnames(conf_matrix)

  # Ensure mu sums to 1 and has no zeros (add small epsilon for stability)
  mu[mu < 1e-6] <- 1e-6
  mu <- mu / sum(mu)

  # Determine concentration parameter alpha
  # Higher alpha = tighter concentration around mean = less noise
  if (!is.null(sample_sizes) && true_class_lower %in% tolower(names(sample_sizes))) {
    n_val <- sample_sizes[tolower(names(sample_sizes)) == true_class_lower]
    alpha <- concentration_scale * n_val
  } else {
    # Default: assume moderate validation sample size
    alpha <- concentration_scale * 20
  }

  # Ensure minimum alpha for numerical stability
  alpha <- max(alpha, 5)

  # Dirichlet parameter vector
  alpha_vec <- alpha * mu

  # Draw from Dirichlet distribution
  # Using the gamma distribution method: X_i ~ Gamma(alpha_i, 1), then normalize
  gamma_draws <- rgamma(length(alpha_vec), shape = alpha_vec, rate = 1)

  # Handle edge case where all draws are zero
  if (sum(gamma_draws) < 1e-10) {
    gamma_draws <- mu  # Fall back to mean
  }

  perturbed <- gamma_draws / sum(gamma_draws)

  # Map back to original fraction names
  # The confusion matrix classes may not exactly match the fraction names
  result <- fractions
  frac_names_lower <- tolower(names(fractions))

  for (i in seq_along(class_names)) {
    class_lower <- tolower(class_names[i])
    # Find matching fraction column (handle "frac_" prefix)
    match_idx <- which(frac_names_lower == class_lower |
                       frac_names_lower == paste0("frac_", class_lower))
    if (length(match_idx) > 0) {
      result[match_idx[1]] <- perturbed[i]
    }
  }

  # Re-normalize to ensure sum = 1 (in case of partial matching)
  if (sum(result) > 0) {
    result <- result / sum(result)
  }

  return(result)
}

# Apply mixture-aware Dirichlet perturbation to fractional covers.
#
# Instead of conditioning on a single "true" class, this builds a Dirichlet mean
# as a weighted mixture of confusion-matrix rows using the current fraction
# vector as weights.
apply_dirichlet_perturbation_mixture <- function(fractions,
                                                 conf_matrix = NULL,
                                                 sample_sizes = NULL,
                                                 concentration_scale = DIRICHLET_CONCENTRATION_SCALE,
                                                 min_weight = 1e-6) {
  if (is.null(conf_matrix)) {
    if (exists(".CONFUSION_MATRIX", envir = globalenv())) {
      conf_matrix <- get(".CONFUSION_MATRIX", envir = globalenv())
    } else {
      stop("[DIRICHLET] No confusion matrix available")
    }
  }

  if (is.null(sample_sizes) && exists(".VALIDATION_SAMPLE_SIZES", envir = globalenv())) {
    sample_sizes <- get(".VALIDATION_SAMPLE_SIZES", envir = globalenv())
  }

  if (is.null(conf_matrix) || !is.matrix(conf_matrix) || nrow(conf_matrix) != ncol(conf_matrix)) {
    stop("[DIRICHLET] Invalid confusion matrix")
  }

  if (is.null(fractions) || length(fractions) == 0) {
    return(fractions)
  }

  fracs <- as.numeric(fractions)
  names(fracs) <- names(fractions)
  fracs[!is.finite(fracs)] <- 0
  fracs[fracs < 0] <- 0

  frac_sum <- sum(fracs)
  if (!is.finite(frac_sum) || frac_sum <= 0) {
    return(fractions)
  }
  fracs <- fracs / frac_sum

  row_names <- tolower(rownames(conf_matrix))
  col_names <- colnames(conf_matrix)
  frac_names_lower <- tolower(names(fracs))

  # Map fraction names to confusion-matrix rows (support 'frac_' prefix)
  mapped_row <- rep(NA_character_, length(fracs))
  for (i in seq_along(fracs)) {
    nm <- frac_names_lower[i]
    nm2 <- sub("^frac_", "", nm)
    if (nm %in% row_names) {
      mapped_row[i] <- nm
    } else if (nm2 %in% row_names) {
      mapped_row[i] <- nm2
    }
  }

  keep <- is.finite(fracs) & fracs > 0 & !is.na(mapped_row) & (fracs >= min_weight)
  if (!any(keep)) {
    # No overlap with confusion matrix; leave unchanged
    return(fracs)
  }

  w <- fracs[keep]
  w <- w / sum(w)
  rows <- mapped_row[keep]
  row_idx <- match(rows, row_names)

  # Weighted mean confusion row
  mu <- rep(0, ncol(conf_matrix))
  for (k in seq_along(row_idx)) {
    mu <- mu + w[k] * as.numeric(conf_matrix[row_idx[k], ])
  }

  mu[!is.finite(mu)] <- 0
  mu[mu < 1e-6] <- 1e-6
  mu <- mu / sum(mu)

  # Concentration based on validation support (fraction-weighted)
  n_eff <- 20
  if (!is.null(sample_sizes) && length(sample_sizes) > 0) {
    ss_names <- tolower(names(sample_sizes))
    ss <- as.numeric(sample_sizes)
    ss[!is.finite(ss) | ss < 0] <- 0
    ss_match <- ss[match(rows, ss_names)]
    ss_match[is.na(ss_match)] <- 0
    if (sum(ss_match) > 0) {
      n_eff <- sum(w * ss_match)
    }
  }

  alpha <- max(concentration_scale * n_eff, 5)
  alpha_vec <- alpha * mu
  gamma_draws <- rgamma(length(alpha_vec), shape = alpha_vec, rate = 1)
  if (sum(gamma_draws) < 1e-10) {
    gamma_draws <- mu
  }
  perturbed <- gamma_draws / sum(gamma_draws)

  # Map back to original fraction names
  result <- fracs
  for (i in seq_along(col_names)) {
    class_lower <- tolower(col_names[i])
    match_idx <- which(frac_names_lower == class_lower |
                         frac_names_lower == paste0("frac_", class_lower))
    if (length(match_idx) > 0) {
      result[match_idx[1]] <- perturbed[i]
    }
  }

  if (sum(result) > 0) {
    result <- result / sum(result)
  }

  return(result)
}

# Batch apply mixture-aware Dirichlet perturbation to a data frame of fractions.
# Useful for resampling location-level results without requiring class labels.
apply_dirichlet_perturbation_batch <- function(coef_df, frac_cols, true_class_col = NULL) {
  if (!isTRUE(ENABLE_CLASSIFICATION_UNCERTAINTY)) {
    return(coef_df)
  }

  # Check if confusion matrix is available
  if (!exists(".CONFUSION_MATRIX", envir = globalenv())) {
    stop("[DIRICHLET] No confusion matrix available")
  }

  result_df <- coef_df

  for (i in seq_len(nrow(coef_df))) {
    # Extract current fractions
    fracs <- as.numeric(coef_df[i, frac_cols])
    names(fracs) <- frac_cols

    # Apply mixture-aware perturbation
    perturbed <- apply_dirichlet_perturbation_mixture(fracs)

    # Store back
    result_df[i, frac_cols] <- perturbed
  }

  return(result_df)
}

# Simple i.i.d. residual bootstrap
# Resamples residuals with replacement (assumes approximate independence)
simple_residual_bootstrap <- function(residuals) {
  n <- length(residuals)
  if (n == 0) return(numeric(0))
  
  # Remove NA values for bootstrapping
  residuals_non_na <- residuals[!is.na(residuals)]
  n_non_na <- length(residuals_non_na)
  
  if (n_non_na == 0) return(residuals) # Return original if all are NA
  
  # Simple random sampling with replacement
  resampled_indices <- sample(1:n_non_na, n_non_na, replace = TRUE)
  
  # Create the bootstrapped residuals vector
  bootstrapped_residuals_non_na <- residuals_non_na[resampled_indices]
  
  # Place bootstrapped residuals back into the original structure with NAs
  residuals_boot <- residuals
  residuals_boot[!is.na(residuals)] <- bootstrapped_residuals_non_na
  
  return(residuals_boot)
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
    
    # For inference/validation data: interpolation is disabled by default to avoid filling gaps in held-out data.
    # It can be enabled by setting `INTERPOLATE_INFERENCE <- TRUE` or `PARAMS$interpolate_inference <- TRUE`.
    base_indices_for_build <- if (!is.null(PARAMS$base_indices)) PARAMS$base_indices else PARAMS$indices
    interpolate_for_inference <- if (!is.null(PARAMS$interpolate_inference)) as.logical(PARAMS$interpolate_inference) else INTERPOLATE_INFERENCE
    if (isTRUE(interpolate_for_inference) && !isTRUE(QUIET_MODE)) {
      cat("[NOTICE] Interpolating missing pentads for inference/validation (INTERPOLATE_INFERENCE=TRUE)\n")
    }
    raw_mat <- build_pentad_matrix(task_data, base_indices_for_build, interpolate = isTRUE(interpolate_for_inference))
    if (is.null(raw_mat)) {
      if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG] loc=%s yr=%d: build_pentad_matrix returned NULL\n", loc, yr))
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

    # Previously we required a minimum fraction of valid observations and would skip tasks
    # with low coverage. Remove that strict filter: proceed whenever there is at least one
    # valid observation, but emit a warning in TESTING_MODE when coverage is very low.
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

      # Apply PCA-LDA transform to the z-scored observation so inference sits in the same
      # PCA-LDA-scored feature space as training. This is applied before masking so column
      # indices remain aligned between obs and library templates.
      if (!is.null(PARAMS$weights) && length(PARAMS$weights) == length(y_norm)) {
        y_norm_full_pca_lda <- y_norm  # PCA-LDA weights applied in solver via feature_weights
      } else {
        y_norm_full_pca_lda <- y_norm
      }

      # If the library was pruned, apply the same column subsetting to the
      # observation vector, validity mask and weights so dimensions match.
      weights_for_mask <- PARAMS$weights
      if (!is.null(.pruned_kept_idx) && length(.pruned_kept_idx) < length(y_norm_full_pca_lda)) {
        y_norm_full_pca_lda <- y_norm_full_pca_lda[.pruned_kept_idx]
        valid_mask <- valid_mask[.pruned_kept_idx]
        if (!is.null(weights_for_mask) && length(weights_for_mask) > length(.pruned_kept_idx)) {
          weights_for_mask <- weights_for_mask[.pruned_kept_idx]
        }
      }

      # Mask the observation (z-scored only, no PCA-LDA weighting applied to data)
      y_norm_masked <- y_norm_full_pca_lda[valid_mask]

      # Extract masked PCA-LDA weights to pass to solver
      # The solver will apply weights directly to both endmembers and observations
      if (!is.null(weights_for_mask) && length(weights_for_mask) == length(y_norm_full_pca_lda)) {
        weights_masked <- weights_for_mask[valid_mask]
      } else {
        weights_masked <- rep(1, length(y_norm_masked))
      }

      # Check if we have sufficient signal (using unweighted norm for validation)
      y_norm_val <- sqrt(sum(y_norm_masked^2, na.rm = TRUE))

      if (is.na(y_norm_val) || y_norm_val < 1e-9) {
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_task] loc=%s yr=%d: insufficient signal (norm=%.6g), skipping\n", loc, yr, y_norm_val))
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
      
      # DEBUG: Check library availability
      if (is.null(OPTIMIZED_LIBRARY) || length(OPTIMIZED_LIBRARY) == 0) {
        cat(sprintf("[ERROR] loc=%s yr=%d: OPTIMIZED_LIBRARY is NULL or empty\n", loc, yr))
        return(NULL)
      }
      if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG LIBRARY] loc=%s yr=%d: OPTIMIZED_LIBRARY contains %d vegetation types: %s\n",
                  loc, yr, length(OPTIMIZED_LIBRARY), paste(names(OPTIMIZED_LIBRARY), collapse=", ")))
      
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

        # Apply PCA-LDA transform to each variant row (if we have matching weights)
        if (!is.null(PARAMS$weights) && length(PARAMS$weights) == ncol(M_full)) {
          M_full_trans <- t(apply(M_full, 1, function(r) {
            as.numeric(r)  # PCA-LDA weights applied in solver via feature_weights
          }))
        } else {
          M_full_trans <- M_full
        }

        # Mask to the valid observation columns
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
          if (isTRUE(TESTING_MODE) && length(masked_vec) != length(y_for_unmixing)) {
            cat(sprintf("[WARN fit_one_task] loc=%s yr=%d veg=%s: template length mismatch (template=%d, y=%d)\n",
                        loc, yr, v, length(masked_vec), length(y_for_unmixing)))
          }
          list(vec = masked_vec, id = lib_ids_kept[i], similarity = sims[i])
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
      
      best_result <- tryCatch({
        evaluate_all_combinations(
          y_for_unmixing,
          top_variants,
          lambda = 0,
          feature_weights = weights_masked
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
      cat(sprintf("[MESMA OUTPUT] loc=%s yr=%d: %s (rmse=%.4f)\n",
                  loc, yr,
                  paste(sprintf("%s=%.4f", names(coefs), coefs), collapse=", "),
                  rmse))

      
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

      # --- Rename high-similarity populus/tamarix variants to "woody_unknown" ---
      # Variants with cross-class similarity > 0.95 between populus and tamarix are indistinguishable
      if (exists("WOODY_UNKNOWN_VARIANTS", envir = globalenv())) {
        woody_unknown_list <- get("WOODY_UNKNOWN_VARIANTS", envir = globalenv())
        if (length(woody_unknown_list) > 0) {
          rename_mask <- coef_df$variant_id %in% woody_unknown_list
          rename_mask[is.na(rename_mask)] <- FALSE
          if (any(rename_mask)) {
            n_renamed <- sum(rename_mask)
            old_vegs <- coef_df$Veg[rename_mask]
            coef_df$Veg[rename_mask] <- "woody_unknown"
          }
        }
      }

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

      # --- Scale vegetation fractions by index-derived total vegetation cover ---
      # MESMA provides relative proportions among vegetation types.
      # We multiply the vegetation fractions by total_veg_cover (from MSAVI) to get absolute fractions.
      # Barren fraction from MESMA is preserved - we do NOT replace it with index-derived estimate.
      # NO equal distribution fallback - if MESMA says 0 vegetation, we keep it as 0.

      veg_coefs_mask <- tolower(coef_agg$Veg) != "barren"
      veg_coefs_mask[is.na(veg_coefs_mask)] <- FALSE
      sum_original_veg_coefs <- sum(coef_agg$coef[veg_coefs_mask], na.rm = TRUE)

      # Scale vegetation fractions: multiply by total_veg_cover / sum_of_veg_fractions
      # This preserves MESMA's relative proportions while scaling to index-derived total veg cover
      scale_factor <- 1.0
      if (is.finite(sum_original_veg_coefs) && sum_original_veg_coefs > 1e-9) {
        scale_factor <- total_veg_cover / sum_original_veg_coefs
        coef_agg$coef[veg_coefs_mask] <- coef_agg$coef[veg_coefs_mask] * scale_factor
      } else {
        # MESMA returned 0 vegetation - keep it as 0, don't fabricate equal fractions
      }

      # Keep MESMA's barren fraction - do NOT replace with index-derived estimate
      # Just ensure barren row exists in coef_agg
      barren_idx <- which(tolower(coef_agg$Veg) == "barren")
      if (length(barren_idx) == 0) {
        # Add barren row from MESMA if it doesn't exist (use 1 - sum of scaled veg fractions)
        mesma_barren <- 1 - sum(coef_agg$coef[veg_coefs_mask], na.rm = TRUE)
        mesma_barren <- pmax(0, pmin(1, mesma_barren))
        barren_row <- data.frame(Veg = "barren", coef = mesma_barren)
        coef_agg <- rbind(coef_agg, barren_row)
      }
      # Note: If barren already exists in coef_agg, we keep MESMA's original value

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
        # Detailed diagnostics to help track down TRUE/FALSE NA issues
        if (exists("coef_df")) {
        }
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
      diag_df$barren_fraction_ppi_based <- barren_fraction # Add PPI-based barren fraction

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

# New function: Process all years for a single location with multi-year bootstrap
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
            if (isTRUE(TESTING_MODE)) cat(sprintf("[ERROR assigning year_results] loc=%s yr=%d: %s\n", loc, yr, e$message))
            stop(e)  # Re-throw to be caught by outer handler
          })

          # Store data needed for multi-year bootstrap - wrapped in tryCatch to isolate errors
          tryCatch({
            # For inference/validation data, do NOT interpolate - use only actual observations.
            # IMPORTANT: build from the same raw feature columns (base_indices) as fit_one_task(),
            # then apply representation + z-scoring to match the library feature space.
            n_bins <- TEMPORAL_BUDGET
            base_indices <- if (!is.null(MESMA_PARAMS$base_indices)) MESMA_PARAMS$base_indices else MESMA_PARAMS$indices
            n_base_idx <- length(base_indices)
            l2_normalize <- isTRUE(MESMA_PARAMS$l2_normalize)

            raw_mat_yr <- build_pentad_matrix(year_data, base_indices, interpolate = FALSE)
            if (!is.null(raw_mat_yr)) {
              vec_raw <- as.numeric(raw_mat_yr)
              y_work <- mesma_apply_representation_vec(vec_raw, n_base_idx, n_bins, l2_normalize)
              y_norm <- mesma_zscore_vec_by_index(y_work, MESMA_PARAMS$indices, MESMA_PARAMS$means, MESMA_PARAMS$sds, n_bins)

              expected_length <- length(MESMA_PARAMS$indices) * n_bins
              if (length(y_norm) != expected_length) {
                if (isTRUE(TESTING_MODE)) {
                  cat(sprintf("[WARN fit_one_location] loc=%s yr=%d: y_norm length mismatch (got %d, expected %d)\n"))
                }
              } else {
                y_norm[!is.finite(y_norm)] <- 0
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
            if (isTRUE(TESTING_MODE)) {
              cat(sprintf("[WARN fit_one_location] loc=%s yr=%d: multi-year processing failed (%s), continuing with single-year results\n",
                loc, yr, e$message))
            }
          })
        } else {
          if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_location] loc=%s yr=%d: fit_one_task returned NULL\n", loc, yr))
        }
      }

      # If we have results for multiple years and ENABLE_UNCERTAINTY, do multi-year bootstrap
      if (length(year_results) == 0) {
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_location] loc=%s: no successful year results\n", loc))
        return(NULL)
      }

      # Run per-year bootstrap uncertainty as long as we have at least one successful year.
      # (Historically this was gated on successfully reconstructing a full feature vector for every year,
      # but the bootstrap below only depends on each year's fitted residuals/endmembers.)
      if (isTRUE(ENABLE_UNCERTAINTY) && isTRUE(ENABLE_MULTI_YEAR_BOOTSTRAP) && length(year_results) >= 1) {
        # Get comp_templates from first successful year result
        top_variants <- NULL

        for (yr in years) {
          res_yr <- year_results[[as.character(yr)]]
          if (!is.null(res_yr) && !is.null(res_yr$chosen_variants)) {
            # Build comp_templates from the MESMA library
            vegs <- names(res_yr$chosen_variants)
            if (is.null(top_variants)) {
              top_variants <- list()
              for (v in vegs) top_variants[[v]] <- list()
            }

            break
          }
        }

        # Build comp_templates from mesma_lib
        if (!is.null(top_variants) && exists("mesma_lib")) {
          comp_templates <- list()
          for (v in names(top_variants)) {
            if (v %in% names(mesma_lib)) {
              comp_templates[[v]] <- mesma_lib[[v]]
            }
          }
        }

        # Perform per-year bootstrap to estimate uncertainty (coef CI, sd, rmse CI, variant frequencies)
        if (isTRUE(ENABLE_UNCERTAINTY)) {
          B_loc <- if (exists("BOOTSTRAP_B")) min(BOOTSTRAP_B, 100L) else 100L
          if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG] Performing local bootstrap (B=%d) for loc=%s\n", B_loc, loc))

          for (yr in years) {
            res_yr <- year_results[[as.character(yr)]]
            if (is.null(res_yr) || is.null(res_yr$residuals)) next

            year_data <- location_data[location_data$pheno_year == yr, , drop = FALSE]
            n_obs <- nrow(year_data)
            if (n_obs < max(6, MIN_OBS_PER_LOC_YEAR)) {
              if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG] Too few observations for bootstrap loc=%s yr=%s (n_obs=%d); skipping\n", loc, yr, n_obs))
              next
            }

            boot_coef_by_veg <- list()
            boot_variant_counts <- list()
            boot_rmse <- numeric(0)

            # Suppress verbose output during bootstraps
            suppress_output_safely({
              for (b in seq_len(B_loc)) {
                residuals_boot <- simple_residual_bootstrap(res_yr$residuals)
                y_boot <- res_yr$y_hat + residuals_boot

                E_boot <- res_yr$E_best
                if (is.null(E_boot) || !is.matrix(E_boot) || ncol(E_boot) < 1) next

                res_b <- tryCatch({
                  solve_weights_fcls(E_boot, y_boot, feature_weights = res_yr$weights_masked)
                }, error = function(e) NULL)

                if (is.null(res_b)) next

                coefs_b <- res_b$w
                chosen_ids_b <- res_yr$coef_df$variant_id
                names(chosen_ids_b) <- res_yr$coef_df$Veg
                names(coefs_b) <- res_yr$coef_df$Veg

                dfb <- data.frame(
                    Veg = sapply(strsplit(chosen_ids_b, "_(v|opt)"), `[`, 1),
                    coef = coefs_b,
                    stringsAsFactors = FALSE
                )

                for (i in seq_len(nrow(dfb))) {
                  vg <- as.character(dfb$Veg[i])
                  if (is.na(vg)) next
                  val <- as.numeric(dfb$coef[i])
                  if (!is.finite(val)) next
                  if (is.null(boot_coef_by_veg[[vg]])) boot_coef_by_veg[[vg]] <- numeric(0)
                  boot_coef_by_veg[[vg]] <- c(boot_coef_by_veg[[vg]], val)
                }

                for (i in seq_along(chosen_ids_b)) {
                  vg <- names(chosen_ids_b)[i]
                  if (is.na(vg) || !is.character(vg)) next
                  vid <- chosen_ids_b[[i]]
                  if (is.na(vid)) next
                  if (is.null(boot_variant_counts[[vg]])) boot_variant_counts[[vg]] <- list()
                  boot_variant_counts[[vg]][[as.character(vid)]] <- (boot_variant_counts[[vg]][[as.character(vid)]] %||% 0) + 1L
                }

                if (!is.null(res_b$rmse) && is.finite(as.numeric(res_b$rmse))) boot_rmse <- c(boot_rmse, as.numeric(res_b$rmse))
              }
            })
            
            # Build coef CI table
            coef_ci_df <- NULL
            if (length(boot_coef_by_veg) > 0) {
              rows <- lapply(names(boot_coef_by_veg), function(vg) {
                vals <- boot_coef_by_veg[[vg]]
                vals <- vals[is.finite(vals)]
                if (length(vals) == 0) return(NULL)
                data.frame(Veg = vg,
                           coef_median = as.numeric(median(vals, na.rm = TRUE)),
                           coef_025 = as.numeric(quantile(vals, 0.025, na.rm = TRUE)),
                           coef_975 = as.numeric(quantile(vals, 0.975, na.rm = TRUE)),
                           coef_sd = as.numeric(stats::sd(vals, na.rm = TRUE)),
                           interval = as.numeric(quantile(vals, 0.975, na.rm = TRUE) - quantile(vals, 0.025, na.rm = TRUE)),
                           stringsAsFactors = FALSE)
              })
              rows <- rows[!sapply(rows, is.null)]
              if (length(rows) > 0) coef_ci_df <- do.call(rbind, rows)
            }

            # Build variant frequency table
            variant_freq_df <- NULL
            if (length(boot_variant_counts) > 0) {
              rows2 <- lapply(names(boot_variant_counts), function(vg) {
                tbl <- boot_variant_counts[[vg]]
                keys <- names(tbl)
                counts <- as.integer(unlist(tbl))
                dfv <- data.frame(Veg = rep(vg, length(keys)), Variant = keys, N = counts, Percent = counts / sum(counts) * 100, stringsAsFactors = FALSE)
                dfv
              })
              rows2 <- rows2[!sapply(rows2, is.null)]
              if (length(rows2) > 0) variant_freq_df <- do.call(rbind, rows2)
            }

            rmse_ci <- NULL
            if (length(boot_rmse) > 0) {
              rmse_ci <- as.numeric(quantile(boot_rmse, c(0.025, 0.975), na.rm = TRUE))
            }

            # Attach uncertainty to year_results
            if (!is.null(coef_ci_df) || !is.null(variant_freq_df) || !is.null(rmse_ci)) {
              if (is.null(year_results[[as.character(yr)]]$uncertainty)) year_results[[as.character(yr)]]$uncertainty <- list()
              year_results[[as.character(yr)]]$uncertainty$coef_ci <- coef_ci_df
              year_results[[as.character(yr)]]$uncertainty$variant_freq <- variant_freq_df
              year_results[[as.character(yr)]]$uncertainty$rmse_ci <- rmse_ci
              
              # Update coef_df with uncertainty metrics if available
              if (!is.null(coef_ci_df)) {
                  curr_df <- year_results[[as.character(yr)]]$coef_df
                  if (!is.null(curr_df)) {
                      # We loop over rows of coef_ci_df
                      for(r_idx in seq_len(nrow(coef_ci_df))) {
                          v_type <- coef_ci_df$Veg[r_idx]
                          # Find matching row(s) in curr_df (should be by Veg)
                          match_idx <- which(curr_df$Veg == v_type)
                          if(length(match_idx) > 0) {
                              curr_df$coef_sd[match_idx] <- coef_ci_df$coef_sd[r_idx]
                              curr_df$coef_025[match_idx] <- coef_ci_df$coef_025[r_idx]
                              curr_df$coef_975[match_idx] <- coef_ci_df$coef_975[r_idx]
                              curr_df$interval[match_idx] <- coef_ci_df$interval[r_idx]
                              # Store median bootstrap estimate for plotting/prediction preference
                              if ("coef_median" %in% names(coef_ci_df)) {
                                curr_df$coef_median[match_idx] <- coef_ci_df$coef_median[r_idx]
                              }
                            }
                        }
                        year_results[[as.character(yr)]]$coef_df <- curr_df
                    }
                }
              }
            }
          }

        if (!is.null(top_variants) && exists("mesma_lib")) {
          comp_templates <- list()
          for (v in names(top_variants)) {
            if (v %in% names(mesma_lib)) {
              comp_templates[[v]] <- mesma_lib[[v]]
            }
          }
        }
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

# Executing main processing steps directly (function removed). This section previously defined
# `main_processing_block <- function() { ... }`. To avoid large function compile stalls, the body
# is now executed at top-level during script run. The original nested helper functions remain in scope.

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
    cat("Preparing locations for batched processing (multi-year bootstrap)...\n")

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

  if (!isTRUE(QUIET_MODE)) cat(sprintf("\r  [Batch %d/%d complete]  ", i, n_batches))
  }
  cat("\n")
  
  end_time <- Sys.time()
  processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(
    "\nTraining processing finished in %.2f seconds (%.2f minutes)\n",
    processing_time, processing_time / 60
  ))
  } # End of n_locs_to_process > 0 conditional

  # Debug log output removed (strict mode)

  # ==========================================================================
  # VALIDATION PROCESSING: unmix held-out validation locations
  # ==========================================================================
  cat("\n=== STARTING VALIDATION PROCESSING ===\n")

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
                                      include_species_plots = TRUE,
                                      include_woody_types_plot = TRUE)
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

    # --- INFERENCE: No-index normalization aggregation ---
    if (!is.null(inference_coefs) && nrow(inference_coefs) > 0) {
      cat("[INFERENCE] Running no-index normalization aggregation (inference mode)\n")
      noindex_inf_full <- tryCatch({ location_bootstrap_noindex(inference_coefs, df_tasks_inference_proc, B = BOOTSTRAP_B, seed = get_mesma_seed(123)) }, error = function(e) { cat(sprintf("[NOINDEX INFERENCE] failed: %s\n", e$message)); NULL })
      if (!is.null(noindex_inf_full) && nrow(noindex_inf_full) > 0) {
        plot_inference_method_results(noindex_inf_full, "No-Index", "noindex",
                                      use_excluded_years_shade = FALSE,
                                      include_species_plots = FALSE)
      } else {
        cat("[INFERENCE] No-index inference aggregation returned no results.\n")
      }
    } else {
      cat("[INFERENCE] No inference coefficients available; skipping no-index aggregation.\n")
    }

    # --- INFERENCE: Dirichlet-perturbed aggregate bootstrap ---
    if (!is.null(inference_coefs) && nrow(inference_coefs) > 0) {
      cat("[INFERENCE] Running Dirichlet-perturbed aggregate bootstrap\n")
      aggregate_inf_full <- tryCatch({
        location_bootstrap_aggregate(inference_coefs, B = BOOTSTRAP_B, seed = get_mesma_seed(123))
      }, error = function(e) { cat(sprintf("[AGGREGATE BOOTSTRAP] failed: %s\n", e$message)); NULL })

      if (!is.null(aggregate_inf_full) && nrow(aggregate_inf_full) > 0) {
        agg_veg_norm <- normalize_veg_name(aggregate_inf_full$Veg)
        agg_inf_veg <- aggregate_inf_full[!agg_veg_norm %in% c("barren"), ]
        agg_inf_barren <- aggregate_inf_full[agg_veg_norm %in% c("barren"), ]

        if (!is.null(agg_inf_veg) && nrow(agg_inf_veg) > 0) {
          p_inf_agg_ts <- ggplot(agg_inf_veg, aes(x = year, y = global_coef, color = Veg, group = Veg)) +
            add_year_lines(is_date = FALSE) +
            geom_line(linewidth = 1) +
            geom_point(show.legend = FALSE) +
            geom_ribbon(aes(ymin = coef_025, ymax = coef_975, fill = Veg), alpha = 0.15, color = NA) +
            labs(title = "Dirichlet-Perturbed Aggregate Vegetation Fractions",
                 x = "Year", y = "Mean Fraction", color = "Veg", fill = "Veg") +
            theme_minimal()
          ggsave(file.path(OUT_DIR, "inference_aggregate_dirichlet_timeseries.png"), p_inf_agg_ts, width = 8, height = 6)
          readr::write_csv(agg_inf_veg, file.path(OUT_DIR, "inference_aggregate_dirichlet_timeseries.csv"))
          cat(sprintf("Saved Dirichlet-perturbed aggregate time series plot to: %s\n", file.path(OUT_DIR, "inference_aggregate_dirichlet_timeseries.png")))

          # Herbs vs Woody
          herbs_woody_agg <- aggregate_inf_full[agg_veg_norm %in% c("herbs", "woody"), ]
          if (nrow(herbs_woody_agg) > 0) {
            herbs_woody_agg$Veg <- ifelse(tolower(herbs_woody_agg$Veg) == "herbs", "Herbs", "Woody")
            p_inf_agg_hw <- ggplot(herbs_woody_agg, aes(x = year, y = global_coef, color = Veg, group = Veg)) +
              add_year_lines(is_date = FALSE) +
              geom_line(linewidth = 1) +
              geom_point(show.legend = FALSE) +
              geom_ribbon(aes(ymin = coef_025, ymax = coef_975, fill = Veg), alpha = 0.15, color = NA) +
              scale_color_manual(values = HERBS_WOODY_COLORS) +
              scale_fill_manual(values = HERBS_WOODY_COLORS) +
              labs(title = "Dirichlet-Perturbed: Herbs vs Woody",
                   x = "Year", y = "Mean Fraction", color = "Type", fill = "Type") +
              theme_minimal()
            ggsave(file.path(OUT_DIR, "inference_aggregate_dirichlet_herbs_vs_woody.png"), p_inf_agg_hw, width = 8, height = 6)
            readr::write_csv(herbs_woody_agg, file.path(OUT_DIR, "inference_aggregate_dirichlet_herbs_vs_woody.csv"))
            cat(sprintf("Saved Dirichlet-perturbed herbs vs woody plot to: %s\n", file.path(OUT_DIR, "inference_aggregate_dirichlet_herbs_vs_woody.png")))
          }
        }
      } else {
        cat("[INFERENCE] Dirichlet-perturbed aggregate bootstrap returned no results.\n")
      }
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

  # Training results CSV saving is disabled by default (user requested removal of training outputs)
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

  # Helper to report validation accuracy (including artificial mixes) to avoid duplication
  report_validation_accuracy <- function(val_coefs, label = "") {
    if (is.null(val_coefs) || nrow(val_coefs) == 0) {
      cat("[NOTICE] No validation coefficients found (validation locations not in results).\n")
      return(invisible(NULL))
    }

    prefix <- if (nzchar(label)) paste0(" (", label, ")") else ""

    # Optional PPI-based truth for vegetation cover (derived from PPI barren fraction)
    ppi_truth <- NULL
    if (exists("all_diagnostics") && !is.null(all_diagnostics) &&
        all(c("location_id", "pheno_year") %in% names(all_diagnostics)) &&
        "barren_fraction_ppi_based" %in% names(all_diagnostics)) {
      ppi_truth <- all_diagnostics |>
        dplyr::select(location_id, pheno_year, barren_fraction_ppi_based) |>
        dplyr::mutate(
          barren_fraction_ppi_based = pmin(pmax(barren_fraction_ppi_based, 0), 1),
          true_veg_cover_ppi = 1 - barren_fraction_ppi_based
        ) |>
        dplyr::distinct()
    }
    
    # --- 1. Identify Valid Locations and Get True Labels from a truly held-out dataset ---
    if (!exists("df_tasks_inference") || is.null(df_tasks_inference) || nrow(df_tasks_inference) == 0) {
      stop("[VALIDATION] No held-out dataset available (df_tasks_inference is empty). Provide a separate INFERENCE_CSV with ground truth labels.")
    }
    if (!"Veg" %in% names(df_tasks_inference)) {
      stop("[VALIDATION] Held-out dataset (df_tasks_inference) is missing ground-truth 'Veg' labels")
    }
    if (exists("df_train") && !is.null(df_train) && nrow(df_train) > 0 && "location_id" %in% names(df_train)) {
      overlap_ids <- intersect(unique(as.character(df_train$location_id)), unique(as.character(df_tasks_inference$location_id)))
      if (length(overlap_ids) > 0) {
        stop(sprintf("[VALIDATION] Held-out dataset overlaps training locations (n_overlap=%d). Refusing to validate on non-held-out data.", length(overlap_ids)))
      }
    }

    labels_df <- df_tasks_inference %>% dplyr::select(location_id, Veg) %>% dplyr::distinct()
    
    # Normalize labels
    labels_df$Veg <- normalize_veg_name(labels_df$Veg)
    labels_df <- labels_df %>% dplyr::rename(true_veg = Veg)
    
    # Filter validation coefficients to only those with labels
    # We join first to filter efficiently
    val_coefs_labeled <- val_coefs %>% 
       dplyr::inner_join(labels_df, by = "location_id")

    # Validation is strictly held-out; do not filter to TRAIN_YEARS.
    
    # Check class distribution in validation set
    if (nrow(val_coefs_labeled) > 0 && "true_veg" %in% names(val_coefs_labeled)) {
      val_class_dist <- table(val_coefs_labeled$true_veg)
      cat(sprintf("[VALIDATION] Class distribution: %s\n", 
                  paste(names(val_class_dist), "=", val_class_dist, collapse=", ")))
    }



    
    cat(sprintf("\n=== VALIDATION ACCURACY ON HELD-OUT SET%s ===\n", prefix))
    cat(sprintf("Validation set: %d locations, %d location-year pairs\n", 
                length(unique(val_coefs_labeled$location_id)), 
                length(unique(paste(val_coefs_labeled$location_id, val_coefs_labeled$pheno_year)))))

    # --- 2. Pivot Long to Wide (Correctly) ---
    # Determine which column holds the value
    measure_col <- dplyr::case_when(
      "coef" %in% names(val_coefs_labeled) ~ "coef",
      "pred_coef_rel" %in% names(val_coefs_labeled) ~ "pred_coef_rel",
      TRUE ~ NA_character_
    )
    
    val_coefs_wide <- NULL
    
    if (is.na(measure_col)) {
       # If no single value column, check if it's already wide (has frac_* columns)
       veg_cols <- grep("^frac_", names(val_coefs_labeled), value = TRUE)
       if(length(veg_cols) > 0) {
          # Already wide, just ensure true_veg is present (it is, from the join)
          val_coefs_wide <- val_coefs_labeled
       } else {
          cat("[ERROR] Could not identify coefficient column (coef) or fraction columns (frac_*) in input.\n")
          return(invisible(NULL))
       }
    } else {
       # Pivot
       # We strictly use location_id, pheno_year, and true_veg as keys. 
       # 'Veg' in val_coefs is the PREDICTED class (if present).
       # Note: pivot_wider will take 'Veg' (predicted) values to make columns.
       
       # Ensure 'Veg' exists for pivoting
       if (!"Veg" %in% names(val_coefs_labeled)) {
          cat("[ERROR] 'Veg' column (predicted class) missing for pivoting.\n")
          return(invisible(NULL))
       }
       
       # AGENT: Aggregate variants to vegetation class level (e.g. "herbs_opt_1" -> "herbs")
       val_coefs_wide <- val_coefs_labeled %>%
         dplyr::mutate(
             # Remove common suffixes to get base class: _opt_N, _single, _ppi
             Veg_class = sub("(_opt_[0-9]+|_single|_ppi)$", "", Veg)
         ) %>%
         dplyr::select(location_id, pheno_year, true_veg, Veg_class, dplyr::all_of(measure_col)) %>%
         tidyr::pivot_wider(names_from = Veg_class, 
                            values_from = dplyr::all_of(measure_col), 
                            values_fill = 0, 
                            values_fn = sum, 
                            names_prefix = "frac_")
    }

    # Identify all fraction columns (now aggregated by class)
    all_veg_cols <- grep("^frac_", names(val_coefs_wide), value = TRUE)
    


    # For PPI analysis, we might exclude barren if strictly analyzing "vegetation" cover
    veg_cols_no_barren <- all_veg_cols[all_veg_cols != "frac_barren"]

    # --- PPI-based total vegetation cover validation (truth from PPI-derived cover) ---
    if (!is.null(ppi_truth) && !is.null(val_coefs_wide) && length(veg_cols_no_barren) > 0) {
      ppi_join <- val_coefs_wide %>%
        dplyr::left_join(ppi_truth, by = c("location_id", "pheno_year"))

      if ("true_veg_cover_ppi" %in% names(ppi_join)) {
        ppi_join$pred_total_veg_cover <- rowSums(ppi_join[, veg_cols_no_barren, drop = FALSE], na.rm = TRUE)
        ppi_join$true_veg_cover_ppi <- pmin(pmax(ppi_join$true_veg_cover_ppi, 0), 1)

        cover_rows <- ppi_join[is.finite(ppi_join$true_veg_cover_ppi), , drop = FALSE]
        if (nrow(cover_rows) > 0) {
          mae_cover <- mean(abs(cover_rows$pred_total_veg_cover - cover_rows$true_veg_cover_ppi), na.rm = TRUE)
          bias_cover <- mean(cover_rows$pred_total_veg_cover - cover_rows$true_veg_cover_ppi, na.rm = TRUE)
          cat("\n--- PPI-BASED VEGETATION COVER VALIDATION ---\n")
          cat(sprintf("Samples with PPI truth: %d\n", nrow(cover_rows)))
          cat(sprintf("Mean absolute error (pred veg cover vs PPI): %.3f\n", mae_cover))
          cat(sprintf("Mean bias (pred - truth): %.3f\n", bias_cover))
        } else {
          cat("[NOTICE] PPI truth found but no matching location-year rows after filtering.\n")
        }
      }
    }
    
    # --- 3. Compute Mean Predicted Fraction ---
    # Use raw fractions directly (MESMA coefficients already sum to 1)
    if (length(all_veg_cols) > 0) {
      for (v in unique(val_coefs_wide$true_veg)) {
        if (tolower(v) == "barren") next
        sub <- val_coefs_wide[val_coefs_wide$true_veg == v, ]
        if (nrow(sub) > 0) {
          frac_col <- paste0("frac_", v) # e.g. frac_populus
          
          if (frac_col %in% names(sub)) {
            # Use raw fraction directly (no rescaling needed - MESMA coeffs sum to 1)
            mean_frac <- mean(sub[[frac_col]], na.rm = TRUE)
            cat(sprintf("  %s: mean predicted fraction = %.3f\n", v, mean_frac))
          } else {
             # If the column doesn't exist, it means the model NEVER predicted this class
             cat(sprintf("  %s: mean predicted fraction = 0.000 (never predicted)\n", v))
          }
        }
      }
    }

    # --- 4. Row-Normalized Confusion Matrix (excluding barren) ---
    if (length(all_veg_cols) > 0 && nrow(val_coefs_wide) > 0) {
      # Filter to non-barren vegetation classes only
      matrix_veg_cols <- all_veg_cols[!grepl("barren", all_veg_cols, ignore.case = TRUE)]
      
      if (length(matrix_veg_cols) == 0) {
        cat("\n[NOTICE] No non-barren vegetation columns - skipping confusion matrix\n")
      } else {
        # Filter rows to non-barren true classes
        val_coefs_veg <- val_coefs_wide %>%
          dplyr::filter(tolower(true_veg) != "barren", !is.na(true_veg), true_veg != "")
        
        n_true_classes <- length(unique(val_coefs_veg$true_veg))
        
        if (n_true_classes < 2) {
          cat(sprintf("\n[NOTICE] Only %d non-barren class in validation - skipping confusion matrix\n", n_true_classes))
        } else {
          cat("\n=== CONFUSION MATRIX (Row-Normalized, Excluding Barren) ===\n")
          cat("Rows: True Class | Columns: Mean Predicted Fraction\n")
          cat("(Each row normalized to sum to 1.0)\n\n")
          
          # For each prediction, normalize vegetation fractions (excluding barren) to sum to 1
          # NOTE: When veg_sum is near zero, preserve original values instead of zeroing them
          val_coefs_norm <- val_coefs_veg %>%
            dplyr::mutate(
              veg_sum = rowSums(dplyr::across(dplyr::all_of(matrix_veg_cols)), na.rm = TRUE)
            ) %>%
            dplyr::mutate(
              dplyr::across(
                dplyr::all_of(matrix_veg_cols),
                ~ ifelse(veg_sum > 1e-9, .x / veg_sum, .x),
                .names = "{.col}_norm"
              )
            )
          
          # Column names after normalization
          norm_cols <- paste0(matrix_veg_cols, "_norm")
          
          # Average by true class
          avg_fractions <- val_coefs_norm %>%
            dplyr::group_by(true_veg) %>%
            dplyr::summarize(
              dplyr::across(dplyr::all_of(norm_cols), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop"
            )
          
          # Clean column names (remove frac_ prefix and _norm suffix)
          clean_names <- sub("_norm$", "", colnames(avg_fractions))
          clean_names <- sub("^frac_", "", clean_names)
          colnames(avg_fractions) <- clean_names
          
          # Convert to matrix
          mat_data <- as.matrix(avg_fractions[, -1])  # Remove true_veg column
          rownames(mat_data) <- avg_fractions$true_veg
          
          # Ensure square matrix with all classes in both dimensions
          all_classes <- sort(unique(c(rownames(mat_data), colnames(mat_data))))
          
          # Add missing columns
          for (cls in all_classes) {
            if (!cls %in% colnames(mat_data)) {
              new_col <- matrix(0, nrow = nrow(mat_data), ncol = 1)
              colnames(new_col) <- cls
              mat_data <- cbind(mat_data, new_col)
            }
          }
          
          # Add missing rows
          for (cls in all_classes) {
            if (!cls %in% rownames(mat_data)) {
              new_row <- matrix(0, nrow = 1, ncol = ncol(mat_data))
              rownames(new_row) <- cls
              mat_data <- rbind(mat_data, new_row)
            }
          }
          
          # Reorder to match all_classes
          mat_data <- mat_data[all_classes, all_classes, drop = FALSE]
          
          # Row-normalize (already done per-prediction, but verify)
          mat_rownorm <- mat_data
          row_sums <- rowSums(mat_rownorm, na.rm = TRUE)
          for (i in seq_len(nrow(mat_rownorm))) {
            if (row_sums[i] > 1e-9) {
              mat_rownorm[i, ] <- mat_rownorm[i, ] / row_sums[i]
            }
          }
          
          # Display matrix
          print(round(mat_rownorm, 3))
          
          # Compute diagonal accuracy
          diag_vals <- diag(mat_rownorm)
          mean_diagonal <- mean(diag_vals, na.rm = TRUE)
          
          cat(sprintf("\nMean diagonal (correctly predicted fraction): %.3f (%.1f%%)\n",
                      mean_diagonal, mean_diagonal * 100))
          cat(sprintf("Per-class accuracy:\n"))
          for (i in seq_along(all_classes)) {
            cat(sprintf("  %s: %.3f (%.1f%%)\n",
                        all_classes[i], diag_vals[i], diag_vals[i] * 100))
          }

          # Store confusion matrix for Dirichlet perturbation in bootstrap
          # Compute validation sample sizes per class
          val_sample_sizes <- table(val_coefs_veg$true_veg)
          val_sample_sizes <- setNames(as.numeric(val_sample_sizes), names(val_sample_sizes))

          store_confusion_matrix(mat_rownorm, sample_sizes = val_sample_sizes)
          cat(sprintf("\n[DIRICHLET] Confusion matrix stored for classification uncertainty propagation\n"))

          # Store OOB fraction residuals for MC uncertainty propagation
          # Residual = predicted_fraction - true_fraction (where true is 1 for own class, 0 for others)
          if (isTRUE(ENABLE_OOB_FRACTION_UNCERTAINTY)) {
            cat("\n[OOB_FRAC] Computing and storing OOB fraction residuals for MC uncertainty...\n")

            # Build residuals by true class
            residuals_by_class <- list()

            for (true_cls in all_classes) {
              # Get rows where true class matches
              cls_data <- val_coefs_norm[val_coefs_norm$true_veg == true_cls, , drop = FALSE]
              if (nrow(cls_data) == 0) next

              # Build true fraction vector: 1 for own class, 0 for others
              # Compute residuals: predicted - true
              n_samples <- nrow(cls_data)
              resid_mat <- matrix(0, nrow = n_samples, ncol = length(all_classes))
              colnames(resid_mat) <- all_classes

              for (i in seq_len(n_samples)) {
                for (pred_cls in all_classes) {
                  # Get predicted fraction (normalized)
                  pred_col <- paste0("frac_", pred_cls, "_norm")
                  pred_frac <- if (pred_col %in% names(cls_data)) cls_data[[pred_col]][i] else 0
                  if (!is.finite(pred_frac)) pred_frac <- 0

                  # True fraction: 1 if pred_cls == true_cls, else 0
                  true_frac <- if (pred_cls == true_cls) 1.0 else 0.0

                  # Residual = predicted - true
                  resid_mat[i, pred_cls] <- pred_frac - true_frac
                }
              }

              residuals_by_class[[true_cls]] <- resid_mat
            }

            if (length(residuals_by_class) > 0) {
              store_oob_fraction_residuals(residuals_by_class)
              cat(sprintf("[OOB_FRAC] Stored residuals for %d classes\n", length(residuals_by_class)))
            }
          }
        }
      }
    }

    # --- 5. Artificial Mix Logic ---
    if (length(all_veg_cols) > 0) {
      cat("DEBUG: Entering artificial mix section\n")
      # --- INTER-CLASS MIXTURE DISCRIMINATION TEST ---
      # Mix different vegetation classes to test if the model coefficients reflect the mix.
      # Note: This averages predictions (coefficients), so it tests linearity of the output space.

      # Only consider true vegetation classes (exclude barren/NA/empty)
      veg_classes_only <- unique(val_coefs_wide$true_veg[
        !is.na(val_coefs_wide$true_veg) &
          trimws(as.character(val_coefs_wide$true_veg)) != "" &
          tolower(val_coefs_wide$true_veg) != "barren"
      ])
      veg_classes_only <- sort(unique(tolower(as.character(veg_classes_only))))

      # Keep only classes that have a corresponding fraction column (model actually produced that column)
      pred_classes <- sub("^frac_", "", all_veg_cols)
      pred_classes <- tolower(pred_classes)
      pred_classes <- setdiff(pred_classes, "barren")
      veg_classes_only <- intersect(veg_classes_only, pred_classes)
      
      if (length(veg_classes_only) >= 2) {
        cat("\n--- INTER-CLASS MIXTURE DISCRIMINATION (Synthetic 50/50 Mixes) ---\n")
        
        # Generate all unique pairs of different vegetation classes
        class_pairs <- utils::combn(veg_classes_only, 2, simplify = FALSE)
        
        for (pair in class_pairs) {
          class_a <- pair[1]
          class_b <- pair[2]
          
          sub_a <- val_coefs_wide[val_coefs_wide$true_veg == class_a, ]
          sub_b <- val_coefs_wide[val_coefs_wide$true_veg == class_b, ]
          
          if (nrow(sub_a) > 0 && nrow(sub_b) > 0) {
             # Sample pairs to create mixes
             n_mix <- min(nrow(sub_a), nrow(sub_b), 100)
             idx_a <- sample(nrow(sub_a), n_mix, replace = TRUE)
             idx_b <- sample(nrow(sub_b), n_mix, replace = TRUE)
             
             # Create mixed coefficients (average of A and B predictions)
             mixed_coefs <- (sub_a[idx_a, all_veg_cols, drop = FALSE] + sub_b[idx_b, all_veg_cols, drop = FALSE]) / 2

             # Report using ROW-NORMALIZED vegetation fractions (exclude barren).
             veg_frac_cols_local <- all_veg_cols[!grepl("barren", all_veg_cols, ignore.case = TRUE)]
             if (length(veg_frac_cols_local) == 0) next

             mixed_norm <- mixed_coefs
             veg_sum <- rowSums(mixed_norm[, veg_frac_cols_local, drop = FALSE], na.rm = TRUE)
             for (cc in veg_frac_cols_local) {
               mixed_norm[[cc]] <- ifelse(veg_sum > 1e-9, mixed_norm[[cc]] / veg_sum, mixed_norm[[cc]])
             }

             frac_col_a <- paste0("frac_", tolower(class_a))
             frac_col_b <- paste0("frac_", tolower(class_b))
             if (!(frac_col_a %in% names(mixed_norm)) || !(frac_col_b %in% names(mixed_norm))) next

             mean_a <- mean(mixed_norm[[frac_col_a]], na.rm = TRUE)
             mean_b <- mean(mixed_norm[[frac_col_b]], na.rm = TRUE)
             mean_sum_ab <- mean(mixed_norm[[frac_col_a]] + mixed_norm[[frac_col_b]], na.rm = TRUE)

             excluded_cols <- setdiff(veg_frac_cols_local, c(frac_col_a, frac_col_b))
             mean_excluded_sum <- if (length(excluded_cols) > 0) {
               mean(rowSums(mixed_norm[, excluded_cols, drop = FALSE], na.rm = TRUE), na.rm = TRUE)
             } else {
               0
             }

             cat(sprintf(
               "  Mix %s + %s: mean(row-norm %s)=%.3f, mean(row-norm %s)=%.3f, mean(sum two)=%.3f (Expected ~1.000), mean(excluded sum)=%.3f (Expected ~0.000)\n",
               class_a, class_b, class_a, mean_a, class_b, mean_b, mean_sum_ab, mean_excluded_sum
             ))
          }
        }
      }
    }

    # --- 6. Compute Per-Class Prediction RMSE (for CI Inflation) ---
    per_class_rmse <- numeric(0)
    if (length(all_veg_cols) > 0 && nrow(val_coefs_wide) > 0) {
      cat("\n--- VALIDATION RMSE (for CI adjustment) ---\n")
      for (v_class in unique(c(val_coefs_wide$true_veg, sub("^frac_", "", all_veg_cols)))) {
        frac_col <- paste0("frac_", tolower(v_class))
        if(frac_col %in% names(val_coefs_wide)) {
           predicted <- val_coefs_wide[[frac_col]]
           # Truth is 1 if true_veg matches class, 0 otherwise
           truth <- ifelse(tolower(val_coefs_wide$true_veg) == tolower(v_class), 1, 0)
           
           # Filter out NAs
           valid_idx <- !is.na(predicted) & !is.na(truth)
           if (sum(valid_idx) > 2) {
             mse <- mean((predicted[valid_idx] - truth[valid_idx])^2)
             rmse_val <- sqrt(mse)
             per_class_rmse[tolower(v_class)] <- rmse_val
             cat(sprintf("  %s: RMSE = %.4f (N=%d)\n", v_class, rmse_val, sum(valid_idx)))
           }
        }
      }
    }

    cat("Validation accuracy computed", prefix, ".\n", sep = "")
    return(per_class_rmse)
  }

  # Compute accuracy on VALIDATION data only (held-out locations, TRAIN_YEARS)
  validation_rmse_adjustment <- NULL
  if (exists("validation_coefs") && !is.null(validation_coefs) && nrow(validation_coefs) > 0) {
    cat("\n=== VALIDATION ACCURACY (held-out locations, TRAIN_YEARS only) ===\n")
    val_coefs <- validation_coefs
    validation_rmse_adjustment <- report_validation_accuracy(val_coefs, "held-out validation set")
  } else {
    cat("[WARNING] No validation coefficients available for accuracy computation\n")
  }

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

  # ==========================================================================
  # ARTIFICIAL 50/50 MIX VALIDATION (Spectral-level mixing)
  # Create synthetic mixed spectra from validation samples and validate unmixing
  # ==========================================================================
  cat("\n=== ARTIFICIAL 50/50 MIX VALIDATION (Spectral-level) ===\n")

  artificial_mix_validation <- function() {
    # Check prerequisites
    if (!exists("df_validation") || is.null(df_validation) || nrow(df_validation) == 0) {
      cat("[SKIP] No df_validation available for artificial mixing\n")
      return(NULL)
    }

    if (!exists("MESMA_PARAMS") || is.null(MESMA_PARAMS) || (!("indices" %in% names(MESMA_PARAMS)) && !("base_indices" %in% names(MESMA_PARAMS)))) {
      cat("[SKIP] MESMA_PARAMS not available - cannot identify feature columns\n")
      return(NULL)
    }

    # Use the same raw feature columns as fit_one_task() uses for build_pentad_matrix()
    spectral_cols <- if (!is.null(MESMA_PARAMS$base_indices) && length(MESMA_PARAMS$base_indices) > 0) MESMA_PARAMS$base_indices else MESMA_PARAMS$indices
    if (length(spectral_cols) == 0) {
      cat("[SKIP] No spectral columns identified\n")
      return(NULL)
    }

    # Check that spectral columns exist in df_validation
    available_cols <- intersect(spectral_cols, names(df_validation))
    if (length(available_cols) < length(spectral_cols) * 0.5) {
      cat(sprintf("[SKIP] Too few spectral columns available in df_validation (%d/%d)\n",
                  length(available_cols), length(spectral_cols)))
      return(NULL)
    }
    spectral_cols <- available_cols

    # Get vegetation classes (excluding barren for mixing)
    df_val_with_veg <- df_validation
    if (!"Veg" %in% names(df_val_with_veg)) {
      cat("[SKIP] No Veg column in df_validation\n")
      return(NULL)
    }
    df_val_with_veg$Veg <- normalize_veg_name(df_val_with_veg$Veg)
    veg_classes <- unique(df_val_with_veg$Veg[df_val_with_veg$Veg != "barren"])

    if (length(veg_classes) < 2) {
      cat(sprintf("[SKIP] Need at least 2 vegetation classes for mixing, found: %d\n", length(veg_classes)))
      return(NULL)
    }

    cat(sprintf("[MIX] Using %d spectral columns and %d vegetation classes: %s\n",
                length(spectral_cols), length(veg_classes), paste(veg_classes, collapse=", ")))

    # Filter to TRAIN_YEARS if defined
    if (exists("TRAIN_YEARS") && !is.null(TRAIN_YEARS) && "pheno_year" %in% names(df_val_with_veg)) {
      df_val_with_veg <- df_val_with_veg[df_val_with_veg$pheno_year %in% TRAIN_YEARS, ]
      cat(sprintf("[MIX] Filtered to TRAIN_YEARS: %d rows\n", nrow(df_val_with_veg)))
    }

    # Aggregate spectra to location-year level (mean across observations)
    # This gives us one "representative spectrum" per location-year
    required_cols <- c("location_id", "pheno_year", "Veg", spectral_cols)
    if ("doy" %in% names(df_val_with_veg)) required_cols <- c(required_cols, "doy")
    if ("date" %in% names(df_val_with_veg)) required_cols <- c(required_cols, "date")
    if ("lat" %in% names(df_val_with_veg)) required_cols <- c(required_cols, "lat")
    if ("lon" %in% names(df_val_with_veg)) required_cols <- c(required_cols, "lon")

    df_val_subset <- df_val_with_veg[, intersect(required_cols, names(df_val_with_veg)), drop = FALSE]

    # Get unique location-years per class
    loc_year_by_class <- list()
    for (vc in veg_classes) {
      class_data <- df_val_subset[df_val_subset$Veg == vc, ]
      if (nrow(class_data) > 0) {
        # Get unique location-year combinations
        class_data$loc_year <- paste(class_data$location_id, class_data$pheno_year, sep = "_")
        loc_year_by_class[[vc]] <- unique(class_data$loc_year)
      }
    }

    # Generate all pairs of vegetation classes
    class_pairs <- utils::combn(veg_classes, 2, simplify = FALSE)

    mix_results <- list()
    n_mixes_total <- 0

    for (pair in class_pairs) {
      class_a <- pair[1]
      class_b <- pair[2]

      loc_years_a <- loc_year_by_class[[class_a]]
      loc_years_b <- loc_year_by_class[[class_b]]

      if (is.null(loc_years_a) || is.null(loc_years_b) ||
          length(loc_years_a) == 0 || length(loc_years_b) == 0) {
        cat(sprintf("[MIX] Skipping %s + %s: insufficient samples\n", class_a, class_b))
        next
      }

      # Create up to 50 mixes per pair
      n_mix <- min(length(loc_years_a), length(loc_years_b), 50)

      cat(sprintf("[MIX] Creating %d synthetic 50/50 mixes: %s + %s\n", n_mix, class_a, class_b))

      set.seed(get_mesma_seed(42) + which(sapply(class_pairs, function(p) identical(p, pair))))
      sampled_a <- sample(loc_years_a, n_mix, replace = TRUE)
      sampled_b <- sample(loc_years_b, n_mix, replace = TRUE)

      mixed_data_list <- list()

      for (i in seq_len(n_mix)) {
        # Get data for each component
        loc_year_a <- sampled_a[i]
        loc_year_b <- sampled_b[i]

        parts_a <- strsplit(loc_year_a, "_")[[1]]
        parts_b <- strsplit(loc_year_b, "_")[[1]]

        loc_a <- paste(parts_a[-length(parts_a)], collapse = "_")
        yr_a <- as.integer(parts_a[length(parts_a)])
        loc_b <- paste(parts_b[-length(parts_b)], collapse = "_")
        yr_b <- as.integer(parts_b[length(parts_b)])

        data_a <- df_val_subset[df_val_subset$location_id == loc_a &
                                 df_val_subset$pheno_year == yr_a &
                                 df_val_subset$Veg == class_a, ]
        data_b <- df_val_subset[df_val_subset$location_id == loc_b &
                                 df_val_subset$pheno_year == yr_b &
                                 df_val_subset$Veg == class_b, ]

        if (nrow(data_a) == 0 || nrow(data_b) == 0) next

        # Match observations by DOY if possible, otherwise use mean spectra
        if ("doy" %in% names(data_a) && "doy" %in% names(data_b)) {
          # Find common or closest DOYs
          common_doys <- intersect(data_a$doy, data_b$doy)

          if (length(common_doys) >= 3) {
            # Use common DOYs
            data_a_matched <- data_a[data_a$doy %in% common_doys, ]
            data_b_matched <- data_b[data_b$doy %in% common_doys, ]

            # Ensure same number of rows by matching DOYs
            data_a_matched <- data_a_matched[order(data_a_matched$doy), ]
            data_b_matched <- data_b_matched[order(data_b_matched$doy), ]

            # Take first occurrence per DOY
            data_a_matched <- data_a_matched[!duplicated(data_a_matched$doy), ]
            data_b_matched <- data_b_matched[!duplicated(data_b_matched$doy), ]

            common_doys_final <- intersect(data_a_matched$doy, data_b_matched$doy)
            data_a_matched <- data_a_matched[data_a_matched$doy %in% common_doys_final, ]
            data_b_matched <- data_b_matched[data_b_matched$doy %in% common_doys_final, ]
          } else {
            # Not enough common DOYs, aggregate to mean
            data_a_matched <- data_a[1, , drop = FALSE]
            data_b_matched <- data_b[1, , drop = FALSE]
            for (sc in spectral_cols) {
              data_a_matched[[sc]] <- mean(data_a[[sc]], na.rm = TRUE)
              data_b_matched[[sc]] <- mean(data_b[[sc]], na.rm = TRUE)
            }
            data_a_matched$doy <- median(data_a$doy, na.rm = TRUE)
            data_b_matched$doy <- median(data_b$doy, na.rm = TRUE)
          }
        } else {
          # No DOY, use first row with mean values
          data_a_matched <- data_a[1, , drop = FALSE]
          data_b_matched <- data_b[1, , drop = FALSE]
          for (sc in spectral_cols) {
            data_a_matched[[sc]] <- mean(data_a[[sc]], na.rm = TRUE)
            data_b_matched[[sc]] <- mean(data_b[[sc]], na.rm = TRUE)
          }
        }

        # Create 50/50 mixed spectra
        n_obs <- min(nrow(data_a_matched), nrow(data_b_matched))
        if (n_obs == 0) next

        mixed_data <- data_a_matched[1:n_obs, , drop = FALSE]
        mixed_data$location_id <- sprintf("MIX_%s_%s_%d", class_a, class_b, i)
        mixed_data$pheno_year <- yr_a  # Use year from A
        mixed_data$Veg <- paste0("mix_", class_a, "_", class_b)
        mixed_data$true_frac_a <- 0.5
        mixed_data$true_frac_b <- 0.5
        mixed_data$class_a <- class_a
        mixed_data$class_b <- class_b

        # Mix spectral values: 0.5 * A + 0.5 * B
        for (sc in spectral_cols) {
          val_a <- data_a_matched[[sc]][1:n_obs]
          val_b <- data_b_matched[[sc]][1:n_obs]
          if (length(val_a) == n_obs && length(val_b) == n_obs) {
            mixed_data[[sc]] <- 0.5 * val_a + 0.5 * val_b
          }
        }

        mixed_data_list[[length(mixed_data_list) + 1]] <- mixed_data
      }

      if (length(mixed_data_list) > 0) {
        pair_key <- paste(class_a, class_b, sep = "_")
        mix_results[[pair_key]] <- list(
          class_a = class_a,
          class_b = class_b,
          mixed_data = do.call(rbind, mixed_data_list)
        )
        n_mixes_total <- n_mixes_total + nrow(mix_results[[pair_key]]$mixed_data)
      }
    }

    if (n_mixes_total == 0) {
      cat("[MIX] No valid mixes created\n")
      return(NULL)
    }

    cat(sprintf("[MIX] Created %d total artificial mix samples across %d class pairs\n",
                n_mixes_total, length(mix_results)))

    # Run unmixing on mixed data
    cat("[MIX] Running unmixing on artificial mixes...\n")

    all_mix_coefs <- list()

    for (pair_key in names(mix_results)) {
      mr <- mix_results[[pair_key]]
      mixed_df <- mr$mixed_data

      if (is.null(mixed_df) || nrow(mixed_df) == 0) next

      # Ensure required columns for unmixing
      if (!"date" %in% names(mixed_df) && "doy" %in% names(mixed_df) && "pheno_year" %in% names(mixed_df)) {
        # Approximate date from DOY and year
        mixed_df$date <- as.Date(paste0(mixed_df$pheno_year, "-01-01")) + mixed_df$doy - 1
      }

      # Run inference using run_inference_silent if available
      if (exists("run_inference_silent")) {
        inf_result <- tryCatch({
          run_inference_silent(mixed_df)
        }, error = function(e) {
          cat(sprintf("[MIX ERROR] %s: %s\n", pair_key, e$message))
          NULL
        })

        if (!is.null(inf_result) && !is.null(inf_result$all_coefs) && nrow(inf_result$all_coefs) > 0) {
          coefs <- inf_result$all_coefs
          coefs$class_a <- mr$class_a
          coefs$class_b <- mr$class_b
          coefs$true_frac_a <- 0.5
          coefs$true_frac_b <- 0.5
          all_mix_coefs[[pair_key]] <- coefs
        }
      }
    }

    if (length(all_mix_coefs) == 0) {
      cat("[MIX] Unmixing produced no results\n")
      return(NULL)
    }

    mix_coefs_combined <- do.call(rbind, all_mix_coefs)

    # Analyze results
    cat("\n--- ARTIFICIAL MIX VALIDATION RESULTS ---\n")
    cat(sprintf("Total unmixed mix samples: %d\n", length(unique(mix_coefs_combined$location_id))))

    # For each class pair, check if predicted fractions match expected 50/50
    for (pair_key in names(all_mix_coefs)) {
      coefs <- all_mix_coefs[[pair_key]]
      class_a <- coefs$class_a[1]
      class_b <- coefs$class_b[1]

      # Pivot to wide format
      coefs_wide <- coefs %>%
        dplyr::mutate(Veg_class = sub("(_opt_[0-9]+|_single|_ppi)$", "", tolower(Veg))) %>%
        dplyr::group_by(location_id, pheno_year, Veg_class) %>%
        dplyr::summarise(coef = sum(coef, na.rm = TRUE), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = Veg_class, values_from = coef, values_fill = 0)

      frac_col_a <- tolower(class_a)
      frac_col_b <- tolower(class_b)

      if (frac_col_a %in% names(coefs_wide) && frac_col_b %in% names(coefs_wide)) {
        pred_frac_a <- coefs_wide[[frac_col_a]]
        pred_frac_b <- coefs_wide[[frac_col_b]]

        # Calculate metrics
        mean_pred_a <- mean(pred_frac_a, na.rm = TRUE)
        mean_pred_b <- mean(pred_frac_b, na.rm = TRUE)
        sd_pred_a <- sd(pred_frac_a, na.rm = TRUE)
        sd_pred_b <- sd(pred_frac_b, na.rm = TRUE)

        # Expected is 0.5 for both
        rmse_a <- sqrt(mean((pred_frac_a - 0.5)^2, na.rm = TRUE))
        rmse_b <- sqrt(mean((pred_frac_b - 0.5)^2, na.rm = TRUE))

        cat(sprintf("\n  Mix %s + %s (N=%d):\n", class_a, class_b, nrow(coefs_wide)))
        cat(sprintf("    %s: predicted=%.3f ± %.3f (expected=0.500, RMSE=%.3f)\n",
                    class_a, mean_pred_a, sd_pred_a, rmse_a))
        cat(sprintf("    %s: predicted=%.3f ± %.3f (expected=0.500, RMSE=%.3f)\n",
                    class_b, mean_pred_b, sd_pred_b, rmse_b))
        cat(sprintf("    Sum: %.3f (expected=1.000)\n", mean_pred_a + mean_pred_b))
      }
    }

    return(mix_coefs_combined)
  }

  # Run artificial mix validation
  artificial_mix_coefs <- tryCatch({
    artificial_mix_validation()
  }, error = function(e) {
    cat(sprintf("[MIX ERROR] Artificial mix validation failed: %s\n", e$message))
    NULL
  })

  if (!is.null(artificial_mix_coefs) && nrow(artificial_mix_coefs) > 0) {
    # Overall accuracy summary across all pairs
    mix_veg_classes <- unique(artificial_mix_coefs$class_a)
    mix_veg_classes <- union(mix_veg_classes, unique(artificial_mix_coefs$class_b))
    overall_errors <- c()
    for (loc_id in unique(artificial_mix_coefs$location_id)) {
      loc_coefs <- artificial_mix_coefs[artificial_mix_coefs$location_id == loc_id, ]
      ca <- loc_coefs$class_a[1]; cb <- loc_coefs$class_b[1]
      frac_a <- sum(loc_coefs$coef[tolower(loc_coefs$Veg) == tolower(ca)], na.rm = TRUE)
      frac_b <- sum(loc_coefs$coef[tolower(loc_coefs$Veg) == tolower(cb)], na.rm = TRUE)
      overall_errors <- c(overall_errors, abs(frac_a - 0.5), abs(frac_b - 0.5))
    }
    cat(sprintf("\n  Overall artificial unmixing MAE: %.3f (across %d mix samples)\n",
                mean(overall_errors, na.rm = TRUE), length(unique(artificial_mix_coefs$location_id))))
  }

  cat("=== END ARTIFICIAL MIX VALIDATION ===\n\n")

  # Results are now in validation_coefs and inference_coefs (combined as all_coefs)
  cat(sprintf("[RESULTS] Validation: %d rows, Inference: %d rows, Combined: %d rows\n",
              if(exists("validation_coefs") && !is.null(validation_coefs)) nrow(validation_coefs) else 0,
              if(exists("inference_coefs") && !is.null(inference_coefs)) nrow(inference_coefs) else 0,
              if(exists("all_coefs") && !is.null(all_coefs)) nrow(all_coefs) else 0))
              
  # --- Apply Validation RMSE Adjustment to CIs ---
  if (!is.null(validation_rmse_adjustment) && length(validation_rmse_adjustment) > 0) {
    cat("\n[UNCERTAINTY ADJUSTMENT] Incorporating validation RMSE into Confidence Intervals...\n")
    
    # Function to adjust CIs
    adjust_ci <- function(df, rmse_vec) {
      if (is.null(df) || nrow(df) == 0) return(df)
      changed_rows <- 0
      
      # Iterate over known classes in RMSE vector
      for (v_class in names(rmse_vec)) {
        rmse_val <- rmse_vec[[v_class]]
        # Find rows for this class (case-insensitive)
        idx <- which(normalize_veg_name(df$Veg) == tolower(v_class))
        
        if (length(idx) > 0) {
           # Extract current intervals
           current_lower <- df$coef_025[idx]
           current_upper <- df$coef_975[idx]
           current_coef <- df$coef[idx]
           
           # Approximate standard error including validation RMSE
           # Current Width ~ 2 * 1.96 * SE_fit
           # New SE ~ sqrt(SE_fit^2 + RMSE_val^2)
           
           # Recover implicit SE from CI width
           # If width is 0 (or NA), assume SE is 0 (or skip)
           width <- current_upper - current_lower
           width[is.na(width)] <- 0
           implicit_se <- width / (2 * 1.96)
           
           # Combine variances
           new_se <- sqrt(implicit_se^2 + rmse_val^2)
           
           # Construct new intervals centered on prediction
           # (We assume symmetric error distribution governed by new_se)
           new_lower <- current_coef - 1.96 * new_se
           new_upper <- current_coef + 1.96 * new_se
           
           # Clamp to [0, 1]
           new_lower <- pmax(0, pmin(1, new_lower))
           new_upper <- pmax(0, pmin(1, new_upper))
           
           # Update DF
           df$coef_025[idx] <- new_lower
           df$coef_975[idx] <- new_upper
           df$coef_sd[idx] <- new_se # Update SD if present
           
           changed_rows <- changed_rows + length(idx)
        }
      }
      return(df)
    }

    # Apply to all_coefs (which combines val & inference)
    if (!is.null(all_coefs) && nrow(all_coefs) > 0) {
       all_coefs <- adjust_ci(all_coefs, validation_rmse_adjustment)
    }
    # Also update constituent dataframes so they stay in sync if used later
    if (exists("inference_coefs") && !is.null(inference_coefs) && nrow(inference_coefs) > 0) {
       inference_coefs <- adjust_ci(inference_coefs, validation_rmse_adjustment)
    }
    if (exists("validation_coefs") && !is.null(validation_coefs) && nrow(validation_coefs) > 0) {
       validation_coefs <- adjust_ci(validation_coefs, validation_rmse_adjustment)
    }
    cat(sprintf("[UNCERTAINTY ADJUSTMENT] Adjusted CIs for ~%d rows across %d classes.\n", 
                nrow(all_coefs), length(validation_rmse_adjustment)))
  }

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
