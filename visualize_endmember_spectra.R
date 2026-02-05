# Visualization of endmember spectral features across pentads
# Each endmember gets its own plot, with each index/band as a separate colored line
#
# USAGE: Run this script AFTER fit_veg_mixture_mesma.R has been executed,
#        so that the required data structures (raw_lib_templates, mesma_lib, etc.)
#        are available in the global environment.
#
# Alternatively, uncomment the source() line below to run the main script first:
# source("fit_veg_mixture_mesma.R")

library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(dplyr)

#' Plot endmember spectral signatures by index
#'
#' Creates one plot per endmember showing all indices/bands across pentads
#' Each index is shown as a separate colored line
#'
#' @param lib The MESMA library (list of vegetation types, each with variants)
#' @param indices Vector of index names
#' @param out_dir Output directory for PNG files
#' @param prefix Filename prefix
#' @param save_png Whether to save plots as PNG
#' @param dpi Resolution for saved images
#' @return List of ggplot objects (invisibly)
plot_endmember_spectra_by_index <- function(lib, indices = NULL, out_dir = "endmember_spectra_plots",
                                             prefix = "endmember_spectra", save_png = TRUE, dpi = 150) {
  if (is.null(lib) || length(lib) == 0) {
    message("Library is empty or NULL")
    return(NULL)
  }

  # Get indices from MESMA_PARAMS if available
  if (is.null(indices)) {
    if (exists("MESMA_PARAMS") && !is.null(MESMA_PARAMS) && !is.null(MESMA_PARAMS$indices)) {
      indices <- MESMA_PARAMS$indices
    } else if (exists("OPTIMAL_INDICES")) {
      indices <- OPTIMAL_INDICES
    } else {
      stop("No indices provided and MESMA_PARAMS/OPTIMAL_INDICES not found")
    }
  }

  # Get temporal budget
  temporal_budget <- if (exists("TEMPORAL_BUDGET")) TEMPORAL_BUDGET else 37

  # Create output directory
  if (save_png && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Build long dataframe for plotting
  rows <- list()
  for (v in names(lib)) {
    variants <- lib[[v]]

    # Handle both list of variants and single variant structures
    if (!is.list(variants)) next
    if (!is.null(variants$vec) && is.numeric(variants$vec)) {
      # Single variant case
      variants <- list(variants)
    }

    for (var in variants) {
      vec <- as.numeric(var$vec)
      vid <- if (!is.null(var$id)) var$id else if (!is.null(var$variant_id)) var$variant_id else paste0(v, "_1")
      n_idx <- length(indices)
      expected_len <- n_idx * temporal_budget

      if (length(vec) < expected_len) {
        warning(sprintf("Skipping variant %s for veg %s: length(vec)=%d != expected=%d",
                        vid, v, length(vec), expected_len))
        next
      }

      for (k in seq_len(n_idx)) {
        idx_name <- indices[k]
        start <- (k-1) * temporal_budget + 1
        end <- k * temporal_budget
        vals <- vec[start:end]
        df_tmp <- data.frame(
          pentad = seq_len(temporal_budget),
          value = vals,
          Veg = v,
          variant_id = vid,
          index = idx_name,
          stringsAsFactors = FALSE
        )
        rows[[length(rows) + 1]] <- df_tmp
      }
    }
  }

  if (length(rows) == 0) {
    message("No valid data found in library")
    return(NULL)
  }

  proto_df <- do.call(rbind, rows)

  # Get unique vegetation types and indices
  veg_types <- unique(proto_df$Veg)
  all_indices <- unique(proto_df$index)

  # Create color palette for indices
  n_colors <- max(3, length(all_indices))
  if (n_colors <= 12) {
    index_colors <- brewer.pal(n = min(n_colors, 12), name = "Set3")
    if (n_colors > 12) {
      index_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_colors)
    }
  } else {
    # For many indices, use a rainbow palette
    index_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_colors)
  }
  names(index_colors) <- all_indices[seq_len(length(index_colors))]

  # Generate one plot per vegetation type (endmember)
  plots <- list()

  for (veg in veg_types) {
    df_veg <- proto_df[proto_df$Veg == veg, , drop = FALSE]

    # Average across variants if multiple exist
    df_avg <- df_veg %>%
      group_by(pentad, index) %>%
      summarize(
        value = mean(value, na.rm = TRUE),
        .groups = "drop"
      )

    # Create the plot
    p <- ggplot(df_avg, aes(x = pentad, y = value, color = index)) +
      geom_line(size = 1, alpha = 0.9) +
      scale_color_manual(values = index_colors, name = "Index/Band") +
      labs(
        title = sprintf("Endmember: %s", tools::toTitleCase(veg)),
        subtitle = "Spectral signature across pentads",
        x = "Pentad (10-day period)",
        y = "Value (z-score normalized)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
        legend.position = "right",
        legend.key.size = unit(0.8, "lines"),
        panel.grid.minor = element_blank()
      ) +
      guides(color = guide_legend(ncol = 1))

    plots[[veg]] <- p

    # Save to file
    if (save_png) {
      fn <- file.path(out_dir, sprintf("%s_%s.png", prefix, tolower(veg)))
      ggsave(filename = fn, plot = p, width = 12, height = 6, dpi = dpi)
      message(sprintf("Saved: %s", fn))
    }
  }

  # Also create a combined faceted plot
  df_avg_all <- proto_df %>%
    group_by(pentad, index, Veg) %>%
    summarize(
      value = mean(value, na.rm = TRUE),
      .groups = "drop"
    )

  p_combined <- ggplot(df_avg_all, aes(x = pentad, y = value, color = index)) +
    geom_line(size = 0.8, alpha = 0.85) +
    facet_wrap(~ Veg, scales = "free_y", ncol = 2) +
    scale_color_manual(values = index_colors, name = "Index/Band") +
    labs(
      title = "Endmember Spectral Signatures by Vegetation Type",
      x = "Pentad (10-day period)",
      y = "Value (z-score normalized)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "bottom",
      legend.key.size = unit(0.6, "lines"),
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(nrow = 3))

  plots[["combined"]] <- p_combined

  if (save_png) {
    fn <- file.path(out_dir, sprintf("%s_combined.png", prefix))
    ggsave(filename = fn, plot = p_combined, width = 14, height = 10, dpi = dpi)
    message(sprintf("Saved: %s", fn))
  }

  invisible(plots)
}

# Alternative: use raw_lib_templates if the optimized library isn't available
plot_raw_templates_by_index <- function(templates = NULL, feature_cols = NULL,
                                         out_dir = "endmember_spectra_plots",
                                         prefix = "raw_template_spectra",
                                         save_png = TRUE, dpi = 150) {

  # Get templates from global environment if not provided
  if (is.null(templates)) {
    if (exists("raw_lib_templates")) {
      templates <- raw_lib_templates
    } else {
      stop("No templates provided and raw_lib_templates not found in environment")
    }
  }

  # Get feature columns
  if (is.null(feature_cols)) {
    if (exists("avail")) {
      feature_cols <- avail
    } else if (exists("OPTIMAL_INDICES")) {
      feature_cols <- OPTIMAL_INDICES
    } else {
      stop("No feature_cols provided and avail/OPTIMAL_INDICES not found")
    }
  }

  # Get temporal budget
  temporal_budget <- if (exists("TEMPORAL_BUDGET")) TEMPORAL_BUDGET else 37

  if (is.null(templates) || length(templates) == 0) {
    message("Templates are empty or NULL")
    return(NULL)
  }

  # Create output directory
  if (save_png && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Build long dataframe
  rows <- list()
  for (veg in names(templates)) {
    t_mat <- templates[[veg]]$T  # This is a matrix: pentads x features

    if (is.null(t_mat)) next

    n_pentads <- nrow(t_mat)
    n_features <- ncol(t_mat)

    for (k in seq_len(min(n_features, length(feature_cols)))) {
      df_tmp <- data.frame(
        pentad = seq_len(n_pentads),
        value = t_mat[, k],
        Veg = veg,
        index = feature_cols[k],
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1]] <- df_tmp
    }
  }

  if (length(rows) == 0) {
    message("No valid data found in templates")
    return(NULL)
  }

  proto_df <- do.call(rbind, rows)

  # Get unique vegetation types and indices
  veg_types <- unique(proto_df$Veg)
  all_indices <- unique(proto_df$index)

  # Create color palette for indices
  n_colors <- length(all_indices)
  if (n_colors <= 12) {
    index_colors <- brewer.pal(n = max(3, n_colors), name = "Set3")[seq_len(n_colors)]
  } else {
    index_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_colors)
  }
  names(index_colors) <- all_indices

  # Generate plots
  plots <- list()

  for (veg in veg_types) {
    df_veg <- proto_df[proto_df$Veg == veg, , drop = FALSE]

    p <- ggplot(df_veg, aes(x = pentad, y = value, color = index)) +
      geom_line(size = 1, alpha = 0.9) +
      scale_color_manual(values = index_colors, name = "Index/Band") +
      labs(
        title = sprintf("Endmember: %s", tools::toTitleCase(veg)),
        subtitle = "Raw template spectral signature across pentads",
        x = "Pentad (10-day period)",
        y = "Value (standardized)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
        legend.position = "right",
        legend.key.size = unit(0.8, "lines"),
        panel.grid.minor = element_blank()
      ) +
      guides(color = guide_legend(ncol = 1))

    plots[[veg]] <- p

    if (save_png) {
      fn <- file.path(out_dir, sprintf("%s_%s.png", prefix, tolower(veg)))
      ggsave(filename = fn, plot = p, width = 12, height = 6, dpi = dpi)
      message(sprintf("Saved: %s", fn))
    }
  }

  # Combined faceted plot
  p_combined <- ggplot(proto_df, aes(x = pentad, y = value, color = index)) +
    geom_line(size = 0.8, alpha = 0.85) +
    facet_wrap(~ Veg, scales = "free_y", ncol = 2) +
    scale_color_manual(values = index_colors, name = "Index/Band") +
    labs(
      title = "Raw Template Spectral Signatures by Vegetation Type",
      x = "Pentad (10-day period)",
      y = "Value (standardized)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "bottom",
      legend.key.size = unit(0.6, "lines"),
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(nrow = 3))

  plots[["combined"]] <- p_combined

  if (save_png) {
    fn <- file.path(out_dir, sprintf("%s_combined.png", prefix))
    ggsave(filename = fn, plot = p_combined, width = 14, height = 10, dpi = dpi)
    message(sprintf("Saved: %s", fn))
  }

  invisible(plots)
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

message("=== Generating Endmember Spectral Visualizations ===\n")

# Check what data structures are available
if (exists("mesma_lib") && !is.null(mesma_lib)) {
  message("Using optimized MESMA library for visualization...")
  plots <- plot_endmember_spectra_by_index(
    lib = mesma_lib,
    out_dir = "endmember_spectra_plots",
    save_png = TRUE
  )
} else if (exists("OPTIMIZED_LIBRARY") && !is.null(OPTIMIZED_LIBRARY)) {
  message("Using OPTIMIZED_LIBRARY for visualization...")
  plots <- plot_endmember_spectra_by_index(
    lib = OPTIMIZED_LIBRARY,
    out_dir = "endmember_spectra_plots",
    save_png = TRUE
  )
} else if (exists("raw_lib_templates") && !is.null(raw_lib_templates)) {
  message("Using raw_lib_templates for visualization...")
  plots <- plot_raw_templates_by_index(
    templates = raw_lib_templates,
    out_dir = "endmember_spectra_plots",
    save_png = TRUE
  )
} else {
  message("No library data available. Please run fit_veg_mixture_mesma.R first.")
  message("Then call one of the plotting functions:")
  message("  - plot_endmember_spectra_by_index(mesma_lib)")
  message("  - plot_raw_templates_by_index(raw_lib_templates)")
}

message("\n=== Done ===")
