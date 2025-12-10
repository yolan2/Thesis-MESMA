#!/usr/bin/env Rscript

# Test script for newly added geometric MESMA functions
# Extract and test cos_angle, geometric_select_best_partner, solve_fclsu, unmix_stage2_geometric, geometric_stage1_with_selection

# Source the main script to load functions
source("fit_veg_mixture_mesma.R")

cat("Testing cos_angle function...\n")
try({
  # Test cos_angle with simple vectors
  v1 <- c(1, 0, 0)
  v2 <- c(0, 1, 0)
  angle <- cos_angle(v1, v2)
  expected <- 0  # 90 degrees, cos(90°) = 0
  cat(sprintf("cos_angle result: %.6f (expected: %.6f)\n", angle, expected))
  if (abs(angle - expected) > 1e-6) stop("cos_angle test failed")
  cat("✓ cos_angle test PASSED\n")
}, error = function(e) {
  cat(sprintf("✗ cos_angle test FAILED: %s\n", e$message))
})

cat("\nTesting geometric_select_best_partner function...\n")
try({
  # Create test data
  y <- c(0.5, 0.4, 0.3)
  m1 <- c(0.6, 0.5, 0.4)
  M2_list <- list(
    list(vec = c(0.4, 0.3, 0.2), id = "partner1"),
    list(vec = c(0.7, 0.6, 0.5), id = "partner2")
  )

  result <- geometric_select_best_partner(y, m1, M2_list)
  cat(sprintf("Selected partner: %s\n", result$m2$id))
  cat("✓ geometric_select_best_partner test PASSED\n")
}, error = function(e) {
  cat(sprintf("✗ geometric_select_best_partner test FAILED: %s\n", e$message))
})

cat("\nTesting solve_fclsu function...\n")
try({
  # Create test data for constrained least squares
  y <- c(0.5, 0.4, 0.3)
  M <- matrix(c(0.6, 0.4, 0.5, 0.3, 0.4, 0.2), nrow = 3, ncol = 2)

  result <- solve_fclsu(y, M)
  cat(sprintf("Weights: %s\n", paste(round(result$weights, 4), collapse = ", ")))
  cat(sprintf("Residual: %.6f\n", result$residual))
  cat("✓ solve_fclsu test PASSED\n")
}, error = function(e) {
  cat(sprintf("✗ solve_fclsu test FAILED: %s\n", e$message))
})

cat("\nTesting unmix_stage2_geometric function...\n")
try({
  # Create test data
  y <- c(0.5, 0.4, 0.3)
  veg_libs <- list(
    veg1 = list(list(vec = c(0.6, 0.5, 0.4), id = "v1")),
    veg2 = list(list(vec = c(0.4, 0.3, 0.2), id = "v2"))
  )

  result <- unmix_stage2_geometric(y, veg_libs)
  cat(sprintf("Fractions: %s\n", paste(sprintf("%s=%.4f", names(result$fractions), result$fractions), collapse = ", ")))
  cat("✓ unmix_stage2_geometric test PASSED\n")
}, error = function(e) {
  cat(sprintf("✗ unmix_stage2_geometric test FAILED: %s\n", e$message))
})

cat("\nTesting geometric_stage1_with_selection function...\n")
try({
  # Create test data
  y <- c(0.4, 0.35, 0.3)
  barren <- c(0.3, 0.3, 0.3)
  veg_agg <- c(0.5, 0.45, 0.4)

  result <- geometric_stage1_with_selection(y, barren, veg_agg)
  cat(sprintf("Fractions: barren=%.4f, veg=%.4f\n", result$barren, result$veg))
  cat("✓ geometric_stage1_with_selection test PASSED\n")
}, error = function(e) {
  cat(sprintf("✗ geometric_stage1_with_selection test FAILED: %s\n", e$message))
})

cat("\nAll tests completed!\n")

# New: Per-year meta-variant visualization
cat("\nGenerating per-year meta-variant visualization...\n")

# Function to analyze endmember stability per year
analyze_endmember_stability_per_year <- function(lib_df, veg_types, raw_lib_templates, avail_idx, years, B = 50, seed = 123) {
  results_per_year <- list()
  
  for (yr in years) {
    cat(sprintf("\n--- Analyzing year %d ---\n", yr))
    
    # Filter data for this year
    lib_df_year <- lib_df[lib_df$year == yr, , drop = FALSE]
    if (nrow(lib_df_year) == 0) {
      cat(sprintf("No data for year %d, skipping\n", yr))
      next
    }
    
    # Get unique location-year pairs for this year
    loc_years_year <- data.frame(location_id = unique(lib_df_year$location_id), year = yr)
    
    # Run bootstrap for this year (using raw templates instead of PCA)
    boot_result <- spatial_bootstrap_library(lib_df_year, veg_types, raw_lib_templates, avail_idx, B = B, seed = seed + yr)
    
    if (is.null(boot_result) || length(boot_result$boot_libs) == 0) {
      cat(sprintf("Bootstrap failed for year %d\n", yr))
      next
    }
    
    # Analyze stability
    stability <- analyze_endmember_stability(boot_result, reference_lib = NULL)
    
    results_per_year[[as.character(yr)]] <- list(
      year = yr,
      stability_results = stability,
      n_boot = boot_result$n_boot
    )
  }
  
  results_per_year
}

# Function to plot meta-variants per year
plot_meta_variants_per_year <- function(results_per_year, output_dir = "phenology_results") {
  if (!requireNamespace('ggplot2', quietly = TRUE)) {
    cat("ggplot2 not available, skipping plot\n")
    return(NULL)
  }
  library(ggplot2)
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Collect data for plotting
  plot_data <- list()
  
  for (yr_str in names(results_per_year)) {
    res <- results_per_year[[yr_str]]
    stability <- res$stability_results
    
    for (veg in names(stability)) {
      metas <- stability[[veg]]$meta_variants
      if (length(metas) == 0) next
      
      for (k in seq_along(metas)) {
        mv <- metas[[k]]
        centroid <- mv$centroid
        feature_idx <- seq_along(centroid)
        
        plot_data[[length(plot_data) + 1]] <- data.frame(
          year = res$year,
          veg = veg,
          meta_id = sprintf("%s_meta%d", veg, k),
          feature = feature_idx,
          value = centroid,
          n_members = mv$n_members,
          cv = mv$coefficient_of_variation
        )
      }
    }
  }
  
  if (length(plot_data) == 0) {
    cat("No meta-variant data to plot\n")
    return(NULL)
  }
  
  all_data <- do.call(rbind, plot_data)
  
  # Create faceted plot
  p <- ggplot(all_data, aes(x = feature, y = value, color = meta_id, group = meta_id)) +
    geom_line(size = 1) +
    facet_grid(veg ~ year, scales = "free_y") +
    labs(title = "Meta-Variants per Year",
         x = "Feature Index",
         y = "Value",
         color = "Meta-Variant") +
    theme_minimal() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 10))
  
  # Save plot
  plot_file <- file.path(output_dir, "meta_variants_per_year.png")
  ggsave(plot_file, p, width = 12, height = 8, dpi = 300)
  cat(sprintf("Per-year meta-variant plot saved to: %s\n", plot_file))
  
  # Also create a summary table
  summary_data <- all_data %>%
    group_by(year, veg, meta_id) %>%
    summarize(
      n_members = first(n_members),
      cv = first(cv),
      .groups = "drop"
    )
  
  summary_file <- file.path(output_dir, "meta_variants_per_year_summary.csv")
  write.csv(summary_data, summary_file, row.names = FALSE)
  cat(sprintf("Summary table saved to: %s\n", summary_file))
  
  list(plot = p, data = all_data, summary = summary_data)
}

# Run the analysis (assuming lib_df, veg_types, etc. are available from sourcing the main script)
if (exists("lib_df") && exists("vegs") && exists("raw_lib_templates") && exists("avail") && !isTRUE(TESTING_MODE)) {
  years <- sort(unique(lib_df$year))
  cat(sprintf("Analyzing stability for years: %s\n", paste(years, collapse = ", ")))
  
  per_year_results <- analyze_endmember_stability_per_year(
    lib_df = lib_df,
    veg_types = vegs,
    raw_lib_templates = raw_lib_templates,
    avail_idx = avail,
    years = years,
    B = 20  # Smaller B for testing
  )
  
  if (length(per_year_results) > 0) {
    plot_result <- plot_meta_variants_per_year(per_year_results)
    if (!is.null(plot_result)) {
      cat("Per-year visualization completed successfully!\n")
    }
  } else {
    cat("No per-year results to visualize\n")
  }
} else {
  if (isTRUE(TESTING_MODE)) {
    cat("Skipping per-year stability analysis in testing mode.\n")
  } else {
    cat("Required data not available (lib_df, vegs, raw_lib_templates, avail). Run the main script first.\n")
  }
}