#!/usr/bin/env Rscript
# Quick check of raw inference data file

# Find the inference file path from the main script
INFERENCE_CSV <- NULL

# Try common paths
possible_paths <- c(
  "inference_data.csv",
  "data/inference_data.csv",
  "../data/inference_data.csv"
)

# Check if INFERENCE_CSV is defined in an env file or config
if (file.exists(".Rprofile")) {
  source(".Rprofile")
}

# Try to find any CSV files that might be the inference data
csv_files <- list.files(".", pattern = ".*inference.*\\.csv$", ignore.case = TRUE, full.names = TRUE, recursive = TRUE)
if (length(csv_files) > 0) {
  cat("Found possible inference files:\n")
  for (f in csv_files) {
    cat(sprintf("  %s (%.2f MB)\n", f, file.size(f) / 1024^2))
  }

  # Use the first one
  INFERENCE_CSV <- csv_files[1]
  cat(sprintf("\nUsing: %s\n", INFERENCE_CSV))

  # Read and check
  df_inf <- read.csv(INFERENCE_CSV, nrows = 1000)  # Read first 1000 rows for speed

  cat(sprintf("\nDimensions (first 1000 rows): %d rows x %d cols\n", nrow(df_inf), ncol(df_inf)))
  cat("\nColumn names:\n")
  print(colnames(df_inf))

  cat("\n=== CHECKING FOR SPECTRAL INDICES ===\n")
  required_indices <- c("NBR", "NDTI", "NDVI", "EVI", "SAVI", "MSAVI", "kNDVI",
                        "NDMI", "NDWI", "TCB", "TCG", "TCW", "PPI")

  present <- required_indices[required_indices %in% colnames(df_inf)]
  missing <- setdiff(required_indices, colnames(df_inf))

  cat(sprintf("\nPresent indices (%d): %s\n", length(present), paste(present, collapse=", ")))
  cat(sprintf("Missing indices (%d): %s\n", length(missing), paste(missing, collapse=", ")))

  if (length(present) > 0) {
    cat("\n=== COMPLETENESS OF PRESENT INDICES ===\n")
    for (idx in present) {
      n_finite <- sum(is.finite(df_inf[[idx]]))
      pct <- 100 * n_finite / nrow(df_inf)
      cat(sprintf("  %-8s: %6.1f%% complete (%d/%d)\n", idx, pct, n_finite, nrow(df_inf)))
    }
  }

  # Check for raw bands
  cat("\n=== CHECKING FOR RAW BANDS ===\n")
  raw_bands <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "NIR", "RED", "SWIR1", "SWIR2")
  present_bands <- raw_bands[raw_bands %in% colnames(df_inf)]
  cat(sprintf("Present bands (%d): %s\n", length(present_bands), paste(present_bands, collapse=", ")))

  if (length(present_bands) > 0) {
    cat("\n=== COMPLETENESS OF PRESENT BANDS ===\n")
    for (band in present_bands) {
      n_finite <- sum(is.finite(df_inf[[band]]))
      pct <- 100 * n_finite / nrow(df_inf)
      cat(sprintf("  %-8s: %6.1f%% complete (%d/%d)\n", band, pct, n_finite, nrow(df_inf)))
    }
  }

} else {
  cat("No inference CSV files found. Please provide the path to your inference data file.\n")
}
