#!/usr/bin/env Rscript
# Check the actual inference input file

INFERENCE_CSV <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/landsat_lower_inference.csv"

if (file.exists(INFERENCE_CSV)) {
  cat(sprintf("File found: %s\n", INFERENCE_CSV))
  cat(sprintf("File size: %.2f MB\n", file.size(INFERENCE_CSV) / 1024^2))

  # Read a sample
  df <- read.csv(INFERENCE_CSV, nrows = 500)
  cat(sprintf("\nDimensions (first 500 rows): %d rows x %d cols\n", nrow(df), ncol(df)))

  cat("\nColumn names:\n")
  print(colnames(df))

  cat("\n=== CHECKING FOR SPECTRAL INDICES ===\n")
  required_indices <- c("NBR", "NDTI", "NDVI", "EVI", "SAVI", "MSAVI", "kNDVI",
                        "NDMI", "NDWI", "TCB", "TCG", "TCW", "PPI")

  present <- intersect(required_indices, colnames(df))
  missing <- setdiff(required_indices, colnames(df))

  cat(sprintf("\nPresent indices (%d): %s\n", length(present), paste(present, collapse=", ")))
  cat(sprintf("Missing indices (%d): %s\n", length(missing), paste(missing, collapse=", ")))

  if (length(present) > 0) {
    cat("\n=== COMPLETENESS OF PRESENT INDICES ===\n")
    for (idx in present) {
      n_finite <- sum(is.finite(df[[idx]]))
      pct <- 100 * n_finite / nrow(df)
      cat(sprintf("  %-8s: %6.1f%% complete (%d/%d)\n", idx, pct, n_finite, nrow(df)))
    }
  }

  # Check for raw bands
  cat("\n=== CHECKING FOR RAW BANDS ===\n")
  raw_bands <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "NIR", "RED", "SWIR1", "SWIR2", "BLUE", "GREEN")
  present_bands <- intersect(raw_bands, colnames(df))
  cat(sprintf("Present bands (%d): %s\n", length(present_bands), paste(present_bands, collapse=", ")))

  if (length(present_bands) > 0) {
    cat("\n=== COMPLETENESS OF PRESENT BANDS ===\n")
    for (band in present_bands) {
      n_finite <- sum(is.finite(df[[band]]))
      pct <- 100 * n_finite / nrow(df)
      cat(sprintf("  %-8s: %6.1f%% complete (%d/%d)\n", band, pct, n_finite, nrow(df)))
    }
  }

  # Check for location_id and date
  cat("\n=== CHECKING KEY COLUMNS ===\n")
  if ("location_id" %in% colnames(df)) {
    cat(sprintf("location_id: %d unique values\n", length(unique(df$location_id))))
  } else {
    cat("location_id: MISSING\n")
  }

  date_cols <- c("date", "prediction_date", "reference_date")
  date_found <- FALSE
  for (dc in date_cols) {
    if (dc %in% colnames(df)) {
      cat(sprintf("%s: found (%d values)\n", dc, sum(!is.na(df[[dc]]))))
      date_found <- TRUE
    }
  }
  if (!date_found) {
    cat("No date column found (checked: date, prediction_date, reference_date)\n")
  }

} else {
  cat(sprintf("ERROR: File not found: %s\n", INFERENCE_CSV))
}
