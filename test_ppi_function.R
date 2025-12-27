# test_ppi_function.R

# Load necessary libraries (assuming they are installed in the R environment)
library(dplyr)
library(lubridate)

# Source ppi_helpers.R to get auto_add_ppi_columns and its dependencies
if (file.exists("ppi_helpers.R")) {
  source("ppi_helpers.R")
  cat("Sourced ppi_helpers.R successfully.\n")
} else {
  stop("ppi_helpers.R not found in the current directory. Please ensure it's accessible.")
}

# --- Simulate an inference dataset (df_inf) ---
# It should contain at least date, lat, nir, red for PPI calculation
dummy_df_inf <- data.frame(
  location_id = rep(c("loc_A", "loc_B"), each = 3),
  date = as.Date(c("2023-06-15", "2023-07-01", "2023-08-10", "2023-06-20", "2023-07-05", "2023-08-15")),
  lat = rep(c(40.0, 41.0), each = 3),
  lon = rep(c(-105.0, -104.0), each = 3),
  nir = c(0.4, 0.5, 0.6, 0.35, 0.45, 0.55),
  red = c(0.1, 0.12, 0.15, 0.08, 0.1, 0.13),
  swir1 = c(0.2, 0.22, 0.25, 0.18, 0.2, 0.23),
  stringsAsFactors = FALSE
)

# Ensure phenological columns are added as fit_veg_mixture_mesma.R would
dummy_df_inf$date <- as.Date(dummy_df_inf$date)
dummy_df_inf$year <- lubridate::year(dummy_df_inf$date)
dummy_df_inf$pheno_year <- ifelse(lubridate::month(dummy_df_inf$date) >= 3, lubridate::year(dummy_df_inf$date), lubridate::year(dummy_df_inf$date) - 1)
dummy_df_inf$month <- lubridate::month(dummy_df_inf$date)
dummy_df_inf$doy <- lubridate::yday(dummy_df_inf$date)

cat("Simulated dummy_df_inf:\n")
print(head(dummy_df_inf))
cat("Columns in dummy_df_inf:", paste(names(dummy_df_inf), collapse = ", "), "\n\n")

# --- Simulate GLOBAL_TRAINING_DVI_SOIL ---
# Test scenario 1: GLOBAL_TRAINING_DVI_SOIL is available
GLOBAL_TRAINING_DVI_SOIL <- 0.05 # A typical DVI soil value
assign("GLOBAL_TRAINING_DVI_SOIL", GLOBAL_TRAINING_DVI_SOIL, envir = globalenv())
cat(sprintf("Simulated GLOBAL_TRAINING_DVI_SOIL: %.4f\n\n", GLOBAL_TRAINING_DVI_SOIL))

# --- Test auto_add_ppi_columns ---
cat("--- Testing auto_add_ppi_columns with GLOBAL_TRAINING_DVI_SOIL ---\n")
ppi_test_res_with_global <- auto_add_ppi_columns(dummy_df_inf, dvi_soil = get("GLOBAL_TRAINING_DVI_SOIL", envir = globalenv()))

if (!is.null(ppi_test_res_with_global) && isTRUE(ppi_test_res_with_global$added)) {
  cat("SUCCESS: PPI added with GLOBAL_TRAINING_DVI_SOIL. Reason:", ppi_test_res_with_global$reason, "\n")
  cat("Head of result df:\n")
  print(head(ppi_test_res_with_global$df))
  cat("Columns in result df:", paste(names(ppi_test_res_with_global$df), collapse = ", "), "\n")
  if ("PPI" %in% names(ppi_test_res_with_global$df)) {
    cat("PPI values (sample):", paste(head(ppi_test_res_with_global$df$PPI), collapse = ", "), "\n")
  }
} else {
  cat("FAILURE: PPI not added with GLOBAL_TRAINING_DVI_SOIL.\n")
  cat("Reason:", ppi_test_res_with_global$reason, "\n")
}

cat("\n--- Testing auto_add_ppi_columns WITHOUT providing dvi_soil explicitly (should try to calculate or fail) ---\n")
# Remove GLOBAL_TRAINING_DVI_SOIL from globalenv for this test
rm(GLOBAL_TRAINING_DVI_SOIL, envir = globalenv())
ppi_test_res_no_global <- auto_add_ppi_columns(dummy_df_inf) # This should try to calculate from barren or fail

if (!is.null(ppi_test_res_no_global) && isTRUE(ppi_test_res_no_global$added)) {
  cat("SUCCESS: PPI added without explicit dvi_soil. Reason:", ppi_test_res_no_global$reason, "\n")
  cat("Head of result df:\n")
  print(head(ppi_test_res_no_global$df))
  cat("Columns in result df:", paste(names(ppi_test_res_no_global$df), collapse = ", "), "\n")
  if ("PPI" %in% names(ppi_test_res_no_global$df)) {
    cat("PPI values (sample):", paste(head(ppi_test_res_no_global$df$PPI), collapse = ", "), "\n")
  }
} else {
  cat("FAILURE: PPI not added without explicit dvi_soil.\n")
  cat("Reason:", ppi_test_res_no_global$reason, "\n")
  cat("This is expected if dummy_df_inf has no 'Veg' == 'barren' rows.\n")
}

cat("\n--- Test script finished ---\n")
