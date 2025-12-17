# Quick check if Stage 2 is working
library(openxlsx)

# Read the results
results_file <- "phenology_results/veg_mixture_fit/mesma_results.xlsx"
if (!file.exists(results_file)) {
  cat("Results file not found yet\n")
  quit()
}

# Read coefficient data
sheets <- getSheetNames(results_file)
if ("Coefficients" %in% sheets) {
  coefs <- read.xlsx(results_file, sheet = "Coefficients")
  
  cat("\n=== STAGE 2 EXECUTION CHECK ===\n")
  cat(sprintf("Total coefficient rows: %d\n", nrow(coefs)))
  
  if ("veg_type" %in% names(coefs)) {
    veg_types <- table(coefs$veg_type, useNA = "ifany")
    cat("\nVegetation types found:\n")
    print(veg_types)
    
    if (any(!is.na(coefs$veg_type) & coefs$veg_type != "")) {
      cat("\n✓ SUCCESS: Stage 2 IS running! Found vegetation type unmixing results.\n")
    } else {
      cat("\n✗ ISSUE: Stage 2 still not running - no vegetation types found.\n")
    }
  } else {
    cat("\n✗ ISSUE: 'veg_type' column not found in results.\n")
  }
  
  # Show sample
  cat("\nSample coefficients (first 10 rows):\n")
  print(head(coefs[, c("location_id", "pheno_year", "endmember", "veg_type", "coefficient")], 10))
} else {
  cat("Coefficients sheet not found\n")
}
