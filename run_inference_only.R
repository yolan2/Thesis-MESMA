# Load the cached model and run inference
cat("Loading MESMA model cache and running inference...\n")

# Source the main script to get function definitions
source("fit_veg_mixture_mesma.R")

# Now manually trigger the main processing
cat("\n=== Running inference manually ===\n")

# Prefer process_inference(); fall back to main_processing_block for legacy behavior
if (exists("process_inference") && is.function(process_inference)) {
  result <- tryCatch({
    ri <- process_inference()
    if (!is.null(ri$all_coefs) && nrow(ri$all_coefs) > 0) {
      cat("\n[SUCCESS] process_inference completed with coefficients\n")
    } else {
      cat("\n[NOTICE] process_inference completed but produced no coefficients\n")
    }
    list(success = TRUE)
  }, error = function(e) {
    cat(sprintf("\nERROR: %s\n", e$message))
    list(success = FALSE, error = e$message)
  })

  if (isTRUE(result$success)) {
    cat("\n[SUCCESS] Inference completed successfully!\n")
  } else {
    cat("\n[FAILED] Inference failed\n")
  }
} else if (exists("main_processing_block")) {
  result <- tryCatch({
    main_processing_block()
    list(success = TRUE)
  }, error = function(e) {
    cat(sprintf("\nERROR: %s\n", e$message))
    list(success = FALSE, error = e$message)
  })

  if (isTRUE(result$success)) {
    cat("\n[SUCCESS] Inference completed successfully!\n")
  } else {
    cat("\n[FAILED] Inference failed\n")
  }
} else {
  cat("\n[ERROR] No inference function found (process_inference or main_processing_block)\n")
  cat("Available functions:\n")
  print(ls(pattern = "process|main"))
}
