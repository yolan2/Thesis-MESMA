#!/usr/bin/env Rscript

# Small test for unmix_stage1_compressed debug

source("fit_veg_mixture_mesma.R")

# Check if COMPRESSED_STAGE1_LIB exists
if (exists("COMPRESSED_STAGE1_LIB") && !is.null(COMPRESSED_STAGE1_LIB)) {
  cat("COMPRESSED_STAGE1_LIB exists\n")
  cat("barren:", paste(head(COMPRESSED_STAGE1_LIB$barren, 5), collapse=", "), "\n")
  cat("vegetation:", paste(head(COMPRESSED_STAGE1_LIB$vegetation, 5), collapse=", "), "\n")

  # Create a dummy y
  y <- COMPRESSED_STAGE1_LIB$barren + 0.1  # slightly different
  cat("y:", paste(head(y, 5), collapse=", "), "\n")

  # Call the function
  result <- unmix_stage1_compressed(y, COMPRESSED_STAGE1_LIB)
  cat("Result:", result, "\n")
} else {
  cat("COMPRESSED_STAGE1_LIB not available\n")
}