#!/usr/bin/env Rscript
# Quick test: ensure the main script contains the enforced check for ENABLE_UNCERTAINTY requiring COMPRESSED_STAGE1_LIB
file_text <- readLines("fit_veg_mixture_mesma.R")
msg <- 'ENABLE_UNCERTAINTY = TRUE requires a built COMPRESSED_STAGE1_LIB'
if (!any(grepl(msg, file_text, fixed = TRUE))) stop(sprintf("Expected enforcement message not found: %s", msg))
cat("Bootstrap enforcement check present in script (PASS)\n")
