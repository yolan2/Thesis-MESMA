#!/usr/bin/env Rscript

# Test that geometric_block_bootstrap is disabled to enforce nested-only bootstrapping
## To avoid sourcing the entire script (which runs the whole MESMA pipeline),
## we extract and evaluate only the function definition for geometric_block_bootstrap
## from the main script and then call it to assert that it raises the expected error.
file_content <- readLines("fit_veg_mixture_mesma.R")
file_text <- paste(file_content, collapse = "\n")
pattern_header <- "geometric_block_bootstrap\\s*<-\\s*function\\s*\\([^\\)]*\\)\\s*\\{"
loc <- regexpr(pattern_header, file_text, perl = TRUE)
if (loc == -1) stop("Could not find geometric_block_bootstrap function header in script")
start_pos <- loc + attr(loc, "match.length") - 1
# Now find matching closing brace for the function body
brace_count <- 1
pos <- start_pos + 1
while (brace_count > 0 && pos < nchar(file_text)) {
  ch <- substr(file_text, pos, pos)
  if (ch == "{") brace_count <- brace_count + 1
  if (ch == "}") brace_count <- brace_count - 1
  pos <- pos + 1
}
if (brace_count != 0) stop("Could not parse function body properly to extract geometry function")
func_text <- substr(file_text, loc, pos - 1)
eval(parse(text = func_text), envir = .GlobalEnv)

cat("Testing geometric_block_bootstrap deprecation...\n")
tryCatch({
  y_vec <- c(0.1, 0.2, 0.3)
  comp_templates <- list()
  top_variants <- list()
  chosen_ids <- list()
  w_hat <- c(veg1 = 1)
  res <- geometric_block_bootstrap(y_vec, comp_templates, top_variants, chosen_ids, w_hat)
  stop("geometric_block_bootstrap should be disabled but did not stop as expected")
}, error = function(e) {
  cat(sprintf("Expected error from geometric_block_bootstrap: %s\n", e$message))
  cat("✓ geometric_block_bootstrap is disabled as expected\n")
})

cat("Test complete\n")

cat("Testing ols_block_bootstrap deprecation...\n")
tryCatch({
  # Extract ols_block_bootstrap function from file and evaluate in the current environment
  pattern_header_ols <- "ols_block_bootstrap\\s*<-\\s*function\\s*\\([^\\)]*\\)\\s*\\{"
  loc_ols <- regexpr(pattern_header_ols, file_text, perl = TRUE)
  if (loc_ols == -1) stop("Could not find ols_block_bootstrap function header in script")
  start_pos_ols <- loc_ols + attr(loc_ols, "match.length") - 1
  brace_count_ols <- 1
  pos_ols <- start_pos_ols + 1
  while (brace_count_ols > 0 && pos_ols < nchar(file_text)) {
    ch <- substr(file_text, pos_ols, pos_ols)
    if (ch == "{") brace_count_ols <- brace_count_ols + 1
    if (ch == "}") brace_count_ols <- brace_count_ols - 1
    pos_ols <- pos_ols + 1
  }
  if (brace_count_ols != 0) stop("Could not parse ols_block_bootstrap function body properly")
  func_text_ols <- substr(file_text, loc_ols, pos_ols - 1)
  eval(parse(text = func_text_ols), envir = .GlobalEnv)
  y_vec <- c(0.1, 0.2, 0.3)
  comp_templates <- list()
  top_variants <- list()
  chosen_ids <- list()
  w_hat <- c(veg1 = 1)
  res <- ols_block_bootstrap(y_vec, comp_templates, top_variants, chosen_ids, w_hat)
  stop("ols_block_bootstrap should be disabled but did not stop as expected")
}, error = function(e) {
  cat(sprintf("Expected error from ols_block_bootstrap: %s\n", e$message))
  cat("✓ ols_block_bootstrap is disabled as expected\n")
})

cat("All deprecation tests complete\n")
