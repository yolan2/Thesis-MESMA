# Smoke-test: ensure NDSI code paths were removed (no functional NDSI usage)
# This performs a static check: no non-comment R source line should reference `NDSI`.
files <- list.files(pattern = "\\.R$", recursive = TRUE)
# Skip this test file (may contain intentional mention of NDSI) and skip generated/docs
files <- setdiff(files, c(basename("check_no_ndsi.R")))
bad <- character(0)
for (f in files) {
  txt <- readLines(f, warn = FALSE)
  for (i in seq_along(txt)) {
    line <- txt[i]
    if (grepl("\\bNDSI\\b", line, perl = TRUE) && !grepl("^\\s*#", line)) {
      bad <- c(bad, sprintf("%s:%d: %s", f, i, trimws(line)))
    }
  }
}
if (length(bad) > 0) {
  cat("[FAIL] Found non-comment references to NDSI:\n")
  cat(paste(bad, collapse = "\n"), "\n")
  quit(status = 2)
}
cat("[PASS] No active NDSI references found in R source (only comments allowed).\n")
quit(status = 0)
