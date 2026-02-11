# Quick reproducibility smoke-test for MESMA seed handling
# Exit code 0 => PASS, non-zero => FAIL
suppressWarnings(suppressMessages({
  if (file.exists("mesma_helpers.R")) source("mesma_helpers.R")
}))

cat("[TEST] Running MESMA seed reproducibility test...\n")
set_mesma_seed(123)
A <- sample.int(1e6, 20)
set_mesma_seed(123)
B <- sample.int(1e6, 20)
if (!identical(A, B)) {
  cat("[FAIL] RNG outputs differ for same MESMA_SEED!\n")
  cat("A:", paste(A, collapse=","), "\n")
  cat("B:", paste(B, collapse=","), "\n")
  quit(status = 2)
}

# Check that derived seeds are deterministic
s1 <- get_mesma_seed(0)
s2 <- get_mesma_seed(0)
s3 <- get_mesma_seed(5)
if (!(identical(s1, s2) && !identical(s1, s3))) {
  cat("[FAIL] get_mesma_seed() behaviour unexpected\n")
  quit(status = 3)
}

cat("[PASS] MESMA_SEED reproducibility smoke-test passed.\n")
quit(status = 0)
