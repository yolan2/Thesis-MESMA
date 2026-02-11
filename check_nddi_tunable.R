# Smoke-test: verify NDDI_DUST_THRESHOLD tunability via MESMA_NDDI_THRESHOLD
suppressWarnings(suppressMessages({
  if (file.exists("mesma_helpers.R")) source("mesma_helpers.R")
}))

make_df <- function(nddis) {
  # construct red/nir values that produce given NDDI values (solve for red given nir=1)
  nir <- rep(1, length(nddis))
  red <- (1 + nddis) / (1 - nddis) * nir  # from (red - nir)/(red + nir) = nddi
  data.frame(red = red, nir = nir)
}

cat(sprintf("[TEST] Default threshold (env unset) -> NDDI_DUST_THRESHOLD= %s\n", .nddi_thresh_fmt()))
# baseline: create three values relative to the runtime threshold (below, just-below, just-above)
th <- as.numeric(NDDI_DUST_THRESHOLD)
vals <- c(th - 0.02, th - 0.001, th + 0.02)
df <- make_df(vals)
if (exists("NDDI_DUST_THRESHOLD")) {
  kept_default <- sum(!((df$red - df$nir)/(df$red + df$nir) > NDDI_DUST_THRESHOLD))
  # Expect the two values strictly below the threshold to be kept
  if (kept_default != 2) {
    cat("[FAIL] Unexpected filtering with default threshold. Kept rows:", kept_default, "expected 2\n")
    quit(status = 2)
  }
} else {
  cat("[FAIL] NDDI_DUST_THRESHOLD not defined (mesma_helpers.R not sourced?)\n")
  quit(status = 3)
}

# Now override via env var and re-source helpers to pick up change
Sys.setenv(MESMA_NDDI_THRESHOLD = "0.10")
source("mesma_helpers.R")
cat(sprintf("[TEST] Overridden threshold -> NDDI_DUST_THRESHOLD= %s\n", .nddi_thresh_fmt()))
df2 <- make_df(c(0.05, 0.12, 0.2))
kept_low_thresh <- sum(!((df2$red - df2$nir)/(df2$red + df2$nir) > NDDI_DUST_THRESHOLD))
if (kept_low_thresh != 1) {
  cat("[FAIL] Unexpected filtering after overriding MESMA_NDDI_THRESHOLD=0.10; kept:", kept_low_thresh, "expected 1\n")
  quit(status = 4)
}

cat("[PASS] NDDI threshold tunability smoke-test passed.\n")
quit(status = 0)
