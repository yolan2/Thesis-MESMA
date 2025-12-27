source("fit_veg_mixture_mesma.R")

# Mock data
# Simulate Z-score normalized PPI (negative values possible)
# And raw PPI (always positive)
df_mock <- data.frame(
  location_id = c("L1", "L1", "L2", "L2"),
  pheno_year = c(2020, 2020, 2020, 2020),
  PPI = c(-1.5, -1.5, 0.5, 0.5), # Normalized values, some negative
  PPI_raw = c(0.1, 0.1, 0.2, 0.2), # Raw values
  Veg = c("shrub", "shrub", "tree", "tree")
)

coefs_mock <- data.frame(
  location_id = c("L1", "L2"),
  pheno_year = c(2020, 2020),
  Veg = c("shrub", "tree"),
  coef = c(0.5, 0.8)
)

cat("Testing location_bootstrap_ppi with PPI_raw...\n")
res <- location_bootstrap_ppi(coefs_mock, df_mock, B = 10, seed = 123)

print(res)

if (!is.null(res) && all(res$global_coef > 0)) {
  cat("SUCCESS: Global coefs are positive using PPI_raw.\n")
} else {
  cat("FAILURE: Global coefs are not positive or NULL using PPI_raw.\n")
}

# Test fallback (simulating the bug)
df_mock_no_raw <- df_mock
df_mock_no_raw$PPI_raw <- NULL

cat("\nTesting location_bootstrap_ppi without PPI_raw (fallback to normalized PPI)...\n")
res_fallback <- location_bootstrap_ppi(coefs_mock, df_mock_no_raw, B = 10, seed = 123)
print(res_fallback)

if (!is.null(res_fallback) && any(res_fallback$global_coef < 0)) {
    cat("CONFIRMED: Fallback produces negative values (reproducing the reported issue).\n")
}