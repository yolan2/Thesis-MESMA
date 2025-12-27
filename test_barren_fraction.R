source("fit_veg_mixture_mesma.R")

# Mock data
# L1: High PPI (full veg) -> Barren ~ 0
# L2: Low PPI (soil) -> Barren ~ 1
df_mock <- data.frame(
  location_id = c("L1", "L1", "L2", "L2"),
  pheno_year = c(2020, 2020, 2020, 2020),
  PPI_raw = c(0.4, 0.4, 0.0, 0.0), # L1=0.4 (max), L2=0.0
  Veg = c("shrub", "shrub", "tree", "tree")
)

coefs_mock <- data.frame(
  location_id = c("L1", "L2"),
  pheno_year = c(2020, 2020),
  Veg = c("shrub", "tree"),
  coef = c(1.0, 1.0)
)

PPI_FULL_VEG_COVER <<- 0.4 # Ensure global is set

cat("Testing location_bootstrap_ppi for Barren fraction...\n")
res <- location_bootstrap_ppi(coefs_mock, df_mock, B = 10, seed = 123)

print(res)

barren_res <- res[res$Veg == "Barren", ]
if (nrow(barren_res) > 0) {
    cat("SUCCESS: Barren class found in results.\n")
    print(barren_res)
} else {
    cat("FAILURE: Barren class NOT found.\n")
}

