# Two-Stage MESMA Implementation

## Overview

The MESMA analysis now uses a two-stage approach to accurately estimate vegetation composition for pixels with mixed barren/vegetation fractions:

### Stage 1: Vegetated Fraction Estimation
- **Endmembers**: Pure barren (`no soil == 0`) and pure vegetation (`no soil == 1`)
- **Goal**: Unmix the vegetated fraction for each pixel using spectral unmixing
- **Method**: Constrained least squares solving `pixel = α·barren + (1-α)·vegetation`
- **Output**: Vegetated fraction for each location-year

### Stage 2: Vegetation Type Decomposition
- **Input**: Vegetated fraction from Stage 1
- **Goal**: Decompose the vegetated portion into vegetation types (e.g., populus, tamarix, phragmites)
- **Method**: Standard MESMA with vegetation endmembers
- **Output**: Scaled coefficients for each vegetation type

## How It Works

### Training Data Requirements

Your training data should include three classes in the `no soil` field:

1. **Barren endmembers** (`no soil = 0`): Pure bare soil/barren locations
2. **Pure vegetation endmembers** (`no soil = 1`): Fully vegetated locations  
3. **Mixed pixels** (`no soil ∈ (0,1)`): Intermediate cases to be unmixed

The `Veg` field should indicate:
- `"barren"` for barren locations
- Vegetation type (e.g., `"grassland"`, `"populus"`) for vegetated locations

### Stage 1 Library Construction

The function `build_barren_veg_library()` creates two endmember libraries:

```r
stage1_lib <- list(
  barren = list(NDVI = list(mu = ..., mv = ...)),      # Barren temporal signature
  vegetation = list(NDVI = list(mu = ..., mv = ...))   # Pure veg temporal signature
)
```

- Aggregates spectral values by DOY for barren and pure vegetation training samples
- Uses median aggregation to reduce outlier influence
- Interpolates missing DOYs for complete temporal coverage (365 days)

### Stage 1 Unmixing

The function `unmix_vegetated_fraction()` estimates vegetated fraction:

```r
# For each observation at DOY d:
# pixel(d) = α·barren(d) + (1-α)·vegetation(d)
# 
# Solve for α (barren fraction):
# α = (pixel - veg) · (barren - veg) / ||barren - veg||²
#
# vegetated_fraction = 1 - α
```

- Performs per-pixel unmixing for each observation
- Returns median vegetated fraction across all observations for the location-year
- Robust to temporal variability and missing data

### Stage 2 Decomposition

After Stage 1 estimates the vegetated fraction (e.g., 0.65), Stage 2:

1. Runs standard MESMA to decompose vegetation types
2. Scales each vegetation coefficient by the vegetated fraction
3. Adds a barren component with coefficient = (1 - vegetated_fraction)

**Example output:**
```r
# Stage 1: vegetated_fraction = 0.65, barren_fraction = 0.35
# Stage 2: grassland = 0.7, shrub = 0.3 (within vegetated portion)
#
# Final output:
coef_df <- data.frame(
  Veg = c("grassland", "shrub", "barren"),
  coef = c(0.455, 0.195, 0.35)  # 0.7*0.65, 0.3*0.65, 0.35
)
```

## Fallback Behavior

**No fallback to metadata-based estimation.** If Stage 1 library cannot be built (insufficient `no soil == 0` or `no soil == 1` samples), the script stops with an error. This ensures all results are based on spectral unmixing rather than metadata assumptions. Stage 1 unmixing must succeed for each location-year, otherwise the location-year is skipped.

### `build_barren_veg_library(df_local, avail_idx, min_samples = 5)`
Builds stage 1 endmember library from training data.

**Parameters:**
- `df_local`: Training dataframe with `no soil`, `Veg`, `doy`, and spectral indices
- `avail_idx`: Vector of spectral index names (e.g., `c("NDVI", "MSAVI")`)
- `min_samples`: Minimum samples required for each endmember class

**Returns:** List with `barren` and `vegetation` endmember libraries, or `NULL` if insufficient data

### `unmix_vegetated_fraction(dly_local, stage1_lib, avail_idx)`
Performs stage 1 unmixing to estimate vegetated fraction.

**Parameters:**
- `dly_local`: Observations for a location-year
- `stage1_lib`: Stage 1 library from `build_barren_veg_library()`
- `avail_idx`: Vector of spectral index names

**Returns:** Estimated vegetated fraction (0 to 1), or `NA` if unmixing fails

## Testing

Run the comprehensive test suite:

```powershell
Rscript tests/test_two_stage_mesma.R
```

The test:
1. Creates synthetic training data with known vegetated fractions
2. Builds stage 1 library from pure endmembers
3. Unmixes mixed pixels and validates against ground truth
4. Verifies pure class recovery (barren → 0, vegetation → 1)

**Expected output:**
```
✓ ALL TESTS PASSED
Two-stage MESMA implementation is working correctly
```

## Benefits

1. **Spectral unmixing**: Uses actual spectral signatures instead of hardcoded metadata
2. **Handles intermediate cases**: Correctly processes `no soil ∈ (0,1)` pixels
3. **Temporal robustness**: Leverages full seasonal patterns for both endmembers
4. **Strict spectral-only**: Requires sufficient stage-1 training data to proceed
5. **Physical constraints**: Ensures fractions sum to 1 and remain in [0,1]

## Configuration

Relevant configuration parameters in `fit_veg_mixture_mesma.R`:

- `ALLOWED_VEG`: Vegetation types for stage 2 decomposition (automatically includes "barren")
- Filtering now preserves barren rows: `tolower(df$Veg) %in% ALLOWED_VEG | tolower(df$Veg) == "barren"`

## Uncertainty Estimation (Bootstrap)

This implementation now supports uncertainty estimation that propagates both Stage 1 (barren vs vegetation) and Stage 2 (vegetation decomposition) through a nested bootstrap when compressed Stage 1 libraries are available.

- By default, the script performs a Stage 2-only bootstrap for vegetation proportions (fixed Stage 1), because we only resampled residuals from the Stage 2 fit and recomputed weights.
- When `COMPRESSED_STAGE1_LIB` exists, the script runs a nested two-stage bootstrap that resamples Stage 2 residuals in compressed space and recomputes both Stage 1 vegetated fractions and Stage 2 weights for each bootstrap replicate (`nested_two_stage_bootstrap`).
- When `COMPRESSED_STAGE1_LIB` exists, the script runs a nested two-stage bootstrap that resamples Stage 2 residuals in compressed space and recomputes both Stage 1 vegetated fractions and Stage 2 weights for each bootstrap replicate (`nested_two_stage_bootstrap`).
  - Important: If the nested bootstrap fails for any reason, there is NO fallback to the old Stage 2-only bootstrap. `uncertainty` will be set to `NULL` and a warning will be logged. This avoids silently mixing two different bootstrap approaches that represent different statistical assumptions.
- The nested bootstrap returns CIs for final vegetation coefficients (already scaled by the vegetated fraction), and also supplies a CI for the vegetated fraction itself (`veg_frac_ci`).

## L1 Regularization (Lasso) option

You can optionally enable Lasso (L1) regularization for Stage 2 to bias toward sparser solutions (fewer active vegetation types):

```r
# Run unmixing with Lasso penalty - this performs L1 regularization on vegetation fractions
mesma_result <- unmix_stage2_compressed(
  y = y_vec,
  grid_type = comp_res$grid_type,
  mesma_lib = mesma_lib[veg_kept],
  topK = TOPK_VARIANTS,
  penalty = "lasso",
  lasso_lambda = NULL # if NULL, glmnet::cv.glmnet selects lambda via CV
)
```

Note: Lasso requires the `glmnet` package to be installed. If it is not available, the solver falls back to the ridge/simplex approach.

### Lambda tuning and visualization

You can inspect the lambda tradeoff between fit quality and solution sparsity:

```r
library(glmnet)
lambda_grid <- 10^seq(-3, 1, length.out = 20)
results <- lapply(lambda_grid, function(lam) {
  fit <- fit_lasso_mesma(y_norm, E_pool, lambda = lam, sum_to_one = TRUE)
  c(rmse = fit$rmse, n_active = sum(fit$w > 0.01))
})
df <- do.call(rbind, results)
df <- data.frame(lambda = lambda_grid, rmse = df[, "rmse"], n_active = df[, "n_active"]) 
library(ggplot2)
ggplot(df, aes(x = n_active, y = rmse)) + geom_point() + geom_line() + labs(x = "Number of Active Vegetation Types", y = "RMSE", title = "Pareto Front: Complexity vs Fit Quality")
```

### Hybrid Monte Carlo Propagation (recommended)

When both bootstraps are available (Stage 1's `veg_frac_boot` and Stage 2's `w_boot`), the script now performs a Monte Carlo propagation step that draws random samples from both bootstrap distributions and multiplies them to get a sample of the final coefficients. This yields percentile-based CI estimates for the final scaled vegetation coefficients that incorporate uncertainty from both stages and their covariance.

- Stage 1 bootstrap: `stage1_block_bootstrap(y, compressed_stage1_lib, B)` returns `veg_frac_boot` and `veg_frac_ci`.
- Stage 2 bootstrap: `gls_block_bootstrap(y, ...)` returns `w_boot` and `w_boot_raw` matrices along with `coef_ci` and other diagnostics.
- Propagation: `propagate_uncertainty_monte_carlo(veg_frac_boot, w_boot, n_samples = 1000)` draws random pairs from both distributions and computes combined CIs from the empirical distribution of the products.

Behavior:
- If Stage 1 compressed library is available the script will attempt to compute `stage1_block_bootstrap` and `gls_block_bootstrap` and then perform MC propagation. If the nested bootstrap (which recomputes both stages) fails, there is NO fallback to the old single-stage bootstrap. In contrast, if `COMPRESSED_STAGE1_LIB` is not present the script keeps the Stage 2-only bootstrap behavior.

You may toggle variant-switching during Stage 2 bootstrapping with the global variables `VARIANT_SWITCH_BOOTSTRAP` (default `FALSE`) and `VARIANT_SWITCH_RE_CENTER` (default `TRUE`).

## Example Workflow

```r
# 1. Prepare training data with no_soil field:
#    - 0 for barren locations
#    - 1 for pure vegetation
#    - intermediate values for mixed pixels

# 2. Run fit_veg_mixture_mesma.R
#    Stage 1 library will be built automatically from no_soil==0 and no_soil==1

# 3. For each location-year:
#    a. Stage 1 unmixes vegetated fraction
#    b. Stage 2 decomposes vegetation types
#    c. Final output includes scaled vegetation + barren fraction

# 4. Outputs include:
#    - coef_df: Vegetation coefficients + barren fraction
#    - diagnostics: vegetated_fraction, barren_fraction per location-year
```

## References

- **MESMA**: Roberts et al. (1998), Multiple Endmember Spectral Mixture Analysis
- **Constrained unmixing**: Heinz & Chang (2001), Fully constrained least squares linear spectral mixture analysis
