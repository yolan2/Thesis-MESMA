
library(data.table)

# ==============================================================================
# 1. Define Endmember Spectra (Reflectance 0-1)
# ==============================================================================

# Healthy Green Vegetation
veg_spec <- c(
  blue  = 0.02,
  green = 0.06,
  red   = 0.03,
  nir   = 0.45,
  swir1 = 0.15,
  swir2 = 0.07
)

# Dry/Barren Soil
soil_spec <- c(
  blue  = 0.08,
  green = 0.12,
  red   = 0.18,
  nir   = 0.24,
  swir1 = 0.35,
  swir2 = 0.30
)

# ==============================================================================
# 2. Define Index Calculation Function
# ==============================================================================

calculate_indices <- function(df) {
  # Expects columns: blue, green, red, nir, swir1, swir2
  eps <- 1e-9
  
  df[, `:=`(
    DVI   = nir - red,
    OSAVI = (nir - red) / (nir + red + 0.16),
    MCARI = ((red - green) - 0.2*(red - blue)) * (red / (green + eps)),
    NIRv  = (nir * ((nir - red) / (nir + red + eps))) * 1.3,
    PSRI  = (red - blue) / (nir + eps),
    NBR   = (nir - swir2) / (nir + swir2 + eps),
    TCW   = (swir1 - swir2) / (swir1 + swir2 + eps),
    NDVI   = (nir - red) / (nir + red + eps),
    MSAVI2 = (2 * nir + 1 - sqrt(pmax(0, (2 * nir + 1)^2 - 8 * (nir - red)))) / 2,
    NDMI   = (nir - swir1) / (nir + swir1 + eps),
    TCB    = 0.3029 * blue + 0.2786 * green + 0.4733 * red + 0.5599 * nir + 0.508 * swir1 + 0.1872 * swir2,
    GVI    = -0.2941 * blue - 0.243 * green - 0.5424 * red + 0.7276 * nir + 0.0713 * swir1 - 0.1608 * swir2
  )]
  
  return(df)
}

# ==============================================================================
# 3. Perform Synthetic Mixing
# ==============================================================================

fractions <- seq(0, 1, by = 0.01)
mixtures <- data.table(fraction_veg = fractions)
bands <- names(veg_spec)
for (b in bands) {
  mixtures[[b]] <- fractions * veg_spec[b] + (1 - fractions) * soil_spec[b]
}
mixtures <- calculate_indices(mixtures)

# ==============================================================================
# 4. Test Transformations
# ==============================================================================

evaluate_linearity <- function(y, x) {
  y0 <- y[1]
  y1 <- y[length(y)]
  linear_pred <- y0 + (y1 - y0) * x
  residuals <- y - linear_pred
  max_dev <- max(abs(residuals))
  range_val <- abs(y1 - y0)
  ifelse(range_val > 1e-6, max_dev / range_val, 0)
}

indices <- c("NDVI", "MSAVI2", "NIRv", "PSRI", "MCARI", "TCW", "NBR", "NDMI")

cat("Original Linearity Scores (Normalized Max Dev):\n")
for (idx in indices) {
  score <- evaluate_linearity(mixtures[[idx]], mixtures$fraction_veg)
  cat(sprintf("%s: %.4f\n", idx, score))
}

cat("\nTesting Transformations...\n")

# Define transforms to test
transforms <- list(
  "Square" = function(x) x^2,
  "Sqrt" = function(x) sqrt(pmax(x, 0)),
  "SignedSqrt" = function(x) sign(x) * sqrt(abs(x)),
  "Log" = function(x) log(pmax(x, 1e-6)),
  "Log1p" = function(x) log(x + 1),
  "Exp" = function(x) exp(x),
  "Fisher" = function(x) 0.5 * log((1+x)/(1-x+1e-6)),
  "M_NDVI" = function(x) 1 - (1-x)^2,
  "Power3" = function(x) x^3,
  "Inv" = function(x) 1/(x+1e-6)
)

best_transforms <- list()

for (idx in indices) {
  vals <- mixtures[[idx]]
  base_score <- evaluate_linearity(vals, mixtures$fraction_veg)
  best_score <- base_score
  best_name <- "None"
  
  for (t_name in names(transforms)) {
    t_func <- transforms[[t_name]]
    # Handle potential domain errors
    t_vals <- tryCatch(t_func(vals), warning=function(w) NULL, error=function(e) NULL)
    if (is.null(t_vals) || any(!is.finite(t_vals))) next
    
    score <- evaluate_linearity(t_vals, mixtures$fraction_veg)
    if (score < best_score) {
      best_score <- score
      best_name <- t_name
    }
  }
  
  cat(sprintf("%s: Best=%s (%.4f -> %.4f)\n", idx, best_name, base_score, best_score))
  best_transforms[[idx]] <- best_name
}

cat("\nDetailed NDVI Analysis:\n")
idx <- "NDVI"
vals <- mixtures[[idx]]
base_score <- evaluate_linearity(vals, mixtures$fraction_veg)
cat(sprintf("Base: %.4f\n", base_score))

for (t_name in names(transforms)) {
    t_func <- transforms[[t_name]]
    t_vals <- tryCatch(t_func(vals), warning=function(w) NULL, error=function(e) NULL)
    if (is.null(t_vals) || any(!is.finite(t_vals))) {
        cat(sprintf("%s: Failed\n", t_name))
        next
    }
    score <- evaluate_linearity(t_vals, mixtures$fraction_veg)
    cat(sprintf("%s: %.4f\n", t_name, score))
}
