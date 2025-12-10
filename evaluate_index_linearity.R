# evaluate_index_linearity.R
# Script to evaluate the linearity of spectral indices through synthetic mixing
# of Vegetation and Soil endmembers, and between different Vegetation types.

library(data.table)
library(ggplot2)

# ==============================================================================
# 1. Define Endmember Spectra (Reflectance 0-1)
# ==============================================================================

# Approximate typical reflectance values
# Band order: Blue, Green, Red, NIR, SWIR1, SWIR2
# Based on typical spectral signatures

# Healthy Green Vegetation (baseline)
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

# Different Vegetation Types (variations of healthy vegetation)
# Based on typical spectral differences between species
veg_endmembers <- list(
  populus = c(  # Deciduous tree - high NIR, moderate red absorption
    blue  = 0.025,
    green = 0.065,
    red   = 0.035,
    nir   = 0.48,
    swir1 = 0.16,
    swir2 = 0.08
  ),
  tamarix = c(  # Shrub - slightly stressed, lower NIR
    blue  = 0.028,
    green = 0.058,
    red   = 0.038,
    nir   = 0.42,
    swir1 = 0.18,
    swir2 = 0.09
  ),
  phragmites = c(  # Wetland grass - high water content, different SWIR
    blue  = 0.022,
    green = 0.062,
    red   = 0.032,
    nir   = 0.46,
    swir1 = 0.12,
    swir2 = 0.06
  ),
  salicornia = c(  # Salt-tolerant - stressed, lower NIR, higher SWIR
    blue  = 0.030,
    green = 0.055,
    red   = 0.040,
    nir   = 0.38,
    swir1 = 0.22,
    swir2 = 0.12
  ),
  halocnemum = c(  # Succulent - very stressed, low NIR
    blue  = 0.035,
    green = 0.052,
    red   = 0.045,
    nir   = 0.35,
    swir1 = 0.25,
    swir2 = 0.15
  ),
  alhagi = c(  # Legume shrub - moderate health
    blue  = 0.026,
    green = 0.060,
    red   = 0.034,
    nir   = 0.44,
    swir1 = 0.17,
    swir2 = 0.08
  )
)

# ==============================================================================
# 2. Define Index Calculation Function
# ==============================================================================

calculate_indices <- function(df) {
  # Expects columns: blue, green, red, nir, swir1, swir2
  eps <- 1e-9
  
  df[, `:=`(
    # Original Set
    DVI   = nir - red,
    OSAVI = (nir - red) / (nir + red + 0.16),
    MCARI = ((red - green) - 0.2*(red - blue)) * (red / (green + eps)),
    #CRI   = (1/(green + eps)) - (1/(red + eps)),
    #PRI   = (green - red) / (green + red + eps),
    NIRv  = (nir * ((nir - red) / (nir + red + eps))) * 1.3, # Including the 1.3 scaling
    PSRI  = (red - blue) / (nir + eps),
    NBR   = (nir - swir2) / (nir + swir2 + eps),
    TCW   = (swir1 - swir2) / (swir1 + swir2 + eps),
    #TCG   = (green - red) / (green + red + eps),
    #MNDWI = (green - swir1) / (green + swir1 + eps),
    
    # New Additions
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

# Create mixing fractions (0 to 100% vegetation)
fractions <- seq(0, 1, by = 0.01)

# Initialize data table for mixtures
mixtures <- data.table(fraction_veg = fractions)

# Linear mixing of reflectance bands
# R_mix = f * R_veg + (1 - f) * R_soil
bands <- names(veg_spec)
for (b in bands) {
  mixtures[[b]] <- fractions * veg_spec[b] + (1 - fractions) * soil_spec[b]
}

# Calculate indices for all mixtures
mixtures <- calculate_indices(mixtures)

# ==============================================================================
# 4. Evaluate Linearity
# ==============================================================================

# Reshape to long format for plotting and analysis
indices <- c("DVI", "OSAVI", "MCARI", "NIRv", "PSRI", "NBR", 
             "TCW", "NDVI", "MSAVI2", "NDMI", "TCB", "GVI")

long_res <- melt(mixtures, id.vars = "fraction_veg", measure.vars = indices, 
                 variable.name = "Index", value.name = "Value")

# Calculate Linearity Metrics
# We compare the actual curve to a perfect line connecting the endpoints (0% and 100%)
linearity_scores <- long_res[, {
  # Endpoints
  y0 <- Value[fraction_veg == 0]
  y1 <- Value[fraction_veg == 1]
  
  # Ideal linear line
  linear_pred <- y0 + (y1 - y0) * fraction_veg
  
  # Residuals
  residuals <- Value - linear_pred
  
  # Metrics
  rmse <- sqrt(mean(residuals^2))
  max_dev <- max(abs(residuals))
  # Normalized max deviation (relative to range)
  range_val <- abs(y1 - y0)
  norm_max_dev <- ifelse(range_val > 1e-6, max_dev / range_val, 0)
  
  list(
    RMSE = rmse,
    Max_Deviation = max_dev,
    Normalized_Max_Dev = norm_max_dev,
    Range = range_val
  )
}, by = Index]

# Sort by linearity (Normalized Max Deviation ascending)
setorder(linearity_scores, Normalized_Max_Dev)

cat("\n=== LINEARITY EVALUATION (Sorted by most linear) ===\n")
print(linearity_scores)

# ==============================================================================
# 5. Plotting
# ==============================================================================

# Add the ideal linear line to the plot data for visualization
long_res[, Linear_Ref := {
  y0 <- Value[fraction_veg == 0]
  y1 <- Value[fraction_veg == 1]
  y0 + (y1 - y0) * fraction_veg
}, by = Index]

p <- ggplot(long_res, aes(x = fraction_veg)) +
  geom_line(aes(y = Linear_Ref), linetype = "dashed", color = "gray50") +
  geom_line(aes(y = Value, color = Index), size = 1) +
  facet_wrap(~Index, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Linearity of Spectral Indices (Vegetation Fraction)",
    subtitle = "Dashed line represents perfect linear mixing. Solid line is actual index behavior.",
    x = "Vegetation Fraction (0 = Soil, 1 = Veg)",
    y = "Index Value"
  ) +
  theme(legend.position = "none")

# Save plot
ggsave("index_linearity_check.png", p, width = 12, height = 10)

# Save metrics
write.csv(linearity_scores, "index_linearity_scores.csv", row.names = FALSE)

cat("\nAnalysis complete. Check 'index_linearity_check.png' and 'index_linearity_scores.csv'.\n")

# ==============================================================================
# 6. Vegetation-Vegetation Mixing Analysis
# ==============================================================================

cat("\n=== VEGETATION-VEGETATION LINEARITY ANALYSIS ===\n")

# Test mixing between different vegetation types
veg_veg_results <- list()

veg_pairs <- combn(names(veg_endmembers), 2, simplify = FALSE)

for (pair in veg_pairs) {
  veg1_name <- pair[1]
  veg2_name <- pair[2]
  veg1_spec <- veg_endmembers[[veg1_name]]
  veg2_spec <- veg_endmembers[[veg2_name]]
  
  cat(sprintf("\nTesting %s vs %s mixing...\n", veg1_name, veg2_name))
  
  # Create mixing fractions
  fractions <- seq(0, 1, by = 0.01)
  veg_mixtures <- data.table(fraction_veg1 = fractions)
  
  # Linear mixing of reflectance bands
  bands <- names(veg1_spec)
  for (b in bands) {
    veg_mixtures[[b]] <- fractions * veg1_spec[b] + (1 - fractions) * veg2_spec[b]
  }
  
  # Calculate indices
  veg_mixtures <- calculate_indices(veg_mixtures)
  
  # Evaluate linearity for each index
  indices <- c("DVI", "OSAVI", "MCARI", "NIRv", "PSRI", "NBR", 
               "TCW", "NDVI", "MSAVI2", "NDMI", "TCB", "GVI")
  
  pair_results <- data.table(
    Index = indices,
    Veg1 = veg1_name,
    Veg2 = veg2_name
  )
  
  for (idx in indices) {
    # Endpoints
    y0 <- veg_mixtures[[idx]][1]  # 100% veg1
    y1 <- veg_mixtures[[idx]][length(veg_mixtures[[idx]])]  # 100% veg2
    
    # Ideal linear line
    linear_pred <- y0 + (y1 - y0) * fractions
    
    # Residuals
    residuals <- veg_mixtures[[idx]] - linear_pred
    
    # Metrics
    rmse <- sqrt(mean(residuals^2))
    max_dev <- max(abs(residuals))
    range_val <- abs(y1 - y0)
    norm_max_dev <- ifelse(range_val > 1e-6, max_dev / range_val, 0)
    
    pair_results[Index == idx, `:=`(
      RMSE = rmse,
      Max_Deviation = max_dev,
      Normalized_Max_Dev = norm_max_dev,
      Range = range_val
    )]
  }
  
  veg_veg_results[[paste(veg1_name, veg2_name, sep = "_vs_")]] <- pair_results
}

# Combine all vegetation-vegetation results
all_veg_veg <- rbindlist(veg_veg_results)

# Summary statistics across all vegetation pairs
veg_veg_summary <- all_veg_veg[, .(
  Mean_Norm_Max_Dev = mean(Normalized_Max_Dev),
  Max_Norm_Max_Dev = max(Normalized_Max_Dev),
  Min_Norm_Max_Dev = min(Normalized_Max_Dev),
  SD_Norm_Max_Dev = sd(Normalized_Max_Dev),
  N_Pairs = .N
), by = Index]

setorder(veg_veg_summary, Mean_Norm_Max_Dev)

cat("\n=== VEGETATION-VEGETATION LINEARITY SUMMARY (Sorted by mean normalized max deviation) ===\n")
print(veg_veg_summary)

# Save vegetation-vegetation results
write.csv(all_veg_veg, "veg_veg_linearity_scores.csv", row.names = FALSE)
write.csv(veg_veg_summary, "veg_veg_linearity_summary.csv", row.names = FALSE)

# Compare soil-veg vs veg-veg linearity
comparison <- merge(
  linearity_scores[, .(Index, Soil_Veg_Norm_Max_Dev = Normalized_Max_Dev)],
  veg_veg_summary[, .(Index, Veg_Veg_Mean_Norm_Max_Dev = Mean_Norm_Max_Dev)],
  by = "Index"
)

setorder(comparison, Veg_Veg_Mean_Norm_Max_Dev)

cat("\n=== SOIL-VEG vs VEG-VEG LINEARITY COMPARISON ===\n")
print(comparison)

write.csv(comparison, "linearity_comparison.csv", row.names = FALSE)

# Create visualization for vegetation-vegetation mixing
# Use the first pair as example
example_pair <- veg_pairs[[1]]
veg1_name <- example_pair[1]
veg2_name <- example_pair[2]

fractions <- seq(0, 1, by = 0.01)
example_mixtures <- data.table(fraction_veg1 = fractions)

veg1_spec <- veg_endmembers[[veg1_name]]
veg2_spec <- veg_endmembers[[veg2_name]]

for (b in names(veg1_spec)) {
  example_mixtures[[b]] <- fractions * veg1_spec[b] + (1 - fractions) * veg2_spec[b]
}

example_mixtures <- calculate_indices(example_mixtures)

# Reshape for plotting
veg_veg_long <- melt(example_mixtures, id.vars = "fraction_veg1", 
                     measure.vars = indices, 
                     variable.name = "Index", value.name = "Value")

# Add linear reference
veg_veg_long[, Linear_Ref := {
  y0 <- Value[fraction_veg1 == 0]
  y1 <- Value[fraction_veg1 == 1]
  y0 + (y1 - y0) * fraction_veg1
}, by = Index]

p_veg <- ggplot(veg_veg_long, aes(x = fraction_veg1)) +
  geom_line(aes(y = Linear_Ref), linetype = "dashed", color = "gray50") +
  geom_line(aes(y = Value, color = Index), size = 1) +
  facet_wrap(~Index, scales = "free_y") +
  theme_minimal() +
  labs(
    title = sprintf("Vegetation-Vegetation Linearity (%s vs %s)", veg1_name, veg2_name),
    subtitle = "Dashed line represents perfect linear mixing. Solid line is actual index behavior.",
    x = sprintf("Fraction of %s (0 = 100%% %s, 1 = 100%% %s)", veg1_name, veg2_name, veg1_name),
    y = "Index Value"
  ) +
  theme(legend.position = "none")

ggsave("veg_veg_linearity_check.png", p_veg, width = 12, height = 10)

cat("\nVegetation-vegetation analysis complete.\n")
cat("Check 'veg_veg_linearity_scores.csv', 'veg_veg_linearity_summary.csv', 'linearity_comparison.csv', and 'veg_veg_linearity_check.png'.\n")
