# evaluate_index_linearity.R
# Script to evaluate the linearity of spectral indices through synthetic mixing
# of TRUE Barren and Pure Vegetation endmembers from actual data, and between different Vegetation types.

library(data.table)
library(ggplot2)
library(dplyr)
library(lubridate)

options(warn = -1)

# Raw spectral bands (optional - included if present)
RAW_BANDS <- c("blue", "green", "red", "nir", "swir1", "swir2")

# Normalize known raw band column names (case and prefix variants) to canonical lower-case names
normalize_band_names <- function(df, bands = RAW_BANDS) {
  if (is.null(df) || nrow(df) == 0) return(df)
  current_names <- names(df)
  for (b in bands) {
    # Candidate variations to match
    candidates <- c(b, toupper(b), tools::toTitleCase(b), paste0('band_', b), toupper(paste0('band_', b)), paste0('Band_', b))
    for (cand in candidates) {
      if (cand %in% current_names && !(b %in% current_names)) {
        names(df)[names(df) == cand] <- b
        cat(sprintf("[NOTICE] Normalized band column '%s' -> '%s'\n", cand, b))
        current_names <- names(df)
        break
      }
    }
  }
  df
}

# ==============================================================================
# 0. Load Real Data and Extract Endmembers
# ==============================================================================

cat("Loading real phenology data to extract true endmembers...\n")

# Try to load from common data file locations
data_file <- NULL
possible_paths <- c(
  "C:\\Users\\yolan\\Downloads\\LS_S2_Harmonized_Timeseries.csv",
  "C:\\Users\\yolan\\OneDrive\\Documenten\\UGENT\\Master\\masterproef\\GIS\\landsat_lower.csv",
  "../masterproef/GIS/landsat_lower.csv",
  "landsat_lower.csv"
)

for (path in possible_paths) {
  if (file.exists(path)) {
    data_file <- path
    break
  }
}

if (is.null(data_file)) {
  stop("Could not find phenology data file. Please ensure landsat_lower.csv is accessible.")
}

cat(sprintf("Loading data from: %s\n", data_file))
df_raw <- fread(data_file)

# Normalize band names to canonical lower-case
df_raw <- normalize_band_names(df_raw)

# GeoJSON support removed: expect Veg to be provided in the CSV
if ("vegetation" %in% names(df_raw) && !"Veg" %in% names(df_raw)) {
  df_raw$Veg <- df_raw$vegetation
  cat("[NOTICE] Renamed 'vegetation' -> 'Veg' in phenology data\n")
}
# Legacy 'no.soil'/'no_soil' columns are ignored and not used by this script

# Require that Veg be present to extract endmembers
if (!("Veg" %in% names(df_raw))) {
  stop("Phenology data must include a 'Veg' column. GeoJSON support has been removed — include this field in landsat_lower.csv.")
}

# Extract TRUE barren endmember (Veg == 'barren')
cat("\nExtracting TRUE barren endmember from data...\n")
barren_data <- df_raw |> 
  dplyr::filter(!is.na(Veg) & tolower(Veg) == "barren") |> 
  dplyr::filter(across(c(blue, green, red, nir, swir1, swir2), is.finite))

if (nrow(barren_data) == 0) {
  stop("No barren data found! Cannot extract true soil endmember.")
}

cat(sprintf("Found %d barren observations\n", nrow(barren_data)))

# Calculate mean soil spectrum
soil_spec <- c(
  blue  = mean(barren_data$blue, na.rm = TRUE),
  green = mean(barren_data$green, na.rm = TRUE),
  red   = mean(barren_data$red, na.rm = TRUE),
  nir   = mean(barren_data$nir, na.rm = TRUE),
  swir1 = mean(barren_data$swir1, na.rm = TRUE),
  swir2 = mean(barren_data$swir2, na.rm = TRUE)
)

cat("True soil endmember spectrum:\n")
print(soil_spec)

# Extract TRUE pure vegetation endmember (identified by Veg != 'barren')
cat("\nExtracting TRUE pure vegetation endmember from data...\n")
pure_veg_data <- df_raw |> 
  dplyr::filter(!is.na(Veg) & tolower(Veg) != 'barren') |> 
  dplyr::filter(across(c(blue, green, red, nir, swir1, swir2), is.finite))

if (nrow(pure_veg_data) == 0) {
  stop("No pure vegetation data found (no vegetation types present)! Cannot extract true vegetation endmember.")
}

cat(sprintf("Found %d pure vegetation observations\n", nrow(pure_veg_data)))

# Calculate mean pure vegetation spectrum
veg_spec <- c(
  blue  = mean(pure_veg_data$blue, na.rm = TRUE),
  green = mean(pure_veg_data$green, na.rm = TRUE),
  red   = mean(pure_veg_data$red, na.rm = TRUE),
  nir   = mean(pure_veg_data$nir, na.rm = TRUE),
  swir1 = mean(pure_veg_data$swir1, na.rm = TRUE),
  swir2 = mean(pure_veg_data$swir2, na.rm = TRUE)
)

cat("True pure vegetation endmember spectrum:\n")
print(veg_spec)

# Extract TRUE vegetation endmembers for each type from pure vegetation data
cat("\nExtracting TRUE vegetation type endmembers from pure vegetation data...\n")
veg_types <- c("populus", "tamarix", "phragmites")  # Focus on main types
veg_endmembers <- list()

for (veg_type in veg_types) {
  veg_subset <- pure_veg_data |> 
    dplyr::filter(!is.na(Veg) & tolower(Veg) == veg_type)

  if (nrow(veg_subset) > 0) {
    veg_endmembers[[veg_type]] <- c(
      blue  = mean(veg_subset$blue, na.rm = TRUE),
      green = mean(veg_subset$green, na.rm = TRUE),
      red   = mean(veg_subset$red, na.rm = TRUE),
      nir   = mean(veg_subset$nir, na.rm = TRUE),
      swir1 = mean(veg_subset$swir1, na.rm = TRUE),
      swir2 = mean(veg_subset$swir2, na.rm = TRUE)
    )
    cat(sprintf("  %s: %d pure samples\n", veg_type, nrow(veg_subset)))
  } else {
    cat(sprintf("  %s: NO pure samples found, skipping\n", veg_type))
  }
}

if (length(veg_endmembers) == 0) {
  stop("No vegetation type endmembers could be extracted!")
}

cat("\nExtracted endmembers for vegetation types:\n")
print(names(veg_endmembers))

# ==============================================================================
# 2. Define Index Calculation Function
# ==============================================================================

# Source PPI helpers for consistent calculation
if (file.exists("ppi_helpers.R")) {
  source("ppi_helpers.R")
  cat("Loaded ppi_helpers.R for PPI calculation\n")
} else {
  warning("ppi_helpers.R not found - PPI will not be calculated")
}

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
    MSAVI  = (2 * nir + 1 - sqrt(pmax(0, (2 * nir + 1)^2 - 8 * (nir - red)))) / 2,
    NDMI   = (nir - swir1) / (nir + swir1 + eps),
    TCB    = 0.3029 * blue + 0.2786 * green + 0.4733 * red + 0.5599 * nir + 0.508 * swir1 + 0.1872 * swir2,
    GVI    = -0.2941 * blue - 0.243 * green - 0.5424 * red + 0.7276 * nir + 0.0713 * swir1 - 0.1608 * swir2,
    SATVI  = (swir1 - red) / (swir1 + red + 0.5) * (1 + 0.5),
    EVI    = 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)
  )]

  # Add PPI if ppi helpers are available
  if (exists("ppi") && exists("calculate_solar_zenith")) {
    # For synthetic mixing: use soil endmember DVI as dvi_soil
    # Assume first row (fraction_veg = 0) is pure soil
    if ("fraction_veg" %in% names(df)) {
      dvi_soil_val <- df$DVI[df$fraction_veg == 0][1]
      if (is.finite(dvi_soil_val)) {
        # Calculate zenith angle (use typical mid-latitude, mid-season value)
        # Assume lat=40°N, DOY=180 (summer solstice), 10:30 AM
        zenith_rad <- calculate_solar_zenith(lat = 40, doy = 180, hour = 10.5)

        # Calculate PPI using ppi() function
        df$PPI <- ppi(dvi = df$DVI, zenith.angle = zenith_rad, M = 0.7, dvi.soil = dvi_soil_val)
        cat(sprintf("Calculated PPI for synthetic mixing (dvi_soil=%.6f, zenith=%.4f rad)\n", dvi_soil_val, zenith_rad))
      }
    }
  }

  return(df)
}

# ==============================================================================
# 2.5. Linearize Indices Function
# ==============================================================================

linearize_indices <- function(df) {
  cat("Applying linearization transformations to indices...\n")
  
  # NIRv: 2*x - x^2 (M_NDVI transform)
  if ("NIRv" %in% names(df)) {
    cat("  Linearizing NIRv (2x - x^2)...\n")
    df$NIRv <- 2 * df$NIRv - df$NIRv^2
  }
  
  # MSAVI2: log(x + 1)
  if ("MSAVI2" %in% names(df)) {
    cat("  Linearizing MSAVI2 (log(x+1))...\n")
    df$MSAVI2 <- log(df$MSAVI2 + 1)
  }

  # MSAVI: same treatment as MSAVI2 for comparability
  if ("MSAVI" %in% names(df)) {
    cat("  Linearizing MSAVI (log(x+1))...\n")
    df$MSAVI <- log(df$MSAVI + 1)
  }
  
  # PSRI: Signed Sqrt
  if ("PSRI" %in% names(df)) {
    cat("  Linearizing PSRI (signed sqrt)...\n")
    df$PSRI <- sign(df$PSRI) * sqrt(abs(df$PSRI))
  }
  
  # TCW: log(x + 1) - handle negatives safely
  if ("TCW" %in% names(df)) {
    cat("  Linearizing TCW (log(x+1))...\n")
    # Ensure x+1 > 0. If TCW < -0.99, clamp?
    # TCW is usually > -1.
    df$TCW <- log(pmax(df$TCW + 1, 1e-6))
  }
  
  df
}

# Raw spectral bands (optional - included if present)
RAW_BANDS <- c("blue", "green", "red", "nir", "swir1", "swir2")

# Compute a set of spectral indices from raw bands (compatible with R_extract_hls.R's formulae)
compute_indices_from_bands <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  eps <- 1e-9
  has_bands <- intersect(RAW_BANDS, names(df))
  if (length(has_bands) == 0) return(df)

  # Safely compute each index only if its required bands are present
  if (all(c('nir','red') %in% names(df))) df$DVI <- as.numeric(df$nir) - as.numeric(df$red)
  if (all(c('nir','red') %in% names(df))) df$OSAVI <- (as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + as.numeric(df$red) + 0.16)
  if (all(c('red','green','blue') %in% names(df))) df$MCARI <- ((as.numeric(df$red) - as.numeric(df$green)) - 0.2*(as.numeric(df$red) - as.numeric(df$blue))) * (as.numeric(df$red) / (as.numeric(df$green) + eps))
  if (all(c('green','red') %in% names(df))) df$PRI <- (as.numeric(df$green) - as.numeric(df$red)) / (as.numeric(df$green) + as.numeric(df$red) + eps)
  if (all(c('nir','red') %in% names(df))) df$NIRv <- as.numeric(df$nir) * ((as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + as.numeric(df$red) + eps))
  if (all(c('red','blue','nir') %in% names(df))) df$PSRI <- (as.numeric(df$red) - as.numeric(df$blue)) / (as.numeric(df$nir) + eps)
  if (all(c('nir','swir2') %in% names(df))) df$NBR <- (as.numeric(df$nir) - as.numeric(df$swir2)) / (as.numeric(df$nir) + as.numeric(df$swir2) + eps)
  if (all(c('swir1','swir2') %in% names(df))) df$TCW <- (as.numeric(df$swir1) - as.numeric(df$swir2)) / (as.numeric(df$swir1) + as.numeric(df$swir2) + eps)
  if (all(c('green','red','nir','swir1','swir2','blue') %in% names(df))) {
    df$TCB <- 0.3029 * as.numeric(df$blue) + 0.2786 * as.numeric(df$green) + 0.4733 * as.numeric(df$red) + 0.5599 * as.numeric(df$nir) + 0.508 * as.numeric(df$swir1) + 0.1872 * as.numeric(df$swir2)
    df$GVI <- -0.2941 * as.numeric(df$blue) - 0.243 * as.numeric(df$green) - 0.5424 * as.numeric(df$red) + 0.7276 * as.numeric(df$nir) + 0.0713 * as.numeric(df$swir1) - 0.1608 * as.numeric(df$swir2)
  }
  if (all(c('nir','red') %in% names(df))) df$NDVI <- (as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + as.numeric(df$red) + eps)
  if (all(c('nir','red') %in% names(df))) df$MSAVI2 <- (2 * as.numeric(df$nir) + 1 - sqrt(pmax(0, (2 * as.numeric(df$nir) + 1)^2 - 8 * (as.numeric(df$nir) - as.numeric(df$red))))) / 2
  if (all(c('nir','red') %in% names(df))) df$MSAVI <- (2 * as.numeric(df$nir) + 1 - sqrt(pmax(0, (2 * as.numeric(df$nir) + 1)^2 - 8 * (as.numeric(df$nir) - as.numeric(df$red))))) / 2
  if (all(c('nir','swir1') %in% names(df))) df$NDMI <- (as.numeric(df$nir) - as.numeric(df$swir1)) / (as.numeric(df$nir) + as.numeric(df$swir1) + eps)

  # Apply kNDVI-like tweak used in extractor (NIRv * 1.3)
  if ('NIRv' %in% names(df)) df$NIRv <- df$NIRv * 1.3

  df
}

# Function to calculate robust variance using STL decomposition + MAD
# Helper: compute MAD^2 with minimal sample requirement (top-level helper, no nested defs)
compute_mad2 <- function(x, min_samples = 3) {
  x <- x[is.finite(x)]
  if (length(x) < min_samples) return(NA_real_)
  m <- mad(x, na.rm = TRUE, constant = 1.4826)
  if (!is.finite(m)) return(NA_real_)
  m^2
}

# Safe multiplication helper: explicitly handle recycling and avoid implicit warnings.
# Returns element-wise product (with explicit replication of the shorter vector when appropriate),
# or throws an informative error when lengths are incompatible.
safe_mul_vec <- function(a, b, allow_recycle = TRUE, caller = NULL) {
  la <- length(a); lb <- length(b)
  if (la == 0 || lb == 0) return(numeric(0))
  if (la == lb) return(a * b)
  if (la == 1) return(rep(a, lb) * b)
  if (lb == 1) return(a * rep(b, la))
  if (allow_recycle && (la %% lb == 0 || lb %% la == 0)) {
    # replicate shorter to the longer length explicitly
    if (la < lb) a <- rep(a, length.out = lb) else b <- rep(b, length.out = la)
    return(a * b)
  }
  caller_text <- if (is.null(caller)) "safe_mul_vec" else paste0(caller, ": ")
  stop(sprintf("%sIncompatible lengths for multiplication: %d vs %d", caller_text, la, lb))
}

# Safe dot product (handles NA removal and length mismatch via replication rules above)
safe_dot <- function(a, b, na.rm = TRUE) {
  if (length(a) == 0 || length(b) == 0) return(0)
  prod <- safe_mul_vec(a, b, allow_recycle = TRUE, caller = "safe_dot")
  if (na.rm) sum(prod, na.rm = TRUE) else sum(prod)
}

# Safe column-weighted average helper: ensures wts length equals number of rows
safe_col_weighted_avg <- function(mat, wts) {
  if (is.null(mat) || nrow(mat) == 0) return(rep(0, ifelse(is.null(mat), 1, ncol(mat))))
  n <- nrow(mat)
  if (length(wts) == 0) wts <- rep(1, n)
  if (length(wts) != n) {
    # replicate shorter vector to length n when compatible, otherwise error
    if (length(wts) == 1 || n %% length(wts) == 0 || length(wts) %% n == 0) {
      wts <- rep(wts, length.out = n)
      warning(sprintf("safe_col_weighted_avg: adjusted weights vector to length %d", n))
    } else {
      stop(sprintf("safe_col_weighted_avg: weights length (%d) incompatible with rows in mat (%d)", length(wts), n))
    }
  }
  if (sum(wts, na.rm = TRUE) == 0) {
    # fallback: unweighted mean
    return(as.numeric(colMeans(mat, na.rm = TRUE)))
  }
  wts <- as.numeric(wts) / sum(wts, na.rm = TRUE)
  as.numeric(colSums(mat * wts, na.rm = TRUE))
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

# Linearize indices
# mixtures <- linearize_indices(mixtures)

# ==============================================================================
# 4. Evaluate Linearity
# ==============================================================================

# Reshape to long format for plotting and analysis
indices <- c("DVI", "OSAVI", "MCARI", "NIRv", "PSRI", "NBR",
             "TCW", "NDVI", "MSAVI2", "MSAVI", "NDMI", "TCB", "GVI", "PPI", "SATVI", "EVI")

long_res <- melt(mixtures, id.vars = "fraction_veg", measure.vars = indices, 
                 variable.name = "Index", value.name = "Value")

target_aic_indices <- c("PPI", "NIRv", "OSAVI", "NDMI", "MSAVI", "NDVI", "SATVI", "EVI")

compute_aic_table <- function(long_dt, response_col = "fraction_veg", value_col = "Value", index_col = "Index", targets = target_aic_indices) {
  if (!data.table::is.data.table(long_dt)) long_dt <- as.data.table(long_dt)

  results <- lapply(targets, function(idx) {
    sub <- long_dt[get(index_col) == idx]
    sub <- sub[is.finite(get(response_col)) & is.finite(get(value_col))]
    n_obs <- nrow(sub)

    if (n_obs < 3) {
      return(data.table(Index = idx, AIC = NA_real_, LogLik = NA_real_, Df = NA_integer_, N = n_obs, Note = "insufficient data"))
    }

    fit <- tryCatch({
      lm_out <- lm(sub[[response_col]] ~ sub[[value_col]])
      list(aic = AIC(lm_out), loglik = as.numeric(logLik(lm_out)), df = length(coef(lm_out)))
    }, error = function(e) NULL)

    if (is.null(fit)) {
      return(data.table(Index = idx, AIC = NA_real_, LogLik = NA_real_, Df = NA_integer_, N = n_obs, Note = "fit_failed"))
    }

    data.table(Index = idx, AIC = fit$aic, LogLik = fit$loglik, Df = fit$df, N = n_obs, Note = NA_character_)
  })

  out <- rbindlist(results, use.names = TRUE, fill = TRUE)
  setorder(out, AIC)
  out
}

# Calculate Linearity Metrics using linear regression
linearity_scores <- long_res[, {
  # Fit linear model fraction ~ index
  lm_fit <- lm(fraction_veg ~ Value)
  r2 <- summary(lm_fit)$r.squared
  residuals <- residuals(lm_fit)
  rmse <- sqrt(mean(residuals^2))
  max_dev <- max(abs(residuals))
  range_val <- diff(range(Value))
  norm_max_dev <- ifelse(range_val > 1e-6, max_dev / range_val, 0)
  
  list(
    R2 = r2,
    RMSE = rmse,
    Max_Deviation = max_dev,
    Normalized_Max_Dev = norm_max_dev,
    Range = range_val
  )
}, by = Index]

# Sort by linearity (R2 descending)
setorder(linearity_scores, -R2)

cat("\n=== LINEARITY EVALUATION (Sorted by R2 descending) ===\n")
print(linearity_scores)

# AIC comparison focused on key indices
aic_table <- compute_aic_table(long_res, response_col = "fraction_veg", targets = target_aic_indices)
cat("\n=== AIC COMPARISON (Lower is better) ===\n")
print(aic_table)
write.csv(aic_table, "target_index_aic.csv", row.names = FALSE)

# Filter linearizable indices
threshold_r2 <- 0.95
linearizable <- linearity_scores[R2 >= threshold_r2]
unlinearizable <- linearity_scores[R2 < threshold_r2]

cat("\n=== LINEARIZABLE INDICES (R2 >= ", threshold_r2, ") ===\n")
print(linearizable)

cat("\n=== UNLINEARIZABLE INDICES ===\n")
print(unlinearizable)

# Save linearizable
write.csv(linearizable, "linearizable_indices.csv", row.names = FALSE)

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
  
  # Linearize indices
  # veg_mixtures <- linearize_indices(veg_mixtures)
  
  # Evaluate linearity for each index
  indices <- c("DVI", "OSAVI", "MCARI", "NIRv", "PSRI", "NBR",
               "TCW", "NDVI", "MSAVI2", "MSAVI", "NDMI", "TCB", "GVI", "PPI", "SATVI", "EVI")
  
  pair_results <- data.table(
    Index = indices,
    Veg1 = veg1_name,
    Veg2 = veg2_name
  )
  
  for (idx in indices) {
    # Check if index exists and has valid data
    if (!idx %in% names(veg_mixtures) || is.null(veg_mixtures[[idx]]) || all(is.na(veg_mixtures[[idx]]))) {
      cat(sprintf("  Skipping %s (not available or all NA)\n", idx))
      pair_results[Index == idx, `:=`(
        R2 = NA_real_,
        RMSE = NA_real_,
        Max_Deviation = NA_real_,
        Normalized_Max_Dev = NA_real_,
        Range = NA_real_
      )]
      next
    }
    
    # Fit linear model
    tryCatch({
      lm_fit <- lm(fractions ~ veg_mixtures[[idx]])
      r2 <- summary(lm_fit)$r.squared
      residuals <- residuals(lm_fit)
      rmse <- sqrt(mean(residuals^2))
      max_dev <- max(abs(residuals))
      range_val <- diff(range(veg_mixtures[[idx]], na.rm = TRUE))
      norm_max_dev <- ifelse(range_val > 1e-6, max_dev / range_val, 0)
      
      pair_results[Index == idx, `:=`(
        R2 = r2,
        RMSE = rmse,
        Max_Deviation = max_dev,
        Normalized_Max_Dev = norm_max_dev,
        Range = range_val
      )]
    }, error = function(e) {
      cat(sprintf("  Error fitting %s: %s\n", idx, e$message))
      pair_results[Index == idx, `:=`(
        R2 = NA_real_,
        RMSE = NA_real_,
        Max_Deviation = NA_real_,
        Normalized_Max_Dev = NA_real_,
        Range = NA_real_
      )]
    })
  }
  
  veg_veg_results[[paste(veg1_name, veg2_name, sep = "_vs_")]] <- pair_results
}

# Combine all vegetation-vegetation results
all_veg_veg <- rbindlist(veg_veg_results)

# Summary statistics across all vegetation pairs
veg_veg_summary <- all_veg_veg[, .(
  Mean_R2 = mean(R2),
  Max_R2 = max(R2),
  Min_R2 = min(R2),
  SD_R2 = sd(R2),
  N_Pairs = .N
), by = Index]

setorder(veg_veg_summary, -Mean_R2)

cat("\n=== VEGETATION-VEGETATION LINEARITY SUMMARY (Sorted by mean R2 descending) ===\n")
print(veg_veg_summary)

# Save vegetation-vegetation results
write.csv(all_veg_veg, "veg_veg_linearity_scores.csv", row.names = FALSE)
write.csv(veg_veg_summary, "veg_veg_linearity_summary.csv", row.names = FALSE)

# Compare soil-veg vs veg-veg linearity
comparison <- merge(
  linearity_scores[, .(Index, Soil_Veg_R2 = R2)],
  veg_veg_summary[, .(Index, Veg_Veg_Mean_R2 = Mean_R2)],
  by = "Index"
)

setorder(comparison, -Veg_Veg_Mean_R2)

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

# Linearize indices
# example_mixtures <- linearize_indices(example_mixtures)

# Reshape for plotting
# Only use indices that are available in example_mixtures
available_indices <- intersect(indices, names(example_mixtures))
veg_veg_long <- melt(example_mixtures, id.vars = "fraction_veg1", 
                     measure.vars = available_indices, 
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

# Vegetation-vegetation analysis complete.
