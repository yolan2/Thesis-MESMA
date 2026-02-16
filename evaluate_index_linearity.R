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
  "C:\\Users\\yolan\\Downloads\\LS_S2_Harmonized_Timeseries.csv"
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
# Support legacy 'no.soil'/'no_soil' columns: populate or fill missing `Veg` values when appropriate
if (!"Veg" %in% names(df_raw)) {
  if ("no.soil" %in% names(df_raw)) {
    df_raw$Veg <- as.character(df_raw$no.soil)
    cat("[NOTICE] Using legacy 'no.soil' column as 'Veg'\n")
  } else if ("no_soil" %in% names(df_raw)) {
    df_raw$Veg <- as.character(df_raw$no_soil)
    cat("[NOTICE] Using legacy 'no_soil' column as 'Veg'\n")
  }
} else {
  # Fill missing entries in Veg from legacy columns if present
  if (any(is.na(df_raw$Veg)) && "no.soil" %in% names(df_raw)) {
    na_idx <- which(is.na(df_raw$Veg) & !is.na(df_raw$no.soil))
    if (length(na_idx) > 0) {
      df_raw$Veg[na_idx] <- as.character(df_raw$no.soil[na_idx])
      cat(sprintf("[NOTICE] Filled %d missing 'Veg' entries from 'no.soil'\n", length(na_idx)))
    }
  }
  if (any(is.na(df_raw$Veg)) && "no_soil" %in% names(df_raw)) {
    na_idx <- which(is.na(df_raw$Veg) & !is.na(df_raw$no_soil))
    if (length(na_idx) > 0) {
      df_raw$Veg[na_idx] <- as.character(df_raw$no_soil[na_idx])
      cat(sprintf("[NOTICE] Filled %d missing 'Veg' entries from 'no_soil'\n", length(na_idx)))
    }
  }
}
# Normalize Veg to trimmed character values for downstream comparisons
if ("Veg" %in% names(df_raw)) df_raw$Veg <- as.character(trimws(df_raw$Veg))

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

# Extract TRUE pure vegetation endmember (use only Veg == 'agri' or 'agriculture')
cat("\nExtracting TRUE pure vegetation endmember from data (Veg in {agri, agriculture})...\n")
pure_veg_data <- df_raw |> 
  dplyr::filter(!is.na(Veg) & tolower(trimws(Veg)) %in% c('agri', 'agric', 'agriculture', 'agricultural')) |> 
  dplyr::filter(across(c(blue, green, red, nir, swir1, swir2), is.finite))

if (nrow(pure_veg_data) == 0) {
  stop("No pure vegetation data found in Veg in {agri, agric, agriculture}! Cannot extract true vegetation endmember.")
}

cat(sprintf("Found %d pure vegetation observations (agri/agriculture)\n", nrow(pure_veg_data)))

# Calculate mean pure vegetation spectrum from agri/agriculture pool
veg_spec <- c(
  blue  = mean(pure_veg_data$blue, na.rm = TRUE),
  green = mean(pure_veg_data$green, na.rm = TRUE),
  red   = mean(pure_veg_data$red, na.rm = TRUE),
  nir   = mean(pure_veg_data$nir, na.rm = TRUE),
  swir1 = mean(pure_veg_data$swir1, na.rm = TRUE),
  swir2 = mean(pure_veg_data$swir2, na.rm = TRUE)
)

cat("True pure vegetation endmember spectrum (agri/agriculture):\n")
print(veg_spec)

# Extract TRUE vegetation endmembers for each type from pure vegetation data
cat("\nExtracting TRUE vegetation type endmembers from the PURE vegetation pool (using legacy labels only when present in that pool)...\n")
veg_types <- c("populus", "tamarix", "herbs")  # Focus on main types
veg_endmembers <- list()

# Operate only on the `pure_veg_data` subset (Veg in {agri, agriculture}) to ensure purity
has_no_soil <- "no.soil" %in% names(pure_veg_data)
has_no_soil_underscore <- "no_soil" %in% names(pure_veg_data)

for (veg_type in veg_types) {
  # Safely aggregate candidate rows from whichever label columns exist, avoiding direct references
  candidates <- df_raw[0, ]
  if ("Veg" %in% names(df_raw)) {
    candidates <- dplyr::bind_rows(candidates, df_raw |> dplyr::filter(!is.na(Veg) & tolower(trimws(as.character(Veg))) == veg_type))
  }
  if ("no.soil" %in% names(df_raw)) {
    candidates <- dplyr::bind_rows(candidates, df_raw |> dplyr::filter(!is.na(no.soil) & tolower(trimws(as.character(no.soil))) == veg_type))
  }
  if ("no_soil" %in% names(df_raw)) {
    candidates <- dplyr::bind_rows(candidates, df_raw |> dplyr::filter(!is.na(no_soil) & tolower(trimws(as.character(no_soil))) == veg_type))
  }

  # Ensure samples have finite reflectance in all required bands and remove duplicates
  candidates <- candidates |> dplyr::filter(across(c(blue, green, red, nir, swir1, swir2), is.finite)) |> dplyr::distinct()

  if (nrow(candidates) == 0) {
    cat(sprintf("  %s: NO labeled candidate samples found, skipping\n", veg_type))
    next
  }

  # Compute NDVI for candidate set (if needed)
  eps <- 1e-9
  if (!"NDVI" %in% names(candidates)) {
    candidates <- candidates |> dplyr::mutate(NDVI = (nir - red) / (nir + red + eps))
  }

  # Select high-NDVI samples to approximate 'pure' vegetation for the species
  ndvi_thresh <- 0.60
  selected <- candidates |> dplyr::filter(NDVI >= ndvi_thresh)

  # If none meet the threshold, take the top 10% by NDVI (ensures we pick the most vegetated samples)
  if (nrow(selected) == 0) {
    top_frac <- 0.10
    cutoff <- quantile(candidates$NDVI, probs = 1 - top_frac, na.rm = TRUE)
    selected <- candidates |> dplyr::filter(NDVI >= cutoff)
    cat(sprintf("[NOTICE] %s: no candidates >= %.2f NDVI; using top %.0f%% (NDVI >= %.3f) — %d samples\n", veg_type, ndvi_thresh, top_frac*100, cutoff, nrow(selected)))
  } else {
    cat(sprintf("  %s: %d candidates selected with NDVI >= %.2f\n", veg_type, nrow(selected), ndvi_thresh))
  }

  # Require a minimal number of samples; if too few, expand selection to top-50 samples by NDVI
  min_samples <- 3
  if (nrow(selected) < min_samples) {
    selected <- candidates |> dplyr::arrange(dplyr::desc(NDVI)) |> dplyr::slice_head(n = 50)
    cat(sprintf("[NOTICE] %s: insufficient selected samples (< %d). Falling back to top 50 NDVI candidates (%d samples)\n", veg_type, min_samples, nrow(selected)))
  }

  # Use the selected set to compute the mean spectrum
  veg_endmembers[[veg_type]] <- c(
    blue  = mean(selected$blue, na.rm = TRUE),
    green = mean(selected$green, na.rm = TRUE),
    red   = mean(selected$red, na.rm = TRUE),
    nir   = mean(selected$nir, na.rm = TRUE),
    swir1 = mean(selected$swir1, na.rm = TRUE),
    swir2 = mean(selected$swir2, na.rm = TRUE)
  )
  cat(sprintf("  %s: computed endmember from %d selected samples (from %d candidates)\n", veg_type, nrow(selected), nrow(candidates)))
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
        M_val <- suppressWarnings(max(df$DVI, na.rm = TRUE))
        if (!is.finite(M_val)) stop("[PPI] Cannot compute finite M for synthetic mixing")
        df$PPI <- ppi(dvi = df$DVI, zenith.angle = zenith_rad, M = M_val, dvi.soil = dvi_soil_val)
        cat(sprintf("Calculated PPI for synthetic mixing (dvi_soil=%.6f, zenith=%.4f rad)\n", dvi_soil_val, zenith_rad))
      }
    }
  }

  return(df)
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

# ==============================================================================
# 4. Linearity Diagnostics & Plots
# ==============================================================================
cat("\nPerforming linearity diagnostics and generating plots for mixtures...\n")

# Prepare output directory
OUTPUT_DIR <- "C:/MAP/linearity_results"
out_dir <- file.path(OUTPUT_DIR, "linearity_plots")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Indices to check (intersect with available columns)
# Include PPI when available (calculated if ppi_helpers.R is present)
candidate_indices <- c("NDVI", "DVI", "NIRv", "OSAVI", "MSAVI2", "EVI", "PSRI", "NBR", "TCW", "TCB", "GVI", "NDMI", "PPI")
indices_to_check <- intersect(candidate_indices, names(mixtures))
if (length(indices_to_check) == 0) {
  warning("No common indices found in mixtures to evaluate linearity.")
} else {
  cat(sprintf("Indices to be checked for linearity: %s\n", paste(indices_to_check, collapse = ", ")))
}

# Helper to analyze a mixture dataframe and save plots + coefficients
analyze_mixture <- function(df, label) {
  results <- list()
  for (idx in indices_to_check) {
    # Skip non-finite values
    dat <- df |> dplyr::select(fraction_veg, all_of(idx)) |> dplyr::filter(is.finite(.data[[idx]]))
    if (nrow(dat) < 3 || var(dat[[idx]], na.rm = TRUE) == 0) {
      cat(sprintf("[WARN] %s - index %s: insufficient data or zero variance, skipping\n", label, idx))
      next
    }

    fit <- tryCatch(lm(as.formula(paste(idx, "~ fraction_veg")), data = dat), error = function(e) return(NULL))
    if (is.null(fit)) {
      cat(sprintf("[WARN] %s - index %s: lm failed, skipping\n", label, idx))
      next
    }
    s <- summary(fit)
    slope <- coef(fit)["fraction_veg"]
    intercept <- coef(fit)["(Intercept)"]
    r2 <- as.numeric(s$r.squared)
    residuals <- resid(fit)
    rmse <- sqrt(mean(residuals^2, na.rm = TRUE))
    pval <- coef(s)["fraction_veg", "Pr(>|t|)"]

    # Build plot
    y_max <- max(dat[[idx]], na.rm = TRUE)
    y_min <- min(dat[[idx]], na.rm = TRUE)
    ann_y <- y_max - 0.05 * (y_max - y_min)
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = fraction_veg, y = .data[[idx]])) +
      ggplot2::geom_point(alpha = 0.6, size = 1) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, color = "#0072B2") +
      ggplot2::labs(title = paste0(label, " — ", idx), x = "Fraction vegetation", y = idx) +
      ggplot2::theme_minimal() +
      ggplot2::annotate("text", x = 0.02, y = ann_y,
                        label = sprintf("R²=%.3f\nRMSE=%.4f\nSlope=%.4f", r2, rmse, slope),
                        hjust = 0, size = 3)

    fname <- file.path(out_dir, sprintf("%s_%s.png", gsub("\\s+", "_", label), idx))
    tryCatch({
      ggplot2::ggsave(filename = fname, plot = p, width = 6, height = 4, dpi = 150)
      cat(sprintf("Saved plot: %s\n", fname))
    }, error = function(e) {
      warning(sprintf("Failed to save plot %s: %s", fname, e$message))
    })

    results[[length(results) + 1]] <- data.frame(
      label = label,
      index = idx,
      slope = as.numeric(slope),
      intercept = as.numeric(intercept),
      r_squared = as.numeric(r2),
      rmse = as.numeric(rmse),
      p_value = as.numeric(pval),
      n = nrow(dat)
    )
  }
  if (length(results) == 0) return(NULL)
  do.call(rbind, results)
}

# Analyze the 'agriculture' mixture (veg_spec) created earlier
scores <- list()
if (exists("mixtures")) {
  sc <- analyze_mixture(mixtures, label = "agriculture")
  if (!is.null(sc)) scores[[length(scores) + 1]] <- sc
}

# Also analyze each species-specific endmember mixture
for (vt in names(veg_endmembers)) {
  veg_s <- veg_endmembers[[vt]]
  m <- data.table(fraction_veg = fractions)
  for (b in names(veg_s)) m[[b]] <- fractions * veg_s[b] + (1 - fractions) * soil_spec[b]
  m <- calculate_indices(m)
  sc <- analyze_mixture(m, label = paste0("species_", vt))
  if (!is.null(sc)) scores[[length(scores) + 1]] <- sc
}

# Combine and write out scores
if (length(scores) > 0) {
  scores_df <- do.call(rbind, scores)
  scores_out <- file.path(OUTPUT_DIR, "index_linearity_scores_with_coefficients.csv")
  write.csv(scores_df, file = scores_out, row.names = FALSE)
  cat(sprintf("Wrote linearity scores to %s\n", scores_out))
} else {
  warning("No linearity scores were produced — no index/species combos met criteria.")
}

