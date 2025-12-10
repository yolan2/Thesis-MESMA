 
library(zoo)
library(dplyr)
library(cluster)
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)
if (!requireNamespace("scales", quietly = TRUE)) {
  install.packages("scales")
}
library(scales)
if (!requireNamespace("nlme", quietly = TRUE)) {
  install.packages("nlme")
}
library(nlme)
if (!requireNamespace("pbapply", quietly = TRUE)) {
  install.packages("pbapply")
}
library(pbapply)
options(warn = 1)  # print warnings as they occur for debugging

# =============================================================================
# USER-MODIFIABLE PARAMETERS
# =============================================================================

# Skip PPI calculation to avoid zenith.angle errors
SKIP_PPI <- FALSE
PPI_DVI_SOIL <- 0.09  # Set to match january_averages.R default

# Input/Output Files
INPUT_CSV <- "phenology_results/hls_phenology_data.csv"
INFERENCE_CSV <- "C:\\Users\\yolan\\OneDrive\\Documenten\\UGENT\\Master\\masterproef\\GIS\\landsat_lower_inference.csv"

# Training Configuration
TRAIN_YEARS <- c(2024)  # Years to use for training (note: script uses all data)

# Parallel Processing
PARALLEL_ENABLE <- TRUE
PARALLEL_WORKERS <- 4
COMBO_PARALLEL_WORKERS <- max(1L, floor(PARALLEL_WORKERS/2))
PERSISTENT_PARALLEL_BACKEND <- TRUE

# MESMA Configuration
ENABLE_DIAGNOSTICS <- TRUE
# ENABLE_GEOMETRIC_MESMA <- TRUE  # Always use geometric MESMA now
TOP_K_CANDIDATES <- 10L
MAX_VARIANTS_PER_VEG <- 10
MIN_VARIANTS_PER_VEG <- 1
MAX_COMBOS_FOR_FULL_SEARCH <- 5000

# Vegetation Classes
ALLOWED_VEG <- c("populus", "tamarix", "phragmites")

# Spectral Indices
OPTIMAL_INDICES <- c("DVI", "GVI", "NIRv", "MSAVI2", "PSRI", "TCW")

# Inference Control (skip when testing the script)
SKIP_INFERENCE <- FALSE

# Source PPI helpers to ensure consistent PPI computation across scripts
if (file.exists("ppi_helpers.R")) {
  source("ppi_helpers.R")
} else {
  warning("ppi_helpers.R not found; using built-in PPI calculation. To standardize PPI, add ppi_helpers.R to the project root.")
}

# Logging
PROGRESS_LOG_TO_FILE <- FALSE
LOG_FILE <- "mesma_progress.log"

# =============================================================================
# FUNCTIONS AND IMPLEMENTATION
# =============================================================================

fit_cost_mkl <- function(obs, weight, t_row) {
  # fit_cost_mkl: vectorized cost calculation helper
  if (length(obs) != length(t_row)) t_row <- t_row[seq_len(length(obs))]
  # Ensure weight length is compatible: allow scalar rep or exact-match length; otherwise repeat with warning
  if (length(weight) == 1L) {
    weight <- rep(weight, length(obs))
  } else if (length(weight) != length(obs)) {
    # If weight length is a multiple of obs (or vice versa) explicitly replicate
    if (length(weight) %% length(obs) == 0 || length(obs) %% length(weight) == 0) {
      warning(sprintf("fit_cost_mkl: duplicating/replicating weight vector (length %d) to match obs length %d", length(weight), length(obs)))
      weight <- rep(weight, length.out = length(obs))
    } else {
      stop(sprintf("fit_cost_mkl: weight length (%d) incompatible with obs length (%d)", length(weight), length(obs)))
    }
      found_inf_raw <- intersect(RAW_BANDS, names(df_inf))
      if (length(found_inf_raw) > 0) {
        cat(sprintf("Found %d RAW_BANDS in inference file: %s\n", length(found_inf_raw), paste(found_inf_raw, collapse = ", ")))
  }
  if (any(!is.finite(c(obs, weight, t_row)))) stop("fit_cost_mkl: inputs contain non-finite values; clean data before calling")
  
  # Use vectorized operations
  # Use explicit safe multiplication to avoid R recycling warnings
  pred <- as.numeric(safe_mul_vec(weight, t_row))
  obs_weighted <- as.numeric(safe_mul_vec(obs, weight))
  residuals <- obs_weighted - pred
  sum(residuals^2)
}

# PPI calculation is provided by `ppi_helpers.R` (sourced earlier) to ensure consistency across scripts.
if (!exists("ppi") && file.exists("ppi_helpers.R")) source("ppi_helpers.R")

## Prefer centralized `add_ppi_columns` from ppi_helpers.R, fallback to local definition
if (!exists("add_ppi_columns") && file.exists("ppi_helpers.R")) source("ppi_helpers.R")
if (!exists("add_ppi_columns")) add_ppi_columns <- function(df, lat) {
  # Ensure df is a data.frame (not tibble) for reliable column access
  df <- as.data.frame(df)
  
  # Calculate DVI from nir and red bands
  df[["DVI"]] <- df[["nir"]] - df[["red"]]

  # Ensure date, doy, and year columns exist
  df[["date"]] <- as.Date(df[["date"]])
  if (!"doy" %in% names(df)) df[["doy"]] <- lubridate::yday(df[["date"]])
  if (!"year" %in% names(df)) df[["year"]] <- lubridate::year(df[["date"]])
  if (!"pheno_year" %in% names(df)) df[["pheno_year"]] <- ifelse(lubridate::month(df[["date"]]) >= 3, lubridate::year(df[["date"]]), lubridate::year(df[["date"]]) - 1)

  # Peak stats per location-year (avoids per-DOY PPI time series)
  peak_df <- df %>%
    group_by(location_id, pheno_year) %>%
    summarise(
      DVI_max = max(DVI, na.rm = TRUE),
      doy_peak = {
        if (all(is.na(DVI))) NA_integer_ else {
          idx <- which.max(DVI)
          doy[idx[1]]
        }
      },
      .groups = "drop"
    )

  peak_df$DVI_max[!is.finite(peak_df$DVI_max)] <- NA_real_
  peak_df$doy_peak[!is.finite(peak_df$doy_peak)] <- NA_integer_

  # Use per-location latitude when available; otherwise fall back to provided lat
  if ("lat" %in% names(df)) {
    lat_lookup <- df %>%
      group_by(location_id) %>%
      summarise(lat_use = mean(lat, na.rm = TRUE), .groups = "drop")
  } else {
    lat_lookup <- data.frame(location_id = unique(df$location_id), lat_use = lat)
  }

  peak_df <- peak_df %>% left_join(lat_lookup, by = "location_id")
  peak_df$lat_use[!is.finite(peak_df$lat_use)] <- lat

  # Single zenith angle per loc-year using DOY of peak DVI
  peak_df$zenith.angle <- NA_real_
  zen_idx <- complete.cases(peak_df$lat_use, peak_df$doy_peak)
  if (any(zen_idx)) {
    peak_df$zenith.angle[zen_idx] <- calculate_solar_zenith(
      lat = peak_df$lat_use[zen_idx],
      doy = peak_df$doy_peak[zen_idx]
    )
  }

  # Single PPI per loc-year using the peak DVI
  peak_df$PPI <- NA_real_
  ppi_idx <- complete.cases(peak_df$DVI_max, peak_df$zenith.angle)
  if (any(ppi_idx)) {
    peak_df$PPI[ppi_idx] <- ppi(
      dvi = peak_df$DVI_max[ppi_idx],
      zenith.angle = peak_df$zenith.angle[ppi_idx],
      M = peak_df$DVI_max[ppi_idx] + 0.005,
      dvi.soil = PPI_DVI_SOIL
    )
  }

  # Drop any existing PPI/DVI_max/zenith columns before join to avoid suffixes
  drop_cols <- intersect(c("DVI_max", "PPI", "zenith.angle"), names(df))
  if (length(drop_cols)) df[drop_cols] <- NULL

  df <- df %>%
    left_join(peak_df %>% select(location_id, year, DVI_max, zenith.angle, PPI),
              by = c("location_id", "year"))

  df <- as.data.frame(df)
  return(df)
}

# Function to calculate solar zenith angle (in radians)
# Assumes 10:30 AM local solar time (typical for Landsat/Sentinel)
## calculate_solar_zenith provided by `ppi_helpers.R` if available
if (!exists("calculate_solar_zenith") && file.exists("ppi_helpers.R")) source("ppi_helpers.R")

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

calc_moving_var <- function(df, index_name, window = 14, span_loess = 0.1, min_obs_loess = 6) {
  # REPLACED: Moving Variance -> Moving Min/Max Range (Linearity-Preserving)
  # This function now calculates the range (Max - Min) within the window.
  # Range is roughly linear for mixing (Range(A+B) approx Range(A) + Range(B) if phases align)
  # Ideally we would use Mean, but this function is used for "variability" features.
  # So we use Range as a proxy for variability that behaves better linearly than Variance.
  
  if (!"date" %in% names(df)) stop("calc_moving_var: df must have a 'date' column")
  if (!index_name %in% names(df)) stop(sprintf("calc_moving_var: index '%s' not found in df", index_name))

  n <- nrow(df)
  out <- rep(NA_real_, n)

  process_series_range <- function(dates, vals) {
    if (length(vals) < 1) return(rep(NA_real_, length(vals)))
    
    # Build regular series
    dts <- as.Date(dates)
    full_days <- seq(from = as.Date(min(dts, na.rm = TRUE)), to = as.Date(max(dts, na.rm = TRUE)), by = "day")
    full_vec <- rep(NA_real_, length(full_days))
    pos_map <- match(dts, full_days)
    full_vec[pos_map] <- as.numeric(vals)
    
    # Interpolate gaps linearly
    idx_finite <- which(is.finite(full_vec))
    if (length(idx_finite) < 2) return(rep(NA_real_, length(vals)))
    
    full_vec_interp <- approx(x = idx_finite, y = full_vec[idx_finite], xout = seq_along(full_vec), rule = 2)$y
    
    # Calculate Rolling Range (Max - Min)
    # Using zoo::rollapply
    r_min <- zoo::rollapply(full_vec_interp, width = window, FUN = min, fill = NA, align = "center")
    r_max <- zoo::rollapply(full_vec_interp, width = window, FUN = max, fill = NA, align = "center")
    r_range <- r_max - r_min
    
    res_for_rows <- r_range[pos_map]
    as.numeric(res_for_rows)
  }

  if ("location_id" %in% names(df)) {
    locs <- unique(df$location_id)
    for (loc in locs) {
      rows <- which(df$location_id == loc)
      if (length(rows) < 1) next
      sub_dates <- df$date[rows]
      ord <- order(sub_dates)
      rrows <- rows[ord]
      vals <- df[[index_name]][rrows]
      out_vals <- process_series_range(sub_dates, vals)
      out[rrows] <- out_vals
    }
  } else {
    out <- process_series_range(df$date, df[[index_name]])
  }

  out[!is.finite(out)] <- NA_real_
  out
}

# Helper: normalize variants of the 'no soil' column name to a single
# canonical name 'no soil' (space). This accepts 'no_soil', 'no.soil',
# '.__no soil__' etc and copies the first matching column into 'no soil'
# if present. This makes downstream access with backticks consistent.
normalize_no_soil_col <- function(tbl) {
  if (is.null(tbl) || !is.data.frame(tbl)) return(tbl)
  nm <- names(tbl)
  candidates <- c("no soil", "no_soil", "no.soil", "__no soil__", "__no_soil__", ".__no soil__", ".__no_soil__")
  # Keep the user's exact 'no soil' when present; otherwise pick the first matching candidate
  # If present, coerce the column to safe numeric form (handles factors/strings)
  if ("no soil" %in% nm) {
    tbl[["no soil"]] <- safe_as_numeric(tbl[["no soil"]])
    return(tbl)
  }
  found <- intersect(candidates, nm)
  if (length(found) > 0) {
    src <- found[1]
    tbl[["no soil"]] <- safe_as_numeric(tbl[[src]])
  }
  tbl
}


# Helper: safely coerce values to numeric. Handles factors, character "TRUE"/"FALSE",
# trims whitespace and avoids direct as.numeric(factor) pitfalls.
safe_as_numeric <- function(x) {
  if (is.null(x)) return(x)
  # Factors -> characters first
  if (is.factor(x)) x <- as.character(x)
  # Characters: trim, handle logical strings
  if (is.character(x)) {
    s <- trimws(x)
    lower <- tolower(s)
    # Map boolean-like strings to 1/0
    lower[lower %in% c("true", "t")] <- "1"
    lower[lower %in% c("false", "f")] <- "0"
    # Attempt numeric conversion
    suppressWarnings(num <- as.numeric(lower))
    return(num)
  }
  # Already numeric (or integer) -> coerce to numeric
  if (is.numeric(x)) return(as.numeric(x))
  # Otherwise try character roundtrip
  suppressWarnings(num <- as.numeric(as.character(x)))
  num
}

# Helper: canonical location_id creation from lon/lat (centralize formatting)
make_location_id <- function(lon, lat) {
  if (length(lon) == 1 && length(lat) == 1) {
    lon_num <- as.numeric(lon); lat_num <- as.numeric(lat)
    if (!is.finite(lon_num) || !is.finite(lat_num)) return(NA_character_)
    sprintf("L_%0.4f_%0.4f", round(lon_num, 4), round(lat_num, 4))
  } else {
    # vectorized form: return character vector
    sapply(seq_along(lon), function(i) {
      lon_i <- as.numeric(lon[i]); lat_i <- as.numeric(lat[i])
      if (!is.finite(lon_i) || !is.finite(lat_i)) NA_character_ else sprintf("L_%0.4f_%0.4f", round(lon_i,4), round(lat_i,4))
    }, USE.NAMES = FALSE)
  }
}

# Persistent parallel backend (set up once)
setup_parallel_backend <- function() {
  if (isTRUE(PARALLEL_ENABLE) && requireNamespace("future", quietly = TRUE)) {
    old_plan <- future::plan()
    options(future.globals.maxSize = max(getOption("future.globals.maxSize", 2e9), 12e9))
    future::plan(future::multisession, workers = PARALLEL_WORKERS, gc = TRUE, earlySignal = TRUE)
    return(function() {
      try(future::plan(old_plan), silent = TRUE)
    })
  }
  function() {}
}

# User-configurable defaults
OUTPUT_DIR <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/phenology_results"
OUT_DIR <- file.path(OUTPUT_DIR, "veg_mixture_fit")

# Allow skipping a minimal number of DOYs per-location when computing sufficiency
# This is useful for small gaps — set to 0 to preserve the strict 50% requirement
MIN_SKIP_DOYS_PER_LOCATION <- 2L

# General algorithm toggles
FAST_VAR <- TRUE
TEMPORAL_BUDGET <- 365L  # Full year resolution (no compression)
TOPK_VARIANTS <- 5L
N_VARIANTS_PER_VEG <- 7L
MIN_CLUSTER_SIZE <- 10L
PCA_VARIANCE_THRESHOLD <- 0.95
LDA_WEIGHT_FLOOR <- 0.01
STAGE1_INDEX_WEIGHT_THRESHOLD <- 0.02
ENABLE_PHASE_ALIGNMENT <- FALSE

# Use raw spectral indices directly (no dimensionality reduction)
SKIP_MOVING_VARIANCE <- TRUE
REFERENCE_PHASE_MARKERS <- c(1, 90, 180, 270, 365)
ENABLE_MULTISCALE <- FALSE
MULTISCALE_WINDOWS <- c(7L, 14L, 30L)
ENABLE_QP_SOLVER <- TRUE
COMBO_PARALLEL_ENABLE <- FALSE
EARLY_STOP_RMSE_THRESHOLD <- 0.0
ENABLE_DIAGNOSTICS <- TRUE

# =============================================================================
# GEOMETRIC MESMA IMPLEMENTATION (Following Tits et al. paper)
# Uses angle-based endmember selection + projection-based unmixing (geometric method)
# =============================================================================

#' Geometric projection-based unmixing for two endmembers
#' This is the core of the paper's approach (Figure 1)
#' 
#' @param y Mixed pixel spectrum (vector)
#' @param m1 First endmember spectrum (vector)
#' @param m2 Second endmember spectrum (vector)
#' @return List with f1, f2 (fractions), y_proj (projection), residual
geometric_unmix_two_endmembers <- function(y, m1, m2) {
  y <- as.numeric(y)
  m1 <- as.numeric(m1)
  m2 <- as.numeric(m2)
  
  # The EM-line direction vector
  em_line <- m2 - m1
  em_norm_sq <- sum(em_line^2)
  
  # Degenerate case: identical endmembers
  if (em_norm_sq < 1e-10) {
    return(list(f1 = 0.5, f2 = 0.5, y_proj = m1, residual = sqrt(sum((y - m1)^2)), t = 0.5))
  }
  
  # Project y onto the EM-line
  # t is the parameter: y' = m1 + t * (m2 - m1)
  # t = 0 means y' = m1, t = 1 means y' = m2
  y_minus_m1 <- y - m1
  t_unclamped <- sum(y_minus_m1 * em_line) / em_norm_sq
  
  # Clamp t to [0, 1] for valid fractions (fully constrained)
  t_clamped <- max(0, min(1, t_unclamped))
  
  # Projected point on EM-line
  y_proj <- m1 + t_clamped * em_line
  
  # Residual (reconstruction error)
  residual <- sqrt(sum((y - y_proj)^2))
  
  # Fractions: f2 = t (fraction of m2), f1 = 1 - t (fraction of m1)
  f2 <- t_clamped
  f1 <- 1 - t_clamped
  
  list(f1 = f1, f2 = f2, y_proj = y_proj, residual = residual, t = t_unclamped)
}


#' Calculate angle between target-line and EM-line (Equation 5 in paper)
#' 
#' @param y Mixed pixel spectrum
#' @param m1 First endmember (anchor point)
#' @param m2 Second endmember (candidate partner)
#' @return Cosine of angle alpha (higher = better fit)
calculate_em_angle_cosine <- function(y, m1, m2) {
  # Target-line: y - m1
  target_line <- y - m1
  
  # EM-line: m2 - m1
  em_line <- m2 - m1
  
  # Norms
  norm_target <- sqrt(sum(target_line^2))
  norm_em <- sqrt(sum(em_line^2))
  
  # Handle degenerate cases
  if (norm_target < 1e-10 || norm_em < 1e-10) {
    return(-1)  # Invalid angle
  }
  
  # Cosine of angle (Equation 5)
  cos_alpha <- sum(target_line * em_line) / (norm_target * norm_em)
  
  # Clamp to valid range
  cos_alpha <- max(-1, min(1, cos_alpha))
  
  return(cos_alpha)
}


#' Select best partner endmember from library using angle criterion (Section 3.2)
#' 
#' For a given m1, find the m2 from library M2 that minimizes the angle
#' between the target-line (y-m1) and the EM-line (m2-m1)
#' 
#' @param y Mixed pixel spectrum
#' @param m1 Fixed first endmember
#' @param M2_library List of candidate second endmembers (each with $vec and $id)
#' @return Best m2 candidate and associated metrics
select_best_partner_by_angle <- function(y, m1, M2_library) {
  if (length(M2_library) == 0) return(NULL)
  
  best_cos_alpha <- -Inf
  best_m2 <- NULL
  best_idx <- NA
  
  for (i in seq_along(M2_library)) {
    m2_candidate <- M2_library[[i]]
    m2_vec <- if (is.list(m2_candidate) && !is.null(m2_candidate$vec)) {
      m2_candidate$vec
    } else {
      m2_candidate
    }
    
    cos_alpha <- calculate_em_angle_cosine(y, m1, m2_vec)
    
    if (cos_alpha > best_cos_alpha) {
      # Verify projection is valid (between m1 and m2)
      unmix_result <- geometric_unmix_two_endmembers(y, m1, m2_vec)
      
      # Only accept if t is in reasonable range (projection on or near EM-line)
      if (unmix_result$t >= -0.1 && unmix_result$t <= 1.1) {
        best_cos_alpha <- cos_alpha
        best_m2 <- m2_candidate
        best_idx <- i
      }
    }
  }
  
  if (is.null(best_m2)) {
    # Fallback: just pick the one with smallest residual
    best_residual <- Inf
    for (i in seq_along(M2_library)) {
      m2_candidate <- M2_library[[i]]
      m2_vec <- if (is.list(m2_candidate) && !is.null(m2_candidate$vec)) {
        m2_candidate$vec
      } else {
        m2_candidate
      }
      unmix_result <- geometric_unmix_two_endmembers(y, m1, m2_vec)
      if (unmix_result$residual < best_residual) {
        best_residual <- unmix_result$residual
        best_m2 <- m2_candidate
        best_idx <- i
      }
    }
  }
  
  list(m2 = best_m2, index = best_idx, cos_alpha = best_cos_alpha)
}


#' Select best endmember pair from two libraries (Algorithm from Section 3.2)
#' 
#' This implements the paper's angle-based selection:
#' 1. For each m1 in M1, find best m2 from M2 using angle criterion
#' 2. Calculate projection distance for each (m1, best_m2) pair
#' 3. Select the pair with smallest projection distance
#' 
#' @param y Mixed pixel spectrum
#' @param M1_library List of first endmember candidates
#' @param M2_library List of second endmember candidates
#' @return Best endmember pair and unmixing result
select_best_em_pair_geometric <- function(y, M1_library, M2_library) {
  if (length(M1_library) == 0) return(NULL)

  if (length(M2_library) == 0) return(NULL)
  
  best_residual <- Inf
  best_result <- NULL
  best_m1 <- NULL
  best_m2 <- NULL
  
  # Iterate over M1 (the smaller library for efficiency)
  for (i in seq_along(M1_library)) {
    m1_candidate <- M1_library[[i]]
    m1_vec <- if (is.list(m1_candidate) && !is.null(m1_candidate$vec)) {
      m1_candidate$vec
    } else {
      m1_candidate
    }
    m1_id <- if (is.list(m1_candidate) && !is.null(m1_candidate$id)) {
      m1_candidate$id
    } else {
      paste0("m1_", i)
    }
    
    # Find best partner from M2 using angle criterion
    partner_result <- select_best_partner_by_angle(y, m1_vec, M2_library)
    
    if (is.null(partner_result$m2)) next
    
    m2_vec <- if (is.list(partner_result$m2) && !is.null(partner_result$m2$vec)) {
      partner_result$m2$vec
    } else {
      partner_result$m2
    }
    m2_id <- if (is.list(partner_result$m2) && !is.null(partner_result$m2$id)) {
      partner_result$m2$id
    } else {
      paste0("m2_", partner_result$index)
    }
    
    # Compute unmixing with this pair
    unmix_result <- geometric_unmix_two_endmembers(y, m1_vec, m2_vec)
    
    # Select based on projection distance (residual)
    if (unmix_result$residual < best_residual) {
      best_residual <- unmix_result$residual
      best_result <- unmix_result
      best_m1 <- list(vec = m1_vec, id = m1_id)
      best_m2 <- list(vec = m2_vec, id = m2_id)
    }
  }
  
  if (is.null(best_result)) return(NULL)
  
  list(
    m1 = best_m1,
    m2 = best_m2,
    f1 = best_result$f1,
    f2 = best_result$f2,
    residual = best_result$residual,
    y_proj = best_result$y_proj
  )
}


#' Project vector onto probability simplex (sum=1, all >= 0)
#' Uses the algorithm from Duchi et al. (2008)
project_to_simplex <- function(v) {
  n <- length(v)
  if (n == 0) return(numeric(0))
  if (n == 1) return(1)
  
  # Sort in descending order
  u <- sort(v, decreasing = TRUE)
  
  # Find rho
  cssv <- cumsum(u)
  rho <- max(which(u + (1 - cssv) / seq_along(u) > 0))
  
  # Compute theta
  theta <- (cssv[rho] - 1) / rho
  
  # Project
  w <- pmax(v - theta, 0)
  
  return(w)
}


#' Geometric unmixing for N endmembers using constrained projection
#' Based on extending the 2-EM case to higher dimensions
#' 
#' @param y Mixed pixel spectrum
#' @param M Matrix of endmembers (columns are endmembers)
#' @return List with fractions and residual
geometric_unmix_n_endmembers <- function(y, M) {
  y <- as.numeric(y)
  N <- ncol(M)  # Number of endmembers
  P <- nrow(M)  # Number of bands/features
  
  if (N == 1) {
    residual <- sqrt(sum((y - M[, 1])^2))
    return(list(f = 1, residual = residual, y_proj = M[, 1]))
  }
  
  if (N == 2) {
    result <- geometric_unmix_two_endmembers(y, M[, 1], M[, 2])
    return(list(f = c(result$f1, result$f2), residual = result$residual, y_proj = result$y_proj))
  }
  
  # For N > 2: Use constrained least squares with simplex projection
  
  # Try quadprog for exact constrained solution (NNLS + Sum-to-One)
  if (requireNamespace("quadprog", quietly = TRUE)) {
    # min || y - M*f ||^2  => min f'M'Mf - 2y'Mf
    # Dmat = M'M
    # dvec = M'y
    # Constraints: sum(f)=1 (equality), f>=0 (inequality)
    
    Dmat <- t(M) %*% M
    # Add small ridge for positive definiteness
    Dmat <- Dmat + diag(N) * 1e-8
    dvec <- t(M) %*% y
    
    # Amat: first col is 1s (equality), then identity (inequality)
    Amat <- cbind(rep(1, N), diag(N))
    bvec <- c(1, rep(0, N))
    
    res <- tryCatch({
      quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
    }, error = function(e) NULL)
    
    if (!is.null(res)) {
      f <- res$solution
      f[f < 0] <- 0 # Clamp small negatives
      if (sum(f) > 0) f <- f / sum(f) # Re-normalize
      y_proj <- as.numeric(M %*% f)
      residual <- sqrt(sum((y - y_proj)^2))
      return(list(f = f, residual = residual, y_proj = y_proj))
    }
  }

  # Fallback: Solve unconstrained least squares: f = (M'M)^(-1) M'y
  MtM <- t(M) %*% M
  Mty <- t(M) %*% y
  
  # Add small ridge for numerical stability
  ridge <- 1e-8 * diag(N)
  
  f_unconstrained <- tryCatch({
    solve(MtM + ridge, Mty)
  }, error = function(e) {
    # Fallback to pseudoinverse
    if (requireNamespace("MASS", quietly = TRUE)) {
      MASS::ginv(MtM + ridge) %*% Mty
    } else {
      rep(1/N, N)  # Equal fractions fallback
    }
  })
  
  # Project onto probability simplex (sum=1, all >= 0)
  f <- project_to_simplex(as.numeric(f_unconstrained))
  
  # Compute projected point and residual
  y_proj <- as.numeric(M %*% f)
  residual <- sqrt(sum((y - y_proj)^2))
  
  list(f = f, residual = residual, y_proj = y_proj)
}


#' Compute LDA-based feature weights for ranking/weighting indices
#' 
#' @param mesma_lib The constructed MESMA library (list of variants)
#' @param compressed_templates_accessor Accessor for compressed feature vectors
#' @param grid_type The temporal grid type to use (default "full")
#' @return A numeric vector of feature weights, or NULL if LDA fails
compute_lda_weights <- function(mesma_lib, compressed_templates_accessor, grid_type = "full") {
  cat("\n=== COMPUTING LDA FEATURE WEIGHTS ===\n")
  
  # Gather training data (variants as samples)
  X_list <- list()
  y_list <- list()
  
  for (veg in names(mesma_lib)) {
    variants <- mesma_lib[[veg]]
    for (variant in variants) {
      vid <- variant$variant_id
      # Get the compressed feature vector
      vec <- tryCatch(compressed_templates_accessor[[veg]][[vid]][[grid_type]], error = function(e) NULL)
      
      if (!is.null(vec)) {
        X_list[[length(X_list) + 1]] <- vec
        y_list[[length(y_list) + 1]] <- veg
      }
    }
  }
  
  if (length(X_list) == 0) return(NULL)
  
  X <- do.call(rbind, X_list)
  y <- unlist(y_list)
  
  # Check if we have enough samples per class
  counts <- table(y)
  if (any(counts < 2) || length(unique(y)) < 2) {
    cat("LDA skipped: insufficient classes or samples per class.\n")
    return(NULL)
  }
  
  # Remove constant columns
  sds <- apply(X, 2, sd)
  keep_cols <- sds > 1e-6
  if (sum(keep_cols) < 2) {
    cat("LDA skipped: too few variable features.\n")
    return(NULL)
  }
  
  X_clean <- X[, keep_cols, drop = FALSE]
  
  # --- PCA Step before LDA ---
  # Perform PCA to handle multicollinearity and reduce dimensionality
  # We use scale. = TRUE to standardize features (indices have different ranges)
  pca_res <- tryCatch({
    prcomp(X_clean, center = TRUE, scale. = TRUE)
  }, error = function(e) {
    cat(sprintf("PCA failed: %s\n", e$message))
    NULL
  })
  
  if (is.null(pca_res)) return(NULL)
  
  # Keep components explaining 99% of variance to remove noise/singularity
  # but retain almost all information (including unique variance)
  cum_var <- cumsum(pca_res$sdev^2) / sum(pca_res$sdev^2)
  n_pcs <- which(cum_var >= 0.99)[1]
  if (is.na(n_pcs)) n_pcs <- length(cum_var)
  
  # Use at least 2 PCs if possible, or all if fewer
  n_pcs <- max(min(n_pcs, ncol(X_clean)), 2)
  if (n_pcs > ncol(X_clean)) n_pcs <- ncol(X_clean)
  
  X_pca <- pca_res$x[, 1:n_pcs, drop = FALSE]
  cat(sprintf("LDA Pre-processing: Reduced %d features to %d PCs (%.1f%% variance)\n", 
              ncol(X_clean), n_pcs, cum_var[n_pcs]*100))
  
  # Run LDA on PCA components
  lda_res <- tryCatch({
    MASS::lda(X_pca, grouping = y)
  }, error = function(e) {
    cat(sprintf("LDA failed: %s\n", e$message))
    NULL
  })
  
  if (is.null(lda_res)) return(NULL)
  
  # Compute weights from loadings
  # lda_res$scaling contains the coefficients for the PCs (W_pc)
  # pca_res$rotation contains the loadings of original features on PCs (R)
  # We project LDA weights back to original feature space: W_std = R * W_pc
  
  W_pc <- lda_res$scaling
  R <- pca_res$rotation[, 1:n_pcs, drop = FALSE]
  
  # Effective weights in standardized feature space
  W_std <- R %*% W_pc
  
  svd <- lda_res$svd
  prop <- svd^2 / sum(svd^2)
  
  if (ncol(W_std) > 1) {
    # Weighted sum of absolute contributions across discriminant axes
    # Ensure prop matches dimensions of W_std columns
    n_dim <- min(length(prop), ncol(W_std))
    weights_clean <- rowSums(abs(W_std[, 1:n_dim, drop=FALSE]) %*% diag(prop[1:n_dim], nrow=n_dim))
  } else {
    weights_clean <- abs(W_std[, 1])
  }
  
  weights <- rep(0, ncol(X))
  weights[keep_cols] <- weights_clean
  
  # Normalize to [0, 1]
  weights <- weights / max(weights)
  
  # --- Per-Index Thresholding ---
  # Determine number of indices (K) from the first variant to aggregate weights per index
  first_veg <- names(mesma_lib)[1]
  if (!is.null(first_veg) && length(mesma_lib[[first_veg]]) > 0) {
    # Try to get K from raw_mat of the first variant
    raw_mat <- mesma_lib[[first_veg]][[1]]$raw_mat
    if (!is.null(raw_mat)) {
      K <- ncol(raw_mat)
      if (K > 0) {
        # Aggregate weights per index (Max across time steps)
        index_max_weights <- numeric(K)
        for (k in 1:K) {
          idx_indices <- seq(k, length(weights), by = K)
          idx_indices <- idx_indices[idx_indices <= length(weights)]
          if (length(idx_indices) > 0) {
             index_max_weights[k] <- max(weights[idx_indices], na.rm = TRUE)
          }
        }
        
        # Identify indices to skip (< STAGE1_INDEX_WEIGHT_THRESHOLD relative weight)
        skip_indices <- which(index_max_weights < STAGE1_INDEX_WEIGHT_THRESHOLD)
        
        if (length(skip_indices) > 0) {
          cat(sprintf("LDA Thresholding: Skipping %d indices with max relative weight < %.3f (indices: %s)\n", 
                      length(skip_indices), STAGE1_INDEX_WEIGHT_THRESHOLD, paste(skip_indices, collapse=", ")))
          
          # Zero out ALL features for these indices
          for (k in skip_indices) {
            idx_indices <- seq(k, length(weights), by = K)
            idx_indices <- idx_indices[idx_indices <= length(weights)]
            weights[idx_indices] <- 0
          }
        }
      }
    }
  }

  # Thresholding: set weights below 1% of max to zero to remove non-contributing features
  weights[weights < 0.01] <- 0
  
  # Re-normalize non-zero weights
  if (max(weights) > 0) {
    weights <- weights / max(weights)
  }
  
  cat("LDA weights computed.\n")
  return(weights)
}


#' Stage 2 geometric unmixing: Separate vegetation types
#' Uses angle-based selection from the paper
#' 
#' @param y Observation vector (in feature space)
#' @param veg_libraries Named list of vegetation endmember libraries
#' @param topK Number of top candidates to consider per vegetation type
#' @param feature_weights Optional vector of weights to apply to features (LDA-based)
#' @return Unmixing result with fractions per vegetation type
stage2_geometric_unmix <- function(y, veg_libraries, topK = 2, feature_weights = NULL) {
  veg_types <- names(veg_libraries)
  n_veg <- length(veg_types)
  
  if (n_veg == 0) return(NULL)
  
  y <- as.numeric(y)
  
  # Prepare weighted versions for selection if weights are provided
  y_select <- y
  veg_libraries_select <- veg_libraries
  
  if (!is.null(feature_weights)) {
    if (length(feature_weights) != length(y)) {
      # Try to recycle or warn
      if (length(feature_weights) > 0) {
         warning(sprintf("stage2_geometric_unmix: feature_weights length (%d) mismatch with y (%d), ignoring weights", length(feature_weights), length(y)))
      }
    } else {
      # Apply weights to y_select
      y_select <- y * feature_weights
      
      # Apply weights to all endmembers in veg_libraries_select (deep copy to avoid side effects)
      new_libs <- list()
      for (veg in veg_types) {
        lib <- veg_libraries[[veg]]
        new_lib <- list()
        for (i in seq_along(lib)) {
          em <- lib[[i]]
          if (is.list(em) && !is.null(em$vec)) {
            em_copy <- em
            em_copy$vec <- em$vec * feature_weights
            new_lib[[i]] <- em_copy
          } else {
            new_lib[[i]] <- em * feature_weights
          }
        }
        new_libs[[veg]] <- new_lib
      }
      veg_libraries_select <- new_libs
    }
  }
  
  # Step 1: For each vegetation type, select top-K candidates by distance
  top_candidates <- list()
  
  for (veg in veg_types) {
    lib <- veg_libraries[[veg]] # Original (unweighted)
    lib_sel <- veg_libraries_select[[veg]] # Weighted (or original)
    
    if (length(lib) == 0) next
    
    # Calculate distance from y_select to each endmember in lib_sel
    scores <- sapply(lib_sel, function(em) {
      em_vec <- if (is.list(em) && !is.null(em$vec)) em$vec else em
      # Use negative distance as score (higher = closer = better)
      -sqrt(sum((y_select - em_vec)^2))
    })
    
    # Select top-K
    n_select <- min(topK, length(lib))
    top_idx <- order(scores, decreasing = TRUE)[1:n_select]
    top_candidates[[veg]] <- lib[top_idx]
  }
  
  if (length(top_candidates) == 0) return(NULL)
  
  # Step 2: If only one veg type, simple case
  if (length(top_candidates) == 1) {
    veg <- names(top_candidates)[1]
    best_em <- top_candidates[[veg]][[1]]
    em_vec <- if (is.list(best_em) && !is.null(best_em$vec)) best_em$vec else best_em
    em_id <- if (is.list(best_em) && !is.null(best_em$id)) best_em$id else paste0(veg, "_1")
    
    residual <- sqrt(sum((y - em_vec)^2))
    return(list(
      fractions = setNames(1.0, veg),
      chosen_variants = setNames(em_id, veg),
      residual = residual
    ))
  }
  
  # Step 3: For two veg types, use the paper's angle-based pair selection
  if (length(top_candidates) == 2) {
    veg_names <- names(top_candidates)
    M1_lib <- top_candidates[[veg_names[1]]]
    M2_lib <- top_candidates[[veg_names[2]]]
    
    # Use angle-based selection
    pair_result <- select_best_em_pair_geometric(y, M1_lib, M2_lib)
    
    if (is.null(pair_result)) return(NULL)
    
    fractions <- c(pair_result$f1, pair_result$f2)
    names(fractions) <- veg_names
    
    chosen <- c(pair_result$m1$id, pair_result$m2$id)
    names(chosen) <- veg_names
    
    return(list(
      fractions = fractions,
      chosen_variants = chosen,
      residual = pair_result$residual
    ))
  }
  
  # Step 4: For 3+ veg types, use combinatorial or iterative approach
  
  # Calculate total combinations
  lib_sizes <- sapply(top_candidates, length)
  total_combos <- prod(lib_sizes)
  MAX_COMBOS_FOR_FULL_SEARCH <- 5000
  
  selected_ems <- NULL
  
  if (total_combos <= MAX_COMBOS_FOR_FULL_SEARCH) {
    # Full combinatorial search
    # Generate all combinations of indices
    # Note: expand.grid columns will correspond to the order of veg_types in top_candidates
    combo_indices <- expand.grid(lapply(lib_sizes, seq_len))
    
    best_residual <- Inf
    best_combo_idx <- 1
    
    for (i in seq_len(nrow(combo_indices))) {
      # Construct M for this combination
      M_cols <- list()
      
      for (j in seq_along(veg_types)) {
        veg <- veg_types[j]
        idx <- combo_indices[i, j]
        em <- top_candidates[[veg]][[idx]]
        M_cols[[j]] <- if (is.list(em) && !is.null(em$vec)) em$vec else em
      }
      
      M <- do.call(cbind, M_cols)
      unmix_result <- geometric_unmix_n_endmembers(y, M)
      
      if (unmix_result$residual < best_residual) {
        best_residual <- unmix_result$residual
        best_combo_idx <- i
      }
    }
    
    # Reconstruct best result
    best_indices <- combo_indices[best_combo_idx, ]
    selected_ems <- list()
    for (j in seq_along(veg_types)) {
      veg <- veg_types[j]
      idx <- best_indices[[j]]
      selected_ems[[veg]] <- top_candidates[[veg]][[idx]]
    }
    
  } else {
    # Fallback to iterative greedy approach for large search spaces
    # Order libraries by size (smallest first for efficiency)
    lib_order <- names(sort(lib_sizes))
    
    # Start with first two libraries
    v1 <- lib_order[1]
    v2 <- lib_order[2]
    
    pair_result <- select_best_em_pair_geometric(y, top_candidates[[v1]], top_candidates[[v2]])
    
    if (is.null(pair_result)) return(NULL)
    
    selected_ems <- list()
    selected_ems[[v1]] <- pair_result$m1
    selected_ems[[v2]] <- pair_result$m2
    
    # Iteratively add remaining vegetation types
    for (k in seq(3, length(lib_order))) {
      vk <- lib_order[k]
      
      best_residual <- Inf
      best_em_k <- NULL
      
      # Try each candidate for vk
      for (em_candidate in top_candidates[[vk]]) {
        em_vec <- if (is.list(em_candidate) && !is.null(em_candidate$vec)) {
          em_candidate$vec
        } else {
          em_candidate
        }
        
        # Build endmember matrix with current selection + candidate
        M_cols <- lapply(selected_ems, function(em) {
          if (is.list(em) && !is.null(em$vec)) em$vec else em
        })
        M_cols[[length(M_cols) + 1]] <- em_vec
        M <- do.call(cbind, M_cols)
        
        # Unmix with all endmembers
        unmix_result <- geometric_unmix_n_endmembers(y, M)
        
        if (unmix_result$residual < best_residual) {
          best_residual <- unmix_result$residual
          best_em_k <- em_candidate
        }
      }
      
      if (!is.null(best_em_k)) {
        selected_ems[[vk]] <- best_em_k
      }
    }
  }
  
  # Final unmixing with all selected endmembers
  M_final_cols <- lapply(selected_ems, function(em) {
    if (is.list(em) && !is.null(em$vec)) em$vec else em
  })
  M_final <- do.call(cbind, M_final_cols)
  colnames(M_final) <- names(selected_ems)
  
  final_result <- geometric_unmix_n_endmembers(y, M_final)
  
  fractions <- final_result$f
  names(fractions) <- names(selected_ems)
  
  chosen <- sapply(selected_ems, function(em) {
    if (is.list(em) && !is.null(em$id)) em$id else "unknown"
  })
  
  list(
    fractions = fractions,
    chosen_variants = chosen,
    residual = final_result$residual
  )
}

# =============================================================================
# END GEOMETRIC MESMA CORE FUNCTIONS
# =============================================================================

# Combination expansion safety thresholds
COMBO_SAFE_EXPAND_LIMIT <- 1e6    # fully expand grid up to this many combos
COMBO_ABORT_LIMIT <- 5e7          # abort if combos exceed this hard limit

# Bootstrap settings 
BOOTSTRAP_B <- 200L
TESTING_MODE <- FALSE
ENABLE_UNCERTAINTY <- TRUE
DEBUG_UNCERTAINTY <- FALSE  # global debug flag for uncertainty diagnostics
# NOTE: Bootstrapping is now "nested two-stage" only when `ENABLE_UNCERTAINTY = TRUE`.
# Other bootstrap mechanisms (geometric_block_bootstrap, ols_block_bootstrap) are disabled as fallbacks.
# Ensure `COMPRESSED_STAGE1_LIB` is built before enabling uncertainty.
# Minimum number of observations per location-year required for MESMA tasks
MIN_OBS_PER_LOC_YEAR <- 1L
# Minimum unique DOYs per task used to avoid trivial predictions; default for training/test
MIN_UNIQUE_DOY_DEFAULT <- 5L
# For inference locations, allow fewer DOYs so we don't skip these locations on that basis
MIN_UNIQUE_DOY_INFERENCE <- 1L

# Maximum number of inference locations to include in the separate inference results
# INFERENCE_MAX_LOCATIONS <- 2000L  # Removed limit
# Bootstrap variant-switching control: keep variants fixed by default (recommended).
VARIANT_SWITCH_BOOTSTRAP <- TRUE
VARIANT_SWITCH_RE_CENTER <- TRUE


# Optimization/solver defaults
VARIANCE_THRESHOLD <- 0.90
MAX_VEG_COMPONENTS <- 8
GAM_K_MAX <- 40
GAM_GAMMA <- 1.0


# Index selection and prefiltering
USE_INDICES_MIN <- 1L
MIN_INDEX_SD <- 0.05

# Sample-balancing and augmentation
ENABLE_SAMPLE_BALANCING <- FALSE

# Soil-preprocessing: run a preliminary soil/no-soil MESMA and subtract soil
# contribution from all spectral time series before building vegetation libs.
# This helps avoid soil contamination when clustering vegetation time-series.
ENABLE_SOIL_PREPROCESS <- TRUE
# Threshold (0-1) to treat a geo-labelled sample as pure soil when building the soil prototype
SOIL_PURE_THRESHOLD <- 0.95
# Minimum number of pure-soil rows required to build a robust soil prototype
SOIL_MIN_SAMPLES <- 3L
# When insufficient pure soil points are present, the soil prototype will not
# be built (returns NULL) unless a minimum number of candidate rows meet
# the SOIL_PURE_THRESHOLD.

# Memory-safe clustering limits (sampling parameters to avoid RAM overload)
MAX_PROJECTIONS_PER_VEG <- 25000L  # subsample before clustering to avoid OOM
SILHOUETTE_SAMPLE_SIZE <- 20000L   # subsample for silhouette distance matrix
MEDOID_SAMPLE_SIZE <- 10000L       # subsample for medoid distance computation

# Numeric safety constants
EPS_SIGMA <- 1e-8
LOWER_BND <- 0

# DOY / index presence tuning
MIN_IDX_PRESENCE <- 0.5

# Strict validation: enforce extractor/predictor output format only
if (!file.exists(INPUT_CSV)) stop(paste0("Required input CSV not found: ", INPUT_CSV))
raw_df <- tryCatch(
  {
    readr::read_csv(INPUT_CSV, show_col_types = FALSE)
  },
  error = function(e) stop(paste0("Failed to read INPUT_CSV: ", e$message))
)
if (nrow(raw_df) < 2) stop("INPUT_CSV contains fewer than 2 rows")
if (!"location_id" %in% names(raw_df) || (!"prediction_date" %in% names(raw_df) && !"date" %in% names(raw_df))) {
  stop("INPUT_CSV must contain 'location_id' and 'prediction_date' (or 'date') columns")
}

# Deduplicate observations based on location_id and date
date_col <- if ("prediction_date" %in% names(raw_df)) "prediction_date" else "date"
raw_df <- raw_df %>% distinct(location_id, !!sym(date_col), .keep_all = TRUE)
cat(sprintf("After deduplication: %d rows remaining from original %d rows.\n", nrow(raw_df), nrow(raw_df)))



df <- raw_df
df <- normalize_band_names(df)

## Index scale detection and correction (runs after linearization)
detect_and_rescale_indices <- function(df, cols = OPTIMAL_INDICES) {
  present <- intersect(cols, names(df))
  if (length(present) == 0) return(list(df = df, BAND_SCALE = 1))

  # Use absolute 99th percentile of indices to detect unusually large magnitudes
  idx99 <- sapply(present, function(c) {
    v <- as.numeric(df[[c]])
    if (all(is.na(v))) return(NA_real_)
    quantile(abs(v), probs = 0.99, na.rm = TRUE)
  })
  global_max <- max(idx99, na.rm = TRUE)
  if (!is.finite(global_max)) global_max <- 0

  if (global_max > 2) {
    # heuristics for scale inference: >1000 -> 10000, >10 -> 100
    scale_factor <- if (global_max > 1000) 10000 else if (global_max > 10) 100 else 1
    if (scale_factor > 1) {
      for (c in present) df[[c]] <- as.numeric(df[[c]]) / scale_factor
      cat(sprintf("[NOTICE] Detected index magnitudes with 99%%ile ~ %f; rescaled indices by 1/%d to expected range.\n", global_max, scale_factor))
      return(list(df = df, BAND_SCALE = 1/scale_factor))
    }
  }
  return(list(df = df, BAND_SCALE = 1))
}

# --- PPI Calculation Start ---
# Ensure date is Date object for DOY calculation
if (!"date" %in% names(df) && "prediction_date" %in% names(df)) df$date <- as.Date(df$prediction_date)
if ("date" %in% names(df)) df$date <- as.Date(df$date)

cat("Checking/Calculating PPI...\n")

if (SKIP_PPI) {
  cat("Skipping PPI calculation\n")
  df$PPI <- NA_real_
  df$zenith.angle <- NA_real_
} else {
  df <- add_ppi_columns(df, 40.2)
  cat("Calculated PPI using add_ppi_columns with lat=40.2\n")
}
# --- PPI Calculation End ---

# Ensure PPI-related columns exist
if (!"zenith.angle" %in% names(df)) df$zenith.angle <- NA_real_
if (!"PPI" %in% names(df)) df$PPI <- NA_real_

# Apply linearization to indices
df <- linearize_indices(df)

# Detect and rescale indices after linearization if they appear scaled (e.g. inputs in 0-10000)
idx_scale_res <- detect_and_rescale_indices(df, cols = OPTIMAL_INDICES)
df <- idx_scale_res$df
BAND_SCALE <- idx_scale_res$BAND_SCALE
if (!exists("BAND_SCALE") || is.null(BAND_SCALE)) BAND_SCALE <- 1

# Normalize indices after linearization to robust [-1..1] range per-index to avoid overly scaled indices
normalize_indices_after_linearization <- function(df, cols = OPTIMAL_INDICES) {
  present <- intersect(cols, names(df))
  idx_scales <- list()
  for (i in present) {
    colv <- df[[i]]
    if (!is.numeric(colv) || all(is.na(colv))) { idx_scales[[i]] <- 1; next }
    scale_val <- as.numeric(stats::quantile(abs(colv), probs = 0.99, na.rm = TRUE))
    if (!is.finite(scale_val) || scale_val <= 0) scale_val <- 1
    df[[i]] <- df[[i]] / scale_val
    idx_scales[[i]] <- scale_val
  }
  attr(df, 'INDEX_SCALES') <- idx_scales
  df
}

# Recompute PPI from scaled DVI values (in case DVI was rescaled post-linearization)
recompute_ppi_from_scaled_dvi <- function(df, lat_default = 40.2) {
  if (!"DVI" %in% names(df) || !"date" %in% names(df)) return(df)
  df$date <- as.Date(df$date)
  if (!"doy" %in% names(df)) df$doy <- lubridate::yday(df$date)
  if (!"year" %in% names(df)) df$year <- lubridate::year(df$date)

  peak_df <- df %>%
    group_by(location_id, year) %>%
    summarise(
      DVI_max = if (all(is.na(DVI))) NA_real_ else max(DVI, na.rm = TRUE),
      doy_peak = { if (all(is.na(DVI))) NA_integer_ else { idx <- which.max(DVI); doy[idx[1]] } },
      .groups = "drop"
    )
  peak_df$DVI_max[!is.finite(peak_df$DVI_max)] <- NA_real_
  peak_df$doy_peak[!is.finite(peak_df$doy_peak)] <- NA_integer_

  if ("lat" %in% names(df)) {
    lat_lookup <- df %>% group_by(location_id) %>% summarise(lat_use = mean(lat, na.rm = TRUE), .groups = "drop")
  } else {
    lat_lookup <- data.frame(location_id = unique(df$location_id), lat_use = lat_default)
  }
  peak_df <- peak_df %>% left_join(lat_lookup, by = "location_id")
  peak_df$lat_use[!is.finite(peak_df$lat_use)] <- lat_default

  peak_df$zenith.angle <- NA_real_
  zen_idx <- complete.cases(peak_df$lat_use, peak_df$doy_peak)
  if (any(zen_idx)) {
    peak_df$zenith.angle[zen_idx] <- calculate_solar_zenith(lat = peak_df$lat_use[zen_idx], doy = peak_df$doy_peak[zen_idx])
  }

  peak_df$PPI <- NA_real_
  ppi_idx <- complete.cases(peak_df$DVI_max, peak_df$zenith.angle)
  if (any(ppi_idx)) {
    peak_df$PPI[ppi_idx] <- ppi(dvi = peak_df$DVI_max[ppi_idx], zenith.angle = peak_df$zenith.angle[ppi_idx], M = peak_df$DVI_max[ppi_idx] + 0.005, dvi.soil = PPI_DVI_SOIL)
  }

  # Replace existing columns in df
  drop_cols <- intersect(c("DVI_max", "PPI", "zenith.angle"), names(df))
  if (length(drop_cols)) df[drop_cols] <- NULL
  df <- df %>% left_join(peak_df %>% select(location_id, pheno_year, DVI_max, zenith.angle, PPI), by = c("location_id", "pheno_year"))
  as.data.frame(df)
}

# recompute PPI now after possible index rescaling
df <- recompute_ppi_from_scaled_dvi(df)

df <- normalize_indices_after_linearization(df)

# Filter indices based on linearity after linearization
linearity_file <- "index_linearity_scores.csv"
if (file.exists(linearity_file)) {
  linearity_scores <- read.csv(linearity_file)
  # Filter indices with Normalized_Max_Dev <= 0.2
  good_indices <- linearity_scores$Index[linearity_scores$Normalized_Max_Dev <= 0.2]
  cat(sprintf("Filtering indices with linearity > 0.2 after linearization. Retained indices: %s\n", paste(good_indices, collapse = ", ")))
  
  # Identify index columns (assuming they are the ones not in meta columns)
  meta_cols <- c("location_id", "lat", "lon", "imagery_lat", "imagery_lon", "target_lat", "target_lon", "date", "year", "doy", "prediction_date", "DUSTI", "PPI", "zenith.angle", "DVI_max")
  index_cols <- setdiff(names(df), meta_cols)
  # Keep only good indices
  keep_cols <- c(meta_cols, intersect(index_cols, good_indices))
  df <- df[, keep_cols, drop = FALSE]
  
  # Update OPTIMAL_INDICES to the filtered list
  OPTIMAL_INDICES <- good_indices
  if ("PPI" %in% names(df)) OPTIMAL_INDICES <- c(OPTIMAL_INDICES, "PPI")
} else {
  warning("Linearity scores file not found, proceeding without filtering")
}

# Ensure all OPTIMAL_INDICES are present in the CSV (or calculated)
missing_idx <- setdiff(OPTIMAL_INDICES, names(df))
if (length(missing_idx) > 0) {
  stop(paste0("INPUT_CSV missing required indices: ", paste(missing_idx, collapse = ", ")))
}

# Assign valid location_ids as row numbers to replace invalid "L_NA_NA" values
# Replace only the rows that are invalid (NA or explicitly "L_NA_NA") rather than
# overwriting the entire column, and coerce factors to character first.
if ("location_id" %in% names(df)) {
  # ensure it is character for safe comparison and assignment
  df$location_id <- as.character(df$location_id)
  bad_idx <- which(is.na(df$location_id) | df$location_id == "L_NA_NA")
  if (length(bad_idx) > 0L) {
    # Use the row numbers for only those invalid rows
    df$location_id[bad_idx] <- as.character(bad_idx)
    cat(sprintf("[NOTICE] Replaced %d invalid location_id entries with row numbers\n", length(bad_idx)))
  }
}
# Normalize any variants of the no-soil column name for consistent downstream access
df <- normalize_no_soil_col(df)

# Read GeoJSON for vegetation ground truth
if (!"date" %in% names(df) && "Date" %in% names(df)) df$date <- df$Date
if ("date" %in% names(df)) {
  df$date <- as.Date(df$date)
  if (!"year" %in% names(df)) df$year <- as.integer(lubridate::year(df$date))
}

geojson_path <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/updated_zuizer_zonder_foto_UTM.geojson"
if (!file.exists(geojson_path)) stop(paste0("GeoJSON points not found at ", geojson_path))
gpts_raw <- sf::st_read(geojson_path, quiet = TRUE)

geojson_names <- names(gpts_raw)
normalized_names <- gsub("[^a-z0-9]+", "_", tolower(geojson_names))
no_soil_col <- geojson_names[normalized_names == "no_soil"]
if (length(no_soil_col) == 0) stop("GeoJSON point data requires a 'no soil' column")
no_soil_raw <- gpts_raw[[no_soil_col[1]]]
no_soil_vals <- safe_as_numeric(no_soil_raw)
if (is.null(no_soil_vals) || (!is.numeric(no_soil_vals) && !is.logical(no_soil_vals))) {
  stop("'no soil' column must be numeric, logical, or character coercible to numeric")
}
  if (any(is.na(no_soil_vals))) {
  # Leave missing 'no soil' values as NA instead of coercing to 0.
  # Coercing to 0 earlier caused many locations to be treated as fully barren/mixed
  # when GeoJSON contained missing entries. Keeping NA helps us detect missing
  # metadata and prevents silently forcing binary fractions.
  cat("[NOTICE] Found missing/non-numeric 'no soil' values in GeoJSON; keeping as NA (do not coerce to 0).\n")
}
## Validate numeric range only for non-missing values; keep NA values allowed
non_na_no_soil <- no_soil_vals[!is.na(no_soil_vals)]
if (length(non_na_no_soil) > 0 && any(non_na_no_soil < 0 | non_na_no_soil > 1)) stop("'no soil' values must lie within [0,1]")
gpts_raw$`.__no soil__` <- no_soil_vals

matched_cols <- names(gpts_raw)[tolower(names(gpts_raw)) %in% c("vegetation", "veg", "class")]
if (length(matched_cols) > 0) {
  veg_col_orig <- matched_cols[1]
  gpts_raw$.__veg__ <- as.character(gpts_raw[[veg_col_orig]])
} else {
  gpts_raw$.__veg__ <- NA_character_
}
coords <- sf::st_coordinates(gpts_raw)
gpts_raw$.__lon__ <- coords[, 1]
gpts_raw$.__lat__ <- coords[, 2]
# Preserve the coordinate-based id under a separate name
gpts_raw$location_id_geo <- make_location_id(gpts_raw$.__lon__, gpts_raw$.__lat__)
# Use row numbers as location_id to match CSV format
gpts_raw$location_id <- as.character(seq_len(nrow(gpts_raw)))

gpts_map <- sf::st_drop_geometry(gpts_raw) %>%
  dplyr::select(location_id_geo = location_id, Veg = .__veg__, `no soil` = `.__no soil__`) %>%
  dplyr::mutate(location_id = as.character(seq_len(dplyr::n()))) %>%  # Use row number as location_id
  dplyr::distinct(location_id, .keep_all = TRUE)

# DEBUG: Print gpts_map summary including barren
cat("\n=== GPTS_MAP LOADING DEBUG ===\n")
cat(sprintf("Total locations in gpts_map: %d\n", nrow(gpts_map)))
cat(sprintf("Veg class distribution:\n"))
print(table(gpts_map$Veg, useNA = "ifany"))
cat(sprintf("Sample location_ids (first 5): %s\n", paste(head(gpts_map$location_id, 5), collapse=", ")))
barren_gpts_check <- gpts_map[tolower(gpts_map$Veg) == "barren", ]
cat(sprintf("Barren locations: %d\n", nrow(barren_gpts_check)))
if (nrow(barren_gpts_check) > 0) {
  cat(sprintf("  Barren location_ids: %s\n", paste(head(barren_gpts_check$location_id, 10), collapse=", ")))
}
cat("==============================\n\n")

if (nrow(gpts_map) == 0) stop("GeoJSON mapping produced no valid points")

if ("location_id" %in% names(df) && "location_id" %in% names(gpts_map)) {
  # Robust left-join: normalize and coerce types so mismatched formats (whitespace/case)
  # do not prevent the mapping.
  if (!is.character(df$location_id)) df$location_id <- as.character(df$location_id)
  if (!is.character(gpts_map$location_id)) gpts_map$location_id <- as.character(gpts_map$location_id)

  # Trim whitespace and coerce empty->NA for stable joining
  df$location_id <- trimws(df$location_id)
  df$location_id[df$location_id == ""] <- NA_character_
  gpts_map$location_id <- trimws(gpts_map$location_id)
  gpts_map$location_id[gpts_map$location_id == ""] <- NA_character_

  # For safety, normalize simple differences (e.g. lowercase/uppercase 'L_' prefix)
  df$location_id <- ifelse(grepl('^l_', df$location_id, ignore.case = TRUE), sub('^l_', 'L_', df$location_id, ignore.case = TRUE), df$location_id)
  gpts_map$location_id <- ifelse(grepl('^l_', gpts_map$location_id, ignore.case = TRUE), sub('^l_', 'L_', gpts_map$location_id, ignore.case = TRUE), gpts_map$location_id)

  # Try joining in two ways:
  #  1) by 'location_id' (preferable), 2) by row-number (if CSV location_id are numeric row indices)
  # Ensure `Veg` exists so subsequent checks don't raise 'unknown column' warnings
  if (!"Veg" %in% names(df)) df$Veg <- NA_character_
  pre_non_na <- sum(!is.na(df$Veg) & df$Veg != "")

  joined <- dplyr::left_join(df, gpts_map, by = "location_id", suffix = c("", ".geo"))
  # DEBUG: Check barren before/after join
  cat("\n=== BARREN JOIN DEBUG ===\n")
  barren_gpts <- gpts_map[tolower(gpts_map$Veg) == "barren", ]
  cat(sprintf("Barren locations in gpts_map: %d\n", length(unique(barren_gpts$location_id))))
  if (nrow(barren_gpts) > 0) {
    sample_barren_ids <- head(unique(barren_gpts$location_id), 5)
    cat(sprintf("  Sample barren location_ids: %s\n", paste(sample_barren_ids, collapse=", ")))
    # Check if these exist in joined df
    for (bid in sample_barren_ids) {
      matched_rows <- joined[joined$location_id == bid, ]
      cat(sprintf("  Location %s: %d rows in joined, Veg values: %s\n", 
                  bid, nrow(matched_rows), 
                  paste(unique(matched_rows$Veg), collapse=",")))
      if ("Veg.geo" %in% names(matched_rows)) {
        cat(sprintf("    Veg.geo values: %s\n", paste(unique(matched_rows$Veg.geo), collapse=",")))
      }
    }
  }
  cat("=========================\n\n")
  
  # Prefer existing Veg values in df; fall back to GeoJSON mapping
  if ("Veg.geo" %in% names(joined)) {
    joined$Veg <- ifelse(is.na(joined$Veg) | joined$Veg == "", joined$Veg.geo, joined$Veg)
    joined$Veg.geo <- NULL
  }
  # Merge 'no soil' from gpts_map into df when df has missing values
  if ("no soil.geo" %in% names(joined)) {
    # If df already has a 'no soil' column, prefer it when non-NA; otherwise use geo value
    if ("no soil" %in% names(joined)) {
      joined$`no soil` <- ifelse(is.na(joined$`no soil`), joined$`no soil.geo`, joined$`no soil`)
    } else {
      joined$`no soil` <- joined$`no soil.geo`
    }
    joined$`no soil.geo` <- NULL
  }

  post_non_na <- sum(!is.na(joined$Veg) & joined$Veg != "")

  # If no Veg found, try matching by geo row-number values (CSV may contain numeric row ids)
  if (post_non_na == pre_non_na && "location_row" %in% names(gpts_map)) {
    # check for intersection between CSV values and location_row
    df_ids <- unique(na.omit(as.character(df$location_id)))
    match_count <- length(intersect(df_ids, unique(na.omit(as.character(gpts_map$location_row)))))
    if (match_count > 0) {
      cat(sprintf("[NOTICE] No matches by 'location_id' — attempting join by geojson row-number mapping (matched ids=%d)\n", match_count))
      joined2 <- dplyr::left_join(df, gpts_map, by = c("location_id" = "location_row"), suffix = c("", ".geo"))
      if ("Veg.geo" %in% names(joined2)) {
        joined2$Veg <- ifelse(is.na(joined2$Veg) | joined2$Veg == "", joined2$Veg.geo, joined2$Veg)
        joined2$Veg.geo <- NULL
      }
      if ("no soil.geo" %in% names(joined2)) {
        if ("no soil" %in% names(joined2)) {
          joined2$`no soil` <- ifelse(is.na(joined2$`no soil`), joined2$`no soil.geo`, joined2$`no soil`)
        } else {
          joined2$`no soil` <- joined2$`no soil.geo`
        }
        joined2$`no soil.geo` <- NULL
      }
      if (sum(!is.na(joined2$Veg) & joined2$Veg != "") > post_non_na) {
        joined <- normalize_no_soil_col(joined2)
        post_non_na <- sum(!is.na(joined$Veg) & joined$Veg != "")
        cat(sprintf("[NOTICE] GeoJSON join (row-number) gained %d Veg rows\n", post_non_na - pre_non_na))
      } else {
        cat("[NOTICE] Row-number join did not increase Veg mapping; keeping original join state.\n")
      }
    }
  }

  df <- normalize_no_soil_col(joined)

  matched_locs <- length(intersect(na.omit(unique(as.character(df$location_id))), na.omit(unique(as.character(gpts_map$location_id)))))
  cat(sprintf("[NOTICE] GeoJSON join results - Veg before=%d after=%d; matched location_id strings=%d\n", pre_non_na, post_non_na, matched_locs))

  if (post_non_na == 0L) {
    sample_df_ids <- unique(na.omit(as.character(head(df$location_id, 20))))
    sample_geo_ids <- unique(na.omit(as.character(head(gpts_map$location_id, 20))))
    sample_geo_rows <- unique(na.omit(as.character(head(gpts_map$location_row, 20))))
    cat("[WARNING] GeoJSON join produced no Veg values. Sample df$location_id (first 20):\n")
    print(sample_df_ids)
    cat("Sample gpts_map$location_id (first 20):\n")
    print(sample_geo_ids)
    cat("Sample gpts_map row-numbers (first 20):\n")
    print(sample_geo_rows)
    cat("Hint: CSV 'location_id' might be numeric row indices; try setting location IDs in the CSV to match the GeoJSON row order or use the row-number mapping.\n")
  }
}

# Add timing for major operations
timing_info <- list()
timing_info$start_time <- Sys.time()

cat("Starting vegetation mixture analysis with MESMA approach...\n")
cat("✓ Multiple endmember spectral mixture analysis\n")
cat("✓ Dynamic programming for year-constant proportions\n")
cat("✓ Variant switching allowed throughout the year\n")
cat("\n")
cat(sprintf("Dataset size: %d rows, %d columns\n", nrow(df), ncol(df)))
cat(sprintf("Number of locations: %d\n", length(unique(df$location_id))))
cat(sprintf("Date range: %s to %s\n", min(df$date, na.rm = TRUE), max(df$date, na.rm = TRUE)))

# Train/Test Split Implementation
cat("\n=== TRAIN/INFERENCE DATA CONFIGURATION ===\n")
cat(sprintf("Training years (config): %s\n", paste(TRAIN_YEARS, collapse = ", ")))
cat("Note: the script will always use all available years for data; TEST_YEARS is not supported.\n")

df$year <- lubridate::year(df$date)

## Use all available data for training by default. Keep a df_train object to
## preserve a clear training/inference distinction for downstream code.
df_train <- df
cat(sprintf(
  "Training dataset: %d rows from %d locations\n",
  nrow(df_train), length(unique(df_train$location_id))
))

## Inference/testing uses the full dataset (all years).
df_test <- df
cat(sprintf(
  "Testing dataset: %d rows from %d locations\n",
  nrow(df_test), length(unique(df_test$location_id))
))

df_full <- df
df_full <- normalize_no_soil_col(df_full)
df <- df_train

# Ensure training subset has at least some samples for each allowed vegetation class
# If a class from ALLOWED_VEG has zero samples in the available training set,
# add a small number of samples for that class from the full dataset.
if (exists("df_full") && "Veg" %in% names(df_full) && length(ALLOWED_VEG) > 0) {
  missing_vegs <- sapply(ALLOWED_VEG, function(v) {
    sum(tolower(df$Veg) == v, na.rm = TRUE)
  })
  missing_names <- names(missing_vegs)[missing_vegs == 0]
  if (length(missing_names) > 0) {
    for (mv in missing_names) {
      # find candidates anywhere in the full dataset
      cand <- df_full[tolower(df_full$Veg) == mv, , drop = FALSE]
      if (nrow(cand) > 0) {
        # Add a modest number of examples so the training set has representation
        add_n <- min(nrow(cand), max(5L, as.integer(floor(nrow(df) / 10))))
        add_n <- max(1L, add_n)
        to_add <- cand[seq_len(add_n), , drop = FALSE]
        df <- rbind(df, to_add)
        cat(sprintf("[NOTICE] Added %d samples from non-training years for Veg='%s' to ensure representation in training set\n", nrow(to_add), mv))
      } else {
        cat(sprintf("[WARNING] No samples found anywhere for Veg='%s'; cannot add examples to training set\n", mv))
      }
    }
    # re-report training dataset size after augmentation
    cat(sprintf("Training dataset after augmentation: %d rows from %d locations\n", nrow(df), length(unique(df$location_id))))
    # keep df_train consistent with df (augmented)
    df_train <- df
  }
}

cat("Using training data for vegetation library construction\n")
cat("=====================================\n\n")

# Initialize persistent parallel backend if enabled
cleanup_parallel <- setup_parallel_backend()
on.exit(cleanup_parallel(), add = TRUE)

# Progress logging helpers

log_msg <- function(...) {
  ts <- format(Sys.time(), "%H:%M:%S")
  msg <- sprintf("[%s] %s\n", ts, sprintf(...))
  cat(msg)
  if (isTRUE(PROGRESS_LOG_TO_FILE)) try(cat(msg, file = LOG_FILE, append = TRUE), silent = TRUE)
}



# Memory helpers: chunked rbind to avoid building huge intermediate objects
chunked_rbind <- function(lst, chunk_size = 50L) {
  if (length(lst) == 0) return(matrix(nrow = 0, ncol = 0))
  if (length(lst) == 1) return(lst[[1]])
  n <- length(lst)
  # Iteratively bind in chunks
  i <- 1L
  acc <- NULL
  while (i <= n) {
    end <- min(n, i + chunk_size - 1L)
    block <- do.call(rbind, lst[i:end])
    if (is.null(acc)) acc <- block else acc <- rbind(acc, block)
    # free memory for block
    rm(block); gc()
    i <- end + 1L
  }
  acc
}

# Debug helper: produce a nil return object for the main processing loop
# Return NULL so failed tasks do not produce placeholder rows that break downstream analysis
dbg_return_null <- function(reason = NULL) {
  # if (!is.null(reason)) cat(sprintf("[DEBUG] abort: %s\n", as.character(reason)))
  invisible(NULL)
}

# Clamp for parallel workers memory usage: limit workers according to available RAM
recommended_workers <- function(max_workers_hint = NULL) {
  # crude heuristic: check memory.limit on Windows or try to estimate via gc
  total_mem_gb <- tryCatch({
    if (.Platform$OS.type == "windows") {
      mem <- memory.limit() / 1024
      mem
    } else {
      NA_real_
    }
  }, error = function(e) stop(sprintf("recommended_workers: memory inspection failed: %s", e$message)))
  # conservative default
  if (!is.na(total_mem_gb) && total_mem_gb > 0) {
    # one worker per ~2 GB
    max(1L, floor(total_mem_gb / 2))
  } else {
    if (!is.null(max_workers_hint)) max(1L, as.integer(max_workers_hint)) else PARALLEL_WORKERS
  }
}

compute_inter_class_similarity_table <- function(mesma_lib, compressed_templates, grid_type = "full") {
  if (length(mesma_lib) < 2) return(NULL)
  
  get_vec <- function(veg, vid) {
    candidate <- NULL
    if (is.list(compressed_templates)) {
      candidate <- tryCatch(compressed_templates[[veg]][[vid]][[grid_type]], error = function(e) NULL)
    } else if (is.environment(compressed_templates)) {
      candidate <- tryCatch(get(vid, envir = compressed_templates[[veg]][[grid_type]]), error = function(e) NULL)
    } else if (is.function(compressed_templates)) {
      candidate <- tryCatch(compressed_templates(veg, vid, grid_type = grid_type), error = function(e) NULL)
    }
    if (!is.null(candidate) && length(candidate) > 0) {
      return(as.numeric(candidate))
    }
    NULL
  }

  veg_types <- names(mesma_lib)
  combo_types <- utils::combn(veg_types, 2, simplify = FALSE)
  rows <- list()

  for (cmb in combo_types) {
    veg1 <- cmb[1]
    veg2 <- cmb[2]
    vars1 <- mesma_lib[[veg1]]
    vars2 <- mesma_lib[[veg2]]
    
    for (v1 in vars1) {
      vec1 <- get_vec(veg1, v1$variant_id)
      if (is.null(vec1)) next
      
      for (v2 in vars2) {
        vec2 <- get_vec(veg2, v2$variant_id)
        if (is.null(vec2)) next
        
        if (length(vec1) != length(vec2)) {
          common_len <- min(length(vec1), length(vec2))
          vec1 <- vec1[seq_len(common_len)]
          vec2 <- vec2[seq_len(common_len)]
        }
        
        cos_val <- cos_sim(vec1, vec2)
        dist_val <- sqrt(sum((vec1 - vec2)^2))
        
        rows[[length(rows) + 1]] <- data.frame(
          veg1 = veg1,
          variant1 = v1$variant_id,
          veg2 = veg2,
          variant2 = v2$variant_id,
          cos_sim = cos_val,
          euclidean_dist = dist_val,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(rows) > 0) do.call(rbind, rows) else NULL
}





# Whitening helper: moved to top-level so other stages can reuse it
whiten_matrix <- function(X, epsilon = 1e-6) {
  if (is.null(X) || nrow(X) < 2 || ncol(X) < 1) return(list(Xw = X, W = diag(ncol(X)), mu = rep(0, ncol(X))))
  mu <- colMeans(X)
  Xc <- sweep(X, 2, mu, "-")
  Sigma <- cov(Xc)
  eig <- eigen(Sigma)
  vals <- eig$values
  vals[vals < epsilon] <- epsilon
  D_inv_sqrt <- diag(1 / sqrt(vals))
  W <- eig$vectors %*% D_inv_sqrt %*% t(eig$vectors)
  Xw <- Xc %*% W
  list(Xw = Xw, W = W, mu = mu)
}

#' Two-stage dimensionality reduction: PCA for decorrelation, LDA for discrimination
#' 
#' Stage 1: PCA removes multicollinearity while preserving most variance
#' Stage 2: LDA finds directions that best separate vegetation classes
#' Compute PCA-LDA weights (weekly resolution for improved dimensionality reduction)
compute_pca_lda_weights <- function(lib_df, avail_idx, pca_variance_threshold = PCA_VARIANCE_THRESHOLD, lda_weight_floor = LDA_WEIGHT_FLOOR) {
  
  cat("\n=== COMPUTING PCA → LDA FEATURE WEIGHTS ===\n")
  
  # Step 1: Build compressed feature matrix from all training traces
  veg_types <- unique(na.omit(lib_df$Veg))
  veg_types <- veg_types[veg_types != "" & tolower(veg_types) != "barren"]
  
  if (length(veg_types) < 2) {
    warning("Need at least 2 vegetation classes for LDA")
    return(NULL)
  }
  
  # Collect all traces
  all_features <- list()
  all_labels <- c()
  
  for (veg in veg_types) {
    veg_data <- lib_df[tolower(lib_df$Veg) == tolower(veg), ]
    traces <- unique(veg_data[, c("location_id", "year")])
    
    for (i in seq_len(nrow(traces))) {
      loc <- traces$location_id[i]
      yr <- traces$year[i]
      
      dly_year <- veg_data[veg_data$location_id == loc & veg_data$year == yr, ]
      n_unique_doys <- length(unique(dly_year$doy))
      if (n_unique_doys < 5) next
      
      # Build weekly 53 × K matrix (reduced temporal resolution)
      raw_mat <- build_weekly_matrix(dly_year, avail_idx)
      if (is.null(raw_mat)) next
      
      # Require at least 5 unique DOYs
      n_unique_doys <- length(unique(dly_year$doy))
      if (n_unique_doys < 5) next
      
      # Use weekly matrix flattened as features (reduced dimensionality)
      compressed <- as.numeric(raw_mat)  # 53 * K vector
      if (any(!is.finite(compressed))) next
      
      all_features[[length(all_features) + 1]] <- compressed
      all_labels <- c(all_labels, veg)
    }
  }
  
  if (length(all_features) < 20) {
    warning("Insufficient traces for PCA-LDA")
    return(NULL)
  }
  
  X <- do.call(rbind, all_features)
  y <- factor(all_labels)
  
  cat(sprintf("Feature matrix: %d samples × %d features (weekly resolution)\n", nrow(X), ncol(X)))
  cat(sprintf("Class distribution: %s\n", 
              paste(sprintf("%s=%d", levels(y), table(y)), collapse=", ")))
  
  # Step 2: PCA for multicollinearity removal
  # Center but don't scale (preserve relative magnitudes across indices)
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  X_centered[!is.finite(X_centered)] <- 0
  
  # Remove zero-variance columns
  col_vars <- apply(X_centered, 2, var, na.rm = TRUE)
  keep_cols <- col_vars > 1e-10
  if (sum(keep_cols) < 10) {
    warning("Too few variable features after filtering")
    return(NULL)
  }
  X_filtered <- X_centered[, keep_cols, drop = FALSE]
  
  cat(sprintf("After variance filtering: %d features\n", ncol(X_filtered)))
  
  # PCA
  pca_result <- prcomp(X_filtered, center = FALSE, scale. = FALSE)
  
  # Determine number of PCs to retain
  cum_var <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
  n_pcs <- which(cum_var >= pca_variance_threshold)[1]
  if (is.na(n_pcs)) n_pcs <- length(cum_var)
  n_pcs <- max(n_pcs, length(veg_types))  # At least as many PCs as classes
  n_pcs <- min(n_pcs, ncol(X_filtered), nrow(X_filtered) - 1)
  
  cat(sprintf("PCA: Retaining %d PCs (%.1f%% variance)\n", n_pcs, 100 * cum_var[n_pcs]))
  
  X_pca <- pca_result$x[, 1:n_pcs, drop = FALSE]
  
  # Step 3: LDA on PCA scores
  # Check class sizes
  class_counts <- table(y)
  if (any(class_counts < 3)) {
    warning("Some classes have fewer than 3 samples")
    # Remove small classes
    keep_classes <- names(class_counts)[class_counts >= 3]
    keep_rows <- y %in% keep_classes
    X_pca <- X_pca[keep_rows, , drop = FALSE]
    y <- factor(y[keep_rows])
  }
  
  if (length(levels(y)) < 2) {
    warning("Fewer than 2 classes after filtering")
    return(NULL)
  }
  
  lda_result <- tryCatch({
    MASS::lda(X_pca, grouping = y)
  }, error = function(e) {
    warning(sprintf("LDA failed: %s", e$message))
    NULL
  })
  
  if (is.null(lda_result)) return(NULL)
  
  # Step 4: Project LDA weights back to original feature space
  # LDA scaling: n_pcs × (n_classes - 1)
  # PCA rotation: n_original × n_pcs
  
  W_lda <- lda_result$scaling  # n_pcs × n_discriminants
  R_pca <- pca_result$rotation[, 1:n_pcs, drop = FALSE]  # n_filtered × n_pcs
  
  # Combined projection: original → PCA → LDA
  # Weight in original space = R_pca %*% W_lda
  W_combined <- R_pca %*% W_lda  # n_filtered × n_discriminants
  
  # Weight importance: sum of absolute contributions across discriminant axes
  # Weighted by proportion of between-class variance explained
  svd_vals <- lda_result$svd
  prop_var <- svd_vals^2 / sum(svd_vals^2)
  
  if (ncol(W_combined) > 1) {
    # Weighted sum across discriminant functions
    feature_importance <- rowSums(abs(W_combined) %*% diag(prop_var))
  } else {
    feature_importance <- abs(W_combined[, 1])
  }
  
  # Map back to full feature set (including zero-variance columns)
  full_weights <- rep(0, ncol(X_centered))
  full_weights[keep_cols] <- feature_importance
  
  # Normalize to [0, 1]
  if (max(full_weights) > 0) {
    full_weights <- full_weights / max(full_weights)
  }
  
  # Apply floor (don't completely zero out features)
  full_weights <- pmax(full_weights, lda_weight_floor * (full_weights > 0))
  
  # Re-normalize
  if (max(full_weights) > 0) {
    full_weights <- full_weights / max(full_weights)
  }
  
  cat(sprintf("LDA weights computed: %d non-zero features\n", sum(full_weights > lda_weight_floor)))
  
  # Return weights along with PCA/LDA objects for potential reuse
  list(
    weights = full_weights,
    pca = pca_result,
    lda = lda_result,
    n_pcs = n_pcs,
    keep_cols = keep_cols,
    feature_names = colnames(X)
  )
}

#' Build weekly aggregated K × 52 index matrix from daily data
#' Reduces temporal resolution from daily to weekly for better LDA performance
build_weekly_matrix <- function(dly_year, avail_idx) {
  if (nrow(dly_year) == 0) return(NULL)

  if (!"doy" %in% names(dly_year)) {
    dly_year$doy <- lubridate::yday(dly_year$date)
  }

  K <- length(avail_idx)
  # 52 weeks + 1 extra day bin = 53 bins total
  raw_mat <- matrix(NA_real_, nrow = 53, ncol = K)
  colnames(raw_mat) <- avail_idx

  for (j in seq_along(avail_idx)) {
    idx <- avail_idx[j]
    if (!idx %in% names(dly_year)) next

    # Aggregate by DOY (median for robustness)
    vals_by_doy <- tapply(dly_year[[idx]], dly_year$doy, function(v) {
      v <- v[is.finite(v)]
      if (length(v) == 0) NA_real_ else median(v)
    })

    doy_values <- as.integer(names(vals_by_doy))
    valid <- doy_values >= 1 & doy_values <= 365
    daily_vals <- rep(NA_real_, 365)
    daily_vals[doy_values[valid]] <- vals_by_doy[valid]

    # Interpolate missing DOYs first
    finite_idx <- which(is.finite(daily_vals))
    if (length(finite_idx) >= 2) {
      missing_idx <- which(!is.finite(daily_vals))
      if (length(missing_idx) > 0) {
        daily_vals[missing_idx] <- approx(
          x = finite_idx,
          y = daily_vals[finite_idx],
          xout = missing_idx,
          rule = 2
        )$y
      }
    } else if (length(finite_idx) == 1) {
      daily_vals[] <- daily_vals[finite_idx]
    } else {
      daily_vals[] <- 0
    }

    # Aggregate into weekly bins (7-day windows)
    for (week in 1:52) {
      start_day <- (week - 1) * 7 + 1
      end_day <- min(week * 7, 365)
      week_vals <- daily_vals[start_day:end_day]
      week_vals <- week_vals[is.finite(week_vals)]
      if (length(week_vals) > 0) {
        raw_mat[week, j] <- mean(week_vals)
      }
    }

    # Handle remaining days (366th day if leap year, or partial last week)
    remaining_vals <- daily_vals[366:min(365, length(daily_vals))]
    remaining_vals <- remaining_vals[is.finite(remaining_vals)]
    if (length(remaining_vals) > 0) {
      raw_mat[53, j] <- mean(remaining_vals)
    } else if (is.finite(raw_mat[52, j])) {
      # Use last complete week if no remaining data
      raw_mat[53, j] <- raw_mat[52, j]
    } else {
      raw_mat[53, j] <- 0
    }
  }

  raw_mat
}

#' Build raw 365 × K index matrix from daily data (NO whitening)
build_raw_365_matrix <- function(dly_year, avail_idx) {
  if (nrow(dly_year) == 0) return(NULL)

  if (!"doy" %in% names(dly_year)) {
    dly_year$doy <- lubridate::yday(dly_year$date)
  }

  K <- length(avail_idx)
  raw_mat <- matrix(NA_real_, nrow = 365, ncol = K)
  colnames(raw_mat) <- avail_idx

  for (j in seq_along(avail_idx)) {
    idx <- avail_idx[j]
    if (!idx %in% names(dly_year)) next

    # Aggregate by DOY (median for robustness)
    vals_by_doy <- tapply(dly_year[[idx]], dly_year$doy, function(v) {
      v <- v[is.finite(v)]
      if (length(v) == 0) NA_real_ else median(v)
    })

    doy_values <- as.integer(names(vals_by_doy))
    valid <- doy_values >= 1 & doy_values <= 365
    raw_mat[doy_values[valid], j] <- vals_by_doy[valid]

    # Interpolate missing DOYs
    finite_idx <- which(is.finite(raw_mat[, j]))
    if (length(finite_idx) >= 2) {
      missing_idx <- which(!is.finite(raw_mat[, j]))
      if (length(missing_idx) > 0) {
        raw_mat[missing_idx, j] <- approx(
          x = finite_idx,
          y = raw_mat[finite_idx, j],
          xout = missing_idx,
          rule = 2
        )$y
      }
    } else if (length(finite_idx) == 1) {
      raw_mat[, j] <- raw_mat[finite_idx, j]
    } else {
      raw_mat[, j] <- 0
    }
  }

  raw_mat
}

#' Simple grid compression (NO whitening, NO PCA on individual trace)
#' Just samples the raw index values at fixed time points
#' Apply PCA-LDA weights during unmixing
#' This replaces the whitening step in the original code
apply_pca_lda_transform <- function(y, pca_lda_result) {
  if (is.null(pca_lda_result)) return(y)
  
  weights <- pca_lda_result$weights
  
  if (length(weights) != length(y)) {
    warning("Weight length mismatch, using unweighted")
    return(y)
  }
  
  # Apply weights (element-wise multiplication)
  y_weighted <- y * weights
  
  y_weighted
}

#' Compute weighted cosine similarity using PCA-LDA weights
weighted_cosine_similarity <- function(a, b, weights = NULL) {
  if (is.null(weights) || length(weights) != length(a)) {
    # Unweighted
    return(sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2))))
  }
  
  # Apply weights
  a_w <- a * sqrt(weights)
  b_w <- b * sqrt(weights)
  
  sum(a_w * b_w) / (sqrt(sum(a_w^2)) * sqrt(sum(b_w^2)))
}

# Helper: index of medoid row (row with minimal total distance)
# For large matrices, sample to avoid O(N^2) memory allocation
medoid_row_index <- function(M) {
  if (is.null(M) || !is.matrix(M) || nrow(M) == 0) return(NA_integer_)
  if (nrow(M) == 1) return(1L)
  X <- M
  X[!is.finite(X)] <- 0
  
  # For large matrices, compute medoid on a random subsample
  n <- nrow(X)
  
  if (n > MEDOID_SAMPLE_SIZE) {
    # Sample rows, find medoid among sample, return index in original matrix
    samp_idx <- sample.int(n, size = min(as.integer(MEDOID_SAMPLE_SIZE), n))
    X_samp <- X[samp_idx, , drop = FALSE]
    D_samp <- tryCatch(as.matrix(stats::dist(X_samp)), error = function(e) stop(sprintf("medoid_row_index: distance computation failed: %s", e$message)))
    rs_samp <- rowSums(D_samp, na.rm = TRUE)
    local_idx <- which.min(rs_samp)
    if (length(local_idx) == 0 || !is.finite(local_idx)) local_idx <- 1L
    idx <- samp_idx[local_idx]
  } else {
    D <- tryCatch(as.matrix(stats::dist(X)), error = function(e) stop(sprintf("medoid_row_index: distance computation failed: %s", e$message)))
    rs <- rowSums(D, na.rm = TRUE)
    idx <- which.min(rs)
    if (length(idx) == 0 || !is.finite(idx)) idx <- 1L
  }
  as.integer(idx)
}

# Cosine similarity helper (promoted to top-level)
cos_sim <- function(a, b) {
  da <- sqrt(safe_dot(a, a)); db <- sqrt(safe_dot(b, b))
  if (da == 0 || db == 0) return(0)
  safe_dot(a, b) / (da * db)
}

# Helper: parallel map
.run_map <- function(X, FUN) {
  f_FUN <- FUN
  
  # Check for progress bar support
  use_pbapply <- requireNamespace("pbapply", quietly = TRUE)

  if (!PARALLEL_ENABLE) {
    if (use_pbapply) {
      pbapply::pblapply(X, function(x) { f_FUN(x) })
    } else {
      lapply(X, function(x) { f_FUN(x) })
    }
  } else {
    if (!requireNamespace("future.apply", quietly = TRUE) || !requireNamespace("future", quietly = TRUE)) {
      stop("future.apply and future packages required for parallel processing")
    }

    if (!isTRUE(PERSISTENT_PARALLEL_BACKEND)) {
      old_plan <- future::plan()
      options(future.globals.maxSize = max(getOption("future.globals.maxSize", 2e9), 12e9))
      future::plan(future::multisession, workers = PARALLEL_WORKERS)
      on.exit(future::plan(old_plan))
    }

    if (use_pbapply) {
      # pbapply supports future backend if plan is set
      pbapply::pblapply(X, function(x) { f_FUN(x) }, cl = "future", future.seed = TRUE)
    } else {
      cat("Processing (no progress bar available)...\n")
      future.apply::future_lapply(X, function(x) {
        f_FUN(x)
      }, future.seed = TRUE)
    }
  }
}

# Read GeoJSON for vegetation ground truth


loc_years <- data.frame(location_id = character(0), year = integer(0), stringsAsFactors = FALSE)

if (!"year" %in% names(df)) {
  if ("date" %in% names(df)) {
    if (!requireNamespace("lubridate", quietly = TRUE)) stop("The package 'lubridate' is required")
    df$year <- as.integer(lubridate::year(as.Date(df$date)))
  }
}

df <- df %>% filter(year >= 1985 & year <= 2025)

if (!"Veg" %in% names(df)) df$Veg <- NA_character_

lon_candidates <- names(df)[grepl("(^|_)(lon|longitude|x)(_|$)", names(df), ignore.case = TRUE)]
lat_candidates <- names(df)[grepl("(^|_)(lat|latitude|y)(_|$)", names(df), ignore.case = TRUE)]
if (length(lon_candidates) > 0 && length(lat_candidates) > 0) {
  if ("location_id" %in% names(df) && "location_id" %in% names(gpts_map) && nrow(gpts_map) > 0) {
    # coerce consistent type before join
    if (!is.character(df$location_id)) {
      df$location_id <- as.character(df$location_id)
      cat("[NOTICE] Coerced df$location_id to character to match gpts_map for joining (second join).\n")
    }
    if (!is.character(gpts_map$location_id)) {
      gpts_map$location_id <- as.character(gpts_map$location_id)
      cat("[NOTICE] Coerced gpts_map$location_id to character for joining (second join).\n")
    }
    df <- dplyr::left_join(df, gpts_map, by = "location_id", suffix = c("", ".y"))
    if ("Veg.y" %in% names(df)) {
      if (any(!is.na(df$Veg.y))) {
        df$Veg <- ifelse(is.na(df$Veg) | df$Veg == "", df$Veg.y, df$Veg)
      }
      df$Veg.y <- NULL
    }
  }
}

df$Veg <- tolower(df$Veg)

# Join df_train with gpts_map to add Veg column
df_train$location_id <- as.character(df_train$location_id)
df_train <- dplyr::left_join(df_train, gpts_map, by = "location_id", suffix = c("", ".geo"))
# Prefer existing Veg values in df_train; fall back to GeoJSON mapping
if ("Veg.geo" %in% names(df_train)) {
  df_train$Veg <- ifelse(is.na(df_train$Veg) | df_train$Veg == "", df_train$Veg.geo, df_train$Veg)
  df_train$Veg.geo <- NULL
}
# Merge 'no soil' from gpts_map into df_train when df_train has missing values
if ("no soil.geo" %in% names(df_train)) {
  # If df_train already has a 'no soil' column, prefer it when non-NA; otherwise use geo value
  if ("no soil" %in% names(df_train)) {
    df_train$`no soil` <- ifelse(is.na(df_train$`no soil`), df_train$`no soil.geo`, df_train$`no soil`)
  } else {
    df_train$`no soil` <- df_train$`no soil.geo`
  }
  df_train$`no soil.geo` <- NULL
}
df_train <- normalize_no_soil_col(df_train)
df_train$Veg <- tolower(df_train$Veg)

if (!"date" %in% names(df)) stop("Input CSV must contain a 'date' column")
df$date <- as.Date(df$date)
if (!"location_id" %in% names(df)) stop("Input CSV must contain a 'location_id' column")

if ("Veg" %in% names(df)) df$Veg <- tolower(as.character(df$Veg))

required_pkgs <- c("future", "future.apply")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    try(
      {
        install.packages(pkg, repos = "https://cloud.r-project.org")
      },
      silent = TRUE
    )
  }
}

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx", repos = "https://cloud.r-project.org")
}

df$doy <- lubridate::yday(df$date)
df$doy[df$doy < 1 | df$doy > 366] <- NA_integer_

# DEBUG: Check barren loading status
cat("\n=== BARREN LOADING DEBUG ===\n")
if (exists("gpts_map") && nrow(gpts_map) > 0) {
  # Check what barren locations exist in gpts_map
  barren_in_gpts <- gpts_map[tolower(gpts_map$Veg) == "barren" | (!is.na(gpts_map$`no soil`) & gpts_map$`no soil` == 0), ]
  cat(sprintf("Barren locations in gpts_map: %d\n", nrow(barren_in_gpts)))
  if (nrow(barren_in_gpts) > 0) {
    cat(sprintf("  Sample location_ids: %s\n", paste(head(unique(barren_in_gpts$location_id), 5), collapse=", ")))
  }
  
  # Check if these locations exist in df with actual phenology data
  barren_locs <- unique(barren_in_gpts$location_id)
  barren_in_df <- df[df$location_id %in% barren_locs, ]
  cat(sprintf("Barren locations found in phenology df: %d rows from %d locations\n", 
              nrow(barren_in_df), length(unique(barren_in_df$location_id))))
  
  # Check if Veg column is properly set for barren
  barren_veg_df <- df[tolower(df$Veg) == "barren", ]
  cat(sprintf("Rows with Veg='barren' in df: %d\n", nrow(barren_veg_df)))
  
  # If barren locations exist in df but Veg is not set, we need to fix the join
  if (nrow(barren_in_df) > 0 && nrow(barren_veg_df) == 0) {
    cat("[WARNING] Barren locations exist in phenology data but Veg column not set to 'barren'!\n")
    cat("  This indicates a join issue. Checking Veg values for barren locations:\n")
    cat(sprintf("  Unique Veg values: %s\n", paste(unique(barren_in_df$Veg), collapse=", ")))
  }
}
cat("=============================\n\n")

# NOTE: We do NOT add synthetic barren rows - barren data must come from actual phenology measurements
# If barren locations are missing from the phenology CSV, they simply won't be used for training

veg_counts <- sort(table(na.omit(df$Veg)), decreasing = TRUE)
cat("Vegetation class counts after loading:\n")
print(veg_counts)

# Define indices to use
meta_cols <- intersect(c(
  "date", "location_id", "Veg", "coverage", "lat", "lon", "latitude", "longitude",
  "target_lon", "target_lat", "imagery_lat", "imagery_lon", "doy", "year"
), names(df))
numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
candidate_indices <- intersect(OPTIMAL_INDICES, numeric_cols)
found_opt <- intersect(OPTIMAL_INDICES, numeric_cols)
missing_opt <- setdiff(OPTIMAL_INDICES, numeric_cols)
found_raw <- intersect(RAW_BANDS, numeric_cols)
missing_raw <- setdiff(RAW_BANDS, numeric_cols)

if (length(found_opt) > 0) {
  cat(sprintf("Found %d OPTIMAL_INDICES in input: %s\n", length(found_opt), paste(found_opt, collapse = ", ")))
}
if (length(missing_opt) > 0) {
  cat(sprintf("Missing OPTIMAL_INDICES in input: %s\n", paste(missing_opt, collapse = ", ")))
}
if (length(found_raw) > 0) {
  cat(sprintf("Found %d RAW_BANDS in input: %s\n", length(found_raw), paste(found_raw, collapse = ", ")))
}
if (length(missing_raw) > 0) {
  cat(sprintf("Missing RAW_BANDS in input: %s\n", paste(missing_raw, collapse = ", ")))
}

candidate_indices <- c(found_opt, found_raw)
if (length(candidate_indices) == 0) {
  stop("No OPTIMAL_INDICES present in input CSV")
}

# If some optimal indices are missing but raw bands exist, compute those indices from raw bands (fusion path)
if (length(missing_opt) > 0 && length(found_raw) >= 2) {
  before_cols <- names(df)
  df <- compute_indices_from_bands(df)
  numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  found_opt <- intersect(OPTIMAL_INDICES, numeric_cols)
  missing_opt <- setdiff(OPTIMAL_INDICES, numeric_cols)
  candidate_indices <- unique(c(found_opt, intersect(RAW_BANDS, numeric_cols)))
  new_cols <- setdiff(names(df), before_cols)
  if (length(new_cols) > 0) cat(sprintf("[NOTICE] Computed indices from raw bands in training data: %s\n", paste(new_cols, collapse = ", ")))
}

# Candidate indices available for processing
avail <- candidate_indices

# Correlation-based filtering REMOVED as per user request
# (Code block removed)


if (length(avail) == 0) {
  stop("No indices remain after correlation filtering; check input candidate indices and correlation threshold")
}

if (length(avail) < USE_INDICES_MIN) {
  stop(sprintf(
    "Only %d indices remain after filtering, minimum required is %d",
    length(avail), USE_INDICES_MIN
  ))
}

cat(sprintf("Selected %d indices: %s\n", length(avail), paste(avail, collapse = ", ")))

  # Data sufficiency check DISABLED per user request - all locations are kept regardless of DOY coverage
  
  if (FALSE && (is.null(TRAIN_YEARS) || length(TRAIN_YEARS) == 0)) stop("TRAIN_YEARS must be defined and non-empty for data sufficiency checks; no fallback permitted.")
  if (FALSE && "location_id" %in% names(df_full)) {
    # For each location and each selected index, count finite observations and error if below threshold
    offenders <- list()
    for (loc in unique(df_full$location_id)) {
      sub <- df_full[df_full$location_id == loc, , drop = FALSE]
      # Determine number of training years present for this location
      loc_years <- unique(sub$year[is.finite(sub$year) & !is.na(sub$year)])
      train_loc_years <- intersect(as.integer(loc_years), as.integer(TRAIN_YEARS))
      n_years <- length(train_loc_years)
      if (n_years == 0) {
        stop(sprintf("Data sufficiency check failed: location '%s' has no observations in TRAIN_YEARS (%s) — cannot compute required DOY threshold; aborting (no fallback permitted).", loc, paste(TRAIN_YEARS, collapse = ", ")))
      }

      # required days = 50% of all possible DOYs across the location's training years
      # Allow skipping a minimal number of DOYs per-location (controlled by MIN_SKIP_DOYS_PER_LOCATION)
      base_required <- as.integer(ceiling(n_years * 365 / 2))
      required_days_for_loc <- as.integer(max(1L, base_required - as.integer(MIN_SKIP_DOYS_PER_LOCATION)))

      # Count unique DOYs in the training years where at least one of the selected indices is finite
      if (nrow(sub) == 0) {
        n_doys_present <- 0L
      } else {
        sub_tr <- sub[sub$year %in% train_loc_years, , drop = FALSE]
        if (nrow(sub_tr) == 0) {
          n_doys_present <- 0L
        } else {
          # Mark rows where any of the selected indices are finite
          has_any_index <- apply(sub_tr[, intersect(avail, names(sub_tr)), drop = FALSE], 1, function(r) any(is.finite(r)))
          # Count unique dates (or doy) present with any index
          if ('date' %in% names(sub_tr)) {
            n_doys_present <- length(unique(as.character(sub_tr$date[which(has_any_index)])))
          } else if ('doy' %in% names(sub_tr)) {
            n_doys_present <- length(unique(sub_tr$doy[which(has_any_index)]))
          } else {
            # fallback to counting rows with any index finite
            n_doys_present <- sum(has_any_index)
          }
        }
      }

      if (n_doys_present < required_days_for_loc) {
        offenders[[length(offenders) + 1]] <- list(location = loc, n = n_doys_present, required = required_days_for_loc, years = paste(train_loc_years, collapse = ","))
      }
    }
    if (length(offenders) > 0) {
      # Build detailed messages per offending location (aggregate across indices)
      msg_lines <- vapply(offenders, function(x) sprintf("loc=%s n=%d required=%d train_years=%s", x$location, x$n, x$required, x$years), character(1))
      bad_locs <- unique(vapply(offenders, function(x) as.character(x$location), character(1)))

      warning(sprintf(
        "Data sufficiency: skipping %d locations because they do not meet the per-location DOY threshold (50%% of DOYs across training years) after allowing up to %d missing DOYs per-location.\nSkipped locations: %s\nDetails (loc,obs,required,train_years):\n%s",
        length(bad_locs), as.integer(MIN_SKIP_DOYS_PER_LOCATION), paste(head(bad_locs, 100), collapse = ", "), paste(msg_lines, collapse = "\n")
      ))

      
      # we flag them instead so downstream code can decide how to handle them.
      flag_col <- ".insufficient_data_by_DOY_threshold"
      if (!flag_col %in% names(df_full)) df_full[[flag_col]] <- FALSE
      df_full[[flag_col]][df_full$location_id %in% bad_locs] <- TRUE
      if (exists("df")) {
        if (!flag_col %in% names(df)) df[[flag_col]] <- FALSE
        df[[flag_col]][df$location_id %in% bad_locs] <- TRUE
      }
      if (exists("df_train")) {
        if (!flag_col %in% names(df_train)) df_train[[flag_col]] <- FALSE
        df_train[[flag_col]][df_train$location_id %in% bad_locs] <- TRUE
      }
      if (exists("df_test")) {
        if (!flag_col %in% names(df_test)) df_test[[flag_col]] <- FALSE
        df_test[[flag_col]][df_test$location_id %in% bad_locs] <- TRUE
      }

      # Keep all locations but notify if this leaves nothing usable
      if (nrow(df_full) == 0 || length(unique(na.omit(df_full$location_id))) == 0) {
        stop("Data sufficiency check left no usable data — aborting.")
      }
    }
  } else {
    # Global time series: require TRAIN_YEARS and enforce proportional check globally
    if (is.null(TRAIN_YEARS) || length(TRAIN_YEARS) == 0) stop("TRAIN_YEARS must be defined for global sufficiency checks; no fallback permitted.")
    base_req_global <- as.integer(ceiling(length(TRAIN_YEARS) * 365 / 2))
    required_days_global <- as.integer(max(1L, base_req_global - as.integer(MIN_SKIP_DOYS_PER_LOCATION)))
    if (nrow(df_full) < required_days_global) stop(sprintf("Data sufficiency check failed: overall data contains fewer than %d rows (50%% of days across TRAIN_YEARS after allowing up to %d missing DOYs)", required_days_global, MIN_SKIP_DOYS_PER_LOCATION))
  }

#
# SKIP MOVING VARIANCE CALCULATION IF FLAG IS SET
if (isTRUE(SKIP_MOVING_VARIANCE)) {
  cat("Skipping moving variance calculation (SKIP_MOVING_VARIANCE = TRUE)...\n")
  seasonal_var_cols <- c()
  timing_info$moving_var_done <- Sys.time()
} else if (!isTRUE(TESTING_MODE)) {
  if (length(avail) > 0) {
    pb <- utils::txtProgressBar(min = 0, max = length(avail), style = 3)
    i <- 0
    for (idx_name in avail) {
      var_col <- paste0(idx_name, "_var14")
      v_full <- calc_moving_var(df_full, idx_name, window = 14)
      df_full[[var_col]] <- v_full

      if (exists("df") && nrow(df) > 0) {
        key_full <- paste0(as.character(df_full$location_id), "__", as.character(df_full$date))
        key_df <- paste0(as.character(df$location_id), "__", as.character(df$date))
        mpos <- match(key_df, key_full)
        df[[var_col]] <- ifelse(is.na(mpos), NA_real_, df_full[[var_col]][mpos])
      }

      if (exists("df_test") && nrow(df_test) > 0) {
        key_test <- paste0(as.character(df_test$location_id), "__", as.character(df_test$date))
        key_full <- paste0(as.character(df_full$location_id), "__", as.character(df_full$date))
        mpos2 <- match(key_test, key_full)
        df_test[[var_col]] <- ifelse(is.na(mpos2), NA_real_, df_full[[var_col]][mpos2])
      }

      i <- i + 1
      try(utils::setTxtProgressBar(pb, i), silent = TRUE)
    }
    try(close(pb), silent = TRUE)
  }

  seasonal_var_cols <- c()
  if (length(avail) > 0) {
    peak_trough_data <- lapply(avail, function(idx_name) {
      var_col <- paste0(idx_name, "_var14")
      v <- df[[var_col]]
      if (all(is.na(v))) {
        return(NULL)
      }
      peaks <- which(diff(sign(diff(v))) == -2)
      troughs <- which(diff(sign(diff(v))) == 2)
      if (length(peaks) >= 1 && length(troughs) >= 1) var_col else NULL
    })
    seasonal_var_cols <- unlist(peak_trough_data[!sapply(peak_trough_data, is.null)])
  }

  df <- df[, c(names(df)[!grepl("_var14$", names(df))], seasonal_var_cols)]

  timing_info$moving_var_done <- Sys.time()
  cat(sprintf(
    "Moving variance calculation completed in %.1f seconds\n",
    as.numeric(difftime(timing_info$moving_var_done, timing_info$start_time, units = "secs"))
  ))
}

# Post-variance cleaning - SKIP IF MOVING VARIANCE WAS SKIPPED
if (!isTRUE(SKIP_MOVING_VARIANCE)) {
var_cols <- names(df)[grepl("(_var14$|_mv$)", names(df))]
if (length(var_cols) > 0) {
  cat(sprintf("Post-variance cleaning: %d variance columns detected\n", length(var_cols)))

  ac1 <- numeric(length(var_cols))
  names(ac1) <- var_cols
  cv <- numeric(length(var_cols))
  names(cv) <- var_cols
  finite_frac <- numeric(length(var_cols))
  names(finite_frac) <- var_cols

  for (i in seq_along(var_cols)) {
    v <- df[[var_cols[i]]]
    v_f <- v[is.finite(v)]
    finite_frac[i] <- if (length(v) == 0) 0 else sum(is.finite(v)) / length(v)
    if (length(v_f) > 3) {
      ac1[i] <- tryCatch(
        {
          cor(v_f[-1], v_f[-length(v_f)], use = "complete.obs")
        },
        error = function(e) stop(sprintf("Post-variance cleaning: failed to compute lag-1 autocorrelation for %s: %s", var_cols[i], e$message))
      )
      mm <- median(v_f, na.rm = TRUE)
      msd <- stats::sd(v_f, na.rm = TRUE)
      cv[i] <- if (is.finite(mm) && mm != 0) msd / abs(mm) else Inf
    } else {
      ac1[i] <- NA_real_
      cv[i] <- NA_real_
    }
  }

  keep_mask <- (finite_frac >= 0.30) & (is.finite(ac1) & ac1 > 0.10 | is.finite(cv) & cv > 0.05)
  remove_rand <- var_cols[!keep_mask]
  if (length(remove_rand) > 0) {
    cat(sprintf("Removing %d random-looking variance cols: %s\n", length(remove_rand), paste(remove_rand, collapse = ", ")))
    df[remove_rand] <- NULL
    var_cols <- setdiff(var_cols, remove_rand)
  }

# Correlation-based filtering of variance columns REMOVED as per user request
# (Code block removed)

  
  # Add surviving variance columns to avail
  if (length(var_cols) > 0) {
    cat(sprintf("Adding %d moving variance features to avail: %s\n", length(var_cols), paste(var_cols, collapse=", ")))
    avail <- unique(c(avail, var_cols))
  }
}

  # Final correlation-based filtering REMOVED as per user request
  # (Code block removed)


} # End of if (!isTRUE(SKIP_MOVING_VARIANCE))

cat(sprintf("Post-processing rows (baseline subtraction disabled): %d\n", nrow(df)))
cat("Data preprocessing complete.\n")
adj_cols <- intersect(avail, names(df))

# Simple feature pruning
if (FALSE && length(avail) > 0) {
  # Originally: compute per-index SD and drop indices below threshold
  idx_sd <- vapply(avail, function(nm) {
    x <- df[[nm]]
    x <- x[is.finite(x)]
    if (length(x) < 5) {
      return(0)
    }
    stats::sd(x, na.rm = TRUE)
  }, numeric(1))
  keep <- idx_sd >= MIN_INDEX_SD
  dropped <- avail[!keep]
  if (length(dropped) > 0) cat("Pruned low-variance indices:", paste(dropped, collapse = ", "), "\n")
  avail <- avail[keep]
}
cat("[NOTICE] SD-based simple feature pruning disabled; keeping all indices (MIN_INDEX_SD not applied)\n")

# Restrict dataset to allowed vegetation classes
if ("Veg" %in% names(df) && length(ALLOWED_VEG) > 0) {
  try(
    {
      for (av in ALLOWED_VEG) {
        sel <- grepl(av, df$Veg, ignore.case = TRUE) & !is.na(df$Veg)
        if (any(sel)) {
          df$Veg[sel] <- av
          cat(sprintf("Normalized %d rows to canonical Veg '%s'\n", sum(sel, na.rm = TRUE), av))
        }
      }
    },
    silent = TRUE
  )
}

if ("Veg" %in% names(df) && length(ALLOWED_VEG) > 0) {
  keep_rows <- tolower(df$Veg) %in% ALLOWED_VEG | tolower(df$Veg) == "barren"
  n_before <- nrow(df)
  df <- df[keep_rows | is.na(df$Veg), , drop = FALSE]
  cat(sprintf(
    "Filtered to allowed classes (%s) + barren: kept %d/%d rows\n",
    paste(ALLOWED_VEG, collapse = ","), nrow(df), n_before
  ))
}

# Per-vegetation quick summary
try(
  {
    cat("Per-veg quick summary:\n")
    all_veg_classes <- c(ALLOWED_VEG, "barren")
    for (av in all_veg_classes) {
      sel <- tolower(df$Veg) == av
      rows <- sum(sel, na.rm = TRUE)
      unique_doys <- length(unique(df$doy[sel & is.finite(df$doy)]))
      unique_locs <- length(unique(df$location_id[sel]))
      cat(sprintf("  %s: rows=%d unique_doys=%d unique_locs=%d\n", av, rows, unique_doys, unique_locs))
    }
  },
  silent = TRUE
)

matched_veg_n <- sum(!is.na(df$Veg))
cat("Non-NA Veg rows:", matched_veg_n, "of", nrow(df), "\n")
if (matched_veg_n == 0) {
  stop("No vegetation classes found after join; cannot build library")
}

loc_years <- df %>%
  dplyr::filter(!is.na(.data$location_id) & .data$location_id != "" & !is.na(.data$year) & .data$year > 0 & !is.na(.data$Veg)) %>%
  dplyr::distinct(.data$location_id, .data$year)
cat(sprintf("Constructed loc_years with %d rows from filtered df\n", nrow(loc_years)))
if (nrow(loc_years) == 0) {
  # Hard fail: do not fall back silently to the full dataset. Training data
  # must contain at least one location-year pair after filtering.
  stop(paste0(
    "No location-year pairs found after filtering. This is a fatal error — library construction cannot continue.\n",
    "Possible causes and suggestions:\n",
    " - Your filtered training dataset has zero valid location/year pairs (check column 'location_id' and 'year').\n",
    " - Verify TRAIN_YEARS and any prior filtering steps do not remove all data (e.g. TRAIN_YEARS <- 2019:2024).\n",
    " - Ensure your transformer produced valid 'location_id' values that match your geojson mapping (transform_phenology.py formats 'L_lon_lat').\n",
    " - If your data are intentionally sparse, reduce the filtering thresholds or increase available training data.\n",
    "Processing cannot continue without at least one location-year pair in filtered training data.")
  )
}

# Construct vegetation library
cat("Constructing lib from TRAINING dataset...\n")
lib <- list()
## Soil-first MESMA preprocessing: build a soil prototype from geojson-labelled
## 'no soil' points and subtract estimated soil fraction from every observation
## in the working data. This runs before vegetation library construction.
build_soil_prototype <- function(df_local, avail_idx, threshold = SOIL_PURE_THRESHOLD, min_samples = SOIL_MIN_SAMPLES) {
  # Prefer using explicit 'barren' Veg-labeled rows as soil (bare ground) references.
  # If no 'barren' rows are present, do not attempt to create a soil prototype.
  if (!"Veg" %in% names(df_local)) return(NULL)
  candidates <- df_local[tolower(as.character(df_local$Veg)) == "barren", , drop = FALSE]

  ## Strict requirement: do not fallback to selecting an arbitrary top-N set.
  ## If there are fewer than min_samples candidates that meet the threshold,
  ## return NULL and let the caller handle the absence of a soil prototype.
  if (nrow(candidates) < min_samples) {
    return(NULL)
  }

  if (nrow(candidates) == 0) return(NULL)

  soil_lib <- list()
  # ensure 'doy' is present
  if (!"doy" %in% names(candidates)) candidates$doy <- lubridate::yday(as.Date(candidates$date))

  for (idx in avail_idx) {
    if (!idx %in% names(candidates)) next
    vals_by_doy <- tryCatch(
      {
        tapply(seq_along(candidates[[idx]]), candidates$doy, function(indices) {
          vals <- candidates[[idx]][indices]
          vals <- vals[is.finite(vals)]
          if (length(vals) == 0) return(NA_real_)
          if (length(vals) == 1) return(as.numeric(vals[1]))
          m <- mean(vals)
          as.numeric(vals[which.min(abs(vals - m))])
        })
      },
      error = function(e) NULL
    )

    mu <- rep(NA_real_, 365)
    if (!is.null(vals_by_doy) && length(vals_by_doy) > 0) {
      doy_values <- as.integer(names(vals_by_doy))
      valid_doy <- doy_values >= 1 & doy_values <= 365
      mu[doy_values[valid_doy]] <- vals_by_doy[valid_doy]
    }

    if (all(!is.finite(mu))) next
    mv <- tryCatch(calc_moving_var(data.frame(date = 1:365, idx = mu), "idx", window = 14), error = function(e) rep(NA_real_, 365))
    soil_lib[[idx]] <- list(mu = mu, mv = mv)
  }

  attr(soil_lib, "rows_used") <- nrow(candidates)
  soil_lib
}

## Build a two-endmember library for stage 1 MESMA: barren (no soil==0 OR Veg=='barren') vs pure vegetation (no soil==1)
## This library is used to unmix the vegetated fraction for all pixels before decomposing vegetation types
build_barren_veg_library <- function(df_local, avail_idx, min_samples = 5) {
  if (!"doy" %in% names(df_local)) df_local$doy <- lubridate::yday(as.Date(df_local$date))
  
  # Extract barren endmember using BOTH methods:
  # 1. no soil ≈ 0 (if available)
  # 2. Veg == 'barren' (fallback/additional)
  barren_by_nosoil <- if ("no soil" %in% names(df_local)) {
    df_local[!is.na(df_local$`no soil`) & {
      val <- df_local$`no soil`
      if (is.character(val)) val <- as.numeric(val)
      abs(val - 0) < 0.01
    }, , drop = FALSE]
  } else {
    df_local[FALSE, , drop = FALSE]  # empty
  }
  
  barren_by_veg <- if ("Veg" %in% names(df_local)) {
    df_local[!is.na(df_local$Veg) & tolower(df_local$Veg) == "barren", , drop = FALSE]
  } else {
    df_local[FALSE, , drop = FALSE]  # empty
  }
  
  # Combine both sources (remove duplicates by row)
  barren_rows <- unique(rbind(barren_by_nosoil, barren_by_veg))
  
  cat(sprintf("[Stage1] Barren rows: %d from no_soil==0, %d from Veg=='barren', %d total unique\n",
              nrow(barren_by_nosoil), nrow(barren_by_veg), nrow(barren_rows)))
  
  # HARD REQUIREMENT: Must have at least min_samples barren observations
  if (nrow(barren_rows) < min_samples) {
    # stop(sprintf("[Stage1] CRITICAL ERROR: Insufficient barren training data. Found %d barren observations, but require at least %d. Cannot proceed with MESMA analysis.", nrow(barren_rows), min_samples))
    cat(sprintf("[Stage1] Insufficient barren training data. Found %d barren observations, but require at least %d. Returning NULL.\n", nrow(barren_rows), min_samples))
    return(NULL)
  }
  
  # Extract pure vegetation endmember (no soil ≈ 1)
  veg_rows <- if ("no soil" %in% names(df_local)) {
    df_local[!is.na(df_local$`no soil`) & {
      val <- df_local$`no soil`
      if (is.character(val)) val <- as.numeric(val)
      abs(val - 1) < 0.01
    }, , drop = FALSE]
  } else {
    # Fallback: use non-barren vegetation rows
    df_local[!is.na(df_local$Veg) & tolower(df_local$Veg) != "barren", , drop = FALSE]
  }
  
  cat(sprintf("[Stage1] After filtering: barren=%d, veg=%d rows\n", 
              nrow(barren_rows), nrow(veg_rows)))
  
  if (nrow(barren_rows) < min_samples) {
     cat(sprintf("[Stage1] Insufficient training data: barren=%d (need >=%d)\n", nrow(barren_rows), min_samples))
     return(NULL)
  }

  if (nrow(veg_rows) < min_samples) {
    cat(sprintf("[Stage1] Insufficient training data: veg=%d (need >=%d)\n", 
                nrow(veg_rows), min_samples))
    return(NULL)
  }
  
  cat(sprintf("[Stage1] Building barren-veg library from barren=%d, veg=%d rows\n", 
              nrow(barren_rows), nrow(veg_rows)))
  
  stage1_lib <- list()
  
  # Build barren endmember
  barren_lib <- list()
  for (idx in avail_idx) {
    if (!idx %in% names(barren_rows)) next
    vals_by_doy <- tapply(seq_along(barren_rows[[idx]]), barren_rows$doy, function(indices) {
      vals <- barren_rows[[idx]][indices]
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) return(NA_real_)
      median(vals)
    })
    mu <- rep(NA_real_, 365)
    if (length(vals_by_doy) > 0) {
      doy_values <- as.integer(names(vals_by_doy))
      valid_doy <- doy_values >= 1 & doy_values <= 365
      mu[doy_values[valid_doy]] <- vals_by_doy[valid_doy]
    }
    if (all(!is.finite(mu))) next
    mv <- tryCatch(calc_moving_var(data.frame(date = 1:365, idx = mu), "idx", window = 14), 
                   error = function(e) rep(NA_real_, 365))
    barren_lib[[idx]] <- list(mu = mu, mv = mv)
  }
  
  # Build vegetation endmember
  veg_lib <- list()
  for (idx in avail_idx) {
    if (!idx %in% names(veg_rows)) next
    vals_by_doy <- tapply(seq_along(veg_rows[[idx]]), veg_rows$doy, function(indices) {
      vals <- veg_rows[[idx]][indices]
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) return(NA_real_)
      median(vals)
    })
    mu <- rep(NA_real_, 365)
    if (length(vals_by_doy) > 0) {
      doy_values <- as.integer(names(vals_by_doy))
      valid_doy <- doy_values >= 1 & doy_values <= 365
      mu[doy_values[valid_doy]] <- vals_by_doy[valid_doy]
    }
    # Interpolate missing DOYs for complete temporal coverage
    if (sum(is.finite(mu)) > 2) {
      mu <- approx(x = which(is.finite(mu)), y = mu[is.finite(mu)], 
                   xout = 1:365, rule = 2)$y
    }
    if (all(!is.finite(mu))) next
    mv <- tryCatch(calc_moving_var(data.frame(date = 1:365, idx = mu), "idx", window = 14), 
                   error = function(e) rep(NA_real_, 365))
    veg_lib[[idx]] <- list(mu = mu, mv = mv)
  }
  
  if (length(barren_lib) == 0 || length(veg_lib) == 0) return(NULL)
  
  stage1_lib$barren <- barren_lib
  stage1_lib$vegetation <- veg_lib
  stage1_lib
}

## Unmix vegetated fraction using stage 1 MESMA (barren vs vegetation endmembers)
## Returns the estimated vegetated fraction for this location-year
## Uses GEOMETRIC PROJECTION (Tits et al. paper method) for all unmixing
unmix_vegetated_fraction <- function(dly_local, stage1_lib, avail_idx) {
  if (is.null(stage1_lib) || length(stage1_lib) == 0) return(NA_real_)
  if (!"doy" %in% names(dly_local)) dly_local$doy <- lubridate::yday(as.Date(dly_local$date))
  
  # For each row, perform geometric projection unmixing
  veg_fractions <- numeric(nrow(dly_local))
  valid_count <- 0
  
  for (r in seq_len(nrow(dly_local))) {
    row <- dly_local[r, , drop = FALSE]
    doy <- as.integer(row$doy)
    if (is.na(doy) || doy < 1 || doy > 365) {
      veg_fractions[r] <- NA_real_
      next
    }
    
    # Build vectors for this DOY
    y <- numeric(0)  # Observation
    B <- numeric(0)  # Barren endmember
    V <- numeric(0)  # Vegetation endmember
    
    for (idx in avail_idx) {
      if (!idx %in% names(row)) next
      
      y_val <- as.numeric(row[[idx]])
      b_val <- if (!is.null(stage1_lib$barren[[idx]]) && 
                   !is.null(stage1_lib$barren[[idx]]$mu)) {
        stage1_lib$barren[[idx]]$mu[doy]
      } else NA_real_
      v_val <- if (!is.null(stage1_lib$vegetation[[idx]]) && 
                   !is.null(stage1_lib$vegetation[[idx]]$mu)) {
        stage1_lib$vegetation[[idx]]$mu[doy]
      } else NA_real_
      
      # Only include if all three values are valid
      if (is.finite(y_val) && is.finite(b_val) && is.finite(v_val)) {
        y <- c(y, y_val)
        B <- c(B, b_val)
        V <- c(V, v_val)
      }
    }
    
    # Need at least 2 dimensions for meaningful unmixing
    if (length(y) < 2) {
      veg_fractions[r] <- NA_real_
      next
    }
    
    # Check endmember separation for this DOY
    separation <- sqrt(sum((B - V)^2))
    if (separation < 0.01) {
      # Endmembers too similar at this DOY - skip
      veg_fractions[r] <- NA_real_
      next
    }
    
    # GEOMETRIC PROJECTION: Project y onto line from B to V (no normalization)
    result <- geometric_unmix_two_endmembers(y, B, V)
    
    # f2 is the fraction of V (vegetation)
    veg_fractions[r] <- result$f2
    valid_count <- valid_count + 1
  }
  
  # Return median vegetated fraction across all valid observations
  valid_fracs <- veg_fractions[is.finite(veg_fractions)]
  
  if (length(valid_fracs) == 0) {
    cat(sprintf("[Stage1 Geometric RAW] No valid unmixing results (checked %d rows)\n", nrow(dly_local)))
    return(NA_real_)
  }
  
  veg_frac <- median(valid_fracs)
  
  cat(sprintf("[Stage1 Geometric RAW] Veg fraction: %.3f (from %d/%d valid obs)\n", 
              veg_frac, length(valid_fracs), nrow(dly_local)))
  
  veg_frac
}

 


# Select indices that maximize barren-vegetation separation for Stage 1
select_stage1_indices <- function(barren_lib, veg_lib, avail_idx, min_separation = 0.1) {
  separations <- sapply(avail_idx, function(idx) {
    if (!idx %in% names(barren_lib) || !idx %in% names(veg_lib)) return(0)
    b_mean <- mean(barren_lib[[idx]]$mu, na.rm = TRUE)
    v_mean <- mean(veg_lib[[idx]]$mu, na.rm = TRUE)
    abs(v_mean - b_mean)
  })
  # Keep indices with separation above threshold
  keep <- avail_idx[separations >= min_separation]
  if (length(keep) < 3) keep <- avail_idx[order(separations, decreasing = TRUE)[1:min(5, length(avail_idx))]]
  keep
}


# Compress Stage 1 Library (Barren/Veg prototypes) using raw spectral indices
build_stage1_lib <- function(STAGE1_LIB, grid_name = "full") {
  if (is.null(STAGE1_LIB)) return(NULL)
  
  barren_raw <- STAGE1_LIB$barren
  vegetation_raw <- STAGE1_LIB$vegetation
  
  # Get available indices
  avail_idx <- intersect(names(barren_raw), names(vegetation_raw))
  if (length(avail_idx) == 0) return(NULL)
  
  # Select indices that maximize barren-vegetation separation
  avail_idx <- select_stage1_indices(barren_raw, vegetation_raw, avail_idx, min_separation = 0.1)
  cat(sprintf("[build_stage1_lib] Selected %d indices with good barren-veg separation\n", length(avail_idx)))
  
  # Helper function to build a 365-day raw index matrix
  build_raw_trace <- function(endmember_lib, avail_idx) {
    # Build raw 365 x n_indices matrix
    raw_mat <- matrix(NA_real_, nrow = 365, ncol = length(avail_idx))
    colnames(raw_mat) <- avail_idx
    
    for (j in seq_along(avail_idx)) {
      idx <- avail_idx[j]
      if (!is.null(endmember_lib[[idx]]) && !is.null(endmember_lib[[idx]]$mu)) {
        raw_mat[, j] <- endmember_lib[[idx]]$mu
      }
    }
    
    # Interpolate missing days
    for (j in seq_len(ncol(raw_mat))) {
      col_vals <- raw_mat[, j]
      finite_idx <- which(is.finite(col_vals))
      if (length(finite_idx) >= 2) {
        raw_mat[is.na(raw_mat[, j]), j] <- approx(finite_idx, col_vals[finite_idx], 
                                                   xout = which(is.na(raw_mat[, j])), rule = 2)$y
      } else if (length(finite_idx) == 1) {
        raw_mat[, j] <- col_vals[finite_idx]
      } else {
        raw_mat[, j] <- 0
      }
    }
    raw_mat[!is.finite(raw_mat)] <- 0
    
    raw_mat
  }
  
  # Build raw index traces for barren and vegetation (full temporal resolution)
  raw_barren <- build_raw_trace(barren_raw, avail_idx)
  raw_veg <- build_raw_trace(vegetation_raw, avail_idx)
  
  # Check similarity between barren and vegetation endmembers
  sim <- cos_sim(as.numeric(raw_barren), as.numeric(raw_veg))
  cat(sprintf("[build_stage1_lib] Full resolution: %d days × %d indices. Cosine Similarity(Barren, Veg) = %.4f\n", nrow(raw_barren), ncol(raw_barren), sim))

  # Warn if the endmembers are too similar (unmixing may be unreliable)
  if (is.finite(sim) && sim > 0.95) {
    warning(sprintf("[build_stage1_lib] Barren and Vegetation endmembers are very similar (cos_sim=%.4f). Stage 1 unmixing may be unreliable.", sim))
  }

  # Also check Euclidean separation
  euclidean_dist <- sqrt(sum((as.numeric(raw_barren) - as.numeric(raw_veg))^2))
  cat(sprintf("[build_stage1_lib] Euclidean distance(Barren, Veg) = %.4f\n", euclidean_dist))
  if (euclidean_dist < 0.1) {
    warning("[build_stage1_lib] Barren and Vegetation endmembers have very small separation.")
  }
  # Check gross scale differences across all features; warn if the overall L2 norms differ a lot
  norm_b <- sqrt(sum(as.numeric(raw_barren)^2, na.rm = TRUE))
  norm_v <- sqrt(sum(as.numeric(raw_veg)^2, na.rm = TRUE))
  if (is.finite(norm_b) && is.finite(norm_v) && (norm_b > 0) && (norm_v > 0)) {
    ratio <- max(norm_b / norm_v, norm_v / norm_b)
    if (ratio > 10) {
      warning(sprintf("[build_stage1_lib] Barren and Vegetation full-trace norms differ substantially (barren=%.4g, veg=%.4g, ratio=%.4g). Consider normalizing stage1 endmembers to avoid projection bias.", norm_b, norm_v, ratio))
    }
  }
  
  list(barren = raw_barren, vegetation = raw_veg)
}

# =============================================================================
# OPTIMIZED LIBRARY PRE-COMPUTATION
# =============================================================================
precompute_optimized_library <- function(mesma_lib, compressed_templates_accessor, grid_type, feature_weights = NULL) {
  opt_lib <- list()
  
  for (veg in names(mesma_lib)) {
    variants <- mesma_lib[[veg]]
    n_vars <- length(variants)
    if (n_vars == 0) next
    
    # Collect all vectors
    vec_list <- list()
    ids <- character(n_vars)
    valid_idx <- integer(0)
    
    for (i in seq_along(variants)) {
      vid <- variants[[i]]$variant_id
      vec <- compressed_templates_accessor[[veg]][[vid]][["full"]]
      if (!is.null(vec)) {
        vec_list[[length(vec_list) + 1]] <- vec
        ids[length(vec_list)] <- vid
        valid_idx <- c(valid_idx, i)
      }
    }
    
    if (length(vec_list) == 0) next
    
    # Create matrix (Rows = Variants, Cols = Features)
    M <- do.call(rbind, vec_list)
    ids <- ids[1:length(vec_list)]
    
    # Pre-compute weighted normalized matrix for fast cosine similarity
    M_weighted <- M
    if (!is.null(feature_weights)) {
      if (length(feature_weights) == ncol(M)) {
        # Apply weights to each row
        M_weighted <- t(t(M) * feature_weights)
      }
    }
    
    # Normalize rows to unit length
    row_norms <- sqrt(rowSums(M_weighted^2))
    # Avoid division by zero
    row_norms[row_norms < 1e-9] <- 1
    M_norm <- M_weighted / row_norms
    
    opt_lib[[veg]] <- list(
      M = M,              # Original vectors (for unmixing)
      M_norm = M_norm,    # Weighted & Normalized vectors (for similarity search)
      ids = ids
    )
  }
  return(opt_lib)
}

# =============================================================================
# GEOMETRIC Stage 2 Unmixing
# Uses angle-based endmember selection + projection as per Tits et al. paper
# =============================================================================
unmix_stage2_compressed <- function(veg_kept, veg_frac, y, grid_type, compressed_templates_accessor, mesma_lib, topK = TOPK_VARIANTS, feature_weights = NULL, optimized_library = NULL) {
  # Track verbose output to first 3 calls only
  if (!exists(".GEOM_DEBUG_COUNTER", envir = globalenv())) assign(".GEOM_DEBUG_COUNTER", 0L, envir = globalenv())
  debug_counter <- get(".GEOM_DEBUG_COUNTER", envir = globalenv())
  verbose_this_call <- (debug_counter < 5)
  if (verbose_this_call) assign(".GEOM_DEBUG_COUNTER", debug_counter + 1L, envir = globalenv())
  
  veg_types <- veg_kept
  if (length(veg_types) == 0) return(NULL)
  
  # Filter veg_kept to valid vegetation types in mesma_lib
  veg_types <- intersect(veg_types, names(mesma_lib))
  if (length(veg_types) == 0) return(NULL)
  
  top_variants <- list()
  
  # FAST PATH: Use Optimized Library if available
  if (!is.null(optimized_library)) {
    y_vec <- as.numeric(y)
    
    # Prepare y for similarity search
    y_weighted <- y_vec
    if (!is.null(feature_weights) && length(feature_weights) == length(y_vec)) {
      y_weighted <- y_vec * feature_weights
    }
    y_norm <- sqrt(sum(y_weighted^2))
    if (y_norm < 1e-9) y_norm <- 1
    y_weighted_norm <- y_weighted / y_norm
    
    for (v in veg_types) {
      lib_v <- optimized_library[[v]]
      if (is.null(lib_v)) next
      
      # Vectorized Cosine Similarity: M_norm %*% y_weighted_norm
      # M_norm is (N_vars x N_feats), y is (N_feats x 1) -> (N_vars x 1)
      sims <- as.numeric(lib_v$M_norm %*% y_weighted_norm)
      
      # Select Top K
      ord <- order(sims, decreasing = TRUE)
      keep_idx <- ord[seq_len(min(topK, length(ord)))]
      
      # Retrieve original vectors for unmixing
      selected_vecs <- lib_v$M[keep_idx, , drop = FALSE]
      selected_ids <- lib_v$ids[keep_idx]
      
      # Format as list for stage2_geometric_unmix
      variant_list <- list()
      for (k in seq_along(keep_idx)) {
        variant_list[[k]] <- list(vec = selected_vecs[k, ], id = selected_ids[k])
      }
      top_variants[[v]] <- variant_list
    }
    
  } else {
    # SLOW PATH: Legacy list-based approach
    # Access compressed templates as nested list
    comp_templates <- list()
    for (v in veg_types) {
      comp_templates[[v]] <- list()
      for (variant in mesma_lib[[v]]) {
        vid <- variant$variant_id
        vec <- compressed_templates_accessor[[v]][[vid]][[grid_type]]
        # if (verbose_this_call) cat(sprintf("[DEBUG Stage2 Geometric] v=%s, vid=%s, vec length=%d\n", v, vid, length(vec)))
        if (!is.null(vec) && length(vec) > 0) {
          comp_templates[[v]][[length(comp_templates[[v]]) + 1]] <- list(vec = vec, id = vid)
        }
      }
      # if (verbose_this_call) cat(sprintf("[DEBUG unmix_stage2_compressed] comp_templates[[%s]] has %d variants\n", v, length(comp_templates[[v]])))
    }
    
    if (length(comp_templates) == 0) return(NULL)
    
    for (v in veg_types) {
      cand <- comp_templates[[v]]
      if (length(cand) == 0) next
      sims <- sapply(cand, function(x) {
        if (!is.null(feature_weights)) {
          # Use weighted cosine similarity
          weighted_cosine_similarity(y, x$vec, feature_weights)
        } else {
          cos_sim(y, x$vec)
        }
      })
      ord <- order(sims, decreasing = TRUE)
      keep <- ord[seq_len(min(topK, length(ord)))]
      top_variants[[v]] <- cand[keep]
      # if (verbose_this_call) cat(sprintf("[DEBUG Stage2 Geometric] top_variants[[%s]] has %d variants\n", v, length(top_variants[[v]])))
    }
  }
  
  if (length(top_variants) == 0) {
    # if (verbose_this_call) cat("[DEBUG Stage2 Geometric] No top variants, returning NULL\n")
    return(NULL)
  }
  
  # ========== GEOMETRIC UNMIXING (Paper's Method) ==========
  # Use angle-based endmember selection + geometric projection
  
  y <- as.numeric(y)
  
  # Build vegetation library for geometric unmixing
  veg_libraries <- list()
  for (v in names(top_variants)) {
    veg_libraries[[v]] <- top_variants[[v]]
  }
  
  # Call the geometric stage 2 unmixing
  geom_result <- stage2_geometric_unmix(y, veg_libraries, topK = topK, feature_weights = feature_weights)
  
  if (is.null(geom_result)) {
    # if (verbose_this_call) cat("[DEBUG Stage2 Geometric] stage2_geometric_unmix returned NULL\n")
    return(NULL)
  }
  
  # Extract results
  w <- geom_result$fractions
  chosen <- geom_result$chosen_variants
  rmse <- geom_result$residual
  
  # Apply minimum fraction threshold (10%)
  # Set fractions below 10% to 0 and renormalize
  MIN_VEG_FRACTION <- 0.10
  w[w < MIN_VEG_FRACTION] <- 0
  if (sum(w) > 0) {
    w <- w / sum(w)
  }
  
  if (verbose_this_call) {
    cat(sprintf("\n[Stage2 Geometric] Fractions: %s\n", 
                paste(sprintf("%s=%.4f", names(w), w), collapse=", ")))
    cat(sprintf("[Stage2 Geometric] Chosen variants: %s\n", 
                paste(sprintf("%s=%s", names(chosen), chosen), collapse=", ")))
    cat(sprintf("[Stage2 Geometric] Residual: %.6f\n", rmse))
  }
  
  best <- list(rmse = rmse, w = w, chosen = as.list(chosen))
  
  diagnostics <- NULL
  uncertainty <- NULL
  
  if (isTRUE(ENABLE_DIAGNOSTICS) && !is.null(best$w) && !is.null(best$chosen)) {
    # Build endmember matrix for diagnostics
    E_best_cols <- list()
    for (v in names(best$chosen)) {
      for (variant in top_variants[[v]]) {
        if (variant$id == best$chosen[[v]]) {
          E_best_cols[[v]] <- variant$vec
          break
        }
      }
    }
    if (length(E_best_cols) > 0) {
      E_best <- do.call(cbind, E_best_cols)
      for (j in seq_len(ncol(E_best))) {
        nj <- sqrt(sum(E_best[, j]^2))
        if (nj > 0) E_best[, j] <- E_best[, j] / nj
      }
      diagnostics <- tryCatch({
        compute_diagnostics(y, E_best, best$w, mesma_result = best)
      }, error = function(e) NULL)
    }
  }

  if (isTRUE(ENABLE_UNCERTAINTY) && !is.null(best$w) && !is.null(best$chosen)) {
    # Enforce nested two-stage bootstrap only
    if (is.null(comp_templates) || is.null(compressed_stage1_lib)) {
      stop("ENABLED UNCERTAINTY requires nested two-stage bootstrap: compressed stage1 lib and comp_templates must be available")
    }
    if (!exists("nested_two_stage_bootstrap")) {
      stop("nested_two_stage_bootstrap function not found - nested bootstrap is required for ENABLE_UNCERTAINTY = TRUE")
    }
    uncertainty <- nested_two_stage_bootstrap(
      y_vec = y,
      compressed_stage1_lib = compressed_stage1_lib,
      comp_templates = comp_templates,
      top_variants = top_variants,
      chosen_ids = best$chosen,
      w_hat = best$w,
      B = BOOTSTRAP_B,
      lambda_star = 1e-6,
      seed = 123,
      variant_switch = VARIANT_SWITCH_BOOTSTRAP,
      re_center = TRUE
    )
  }
  
  list(vegetation_proportions = best$w, chosen_variants = best$chosen, rmse = best$rmse, diagnostics = diagnostics, uncertainty = uncertainty)
}

#' Geometric block bootstrap for uncertainty estimation
#' Uses geometric unmixing for weight estimation
geometric_block_bootstrap <- function(y_vec, comp_templates, top_variants, chosen_ids, w_hat, 
                                        B = 100, seed = 123, variant_switch = FALSE) {
  stop("geometric_block_bootstrap is disabled: nested two-stage bootstrap is the only supported uncertainty method when ENABLE_UNCERTAINTY = TRUE")
  set.seed(seed)
  
  veg_names <- names(w_hat)
  n_veg <- length(veg_names)
  
  if (n_veg == 0) return(NULL)
  
  # Build vegetation library
  veg_libraries <- list()
  for (v in veg_names) {
    veg_libraries[[v]] <- top_variants[[v]]
  }
  
  # Bootstrap results storage
  boot_weights <- matrix(NA_real_, nrow = B, ncol = n_veg)
  colnames(boot_weights) <- veg_names
  boot_residuals <- numeric(B)
  
  n <- length(y_vec)
  block_size <- max(3, floor(sqrt(n)))
  
  for (b in seq_len(B)) {
    # Generate bootstrap sample
    n_blocks <- ceiling(n / block_size)
    block_starts <- sample(1:(n - block_size + 1), n_blocks, replace = TRUE)
    boot_idx <- unlist(lapply(block_starts, function(s) s:(s + block_size - 1)))
    boot_idx <- boot_idx[boot_idx <= n][1:n]
    
    y_boot <- y_vec[boot_idx]
    
    # Geometric unmixing on bootstrap sample
    if (variant_switch) {
      # Allow variant switching in bootstrap
      result <- stage2_geometric_unmix(y_boot, veg_libraries, topK = 2)
    } else {
      # Fixed variants - just re-unmix with same endmembers
      M_cols <- list()
      for (v in veg_names) {
        for (variant in top_variants[[v]]) {
          if (variant$id == chosen_ids[[v]]) {
            M_cols[[v]] <- variant$vec[boot_idx]
            break
          }
        }
      }
      if (length(M_cols) < n_veg) next
      M <- do.call(cbind, M_cols)
      result <- geometric_unmix_n_endmembers(y_boot, M)
      result <- list(fractions = setNames(result$f, veg_names), residual = result$residual)
    }
    
    if (!is.null(result) && !is.null(result$fractions)) {
      for (v in veg_names) {
        if (v %in% names(result$fractions)) {
          boot_weights[b, v] <- result$fractions[v]
        }
      }
      boot_residuals[b] <- result$residual
    }
  }
  
  # Compute CIs
  coef_ci <- data.frame(
    Veg = veg_names,
    coef_025 = apply(boot_weights, 2, quantile, 0.025, na.rm = TRUE),
    coef_975 = apply(boot_weights, 2, quantile, 0.975, na.rm = TRUE),
    coef_sd = apply(boot_weights, 2, sd, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  rmse_ci <- quantile(boot_residuals, c(0.025, 0.975), na.rm = TRUE)
  
  w_median <- apply(boot_weights, 2, median, na.rm = TRUE)
  names(w_median) <- veg_names
  
  list(
    coef_ci = coef_ci,
    rmse_ci = rmse_ci,
    w_median = w_median,
    block_size = block_size,
    n_valid = sum(!is.na(boot_weights[, 1]))
  )
}

# Build stage 1 library for vegetated fraction unmixing (barren vs pure vegetation)
STAGE1_LIB <- NULL
try({
  STAGE1_LIB <- build_barren_veg_library(df, avail, min_samples = 5)
}, silent = TRUE)

if (is.null(STAGE1_LIB)) {
  stop("[Stage1] Could not build barren-veg library (insufficient no soil==0 or no soil==1 rows). Two-stage MESMA requires sufficient stage-1 training data. Cannot proceed with spectral-only unmixing.")
} else {
  cat("[Stage1] Barren-vegetation library built successfully for stage 1 unmixing\n")
  
  # ========== STAGE 1 DIAGNOSTIC 6: Check Stage 1 Endmember Quality ==========
  cat("\n=== STAGE 1 LIBRARY CHECK ===\n")
  
  # Check barren endmember
  barren_lib <- STAGE1_LIB$barren
  barren_vals <- sapply(barren_lib, function(x) mean(x$mu, na.rm = TRUE))
  cat(sprintf("Barren: mean across indices = %.4f (range: [%.4f, %.4f])\n",
              mean(barren_vals), min(barren_vals), max(barren_vals)))
  
  # Check vegetation endmember
  veg_lib <- STAGE1_LIB$vegetation
  veg_vals <- sapply(veg_lib, function(x) mean(x$mu, na.rm = TRUE))
  cat(sprintf("Vegetation: mean across indices = %.4f (range: [%.4f, %.4f])\n",
              mean(veg_vals), min(veg_vals), max(veg_vals)))
  
  # Check separation per index
  cat("\nPer-index separation:\n")
  for (idx in names(barren_lib)) {
    if (idx %in% names(veg_lib)) {
      b_mean <- mean(barren_lib[[idx]]$mu, na.rm = TRUE)
      v_mean <- mean(veg_lib[[idx]]$mu, na.rm = TRUE)
      separation <- abs(v_mean - b_mean)
      cat(sprintf("  %s: barren=%.4f, veg=%.4f, sep=%.4f\n", 
                  idx, b_mean, v_mean, separation))
      
      if (separation < 0.05) {
        warning(sprintf("Index %s has low barren-veg separation!", idx))
      }
    }
  }
  cat("=============================\n\n")
  # ========== END STAGE 1 DIAGNOSTICS ==========
}

lib_df <- df
vegs <- unique(na.omit(lib_df$Veg))
vegs <- vegs[vegs != ""]
vegs <- vegs[vegs %in% c("phragmites", "populus", "tamarix")]

if (ENABLE_SAMPLE_BALANCING) {
  # Sample Size Balancing - Aggressive Equalization
  lib_df <- lib_df[!is.na(lib_df$Veg), ]
  cat("Applying aggressive sample size balancing to training data...\n")
  lib_df$doy_for_lib <- lubridate::yday(lib_df$date)

  veg_sample_sizes <- sapply(vegs, function(v) {
    dveg <- lib_df[lib_df$Veg == v & is.finite(lib_df$doy_for_lib), , drop = FALSE]
    nrow(dveg)
  })

  min_samples <- min(veg_sample_sizes)
  max_samples <- max(veg_sample_sizes)
  median_samples <- median(veg_sample_sizes)

  cat(sprintf(
    "Original sample sizes - Min: %.0f, Max: %.0f, Median: %.1f\n",
    min_samples, max_samples, median_samples
  ))

  # Calculate the final balancing size as 2x the minimum original sample size
  final_balancing_size <- 2 * min(veg_sample_sizes)

  cat(sprintf("Final balancing size (2x minimum original sample size): %d\n", final_balancing_size))

# SMOTE-like augmentation for minority classes
augment_minority_class <- function(df_class, target_n, seed = NULL, alpha_range = c(0.3, 0.7), jitter_frac = 1e-6) {
  if (!is.data.frame(df_class)) stop("augment_minority_class requires a data.frame")
  if (nrow(df_class) >= target_n) {
    return(df_class)
  }
  if (nrow(df_class) <= 1) {
    reps <- target_n - nrow(df_class)
    base <- df_class[rep(1, reps), , drop = FALSE]
    num_cols <- sapply(base, is.numeric)
    if (any(num_cols)) {
      rng <- apply(df_class[, num_cols, drop = FALSE], 2, function(x) {
        finite_x <- x[is.finite(x)]
        if (length(finite_x) < 2) 0 else diff(range(finite_x))
      })
      rng[rng == 0] <- 1.0
      base[, num_cols] <- base[, num_cols, drop = FALSE] +
        matrix(rnorm(reps * sum(num_cols), sd = sqrt(mean(rng)) * jitter_frac), nrow = reps)
    }
    out <- rbind(df_class, base)
    rownames(out) <- NULL
    return(out)
  }

  if (!is.null(seed)) set.seed(seed)
  n_add <- target_n - nrow(df_class)
  synth_rows <- vector("list", n_add)

  num_cols <- names(df_class)[sapply(df_class, is.numeric)]
  other_cols <- setdiff(names(df_class), num_cols)

  for (i in seq_len(n_add)) {
    ids <- sample.int(nrow(df_class), 2, replace = FALSE)
    a <- df_class[ids[1], , drop = FALSE]
    b <- df_class[ids[2], , drop = FALSE]
    alpha <- runif(1, min(alpha_range), max(alpha_range))
    new_row <- a

    if (length(num_cols) > 0) {
      va <- as.numeric(a[1, num_cols, drop = FALSE])
      vb <- as.numeric(b[1, num_cols, drop = FALSE])
      interp <- (1 - alpha) * va + alpha * vb
      rng <- apply(df_class[, num_cols, drop = FALSE], 2, function(x) {
        finite_x <- x[is.finite(x)]
        if (length(finite_x) < 2) 0 else diff(range(finite_x))
      })
      jitter <- rnorm(length(interp), sd = pmax(abs(rng), 1e-8) * jitter_frac)
      interp <- interp + jitter
      # Assign numeric interpolated values element-wise to avoid tibble/vctrs
      # assignment errors when new_row is a tibble-like object
      if (length(num_cols) == length(interp)) {
        for (k in seq_along(num_cols)) {
          colname <- num_cols[k]
          new_row[[colname]] <- interp[k]
        }
      } else {
        # Do not silently recycle/interpolate when shapes don't match: this
        # hides subtle bugs in numeric column handling. Fail loudly so calling
        # code can be fixed instead of producing silently corrupted samples.
        stop(sprintf("augment_minority_class: numeric assignment length mismatch (num_cols=%d, interp_len=%d)", length(num_cols), length(interp)))
      }
    }

    for (col in other_cols) {
      if (col %in% c("date", "doy", "doy_for_lib")) next
      # Preserve Veg column (do not interpolate label)
      if (col == "Veg") {
        new_row[[col]] <- a[[col]]
        next
      }
      val <- if (runif(1) < 0.5) a[[col]] else b[[col]]
      new_row[[col]] <- val
    }

    if ("doy_for_lib" %in% names(df_class)) {
      da <- as.integer(a$doy_for_lib)
      db <- as.integer(b$doy_for_lib)
      if (is.finite(da) && is.finite(db)) {
        diff_ab <- ((db - da + 365) %% 365)
        pick_doy <- ((da + round(alpha * diff_ab) - 1) %% 365) + 1
        new_row$doy_for_lib <- as.integer(pick_doy)
      } else {
        new_row$doy_for_lib <- if (!is.na(da) && is.finite(da)) da else if (!is.na(db) && is.finite(db)) db else NA_integer_
      }
    }

    if ("date" %in% names(df_class)) {
      if (!is.null(new_row$doy_for_lib) && is.finite(new_row$doy_for_lib)) {
        new_row$date <- as.Date(sprintf("1970-01-01")) + (as.integer(new_row$doy_for_lib) - 1)
      } else {
        new_row$date <- NA
      }
    }

    synth_rows[[i]] <- new_row[1, , drop = FALSE]
  }

  synth_df <- chunked_rbind(synth_rows, chunk_size = 50L)
  gc()
  out <- rbind(df_class, synth_df)
  for (col in other_cols) {
    if (is.factor(df_class[[col]])) out[[col]] <- factor(out[[col]], levels = levels(df_class[[col]]))
  }
  rownames(out) <- NULL
  out
}

balanced_dfs <- list()
for (v in vegs) {
  dveg <- lib_df[lib_df$Veg == v & is.finite(lib_df$doy_for_lib), , drop = FALSE]
  if (nrow(dveg) == 0) next
  if (isTRUE(TESTING_MODE)) {
    # Cap samples at 5,000 per class for testing
    n_keep <- min(5000L, nrow(dveg))
    if (n_keep < nrow(dveg)) {
      dveg <- dveg[sample(1:nrow(dveg), n_keep), , drop = FALSE]
    }
    balanced_dfs[[v]] <- dveg
    cat(sprintf("Capped %s to %d samples for testing\n", v, nrow(dveg)))
  } else {
    # Aggressive augmentation to exactly match final balancing size
    orig_n <- nrow(dveg)
    target_for_class <- final_balancing_size

    if (target_for_class > orig_n) {
      # Augment up to the target
      augmented <- tryCatch(
        {
          augment_minority_class(dveg, target_for_class)
        },
        error = function(e) stop(sprintf("augment_minority_class failed for %s: %s", v, e$message))
      )
      balanced_dfs[[v]] <- augmented
      cat(sprintf("Augmented %s from %d to %d samples\n", v, orig_n, nrow(augmented)))
    } else {
      # Downsample if needed to reach target
      if (orig_n > target_for_class) {
        set.seed(123)
        dveg <- dveg[sample(1:nrow(dveg), target_for_class), , drop = FALSE]
        cat(sprintf("Downsampled %s from %d to %d samples\n", v, orig_n, target_for_class))
      }
      balanced_dfs[[v]] <- dveg
      cat(sprintf("Using %s with %d samples\n", v, nrow(dveg)))
    }
  }
}

if (length(balanced_dfs) > 0) {
  # All classes are now balanced to exactly final_balancing_size
  lib_df <- chunked_rbind(balanced_dfs, chunk_size = 50L)
  # Safeguard: chunked_rbind may return a matrix; convert to data.frame to preserve column access
  if (is.matrix(lib_df)) {
    lib_df <- as.data.frame(lib_df, stringsAsFactors = FALSE)
  }
  gc()
  cat(sprintf(
    "Training data balanced: %d total samples across %d vegetation types (exactly equal)\n",
    nrow(lib_df), length(balanced_dfs)
  ))
}
} else {
  cat("Sample size balancing disabled.\n")
}

# ===========================================================================
# Stage 1 balancing: ensure barren vs pure-vegetation classes are balanced
# ===========================================================================
if (ENABLE_SAMPLE_BALANCING) {
  cat("Applying Stage-1 (barren vs vegetation pure) sample size balancing...\n")

  # Use the original full dataframe (df) for stage 1 selection because some
  # barren observations rely on 'no soil' flag and may not have a 'Veg' label.
  if (!exists("df")) df <- lib_df

  # Identify barren rows consistent with build_barren_veg_library
  barren_by_nosoil <- if ("no soil" %in% names(df)) {
    df[!is.na(df$`no soil`) & { val <- df$`no soil`; if (is.character(val)) val <- as.numeric(val); abs(val - 0) < 0.01 }, , drop = FALSE]
  } else df[FALSE, , drop = FALSE]
  barren_by_veg <- if ("Veg" %in% names(df)) {
    df[!is.na(df$Veg) & tolower(df$Veg) == "barren", , drop = FALSE]
  } else df[FALSE, , drop = FALSE]
  barren_rows <- unique(rbind(barren_by_nosoil, barren_by_veg))

  # Vegetation pure rows (no soil == 1 or non-barren Veg)
  veg_by_nosoil <- if ("no soil" %in% names(df)) {
    df[!is.na(df$`no soil`) & { val <- df$`no soil`; if (is.character(val)) val <- as.numeric(val); abs(val - 1) < 0.01 }, , drop = FALSE]
  } else df[FALSE, , drop = FALSE]
  veg_by_veg <- if ("Veg" %in% names(df)) {
    df[!is.na(df$Veg) & tolower(df$Veg) != "barren", , drop = FALSE]
  } else df[FALSE, , drop = FALSE]
  veg_rows <- unique(rbind(veg_by_nosoil, veg_by_veg))

  n_barren <- nrow(barren_rows)
  n_veg_pure <- nrow(veg_rows)
  cat(sprintf("Stage1 original sample sizes - Barren: %d, VegPure: %d\n", n_barren, n_veg_pure))

  # Use the same final_balancing_size computed for veg types if available
  stage1_target <- if (exists("final_balancing_size") && is.numeric(final_balancing_size) && final_balancing_size > 0) final_balancing_size else NA_integer_
  if (is.na(stage1_target)) {
    # Fallback: 2x min(barren, veg) to avoid tiny classes
    stage1_target <- 2 * min(c(max(1L, n_barren), max(1L, n_veg_pure)))
  }

  # Augment or downsample to match exactly stage1_target for both classes
  balance_stage1_class <- function(df_class, target_n) {
    if (nrow(df_class) == 0) return(df_class)
    if (nrow(df_class) >= target_n) {
      # Downsample to target size
      set.seed(123)
      if (nrow(df_class) > target_n) df_class <- df_class[sample.int(nrow(df_class), target_n), , drop = FALSE]
      return(df_class)
    } else {
      # Augment minority up to the target using previously defined augment_minority_class
      return(augment_minority_class(df_class, target_n))
    }
  }

  balanced_barren <- balance_stage1_class(barren_rows, stage1_target)
  balanced_veg <- balance_stage1_class(veg_rows, stage1_target)
  cat(sprintf("Stage1 balanced sample sizes - Barren: %d, VegPure: %d (target=%d)\n", nrow(balanced_barren), nrow(balanced_veg), stage1_target))

  # Combine into stage1 training dataframe
  stage1_df <- if (nrow(balanced_barren) > 0 && nrow(balanced_veg) > 0) {
    chunked_rbind(list(balanced_barren, balanced_veg), chunk_size = 50L)
  } else {
    # If either class empty, fallback to original df
    df
  }

  # Rebuild STAGE1_LIB using the balanced stage1_df
  STAGE1_LIB <- build_barren_veg_library(stage1_df, avail, min_samples = 5)
  if (is.null(STAGE1_LIB)) {
    warning("[Stage1 Balancing] Rebuilt STAGE1_LIB is NULL after balancing - using original STAGE1_LIB (if available)")
  } else {
    cat("[Stage1 Balancing] Rebuilt STAGE1_LIB from balanced stage-1 data\n")
  }
}

# Build basic library
# REVISED: We skip the per-DOY averaging step here.
# The 'lib' object will be constructed later from the variant centroids.
# We initialize an empty list for now to satisfy downstream checks until we rebuild it.
lib <- list()
# Placeholder for legacy code compatibility
for (v in vegs) {
  lib[[v]] <- list(n_samples = 0)
}

cat("Skipping per-DOY averaging (will build variants from raw traces)...\n")

timing_info$lib_construction_done <- Sys.time()
cat(sprintf(
  "Library construction (skipped) completed in %.1f seconds\n",
  as.numeric(difftime(timing_info$lib_construction_done, timing_info$moving_var_done, units = "secs"))
))

# Build raw index library for vegetation types
# Using raw spectral indices directly (no PCA)
cat("=== Building raw index library ===\n")

feature_cols <- avail
cat(sprintf("Raw index features: %s\n", paste(feature_cols, collapse=", ")))

# Compute column means and SDs for standardization
X_all <- as.matrix(lib_df[, feature_cols, drop = FALSE])
for (j in seq_len(ncol(X_all))) {
  col_vals <- X_all[, j]
  if (any(!is.finite(col_vals))) {
    mu_j <- mean(col_vals[is.finite(col_vals)], na.rm = TRUE)
    if (!is.finite(mu_j)) mu_j <- 0
    X_all[!is.finite(col_vals), j] <- mu_j
  }
}
mu_all <- colMeans(X_all)
feature_sds <- apply(X_all, 2, sd)
feature_sds[feature_sds <= 1e-10] <- 1.0

n_features <- length(feature_cols)
cat(sprintf("Raw index library: %d features from %d training samples\n", n_features, nrow(X_all)))

# Build raw_lib_templates: 365 x n_features matrix per vegetation type
cat("Computing raw index templates per vegetation type...\n")
raw_lib_templates <- list()

for (vname in vegs) {
  dveg <- lib_df[lib_df$Veg == vname, , drop = FALSE]
  if (nrow(dveg) < 10) {
    raw_lib_templates[[vname]] <- list(T = matrix(0, nrow = 365, ncol = n_features), n_samples = 0)
    next
  }
  
  # Extract raw features
  X_v <- as.matrix(dveg[, feature_cols, drop = FALSE])
  for (j in seq_len(ncol(X_v))) {
    col_vals <- X_v[, j]
    if (any(!is.finite(col_vals))) {
      X_v[!is.finite(col_vals), j] <- mu_all[j]
    }
  }
  
  # Standardize (center and scale)
  X_v_c <- sweep(X_v, 2, mu_all, "-")
  X_v_std <- sweep(X_v_c, 2, feature_sds, "/")
  
  # Group by DOY and compute medoid
  doy_vec <- lubridate::yday(dveg$date)
  T_medoid <- matrix(NA_real_, nrow = 365, ncol = n_features)
  
  for (d in 1:365) {
    rows_d <- which(doy_vec == d)
    if (length(rows_d) > 0) {
      sub <- X_v_std[rows_d, , drop = FALSE]
      if (nrow(sub) == 1) {
        T_medoid[d, ] <- sub[1, ]
      } else {
        center <- colMeans(sub)
        dists <- rowSums(sweep(sub, 2, center, "-")^2)
        T_medoid[d, ] <- sub[which.min(dists), ]
      }
    }
  }
  
  # Interpolate missing days
  for (j in 1:n_features) {
    y <- T_medoid[, j]
    if (any(is.na(y))) {
      x <- which(!is.na(y))
      if (length(x) >= 2) {
        T_medoid[is.na(y), j] <- approx(x, y[x], xout = which(is.na(y)), rule = 2)$y
      } else if (length(x) == 1) {
        T_medoid[, j] <- y[x[1]]
      } else {
        T_medoid[, j] <- 0
      }
    }
  }
  
  raw_lib_templates[[vname]] <- list(T = T_medoid, n_samples = nrow(dveg))
  cat(sprintf("  %s: T_medoid range [%.4f, %.4f], mean=%.4f\n", 
              vname, min(T_medoid), max(T_medoid), mean(T_medoid)))
}

cat("Raw index library templates computed.\n")

# Build Stage 1 Library for inference
cat("Building full-resolution Stage 1 Library for inference...\n")
COMPRESSED_STAGE1_LIB <- build_stage1_lib(STAGE1_LIB, grid_name = "full")

if (is.null(COMPRESSED_STAGE1_LIB)) {
  warning("[Stage1] Failed to build Stage 1 library")
} else {
  cat(sprintf("[Stage1] Full-resolution Stage 1 library: %d features per endmember\n",
              length(COMPRESSED_STAGE1_LIB$barren)))
  
}
if (isTRUE(ENABLE_UNCERTAINTY) && is.null(COMPRESSED_STAGE1_LIB)) {
  stop("ENABLE_UNCERTAINTY = TRUE requires a built COMPRESSED_STAGE1_LIB. Please build a compressed Stage 1 library before enabling uncertainty.")
}

timing_info$lib_templates_done <- Sys.time()
cat(sprintf(
  "Library templates computed in %.1f seconds\n",
  as.numeric(difftime(timing_info$lib_templates_done,
    timing_info$lib_construction_done, units = "secs"))
))

# ================== MESMA-SPECIFIC FUNCTIONS ==================

# Optional: seasonal phase alignment before temporal compression
safe_diff <- function(x) {
  if (length(x) <= 1) return(numeric(0))
  dx <- diff(x)
  dx[!is.finite(dx)] <- 0
  dx
}

detect_phase_markers <- function(dvi_proxy) {
  n <- length(dvi_proxy)
  if (n < 5) return(c(1L, round(n/4), round(n/2), round(3*n/4), n))
  dvi <- as.numeric(dvi_proxy)
  dvi[!is.finite(dvi)] <- 0
  d1 <- safe_diff(dvi)
  d2 <- safe_diff(d1)
  # Windows for green-up and senescence
  gup_win <- 1:max(2L, floor(n * 0.5))
  sen_win <- (max(1L, floor(n * 0.5))):n
  green_up <- if (length(gup_win) > 1 && length(d1) >= max(gup_win)) which.max(d1[gup_win]) else floor(n * 0.25)
  green_up <- max(2L, min(green_up, n - 3L))
  peak <- suppressWarnings(which.max(dvi)); if (!is.finite(peak) || length(peak) == 0) peak <- floor(n * 0.5)
  senescence <- if (length(sen_win) > 1 && length(d1) >= max(sen_win) && min(sen_win) <= length(d1)) {
    which.max(-d1[sen_win]) + (min(sen_win) - 1)
  } else floor(n * 0.75)
  senescence <- max(green_up + 1L, min(senescence, n - 1L))
  sort(unique(as.integer(c(1L, green_up, peak, senescence, n))))
}

align_to_phenological_phase <- function(pca_mat, reference_phase = REFERENCE_PHASE_MARKERS) {
  # pca_mat: 365 x k matrix of factor time series
  if (is.null(pca_mat) || !is.matrix(pca_mat) || nrow(pca_mat) < 5) return(pca_mat)
  k <- ncol(pca_mat)
  # DVI proxy from first up to 3 components
  use_cols <- seq_len(min(3L, k))
  dvi_proxy <- rowSums(abs(pca_mat[, use_cols, drop = FALSE]))
  dvi_proxy[!is.finite(dvi_proxy)] <- 0
  phase_markers <- detect_phase_markers(dvi_proxy)
  ref <- reference_phase
  if (is.null(ref) || length(ref) < 2) ref <- c(1, 90, 180, 270, 365)
  ref <- sort(unique(as.integer(pmax(1, pmin(nrow(pca_mat), ref)))))
  # Build piecewise-linear mapping from phase_markers -> ref
  aligned_mat <- matrix(NA_real_, nrow = nrow(pca_mat), ncol = k)
  for (j in seq_len(k)) {
    y <- pca_mat[, j]
    y[!is.finite(y)] <- 0
    # Sample values at phase markers
    pm <- phase_markers
    rm <- ref
    # Interpolate segment-wise by mapping original day indices to reference markers
    # Construct a monotone mapping function f: [1..n] -> [1..n]
    x_src <- pm
    x_dst <- rm[seq_len(min(length(rm), length(pm)))]
    # Ensure strictly increasing
    x_src <- unique(sort(x_src))
    x_dst <- unique(sort(x_dst))
    len_map <- min(length(x_src), length(x_dst))
    x_src <- x_src[seq_len(len_map)]
    x_dst <- x_dst[seq_len(len_map)]
    if (length(x_src) < 2) {
      aligned_mat[, j] <- y
    } else {
      # Map destination grid 1..n back to source via inverse approx
      fwd <- stats::approx(x = x_dst, y = x_src, xout = seq_len(nrow(pca_mat)), rule = 2)$y
      aligned_mat[, j] <- stats::approx(x = seq_len(nrow(pca_mat)), y = y, xout = fwd, rule = 2)$y
    }
  }
  aligned_mat[!is.finite(aligned_mat)] <- 0
  aligned_mat
}

# ===== Uncertainty: OLS + Block Bootstrap Helpers =====
estimate_block_size <- function(resid_series, min_block = 5, max_block = 60) {
  x <- as.numeric(resid_series)
  x <- x[is.finite(x)]
  if (length(x) < 10) return(max(min_block, 5L))
  acf_vals <- tryCatch(as.numeric(stats::acf(x, plot = FALSE, na.action = na.pass)$acf)[-1], error = function(e) stop(sprintf("estimate_block_size: acf computation failed: %s", e$message)))
  if (length(acf_vals) == 0) stop("estimate_block_size: acf returned no values; cannot estimate block size")
  thr <- 0.1
  k <- which(abs(acf_vals) < thr)[1]
  if (is.na(k)) k <- length(acf_vals)
  b <- max(min_block, min(max_block, 2L * k + 1L))
  as.integer(b)
}

## Bootstrap for Stage 1 (vegetated fraction) on compressed traces
stage1_block_bootstrap <- function(y_vec, compressed_stage1_lib, B = BOOTSTRAP_B, seed = 123) {
  if (is.null(y_vec) || is.null(compressed_stage1_lib)) return(NULL)
  set.seed(seed)

  B_vec <- compressed_stage1_lib$barren
  V_vec <- compressed_stage1_lib$vegetation
  len <- min(length(y_vec), length(B_vec), length(V_vec))
  y <- y_vec[1:len]; Bv <- B_vec[1:len]; Vv <- V_vec[1:len]

  # Compute alpha and predicted y (Stage 1 model)
  diff_BV <- Bv - Vv
  diff_XV <- y - Vv
  denom <- sum(diff_BV^2)
  if (denom < 1e-10) return(NULL)
  alpha <- safe_dot(diff_XV, diff_BV) / denom
  alpha <- max(0, min(1, alpha))
  veg_frac_hat <- 1 - alpha
  y_pred <- safe_mul_vec(alpha, Bv) + safe_mul_vec(1 - alpha, Vv)
  r <- y - y_pred

  bsz <- estimate_block_size(r)
  n <- length(y)

  veg_frac_boot <- rep(NA_real_, B)
  starts <- seq_len(n)
  for (b in seq_len(B)) {
    idx <- integer(0)
    while (length(idx) < n) {
      s <- sample(starts, 1)
      block <- ((s - 1) + seq_len(bsz) - 1) %% n + 1
      idx <- c(idx, block)
    }
    idx <- idx[seq_len(n)]
    r_b <- r[idx]
    yb <- y_pred + r_b
    veg_frac_b <- geometric_stage1_unmix(yb, compressed_stage1_lib$barren, compressed_stage1_lib$vegetation)$veg_frac
    if (length(veg_frac_b) > 1) veg_frac_b <- as.numeric(veg_frac_b[1])
    if (!is.finite(veg_frac_b)) veg_frac_b <- NA_real_
    veg_frac_boot[b] <- veg_frac_b
  }

  veg_frac_boot <- veg_frac_boot[is.finite(veg_frac_boot)]
  veg_frac_ci <- if (length(veg_frac_boot) >= 10) as.numeric(stats::quantile(veg_frac_boot, c(0.025, 0.975), na.rm = TRUE)) else c(NA_real_, NA_real_)
  if (length(veg_frac_boot) == 0) {
    warning(sprintf("stage1_block_bootstrap: no finite bootstrap samples returned (B=%d)" , B))
  } else {
    cat(sprintf("stage1_block_bootstrap: %d finite bootstrap samples (out of B=%d)\n", length(veg_frac_boot), B))
    # Debugging message removed: suppressed printing of veg_frac_boot summary to avoid excessive logs
  }
  list(veg_frac_boot = veg_frac_boot, veg_frac_ci = veg_frac_ci, veg_frac_hat = veg_frac_hat, n = length(veg_frac_boot), block_size = bsz)
}

## NOTE: Weight estimation now uses geometric projection (simplex-constrained least squares)
## Use `solve_weights_ols` (geometric projection) for all weight estimation.

ols_block_bootstrap <- function(y_vec, comp_templates, top_variants, chosen_ids, w_hat, B = BOOTSTRAP_B, lambda_star = 1e-6, seed = 123, variant_switch = FALSE, re_center = TRUE) {
  stop("ols_block_bootstrap is disabled: nested two-stage bootstrap is the only supported uncertainty method when ENABLE_UNCERTAINTY = TRUE")
  # Block bootstrap for OLS-weight uncertainty estimation.
  # By default (variant_switch = FALSE), the bootstrap keeps the fitted model fixed
  # (the same selected endmembers) and only perturbs observations via circular block bootstrap
  # of residuals. This ensures the bootstrap is assessing uncertainty for the same model.
  # If variant_switch = TRUE the function will reselect per-veg variants during each
  # iteration; if re_center = TRUE the bootstrap coefficient distribution will be
  # re-centered around the original point estimate to avoid the point estimate falling
  # outside the computed bootstrap CIs (a consequence of model changes between
  # iterations).
  # Returns list with coef_ci (data.frame), rmse_ci (numeric vector), variant_freq (data.frame)
  if (length(y_vec) < 5) return(NULL)
  set.seed(seed)
  # Rebuild E_best with columns per veg in names(top_variants) in fixed order
  vegs <- names(top_variants)
  cols <- list(); id_map <- list()
  for (v in vegs) {
    cid <- chosen_ids[[v]]
    vec <- NULL
    for (cand in comp_templates[[v]]) { if (!is.null(cand$id) && cand$id == cid) { vec <- cand$vec; break } }
    if (is.null(vec)) return(NULL)
    cols[[length(cols) + 1]] <- as.numeric(vec)
    id_map[[v]] <- sapply(comp_templates[[v]], function(c) c$id)
  }
  E_best <- do.call(cbind, cols)
  for (j in seq_len(ncol(E_best))) { nj <- sqrt(sum(E_best[, j]^2)); if (nj > 0) E_best[, j] <- E_best[, j] / nj }
  # Residuals from current fit to size block
  r <- as.numeric(y_vec - E_best %*% as.numeric(w_hat))
  bsz <- estimate_block_size(r)
  n <- length(y_vec)
  # Precompute blocks (circular to preserve ends)
  starts <- seq_len(n)
  res_mat <- matrix(NA_real_, nrow = n, ncol = B)
  coef_mat <- matrix(NA_real_, nrow = length(vegs), ncol = B, dimnames = list(vegs, NULL))
  # Keep a copy of raw bootstrap weights (before optional re-centering)
  w_boot_raw <- matrix(NA_real_, nrow = length(vegs), ncol = B, dimnames = list(vegs, NULL))
  rmse_vec <- rep(NA_real_, B)
  dom_counts <- lapply(vegs, function(.) list())
  names(dom_counts) <- vegs

  for (b in seq_len(B)) {
    # Circular block bootstrap of residuals
    idx <- integer(0)
    while (length(idx) < n) {
      s <- sample(starts, 1)
      block <- ((s - 1) + seq_len(bsz) - 1) %% n + 1
      idx <- c(idx, block)
    }
    idx <- idx[seq_len(n)]
    r_b <- r[idx]
    # Perturb y: y* = E_best w_hat + r_b
    yb <- as.numeric(E_best %*% as.numeric(w_hat)) + r_b
    if (isTRUE(variant_switch)) {
      # Variant reselection on perturbed yb: pick top per-veg by cosine, then solve with OLS
      chosen_b <- list(); cols_b <- list()
      for (v in vegs) {
        cand <- comp_templates[[v]]
        sims <- sapply(cand, function(x) cos_sim(yb, x$vec))
        i_best <- which.max(sims)
        chosen_b[[v]] <- cand[[i_best]]$id
        cols_b[[length(cols_b) + 1]] <- cand[[i_best]]$vec
      }
      E_b <- do.call(cbind, cols_b)
      for (j in seq_len(ncol(E_b))) { nj <- sqrt(sum(E_b[, j]^2)); if (nj > 0) E_b[, j] <- E_b[, j] / nj }
      w_b_res <- solve_weights_ols(E_b, yb, enforce_constraints = TRUE)
      w_b <- if (is.list(w_b_res) && !is.null(w_b_res$w)) as.numeric(w_b_res$w) else as.numeric(w_b_res)
    } else {
      # Fixed: use E_best (the original endmembers) for ALL bootstrap iterations
      E_b <- E_best
      chosen_b <- chosen_ids
      w_b_res <- solve_weights_ols(E_b, yb, enforce_constraints = TRUE)
      w_b <- if (is.list(w_b_res) && !is.null(w_b_res$w)) as.numeric(w_b_res$w) else as.numeric(w_b_res)
    }
    # If a solver returns a vector of different length than expected, put NA
    wb_num <- as.numeric(w_b)
    if (length(wb_num) == nrow(coef_mat)) {
      coef_mat[, b] <- wb_num
    } else {
      coef_mat[, b] <- rep(NA_real_, nrow(coef_mat))
    }
    w_boot_raw[, b] <- coef_mat[, b]
    pred_b <- as.numeric(E_b %*% w_b)
    rmse_vec[b] <- sqrt(mean((yb - pred_b)^2, na.rm = TRUE))
    # Track variant dominance
    if (isTRUE(variant_switch)) {
      for (v in vegs) {
        key <- chosen_b[[v]]
        prev <- dom_counts[[v]][[key]]
        if (is.null(prev)) prev <- 0
        dom_counts[[v]][[key]] <- prev + 1
      }
    }
  }
  # Compute CIs (2.5/97.5) for coefficients and RMSE
  ci_low <- ci_high <- rep(NA_real_, length(vegs)); names(ci_low) <- vegs; names(ci_high) <- vegs
  for (i in seq_along(vegs)) {
    x <- coef_mat[i, ]; x <- x[is.finite(x)]
    if (length(x) >= 1) { ci <- stats::quantile(x, c(0.025, 0.975), na.rm = TRUE); ci_low[i] <- ci[1]; ci_high[i] <- ci[2] } else {
      warning(sprintf("[BOOTSTRAP] Inflated uncertainty: insufficient finite bootstrap samples for %s (need %d, got %d)", vegs[i], 1, length(x)))
    }
  }
  if (isTRUE(DEBUG_UNCERTAINTY)) {
    for (i in seq_along(vegs)) {
      x <- as.numeric(coef_mat[i, ]); nfinite <- sum(is.finite(x));
      cat(sprintf("[BOOTSTRAP OLS DEBUG] Veg=%s: %d/%d finite samples (B=%d)\n", vegs[i], nfinite, B, B))
    }
    cat(sprintf("[BOOTSTRAP OLS DEBUG] Coef CI (first 10 rows):\n"))
    print(head(data.frame(Veg=vegs, coef_025=ci_low, coef_975=ci_high), 10))
  }
  rmse_ci <- if (sum(is.finite(rmse_vec)) >= 10) as.numeric(stats::quantile(rmse_vec, c(0.025, 0.975), na.rm = TRUE)) else c(NA_real_, NA_real_)
  coef_ci_df <- data.frame(Veg = vegs, coef_025 = ci_low, coef_975 = ci_high, stringsAsFactors = FALSE)
  if (isTRUE(DEBUG_UNCERTAINTY)) cat(sprintf("[BOOTSTRAP OLS DEBUG] Completed OLS bootstrap: coef_ci rows = %d\n", nrow(coef_ci_df)))
  # Variant frequencies
  var_rows <- NULL
  if (isTRUE(variant_switch)) {
    var_rows <- chunked_rbind(lapply(names(dom_counts), function(v) {
    freqs <- dom_counts[[v]]
    if (length(freqs) == 0) return(NULL)
    data.frame(Veg = v, Variant = names(freqs), Freq = as.numeric(unlist(freqs)), stringsAsFactors = FALSE)
  }))
    if (!is.null(var_rows) && nrow(var_rows) > 0) {
    var_rows <- var_rows %>% dplyr::group_by(.data$Veg, .data$Variant) %>% dplyr::summarize(N = sum(.data$Freq), .groups = "drop")
    var_rows$Percent <- 100 * var_rows$N / sum(var_rows$N)
  }
  }
  
  # If variant switching was used and re-centering requested, adjust the bootstrap distribution
  if (isTRUE(variant_switch) && isTRUE(re_center)) {
    # Shift each vegetation's bootstrap coefficients so that their mean equals w_hat
    for (i in seq_along(vegs)) {
      xi <- coef_mat[i, ]
      finite_mask <- is.finite(xi)
      if (sum(finite_mask) >= 3) {
        m <- mean(xi[finite_mask], na.rm = TRUE)
        shift <- as.numeric(w_hat[i] - m)
        coef_mat[i, finite_mask] <- coef_mat[i, finite_mask] + shift
      }
    }
    # Recompute the CIs from the re-centered distribution
    ci_low <- ci_high <- rep(NA_real_, length(vegs)); names(ci_low) <- vegs; names(ci_high) <- vegs
    for (i in seq_along(vegs)) {
      x <- coef_mat[i, ]; x <- x[is.finite(x)]
      if (length(x) >= 10) { ci <- stats::quantile(x, c(0.025, 0.975), na.rm = TRUE); ci_low[i] <- ci[1]; ci_high[i] <- ci[2] }
    }
  }
  # Provide raw and optionally centered bootstrap samples for potential future propagation methods (deprecated)
  w_median <- apply(coef_mat, 1, function(x) median(x, na.rm = TRUE))
  names(w_median) <- vegs
  list(coef_ci = coef_ci_df, rmse_ci = rmse_ci, variant_freq = var_rows, block_size = bsz,
       w_boot = coef_mat, # possibly re-centered bootstrap weights
       w_boot_raw = w_boot_raw, w_median = w_median)
}

## Nested two-stage bootstrap: propagate Stage 1 (barren vs vegetation fraction) uncertainty
# into the final vegetation coefficients by resampling residuals in compressed space and
# recomputing both Stage 1 and Stage 2 on each perturbed observation.
nested_two_stage_bootstrap <- function(y_vec, compressed_stage1_lib, comp_templates, top_variants, chosen_ids, w_hat, B = BOOTSTRAP_B, lambda_star = 1e-6, seed = 123, variant_switch = FALSE, re_center = TRUE) {
  if (is.null(y_vec) || length(y_vec) < 5) return(NULL)
  if (is.null(compressed_stage1_lib)) return(NULL)
  set.seed(seed)

  # Build E_best (fixed endmembers) from chosen_ids unless variant_switch requested
  vegs <- names(top_variants)
  cols <- list()
  for (v in vegs) {
    cid <- chosen_ids[[v]]
    vec <- NULL
    for (cand in comp_templates[[v]]) { if (!is.null(cand$id) && cand$id == cid) { vec <- cand$vec; break } }
    if (is.null(vec)) return(NULL)
    cols[[length(cols) + 1]] <- as.numeric(vec)
  }
  E_best <- do.call(cbind, cols)
  for (j in seq_len(ncol(E_best))) { nj <- sqrt(sum(E_best[, j]^2)); if (nj > 0) E_best[, j] <- E_best[, j] / nj }

  # Compute residuals from the stage 2 point fit
  r <- as.numeric(y_vec - E_best %*% as.numeric(w_hat))
  bsz <- estimate_block_size(r)
  n <- length(y_vec)

  coef_final_mat <- matrix(NA_real_, nrow = length(vegs), ncol = B, dimnames = list(vegs, NULL))
  veg_frac_vec <- rep(NA_real_, B)
  rmse_vec <- rep(NA_real_, B)
  # optionally collect stage2 raw bootstrap weights for compatibility with MC approach
  stage2_w_boot_raw <- matrix(NA_real_, nrow = length(vegs), ncol = B, dimnames = list(vegs, NULL))

  # Perform outer loop bootstrap: resample residuals and recompute both stages
  starts <- seq_len(n)
  for (b in seq_len(B)) {
    # Circular block bootstrap of residuals
    idx <- integer(0)
    while (length(idx) < n) {
      s <- sample(starts, 1)
      block <- ((s - 1) + seq_len(bsz) - 1) %% n + 1
      idx <- c(idx, block)
    }
    idx <- idx[seq_len(n)]
    r_b <- r[idx]

    # Perturb y
    yb <- as.numeric(E_best %*% as.numeric(w_hat)) + r_b

    # Stage 1: estimate vegetated fraction from compressed stage1 lib (using the perturbed yb)
    veg_frac_b <- geometric_stage1_unmix(yb, compressed_stage1_lib$barren, compressed_stage1_lib$vegetation)$veg_frac
    if (length(veg_frac_b) > 1) veg_frac_b <- as.numeric(veg_frac_b[1])
    if (!is.finite(veg_frac_b)) veg_frac_b <- NA_real_
    # Ensure veg_frac_b is scalar numeric
    if (length(veg_frac_b) > 1) veg_frac_b <- as.numeric(veg_frac_b[1])
    if (!is.finite(veg_frac_b)) veg_frac_b <- NA_real_
    if (!is.finite(veg_frac_b)) veg_frac_b <- NA_real_
    veg_frac_vec[b] <- veg_frac_b

    # Stage 2: fit weights on the same E_best (fixed endmembers) using yb
    w_b_res <- solve_weights_ols(E_best, yb, enforce_constraints = TRUE)
    w_b <- if (is.list(w_b_res) && !is.null(w_b_res$w)) as.numeric(w_b_res$w) else as.numeric(w_b_res)
    w_b_num <- as.numeric(w_b)
    if (length(w_b_num) == nrow(coef_final_mat)) {
      stage2_w_boot_raw[, b] <- w_b_num
    }
    if (length(w_b_num) == nrow(coef_final_mat)) {
      coef_final_mat[, b] <- as.numeric(safe_mul_vec(w_b_num, veg_frac_b))
    } else {
      coef_final_mat[, b] <- rep(NA_real_, nrow(coef_final_mat))
    }
    # Pred and RMSE for information
    pred_b <- as.numeric(E_best %*% w_b)
    rmse_vec[b] <- sqrt(mean((yb - pred_b)^2, na.rm = TRUE))
  }

  # Compute final CIs on the product (veg fraction * stage 2 coefficients)
  ci_low <- ci_high <- rep(NA_real_, length(vegs)); names(ci_low) <- vegs; names(ci_high) <- vegs
  for (i in seq_along(vegs)) {
    x <- coef_final_mat[i, ]; x <- x[is.finite(x)]
    if (length(x) >= 1) { ci <- stats::quantile(x, c(0.025, 0.975), na.rm = TRUE); ci_low[i] <- ci[1]; ci_high[i] <- ci[2] } else {
      warning(sprintf("[BOOTSTRAP NESTED] Insufficient finite bootstrap samples for %s (need %d, got %d)", vegs[i], 1, length(x)))
    }
  }
  if (isTRUE(DEBUG_UNCERTAINTY)) {
    for (i in seq_along(vegs)) {
      x_all <- as.numeric(coef_final_mat[i, ]); nfinite <- sum(is.finite(x_all));
      cat(sprintf("[BOOTSTRAP NESTED DEBUG] Veg=%s: %d/%d finite final samples (B=%d)\n", vegs[i], nfinite, B, B))
    }
    cat(sprintf("[BOOTSTRAP NESTED DEBUG] Coef CI computed for vegs: %s\n", paste(vegs, collapse=", ")))
    print(head(data.frame(Veg=vegs, coef_025=ci_low, coef_975=ci_high), 10))
  }

  veg_frac_ci <- if (sum(is.finite(veg_frac_vec)) >= 10) as.numeric(stats::quantile(veg_frac_vec, c(0.025, 0.975), na.rm = TRUE)) else c(NA_real_, NA_real_)

  coef_ci_df <- data.frame(Veg = vegs, coef_025 = ci_low, coef_975 = ci_high, stringsAsFactors = FALSE)
  if (isTRUE(DEBUG_UNCERTAINTY)) cat(sprintf("[BOOTSTRAP NESTED DEBUG] Completed nested bootstrap: coef_ci rows = %d\n", nrow(coef_ci_df)))
  w_median <- apply(coef_final_mat, 1, function(x) median(x, na.rm = TRUE))
  names(w_median) <- vegs
  list(coef_ci = coef_ci_df, veg_frac_ci = veg_frac_ci, rmse_ci = if (sum(is.finite(rmse_vec)) >= 10) as.numeric(stats::quantile(rmse_vec, c(0.025, 0.975), na.rm = TRUE)) else c(NA_real_, NA_real_), block_size = bsz, w_boot_raw = stage2_w_boot_raw, w_median = w_median)
}

# Perform full-resolution processing on all traces BEFORE variant construction
# Uses raw spectral indices directly (no PCA projection)
reduce_all_traces <- function(lib_df, veg_types, avail_idx, fixed_grid_size = TEMPORAL_BUDGET, 
                               enable_phase_alignment = FALSE, reference_phase = NULL,
                               enable_multiscale = FALSE, multiscale_windows = NULL) {
  cat("Performing full-resolution processing on all traces using raw indices...\n")
  
  # Fixed grid for all traces (ensures comparability)
  fixed_grid <- unique(round(seq(1, 365, length.out = fixed_grid_size)))
  
  reduced_data <- list()
  
  # Ensure lib_df is a data.frame (handle matrix or tibble inputs gracefully)
  if (is.matrix(lib_df)) lib_df <- as.data.frame(lib_df, stringsAsFactors = FALSE)
  # Resolve which column in lib_df holds vegetation labels
  veg_col_name <- NULL
  if ("Veg" %in% names(lib_df)) veg_col_name <- "Veg"
  else if ("veg" %in% names(lib_df)) veg_col_name <- "veg"
  else {
    # Try a case-insensitive search for a likely vegetation column
    lcnames <- tolower(names(lib_df))
    idx <- which(lcnames %in% c("veg", "vegetation", "vegetation_type", "type"))
    if (length(idx) > 0) veg_col_name <- names(lib_df)[idx[1]]
  }

  if (is.null(veg_col_name)) {
    cat("[WARN] reduce_all_traces: training data has no vegetation label column ('Veg' or similar). No per-veg reduction will be performed.\n")
    # cat(sprintf("[DEBUG] lib_df class: %s\n", paste(class(lib_df), collapse = ", ")))
    if (is.data.frame(lib_df) || is.matrix(lib_df)) {
      # cat(sprintf("[DEBUG] lib_df dimensions: %s x %s\n", if (!is.null(nrow(lib_df))) nrow(lib_df) else 0, if (!is.null(ncol(lib_df))) ncol(lib_df) else 0))
      # cat("[DEBUG] lib_df columns: ", paste(names(lib_df), collapse = ", "), "\n")
    }
    return(list())
  }

  for (veg in veg_types) {
    # subset by resolved vegetation column
    veg_data <- lib_df[tolower(as.character(lib_df[[veg_col_name]])) == tolower(as.character(veg)), , drop = FALSE]
    
    if (!"location_id" %in% names(veg_data) || !"year" %in% names(veg_data)) {
      cat(sprintf("  [%s] Missing location_id or year, skipping\n", veg))
      next
    }
    
    traces <- unique(veg_data[, c("location_id", "year")])
    cat(sprintf("  [%s] Reducing %d traces to %d time points...\n", veg, nrow(traces), length(fixed_grid)))
    
    feature_list <- list()
    Z_list <- list()
    trace_info <- list()
    
    for (i in seq_len(nrow(traces))) {
      loc <- traces$location_id[i]
      yr <- traces$year[i]
      
      dly_year <- veg_data[veg_data$location_id == loc & veg_data$year == yr, , drop = FALSE]
      n_unique_doys <- length(unique(dly_year$doy))
      if (n_unique_doys < 5) next
      
      # Build a full 365 x K raw-index matrix for this trace
      # For each index in avail_idx we choose the medoid per DOY and interpolate missing days
      idxs <- avail_idx
      K <- length(idxs)
      raw_mat <- matrix(NA_real_, nrow = 365, ncol = K)
      colnames(raw_mat) <- idxs

      for (j in seq_along(idxs)) {
        idn <- idxs[j]
        if (!idn %in% names(dly_year)) next
        vals_by_doy <- tryCatch(
          {
            tapply(seq_along(dly_year[[idn]]), dly_year$doy, function(indices) {
              vals <- dly_year[[idn]][indices]
              vals <- vals[is.finite(vals)]
              if (length(vals) == 0) return(NA_real_)
              if (length(vals) == 1) return(as.numeric(vals[1]))
              m <- mean(vals)
              as.numeric(vals[which.min(abs(vals - m))])
            })
          }, error = function(e) NULL)

        if (!is.null(vals_by_doy) && length(vals_by_doy) > 0) {
          doy_values <- as.integer(names(vals_by_doy))
          valid_doy <- doy_values >= 1 & doy_values <= 365
          raw_mat[doy_values[valid_doy], j] <- vals_by_doy[valid_doy]
        }

        # Interpolate missing days per index column
        colv <- raw_mat[, j]
        finite_idx <- which(is.finite(colv))
        if (length(finite_idx) == 0) {
          raw_mat[, j] <- 0
        } else if (length(finite_idx) == 1) {
          raw_mat[, j] <- colv[finite_idx]
        } else {
          x <- c(finite_idx - 365, finite_idx, finite_idx + 365)
          y <- rep(colv[finite_idx], 3)
          interp <- tryCatch(stats::approx(x = x, y = y, xout = seq_len(365), rule = 2)$y, error = function(e) rep(colv[finite_idx[1]], 365))
          interp <- as.numeric(interp)
          interp[!is.finite(interp)] <- median(colv[finite_idx], na.rm = TRUE)
          raw_mat[, j] <- interp
        }
      }

      # Full-resolution processing performed ON RAW INDEX MATRIX (features = 365 × indices)
      feat <- as.numeric(raw_mat)
      
      if (any(!is.finite(feat))) next
      
      # Ensure feature vector has consistent names: index_t1..index_tN
      if (is.numeric(feat) && length(feat) == length(idxs) * length(fixed_grid)) {
        grid_count <- length(fixed_grid)
        nm <- unlist(lapply(idxs, function(x) paste0(x, "_t", seq_len(grid_count))))
        names(feat) <- nm
      }
      feature_list[[length(feature_list) + 1]] <- feat
      Z_list[[length(Z_list) + 1]] <- raw_mat  # Store raw_mat for discriminative templates
      trace_info[[length(trace_info) + 1]] <- list(
        location_id = loc,
        year = yr,
        trace_index = i
      )
    }
    
    if (length(feature_list) > 0) {
      reduced_data[[veg]] <- list(
        features = do.call(rbind, feature_list),  # Matrix: n_traces x n_features
        Z_matrices = Z_list,                       # List of full 365 x K RAW INDEX matrices
        trace_info = trace_info,                   # Metadata for each trace
        n_samples = nrow(veg_data)                 # Total samples for this veg
      )
      cat(sprintf("  [%s] Reduced %d traces successfully\n", veg, length(feature_list)))
      # Sanity checks: expected feature width = TEMPORAL_BUDGET * n_indices
      feat_mat <- reduced_data[[veg]]$features
      grid_count <- length(fixed_grid)
      expected_cols <- grid_count * length(idxs)
      if (!is.null(dim(feat_mat)) && ncol(feat_mat) != expected_cols) {
        warning(sprintf("reduce_all_traces: feature column count mismatch for '%s' (got %d, expected %d).", veg, ncol(feat_mat), expected_cols))
      }
      keys <- sapply(reduced_data[[veg]]$trace_info, function(x) paste0(as.character(x$location_id), "__", as.character(x$year)))
      if (any(duplicated(keys))) {
        warning(sprintf("reduce_all_traces: duplicate location-year traces found for '%s' (%d duplicates) — expected exactly one row per loc/year.", veg, sum(duplicated(keys))))
      }
    }
  }
  
  cat(sprintf("Full-resolution processing complete: %d vegetation types processed\n", length(reduced_data)))
  return(reduced_data)
}

analyze_library_similarity <- function(mesma_lib, compressed_templates_accessor, grid_type = "full") {
  cat("\n=== INTER-CLASS VARIANT SIMILARITY ANALYSIS ===\n")
  
  # Collect all variant vectors
  all_variants <- list()
  variant_names <- c()
  variant_vegs <- c()
  
  for (veg in names(mesma_lib)) {
    for (variant in mesma_lib[[veg]]) {
      vid <- variant$variant_id
      vec <- compressed_templates_accessor[[veg]][[vid]][["full"]]
      if (!is.null(vec)) {
        all_variants[[length(all_variants) + 1]] <- list(vec = vec, variant = variant, veg = veg)
        variant_names <- c(variant_names, vid)
        variant_vegs <- c(variant_vegs, veg)
      }
    }
  }
  
  if (length(all_variants) < 2) {
    cat("Not enough variants for similarity analysis.\n")
    return(NULL)
  }
  
  # Report small variants
  cat("\n=== SMALL VARIANTS (potential outliers) ===\n")
  for (i in seq_along(all_variants)) {
    var <- all_variants[[i]]
    if (var$variant$n_samples < 10) {
      cat(sprintf("Variant %s (%s): %d samples\n", var$variant$variant_id, var$veg, var$variant$n_samples))
      if (!is.null(var$variant$sample_ids)) {
        cat(sprintf("  Samples: %s\n", paste(var$variant$sample_ids, collapse = ", ")))
      }
    }
  }
  cat("=======================================\n\n")
  
  # Compute pairwise similarity
  n_vars <- length(all_variants)
  
  high_sim_pairs <- data.frame(
    Variant1 = character(),
    Veg1 = character(),
    Variant2 = character(),
    Veg2 = character(),
    Similarity = numeric(),
    threshold_tier = character(),
    euclidean_dist = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:(n_vars - 1)) {
    for (j in (i + 1):n_vars) {
      # Cosine Similarity
      sim <- cos_sim(all_variants[[i]]$vec, all_variants[[j]]$vec)
      
      # DEBUG: Check for zero similarity
      if (sim == 0) {
        cat(sprintf("DEBUG: Zero similarity between %s (%s) and %s (%s)\n", 
                    variant_names[i], variant_vegs[i], variant_names[j], variant_vegs[j]))
        
        # Check if vectors are degenerate
        vec1 <- all_variants[[i]]$vec
        vec2 <- all_variants[[j]]$vec
        vec1_norm <- sqrt(sum(vec1^2))
        vec2_norm <- sqrt(sum(vec2^2))
        vec1_zero <- all(vec1 == 0)
        vec2_zero <- all(vec2 == 0)
        
        cat(sprintf("  Vec1 (%s): norm=%.6f, all_zero=%s, range=[%.6f, %.6f]\n", 
                    variant_names[i], vec1_norm, vec1_zero, min(vec1), max(vec1)))
        cat(sprintf("  Vec2 (%s): norm=%.6f, all_zero=%s, range=[%.6f, %.6f]\n", 
                    variant_names[j], vec2_norm, vec2_zero, min(vec2), max(vec2)))
      }
      
      # Euclidean Distance (on normalized vectors if possible, but here on compressed features)
      dist <- sqrt(sum((all_variants[[i]]$vec - all_variants[[j]]$vec)^2))
      
      # Check for high similarity between DIFFERENT vegetation types
      if (variant_vegs[i] != variant_vegs[j]) {
        tier <- NA_character_
        if (sim > 0.98) tier <- "CRITICAL (>0.98)"
        else if (sim > 0.95) tier <- "HIGH (>0.95)"
        else if (sim > 0.90) tier <- "WARN (>0.90)"
        
        if (!is.na(tier)) {
           high_sim_pairs[nrow(high_sim_pairs) + 1, ] <- list(
             variant_names[i], variant_vegs[i],
             variant_names[j], variant_vegs[j],
             sim, tier, dist
           )
        }
      }
    }
  }
  
  # Print high similarity pairs
  if (nrow(high_sim_pairs) > 0) {
    high_sim_pairs <- high_sim_pairs[order(-high_sim_pairs$Similarity), ]
    cat(sprintf("Found %d pairs of variants from DIFFERENT vegetation types with similarity > 0.90:\n", nrow(high_sim_pairs)))
    print(head(high_sim_pairs, 20))
    if (nrow(high_sim_pairs) > 20) cat(sprintf("... and %d more.\n", nrow(high_sim_pairs) - 20))
    
    # Store for later use if needed
    assign("INSEPARABLE_VARIANT_INFO", list(pairs = high_sim_pairs), envir = globalenv())
    
    # Create a set of inseparable variants to flag
    inseparable_vars <- list()
    for (i in 1:nrow(high_sim_pairs)) {
       # Flag anything > 0.95 as potentially problematic
       if (high_sim_pairs$Similarity[i] > 0.95) {
         v1 <- high_sim_pairs$Veg1[i]
         v2 <- high_sim_pairs$Veg2[i]
         vid1 <- high_sim_pairs$Variant1[i]
         vid2 <- high_sim_pairs$Variant2[i]
         
         if (is.null(inseparable_vars[[v1]])) inseparable_vars[[v1]] <- c()
         inseparable_vars[[v1]] <- unique(c(inseparable_vars[[v1]], vid1))
         
         if (is.null(inseparable_vars[[v2]])) inseparable_vars[[v2]] <- c()
         inseparable_vars[[v2]] <- unique(c(inseparable_vars[[v2]], vid2))
       }
    }
    assign("INSEPARABLE_VARIANTS", inseparable_vars, envir = globalenv())
    
  } else {
    cat("No high similarity (> 0.90) pairs found between different vegetation types.\n")
    assign("INSEPARABLE_VARIANTS", NULL, envir = globalenv())
  }
  
  cat("==============================================\n\n")
}

#' Simplified reduce_all_traces (NO compression - full year resolution)
reduce_all_traces_simple <- function(lib_df, veg_types, avail_idx) {
  
  reduced_data <- list()
  
  for (veg in veg_types) {
    veg_data <- lib_df[tolower(lib_df$Veg) == tolower(veg), ]
    if (nrow(veg_data) == 0) next
    
    traces <- unique(veg_data[, c("location_id", "year")])
    
    feature_list <- list()
    Z_list <- list()
    trace_info <- list()
    
    for (i in seq_len(nrow(traces))) {
      loc <- traces$location_id[i]
      yr <- traces$year[i]
      
      dly_year <- veg_data[veg_data$location_id == loc & veg_data$year == yr, ]
      n_unique_doys <- length(unique(dly_year$doy))
      if (n_unique_doys < 5) next
      
      # Build raw matrix
      raw_mat <- build_raw_365_matrix(dly_year, avail_idx)
      if (is.null(raw_mat)) next
      
      # Require at least 5 unique DOYs
      n_unique_doys <- length(unique(dly_year$doy))
      if (n_unique_doys < 5) next
      
      # Use full raw matrix flattened as features (no compression)
      full_features <- as.numeric(raw_mat)  # 365 * K vector
      if (any(!is.finite(full_features))) next
      
      feature_list[[length(feature_list) + 1]] <- full_features
      Z_list[[length(Z_list) + 1]] <- raw_mat
      trace_info[[length(trace_info) + 1]] <- list(
        location_id = loc, year = yr, trace_index = i
      )
    }
    
    if (length(feature_list) > 0) {
      reduced_data[[veg]] <- list(
        features = do.call(rbind, feature_list),
        Z_matrices = Z_list,
        trace_info = trace_info,
        n_samples = nrow(veg_data)
      )
      cat(sprintf("  [%s] Processed %d traces\n", veg, length(feature_list)))
    }
  }
  
  reduced_data
}

#' Modified build_mesma_variants (uses PCA-LDA weights for clustering)
build_mesma_variants_weighted <- function(reduced_data, raw_lib_templates, 
                                           pca_lda_weights = NULL,
                                           n_variants = 10,
                                           min_cluster_size = 10) {
  mesma_lib <- list()
  
  for (veg in names(reduced_data)) {
    veg_info <- reduced_data[[veg]]
    X_feat <- veg_info$features
    Z_list <- veg_info$Z_matrices
    
    if (nrow(X_feat) < min_cluster_size) {
      # Single variant fallback
      mesma_lib[[veg]] <- list(list(
        raw_mat = raw_lib_templates[[veg]]$T,
        variant_id = paste0(veg, "_single"),
        n_samples = veg_info$n_samples,
        sample_ids = sapply(veg_info$trace_info, function(x) paste0(as.character(x$location_id), "__", as.character(x$year)))
      ))
      next
    }
    
    # Per-index standardization across traces and time to equalize index importance
    if (!is.null(Z_list) && length(Z_list) > 0) {
      n_timepoints <- nrow(Z_list[[1]])
      K <- ncol(Z_list[[1]])
    } else {
      n_timepoints <- 365
      K <- as.integer(ncol(X_feat) / n_timepoints)
    }
    X_std <- matrix(NA_real_, nrow = nrow(X_feat), ncol = ncol(X_feat))
    for (k_idx in seq_len(K)) {
      col_start <- (k_idx - 1) * n_timepoints + 1
      col_end <- k_idx * n_timepoints
      cols <- seq(col_start, col_end)
      if (max(cols) > ncol(X_feat)) cols <- cols[cols <= ncol(X_feat)]
      index_data <- X_feat[, cols, drop = FALSE]
      global_mean <- mean(index_data, na.rm = TRUE)
      global_sd <- sd(as.vector(index_data), na.rm = TRUE)
      if (!is.finite(global_sd) || global_sd < 1e-10) global_sd <- 1.0
      X_std[, cols] <- (index_data - global_mean) / global_sd
    }
    X_std[!is.finite(X_std)] <- 0
    # Apply PCA-LDA weights for clustering (if available)
    # NOTE: We apply weights after performing per-index standardization
    X_weighted <- X_std
    if (!is.null(pca_lda_weights) && length(pca_lda_weights) == ncol(X_feat)) {
      # Weight features by discriminative importance
      X_weighted <- sweep(X_std, 2, sqrt(pca_lda_weights), "*")
    }
    
    # Remove constant columns before scaling to avoid scale() error
    col_vars <- apply(X_weighted, 2, var, na.rm = TRUE)
    keep_cols <- col_vars > 1e-10
    if (sum(keep_cols) < 10) {
      warning(sprintf("[%s] Too few variable features, using single variant", veg))
      mesma_lib[[veg]] <- list(list(
        raw_mat = raw_lib_templates[[veg]]$T,
        variant_id = paste0(veg, "_single"),
        n_samples = veg_info$n_samples
      ))
      next
    }
    X_weighted <- X_weighted[, keep_cols, drop = FALSE]
    
    # Standardize for clustering (center and scale)
    X_std <- scale(X_weighted, center = TRUE, scale = TRUE)
    X_std[!is.finite(X_std)] <- 0
    
    # Optional: PCA for dimensionality reduction before clustering
    # This addresses multicollinearity in clustering space
    if (ncol(X_std) > 50) {
      pca_clust <- prcomp(X_std, center = FALSE, scale. = FALSE)
      cum_var <- cumsum(pca_clust$sdev^2) / sum(pca_clust$sdev^2)
      n_pcs <- min(50, which(cum_var >= 0.95)[1])
      if (is.na(n_pcs)) n_pcs <- min(50, ncol(pca_clust$x))
      X_clust <- pca_clust$x[, 1:n_pcs, drop = FALSE]
      cat(sprintf("  [%s] Clustering in %d-PC space\n", veg, n_pcs))
    } else {
      X_clust <- X_std
    }
    
    # K-means clustering
    k <- min(n_variants, floor(nrow(X_clust) / 3))
    k <- max(k, 1)
    k <- min(k, 7)  # Cap at 7
    
    km <- kmeans(X_clust, centers = k, nstart = 25, iter.max = 100)
    
    # Build variants from cluster medoids
    variants <- list()
    for (clust in seq_len(k)) {
      members <- which(km$cluster == clust)
      if (length(members) == 0) next
      
      # Find medoid in ORIGINAL (unweighted) feature space
      X_members <- X_feat[members, , drop = FALSE]
      center <- colMeans(X_members)
      dists <- rowSums(sweep(X_members, 2, center, "-")^2)
      medoid_idx <- members[which.min(dists)]
      
      raw_mat <- Z_list[[medoid_idx]]
      
      # Check if raw_mat is degenerate
      raw_vec <- as.numeric(raw_mat)
      raw_vec[is.na(raw_vec)] <- 0
      raw_norm <- sqrt(sum(raw_vec^2))
      if (raw_norm < 1e-10) {
        warning(sprintf("Skipping zero norm variant %s_v%d (norm=%.2e)", veg, clust, raw_norm))
        next
      }
      raw_sd <- sd(raw_vec, na.rm = TRUE)
      raw_range <- diff(range(raw_vec, na.rm = TRUE))
      raw_unique <- length(unique(raw_vec[is.finite(raw_vec)]))
      
      if (!is.finite(raw_sd) || raw_sd < 1e-10) {
        warning(sprintf("Skipping degenerate variant %s_v%d (SD=%.2e)", veg, clust, raw_sd))
        next
      }
      
      if (!is.finite(raw_range) || raw_range < 1e-10) {
        warning(sprintf("Skipping degenerate variant %s_v%d (range=%.2e)", veg, clust, raw_range))
        next
      }
      
      if (raw_unique < 3) {
        warning(sprintf("Skipping degenerate variant %s_v%d (only %d unique values)", veg, clust, raw_unique))
        next
      }
      
      variants[[length(variants) + 1]] <- list(
        raw_mat = raw_mat,
        variant_id = paste0(veg, "_v", clust),
        n_samples = length(members),
        cluster_center = km$centers[clust, ],
        sample_ids = sapply(veg_info$trace_info[members], function(x) paste0(as.character(x$location_id), "__", as.character(x$year)))
      )
    }

    # Post-processing: Remove outlier variants based on similarity to other variants
    # Addresses issue where one prototype might be "off zero similar" (outlier)
    if (length(variants) >= 3) {
      var_vecs <- do.call(rbind, lapply(variants, function(v) {
        vec <- as.numeric(v$raw_mat)
        vec[is.na(vec)] <- 0
        vec
      }))
      n_v <- nrow(var_vecs)
      
      # Compute pairwise cosine similarity matrix
      sim_mat <- matrix(1, n_v, n_v)
      for (i in 1:(n_v-1)) {
        for (j in (i+1):n_v) {
          # Cosine similarity
          dot_prod <- sum(var_vecs[i,] * var_vecs[j,])
          norm_i <- sqrt(sum(var_vecs[i,]^2))
          norm_j <- sqrt(sum(var_vecs[j,]^2))
          
          if (norm_i > 0 && norm_j > 0) {
            s <- dot_prod / (norm_i * norm_j)
          } else {
            s <- 0
          }
          sim_mat[i,j] <- s
          sim_mat[j,i] <- s
        }
      }
      
      # Average similarity to others (excluding self)
      avg_sim <- rowSums(sim_mat - diag(n_v)) / (n_v - 1)
      
      # Identify outliers: avg_sim < 0.5 is very low for within-class variants
      outlier_mask <- rep(FALSE, length(avg_sim))
      
      if (any(outlier_mask)) {
         cat(sprintf("  [%s] Removed %d outlier variants (avg pairwise sim < 0.5)\n", veg, sum(outlier_mask)))
         variants <- variants[!outlier_mask]
      }
    }
    
    mesma_lib[[veg]] <- variants
    cat(sprintf("  [%s] Created %d variants\n", veg, length(variants)))
  }
  
  mesma_lib
}

# Build MESMA library with endmember variants
# Uses raw index matrices (365 x n_features) for templates
build_mesma_variants <- function(reduced_data, raw_lib_templates, min_cluster_size = 10) {
  mesma_lib <- list()

  # Identify barren prototype from reduced data if available
  barren_proto_raw <- NULL
  if ("barren" %in% names(reduced_data)) {
    b_feat <- reduced_data[["barren"]]$features
    if (!is.null(b_feat) && nrow(b_feat) > 0) {
      barren_proto_raw <- apply(b_feat, 2, median, na.rm = TRUE)
      cat(sprintf("[Soil Correction] Identified barren prototype from %d traces\n", nrow(b_feat)))
    }
  }

  for (veg in names(reduced_data)) {
    veg_info <- reduced_data[[veg]]
    
    # Extract pre-reduced features and full raw index matrices
    X_feat <- veg_info$features        # n_traces x n_features (time-reduced)
    Z_list <- veg_info$Z_matrices      # List of 365 x n_indices raw matrices
    n_samples <- veg_info$n_samples
    
    if (nrow(X_feat) < min_cluster_size) {
      # Fallback to single variant using raw template
      mesma_lib[[veg]] <- list(
        list(
          raw_mat = raw_lib_templates[[veg]]$T,
          variant_id = paste0(veg, "_single"),
          n_samples = n_samples
        )
      )
      cat(sprintf("  [%s] Insufficient traces (%d), using single variant\n", veg, nrow(X_feat)))
      next
    }
    
    cat(sprintf("  [%s] Building variants from %d reduced traces...\n", veg, nrow(X_feat)))
    
    # Subsample if too many traces
    if (nrow(X_feat) > MAX_PROJECTIONS_PER_VEG) {
      set.seed(123)
      idx <- sample(nrow(X_feat), MAX_PROJECTIONS_PER_VEG)
      X_feat <- X_feat[idx, , drop = FALSE]
      Z_list <- Z_list[idx]
      cat(sprintf("  [%s] Subsampled to %d traces\n", veg, nrow(X_feat)))
    }
    # Whitening on the reduced features (if Xw precomputed in the reduced_data, use that)
    if (!is.null(veg_info$Xw)) {
      X_w <- as.matrix(veg_info$Xw)
      whitened <- list(Xw = X_w, W = if (!is.null(veg_info$W)) veg_info$W else diag(ncol(X_w)), mu = if (!is.null(veg_info$mu)) veg_info$mu else rep(0, ncol(X_w)))
    } else {
      # ========== HYBRID APPROACH: Standardization + Optional PCA instead of Whitening ==========
      cat(sprintf("\n[CLUSTERING PREP in build_mesma_variants] %s: %d traces × %d features\n", 
                  veg, nrow(X_feat), ncol(X_feat)))
      
      # 1. Standardize per index (z-scores) across all traces and all time points
      # This preserves temporal dependence while equalizing index importance.
      # Use Z_list[[]] to determine number of timepoints per trace (n_timepoints) and K indices
      if (!is.null(Z_list) && length(Z_list) > 0) {
        n_timepoints <- nrow(Z_list[[1]])
        K <- ncol(Z_list[[1]])
      } else {
        # Fallback: We assume 365 DOYs—this should rarely be needed
        n_timepoints <- 365
        K <- as.integer(ncol(X_feat) / n_timepoints)
      }
      # Pre-allocate standardized matrix
      X_std <- matrix(NA_real_, nrow = nrow(X_feat), ncol = ncol(X_feat))
      for (k_idx in seq_len(K)) {
        col_start <- (k_idx - 1) * n_timepoints + 1
        col_end <- k_idx * n_timepoints
        cols <- seq(col_start, col_end)
        if (max(cols) > ncol(X_feat)) {
          cols <- cols[cols <= ncol(X_feat)]
        }
        index_data <- X_feat[, cols, drop = FALSE]
        global_mean <- mean(index_data, na.rm = TRUE)
        global_sd <- sd(as.vector(index_data), na.rm = TRUE)
        if (!is.finite(global_sd) || global_sd < 1e-10) global_sd <- 1.0
        X_std[, cols] <- (index_data - global_mean) / global_sd
      }
      X_std[!is.finite(X_std)] <- 0
      
      # 2. Dimensionality reduction ONLY if features >> samples
      if (ncol(X_std) > 50 && nrow(X_std) < ncol(X_std) * 3) {
        # Use PCA to avoid curse of dimensionality in clustering
        pca_result <- prcomp(X_std, center = FALSE, scale. = FALSE)
        var_explained <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
        n_pcs <- min(50, which(var_explained > 0.95)[1])
        if (is.na(n_pcs)) n_pcs <- min(50, ncol(pca_result$x))
        X_w <- pca_result$x[, 1:n_pcs, drop = FALSE]
        cat(sprintf("  Reduced to %d PCs (%.1f%% variance) for clustering\n", 
                    n_pcs, 100 * var_explained[n_pcs]))
      } else {
        X_w <- X_std
        cat(sprintf("  Using all %d standardized features for clustering\n", ncol(X_std)))
      }
      
      whitened <- list(Xw = X_w, W = diag(ncol(X_feat)), mu = X_means)
      # ========== END HYBRID APPROACH ==========
    }
    
    # Soil subtraction is only performed here when data are NOT pre-whitened and we have a raw barren prototype
    if (is.null(veg_info$Xw) && !is.null(barren_proto_raw) && veg != "barren") {
      if (length(barren_proto_raw) == ncol(X_w)) {
        b_centered <- barren_proto_raw - whitened$mu
        b_proj <- as.numeric(b_centered %*% whitened$W)
        b_norm2 <- sum(b_proj^2)
        if (b_norm2 > 1e-9) {
          alphas <- (X_w %*% b_proj) / b_norm2
          # Orthogonalization: remove component parallel to soil
          X_w <- X_w - (alphas %*% t(b_proj))
          cat(sprintf("  [%s] Applied soil subtraction (orthogonalization) in whitened space\n", veg))
        }
      }
    }
    
    # Special case for barren: use single prototype instead of clustering
    if (veg == "barren") {
      clust_members <- seq_len(nrow(X_w))
      median_center <- apply(X_w, 2, median, na.rm = TRUE)
      dists <- rowSums(sweep(X_w[clust_members, , drop=FALSE], 2, median_center, "-")^2)
      best_idx_local <- which.min(dists)
      best_idx_global <- clust_members[best_idx_local]
      medoid_Z <- Z_list[[best_idx_global]]
      
      # Check if medoid_Z is degenerate
      raw_vec <- as.numeric(medoid_Z)
      raw_sd <- sd(raw_vec, na.rm = TRUE)
      raw_range <- diff(range(raw_vec, na.rm = TRUE))
      raw_unique <- length(unique(raw_vec[is.finite(raw_vec)]))
      
      if (!is.finite(raw_sd) || raw_sd < 1e-10 ||
          !is.finite(raw_range) || raw_range < 1e-10 ||
          raw_unique < 3) {
        warning(sprintf("Barren prototype is degenerate, using fallback"))
        medoid_Z <- raw_lib_templates[["barren"]]$T
      }
      
      variants <- list(list(
        raw_mat = medoid_Z,
        variant_id = paste0(veg, "_single"),
        n_samples = length(clust_members),
        cluster_center = median_center
      ))
      mesma_lib[[veg]] <- variants
      next
    }
    
    # Endmember Construction: Cluster the time-reduced features
    
    # User requested exactly 7 variants per vegtype
    best_k <- 7
    
    # Ensure we don't exceed number of samples (though min_cluster_size check above should prevent this)
    if (best_k > nrow(X_w)) {
       best_k <- nrow(X_w)
       cat(sprintf("  [%s] Warning: Requested 10 variants but only %d samples available. Using %d variants.\n", veg, nrow(X_w), best_k))
    }
    
    best_sil <- NA_real_ # Silhouette not calculated
    
    km_final <- kmeans(X_w, centers = best_k, nstart = 25, iter.max = 100)
    
    # Clustering diagnostics
    cat(sprintf("\n[CLUSTERING DEBUG] %s: best_k=%d (silhouette=%.4f)\n", 
                veg, best_k, best_sil))
    
    cluster_sizes <- table(km_final$cluster)
    cat(sprintf("  Cluster sizes: %s\n", 
                paste(sprintf("C%d=%d", 1:best_k, cluster_sizes), collapse=", ")))
    
    bss_tss <- km_final$betweenss / km_final$totss
    cat(sprintf("  Between-cluster variance: %.1f%%\n", 100 * bss_tss))
    
    if (bss_tss < 0.3) {
      warning(sprintf("%s: Clusters are poorly separated (%.1f%% variance explained)", 
                      veg, 100 * bss_tss))
    }
    
    # Create variants
    variants <- list()
    for (clust in seq_len(best_k)) {
      clust_members <- which(km_final$cluster == clust)
      if (length(clust_members) < 1) next
      
      # Find medoid: sample closest to cluster centroid
      dists <- rowSums(sweep(X_w[clust_members, , drop=FALSE], 2, km_final$centers[clust, ], "-")^2)
      best_idx_local <- which.min(dists)
      best_idx_global <- clust_members[best_idx_local]
      
      # The variant is defined by the raw index matrix of the medoid trace
      medoid_Z <- Z_list[[best_idx_global]]
      
      # Check if medoid_Z is degenerate
      raw_vec <- as.numeric(medoid_Z)
      raw_sd <- sd(raw_vec, na.rm = TRUE)
      raw_range <- diff(range(raw_vec, na.rm = TRUE))
      raw_unique <- length(unique(raw_vec[is.finite(raw_vec)]))
      
      if (!is.finite(raw_sd) || raw_sd < 1e-10) {
        warning(sprintf("Skipping degenerate variant %s_v%d (SD=%.2e)", veg, clust, raw_sd))
        next
      }
      
      if (!is.finite(raw_range) || raw_range < 1e-10) {
        warning(sprintf("Skipping degenerate variant %s_v%d (range=%.2e)", veg, clust, raw_range))
        next
      }
      
      if (raw_unique < 3) {
        warning(sprintf("Skipping degenerate variant %s_v%d (only %d unique values)", veg, clust, raw_unique))
        next
      }
      
      variants[[length(variants) + 1]] <- list(
        raw_mat = medoid_Z,
        variant_id = paste0(veg, "_v", clust),
        n_samples = length(clust_members),
        cluster_center = km_final$centers[clust, ]
      )
    }
    
    mesma_lib[[veg]] <- variants
  }
  
  return(mesma_lib)
}

  # Geometric solver using simplex projection for constrained least squares unmixing
  # Uses pure geometric projection with sum-to-one and non-negativity constraints
  solve_weights_ols <- function(E, y, enforce_constraints = TRUE) {
    # E: n x p matrix (endmembers), y: n-vector (spectrum)
    # Returns: list(w = weights, rmse = RMSE, active_count = number of non-zero weights)
    
    n <- nrow(E)
    p <- ncol(E)
    
    if (n < p) stop("solve_weights_ols: more endmembers than bands")
    if (any(!is.finite(E)) || any(!is.finite(y))) stop("solve_weights_ols: non-finite values in E or y")
    
    # Use geometric projection with simplex constraint
    # Solve unconstrained least squares: w = (E'E)^(-1) E'y
    EtE <- t(E) %*% E
    Ety <- t(E) %*% y
    
    # Add small ridge for numerical stability
    ridge <- 1e-8 * diag(p)
    
    w_unconstrained <- tryCatch({
      solve(EtE + ridge, Ety)
    }, error = function(e) {
      # Fallback to pseudoinverse
      if (requireNamespace("MASS", quietly = TRUE)) {
        MASS::ginv(EtE + ridge) %*% Ety
      } else {
        rep(1/p, p)  # Equal weights fallback
      }
    })
    
    w <- as.numeric(w_unconstrained)
    
    # Enforce constraints if requested
    if (enforce_constraints) {
      # Project onto probability simplex (sum=1, all >= 0)
      w <- project_to_simplex(w)
    }
    
    # Compute final RMSE after constraints
    pred <- E %*% w
    rmse <- sqrt(mean((y - pred)^2))
    
    list(
      w = w,
      rmse = rmse,
      active_count = sum(w > 1e-6)  # Number of "active" endmembers
    )
  }

  # Spectral angle function for geometric unmixing
  spectral_angle <- function(a, b) {
    cs <- sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
    acos(pmax(-1, pmin(1, cs)))
  }

  cos_angle <- function(a, b) {
    sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
  }

  geometric_select_best_partner <- function(y, M1, M2) {
    best_dist <- Inf
    best_m1 <- NULL
    best_m2 <- NULL
    for (i in seq_along(M1)) {
      m1 <- M1[[i]]$vec
      best_angle <- Inf
      best_m2_idx <- NULL
      for (j in seq_along(M2)) {
        m2 <- M2[[j]]$vec
        diff <- y - m1
        diff2 <- m1 - m2
        angle <- cos_angle(diff, diff2)
        if (angle < best_angle) {
          best_angle <- angle
          best_m2_idx <- j
        }
      }
      if (!is.null(best_m2_idx)) {
        m2 <- M2[[best_m2_idx]]$vec
        # Compute projection distance
        E <- cbind(m1, m2)
        w <- tryCatch(solve(E, y), error = function(e) rep(0, 2))
        w <- pmax(0, w)
        if (sum(w) > 0) w <- w / sum(w)
        dist <- sqrt(sum((y - E %*% w)^2))
        if (dist < best_dist) {
          best_dist <- dist
          best_m1 <- M1[[i]]
          best_m2 <- M2[[best_m2_idx]]
        }
      }
    }
    list(m1 = best_m1, m2 = best_m2, dist = best_dist)
  }

  # Fully Constrained Least Squares Unmixing (using geometric projection)
  # Minimize ||E w - y||^2 subject to sum(w) = 1, w >= 0
  solve_fclsu <- function(E, y, lambda = 1e-6) {
    # Use geometric simplex projection for constraints
    result <- solve_weights_ols(E, y, enforce_constraints = TRUE)
    list(w = result$w, residual = result$rmse * sqrt(nrow(E)))
  }

  unmix_stage2_geometric <- function(y, vegetation_libraries, topK = 2) {
    veg_types <- names(vegetation_libraries)
    if (length(veg_types) == 0) return(NULL)
    X <- as.numeric(y); ny <- sqrt(sum(X^2)); X_norm <- if (ny > 0) X / ny else X
    top_variants <- list()
    for (v in veg_types) {
      cand <- vegetation_libraries[[v]]; if (length(cand) == 0) next
      sims <- sapply(cand, function(x) { vec <- if (is.list(x) && !is.null(x$vec)) x$vec else x; norm_x <- sqrt(sum(vec^2)); if (norm_x == 0 || ny == 0) return(0); sum(vec * X) / (norm_x * ny) })
      ord <- order(sims, decreasing = TRUE); keep <- ord[seq_len(min(topK, length(ord)))]; top_variants[[v]] <- lapply(keep, function(i) { x <- cand[[i]]; if (is.list(x) && !is.null(x$vec)) x else list(vec = x, id = paste0(v, "_", i)) })
    }
    if (length(top_variants) == 0) return(NULL)
    
    # Use exhaustive combinatorial search instead of greedy angle-based search
    # This guarantees finding the best combination among the topK variants
    result <- evaluate_all_combinations(X_norm, top_variants, lambda = 0)
    
    if (is.null(result)) return(NULL)
    list(proportions = result$w, chosen_variants = result$ids, residual = result$rmse)
  }

  geometric_stage1_with_selection <- function(y, barren_endmember, vegetation_endmember) {
    if (is.null(barren_endmember) || is.null(vegetation_endmember)) return(list(veg_frac = NA_real_, barren_frac = NA_real_, residual = NA_real_))
    B <- as.numeric(barren_endmember); V <- as.numeric(vegetation_endmember); X <- as.numeric(y)
    len <- min(length(X), length(B), length(V)); X <- X[1:len]; B <- B[1:len]; V <- V[1:len]
    res <- geometric_project_and_unmix(X, B, V)
    barren_frac <- res$f1; veg_frac <- res$f2
    veg_frac <- max(0, min(1, veg_frac)); barren_frac <- 1 - veg_frac
    list(veg_frac = veg_frac, barren_frac = barren_frac, residual = res$residual, projection = res$y_proj)
  }

  # Geometric selection for two endmember sets
  geometric_select_pair <- function(y, M1, M2) {
    best_dist <- Inf
    best_m1 <- NULL
    best_m2 <- NULL
    for (i in seq_along(M1)) {
      m1 <- M1[[i]]$vec
      best_angle <- Inf
      best_m2_idx <- NULL
      for (j in seq_along(M2)) {
        m2 <- M2[[j]]$vec
        diff <- y - m1
        diff2 <- m1 - m2
        angle <- spectral_angle(diff, diff2)
        if (angle < best_angle) {
          best_angle <- angle
          best_m2_idx <- j
        }
      }
      if (!is.null(best_m2_idx)) {
        m2 <- M2[[best_m2_idx]]$vec
        # Compute projection distance
        E <- cbind(m1, m2)
        w <- tryCatch(solve(E, y), error = function(e) rep(0, 2))
        w <- pmax(0, w)
        if (sum(w) > 0) w <- w / sum(w)
        dist <- sqrt(sum((y - E %*% w)^2))
        if (dist < best_dist) {
          best_dist <- dist
          best_m1 <- M1[[i]]
          best_m2 <- M2[[best_m2_idx]]
        }
      }
    }
    list(m1 = best_m1, m2 = best_m2, dist = best_dist)
  }

  # ----------------------------------------------------------------------------
  # GEOMETRIC PROJECTION-BASED UNMIXING (Two Endmembers)
  # ----------------------------------------------------------------------------
  geometric_project_and_unmix <- function(y, m1, m2) {
    em_line <- m2 - m1
    em_norm_sq <- sum(em_line^2)
    if (em_norm_sq < 1e-10) {
      # Degenerate: identical endmembers
      return(list(f1 = 0.5, f2 = 0.5, y_proj = m1, residual = sqrt(sum((y - m1)^2))))
    }
    y_minus_m1 <- y - m1
    t <- sum(y_minus_m1 * em_line) / em_norm_sq
    t_clamped <- max(0, min(1, t))
    y_proj <- m1 + t_clamped * em_line
    residual <- sqrt(sum((y - y_proj)^2))
    f2 <- t_clamped
    f1 <- 1 - f2
    list(f1 = f1, f2 = f2, y_proj = y_proj, residual = residual, t = t)
  }

  # ----------------------------------------------------------------------------
  # GEOMETRIC UNMIXING FOR MULTIPLE ENDMEMBERS (N > 2)
  # ----------------------------------------------------------------------------
  geometric_unmix_simplex <- function(y, M, enforce_sum_to_one = TRUE) {
    N <- ncol(M)
    P <- nrow(M)
    if (N == 1) {
      residual <- sqrt(sum((y - M[, 1])^2))
      return(list(f = 1, residual = residual, y_proj = M[, 1]))
    }
    if (N == 2) {
      res <- geometric_project_and_unmix(y, M[,1], M[,2])
      return(list(f = c(res$f1, res$f2), residual = res$residual, y_proj = res$y_proj))
    }
    MtM <- t(M) %*% M
    Mty <- t(M) %*% y
    ridge <- 1e-8 * diag(N)
    f_unconstrained <- tryCatch({
      solve(MtM + ridge, Mty)
    }, error = function(e) {
      if (requireNamespace("MASS", quietly = TRUE)) MASS::ginv(MtM + ridge) %*% Mty else rep(1 / N, N)
    })
    f <- project_to_simplex(as.numeric(f_unconstrained))
    y_proj <- as.numeric(M %*% f)
    residual <- sqrt(sum((y - y_proj)^2))
    list(f = f, residual = residual, y_proj = y_proj)
  }

  # ----------------------------------------------------------------------------
  # ANGLE-BASED MESMA (Two or N libraries)
  # ----------------------------------------------------------------------------
  angle_based_mesma <- function(y, library_list, max_components = NULL) {
    veg_names <- names(library_list)
    n_libs <- length(veg_names)
    if (n_libs == 0) return(NULL)
    if (n_libs == 1) {
      lib <- library_list[[1]]
      best_dist <- Inf
      best_idx <- 1
      for (i in seq_along(lib)) {
        dist <- sqrt(sum((y - lib[[i]]$vec)^2))
        if (dist < best_dist) {
          best_dist <- dist; best_idx <- i
        }
      }
      return(list(fractions = setNames(1.0, veg_names[1]), chosen = setNames(lib[[best_idx]]$id, veg_names[1]), residual = best_dist))
    }
    if (n_libs == 2) {
      res <- geometric_select_pair(y, library_list[[1]], library_list[[2]])
      if (is.null(res$m1) || is.null(res$m2)) return(NULL)
      unmix_res <- geometric_project_and_unmix(y, res$m1$vec, res$m2$vec)
      fractions <- c(unmix_res$f1, unmix_res$f2)
      names(fractions) <- veg_names
      chosen <- c(res$m1$id, res$m2$id); names(chosen) <- veg_names
      return(list(fractions = fractions, chosen = chosen, residual = unmix_res$residual))
    }
    # Hierarchical for >2 libs: iterative selection
    lib_sizes <- sapply(library_list, length)
    lib_order <- order(lib_sizes)
    v1 <- veg_names[lib_order[1]]; v2 <- veg_names[lib_order[2]]
    pair_result <- geometric_select_pair(y, library_list[[v1]], library_list[[v2]])
    if (is.null(pair_result$m1)) return(NULL)
    selected_endmembers <- list(); selected_endmembers[[v1]] <- pair_result$m1; selected_endmembers[[v2]] <- pair_result$m2
    for (k in seq(3, n_libs)) {
      vk <- veg_names[lib_order[k]]
      best_dist <- Inf; best_em <- NULL
      for (candidate in library_list[[vk]]) {
        M_cols <- lapply(selected_endmembers, function(em) em$vec)
        M_cols[[length(M_cols) + 1]] <- candidate$vec
        M <- do.call(cbind, M_cols)
        unmix_result <- geometric_unmix_simplex(y, M)
        if (unmix_result$residual < best_dist) {
          best_dist <- unmix_result$residual; best_em <- candidate
        }
      }
      if (!is.null(best_em)) selected_endmembers[[vk]] <- best_em
    }
    M_final <- do.call(cbind, lapply(selected_endmembers, function(em) em$vec))
    colnames(M_final) <- names(selected_endmembers)
    final_unmix <- geometric_unmix_simplex(y, M_final)
    fractions <- final_unmix$f; names(fractions) <- names(selected_endmembers)
    chosen <- sapply(selected_endmembers, function(em) em$id)
    list(fractions = fractions, chosen = chosen, residual = final_unmix$residual)
  }

  # ----------------------------------------------------------------------------
  # STAGE 1 AND STAGE 2 GEOMETRIC HIERACHICAL WRAPPERS
  # ----------------------------------------------------------------------------
  geometric_stage1_unmix <- function(y, barren_endmember, vegetation_endmember) {
      if (is.null(barren_endmember) || is.null(vegetation_endmember)) return(list(veg_frac = NA_real_, barren_frac = NA_real_, residual = NA_real_))
      # Convert to numeric and trim to aligned length
      B <- as.numeric(barren_endmember); V <- as.numeric(vegetation_endmember); X <- as.numeric(y)
      len <- min(length(X), length(B), length(V)); X <- X[1:len]; B <- B[1:len]; V <- V[1:len]
      res <- geometric_project_and_unmix(X, B, V)
    barren_frac <- res$f1; veg_frac <- res$f2
    veg_frac <- max(0, min(1, veg_frac)); barren_frac <- 1 - veg_frac
    list(veg_frac = veg_frac, barren_frac = barren_frac, residual = res$residual, projection = res$y_proj)
  }

  geometric_stage2_unmix <- function(y, vegetation_libraries, topK = 2) {
    veg_types <- names(vegetation_libraries)
    if (length(veg_types) == 0) return(NULL)
    X <- as.numeric(y); ny <- sqrt(sum(X^2)); X_norm <- if (ny > 0) X / ny else X
    top_variants <- list()
    for (v in veg_types) {
      cand <- vegetation_libraries[[v]]; if (length(cand) == 0) next
      sims <- sapply(cand, function(x) { vec <- if (is.list(x) && !is.null(x$vec)) x$vec else x; norm_x <- sqrt(sum(vec^2)); if (norm_x == 0 || ny == 0) return(0); sum(vec * X) / (norm_x * ny) })
      ord <- order(sims, decreasing = TRUE); keep <- ord[seq_len(min(topK, length(ord)))]; top_variants[[v]] <- lapply(keep, function(i) { x <- cand[[i]]; if (is.list(x) && !is.null(x$vec)) x else list(vec = x, id = paste0(v, "_", i)) })
    }
    if (length(top_variants) == 0) return(NULL)
    result <- angle_based_mesma(X_norm, top_variants)
    if (is.null(result)) return(NULL)
    list(proportions = result$fractions, chosen_variants = result$chosen, residual = result$residual)
  }

  hierarchical_geometric_mesma <- function(y, barren_endmember, vegetation_endmember, vegetation_libraries, topK = 2) {
    X <- as.numeric(y)
    # Stage 1
    if (is.list(barren_endmember) || is.list(vegetation_endmember)) {
      B_vars <- if (is.list(barren_endmember)) barren_endmember else list(list(vec = as.numeric(barren_endmember), id = 'B1'))
      V_vars <- if (is.list(vegetation_endmember)) vegetation_endmember else list(list(vec = as.numeric(vegetation_endmember), id = 'V1'))
      stage1_result <- geometric_select_pair(X, B_vars, V_vars)
      if (is.null(stage1_result$m1) || is.null(stage1_result$m2)) stage1_result <- list(m1 = B_vars[[1]], m2 = V_vars[[1]], dist = NA_real_)
      unmix_stage1 <- geometric_project_and_unmix(X, stage1_result$m1$vec, stage1_result$m2$vec)
      veg_fraction_total <- unmix_stage1$f2; barren_fraction_total <- unmix_stage1$f1
    } else {
      stage1_unmix <- geometric_stage1_unmix(X, barren_endmember, vegetation_endmember)
      veg_fraction_total <- stage1_unmix$veg_frac; barren_fraction_total <- stage1_unmix$barren_frac
    }
    if (!is.finite(veg_fraction_total)) return(list(fractions = NULL, stage1 = stage1_result, stage2 = NULL, error = 'Stage 1 failed'))
    if (veg_fraction_total <= 0.01) return(list(fractions = c(barren = 1.0), stage1 = stage1_result, stage2 = NULL, veg_fraction = 0, barren_fraction = 1))
    stage2_result <- geometric_stage2_unmix(X, vegetation_libraries, topK = topK)
    if (is.null(stage2_result)) return(list(fractions = c(barren = barren_fraction_total, vegetation = veg_fraction_total), stage1 = stage1_result, stage2 = NULL, veg_fraction = veg_fraction_total, barren_fraction = barren_fraction_total))
    stage2_proportions <- stage2_result$proportions
    final_fractions <- stage2_proportions * veg_fraction_total
    final_fractions <- c(final_fractions, barren = barren_fraction_total)
    total <- sum(final_fractions); if (abs(total - 1) > 1e-6 && total > 0) final_fractions <- final_fractions / total
    list(fractions = final_fractions, stage1 = stage1_result, stage2 = stage2_result, veg_fraction = veg_fraction_total, barren_fraction = barren_fraction_total, chosen_variants = stage2_result$chosen_variants, residual = list(stage1 = if (exists('unmix_stage1')) unmix_stage1$residual else NA_real_, stage2 = stage2_result$residual))
  }

  # Retrieve precomputed templates into the structure used by downstream code
  get_variant_templates <- function(veg_types, grid_type, compressed_templates, mesma_lib) {
    templates <- list()
    for (veg in veg_types) {
      veg_templates <- list()
      for (variant in mesma_lib[[veg]]) {
        vec <- compressed_templates[[veg]][[variant$variant_id]][[grid_type]]
        if (!is.null(vec) && length(vec) > 0) {
          veg_templates[[length(veg_templates) + 1]] <- list(vec = vec, id = variant$variant_id)
        }
      }
      templates[[veg]] <- veg_templates
    }
    templates
  }

  # Vectorized evaluation of combinations with chunking
  evaluate_all_combinations <- function(y, templates, lambda, early_stop_rmse = EARLY_STOP_RMSE_THRESHOLD) {
    veg_names <- names(templates)
    n_veg <- length(veg_names)
    if (n_veg == 0) return(NULL)
    # all indices per veg
    idx_lists <- lapply(templates, function(lst) seq_along(lst))
    if (any(lengths(idx_lists) == 0)) return(NULL)

    # Avoid materializing the full Cartesian product in memory (expand.grid)
    # which can easily explode (OOM) when many veg types * candidates per veg.
    radices <- as.integer(lengths(idx_lists))
    n_veg <- length(radices)
    n_combos <- as.numeric(Reduce(`*`, as.list(radices), init = 1L))
    chunk_size <- 100L

    # Safety thresholds
    safe_expand_limit <- COMBO_SAFE_EXPAND_LIMIT    # small enough to expand safely in memory
    abort_limit <- COMBO_ABORT_LIMIT                # too many combinations — abort rather than OOM

    use_expanded_grid <- !is.na(n_combos) && n_combos <= safe_expand_limit
    if (use_expanded_grid) combos <- do.call(expand.grid, idx_lists)
    if (!is.na(n_combos) && n_combos > abort_limit) {
      stop(sprintf("Combination space is too large (%d combinations). Reduce TOPK_VARIANTS, MAX_VEG_COMPONENTS, or limit vegetation types to avoid huge memory allocations.", as.integer(n_combos)))
    }
    best_rmse <- Inf
    best_result <- NULL
    for (start in seq(1, n_combos, by = chunk_size)) {
      end <- min(start + chunk_size - 1L, n_combos)
      if (use_expanded_grid) {
        # small search space — use the materialized combos directly
        chunk <- combos[start:end, , drop = FALSE]
      } else {
        # lazily generate the mixed-radix indices for this chunk to avoid full expand.grid
        nrows <- end - start + 1L
        chunk <- matrix(NA_integer_, nrow = nrows, ncol = n_veg)
        colnames(chunk) <- names(idx_lists)
        # 0-based start for easier mixed-radix math
        base_start <- as.numeric(start - 1L)
        for (r in seq_len(nrows)) {
          x <- base_start + (r - 1L)
          # compute indices in each digit position
          rem <- x
          for (j in seq_len(n_veg)) {
            base <- radices[j]
            # mixed-radix digit (1-based index into candidate list)
            digit <- (rem %% base) + 1L
            chunk[r, j] <- as.integer(digit)
            rem <- rem %/% base
          }
        }
        # convert to data.frame-like access used later
        chunk <- as.data.frame(chunk, stringsAsFactors = FALSE, check.names = FALSE)
      }
      # build E matrices for the chunk
      results <- vector("list", nrow(chunk))
      for (i in seq_len(nrow(chunk))) {
        cols <- list(); ids <- character(0)
        for (v in veg_names) {
          idx <- as.integer(chunk[i, v])
          cols[[length(cols) + 1]] <- templates[[v]][[idx]]$vec
          ids <- c(ids, templates[[v]][[idx]]$id)
        }
        E <- do.call(cbind, cols)
        for (j in seq_len(ncol(E))) { nj <- sqrt(sum(E[, j]^2)); if (nj > 0) E[, j] <- E[, j] / nj }
        # Always use Ordinary Least Squares (OLS) for consistency
        w_res <- solve_weights_ols(E, y, enforce_constraints = TRUE)
        w <- if (is.list(w_res) && !is.null(w_res$w)) as.numeric(w_res$w) else as.numeric(w_res)
        pred <- as.numeric(E %*% w)
        rmse <- sqrt(mean((y - pred)^2))
        results[[i]] <- list(w = w, rmse = rmse, ids = ids)
      }
      # pick best in chunk
      rmses <- vapply(results, function(r) r$rmse, numeric(1))
      i_best <- which.min(rmses)
      if (is.finite(rmses[i_best]) && rmses[i_best] < best_rmse) {
        best_rmse <- rmses[i_best]
        best_result <- results[[i_best]]
        names(best_result$w) <- veg_names
        names(best_result$ids) <- veg_names
      }
      if (is.finite(early_stop_rmse) && early_stop_rmse > 0 && best_rmse <= early_stop_rmse) break
    }
    best_result
  }

# Compress a single trace to a feature vector (using raw spectral indices)
compress_trace <- function(dly_year, avail_idx, budget = 365L) {
  # Simply build raw 365 × K matrix and flatten
  raw_mat <- build_raw_365_matrix(dly_year, avail_idx)
  if (is.null(raw_mat)) return(NULL)
  
  y <- as.numeric(raw_mat)  # Flatten to vector
  if (any(!is.finite(y))) return(NULL)
  
  list(y = y, grid_type = "full", raw_mat = raw_mat)
}



compress_and_unmix_year <- function(dly_year, mesma_lib, budget = TEMPORAL_BUDGET, topK = TOPK_VARIANTS) {
  res <- compress_trace(dly_year, avail, budget)
  if (is.null(res)) return(NULL)
  unmix_stage2_compressed(res$y, res$grid_type, mesma_lib, topK, feature_weights = LDA_FEATURE_WEIGHTS)
}

# Integration function: use geometric hierarchical unmix on compressed templates
geometric_hierarchical_unmix_compressed <- function(y, compressed_stage1_lib, mesma_lib, grid_type = 'medium', topK = TOPK_VARIANTS) {
  
  # --- Always use Geometric MESMA ---
  if (is.null(compressed_stage1_lib) || is.null(compressed_stage1_lib$barren) || is.null(compressed_stage1_lib$vegetation)) {
    return(list(error = 'Missing Stage 1 compressed library'))
  }
  stage1_res <- geometric_stage1_unmix(y, compressed_stage1_lib$barren, compressed_stage1_lib$vegetation)
  veg_frac <- stage1_res$veg_frac; barren_frac <- stage1_res$barren_frac
  
  if (!is.finite(veg_frac)) return(list(error = 'Stage 1 failed'))
  if (veg_frac <= 0.01) return(list(vegetation_proportions = c(barren = 1.0), chosen_variants = NULL, rmse = NA_real_, veg_fraction = 0, barren_fraction = 1))

  # Build vegetation libraries for stage 2 using compressed templates accessor
  if (!exists('.COMPRESSED_TEMPLATES_ACCESSOR', envir = globalenv())) return(list(error = 'Missing compressed templates accessor'))
  compressed_templates <- get('.COMPRESSED_TEMPLATES_ACCESSOR', envir = globalenv())
  veg_types <- names(mesma_lib)
  veg_libs <- list()
  for (veg in veg_types) {
    veg_templates <- list()
    for (variant in mesma_lib[[veg]]) {
      vec <- compressed_templates[[veg]][[variant$variant_id]][[grid_type]]
      if (!is.null(vec)) veg_templates[[length(veg_templates) + 1]] <- list(vec = vec, id = variant$variant_id)
    }
    if (length(veg_templates) > 0) veg_libs[[veg]] <- veg_templates
  }
  if (length(veg_libs) == 0) return(list(vegetation_proportions = c(vegetation = veg_frac, barren = barren_frac), chosen_variants = NULL, rmse = NA_real_, veg_fraction = veg_frac, barren_fraction = barren_frac))
  
  # Stage 2: geometric
  # Note: geometric_stage2_unmix normalizes inputs to unit length for shape matching
  ny <- sqrt(sum(y^2)); y_norm <- if (ny > 0) y / ny else y
  stage2_res <- geometric_stage2_unmix(y_norm, veg_libs, topK = topK)
  
  if (is.null(stage2_res)) return(list(vegetation_proportions = c(vegetation = veg_frac, barren = barren_frac), chosen_variants = NULL, rmse = NA_real_, veg_fraction = veg_frac, barren_fraction = barren_frac))
  
  final_veg_fractions <- stage2_res$proportions * veg_frac
  final_fractions <- c(final_veg_fractions, barren = barren_frac)
  
  # Normalize if needed (should sum to 1)
  total <- sum(final_fractions)
  if (abs(total - 1) > 1e-6 && total > 0) final_fractions <- final_fractions / total
  
  list(vegetation_proportions = final_fractions, chosen_variants = stage2_res$chosen_variants, rmse = stage2_res$residual, veg_fraction = veg_frac, barren_fraction = barren_frac)
}

  # End of mesma_dynamic_programming inner helpers
  # Helper to prepare factor data
  prepare_factor_data <- function(dly, gpca, avail_idx, veg_type) {
    date_list <- list()
    dly$doy <- lubridate::yday(dly$date)
    dts <- sort(unique(dly$date))

    for (i_dt in seq_along(dts)) {
      dt <- dts[i_dt]
      sub <- dly[dly$date == dt, , drop = FALSE]
      sub <- sub[is.finite(sub$doy), , drop = FALSE]
      if (nrow(sub) == 0) next

      idx_present_orig <- intersect(avail_idx, names(sub))
      min_required <- ceiling(length(avail_idx) * MIN_IDX_PRESENCE)
      if (length(idx_present_orig) < min_required) {
        next
      }

      vals_aug <- c()
      vrows_list <- list()

      for (idx in idx_present_orig) {
        vals_idx <- sub[[idx]]
        vals_idx <- vals_idx[is.finite(vals_idx)]
        if (length(vals_idx) == 0) {
          vals_aug <- NULL
          break
        }

        # Medoid selection for duplicate observations (closest to mean)
        yy <- if (length(vals_idx) == 1) {
          as.numeric(vals_idx[1])
        } else {
          m <- mean(vals_idx, na.rm = TRUE)
          as.numeric(vals_idx[which.min(abs(vals_idx - m))])
        }
        if (!is.finite(yy)) next

        kpos <- match(idx, gpca$idx_order)
        if (is.na(kpos)) next

        mi <- gpca$col_means[kpos]
        # Do not apply per-band scaling — assume identity scaling (1.0)
        valid_veg <- !is.null(veg_type) && !is.na(veg_type) && nzchar(as.character(veg_type))
        veg_type_chr <- if (valid_veg) as.character(veg_type) else NA_character_
        veg_type_lc <- if (valid_veg) tolower(veg_type_chr) else NA_character_

        si <- 1.0
        vals_aug <- c(vals_aug, as.numeric((yy - mi) / si))
        vrows_list[[length(vrows_list) + 1]] <- gpca$V[kpos, , drop = FALSE]
      }
      if (is.null(vals_aug)) next

      for (idx in idx_present_orig) {
        mv_col <- paste0(idx, "_var14")
        if (!mv_col %in% names(sub)) next

        vals_mv <- sub[[mv_col]]
        vals_mv <- vals_mv[is.finite(vals_mv)]
        if (length(vals_mv) == 0) {
          vals_aug <- NULL
          break
        }

        # Medoid selection for duplicate moving-variance observations
        yy_mv <- if (length(vals_mv) == 1) {
          as.numeric(vals_mv[1])
        } else {
          m_mv <- mean(vals_mv, na.rm = TRUE)
          as.numeric(vals_mv[which.min(abs(vals_mv - m_mv))])
        }
        if (!is.finite(yy_mv)) next

        mv_name <- paste0(idx, "_mv")
        kpos_mv <- match(mv_name, gpca$idx_order)
        if (is.na(kpos_mv)) next

        mi_mv <- gpca$col_means[kpos_mv]
        # No per-band scaling: assume scale = 1
        si_mv <- 1.0
        vals_aug <- c(vals_aug, as.numeric((yy_mv - mi_mv) / si_mv))
        vrows_list[[length(vrows_list) + 1]] <- gpca$V[kpos_mv, , drop = FALSE]
      }
      if (is.null(vals_aug)) next

      if (length(vrows_list) == 0) next
      vsub <- do.call(rbind, vrows_list)
      if (ncol(vsub) < 1) next

      xtx <- crossprod(vsub)
      xty <- as.numeric(crossprod(vsub, matrix(vals_aug, ncol = 1)))

      ridge_small <- 0.01
      xtx_reg <- xtx + ridge_small * diag(ncol(xtx))
      z <- solve(xtx_reg, xty)

      if (any(!is.finite(z))) {
        stop("All matrix inversion methods failed")
      }

      date_list[[length(date_list) + 1]] <- list(doy = as.integer(sub$doy[1]), z = as.numeric(z))
    }
    date_list
  }

  # Missing utility functions for MESMA
  build_Z365 <- function(date_list, k) {
    # Build a 365 x k matrix from the date_list factor projections
    Z <- matrix(NA_real_, nrow = 365, ncol = k)
    for (entry in date_list) {
      if (!is.null(entry$doy) && !is.null(entry$z) && 
          entry$doy >= 1 && entry$doy <= 365 && 
          length(entry$z) == k) {
        Z[entry$doy, ] <- entry$z
      }
    }
    # Interpolate missing DOYs per component
    for (col_idx in seq_len(k)) {
      finite_mask <- is.finite(Z[, col_idx])
      if (sum(finite_mask) >= 2) {
        Z[!finite_mask, col_idx] <- approx(
          x = which(finite_mask),
          y = Z[finite_mask, col_idx],
          xout = which(!finite_mask),
          rule = 2
        )$y
      } else if (sum(finite_mask) == 1) {
        Z[, col_idx] <- Z[finite_mask, col_idx][1]
      } else {
        Z[, col_idx] <- 0
      }
    }
    Z[!is.finite(Z)] <- 0
    Z
  }

  precompute_whitening <- function(mesma_lib, compressed_templates, combinations) {
    # Placeholder for whitening precomputation - returns a simple accessor
    function(combo_key) {
      # Simple identity whitening for now
      list(whiten = function(x) x, unwhiten = function(x) x)
    }
  }

  compute_diagnostics <- function(y, E, w, mesma_result = NULL) {
    # Compute basic diagnostic statistics
    if (!is.matrix(E) || ncol(E) == 0 || length(w) != ncol(E)) {
      return(data.frame(
        condition_number = NA_real_,
        residual_sum_of_squares = NA_real_,
        r_squared = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    
    pred <- as.numeric(E %*% w)
    residuals <- y - pred
    rss <- sum(residuals^2)
    tss <- sum((y - mean(y))^2)
    r_squared <- if (tss > 0) 1 - rss/tss else 0
    
    # Condition number of E
    cond_num <- tryCatch({
      s <- svd(E)$d
      if (length(s) > 0 && min(s) > 0) max(s) / min(s) else stop("compute_diagnostics: singular values invalid for condition number calculation")
    }, error = function(e) stop(sprintf("compute_diagnostics: failed computing SVD for condition number: %s", e$message)))
    
    data.frame(
      condition_number = cond_num,
      residual_sum_of_squares = rss,
      r_squared = r_squared,
      stringsAsFactors = FALSE
    )
  }

  # Build MESMA library
  cat("Building MESMA endmember library...\n")
  
  # STEP 1: Compute global PCA-LDA weights for feature selection
  cat("\n=== COMPUTING GLOBAL PCA-LDA WEIGHTS ===\n")
  
  PCA_LDA_RESULT <- compute_pca_lda_weights(
    lib_df = lib_df,
    avail_idx = avail,
    pca_variance_threshold = 0.95,
    lda_weight_floor = 0.01
  )
  
  if (!is.null(PCA_LDA_RESULT)) {
    LDA_FEATURE_WEIGHTS <- PCA_LDA_RESULT$weights
    cat(sprintf("PCA-LDA weights computed: %d features, %d PCs retained\n", 
                length(LDA_FEATURE_WEIGHTS), PCA_LDA_RESULT$n_pcs))
  } else {
    LDA_FEATURE_WEIGHTS <- NULL
    warning("PCA-LDA weight computation failed, using unweighted features")
  }
  
  # Store globally for use in unmixing
  assign("LDA_FEATURE_WEIGHTS", LDA_FEATURE_WEIGHTS, envir = globalenv())
  assign("PCA_LDA_RESULT", PCA_LDA_RESULT, envir = globalenv())
  
  # STEP 2: Process traces at full resolution (no time compression)
  cat("\n=== PROCESSING TRACES AT FULL YEAR RESOLUTION ===\n")
  
  reduced_traces <- reduce_all_traces_simple(
    lib_df = lib_df,
    veg_types = names(lib),
    avail_idx = avail
  )
  
  # Print rows per class after processing
  cat("\n=== Rows per class after full-resolution processing ===\n")
  for (v in names(reduced_traces)) {
    cat(sprintf("  %s: %d rows\n", v, nrow(reduced_traces[[v]]$features)))
  }
  cat("=============================================\n\n")
  
  # STEP 3: Build variants using weighted clustering
  mesma_lib <- build_mesma_variants_weighted(
    reduced_data = reduced_traces, 
    raw_lib_templates = lib,
    pca_lda_weights = LDA_FEATURE_WEIGHTS,
    n_variants = N_VARIANTS_PER_VEG,
    min_cluster_size = MIN_CLUSTER_SIZE
  )

  # Precompute compressed templates accessor
  precompute_compressed_templates <- function(mesma_lib, grid_name = "full") {
    template_db <- list()
    
    # Debug: Check raw_mat matrices before compression
    cat("\n=== RAW INDEX MATRIX DEBUG (before compression) ===\n")
    for (veg in names(mesma_lib)) {
      for (variant in mesma_lib[[veg]]) {
        raw_mat <- variant$raw_mat
        if (!is.null(raw_mat)) {
          cat(sprintf("  %s/%s: raw_mat dim=%dx%d, range=[%.4f, %.4f], mean=%.4f\n",
                      veg, variant$variant_id, nrow(raw_mat), ncol(raw_mat),
                      min(raw_mat, na.rm=TRUE), max(raw_mat, na.rm=TRUE), mean(raw_mat, na.rm=TRUE)))
        } else {
          cat(sprintf("  %s/%s: raw_mat is NULL!\n", veg, variant$variant_id))
        }
      }
    }
    cat("==============================================\n\n")
    
    for (veg in names(mesma_lib)) {
      template_db[[veg]] <- list()
      for (variant in mesma_lib[[veg]]) {
        template_db[[veg]][[variant$variant_id]] <- list()
        
        # Since no temporal compression, use full raw_mat as features
        compressed_vec <- as.numeric(variant$raw_mat)
        
        # DEBUG: Check for degenerate compressed vectors
        vec_sd <- sd(compressed_vec, na.rm = TRUE)
        vec_range <- diff(range(compressed_vec, na.rm = TRUE))
        vec_unique <- length(unique(compressed_vec[is.finite(compressed_vec)]))
        vec_norm <- sqrt(sum(compressed_vec^2))
        
        if (!is.finite(vec_sd) || vec_sd < 1e-10 || 
            !is.finite(vec_range) || vec_range < 1e-10 || 
            vec_unique < 3 || vec_norm < 1e-10) {
          cat(sprintf("WARNING: Degenerate compressed vector for %s/%s: sd=%.2e, range=%.2e, unique=%d, norm=%.2e\n",
                      veg, variant$variant_id, vec_sd, vec_range, vec_unique, vec_norm))
        }
        
        # Debug: Show compressed vector stats
        # cat(sprintf("  [COMPRESS] %s/%s: compressed length=%d, range=[%.4f, %.4f]\n",
        #             veg, variant$variant_id, length(compressed_vec),
        #             min(compressed_vec), max(compressed_vec)))
        
        template_db[[veg]][[variant$variant_id]][["full"]] <- compressed_vec
      }
    }
    template_db
  }
  compressed_templates_accessor <- precompute_compressed_templates(mesma_lib, "full")
  assign(".COMPRESSED_TEMPLATES_ACCESSOR", compressed_templates_accessor, envir = globalenv())

  # Analyze inter-class similarity
  analyze_library_similarity(mesma_lib, compressed_templates_accessor, grid_type = "full")

  # Compute LDA feature weights
  LDA_FEATURE_WEIGHTS <- compute_lda_weights(mesma_lib, compressed_templates_accessor, grid_type = "full")
  assign("LDA_FEATURE_WEIGHTS", LDA_FEATURE_WEIGHTS, envir = globalenv())

  # Diagnostic: Log compressed template dimensions for each veg/vid/gt
  cat("Compressed templates accessor dimensions:\n")
  for (veg in names(compressed_templates_accessor)) {
    for (vid in names(compressed_templates_accessor[[veg]])) {
      for (gt in names(compressed_templates_accessor[[veg]][[vid]])) {
        vec <- compressed_templates_accessor[[veg]][[vid]][[gt]]
        cat(sprintf("  %s/%s/%s: length=%d\n", veg, vid, gt, length(vec)))
      }
    }
  }

  # Inter-class similarity check
  INTER_CLASS_SIMILARITY <- compute_inter_class_similarity_table(mesma_lib, compressed_templates_accessor, grid_type = "full")
  if (!is.null(INTER_CLASS_SIMILARITY)) {
    cat("\n=== INTER-CLASS VARIANT SIMILARITY CHECK ===\n")
    # Filter for high similarity
    high_sim <- INTER_CLASS_SIMILARITY[INTER_CLASS_SIMILARITY$cos_sim > 0.95, ]
    if (nrow(high_sim) > 0) {
      cat(sprintf("Found %d pairs of variants from DIFFERENT classes with cos_sim > 0.95:\n", nrow(high_sim)))
      print(high_sim[order(-high_sim$cos_sim), ])
    } else {
      cat("No high similarity (>0.95) detected between variants of different classes.\n")
    }
    assign("INTER_CLASS_SIMILARITY", INTER_CLASS_SIMILARITY, envir = globalenv())
    
    # Reclassify variants that are too similar to other classes as "unknown"
    if (!is.null(INTER_CLASS_SIMILARITY) && nrow(INTER_CLASS_SIMILARITY) > 0) {
      threshold <- 0.95
      unknown_threshold <- threshold
      
      # Collect all variants
      all_variants <- list()
      for (veg in names(mesma_lib)) {
        for (variant in mesma_lib[[veg]]) {
          all_variants[[variant$variant_id]] <- list(veg = veg, variant = variant)
        }
      }
      
      # For each variant, check max inter-class similarity
      variants_to_unknown <- c()
      for (vid in names(all_variants)) {
        var_info <- all_variants[[vid]]
        current_veg <- var_info$veg
        
        # Find similarities to other classes
        inter_sims <- INTER_CLASS_SIMILARITY[
          (INTER_CLASS_SIMILARITY$variant_a == vid | INTER_CLASS_SIMILARITY$variant_b == vid) &
          INTER_CLASS_SIMILARITY$veg_a != INTER_CLASS_SIMILARITY$veg_b,
        ]
        
        if (nrow(inter_sims) > 0) {
          max_inter_sim <- max(inter_sims$cos_sim, na.rm = TRUE)
          if (max_inter_sim > unknown_threshold) {
            variants_to_unknown <- c(variants_to_unknown, vid)
            cat(sprintf("Reclassifying %s from %s to unknown (max inter-class sim=%.3f)\n", 
                        vid, current_veg, max_inter_sim))
          }
        }
      }
      
      # Move variants to unknown
      if (length(variants_to_unknown) > 0) {
        mesma_lib[["unknown"]] <- list()
        for (vid in variants_to_unknown) {
          var_info <- all_variants[[vid]]
          # Remove from original veg
          mesma_lib[[var_info$veg]] <- mesma_lib[[var_info$veg]][
            sapply(mesma_lib[[var_info$veg]], function(x) x$variant_id != vid)
          ]
          # Add to unknown
          var_info$variant$variant_id <- paste0("unknown_", sub(".*_v", "v", vid))
          mesma_lib[["unknown"]] <- c(mesma_lib[["unknown"]], list(var_info$variant))
        }
        cat(sprintf("Moved %d variants to 'unknown' class\n", length(variants_to_unknown)))
      }
    }
    
    # --- FULL SIMILARITY HEATMAP ---
    cat("Generating full variant similarity heatmap...\n")
    
    # 1. Collect all vectors
    all_vecs <- list()
    all_ids <- c()
    all_vegs <- c()
    
    for (veg in names(mesma_lib)) {
      for (variant in mesma_lib[[veg]]) {
        vid <- variant$variant_id
        # Use the accessor to get the compressed vector
        vec <- NULL
        if (!is.null(compressed_templates_accessor[[veg]]) && 
            !is.null(compressed_templates_accessor[[veg]][[vid]]) &&
            !is.null(compressed_templates_accessor[[veg]][[vid]][["full"]])) {
           vec <- compressed_templates_accessor[[veg]][[vid]][["full"]]
        }
        
        if (!is.null(vec)) {
           all_vecs[[length(all_vecs)+1]] <- vec
           all_ids <- c(all_ids, vid)
           all_vegs <- c(all_vegs, veg)
        }
      }
    }
    
    if (length(all_vecs) > 1) {
       n_v <- length(all_vecs)
       sim_mat <- matrix(NA, n_v, n_v)
       rownames(sim_mat) <- all_ids
       colnames(sim_mat) <- all_ids
       
       for (i in 1:n_v) {
         for (j in 1:n_v) {
            sim_mat[i,j] <- pmax(pmin(cos_sim(all_vecs[[i]], all_vecs[[j]]), 1.01), -0.01)
         }
       }
       
       # Convert to long format for ggplot
       sim_df <- as.data.frame(as.table(sim_mat))
       colnames(sim_df) <- c("Var1", "Var2", "Similarity")
       
       # Order factors to group by vegetation type
       # Sort by veg type then variant ID
       ord_idx <- order(all_vegs, all_ids)
       ordered_ids <- all_ids[ord_idx]
       ordered_vegs <- all_vegs[ord_idx]
       
       sim_df$Var1 <- factor(sim_df$Var1, levels = ordered_ids)
       sim_df$Var2 <- factor(sim_df$Var2, levels = rev(ordered_ids)) # Reverse for y-axis to match matrix layout
       
       # Add lines between different species
       veg_changes <- which(diff(as.numeric(factor(ordered_vegs))) != 0)
       if (length(veg_changes) > 0) {
         vline_positions <- veg_changes + 0.5
         hline_positions <- length(ordered_ids) - veg_changes + 0.5
       } else {
         vline_positions <- numeric(0)
         hline_positions <- numeric(0)
       }
       
       p_heat <- ggplot(sim_df, aes(Var1, Var2, fill = Similarity)) +
         geom_tile() +
         scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.9, limits = c(-0.01, 1.01), na.value = "white") +
         theme_minimal() +
         theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
               axis.text.y = element_text(size = 6)) +
         labs(title = "Pairwise Cosine Similarity of All Variants", x = NULL, y = NULL)
         
       # Add vertical and horizontal lines between species
       if (length(vline_positions) > 0) {
         p_heat <- p_heat + 
           geom_vline(xintercept = vline_positions, color = "black", linewidth = 0.5) +
           geom_hline(yintercept = hline_positions, color = "black", linewidth = 0.5)
       }
         
       ggsave(file.path(OUT_DIR, "variant_similarity_heatmap.png"), p_heat, width = 12, height = 10)
       cat(sprintf("Saved similarity heatmap to: %s\n", file.path(OUT_DIR, "variant_similarity_heatmap.png")))
       
       # Also save the matrix
       write.csv(sim_mat, file.path(OUT_DIR, "variant_similarity_matrix.csv"))
    }
  }

  # ========== UNMIXING SPACE MULTICOLLINEARITY CHECK ==========
  cat("\n=== UNMIXING SPACE MULTICOLLINEARITY CHECK ===\n")
  
  # Build endmember matrix from compressed templates
  for (grid_type in c("full")) {
    cat(sprintf("\nGrid type: %s\n", grid_type))
    
    # Get one variant per vegetation type
    E_cols <- list()
    for (veg in names(mesma_lib)) {
      if (length(mesma_lib[[veg]]) > 0) {
        vid <- mesma_lib[[veg]][[1]]$variant_id
        if (!is.null(compressed_templates_accessor[[veg]]) && 
            !is.null(compressed_templates_accessor[[veg]][[vid]]) &&
            !is.null(compressed_templates_accessor[[veg]][[vid]][[grid_type]])) {
          vec <- compressed_templates_accessor[[veg]][[vid]][[grid_type]]
          if (!is.null(vec) && length(vec) > 0) {
            E_cols[[veg]] <- vec
          }
        }
      }
    }
    
    if (length(E_cols) >= 2) {
      E <- do.call(cbind, E_cols)
      
      # Check condition number
      cond_E <- tryCatch(kappa(t(E) %*% E), error = function(e) NA)
      cat(sprintf("  Endmember matrix condition number: %.2e\n", cond_E))
      
      # Check pairwise correlations
      E_cor <- cor(E)
      max_cor <- max(abs(E_cor[upper.tri(E_cor)]), na.rm = TRUE)
      cat(sprintf("  Max abs correlation between endmembers: %.4f\n", max_cor))
      
      # Print correlation matrix
      cat("  Endmember correlation matrix:\n")
      print(round(E_cor, 3))
      
      if (!is.na(cond_E) && cond_E > 1e10) {
        warning("CRITICAL: Endmember matrix is severely ill-conditioned for unmixing!")
      }
      if (!is.na(max_cor) && max_cor > 0.95) {
        warning("CRITICAL: Some endmembers are nearly identical (correlation > 0.95)!")
      }
    } else {
      cat("  Not enough endmembers for multicollinearity check\n")
    }
  }
  cat("==============================================\n\n")
  # ========== END MULTICOLLINEARITY CHECK ==========

  # Optionally precompute whitening for common combinations (3-way if available)
  if (length(mesma_lib) >= 3) {
    combos <- utils::combn(names(mesma_lib), min(3, length(mesma_lib)), simplify = FALSE)
    .WHITENING_DB_ACCESSOR <- precompute_whitening(mesma_lib, .COMPRESSED_TEMPLATES_ACCESSOR, combos)
    assign(".WHITENING_DB_ACCESSOR", .WHITENING_DB_ACCESSOR, envir = globalenv())
  }

  # Report library statistics
  for (veg in names(mesma_lib)) {
    n_variants <- length(mesma_lib[[veg]])
    cat(sprintf("%s: %d variants\n", veg, n_variants))
    for (var in mesma_lib[[veg]]) {
      cat(sprintf("  %s: %d samples\n", var$variant_id, var$n_samples))
    }
  }

  # Spatial bootstrap for library stability analysis (resample training traces)
  stability_results <- NULL
  # if (isTRUE(ENABLE_UNCERTAINTY) && !isTRUE(TESTING_MODE)) {


  # Save a MESMA model cache so the trained model (library + templates)
  # can be re-used later for testing/inference without re-training.
  save_mesma_cache <- function(cache_dir = file.path(OUT_DIR, "mesma_cache")) {
    cat("\n=== SAVING MESMA MODEL CACHE ===\n")
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    # 1. Core library components
    core <- list(
      lib = lib,
      mesma_lib = mesma_lib,
      raw_lib_templates = raw_lib_templates,
      veg_counts = veg_counts,
      avail = avail,
      ALLOWED_VEG = ALLOWED_VEG,
      BAND_SCALE = if (exists("BAND_SCALE")) BAND_SCALE else NULL,
      COMPRESSED_STAGE1_LIB = if (exists("COMPRESSED_STAGE1_LIB")) COMPRESSED_STAGE1_LIB else NULL,
      LDA_FEATURE_WEIGHTS = if (exists("LDA_FEATURE_WEIGHTS")) LDA_FEATURE_WEIGHTS else NULL
    )
    saveRDS(core, file = file.path(cache_dir, "mesma_library.rds"))

    # 2. Raw templates (365 x K index matrices per variant)
    raw_templates <- list(
      raw_lib_templates = if (exists("raw_lib_templates")) raw_lib_templates else NULL
    )
    saveRDS(raw_templates, file = file.path(cache_dir, "raw_templates.rds"))

    # 3. Compressed templates (precomputed feature vectors) accessed via accessor
    if (exists(".COMPRESSED_TEMPLATES_ACCESSOR", envir = globalenv())) {
      template_data <- list()
      for (veg in names(mesma_lib)) {
        template_data[[veg]] <- list()
        for (variant in mesma_lib[[veg]]) {
          for (grid_type in c("sparse", "full", "dense")) {
            key <- paste(veg, variant$variant_id, grid_type, sep = "|")
            vec <- tryCatch(.COMPRESSED_TEMPLATES_ACCESSOR[[veg]][[variant$variant_id]][[grid_type]], error = function(e) stop(sprintf("Compressed templates accessor failed for %s|%s|%s: %s", veg, variant$variant_id, grid_type, e$message)))
            if (!is.null(vec)) template_data[[veg]][[key]] <- vec
          }
        }
      }
      saveRDS(template_data, file = file.path(cache_dir, "compressed_templates.rds"))
    }

    # 4. Configuration parameters
    cfg <- list(
      TEMPORAL_BUDGET = if (exists("TEMPORAL_BUDGET")) TEMPORAL_BUDGET else NULL,
      TOPK_VARIANTS = if (exists("TOPK_VARIANTS")) TOPK_VARIANTS else NULL,
      ENABLE_PHASE_ALIGNMENT = if (exists("ENABLE_PHASE_ALIGNMENT")) ENABLE_PHASE_ALIGNMENT else NULL,
      REFERENCE_PHASE_MARKERS = if (exists("REFERENCE_PHASE_MARKERS")) REFERENCE_PHASE_MARKERS else NULL,
      ENABLE_MULTISCALE = if (exists("ENABLE_MULTISCALE")) ENABLE_MULTISCALE else NULL,
      MULTISCALE_WINDOWS = if (exists("MULTISCALE_WINDOWS")) MULTISCALE_WINDOWS else NULL,
      ENABLE_QP_SOLVER = if (exists("ENABLE_QP_SOLVER")) ENABLE_QP_SOLVER else NULL,
      ENABLE_DIAGNOSTICS = if (exists("ENABLE_DIAGNOSTICS")) ENABLE_DIAGNOSTICS else NULL,
      ENABLE_UNCERTAINTY = if (exists("ENABLE_UNCERTAINTY")) ENABLE_UNCERTAINTY else NULL,
      DEBUG_UNCERTAINTY = if (exists("DEBUG_UNCERTAINTY")) DEBUG_UNCERTAINTY else NULL,
      BOOTSTRAP_B = if (exists("BOOTSTRAP_B")) BOOTSTRAP_B else NULL,
      MAX_VEG_COMPONENTS = if (exists("MAX_VEG_COMPONENTS")) MAX_VEG_COMPONENTS else NULL,
      MIN_IDX_PRESENCE = if (exists("MIN_IDX_PRESENCE")) MIN_IDX_PRESENCE else NULL,
      EPS_SIGMA = if (exists("EPS_SIGMA")) EPS_SIGMA else NULL,
      LOWER_BND = if (exists("LOWER_BND")) LOWER_BND else NULL,
      USE_INDICES_MIN = if (exists("USE_INDICES_MIN")) USE_INDICES_MIN else NULL,
      MIN_INDEX_SD = if (exists("MIN_INDEX_SD")) MIN_INDEX_SD else NULL
    )
    saveRDS(cfg, file = file.path(cache_dir, "config_params.rds"))

    # 5. Training metadata
    meta <- list(
      training_years = if (exists("TRAIN_YEARS")) TRAIN_YEARS else NULL,
      training_date_range = if (exists("df_train")) range(df_train$date, na.rm = TRUE) else NULL,
      n_training_samples = if (exists("df_train")) nrow(df_train) else NULL,
      n_locations_trained = if (exists("df_train")) length(unique(df_train$location_id)) else NULL,
      indices_used = avail,
      vegetation_types = names(mesma_lib),
      model_creation_time = Sys.time()
    )
    saveRDS(meta, file = file.path(cache_dir, "training_metadata.rds"))

    # 6. Manifest with checksums
    manifest <- list(
      version = "1.0",
      created = Sys.time(),
      files = c("mesma_library.rds", "raw_templates.rds", "compressed_templates.rds", "config_params.rds", "training_metadata.rds"),
      checksum = list()
    )
    for (f in manifest$files) {
      fpath <- file.path(cache_dir, f)
      if (file.exists(fpath)) manifest$checksum[[f]] <- tryCatch(tools::md5sum(fpath), error = function(e) stop(sprintf("Failed to compute checksum for %s: %s", fpath, e$message)))
    }
    saveRDS(manifest, file = file.path(cache_dir, "manifest.rds"))

    cat(sprintf("Model cache saved to: %s\n", cache_dir))
    cat(sprintf("Total cache size: %.2f MB\n", sum(file.info(list.files(cache_dir, full.names = TRUE))$size) / 1024^2))
    invisible(cache_dir)
  }

  # Persist cache now so the trained model can be used later without re-training
  tryCatch({
    cache_dir <- save_mesma_cache()
  }, error = function(e) {
    cat(sprintf("Failed to write MESMA cache: %s\n", e$message))
  })

  # Prepare inference/test tasks — always use the full dataset (all years)
  cat("Constructing task list from inference dataset (all years)...\n")
  df_tasks <- df_full

# --- INFERENCE DATA LOADING ---
if (!isTRUE(SKIP_INFERENCE)) {
  df_inf <- NULL
  if (exists("INFERENCE_CSV")) {
    cat(sprintf("Checking inference file at: %s\n", INFERENCE_CSV))
    if (file.exists(INFERENCE_CSV)) {
      cat(sprintf("Loading inference data from %s...\n", INFERENCE_CSV))
      
      # Check extension to decide read method
      if (grepl("\\.csv$", INFERENCE_CSV, ignore.case = TRUE)) {
         df_inf <- tryCatch(read.csv(INFERENCE_CSV), error = function(e) {
            cat(sprintf("[WARNING] Error reading inference CSV: %s\n", e$message))
            NULL
         })
      } else {
         df_inf <- tryCatch(openxlsx::read.xlsx(INFERENCE_CSV), error = function(e) {
            cat(sprintf("[WARNING] Error reading inference XLSX: %s\n", e$message))
            NULL
         })
      }
      
      if (!is.null(df_inf)) {
         cat(sprintf("Loaded %d rows from inference file.\n", nrow(df_inf)))
         if ("location_id" %in% names(df_inf)) {
            cat(sprintf("Found %d unique inference location IDs.\n", length(unique(df_inf$location_id))))
         } else {
            cat("WARNING: 'location_id' column missing from inference file.\n")
         }
      }
    } else {
      cat(sprintf("Inference file not found at: %s\n", INFERENCE_CSV))
    }
  } else {
    cat("INFERENCE_CSV variable not defined.\n")
  }
    
  if (!is.null(df_inf) && nrow(df_inf) > 0) {
      # Normalize band column names to our canonical lower-case names (accept Blue/NIR/SWIR uppercase variants)
      df_inf <- normalize_band_names(df_inf)
    # Normalize columns
    # If the first column is unnamed or empty in Excel (row names), read.xlsx might name it "...1"
    if ("...1" %in% names(df_inf) && !"location_id" %in% names(df_inf)) {
       # Check if it looks like location_id
       if (is.character(df_inf$...1) || is.numeric(df_inf$...1)) {
         names(df_inf)[names(df_inf) == "...1"] <- "location_id"
       }
    }
    
    # Handle date
    if ("prediction_date" %in% names(df_inf)) {
      df_inf$date <- as.Date(df_inf$prediction_date)
    } else if ("date" %in% names(df_inf)) {
      df_inf$date <- as.Date(df_inf$date)
    } else {
      # Try to find a date-like column
      # Look for columns that parse as date
      for (col in names(df_inf)) {
         if (inherits(df_inf[[col]], "Date")) {
           df_inf$date <- df_inf[[col]]
           break
         }
         # Check string format YYYY-MM-DD
         if (is.character(df_inf[[col]]) && all(grepl("^\\d{4}-\\d{2}-\\d{2}", na.omit(df_inf[[col]][1:min(10, nrow(df_inf))])))) {
           df_inf$date <- as.Date(df_inf[[col]])
           break
         }
      }
    }
    
    if ("location_id" %in% names(df_inf) && "date" %in% names(df_inf)) {
      df_inf$location_id <- as.character(df_inf$location_id)
      
      # Deduplicate inference observations based on location_id and date
      original_rows <- nrow(df_inf)
      df_inf <- df_inf %>% distinct(location_id, date, .keep_all = TRUE)
      if (nrow(df_inf) < original_rows) {
        cat(sprintf("Inference data deduplicated: %d rows remaining from %d original rows.\n", nrow(df_inf), original_rows))
      }
      
      # If we have raw bands, compute missing indices from them
      if (length(intersect(RAW_BANDS, names(df_inf))) >= 2) {
        before_cols <- names(df_inf)
        df_inf <- compute_indices_from_bands(df_inf)
        new_cols <- setdiff(names(df_inf), before_cols)
        if (length(new_cols) > 0) cat(sprintf("[NOTICE] Computed indices from raw bands in inference data: %s\n", paste(new_cols, collapse=", ")))
      }
      # Ensure indices are present
      missing_idx <- setdiff(avail, names(df_inf))
      if (length(missing_idx) > 0) {
        cat(sprintf("[WARNING] Inference data missing indices: %s. Filling with NA (will likely fail unmixing).\n", paste(missing_idx, collapse=", ")))
        for (col in missing_idx) df_inf[[col]] <- NA_real_
      }
      
      # Add missing metadata columns
      if (!"Veg" %in% names(df_inf)) df_inf$Veg <- NA_character_
      if (!"year" %in% names(df_inf)) df_inf$year <- lubridate::year(df_inf$date)
      if (!"doy" %in% names(df_inf)) df_inf$doy <- lubridate::yday(df_inf$date)
      if (!"PPI" %in% names(df_inf)) df_inf$PPI <- NA_real_
      if (!"zenith.angle" %in% names(df_inf)) df_inf$zenith.angle <- NA_real_
      if (!"DVI_max" %in% names(df_inf)) df_inf$DVI_max <- NA_real_
      
      # Ensure prediction_date is Date type
      if ("prediction_date" %in% names(df_inf)) {
        df_inf$prediction_date <- as.Date(df_inf$prediction_date)
      }
      # Ensure reference_date is Date type
      if ("reference_date" %in% names(df_inf)) {
        df_inf$reference_date <- as.Date(df_inf$reference_date)
      }

      # Recompute PPI from DVI if possible
      if ("DVI" %in% names(df_inf)) {
        df_inf <- recompute_ppi_from_scaled_dvi(df_inf)
      }
      
      # Filter to 1985-2025
      df_inf <- df_inf %>% filter(year >= 1985 & year <= 2025)
      
      # Ensure numeric types for indices
      for (col in avail) {
        if (col %in% names(df_inf)) df_inf[[col]] <- as.numeric(df_inf[[col]])
      }

      # Filter inference data to 2020-2025
      df_inf <- df_inf %>% filter(year >= 2020 & year <= 2025)

      # Create separate df_tasks for inference
      df_tasks_inference <- df_inf
      cat(sprintf("Created separate inference task list with %d rows from %d locations.\n", nrow(df_tasks_inference), length(unique(df_tasks_inference$location_id))))
      
      # Deduplicate within inference data
      before_dedup <- nrow(df_tasks_inference)
      df_tasks_inference <- df_tasks_inference %>% dplyr::distinct(location_id, date, .keep_all = TRUE)
      after_dedup <- nrow(df_tasks_inference)
      if (before_dedup > after_dedup) {
        cat(sprintf("Removed %d duplicate rows from inference data.\n", before_dedup - after_dedup))
      }

      # Report inference dataset in terms of unique location-year pairs and compare with training
      if (!"year" %in% names(df_tasks_inference) && "date" %in% names(df_tasks_inference)) df_tasks_inference$year <- lubridate::year(df_tasks_inference$date)
      n_infer_loc_years <- nrow(unique(df_tasks_inference[c("location_id", "year")]))
      if (exists("df_train") && !is.null(df_train) && nrow(df_train) > 0) {
        if (!"year" %in% names(df_train) && "date" %in% names(df_train)) df_train$year <- lubridate::year(df_train$date)
        n_train_loc_years <- nrow(unique(df_train[c("location_id", "year")]))
      } else {
        n_train_loc_years <- 0
      }
      cat(sprintf("(NOTICE) Inference dataset location-years: %d\n", n_infer_loc_years))
      cat(sprintf("(NOTICE) Training dataset location-years: %d\n", n_train_loc_years))
      if (n_train_loc_years > 0 && n_infer_loc_years == n_train_loc_years) {
        stop(sprintf("ERROR: Training and inference datasets appear to have the same number of location-years (%d). This may indicate you passed the same data for training and inference — aborting to avoid accidental overlap.", n_train_loc_years))
      }

      # Keep df_tasks as training tasks - do not overwrite
      cat("Keeping df_tasks as training data for separate processing.\n")
    } else {
      cat("[WARNING] Inference data missing 'location_id' or 'date' column. Skipping.\n")
      cat(sprintf("Columns found: %s\n", paste(names(df_inf), collapse=", ")))
    }
  }
  # Capture inference location IDs for separate output
  inference_location_ids <- if (!is.null(df_inf) && nrow(df_inf) > 0 && "location_id" %in% names(df_inf)) unique(df_inf$location_id) else character(0)
  # No limit on inference locations
} else {
  cat("Skipping inference data loading (SKIP_INFERENCE = TRUE).\n")
  df_inf <- NULL
  inference_location_ids <- character(0)
}

  if (exists("gpts_map")) {
    if ("location_id" %in% names(df_tasks) && "location_id" %in% names(gpts_map)) {
      if (!is.character(df_tasks$location_id)) {
        df_tasks$location_id <- as.character(df_tasks$location_id)
        cat("[NOTICE] Coerced df_tasks$location_id to character to match gpts_map for joining.\n")
      }
      if (!is.character(gpts_map$location_id)) {
        gpts_map$location_id <- as.character(gpts_map$location_id)
        cat("[NOTICE] Coerced gpts_map$location_id to character for joining with df_tasks.\n")
      }
      # warn if formats look mismatched
      s1 <- unique(na.omit(df_tasks$location_id))
      s2 <- unique(na.omit(gpts_map$location_id))
      if (length(s1) && length(s2) && all(grepl("^[0-9]+$", s1)) && any(grepl("^L_", s2))) {
        cat("[WARNING] df_tasks$location_id looks numeric while gpts_map$location_id looks like 'L_x_y' strings — matching will likely fail.\n")
      }
      # Use suffix to avoid overwriting existing columns and prefer existing task values
      df_tasks <- dplyr::left_join(df_tasks, gpts_map, by = "location_id", suffix = c("", ".geo"))
      if ("Veg.geo" %in% names(df_tasks)) {
        df_tasks$Veg <- ifelse(is.na(df_tasks$Veg) | df_tasks$Veg == "", df_tasks$Veg.geo, df_tasks$Veg)
        df_tasks$Veg.geo <- NULL
      }
      if ("no soil.geo" %in% names(df_tasks)) {
        if ("no soil" %in% names(df_tasks)) {
          df_tasks$`no soil` <- ifelse(is.na(df_tasks$`no soil`), df_tasks$`no soil.geo`, df_tasks$`no soil`)
        } else {
          df_tasks$`no soil` <- df_tasks$`no soil.geo`
        }
        df_tasks$`no soil.geo` <- NULL
      }
      # Normalize possible alternate column names (e.g. no.soil, no_soil)
      df_tasks <- normalize_no_soil_col(df_tasks)
    } else {
      # perform a straight join attempt (this will error if keys missing)
      df_tasks <- dplyr::left_join(df_tasks, gpts_map, by = "location_id", suffix = c("", ".geo"))
      if ("Veg.geo" %in% names(df_tasks)) {
        df_tasks$Veg <- ifelse(is.na(df_tasks$Veg) | df_tasks$Veg == "", df_tasks$Veg.geo, df_tasks$Veg)
        df_tasks$Veg.geo <- NULL
      }
      if ("no soil.geo" %in% names(df_tasks)) {
        if ("no soil" %in% names(df_tasks)) {
          df_tasks$`no soil` <- ifelse(is.na(df_tasks$`no soil`), df_tasks$`no soil.geo`, df_tasks$`no soil`)
        } else {
          df_tasks$`no soil` <- df_tasks$`no soil.geo`
        }
        df_tasks$`no soil.geo` <- NULL
      }
      df_tasks <- normalize_no_soil_col(df_tasks)
    }
  }

  if ("Veg" %in% names(df_tasks)) df_tasks$Veg <- tolower(as.character(df_tasks$Veg))
  if ("date" %in% names(df_tasks)) {
    df_tasks$date <- as.Date(df_tasks$date)
    df_tasks$year <- as.integer(lubridate::year(df_tasks$date))
    df_tasks$doy <- lubridate::yday(df_tasks$date)
    df_tasks$doy[df_tasks$doy < 1 | df_tasks$doy > 366] <- NA_integer_
  }

  test_loc_years <- df_tasks %>%
    dplyr::filter(!is.na(.data$location_id) & trimws(.data$location_id) != "" & !is.na(.data$year) & .data$year > 0) %>%
    dplyr::distinct(.data$location_id, .data$year)

  test_loc_years$location_id <- trimws(as.character(test_loc_years$location_id))

  cat(sprintf(
    "Final test_loc_years: %d rows from %d unique locations\n",
    nrow(test_loc_years), length(unique(test_loc_years$location_id))
  ))

  n_loc_years <- nrow(test_loc_years)
  task_list <- lapply(seq_len(n_loc_years), function(i) {
    list(loc = test_loc_years$location_id[i], yr = test_loc_years$year[i])
  })

  # Before task processing, validate one location
cat("\n=== PRE-FLIGHT CHECK ===\n")
test_task <- task_list[[1]]
test_dly <- df_tasks[df_tasks$location_id == test_task$loc, , drop = FALSE]
test_dly_year <- test_dly[lubridate::year(test_dly$date) == test_task$yr, , drop = FALSE]

cat(sprintf("Test location: %s, year: %d\n", test_task$loc, test_task$yr))
cat(sprintf("Rows in test_dly_year: %d\n", nrow(test_dly_year)))
cat(sprintf("Columns available: %s\n", paste(names(test_dly_year), collapse=", ")))
cat(sprintf("Var14 columns: %s\n", paste(grep("_var14$", names(test_dly_year), value=TRUE), collapse=", ")))
  cat("======================\n\n")

  # ===== ADD COMPREHENSIVE DATA DISTRIBUTION CHECK =====
  cat("\n=== DATA DISTRIBUTION ANALYSIS ===\n")
  if (exists("df_tasks") && nrow(df_tasks) > 0) {
    # Check overall distribution
    cat(sprintf("Total locations in df_tasks: %d\n", length(unique(df_tasks$location_id))))
    cat(sprintf("Total location-years: %d\n", nrow(test_loc_years)))

    sample_sizes <- df_tasks %>%
      dplyr::group_by(location_id, year) %>%
      dplyr::summarize(n_obs = dplyr::n(), .groups = "drop")

    if (nrow(sample_sizes) > 0) {
      # Ensure numeric (handle possible factors)
      sample_sizes$n_obs <- as.numeric(sample_sizes$n_obs)
      cat("\nObservations per location-year distribution:\n")
      cat(sprintf("  Min:    %d\n", as.integer(min(sample_sizes$n_obs, na.rm = TRUE))))
      cat(sprintf("  Q1:     %d\n", as.integer(stats::quantile(sample_sizes$n_obs, 0.25, na.rm = TRUE))))
      cat(sprintf("  Median: %d\n", as.integer(median(sample_sizes$n_obs, na.rm = TRUE))))
      cat(sprintf("  Q3:     %d\n", as.integer(stats::quantile(sample_sizes$n_obs, 0.75, na.rm = TRUE))))
      cat(sprintf("  Max:    %d\n", as.integer(max(sample_sizes$n_obs, na.rm = TRUE))))
      cat(sprintf("  Mean:   %.1f\n", mean(sample_sizes$n_obs, na.rm = TRUE)))

      # Identify problematic location-years
      problem_threshold <- MIN_OBS_PER_LOC_YEAR  # Minimum observations for meaningful MESMA
      n_problem <- sum(sample_sizes$n_obs < problem_threshold, na.rm = TRUE)
      cat(sprintf("\nLocation-years with < %d observations: %d (%.1f%%)\n",
                  problem_threshold, n_problem, 100 * n_problem / nrow(sample_sizes)))
      if (n_problem > 0) {
        cat("\nSample of problematic location-years:\n")
        problem_cases <- sample_sizes[sample_sizes$n_obs < problem_threshold, , drop = FALSE]
        print(utils::head(problem_cases, 10))
      }

      # Check DOY coverage
      if ("doy" %in% names(df_tasks)) {
        doy_coverage <- df_tasks %>%
          dplyr::group_by(location_id, year) %>%
          dplyr::summarize(
            n_unique_doys = length(unique(doy[!is.na(doy)])),
            doy_span = if (any(!is.na(doy))) as.integer(max(doy, na.rm = TRUE) - min(doy, na.rm = TRUE)) else 0,
            .groups = "drop"
          )
        cat("\nDOY coverage per location-year:\n")
        cat(sprintf("  Median unique DOYs: %d\n", as.integer(median(doy_coverage$n_unique_doys, na.rm = TRUE))))
        cat(sprintf("  Median DOY span:    %d days\n", as.integer(median(doy_coverage$doy_span, na.rm = TRUE))))
        n_poor_coverage <- sum(doy_coverage$n_unique_doys < 30, na.rm = TRUE)
        cat(sprintf("  Location-years with < 30 unique DOYs: %d (%.1f%%)\n",
                    n_poor_coverage, 100 * n_poor_coverage / nrow(doy_coverage)))
      }
    } else {
      cat("No location-year sample sizes available to summarize.\n")
    }

    # Check test_dly_year specifically
    cat(sprintf("\n=== SPECIFIC TEST LOCATION ===\n"))
    cat(sprintf("Location %s, Year %d:\n", test_task$loc, test_task$yr))
    cat(sprintf("  Total observations: %d\n", nrow(test_dly_year)))
    if (nrow(test_dly_year) > 0) {
      if ("date" %in% names(test_dly_year)) {
        cat(sprintf("  Date range: %s to %s\n", min(test_dly_year$date, na.rm = TRUE), max(test_dly_year$date, na.rm = TRUE)))
      }
      cat(sprintf("  Unique DOYs: %d\n", length(unique(test_dly_year$doy))))
      if (any(!is.na(test_dly_year$doy))) {
        cat(sprintf("  DOY range: %d to %d\n", min(test_dly_year$doy, na.rm = TRUE), max(test_dly_year$doy, na.rm = TRUE)))
      }
    }
    cat("======================\n\n")
  } else {
    cat("df_tasks not present or empty — skipping data distribution analysis.\n")
  }

  # Filter out location-years with insufficient observations (inference tasks)
  # DISABLED: df_tasks may be a task list (1 row per task) and not contain full observations.
  # Filtering for sufficient data (e.g. >5 unique DOYs) is handled inside the processing loop.
  # if (exists("df_tasks") && nrow(df_tasks) > 0) {
  #   # Recompute sample sizes reliably (in case diagnostics ran in different scope)
  #   sample_sizes <- df_tasks %>%
  #     dplyr::mutate(location_id = trimws(as.character(location_id))) %>%
  #     dplyr::group_by(location_id, year) %>%
  #     dplyr::summarize(n_obs = dplyr::n(), .groups = "drop")
  #   keep_pairs <- sample_sizes[sample_sizes$n_obs >= MIN_OBS_PER_LOC_YEAR, , drop = FALSE]
  #   # Coerce join columns to consistent types in case they are factors
  #   keep_pairs$location_id <- trimws(as.character(keep_pairs$location_id))
  #   keep_pairs$year <- as.integer(keep_pairs$year)
  #   test_loc_years$location_id <- trimws(as.character(test_loc_years$location_id))
  #   test_loc_years$year <- as.integer(test_loc_years$year)
  #   n_before <- nrow(test_loc_years)
  #   # Perform an inner join to keep only the pairs that meet the threshold
  #   test_loc_years <- dplyr::inner_join(test_loc_years, keep_pairs[, c('location_id', 'year')], by = c('location_id','year'))
  #   n_after <- nrow(test_loc_years)
  #   cat(sprintf("\nFiltered tasks: kept %d/%d location-years with >= %d observations\n", n_after, n_before, MIN_OBS_PER_LOC_YEAR))
  #   if (n_after < n_before) {
  #     cat("Sample of excluded location-years:\n")
  #     excluded <- setdiff(sample_sizes[, c('location_id','year')], keep_pairs[, c('location_id','year')])
  #     print(utils::head(excluded, 10))
  #   }
  #   if (n_after == 0) {
  #     stop(sprintf("No testing tasks remain after filtering location-years with < %d observations; aborting.", MIN_OBS_PER_LOC_YEAR))
  #   }
  #
  #   # Add debugging before task list creation
  #   cat(sprintf("DEBUG: test_loc_years has %d rows\n", nrow(test_loc_years)))
  #   cat(sprintf("DEBUG: NA location_ids in test_loc_years: %d\n", sum(is.na(test_loc_years$location_id))))
  #   cat(sprintf("DEBUG: Sample location_ids: %s\n", paste(head(test_loc_years$location_id, 10), collapse=", ")))
  #
  #   # Rebuild the task_list to reflect the filtered location-years
  #   n_loc_years <- nrow(test_loc_years)
  #   task_list <- lapply(seq_len(n_loc_years), function(i) {
  #     list(loc = test_loc_years$location_id[i], yr = test_loc_years$year[i])
  #   })
  # }

  if (length(task_list) == 0) {
    stop("No testing tasks found")
  }

  # ============================================================================
  # UNIFIED COMPRESSION PIPELINE FUNCTIONS
  # ============================================================================

  #' Build temporal matrix from daily data
  #' @param dly_year Data frame with daily observations including 'doy' column
  #' @param avail Vector of available indices
  #' @return Matrix with 365 rows (DOY 1-365) and columns for each index in avail
  build_temporal_matrix <- function(dly_year, avail) {
    if (nrow(dly_year) == 0 || length(avail) == 0) return(NULL)
    
    # Ensure doy column exists

    if (!"doy" %in% names(dly_year)) {
      if ("date" %in% names(dly_year)) {
        dly_year$doy <- lubridate::yday(dly_year$date)
      } else {
        return(NULL)
      }
    }

    # Initialize matrix: 365 rows (DOY 1-365), columns for avail indices
    temporal_matrix <- matrix(NA_real_, nrow = 365, ncol = length(avail))
    colnames(temporal_matrix) <- avail

    # Get indices that exist in the data
    avail_present <- intersect(avail, names(dly_year))
    if (length(avail_present) == 0) return(NULL)

    # Fill matrix with values from dly_year based on doy
    for (i in seq_len(nrow(dly_year))) {
      doy <- dly_year$doy[i]
      if (!is.na(doy) && doy >= 1 && doy <= 365) {
        for (idx in avail_present) {
          val <- dly_year[[idx]][i]
          if (is.finite(val)) {
            # If multiple observations for same DOY, keep the one closest to existing mean
            existing <- temporal_matrix[doy, idx]
            if (is.na(existing)) {
              temporal_matrix[doy, idx] <- val
            } else {
              # Average with existing value
              temporal_matrix[doy, idx] <- (existing + val) / 2
            }
          }
        }
      }
    }

    temporal_matrix
  }

  #' Compress temporal matrix using adaptive grids and information content weighting
  #' @param temporal_matrix Matrix with rows as time points, columns as indices
  #' @param temporal_budget Maximum number of time points to retain
  #' @param avail Vector of available indices (column names or indices)
  #' @return List with compressed matrix, grid information, and metadata
  #' Build raw index matrix from temporal data
  #' @param temporal_matrix Matrix with spectral indices over time
  #' @param avail Vector of available index names
  #' @return Matrix ready for information content computation
  build_raw_index_matrix <- function(temporal_matrix, avail) {
    # Use whatever indices are available - don't restrict to a hardcoded list
    if (is.null(temporal_matrix) || length(avail) == 0) return(NULL)
    
    # If temporal_matrix has column names, use them
    if (!is.null(colnames(temporal_matrix))) {
      available_indices <- intersect(avail, colnames(temporal_matrix))
      if (length(available_indices) > 0) {
        raw_matrix <- temporal_matrix[, available_indices, drop = FALSE]
      } else {
        # No matching names - assume columns correspond to avail order
        if (ncol(temporal_matrix) >= length(avail)) {
          raw_matrix <- temporal_matrix[, seq_along(avail), drop = FALSE]
          colnames(raw_matrix) <- avail
        } else {
          raw_matrix <- temporal_matrix
        }
      }
    } else {
      # No column names - use as-is
      raw_matrix <- temporal_matrix
      if (ncol(raw_matrix) == length(avail)) {
        colnames(raw_matrix) <- avail
      }
    }
    
    # Ensure finite values
    raw_matrix[!is.finite(raw_matrix)] <- NA
    
    raw_matrix
  }

  #' Unified trace compression function - wrapper around original compress_trace
  #' @param dly_year Data frame with daily observations
  #' @param avail Available indices
  #' @param temporal_budget Maximum temporal points
  #' @return Compressed trace result
  compress_trace_unified <- function(dly_year, avail, temporal_budget) {
    compress_trace(dly_year, avail, temporal_budget)
  }

  #' Unified Stage 1 library compression
  #' @param stage1_lib Raw Stage 1 endmember library
  #' @param avail Available indices
  #' @param temporal_budget Temporal compression budget
  #' @return Compressed Stage 1 library

# Dummy compress_temporal_matrix for compatibility (no compression)
compress_temporal_matrix <- function(data, temporal_budget, avail) {
  # No compression - return data as-is
  list(
    compressed_matrix = data,
    grid = 1:365,
    grid_type = "full"
  )
}
  compress_stage1_lib_unified <- function(stage1_lib, avail, temporal_budget) {
    if (is.null(stage1_lib)) return(NULL)

    compressed_lib <- list()

    for (endmember_name in names(stage1_lib)) {
      endmember_data <- stage1_lib[[endmember_name]]

      # Compress each endmember's temporal profile
      compressed <- compress_temporal_matrix(endmember_data, temporal_budget, avail)

      if (!is.null(compressed)) {
        compressed_lib[[endmember_name]] <- list(
          compressed = compressed$compressed_matrix,
          grid = compressed$grid,
          grid_type = compressed$grid_type,
          veg_frac = compressed$veg_frac
        )
      }
    }

    compressed_lib
  }

  #' Precompute compressed templates using unified pipeline
  #' @param mesma_lib MESMA library
  #' @param avail Available indices
  #' @param temporal_budget Temporal budget
  #' @return Precomputed compressed templates
  precompute_compressed_templates_unified <- function(mesma_lib, avail, temporal_budget) {
    if (is.null(mesma_lib)) return(NULL)

    compressed_templates <- list()

    for (veg_type in names(mesma_lib)) {
      veg_variants <- mesma_lib[[veg_type]]

      compressed_variants <- list()

      for (variant_name in names(veg_variants)) {
        variant_data <- veg_variants[[variant_name]]

        # Compress variant
        compressed <- compress_temporal_matrix(variant_data, temporal_budget, avail)

        if (!is.null(compressed)) {
          compressed_variants[[variant_name]] <- list(
            compressed = compressed$compressed_matrix,
            grid = compressed$grid,
            grid_type = compressed$grid_type,
            veg_frac = compressed$veg_frac
          )
        }
      }

      compressed_templates[[veg_type]] <- compressed_variants
    }

    compressed_templates
  }

  #' Unified trace reduction for all traces
  #' @param traces List of trace data frames
  #' @param avail Available indices
  #' @param temporal_budget Temporal budget
  #' @return List of compressed traces
  reduce_all_traces_unified <- function(traces, avail, temporal_budget) {
    if (is.null(traces) || length(traces) == 0) return(list())

    compressed_traces <- list()

    for (trace_name in names(traces)) {
      trace_data <- traces[[trace_name]]

      compressed <- compress_trace_unified(trace_data, avail, temporal_budget)

      if (!is.null(compressed)) {
        compressed_traces[[trace_name]] <- compressed
      }
    }

    compressed_traces
  }

  #' Build MESMA variants using unified compression
  #' @param lib Library data
  #' @param avail Available indices
  #' @param temporal_budget Temporal budget
  #' @return MESMA variants
  build_mesma_variants_unified <- function(lib, avail, temporal_budget) {
    if (is.null(lib)) return(NULL)

    mesma_variants <- list()

    for (veg_type in names(lib)) {
      veg_data <- lib[[veg_type]]

      # Compress using unified pipeline
      compressed <- compress_temporal_matrix(veg_data, temporal_budget, avail)

      if (!is.null(compressed)) {
        mesma_variants[[veg_type]] <- list(
          compressed = compressed$compressed_matrix,
          grid = compressed$grid,
          grid_type = compressed$grid_type,
          veg_frac = compressed$veg_frac
        )
      }
    }

    mesma_variants
  }

  #' Precompute templates from unified library
  #' @param unified_lib Unified compressed library
  #' @param veg_types Vegetation types to include
  #' @return Template accessor function
  precompute_templates_from_unified_lib <- function(unified_lib, veg_types) {
    if (is.null(unified_lib)) return(NULL)

    # Create accessor function for compressed templates
    function(veg_type, variant_id = NULL) {
      if (!veg_type %in% names(unified_lib)) return(NULL)

      veg_templates <- unified_lib[[veg_type]]

      if (is.null(variant_id)) {
        # Return all variants for this vegetation type
        return(veg_templates)
      } else {
        # Return specific variant
        return(veg_templates[[variant_id]])
      }
    }
  }

  # Define fit_one_task function for MESMA
  # OPTIMIZED: Accepts pre-filtered data frame (task_data) instead of looking up in global df_tasks
  fit_one_task <- function(task_data) {
    if (is.null(task_data) || nrow(task_data) == 0) return(NULL)
    
    # Extract metadata from first row
    loc <- as.character(task_data$location_id[1])
    yr <- as.integer(task_data$year[1])
    
    loc <- trimws(loc)
    
    # Early exit if loc is NA or empty
    if (is.na(loc) || loc == "") {
      return(NULL)
    }

    res_safe <- tryCatch(
      {
        # Use passed data directly
        dly_year <- task_data
        
        # Ensure DOY column exists and check sufficiency
        if (!"doy" %in% names(dly_year)) dly_year$doy <- lubridate::yday(dly_year$date)
        # Determine whether this task belongs to inference (allow lower DOY counts)
        min_doy_needed <- if (exists("inference_location_ids") && loc %in% inference_location_ids) MIN_UNIQUE_DOY_INFERENCE else MIN_UNIQUE_DOY_DEFAULT
        if (length(unique(dly_year$doy)) < min_doy_needed) {
          # For inference locations we do not skip tasks solely due to low observation count
          if (min_doy_needed > MIN_UNIQUE_DOY_INFERENCE) {
            return(NULL)
          } else {
            # DEBUG: still processing but with few DOYs
            if (exists("DEBUG") && isTRUE(DEBUG)) cat(sprintf("[DEBUG] Inference location %s/%d has only %d unique DOYs — still processing.\n", loc, yr, length(unique(dly_year$doy))))
          }
        }
        
        # Calculate Q10 and Q90 DVI
        # Note: We use the passed data (dly_year) for this. 
        # If Q10/Q90 was intended to be over ALL years for that location, we'd need the full history.
        # Assuming per-year stats or that task_data is sufficient.
        # If full history is needed, we'd need to pass it or accept the trade-off for performance.
        # Given the previous code: dly_train <- df_tasks[... & year == yr], it was per-year.
        dvi_vals <- dly_year$DVI[is.finite(dly_year$DVI)]
        q10_dvi <- if (length(dvi_vals) > 0) stats::quantile(dvi_vals, 0.10, na.rm = TRUE) else NA
        q90_dvi <- if (length(dvi_vals) > 0) stats::quantile(dvi_vals, 0.90, na.rm = TRUE) else NA

        # 1. Build Full Resolution Trace (No Compression)
        # Build raw 365 x K matrix
        raw_mat <- build_raw_365_matrix(dly_year, avail)
        if (is.null(raw_mat)) return(NULL)
        
        # Flatten to vector
        y <- as.numeric(raw_mat)
        if (length(y) == 0 || any(!is.finite(y))) return(NULL)
        
        # 2. STAGE 1 MESMA: Unmix vegetated fraction using FULL RESOLUTION endmembers (geometric)
        # Must use the same feature space as the observation (full resolution)
        if (!is.null(COMPRESSED_STAGE1_LIB) && length(y) > 0) {
          stage1_res <- geometric_stage1_unmix(y, COMPRESSED_STAGE1_LIB$barren, COMPRESSED_STAGE1_LIB$vegetation)
          veg_fraction_total <- stage1_res$veg_frac
        } else {
          # Fallback to raw unmixing (but this operates in different feature space!)
          veg_fraction_total <- unmix_vegetated_fraction(dly_year, STAGE1_LIB, avail)
        }
        
        if (!is.finite(veg_fraction_total)) {
          # cat(sprintf("[DEBUG fit_one_task] Stage 1 unmixing returned non-finite for %s/%d\n", loc, yr))
          veg_fraction_total <- 0.5  # Default to 50/50 if unmixing fails
        }
        veg_fraction_total <- pmax(0, pmin(1, veg_fraction_total))
        barren_fraction_total <- 1 - veg_fraction_total
        
        # Apply minimum fraction threshold (10%) like in stage 2
        MIN_FRACTION <- 0.10
        if (veg_fraction_total < MIN_FRACTION) {
          veg_fraction_total <- 0
          barren_fraction_total <- 1
        } else if (barren_fraction_total < MIN_FRACTION) {
          barren_fraction_total <- 0
          veg_fraction_total <- 1
        }
        
        if (!is.finite(veg_fraction_total)) return(NULL)
        # If veg_fraction_total is exactly zero, this is a pure barren observation — proceed
        if (veg_fraction_total <= 0) {
          veg_fraction_total <- 0
          barren_fraction_total <- 1
          # Prepare output with only barren fraction
          coef_df <- data.frame(
            location_id = loc,
            year = yr,
            Veg = "barren",
            coef = barren_fraction_total,
            rmse = NA_real_,
            stringsAsFactors = FALSE
          )
          coef_df$coef_025 <- NA_real_
          coef_df$coef_975 <- NA_real_
          coef_df$coef_sd <- NA_real_
          coef_df$interval <- NA_real_
          coef_df$inseparable_variant_flag <- FALSE
          coef_df$inseparable_variant_details <- NA_character_
          diag_df <- data.frame(location_id = loc, year = yr, vegetated_fraction = veg_fraction_total, barren_fraction = barren_fraction_total, stringsAsFactors = FALSE)
          return(list(coef_df = coef_df, variant_trajectory = NULL, diagnostics = diag_df, uncertainty = NULL, q10_dvi = q10_dvi, q90_dvi = q90_dvi, vegetated_fraction = veg_fraction_total, barren_fraction = barren_fraction_total, inseparable_flag = FALSE, inseparable_details = NA_character_))
        }

        veg_names <- names(lib)
        veg_kept <- veg_names[tolower(veg_names) %in% ALLOWED_VEG]

        # (Logging removed for speed in optimized loop)
        # cat(sprintf("veg_names: %s\n", paste(veg_names, collapse = ",")))

        if (length(veg_kept) < MAX_VEG_COMPONENTS && length(veg_counts) > 0) {
          global_order <- names(veg_counts)
          global_order <- intersect(global_order, veg_names)
          add <- setdiff(global_order, veg_kept)
          add <- add[tolower(add) != "barren"]
          if (length(add) > 0) {
            veg_kept <- c(veg_kept, add[seq_len(min(length(add), MAX_VEG_COMPONENTS - length(veg_kept)))])
          }
        }

        if (length(veg_kept) > MAX_VEG_COMPONENTS) {
          veg_kept <- veg_kept[seq_len(MAX_VEG_COMPONENTS)]
        }

        if (length(veg_kept) == 0) {
          return(dbg_return_null("no_veg_kept"))
        }

        # 3. STAGE 2 MESMA: Unmix vegetation using full resolution data
        # OPTIMIZED: Pass OPTIMIZED_LIBRARY for fast vectorized similarity search
        mesma_result <- unmix_stage2_compressed(
          veg_kept, veg_fraction_total, y, "full", 
          compressed_templates_accessor, mesma_lib, 
          topK = TOPK_VARIANTS, feature_weights = LDA_FEATURE_WEIGHTS,
          optimized_library = OPTIMIZED_LIBRARY
        )
        if (isTRUE(DEBUG_UNCERTAINTY)) {
          # cat(sprintf("[DEBUG] Location %s, Year %d: ENABLE_UNCERTAINTY = %s\n", loc, yr, as.character(ENABLE_UNCERTAINTY)))
          # cat(sprintf("[DEBUG] mesma_result$uncertainty is NULL: %s\n", as.character(is.null(mesma_result$uncertainty))))
          if (!is.null(mesma_result$uncertainty)) {
            # cat(sprintf("[DEBUG] coef_ci exists: %s\n", as.character(!is.null(mesma_result$uncertainty$coef_ci))))
            if (!is.null(mesma_result$uncertainty$coef_ci)) {
              # cat(sprintf("[DEBUG] coef_ci rows: %d\n", nrow(mesma_result$uncertainty$coef_ci)))
              print(mesma_result$uncertainty$coef_ci)
            }
          }
        }
        if (is.null(mesma_result)) return(dbg_return_null("mesma_pca_failed"))

        inseparable_flag <- FALSE
        inseparable_details <- character(0)

        coef_stage2 <- mesma_result$vegetation_proportions
        use_median <- !is.null(mesma_result$uncertainty$w_median)
        if (use_median) {
          coef_stage2 <- mesma_result$uncertainty$w_median
          is_nested <- !is.null(mesma_result$uncertainty$veg_frac_ci)
          if (is_nested) {
            # w_median is final coef
            coef_values <- coef_stage2[names(coef_stage2)]
          } else {
            # w_median is stage2 weights
            coef_values <- coef_stage2[names(coef_stage2)]  # Keep as proportions within vegetation
          }
        } else {
          coef_values <- as.numeric(coef_stage2)  # Keep as proportions within vegetation
        }
        
        # Multiply vegetation coefficients by vegetated fraction for correct response
        if (!(use_median && is_nested)) {
          coef_values <- coef_values * veg_fraction_total
        }
        
        coef_df <- data.frame(
          location_id = loc,
          year = yr,
          Veg = names(coef_stage2),
          coef = coef_values,
          rmse = mesma_result$rmse,
          stringsAsFactors = FALSE
        )

        inseparable_text <- if (inseparable_flag && length(inseparable_details) > 0) paste(unique(inseparable_details), collapse = "; ") else NA_character_

        # Initialize CI columns
        coef_df$coef_025 <- NA_real_
        coef_df$coef_975 <- NA_real_
        coef_df$coef_sd <- NA_real_
        coef_df$interval <- NA_real_
        coef_df$inseparable_variant_flag <- inseparable_flag
        coef_df$inseparable_variant_details <- inseparable_text

        # Populate CIs if uncertainty was computed
        if (!is.null(mesma_result$uncertainty) && !is.null(mesma_result$uncertainty$coef_ci)) {
          ci_tbl <- mesma_result$uncertainty$coef_ci
          if (isTRUE(DEBUG_UNCERTAINTY)) {
            # cat(sprintf("[DEBUG] CI table vegs: %s\n", paste(ci_tbl$Veg, collapse=", ")))
            # cat(sprintf("[DEBUG] coef_df vegs: %s\n", paste(coef_df$Veg, collapse=", ")))
          }
          
          # Determine if this is nested bootstrap (already includes veg_fraction)
          is_nested <- !is.null(mesma_result$uncertainty$veg_frac_ci)
          
          for (i in seq_len(nrow(coef_df))) {
            vname <- coef_df$Veg[i]
            row_ci <- ci_tbl[ci_tbl$Veg == vname, , drop = FALSE]
            
            if (isTRUE(DEBUG_UNCERTAINTY)) {
              # cat(sprintf("[DEBUG] Looking for %s: found %d matching rows\n", vname, nrow(row_ci)))
            }
            
            if (nrow(row_ci) >= 1) {  # Changed from == 1 to >= 1
              # Take first match if multiple (shouldn't happen but be safe)
              ci_025 <- row_ci$coef_025[1]
              ci_975 <- row_ci$coef_975[1]
              ci_sd <- if ("coef_sd" %in% names(row_ci)) row_ci$coef_sd[1] else NA_real_
              
              if (is.finite(ci_025) && is.finite(ci_975)) {
                if (is_nested) {
                  # Already includes veg_fraction
                  coef_df$coef_025[i] <- ci_025
                  coef_df$coef_975[i] <- ci_975
                  coef_df$coef_sd[i] <- ci_sd
                } else {
                  # Need to multiply by veg_fraction
                  coef_df$coef_025[i] <- veg_fraction_total * ci_025
                  coef_df$coef_975[i] <- veg_fraction_total * ci_975
                  coef_df$coef_sd[i] <- veg_fraction_total * ci_sd
                }
                # Calculate interval
                coef_df$interval[i] <- coef_df$coef_975[i] - coef_df$coef_025[i]
              } else if (isTRUE(DEBUG_UNCERTAINTY)) {
                # cat(sprintf("[DEBUG] Non-finite CIs for %s: ci_025=%s, ci_975=%s\n",
                #             vname, ci_025, ci_975))
              }
            } else if (isTRUE(DEBUG_UNCERTAINTY)) {
              cat(sprintf("[DEBUG] No CI rows found for vegetation: %s\n", vname))
            }
          }
        }
        variant_info_pca <- if (!is.null(mesma_result$chosen_variants)) {
          vi <- as.list(mesma_result$chosen_variants)
          names(vi) <- paste0(names(vi), "_variant")
          data.frame(location_id = loc, year = yr, vi, stringsAsFactors = FALSE, check.names = FALSE)
        } else NULL
        unc <- mesma_result$uncertainty
        diag_df <- if (!is.null(mesma_result$diagnostics)) {
          dd <- as.data.frame(mesma_result$diagnostics, stringsAsFactors = FALSE)
          dd$location_id <- loc; dd$year <- yr
          dd <- dd[, c("location_id", "year", setdiff(names(dd), c("location_id", "year"))) , drop = FALSE]
          dd
        } else {
          data.frame(location_id = loc, year = yr, stringsAsFactors = FALSE)
        }
        diag_df$vegetated_fraction <- veg_fraction_total
        diag_df$barren_fraction <- barren_fraction_total
        diag_df$inseparable_variant_flag <- inseparable_flag
        diag_df$inseparable_variant_details <- inseparable_text
        # If we have a Stage 1 bootstrap result, show its CI
        # Note: stage1_boot not defined, using unc instead
        if (!is.null(unc) && !is.null(unc$veg_frac_ci)) {
          diag_df$vegetated_frac_lo <- unc$veg_frac_ci[1]
          diag_df$vegetated_frac_hi <- unc$veg_frac_ci[2]
        }

        if (barren_fraction_total >= 0) {
          barren_row <- data.frame(
            location_id = loc,
            year = yr,
            Veg = "barren",
            coef = barren_fraction_total,
            rmse = mesma_result$rmse,
            coef_025 = NA_real_,
            coef_975 = NA_real_,
            coef_sd = NA_real_,
            interval = NA_real_,
            inseparable_variant_flag = FALSE,
            inseparable_variant_details = NA_character_,
            stringsAsFactors = FALSE
          )
          coef_df <- rbind(coef_df, barren_row)
        }
        if (!is.null(unc) && !is.null(unc$coef_ci) && is.data.frame(unc$coef_ci)) {
          is_nested_unc <- !is.null(unc$veg_frac_ci)
          if (!is_nested_unc) {
            # Stage 2-only uncertainty: coef_ci refers to stage2 weights, so multiply by overall vegetated fraction
            unc$coef_ci$coef_025 <- veg_fraction_total * unc$coef_ci$coef_025
            unc$coef_ci$coef_975 <- veg_fraction_total * unc$coef_ci$coef_975
            if ("coef_sd" %in% names(unc$coef_ci)) {
              unc$coef_ci$coef_sd <- veg_fraction_total * unc$coef_ci$coef_sd
            }
          }
          # If is_nested_unc then coef_ci already contains final coefficients
          # Ensure returned coef_ci contains location_id and year for later aggregation
          if (!"location_id" %in% names(unc$coef_ci)) unc$coef_ci$location_id <- loc
          if (!"year" %in% names(unc$coef_ci)) unc$coef_ci$year <- yr
          # Reorder columns to predictable order: location_id, year, Veg, coef_025, coef_975, coef_sd, others...
          col_order <- c("location_id", "year", "Veg", "coef_025", "coef_975", "coef_sd")
          extra <- setdiff(names(unc$coef_ci), col_order)
          unc$coef_ci <- unc$coef_ci[, c(intersect(col_order, names(unc$coef_ci)), extra), drop = FALSE]
        }
        # Ensure uncertainty has a coefficient CI table (fallback with NA rows if not computed)
        if (!is.null(unc)) {
          if (is.null(unc$coef_ci) || !is.data.frame(unc$coef_ci) || nrow(unc$coef_ci) == 0) {
            # Fallback to NA-filled CI rows matching the coefficients in coef_df (non-barren vegs)
            vegs_list <- coef_df$Veg[coef_df$Veg != "barren"]
            if (length(vegs_list) > 0) {
              ci_fb <- data.frame(Veg = vegs_list, coef_025 = NA_real_, coef_975 = NA_real_, coef_sd = NA_real_, stringsAsFactors = FALSE)
              ci_fb$location_id <- loc; ci_fb$year <- yr
              unc$coef_ci <- ci_fb[, c("location_id","year","Veg","coef_025","coef_975","coef_sd")]
              warning(sprintf("[UNCERTAINTY] No coefficient bootstrap CI computed for %s/%s; filled with NA rows", loc, yr))
            }
          }
        }

        # Hybrid propagation removed: the nested two-stage bootstrap provides combined CIs
        # by recomputing both Stage 1 and Stage 2 per bootstrap replicate. Reintroduce a
        # dedicated helper and restore the propagation block only if you need to compare
        # alternative propagation approaches (not recommended for default runs).

        # Add 95% CI interval to coef_df (difference between upper and lower CI bounds)
        coef_df$interval <- ifelse(is.finite(coef_df$coef_975) & is.finite(coef_df$coef_025), coef_df$coef_975 - coef_df$coef_025, NA_real_)

        return(list(
          coef_df = coef_df,
          variant_trajectory = variant_info_pca,
          diagnostics = diag_df,
          uncertainty = unc,
          q10_dvi = q10_dvi,
          q90_dvi = q90_dvi,
          vegetated_fraction = veg_fraction_total,
          barren_fraction = barren_fraction_total,
          inseparable_flag = inseparable_flag,
          inseparable_details = inseparable_text
        ))
      },
      error = function(e) {
        dbg_return_null(paste0("error:", as.character(e$message)))
      }
    )
    res_safe
  }

  # Pre-compute Optimized Library for Vectorized Similarity Search
  cat("Pre-computing optimized library for vectorized similarity search...\n")
  OPTIMIZED_LIBRARY <- precompute_optimized_library(mesma_lib, compressed_templates_accessor, grid_type = "full", feature_weights = LDA_FEATURE_WEIGHTS)
  assign("OPTIMIZED_LIBRARY", OPTIMIZED_LIBRARY, envir = globalenv())

  # Prepare task function environment
  # OPTIMIZED: Removed 'df' and 'df_tasks' from required_globals to reduce serialization overhead
  required_globals <- list(
    lib = lib,
    raw_lib_templates = raw_lib_templates,
    mesma_lib = mesma_lib,
    avail = avail,
    LOWER_BND = LOWER_BND,
    EPS_SIGMA = EPS_SIGMA,
    MAX_VEG_COMPONENTS = MAX_VEG_COMPONENTS,
    veg_counts = veg_counts,
    OUT_DIR = OUT_DIR,
    # df_tasks = df_tasks, # REMOVED
    prepare_factor_data = prepare_factor_data,
    compress_trace = compress_trace_unified,
    compress_temporal_matrix = compress_temporal_matrix,
    build_temporal_matrix = build_temporal_matrix,
    build_raw_index_matrix = build_raw_index_matrix,
    compress_trace_unified = compress_trace_unified,
    compress_stage1_lib_unified = compress_stage1_lib_unified,
    precompute_compressed_templates_unified = precompute_compressed_templates_unified,
    reduce_all_traces_unified = reduce_all_traces_unified,
    build_mesma_variants_unified = build_mesma_variants_unified,
    precompute_templates_from_unified_lib = precompute_templates_from_unified_lib,
    unmix_stage2_compressed = unmix_stage2_compressed,
    COMPRESSED_STAGE1_LIB = if (exists("COMPRESSED_STAGE1_LIB")) COMPRESSED_STAGE1_LIB else NULL,
    STAGE1_LIB = if (exists("STAGE1_LIB")) STAGE1_LIB else NULL,
    compressed_templates_accessor = compressed_templates_accessor,
    # bootstrap helpers (exposed for worker env)
    ols_block_bootstrap = ols_block_bootstrap,
    estimate_block_size = estimate_block_size,
    cos_sim = cos_sim,
    BOOTSTRAP_B = BOOTSTRAP_B,
    ENABLE_UNCERTAINTY = ENABLE_UNCERTAINTY,
    ENABLE_QP_SOLVER = ENABLE_QP_SOLVER,
    compute_diagnostics = compute_diagnostics,
    project_to_simplex = project_to_simplex,
    dbg_return_null = dbg_return_null,
    inference_location_ids = if (exists("inference_location_ids")) inference_location_ids else character(0),
    MIN_UNIQUE_DOY_DEFAULT = MIN_UNIQUE_DOY_DEFAULT,
    MIN_UNIQUE_DOY_INFERENCE = MIN_UNIQUE_DOY_INFERENCE,
    nested_two_stage_bootstrap = if (exists("nested_two_stage_bootstrap")) nested_two_stage_bootstrap else NULL,
    stage1_block_bootstrap = if (exists("stage1_block_bootstrap")) stage1_block_bootstrap else NULL,
    build_Z365 = if (exists("build_Z365")) build_Z365 else NULL,
    DEBUG_UNCERTAINTY = DEBUG_UNCERTAINTY,
    VARIANT_SWITCH_BOOTSTRAP = VARIANT_SWITCH_BOOTSTRAP,
    INSEPARABLE_VARIANT_INFO = if (exists("INSEPARABLE_VARIANT_INFO")) INSEPARABLE_VARIANT_INFO else NULL,
    INSEPARABLE_VARIANTS = if (exists("INSEPARABLE_VARIANTS")) INSEPARABLE_VARIANTS else NULL,
    VARIANT_SIMILARITY_TABLE = if (exists("VARIANT_SIMILARITY_TABLE")) VARIANT_SIMILARITY_TABLE else NULL,
    LDA_FEATURE_WEIGHTS = if (exists("LDA_FEATURE_WEIGHTS")) LDA_FEATURE_WEIGHTS else NULL,
    OPTIMIZED_LIBRARY = OPTIMIZED_LIBRARY
  )

  env_task <- list2env(required_globals, parent = globalenv())
  environment(fit_one_task) <- env_task

  # Main processing loop
  # Report training/inference dataset sizes (in # location-years) and check for accidental duplication
  cat("Starting main processing loop...\n")

  # Compute unique location-year counts for training (df_train) and inference (df_tasks_inference or df_tasks)
  n_train_loc_years <- 0L
  n_infer_loc_years <- 0L
  if (exists("df_train") && !is.null(df_train) && nrow(df_train) > 0) {
    if (!"year" %in% names(df_train) && "date" %in% names(df_train)) df_train$year <- lubridate::year(df_train$date)
    n_train_loc_years <- nrow(unique(df_train[c("location_id", "year")]))
  }
  if (exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0) {
    if (!"year" %in% names(df_tasks_inference) && "date" %in% names(df_tasks_inference)) df_tasks_inference$year <- lubridate::year(df_tasks_inference$date)
    n_infer_loc_years <- nrow(unique(df_tasks_inference[c("location_id", "year")]))
  } else if (exists("df_tasks") && !is.null(df_tasks) && nrow(df_tasks) > 0) {
    if (!"year" %in% names(df_tasks) && "date" %in% names(df_tasks)) df_tasks$year <- lubridate::year(df_tasks$date)
    n_infer_loc_years <- nrow(unique(df_tasks[c("location_id", "year")]))
  }

  cat(sprintf("Training dataset location-years: %d\n", n_train_loc_years))
  cat(sprintf("Inference dataset location-years: %d\n", n_infer_loc_years))

  # If separate inference data exists, ensure we are not accidentally running on identically sized loc-year sets
  if ((exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0) && n_train_loc_years == n_infer_loc_years) {
    stop(sprintf("ERROR: Training and inference datasets appear to have the same number of location-years (%d). This may indicate you passed the same data for training and inference — aborting to avoid accidental overlap.", n_train_loc_years))
  }
  
  # OPTIMIZED: Batched processing to minimize RAM usage
  # Instead of splitting the entire dataset at once (which duplicates it in memory),
  # we process in batches.
  cat("Preparing task keys for batched processing...\n")
  
  # Ensure task key exists
  if (!"task_key" %in% names(df_tasks)) {
    df_tasks$task_key <- paste(df_tasks$location_id, df_tasks$year, sep = "_")
  }
  
  # Identify target tasks
  target_keys <- paste(test_loc_years$location_id, test_loc_years$year, sep = "_")
  # Only keep keys that actually exist in the data
  available_keys <- unique(df_tasks$task_key)
  target_keys <- intersect(target_keys, available_keys)
  
  n_keys <- length(target_keys)
  BATCH_SIZE <- 200 # Adjust based on available RAM
  
  # Create batches of keys
  key_batches <- split(target_keys, ceiling(seq_along(target_keys) / BATCH_SIZE))
  
  cat(sprintf("Processing %d tasks in %d batches (approx %d tasks/batch)...\n", n_keys, length(key_batches), BATCH_SIZE))
  
  # Pre-allocate results list to avoid growing it
  results_list <- vector("list", n_keys)
  names(results_list) <- target_keys
  
  start_time <- Sys.time()
  
  # Process batches
  for (i in seq_along(key_batches)) {
    batch_keys <- key_batches[[i]]
    cat(sprintf("  Batch %d/%d: Processing %d tasks...\n", i, length(key_batches), length(batch_keys)))
    
    # Subset data for this batch ONLY
    # This avoids creating a copy of the full dataset
    batch_df <- df_tasks[df_tasks$task_key %in% batch_keys, ]
    
    # Split into individual tasks for this batch
    batch_task_list <- split(batch_df, batch_df$task_key)
    
    # Run parallel map on this batch
    batch_results <- .run_map(batch_task_list, fit_one_task)
    
    # Store results
    # We match by name to ensure correct placement
    results_list[names(batch_results)] <- batch_results
    
    # Explicitly clean up to free RAM before next batch
    rm(batch_df, batch_task_list, batch_results)
    gc(verbose = FALSE)
  }
  
  end_time <- Sys.time()
  processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(
    "Main processing loop finished in %.2f seconds (%.2f minutes)\n",
    processing_time, processing_time / 60
  ))
  cat(sprintf("Average time per task: %.2f seconds\n", processing_time / n_keys))

  # Process results and write to Excel
  cat("Processing results and writing to Excel files...\n")

  # Filter out null results and any results with missing/invalid location IDs
  results_list <- results_list[!sapply(results_list, is.null)]
  # Exclude results where coef_df is missing or location_id is entirely NA/empty
  is_invalid_res <- sapply(results_list, function(res) {
    if (is.null(res)) return(TRUE)
    if (is.null(res$coef_df) || !is.data.frame(res$coef_df)) return(TRUE)
    lid <- res$coef_df$location_id
    if (is.null(lid) || length(lid) == 0) return(TRUE)
    # consider empty string or all NA characters as invalid
    all(is.na(lid) | trimws(as.character(lid)) == "")
  })
  if (any(is_invalid_res)) {
    cat(sprintf("[NOTICE] Filtering out %d invalid result(s) with missing location_id\n", sum(is_invalid_res)))
    results_list <- results_list[!is_invalid_res]
  }
  cat(sprintf("After filtering NULL results: %d results remaining\n", length(results_list)))

  # Check for barren fraction = 1 in more than 10% of predictions
  if (length(results_list) > 0) {
    barren_one_count <- sum(sapply(results_list, function(res) {
      if (!is.null(res$barren_fraction) && is.finite(res$barren_fraction)) {
        res$barren_fraction == 1
      } else {
        FALSE
      }
    }))
    barren_one_pct <- barren_one_count / length(results_list) * 100
    if (barren_one_pct > 50) {
      stop(sprintf("ERROR: Barren fraction = 1 in %.1f%% of predictions (>%d/%d). This indicates unnormalized endmembers causing bias toward barren. Check endmember normalization in geometric stage 1 unmixing.", barren_one_pct, barren_one_count, length(results_list)))
    }
    cat(sprintf("Barren fraction = 1 in %.1f%% of predictions (%d/%d) - within acceptable limits\n", barren_one_pct, barren_one_count, length(results_list)))
  } else {
    cat("No results to check for barren fraction\n")
  }

  if (length(results_list) == 0) {
    cat("ERROR: All tasks returned NULL results!\n")
    cat("Most likely causes:\n")
    cat("1. No valid testing/inference data available (no location-year pairs found)\n")
    cat("2. Data filtering issues or missing indices\n")
    cat("3. All locations have insufficient data for fitting\n")
    stop("No valid results to process")
  }

  if (length(results_list) > 0) {
    # Combine all coefficient data frames
    cat("Combining coefficient data frames...\n")
    coef_list <- lapply(results_list, function(res) {
      if (is.null(res$coef_df) || nrow(res$coef_df) == 0) {
        NULL
      } else {
        res$coef_df
      }
    })
    coef_list <- coef_list[!sapply(coef_list, is.null)]

    if (length(coef_list) == 0) {
      cat("ERROR: No valid coefficient data frames found!\n")
      stop("No coefficient data to process")
    }

    # Attempt to combine coefficient data frames robustly. Use dplyr::bind_rows
    # to gracefully handle differing column sets (fills missing columns with NA).
    all_coefs <- tryCatch({
      if (requireNamespace("dplyr", quietly = TRUE)) {
        # Debug: print unique column sets if DEBUG_UNCERTAINTY
        if (exists("DEBUG_UNCERTAINTY") && isTRUE(DEBUG_UNCERTAINTY)) {
          col_sets <- lapply(seq_along(coef_list), function(i) { list(i = i, cols = names(coef_list[[i]]), ncol = ncol(coef_list[[i]])) })
          uniq <- unique(sapply(col_sets, function(x) paste(sort(x$cols), collapse = ",")))
          cat(sprintf("[DEBUG] combining %d coefficient frames; unique column sets: %d\n", length(col_sets), length(uniq)))
          if (length(uniq) > 1) {
            cat("[DEBUG] Column sets per frame (first 10 shown):\n")
            for (i in seq_len(min(length(col_sets), 10))) {
              cat(sprintf("  frame %d: %d cols: %s\n", col_sets[[i]]$i, col_sets[[i]]$ncol, paste(col_sets[[i]]$cols, collapse = ",")))
            }
          }
        }
        dplyr::bind_rows(coef_list)
      } else {
        # Fallback: attempt rbind and provide a clearer error if it fails
        do.call(rbind, coef_list)
      }
    }, error = function(e) stop(sprintf("ERROR combining coef_df: %s", e$message))
    )

    if (is.null(all_coefs)) {
      cat("Failed to combine coefficient data frames\n")
      stop("Cannot proceed without coefficient data")
    }

    # Ensure required coef columns exist and are consistent across rows
    required_coef_cols <- c("location_id", "year", "Veg", "coef", "rmse", "coef_025", "coef_975", "interval")
    missing <- setdiff(required_coef_cols, names(all_coefs))
    if (length(missing) > 0) {
      cat(sprintf("[NOTICE] Filling missing coefficient columns with NA: %s\n", paste(missing, collapse = ", ")))
      for (col in missing) all_coefs[[col]] <- NA
    }
    # Type coercion for predictable downstream joins/aggregations
    all_coefs$location_id <- as.character(all_coefs$location_id)
    all_coefs$year <- as.integer(all_coefs$year)
    all_coefs$Veg <- as.character(all_coefs$Veg)
    all_coefs$coef <- as.numeric(all_coefs$coef)
    all_coefs$rmse <- as.numeric(all_coefs$rmse)
    all_coefs$coef_025 <- as.numeric(all_coefs$coef_025)
    all_coefs$coef_975 <- as.numeric(all_coefs$coef_975)
    all_coefs$interval <- as.numeric(all_coefs$interval)

    # Trim whitespace from location_id to handle any remaining issues
    all_coefs$location_id <- trimws(all_coefs$location_id)
    cat(sprintf("Combined coefficients: %d rows\n", nrow(all_coefs)))

  # Combine chosen variant summaries (no per-DOY trajectories in new approach)
    cat("Combining chosen variant summaries...\n")
    variant_list_pca <- lapply(results_list, function(res) if (!is.null(res$variant_trajectory)) res$variant_trajectory else NULL)
    variant_list_pca <- variant_list_pca[!sapply(variant_list_pca, is.null)]
    all_variants_pca <- NULL
    if (length(variant_list_pca) > 0) {
      all_variants_pca <- tryCatch({
        if (requireNamespace("dplyr", quietly = TRUE)) {
          dplyr::bind_rows(variant_list_pca)
        } else {
          do.call(rbind, variant_list_pca)
        }
      }, error = function(e) stop(sprintf("Failed to combine variant PCA summaries: %s", e$message)))
    }

    # Collect diagnostics
    diag_list <- lapply(results_list, function(res) {
      if (!is.null(res$diagnostics)) res$diagnostics else NULL
    })
    diag_list <- diag_list[!sapply(diag_list, is.null)]
    # Use dplyr::bind_rows so differing diagnostic column sets are safely merged (missing columns filled with NA)
    all_diagnostics <- if (length(diag_list) > 0) tryCatch({
      if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required to combine diagnostics robustly")
      dplyr::bind_rows(diag_list)
    }, error = function(e) stop(sprintf("Failed to combine diagnostics: %s", e$message))) else NULL

    # Collect Q10 and Q90 DVI values
    q_dvi_data <- tryCatch({
      if (requireNamespace("dplyr", quietly = TRUE)) {
        # Use bind_rows to allow some results to be NULL
        dplyr::bind_rows(lapply(results_list, function(res) {
          if ((!is.null(res$q10_dvi) && !is.na(res$q10_dvi)) || (!is.null(res$q90_dvi) && !is.na(res$q90_dvi))) {
            data.frame(
              location_id = res$coef_df$location_id[1],
              year = res$coef_df$year[1],
              q10_dvi = res$q10_dvi,
              q90_dvi = res$q90_dvi,
              stringsAsFactors = FALSE
            )
          } else {
            NULL
          }
        }))
      } else {
        do.call(rbind, lapply(results_list, function(res) {
          if ((!is.null(res$q10_dvi) && !is.na(res$q10_dvi)) || (!is.null(res$q90_dvi) && !is.na(res$q90_dvi))) {
            data.frame(
              location_id = res$coef_df$location_id[1],
              year = res$coef_df$year[1],
              q10_dvi = res$q10_dvi,
              q90_dvi = res$q90_dvi,
              stringsAsFactors = FALSE
            )
          } else {
            NULL
          }
        }))
      }
    }, error = function(e) {
      warning(sprintf("q_dvi_data assembly failed: %s", e$message)); NULL
    })

    variant_similarity_table <- NULL
    variant_similarity_summary <- NULL
    if (exists("VARIANT_SIMILARITY_TABLE")) {
      variant_similarity_table <- VARIANT_SIMILARITY_TABLE
    } else if (exists("INSEPARABLE_VARIANT_INFO") && !is.null(INSEPARABLE_VARIANT_INFO$similarity_table)) {
      variant_similarity_table <- INSEPARABLE_VARIANT_INFO$similarity_table
    }
    if (!is.null(variant_similarity_table) && nrow(variant_similarity_table) > 0) {
      variant_similarity_table <- variant_similarity_table[order(variant_similarity_table$veg, -variant_similarity_table$cos_sim, variant_similarity_table$euclidean_dist), , drop = FALSE]
      if (requireNamespace("dplyr", quietly = TRUE)) {
        variant_similarity_summary <- variant_similarity_table %>%
          dplyr::group_by(.data$veg) %>%
          dplyr::summarise(
            pair_count = dplyr::n(),
            max_cos_sim = if (all(is.na(.data$cos_sim))) NA_real_ else max(.data$cos_sim, na.rm = TRUE),
            min_cos_sim = if (all(is.na(.data$cos_sim))) NA_real_ else min(.data$cos_sim, na.rm = TRUE),
            median_cos_sim = if (all(is.na(.data$cos_sim))) NA_real_ else stats::median(.data$cos_sim, na.rm = TRUE),
            min_euclidean_dist = if (all(is.na(.data$euclidean_dist))) NA_real_ else min(.data$euclidean_dist, na.rm = TRUE),
            median_euclidean_dist = if (all(is.na(.data$euclidean_dist))) NA_real_ else stats::median(.data$euclidean_dist, na.rm = TRUE),
            .groups = "drop"
          )
      } else {
        by_veg <- split(variant_similarity_table, variant_similarity_table$veg)
        summary_rows <- lapply(names(by_veg), function(veg_name) {
          tbl <- by_veg[[veg_name]]
          data.frame(
            veg = veg_name,
            pair_count = nrow(tbl),
            max_cos_sim = if (all(is.na(tbl$cos_sim))) NA_real_ else max(tbl$cos_sim, na.rm = TRUE),
            min_cos_sim = if (all(is.na(tbl$cos_sim))) NA_real_ else min(tbl$cos_sim, na.rm = TRUE),
            median_cos_sim = if (all(is.na(tbl$cos_sim))) NA_real_ else stats::median(tbl$cos_sim, na.rm = TRUE),
            min_euclidean_dist = if (all(is.na(tbl$euclidean_dist))) NA_real_ else min(tbl$euclidean_dist, na.rm = TRUE),
            median_euclidean_dist = if (all(is.na(tbl$euclidean_dist))) NA_real_ else stats::median(tbl$euclidean_dist, na.rm = TRUE),
            stringsAsFactors = FALSE
          )
        })
        variant_similarity_summary <- do.call(rbind, summary_rows)
      }
    }

    # Get unique locations (remove NA and empty strings to avoid Excel row-name errors)
    cat(sprintf("[DEBUG] all_coefs$location_id class: %s\n", class(all_coefs$location_id)))
    cat(sprintf("[DEBUG] all_coefs$location_id sample: %s\n", paste(head(all_coefs$location_id, 10), collapse = ", ")))
    cat(sprintf("[DEBUG] all_coefs$location_id is.na sum: %d\n", sum(is.na(all_coefs$location_id))))
    cat(sprintf("[DEBUG] all_coefs$location_id == '' sum: %d\n", sum(all_coefs$location_id == "", na.rm = TRUE)))
    unique_locations <- unique(trimws(as.character(all_coefs$location_id)))
    unique_locations <- unique_locations[!is.na(unique_locations) & unique_locations != ""]
    if (length(unique_locations) == 0) {
      stop("No valid location IDs found in results")
    }
    cat(sprintf("Creating Excel file with %d valid locations\n", length(unique_locations)))

    # Join with ground truth vegetation data
    true_veg_map <- gpts_map %>% dplyr::select(location_id, true_veg = Veg)
    # Coerce types to avoid incompatible-type join errors
    if (!is.character(all_coefs$location_id)) {
      all_coefs$location_id <- as.character(all_coefs$location_id)
      cat("[NOTICE] Coerced all_coefs$location_id to character to match true_veg_map for joining.\n")
    }
    if (!is.character(true_veg_map$location_id)) {
      true_veg_map$location_id <- as.character(true_veg_map$location_id)
      cat("[NOTICE] Coerced true_veg_map$location_id to character for joining.\n")
    }
    if (length(all_coefs$location_id) > 0 && length(true_veg_map$location_id) > 0) {
      s1 <- unique(na.omit(all_coefs$location_id))
      s2 <- unique(na.omit(true_veg_map$location_id))
      if (length(s1) && length(s2) && all(grepl("^[0-9]+$", s1)) && any(grepl("^L_", s2))) {
        cat("[WARNING] all_coefs$location_id looks numeric while true_veg_map$location_id looks like 'L_x_y' strings — matching will likely fail.\n")
      }
    }
    all_coefs <- dplyr::left_join(all_coefs, true_veg_map, by = "location_id")
    # Ensure CI columns exist and compute interval (upper - lower) for each row
    if (!"coef_025" %in% names(all_coefs)) all_coefs$coef_025 <- NA_real_
    if (!"coef_975" %in% names(all_coefs)) all_coefs$coef_975 <- NA_real_
    # Interval is already computed in fit_one_task
    if (DEBUG_UNCERTAINTY) {
      cat("After combining all_coefs with true_veg_map:\n")
      cat("Number of rows with finite coef_025:", sum(is.finite(all_coefs$coef_025)), "\n")
      cat("Number of rows with finite coef_975:", sum(is.finite(all_coefs$coef_975)), "\n")
      cat("Number of rows with finite interval:", sum(is.finite(all_coefs$interval)), "\n")
    }

    inseparable_flags_detail <- NULL
    inseparable_flags_summary <- NULL
    if ("inseparable_variant_flag" %in% names(all_coefs)) {
      flagged_rows <- which(!is.na(all_coefs$inseparable_variant_flag) & all_coefs$inseparable_variant_flag)
      if (length(flagged_rows) > 0) {
        cols_keep <- intersect(c("location_id", "year", "Veg", "coef", "rmse", "true_veg", "inseparable_variant_details"), names(all_coefs))
        inseparable_flags_detail <- all_coefs[flagged_rows, cols_keep, drop = FALSE]
        inseparable_flags_detail <- inseparable_flags_detail[order(inseparable_flags_detail$location_id, inseparable_flags_detail$year, inseparable_flags_detail$Veg), , drop = FALSE]

        if (requireNamespace("dplyr", quietly = TRUE)) {
          inseparable_flags_summary <- inseparable_flags_detail %>%
            dplyr::group_by(.data$location_id, .data$year) %>%
            dplyr::summarise(
              veg_components = paste(sort(unique(.data$Veg)), collapse = ", "),
              max_coef = if (all(is.na(.data$coef))) NA_real_ else max(.data$coef, na.rm = TRUE),
              details = paste(unique(na.omit(.data$inseparable_variant_details)), collapse = " | "),
              .groups = "drop"
            ) %>%
            dplyr::arrange(.data$location_id, .data$year)
        } else {
          summary_keys <- unique(inseparable_flags_detail[, c("location_id", "year"), drop = FALSE])
          summary_keys <- summary_keys[order(summary_keys$location_id, summary_keys$year), , drop = FALSE]
          build_summary <- function(loc, yr) {
            rows <- inseparable_flags_detail$location_id == loc & inseparable_flags_detail$year == yr
            vegs <- sort(unique(inseparable_flags_detail$Veg[rows]))
            coef_vals <- inseparable_flags_detail$coef[rows]
            details <- unique(na.omit(inseparable_flags_detail$inseparable_variant_details[rows]))
            data.frame(
              location_id = loc,
              year = yr,
              veg_components = paste(vegs, collapse = ", "),
              max_coef = if (all(is.na(coef_vals))) NA_real_ else max(coef_vals, na.rm = TRUE),
              details = paste(details, collapse = " | "),
              stringsAsFactors = FALSE
            )
          }
          if (nrow(summary_keys) > 0) {
            summary_list <- apply(summary_keys, 1, function(row) build_summary(row[["location_id"]], as.integer(row[["year"]])))
            inseparable_flags_summary <- do.call(rbind, summary_list)
          }
        }
      }
    }

    # Create Excel file with all locations as separate sheets
    cat("Creating single Excel file with all location results...\n")

    wb <- tryCatch(
      {
        openxlsx::createWorkbook()
      },
      error = function(e) {
        cat(sprintf("ERROR creating workbook: %s\n", e$message))
        return(NULL)
      }
    )

    if (is.null(wb)) {
      cat("Failed to create workbook\n")
      stop("Cannot create Excel output")
    }

    # Add summary sheet
    summary_data <- data.frame(
      Location_ID = unique_locations,
      Total_Years = sapply(unique_locations, function(loc) {
        length(unique(all_coefs$pheno_year[all_coefs$location_id == loc]))
      }),
      Total_Observations = sapply(unique_locations, function(loc) {
        nrow(all_coefs[all_coefs$location_id == loc, ])
      }),
      # Projection method column not included: a single method is used per run
      stringsAsFactors = FALSE
    )

    openxlsx::addWorksheet(wb, "Summary")
    openxlsx::writeData(wb, "Summary", summary_data)

    # Diagnostics sheet
    if (!is.null(all_diagnostics) && nrow(all_diagnostics) > 0) {
      openxlsx::addWorksheet(wb, "Diagnostics")
      openxlsx::writeData(wb, "Diagnostics", all_diagnostics)
    }

    # Add MESMA Variant Summary sheet if variant data exists (chosen variants only)
            if (!is.null(all_variants_pca) && nrow(all_variants_pca) > 0) {
      openxlsx::addWorksheet(wb, "Variant_Summary")
      openxlsx::writeData(wb, "Variant_Summary", all_variants_pca)
    }

    if (!is.null(variant_similarity_table) && nrow(variant_similarity_table) > 0) {
      openxlsx::addWorksheet(wb, "Variant_Similarity")
      openxlsx::writeData(wb, "Variant_Similarity", "Pairwise similarity across variants (cosine similarity and Euclidean distance)", startRow = 1, startCol = 1)
      start_row <- 2
      if (!is.null(variant_similarity_summary) && nrow(variant_similarity_summary) > 0) {
        openxlsx::writeData(wb, "Variant_Similarity", variant_similarity_summary, startRow = start_row, startCol = 1)
        start_row <- start_row + nrow(variant_similarity_summary) + 2
      }
      openxlsx::writeData(wb, "Variant_Similarity", variant_similarity_table, startRow = start_row, startCol = 1)
    }

    # Add Inter-Class Similarity sheet
    if (exists("INTER_CLASS_SIMILARITY") && !is.null(INTER_CLASS_SIMILARITY) && nrow(INTER_CLASS_SIMILARITY) > 0) {
      openxlsx::addWorksheet(wb, "Inter_Veg_Similarity")
      openxlsx::writeData(wb, "Inter_Veg_Similarity", "Pairwise similarity between variants of different vegetation types", startRow = 1, startCol = 1)
      openxlsx::writeData(wb, "Inter_Veg_Similarity", INTER_CLASS_SIMILARITY, startRow = 2, startCol = 1)
    }

    if (!is.null(inseparable_flags_summary) && nrow(inseparable_flags_summary) > 0) {
      openxlsx::addWorksheet(wb, "Inseparable_Flags")
      openxlsx::writeData(wb, "Inseparable_Flags", "Location-years with inseparable MESMA variants detected", startRow = 1, startCol = 1)
      openxlsx::writeData(wb, "Inseparable_Flags", inseparable_flags_summary, startRow = 2, startCol = 1)
      next_row <- nrow(inseparable_flags_summary) + 4
      if (!is.null(inseparable_flags_detail) && nrow(inseparable_flags_detail) > 0) {
        openxlsx::writeData(wb, "Inseparable_Flags", "Detailed component rows", startRow = next_row, startCol = 1)
        openxlsx::writeData(wb, "Inseparable_Flags", inseparable_flags_detail, startRow = next_row + 1, startCol = 1)
      }
    }

    # Collect uncertainty across locations if available
    unc_coef_rows <- list(); unc_var_rows <- list(); unc_rmse_rows <- list(); unc_meta_rows <- list()
    for (res in results_list) {
      if (is.null(res$uncertainty)) next
      loc <- if (!is.null(res$coef_df$location_id)) res$coef_df$location_id[1] else NA_character_
      yr <- if (!is.null(res$coef_df$year)) res$coef_df$year[1] else NA_integer_
      ci <- res$uncertainty$coef_ci
      if (!is.null(ci) && nrow(ci) > 0) {
        ci$location_id <- loc; ci$year <- yr
        unc_coef_rows[[length(unc_coef_rows) + 1]] <- ci[, c("location_id","year","Veg","coef_025","coef_975"), drop = FALSE]
      }
      vr <- res$uncertainty$variant_freq
      if (!is.null(vr) && nrow(vr) > 0) {
        vr$location_id <- loc; vr$year <- yr
        unc_var_rows[[length(unc_var_rows) + 1]] <- vr[, c("location_id","year","Veg","Variant","N","Percent"), drop = FALSE]
      }
      rc <- res$uncertainty$rmse_ci
      if (!is.null(rc) && length(rc) == 2) {
        unc_rmse_rows[[length(unc_rmse_rows) + 1]] <- data.frame(location_id = loc, year = yr, rmse_lo = rc[1], rmse_hi = rc[2], stringsAsFactors = FALSE)
      }
      if (!is.null(res$uncertainty$block_size)) {
        unc_meta_rows[[length(unc_meta_rows) + 1]] <- data.frame(location_id = loc, year = yr, block_size = res$uncertainty$block_size, B = BOOTSTRAP_B, stringsAsFactors = FALSE)
      }
    }
    if (length(unc_coef_rows) > 0) {
      if (requireNamespace("dplyr", quietly = TRUE)) all_unc_coef <- dplyr::bind_rows(unc_coef_rows) else all_unc_coef <- do.call(rbind, unc_coef_rows)
    } else all_unc_coef <- NULL
    if (!is.null(all_unc_coef)) {
      n_total_unc <- nrow(all_unc_coef)
      n_non_na_ci <- sum(is.finite(all_unc_coef$coef_025) | is.finite(all_unc_coef$coef_975))
      cat(sprintf("Uncertainty CIs: total rows = %d; rows with at least one finite CI = %d\n", n_total_unc, n_non_na_ci))
    }
    if (length(unc_var_rows) > 0) {
      if (requireNamespace("dplyr", quietly = TRUE)) all_unc_var <- dplyr::bind_rows(unc_var_rows) else all_unc_var <- do.call(rbind, unc_var_rows)
    } else all_unc_var <- NULL
    if (length(unc_rmse_rows) > 0) {
      if (requireNamespace("dplyr", quietly = TRUE)) all_unc_rmse <- dplyr::bind_rows(unc_rmse_rows) else all_unc_rmse <- do.call(rbind, unc_rmse_rows)
    } else all_unc_rmse <- NULL
    if (length(unc_meta_rows) > 0) {
      if (requireNamespace("dplyr", quietly = TRUE)) all_unc_meta <- dplyr::bind_rows(unc_meta_rows) else all_unc_meta <- do.call(rbind, unc_meta_rows)
    } else all_unc_meta <- NULL

    if (!is.null(all_unc_coef) || !is.null(all_unc_var) || !is.null(all_unc_rmse)) {
      openxlsx::addWorksheet(wb, "Uncertainty")
      start_row <- 1
      openxlsx::writeData(wb, "Uncertainty", data.frame(Setting = c("ENABLE_UNCERTAINTY","BOOTSTRAP_B","BOOTSTRAP_METHOD"), Value = c(ENABLE_UNCERTAINTY, BOOTSTRAP_B, if (isTRUE(ENABLE_UNCERTAINTY)) "nested" else "none")), startRow = start_row, startCol = 1)
      start_row <- start_row + 3
      # Coefficient CI table suppressed: coefficients CIs are now included in the COEFFICIENTS sheet as columns
      if (!is.null(all_unc_rmse)) {
        openxlsx::writeData(wb, "Uncertainty", "RMSE CI (2.5%/97.5%)", startRow = start_row, startCol = 1)
        openxlsx::writeData(wb, "Uncertainty", all_unc_rmse, startRow = start_row + 1, startCol = 1)
        start_row <- start_row + nrow(all_unc_rmse) + 3
      }
      if (!is.null(all_unc_var)) {
        openxlsx::writeData(wb, "Uncertainty", "Variant Dominance Frequencies (%)", startRow = start_row, startCol = 1)
        openxlsx::writeData(wb, "Uncertainty", all_unc_var, startRow = start_row + 1, startCol = 1)
      }
      if (!is.null(all_unc_meta)) {
        start_row <- start_row + ifelse(!is.null(all_unc_var), nrow(all_unc_var) + 3, 3)
        openxlsx::writeData(wb, "Uncertainty", "Bootstrap Meta", startRow = start_row, startCol = 1)
        openxlsx::writeData(wb, "Uncertainty", all_unc_meta, startRow = start_row + 1, startCol = 1)
      }
    }

    # Endmember Stability sheet
    if (exists('stability_results') && !is.null(stability_results)) {
      openxlsx::addWorksheet(wb, 'Endmember_Stability')
      current_row <- 1
      for (veg in names(stability_results)) {
        res <- stability_results[[veg]]
        openxlsx::writeData(wb, 'Endmember_Stability', sprintf('=== %s ENDMEMBER STABILITY ===', toupper(veg)), startRow = current_row, startCol = 1)
        current_row <- current_row + 2
        openxlsx::writeData(wb, 'Endmember_Stability', 'Variant Count Distribution:', startRow = current_row, startCol = 1)
        vcd <- as.data.frame(res$variant_count_distribution)
        names(vcd) <- c('N_Variants', 'Frequency')
        openxlsx::writeData(wb, 'Endmember_Stability', vcd, startRow = current_row + 1, startCol = 1)
        current_row <- current_row + nrow(vcd) + 3
        openxlsx::writeData(wb, 'Endmember_Stability', 'Meta-Variant Stability Metrics:', startRow = current_row, startCol = 1)
        mv_stats <- if (length(res$meta_variants) > 0) do.call(rbind, lapply(res$meta_variants, function(mv) {
          data.frame(Meta_Variant = mv$meta_variant_id, N_Bootstrap_Members = mv$n_members, Coefficient_of_Variation = mv$coefficient_of_variation, Spectral_Angle_Mean_deg = mv$mean_spectral_angle, Spectral_Angle_SD_deg = sqrt(mv$spectral_angle_variance), Bootstrap_Frequency_pct = length(unique(mv$bootstrap_iters)) / res$n_bootstrap_iters * 100, stringsAsFactors = FALSE)
        })) else data.frame()
        if (nrow(mv_stats) > 0) {
          openxlsx::writeData(wb, 'Endmember_Stability', mv_stats, startRow = current_row + 1, startCol = 1)
          current_row <- current_row + nrow(mv_stats) + 4
        } else {
          current_row <- current_row + 4
        }
      }
      openxlsx::writeData(wb, 'Endmember_Stability', 'INTERPRETATION GUIDE', startRow = current_row, startCol = 1)
      current_row <- current_row + 1
      guide <- data.frame(Metric = c('Variant Count Distribution', 'Coefficient of Variation (CV)', 'Spectral Angle SD', 'Bootstrap Frequency'), Interpretation = c('How many variants are identified across bootstrap iterations. Stable if concentrated.', 'Variability of signature features. Lower CV = more stable endmember.', 'Directional variability in spectral space (degrees). Lower = more stable.', '% of bootstrap iterations where this pattern appears. Higher = more robust.'), Good_Value = c('Narrow distribution (e.g., 90% have same count)', '< 0.10 (very stable), 0.10-0.20 (stable), > 0.20 (unstable)', '< 5° (very stable), 5-10° (stable), > 10° (unstable)', '> 80% (very robust), 50-80% (robust), < 50% (unreliable)'), stringsAsFactors = FALSE)
      openxlsx::writeData(wb, 'Endmember_Stability', guide, startRow = current_row, startCol = 1)
    }

    # Build best fit summary per location-year (single projection per run)
    # NOTE: We now compute both absolute fit (original) and relative fit among
    # vegetation classes (ignoring 'barren'). The relative fit compares
    # predicted fraction of the true vegetation class to the sum of all
    # predicted vegetation fractions (i.e., normalized among vegetation only).
    best_fit_summary <- do.call(rbind, lapply(unique_locations, function(loc) {
      yrs <- unique(all_coefs$pheno_year[all_coefs$location_id == loc])
      tv <- true_veg_map$true_veg[true_veg_map$location_id == loc]
      if (length(tv) == 0) tv <- NA_character_

        do.call(rbind, lapply(yrs, function(yr) {
        row <- all_coefs[all_coefs$location_id == loc & all_coefs$pheno_year == yr & tolower(all_coefs$Veg) == tolower(tv), , drop = FALSE]
        pred <- if (nrow(row) == 1) row$coef else NA_real_
        pred_abs <- pred
        rmse_val <- if (nrow(row) == 1 && "rmse" %in% names(row)) row$rmse else NA_real_
          # Compute sum of predicted vegetated fractions (exclude barren)
          sum_veg_coef <- sum(all_coefs$coef[all_coefs$location_id == loc & all_coefs$pheno_year == yr & tolower(all_coefs$Veg) != "barren"], na.rm = TRUE)
          # If there is a valid predicted fraction and vegetated sum > 0, compute
          # prediction relative to vegetation (ignore barren). Otherwise NA.
          pred_rel <- if (!is.na(pred) && is.finite(sum_veg_coef) && sum_veg_coef > 0) pred / sum_veg_coef else NA_real_
          pred_rel <- pmin(pred_rel, 1)  # Clamp to 1 to prevent >100%
          abs_pct <- if (!is.na(pred_abs) && !is.na(tv)) abs(1 - pred_abs) * 100 else NA_real_
          abs_pct_rel <- if (!is.na(pred_rel) && !is.na(tv)) abs(1 - pred_rel) * 100 else NA_real_

        data.frame(
          location_id = loc,
          year = yr,
          true_veg = tv,
          # pred_coef now stores relative fraction among vegetation only (excludes barren)
          pred_coef = pred_rel,
          pred_coef_abs = pred_abs,
          pred_coef_rel = pred_rel,
          rmse = rmse_val,
          # abs_pct_diff now represents percentage error relative to vegetated-only normalization
          abs_pct_diff = abs_pct_rel,
          abs_pct_diff_abs = abs_pct,
          abs_pct_diff_rel = abs_pct_rel,
          stringsAsFactors = FALSE
        )
      }))
    }))

    eval_years <- sort(unique(c(TRAIN_YEARS, TRAIN_YEARS - 1L, TRAIN_YEARS + 1L)))
    eval_years <- eval_years[is.finite(eval_years)]
    best_fit_summary$eval_window <- best_fit_summary$year %in% eval_years
    best_fit_eval <- best_fit_summary[best_fit_summary$eval_window, , drop = FALSE]
    if (nrow(best_fit_eval) == 0) {
      warning("No location-years fall within the TRAIN_YEARS +/- 1 evaluation window; overall fit cannot be computed")
    }

    # Compute overall fit score (single projection): vegetation-relative only
    # Definition: per location-year, compute predicted fraction for the true Veg class
    # normalized by the sum of predicted vegetation fractions (exclude barren).
    # This represents the average "purity" or accuracy of the vegetation classification.
    overall_fit <- suppressWarnings(as.numeric(mean(best_fit_eval$pred_coef_rel * 100, na.rm = TRUE)))
    if (!is.finite(overall_fit)) overall_fit <- NA_real_

    # Write overall fit score at top of Summary sheet
    # Write the new Overall_Fit_pct metric (vegetation-relative accuracy)
    openxlsx::writeData(wb, "Summary", data.frame(
      Overall_Fit_pct = overall_fit
    ), startRow = 1, startCol = ncol(summary_data) + 2)

    # Add location sheets
    for (i in seq_along(unique_locations)) {
      loc_id <- unique_locations[i]

          # Create sheet name (Excel limits sheet names to 31 characters)
          sheet_name <- substr(gsub("[^A-Za-z0-9]", "_", loc_id), 1, 31)

          # Add worksheet for this location
          openxlsx::addWorksheet(wb, sheet_name)

          # Quality Metrics Calculation
          loc_coefs <- all_coefs[all_coefs$location_id == loc_id, ]

          quality_metrics <- loc_coefs %>%
            dplyr::group_by(.data$year) %>%
            dplyr::summarize(
              deviation = sum(abs(.data$coef - (tolower(.data$Veg) == tolower(true_veg[1])))),
              avg_rmse = mean(.data$rmse, na.rm = TRUE),
              .groups = "drop"
            ) %>%
            dplyr::summarize(
              avg_pct_deviation = mean(.data$deviation, na.rm = TRUE) * 100,
              avg_rmse = mean(.data$avg_rmse, na.rm = TRUE),
              .groups = "drop"
            )

          # Calculate peak Q10 and Q90 DVI for this location
          loc_q_data <- if (!is.null(q_dvi_data)) {
            q_dvi_data[q_dvi_data$location_id == loc_id, ]
          } else {
            NULL
          }

          peak_q10_dvi <- if (!is.null(loc_q_data) && nrow(loc_q_data) > 0) {
            max(loc_q_data$q10_dvi, na.rm = TRUE)
          } else {
            NA
          }

          peak_q90_dvi <- if (!is.null(loc_q_data) && nrow(loc_q_data) > 0) {
            max(loc_q_data$q90_dvi, na.rm = TRUE)
          } else {
            NA
          }

          quality_metrics$peak_q10_dvi <- peak_q10_dvi
          quality_metrics$peak_q90_dvi <- peak_q90_dvi

          # Write quality metrics to the sheet
          openxlsx::writeData(wb, sheet_name, "QUALITY METRICS", startRow = 1, startCol = 1)
          openxlsx::writeData(wb, sheet_name, quality_metrics, startRow = 2, startCol = 1)

          # Data Sections
          current_row <- nrow(quality_metrics) + 4

          # Optional: per-location diagnostics
          loc_diag <- if (exists("all_diagnostics") && !is.null(all_diagnostics)) all_diagnostics[all_diagnostics$location_id == loc_id, , drop = FALSE] else NULL
          if (!is.null(loc_diag) && nrow(loc_diag) > 0) {
            openxlsx::writeData(wb, sheet_name, "DIAGNOSTICS", startRow = current_row, startCol = 1)
            openxlsx::writeData(wb, sheet_name, loc_diag, startRow = current_row + 1, startCol = 1)
            current_row <- current_row + nrow(loc_diag) + 3
          }

          # BEST FIT SUMMARY for this location
          loc_best <- best_fit_summary[best_fit_summary$location_id == loc_id, , drop = FALSE]
          if (!is.null(loc_best) && nrow(loc_best) > 0) {
            desired_cols <- c("year", "true_veg", "pred_coef", "pred_coef_abs", "rmse", "abs_pct_diff", "abs_pct_diff_abs", "eval_window")
            write_tbl <- loc_best[, c("location_id", intersect(desired_cols, names(loc_best))), drop = FALSE]
            openxlsx::writeData(wb, sheet_name, "BEST FIT SUMMARY (per-year) — pred_coef is proportion relative to vegetated area (barren excluded)",
              startRow = current_row, startCol = 1
            )
            openxlsx::writeData(wb, sheet_name, write_tbl, startRow = current_row + 1, startCol = 1)
            current_row <- current_row + nrow(write_tbl) + 3
          }

          # Filter data for the current location
          # Unified coefficients for this location (single projection per run)
          loc_coefs_unified <- all_coefs[all_coefs$location_id == loc_id, ]

          loc_variants_pca <- if (!is.null(all_variants_pca)) {
            all_variants_pca[all_variants_pca$location_id == loc_id, ]
          } else {
            NULL
          }

          # Add MESMA PCA coefficients
          if (nrow(loc_coefs_unified) > 0) {
            openxlsx::writeData(wb, sheet_name, "COEFFICIENTS",
              startRow = current_row, startCol = 1
            )
            openxlsx::writeData(wb, sheet_name, loc_coefs_unified,
              startRow = current_row + 1, startCol = 1
            )
            current_row <- current_row + nrow(loc_coefs_unified) + 3
          }

          # If uncertainty coefficients exist, write them next
          # COEFFICIENT CI table suppressed: CI bounds have been integrated into the COEFFICIENTS sheet as columns

          # Add MESMA PCA variant trajectory summary
          if (!is.null(loc_variants_pca) && nrow(loc_variants_pca) > 0) {
            # Create summary of variant usage
            variant_usage <- data.frame(
              year = unique(loc_variants_pca$year),
              stringsAsFactors = FALSE
            )

            # Determine veg variants to summarise for this location
            # Prefer the veg list used during fitting but fall back to
            # available Vegetation columns or ALLOWED_VEG if not present.
            veg_candidates <- unique(c(ALLOWED_VEG, na.omit(unique(as.character(all_coefs$Veg)))))
            veg_candidates <- veg_candidates[!is.na(veg_candidates) & veg_candidates != ""]
            veg_kept_local <- veg_candidates
            if (length(veg_kept_local) == 0) veg_kept_local <- ALLOWED_VEG
            for (veg in veg_kept_local) {
              var_col <- paste0(veg, "_variant")
              if (var_col %in% names(loc_variants_pca)) {
                variant_usage[[paste0(veg, "_most_common")]] <- sapply(variant_usage$year, function(y) {
                  year_data <- loc_variants_pca[loc_variants_pca$year == y, var_col]
                  if (length(year_data) > 0) {
                    names(sort(table(year_data), decreasing = TRUE))[1]
                  } else {
                    NA
                  }
                })
              }
            }

            openxlsx::writeData(wb, sheet_name, "VARIANT USAGE",
              startRow = current_row, startCol = 1
            )
            openxlsx::writeData(wb, sheet_name, variant_usage,
              startRow = current_row + 1, startCol = 1
            )
            current_row <- current_row + nrow(variant_usage) + 3
          }

          # If uncertainty variant dominance exists, write it
          if (exists("all_unc_var") && !is.null(all_unc_var)) {
            loc_unc_var <- all_unc_var[all_unc_var$location_id == loc_id, , drop = FALSE]
            if (nrow(loc_unc_var) > 0) {
              openxlsx::writeData(wb, sheet_name, "VARIANT DOMINANCE (%)",
                startRow = current_row, startCol = 1
              )
              openxlsx::writeData(wb, sheet_name, loc_unc_var,
                startRow = current_row + 1, startCol = 1
              )
              current_row <- current_row + nrow(loc_unc_var) + 3
            }
          }

          # Add Q10/Q90 DVI data if available
          if (!is.null(loc_q_data) && nrow(loc_q_data) > 0) {
            openxlsx::writeData(wb, sheet_name, "Q10/Q90 DVI TREND",
              startRow = current_row, startCol = 1
            )
            openxlsx::writeData(wb, sheet_name, loc_q_data,
              startRow = current_row + 1, startCol = 1
            )
            current_row <- current_row + nrow(loc_q_data) + 3
          }

          # If RMSE CI exists for this location
          if (exists("all_unc_rmse") && !is.null(all_unc_rmse)) {
            loc_unc_rmse <- all_unc_rmse[all_unc_rmse$location_id == loc_id, , drop = FALSE]
            if (nrow(loc_unc_rmse) > 0) {
              openxlsx::writeData(wb, sheet_name, "RMSE CI (2.5%/97.5%)",
                startRow = current_row, startCol = 1
              )
              openxlsx::writeData(wb, sheet_name, loc_unc_rmse,
                startRow = current_row + 1, startCol = 1
              )
              current_row <- current_row + nrow(loc_unc_rmse) + 3
            }
          }

          if (!is.null(inseparable_flags_summary) && nrow(inseparable_flags_summary) > 0) {
            loc_flagged <- inseparable_flags_summary[inseparable_flags_summary$location_id == loc_id, , drop = FALSE]
            if (nrow(loc_flagged) > 0) {
              openxlsx::writeData(wb, sheet_name, "INSEPARABLE VARIANT ALERTS (summary)", startRow = current_row, startCol = 1)
              openxlsx::writeData(wb, sheet_name, loc_flagged, startRow = current_row + 1, startCol = 1)
              current_row <- current_row + nrow(loc_flagged) + 3
            }
          }

          if (!is.null(inseparable_flags_detail) && nrow(inseparable_flags_detail) > 0) {
            loc_flagged_detail <- inseparable_flags_detail[inseparable_flags_detail$location_id == loc_id, , drop = FALSE]
            if (nrow(loc_flagged_detail) > 0) {
              openxlsx::writeData(wb, sheet_name, "INSEPARABLE VARIANT COMPONENTS", startRow = current_row, startCol = 1)
              openxlsx::writeData(wb, sheet_name, loc_flagged_detail, startRow = current_row + 1, startCol = 1)
              current_row <- current_row + nrow(loc_flagged_detail) + 3
            }
          }

        }

    # Save the workbook
    output_filename <- file.path(OUT_DIR, "mesma_results.xlsx")

    if (!dir.exists(OUT_DIR)) {
      cat(sprintf("Creating output directory: %s\n", OUT_DIR))
      dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
    }

    if (!dir.exists(OUT_DIR)) {
      cat(sprintf("ERROR: Cannot create output directory: %s\n", OUT_DIR))
      stop("Cannot create output directory")
    }

    cat(sprintf("Saving workbook to: %s\n", output_filename))

    save_result <- tryCatch(
      {
        openxlsx::saveWorkbook(wb, output_filename, overwrite = TRUE)
        TRUE
      },
      error = function(e) {
        cat(sprintf("ERROR saving workbook: %s\n", e$message))
        FALSE
      }
    )

    if (save_result) {
      cat(sprintf(
        "Created Excel file '%s' with %d location sheets\n",
        basename(output_filename), length(unique_locations)
      ))
    } else {
      cat("Failed to save Excel file\n")
    }
  }

  # ============================================================================
  # GLOBAL VEGETATION PATTERN AGGREGATION WITH UNCERTAINTY
  # ============================================================================

  #' Aggregate location-level coefficients to global patterns with uncertainty
  #' 
  #' @param all_coefs Data frame with columns: location_id, year, Veg, coef, coef_025, coef_975
  #' @param method Aggregation method: "hierarchical", "bootstrap" (default: "bootstrap")
  #' @return Data frame with global patterns and confidence intervals per year and vegetation type
  aggregate_to_global_pattern <- function(all_coefs, method = "bootstrap") {
    
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")
    
    # Ensure required columns exist
    required_cols <- c("location_id", "year", "Veg", "coef")
    missing <- setdiff(required_cols, names(all_coefs))
    if (length(missing) > 0) stop(paste("Missing columns:", paste(missing, collapse = ", ")))
    
    # Add CI columns if missing
    if (!"coef_025" %in% names(all_coefs)) all_coefs$coef_025 <- NA_real_
    if (!"coef_975" %in% names(all_coefs)) all_coefs$coef_975 <- NA_real_
    
    # Calculate interval width (proxy for precision)
    all_coefs$interval <- all_coefs$coef_975 - all_coefs$coef_025
    all_coefs$interval[!is.finite(all_coefs$interval)] <- NA_real_
    
    # Calculate precision weight (inverse variance proxy)
    # Use interval width as proxy for standard error: SE ≈ interval / (2 * 1.96)
    all_coefs$se_proxy <- all_coefs$interval / 3.92
    all_coefs$se_proxy[all_coefs$se_proxy <= 0 | !is.finite(all_coefs$se_proxy)] <- NA_real_
    all_coefs$weight <- 1 / (all_coefs$se_proxy^2)
    all_coefs$weight[!is.finite(all_coefs$weight)] <- 1  # Default weight if no CI
    
    # Prefer more robust aggregations: hierarchical (mixed-effects) or bootstrap.
    if (method == "hierarchical") {
      result <- aggregate_hierarchical(all_coefs)
    } else if (method == "bootstrap") {
      result <- aggregate_bootstrap(all_coefs)
    } else {
      stop("Unknown method. Use 'hierarchical' or 'bootstrap'")
    }
    
    result
  }

  # NOTE: 'weighted_mean' aggregation removed by request - use 'hierarchical' or 'bootstrap'.

  #' Hierarchical/Mixed-effects aggregation
  aggregate_hierarchical <- function(all_coefs) {
    
    if (!requireNamespace("lme4", quietly = TRUE)) {
      warning("lme4 not available, falling back to bootstrap aggregation")
      return(aggregate_bootstrap(all_coefs))
    }
    
    veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
    results_list <- list()
    
    for (veg in veg_types) {
      veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), ]
      
      if (nrow(veg_data) < 10) {
        # Not enough data for mixed model, use simple mean
        simple_result <- veg_data %>%
          dplyr::group_by(pheno_year) %>%
          dplyr::summarize(
            n_locations = dplyr::n(),
            global_coef = mean(coef, na.rm = TRUE),
            se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
            ci_lower = global_coef - 1.96 * se,
            ci_upper = global_coef + 1.96 * se,
            .groups = "drop"
          ) %>%
          dplyr::mutate(Veg = veg, method = "simple_mean")
        
        results_list[[veg]] <- simple_result
        next
      }
      
      # Fit mixed-effects model: coef ~ year + (1|location_id)
      # This accounts for repeated measures within locations
      tryCatch({
        # Convert year to factor for categorical treatment
        veg_data$year_factor <- as.factor(veg_data$pheno_year)
        
        model <- lme4::lmer(coef ~ year_factor + (1|location_id), data = veg_data)
        
        # Extract fixed effects (year estimates)
        fe <- lme4::fixef(model)
        vcov_fe <- as.matrix(vcov(model))
        
        # Get predictions for each year
        years <- sort(unique(veg_data$year))
        pred_data <- data.frame(year_factor = as.factor(years))
        
        # Use predict with se.fit if available, otherwise bootstrap
        preds <- predict(model, newdata = pred_data, re.form = NA)
        
        # Bootstrap for confidence intervals
        boot_preds <- lme4::bootMer(model, function(m) {
          predict(m, newdata = pred_data, re.form = NA)
        }, nsim = 100, type = "parametric")
        
        ci_lower <- apply(boot_preds$t, 2, quantile, 0.025)
        ci_upper <- apply(boot_preds$t, 2, quantile, 0.975)
        
        hier_result <- data.frame(
          year = years,
          Veg = veg,
          n_locations = sapply(years, function(y) sum(veg_data$pheno_year == y)),
          global_coef = preds,
          se = apply(boot_preds$t, 2, sd),
          ci_lower = pmax(0, ci_lower),
          ci_upper = pmin(1, ci_upper),
          method = "hierarchical"
        )
        
        results_list[[veg]] <- hier_result
        
      }, error = function(e) {
        warning(sprintf("Hierarchical model failed for %s: %s. Using simple mean.", veg, e$message))
        simple_result <- veg_data %>%
          dplyr::group_by(pheno_year) %>%
          dplyr::summarize(
            n_locations = dplyr::n(),
            global_coef = mean(coef, na.rm = TRUE),
            se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
            ci_lower = pmax(0, global_coef - 1.96 * se),
            ci_upper = pmin(1, global_coef + 1.96 * se),
            .groups = "drop"
          ) %>%
          dplyr::mutate(Veg = veg, method = "simple_mean")
        
        results_list[[veg]] <- simple_result
      })
    }
    
    dplyr::bind_rows(results_list)
  }

  #' Bootstrap aggregation - resample locations
  aggregate_bootstrap <- function(all_coefs, B = 1000, seed = 123) {
    set.seed(seed)
    
    veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
    years <- sort(unique(all_coefs$pheno_year[!is.na(all_coefs$pheno_year)]))
    locations <- unique(all_coefs$location_id)
    n_locs <- length(locations)
    
    results_list <- list()
    
    for (veg in veg_types) {
      veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), ]
      
      # Bootstrap: resample locations with replacement
      boot_means <- matrix(NA_real_, nrow = B, ncol = length(years))
      colnames(boot_means) <- as.character(years)
      
      for (b in seq_len(B)) {
        # Resample locations
        boot_locs <- sample(locations, n_locs, replace = TRUE)
        
        # Get data for resampled locations
        boot_data <- veg_data[veg_data$location_id %in% boot_locs, ]
        
        # Calculate mean per year
        for (i in seq_along(years)) {
          yr <- years[i]
          yr_data <- boot_data[boot_data$year == yr, ]
          if (nrow(yr_data) > 0) {
            boot_means[b, i] <- mean(yr_data$coef, na.rm = TRUE)
          }
        }
      }
      
      # Calculate statistics from bootstrap distribution
      boot_result <- data.frame(
        year = years,
        Veg = veg,
        n_locations = sapply(years, function(y) sum(veg_data$year == y & !is.na(veg_data$coef))),
        global_coef = apply(boot_means, 2, median, na.rm = TRUE),
        se = apply(boot_means, 2, sd, na.rm = TRUE),
        ci_lower = apply(boot_means, 2, quantile, 0.025, na.rm = TRUE),
        ci_upper = apply(boot_means, 2, quantile, 0.975, na.rm = TRUE),
        method = "bootstrap"
      )
      
      # Clamp to [0, 1]
      boot_result$ci_lower <- pmax(0, boot_result$ci_lower)
      boot_result$ci_upper <- pmin(1, boot_result$ci_upper)
      
      results_list[[veg]] <- boot_result
    }
    
    dplyr::bind_rows(results_list)
  }

  # ============================================================================
  # VISUALIZATION OF GLOBAL VEGETATION PATTERNS
  # ============================================================================

  #' Plot global vegetation pattern over time
  #' 
  #' @param global_pattern Data frame from aggregate_to_global_pattern()
  #' @param title Plot title
  #' @param show_ci Whether to show confidence intervals
  plot_global_vegetation_pattern <- function(global_pattern, 
                                              title = "Global Vegetation Composition Over Time",
                                              show_ci = TRUE,
                                              ci_type = "auto") {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    
    # Determine which CI/coef columns to use
    # Prefer 'global_coef' (bootstrap/hierarchical), then 'mean_coef', finally 'weighted_mean_coef' if present
    if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
      if ("ci_lower" %in% names(global_pattern) && "ci_upper" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower
        global_pattern$ci_upper <- global_pattern$ci_upper
      }
    } else if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
      if ("ci_lower_simple" %in% names(global_pattern) && "ci_upper_simple" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower_simple
        global_pattern$ci_upper <- global_pattern$ci_upper_simple
      }
    } else if ("weighted_mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$weighted_mean_coef
      if ("ci_lower_pooled" %in% names(global_pattern) && "ci_upper_pooled" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower_pooled
        global_pattern$ci_upper <- global_pattern$ci_upper_pooled
      }
    }
    
    # Separate barren and vegetation
    global_pattern_veg <- global_pattern[tolower(global_pattern$Veg) != "barren", ]
    global_pattern_barren <- global_pattern[tolower(global_pattern$Veg) == "barren", ]
    
    # Create main vegetation plot
    p <- ggplot2::ggplot(global_pattern_veg, ggplot2::aes(x = year, y = coef, color = Veg, fill = Veg)) +
      ggplot2::geom_line(size = 1.2) +
      ggplot2::geom_point(size = 2)
    
    if (show_ci && "ci_lower" %in% names(global_pattern_veg) && "ci_upper" %in% names(global_pattern_veg)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
        alpha = 0.2,
        color = NA
      )
    }
    
    # Add barren line with adaptive secondary axis
    if (nrow(global_pattern_barren) > 0) {
      # Compute adaptive scale factor based on the visible ranges of vegetation and barren
      veg_max <- suppressWarnings(max(global_pattern_veg$coef, global_pattern_veg$ci_upper, na.rm = TRUE))
      barren_max <- suppressWarnings(max(global_pattern_barren$coef, global_pattern_barren$ci_upper, na.rm = TRUE))
      if (is.na(veg_max) || veg_max <= 0 || is.na(barren_max) || barren_max <= 0) {
        barren_scale_factor <- 1
      } else {
        barren_scale_factor <- veg_max / barren_max
      }
      # Create scaled values for plotting on the primary (vegetation) axis
      global_pattern_barren$coef_scaled <- global_pattern_barren$coef * barren_scale_factor
      
      p <- p + 
        ggplot2::geom_line(data = global_pattern_barren, 
              ggplot2::aes(x = year, y = coef_scaled), 
                          color = "brown", size = 1.2, linetype = "dashed") +
        ggplot2::geom_point(data = global_pattern_barren, 
               ggplot2::aes(x = year, y = coef_scaled), 
                           color = "brown", size = 2)
      
      if (show_ci && "ci_lower" %in% names(global_pattern_barren) && "ci_upper" %in% names(global_pattern_barren)) {
        global_pattern_barren$ci_lower_scaled <- global_pattern_barren$ci_lower * barren_scale_factor
        global_pattern_barren$ci_upper_scaled <- global_pattern_barren$ci_upper * barren_scale_factor
        
        p <- p + ggplot2::geom_ribbon(data = global_pattern_barren,
                                     ggplot2::aes(x = year, ymin = ci_lower_scaled, ymax = ci_upper_scaled),
                                     fill = "brown", alpha = 0.1, color = NA)
      }
    }
    
    p <- p +
      ggplot2::labs(
        title = title,
        subtitle = sprintf("Based on %d locations", max(global_pattern$n_locations, na.rm = TRUE)),
        x = "Year",
        y = "Vegetation Fraction",
        color = "Vegetation Type",
        fill = "Vegetation Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)
      ) +
          ggplot2::scale_y_continuous(
        labels = scales::percent_format(),
        sec.axis = ggplot2::sec_axis(~ . / barren_scale_factor, name = "Barren Fraction", labels = scales::percent_format())
      ) +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::scale_fill_brewer(palette = "Set1")
    
    p
  }

  #' Plot stacked area chart of vegetation composition
  plot_vegetation_stacked_area <- function(global_pattern, 
                                            title = "Global Vegetation Composition Over Time") {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    
    # Use mean_coef or global_coef
    if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
    } else if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
    }
    
    # Normalize to sum to 1 per year
    global_pattern <- global_pattern %>%
      dplyr::group_by(pheno_year) %>%
      dplyr::mutate(coef_normalized = coef / sum(coef, na.rm = TRUE)) %>%
      dplyr::ungroup()
    
    p <- ggplot2::ggplot(global_pattern, 
                ggplot2::aes(x = pheno_year, y = coef_normalized, fill = Veg)) +
      ggplot2::geom_area(alpha = 0.8, position = "stack") +
      ggplot2::labs(
        title = title,
        x = "Year",
        y = "Relative Vegetation Fraction",
        fill = "Vegetation Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold")
      ) +
      ggplot2::scale_y_continuous(labels = scales::percent_format()) +
      ggplot2::scale_fill_brewer(palette = "Set2")
    
    p
  }

  #' Plot heatmap of vegetation changes
  plot_vegetation_heatmap <- function(global_pattern, 
                                       title = "Vegetation Fraction by Year") {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    
    if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
    } else if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
    }
    
    p <- ggplot2::ggplot(global_pattern, 
                          ggplot2::aes(x = year, y = Veg, fill = coef)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", coef * 100)), 
                         color = "white", size = 3) +
      ggplot2::scale_fill_viridis_c(option = "plasma", 
                                     labels = scales::percent_format()) +
      ggplot2::labs(
        title = title,
        x = "Year",
        y = "Vegetation Type",
        fill = "Fraction"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    
    p
  }

  # ============================================================================
  # INFERENCE RESULTS - SEPARATE EXCEL AND PLOT
  # ============================================================================
  if (length(inference_location_ids) > 0) {
    cat("\nProcessing inference results for separate output...\n")
    
    # Filter all_coefs for inference locations
    inference_coefs <- all_coefs[all_coefs$location_id %in% inference_location_ids, ]
    
    if (nrow(inference_coefs) > 0) {
      # Create separate workbook for inference
      wb_inference <- openxlsx::createWorkbook()
      
      # Add sheets for each inference location
      unique_inference_locations <- unique(inference_coefs$location_id)
      for (loc_id in unique_inference_locations) {
        loc_data <- inference_coefs[inference_coefs$location_id == loc_id, ]
        sheet_name <- paste0("Loc_", loc_id)
        openxlsx::addWorksheet(wb_inference, sheet_name)
        openxlsx::writeData(wb_inference, sheet_name, loc_data)
      }
      
      # Save inference workbook
      inference_output_filename <- file.path(OUT_DIR, "inference_results.xlsx")
      openxlsx::saveWorkbook(wb_inference, inference_output_filename, overwrite = TRUE)
      cat(sprintf("Saved inference Excel file to: %s\n", inference_output_filename))
      
      # Create inference plot with bootstrapping
      # Aggregate inference data with bootstrapping for uncertainty
      inference_global_pattern <- aggregate_to_global_pattern(inference_coefs, method = "bootstrap")
      
      if (nrow(inference_global_pattern) > 0) {
        # Ensure plotting packages are available
        if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
        if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")
        
        # Separate barren and vegetation for inference
        inference_global_pattern_veg <- inference_global_pattern[tolower(inference_global_pattern$Veg) != "barren", ]
        inference_global_pattern_barren <- inference_global_pattern[tolower(inference_global_pattern$Veg) == "barren", ]
        
        # Compute adaptive scale between vegetation and barren for the inference plot
        veg_max <- suppressWarnings(max(inference_global_pattern_veg$global_coef, inference_global_pattern_veg$ci_upper, na.rm = TRUE))
        barren_max <- suppressWarnings(max(inference_global_pattern_barren$global_coef, inference_global_pattern_barren$ci_upper, na.rm = TRUE))
        if (is.na(veg_max) || veg_max <= 0 || is.na(barren_max) || barren_max <= 0) {
          inference_barren_scale <- 1
        } else {
          inference_barren_scale <- veg_max / barren_max
        }
        # keep vegetation as-is, scale barren for visualization on primary axis
        if (nrow(inference_global_pattern_barren) > 0) {
          inference_global_pattern_barren$global_coef_scaled <- inference_global_pattern_barren$global_coef * inference_barren_scale
          inference_global_pattern_barren$ci_lower_scaled <- inference_global_pattern_barren$ci_lower * inference_barren_scale
          inference_global_pattern_barren$ci_upper_scaled <- inference_global_pattern_barren$ci_upper * inference_barren_scale
        }
        
        # Create main vegetation plot with CI
        p_inference <- ggplot(inference_global_pattern_veg, aes(x = year, y = global_coef, color = Veg, fill = Veg)) +
          geom_line(size = 1.2) +
          geom_point(size = 2, show.legend = FALSE) +
          geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, color = NA)
        
        # Add barren line on the same axis if present (scaled to primary veg axis)
        if (nrow(inference_global_pattern_barren) > 0) {
          p_inference <- p_inference + 
            geom_line(data = inference_global_pattern_barren, 
                     aes(x = year, y = global_coef_scaled), 
                     color = "brown", size = 1.2, linetype = "dashed") +
            geom_point(data = inference_global_pattern_barren, 
                      aes(x = year, y = global_coef_scaled), 
                      color = "brown", size = 2, show.legend = FALSE) +
            geom_ribbon(data = inference_global_pattern_barren,
                       aes(x = year, ymin = ci_lower_scaled, ymax = ci_upper_scaled),
                       fill = "brown", alpha = 0.1, color = NA)
        }
        
        # Calculate adaptive limits for y-axis
        all_y_values <- c(inference_global_pattern_veg$global_coef, inference_global_pattern_veg$ci_lower, inference_global_pattern_veg$ci_upper)
        if (nrow(inference_global_pattern_barren) > 0) {
          all_y_values <- c(all_y_values, inference_global_pattern_barren$global_coef_scaled, inference_global_pattern_barren$ci_lower_scaled, inference_global_pattern_barren$ci_upper_scaled)
        }
        min_y <- min(all_y_values, na.rm = TRUE)
        max_y <- max(all_y_values, na.rm = TRUE)
        
        p_inference <- p_inference +
          labs(
            title = "Inference Locations: Average Coverage Percentage per Vegetation Type (2020-2025)",
            subtitle = sprintf("Based on %d locations with bootstrap uncertainty", max(inference_global_pattern$n_locations, na.rm = TRUE)),
            x = "Year",
            y = "Vegetation Fraction",
            color = "Vegetation Type",
            fill = "Vegetation Type"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 10),
            axis.title = element_text(size = 12),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10)
          ) +
          scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                             limits = c(min_y, max_y),
                             expand = c(0, 0),
                             sec.axis = ggplot2::sec_axis(~ . / inference_barren_scale, name = "Barren Fraction", labels = scales::percent_format(accuracy = 1))) +
          scale_x_continuous(limits = c(2020, 2025)) +
          scale_color_brewer(palette = "Set1") +
          scale_fill_brewer(palette = "Set1")
        
        # Save the inference plot
        inference_plot_filename <- file.path(OUT_DIR, "inference_average_coverage_plot.png")
        ggsave(inference_plot_filename, p_inference, width = 10, height = 6, dpi = 300)
        cat(sprintf("Saved inference plot to: %s\n", inference_plot_filename))
        
        # PPI-normalized cumulative for inference
        if ("PPI" %in% names(df) && any(!is.na(df$PPI))) {
          ppi_per_loc_year_inf <- df %>%
            filter(location_id %in% inference_location_ids, doy >= 182, doy <= 243) %>%
            group_by(location_id, year) %>%
            summarize(median_ppi = median(PPI, na.rm = TRUE), .groups = "drop")
          
          inference_coefs_ppi <- inference_coefs %>%
            left_join(ppi_per_loc_year_inf, by = c("location_id", "year"))
          
          inference_coefs_ppi <- inference_coefs_ppi %>%
            mutate(normalized_coef = coef * median_ppi)
          
          time_series_per_veg_inf <- inference_coefs_ppi %>%
            group_by(Veg, year) %>%
            summarize(total_normalized = sum(normalized_coef, na.rm = TRUE), .groups = "drop")
          
          # Plot time series for inference
          p_inf_ppi_ts <- ggplot(time_series_per_veg_inf, aes(x = year, y = total_normalized, color = Veg, group = Veg)) +
            geom_line(size = 1) +
            geom_point(show.legend = FALSE) +
            labs(title = "Inference: PPI-Normalized Vegetation Fractions Over Time",
                 x = "Year", y = "Total Normalized Fraction") +
            theme_minimal()
          
          inf_ppi_plot_filename <- file.path(OUT_DIR, "inference_ppi_normalized_timeseries.png")
          ggsave(inf_ppi_plot_filename, p_inf_ppi_ts, width = 8, height = 6)
          cat(sprintf("Saved inference PPI-normalized time series plot to: %s\n", inf_ppi_plot_filename))
        }
      }
    } else {
      cat("No inference data found for specified location IDs.\n")
    }
  } else {
    cat("No inference location IDs specified.\n")
  }

  # ============================================================================
  # GENERATE AVERAGE COVERAGE PLOT (2000-2024)
  # ============================================================================
  cat("\nGenerating average coverage plot...\n")

  # Load ggplot2 if not already loaded
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  library(ggplot2)

  # Aggregate average coverage with bootstrapping for uncertainty across all locations
  global_pattern_all <- aggregate_to_global_pattern(all_coefs, method = "bootstrap")
  
  # Filter to 1985-2025 and separate barren/vegetation
  global_pattern_all <- global_pattern_all %>% filter(year >= 1985 & year <= 2025)
  global_pattern_all_veg <- global_pattern_all[tolower(global_pattern_all$Veg) != "barren", ]
  global_pattern_all_barren <- global_pattern_all[tolower(global_pattern_all$Veg) == "barren", ]

  # Create the plot with bootstrapping
  p <- ggplot(global_pattern_all_veg, aes(x = year, y = global_coef, color = Veg, fill = Veg)) +
    geom_line(size = 1.2) +
    geom_point(size = 2, show.legend = FALSE) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, color = NA)
  
  # Add barren line with adaptive secondary axis if present
  barren_scale_factor <- 1
  if (nrow(global_pattern_all_barren) > 0) {
    veg_max <- suppressWarnings(max(global_pattern_all_veg$global_coef, global_pattern_all_veg$ci_upper, na.rm = TRUE))
    barren_max <- suppressWarnings(max(global_pattern_all_barren$global_coef, global_pattern_all_barren$ci_upper, na.rm = TRUE))
    if (is.na(veg_max) || veg_max <= 0 || is.na(barren_max) || barren_max <= 0) {
      barren_scale_factor <- 1
    } else {
      barren_scale_factor <- veg_max / barren_max
    }
    global_pattern_all_barren$coef_scaled <- global_pattern_all_barren$global_coef * barren_scale_factor
    global_pattern_all_barren$ci_lower_scaled <- global_pattern_all_barren$ci_lower * barren_scale_factor
    global_pattern_all_barren$ci_upper_scaled <- global_pattern_all_barren$ci_upper * barren_scale_factor
    
    p <- p + 
      geom_line(data = global_pattern_all_barren, 
           aes(x = year, y = coef_scaled), 
               color = "brown", size = 1.2, linetype = "dashed") +
      geom_point(data = global_pattern_all_barren, 
            aes(x = year, y = coef_scaled), 
                color = "brown", size = 2, show.legend = FALSE) +
      geom_ribbon(data = global_pattern_all_barren,
                 aes(x = year, ymin = ci_lower_scaled, ymax = ci_upper_scaled),
                 fill = "brown", alpha = 0.1, color = NA)
  }
  
  # Calculate adaptive limits for y-axis
  all_y_values <- c(global_pattern_all_veg$global_coef, global_pattern_all_veg$ci_lower, global_pattern_all_veg$ci_upper)
  if (exists("global_pattern_all_barren") && nrow(global_pattern_all_barren) > 0) {
    all_y_values <- c(all_y_values, global_pattern_all_barren$coef_scaled, global_pattern_all_barren$ci_lower_scaled, global_pattern_all_barren$ci_upper_scaled)
  }
  min_y <- min(all_y_values, na.rm = TRUE)
  max_y <- max(all_y_values, na.rm = TRUE)
  
  p <- p +
    labs(
      title = "Average Coverage Percentage per Vegetation Type (2000-2024)",
      subtitle = sprintf("Based on %d locations with bootstrap uncertainty", max(global_pattern_all$n_locations, na.rm = TRUE)),
      x = "Year",
      y = "Vegetation Fraction",
      color = "Vegetation Type",
      fill = "Vegetation Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(min_y, max_y),
      expand = c(0, 0),
      sec.axis = sec_axis(~ . / barren_scale_factor, name = "Barren Fraction", labels = scales::percent_format(accuracy = 1))
    ) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1")

  # Save the plot
  plot_filename <- file.path(OUT_DIR, "average_coverage_plot.png")
  ggsave(plot_filename, p, width = 10, height = 6, dpi = 300)
  cat(sprintf("Saved average coverage plot to: %s\n", plot_filename))

  # ============================================================================
  # PLOT: #OBSERVATIONS vs PERCENTAGE CORRECTLY PREDICTED
  # ============================================================================
  cat("\nGenerating Observations vs Accuracy plot...\n")
  
  # Prepare data: Join summary_data (Total_Observations) with best_fit_summary (accuracy)
  # best_fit_summary has one row per location-year. We average accuracy per location.
  
  loc_accuracy <- best_fit_summary %>%
    group_by(location_id) %>%
    summarize(
      mean_pred_coef_rel = mean(pred_coef_rel, na.rm = TRUE),
      mean_pred_coef_abs = mean(pred_coef_abs, na.rm = TRUE),
      .groups = "drop"
    )
  
  obs_vs_acc_data <- summary_data %>%
    left_join(loc_accuracy, by = c("Location_ID" = "location_id"))
  
  # Filter out locations with no accuracy data (e.g. barren locations if we only look at veg accuracy)
  obs_vs_acc_data <- obs_vs_acc_data %>% filter(!is.na(mean_pred_coef_rel))
  
  if (nrow(obs_vs_acc_data) > 0) {
    p_obs_acc <- ggplot(obs_vs_acc_data, aes(x = Total_Observations, y = mean_pred_coef_rel)) +
      geom_point(alpha = 0.6, color = "darkblue") +
      geom_smooth(method = "lm", color = "red", se = FALSE) +
      labs(
        title = "Number of Observations vs. Vegetation Prediction Accuracy",
        subtitle = "Accuracy = Mean Predicted Fraction of True Vegetation Class (Relative to Total Vegetation)",
        x = "Total Observations (All Years)",
        y = "Mean Correct Prediction Fraction"
      ) +
      theme_minimal() +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1))
    
    ggsave(file.path(OUT_DIR, "observations_vs_accuracy.png"), p_obs_acc, width = 8, height = 6)
    cat(sprintf("Saved Observations vs Accuracy plot to: %s\n", file.path(OUT_DIR, "observations_vs_accuracy.png")))
  } else {
    cat("No data available for Observations vs Accuracy plot.\n")
  }

  # ============================================================================
  # CLUSTER SENSITIVITY ANALYSIS (REMOVED)
  # ============================================================================
  # Cluster sensitivity analysis has been removed as per user request.
  # Default cluster size is now fixed at 10 variants per vegetation type.


  # ============================================================================
  # GLOBAL VEGETATION PATTERN AGGREGATION WITH UNCERTAINTY
  # ============================================================================

  #' Aggregate location-level coefficients to global patterns with uncertainty
  #' 
  #' @param all_coefs Data frame with columns: location_id, year, Veg, coef, coef_025, coef_975
  #' @param method Aggregation method: "hierarchical", "bootstrap" (default: "bootstrap")
  #' @return Data frame with global patterns and confidence intervals per year and vegetation type
  aggregate_to_global_pattern <- function(all_coefs, method = "bootstrap") {
    
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")
    
    # Ensure required columns exist
    required_cols <- c("location_id", "year", "Veg", "coef")
    missing <- setdiff(required_cols, names(all_coefs))
    if (length(missing) > 0) stop(paste("Missing columns:", paste(missing, collapse = ", ")))
    
    # Add CI columns if missing
    if (!"coef_025" %in% names(all_coefs)) all_coefs$coef_025 <- NA_real_
    if (!"coef_975" %in% names(all_coefs)) all_coefs$coef_975 <- NA_real_
    
    # Calculate interval width (proxy for precision)
    all_coefs$interval <- all_coefs$coef_975 - all_coefs$coef_025
    all_coefs$interval[!is.finite(all_coefs$interval)] <- NA_real_
    
    # Calculate precision weight (inverse variance proxy)
    # Use interval width as proxy for standard error: SE ≈ interval / (2 * 1.96)
    all_coefs$se_proxy <- all_coefs$interval / 3.92
    all_coefs$se_proxy[all_coefs$se_proxy <= 0 | !is.finite(all_coefs$se_proxy)] <- NA_real_
    all_coefs$weight <- 1 / (all_coefs$se_proxy^2)
    all_coefs$weight[!is.finite(all_coefs$weight)] <- 1  # Default weight if no CI
    
    # Prefer more robust aggregations: hierarchical (mixed-effects) or bootstrap.
    if (method == "hierarchical") {
      result <- aggregate_hierarchical(all_coefs)
    } else if (method == "bootstrap") {
      result <- aggregate_bootstrap(all_coefs)
    } else {
      stop("Unknown method. Use 'hierarchical' or 'bootstrap'")
    }
    
    result
  }

  # NOTE: 'weighted_mean' aggregation removed by request - use 'hierarchical' or 'bootstrap'.

  #' Hierarchical/Mixed-effects aggregation
  aggregate_hierarchical <- function(all_coefs) {
    
    if (!requireNamespace("lme4", quietly = TRUE)) {
      warning("lme4 not available, falling back to bootstrap aggregation")
      return(aggregate_bootstrap(all_coefs))
    }
    
    veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
    results_list <- list()
    
    for (veg in veg_types) {
      veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), ]
      
      if (nrow(veg_data) < 10) {
        # Not enough data for mixed model, use simple mean
        simple_result <- veg_data %>%
          dplyr::group_by(pheno_year) %>%
          dplyr::summarize(
            n_locations = dplyr::n(),
            global_coef = mean(coef, na.rm = TRUE),
            se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
            ci_lower = global_coef - 1.96 * se,
            ci_upper = global_coef + 1.96 * se,
            .groups = "drop"
          ) %>%
          dplyr::mutate(Veg = veg, method = "simple_mean")
        
        results_list[[veg]] <- simple_result
        next
      }
      
      # Fit mixed-effects model: coef ~ year + (1|location_id)
      # This accounts for repeated measures within locations
      tryCatch({
        # Convert year to factor for categorical treatment
        veg_data$year_factor <- as.factor(veg_data$pheno_year)
        
        model <- lme4::lmer(coef ~ year_factor + (1|location_id), data = veg_data)
        
        # Extract fixed effects (year estimates)
        fe <- lme4::fixef(model)
        vcov_fe <- as.matrix(vcov(model))
        
        # Get predictions for each year
        years <- sort(unique(veg_data$year))
        pred_data <- data.frame(year_factor = as.factor(years))
        
        # Use predict with se.fit if available, otherwise bootstrap
        preds <- predict(model, newdata = pred_data, re.form = NA)
        
        # Bootstrap for confidence intervals
        boot_preds <- lme4::bootMer(model, function(m) {
          predict(m, newdata = pred_data, re.form = NA)
        }, nsim = 100, type = "parametric")
        
        ci_lower <- apply(boot_preds$t, 2, quantile, 0.025)
        ci_upper <- apply(boot_preds$t, 2, quantile, 0.975)
        
        hier_result <- data.frame(
          year = years,
          Veg = veg,
          n_locations = sapply(years, function(y) sum(veg_data$pheno_year == y)),
          global_coef = preds,
          se = apply(boot_preds$t, 2, sd),
          ci_lower = pmax(0, ci_lower),
          ci_upper = pmin(1, ci_upper),
          method = "hierarchical"
        )
        
        results_list[[veg]] <- hier_result
        
      }, error = function(e) {
        warning(sprintf("Hierarchical model failed for %s: %s. Using simple mean.", veg, e$message))
        simple_result <- veg_data %>%
          dplyr::group_by(pheno_year) %>%
          dplyr::summarize(
            n_locations = dplyr::n(),
            global_coef = mean(coef, na.rm = TRUE),
            se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
            ci_lower = pmax(0, global_coef - 1.96 * se),
            ci_upper = pmin(1, global_coef + 1.96 * se),
            .groups = "drop"
          ) %>%
          dplyr::mutate(Veg = veg, method = "simple_mean")
        
        results_list[[veg]] <- simple_result
      })
    }
    
    dplyr::bind_rows(results_list)
  }

  #' Bootstrap aggregation - resample locations
  aggregate_bootstrap <- function(all_coefs, B = 1000, seed = 123) {
    set.seed(seed)
    
    veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
    years <- sort(unique(all_coefs$year[!is.na(all_coefs$year)]))
    locations <- unique(all_coefs$location_id)
    n_locs <- length(locations)
    
    results_list <- list()
    
    for (veg in veg_types) {
      veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), ]
      
      # Bootstrap: resample locations with replacement
      boot_means <- matrix(NA_real_, nrow = B, ncol = length(years))
      colnames(boot_means) <- as.character(years)
      
      for (b in seq_len(B)) {
        # Resample locations
        boot_locs <- sample(locations, n_locs, replace = TRUE)
        
        # Get data for resampled locations
        boot_data <- veg_data[veg_data$location_id %in% boot_locs, ]
        
        # Calculate mean per year
        for (i in seq_along(years)) {
          yr <- years[i]
          yr_data <- boot_data[boot_data$year == yr, ]
          if (nrow(yr_data) > 0) {
            boot_means[b, i] <- mean(yr_data$coef, na.rm = TRUE)
          }
        }
      }
      
      # Calculate statistics from bootstrap distribution
      boot_result <- data.frame(
        year = years,
        Veg = veg,
        n_locations = sapply(years, function(y) sum(veg_data$year == y & !is.na(veg_data$coef))),
        global_coef = apply(boot_means, 2, median, na.rm = TRUE),
        se = apply(boot_means, 2, sd, na.rm = TRUE),
        ci_lower = apply(boot_means, 2, quantile, 0.025, na.rm = TRUE),
        ci_upper = apply(boot_means, 2, quantile, 0.975, na.rm = TRUE),
        method = "bootstrap"
      )
      
      # Clamp to [0, 1]
      boot_result$ci_lower <- pmax(0, boot_result$ci_lower)
      boot_result$ci_upper <- pmin(1, boot_result$ci_upper)
      
      results_list[[veg]] <- boot_result
    }
    
    dplyr::bind_rows(results_list)
  }

  # ============================================================================
  # VISUALIZATION OF GLOBAL VEGETATION PATTERNS
  # ============================================================================

  #' Plot global vegetation pattern over time
  #' 
  #' @param global_pattern Data frame from aggregate_to_global_pattern()
  #' @param title Plot title
  #' @param show_ci Whether to show confidence intervals
  plot_global_vegetation_pattern <- function(global_pattern, 
                                              title = "Global Vegetation Composition Over Time",
                                              show_ci = TRUE,
                                              ci_type = "auto") {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    
    # Determine which CI/coef columns to use
    # Prefer 'global_coef' (bootstrap/hierarchical), then 'mean_coef', finally 'weighted_mean_coef' if present
    if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
      if ("ci_lower" %in% names(global_pattern) && "ci_upper" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower
        global_pattern$ci_upper <- global_pattern$ci_upper
      }
    } else if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
      if ("ci_lower_simple" %in% names(global_pattern) && "ci_upper_simple" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower_simple
        global_pattern$ci_upper <- global_pattern$ci_upper_simple
      }
    } else if ("weighted_mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$weighted_mean_coef
      if ("ci_lower_pooled" %in% names(global_pattern) && "ci_upper_pooled" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower_pooled
        global_pattern$ci_upper <- global_pattern$ci_upper_pooled
      }
    }
    
    # Separate barren and vegetation
    global_pattern_veg <- global_pattern[tolower(global_pattern$Veg) != "barren", ]
    global_pattern_barren <- global_pattern[tolower(global_pattern$Veg) == "barren", ]
    
    # Create main vegetation plot
    p <- ggplot2::ggplot(global_pattern_veg, ggplot2::aes(x = year, y = coef, color = Veg, fill = Veg)) +
      ggplot2::geom_line(size = 1.2) +
      ggplot2::geom_point(size = 2)
    
    if (show_ci && "ci_lower" %in% names(global_pattern_veg) && "ci_upper" %in% names(global_pattern_veg)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
        alpha = 0.2,
        color = NA
      )
    }
    
    # Add barren line with adaptive secondary axis
    barren_scale_factor <- 1
    if (nrow(global_pattern_barren) > 0) {
      veg_max <- suppressWarnings(max(global_pattern_veg$coef, global_pattern_veg$ci_upper, na.rm = TRUE))
      barren_max <- suppressWarnings(max(global_pattern_barren$coef, global_pattern_barren$ci_upper, na.rm = TRUE))
      if (is.na(veg_max) || veg_max <= 0 || is.na(barren_max) || barren_max <= 0) {
        barren_scale_factor <- 1
      } else {
        barren_scale_factor <- veg_max / barren_max
      }
      global_pattern_barren$coef_scaled <- global_pattern_barren$coef * barren_scale_factor
      
      p <- p + 
        ggplot2::geom_line(data = global_pattern_barren, 
                          ggplot2::aes(x = year, y = coef_scaled), 
                          color = "brown", size = 1.2, linetype = "dashed") +
        ggplot2::geom_point(data = global_pattern_barren, 
                           ggplot2::aes(x = year, y = coef_scaled), 
                           color = "brown", size = 2)
      
      if (show_ci && "ci_lower" %in% names(global_pattern_barren) && "ci_upper" %in% names(global_pattern_barren)) {
        global_pattern_barren$ci_lower_scaled <- global_pattern_barren$ci_lower * barren_scale_factor
        global_pattern_barren$ci_upper_scaled <- global_pattern_barren$ci_upper * barren_scale_factor
        
        p <- p + ggplot2::geom_ribbon(data = global_pattern_barren,
                                     ggplot2::aes(x = year, ymin = ci_lower_scaled, ymax = ci_upper_scaled),
                                     fill = "brown", alpha = 0.1, color = NA)
      }
    }
    
    p <- p +
      ggplot2::labs(
        title = title,
        subtitle = sprintf("Based on %d locations", max(global_pattern$n_locations, na.rm = TRUE)),
        x = "Year",
        y = "Vegetation Fraction",
        color = "Vegetation Type",
        fill = "Vegetation Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)
      ) +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(),
        sec.axis = ggplot2::sec_axis(~ . / barren_scale_factor, name = "Barren Fraction", labels = scales::percent_format())
      ) +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::scale_fill_brewer(palette = "Set1")
    
    p
  }

  #' Plot stacked area chart of vegetation composition
  plot_vegetation_stacked_area <- function(global_pattern, 
                                            title = "Global Vegetation Composition Over Time") {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    
    # Use mean_coef or global_coef
    if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
    } else if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
    }
    
    # Normalize to sum to 1 per year
    global_pattern <- global_pattern %>%
      dplyr::group_by(pheno_year) %>%
      dplyr::mutate(coef_normalized = coef / sum(coef, na.rm = TRUE)) %>%
      dplyr::ungroup()
    
    p <- ggplot2::ggplot(global_pattern, 
                ggplot2::aes(x = pheno_year, y = coef_normalized, fill = Veg)) +
      ggplot2::geom_area(alpha = 0.8, position = "stack") +
      ggplot2::labs(
        title = title,
        x = "Year",
        y = "Relative Vegetation Fraction",
        fill = "Vegetation Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold")
      ) +
      ggplot2::scale_y_continuous(labels = scales::percent_format()) +
      ggplot2::scale_fill_brewer(palette = "Set2")
    
    p
  }

  #' Plot heatmap of vegetation changes
  plot_vegetation_heatmap <- function(global_pattern, 
                                       title = "Vegetation Fraction by Year") {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    
    if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
    } else if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
    }
    
    p <- ggplot2::ggplot(global_pattern, 
                          ggplot2::aes(x = year, y = Veg, fill = coef)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", coef * 100)), 
                         color = "white", size = 3) +
      ggplot2::scale_fill_viridis_c(option = "plasma", 
                                     labels = scales::percent_format()) +
      ggplot2::labs(
        title = title,
        x = "Year",
        y = "Vegetation Type",
        fill = "Fraction"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    
    p
  }

  # ============================================================================
  # TREND ANALYSIS FOR VEGETATION CHANGES
  # ============================================================================

  #' Analyze trends in global vegetation patterns with uncertainty
  #' 
  #' @param all_coefs Data frame with location-level coefficients
  #' @return Data frame with trend statistics per vegetation type
  analyze_vegetation_trends <- function(all_coefs) {
    
    veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
    trend_results <- list()
    
    for (veg in veg_types) {
      veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), ]
      
      # Get unique locations
      locations <- unique(veg_data$location_id)
      
      # Compute slope for each location
      loc_slopes <- sapply(locations, function(loc) {
      
      loc_data <- veg_data[veg_data$location_id == loc, ]
      
      if (nrow(loc_data) < 3) return(NA_real_)
      
      # Expand yearly coefficients to weekly timeseries (52 weeks per year)
      n_years <- nrow(loc_data)
      weekly_data <- data.frame(
        week = 1:(n_years * 52),
        coef = rep(loc_data$coef, each = 52)
      )
      
      # Fit GLS trend on weekly data
      model <- tryCatch(nlme::gls(coef ~ week, data = weekly_data, correlation = corAR1()), error = function(e) lm(coef ~ week, data = weekly_data))
      
      fitted_vals <- fitted(model)
      resids <- residuals(model)
      n <- length(resids)
      block_length <- min(52, n)
      if (block_length >= 2) {
        n_blocks <- n - block_length + 1
        blocks <- lapply(1:n_blocks, function(i) resids[i:(i+block_length-1)])
        boot_slopes <- replicate(200, {
          sampled_blocks <- sample(blocks, n_blocks, replace = TRUE)
          boot_resids <- unlist(sampled_blocks)[1:n]
          boot_coef <- fitted_vals + boot_resids
          boot_model <- tryCatch(nlme::gls(boot_coef ~ weekly_data$week, correlation = corAR1()), error = function(e) lm(boot_coef ~ weekly_data$week))
          coef(boot_model)["week"]
        })
      } else {
        boot_slopes <- replicate(200, {
          boot_idx <- sample(n, replace = TRUE)
          boot_coef <- fitted_vals + resids[boot_idx]
          boot_model <- tryCatch(nlme::gls(boot_coef ~ weekly_data$week, correlation = corAR1()), error = function(e) lm(boot_coef ~ weekly_data$week))
          coef(boot_model)["week"]
        })
      }
      median(boot_slopes, na.rm = TRUE)
    })

    # Aggregate location-level slopes for this vegetation type
    slope_mean <- mean(loc_slopes, na.rm = TRUE)
    slope_median <- median(loc_slopes, na.rm = TRUE)
    slope_ci_lower <- quantile(loc_slopes, 0.025, na.rm = TRUE)
    slope_ci_upper <- quantile(loc_slopes, 0.975, na.rm = TRUE)
    prob_positive <- mean(loc_slopes > 0, na.rm = TRUE)
    prob_negative <- mean(loc_slopes < 0, na.rm = TRUE)

    trend_results[[veg]] <- data.frame(
      Veg = veg,
      slope_mean = slope_mean,
      slope_median = slope_median,
      slope_ci_lower = slope_ci_lower,
      slope_ci_upper = slope_ci_upper,
      prob_positive = prob_positive,
      prob_negative = prob_negative,
      stringsAsFactors = FALSE
    )
  }
  
  dplyr::bind_rows(trend_results)
}
  

  #' GLS Moving Block Bootstrap trend uncertainty
  bootstrap_trend_ci <- function(all_coefs, B = 1000, block_length = 3, seed = 123) {
    set.seed(seed)
    
    veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
    
    results_list <- list()
    
    for (veg in veg_types) {
      veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), ]
      
      yearly_data <- veg_data %>%
        group_by(pheno_year) %>%
        summarize(global_coef = mean(coef, na.rm = TRUE), .groups = "drop") %>%
        dplyr::rename(year = pheno_year)
      
      if (nrow(yearly_data) < 5) {
        results_list[[veg]] <- data.frame(Veg = veg, slope_mean = NA_real_, slope_median = NA_real_, slope_ci_lower = NA_real_, slope_ci_upper = NA_real_, prob_positive = NA_real_, prob_negative = NA_real_)
        next
      }
      
      model <- tryCatch(nlme::gls(global_coef ~ year, data = yearly_data, correlation = corAR1()), error = function(e) lm(global_coef ~ year, data = yearly_data))
      
      fitted_vals <- fitted(model)
      resids <- residuals(model)
      
      n <- length(resids)
      if (block_length >= n) block_length <- n
      n_blocks <- n - block_length + 1
      blocks <- lapply(1:n_blocks, function(i) resids[i:(i+block_length-1)])
      
      boot_slopes <- rep(NA_real_, B)
      for (b in 1:B) {
        sampled_blocks <- sample(blocks, n_blocks, replace = TRUE)
        boot_resids <- unlist(sampled_blocks)[1:n]
        boot_coef <- fitted_vals + boot_resids
        boot_model <- tryCatch(nlme::gls(boot_coef ~ yearly_data$year, correlation = corAR1()), error = function(e) lm(boot_coef ~ yearly_data$year))
        boot_slopes[b] <- coef(boot_model)[if (inherits(boot_model, "gls")) "yearly_data$year" else "yearly_data$year"]
      }
      
      results_list[[veg]] <- data.frame(Veg = veg, slope_mean = mean(boot_slopes, na.rm = TRUE), slope_median = median(boot_slopes, na.rm = TRUE), slope_ci_lower = quantile(boot_slopes, 0.025, na.rm = TRUE), slope_ci_upper = quantile(boot_slopes, 0.975, na.rm = TRUE), prob_positive = mean(boot_slopes > 0, na.rm = TRUE), prob_negative = mean(boot_slopes < 0, na.rm = TRUE))
    }
    
    dplyr::bind_rows(results_list)
  }

  # ============================================================================
  # EXAMPLE USAGE
  # ============================================================================

  # After your MESMA fitting completes and you have all_coefs:

  # 1. Aggregate to global pattern (use Bootstrap by default; hierarchical also supported)
  global_pattern <- aggregate_to_global_pattern(all_coefs, method = "bootstrap")
  global_pattern_bootstrap <- global_pattern

  # 2. Visualize - SINGLE FIGURE WITH ALL VEGETATION TYPES ON SAME SCALE
  # Create single plot with all vegetation types
  p1_all <- plot_global_vegetation_pattern(global_pattern, 
                                          title = "All Vegetation Trends",
                                          show_ci = TRUE)
  p1_all <- p1_all + scale_y_continuous(labels = scales::percent_format()) + labs(y = "Fraction")
  ggsave(file.path(OUT_DIR, "all_vegetation_trends.png"), p1_all, width = 10, height = 6)



  # Stacked area and heatmap (only for vegetation, excluding barren)
  p3 <- plot_vegetation_stacked_area(global_pattern)
  ggsave(file.path(OUT_DIR, "vegetation_stacked_area.png"), p3, width = 10, height = 6)

  p4 <- plot_vegetation_heatmap(global_pattern)
  ggsave(file.path(OUT_DIR, "vegetation_heatmap.png"), p4, width = 10, height = 6)

  # 3. Analyze trends
  trends <- analyze_vegetation_trends(all_coefs)
  print(trends)

  trend_ci <- bootstrap_trend_ci(all_coefs, B = 1000, block_length = 3)
  print(trend_ci)

  # PPI-normalized cumulative per vegetation type
  cat("\nGenerating PPI-normalized cumulative plot...\n")
  if ("PPI" %in% names(df) && any(!is.na(df$PPI))) {
    ppi_per_loc_year <- df %>%
      filter(doy >= 182, doy <= 243) %>%
      mutate(pheno_year = ifelse(lubridate::month(date) >= 3, lubridate::year(date), lubridate::year(date) - 1)) %>%
      group_by(location_id, pheno_year) %>%
      summarize(median_ppi = median(PPI, na.rm = TRUE), .groups = "drop")
    
    all_coefs_ppi <- all_coefs %>%
      left_join(ppi_per_loc_year, by = c("location_id", "pheno_year"))
    
    all_coefs_ppi <- all_coefs_ppi %>%
      mutate(normalized_coef = coef * median_ppi)
    
    time_series_per_veg <- all_coefs_ppi %>%
      group_by(Veg, pheno_year) %>%
      summarize(total_normalized = sum(normalized_coef, na.rm = TRUE), .groups = "drop")
    
    # Plot time series
    p_ppi_ts <- ggplot(time_series_per_veg, aes(x = pheno_year, y = total_normalized, color = Veg, group = Veg)) +
      geom_line(size = 1) +
      geom_point(show.legend = FALSE) +
      labs(title = "PPI-Normalized Vegetation Fractions Over Time",
           x = "Year", y = "Total Normalized Fraction") +
      theme_minimal()
    
    ggsave(file.path(OUT_DIR, "ppi_normalized_timeseries.png"), p_ppi_ts, width = 8, height = 6)
    cat(sprintf("Saved PPI-normalized time series plot to: %s\n", file.path(OUT_DIR, "ppi_normalized_timeseries.png")))
  } else {
    cat("PPI data not available, skipping PPI-normalized plot.\n")
  }

  # 4. Save results to Excel
  openxlsx::addWorksheet(wb, "Global_Pattern")
  openxlsx::writeData(wb, "Global_Pattern", global_pattern)

  openxlsx::addWorksheet(wb, "Vegetation_Trends")
  openxlsx::writeData(wb, "Vegetation_Trends", trends)

  openxlsx::addWorksheet(wb, "Trend_Bootstrap_CI")
  openxlsx::writeData(wb, "Trend_Bootstrap_CI", trend_ci)

  # Final timing summary
  timing_info$end_time <- Sys.time()
  total_time <- as.numeric(difftime(timing_info$end_time, timing_info$start_time, units = "secs"))

  cat(sprintf("\nTotal execution time: %.1f seconds (%.1f minutes)\n", total_time, total_time / 60))
  if (!is.null(timing_info$moving_var_done)) {
    cat(sprintf(
      "Moving variance: %.1f seconds\n",
      as.numeric(difftime(timing_info$moving_var_done, timing_info$start_time, units = "secs"))
    ))
  }
  if (!is.null(timing_info$lib_construction_done)) {
    cat(sprintf(
      "Library construction: %.1f seconds\n",
      as.numeric(difftime(timing_info$lib_construction_done, timing_info$moving_var_done, units = "secs"))
    ))
  }
  if (!is.null(timing_info$pca_computation_done)) {
    cat(sprintf(
      "PCA computation: %.1f seconds\n",
      as.numeric(difftime(timing_info$pca_computation_done, timing_info$lib_construction_done, units = "secs"))
    ))
  }
  cat(sprintf(
    "Main processing + Excel: %.1f seconds\n",
    as.numeric(difftime(timing_info$end_time, timing_info$pca_computation_done, units = "secs"))
  ))

  cat("\nMESMA fitting completed successfully!\n")

print("Script fit_veg_mixture_mesma.R execution finished.")






}
