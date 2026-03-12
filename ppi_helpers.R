## Centralized PPI helpers for consistent calculation
## Provides: safe_as_numeric, calculate_solar_zenith, ppi, add_ppi_columns, auto_add_ppi_columns
## Also provides shared utility functions: make_location_id, assign_pheno_year, pheno_doy

# Shared utilities moved to `mesma_helpers.R` — canonical definitions live there.
if (!exists("safe_as_numeric") && file.exists("mesma_helpers.R")) source("mesma_helpers.R")

# Provide basic phenology and location helpers for all MESMA scripts.  Some of the
# older lightweight pipelines (e.g. fit_veg_svm_pca_lda.R or
# january_averages.R) re‑implement these locally; defining them here ensures
# they are available as soon as ppi_helpers.R is sourced.  We also guard with
# exists() checks so callers can override or re‑define for testing.
if (!exists("assign_pheno_year")) {
  # phenological year starts March 1st; Jan/Feb belong to the previous year.
  assign_pheno_year <- function(d) {
    d <- tryCatch(as.Date(d), error = function(e) NA)
    ifelse(is.na(d), NA_integer_,
           ifelse(lubridate::month(d) >= 3,
                  lubridate::year(d),
                  lubridate::year(d) - 1))
  }
}

if (!exists("pheno_doy")) {
  # day-of-year adjusted so that March 1 = 1
  pheno_doy <- function(d) {
    d <- tryCatch(as.Date(d), error = function(e) NA)
    month <- lubridate::month(d)
    ifelse(is.na(d), NA_integer_,
           ifelse(month >= 3,
                  as.integer(d - as.Date(paste0(lubridate::year(d), "-03-01"))) + 1L,
                  as.integer(d - as.Date(paste0(lubridate::year(d) - 1, "-03-01"))) + 1L))
  }
}

if (!exists("make_location_id")) {
  # create a deterministic string ID from lon/lat pairs
  make_location_id <- function(lon, lat) {
    lon <- as.numeric(lon)
    lat <- as.numeric(lat)
    if (length(lon) == 1 && length(lat) == 1) {
      if (!is.finite(lon) || !is.finite(lat)) return(NA_character_)
      sprintf("L_%0.6f_%0.6f", round(lat, 6), round(lon, 6))
    } else {
      res <- rep(NA_character_, length(lon))
      valid <- is.finite(lon) & is.finite(lat)
      if (any(valid)) {
        res[valid] <- sprintf("L_%0.6f_%0.6f", round(lat[valid], 6), round(lon[valid], 6))
      }
      res
    }
  }
}

# PPI-specific helpers continue below.
calculate_solar_zenith <- function(lat, doy, hour = 10.5) {
  # This function is CORRECT. It addresses the SZA issue.
  lat_rad <- lat * pi / 180
  dec_rad <- 23.45 * sin(2 * pi * (284 + doy) / 365) * pi / 180
  h_rad <- (hour - 12) * 15 * pi / 180
  cos_z <- sin(lat_rad) * sin(dec_rad) + cos(lat_rad) * cos(dec_rad) * cos(h_rad)
  acos(pmin(pmax(cos_z, -1), 1))
}

# Default DVI soil baseline to use when no barren-derived baseline or explicit
# parameter is provided. This value is a stable choice used in several experiments.
DEFAULT_DVI_SOIL <- 0.0308

ppi <- function(dvi, zenith.angle, M, dvi.soil, G = 0.5){
  # Strict parameter requirements: M and dvi.soil must be provided and finite.
  if (missing(M) || is.null(M) || any(!is.finite(M))) {
    stop("[PPI ERROR] M parameter must be provided and finite. No automatic fallback allowed.")
  }
  if (missing(dvi.soil) || is.null(dvi.soil) || any(!is.finite(dvi.soil))) {
    stop("[PPI ERROR] dvi.soil parameter is required and must be finite.")
  }

  # Ensure M sits above dvi.soil to avoid zero/negative denominators
  min_dsoil <- min(dvi.soil, na.rm = TRUE)
  if (is.finite(min_dsoil) && any(is.finite(M)) && any(M <= min_dsoil)) {
    stop(sprintf("[PPI ERROR] Provided M (<= %.6f) must be greater than dvi.soil baseline (%.6f).", min(M, na.rm = TRUE), min_dsoil))
  }
  
  d_c <- 0.0336 + 0.0477/cos(zenith.angle)
  Q_E <- d_c + (1 - d_c) * G / cos(zenith.angle)
  K <- 1/(4*Q_E) * (1 + M)/(1 - M)

  denom <- M - dvi
  numer <- M - dvi.soil
  
  # Avoid division by zero
  invalid_denom <- (abs(denom) <= 1e-12)
  ratio <- numer / denom
  ratio[invalid_denom] <- NA_real_
  
  # CRITICAL: Clamp ratio to handle cases where DVI slightly dips below dvi.soil due to noise
  # This prevents log(negative) generation
  ratio <- pmax(ratio, 1e-6)

  res <- K * log( ratio )
  res[!is.finite(res)] <- NA_real_
  # Allow negative PPI values (detrending may produce values below zero); do NOT clamp to 0
  res
}

add_ppi_columns <- function(df, dvi_soil = NULL) {
  df <- as.data.frame(df)

  # Ensure necessary columns exist
  if (!"DVI" %in% names(df) && all(c("nir", "red") %in% names(df))) df$DVI <- df$nir - df$red
  if (!inherits(df[["date"]], "Date")) df[["date"]] <- as.Date(df[["date"]])
  if (!"doy" %in% names(df)) df[["doy"]] <- lubridate::yday(df[["date"]])

  valid_dvi <- is.finite(df$DVI)

  # --- CRITICAL: Require explicit dvi_soil baseline ---
  # dvi_soil must be provided (scalar or a vector matching rows). No automatic computation from 'Veg' or DJF is allowed.
  if (missing(dvi_soil) || is.null(dvi_soil) || all(!is.finite(dvi_soil))) {
    stop("[PPI ERROR] dvi_soil parameter must be provided (scalar or vector) and contain finite values for rows with finite DVI.")
  }

  # Assign per-row dvi_soil
  df$dvi_soil <- rep(NA_real_, nrow(df))
  if (length(dvi_soil) == 1) {
    if (!is.finite(dvi_soil)) stop("[PPI ERROR] provided scalar dvi_soil is not finite.")
    df$dvi_soil[] <- as.numeric(dvi_soil)
  } else if (length(dvi_soil) == nrow(df)) {
    df$dvi_soil <- as.numeric(dvi_soil)
  } else {
    stop("[PPI ERROR] dvi_soil must be a scalar or a vector with length equal to number of rows in df.")
  }

  # Require that every row with a finite DVI has a finite dvi_soil baseline; fail fast if not
  need_idx <- valid_dvi & !is.finite(df$dvi_soil)
  if (any(need_idx)) {
    stop(sprintf("[PPI ERROR] dvi_soil baseline not established for %d rows; provide explicit dvi_soil for these rows.\n", sum(need_idx)))
  }

  # NOTE: Per the current experimental setting, the canopy-maximum M is fixed
  # to PPI_FULL_VEG_COVER (default 0.7). We DO NOT compute per-location max(DVI)
  # for M — use the fixed value and nudge above dvi_soil when required.


  # Latitude Handling for SZA
  lat_fallback <- 40.2
  lat_col <- intersect(c("lat", "latitude", "target_lat"), names(df))[1]
  if (!is.na(lat_col)) {
    lat_vals <- as.numeric(as.character(df[[lat_col]]))
    lat_vals[!is.finite(lat_vals)] <- lat_fallback
    df$lat_use <- lat_vals
  } else {
    df$lat_use <- lat_fallback
  }

  # Calculate SZA
  df$zenith.angle <- calculate_solar_zenith(df$lat_use, df$doy)
  
  # Calculate PPI
  df$PPI <- NA_real_
  calc_idx <- complete.cases(df$DVI, df$zenith.angle, df$dvi_soil)

  # Debug: Check if we have any valid data
  if (!any(calc_idx)) {
    n_missing_dvi <- sum(is.na(df$DVI) | !is.finite(df$DVI))
    n_missing_zenith <- sum(is.na(df$zenith.angle) | !is.finite(df$zenith.angle))
    n_missing_dvi_soil <- sum(is.na(df$dvi_soil) | !is.finite(df$dvi_soil))
    warning(sprintf("No complete cases for PPI calculation (missing: %d DVI, %d zenith.angle, %d dvi_soil out of %d rows)",
                n_missing_dvi, n_missing_zenith, n_missing_dvi_soil, nrow(df)))
  }

  # Only run if we have data
  if (any(calc_idx)) {
    # Prefer per-location M: if location_id present, compute M as max DVI per location and call ppi per-location
    if ("location_id" %in% names(df)) {
      locs <- unique(df$location_id[calc_idx])
      for (loc in locs) {
        idx_loc <- calc_idx & df$location_id == loc
        if (!any(idx_loc)) next
        dvi_loc <- df$DVI[idx_loc]
        zen_loc <- df$zenith.angle[idx_loc]
        dsoil_loc <- df$dvi_soil[idx_loc]
        # pick a single dvi_soil per location (use first finite)
        dsoil_val <- dsoil_loc[is.finite(dsoil_loc)][1]
        if (!is.finite(dsoil_val)) next
        if (!any(is.finite(dvi_loc) & is.finite(zen_loc))) next
        # Use fixed canopy-maximum M (do NOT compute max(DVI) dynamically)
        M_loc <- if (exists("PPI_FULL_VEG_COVER")) get("PPI_FULL_VEG_COVER") else 0.7
        # Ensure M is strictly above the local soil baseline to avoid invalid denominator
        if (!is.finite(M_loc) || M_loc <= dsoil_val) {
          M_loc <- dsoil_val + 1e-3
          warning(sprintf("[PPI] Fixed M <= dvi_soil for location '%s' — using M = dvi_soil + 1e-3 (%.6f)", loc, M_loc))
        }
        df$PPI[idx_loc] <- ppi(dvi_loc, zen_loc, M = M_loc, dvi.soil = dsoil_val)
      }
    } else {
      # No location info: compute a single M from the available rows (no fallback constants)
      # Use fixed global M (do NOT compute max(DVI) dynamically)
      M_global <- if (exists("PPI_FULL_VEG_COVER")) get("PPI_FULL_VEG_COVER") else 0.7
      min_dsoil <- suppressWarnings(min(df$dvi_soil[calc_idx], na.rm = TRUE))
      if (is.finite(min_dsoil) && M_global <= min_dsoil) {
        # Nudge M above the soil baseline to avoid errors
        M_global <- min_dsoil + 1e-3
        warning(sprintf("[PPI] Fixed global M <= min dvi_soil; using M = dvi_soil + 1e-3 (%.6f)", M_global))
      }
      df$PPI[calc_idx] <- ppi(df$DVI[calc_idx], df$zenith.angle[calc_idx], M = M_global, dvi.soil = df$dvi_soil[calc_idx])
    }

  }

  df$lat_use <- NULL
  return(df)
}

auto_add_ppi_columns <- function(df, dvi_soil = NULL, env_var = "MESMA_DVI_SOIL") {
  # Strict behavior: require caller to pass explicit dvi_soil. Do not attempt environment or automatic fallbacks.
  if (!is.null(dvi_soil) && (length(dvi_soil) == 1 || length(dvi_soil) == nrow(df))) {
    df_out <- add_ppi_columns(df, dvi_soil = dvi_soil)
    return(list(df = df_out, added = TRUE, reason = "passed_baseline"))
  }
  return(list(df = df, added = FALSE, reason = "dvi_soil_required"))
}

compute_indices_from_bands <- function(df,
                                      raw_bands = NULL) {
  if (is.null(df)) stop("compute_indices_from_bands: df is NULL")
  if (nrow(df) == 0) return(df)
  eps <- 1e-9

  # Determine raw bands to consider.  We purposely avoid referencing RAW_BANDS
  # here to keep the function self-contained for parallel workers.
  if (is.null(raw_bands)) {
    raw_bands <- c("blue", "green", "red", "nir", "swir1", "swir2")
  }

  has_bands <- intersect(raw_bands, names(df))
  if (length(has_bands) == 0) return(df)

  # Convert bands to numeric once (avoids 50+ repeated as.numeric() calls)
  b <- list()
  for (bn in has_bands) b[[bn]] <- as.numeric(df[[bn]])

  has <- function(...) all(c(...) %in% names(b))

  # Compute all required indices
  if (has('nir','red')) df$DVI <- b$nir - b$red
  if (has('red','nir')) df$NDDI <- (b$red - b$nir) / (b$red + b$nir + eps)
  if (has('nir','red','blue')) df$EVI <- 2.5 * (b$nir - b$red) / (b$nir + 6 * b$red - 7.5 * b$blue + 1)
  if (has('red','blue','nir')) df$PSRI <- (b$red - b$blue) / (b$nir + eps)
  if (has('nir','swir1')) df$NDMI <- (b$nir - b$swir1) / (b$nir + b$swir1 + eps)
  if (has('swir1','swir2')) df$NDTI <- (b$swir1 - b$swir2) / (b$swir1 + b$swir2 + eps)
  if (has('swir1','nir')) df$MSI <- b$swir1 / (b$nir + eps)
  if (has('nir','red')) df$MSAVI <- (2 * b$nir + 1 - sqrt(pmax(0, (2 * b$nir + 1)^2 - 8 * (b$nir - b$red)))) / 2
  if (has('nir','red')) df$NDVI <- (b$nir - b$red) / (b$nir + b$red + eps)
  if (has('nir','red')) df$OSAVI <- (b$nir - b$red) / (b$nir + b$red + 0.16)
  if (has('nir','red')) df$NIRv <- b$nir * ((b$nir - b$red) / (b$nir + b$red + eps)) * 1.3
  if (has('nir','swir2')) df$NBR <- (b$nir - b$swir2) / (b$nir + b$swir2 + eps)
  if (has('swir1','swir2')) df$TCW <- (b$swir1 - b$swir2) / (b$swir1 + b$swir2 + eps)
  if (has('green','red')) df$PRI <- (b$green - b$red) / (b$green + b$red + eps)
  if (has('red','green','blue')) df$MCARI <- ((b$red - b$green) - 0.2*(b$red - b$blue)) * (b$red / (b$green + eps))
  if (has('blue','green','red','nir','swir1','swir2')) {
    df$TCB <- 0.3029 * b$blue + 0.2786 * b$green + 0.4733 * b$red + 0.5599 * b$nir + 0.508 * b$swir1 + 0.1872 * b$swir2
    df$GVI <- -0.2941 * b$blue - 0.243 * b$green - 0.5424 * b$red + 0.7276 * b$nir + 0.0713 * b$swir1 - 0.1608 * b$swir2
  }

  df
}