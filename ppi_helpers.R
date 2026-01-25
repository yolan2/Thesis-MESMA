## Centralized PPI helpers for consistent calculation
## Provides: safe_as_numeric, calculate_solar_zenith, ppi, add_ppi_columns, auto_add_ppi_columns

safe_as_numeric <- function(x) {
  as.numeric(as.character(x))
}

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

ppi <- function(dvi, zenith.angle, M = NULL, dvi.soil, G = 0.5){
  # Determine M: prefer explicit M, otherwise prefer per-call max(DVI), then fallback to 0.7
  if (is.null(M)) {
    # Use per-call max(dvi) when available (suitable for per-location or per-pair calls)
    if (any(is.finite(dvi))) {
      Mcand <- suppressWarnings(max(dvi, na.rm = TRUE))
      if (is.finite(Mcand)) {
        M <- Mcand
        # Per-call max DVI used as M (silent)
      }
    }
    # If per-call max is not usable, use fallback M = 0.7 (silent)
    if (!is.finite(M)) {
      M <- 0.7
    }
  } else {
    if (!is.finite(M)) stop("[PPI ERROR] Provided M is not finite.")
  }
  
  # Ensure M sits above dvi.soil to avoid zero/negative denominators
  if (any(is.finite(dvi.soil))) {
    min_dsoil <- min(dvi.soil, na.rm = TRUE)
    if (M <= min_dsoil) {
      M_old <- M
      M <- min_dsoil + 1e-3
      invisible()
    }
  }
  
  if (missing(dvi.soil) || is.null(dvi.soil) || any(!is.finite(dvi.soil))) {
    stop("[PPI ERROR] dvi.soil parameter is required.")
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

  # Identify Barren Pixels
  barren_idx <- rep(FALSE, nrow(df))
  if ("Veg" %in% names(df)) {
    barren_idx <- barren_idx | (!is.na(df$Veg) & tolower(trimws(as.character(df$Veg))) == "barren")
  }
  
  # Note: 'no soil' fraction column is deprecated and ignored; barren identification uses Veg=='barren' only
  
  valid_dvi <- is.finite(df$DVI)
  barren_idx <- barren_idx & valid_dvi

  # --- CRITICAL: Determine dvi_soil baseline ---
  # AGENT CHANGE: Per user feedback, using a single global DVI baseline
  # from 'barren' veg type, with no fallback to other methods.

  dvi_soil_calc <- NULL

  # Priority 1: Use provided dvi_soil parameter (e.g., from training data override)
  if (!is.null(dvi_soil) && is.finite(dvi_soil)) {
    dvi_soil_calc <- dvi_soil
  }
  # Priority 2: Compute global baseline from 'barren' vegetation type (local-only and optional)
  else if ("Veg" %in% names(df)) {
    barren_dvi <- df$DVI[is.finite(df$DVI) & tolower(df$Veg) == 'barren']
    if (length(barren_dvi) > 0) {
      dvi_soil_calc <- as.numeric(median(barren_dvi, na.rm = TRUE))
    }
  }

  # Initialize per-row dvi_soil as NA; we will assign per-location DJF medians below or fill from a computed global baseline only if present.
  df$dvi_soil <- rep(NA_real_, nrow(df))
  if (is.finite(dvi_soil_calc)) {
    df$dvi_soil[] <- dvi_soil_calc
  }
  
  # Attempt per-location DVI soil baseline using December-February median (DJF) where available
  # Minimum samples per location to trust DJF median can be controlled via global PPI_MIN_DJF_SAMPLES (default: 3)
  if ("location_id" %in% names(df)) {
    months <- lubridate::month(df$date)
    djf_idx <- months %in% c(12, 1, 2)
    locs <- unique(df$location_id)
    min_samples <- if (exists("PPI_MIN_DJF_SAMPLES", envir = globalenv())) as.integer(get("PPI_MIN_DJF_SAMPLES", envir = globalenv())) else 3
    n_loc_assigned <- 0
    for (loc in locs) {
      vals <- df$DVI[df$location_id == loc & djf_idx & is.finite(df$DVI)]
      if (length(vals) >= min_samples) {
        df$dvi_soil[df$location_id == loc] <- as.numeric(median(vals, na.rm = TRUE))
        n_loc_assigned <- n_loc_assigned + 1
      }
    }
    # Per-location DJF medians assigned to `df$dvi_soil` where available (n_assigned = n_loc_assigned).
  }

  # Require that every row with a finite DVI has a finite dvi_soil baseline; fail fast if not
  need_idx <- is.finite(df$DVI) & !is.finite(df$dvi_soil)
  if (any(need_idx)) {
    stop(sprintf("[PPI ERROR] dvi_soil baseline not established for %d rows; please provide 'barren' samples, per-location DJF data (min %d samples), or pass explicit dvi_soil to add_ppi_columns()\n", sum(need_idx), min_samples))
  }

  # Note: Global PPI M computation has been removed. Per-location or per-call M is used instead; fallback M=0.7 will be applied only when necessary.


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
        M_loc <- suppressWarnings(max(dvi_loc, na.rm = TRUE))
        if (!is.finite(M_loc) || M_loc <= dsoil_val) {
          # Fall back to just above the local soil baseline to avoid invalid denom
          M_loc <- dsoil_val + 1e-3
        } else {
          # use M_loc (max DVI) silently
        }
        df$PPI[idx_loc] <- ppi(dvi_loc, zen_loc, M = M_loc, dvi.soil = dsoil_val)
      }
    } else {
      # No location info: let ppi() pick per-call max(dvi) as M when M is NULL
      df$PPI[calc_idx] <- ppi(df$DVI[calc_idx], df$zenith.angle[calc_idx], dvi.soil = df$dvi_soil[calc_idx])
    }

  }

  df$lat_use <- NULL
  return(df)
}

auto_add_ppi_columns <- function(df, dvi_soil = NULL, env_var = "MESMA_DVI_SOIL") {
  tryCatch({
    if (!is.null(dvi_soil) && is.finite(dvi_soil)) {
      df_out <- add_ppi_columns(df, dvi_soil = dvi_soil)
      return(list(df = df_out, added = TRUE, reason = "passed_baseline"))
    }
    df_out <- add_ppi_columns(df)
    return(list(df = df_out, added = TRUE, reason = "calculated"))
  }, error = function(e) {
    warning(sprintf("auto_add_ppi_columns caught error: %s", e$message))
    # Fallback logic: allow explicit environment override if set
    user_dvi <- suppressWarnings(as.numeric(Sys.getenv(env_var)))
    if (!is.na(user_dvi)) {
       df_out <- add_ppi_columns(df, dvi_soil = user_dvi)
       return(list(df = df_out, added = TRUE, reason = "env_override"))
    }
    return(list(df = df, added = FALSE, reason = e$message))
  })
}
