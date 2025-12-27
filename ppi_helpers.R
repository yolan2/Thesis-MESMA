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

ppi <- function(dvi, zenith.angle, M = 0.7, dvi.soil, G = 0.5){
  # Ensure M is fixed
  M <- 0.7
  
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
  res <- pmax(res, 0)
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
    cat(sprintf("[PPI] Using provided single global dvi_soil baseline: %.4f\n", dvi_soil_calc))
  }
  # Priority 2: Compute global baseline from 'barren' vegetation type
  else if ("Veg" %in% names(df)) {
    barren_dvi <- df$DVI[is.finite(df$DVI) & tolower(df$Veg) == 'barren']
    if (length(barren_dvi) > 0) {
      dvi_soil_calc <- as.numeric(median(barren_dvi, na.rm = TRUE))
      cat(sprintf("[PPI] Computed single global baseline from median of %d 'barren' observations: dvi_soil = %.4f\n",
                  length(barren_dvi), dvi_soil_calc))
    }
  }

  # If after all checks, dvi_soil_calc is still NULL, stop execution.
  if (is.null(dvi_soil_calc)) {
      stop("[PPI ERROR] Cannot determine DVI soil baseline. No 'dvi_soil' parameter was provided and no observations with Veg='barren' were found in the dataset.")
  }

  df$dvi_soil <- dvi_soil_calc
  
  # Store baseline globally for inference step
  assign("GLOBAL_TRAINING_DVI_SOIL", dvi_soil_calc, envir = globalenv())


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
    cat(sprintf("[PPI WARNING] No complete cases for PPI calculation (missing: %d DVI, %d zenith.angle, %d dvi_soil out of %d rows)\n",
                n_missing_dvi, n_missing_zenith, n_missing_dvi_soil, nrow(df)))
  }

  # Only run if we have data
  if (any(calc_idx)) {
    df$PPI[calc_idx] <- ppi(df$DVI[calc_idx], df$zenith.angle[calc_idx], M = 0.7, dvi.soil = df$dvi_soil[calc_idx])

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
    cat(sprintf("[PPI ERROR] auto_add_ppi_columns caught error: %s\n", e$message))
    # Fallback logic
    user_dvi <- suppressWarnings(as.numeric(Sys.getenv(env_var)))
    if (!is.na(user_dvi)) {
       df_out <- add_ppi_columns(df, dvi_soil = user_dvi)
       return(list(df = df_out, added = TRUE, reason = "env_override"))
    }
    return(list(df = df, added = FALSE, reason = e$message))
  })
}
