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
  
  # Identify based on "no soil" fraction if available
  ns_cols <- c("no soil", "no.soil", "no_soil")
  ns_col <- intersect(ns_cols, names(df))[1]
  if (!is.na(ns_col)) {
    vals <- as.numeric(as.character(df[[ns_col]]))
    barren_idx <- barren_idx | (is.finite(vals) & vals > 0.5)
  }
  
  valid_dvi <- is.finite(df$DVI)
  barren_idx <- barren_idx & valid_dvi

  # --- CRITICAL: Determine dvi_soil baseline ---
  if ("location_id" %in% names(df)) {
    # Compute 5th percentile DVI per location as baseline
    df <- df %>%
      dplyr::group_by(location_id) %>%
      dplyr::mutate(dvi_soil = quantile(DVI[is.finite(DVI)], 0.05, na.rm = TRUE)) %>%
      dplyr::ungroup()
    cat("[PPI] Using per-location 5th percentile DVI as baseline\n")
  } else {
    # Fallback: single baseline
    dvi_soil_calc <- NULL
    
    # Priority 1: Use provided dvi_soil parameter (e.g., from training data)
    if (!is.null(dvi_soil) && is.finite(dvi_soil)) {
      dvi_soil_calc <- dvi_soil
      cat(sprintf("[PPI] Using provided dvi_soil baseline: %.4f\n", dvi_soil_calc))
    }
    # Priority 2: Compute from barren pixels in current dataset (use mean DVI of barren rows)
    else if (any(barren_idx)) {
      dvi_soil_calc <- as.numeric(mean(df$DVI[barren_idx], na.rm = TRUE))
      cat(sprintf("[PPI] Computed baseline from %d barren observations (mean DVI of barren): dvi_soil = %.4f\n",
                  sum(barren_idx), dvi_soil_calc))
    }
    # Priority 3: Error - no baseline available
    else {
      stop("[PPI ERROR] Cannot compute PPI: No dvi_soil provided and no barren pixels found in data. ",
           "Pass dvi_soil parameter from training data or ensure barren pixels are present.")
    }
    
    df$dvi_soil <- dvi_soil_calc
  }

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
    # Fallback logic
    user_dvi <- suppressWarnings(as.numeric(Sys.getenv(env_var)))
    if (!is.na(user_dvi)) {
       df_out <- add_ppi_columns(df, dvi_soil = user_dvi)
       return(list(df = df_out, added = TRUE, reason = "env_override"))
    }
    return(list(df = df, added = FALSE, reason = e$message))
  })
}
