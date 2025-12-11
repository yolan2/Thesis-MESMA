## Centralized PPI helpers for consistent calculation across scripts
## Provides: PPI_DVI_SOIL, calculate_solar_zenith, ppi, add_ppi_columns

calculate_solar_zenith <- function(lat, doy, hour = 10.5) {
  lat_rad <- lat * pi / 180
  dec_rad <- 23.45 * sin(2 * pi * (284 + doy) / 365) * pi / 180
  h_rad <- (hour - 12) * 15 * pi / 180
  cos_z <- sin(lat_rad) * sin(dec_rad) + cos(lat_rad) * cos(dec_rad) * cos(h_rad)
  acos(pmin(pmax(cos_z, -1), 1))
}

ppi <- function(dvi, zenith.angle, M = max(dvi) + 0.005, G = 0.5){
  stopifnot(!anyNA(dvi), !anyNA(zenith.angle))
  if(any(zenith.angle > pi)) warning("zenith.angle must be in radians, but is most probably in degrees")
  d_c <- 0.0336 + 0.0477/cos(zenith.angle)
  Q_E <- d_c + (1 - d_c) * G / cos(zenith.angle)
  K <- 1/(4*Q_E) * (1 + M)/(1 - M)

  denom <- M
  numer <- M - dvi
  ratio <- numer / denom
  invalid_denom <- (denom <= 1e-6)
  if (any(invalid_denom)) ratio[invalid_denom] <- 1
  invalid_ratio <- (ratio <= 0)
  if (any(invalid_ratio)) ratio[invalid_ratio] <- 1

  res <- - K * log( ratio )
  res[!is.finite(res)] <- NA_real_
  res <- pmin(pmax(res, 0), 3)
  res
}

add_ppi_columns <- function(df, lat) {
  df <- as.data.frame(df)
  df[["DVI"]] <- df[["nir"]] - df[["red"]]
  df[["date"]] <- as.Date(df[["date"]])
  if (!"doy" %in% names(df)) df[["doy"]] <- lubridate::yday(df[["date"]])
  if (!"year" %in% names(df)) df[["year"]] <- lubridate::year(df[["date"]])
  if (!"pheno_year" %in% names(df)) df[["pheno_year"]] <- ifelse(lubridate::month(df[["date"]]) >= 3, lubridate::year(df[["date"]]), lubridate::year(df[["date"]]) - 1)

  peak_df <- df %>%
    group_by(location_id, pheno_year) %>%
    summarise(
      DVI_max = max(DVI, na.rm = TRUE),
      doy_peak = { if (all(is.na(DVI))) NA_integer_ else { idx <- which.max(DVI); doy[idx[1]] } },
      .groups = "drop"
    )

  peak_df$DVI_max[!is.finite(peak_df$DVI_max)] <- NA_real_
  peak_df$doy_peak[!is.finite(peak_df$doy_peak)] <- NA_integer_

  if ("lat" %in% names(df)) {
    lat_lookup <- df %>% group_by(location_id) %>% summarise(lat_use = mean(lat, na.rm = TRUE), .groups = "drop")
  } else {
    lat_lookup <- data.frame(location_id = unique(df$location_id), lat_use = lat)
  }

  peak_df <- peak_df %>% left_join(lat_lookup, by = "location_id")
  peak_df$lat_use[!is.finite(peak_df$lat_use)] <- lat

  peak_df$zenith.angle <- NA_real_
  zen_idx <- complete.cases(peak_df$lat_use, peak_df$doy_peak)
  if (any(zen_idx)) {
    peak_df$zenith.angle[zen_idx] <- calculate_solar_zenith(lat = peak_df$lat_use[zen_idx], doy = peak_df$doy_peak[zen_idx])
  }

  peak_df$PPI <- NA_real_
  ppi_idx <- complete.cases(peak_df$DVI_max, peak_df$zenith.angle)
  if (any(ppi_idx)) {
    peak_df$PPI[ppi_idx] <- ppi(dvi = peak_df$DVI_max[ppi_idx], zenith.angle = peak_df$zenith.angle[ppi_idx], M = peak_df$DVI_max[ppi_idx] + 0.005)
  }

  drop_cols <- intersect(c("DVI_max", "PPI", "zenith.angle"), names(df))
  if (length(drop_cols)) df[drop_cols] <- NULL

  df <- df %>% left_join(peak_df %>% select(location_id, pheno_year, DVI_max, zenith.angle, PPI), by = c("location_id", "pheno_year"))
  df <- as.data.frame(df)
  return(df)
}

## End PPI helpers
