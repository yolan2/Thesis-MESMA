# Script to plot average PPI over time for each vegetation type

# Load required libraries
library(dplyr)
library(lubridate)
library(ggplot2)
library(readr)
library(sf)

# --- Helper Functions (Copied from january_averages.R) ---

make_location_id <- function(lon, lat) {
  # Ensure inputs are numeric
  lon <- as.numeric(lon)
  lat <- as.numeric(lat)
  
  if (length(lon) == 1 && length(lat) == 1) {
    if (!is.finite(lon) || !is.finite(lat)) return(NA_character_)
    sprintf("L_%0.6f_%0.6f", round(lat, 6), round(lon, 6))
  } else {
    # Vectorized approach
    res <- rep(NA_character_, length(lon))
    valid <- is.finite(lon) & is.finite(lat)
    if (any(valid)) {
      res[valid] <- sprintf("L_%0.6f_%0.6f", round(lat[valid], 6), round(lon[valid], 6))
    }
    res
  }
}

safe_as_numeric <- function(x) {
  if (is.null(x)) return(x)
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    s <- trimws(x)
    lower <- tolower(s)
    lower[lower %in% c("true", "t")] <- "1"
    lower[lower %in% c("false", "f")] <- "0"
    suppressWarnings(num <- as.numeric(lower))
    return(num)
  }
  if (is.numeric(x)) return(as.numeric(x))
  suppressWarnings(num <- as.numeric(as.character(x)))
  num
}

# Legacy normalize_no_soil_col removed — 'no soil' / 'no_soil' columns are not used in PPI plotting any more.

## Centralized PPI helpers
if (file.exists("ppi_helpers.R")) {
  source("ppi_helpers.R")
} else {
  warning("ppi_helpers.R not found; using local inline PPI logic (ensure ppi_helpers.R in project root for consistent PPI calculation)")
}
# Optional visualization helpers (provide shading for excluded years)
if (file.exists("mesma_helpers.R")) {
  source("mesma_helpers.R")
} else {
  # mesma_helpers.R not found; shading helper will be unavailable
} 

# --- Main Script ---

# Define input file paths
INPUT_CSV <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/phenology_results/hls_phenology_data.csv"
GEOJSON_PATH <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/updated_zuizer_zonder_foto_UTM.geojson"
OUTPUT_DIR <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/phenology_results"

# Load the phenology data
cat("Loading data from:", INPUT_CSV, "\n")
df <- readr::read_csv(INPUT_CSV, show_col_types = FALSE)
# Exclude years 1992-1999 from plotting and analysis (if present)
if ("date" %in% names(df)) {
  if (!lubridate::is.Date(df$date)) df$date <- as.Date(df$date)
  n_before <- nrow(df)
  years_to_drop <- 1992:1999
  df <- df[!(lubridate::year(df$date) %in% years_to_drop), , drop = FALSE]
  if (n_before != nrow(df)) cat(sprintf("[DATA FILTER] Dropped %d rows from years %d-%d from input in plotting script\n", n_before - nrow(df), min(years_to_drop), max(years_to_drop)))
} else if ("year" %in% names(df)) {
  n_before <- nrow(df)
  years_to_drop <- 1992:1999
  df <- df[!(as.integer(df$year) %in% years_to_drop), , drop = FALSE]
  if (n_before != nrow(df)) cat(sprintf("[DATA FILTER] Dropped %d rows from years %d-%d from input in plotting script (using 'year' column)\n", n_before - nrow(df), min(years_to_drop), max(years_to_drop)))
} else {
  # No date/year column available; nothing to drop
}

# --- Location ID and Join Logic (Copied from january_averages.R) ---
# Check if location_id already exists and is useful
use_existing_id <- FALSE
if ("location_id" %in% names(df)) {
  n_unique_id <- length(unique(df$location_id))
  if (n_unique_id > 1) {
    cat("Using existing 'location_id' from CSV.\n")
    use_existing_id <- TRUE
  }
}

if (!use_existing_id) {
  # Check if target coordinates are available and have more variance
  use_target_coords <- FALSE
  if (all(c("target_lat", "target_lon") %in% names(df))) {
    n_unique_tlat <- length(unique(df$target_lat))
    if (exists("n_unique_lat") && n_unique_lat == 1 && n_unique_tlat > 1) {
      use_target_coords <- TRUE
    }
  }

  # Set location_id from lat lon (or target_lat/target_lon)
  if (use_target_coords) {
    df$location_id <- make_location_id(df$target_lon, df$target_lat)
    df$lat <- df$target_lat
    df$lon <- df$target_lon
  } else {
    df$location_id <- make_location_id(df$lon, df$lat) 
  }
}

# Load GeoJSON for vegetation types
cat("Loading GeoJSON from:", GEOJSON_PATH, "\n")
gpts_raw <- tryCatch({
  sf::st_read(GEOJSON_PATH, quiet = TRUE)
}, error = function(e) {
  cat("Warning: Failed to read GeoJSON (", e$message, "). Skipping Veg mapping.\n")
  NULL
})

# Calculate mean latitude of the study area for zenith angle fallback
if (!is.null(gpts_raw)) {
  tryCatch({
    gpts_wgs84 <- sf::st_transform(gpts_raw, 4326)
    coords_wgs84 <- sf::st_coordinates(gpts_wgs84)
  }, error = function(e) {
    cat("Warning: Could not calculate mean latitude from GeoJSON:", e$message, "\n")
  })
}

# Load vegetation mapping from CSV instead of GeoJSON
MAPPING_CSV <- "C:/Users/yolan/Downloads/landsat_timeseries_vegetation_filtered (4).csv"
if (file.exists(MAPPING_CSV)) {
    cat("Loading vegetation mapping from:", MAPPING_CSV, "\n")
    map_df <- readr::read_csv(MAPPING_CSV, show_col_types = FALSE)
    veg_cols <- names(map_df)[tolower(names(map_df)) %in% c("vegetation", "veg", "class")]
    if (length(veg_cols) > 0) {
      map_df$Veg <- as.character(map_df[[veg_cols[1]]])
      map_df$Veg <- tolower(trimws(as.character(map_df$Veg)))
      map_df$Veg <- ifelse(grepl("phragmites", map_df$Veg, ignore.case = TRUE) |
                           map_df$Veg %in% c("herbs", "alhagi", "salicornia", "halocnemum"),
                           "herbs", map_df$Veg)
    }

    if (!"location_id" %in% names(map_df) && all(c("lon", "lat") %in% names(map_df))) {
      map_df$location_id <- make_location_id(map_df$lon, map_df$lat)
    }

    if ("location_id" %in% names(map_df)) {
      if (!"Veg" %in% names(map_df)) map_df$Veg <- NA_character_
      cat("Mapping CSV will NOT be joined into the main observations data; it will only be used for Veg mapping and summaries.\n")
      cat("DEBUG: Unique Veg types in mapping CSV:\n")
      print(table(map_df$Veg, useNA = "ifany"))
    } else {
      cat("Warning: mapping CSV lacks 'location_id' and 'lon'/'lat'; cannot use location-level mapping.\n")
    }
} else {
    cat("Mapping CSV not found at:", MAPPING_CSV, "- proceeding without Veg/no-soil mapping.\n")
}

# --- PPI Calculation ---
cat("Calculating PPI...\n")

# Ensure date column is Date type
df$date <- as.Date(df$date)
df$doy <- lubridate::yday(df$date)

# Calculate DVI if missing
if (!"DVI" %in% names(df) && all(c("nir", "red") %in% names(df))) {
  df$DVI <- df$nir - df$red
  cat("Calculated DVI from nir and red bands\n")
}

# Calculate DVI_max per location
  if ("DVI" %in% names(df)) {
  df <- df |> dplyr::group_by(location_id) |> dplyr::mutate(DVI_max = max(DVI, na.rm = TRUE)) |> dplyr::ungroup()
  df$DVI_max[!is.finite(df$DVI_max)] <- NA_real_
}

# Zenith Angle
zenith_candidates <- c("zenith.angle", "zenith", "sun_zenith", "sun_zenith_angle", "sun_zenith_deg", "sun_zenith_angle_deg")
elev_candidates <- c("sun_elevation", "sun_elev", "sun_elevation_deg", "sun_elevation_angle")

if (any(zenith_candidates %in% names(df))) {
  zen_col <- zenith_candidates[zenith_candidates %in% names(df)][1]
  df$zenith.angle <- df[[zen_col]]
} else if (any(elev_candidates %in% names(df))) {
  elev_col <- elev_candidates[elev_candidates %in% names(df)][1]
  df$zenith.angle <- (90 - df[[elev_col]]) * pi / 180
} else {
    # Use per-row latitude when available; otherwise leave zenith angle NA
    if ("lat" %in% names(df) && any(is.finite(df$lat))) {
      use_lat <- as.numeric(df$lat)
      df$zenith.angle <- calculate_solar_zenith(use_lat, df$doy)
    } else {
      cat("No per-row latitude available; leaving 'zenith.angle' as NA. PPI will not be computed without SZA.\n")
      df$zenith.angle <- NA_real_
    }
}

# Ensure zenith is in radians
if ("zenith.angle" %in% names(df) && any(df$zenith.angle > pi, na.rm = TRUE)) {
  df$zenith.angle <- df$zenith.angle * pi / 180
}

# Calculate PPI
df$PPI <- NA_real_
if ("DVI" %in% names(df) && "zenith.angle" %in% names(df) && "DVI_max" %in% names(df)) {
    complete_idx <- complete.cases(df$DVI, df$zenith.angle, df$DVI_max)
    if (any(complete_idx)) {
        dsoil <- if (exists("PPI_DVI_SOIL", envir = .GlobalEnv)) get("PPI_DVI_SOIL", envir = .GlobalEnv) else 0.09
        # Use fixed M=0.7 per standardized PPI implementation
        df$PPI[complete_idx] <- ppi(df$DVI[complete_idx], df$zenith.angle[complete_idx], dvi.soil = dsoil)
    }
}

# --- Aggregation and Plotting ---

cat("Aggregating data...\n")

# Filter out rows with missing PPI or Veg
# And filter for specific vegetation types: barren, herbs, populus, tamarix
target_veg <- c("barren", "herbs", "populus", "tamarix", "woody_unknown")

plot_data <- df |> 
  dplyr::filter(!is.na(PPI), !is.na(Veg)) |> 
  dplyr::filter(Veg %in% target_veg) |> 
  dplyr::mutate(year = year(date),
         doy = yday(date))

# Calculate average PPI per Veg type and Date
ppi_summary <- plot_data |> 
  dplyr::group_by(Veg, date) |> 
  dplyr::summarize(
    mean_PPI = mean(PPI, na.rm = TRUE),
    sd_PPI = sd(PPI, na.rm = TRUE),
    n = dplyr::n(),
    se_PPI = sd_PPI / sqrt(n),
    .groups = "drop"
  )

cat("Generating plot...\n")

p <- ggplot(ppi_summary, aes(x = date, y = mean_PPI, color = Veg, fill = Veg)) +
  add_excluded_years_shade(is_date = TRUE) + add_year_lines(is_date = TRUE) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean_PPI - se_PPI, ymax = mean_PPI + se_PPI), alpha = 0.2, color = NA) +
  labs(title = "Average Plant Phenology Index (PPI) over Time by Vegetation Type",
       subtitle = "Shaded area represents Standard Error",
       x = "Date",
       y = "PPI") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Save plot
output_file <- file.path(OUTPUT_DIR, "PPI_timeseries_by_veg.png")
ggsave(output_file, plot = p, width = 12, height = 8, bg = "white")

cat("Plot saved to:", output_file, "\n")

# Also plot by DOY (Day of Year) to see seasonality
p_doy <- ggplot(plot_data, aes(x = doy, y = PPI, color = Veg)) +
  geom_smooth(se = FALSE) + 
  labs(title = "Seasonal PPI Profile by Vegetation Type (Smoothed)",
       x = "Day of Year",
       y = "PPI") +
  theme_minimal() +
  facet_wrap(~year)

output_file_doy <- file.path(OUTPUT_DIR, "PPI_seasonality_by_veg.png")
ggsave(output_file_doy, plot = p_doy, width = 12, height = 8, bg = "white")
cat("Seasonality plot saved to:", output_file_doy, "\n")
