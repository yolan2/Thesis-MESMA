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
    sprintf("L_%0.4f_%0.4f", round(lon, 4), round(lat, 4))
  } else {
    # Vectorized approach
    res <- rep(NA_character_, length(lon))
    valid <- is.finite(lon) & is.finite(lat)
    if (any(valid)) {
      res[valid] <- sprintf("L_%0.4f_%0.4f", round(lon[valid], 4), round(lat[valid], 4))
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

normalize_no_soil_col <- function(tbl) {
  if (is.null(tbl) || !is.data.frame(tbl)) return(tbl)
  nm <- names(tbl)
  candidates <- c("no soil", "no_soil", "no.soil", "__no soil__", "__no_soil__", ".__no soil__", ".__no_soil__")
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

## Centralized PPI helpers
if (file.exists("ppi_helpers.R")) {
  source("ppi_helpers.R")
} else {
  warning("ppi_helpers.R not found; using local inline PPI logic (ensure ppi_helpers.R in project root for consistent PPI calculation)")
}

# --- Main Script ---

# Define input file paths
INPUT_CSV <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/phenology_results/hls_phenology_data.csv"
GEOJSON_PATH <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/updated_zuizer_zonder_foto_UTM.geojson"
OUTPUT_DIR <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/phenology_results"

# Load the phenology data
cat("Loading data from:", INPUT_CSV, "\n")
df <- readr::read_csv(INPUT_CSV, show_col_types = FALSE)

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
    mean_lat_study_area <- mean(coords_wgs84[, 2], na.rm = TRUE)
  }, error = function(e) {
    cat("Warning: Could not calculate mean latitude from GeoJSON:", e$message, "\n")
  })
}

# Process GeoJSON to get Veg mapping
if (!is.null(gpts_raw)) {
    geojson_names <- names(gpts_raw)
    normalized_names <- gsub("[^a-z0-9]+", "_", tolower(geojson_names))
    no_soil_col <- geojson_names[normalized_names == "no_soil"]
    if (length(no_soil_col) > 0) {
        no_soil_raw <- gpts_raw[[no_soil_col[1]]]
        no_soil_vals <- safe_as_numeric(no_soil_raw)
        gpts_raw$`.__no soil__` <- no_soil_vals
    }

    matched_cols <- names(gpts_raw)[tolower(names(gpts_raw)) %in% c("vegetation", "veg", "class")]
    if (length(matched_cols) > 0) {
      veg_col_orig <- matched_cols[1]
      gpts_raw$.__veg__ <- as.character(gpts_raw[[veg_col_orig]])
    } else {
      gpts_raw$.__veg__ <- NA_character_
    }
    
    gpts_wgs84 <- sf::st_transform(gpts_raw, 4326)
    coords_wgs84 <- sf::st_coordinates(gpts_wgs84)
    gpts_raw$.__lon__ <- coords_wgs84[, 1]
    gpts_raw$.__lat__ <- coords_wgs84[, 2]
    gpts_raw$location_id_geo <- make_location_id(gpts_raw$.__lon__, gpts_raw$.__lat__)
    gpts_raw$location_id <- as.character(seq_len(nrow(gpts_raw)))
    
    gdf <- sf::st_drop_geometry(gpts_raw)
    
    # Determine mapping strategy
    use_seq_id <- FALSE
    if (exists("df") && "location_id" %in% names(df)) {
        sample_id <- as.character(df$location_id[1])
        if (!grepl("^L_", sample_id)) {
            use_seq_id <- TRUE
        }
    }
    
    if (use_seq_id) {
        gpts_map <- gdf %>%
        dplyr::select(location_id, Veg = .__veg__, `no soil` = `.__no soil__`) %>%
        dplyr::mutate(location_row = as.character(seq_len(dplyr::n()))) %>%
        dplyr::distinct(location_id, .keep_all = TRUE)
    } else {
        gpts_map <- gdf %>%
        dplyr::select(location_id = location_id_geo, Veg = .__veg__, `no soil` = `.__no soil__`) %>%
        dplyr::mutate(location_row = as.character(seq_len(dplyr::n()))) %>%
        dplyr::distinct(location_id, .keep_all = TRUE)
    }
    
    # Join
    if ("location_id" %in% names(df)) df$location_id <- as.character(df$location_id)
    if ("location_id" %in% names(gpts_map)) gpts_map$location_id <- as.character(gpts_map$location_id)
    
    df <- dplyr::left_join(df, gpts_map, by = "location_id", suffix = c("", ".geo"))
    
    if ("Veg.geo" %in% names(df)) {
      if (!"Veg" %in% names(df)) df$Veg <- NA_character_
      df$Veg <- ifelse(is.na(df$Veg) | df$Veg == "", df$Veg.geo, df$Veg)
      df$Veg.geo <- NULL
    }
    df <- normalize_no_soil_col(df)
    if ("Veg" %in% names(df)) {
      df$Veg <- tolower(df$Veg)
    }
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
  df <- df %>%
    group_by(location_id) %>%
    mutate(DVI_max = max(DVI, na.rm = TRUE)) %>%
    ungroup()
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
    # Fallback calculation
    use_lat <- NULL
    if ("lat" %in% names(df) && !all(is.na(df$lat))) {
      use_lat <- df$lat
    } else if (exists("mean_lat_study_area")) {
      use_lat <- rep(mean_lat_study_area, nrow(df))
    }
    
    if (!is.null(use_lat)) {
        df$zenith.angle <- calculate_solar_zenith(use_lat, df$doy)
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
        df$PPI[complete_idx] <- ppi(df$DVI[complete_idx], df$zenith.angle[complete_idx], M = df$DVI_max[complete_idx] + 0.005)
    }
}

# --- Aggregation and Plotting ---

cat("Aggregating data...\n")

# Filter out rows with missing PPI or Veg
# And filter for specific vegetation types: barren, phragmites, populus, tamarix
target_veg <- c("barren", "phragmites", "populus", "tamarix")

plot_data <- df %>%
  filter(!is.na(PPI), !is.na(Veg)) %>%
  filter(Veg %in% target_veg) %>%
  mutate(year = year(date),
         doy = yday(date))

# Calculate average PPI per Veg type and Date
ppi_summary <- plot_data %>%
  group_by(Veg, date) %>%
  summarize(
    mean_PPI = mean(PPI, na.rm = TRUE),
    sd_PPI = sd(PPI, na.rm = TRUE),
    n = n(),
    se_PPI = sd_PPI / sqrt(n),
    .groups = "drop"
  )

cat("Generating plot...\n")

p <- ggplot(ppi_summary, aes(x = date, y = mean_PPI, color = Veg, fill = Veg)) +
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
