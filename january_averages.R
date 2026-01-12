# Script to calculate January averages of MSAVI, NDVI, and PPI for each vegetation type
# with bootstrapping for uncertainty estimation

# Load required libraries
library(dplyr)
library(lubridate)
library(sf)
library(boot)
library(openxlsx)
library(ggplot2)
library(tidyr)

# Helper functions
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

# Define indices of interest globally
INDICES_OF_INTEREST <- c("MSAVI", "NDVI", "PPI", "OSAVI", "MCARI", "NIRv", "PSRI", "NBR", "TCW", "NDMI", "TCB", "GVI")

calculate_indices <- function(df) {
  # Ensure bands are present
  req_bands <- c("blue", "green", "red", "nir", "swir1", "swir2")
  missing_bands <- setdiff(req_bands, names(df))
  if (length(missing_bands) > 0) {
    cat("Warning: Missing bands for index calculation:", paste(missing_bands, collapse=", "), ". Skipping index calculation.\n")
    return(df)
  }
  
  eps <- 1e-9
  
  # Calculate DVI if missing
  if (!"DVI" %in% names(df)) df$DVI <- df$nir - df$red
  
  # Calculate additional indices
  df$OSAVI <- (df$nir - df$red) / (df$nir + df$red + 0.16)
  df$MCARI <- ((df$red - df$green) - 0.2*(df$red - df$blue)) * (df$red / (df$green + eps))
  df$NIRv  <- (df$nir * ((df$nir - df$red) / (df$nir + df$red + eps))) * 1.3
  df$PSRI  <- (df$red - df$blue) / (df$nir + eps)
  df$NBR   <- (df$nir - df$swir2) / (df$nir + df$swir2 + eps)
  # Project-specific TCW (Normalized Difference) definition
  df$TCW   <- (df$swir1 - df$swir2) / (df$swir1 + df$swir2 + eps) 
  df$NDMI  <- (df$nir - df$swir1) / (df$nir + df$swir1 + eps)
  
  # Tasseled Cap (using project-specific coefficients)
  df$TCB   <- 0.3029 * df$blue + 0.2786 * df$green + 0.4733 * df$red + 0.5599 * df$nir + 0.508 * df$swir1 + 0.1872 * df$swir2
  df$GVI   <- -0.2941 * df$blue - 0.243 * df$green - 0.5424 * df$red + 0.7276 * df$nir + 0.0713 * df$swir1 - 0.1608 * df$swir2
  
  # Ensure NDVI/MSAVI are present and consistent
  df$NDVI <- (df$nir - df$red) / (df$nir + df$red + eps)
  df$MSAVI <- (2 * df$nir + 1 - sqrt(pmax(0, (2 * df$nir + 1)^2 - 8 * (df$nir - df$red)))) / 2
  
  return(df)
}

## Outlier filtering intentionally removed per user request.
## The previous implementation used a MAD-based high-outlier filter
## (remove_outliers_mad) which removed large positive PPI values.
## Per the user's instruction we no longer perform automatic outlier
## removal; therefore that helper has been removed and PPI values are
## left as-is.

# normalize_no_soil_col removed — 'no soil' column handling is deprecated and no longer used.
# If migration of legacy datasets is required, normalize the column externally before running this script.

## Use centralized PPI helper to ensure identical PPI calculation
if (file.exists("ppi_helpers.R")) {
  source("ppi_helpers.R")
} else {
  warning("ppi_helpers.R not found; falling back to local inline PPI logic. Please ensure ppi_helpers.R is present in project root for consistent PPI calculations.")
}

# Define input file paths
INPUT_CSV <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/landsat_lower_inference.csv"
# GeoJSON support removed: this script uses CSV-provided Veg / "no soil" columns and lat/lon when available

# Load the phenology data
cat("Loading data from:", INPUT_CSV, "\n")
df <- readr::read_csv(INPUT_CSV, show_col_types = FALSE)

# Normalize and prefer CSV-provided Veg / no soil values when present
if ("vegetation" %in% names(df) && !"Veg" %in% names(df)) {
  df$Veg <- df$vegetation
  cat("[NOTICE] Renamed 'vegetation' column to 'Veg' (from CSV)\n")
}
# Use CSV-provided Veg values if present (no 'no soil' processing)
if ("Veg" %in% names(df)) {
  cat("[NOTICE] Found Veg values in input CSV; skipping GeoJSON join and using CSV-provided Veg values.\n")
  skip_geojson <- TRUE
} else {
  skip_geojson <- FALSE
}

# Normalize column names for bands to lowercase to match script expectations
# The script expects 'nir' and 'red' for calculations if indices are missing
band_mapping <- c("NIR" = "nir", "Red" = "red", "Blue" = "blue", "Green" = "green", "SWIR1" = "swir1", "SWIR2" = "swir2")
for (orig in names(band_mapping)) {
  if (orig %in% names(df) && !band_mapping[orig] %in% names(df)) {
    names(df)[names(df) == orig] <- band_mapping[orig]
    cat(sprintf("Renamed column '%s' to '%s'\n", orig, band_mapping[orig]))
  }
}

# Check if location_id already exists and is useful
use_existing_id <- FALSE
if ("location_id" %in% names(df)) {
  n_unique_id <- length(unique(df$location_id))
  cat("DEBUG: Unique location_id values in CSV:", n_unique_id, "\n")
  cat("DEBUG: First 5 location_id values in CSV:", paste(head(unique(df$location_id), 5), collapse=", "), "\n")
  if (n_unique_id > 1) {
    cat("Using existing 'location_id' from CSV.\n")
    use_existing_id <- TRUE
  }
}

if (!use_existing_id) {
  # DEBUG: Check coordinate variance
  if (all(c("lat", "lon") %in% names(df))) {
    n_unique_lat <- length(unique(df$lat))
    n_unique_lon <- length(unique(df$lon))
    cat("DEBUG: Unique lat values in CSV:", n_unique_lat, "\n")
    cat("DEBUG: Unique lon values in CSV:", n_unique_lon, "\n")
  }

  # Check if target coordinates are available and have more variance
  use_target_coords <- FALSE
  if (all(c("target_lat", "target_lon") %in% names(df))) {
    n_unique_tlat <- length(unique(df$target_lat))
    n_unique_tlon <- length(unique(df$target_lon))
    cat("DEBUG: Unique target_lat values in CSV:", n_unique_tlat, "\n")
    cat("DEBUG: Unique target_lon values in CSV:", n_unique_tlon, "\n")
    
    if (exists("n_unique_lat") && n_unique_lat == 1 && n_unique_tlat > 1) {
      cat("Detected constant 'lat'/'lon' but varying 'target_lat'/'target_lon'. Switching to target coordinates for location ID.\n")
      use_target_coords <- TRUE
    }
  }

  # Set location_id from lat lon (or target_lat/target_lon)
  if (use_target_coords) {
    df$location_id <- make_location_id(df$target_lon, df$target_lat)
    # Also update lat/lon columns to be target coordinates for zenith calculation later?
    # Or keep them as is? Zenith calculation uses 'lat'. If 'lat' is constant, it might be wrong for other points.
    # Let's update 'lat' and 'lon' to be safe, as they are used for zenith angle.
    df$lat <- df$target_lat
    df$lon <- df$target_lon
  } else {
    df$location_id <- make_location_id(df$lon, df$lat) # Note: make_location_id takes (lon, lat)
  }
}

# Per-user request: remove computation and use of mean latitude. SZA is
# computed per-row only when per-row 'lat' is available. If not available,
# 'zenith.angle' will be left NA so that PPI is not computed with an invalid
# assumed study-area mean.

# Load vegetation mapping from CSV instead of GeoJSON
MAPPING_CSV <- "C:/Users/yolan/Downloads/landsat_timeseries_vegetation_filtered (4).csv"
if (file.exists(MAPPING_CSV)) {
  cat("Loading vegetation mapping from:", MAPPING_CSV, "\n")
  map_df <- readr::read_csv(MAPPING_CSV, show_col_types = FALSE)

  # Normalize common column names (detect veg column robustly)
  veg_cols <- names(map_df)[tolower(names(map_df)) %in% c("vegetation", "veg", "class")]
  if (length(veg_cols) > 0) map_df$Veg <- as.character(map_df[[veg_cols[1]]])
  # 'no soil' normalization removed; mapping CSV will only be used for Veg mapping and barren DVI estimation based on Veg=='barren'

  # Construct location_id from lon/lat if necessary
  if (!"location_id" %in% names(map_df) && all(c("lon", "lat") %in% names(map_df))) {
    map_df$location_id <- make_location_id(map_df$lon, map_df$lat)
  }

  if ("location_id" %in% names(map_df)) {
    # Ensure Veg exists or is NA
    if (!"Veg" %in% names(map_df)) map_df$Veg <- NA_character_
    gpts_map <- map_df |> dplyr::select(location_id, Veg) |> dplyr::mutate(location_id = as.character(location_id)) |> dplyr::distinct(location_id, .keep_all = TRUE)
  } else {
    cat("Warning: mapping CSV lacks 'location_id' and 'lon'/'lat'; cannot join Veg mapping.\n")
    gpts_map <- NULL
  }
} else {
  cat("Mapping CSV not found at:", MAPPING_CSV, "- skipping mapping join.\n")
  gpts_map <- NULL
}

# If mapping provided, show summary and attempt a robust join
if (exists("map_df")) {
    cat("DEBUG: Unique Veg types in mapping CSV:\n")
    print(table(map_df$Veg, useNA = "ifany"))
    cat("Mapping CSV will NOT be joined to main CSV; it will only be used to estimate soil DVI from Veg=='barren'.\n")
}

# --- TRAINING DATA (used for statistics) ---
# Use the training timeseries CSV (if available) rather than the inference CSV
TRAINING_CSV <- "C:/Users/yolan/Downloads/landsat_timeseries_vegetation_filtered (4)_fixed.csv"
if (!file.exists(TRAINING_CSV)) {
  alt_candidates <- c(
    sub("_fixed\\.csv$", ".csv", TRAINING_CSV),
    sub("\\(4\\)_fixed", "(4)", TRAINING_CSV),
    sub("\\(4\\)_fixed", "", TRAINING_CSV)
  )
  alt_found <- alt_candidates[file.exists(alt_candidates)]
  if (length(alt_found) > 0) TRAINING_CSV <- alt_found[1]
}

if (file.exists(TRAINING_CSV)) {
  cat("Loading training data from:", TRAINING_CSV, "\n")
  training_df <- readr::read_csv(TRAINING_CSV, show_col_types = FALSE)

  # Normalize and compute indices as for the main df
  if ("vegetation" %in% names(training_df) && !"Veg" %in% names(training_df)) training_df$Veg <- training_df$vegetation
  # 'no soil' normalization removed — training_df no longer uses 'no soil' column
  for (orig in names(band_mapping)) if (orig %in% names(training_df) && !band_mapping[orig] %in% names(training_df)) names(training_df)[names(training_df) == orig] <- band_mapping[orig]
  if (!"date" %in% names(training_df) && "prediction_date" %in% names(training_df)) training_df$date <- as.Date(training_df$prediction_date)
  if ("date" %in% names(training_df)) training_df$date <- as.Date(training_df$date)
  # Ensure temporal columns exist for filtering/aggregation
  if (!"year" %in% names(training_df) && "date" %in% names(training_df)) training_df$year <- lubridate::year(training_df$date)
  if (!"pheno_year" %in% names(training_df) && "date" %in% names(training_df)) training_df$pheno_year <- ifelse(lubridate::month(training_df$date) >= 3, lubridate::year(training_df$date), lubridate::year(training_df$date) - 1)
  if (!"month" %in% names(training_df) && "date" %in% names(training_df)) training_df$month <- lubridate::month(training_df$date)
  
  # Calculate indices for training data
  training_df <- calculate_indices(training_df)
  
  if (!"location_id" %in% names(training_df) && all(c("lon","lat") %in% names(training_df))) training_df$location_id <- make_location_id(training_df$lon, training_df$lat)

  # Ensure PPI is present for training data (uses mapping or internal barren rows)
  if (exists("auto_add_ppi_columns")) {
    ppi_res <- tryCatch(auto_add_ppi_columns(training_df), error = function(e) NULL)
    if (!is.null(ppi_res) && is.list(ppi_res) && !is.null(ppi_res$df)) training_df <- ppi_res$df
  }
  # Filter training years to the analysis window
  cat("Filtering training data for years 1985-2025...\n")
  training_df <- training_df |> dplyr::filter(year >= 1985 & year <= 2025)
  cat("Training rows after year filtering:", nrow(training_df), "\n")
} else {
  cat("Training CSV not found at:", TRAINING_CSV, "- falling back to using inference CSV for statistics (not recommended).\n")
  training_df <- df
}

# Join df with gpts_map to add Veg column
if ("location_id" %in% names(df)) df$location_id <- as.character(df$location_id) else cat("Warning: input CSV missing 'location_id' column; skipping cast.\n")
cat("Note: No Veg/no-soil joining performed; mapping CSV only used to estimate barren DVI when present.\n")
# df <- normalize_no_soil_col(df)  # Function removed - deprecated
## Enforce policy: if a no_soil/no.soil/no\ssoi l column exists it must contain at least
## one non-NA value. A present-but-all-NA column is a data error and should be
## fixed by the data provider (do not silently fall back to index inference).
if ("no soil" %in% names(df)) {
  n_non_na_no_soil <- sum(!is.na(df$`no soil`))
  if (n_non_na_no_soil == 0) stop("CRITICAL: 'no soil' column is present but contains no non-NA values. Provide valid 'no soil' (or 'no_soil') values, or remove the column to allow index-based inference.")
}
if ("Veg" %in% names(df)) {
  df$Veg <- tolower(df$Veg)
  cat("DEBUG: Unique Veg types in df:\n")
  print(table(df$Veg, useNA = "ifany"))
}

# Ensure date column is Date type
df$date <- as.Date(df$date)

# Extract year and month
df$year <- lubridate::year(df$date)
if (!"pheno_year" %in% names(df)) df$pheno_year <- ifelse(lubridate::month(df$date) >= 3, lubridate::year(df$date), lubridate::year(df$date) - 1)
df$month <- lubridate::month(df$date)

# Filter for years 1985-2025
cat("Filtering data for years 1985-2025...\n")
df <- df |> dplyr::filter(year >= 1985 & year <= 2025)
cat("Data rows after year filtering:", nrow(df), "\n")

# Filter for years 1985-2025
cat("Filtering data for years 1985-2025...\n")
df <- df |> dplyr::filter(year >= 1985 & year <= 2025)
cat("Data rows after year filtering:", nrow(df), "\n")

# Calculate all indices of interest
df <- calculate_indices(df)

# Optionally prefer mitmat/ppi package if installed
if (!exists("ppi", mode = "function") && requireNamespace("ppi", quietly = TRUE)) {
  # Use upstream implementation from the 'ppi' package if available
  assign("ppi", ppi::ppi, envir = globalenv())
  cat("Using 'ppi' function from installed package 'ppi'\n")
}

# Compute DVI if not present so we can estimate soil DVI from barren rows
if (!"DVI" %in% names(df) && all(c("nir", "red") %in% names(df))) {
  df$DVI <- df$nir - df$red
  cat("Computed DVI in january_averages.R for soil DVI estimation\n")
}

# Calculate Solar Zenith Angle (SZA) per-row if missing, using per-location latitude
# PPI relies on correct zenith angles (radians). We prefer per-row 'lat'.
if (!"zenith.angle" %in% names(df)) {
  if ("lat" %in% names(df) && any(is.finite(df$lat))) {
    cat("Calculating Solar Zenith Angle (SZA) per row using 'lat' and 'date'...\n")
    if (!"doy" %in% names(df)) df$doy <- lubridate::yday(df$date)
    # calculate_solar_zenith returns radians
    df$zenith.angle <- calculate_solar_zenith(as.numeric(df$lat), df$doy)
  } else {
    cat("No per-row 'lat' available; leaving 'zenith.angle' NA. PPI will not be computed without SZA.\n")
    df$zenith.angle <- NA_real_
  }
} else {
  # If zenith.angle exists but appears in degrees (> pi), convert to radians
  if (any(is.finite(df$zenith.angle) & df$zenith.angle > pi, na.rm = TRUE)) {
    cat("Converting provided 'zenith.angle' from degrees to radians...\n")
    df$zenith.angle <- df$zenith.angle * pi / 180
  }
}

# Estimate DVI soil from barren plots in the input CSV if available.
## NOTE: We prefer to estimate soil DVI from empirical barren observations.
## First, try to get barren observations from the mapping CSV (if loaded as map_df),
## using Veg == 'barren'. If that fails, fall back to barren rows in the main CSV.
dvi_soil_est <- NA_real_

# Helper to compute DVI column robustly (case-insensitive)
compute_dvi_col <- function(tbl) {
  lc <- tolower(names(tbl))
  if ("dvi" %in% lc) return(tbl[[which(lc == "dvi")[1]]])
  # prefer NIR - Red (case-insensitive)
  if ("nir" %in% lc && "red" %in% lc) {
    nir_col <- which(lc == "nir")[1]
    red_col <- which(lc == "red")[1]
    return(tbl[[nir_col]] - tbl[[red_col]])
  }
  return(NULL)
}

# Try mapping CSV first
if (exists("map_df")) {
  map_dvi <- compute_dvi_col(map_df)
  if (!is.null(map_dvi) && "Veg" %in% names(map_df)) {
    barren_idx_map <- tolower(as.character(map_df$Veg)) == "barren"
    if (any(barren_idx_map, na.rm = TRUE)) {
      mean_map_dvi <- mean(map_dvi[barren_idx_map], na.rm = TRUE)
      if (is.finite(mean_map_dvi)) {
        dvi_soil_est <- mean_map_dvi
        assign("PPI_DVI_SOIL", dvi_soil_est, envir = globalenv())
        cat(sprintf("Estimated PPI_DVI_SOIL from mapping CSV (veg=='barren') = %.4f\n", dvi_soil_est))
        if (!"DVI" %in% names(df)) df$DVI <- compute_dvi_col(df)
        df <- add_ppi_columns(df, dvi_soil = dvi_soil_est)
      }
    }
  }
}

# If mapping CSV provided no barren estimate, fall back to main CSV barren rows
if (is.na(dvi_soil_est) && "DVI" %in% names(df)) {
  barren_idx <- FALSE
  if ("Veg" %in% names(df)) barren_idx <- tolower(as.character(df$Veg)) == "barren"


  if (any(barren_idx, na.rm = TRUE)) {
    dvi_soil_est <- mean(df$DVI[barren_idx], na.rm = TRUE)
    assign("PPI_DVI_SOIL", dvi_soil_est, envir = globalenv())
    cat(sprintf("Estimated PPI_DVI_SOIL (from main CSV barren rows) = %.4f\n", dvi_soil_est))
    df <- add_ppi_columns(df, dvi_soil = dvi_soil_est)
  } else {
    cat("[PPI] No barren observations found in main CSV and mapping CSV — skipping PPI calculation (no empirical baseline).\n")
    if (!"zenith.angle" %in% names(df)) df$zenith.angle <- NA_real_
    df$PPI <- NA_real_
  }
} else if (is.na(dvi_soil_est)) {
  cat("[PPI] DVI not computed; skipping PPI calculation. Ensure 'nir' and 'red' bands are present or mapping CSV contains barren samples.\n")
  if (!"zenith.angle" %in% names(df)) df$zenith.angle <- NA_real_
  df$PPI <- NA_real_
}

# Outlier detection and removal has been removed per user request; PPI values are kept unchanged.
if ("PPI" %in% names(df)) {
  cat("Note: PPI outlier removal disabled per user request; keeping raw PPI values.\n")
}

# Filter for January, July, September data from TRAINING data (not inference)
january_data <- training_df |> dplyr::filter(month == 1)
july_data <- training_df |> dplyr::filter(month == 7)
september_data <- training_df |> dplyr::filter(month == 9)

cat("Number of January observations:", nrow(january_data), "\n")
cat("Number of July observations:", nrow(july_data), "\n")
cat("Number of September observations:", nrow(september_data), "\n")

# Compute dataset-level SNR for indices (per-index SNR across locations)
compute_global_index_snr <- function(df, indices, group_col = "location_id", eps = 1e-8) {
  # Compute SNR per location-year as (July_mean - January_mean) / MAD(within-year),
  # then take median across all location-years to robustly estimate seasonal SNR.
  out <- numeric(length(indices)); names(out) <- indices
  if (!"month" %in% names(df) && "date" %in% names(df)) df$month <- lubridate::month(df$date)
  if (!"pheno_year" %in% names(df) && "date" %in% names(df)) df$pheno_year <- assign_pheno_year(df$date)
  for (idx in indices) {
    if (!idx %in% names(df)) { out[idx] <- NA_real_; next }
    s_vals <- numeric(0)
    locs <- unique(na.omit(as.character(df[[group_col]])))
    for (loc in locs) {
      sub_loc <- df[df[[group_col]] == loc, , drop = FALSE]
      if (nrow(sub_loc) == 0) next
      years <- split(sub_loc, sub_loc$pheno_year)
      for (yr in years) {
        if (!idx %in% names(yr) || !"month" %in% names(yr)) next
        jan_vals <- yr[[idx]][yr$month == 1]; jul_vals <- yr[[idx]][yr$month == 7]
        jan_vals <- jan_vals[is.finite(jan_vals)]; jul_vals <- jul_vals[is.finite(jul_vals)]
        if (length(jan_vals) < 1 || length(jul_vals) < 1) next
        amp <- mean(jul_vals, na.rm = TRUE) - mean(jan_vals, na.rm = TRUE)
        vals_year <- yr[[idx]]; vals_year <- vals_year[is.finite(vals_year)]
        if (length(vals_year) < 2) next
        noise <- tryCatch(mad(vals_year, na.rm = TRUE), error = function(e) NA_real_)
        if (!is.finite(amp) || !is.finite(noise) || noise <= 0) next
        s_vals <- c(s_vals, amp / (noise + eps))
      }
    }
    if (length(s_vals) == 0) { out[idx] <- NA_real_; next }
    out[idx] <- median(s_vals, na.rm = TRUE)
  }
  out
}

# Calculate SNR for the indices of interest using training data
indices_to_snr <- INDICES_OF_INTEREST
index_snr <- compute_global_index_snr(training_df, indices_to_snr, group_col = "location_id")
## For backward-compatible SNR reporting (and to match prior reported values),
## compute SNR as location-wise amplitude (range) / noise (MAD) on non-barren
## training samples and report that value in the summary table.
compute_global_index_snr_simple <- function(df, indices, group_col = "location_id", eps = 1e-8) {
  out <- numeric(length(indices)); names(out) <- indices
  for (idx in indices) {
    if (!idx %in% names(df)) { out[idx] <- NA_real_; next }
    locs <- unique(na.omit(as.character(df[[group_col]])))
    amps <- numeric(0); noises <- numeric(0)
    for (loc in locs) {
      sub <- df[df[[group_col]] == loc, , drop = FALSE]
      vals <- sub[[idx]]
      vals <- vals[is.finite(vals)]
      if (length(vals) < 2) next
      amp_loc <- diff(range(vals, na.rm = TRUE))
      noise_loc <- tryCatch(mad(vals, na.rm = TRUE), error = function(e) NA_real_)
      if (is.finite(amp_loc)) amps <- c(amps, amp_loc)
      if (is.finite(noise_loc)) noises <- c(noises, noise_loc)
    }
    if (length(amps) == 0 || length(noises) == 0) { out[idx] <- NA_real_; next }
    out[idx] <- median(amps, na.rm = TRUE) / (median(noises, na.rm = TRUE) + eps)
  }
  out
}

training_nonbarren <- training_df
if ("Veg" %in% names(training_df)) training_nonbarren <- training_df[!(tolower(training_df$Veg) == "barren"), , drop = FALSE]
index_snr <- compute_global_index_snr_simple(training_nonbarren, indices_to_snr, group_col = "location_id")
cat("Index SNR (amplitude/MAD on non-barren training data):\n")
print(index_snr)
cat("Index SNR computed — will be written to the results folder if available when saving outputs\n")

# Helper function to convert day-of-year to pentad (5-day intervals)
doy_to_pentad <- function(doy) {
  pmin(ceiling(as.numeric(doy) / 5), 73)
}

# Function to bootstrap averages using location resampling
# Resamples unique location_ids, then takes all observations for those locations
bootstrap_hierarchical_means <- function(df, metrics = c("MSAVI", "NDVI", "PPI"), group_col = "location_id", B = 1000) {
  # Validation
  if (!group_col %in% names(df)) return(rep(NA, length(metrics) * 2))

  # Work only with requested columns for speed
  df <- df |> dplyr::select(all_of(c(group_col, metrics)))
  df[[group_col]] <- as.character(df[[group_col]])
  ids <- unique(df[[group_col]])
  n_ids <- length(ids)
  if (n_ids < 2) return(rep(NA, length(metrics) * 2))

  # Pre-split into matrices per location (faster than repeated dataframe subsetting)
  data_map <- lapply(ids, function(id) {
    mat <- as.matrix(df[df[[group_col]] == id, metrics, drop = FALSE])
    # Ensure matrix has correct number of columns even if single column was present
    if (is.null(dim(mat))) mat <- matrix(mat, ncol = length(metrics))
    colnames(mat) <- metrics
    mat
  })
  names(data_map) <- ids

  # Run hierarchical bootstrap: Resample locations, then resample rows within each selected location
  boot_replicates <- replicate(B, {
    # Stage 1: resample locations (with replacement)
    sel_ids <- sample(ids, n_ids, replace = TRUE)
    selected_mats <- data_map[sel_ids]

    # Stage 2: for each selected location matrix, resample rows and compute column means
    cluster_means <- vapply(selected_mats, function(mat) {
      n_obs <- nrow(mat)
      if (is.na(n_obs) || n_obs <= 0) return(rep(NA_real_, length(metrics)))
      rows <- sample.int(n_obs, n_obs, replace = TRUE)
      colMeans(mat[rows, , drop = FALSE], na.rm = TRUE)
    }, numeric(length(metrics)))
    # Ensure cluster_means is a matrix (metrics x n_clusters)
    if (is.null(dim(cluster_means))) cluster_means <- matrix(cluster_means, nrow = length(metrics))

    # Final statistic: mean across clusters (metrics x 1)
    rowMeans(cluster_means, na.rm = TRUE)
  })

  # boot_replicates is [num_metrics x B] (transpose if needed)
  if (is.null(dim(boot_replicates))) {
    boot_replicates <- matrix(boot_replicates, nrow = length(metrics))
  }

  cis <- apply(boot_replicates, 1, function(x) stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))

  # Format output: metric_ci_lower, metric_ci_upper
  output <- as.vector(cis)
  metric_names <- rep(metrics, each = 2)
  ci_types <- rep(c("ci_lower", "ci_upper"), times = length(metrics))
  names(output) <- paste(metric_names, ci_types, sep = "_")

  return(output)
}

# Function to bootstrap medians using location resampling
# Resamples unique location_ids, then takes all observations for those locations
bootstrap_hierarchical_medians <- function(df, metrics = c("PPI"), group_col = "location_id", B = 1000) {
  if (!group_col %in% names(df)) return(rep(NA, length(metrics) * 2))

  df <- df |> dplyr::select(all_of(c(group_col, metrics)))
  df[[group_col]] <- as.character(df[[group_col]])
  ids <- unique(df[[group_col]])
  n_ids <- length(ids)
  if (n_ids < 2) return(rep(NA, length(metrics) * 2))

  data_map <- lapply(ids, function(id) {
    mat <- as.matrix(df[df[[group_col]] == id, metrics, drop = FALSE])
    if (is.null(dim(mat))) mat <- matrix(mat, ncol = length(metrics))
    colnames(mat) <- metrics
    mat
  })
  names(data_map) <- ids

  boot_replicates <- replicate(B, {
    sel_ids <- sample(ids, n_ids, replace = TRUE)
    selected_mats <- data_map[sel_ids]

    cluster_stats <- vapply(selected_mats, function(mat) {
      n_obs <- nrow(mat)
      if (is.na(n_obs) || n_obs <= 0) return(rep(NA_real_, length(metrics)))
      rows <- sample.int(n_obs, n_obs, replace = TRUE)
      apply(mat[rows, , drop = FALSE], 2, median, na.rm = TRUE)
    }, numeric(length(metrics)))
    # Ensure cluster_stats is a matrix (metrics x n_clusters)
    if (is.null(dim(cluster_stats))) cluster_stats <- matrix(cluster_stats, nrow = length(metrics))

    # Aggregate across clusters by median
    apply(cluster_stats, 1, median, na.rm = TRUE)
  })

  if (is.null(dim(boot_replicates))) boot_replicates <- matrix(boot_replicates, nrow = length(metrics))

  cis <- apply(boot_replicates, 1, function(x) stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
  output <- as.vector(cis)
  metric_names <- rep(metrics, each = 2)
  ci_types <- rep(c("lower", "upper"), times = length(metrics))
  names(output) <- paste(metric_names, ci_types, sep = "_")
  return(output)
}

# Helper to calculate mean of location means
mean_of_means <- function(vals, ids) {
  if (length(vals) == 0) return(NA_real_)
  # Calculate mean for each unique ID
  loc_means <- tapply(vals, ids, mean, na.rm = TRUE)
  # Calculate mean of those means
  mean(loc_means, na.rm = TRUE)
}

# Add bootstrapping for uncertainty
set.seed(123)
B <- 1000  # Number of bootstrap replicates

# Helper to process global averages
process_global_averages <- function(data, month_name) {
  cat(sprintf("Calculating Global %s Averages with Hierarchical Bootstrapping...\n", month_name))
  global_boot <- bootstrap_hierarchical_means(data, metrics = INDICES_OF_INTEREST, B = B)
  
  avg_stats <- data |> 
    dplyr::summarize(
      n_observations = dplyr::n(),
      n_locations = dplyr::n_distinct(location_id),
      across(all_of(INDICES_OF_INTEREST), ~ mean_of_means(.x, location_id), .names = "avg_{.col}")
    )
  
  # Bind CIs
  # global_boot is a named vector. Convert to 1-row DF.
  ci_df <- as.data.frame(t(global_boot))
  dplyr::bind_cols(avg_stats, ci_df)
}

# Calculate Global Averages (Aggregated 2020-2024) for January
jan_global_avg <- process_global_averages(january_data, "January")

# Calculate Global Averages (Aggregated 2020-2024) for July
july_global_avg <- process_global_averages(july_data, "July")

# Calculate Global Averages (Aggregated 2020-2024) for September
sept_global_avg <- process_global_averages(september_data, "September")

# --- Averages by Vegetation Type (using mapping CSV at location-level) ---

# Build per-location summaries from observation-level data (do not join mapping to observations)
location_summary <- function(data, metrics = c("MSAVI", "NDVI", "PPI")) {
  data |> dplyr::group_by(location_id) |> dplyr::summarize(
    n_observations = dplyr::n(),
    across(any_of(metrics), ~ mean(.x, na.rm = TRUE), .names = "avg_{.col}"),
    .groups = "drop"
  )
}

loc_jan <- location_summary(january_data, metrics = INDICES_OF_INTEREST)
loc_july <- location_summary(july_data, metrics = INDICES_OF_INTEREST)
loc_sept <- location_summary(september_data, metrics = INDICES_OF_INTEREST)

# Prepare mapping at location level (map_df should already be loaded if available)
map_loc <- NULL
if (exists("map_df")) {
  map_loc <- map_df
  if (!"location_id" %in% names(map_loc) && all(c("lon", "lat") %in% names(map_loc))) {
    map_loc$location_id <- make_location_id(map_loc$lon, map_loc$lat)
  }
  if ("location_id" %in% names(map_loc) && "Veg" %in% names(map_loc)) {
    map_loc$location_id <- as.character(map_loc$location_id)
    map_loc$Veg <- tolower(as.character(map_loc$Veg))
    map_loc <- map_loc |> dplyr::select(location_id, Veg)
  } else {
    map_loc <- NULL
  }
}

# Helper to compute veg-level stats from location summaries
veg_stats_from_locs <- function(loc_tbl, data_tbl, period_name) {
  if (is.null(map_loc)) {
    cat(sprintf("Mapping CSV not available — skipping %s veg-type summaries.\n", period_name))
    return(data.frame())
  }

  loc_joined <- loc_tbl |> dplyr::left_join(map_loc, by = "location_id")
  loc_joined <- loc_joined |> dplyr::filter(!is.na(Veg))

  # Per-veg location-level means and counts
  veg_loc_stats <- loc_joined |> dplyr::group_by(Veg) |> dplyr::summarize(
    n_locations = dplyr::n(),
    n_observations = sum(n_observations, na.rm = TRUE),
    across(starts_with("avg_"), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

  # Bootstrap CIs per veg using the original observation-level data restricted to locations in each veg
  veg_boot <- lapply(unique(loc_joined$Veg), function(v) {
    locs <- loc_joined |> dplyr::filter(Veg == v) |> dplyr::pull(location_id)
    sub <- data_tbl |> dplyr::filter(location_id %in% locs)
    
    # Run bootstrap (returns named vector)
    res <- bootstrap_hierarchical_means(sub, metrics = INDICES_OF_INTEREST, B = B)
    
    # Convert res vector to 1-row DF
    res_df <- as.data.frame(t(res))
    res_df$Veg <- v
    res_df
  })
  
  if (length(veg_boot) > 0) {
    veg_boot_df <- dplyr::bind_rows(veg_boot)
  } else {
    veg_boot_df <- data.frame(Veg = character(0))
  }

  dplyr::left_join(veg_loc_stats, veg_boot_df, by = "Veg")
}

cat("Calculating Vegetation Type January Averages using mapping CSV and location summaries...\n")
veg_type_averages <- veg_stats_from_locs(loc_jan, january_data, "January")

cat("Calculating Vegetation Type July Averages using mapping CSV and location summaries...\n")
veg_type_averages_july <- veg_stats_from_locs(loc_july, july_data, "July")

cat("Calculating Vegetation Type September Averages using mapping CSV and location summaries...\n")
veg_type_averages_sept <- veg_stats_from_locs(loc_sept, september_data, "September")

# Print results
cat("\n=== GLOBAL JANUARY AVERAGES (1985-2025) ===\n")
print(jan_global_avg)

cat("\n=== GLOBAL JULY AVERAGES (1985-2025) ===\n")
print(july_global_avg)

cat("\n=== GLOBAL SEPTEMBER AVERAGES (1985-2025) ===\n")
print(sept_global_avg)

cat("\n=== JANUARY AVERAGES BY VEGETATION TYPE (1985-2025) ===\n")
print(veg_type_averages)

cat("\n=== JULY AVERAGES BY VEGETATION TYPE (1985-2025) ===\n")
print(veg_type_averages_july)

cat("\n=== SEPTEMBER AVERAGES BY VEGETATION TYPE (1985-2025) ===\n")
print(veg_type_averages_sept)

# --- Create Summary Table as requested ---
# Create global averages excluding 'barren' vegetation (for summary table)
exclude_barren_filter <- function(d) {
  if (is.null(d) || nrow(d) == 0) return(d)
  if (!"Veg" %in% names(d)) return(d)
  idx <- tolower(as.character(d$Veg)) == "barren"
  d[is.na(idx) | !idx, , drop = FALSE]
}

nonbarren_jan <- exclude_barren_filter(january_data)
nonbarren_july <- exclude_barren_filter(july_data)
nonbarren_sept <- exclude_barren_filter(september_data)

jan_global_avg_nonbarren <- if (nrow(nonbarren_jan) > 0) process_global_averages(nonbarren_jan, "January (Non-Barren)") else process_global_averages(nonbarren_jan[0,], "January (Empty)")
july_global_avg_nonbarren <- if (nrow(nonbarren_july) > 0) process_global_averages(nonbarren_july, "July (Non-Barren)") else process_global_averages(nonbarren_july[0,], "July (Empty)")
sept_global_avg_nonbarren <- if (nrow(nonbarren_sept) > 0) process_global_averages(nonbarren_sept, "September (Non-Barren)") else process_global_averages(nonbarren_sept[0,], "September (Empty)")
create_summary_entry <- function(idx_name, snr_val, jan_row, july_row, sept_row) {
  # Helper to format mean +/- margin
  fmt_val <- function(mean_val, lower, upper) {
    if (is.na(mean_val)) return("NA")
    margin <- (upper - lower) / 2
    # Format with comma decimal
    m_str <- format(round(mean_val, 3), nsmall=3, decimal.mark=",")
    e_str <- format(round(margin, 3), nsmall=3, decimal.mark=",")
    sprintf("%s (+/-%s)*", m_str, e_str)
  }
  
  jan_str <- fmt_val(jan_row[[paste0("avg_", idx_name)]], jan_row[[paste0(idx_name, "_ci_lower")]], jan_row[[paste0(idx_name, "_ci_upper")]])
  july_str <- fmt_val(july_row[[paste0("avg_", idx_name)]], july_row[[paste0(idx_name, "_ci_lower")]], july_row[[paste0(idx_name, "_ci_upper")]])
  sept_str <- fmt_val(sept_row[[paste0("avg_", idx_name)]], sept_row[[paste0(idx_name, "_ci_lower")]], sept_row[[paste0(idx_name, "_ci_upper")]])
  
  # Format SNR
  if (!is.na(snr_val) && snr_val > 10000) {
     exp_val <- floor(log10(snr_val))
     base_val <- snr_val / (10^exp_val)
     snr_str <- sprintf("%.2f * 10^%d", base_val, exp_val)
     snr_str <- gsub("\\.", ",", snr_str)
  } else if (!is.na(snr_val)) {
     snr_str <- format(round(snr_val, 2), nsmall=2, decimal.mark=",")
  } else {
     snr_str <- "NA"
  }
  
  c(idx_name, snr_str, jan_str, july_str, sept_str)
}

summary_rows <- lapply(indices_to_snr, function(idx) {
  create_summary_entry(idx, index_snr[idx], jan_global_avg_nonbarren, july_global_avg_nonbarren, sept_global_avg_nonbarren)
})

summary_table <- do.call(rbind, summary_rows)
colnames(summary_table) <- c("INDEX", "SNR", "January Average (All areas, excl. barren)", "July average (All areas, excl. barren)", "September average (All areas, excl. barren)")
summary_table <- as.data.frame(summary_table)

cat("\n=== SUMMARY TABLE (excluding Veg=='barren') ===\n")
print(summary_table)

# Optional: Save results to Excel
output_dir <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/phenology_results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

results_list <- list(
  "Jan_Global_Avg" = jan_global_avg,
  "July_Global_Avg" = july_global_avg,
  "Sept_Global_Avg" = sept_global_avg,
  "Jan_VegType_Avg" = veg_type_averages,
  "July_VegType_Avg" = veg_type_averages_july,
  "Sept_VegType_Avg" = veg_type_averages_sept,
  "Summary_Table" = summary_table
)
# Add index SNR results to the output bundle
if (exists("index_snr") && length(index_snr) > 0) {
  results_list$Index_SNR <- data.frame(index = names(index_snr), snr = as.numeric(index_snr), stringsAsFactors = FALSE)
}

openxlsx::write.xlsx(results_list, file = file.path(output_dir, "phenology_averages.xlsx"))
if (exists("index_snr") && length(index_snr) > 0) {
  # Save the SNR CSV alongside the Excel file
  try({
    write.csv(data.frame(index = names(index_snr), snr = as.numeric(index_snr)), file = file.path(output_dir, "index_snr.csv"), row.names = FALSE)
    cat("Index SNR saved to:", file.path(output_dir, "index_snr.csv"), "\n")
  }, silent = TRUE)
}

cat("\nResults saved to:", file.path(output_dir, "phenology_averages.xlsx"), "\n")

# Plot timeseries for all indices with seasonal normalization and bootstrapping
indices_to_plot <- INDICES_OF_INTEREST

plot_data_mean <- list()
plot_data_median <- list()

for (idx in indices_to_plot) {
  cat(sprintf("Processing summer trend data for %s...\n", idx))
  
  # Filter summer data
  if (idx == "PPI") {
    summer_data <- df |> dplyr::filter(month %in% 6:9, .data[[idx]] > 0)
  } else {
    summer_data <- df |> dplyr::filter(month %in% 6:9, is.finite(.data[[idx]]))
  }
  
  if (nrow(summer_data) == 0) {
    cat(sprintf("Warning: No summer %s data found. Skipping.\n", idx))
    next
  }
  
  # --- Normalization Step ---
  # Calculate Day of Year (DOY) if not present
  if (!"doy" %in% names(summer_data)) {
    summer_data$doy <- lubridate::yday(summer_data$date)
  }
  
  # Fit a 3rd degree polynomial to capture the seasonal curve across the summer months
  seasonal_model <- lm(as.formula(paste(idx, "~ poly(doy, 3)")), data = summer_data)
  
  # Predict the seasonal trend for each observation
  summer_data$seasonal_trend <- predict(seasonal_model, newdata = summer_data)
  
  # Calculate the global mean of the seasonal trend to maintain the overall magnitude
  global_seasonal_mean <- mean(summer_data$seasonal_trend, na.rm = TRUE)
  
  # Normalize: Subtract the seasonal deviation from the global mean
  summer_data[[paste0(idx, "_norm")]] <- summer_data[[idx]] - (summer_data$seasonal_trend - global_seasonal_mean)
  
  # --- Bootstrapping with Normalized Data ---
  # Bootstrap means for each year using the normalized index
  summer_yearly_boot <- summer_data |> 
    dplyr::group_by(pheno_year) |> 
    dplyr::do({
      # Note: We pass the normalized metric to bootstrap
      norm_metric <- paste0(idx, "_norm")
      res <- bootstrap_hierarchical_means(., metrics = c(norm_metric), B = B)
      
      # Compute point estimate as mean of location-level means using normalized
      loc_means <- sapply(unique(.$location_id), function(id) {
        sub <- .[.$location_id == id, ]
        mean(sub[[norm_metric]], na.rm = TRUE)
      })
      mean_val <- mean(loc_means, na.rm = TRUE)
      
      data.frame(
        mean_val = mean_val,
        ci_lower = res[1],
        ci_upper = res[2]
      )
    }) |> 
    dplyr::ungroup()
  
  summer_yearly_boot$index <- idx
  plot_data_mean[[idx]] <- summer_yearly_boot
  
  # --- Median Bootstrapping ---
  summer_yearly_boot_median <- summer_data |> 
    dplyr::group_by(pheno_year) |> 
    dplyr::do({
      norm_metric <- paste0(idx, "_norm")
      res <- bootstrap_hierarchical_medians(., metrics = c(norm_metric), B = B)
      
      # Compute point estimate as median of location-level medians using normalized
      loc_medians <- sapply(unique(.$location_id), function(id) {
        sub <- .[.$location_id == id, ]
        median(sub[[norm_metric]], na.rm = TRUE)
      })
      median_val <- median(loc_medians, na.rm = TRUE)
      
      data.frame(
        median_val = median_val,
        ci_lower = res[1],
        ci_upper = res[2]
      )
    }) |> 
    dplyr::ungroup()
  
  summer_yearly_boot_median$index <- idx
  plot_data_median[[idx]] <- summer_yearly_boot_median
}

# Combine data
plot_data_mean <- do.call(rbind, plot_data_mean)
plot_data_median <- do.call(rbind, plot_data_median)

# Plot mean
p <- ggplot(plot_data_mean, aes(x = pheno_year, y = mean_val)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "blue") +
  facet_wrap(~ index, scales = "free_y") +
  labs(title = "Mean Seasonally-Normalized June-Sept Vegetation Indices (1985-2025)",
       subtitle = "Normalization: Index - SeasonalTrend(doy) + Mean(SeasonalTrend)",
       x = "Year",
       y = "Mean Normalized Value") +
  theme_minimal()

print(p)

# Save the plot
ggsave(file.path(output_dir, "all_indices_summer_trend_normalized_mean.png"), plot = p, width = 16, height = 12)
cat("Saved mean plot to:", file.path(output_dir, "all_indices_summer_trend_normalized_mean.png"), "\n")

# Plot median
p_med <- ggplot(plot_data_median, aes(x = pheno_year, y = median_val)) +
  geom_line(color = "darkgreen") +
  geom_point(color = "orange") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "darkgreen") +
  facet_wrap(~ index, scales = "free_y") +
  labs(title = "Median Seasonally-Normalized June-Sept Vegetation Indices (1985-2025)",
       subtitle = "Normalization: Index - SeasonalTrend(doy) + Mean(SeasonalTrend)",
       x = "Year",
       y = "Median Normalized Value") +
  theme_minimal()

print(p_med)

# Save the plot
ggsave(file.path(output_dir, "all_indices_summer_trend_normalized_median.png"), plot = p_med, width = 16, height = 12)
cat("Saved median plot to:", file.path(output_dir, "all_indices_summer_trend_normalized_median.png"), "\n")