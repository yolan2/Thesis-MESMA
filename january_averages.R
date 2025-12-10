# Script to calculate January averages of MSAVI, NDVI, and PPI for each vegetation type
# with bootstrapping for uncertainty estimation

# Load required libraries
library(dplyr)
library(lubridate)
library(sf)
library(boot)
library(openxlsx)
library(ggplot2)

# Helper functions
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

remove_outliers_mad <- function(x, k = 3) {
  if (all(is.na(x))) return(x)
  
  # Only consider positive values for statistics to avoid zero-inflation skewing the MAD to 0
  # This prevents removing all vegetation signals when the majority of the area is barren (PPI=0)
  pos_x <- x[x > 0 & !is.na(x)]
  
  # If not enough positive data points, return original vector
  if (length(pos_x) < 5) return(x)
  
  med <- median(pos_x, na.rm = TRUE)
  mad_val <- mad(pos_x, na.rm = TRUE)
  
  # If MAD is 0 (e.g. all positive values are identical), skip filtering
  if (mad_val == 0) return(x)
  
  # Calculate upper bound based on positive values distribution
  upper <- med + k * mad_val
  
  # Only filter high outliers. We keep 0s and low values.
  x[!is.na(x) & x > upper] <- NA
  x
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

## Use centralized PPI helper to ensure identical PPI calculation
if (file.exists("ppi_helpers.R")) {
  source("ppi_helpers.R")
} else {
  warning("ppi_helpers.R not found; falling back to local inline PPI logic. Please ensure ppi_helpers.R is present in project root for consistent PPI calculations.")
}

# Define input file paths
INPUT_CSV <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/landsat_lower_inference.csv"
GEOJSON_PATH <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/updated_zuizer_zonder_foto_UTM.geojson"

# Load the phenology data
cat("Loading data from:", INPUT_CSV, "\n")
df <- readr::read_csv(INPUT_CSV, show_col_types = FALSE)

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

# Load GeoJSON for vegetation types
cat("Loading GeoJSON from:", GEOJSON_PATH, "\n")
gpts_raw <- tryCatch({
  sf::st_read(GEOJSON_PATH, quiet = TRUE)
}, error = function(e) {
  cat("Warning: Failed to read GeoJSON (", e$message, "). Skipping Veg mapping.\n")
  NULL
})

# Calculate mean latitude of the study area for zenith angle fallback
# Transform to WGS84 to ensure we have lat/lon
if (!is.null(gpts_raw)) {
  tryCatch({
    gpts_wgs84 <- sf::st_transform(gpts_raw, 4326)
    coords_wgs84 <- sf::st_coordinates(gpts_wgs84)
    mean_lat_study_area <- mean(coords_wgs84[, 2], na.rm = TRUE)
    cat("Mean latitude of study area:", mean_lat_study_area, "\n")
  }, error = function(e) {
  cat("Warning: Could not calculate mean latitude from GeoJSON:", e$message, "\n")
  })
} else {
  cat("GeoJSON unavailable; mean latitude fallback not computed\n")
}

# Process GeoJSON to get Veg mapping
geojson_names <- names(gpts_raw)
normalized_names <- gsub("[^a-z0-9]+", "_", tolower(geojson_names))
no_soil_col <- geojson_names[normalized_names == "no_soil"]
if (length(no_soil_col) == 0) stop("GeoJSON point data requires a 'no soil' column")
no_soil_raw <- gpts_raw[[no_soil_col[1]]]
no_soil_vals <- safe_as_numeric(no_soil_raw)
gpts_raw$`.__no soil__` <- no_soil_vals

matched_cols <- names(gpts_raw)[tolower(names(gpts_raw)) %in% c("vegetation", "veg", "class")]
if (length(matched_cols) > 0) {
  veg_col_orig <- matched_cols[1]
  gpts_raw$.__veg__ <- as.character(gpts_raw[[veg_col_orig]])
} else {
  gpts_raw$.__veg__ <- NA_character_
}
gpts_wgs84 <- sf::st_transform(gpts_raw, 4326)
coords_wgs84 <- sf::st_coordinates(gpts_wgs84)
if (!is.null(gpts_raw)) {
  gpts_raw$.__lon__ <- coords_wgs84[, 1]
  gpts_raw$.__lat__ <- coords_wgs84[, 2]
  gpts_raw$location_id_geo <- make_location_id(gpts_raw$.__lon__, gpts_raw$.__lat__)
  gpts_raw$location_id <- as.character(seq_len(nrow(gpts_raw)))
}

gdf <- if (!is.null(gpts_raw)) sf::st_drop_geometry(gpts_raw) else NULL
# Ensure we have simple lon/lat and a location ID in the gdf
if (!is.null(gdf) && (!".__lon__" %in% names(gdf) || !".__lat__" %in% names(gdf))) {
  # try transforming geometry and reading coordinates if __lon__/__lat__ missing
  tryCatch({
    wgs <- sf::st_transform(gpts_raw, 4326)
    coords_tmp <- sf::st_coordinates(wgs)
    gdf$.__lon__ <- coords_tmp[,1]
    gdf$.__lat__ <- coords_tmp[,2]
  }, error = function(e) {
    # leave as-is if fails
  })
}
if (!is.null(gdf)) {
  gdf$location_id_geo <- make_location_id(gdf$.__lon__, gdf$.__lat__)
  gdf$location_id <- as.character(seq_len(nrow(gdf)))
  cat("DEBUG: First 5 location_id values in GeoJSON (sequential):", paste(head(gdf$location_id, 5), collapse=", "), "\n")
  cat("DEBUG: First 5 location_id_geo values in GeoJSON:", paste(head(gdf$location_id_geo, 5), collapse=", "), "\n")
}
gpts_map <- if (!is.null(gdf)) {
  # Try to determine which ID to use for mapping
  # If CSV has "1", "2", etc., we should use the sequential location_id from GeoJSON
  # If CSV has "L_...", we should use location_id_geo
  
  # Default to location_id_geo as before, but we might change this based on CSV inspection
  # But wait, we can't see CSV here easily inside this block if we want to be dynamic.
  # Let's just include BOTH and let the join decide or prepare multiple maps?
  # Or better: check what df$location_id looks like.
  
  # For now, let's assume we want to match the CSV.
  # If df$location_id looks like integers, use sequential.
  
  use_seq_id <- FALSE
  if (exists("df") && "location_id" %in% names(df)) {
     sample_id <- as.character(df$location_id[1])
     if (!grepl("^L_", sample_id)) {
        use_seq_id <- TRUE
        cat("DEBUG: CSV location_id does not start with 'L_', assuming sequential/integer IDs. Using sequential ID from GeoJSON for join.\n")
     }
  }
  
  if (use_seq_id) {
      gdf %>%
      dplyr::select(location_id, Veg = .__veg__, `no soil` = `.__no soil__`) %>%
      dplyr::mutate(location_row = as.character(seq_len(dplyr::n()))) %>%
      dplyr::distinct(location_id, .keep_all = TRUE)
  } else {
      gdf %>%
      dplyr::select(location_id = location_id_geo, Veg = .__veg__, `no soil` = `.__no soil__`) %>%
      dplyr::mutate(location_row = as.character(seq_len(dplyr::n()))) %>%
      dplyr::distinct(location_id, .keep_all = TRUE)
  }
} else {
  NULL
}

if (!is.null(gpts_map)) {
    cat("DEBUG: Unique Veg types in GeoJSON map:\n")
    print(table(gpts_map$Veg, useNA = "ifany"))

    # Check overlap
    ids_in_csv <- unique(df$location_id)
    ids_in_geo <- unique(gpts_map$location_id)
    
    cat("DEBUG: Total unique locations in CSV:", length(ids_in_csv), "\n")
    cat("DEBUG: Total unique locations in GeoJSON:", length(ids_in_geo), "\n")
    cat("DEBUG: Intersection of location IDs:", length(intersect(ids_in_csv, ids_in_geo)), "\n")
    
    # Check which veg types are in the intersection
    matched_ids <- intersect(ids_in_csv, ids_in_geo)
    matched_map <- gpts_map %>% filter(location_id %in% matched_ids)
    cat("DEBUG: Veg types in matched locations:\n")
    print(table(matched_map$Veg))
}

# Join df with gpts_map to add Veg column
if ("location_id" %in% names(df)) df$location_id <- as.character(df$location_id) else cat("Warning: input CSV missing 'location_id' column; skipping cast.\n")
join_possible <- (!is.null(gpts_map) && "location_id" %in% names(df) && "location_id" %in% names(gpts_map))
if (join_possible) {
  gpts_map$location_id <- as.character(gpts_map$location_id)
  df <- dplyr::left_join(df, gpts_map, by = "location_id", suffix = c("", ".geo"))
} else if (!is.null(gpts_map) && !"location_id" %in% names(df)) {
  cat("Skipping Veg join because 'location_id' not present in CSV\n")
} else {
  cat("Skipping Veg join because GeoJSON mapping is not available or incomplete\n")
}
if ("Veg.geo" %in% names(df)) {
  if (!"Veg" %in% names(df)) df$Veg <- NA_character_
  df$Veg <- ifelse(is.na(df$Veg) | df$Veg == "", df$Veg.geo, df$Veg)
  df$Veg.geo <- NULL
}
df <- normalize_no_soil_col(df)
if ("Veg" %in% names(df)) {
  df$Veg <- tolower(df$Veg)
  cat("DEBUG: Unique Veg types in df after join:\n")
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
df <- df %>% filter(year >= 1985 & year <= 2025)
cat("Data rows after year filtering:", nrow(df), "\n")

# Filter for years 1985-2025
cat("Filtering data for years 1985-2025...\n")
df <- df %>% filter(year >= 1985 & year <= 2025)
cat("Data rows after year filtering:", nrow(df), "\n")

# Check if MSAVI and NDVI are present
indices_to_check <- c("MSAVI", "NDVI")
missing_indices <- setdiff(indices_to_check, names(df))

if (length(missing_indices) > 0) {
  cat("Warning: The following indices are missing from the dataset:", paste(missing_indices, collapse = ", "), "\n")
  cat("Available indices:", paste(intersect(indices_to_check, names(df)), collapse = ", "), "\n")

  # If NDVI is missing, calculate it
  if ("NDVI" %in% missing_indices && all(c("nir", "red") %in% names(df))) {
    df$NDVI <- (df$nir - df$red) / (df$nir + df$red)
    cat("Calculated NDVI from nir and red bands\n")
  }

  # If MSAVI is missing, calculate it (MSAVI = (2*nir + 1 - sqrt((2*nir + 1)^2 - 8*(nir - red)))/2 )
  if ("MSAVI" %in% missing_indices && all(c("nir", "red") %in% names(df))) {
    df$MSAVI <- (2 * df$nir + 1 - sqrt((2 * df$nir + 1)^2 - 8 * (df$nir - df$red))) / 2
    cat("Calculated MSAVI from nir and red bands\n")
  }
}

# Calculate PPI using the new function
df <- add_ppi_columns(df, mean_lat_study_area)

# Remove PPI outliers using MAD (grouped by month to account for seasonality)
if ("PPI" %in% names(df)) {
  cat("Removing PPI outliers using MAD (k=3) per month...\n")
  df <- df %>%
    group_by(month) %>%
    mutate(PPI = remove_outliers_mad(PPI)) %>%
    ungroup()
}

# Filter for January and July data
january_data <- df %>% filter(month == 1)
july_data <- df %>% filter(month == 7)
september_data <- df %>% filter(month == 9)

cat("Number of January observations:", nrow(january_data), "\n")
cat("Number of July observations:", nrow(july_data), "\n")
cat("Number of September observations:", nrow(september_data), "\n")

# Compute dataset-level SNR for indices (per-index SNR across locations)
compute_global_index_snr <- function(df, indices, group_col = "location_id", eps = 1e-8) {
  out <- numeric(length(indices)); names(out) <- indices
  for (idx in indices) {
    if (!idx %in% names(df)) { out[idx] <- NA_real_; next }
    # per-location amplitude and noise
    locs <- unique(na.omit(as.character(df[[group_col]])))
    amps <- numeric(0); noises <- numeric(0)
    for (loc in locs) {
      sub <- df[df[[group_col]] == loc, , drop = FALSE]
      vals <- sub[[idx]]
      vals <- vals[is.finite(vals)]
      if (length(vals) < 2) next
      amp_loc <- diff(range(vals, na.rm = TRUE))
      noise_loc <- tryCatch({ mad(vals, na.rm = TRUE) }, error = function(e) NA_real_)
      if (is.finite(amp_loc)) amps <- c(amps, amp_loc)
      if (is.finite(noise_loc)) noises <- c(noises, noise_loc)
    }
    if (length(amps) == 0 || length(noises) == 0) { out[idx] <- NA_real_; next }
    sig <- median(amps, na.rm = TRUE)
    nois <- median(noises, na.rm = TRUE)
    snr_val <- sig / (nois + eps)
    out[idx] <- snr_val
  }
  out
}

# Calculate SNR for the three indices of interest and save results
indices_to_snr <- c("MSAVI", "NDVI", "PPI")
index_snr <- compute_global_index_snr(df, indices_to_snr, group_col = "location_id")
cat("Index SNR (raw):\n")
print(index_snr)
cat("Index SNR computed — will be written to the results folder if available when saving outputs\n")

# Function to bootstrap averages using location resampling
# Resamples unique location_ids, then takes all observations for those locations
bootstrap_location_means <- function(df, metrics = c("MSAVI", "NDVI", "PPI"), B = 1000) {
  # Ensure location_id is character
  if (!"location_id" %in% names(df)) return(rep(NA, length(metrics) * 2))
  df$location_id <- as.character(df$location_id)
  ids <- unique(df$location_id)
  
  if (length(ids) < 2) return(rep(NA, length(metrics) * 2))
  
  # Pre-split data by location_id to speed up access
  # Create a list where each element is a list of vectors for the metrics
  data_map <- lapply(ids, function(id) {
    sub <- df[df$location_id == id, ]
    lapply(metrics, function(m) if(m %in% names(sub)) sub[[m]] else numeric(0))
  })
  names(data_map) <- ids
  
  # Statistic function for boot
  boot_stat <- function(original_ids, indices) {
    selected_ids <- original_ids[indices]
    selected_data <- data_map[selected_ids]
    
    res <- vapply(seq_along(metrics), function(i) {
      # Compute mean per location, then mean of those means
      # This ensures the statistic matches the "mean of means" point estimate
      loc_means <- vapply(selected_data, function(sd) mean(sd[[i]], na.rm = TRUE), numeric(1))
      mean(loc_means, na.rm = TRUE)
    }, numeric(1))
    return(res)
  }
  
  # Run bootstrap
  b_res <- tryCatch(boot::boot(data = ids, statistic = boot_stat, R = B), error = function(e) NULL)
  if (is.null(b_res)) return(rep(NA, length(metrics) * 2))
  
  # Calculate CIs
  cis <- lapply(seq_along(metrics), function(i) {
    tryCatch({
      boot::boot.ci(b_res, type = "perc", index = i)$percent[4:5]
    }, error = function(e) c(NA, NA))
  })
  
  # Flatten result: MSAVI_lower, MSAVI_upper, NDVI_lower, ...
  unlist(cis)
}

# Function to bootstrap medians using location resampling
# Resamples unique location_ids, then takes all observations for those locations
bootstrap_location_medians <- function(df, metrics = c("PPI"), B = 1000) {
  # Ensure location_id is character
  if (!"location_id" %in% names(df)) return(rep(NA, length(metrics) * 2))
  df$location_id <- as.character(df$location_id)
  ids <- unique(df$location_id)
  
  if (length(ids) < 2) return(rep(NA, length(metrics) * 2))
  
  # Pre-split data by location_id to speed up access
  # Create a list where each element is a list of vectors for the metrics
  data_map <- lapply(ids, function(id) {
    sub <- df[df$location_id == id, ]
    lapply(metrics, function(m) if(m %in% names(sub)) sub[[m]] else numeric(0))
  })
  names(data_map) <- ids
  
  # Statistic function for boot
  boot_stat <- function(original_ids, indices) {
    selected_ids <- original_ids[indices]
    selected_data <- data_map[selected_ids]
    
    res <- vapply(seq_along(metrics), function(i) {
      # Compute median per location, then median of those medians
      loc_medians <- sapply(selected_data, function(sd) median(sd[[i]], na.rm = TRUE))
      median(loc_medians, na.rm = TRUE)
    }, numeric(1))
    return(res)
  }
  
  # Run bootstrap
  b_res <- tryCatch(boot::boot(data = ids, statistic = boot_stat, R = B), error = function(e) NULL)
  if (is.null(b_res)) return(rep(NA, length(metrics) * 2))
  
  # Calculate CIs
  cis <- lapply(seq_along(metrics), function(i) {
    tryCatch({
      boot::boot.ci(b_res, type = "perc", index = i)$percent[4:5]
    }, error = function(e) c(NA, NA))
  })
  
  # Flatten result: PPI_lower, PPI_upper, etc.
  unlist(cis)
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

# Calculate Global Averages (Aggregated 2020-2024) for January
cat("Calculating Global January Averages with Location Bootstrapping...\n")
jan_global_boot <- bootstrap_location_means(january_data, metrics = c("MSAVI", "NDVI", "PPI"), B = B)
jan_global_avg <- january_data %>%
  summarize(
    n_observations = n(),
    n_locations = n_distinct(location_id),
    avg_MSAVI = mean_of_means(MSAVI, location_id),
    avg_NDVI = mean_of_means(NDVI, location_id),
    avg_PPI = mean_of_means(PPI, location_id)
  ) %>%
  mutate(
    MSAVI_ci_lower = jan_global_boot[1], MSAVI_ci_upper = jan_global_boot[2],
    NDVI_ci_lower = jan_global_boot[3], NDVI_ci_upper = jan_global_boot[4],
    PPI_ci_lower = jan_global_boot[5], PPI_ci_upper = jan_global_boot[6]
  )

# Calculate Global Averages (Aggregated 2020-2024) for July
cat("Calculating Global July Averages with Location Bootstrapping...\n")
july_global_boot <- bootstrap_location_means(july_data, metrics = c("MSAVI", "NDVI", "PPI"), B = B)
july_global_avg <- july_data %>%
  summarize(
    n_observations = n(),
    n_locations = n_distinct(location_id),
    avg_MSAVI = mean_of_means(MSAVI, location_id),
    avg_NDVI = mean_of_means(NDVI, location_id),
    avg_PPI = mean_of_means(PPI, location_id)
  ) %>%
  mutate(
    MSAVI_ci_lower = july_global_boot[1], MSAVI_ci_upper = july_global_boot[2],
    NDVI_ci_lower = july_global_boot[3], NDVI_ci_upper = july_global_boot[4],
    PPI_ci_lower = july_global_boot[5], PPI_ci_upper = july_global_boot[6]
  )

# Calculate Global Averages (Aggregated 2020-2024) for September
cat("Calculating Global September Averages with Location Bootstrapping...\n")
sept_global_boot <- bootstrap_location_means(september_data, metrics = c("MSAVI", "NDVI", "PPI"), B = B)
sept_global_avg <- september_data %>%
  summarize(
    n_observations = n(),
    n_locations = n_distinct(location_id),
    avg_MSAVI = mean_of_means(MSAVI, location_id),
    avg_NDVI = mean_of_means(NDVI, location_id),
    avg_PPI = mean_of_means(PPI, location_id)
  ) %>%
  mutate(
    MSAVI_ci_lower = sept_global_boot[1], MSAVI_ci_upper = sept_global_boot[2],
    NDVI_ci_lower = sept_global_boot[3], NDVI_ci_upper = sept_global_boot[4],
    PPI_ci_lower = sept_global_boot[5], PPI_ci_upper = sept_global_boot[6]
  )

# --- Averages by Vegetation Type (Aggregated 2020-2024) ---

# Bootstrap for Veg type (January)
cat("Calculating Vegetation Type January Averages with Location Bootstrapping...\n")
veg_boot <- january_data %>%
  filter(!is.na(Veg)) %>%
  group_by(Veg) %>%
  do({
    res <- bootstrap_location_means(., metrics = c("MSAVI", "NDVI", "PPI"), B = B)
    data.frame(
      MSAVI_ci_lower = res[1], MSAVI_ci_upper = res[2],
      NDVI_ci_lower = res[3], NDVI_ci_upper = res[4],
      PPI_ci_lower = res[5], PPI_ci_upper = res[6]
    )
  }) %>%
  ungroup()

veg_type_averages <- january_data %>%
  filter(!is.na(Veg)) %>%
  group_by(Veg) %>%
  summarize(
    n_observations = n(),
    n_locations = n_distinct(location_id),
    avg_MSAVI = mean_of_means(MSAVI, location_id),
    avg_NDVI = mean_of_means(NDVI, location_id),
    avg_PPI = mean_of_means(PPI, location_id),
    .groups = "drop"
  )

veg_type_averages <- left_join(veg_type_averages, veg_boot, by = "Veg")

# Bootstrap for Veg type (July)
cat("Calculating Vegetation Type July Averages with Location Bootstrapping...\n")
veg_boot_july <- july_data %>%
  filter(!is.na(Veg)) %>%
  group_by(Veg) %>%
  do({
    res <- bootstrap_location_means(., metrics = c("MSAVI", "NDVI", "PPI"), B = B)
    data.frame(
      MSAVI_ci_lower = res[1], MSAVI_ci_upper = res[2],
      NDVI_ci_lower = res[3], NDVI_ci_upper = res[4],
      PPI_ci_lower = res[5], PPI_ci_upper = res[6]
    )
  }) %>%
  ungroup()

veg_type_averages_july <- july_data %>%
  filter(!is.na(Veg)) %>%
  group_by(Veg) %>%
  summarize(
    n_observations = n(),
    n_locations = n_distinct(location_id),
    avg_MSAVI = mean_of_means(MSAVI, location_id),
    avg_NDVI = mean_of_means(NDVI, location_id),
    avg_PPI = mean_of_means(PPI, location_id),
    .groups = "drop"
  )

veg_type_averages_july <- left_join(veg_type_averages_july, veg_boot_july, by = "Veg")

# Bootstrap for Veg type (September)
cat("Calculating Vegetation Type September Averages with Location Bootstrapping...\n")
veg_boot_sept <- september_data %>%
  filter(!is.na(Veg)) %>%
  group_by(Veg) %>%
  do({
    res <- bootstrap_location_means(., metrics = c("MSAVI", "NDVI", "PPI"), B = B)
    data.frame(
      MSAVI_ci_lower = res[1], MSAVI_ci_upper = res[2],
      NDVI_ci_lower = res[3], NDVI_ci_upper = res[4],
      PPI_ci_lower = res[5], PPI_ci_upper = res[6]
    )
  }) %>%
  ungroup()

veg_type_averages_sept <- september_data %>%
  filter(!is.na(Veg)) %>%
  group_by(Veg) %>%
  summarize(
    n_observations = n(),
    n_locations = n_distinct(location_id),
    avg_MSAVI = mean_of_means(MSAVI, location_id),
    avg_NDVI = mean_of_means(NDVI, location_id),
    avg_PPI = mean_of_means(PPI, location_id),
    .groups = "drop"
  )

veg_type_averages_sept <- left_join(veg_type_averages_sept, veg_boot_sept, by = "Veg")

# Print results
cat("\n=== GLOBAL JANUARY AVERAGES (1985-2025) ===\n")
print(jan_global_avg, width = Inf)

cat("\n=== GLOBAL JULY AVERAGES (1985-2025) ===\n")
print(july_global_avg, width = Inf)

cat("\n=== GLOBAL SEPTEMBER AVERAGES (1985-2025) ===\n")
print(sept_global_avg, width = Inf)

cat("\n=== JANUARY AVERAGES BY VEGETATION TYPE (1985-2025) ===\n")
print(veg_type_averages, width = Inf)

cat("\n=== JULY AVERAGES BY VEGETATION TYPE (1985-2025) ===\n")
print(veg_type_averages_july, width = Inf)

cat("\n=== SEPTEMBER AVERAGES BY VEGETATION TYPE (1985-2025) ===\n")
print(veg_type_averages_sept, width = Inf)

# --- Create Summary Table as requested ---
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
  create_summary_entry(idx, index_snr[idx], jan_global_avg, july_global_avg, sept_global_avg)
})

summary_table <- do.call(rbind, summary_rows)
colnames(summary_table) <- c("INDEX", "SNR", "January Average (All areas)", "July average (All areas)", "September average (All areas)")
summary_table <- as.data.frame(summary_table)

cat("\n=== SUMMARY TABLE ===\n")
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

# Plot mean June-September PPI from 2000 to 2025 across all locations with bootstrapping
summer_data <- df %>% filter(month %in% 6:9, PPI > 0)

if (nrow(summer_data) > 0) {
  cat("Generating summer trend plot with", nrow(summer_data), "observations...\n")
  
  # --- Normalization Step ---
  # Calculate Day of Year (DOY) if not present
  if (!"doy" %in% names(summer_data)) {
    summer_data$doy <- lubridate::yday(summer_data$date)
  }
  
  cat("Normalizing summer PPI for seasonal trend (June-Sept) using polynomial fit...\n")
  # Fit a 3rd degree polynomial to capture the seasonal curve across the summer months
  # This fits a curve through the data to model the typical seasonal pattern
  seasonal_model <- lm(PPI ~ poly(doy, 3), data = summer_data)
  
  # Predict the seasonal trend for each observation
  summer_data$seasonal_trend <- predict(seasonal_model, newdata = summer_data)
  
  # Calculate the global mean of the seasonal trend to maintain the overall magnitude
  global_seasonal_mean <- mean(summer_data$seasonal_trend, na.rm = TRUE)
  
  # Normalize: Subtract the seasonal deviation from the global mean
  # If a day typically has low PPI (trend < mean), we add the difference (correct upwards).
  # This "warps" the data to make the seasonal curve flat.
  summer_data$PPI_norm <- summer_data$PPI - (summer_data$seasonal_trend - global_seasonal_mean)
  
  # --- Bootstrapping with Normalized Data ---
  # Bootstrap means for each year using PPI_norm
  summer_yearly_boot <- summer_data %>%
    group_by(pheno_year) %>%
    do({
      # Note: We pass PPI_norm as the metric to bootstrap
      res <- bootstrap_location_means(., metrics = c("PPI_norm"), B = B)
      
      # Compute point estimate as mean of location-level means using PPI_norm
      loc_means <- sapply(unique(.$location_id), function(id) {
        sub <- .[.$location_id == id, ]
        mean(sub$PPI_norm, na.rm = TRUE)
      })
      mean_PPI <- mean(loc_means, na.rm = TRUE)
      
      data.frame(
        mean_PPI = mean_PPI,
        PPI_ci_lower = res[1],
        PPI_ci_upper = res[2]
      )
    }) %>%
    ungroup()
  
  p <- ggplot(summer_yearly_boot, aes(x = year, y = mean_PPI)) +
    geom_line(color = "blue") +
    geom_point(color = "red") +
    geom_ribbon(aes(ymin = PPI_ci_lower, ymax = PPI_ci_upper), alpha = 0.2, fill = "blue") +
    labs(title = "Mean Seasonally-Normalized June-Sept PPI (1985-2025)",
         subtitle = "Normalization: PPI - SeasonalTrend(doy) + Mean(SeasonalTrend)",
         x = "Year",
         y = "Mean Normalized PPI") +
    theme_minimal()
  
  print(p)
  
  # Save the plot
  ggsave(file.path(output_dir, "june_september_ppi_trend_normalized.png"), plot = p, width = 8, height = 6)
  cat("Saved plot to:", file.path(output_dir, "june_september_ppi_trend_normalized.png"), "\n")
  
  # --- Median Graph ---
  cat("Generating summer trend plot (MEDIAN) with", nrow(summer_data), "observations...\n")
  
  # Bootstrap medians for each year using PPI_norm
  summer_yearly_boot_median <- summer_data %>%
    group_by(pheno_year) %>%
    do({
      # Use bootstrap_location_medians for median CI
      res <- bootstrap_location_medians(., metrics = c("PPI_norm"), B = B)
      
      # Compute point estimate as median of location-level medians using PPI_norm
      loc_medians <- sapply(unique(.$location_id), function(id) {
        sub <- .[.$location_id == id, ]
        median(sub$PPI_norm, na.rm = TRUE)
      })
      median_PPI <- median(loc_medians, na.rm = TRUE)
      
      data.frame(
        median_PPI = median_PPI,
        PPI_ci_lower = res[1],
        PPI_ci_upper = res[2]
      )
    }) %>%
    ungroup()
  
  p_med <- ggplot(summer_yearly_boot_median, aes(x = year, y = median_PPI)) +
    geom_line(color = "darkgreen") +
    geom_point(color = "orange") +
    geom_ribbon(aes(ymin = PPI_ci_lower, ymax = PPI_ci_upper), alpha = 0.2, fill = "darkgreen") +
    labs(title = "Median Seasonally-Normalized June-Sept PPI (1985-2025)",
         subtitle = "Normalization: PPI - SeasonalTrend(doy) + Mean(SeasonalTrend)",
         x = "Year",
         y = "Median Normalized PPI") +
    theme_minimal()
  
  print(p_med)
  
  # Save the plot
  ggsave(file.path(output_dir, "june_september_ppi_trend_normalized_median.png"), plot = p_med, width = 8, height = 6)
  cat("Saved median plot to:", file.path(output_dir, "june_september_ppi_trend_normalized_median.png"), "\n")
  
} else {
  cat("Warning: No summer PPI data > 0 found. Skipping trend plot.\n")
}