#!/usr/bin/env Rscript
## R_extract_hls.R
# Aggregate per-location ESTARFM predictions, compute spectral indices, and
# produce the `phenology_results/hls_phenology_data.csv` input expected by
# `fit_veg_mixture_weighted.R`.

opts <- list()
opts$pred_dir <- Sys.getenv('ESTARFM_PRED_DIR', unset = 'C:/MAP/estarfm_results_gee')
opts$out_dir <- Sys.getenv('PHENO_OUT_DIR', unset = 'phenology_results')
dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

# Only accept predictor output files with the exact pattern:
# location_<id>_timeseries.csv (one file per location with all years)
files_all <- list.files(path = opts$pred_dir, pattern = '\\.(csv|CSV)$', recursive = TRUE, full.names = TRUE)
files <- files_all[grepl('location_\\d+_timeseries\\.csv$', basename(files_all), perl=TRUE)]
if(length(files) == 0) stop(sprintf('No predictor-style per-location timeseries CSVs found under %s. Run the predictor to generate files named like location_<id>_timeseries.csv', opts$pred_dir))

message(sprintf('Found %d location CSV files - reading...', length(files)))
library(data.table)
dt_list <- vector('list', length(files))
for(i in seq_along(files)){
  f <- files[i]
  dt <- tryCatch(fread(f), error = function(e) stop(sprintf('Failed reading %s: %s', f, e$message)))
  # normalize column names
  setnames(dt, tolower(names(dt)))
  # The predictor produces a per-location timeseries CSV with all years. Ensure the
  # file contains at least 2 observations (preferably more); fail loudly if
  # the file looks like a legacy single-observation export.
  if(nrow(dt) < 2) stop(sprintf('File %s appears to contain fewer than 2 observations. The extractor requires per-location timeseries files (multiple rows).', f))
  dt_list[[i]] <- dt
}
dt <- rbindlist(dt_list, fill = TRUE)

# Ensure date column exists and is Date
if(!'prediction_date' %in% names(dt) && 'date' %in% names(dt)) dt[, prediction_date := date]
if(!'prediction_date' %in% names(dt)) stop('Input CSVs must include a prediction_date (or date) column')
dt[, date := as.IDate(prediction_date, format = '%Y-%m-%d')]
if(any(is.na(dt$date))) stop('Some prediction_date values could not be parsed as YYYY-MM-DD')
dt[, year := year(date)]
dt[, doy := as.integer(strftime(date, '%j'))]

# Create a standardized location_id matching the fit script's expectation:
# "L_<lon>_<lat>" rounded to 4 decimals (so geojson mapping will match)
if(!('lat' %in% names(dt) && 'lon' %in% names(dt))) stop('Per-location CSVs must include lat and lon columns')
dt[, location_id := sprintf('L_%0.4f_%0.4f', round(lon,4), round(lat,4))]

# Ensure expected band columns exist; fail loudly if not
expected_bands <- c('band_blue','band_green','band_red','band_nir','band_swir1','band_swir2')
missing_bands <- setdiff(expected_bands, names(dt))
if(length(missing_bands) > 0) stop(sprintf('Missing required band columns: %s', paste(missing_bands, collapse = ', ')))

# Convert to numeric and coerce NA if non-numeric
for(b in expected_bands) set(dt, j = b, value = as.numeric(dt[[b]]))

# Rename band columns to remove 'band_' prefix for consistency
setnames(dt, 
         old = c('band_blue','band_green','band_red','band_nir','band_swir1','band_swir2'),
         new = c('blue','green','red','nir','swir1','swir2'))

# Compute spectral indices (vectorized). Use a small epsilon to avoid division by zero
eps <- 1e-9
with(dt, {
  DVI <- nir - red
  OSAVI <- (nir - red) / (nir + red + 0.16)
  MCARI <- ((red - green) - 0.2*(red - blue)) * (red / (green + eps))
  CRI <- (1/(green + eps)) - (1/(red + eps))
  PRI <- (green - red) / (green + red + eps)
  NIRv <- nir * ((nir - red) / (nir + red + eps))
  PSRI <- (red - blue) / (nir + eps)
  NBR <- (nir - swir2) / (nir + swir2 + eps)
  TCW <- (swir1 - swir2) / (swir1 + swir2 + eps)  # Tasseled Cap Wetness approximation
  TCG <- (green - red) / (green + red + eps)  # Tasseled Cap Greenness approximation
  MNDWI <- (green - swir1) / (green + swir1 + eps)
  DUSTI <- (red - blue) / (red + blue + eps)
  # Attach computed indices
  dt[, DVI := DVI]
  dt[, OSAVI := OSAVI]
  dt[, MCARI := MCARI]
  dt[, CRI := CRI]
  dt[, PRI := PRI]
  dt[, NIRv := NIRv]
  dt[, PSRI := PSRI]
  dt[, NBR := NBR]
  dt[, TCW := TCW]
  dt[, TCG := TCG]
  dt[, MNDWI := MNDWI]
  dt[, DUSTI := DUSTI]
})

# Apply transformations to improve linearity
# NIRv: kNDVI transform with α ≈ 1.2-1.5
dt[, NIRv := NIRv * 1.3]  # Using midpoint of recommended range

# Validate we produced all OPTIMAL_INDICES expected by the fitter
OPTIMAL_INDICES <- c(
  'DVI','OSAVI','MCARI','CRI','PRI','NIRv','PSRI','NBR','TCW','TCG','MNDWI'
)
# DUSTI is intentionally not in OPTIMAL_INDICES, it is for filtering only
missing_idx <- setdiff(OPTIMAL_INDICES, names(dt))
if(length(missing_idx) > 0) stop(sprintf('Extractor did not produce required indices: %s', paste(missing_idx, collapse = ', ')))

# --- DUSTI filtering: remove observations with high DUSTI values ---
if("DUSTI" %in% names(dt)) {
  # Define a threshold for DUSTI. This might need tuning.
  DUSTI_THRESHOLD <- 0.1 
  n_before <- nrow(dt)
  dt <- dt[dt$DUSTI <= DUSTI_THRESHOLD, ]
  n_after <- nrow(dt)
  if (n_before > n_after) {
    cat(sprintf("Applied DUSTI filter: removed %d observations with DUSTI > %.2f\n", n_before - n_after, DUSTI_THRESHOLD))
  }
}

# --- Spline-based outlier detection: remove extreme values using moving splines ---
# Apply spline outlier detection to all spectral indices per location-year combination
if(all(c("DVI", "date", "doy") %in% names(dt))) {
  # Define the detect_spline_outliers function
  detect_spline_outliers <- function(values, doys, window_size = 45, outlier_threshold = 4.0) {
    if(length(values) < 10 || length(unique(doys)) < 5) {
      return(rep(FALSE, length(values)))  # Not enough data for reliable spline fitting
    }
    
    # Sort by DOY for proper spline fitting
    ord <- order(doys)
    vals_sorted <- values[ord]
    doys_sorted <- doys[ord]
    
    # Remove any remaining NAs for spline fitting
    valid_idx <- is.finite(vals_sorted) & is.finite(doys_sorted)
    vals_clean <- vals_sorted[valid_idx]
    doys_clean <- doys_sorted[valid_idx]
    
    if(length(vals_clean) < 5) {
      return(rep(FALSE, length(values)))
    }
    
    tryCatch({
      # Fit a smooth spline to the time series
      spline_fit <- stats::smooth.spline(doys_clean, vals_clean, df = min(length(vals_clean)/3, 10))
      
      # Predict values for all DOYs
      predicted <- stats::predict(spline_fit, doys_clean)$y
      
      # Calculate residuals
      residuals <- vals_clean - predicted
      
      # Calculate robust scale estimate (MAD)
      mad_residuals <- stats::mad(residuals, na.rm = TRUE, constant = 1.4826)
      
      if(!is.finite(mad_residuals) || mad_residuals <= 0) {
        return(rep(FALSE, length(values)))
      }
      
      # Identify outliers
      outlier_flags_clean <- abs(residuals) > outlier_threshold * mad_residuals
      
      # Map back to original indices
      outlier_flags <- rep(FALSE, length(vals_sorted))
      outlier_flags[valid_idx] <- outlier_flags_clean
      
      # Map back to original order
      result <- rep(FALSE, length(values))
      result[ord] <- outlier_flags
      
      return(result)
    }, error = function(e) {
      # If spline fitting fails, return no outliers
      cat(sprintf("Warning: Spline fitting failed for outlier detection: %s\n", e$message))
      return(rep(FALSE, length(values)))
    })
  }
  
  # Get list of spectral indices to check for outliers
  spectral_indices <- intersect(OPTIMAL_INDICES, names(dt))
  
  # Apply outlier detection to each spectral index per location-year
  dt[, outlier_mask := FALSE]  # Initialize outlier mask
  
  for(idx in spectral_indices) {
    if(idx %in% c("location_id", "year", "date", "doy")) next
    
    dt[, {
      if(.N >= 10) {  # Need minimum observations for spline fitting
        values <- get(idx)
        doys <- doy
        
        # Detect outliers for this index
        idx_outliers <- detect_spline_outliers(values, doys)
        
        # Update overall outlier mask
        outlier_mask <<- outlier_mask | idx_outliers
      }
    }, by = .(location_id, year)]
  }
  
  # Remove detected outliers from the dataset
  n_total_outliers <- sum(dt$outlier_mask, na.rm = TRUE)
  if(n_total_outliers > 0) {
    cat(sprintf("Spline-based outlier detection removed %d observations\n", n_total_outliers))
    dt <- dt[!outlier_mask, ]
  }
  
  # Clean up temporary column
  dt[, outlier_mask := NULL]
}

# --- Dynamic baseline computation and subtraction ---
# Apply baseline correction to normalize spectral indices per location-year
if(all(c("DVI", "date", "doy") %in% names(dt))) {
  # Process each location-year combination
  loc_years <- unique(dt[, .(location_id, year)])
  
  for(i in seq_len(nrow(loc_years))) {
    loc <- loc_years$location_id[i]
    yr <- loc_years$year[i]
    
    # Subset data for this location-year
    idx_rows <- dt$location_id == loc & dt$year == yr
    loc_yr_data <- dt[idx_rows, ]
    
    if(nrow(loc_yr_data) < 30) {
      cat(sprintf("Warning: Insufficient data for baseline computation in %s_%d: only %d observations (minimum 30 required)\n", 
                  loc, yr, nrow(loc_yr_data)))
      next
    }
    
    # Sort by day of year
    loc_yr_data <- loc_yr_data[order(loc_yr_data$doy), ]
    
    # Find the lowest continuous 30-day period
    window_size <- 30
    n_windows <- nrow(loc_yr_data) - window_size + 1
    
    if(n_windows > 0) {
      window_medians <- rep(NA_real_, n_windows)
      
      for(w in 1:n_windows) {
        window_data <- loc_yr_data[w:(w + window_size - 1), ]
        window_dvi <- window_data$DVI[is.finite(window_data$DVI)]
        if(length(window_dvi) >= 10) {  # Require at least 10 valid observations
          window_medians[w] <- median(window_dvi, na.rm = TRUE)
        }
      }
      
      # Find window with lowest median MSAVI
      if(any(is.finite(window_medians))) {
        best_window_idx <- which.min(window_medians)
        baseline_period <- loc_yr_data[best_window_idx:(best_window_idx + window_size - 1), ]
        
        # Calculate baseline as median of DVI in this 30-day period
        baseline_median <- median(baseline_period$DVI[is.finite(baseline_period$DVI)], na.rm = TRUE)
        
        if(is.finite(baseline_median)) {
          # Apply baseline subtraction to all spectral indices
          for(idx in OPTIMAL_INDICES) {
            if(idx %in% names(dt)) {
              sub_vals <- dt[[idx]][idx_rows]
              # Use median from the 20-day baseline period
              baseline_vals <- baseline_period[[idx]]
              base_med <- suppressWarnings(median(baseline_vals[is.finite(baseline_vals)], na.rm = TRUE))
              
              if(is.finite(base_med)) {
                dt[[idx]][idx_rows] <- sub_vals - base_med
              } else {
                cat(sprintf("Warning: Could not compute baseline for index %s in %s_%d\n", idx, loc, yr))
              }
            }
          }
          
          # Diagnostic output
          cat(sprintf("Baseline for %s_%d: 30-day window starting DOY %d, baseline=%.4f\n", 
                      loc, yr, baseline_period$doy[1], baseline_median))
        } else {
          cat(sprintf("Warning: Baseline computation failed for %s_%d: computed baseline is not finite\n", loc, yr))
        }
      } else {
        cat(sprintf("Warning: Failed to find valid 30-day baseline window for %s_%d: no windows with sufficient observations\n", loc, yr))
      }
    } else {
      cat(sprintf("Warning: Failed to find valid 30-day baseline window for %s_%d: insufficient finite DVI values\n", loc, yr))
    }
  }
  
  cat("Baseline normalization completed.\n")
}

# Minimal metadata expected by the fit script
dt[, imagery_lon := lon]
dt[, imagery_lat := lat]
dt[, target_lon := lon]
dt[, target_lat := lat]

# Reorder columns: metadata then indices
meta_cols <- c('location_id','lat','lon','imagery_lat','imagery_lon','target_lat','target_lon','date','year','doy')
idx_cols <- c(OPTIMAL_INDICES, 'DUSTI')
other_cols <- setdiff(names(dt), c(meta_cols, idx_cols))
out_cols <- c(meta_cols, idx_cols, other_cols)
out_cols <- out_cols[out_cols %in% names(dt)]
out_dt <- dt[, ..out_cols]

out_file <- file.path(opts$out_dir, 'hls_phenology_data.csv')
# Always overwrite the phenology data file - no existence checks
fwrite(out_dt, out_file)
message(sprintf('Wrote combined phenology CSV: %s (%d rows)', out_file, nrow(out_dt)))

invisible(NULL)
