#!/usr/bin/env Rscript
## R_extract_hls.R
# Aggregate per-location ESTARFM predictions, compute spectral indices, and
# produce the `phenology_results/hls_phenology_data.csv` input expected by
# `fit_veg_mixture_weighted.R`.

opts <- list()
# Path to a single combined CSV containing all incremental predictions and indices.
# If present, this file will be used and the per-location directory scanning will be skipped.
opts$input_file <- Sys.getenv('PHENO_INPUT_FILE', unset = 'C:/MAP/incremental_results.csv')
opts$pred_dir <- Sys.getenv('ESTARFM_PRED_DIR', unset = 'C:/MAP/estarfm_results_gee')
opts$out_dir <- Sys.getenv('PHENO_OUT_DIR', unset = 'phenology_results')
dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)


library(data.table)

# If a single combined input CSV exists, use it directly; otherwise scan for a
# directory containing per-location timeseries CSVs as before.
if (file.exists(opts$input_file)) {
  message(sprintf('Found combined input CSV: %s - reading...', opts$input_file))
  dt <- tryCatch(fread(opts$input_file), error = function(e) stop(sprintf('Failed reading input file %s: %s', opts$input_file, e$message)))
  setnames(dt, tolower(names(dt)))
  # For combined file, at least 2 rows is required
  if (nrow(dt) < 2) stop(sprintf('Input CSV %s contains fewer than 2 rows; need at least two observations', opts$input_file))
} else {
  stop(sprintf('No combined input CSV found at %s. This extractor now requires a single combined incremental CSV; set the environment variable PHENO_INPUT_FILE to point to it (default: %s).', opts$input_file, opts$input_file))
}

# Ensure date column exists and is Date
if(!'prediction_date' %in% names(dt) && 'date' %in% names(dt)) dt[, prediction_date := date]
if(!'prediction_date' %in% names(dt)) stop('Input CSVs must include a prediction_date (or date) column')
dt[, date := as.IDate(prediction_date, format = '%Y-%m-%d')]
if(any(is.na(dt$date))) stop('Some prediction_date values could not be parsed as YYYY-MM-DD')
dt[, year := year(date)]
dt[, doy := as.integer(strftime(date, '%j'))]

# Ensure expected band columns exist; fail loudly if not
expected_bands <- c('band_blue','band_green','band_red','band_nir','band_swir1','band_swir2')
missing_bands <- setdiff(expected_bands, names(dt))
if(length(missing_bands) > 0) stop(sprintf('Missing required band columns: %s', paste(missing_bands, collapse = ', ')))

# Convert to numeric and coerce NA if non-numeric
for(b in expected_bands) set(dt, j = b, value = as.numeric(dt[[b]]))

#

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
  #CRI <- (1/(green + eps)) - (1/(red + eps))
  PRI <- (green - red) / (green + red + eps)
  NIRv <- nir * ((nir - red) / (nir + red + eps))
  PSRI <- (red - blue) / (nir + eps)
  NBR <- (nir - swir2) / (nir + swir2 + eps)
  TCW <- (swir1 - swir2) / (swir1 + swir2 + eps)  # Tasseled Cap Wetness approximation
  TCG <- (green - red) / (green + red + eps)  # Tasseled Cap Greenness approximation
  MNDWI <- (green - swir1) / (green + swir1 + eps)
  DUSTI <- (red - blue) / (red + blue + eps)
  
  # --- New Indices (Linearity & Soil/Moisture) ---
  # NDVI: Normalized Difference Vegetation Index
  NDVI <- (nir - red) / (nir + red + eps)
  
  # MSAVI2: Modified Soil Adjusted Vegetation Index 2
  # Minimizes soil background influence
  MSAVI2 <- (2 * nir + 1 - sqrt(pmax(0, (2 * nir + 1)^2 - 8 * (nir - red)))) / 2
  
  # NDMI: Normalized Difference Moisture Index
  NDMI <- (nir - swir1) / (nir + swir1 + eps)
  
  # TCB: Tasseled Cap Brightness (Landsat 8 coefficients)
  # 0.3029*Blue + 0.2786*Green + 0.4733*Red + 0.5599*NIR + 0.508*SWIR1 + 0.1872*SWIR2
  TCB <- 0.3029 * blue + 0.2786 * green + 0.4733 * red + 0.5599 * nir + 0.508 * swir1 + 0.1872 * swir2
  
  # GVI: Green Vegetation Index (Tasseled Cap Greenness - Linear Combination)
  # -0.2941*Blue - 0.243*Green - 0.5424*Red + 0.7276*NIR + 0.0713*SWIR1 - 0.1608*SWIR2
  GVI <- -0.2941 * blue - 0.243 * green - 0.5424 * red + 0.7276 * nir + 0.0713 * swir1 - 0.1608 * swir2

  # Attach computed indices
  dt[, DVI := DVI]
  dt[, OSAVI := OSAVI]
  dt[, MCARI := MCARI]
  #dt[, CRI := CRI]
  dt[, PRI := PRI]
  dt[, NIRv := NIRv]
  dt[, PSRI := PSRI]
  dt[, NBR := NBR]
  dt[, TCW := TCW]
  dt[, TCG := TCG]
  dt[, MNDWI := MNDWI]
  dt[, DUSTI := DUSTI]
  dt[, NDVI := NDVI]
  dt[, MSAVI2 := MSAVI2]
  dt[, NDMI := NDMI]
  dt[, TCB := TCB]
  dt[, GVI := GVI]
})

# Apply transformations to improve linearity
# NIRv: kNDVI transform with α ≈ 1.2-1.5
dt[, NIRv := NIRv * 1.3]  # Using midpoint of recommended range

# Validate we produced all OPTIMAL_INDICES expected by the fitter
OPTIMAL_INDICES <- c(
  'DVI','OSAVI','MCARI','PRI','NIRv','PSRI','NBR','TCW','TCG','MNDWI',
  'NDVI', 'MSAVI2', 'NDMI', 'TCB', 'GVI'
)
# DUSTI is intentionally not in OPTIMAL_INDICES, it is for filtering only
missing_idx <- setdiff(OPTIMAL_INDICES, names(dt))
if(length(missing_idx) > 0) stop(sprintf('Extractor did not produce required indices: %s', paste(missing_idx, collapse = ', ')))

# --- DUSTI filtering: remove observations with high DUSTI values ---
if("DUSTI" %in% names(dt)) {
  # Define a threshold for DUSTI. This might need tuning.
  DUSTI_THRESHOLD <- 0.5 
  n_before <- nrow(dt)
  dt <- dt[dt$DUSTI <= DUSTI_THRESHOLD, ]
  n_after <- nrow(dt)
  if (n_before > n_after) {
    cat(sprintf("Applied DUSTI filter: removed %d observations with DUSTI > %.2f\n", n_before - n_after, DUSTI_THRESHOLD))
  }
}

# --- Spline-based outlier detection (robust 2-pass with MAD fallback) ---
# Unified behavior: use the same iterative 2-pass spline + MAD logic as the training pipeline.
# Requires at least 10 observations per group and >=5 finite points for fitting.
if (all(c("DVI", "date", "doy") %in% names(dt))) {
  if (!exists("OUTLIER_MAD_THRESHOLD")) OUTLIER_MAD_THRESHOLD <- 3.5

  detect_spline_outliers <- function(values, doys, mad_thresh = OUTLIER_MAD_THRESHOLD) {
    # Quick checks for sufficient data
    if (length(values) < 5 || length(unique(doys)) < 3) return(rep(FALSE, length(values)))

    ord <- order(doys)
    vals_sorted <- values[ord]
    doys_sorted <- doys[ord]

    # Remove non-finite pairs
    valid_idx <- is.finite(vals_sorted) & is.finite(doys_sorted)
    if (sum(valid_idx) < 5) return(rep(FALSE, length(values)))

    x <- doys_sorted[valid_idx]
    y <- vals_sorted[valid_idx]

    tryCatch({
      # Pass 1: initial fit
      fit1 <- stats::smooth.spline(x, y, df = min(5, length(x)/2))
      pred1 <- predict(fit1, x)$y
      res1 <- y - pred1
      mad1 <- stats::mad(res1, na.rm = TRUE)
      if (!is.finite(mad1) || mad1 <= 1e-6) stop("Invalid MAD in Pass 1")

      # Exclude gross outliers (1.5x scaled threshold) for robust refit
      keep_mask <- abs(res1 - stats::median(res1, na.rm = TRUE)) <= (mad_thresh * 1.5 * mad1)

      # Pass 2: refit if enough points left, otherwise fallback to pass1 predictions
      if (sum(keep_mask) >= 5) {
        fit2 <- stats::smooth.spline(x[keep_mask], y[keep_mask], df = min(5, sum(keep_mask)/2))
        pred_final <- predict(fit2, x)$y
      } else {
        pred_final <- pred1
      }

      # Final residual-based outlier detection using MAD
      residuals <- y - pred_final
      med_res <- stats::median(residuals, na.rm = TRUE)
      mad_res <- stats::mad(residuals, na.rm = TRUE)
      if (!is.finite(mad_res) || mad_res <= 0) stop("Invalid Final MAD")

      this_mask <- rep(FALSE, length(vals_sorted))
      this_mask[valid_idx] <- abs(residuals - med_res) > mad_thresh * mad_res

      result <- rep(FALSE, length(values))
      result[ord] <- this_mask
      return(result)
    }, error = function(e) {
      # Fallback: use simple MAD-based outlier detection on raw values
      med <- stats::median(vals_sorted, na.rm = TRUE)
      m <- stats::mad(vals_sorted, na.rm = TRUE)
      if (is.finite(m) && m > 0) {
        out_flags <- abs(vals_sorted - med) > mad_thresh * m
        res <- rep(FALSE, length(values))
        res[ord] <- out_flags
        return(res)
      }
      return(rep(FALSE, length(values)))
    })
  }

  # Spectral indices to check (intersection with available columns)
  spectral_indices <- intersect(OPTIMAL_INDICES, names(dt))

  # Initialize mask and run detection per (location_id, year)
  dt[, outlier_mask := FALSE]
  for (idx in spectral_indices) {
    if (idx %in% c("location_id", "year", "date", "doy")) next
    dt[, {
      if (.N >= 10) {
        values <- get(idx)
        doys <- doy
        idx_outliers <- detect_spline_outliers(values, doys)
        outlier_mask <<- outlier_mask | idx_outliers
      }
    }, by = .(location_id, year)]
  }

  # Remove detected outliers
  n_total_outliers <- sum(dt$outlier_mask, na.rm = TRUE)
  if (n_total_outliers > 0) {
    cat(sprintf("Spline-based outlier detection removed %d observations\n", n_total_outliers))
    dt <- dt[!outlier_mask, ]
  }
  dt[, outlier_mask := NULL]
} else {
  # Required columns for spline-based detection missing; skipping
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
