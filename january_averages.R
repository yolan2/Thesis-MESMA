library(data.table)
library(ggplot2)
library(dplyr)
library(lubridate)
# Robust regression (Huber weighting)
library(MASS)



# Ensure PPI helpers are loaded early so functions like calculate_solar_zenith are available
if (!exists("calculate_solar_zenith") && file.exists("ppi_helpers.R")) source("ppi_helpers.R")
# Optional visualization helpers (provide shading for excluded years)
if (file.exists("mesma_helpers.R")) {
  source("mesma_helpers.R")
  # Ensure reproducible sampling when running this script standalone
  set_mesma_seed()
}

# Provide a safe fallback logger so this script can run standalone
if (!exists("write_debug", mode = "function")) {
  write_debug <- function(msg) {
    tryCatch({ cat(paste0("[DEBUG] ", as.character(msg), "\n")) }, error = function(e) {})
  }
}

# Raw spectral bands (optional - included if present)
RAW_BANDS <- c("blue", "green", "red", "nir", "swir1", "swir2")
# Indices we consider for FVC calibration & diagnostics
INDICES_OF_INTEREST <- c("MSAVI", "NDVI", "PPI", "OSAVI", "NIRv", "NBR", "TCW", "NDMI", "TCB", "GVI", "EVI")

# Detrending configuration: which months constitute "summer" for seasonal model fitting
# and the minimum number of finite samples required per index to build a seasonal model.
SUMMER_DETREND_MONTHS <- c(6,7,8)  # June - August (use median across these months)
MIN_SEASONAL_SAMPLES <- 50

# Normalize known raw band column names (case and prefix variants) to canonical lower-case names
normalize_band_names <- function(df, bands = RAW_BANDS) {
  if (is.null(df) || nrow(df) == 0) return(df)
  current_names <- names(df)
  for (b in bands) {
    # Candidate variations to match
    candidates <- c(b, toupper(b), tools::toTitleCase(b), paste0('band_', b), toupper(paste0('band_', b)), paste0('Band_', b))
    for (cand in candidates) {
      if (cand %in% current_names && !(b %in% current_names)) {
        names(df)[names(df) == cand] <- b
        cat(sprintf("[NOTICE] Normalized band column '%s' -> '%s'\n", cand, b))
        current_names <- names(df)
        break
      }
    }
  }
  df
}

# Compute a set of spectral indices from raw bands (compatible with R_extract_hls.R's formulae)
compute_indices_from_bands <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  eps <- 1e-9
  has_bands <- intersect(RAW_BANDS, names(df))
  if (length(has_bands) == 0) return(df)

  # Safely compute each index only if its required bands are present
  if (all(c('nir','red') %in% names(df))) df$DVI <- as.numeric(df$nir) - as.numeric(df$red)
  if (all(c('nir','red') %in% names(df))) df$OSAVI <- (as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + as.numeric(df$red) + 0.16)
  if (all(c('red','green','blue') %in% names(df))) df$MCARI <- ((as.numeric(df$red) - as.numeric(df$green)) - 0.2*(as.numeric(df$red) - as.numeric(df$blue))) * (as.numeric(df$red) / (as.numeric(df$green) + eps))
  if (all(c('green','red') %in% names(df))) df$PRI <- (as.numeric(df$green) - as.numeric(df$red)) / (as.numeric(df$green) + as.numeric(df$red) + eps)
  if (all(c('nir','red') %in% names(df))) df$NIRv <- as.numeric(df$nir) * ((as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + as.numeric(df$red) + eps))
  if (all(c('red','blue','nir') %in% names(df))) df$PSRI <- (as.numeric(df$red) - as.numeric(df$blue)) / (as.numeric(df$nir) + eps)
  if (all(c('nir','swir2') %in% names(df))) df$NBR <- (as.numeric(df$nir) - as.numeric(df$swir2)) / (as.numeric(df$nir) + as.numeric(df$swir2) + eps)
  if (all(c('swir1','swir2') %in% names(df))) df$TCW <- (as.numeric(df$swir1) - as.numeric(df$swir2)) / (as.numeric(df$swir1) + as.numeric(df$swir2) + eps)
  if (all(c('green','red','nir','swir1','swir2','blue') %in% names(df))) {
    df$TCB <- 0.3029 * as.numeric(df$blue) + 0.2786 * as.numeric(df$green) + 0.4733 * as.numeric(df$red) + 0.5599 * as.numeric(df$nir) + 0.508 * as.numeric(df$swir1) + 0.1872 * as.numeric(df$swir2)
    df$GVI <- -0.2941 * as.numeric(df$blue) - 0.243 * as.numeric(df$green) - 0.5424 * as.numeric(df$red) + 0.7276 * as.numeric(df$nir) + 0.0713 * as.numeric(df$swir1) - 0.1608 * as.numeric(df$swir2)
  }
  if (all(c('nir','red') %in% names(df))) df$NDVI <- (as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + as.numeric(df$red) + eps)
  if (all(c('nir','red') %in% names(df))) df$MSAVI2 <- (2 * as.numeric(df$nir) + 1 - sqrt(pmax(0, (2 * as.numeric(df$nir) + 1)^2 - 8 * (as.numeric(df$nir) - as.numeric(df$red))))) / 2
  if (all(c('nir','red') %in% names(df))) df$MSAVI <- (2 * as.numeric(df$nir) + 1 - sqrt(pmax(0, (2 * as.numeric(df$nir) + 1)^2 - 8 * (as.numeric(df$nir) - as.numeric(df$red))))) / 2
  if (all(c('nir','swir1') %in% names(df))) df$NDMI <- (as.numeric(df$nir) - as.numeric(df$swir1)) / (as.numeric(df$nir) + as.numeric(df$swir1) + eps)

  # Apply kNDVI-like tweak used in extractor (NIRv * 1.3)
  if ('NIRv' %in% names(df)) df$NIRv <- df$NIRv * 1.3

  df
}

# Outlier configuration defined early so helper can be used during initial filtering
ENABLE_OUTLIER_REMOVAL <- FALSE
OUTLIER_MAD_THRESHOLD <- 3.5

# Remove large outliers robustly per (location_id, pheno_year) where possible, otherwise per-location.
# Uses spline-based outlier detection for groups with sufficient data, otherwise falls back to MAD.
remove_large_outliers <- function(df, candidates = NULL, mad_thresh = OUTLIER_MAD_THRESHOLD) {
  if (!exists("OUTLIER_SPLINE_MAX_DF", inherits = TRUE)) OUTLIER_SPLINE_MAX_DF <- 10L
  if (!isTRUE(ENABLE_OUTLIER_REMOVAL)) return(df)
  if (is.null(candidates)) {
    # Use indices we know are meaningful: INDICES_OF_INTEREST + RAW_BANDS, if present
    candidates <- intersect(unique(c(INDICES_OF_INTEREST, RAW_BANDS)), names(df))
  } else {
    candidates <- intersect(candidates, names(df))
  }
  if (length(candidates) == 0) {
    cat("[OUTLIER] No candidate indices found for outlier detection; skipping\n")
    return(df)
  }
  if (!"location_id" %in% names(df)) {
    cat("[OUTLIER] 'location_id' missing from data; skipping outlier removal\n")
    return(df)
  }

  # Assign pheno_year if possible for better grouping
  if (!"pheno_year" %in% names(df) && "date" %in% names(df)) {
    df$pheno_year <- assign_pheno_year(df$date)
  }

  grp <- interaction(df$location_id, ifelse(is.na(df$pheno_year), "NA", as.character(df$pheno_year)), drop = TRUE)
  groups <- split(seq_len(nrow(df)), grp)
  removed_idx <- logical(nrow(df))
  n_groups <- length(groups)

  for (g in seq_along(groups)) {
    rows <- groups[[g]]
    sub <- df[rows, , drop = FALSE]
    if (length(rows) < 3) next  # not enough data
    out_mask <- rep(FALSE, nrow(sub))

    # Check if we have date for spline
    has_date <- "date" %in% names(sub) && any(!is.na(sub$date))
    use_spline <- has_date && length(rows) >= 10  # Use spline if enough data and date available

    if (use_spline) {
      # Compute DOY
      sub$doy <- as.numeric(format(sub$date, "%j"))
      for (col in candidates) {
        if (!is.numeric(sub[[col]])) next
        colv <- sub[[col]]
        finite_idx <- is.finite(colv) & is.finite(sub$doy)
        if (sum(finite_idx) < 5) next  # not enough for spline
        tryCatch({
          x <- sub$doy[finite_idx]
          y <- colv[finite_idx]
          n_unique <- length(unique(x))
          fit1 <- stats::smooth.spline(x, y, df = min(OUTLIER_SPLINE_MAX_DF, length(x)/2, n_unique - 1))
          pred1 <- predict(fit1, x)$y
          res1 <- y - pred1
          mad1 <- stats::mad(res1, na.rm = TRUE)
          if (!is.finite(mad1) || mad1 <= 1e-6) stop("Invalid MAD")
          keep_mask <- abs(res1 - stats::median(res1, na.rm = TRUE)) <= (mad_thresh * 1.5 * mad1)
          if (sum(keep_mask) >= 5) {
            n_unique2 <- length(unique(x[keep_mask]))
            fit2 <- stats::smooth.spline(x[keep_mask], y[keep_mask], df = min(OUTLIER_SPLINE_MAX_DF, sum(keep_mask)/2, n_unique2 - 1))
            pred_final <- predict(fit2, x)$y
          } else {
            pred_final <- pred1
          }
          residuals <- y - pred_final
          med_res <- stats::median(residuals, na.rm = TRUE)
          mad_res <- stats::mad(residuals, na.rm = TRUE)
          if (!is.finite(mad_res) || mad_res <= 0) stop("Invalid Final MAD")
          this_mask <- rep(FALSE, length(colv))
          this_mask[finite_idx] <- abs(residuals - med_res) > mad_thresh * mad_res
          out_mask <- out_mask | this_mask
        }, error = function(e) {
          med <- stats::median(colv, na.rm = TRUE)
          m <- stats::mad(colv, na.rm = TRUE)
          if (is.finite(m) && m > 0) {
            this_mask <- is.finite(colv) & (abs(colv - med) > mad_thresh * m)
            out_mask <<- out_mask | this_mask
          }
        })
      }
    } else {
      for (col in candidates) {
        if (!is.numeric(sub[[col]])) next
        colv <- sub[[col]]
        if (all(is.na(colv))) next
        med <- stats::median(colv, na.rm = TRUE)
        m <- stats::mad(colv, na.rm = TRUE)
        if (!is.finite(m) || m <= 0) next
        this_mask <- is.finite(colv) & (abs(colv - med) > mad_thresh * m)
        this_mask[is.na(this_mask)] <- FALSE
        out_mask <- out_mask | this_mask
      }
    }
    if (any(out_mask, na.rm = TRUE)) removed_idx[rows[which(out_mask)]] <- TRUE
  }
  if (any(removed_idx, na.rm = TRUE)) {
    n_removed <- sum(removed_idx, na.rm = TRUE)
    cat(sprintf("[OUTLIER] Removed %d observations across %d groups\n", n_removed, n_groups))
    df <- df[!removed_idx, , drop = FALSE]
  }
  df
}

# ==============================================================================
# 0. Load Real Data and Extract Endmembers
# ==============================================================================

cat("Loading real phenology data to extract true endmembers...\n")

# Try to load from common data file locations
data_file <- NULL
possible_paths <- c(
  "C:\\Users\\yolan\\Downloads\\LS_S2_Harmonized_Timeseries.csv"
)

for (path in possible_paths) {
  if (file.exists(path)) {
    data_file <- path
    break
  }
}

if (is.null(data_file)) {
  stop("Could not find phenology data file. Please ensure landsat_lower.csv is accessible.")
}

cat(sprintf("Loading data from: %s\n", data_file))
df_raw <- fread(data_file)

# Normalize band names to canonical lower-case
df_raw <- normalize_band_names(df_raw)

if ("vegetation" %in% names(df_raw) && !"Veg" %in% names(df_raw)) {
  df_raw$Veg <- df_raw$vegetation
  cat("[NOTICE] Renamed 'vegetation' -> 'Veg' in phenology data\n")
}
# If Veg already present, be explicit about skipping GeoJSON join
if ("Veg" %in% names(df_raw)) {
  cat("[NOTICE] Found Veg values in input CSV; skipping GeoJSON join and using CSV-provided Veg values.\n")
}

# -----------------------------------------------------------------------------
# PREPARE INDICES & DETERREND (ensure this runs before extracting endmembers)
# -----------------------------------------------------------------------------
# Ensure indices are computed on df_raw before detrending
missing_indices <- setdiff(setdiff(INDICES_OF_INTEREST, "PPI"), names(df_raw))
if (length(missing_indices) > 0 && all(c('nir','red') %in% names(df_raw))) {
  cat("[INDEX SETUP] Computing missing indices from raw bands on df_raw\n")
  df_raw <- compute_indices_from_bands(df_raw)
}

# Helper: assign phenology year from a date (used by outlier grouping)
assign_pheno_year <- function(d) {
  d <- as.Date(d)
  ifelse(is.na(d), NA_integer_, ifelse(lubridate::month(d) >= 3, lubridate::year(d), lubridate::year(d) - 1))
}

# Early dust filtering (NDDI) to remove contaminated observations before processing
# NOTE: snow-index (previously NDSI) removed; pipeline uses dust-only filtering (NDDI).
eps <- 1e-9
if (all(c('red','nir') %in% names(df_raw))) df_raw$NDDI <- (as.numeric(df_raw$red) - as.numeric(df_raw$nir)) / (as.numeric(df_raw$red) + as.numeric(df_raw$nir) + eps)

if ("NDDI" %in% names(df_raw)) {
  dust_count <- sum(df_raw$NDDI > NDDI_DUST_THRESHOLD, na.rm = TRUE)
  total_before <- nrow(df_raw)
  df_raw <- df_raw[!(df_raw$NDDI > NDDI_DUST_THRESHOLD), , drop = FALSE]
  cat(sprintf("[FILTERING] Filtered out %d observations with dust contamination (NDDI > %s)\n", total_before - nrow(df_raw), .nddi_thresh_fmt()))

  # Also remove large outliers early using the spectral outlier helper
  df_raw <- remove_large_outliers(df_raw)

  # Recalculate DVI soil baseline from the *filtered* training data so PPI computations use a post-filter baseline
  # This ensures that seasonal baselines and PPI soil references are not biased by dust or gross outliers.
  tryCatch({
    compute_soil_line_slope(df_raw, assign_global_dvi = TRUE)
  }, error = function(e) {
    cat(sprintf("[SOIL BASELINE] compute_soil_line_slope failed: %s\n", e$message))
  })
}

# Compute a per-location DVI soil baseline for PPI.
# Baseline is the median of the lowest `quantile_p` fraction of DVI within each location.
compute_dvi_soil_per_location <- function(df, quantile_p = 0.10, min_samples = 5L) {
  if (is.null(df) || nrow(df) == 0) stop("[PPI] compute_dvi_soil_per_location: empty df")
  if (!"location_id" %in% names(df)) stop("[PPI] compute_dvi_soil_per_location: missing location_id")
  if (!"DVI" %in% names(df) && all(c("nir", "red") %in% names(df))) df$DVI <- as.numeric(df$nir) - as.numeric(df$red)
  if (!"DVI" %in% names(df)) stop("[PPI] compute_dvi_soil_per_location: missing DVI (and nir/red not available)")

  locs <- unique(as.character(df$location_id))
  dvi_soil_vec <- rep(NA_real_, nrow(df))
  for (loc in locs) {
    idx <- which(as.character(df$location_id) == loc)
    vals <- df$DVI[idx]
    vals <- vals[is.finite(vals)]
    if (length(vals) < as.integer(min_samples)) next
    q <- suppressWarnings(as.numeric(stats::quantile(vals, probs = quantile_p, na.rm = TRUE, names = FALSE, type = 7)))
    if (!is.finite(q)) next
    low_vals <- vals[vals <= q]
    if (length(low_vals) < 2L) next
    soil <- suppressWarnings(as.numeric(stats::median(low_vals, na.rm = TRUE)))
    if (!is.finite(soil)) next
    dvi_soil_vec[idx] <- soil
  }

  need_idx <- is.finite(df$DVI) & !is.finite(dvi_soil_vec)
  if (any(need_idx)) {
    bad_locs <- unique(as.character(df$location_id[need_idx]))
    stop(sprintf("[PPI] Cannot compute per-location dvi_soil for %d rows across %d locations (example locs: %s)",
                 sum(need_idx), length(bad_locs), paste(head(bad_locs, 10), collapse = ", ")))
  }
  dvi_soil_vec
}

# Ensure PPI exists; require per-location dvi_soil and per-location M (no auto fallbacks)
if (!"PPI" %in% names(df_raw) || all(!is.finite(df_raw$PPI))) {
  if (!"location_id" %in% names(df_raw)) stop("[INDEX SETUP] Cannot compute per-location PPI: missing location_id")
  if (!"DVI" %in% names(df_raw) && all(c("nir", "red") %in% names(df_raw))) df_raw$DVI <- as.numeric(df_raw$nir) - as.numeric(df_raw$red)
  if (!"DVI" %in% names(df_raw) || !any(is.finite(df_raw$DVI))) {
    stop("[INDEX SETUP] Cannot compute PPI because DVI is missing or invalid")
  }
  dvi_soil_vec <- compute_dvi_soil_per_location(df_raw, quantile_p = 0.10)
  df_raw <- add_ppi_columns(df_raw, dvi_soil = dvi_soil_vec)
  cat("[INDEX SETUP] PPI computed via add_ppi_columns() using per-location dvi_soil baseline and per-location M\n")
}

# Instead of fitting seasonal models, use June-August medians per user request
if (!"date" %in% names(df_raw)) stop("Phenology data must include a 'date' column")
if (!lubridate::is.Date(df_raw$date)) df_raw$date <- as.Date(df_raw$date)
cat(sprintf("[INDEX DETREND] Skipping seasonal detrending; using months %s (June-Aug) median per location instead\n", paste(SUMMER_DETREND_MONTHS, collapse = ",")))

# Keep original index to assign back robustly (create before filtering)
df_raw$.orig_row <- seq_len(nrow(df_raw))

# Prepare masks and containers
summer_mask <- lubridate::month(df_raw$date) %in% SUMMER_DETREND_MONTHS
indices_available <- intersect(INDICES_OF_INTEREST, names(df_raw))
if (length(indices_available) == 0) stop("[INDEX DETREND] No indices available to compute June-Aug medians; ensure indices are computed on raw data.")

INDEX_SUMMER_MEDIANS <- list()

# Create "_median" columns but populate them with raw index values only for June-August rows
for (idx in indices_available) {
  df_raw[[paste0(idx, "_median")]] <- NA_real_
  valid_idx <- which(summer_mask & is.finite(df_raw[[idx]]))
  if (length(valid_idx) > 0) {
    df_raw[[paste0(idx, "_median")]][valid_idx] <- df_raw[[idx]][valid_idx]
    INDEX_SUMMER_MEDIANS[[idx]] <- median(df_raw[[idx]][valid_idx], na.rm = TRUE)
  } else {
    INDEX_SUMMER_MEDIANS[[idx]] <- NA_real_
  }
}

# Make sure globals exist for backward compatibility
assign("INDEX_SUMMER_MEDIANS", INDEX_SUMMER_MEDIANS, envir = globalenv())

cat(sprintf("[INDEX DETREND] Created %d '_median' columns populated with June-Aug values\n", length(indices_available)))

# Report location_id diagnostics if provided in CSV
if ("location_id" %in% names(df_raw)) {
  unique_loc <- length(unique(df_raw$location_id))
  cat(sprintf("DEBUG: Unique location_id values in CSV: %d\n", unique_loc))
  first5 <- paste(head(unique(df_raw$location_id), 5), collapse = ", ")
  cat(sprintf("DEBUG: First 5 location_id values in CSV: %s \n", first5))
  cat("Using existing 'location_id' from CSV.\n")
}
# Legacy 'no.soil'/'no_soil' columns are ignored and not used by this script

# Require that Veg be present to extract endmembers
if (!("Veg" %in% names(df_raw))) {
  stop("Phenology data must include a 'Veg' column. GeoJSON support has been removed — include this field in landsat_lower.csv.")
}

# Aggregate summer medians per location to get one representative sample per location (ensures ONE pure sample per location)
cat("\nAggregating summer medians per location to extract representative endmembers...\n")
if (!"date" %in% names(df_raw)) stop("Phenology data must include a 'date' column for aggregation")
if (!lubridate::is.Date(df_raw$date)) df_raw$date <- as.Date(df_raw$date)

# Ensure summer-median indices are available
median_cols <- intersect(paste0(INDICES_OF_INTEREST, "_median"), names(df_raw))
if (length(median_cols) == 0) stop("[INDEX DETREND] No '_median' index columns found in data; cannot build location representatives.")

# Select summer observations that have at least one finite median index value
summer_rows <- df_raw |> dplyr::filter(lubridate::month(date) %in% SUMMER_DETREND_MONTHS)
summer_rows <- summer_rows |> dplyr::filter(dplyr::if_any(dplyr::all_of(median_cols), is.finite))
if (nrow(summer_rows) == 0) stop("No summer observations with median indices available to aggregate endmembers")

# Group by location_id and Veg and compute medians of raw bands (one row per location-Veg pair)
# ALSO: compute per-location median of summer indices so we have ONE representative per location
agg_group <- local({
  d_cols <- intersect(paste0(INDICES_OF_INTEREST, "_median"), names(summer_rows))
  summer_rows |> dplyr::group_by(location_id, Veg) |>
    dplyr::summarise(
      n = dplyr::n(),
      dplyr::across(all_of(RAW_BANDS), ~ median(.x, na.rm = TRUE)),
      dplyr::across(all_of(d_cols), ~ median(.x, na.rm = TRUE)),
      .groups = "drop"
    )
})

# For barren endmember: select aggregated rows where Veg == 'barren'
barren_rows <- agg_group |> dplyr::filter(tolower(Veg) == 'barren')
if (nrow(barren_rows) == 0) stop("No barren location representatives found! Cannot extract true soil endmember.")
cat(sprintf("Found %d barren location representatives\n", nrow(barren_rows)))
# Use a random subset of barren representatives as soil endmembers (do NOT average across all barren locations)
# This allows variability in soil endmembers per run; configure sample size with FVC_SOIL_SAMPLE_SIZE (global) and seed with FVC_SAMPLING_SEED
n_barren <- nrow(barren_rows)
# Default to up to 5 sampled soils, bounded by available rows
n_sample <- min(if (exists("FVC_SOIL_SAMPLE_SIZE", envir = globalenv())) as.integer(get("FVC_SOIL_SAMPLE_SIZE", envir = globalenv())) else 5, n_barren)
if (n_sample < 1) n_sample <- 1
# Apply sampling seed if present to allow reproducibility
if (exists("FVC_SAMPLING_SEED", envir = globalenv())) set.seed(as.integer(get("FVC_SAMPLING_SEED", envir = globalenv())))
sel_idx <- sample(seq_len(n_barren), size = n_sample, replace = FALSE)
# Keep only sampled barren representatives for downstream mixing
barren_rows <- barren_rows[sel_idx, , drop = FALSE]
cat(sprintf("Selected %d random barren location representatives (rows: %s)\n", n_sample, paste(sel_idx, collapse = ",")))
# Compute soil spectrum as mean across the sampled representatives (for display only)
soil_spec <- sapply(RAW_BANDS, function(b) mean(barren_rows[[b]], na.rm = TRUE))
cat("True soil endmember spectrum (from sampled barren representatives):\n")
print(soil_spec)

# For pure vegetation: require 'no_soil' == 1 and restrict Veg to 'populus'
ALLOWED_VEG <- c("populus")
no_soil_col <- NULL
if ("no_soil" %in% names(df_raw)) no_soil_col <- "no_soil"

# Normalize Veg for matching
summer_rows$Veg <- tolower(as.character(summer_rows$Veg))

if (!is.null(no_soil_col)) {
  # Identify locations that have truthy no_soil flag during summer (only these are allowed)
  truthy <- tolower(trimws(as.character(summer_rows[[no_soil_col]]))) %in% c("1","true","t","yes","y")
  locs_with_flag <- unique(summer_rows$location_id[truthy])
  if (length(locs_with_flag) == 0) stop(sprintf("CRITICAL: '%s' column present but no truthy values found in summer observations; aborting. Only observations with %s==1 are permitted for '%s' pure vegetation selection.", no_soil_col, no_soil_col, paste(ALLOWED_VEG, collapse = ",")))
  # Select only locations flagged as no_soil and with allowed Veg
  veg_candidates <- agg_group |> dplyr::filter(location_id %in% locs_with_flag & tolower(Veg) %in% ALLOWED_VEG)
  cat(sprintf("Using '%s' column: selected %d location representatives for pure vegetation (restricted to Veg in %s and %s==1)\n", no_soil_col, nrow(veg_candidates), paste(ALLOWED_VEG, collapse = ","), no_soil_col))
} else {
  stop("Script requires a 'no_soil' column. Only observations with no_soil==1 are permitted for populus pure vegetation selection.")
}

if (nrow(veg_candidates) == 0) stop(sprintf("No pure vegetation location representatives found for Veg in %s with %s==1; cannot extract vegetation endmember.", paste(ALLOWED_VEG, collapse = ","), no_soil_col))

# Ensure only one representative per location: pick the Veg with most summer samples per location
# Restrict candidate rows to allowed Veg and (when available) locations flagged with no_soil
veg_rows_for_choice <- summer_rows |> dplyr::filter(tolower(Veg) %in% ALLOWED_VEG)
if (!is.null(no_soil_col)) {
  veg_rows_for_choice <- veg_rows_for_choice |> dplyr::filter(location_id %in% locs_with_flag)
}
veg_by_loc <- veg_rows_for_choice |>
  dplyr::group_by(location_id, Veg) |>
  dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
  dplyr::group_by(location_id) |>
  dplyr::slice_max(n, n = 1, with_ties = FALSE) |>
  dplyr::ungroup()
# Join to get band medians for selected location-veg pairs
veg_rep <- dplyr::inner_join(veg_by_loc, agg_group, by = c("location_id", "Veg"))
if (nrow(veg_rep) == 0) stop("Failed to derive per-location pure vegetation representatives")

cat(sprintf("Selected %d unique location representatives for pure vegetation\n", nrow(veg_rep)))
# Compute mean vegetation spectrum across selected representatives
veg_spec <- sapply(RAW_BANDS, function(b) mean(veg_rep[[b]], na.rm = TRUE))
cat("True pure vegetation endmember spectrum (from location representatives):\n")
print(veg_spec)

# Also set pure_veg_data and barren_data structures for downstream code compatibility
pure_veg_data <- as.data.frame(veg_rep)
barren_data <- as.data.frame(barren_rows)

# Extract TRUE vegetation endmembers for each type from pure vegetation data
cat("\nExtracting TRUE vegetation type endmembers from pure vegetation data...\n")
veg_types <- c("populus")  # Restrict to 'populus' per user request
veg_endmembers <- list()

for (veg_type in veg_types) {
  veg_subset <- pure_veg_data |>
    dplyr::filter(!is.na(Veg) & tolower(Veg) == veg_type)

  if (nrow(veg_subset) > 0) {
    veg_endmembers[[veg_type]] <- c(
      blue  = mean(veg_subset$blue, na.rm = TRUE),
      green = mean(veg_subset$green, na.rm = TRUE),
      red   = mean(veg_subset$red, na.rm = TRUE),
      nir   = mean(veg_subset$nir, na.rm = TRUE),
      swir1 = mean(veg_subset$swir1, na.rm = TRUE),
      swir2 = mean(veg_subset$swir2, na.rm = TRUE)
    )
    cat(sprintf("  %s: %d pure samples (location representatives)\n", veg_type, nrow(veg_subset)))
  } else {
    cat(sprintf("  %s: NO pure samples found, skipping\n", veg_type))
  }
}

if (length(veg_endmembers) == 0) {
  stop("No vegetation type endmembers could be extracted!")
}

cat("\nExtracted endmembers for vegetation types:\n")
print(names(veg_endmembers))

# ==============================================================================
# 2. Define Index Calculation Function
# ==============================================================================

# Source PPI helpers for consistent calculation
if (file.exists("ppi_helpers.R")) {
  source("ppi_helpers.R")
  cat("Loaded ppi_helpers.R for PPI calculation\n")
} else {
  warning("ppi_helpers.R not found - PPI will not be calculated")
}

calculate_indices <- function(df) {
  # Expects columns: blue, green, red, nir, swir1, swir2
  eps <- 1e-9

  df[, `:=`(
    # Original Set
    DVI   = nir - red,
    OSAVI = (nir - red) / (nir + red + 0.16),
    MCARI = ((red - green) - 0.2*(red - blue)) * (red / (green + eps)),
    #CRI   = (1/(green + eps)) - (1/(red + eps)),
    #PRI   = (green - red) / (green + red + eps),
    NIRv  = (nir * ((nir - red) / (nir + red + eps))) * 1.3, # Including the 1.3 scaling
    PSRI  = (red - blue) / (nir + eps),
    NBR   = (nir - swir2) / (nir + swir2 + eps),
    TCW   = (swir1 - swir2) / (swir1 + swir2 + eps),
    #TCG   = (green - red) / (green + red + eps),
    #MNDWI = (green - swir1) / (green + swir1 + eps),

    # New Additions
    NDVI   = (nir - red) / (nir + red + eps),
    MSAVI2 = (2 * nir + 1 - sqrt(pmax(0, (2 * nir + 1)^2 - 8 * (nir - red)))) / 2,
    MSAVI  = (2 * nir + 1 - sqrt(pmax(0, (2 * nir + 1)^2 - 8 * (nir - red)))) / 2,
    NDMI   = (nir - swir1) / (nir + swir1 + eps),
    TCB    = 0.3029 * blue + 0.2786 * green + 0.4733 * red + 0.5599 * nir + 0.508 * swir1 + 0.1872 * swir2,
    GVI    = -0.2941 * blue - 0.243 * green - 0.5424 * red + 0.7276 * nir + 0.0713 * swir1 - 0.1608 * swir2,
    SATVI  = (swir1 - red) / (swir1 + red + 0.5) * (1 + 0.5),
    EVI    = 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)
  )]

  # Add PPI if ppi helpers are available
  if (exists("ppi") && exists("calculate_solar_zenith")) {
    # For synthetic mixing: use soil endmember DVI as dvi_soil
    # Assume first row (fraction_veg = 0) is pure soil
    if ("fraction_veg" %in% names(df)) {
      dvi_soil_val <- df$DVI[df$fraction_veg == 0][1]
      if (is.finite(dvi_soil_val)) {
        # Calculate zenith angle (use typical mid-latitude, mid-season value)
        # Assume lat=40°N, DOY=180 (summer solstice), 10:30 AM
        zenith_rad <- calculate_solar_zenith(lat = 40, doy = 180, hour = 10.5)

        # Calculate PPI using ppi() function
        M_val <- suppressWarnings(max(df$DVI, na.rm = TRUE))
        if (!is.finite(M_val)) stop("[PPI] Synthetic mixing: cannot compute finite M")
        df$PPI <- ppi(dvi = df$DVI, zenith.angle = zenith_rad, M = M_val, dvi.soil = dvi_soil_val)
        cat(sprintf("Calculated PPI for synthetic mixing (dvi_soil=%.6f, zenith=%.4f rad)\n", dvi_soil_val, zenith_rad))
      }
    }
  }

  return(df)
}


# Raw spectral bands (optional - included if present)
RAW_BANDS <- c("blue", "green", "red", "nir", "swir1", "swir2")

# Compute a set of spectral indices from raw bands (compatible with R_extract_hls.R's formulae)
compute_indices_from_bands <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  eps <- 1e-9
  has_bands <- intersect(RAW_BANDS, names(df))
  if (length(has_bands) == 0) return(df)

  # Safely compute each index only if its required bands are present
  if (all(c('nir','red') %in% names(df))) df$DVI <- as.numeric(df$nir) - as.numeric(df$red)
  if (all(c('nir','red') %in% names(df))) df$OSAVI <- (as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + as.numeric(df$red) + 0.16)
  if (all(c('red','green','blue') %in% names(df))) df$MCARI <- ((as.numeric(df$red) - as.numeric(df$green)) - 0.2*(as.numeric(df$red) - as.numeric(df$blue))) * (as.numeric(df$red) / (as.numeric(df$green) + eps))
  if (all(c('green','red') %in% names(df))) df$PRI <- (as.numeric(df$green) - as.numeric(df$red)) / (as.numeric(df$green) + as.numeric(df$red) + eps)
  if (all(c('nir','red') %in% names(df))) df$NIRv <- as.numeric(df$nir) * ((as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + as.numeric(df$red) + eps))
  if (all(c('red','blue','nir') %in% names(df))) df$PSRI <- (as.numeric(df$red) - as.numeric(df$blue)) / (as.numeric(df$nir) + eps)
  if (all(c('nir','swir2') %in% names(df))) df$NBR <- (as.numeric(df$nir) - as.numeric(df$swir2)) / (as.numeric(df$nir) + as.numeric(df$swir2) + eps)
  if (all(c('swir1','swir2') %in% names(df))) df$TCW <- (as.numeric(df$swir1) - as.numeric(df$swir2)) / (as.numeric(df$swir1) + as.numeric(df$swir2) + eps)
  if (all(c('green','red','nir','swir1','swir2','blue') %in% names(df))) {
    df$TCB <- 0.3029 * as.numeric(df$blue) + 0.2786 * as.numeric(df$green) + 0.4733 * as.numeric(df$red) + 0.5599 * as.numeric(df$nir) + 0.508 * as.numeric(df$swir1) + 0.1872 * as.numeric(df$swir2)
    df$GVI <- -0.2941 * as.numeric(df$blue) - 0.243 * as.numeric(df$green) - 0.5424 * as.numeric(df$red) + 0.7276 * as.numeric(df$nir) + 0.0713 * as.numeric(df$swir1) - 0.1608 * as.numeric(df$swir2)
  }
  if (all(c('nir','red') %in% names(df))) df$NDVI <- (as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + as.numeric(df$red) + eps)
  if (all(c('nir','red') %in% names(df))) df$MSAVI2 <- (2 * as.numeric(df$nir) + 1 - sqrt(pmax(0, (2 * as.numeric(df$nir) + 1)^2 - 8 * (as.numeric(df$nir) - as.numeric(df$red))))) / 2
  if (all(c('nir','red') %in% names(df))) df$MSAVI <- (2 * as.numeric(df$nir) + 1 - sqrt(pmax(0, (2 * as.numeric(df$nir) + 1)^2 - 8 * (as.numeric(df$nir) - as.numeric(df$red))))) / 2
  if (all(c('nir','swir1') %in% names(df))) df$NDMI <- (as.numeric(df$nir) - as.numeric(df$swir1)) / (as.numeric(df$nir) + as.numeric(df$swir1) + eps)

  # Apply kNDVI-like tweak used in extractor (NIRv * 1.3)
  if ('NIRv' %in% names(df)) df$NIRv <- df$NIRv * 1.3

  df
}

# Function to calculate robust variance using STL decomposition + MAD
# Helper: compute MAD^2 with minimal sample requirement (top-level helper, no nested defs)
compute_mad2 <- function(x, min_samples = 3) {
  x <- x[is.finite(x)]
  if (length(x) < min_samples) return(NA_real_)
  m <- mad(x, na.rm = TRUE, constant = 1.4826)
  if (!is.finite(m)) return(NA_real_)
  m^2
}

# Safe multiplication helper: explicitly handle recycling and avoid implicit warnings.
# Returns element-wise product (with explicit replication of the shorter vector when appropriate),
# or throws an informative error when lengths are incompatible.
safe_mul_vec <- function(a, b, allow_recycle = TRUE, caller = NULL) {
  la <- length(a); lb <- length(b)
  if (la == 0 || lb == 0) return(numeric(0))
  if (la == lb) return(a * b)
  if (la == 1) return(rep(a, lb) * b)
  if (lb == 1) return(a * rep(b, la))
  if (allow_recycle && (la %% lb == 0 || lb %% la == 0)) {
    # replicate shorter to the longer length explicitly
    if (la < lb) a <- rep(a, length.out = lb) else b <- rep(b, length.out = la)
    return(a * b)
  }
  caller_text <- if (is.null(caller)) "safe_mul_vec" else paste0(caller, ": ")
  stop(sprintf("%sIncompatible lengths for multiplication: %d vs %d", caller_text, la, lb))
}

# Safe dot product (handles NA removal and length mismatch via replication rules above)
safe_dot <- function(a, b, na.rm = TRUE) {
  if (length(a) == 0 || length(b) == 0) return(0)
  prod <- safe_mul_vec(a, b, allow_recycle = TRUE, caller = "safe_dot")
  if (na.rm) sum(prod, na.rm = TRUE) else sum(prod)
}

# Safe column-weighted average helper: ensures wts length equals number of rows
safe_col_weighted_avg <- function(mat, wts) {
  if (is.null(mat) || nrow(mat) == 0) return(rep(0, ifelse(is.null(mat), 1, ncol(mat))))
  n <- nrow(mat)
  if (length(wts) == 0) wts <- rep(1, n)
  if (length(wts) != n) {
    # replicate shorter vector to length n when compatible, otherwise error
    if (length(wts) == 1 || n %% length(wts) == 0 || length(wts) %% n == 0) {
      wts <- rep(wts, length.out = n)
      warning(sprintf("safe_col_weighted_avg: adjusted weights vector to length %d", n))
    } else {
      stop(sprintf("safe_col_weighted_avg: weights length (%d) incompatible with rows in mat (%d)", length(wts), n))
    }
  }
  if (sum(wts, na.rm = TRUE) == 0) {
    # fallback: unweighted mean
    return(as.numeric(colMeans(mat, na.rm = TRUE)))
  }
  wts <- as.numeric(wts) / sum(wts, na.rm = TRUE)
  as.numeric(colSums(mat * wts, na.rm = TRUE))
}

# ==============================================================================
# 3. Perform Synthetic Mixing
# ==============================================================================

# Create mixing fractions (0 to 100% vegetation)
fractions <- seq(0, 1, by = 0.01)

# -----------------------------------------------------------------------------
# Compute PPI first (summer medians used as representative values) so we have one stable index value per loc-year
# -----------------------------------------------------------------------------
# Compute indices (DVI etc.) on raw data if not present
if (!"DVI" %in% names(df_raw) || !all(c("nir","red") %in% names(df_raw))) {
  df_raw <- compute_indices_from_bands(df_raw)
}

# Ensure PPI exists; try to compute using median barren DVI if no global baseline
if (!"PPI" %in% names(df_raw) || all(!is.finite(df_raw$PPI))) {
  dvi_soil_val <- NA_real_
  if ("Veg" %in% names(df_raw) && any(tolower(df_raw$Veg) == 'barren', na.rm = TRUE)) {
    dvi_soil_val <- median(df_raw$DVI[tolower(df_raw$Veg) == 'barren'], na.rm = TRUE)
  }
  if (!is.finite(dvi_soil_val)) stop("Cannot determine dvi_soil baseline for PPI detrending; provide barren samples or set 'dvi_soil' explicitly")
  zenith_rad <- calculate_solar_zenith(lat = 40, doy = 180, hour = 10.5)
  M_val <- suppressWarnings(max(df_raw$DVI, na.rm = TRUE))
  if (!is.finite(M_val)) stop("[PPI] Cannot compute finite M for df_raw")
  df_raw$PPI <- ppi(df_raw$DVI, zenith.angle = zenith_rad, M = M_val, dvi.soil = dvi_soil_val)
}

# Detrend summer indices (June-Sep) using polynomial fit of DOY for each index
if (!"date" %in% names(df_raw)) stop("Phenology data must include a 'date' column for detrending")
if (!lubridate::is.Date(df_raw$date)) df_raw$date <- as.Date(df_raw$date)

# Reset .orig_row to match the current filtered df_raw row indices
df_raw$.orig_row <- seq_len(nrow(df_raw))

cat(sprintf("[INDEX DETREND] Using months %s for summer median representative (no seasonal detrending)\n", paste(SUMMER_DETREND_MONTHS, collapse = ",")))

# Ensure median columns exist (created earlier). If not, create them now using June-Aug medians
if (!all(paste0(INDICES_OF_INTEREST, "_median") %in% names(df_raw))) {
  cat("[INDEX DETREND] Median columns not found — creating June-Aug '_median' columns now\n")
  INDEX_SUMMER_MEDIANS <- list()
  summer_mask <- lubridate::month(df_raw$date) %in% SUMMER_DETREND_MONTHS
  indices_available <- intersect(INDICES_OF_INTEREST, names(df_raw))
  for (idx in indices_available) {
    df_raw[[paste0(idx, "_median")]] <- NA_real_
    valid_idx <- which(summer_mask & is.finite(df_raw[[idx]]))
    if (length(valid_idx) > 0) {
      df_raw[[paste0(idx, "_median")]][valid_idx] <- df_raw[[idx]][valid_idx]
      INDEX_SUMMER_MEDIANS[[idx]] <- median(df_raw[[idx]][valid_idx], na.rm = TRUE)
    } else {
      INDEX_SUMMER_MEDIANS[[idx]] <- NA_real_
    }
  }
  assign("INDEX_SUMMER_MEDIANS", INDEX_SUMMER_MEDIANS, envir = globalenv())
} else {
  cat("[INDEX DETREND] Median columns present — proceeding with June-Aug medians\n")
  if (!exists("INDEX_SUMMER_MEDIANS", envir = globalenv())) {
    idx_means <- lapply(intersect(INDICES_OF_INTEREST, names(df_raw)), function(idx) median(df_raw[[idx]][lubridate::month(df_raw$date) %in% SUMMER_DETREND_MONTHS], na.rm = TRUE))
    names(idx_means) <- intersect(INDICES_OF_INTEREST, names(df_raw))
    assign("INDEX_SUMMER_MEDIANS", idx_means, envir = globalenv())
  }
}

# Build synthetic mixtures for ALL possible pure combinations (every soil x veg pair)
# Each soil and veg are taken from aggregated location-year endmembers (one per loc-year)
# Convert to plain data.frames to avoid data.table column-selection quirks
soil_rows <- as.data.frame(barren_data)
veg_rows <- as.data.frame(pure_veg_data)

n_soil <- as.numeric(nrow(soil_rows))
n_veg <- as.numeric(nrow(veg_rows))
n_frac <- length(fractions)

npairs <- n_soil * n_veg
# Allow overriding sample size via global FVC_PAIR_SAMPLE_SIZE, otherwise default to 500
max_pairs_to_sample <- if (exists("FVC_PAIR_SAMPLE_SIZE", envir = globalenv())) as.integer(get("FVC_PAIR_SAMPLE_SIZE", envir = globalenv())) else 500L
# Optional: allow reproducible sampling via global FVC_SAMPLING_SEED
if (exists("FVC_SAMPLING_SEED", envir = globalenv())) {
  set.seed(as.integer(get("FVC_SAMPLING_SEED", envir = globalenv())))
}

# If the full cartesian product is small, use all pairs; otherwise sample up to max_pairs_to_sample unique pairs
if (!is.finite(npairs) || npairs <= 0) stop("Invalid counts for soil/veg endmembers")

if (npairs <= max_pairs_to_sample) {
  total_samples <- npairs * n_frac
  if (!is.finite(total_samples)) total_samples <- Inf
  cat(sprintf("Generating synthetic mixtures for full set: %.0f soil x %.0f veg pairs = %.0f pairs; total samples = %.0f\n", n_soil, n_veg, npairs, total_samples))
  
  mixtures_list <- vector("list", npairs)
  pair_idx <- 1
  bands <- names(veg_spec)
  for (i in seq_len(n_soil)) {
    soil_vec <- as.numeric(soil_rows[i, bands])
    names(soil_vec) <- bands
    for (j in seq_len(n_veg)) {
      veg_vec <- as.numeric(veg_rows[j, bands])
      names(veg_vec) <- bands
      df_pair <- data.frame(fraction_veg = fractions)
      for (b in bands) {
        df_pair[[b]] <- fractions * veg_vec[b] + (1 - fractions) * soil_vec[b]
      }
      # Preserve veg type in synthetic mixtures so we can color points by Veg in diagnostics
      df_pair$Veg <- as.character(veg_rows$Veg[j])
      df_pair$soil_idx <- i
      df_pair$veg_idx <- j
      df_pair$pair_id <- pair_idx
      mixtures_list[[pair_idx]] <- df_pair
      pair_idx <- pair_idx + 1
    }
  }
  mixtures <- do.call(rbind, mixtures_list)
  rm(mixtures_list)
} else {
  desired_pairs <- as.integer(max_pairs_to_sample)
  cat(sprintf("Sampling %d random unique soil x veg pairs out of %.0f possible pairs (this limits total samples to %.0f)\n", desired_pairs, npairs, desired_pairs * n_frac))
  
  # Random unique pair sampler (avoids full-cartesian generation)
  pair_keys <- character(0)
  # Loop to collect unique pairs; sample in batches for efficiency
  while (length(pair_keys) < desired_pairs) {
    need <- desired_pairs - length(pair_keys)
    s <- sample.int(n_soil, size = need * 3, replace = TRUE)
    v <- sample.int(n_veg, size = need * 3, replace = TRUE)
    new_keys <- paste(s, v, sep = "_")
    unique_new <- setdiff(unique(new_keys), pair_keys)
    if (length(unique_new) > 0) {
      take <- min(length(unique_new), need)
      pair_keys <- c(pair_keys, unique_new[1:take])
    }
  }
  
  # Parse pair keys into numeric indices
  pair_idx_mat <- do.call(rbind, strsplit(pair_keys, "_"))
  pair_idx_mat <- apply(pair_idx_mat, 2, as.integer)
  sel_soils <- pair_idx_mat[,1]
  sel_vegs <- pair_idx_mat[,2]
  
  bands <- names(veg_spec)
  mixtures_list <- vector("list", length(sel_soils))
  for (k in seq_along(sel_soils)) {
    i <- sel_soils[k]
    j <- sel_vegs[k]
    soil_vec <- as.numeric(soil_rows[i, bands])
    names(soil_vec) <- bands
    veg_vec <- as.numeric(veg_rows[j, bands])
    names(veg_vec) <- bands
    df_pair <- data.frame(fraction_veg = fractions)
    for (b in bands) {
      df_pair[[b]] <- fractions * veg_vec[b] + (1 - fractions) * soil_vec[b]
    }
    # Preserve veg type in synthetic mixtures so we can color points by Veg in diagnostics
    df_pair$Veg <- as.character(veg_rows$Veg[j])
    df_pair$soil_idx <- i
    df_pair$veg_idx <- j
    df_pair$pair_id <- k
    mixtures_list[[k]] <- df_pair
  }
  mixtures <- do.call(rbind, mixtures_list)
  rm(mixtures_list)
}

# Created synthetic mixtures: log and confirm method
cat(sprintf("Created %d synthetic mixtures (pairwise samples). Synthetic mixing uses raw bands %s in 1%% increments producing fraction_veg 0..1.\n",
            nrow(mixtures), paste(bands, collapse = ",")))

# Calculate indices for all mixtures (uses raw bands)
# Use compute_indices_from_bands (base-R variant) to avoid data.table-specific behavior
mixtures <- compute_indices_from_bands(as.data.frame(mixtures))
# Ensure EVI is present on synthetic mixtures (compute it if missing)
if (all(c('nir','red','blue') %in% names(mixtures)) && !('EVI' %in% names(mixtures))) {
  mixtures$EVI <- 2.5 * (as.numeric(mixtures$nir) - as.numeric(mixtures$red)) / (as.numeric(mixtures$nir) + 6 * as.numeric(mixtures$red) - 7.5 * as.numeric(mixtures$blue) + 1)
}

# Diagnostic: report which indices are available and how many finite samples each has
FVC_DIAGNOSTICS <- setdiff(names(mixtures), c("fraction_veg", RAW_BANDS, "soil_idx", "veg_idx", "pair_id", "seasonal_trend", "doy"))
cat(sprintf("[FVC DIAGNOSTIC] Candidate indices: %s\n", paste(FVC_DIAGNOSTICS, collapse = ", ")))
for (idx in FVC_DIAGNOSTICS) {
  vals <- mixtures[[idx]]
  n_finite <- sum(is.finite(vals))
  n_unique <- length(unique(na.omit(vals)))
  cat(sprintf("  [%s] finite=%d unique_vals=%d\n", idx, n_finite, n_unique))
}
if ("PPI_median" %in% names(mixtures)) {
  cat(sprintf("[FVC DIAGNOSTIC] PPI_median finite=%d total_rows=%d\n", sum(is.finite(mixtures$PPI_median)), nrow(mixtures)))
} 

# Extra diagnostic: check raw band columns in mixtures
cat("[FVC DIAGNOSTIC] Raw band overview:\n")
for (b in RAW_BANDS) {
  if (b %in% names(mixtures)) {
    colv <- mixtures[[b]]
    n_finite <- sum(is.finite(colv))
    n_na <- sum(is.na(colv))
    n_unique <- length(unique(na.omit(colv)))
    ctype <- class(colv)
    cat(sprintf("  [%s] class=%s finite=%d na=%d unique=%d\n", b, paste(ctype, collapse = "/"), n_finite, n_na, n_unique))
  } else {
    cat(sprintf("  [%s] MISSING\n", b))
  }
}
cat("[FVC DIAGNOSTIC] Head of mixtures (first 5 rows):\n")
print(utils::head(mixtures, 5))


# Calculate PPI per pair using pair-specific soil baseline (DVI at fraction_veg==0)
if ("DVI" %in% names(mixtures)) {
  zenith_rad <- calculate_solar_zenith(lat = 40, doy = 180, hour = 10.5)
  mixtures$PPI <- NA_real_
  # Defensive: ensure median PPI column exists to avoid warnings
  if (!"PPI_median" %in% names(mixtures)) mixtures$PPI_median <- NA_real_
  pair_ids <- unique(mixtures$pair_id)
  for (pid in pair_ids) {
    idxs <- mixtures$pair_id == pid
    dvi_soil_val <- mixtures$DVI[idxs & mixtures$fraction_veg == 0][1]
    if (!is.finite(dvi_soil_val)) next
    M_val <- suppressWarnings(max(mixtures$DVI[idxs], na.rm = TRUE))
    if (!is.finite(M_val)) next
    mixtures$PPI[idxs] <- ppi(dvi = mixtures$DVI[idxs], zenith.angle = zenith_rad, M = M_val, dvi.soil = dvi_soil_val)
  }
  cat(sprintf("Calculated PPI for synthetic mixtures across %d pairs\n", length(pair_ids)))

  # For synthetic mixtures, create '*_median' columns from raw index values (no seasonal models used)
  mixtures$doy <- 180  # assume mid-summer reference for synthetic mixes
  idx_means <- if (exists("INDEX_SUMMER_MEDIANS", envir = globalenv())) get("INDEX_SUMMER_MEDIANS", envir = globalenv()) else list()
  for (idx in intersect(names(idx_means), names(mixtures))) {
    mixtures[[paste0(idx, "_median")]] <- mixtures[[idx]]
  }
}

# ==============================================================================
# Build Linear FVC Regression Models from Synthetic Mixing
# ONE GLOBAL MODEL FOR ALL LOCATIONS using summer median values
# ==============================================================================
cat("\n=== BUILDING GLOBAL FVC CALIBRATION MODELS ===\n")
cat("Using synthetic mixing of summer-aggregated endmembers\n")
cat("Building ONE linear model per index for ALL locations\n\n")

# Fit linear regression: FVC ~ index_value for each spectral index
# This creates ONE global calibration curve to estimate FVC from index values
FVC_MODELS <- list()
# Build list of candidate indices: prefer median versions if present
FVC_INDICES <- character(0)
for (idx in INDICES_OF_INTEREST) {
  dcol <- paste0(idx, "_median")
  if (dcol %in% names(mixtures)) {
    FVC_INDICES <- c(FVC_INDICES, dcol)
  } else if (idx %in% names(mixtures)) {
    FVC_INDICES <- c(FVC_INDICES, idx)
  }
}
FVC_INDICES <- unique(FVC_INDICES)
if (length(FVC_INDICES) == 0) {
  # Fallback: any numeric index columns (exclude raw bands and metadata)
  cand <- setdiff(names(mixtures), c("fraction_veg", RAW_BANDS, "soil_idx", "veg_idx", "pair_id", "doy"))
  FVC_INDICES <- cand[sapply(cand, function(x) is.numeric(mixtures[[x]]))]
}
# Ensure EVI is included when available (prefer median form)
if ("EVI_median" %in% names(mixtures) && !("EVI_median" %in% FVC_INDICES)) FVC_INDICES <- c(FVC_INDICES, "EVI_median")
if ("EVI" %in% names(mixtures) && !("EVI" %in% FVC_INDICES) && !("EVI_median" %in% FVC_INDICES)) FVC_INDICES <- c(FVC_INDICES, "EVI")
cat(sprintf("[FVC] Using indices for calibration: %s\n", paste(FVC_INDICES, collapse = ", ")))

if ("PPI_median" %in% FVC_INDICES) FVC_INDICES <- c(setdiff(FVC_INDICES, "PPI"), "PPI_median")

for (idx in FVC_INDICES) {
  if (!idx %in% names(mixtures)) next

  # Prepare data: known FVC (fraction_veg) and observed index value
  model_data <- data.frame(
    fvc = mixtures$fraction_veg,
    index_val = mixtures[[idx]]
  )
  model_data <- model_data[is.finite(model_data$index_val), ]

  if (nrow(model_data) < 10) {
    cat(sprintf("  [%s] Skipping: insufficient finite samples (%d)\n", idx, nrow(model_data)))
    next
  }

  # 1) Linear (robust) fit
  lin_fit <- tryCatch({
    MASS::rlm(fvc ~ index_val, data = model_data, psi = MASS::psi.huber, maxit = 100)
  }, error = function(e) { cat(sprintf("  [%s] Linear fit failed: %s\n", idx, e$message)); NULL })

  # Compute linear diagnostics (unweighted RMSE/R² for comparability)
  rmse_lin <- NA_real_; r2_lin <- NA_real_
  if (!is.null(lin_fit)) {
    pred_lin <- predict(lin_fit, newdata = model_data)
    res_lin <- model_data$fvc - pred_lin
    rmse_lin <- sqrt(mean(res_lin^2, na.rm = TRUE))
    ss_res_lin <- sum(res_lin^2, na.rm = TRUE)
    ss_tot_lin <- sum((model_data$fvc - mean(model_data$fvc, na.rm = TRUE))^2, na.rm = TRUE)
    r2_lin <- if (ss_tot_lin > 0) 1 - ss_res_lin / ss_tot_lin else NA_real_
  }

  # 2) Exponential fit: FVC = A * (1 - exp(B * index_val)) (saturating form)
  exp_fit <- NULL; rmse_exp <- NA_real_; r2_exp <- NA_real_
  if (length(unique(model_data$index_val)) >= 4) {
    # sensible starting values
    stA <- min(1, max(0.5, quantile(model_data$fvc, 0.9, na.rm = TRUE)))
    range_idx <- max(model_data$index_val, na.rm = TRUE) - min(model_data$index_val, na.rm = TRUE)
    stB <- ifelse(abs(range_idx) < 1e-6, -1, -1 / (range_idx + 1e-6))

    exp_fit <- tryCatch({
      nls(fvc ~ A * (1 - exp(B * index_val)), data = model_data, start = list(A = stA, B = stB), control = nls.control(maxiter = 50, warnOnly = TRUE))
    }, error = function(e) { cat(sprintf("  [%s] Exponential fit failed: %s\n", idx, e$message)); NULL })

    if (!is.null(exp_fit)) {
      pred_exp <- predict(exp_fit, newdata = model_data)
      pred_exp[!is.finite(pred_exp)] <- NA_real_
      res_exp <- model_data$fvc - pred_exp
      rmse_exp <- sqrt(mean(res_exp^2, na.rm = TRUE))
      ss_res_exp <- sum(res_exp^2, na.rm = TRUE)
      ss_tot_exp <- sum((model_data$fvc - mean(model_data$fvc, na.rm = TRUE))^2, na.rm = TRUE)
      r2_exp <- if (ss_tot_exp > 0) 1 - ss_res_exp / ss_tot_exp else NA_real_
    }
  }

  # Additional: test Menten-style transform for ALL indices: index' = (a + index) / (b + index)
  # Keep a simple transform (a = 1, b = a/2) but fit transformed-model candidates for every index
  menten_a <- 1; menten_b <- menten_a / 2
  menten_best <- list(type = NULL, model = NULL, rmse = NA_real_, r2 = NA_real_, transform = c(a = menten_a, b = menten_b))

  # build transformed predictor for all indices
  model_data_trans <- data.frame(fvc = model_data$fvc, index_val = (menten_a + model_data$index_val) / (menten_b + model_data$index_val))

  # Linear on transformed predictor
  lin_t <- tryCatch({ MASS::rlm(fvc ~ index_val, data = model_data_trans, psi = MASS::psi.huber, maxit = 100) }, error = function(e) NULL)
  if (!is.null(lin_t)) {
    pred_lin_t <- predict(lin_t, newdata = model_data_trans)
    res_lin_t <- model_data_trans$fvc - pred_lin_t
    rmse_lin_t <- sqrt(mean(res_lin_t^2, na.rm = TRUE))
    ss_res_lin_t <- sum(res_lin_t^2, na.rm = TRUE)
    ss_tot_lin_t <- sum((model_data_trans$fvc - mean(model_data_trans$fvc, na.rm = TRUE))^2, na.rm = TRUE)
    r2_lin_t <- if (ss_tot_lin_t > 0) 1 - ss_res_lin_t / ss_tot_lin_t else NA_real_
  } else { rmse_lin_t <- NA_real_; r2_lin_t <- NA_real_ }

  # Exponential on transformed predictor
  exp_t <- NULL; rmse_exp_t <- NA_real_; r2_exp_t <- NA_real_
  if (length(unique(model_data_trans$index_val)) >= 4) {
    stA_t <- min(1, max(0.5, quantile(model_data_trans$fvc, 0.9, na.rm = TRUE)))
    range_idx_t <- max(model_data_trans$index_val, na.rm = TRUE) - min(model_data_trans$index_val, na.rm = TRUE)
    stB_t <- ifelse(abs(range_idx_t) < 1e-6, -1, -1 / (range_idx_t + 1e-6))
    exp_t <- tryCatch({ nls(fvc ~ A * (1 - exp(B * index_val)), data = model_data_trans, start = list(A = stA_t, B = stB_t), control = nls.control(maxiter = 50, warnOnly = TRUE)) }, error = function(e) NULL)
    if (!is.null(exp_t)) {
      pred_exp_t <- predict(exp_t, newdata = model_data_trans); pred_exp_t[!is.finite(pred_exp_t)] <- NA_real_
      res_exp_t <- model_data_trans$fvc - pred_exp_t
      rmse_exp_t <- sqrt(mean(res_exp_t^2, na.rm = TRUE))
      ss_res_exp_t <- sum(res_exp_t^2, na.rm = TRUE)
      ss_tot_exp_t <- sum((model_data_trans$fvc - mean(model_data_trans$fvc, na.rm = TRUE))^2, na.rm = TRUE)
      r2_exp_t <- if (ss_tot_exp_t > 0) 1 - ss_res_exp_t / ss_tot_exp_t else NA_real_
    }
  }

  # Michaelis-like on transformed predictor (fit Vmax and Km; a_hill fixed)
  mm_t <- NULL; rmse_mm_t <- NA_real_; r2_mm_t <- NA_real_
  if (length(unique(model_data_trans$index_val)) >= 4) {
    a_hill <- 2

    # Basic data guards
    model_df <- model_data_trans[is.finite(model_data_trans$index_val) & is.finite(model_data_trans$fvc), , drop = FALSE]
    if (nrow(model_df) < 4 || var(model_df$fvc, na.rm = TRUE) < 1e-8) {
      write_debug(sprintf("[JANUARY] mm_t: insufficient or low-variance data (n=%d, var=%.3g) - skipping", nrow(model_df), var(model_df$fvc, na.rm = TRUE)))
    } else {
      # sensible starting grids (tries a few nearby starts to avoid singular gradient failures)
      stKm_t <- as.numeric(quantile(model_df$index_val, 0.5, na.rm = TRUE))
      stV_t <- min(1, max(0.5, quantile(model_df$fvc, 0.9, na.rm = TRUE)))
      km_candidates <- unique(pmax(1e-6, as.numeric(quantile(model_df$index_val, c(0.25, 0.5, 0.75), na.rm = TRUE))))
      v_candidates <- unique(pmax(1e-6, c(stV_t, as.numeric(quantile(model_df$fvc, c(0.8, 0.9), na.rm = TRUE)))))

      mm_local <- NULL
      # Try multiple starts with nlsLM (preferred) and log warnings but don't convert them to errors
      for (sk in km_candidates) {
        for (sv in v_candidates) {
          if (!is.null(mm_local)) break
          if (requireNamespace("minpack.lm", quietly = TRUE)) {
            try({
              withCallingHandlers({
                tmp <- minpack.lm::nlsLM(
                  fvc ~ Vmax * (index_val^a_hill) / (Km^a_hill + index_val^a_hill),
                  data = model_df,
                  start = list(Vmax = sv, Km = sk),
                  lower = c(Vmax = 1e-6, Km = 1e-6),
                  upper = c(Vmax = max(1, max(model_df$fvc, na.rm = TRUE) * 2), Km = max(model_df$index_val, na.rm = TRUE) * 10),
                  control = minpack.lm::nls.lm.control(maxiter = 500)
                )
                mm_local <- tmp
              }, warning = function(w) { write_debug(sprintf("[JANUARY] mm_t (nlsLM) warning: %s (start V=%.3g Km=%.3g)", w$message, sv, sk)); invokeRestart("muffleWarning") })
            }, silent = TRUE)
          }
        }
        if (!is.null(mm_local)) break
      }

      # If still NULL try bounded nls (algorithm='port') with same start grid
      if (is.null(mm_local)) {
        for (sk in km_candidates) {
          for (sv in v_candidates) {
            if (!is.null(mm_local)) break
            try({
              withCallingHandlers({
                tmp <- tryCatch({
                  nls(
                    fvc ~ Vmax * (index_val^a_hill) / (Km^a_hill + index_val^a_hill),
                    data = model_df,
                    start = list(Vmax = sv, Km = sk),
                    algorithm = "port",
                    lower = c(Vmax = 1e-6, Km = 1e-6),
                    upper = c(Vmax = max(1, max(model_df$fvc, na.rm = TRUE) * 2), Km = max(model_df$index_val, na.rm = TRUE) * 10),
                    control = nls.control(maxiter = 500, warnOnly = TRUE, minFactor = 1e-12)
                  )
                }, error = function(e) { write_debug(sprintf("[JANUARY] mm_t (nls port) failed: %s (start V=%.3g Km=%.3g)", e$message, sv, sk)); NULL })
                if (!is.null(tmp)) mm_local <- tmp
              }, warning = function(w) { write_debug(sprintf("[JANUARY] mm_t (nls port) warning: %s (start V=%.3g Km=%.3g)", w$message, sv, sk)); invokeRestart("muffleWarning") })
            }, silent = TRUE)
          }
          if (!is.null(mm_local)) break
        }
      }

      if (is.null(mm_local)) write_debug("[JANUARY] mm_t: all mm attempts failed (falling back)")
      mm_t <- mm_local
    }

    if (!is.null(mm_t)) {
      pred_mm_t <- predict(mm_t, newdata = model_data_trans); pred_mm_t[!is.finite(pred_mm_t)] <- NA_real_
      res_mm_t <- model_data_trans$fvc - pred_mm_t
      rmse_mm_t <- sqrt(mean(res_mm_t^2, na.rm = TRUE))
      ss_res_mm_t <- sum(res_mm_t^2, na.rm = TRUE)
      ss_tot_mm_t <- sum((model_data_trans$fvc - mean(model_data_trans$fvc, na.rm = TRUE))^2, na.rm = TRUE)
      r2_mm_t <- if (ss_tot_mm_t > 0) 1 - ss_res_mm_t / ss_tot_mm_t else NA_real_
    }
  }

  # Pick best of transformed fits
  transformed_candidates <- list(linear = list(model = lin_t, rmse = rmse_lin_t, r2 = r2_lin_t), exponential = list(model = exp_t, rmse = rmse_exp_t, r2 = r2_exp_t), michaelis = list(model = mm_t, rmse = rmse_mm_t, r2 = r2_mm_t))
  best_t <- NULL; best_rmse_t <- Inf; best_type_t <- NULL; best_r2_t <- NA_real_
  for (nm in names(transformed_candidates)) {
    ent <- transformed_candidates[[nm]]
    if (!is.null(ent$model) && is.finite(ent$rmse) && ent$rmse < best_rmse_t) {
      best_rmse_t <- ent$rmse; best_t <- ent$model; best_type_t <- nm; best_r2_t <- ent$r2
    }
  }
  if (!is.null(best_t)) {
    menten_best$type <- best_type_t; menten_best$model <- best_t; menten_best$rmse <- best_rmse_t; menten_best$r2 <- best_r2_t
  }

  # 3) Michaelis-Menten / Hill fit on raw predictor: fit Vmax and Km (a_hill fixed)
  mm_fit <- NULL; rmse_mm <- NA_real_; r2_mm <- NA_real_
  if (length(unique(model_data$index_val)) >= 4) {
    a_hill <- 2
    # sensible starting values
    stKm <- as.numeric(quantile(model_data$index_val, 0.5, na.rm = TRUE))
    stV <- min(1, max(0.5, quantile(model_data$fvc, 0.9, na.rm = TRUE)))

    mm_fit <- tryCatch({
      mm_local <- NULL
      if (requireNamespace("minpack.lm", quietly = TRUE)) {
        mm_local <- tryCatch({
          minpack.lm::nlsLM(
            fvc ~ Vmax * (index_val^a_hill) / (Km^a_hill + index_val^a_hill),
            data = model_data,
            start = list(Vmax = stV, Km = stKm),
            lower = c(Vmax = 1e-6, Km = 1e-6),
            upper = c(Vmax = max(1, max(model_data$fvc, na.rm = TRUE) * 2), Km = max(model_data$index_val, na.rm = TRUE) * 10),
            control = minpack.lm::nls.lm.control(maxiter = 500)
          )
        }, error = function(e) { cat(sprintf("  [%s] nlsLM mm_fit failed: %s\n", idx, e$message)); NULL })
      }
      if (is.null(mm_local)) {
        mm_local <- tryCatch({
          nls(fvc ~ Vmax * (index_val^a_hill) / (Km^a_hill + index_val^a_hill), data = model_data, start = list(Vmax = stV, Km = stKm), control = nls.control(maxiter = 500, warnOnly = TRUE, minFactor = 1e-12))
        }, error = function(e) { cat(sprintf("  [%s] Michaelis-Menten fit failed unexpectedly: %s\n", idx, e$message)); NULL })
      }
      mm_local
    }, error = function(e) { cat(sprintf("  [%s] Michaelis-Menten fit failed unexpectedly: %s\n", idx, e$message)); NULL })

    if (!is.null(mm_fit)) {
      pred_mm <- predict(mm_fit, newdata = model_data)
      pred_mm[!is.finite(pred_mm)] <- NA_real_
      res_mm <- model_data$fvc - pred_mm
      rmse_mm <- sqrt(mean(res_mm^2, na.rm = TRUE))
      ss_res_mm <- sum(res_mm^2, na.rm = TRUE)
      ss_tot_mm <- sum((model_data$fvc - mean(model_data$fvc, na.rm = TRUE))^2, na.rm = TRUE)
      r2_mm <- if (ss_tot_mm > 0) 1 - ss_res_mm / ss_tot_mm else NA_real_
      # store estimated Vmax info in the model attributes for plotting
      if (exists("Vmax", where = as.list(coef(mm_fit)))) attr(mm_fit, "Vmax_est") <- coef(mm_fit)["Vmax"]
    }
  }

  # Choose best model by RMSE (prefer lower RMSE)
  chosen_type <- NA_character_; chosen_model <- NULL; chosen_rmse <- NA_real_; chosen_r2 <- NA_real_; chosen_use_menten <- FALSE

  # Determine best among raw fits
  raw_candidates <- list(linear = list(model = lin_fit, rmse = rmse_lin, r2 = r2_lin), exponential = list(model = exp_fit, rmse = rmse_exp, r2 = r2_exp), michaelis = list(model = mm_fit, rmse = rmse_mm, r2 = r2_mm))
  best_raw <- NULL; best_raw_rmse <- Inf; best_raw_type <- NULL; best_raw_r2 <- NA_real_
  for (nm in names(raw_candidates)) {
    ent <- raw_candidates[[nm]]
    if (!is.null(ent$model) && is.finite(ent$rmse) && ent$rmse < best_raw_rmse) {
      best_raw_rmse <- ent$rmse; best_raw <- ent$model; best_raw_type <- nm; best_raw_r2 <- ent$r2
    }
  }

  # Compare to transformed (menten) best if available
  best_trans_rmse <- if (!is.null(menten_best$model)) menten_best$rmse else Inf

  if (!is.null(best_raw) && best_raw_rmse <= best_trans_rmse) {
    chosen_type <- best_raw_type; chosen_model <- best_raw; chosen_rmse <- best_raw_rmse; chosen_r2 <- best_raw_r2; chosen_use_menten <- FALSE
  } else if (!is.null(menten_best$model) && is.finite(menten_best$rmse)) {
    chosen_type <- menten_best$type; chosen_model <- menten_best$model; chosen_rmse <- menten_best$rmse; chosen_r2 <- menten_best$r2; chosen_use_menten <- TRUE
  }

  if (!is.null(chosen_model)) {
    # Store as a structured entry so callers can inspect model type and whether it used the Menten transform
    FVC_MODELS[[idx]] <- list(type = chosen_type, model = chosen_model, metrics = list(rmse = chosen_rmse, r2 = chosen_r2), use_menten = chosen_use_menten, transform = if (chosen_use_menten) list(a = menten_a, b = menten_b) else NULL)
    cat(sprintf("  [%s] FVC model chosen: %s (RMSE=%.4f R²=%.4f) use_menten=%s — other: linear RMSE=%.4f R²=%.4f, exp RMSE=%.4f R²=%.4f, mm RMSE=%.4f R²=%.4f\n",
                idx, chosen_type, chosen_rmse, chosen_r2,
                ifelse(chosen_use_menten, "TRUE", "FALSE"),
                ifelse(is.finite(rmse_lin), rmse_lin, NA_real_), ifelse(is.finite(r2_lin), r2_lin, NA_real_),
                ifelse(is.finite(rmse_exp), rmse_exp, NA_real_), ifelse(is.finite(r2_exp), r2_exp, NA_real_),
                ifelse(is.finite(rmse_mm), rmse_mm, NA_real_), ifelse(is.finite(r2_mm), r2_mm, NA_real_)))
  } else {
    cat(sprintf("  [%s] All fitting attempts failed — skipping\n", idx))
  }
}

cat(sprintf("Built %d FVC calibration models from synthetic mixing\n", length(FVC_MODELS)))

# Function to estimate FVC from index values using the calibration models
estimate_fvc_from_index <- function(df, index_name, model_list = FVC_MODELS) {
  if (!index_name %in% names(model_list)) {
    cat(sprintf("[WARN] No FVC model available for index '%s'\n", index_name))
    return(rep(NA_real_, nrow(df)))
  }
  
  if (!index_name %in% names(df)) {
    cat(sprintf("[WARN] Index '%s' not found in data\n", index_name))
    return(rep(NA_real_, nrow(df)))
  }

  model_entry <- model_list[[index_name]]

  # Legacy: model_entry might be a bare model object (rlm) — treat as linear
  if (!is.list(model_entry) && (inherits(model_entry, "rlm") || inherits(model_entry, "lm"))) {
    pred_data <- data.frame(index_val = df[[index_name]])
    return(pmax(0, pmin(1, predict(model_entry, newdata = pred_data))))
  }

  # Structured entry: contains type and model
  if (is.list(model_entry) && !is.null(model_entry$type) && !is.null(model_entry$model)) {
    pred_data <- data.frame(index_val = df[[index_name]])
    if (model_entry$use_menten) {
      a <- model_entry$transform$a; b <- model_entry$transform$b
      pred_data$index_val <- (a + pred_data$index_val) / (b + pred_data$index_val)
    }
    if (model_entry$type == "linear") {
      res <- tryCatch(predict(model_entry$model, newdata = pred_data), error = function(e) NA_real_)
      return(pmax(0, pmin(1, res)))
    } else if (model_entry$type %in% c("exponential", "michaelis")) {
      # exponential and Michaelis-Menten are both non-linear nls objects; predict and sanitize
      res <- tryCatch(predict(model_entry$model, newdata = pred_data), error = function(e) NA_real_)
      res[!is.finite(res)] <- NA_real_
      return(pmax(0, pmin(1, res)))
    }
  }

  cat(sprintf("[WARN] Unrecognized model entry for index '%s'\n", index_name))
  return(rep(NA_real_, nrow(df)))
}

# Store models globally for use in other scripts
assign("FVC_CALIBRATION_MODELS", FVC_MODELS, envir = globalenv())
cat("FVC calibration models stored in global environment as 'FVC_CALIBRATION_MODELS'\n")
# Persist models to disk so other scripts (e.g., fit_veg_mixture_mesma.R) can load them without re-fitting
model_file <- file.path(output_dir, "fvc_calibration_models.rds")
tryCatch({
  # Before saving, ensure MSAVI canonical key points to the best-performing MSAVI variant
  msavi_candidates <- intersect(names(FVC_MODELS), c("MSAVI_median", "MSAVI2_median", "MSAVI", "MSAVI2"))
  if (length(msavi_candidates) > 0) {
    rmse_vals <- sapply(msavi_candidates, function(k) {
      ent <- FVC_MODELS[[k]]
      if (!is.null(ent) && is.list(ent) && !is.null(ent$metrics$rmse)) return(as.numeric(ent$metrics$rmse))
      NA_real_
    })
    if (all(is.na(rmse_vals))) {
      cat("[FVC] MSAVI variants present but RMSE not available; leaving keys as-is.\n")
    } else {
      chosen_idx <- which.min(rmse_vals)
      chosen_key <- msavi_candidates[chosen_idx]
      cat(sprintf("[FVC] Selected MSAVI variant for canonical 'MSAVI' key: %s (RMSE=%.4f)\n", chosen_key, rmse_vals[chosen_idx]))
      # Copy chosen model to canonical 'MSAVI' key so MESMA/fitting code can rely on 'MSAVI'
      FVC_MODELS[["MSAVI"]] <- FVC_MODELS[[chosen_key]]
      # record which variant was chosen for traceability
      attr(FVC_MODELS[["MSAVI"]], "source_variant") <- chosen_key
    }
  }
  cat(sprintf("Saved FVC calibration models to: %s\n", model_file))
}, error = function(e) {
  warning(sprintf("Failed to save FVC calibration models to disk: %s", e$message))
})

# === BUILDING VEG-TYPE-SPECIFIC FVC CALIBRATION MODELS ===
cat("\n=== BUILDING VEG-TYPE-SPECIFIC FVC CALIBRATION MODELS ===\n")

# Treat agri/agriculture as the same vegtype (aliases)
AGRI_ALIASES <- c("agri", "agric", "agriculture", "agricultural")
# Veg types to build per-type models for (agri ignores no_soil flag; others require no_soil==1)
NEW_VEG_TYPES <- c("agri", "tamarix", "populus", "herbs")
# Per-veg color mapping (named vector)
# woody_unknown is a blend color between tamarix and populus for indistinguishable variants
veg_colors <- c("agri" = "#1b9e77", "tamarix" = "#d95f02", "populus" = "#7570b3", "herbs" = "#e7298a", "woody_unknown" = "#a0522d")
assign("VEG_CALIBRATION_COLORS", veg_colors, envir = globalenv())

FVC_MODELS_BY_VEG <- list()

# Helper: build mixtures and calibration models for a given set of soil/veg representatives
build_models_for_soil_veg <- function(soil_rows, veg_rows, veg_label) {
  if (nrow(veg_rows) == 0) {
    cat(sprintf("[FVC][%s] No veg representatives found; skipping\n", veg_label))
    return(NULL)
  }
  soils <- as.data.frame(soil_rows)
  vegs <- as.data.frame(veg_rows)
  fractions <- seq(0, 1, by = 0.01)

  n_soil <- nrow(soils); n_veg <- nrow(vegs); npairs <- n_soil * n_veg
  max_pairs_to_sample <- if (exists("FVC_PAIR_SAMPLE_SIZE", envir = globalenv())) as.integer(get("FVC_PAIR_SAMPLE_SIZE", envir = globalenv())) else 500L

  mixtures_list <- list()
  bands <- intersect(RAW_BANDS, names(soils))

  if (npairs <= max_pairs_to_sample) {
    pair_idx <- 1
    for (i in seq_len(n_soil)) {
      for (j in seq_len(n_veg)) {
        soil_spec <- as.numeric(soils[i, bands, drop = TRUE])
        veg_spec <- as.numeric(vegs[j, bands, drop = TRUE])
        mat <- sapply(fractions, function(f) veg_spec * f + soil_spec * (1 - f))
        dfm <- as.data.frame(t(mat))
        names(dfm) <- bands
        dfm$fraction_veg <- fractions
        dfm$soil_idx <- i
        dfm$veg_idx <- j
        dfm$pair_id <- sprintf("%d_%d", i, j)
        mixtures_list[[pair_idx]] <- dfm
        pair_idx <- pair_idx + 1
      }
    }
  } else {
    desired_pairs <- as.integer(max_pairs_to_sample)
    sel_pairs <- sample(npairs, size = desired_pairs)
    k <- 1
    for (pnum in sel_pairs) {
      p <- pnum - 1
      i <- (p %% n_soil) + 1
      j <- (p %/% n_soil) + 1
      soil_spec <- as.numeric(soils[i, bands, drop = TRUE])
      veg_spec <- as.numeric(vegs[j, bands, drop = TRUE])
      mat <- sapply(fractions, function(f) veg_spec * f + soil_spec * (1 - f))
      dfm <- as.data.frame(t(mat))
      names(dfm) <- bands
      dfm$fraction_veg <- fractions
      dfm$soil_idx <- i
      dfm$veg_idx <- j
      dfm$pair_id <- sprintf("%d_%d", i, j)
      mixtures_list[[k]] <- dfm
      k <- k + 1
    }
  }

  mixtures <- do.call(rbind, mixtures_list)
  mixtures <- compute_indices_from_bands(as.data.frame(mixtures))

  # Ensure EVI present
  if (all(c('nir','red','blue') %in% names(mixtures)) && !('EVI' %in% names(mixtures))) {
    mixtures$EVI <- 2.5 * (as.numeric(mixtures$nir) - as.numeric(mixtures$red)) / (as.numeric(mixtures$nir) + 6 * as.numeric(mixtures$red) - 7.5 * as.numeric(mixtures$blue) + 1)
  }

  # Calculate PPI per pair using the same logic as main pipeline
  if ("DVI" %in% names(mixtures)) {
    zenith_rad <- calculate_solar_zenith(lat = 40, doy = 180, hour = 10.5)
    mixtures$PPI <- NA_real_
    if (!"PPI_median" %in% names(mixtures)) mixtures$PPI_median <- NA_real_
    pair_ids <- unique(mixtures$pair_id)
    for (pid in pair_ids) {
      psub <- mixtures[mixtures$pair_id == pid, , drop = FALSE]
      soil_dvi <- psub$DVI[which(psub$fraction_veg == 0)[1]]
      if (!is.finite(soil_dvi)) soil_dvi <- NA_real_
      M_val <- suppressWarnings(max(psub$DVI, na.rm = TRUE))
      if (!is.finite(M_val) || !is.finite(soil_dvi)) {
        mixtures$PPI[mixtures$pair_id == pid] <- NA_real_
      } else {
        mixtures$PPI[mixtures$pair_id == pid] <- ppi(psub$DVI, zenith.angle = zenith_rad, M = M_val, dvi.soil = soil_dvi)
      }
    }
    mixtures$doy <- 180
    for (idx in INDICES_OF_INTEREST) {
      col_median <- paste0(idx, "_median")
      if (idx %in% names(mixtures)) mixtures[[col_median]] <- mixtures[[idx]]
    }
  }

  # Fit linear regression FVC ~ index_value for all candidate indices
  fmodels <- list()
  FVC_INDICES_local <- character(0)
  for (idx in INDICES_OF_INTEREST) {
    if (paste0(idx, "_median") %in% names(mixtures)) FVC_INDICES_local <- c(FVC_INDICES_local, paste0(idx, "_median"))
    if (idx %in% names(mixtures)) FVC_INDICES_local <- c(FVC_INDICES_local, idx)
  }
  FVC_INDICES_local <- unique(FVC_INDICES_local)
  if ("PPI_median" %in% FVC_INDICES_local) FVC_INDICES_local <- c(setdiff(FVC_INDICES_local, "PPI"), "PPI_median")

  for (idx in FVC_INDICES_local) {
    vals <- mixtures[[idx]]
    ok <- is.finite(vals)
    if (sum(ok) < 10) next
    dfm <- data.frame(fvc = mixtures$fraction_veg[ok], x = vals[ok])
    m <- lm(fvc ~ x, data = dfm)
    preds <- predict(m, newdata = dfm)
    rmse <- sqrt(mean((preds - dfm$fvc)^2))
    r2 <- 1 - sum((preds - dfm$fvc)^2)/sum((dfm$fvc - mean(dfm$fvc))^2)
    fmodels[[idx]] <- list(model = m, metrics = list(rmse = rmse, r2 = r2, n = nrow(dfm)))
  }
  return(list(models = fmodels, mixtures = mixtures))
}

# Build a full set of synthetic mixtures using all vegetation representatives (agg_group)
# for the purpose of coloring combined plots by Veg type. Keep the existing 'mixtures'
# variable (populus-only) for the populus pipeline/plots.
if (exists("agg_group") && nrow(agg_group) > 0) {
  veg_pool <- as.data.frame(agg_group)
  veg_pool$Veg <- tolower(trimws(as.character(veg_pool$Veg)))
  cat(sprintf("Building full mixtures using %d veg representatives (for coloring only)\n", nrow(veg_pool)))
  res_allveg <- build_models_for_soil_veg(barren_data, veg_pool, "allveg")
  mixtures_allveg <- if (!is.null(res_allveg) && is.list(res_allveg) && !is.null(res_allveg$mixtures)) res_allveg$mixtures else NULL
  if (!is.null(mixtures_allveg)) {
    mixtures_allveg$Veg <- as.character(veg_pool$Veg[mixtures_allveg$veg_idx])
    mixtures_allveg$Veg[is.na(mixtures_allveg$Veg)] <- NA_character_
    mixtures_allveg$Veg <- tolower(trimws(as.character(mixtures_allveg$Veg)))
    cat(sprintf("Annotated full mixtures with Veg labels (unique: %s)\n", paste(sort(unique(na.omit(mixtures_allveg$Veg))), collapse = ", ")))
  } else {
    cat("Warning: failed to build full mixtures for all veg types; combined plots will use the current mixtures variable (populus-only).\n")
    mixtures_allveg <- NULL
  }
} else {
  mixtures_allveg <- NULL
  cat("agg_group not available; cannot build full veg mixtures for coloring.\n")
}

# Decide which indices to plot: use FVC_DIAGNOSTICS (already computed) but skip PPI/NDVI/EVI
skip_idxs <- c(if (exists("ppi_col")) ppi_col else NA_character_, if (exists("ndvi_col")) ndvi_col else NA_character_, if (exists("evi_col")) evi_col else NA_character_)
skip_idxs <- na.omit(skip_idxs)
plot_candidates <- setdiff(FVC_DIAGNOSTICS, skip_idxs)

# Use mixtures_allveg if available for combined plots; else fall back to mixtures
mixtures_for_plots <- if (!is.null(mixtures_allveg)) mixtures_allveg else mixtures

# Prepare color mapping from mixtures_for_plots
if (exists("VEG_CALIBRATION_COLORS", envir = globalenv())) {
  veg_palette <- get("VEG_CALIBRATION_COLORS", envir = globalenv())
} else {
  veg_levels <- sort(unique(na.omit(as.character(mixtures_for_plots$Veg))))
  pal <- grDevices::rainbow(max(1, length(veg_levels)))
  veg_palette <- setNames(pal, veg_levels)
}

for (idx in plot_candidates) {
  if (!idx %in% names(mixtures_for_plots)) next
  plot_df <- mixtures_for_plots[is.finite(mixtures_for_plots[[idx]]) & is.finite(mixtures_for_plots$fraction_veg), , drop = FALSE]
  if (nrow(plot_df) <= 10) next

  # Determine point colors by Veg when available
  point_cols <- rep(rgb(0,0,0,0.5), nrow(plot_df))
  if ("Veg" %in% names(plot_df) && any(!is.na(plot_df$Veg))) {
    veg_levels <- sort(unique(na.omit(as.character(plot_df$Veg))))
    for (vl in veg_levels) {
      sel <- which(as.character(plot_df$Veg) == vl)
      col <- if (vl %in% names(veg_palette)) veg_palette[vl] else rgb(0,0,0,0.5)
      point_cols[sel] <- col
    }
  }

  # Determine best model to overlay when available (use global FVC_MODELS or median variant)
  model_entry <- NULL
  if (exists("FVC_MODELS", envir = globalenv())) {
    models_all <- get("FVC_MODELS", envir = globalenv())
    if (idx %in% names(models_all)) model_entry <- models_all[[idx]]
    else if (paste0(idx, "_median") %in% names(models_all)) model_entry <- models_all[[paste0(idx, "_median")]]
  }
  x_seq <- seq(min(plot_df[[idx]], na.rm = TRUE), max(plot_df[[idx]], na.rm = TRUE), length.out = 200)
  y_pred_grid <- rep(NA_real_, length(x_seq))
  rmse <- NA_real_; r2 <- NA_real_
  if (!is.null(model_entry) && !is.null(model_entry$model)) {
    # Use stored calibration model
    y_pred_grid <- tryCatch(predict(model_entry$model, newdata = data.frame(index_val = x_seq)), error = function(e) rep(NA_real_, length(x_seq)))
    y_pred_grid[!is.finite(y_pred_grid)] <- NA_real_
    if (!is.null(model_entry$metrics)) {
      rmse <- as.numeric(model_entry$metrics$rmse)
      r2 <- as.numeric(model_entry$metrics$r2)
    }
    legend_txt <- c(sprintf("Model: %s", ifelse(!is.null(model_entry$type), model_entry$type, "linear")), sprintf("RMSE=%.4f", rmse), sprintf("R²=%.4f", r2))
  } else {
    # Fallback: robust linear fit to the pooled data and compute metrics for legend
    fit_rlm <- tryCatch({ MASS::rlm(as.formula(paste0("fraction_veg ~ ", idx)), data = plot_df, psi = MASS::psi.huber, maxit = 100) }, error = function(e) NULL)
    if (!is.null(fit_rlm)) {
      y_pred_grid <- tryCatch(predict(fit_rlm, newdata = data.frame(x = x_seq)), error = function(e) rep(NA_real_, length(x_seq)))
      res <- residuals(fit_rlm)
      wts <- if (!is.null(fit_rlm$w)) fit_rlm$w else rep(1, length(res))
      rmse <- sqrt(sum(wts * (res^2)) / sum(wts))
      y <- plot_df$fraction_veg
      wy_mean <- sum(wts * y) / sum(wts)
      ss_res <- sum(wts * (res^2))
      ss_tot <- sum(wts * ((y - wy_mean)^2))
      r2 <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
      legend_txt <- c(sprintf("Linear fit"), sprintf("RMSE=%.4f", rmse), sprintf("R²=%.4f", r2))
    } else {
      legend_txt <- character(0)
    }
  }
  if (length(veg_palette) > 0 && any(nzchar(names(veg_palette)))) {
    legend_colors <- veg_palette
    legend_labels <- names(legend_colors)
  } else {
    legend_colors <- character(0); legend_labels <- character(0)
  }
  # Consolidated legend in bottom-right: model metrics followed by veg color keys
  all_leg <- c(legend_txt, legend_labels)
  leg_col <- c(rep(NA_character_, length(legend_txt)), as.character(legend_colors))
  leg_pch <- c(rep(NA_integer_, length(legend_txt)), rep(20L, length(legend_labels)))

  # Open PNG device and create the plot
  plot_file <- file.path(output_dir, sprintf("fvc_vs_%s_allveg.png", idx))
  png(plot_file, width = 800, height = 600)
  plot(plot_df[[idx]], plot_df$fraction_veg, pch = 20, col = point_cols, xlab = idx, ylab = "Fraction vegetation (ground truth)", main = sprintf("FVC vs %s (all veg types)", idx), ylim = c(0,1), xlim = c(min(x_seq, na.rm=TRUE), max(x_seq, na.rm=TRUE)))
  lines(x_seq, y_pred_grid, col = "blue", lwd = 2)

  if (length(all_leg) > 0) {
    if (dev.cur() > 1) {
      tryCatch(legend("bottomright", legend = all_leg, col = leg_col, pch = leg_pch, bty = "n"), error = function(e) warning(sprintf("Failed to draw legend for %s: %s", idx, e$message)))
    } else {
      warning(sprintf("No active plotting device for legend of %s", idx))
    }
  }
  if (dev.cur() > 1) dev.off() else warning(sprintf("No active device to close for %s plot", idx))
  cat(sprintf("[FVC PLOT][allveg][%s] Saved combined plot: %s\n", idx, plot_file))
}

# --- PRINT: Location-level median indices + Veg + derived FVC (first 10 rows only) ---
# Build list of median index column names available
median_cols <- intersect(paste0(INDICES_OF_INTEREST, "_median"), names(df_raw))
if (length(median_cols) == 0) {
  cat("[PRINT] No _median index columns available; skipping location-level print.\n")
} else {
  # Compute per-location median of summer indices using the median window (one representative per location)
  loc_median <- df_raw |>
    dplyr::filter(lubridate::month(date) %in% SUMMER_DETREND_MONTHS) |>
    dplyr::group_by(location_id) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(median_cols), ~ median(.x, na.rm = TRUE)), .groups = "drop")

  # Join with selected location representatives (veg_rep) to get Veg per location
  if (!exists("veg_rep")) {
    cat("[PRINT] 'veg_rep' not found; cannot assemble location representatives for print.\n")
  } else {
    df_print <- dplyr::inner_join(veg_rep |> dplyr::select(location_id, Veg), loc_median, by = "location_id")

    # Derive FVC estimate using available FVC models; prefer PPI_median then NDVI_median then any other available
    fvc_est <- rep(NA_real_, nrow(df_print))
    pref_idxs <- c("PPI", "NDVI", setdiff(INDICES_OF_INTEREST, c("PPI", "NDVI")))
    for (idx in pref_idxs) {
      colname <- paste0(idx, "_median")
      if (colname %in% names(df_print)) {
        preds <- estimate_fvc_from_index(df_print, colname)
        # Fill missing FVC estimates where we don't have one yet
        need <- is.na(fvc_est) & is.finite(preds)
        if (any(need)) fvc_est[need] <- preds[need]
      }
    }
    df_print$FVC_est <- fvc_est

    # Only keep asked columns: location_id, Veg, median indices, FVC_est
    keep_cols <- c("location_id", "Veg", median_cols, "FVC_est")
    df_print <- df_print[, keep_cols, drop = FALSE]

    cat("[PRINT] First 10 rows (location_id, Veg, _median indices, FVC_est):\n")
    print(utils::head(df_print, 10))
  }
}

# --- Plot FVC vs PPI (prefer median variant if available) and save R² and RMSE ---
plot_file <- file.path(output_dir, sprintf("fvc_vs_ppi_median.png"))
ppi_col <- if ("PPI_median" %in% names(mixtures)) "PPI_median" else if ("PPI" %in% names(mixtures)) "PPI" else NA_character_
if (!is.na(ppi_col)) {
  plot_df <- mixtures[is.finite(mixtures[[ppi_col]]) & is.finite(mixtures$fraction_veg), , drop = FALSE]
  if (nrow(plot_df) > 10) {
    model_entry <- if (exists("FVC_MODELS", envir = globalenv()) && ppi_col %in% names(get("FVC_MODELS", envir = globalenv())) ) get("FVC_MODELS", envir = globalenv())[[ppi_col]] else NULL

    # Plotting metrics and predicted curve: prefer stored calibration model when available
    if (!is.null(model_entry)) {
      rmse <- model_entry$metrics$rmse
      r2 <- model_entry$metrics$r2
      # Prediction grid
      x_seq <- seq(min(plot_df[[ppi_col]], na.rm = TRUE), max(plot_df[[ppi_col]], na.rm = TRUE), length.out = 200)
      if (isTRUE(model_entry$use_menten) && !is.null(model_entry$transform)) {
        a <- model_entry$transform$a; b <- model_entry$transform$b
        x_pred_in <- (a + x_seq) / (b + x_seq)
        x_lab <- sprintf("%s (Menten a=%.2f b=%.2f)", ppi_col, a, b)
      } else {
        x_pred_in <- x_seq
        x_lab <- ppi_col
      }
      y_pred_grid <- tryCatch(predict(model_entry$model, newdata = data.frame(index_val = x_pred_in)), error = function(e) rep(NA_real_, length(x_seq)))
      y_pred_grid[!is.finite(y_pred_grid)] <- NA_real_

      # Color points by Veg type if available
      if ("Veg" %in% names(plot_df)) {
        veg_levels <- unique(na.omit(as.character(plot_df$Veg)))
        cols <- grDevices::rainbow(length(veg_levels))
        cols <- grDevices::adjustcolor(cols, alpha.f = 0.5)
        col_map <- setNames(cols, veg_levels)
        point_cols <- col_map[as.character(plot_df$Veg)]
        point_cols[is.na(point_cols)] <- rgb(0,0,0,0.5)
      } else {
        point_cols <- rgb(0,0,0,0.5)
        veg_levels <- character(0)
        cols <- character(0)
      }

      ylim <- c(0, 1)

      png(plot_file, width = 800, height = 600)
      plot(plot_df[[ppi_col]], plot_df$fraction_veg, pch = 20, col = point_cols, xlab = x_lab, ylab = "Fraction vegetation (ground truth)", main = sprintf("FVC vs %s (synthetic mixtures)", ppi_col), ylim = ylim, xlim = c(0, 1))
      lines(x_seq, y_pred_grid, col = "red", lwd = 2)
      extra_info <- character(0)
      if (is.list(model_entry) && !is.null(model_entry$type) && model_entry$type == "michaelis") {
        vmax_val <- NA_real_
        km_val <- NA_real_
        if (!is.null(attr(model_entry$model, "Vmax_fixed"))) vmax_val <- attr(model_entry$model, "Vmax_fixed")
        cf <- tryCatch(coef(model_entry$model), error = function(e) NULL)
        if (!is.null(cf) && "Km" %in% names(cf)) km_val <- as.numeric(cf["Km"])
        if (is.finite(vmax_val)) extra_info <- c(extra_info, sprintf("Vmax=%.2f", vmax_val))
        if (is.finite(km_val)) extra_info <- c(extra_info, sprintf("Km=%.4f", km_val))
      }
      legend_txt <- c(sprintf("Model: %s", model_entry$type), sprintf("use_menten=%s", ifelse(isTRUE(model_entry$use_menten), "TRUE", "FALSE")), sprintf("RMSE=%.4f", rmse), sprintf("R²=%.4f", r2), extra_info)
      all_leg <- c(legend_txt, veg_levels)
      leg_col <- c(rep(NA_character_, length(legend_txt)), as.character(cols))
      leg_pch <- c(rep(NA_integer_, length(legend_txt)), rep(20L, length(veg_levels)))
      if (length(all_leg) > 0) {
        if (dev.cur() > 1) {
          tryCatch(legend("bottomright", legend = all_leg, col = leg_col, pch = leg_pch, bty = "n"), error = function(e) warning(sprintf("Failed to draw PPI legend: %s", e$message)))
        } else {
          warning("No active plotting device for PPI legend")
        }
      }
      dev.off()
      cat(sprintf("[FVC PLOT] Saved FVC vs %s plot using chosen model: %s (R²=%.4f, RMSE=%.4f)\n", ppi_col, model_entry$type, r2, rmse))

    } else {
      # fallback: simple robust linear fit for plotting
      fit_ppi <- tryCatch({ MASS::rlm(as.formula(paste0("fraction_veg ~ ", ppi_col)), data = plot_df, psi = MASS::psi.huber, maxit = 100) }, error = function(e) NULL)
      if (!is.null(fit_ppi)) {
        res <- residuals(fit_ppi)
        wts <- if (!is.null(fit_ppi$w)) fit_ppi$w else rep(1, length(res))
        rmse <- sqrt(sum(wts * (res^2)) / sum(wts))
        y <- plot_df$fraction_veg
        wy_mean <- sum(wts * y) / sum(wts)
        ss_res <- sum(wts * (res^2))
        ss_tot <- sum(wts * ((y - wy_mean)^2))
        r2 <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
        coefs <- coef(fit_ppi)
        y_pred <- coefs[1] + coefs[2] * plot_df[[ppi_col]]
        # Restrict FVC axis to [0,1] for clarity
        ylim <- c(0, 1)

        png(plot_file, width = 800, height = 600)
        # Color points by Veg type if available
        if ("Veg" %in% names(plot_df)) {
          veg_levels <- unique(na.omit(as.character(plot_df$Veg)))
          cols <- grDevices::rainbow(length(veg_levels))
          cols <- grDevices::adjustcolor(cols, alpha.f = 0.5)
          col_map <- setNames(cols, veg_levels)
          point_cols <- col_map[as.character(plot_df$Veg)]
          point_cols[is.na(point_cols)] <- rgb(0,0,0,0.5)
        } else {
          point_cols <- rgb(0,0,0,0.5)
          veg_levels <- character(0)
          cols <- character(0)
        }

        plot(plot_df[[ppi_col]], plot_df$fraction_veg, pch = 20, col = point_cols, xlab = ppi_col, ylab = "Fraction vegetation (ground truth)", main = sprintf("FVC vs %s (fallback linear)", ppi_col), ylim = ylim, xlim = c(0, 1))
        abline(a = coefs[1], b = coefs[2], col = "red", lwd = 2)
        legend_txt <- c(sprintf("Linear fit"), sprintf("RMSE=%.4f", rmse), sprintf("R²=%.4f", r2))
        all_leg <- c(legend_txt, veg_levels)
        leg_col <- c(rep(NA_character_, length(legend_txt)), as.character(cols))
        leg_pch <- c(rep(NA_integer_, length(legend_txt)), rep(20L, length(veg_levels)))
        if (length(all_leg) > 0) {
          if (dev.cur() > 1) {
            tryCatch(legend("bottomright", legend = all_leg, col = leg_col, pch = leg_pch, bty = "n"), error = function(e) warning(sprintf("Failed to draw fallback PPI legend: %s", e$message)))
          } else {
            warning("No active plotting device for fallback PPI legend")
          }
        }
        if (dev.cur() > 1) dev.off() else warning(sprintf("No active device to close for PPI plot"))
        cat(sprintf("[FVC PLOT] Saved FVC vs %s plot using fallback linear model: R²=%.4f, RMSE=%.4f\n", ppi_col, r2, rmse))
      }
    }
  } else {
      cat(sprintf("[FVC PLOT] Not enough points to plot FVC vs %s\n", ppi_col))
  }
}  # close if (!is.na(ppi_col))
  ndvi_col <- if ("NDVI_median" %in% names(mixtures)) "NDVI_median" else if ("NDVI" %in% names(mixtures)) "NDVI" else NA_character_
  if (!is.na(ndvi_col)) {
    plot_df <- mixtures[is.finite(mixtures[[ndvi_col]]) & is.finite(mixtures$fraction_veg), , drop = FALSE]
  if (nrow(plot_df) > 10) {
    model_entry <- if (exists("FVC_MODELS", envir = globalenv()) && ndvi_col %in% names(get("FVC_MODELS", envir = globalenv())) ) get("FVC_MODELS", envir = globalenv())[[ndvi_col]] else NULL
    if (!is.null(model_entry)) {
      rmse <- model_entry$metrics$rmse
      r2 <- model_entry$metrics$r2
      x_seq <- seq(min(plot_df[[ndvi_col]], na.rm = TRUE), max(plot_df[[ndvi_col]], na.rm = TRUE), length.out = 200)
      if (isTRUE(model_entry$use_menten) && !is.null(model_entry$transform)) {
        a <- model_entry$transform$a; b <- model_entry$transform$b
        x_pred_in <- (a + x_seq) / (b + x_seq)
        x_lab <- sprintf("%s (Menten a=%.2f b=%.2f)", ndvi_col, a, b)
      } else {
        x_pred_in <- x_seq
        x_lab <- ndvi_col
      }
      y_pred_grid <- tryCatch(predict(model_entry$model, newdata = data.frame(index_val = x_pred_in)), error = function(e) rep(NA_real_, length(x_seq)))
      y_pred_grid[!is.finite(y_pred_grid)] <- NA_real_

      # Use fixed filename (no date/time) to avoid embedding system timestamps
      ndvi_plot_file <- file.path(output_dir, "fvc_vs_ndvi.png")
      png(ndvi_plot_file, width = 800, height = 600)
      if ("Veg" %in% names(plot_df)) {
        veg_levels <- unique(na.omit(as.character(plot_df$Veg)))
        cols <- grDevices::rainbow(length(veg_levels))
        cols <- grDevices::adjustcolor(cols, alpha.f = 0.5)
        col_map <- setNames(cols, veg_levels)
        point_cols <- col_map[as.character(plot_df$Veg)]
        point_cols[is.na(point_cols)] <- rgb(0,0,0,0.5)
      } else {
        point_cols <- rgb(0,0,0,0.5)
        veg_levels <- character(0)
        cols <- character(0)
      }
      plot(plot_df[[ndvi_col]], plot_df$fraction_veg, pch = 20, col = point_cols, xlab = x_lab, ylab = "Fraction vegetation (ground truth)", main = sprintf("FVC vs %s (synthetic mixtures)", ndvi_col), ylim = c(0,1), xlim = c(0, 1))
      lines(x_seq, y_pred_grid, col = "blue", lwd = 2)
      extra_info <- character(0)
      if (is.list(model_entry) && !is.null(model_entry$type) && model_entry$type == "michaelis") {
        vmax_val <- NA_real_
        km_val <- NA_real_
        if (!is.null(attr(model_entry$model, "Vmax_fixed"))) vmax_val <- attr(model_entry$model, "Vmax_fixed")
        cf <- tryCatch(coef(model_entry$model), error = function(e) NULL)
        if (!is.null(cf) && "Km" %in% names(cf)) km_val <- as.numeric(cf["Km"])
        if (is.finite(vmax_val)) extra_info <- c(extra_info, sprintf("Vmax=%.2f", vmax_val))
        if (is.finite(km_val)) extra_info <- c(extra_info, sprintf("Km=%.4f", km_val))
      }
      legend_txt <- c(sprintf("Model: %s", model_entry$type), sprintf("use_menten=%s", ifelse(isTRUE(model_entry$use_menten), "TRUE", "FALSE")), sprintf("RMSE=%.4f", rmse), sprintf("R²=%.4f", r2), extra_info)
      all_leg <- c(legend_txt, veg_levels)
      leg_col <- c(rep(NA_character_, length(legend_txt)), as.character(cols))
      leg_pch <- c(rep(NA_integer_, length(legend_txt)), rep(20L, length(veg_levels)))
      if (length(all_leg) > 0) {
        if (dev.cur() > 1) {
          tryCatch(legend("bottomright", legend = all_leg, col = leg_col, pch = leg_pch, bty = "n"), error = function(e) warning(sprintf("Failed to draw NDVI legend: %s", e$message)))
        } else {
          warning("No active plotting device for NDVI legend")
        }
      }
      if (dev.cur() > 1) dev.off() else warning(sprintf("No active device to close for NDVI plot"))
      cat(sprintf("[FVC PLOT] Saved FVC vs %s plot using chosen model: %s (R²=%.4f, RMSE=%.4f)\n", ndvi_col, model_entry$type, r2, rmse))
    } else {
      cat(sprintf("[FVC PLOT] Not enough points to plot FVC vs %s\n", ndvi_col))
    }
  } else {
    cat(sprintf("[FVC PLOT] Not enough points to plot FVC vs %s\n", ndvi_col))
  }
} else {
    cat("[FVC PLOT] No NDVI/NDVI_median available in synthetic mixtures; skipping NDVI plot\n")
}

# --- EVI: Save plot and estimator (if available) ---
evi_col <- if ("EVI_median" %in% names(mixtures)) "EVI_median" else if ("EVI" %in% names(mixtures)) "EVI" else NA_character_
if (!is.na(evi_col)) {
  plot_df <- mixtures[is.finite(mixtures[[evi_col]]) & is.finite(mixtures$fraction_veg), , drop = FALSE]
  if (nrow(plot_df) > 10) {
    model_entry <- if (exists("FVC_MODELS", envir = globalenv()) && evi_col %in% names(get("FVC_MODELS", envir = globalenv())) ) get("FVC_MODELS", envir = globalenv())[[evi_col]] else NULL
    if (!is.null(model_entry)) {
      rmse <- model_entry$metrics$rmse
      r2 <- model_entry$metrics$r2
      x_seq <- seq(min(plot_df[[evi_col]], na.rm = TRUE), max(plot_df[[evi_col]], na.rm = TRUE), length.out = 200)
      if (isTRUE(model_entry$use_menten) && !is.null(model_entry$transform)) {
        a <- model_entry$transform$a; b <- model_entry$transform$b
        x_pred_in <- (a + x_seq) / (b + x_seq)
        x_lab <- sprintf("%s (Menten a=%.2f b=%.2f)", evi_col, a, b)
      } else {
        x_pred_in <- x_seq
        x_lab <- evi_col
      }
      y_pred_grid <- tryCatch(predict(model_entry$model, newdata = data.frame(index_val = x_pred_in)), error = function(e) rep(NA_real_, length(x_seq)))
      y_pred_grid[!is.finite(y_pred_grid)] <- NA_real_

      # Use fixed filename (no date/time) to avoid embedding system timestamps
      evi_plot_file <- file.path(output_dir, "fvc_vs_evi.png")
      png(evi_plot_file, width = 800, height = 600)
      if ("Veg" %in% names(plot_df)) {
        veg_levels <- unique(na.omit(as.character(plot_df$Veg)))
        cols <- grDevices::rainbow(length(veg_levels))
        cols <- grDevices::adjustcolor(cols, alpha.f = 0.5)
        col_map <- setNames(cols, veg_levels)
        point_cols <- col_map[as.character(plot_df$Veg)]
        point_cols[is.na(point_cols)] <- rgb(0,0,0,0.5)
      } else {
        point_cols <- rgb(0,0,0,0.5)
        veg_levels <- character(0)
        cols <- character(0)
      }
      plot(plot_df[[evi_col]], plot_df$fraction_veg, pch = 20, col = point_cols, xlab = x_lab, ylab = "Fraction vegetation (ground truth)", main = sprintf("FVC vs %s (synthetic mixtures)", evi_col), ylim = c(0,1), xlim = c(0, 1))
      lines(x_seq, y_pred_grid, col = "blue", lwd = 2)
      legend_txt <- c(sprintf("Model: %s", model_entry$type), sprintf("use_menten=%s", ifelse(isTRUE(model_entry$use_menten), "TRUE", "FALSE")), sprintf("RMSE=%.4f", rmse), sprintf("R²=%.4f", r2))
      all_leg <- c(legend_txt, veg_levels)
      leg_col <- c(rep(NA_character_, length(legend_txt)), as.character(cols))
      leg_pch <- c(rep(NA_integer_, length(legend_txt)), rep(20L, length(veg_levels)))
      if (length(all_leg) > 0) {
        if (dev.cur() > 1) {
          tryCatch(legend("bottomright", legend = all_leg, col = leg_col, pch = leg_pch, bty = "n"), error = function(e) warning(sprintf("Failed to draw EVI legend: %s", e$message)))
        } else {
          warning("No active plotting device for EVI legend")
        }
      }
      if (dev.cur() > 1) dev.off() else warning(sprintf("No active device to close for EVI plot"))
      cat(sprintf("[FVC PLOT] Saved FVC vs %s plot using chosen model: %s (R²=%.4f, RMSE=%.4f)\n", evi_col, model_entry$type, r2, rmse))

      # Save EVI estimator to disk alongside other models (FVC_MODELS already saved earlier), also save individual plot path to model_entry if useful
      # Add model_entry to FVC_MODELS (should already be present), then resave models file
      tryCatch({
        model_file <- file.path(output_dir, "FVC_MODELS.rds")
        if (exists("FVC_MODELS", envir = globalenv())) saveRDS(get("FVC_MODELS", envir = globalenv()), file = model_file)
        cat(sprintf("[FVC] Saved FVC_MODELS (including EVI) to: %s\n", model_file))
      }, error = function(e) cat(sprintf("[FVC] Failed to save FVC_MODELS: %s\n", e$message)))

    } else {
      cat(sprintf("[FVC PLOT] Not enough points to plot FVC vs %s\n", evi_col))
    }
  } else {
    cat(sprintf("[FVC PLOT] Not enough points to plot FVC vs %s\n", evi_col))
  }
} else {
  cat("[FVC PLOT] No EVI/EVI_median available in synthetic mixtures; skipping EVI plot\n")
}

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

# Indices and RAW_BANDS are defined earlier in the file; outlier configuration was moved earlier to make the helper available at load time.

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

## Use centralized PPI helper to ensure identical PPI calculation
if (file.exists("ppi_helpers.R")) {
  source("ppi_helpers.R")
} else {
  warning("ppi_helpers.R not found; falling back to local inline PPI logic. Please ensure ppi_helpers.R is present in project root for consistent PPI calculations.")
}

# Define input file paths
INPUT_CSV <- "C:\\Users\\yolan\\Downloads\\LS_S2_Harmonized_Timeseries.csv"

# Load the phenology data
cat("Loading data from:", INPUT_CSV, "\n")
df <- readr::read_csv(INPUT_CSV, show_col_types = FALSE)

# Drop observations from years 1992-1999 if present
if ("date" %in% names(df)) {
  # Ensure date column is Date
  if (!lubridate::is.Date(df$date)) df$date <- as.Date(df$date)
  n_before <- nrow(df)
  years_to_drop <- 1992:1999
  df <- df[!(lubridate::year(df$date) %in% years_to_drop), , drop = FALSE]
  cat(sprintf("[DATA FILTER] Dropped %d observations from years %d-%d\n", n_before - nrow(df), min(years_to_drop), max(years_to_drop)))
} else if ("year" %in% names(df)) {
  n_before <- nrow(df)
  years_to_drop <- 1992:1999
  df <- df[!(as.integer(df$year) %in% years_to_drop), , drop = FALSE]
  cat(sprintf("[DATA FILTER] Dropped %d observations from years %d-%d (using 'year' column)\n", n_before - nrow(df), min(years_to_drop), max(years_to_drop)))
} else {
  cat("[DATA FILTER] No 'date' or 'year' column found; cannot drop 1992-1999 observations automatically\n")
}

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

# Early dust filtering (NDDI) on the raw input CSV to drop contaminated rows early
eps <- 1e-9
# compute dust-only index (NDSI removed)
if (all(c('red','nir') %in% names(df))) df$NDDI <- (as.numeric(df$red) - as.numeric(df$nir)) / (as.numeric(df$red) + as.numeric(df$nir) + eps)

if ("NDDI" %in% names(df)) {
  dust_count <- sum(df$NDDI > NDDI_DUST_THRESHOLD, na.rm = TRUE)
  total_before <- nrow(df)
  df <- df[!(df$NDDI > NDDI_DUST_THRESHOLD), , drop = FALSE]
  cat(sprintf("[FILTERING] Filtered out %d observations with dust contamination (NDDI > %s)\n", total_before - nrow(df), .nddi_thresh_fmt()))

  # Remove large spectral outliers early
  df <- remove_large_outliers(df)
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
  if (length(veg_cols) > 0) {
    map_df$Veg <- as.character(map_df[[veg_cols[1]]])
    map_df$Veg <- tolower(trimws(as.character(map_df$Veg)))
    map_df$Veg <- ifelse(grepl("phragmites", map_df$Veg, ignore.case = TRUE) |
                         map_df$Veg %in% c("herbs", "alhagi", "salicornia", "halocnemum"),
                         "herbs", map_df$Veg)
  }
  
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
        if ("Veg" %in% names(training_df)) {
            training_df$Veg <- tolower(trimws(as.character(training_df$Veg)))
            training_df$Veg <- ifelse(grepl("phragmites", training_df$Veg, ignore.case = TRUE) |
                                      training_df$Veg %in% c("herbs", "alhagi", "salicornia", "halocnemum"),
                                      "herbs", training_df$Veg)
        }
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

  # Apply dust (NDDI) and MAD filtering to training data
  if (all(c('red','nir') %in% names(training_df))) training_df$NDDI <- (as.numeric(training_df$red) - as.numeric(training_df$nir)) / (as.numeric(training_df$red) + as.numeric(training_df$nir) + eps)

  if ("NDDI" %in% names(training_df)) {
    dust_count <- sum(training_df$NDDI > NDDI_DUST_THRESHOLD, na.rm = TRUE)
    total_before <- nrow(training_df)
    training_df <- training_df[!(training_df$NDDI > NDDI_DUST_THRESHOLD), , drop = FALSE]
    cat(sprintf("[TRAINING FILTERING] Filtered out %d observations with dust contamination (NDDI > %s)\n", total_before - nrow(training_df), .nddi_thresh_fmt()))
    
    # Outlier removal already applied earlier to the raw CSV; skipping duplicate removal here
  }

  # Ensure PPI is present for training data (strict per-location dvi_soil + per-location M)
  if (!"PPI" %in% names(training_df) || all(!is.finite(training_df$PPI))) {
    if (!"location_id" %in% names(training_df)) stop("[TRAINING] Cannot compute PPI: missing location_id")
    if (!"DVI" %in% names(training_df) && all(c("nir", "red") %in% names(training_df))) training_df$DVI <- as.numeric(training_df$nir) - as.numeric(training_df$red)
    dvi_soil_vec <- compute_dvi_soil_per_location(training_df, quantile_p = 0.10)
    training_df <- add_ppi_columns(training_df, dvi_soil = dvi_soil_vec)
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
if (!"pheno_year" %in% names(df)) df$pheno_year <- assign_pheno_year(df$date)
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
    # Use linear FVC regression instead of simple mean
    # Estimate FVC from DVI using calibration model, then use barren observations (FVC~0) to get soil baseline
    if (exists("FVC_CALIBRATION_MODELS") && "DVI" %in% names(FVC_CALIBRATION_MODELS)) {
      # For barren pixels, we expect FVC~0, so DVI~dvi_soil
      # Prefer linear FVC model for deriving DVI_soil (intercept/slope). If an exponential
      # model was chosen for DVI, fall back to simple barren mean.
      dvi_entry <- FVC_CALIBRATION_MODELS[["DVI"]]
      coefs <- NULL

      # Legacy: entry may be a bare model
      if (!is.list(dvi_entry) && (inherits(dvi_entry, "rlm") || inherits(dvi_entry, "lm"))) {
        coefs <- tryCatch(coef(dvi_entry), error = function(e) NULL)
      } else if (is.list(dvi_entry) && !is.null(dvi_entry$type) && dvi_entry$type == "linear") {
        coefs <- tryCatch(coef(dvi_entry$model), error = function(e) NULL)
      }

      if (!is.null(coefs) && length(coefs) == 2 && is.finite(coefs[2]) && abs(coefs[2]) > 1e-6) {
        dvi_soil_est <- -coefs[1] / coefs[2]  # FVC=0 -> DVI_soil
        cat(sprintf("[FVC MODEL] Estimated PPI_DVI_SOIL from linear FVC regression: %.4f\n", dvi_soil_est))
      } else {
        # Fallback to simple mean if model is degenerate or non-linear
        dvi_soil_est <- mean(df$DVI[barren_idx], na.rm = TRUE)
        cat(sprintf("[FALLBACK] Estimated PPI_DVI_SOIL from simple mean: %.4f\n", dvi_soil_est))
      }
    } else {
      # No model available, use simple mean
      dvi_soil_est <- mean(df$DVI[barren_idx], na.rm = TRUE)
      cat(sprintf("Estimated PPI_DVI_SOIL (from main CSV barren rows) = %.4f\n", dvi_soil_est))
    }
    assign("PPI_DVI_SOIL", dvi_soil_est, envir = globalenv())
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

## Compute SNR using estimated FVC fractions (Dynamic Range / Noise)
## Signal = Dynamic Range (99th - 1st percentile of all data) -> represents full FVC range (0-1)
## Noise = Median of location-wise MADs (on non-barren data)
compute_fvc_snr <- function(df_signal, df_noise, indices, group_col = "location_id", eps = 1e-8) {
  out <- numeric(length(indices)); names(out) <- indices
  for (idx in indices) {
    if (!idx %in% names(df_signal) || !idx %in% names(df_noise)) { out[idx] <- NA_real_; next }
    
    # 1. Estimate Dynamic Range (Signal) from full dataset (Soil to Veg)
    vals_sig <- df_signal[[idx]]
    vals_sig <- vals_sig[is.finite(vals_sig)]
    if (length(vals_sig) < 10) { out[idx] <- NA_real_; next }
    
    # Use 1st and 99th percentiles to estimate the full dynamic range (Soil to Dense Veg)
    q <- quantile(vals_sig, probs = c(0.01, 0.99), na.rm = TRUE)
    dynamic_range <- as.numeric(diff(q))
    if (dynamic_range <= eps) dynamic_range <- eps
    
    # 2. Estimate Noise from location-wise MADs
    locs <- unique(na.omit(as.character(df_noise[[group_col]])))
    noises <- numeric(0)
    
    for (loc in locs) {
      sub <- df_noise[df_noise[[group_col]] == loc, , drop = FALSE]
      vals <- sub[[idx]]
      vals <- vals[is.finite(vals)]
      if (length(vals) < 2) next
      
      # Use MAD as robust noise estimator
      noise_loc <- tryCatch(mad(vals, na.rm = TRUE), error = function(e) NA_real_)
      if (is.finite(noise_loc)) noises <- c(noises, noise_loc)
    }
    
    if (length(noises) == 0) { out[idx] <- NA_real_; next }
    median_noise <- median(noises, na.rm = TRUE)
    
    # SNR = Dynamic Range / Median Noise
    out[idx] <- dynamic_range / (median_noise + eps)
  }
  out
}

training_nonbarren <- training_df
if ("Veg" %in% names(training_df)) training_nonbarren <- training_df[!(tolower(training_df$Veg) == "barren"), , drop = FALSE]

# Use training_df (all classes) for Signal (Dynamic Range)
# Use training_nonbarren for Noise
index_snr <- compute_fvc_snr(training_df, training_nonbarren, indices_to_snr, group_col = "location_id")

cat("Index SNR (Dynamic Range / MAD):\n")
print(index_snr)
cat("Index SNR computed using estimated FVC fractions logic (Dynamic Range / Noise)\n")

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
set.seed(get_mesma_seed(123))
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
create_summary_entry <- function(idx_name, bias_info, jan_row, july_row, sept_row) {
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
  
  # Format Bias with CI
  if (!is.null(bias_info) && is.finite(bias_info$bias)) {
     margin <- (bias_info$upper - bias_info$lower) / 2
     b_str <- format(round(bias_info$bias, 4), nsmall=4, decimal.mark=",")
     e_str <- format(round(margin, 4), nsmall=4, decimal.mark=",")
     bias_str <- sprintf("%s (+/-%s)*", b_str, e_str)
     rmse_str <- format(round(bias_info$rmse, 4), nsmall=4, decimal.mark=",")
  } else {
     bias_str <- "NA"
     rmse_str <- "NA"
  }
  
  c(idx_name, bias_str, rmse_str, jan_str, july_str, sept_str)
}

# Compute per-index winter bias (Jan vs zero) using FVC calibration models and hierarchical bootstrapping
index_bias_list <- list()
for (idx in indices_to_snr) {
  # Select the best available FVC model for this index from FVC_MODELS (prefer any '_median' variant, else lowest RMSE among candidates)
  model_key <- NULL
  if (exists("FVC_MODELS")) {
    candidate_keys <- names(FVC_MODELS)[grepl(paste0("^", idx, "(?:_|$)"), names(FVC_MODELS))]
    if (length(candidate_keys) > 0) {
      rmse_vals <- sapply(candidate_keys, function(k) {
        ent <- FVC_MODELS[[k]]
        if (!is.null(ent) && is.list(ent) && !is.null(ent$metrics$rmse)) return(as.numeric(ent$metrics$rmse))
        NA_real_
      })
      if (!all(is.na(rmse_vals))) {
        model_key <- candidate_keys[which.min(rmse_vals)]
      } else {
        model_key <- candidate_keys[1]
      }
      cat(sprintf("[BIAS] Using FVC model '%s' for index '%s' (candidates: %s)\n", model_key, idx, paste(candidate_keys, collapse = ",")))
    }
  }
  if (is.null(model_key)) {
    index_bias_list[[idx]] <- NULL
    next
  }

  # Prepare a temporary copy of January data and ensure the model column exists
  tmp_jan <- january_data
  base <- sub("_median$", "", model_key)
  if (!(model_key %in% names(tmp_jan)) && base %in% names(tmp_jan)) tmp_jan[[model_key]] <- tmp_jan[[base]]

  if (!(model_key %in% names(tmp_jan))) {
    index_bias_list[[idx]] <- NULL
    next
  }

  # Predict FVC per-observation for January only using the chosen (best) model
  preds_jan <- estimate_fvc_from_index(tmp_jan, model_key)
  tmp_jan$.pred_fvc <- preds_jan

  if (!"location_id" %in% names(tmp_jan)) {
    index_bias_list[[idx]] <- NULL
    next
  }

  jan_loc <- tmp_jan |> dplyr::group_by(location_id) |> dplyr::summarise(fvc = mean(.pred_fvc, na.rm = TRUE), .groups = "drop")

  # Need at least 3 locations to estimate bias reliably
  if (nrow(jan_loc) < 3) {
    index_bias_list[[idx]] <- NULL
    next
  }

  # Point estimates: mean predicted FVC in January (bias vs zero) and RMSE vs zero
  bias_point <- mean(jan_loc$fvc, na.rm = TRUE)
  rmse_point <- sqrt(mean((jan_loc$fvc)^2, na.rm = TRUE))

  # Hierarchical bootstrap over locations to get CI (mean and RMSE vs zero)
  B <- 1000
  set.seed(get_mesma_seed(123))
  boot_bias <- numeric(B)
  boot_rmse <- numeric(B)
  locs <- jan_loc$location_id
  nloc <- length(locs)
  for (b in seq_len(B)) {
    sel <- sample(locs, nloc, replace = TRUE)
    sel_d <- jan_loc[jan_loc$location_id %in% sel, , drop = FALSE]
    boot_bias[b] <- mean(sel_d$fvc, na.rm = TRUE)
    boot_rmse[b] <- sqrt(mean((sel_d$fvc)^2, na.rm = TRUE))
  }

  index_bias_list[[idx]] <- list(
    bias = median(boot_bias, na.rm = TRUE),
    lower = as.numeric(quantile(boot_bias, 0.025, na.rm = TRUE)),
    upper = as.numeric(quantile(boot_bias, 0.975, na.rm = TRUE)),
    rmse = median(boot_rmse, na.rm = TRUE)
  )
}

summary_rows <- lapply(indices_to_snr, function(idx) {
  bias_val <- if (exists("index_bias_list")) index_bias_list[[idx]] else NULL
  create_summary_entry(idx, bias_val, jan_global_avg_nonbarren, july_global_avg_nonbarren, sept_global_avg_nonbarren)
})

summary_table <- do.call(rbind, summary_rows)
colnames(summary_table) <- c("INDEX", "Winter Bias (Est. FVC)", "Winter RMSE (Est. FVC)", "January Average (All areas, excl. barren)", "July average (All areas, excl. barren)", "September average (All areas, excl. barren)")
summary_table <- as.data.frame(summary_table)

cat("\n=== SUMMARY TABLE (excluding Veg=='barren') ===\n")
print(summary_table)

# Optional: Save results to Excel
output_dir <- "C:/MAP/january_averages_results"
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
# Add index Bias results to the output bundle
if (exists("index_bias_list") && length(index_bias_list) > 0) {
  bias_df <- do.call(rbind, lapply(names(index_bias_list), function(n) {
    b <- index_bias_list[[n]]
    data.frame(index = n, bias = b$bias, lower = b$lower, upper = b$upper, rmse = b$rmse, stringsAsFactors = FALSE)
  }))
  results_list$Index_Bias <- bias_df
}

openxlsx::write.xlsx(results_list, file = file.path(output_dir, "phenology_averages.xlsx"))
if (exists("bias_df") && nrow(bias_df) > 0) {
  # Save the Bias CSV alongside the Excel file
  try({
    write.csv(bias_df, file = file.path(output_dir, "index_bias.csv"), row.names = FALSE)
    cat("Index Bias saved to:", file.path(output_dir, "index_bias.csv"), "\n")
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
  add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
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
  add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
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

# -----------------------------------------------------------------------------
# NDVI vs PPI comparison: ensure both are plotted on the SAME y-scale for easier comparison
# Create mean-based NDVI/PPI plot with fixed y-axis across the two facets
np_mean <- plot_data_mean |> dplyr::filter(index %in% c("NDVI", "PPI"))
if (nrow(np_mean) > 0) {
  p_np_mean <- ggplot(np_mean, aes(x = pheno_year, y = mean_val, color = index)) +
    add_excluded_years_shade(is_date = FALSE) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = index), alpha = 0.15, inherit.aes = TRUE, colour = NA) +
    geom_line(size = 1.0) +
    geom_point(size = 1.5) +
    facet_wrap(~ index, scales = "fixed") +
    labs(title = "NDVI and PPI: Mean Seasonally-Normalized Trends (same y-scale)",
         x = "Year", y = "Normalized Value") +
    theme_minimal()

  ggsave(file.path(output_dir, "ndvi_ppi_summer_trend_normalized_mean.png"), plot = p_np_mean, width = 10, height = 4.5)
  cat("Saved NDVI+PPI mean comparison plot to:", file.path(output_dir, "ndvi_ppi_summer_trend_normalized_mean.png"), "\n")
} else {
  cat("No NDVI/PPI mean data available to plot comparison\n")
}

# Create median-based NDVI/PPI plot with fixed y-axis across the two facets
np_med <- plot_data_median |> dplyr::filter(index %in% c("NDVI", "PPI"))
if (nrow(np_med) > 0) {
  p_np_med <- ggplot(np_med, aes(x = pheno_year, y = median_val, color = index)) +
    add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = index), alpha = 0.15, inherit.aes = TRUE, colour = NA) +
    geom_line(size = 1.0) +
    geom_point(size = 1.5) +
    facet_wrap(~ index, scales = "fixed") +
    labs(title = "NDVI and PPI: Median Seasonally-Normalized Trends (same y-scale)",
         x = "Year", y = "Normalized Value") +
    theme_minimal()

  ggsave(file.path(output_dir, "ndvi_ppi_summer_trend_normalized_median.png"), plot = p_np_med, width = 10, height = 4.5)
  cat("Saved NDVI+PPI median comparison plot to:", file.path(output_dir, "ndvi_ppi_summer_trend_normalized_median.png"), "\n")
} else {
  cat("No NDVI/PPI median data available to plot comparison\n")
}
# -----------------------------------------------------------------------------