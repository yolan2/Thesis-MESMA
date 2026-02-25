library(data.table)
library(ggplot2)
library(dplyr)
library(lubridate)
# Robust regression (Huber weighting)
library(MASS)

# ---------------------------------------------------------------------------
# configuration flags that users can override before sourcing this script
# ---------------------------------------------------------------------------
# When TRUE we exclude the historically problematic interval 1992-1999 from
# analysis and add shading to plots.  **The default is now FALSE so that all
# years are kept**; set EXCLUDE_PRE2000 <- TRUE before sourcing if you want to
# reproduce earlier behaviour.
if (!exists("EXCLUDE_PRE2000", inherits = TRUE)) EXCLUDE_PRE2000 <- FALSE




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

# === output directory setup ===
# many later operations write files; define a sensible default up front so callers
# can override by setting the variable before sourcing this script.
if (!exists("output_dir")) {
  output_dir <- "C:/MAP/january_averages_results"
}
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Raw spectral bands (optional - included if present)
RAW_BANDS <- c("blue", "green", "red", "nir", "swir1", "swir2")
# Indices we consider for FVC calibration & diagnostics
INDICES_OF_INTEREST <- c("MSAVI", "NDVI", "PPI", "OSAVI", "NIRv", "NBR", "TCW", "NDMI", "TCB", "GVI", "EVI")

# Detrending configuration: which months constitute "summer" for seasonal model fitting
# and the minimum number of finite samples required per index to build a seasonal model.
SUMMER_DETREND_MONTHS <- c(7,8,9)  # July - September (use median across these months)
MIN_SEASONAL_SAMPLES <- 50

# Polynomial degree used when fitting seasonal detrend curves on summer data
# (previously a smooth spline was used across the full year).  Users may override
# by setting DETREND_POLY_DEGREE before sourcing the script.
if (!exists("DETREND_POLY_DEGREE", inherits = TRUE)) DETREND_POLY_DEGREE <- 3

# PPI configuration: previously a fixed M fallback was allowed, but we
# now require a valid M be computed from any samples with no_soil == 1.
# constant defined (for backwards compatibility) but mark it as unused.
PPI_FIXED_M <- NA_real_  # deprecated; no fallback allowed
# Soil baseline months for dvi_soil: use January-March per location mean.
PPI_SOIL_BASELINE_MONTHS <- c(1L, 2L, 3L)
# PPI M calibration: M must be strictly ABOVE the maximum observed DVI so that
# ppi() never receives dvi >= M (avoids negative-log issues).  We simply take
# the maximum DVI among all rows with no_soil==1 over the target months.
compute_ppi_m_from_populus_q90 <- function(df, months = SUMMER_DETREND_MONTHS) {
  # NOTE: despite the name, we now ignore Veg and simply take the maximum DVI
  # among all rows where no_soil == 1 within the requested months.
  if (is.null(df) || nrow(df) == 0) return(NA_real_)
  if (!"DVI" %in% names(df) && all(c("nir", "red") %in% names(df))) df$DVI <- as.numeric(df$nir) - as.numeric(df$red)
  if (!"DVI" %in% names(df)) return(NA_real_)
  if (!"no_soil" %in% names(df)) return(NA_real_)

  no_soil_int <- df$no_soil
  if (is.logical(no_soil_int)) {
    no_soil_int <- as.integer(no_soil_int)
  } else {
    no_soil_int <- suppressWarnings(as.integer(as.character(no_soil_int)))
  }

  veg_norm <- tolower(trimws(as.character(df$Veg)))
  month_mask <- rep(TRUE, nrow(df))
  if (!is.null(months) && length(months) > 0 && ("date" %in% names(df) || "prediction_date" %in% names(df))) {
    dcol <- if ("date" %in% names(df)) "date" else "prediction_date"
    dd <- suppressWarnings(as.Date(df[[dcol]]))
    month_mask <- lubridate::month(dd) %in% months
  }

  keep <- month_mask & is.finite(df$DVI) & (no_soil_int == 1)
  vals <- as.numeric(df$DVI[keep])
  vals <- vals[is.finite(vals)]
  if (length(vals) < 5L) return(NA_real_)

  # M = observed peak DVI (no additional scaling)
  dvi_max <- max(vals)
  if (!is.finite(dvi_max) || dvi_max <= 0) return(NA_real_)
  as.numeric(dvi_max)
}

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
if (!exists("ENABLE_OUTLIER_REMOVAL", inherits = TRUE)) ENABLE_OUTLIER_REMOVAL <- TRUE
if (!exists("OUTLIER_MAD_THRESHOLD", inherits = TRUE)) OUTLIER_MAD_THRESHOLD <- 3.5

# global configuration file (mesma_config.R) defines many constants such as
# OUTLIER_SPLINE_MAX_DF; source it if available so users can override defaults
if (file.exists("mesma_config.R")) {
  source("mesma_config.R")
}

# maximum degrees-of-freedom allowed when fitting splines for outlier detection
# and seasonal normalization.  The value in mesma_config.R takes precedence.
if (!exists("OUTLIER_SPLINE_MAX_DF", inherits = TRUE)) OUTLIER_SPLINE_MAX_DF <- 10L

# Ensure an EVI helper exists before any index logic
ensure_evi <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  # create column if missing
  if (!"EVI" %in% names(df)) df$EVI <- NA_real_
  # if all values NA and nir/red exist, attempt to compute
  if (all(!is.finite(df$EVI)) && all(c("nir","red") %in% names(df))) {
    if ("blue" %in% names(df)) {
      df$EVI <- 2.5 * (as.numeric(df$nir) - as.numeric(df$red)) /
                (as.numeric(df$nir) + 6*as.numeric(df$red) - 7.5*as.numeric(df$blue) + 1)
    } else {
      warning("Blue band missing - approximating EVI with red for blue values")
      df$EVI <- 2.5 * (as.numeric(df$nir) - as.numeric(df$red)) /
                (as.numeric(df$nir) + 6*as.numeric(df$red) - 7.5*as.numeric(df$red) + 1)
    }
  }
  df
}

# Remove large outliers robustly per (location_id, pheno_year).
# Matches fit_veg_mixture_mesma.R behavior: spline-only outlier detection for groups with
# sufficient data (>=10 obs with dates), and drop location-years with <5 observations.
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
    # Match fit script: drop location-years with fewer than 5 observations entirely.
    if (length(rows) < 5) {
      removed_idx[rows] <- TRUE
      next
    }
    out_mask <- rep(FALSE, nrow(sub))

    # Match fit script: only run spline outlier removal if date exists and we have
    # at least 10 observations; otherwise skip outlier detection for this group.
    has_date <- "date" %in% names(sub) && any(!is.na(sub$date))
    if (!has_date || length(rows) < 10) next

    # Compute DOY
    sub$doy <- as.numeric(format(sub$date, "%j"))
    for (col in candidates) {
      if (!is.numeric(sub[[col]])) next
      colv <- sub[[col]]
      finite_idx <- is.finite(colv) & is.finite(sub$doy)
      if (sum(finite_idx) < 5) next

      tryCatch({
        # Pass 1: initial spline fit
        x <- sub$doy[finite_idx]
        y <- colv[finite_idx]

        n_unique <- length(unique(x))
        fit1 <- stats::smooth.spline(x, y, df = min(OUTLIER_SPLINE_MAX_DF, length(x)/2, n_unique - 1))
        pred1 <- predict(fit1, x)$y
        res1 <- y - pred1
        mad1 <- stats::mad(res1, na.rm = TRUE)
        if (!is.finite(mad1) || mad1 <= 1e-6) next

        keep_mask <- abs(res1 - stats::median(res1, na.rm = TRUE)) <= (mad_thresh * 1.5 * mad1)

        # Pass 2: refit on cleaner points if enough remain
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
        if (!is.finite(mad_res) || mad_res <= 0) next

        this_mask <- rep(FALSE, length(colv))
        this_mask[finite_idx] <- abs(residuals - med_res) > mad_thresh * mad_res
        out_mask <- out_mask | this_mask
      }, error = function(e) {
        # Match fit script: skip this column if spline fitting fails.
      })
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
# 0. Load Real Data
# ==============================================================================

cat("Loading real phenology data...\n")

# Try to load from common data file locations
data_file <- NULL
possible_paths <- c(
  # labelled training splits (preferred)
  "C:\\Users\\yolan\\Downloads\\Landsat_Harmonized_Bands_1985_2025_train (2).csv",
  "C:\\Users\\yolan\\Downloads\\Landsat_Harmonized_Bands_1985_2025_train.csv"
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

# Reconstruct location_id from coordinates (same policy as fit_veg_mixture_mesma.R):
# do not trust existing location_id values from CSV.
if ("location_id" %in% names(df_raw)) {
  df_raw$location_id_orig <- df_raw$location_id
  df_raw$location_id <- NULL
  cat("[NOTICE] Removed existing 'location_id' from input; reconstructing from coordinates\n")
}
if (all(c("lon", "lat") %in% names(df_raw))) {
  df_raw$location_id <- make_location_id(df_raw$lon, df_raw$lat)
} else if (all(c("target_lon", "target_lat") %in% names(df_raw))) {
  df_raw$location_id <- make_location_id(df_raw$target_lon, df_raw$target_lat)
} else {
  stop("[INDEX SETUP] Cannot reconstruct location_id: missing lon/lat (or target_lon/target_lat)")
}

if ("vegetation" %in% names(df_raw) && !"Veg" %in% names(df_raw)) {
  df_raw$Veg <- df_raw$vegetation
  cat("[NOTICE] Renamed 'vegetation' -> 'Veg' in phenology data\n")
}
# If Veg already present, be explicit about skipping GeoJSON join
if ("Veg" %in% names(df_raw)) {
  df_raw$Veg <- tolower(trimws(as.character(df_raw$Veg)))
  df_raw$Veg[df_raw$Veg %in% c("", "na", "null")] <- NA_character_
  n_non_missing_veg <- sum(!is.na(df_raw$Veg))
  n_barren_veg <- sum(df_raw$Veg == "barren", na.rm = TRUE)
  cat(sprintf("[VEG DIAG] Non-missing Veg rows: %d/%d; barren rows: %d\n", n_non_missing_veg, nrow(df_raw), n_barren_veg))
  if (n_non_missing_veg > 0) {
    veg_tab <- sort(table(df_raw$Veg), decreasing = TRUE)
    top_n <- min(10L, length(veg_tab))
    cat(sprintf("[VEG DIAG] Top Veg labels: %s\n", paste(sprintf("%s=%d", names(veg_tab)[seq_len(top_n)], as.integer(veg_tab[seq_len(top_n)])), collapse = ", ")))
  } else {
    cat("[VEG DIAG] Veg column exists but is entirely empty/NA in the input CSV used for endmember extraction. Soil-line estimation from barren rows cannot run at this stage.\n")
    # try loading a labelled training split if it's available and different from the current file
    alt_path <- "C:\\Users\\yolan\\Downloads\\Landsat_Harmonized_Bands_1985_2025_train (2).csv"
    if (file.exists(alt_path) && alt_path != data_file) {
      cat(sprintf("[NOTICE] Reloading df_raw from alternate labelled file: %s\n", alt_path))
      df_raw <- fread(alt_path)
      df_raw <- normalize_band_names(df_raw)
      if ("vegetation" %in% names(df_raw) && !"Veg" %in% names(df_raw)) df_raw$Veg <- df_raw$vegetation
      df_raw$Veg <- tolower(trimws(as.character(df_raw$Veg)))
      df_raw$Veg[df_raw$Veg %in% c("", "na", "null")] <- NA_character_
      n_non_missing_veg <- sum(!is.na(df_raw$Veg))
      n_barren_veg <- sum(df_raw$Veg == "barren", na.rm = TRUE)
      cat(sprintf("[VEG DIAG] After reload: Non-missing Veg rows: %d/%d; barren rows: %d\n", n_non_missing_veg, nrow(df_raw), n_barren_veg))
    }
  }
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
  # ensure EVI exists for this preliminary index set
  df_raw <- ensure_evi(df_raw)
  cat("[INDEX SETUP] ensured EVI in df_raw after initial compute_indices_from_bands\n")
}

# Helper: assign phenology year from a date (used by outlier grouping)
assign_pheno_year <- function(d) {
  d <- as.Date(d)
  ifelse(is.na(d), NA_integer_, ifelse(lubridate::month(d) >= 3, lubridate::year(d), lubridate::year(d) - 1))
}

# Early dust filtering (NDDI) to remove contaminated observations before processing
eps <- 1e-9
if (all(c('red','nir') %in% names(df_raw))) df_raw$NDDI <- (as.numeric(df_raw$red) - as.numeric(df_raw$nir)) / (as.numeric(df_raw$red) + as.numeric(df_raw$nir) + eps)

if ("NDDI" %in% names(df_raw)) {
  dust_count <- sum(df_raw$NDDI > NDDI_DUST_THRESHOLD, na.rm = TRUE)
  total_before <- nrow(df_raw)
  df_raw <- df_raw[!(df_raw$NDDI > NDDI_DUST_THRESHOLD), , drop = FALSE]
  cat(sprintf("[FILTERING] Filtered out %d observations with dust contamination (NDDI > %s)\n", total_before - nrow(df_raw), .nddi_thresh_fmt()))
}

# Remove large outliers (if enabled) after any early dust filtering,
# and BEFORE computing per-location Jan-Mar mean dvi_soil / PPI / separability.
df_raw <- remove_large_outliers(df_raw)

# Compute a per-location DVI soil baseline for PPI using January–March mean.
# The soil value for each location is now the average of all DVI observations
# occurring in months 1–3 (calendar year), ignoring any quantile-based logic.
# `min_samples` still guards against extremely small sample sizes; locations
# that do not meet this requirement are left as NA and should be skipped by caller.
compute_dvi_soil_per_location <- function(df, min_samples = 1L, months = PPI_SOIL_BASELINE_MONTHS) {
  if (is.null(df) || nrow(df) == 0) stop("[PPI] compute_dvi_soil_per_location: empty df")
  if (!"location_id" %in% names(df)) stop("[PPI] compute_dvi_soil_per_location: missing location_id")
  if (!"DVI" %in% names(df) && all(c("nir", "red") %in% names(df))) df$DVI <- as.numeric(df$nir) - as.numeric(df$red)
  if (!"DVI" %in% names(df)) stop("[PPI] compute_dvi_soil_per_location: missing DVI (and nir/red not available)")

  dcol <- NULL
  if (!is.null(months) && length(months) > 0 && ("date" %in% names(df) || "prediction_date" %in% names(df))) {
    dcol <- if ("date" %in% names(df)) "date" else "prediction_date"
  }

  locs <- unique(as.character(df$location_id))
  dvi_soil_vec <- rep(NA_real_, nrow(df))
  for (loc in locs) {
    idx <- which(as.character(df$location_id) == loc)
    vals <- df$DVI[idx]
    # restrict to requested months if date info is available
    if (!is.null(dcol)) {
      dd <- suppressWarnings(as.Date(df[[dcol]][idx]))
      msk <- lubridate::month(dd) %in% months
      vals <- vals[msk]
    }
    vals <- vals[is.finite(vals)]
    if (length(vals) < as.integer(min_samples)) next
    soil <- suppressWarnings(as.numeric(mean(vals, na.rm = TRUE)))
    if (!is.finite(soil)) next
    dvi_soil_vec[idx] <- soil
  }
  dvi_soil_vec
}

# Compute PPI using per-location Jan-Mar mean DVI baseline (no fallback).
if (!"location_id" %in% names(df_raw)) stop("[INDEX SETUP] Cannot compute per-location PPI: missing location_id")
if (!"DVI" %in% names(df_raw) && all(c("nir", "red") %in% names(df_raw))) df_raw$DVI <- as.numeric(df_raw$nir) - as.numeric(df_raw$red)
if (!"DVI" %in% names(df_raw) || !any(is.finite(df_raw$DVI))) {
  stop("[INDEX SETUP] Cannot compute per-location PPI because DVI is missing or invalid")
}
dvi_soil_vec <- compute_dvi_soil_per_location(df_raw, min_samples = 1L)
insufficient_idx <- is.finite(df_raw$DVI) & !is.finite(dvi_soil_vec)
if (any(insufficient_idx)) {
  insufficient_locs <- unique(as.character(df_raw$location_id[insufficient_idx]))
  keep_idx <- !(as.character(df_raw$location_id) %in% insufficient_locs)
  cat(sprintf("[PPI] Skipping %d locations with insufficient Jan-Mar rows for dvi_soil; removed %d rows (example locs: %s)\n",
              length(insufficient_locs), sum(!keep_idx), paste(head(insufficient_locs, 10), collapse = ", ")))
  df_raw <- df_raw[keep_idx, , drop = FALSE]
  dvi_soil_vec <- dvi_soil_vec[keep_idx]
}
if (!any(is.finite(dvi_soil_vec))) {
  stop("[PPI] No locations with sufficient Jan-Mar rows remained after skipping")
}
cat(sprintf("[PPI] per-location dvi_soil(Jan-Mar mean) summary: finite=%d/%d, min=%.6f, median=%.6f, max=%.6f\n",
            sum(is.finite(dvi_soil_vec)), length(dvi_soil_vec),
            min(dvi_soil_vec, na.rm = TRUE), stats::median(dvi_soil_vec, na.rm = TRUE), max(dvi_soil_vec, na.rm = TRUE)))

# --- M diagnostics: understand exactly what rows feed into M calculation ---
{
  .dv  <- tolower(trimws(as.character(df_raw$Veg)))
  .ns_raw <- df_raw$no_soil
  .ns  <- suppressWarnings(as.integer(as.character(.ns_raw)))
  .dvi <- if ("DVI" %in% names(df_raw)) as.numeric(df_raw$DVI) else rep(NA_real_, nrow(df_raw))
  .dd  <- suppressWarnings(as.Date(df_raw$date))
  .mon <- lubridate::month(.dd)

  cat("[PPI M DIAG] no_soil column class  :", class(df_raw$no_soil), "\n")
  cat("[PPI M DIAG] no_soil unique values  :", paste(sort(unique(.ns_raw)), collapse = ", "), "\n")
  cat("[PPI M DIAG] no_soil as int uniq    :", paste(sort(unique(.ns)), collapse = ", "), "\n")
  cat("[PPI M DIAG] Veg unique vals        :", paste(sort(unique(.dv[!is.na(.dv)])), collapse = ", "), "\n")

  .n_rows    <- nrow(df)
  .n_no_soil <- sum(.ns == 1L, na.rm = TRUE)
  .n_no_soil_sum <- sum(.ns == 1L & .mon %in% SUMMER_DETREND_MONTHS, na.rm = TRUE)
  .n_no_soil_fin <- sum(.ns == 1L & .mon %in% SUMMER_DETREND_MONTHS & is.finite(.dvi), na.rm = TRUE)

  cat(sprintf("[PPI M DIAG] total rows              : %d\n", .n_rows))
  cat(sprintf("[PPI M DIAG] rows with no_soil==1     : %d\n", .n_no_soil))
  cat(sprintf("[PPI M DIAG] no_soil==1 + months %s : %d\n", paste(SUMMER_DETREND_MONTHS, collapse=","), .n_no_soil_sum))
  cat(sprintf("[PPI M DIAG] ... + finite DVI        : %d\n", .n_no_soil_fin))

  if (.n_no_soil_fin > 0) {
    .vals <- .dvi[.ns == 1L & .mon %in% SUMMER_DETREND_MONTHS & is.finite(.dvi)]
    # use na.rm=TRUE to avoid crashing when unexpected NA/NaN slip through
    cat(sprintf("[PPI M DIAG] DVI of qualifying rows: min=%.4f  Q10=%.4f  Q50=%.4f  Q90=%.4f  Q95=%.4f  max=%.4f  -> M=max=%.4f\n",
                min(.vals, na.rm = TRUE), quantile(.vals, .10, na.rm = TRUE), quantile(.vals, .50, na.rm = TRUE),
                quantile(.vals, .90, na.rm = TRUE), quantile(.vals, .95, na.rm = TRUE),
                max(.vals, na.rm = TRUE), max(.vals, na.rm = TRUE)))
  }
  rm(.dv, .ns_raw, .ns, .dvi, .dd, .mon,
     .n_pop, .n_pop_ns1, .n_pop_ns0, .n_pop_nsNA, .n_pop_sum, .n_pop_fin)
}
# ---- end M diagnostics ----

ppi_m <- compute_ppi_m_from_populus_q90(df_raw, months = SUMMER_DETREND_MONTHS)
if (!is.finite(ppi_m)) {
  # Emit detailed diagnostics before aborting so the user knows exactly what is missing.
  .diag_veg  <- tolower(trimws(as.character(df_raw$Veg)))
  .diag_ns   <- suppressWarnings(as.integer(as.character(df_raw$no_soil)))
  .diag_dvi  <- if ("DVI" %in% names(df_raw)) as.numeric(df_raw$DVI) else rep(NA_real_, nrow(df_raw))
  .diag_date <- if ("date" %in% names(df_raw)) suppressWarnings(as.Date(df_raw$date)) else as.Date(rep(NA_character_, nrow(df_raw)))
  .diag_mon  <- lubridate::month(.diag_date)
  .n_rows     <- nrow(df_raw)
  .n_no_soil  <- sum(.diag_ns == 1L, na.rm = TRUE)
  .n_no_soil_sum <- sum(.diag_ns == 1L & .diag_mon %in% SUMMER_DETREND_MONTHS, na.rm = TRUE)
  .n_no_soil_fin <- sum(.diag_ns == 1L & .diag_mon %in% SUMMER_DETREND_MONTHS & is.finite(.diag_dvi), na.rm = TRUE)
  stop(sprintf(
    paste0("[PPI] Could not compute no_soil==1 summer DVI max for M on df_raw from months %s (need >=5 finite rows).\n",
           "  total rows                   : %d\n",
           "  rows with no_soil==1         : %d\n",
           "  ... and in summer months (%s): %d\n",
           "  ... with finite DVI           : %d\n",
           "Hard requirement violated; aborting."),
    paste(SUMMER_DETREND_MONTHS, collapse = ","),
    .n_rows, .n_no_soil,
    paste(SUMMER_DETREND_MONTHS, collapse = ","),
    .n_no_soil_sum, .n_no_soil_fin
  ))
}
PPI_DVI_MAX_DEFAULT <- ppi_m
df_raw <- add_ppi_columns(df_raw, dvi_soil = dvi_soil_vec)
cat(sprintf("[INDEX SETUP] PPI computed/recomputed using per-location Jan-Mar mean DVI baseline with M=%.6f (no_soil==1 Jul-Sep max DVI); soil median=%.6f\n", ppi_m, stats::median(dvi_soil_vec, na.rm=TRUE)))

# Persist per-location Jan-Mar mean dvi_soil lookup computed from full non-detrended raw data.
# This lookup must be reused for all downstream subsets to avoid recalculating on filtered/detrended slices.
DVI_SOIL_JFM_MEAN_BY_LOCATION <- tapply(
  dvi_soil_vec,
  as.character(df_raw$location_id),
  function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) NA_real_ else as.numeric(x[1])
  }
)

map_dvi_soil_from_full_baseline <- function(df, lookup, context = "data") {
  if (is.null(df) || nrow(df) == 0) stop(sprintf("[PPI] %s: empty df", context))
  if (!"location_id" %in% names(df)) stop(sprintf("[PPI] %s: missing location_id", context))
  if (!"DVI" %in% names(df) && all(c("nir", "red") %in% names(df))) df$DVI <- as.numeric(df$nir) - as.numeric(df$red)
  if (!"DVI" %in% names(df)) stop(sprintf("[PPI] %s: missing DVI (and nir/red not available)", context))

  mapped <- as.numeric(lookup[as.character(df$location_id)])
  need_idx <- is.finite(df$DVI) & !is.finite(mapped)
  if (any(need_idx)) {
    bad_locs <- unique(as.character(df$location_id[need_idx]))
    stop(sprintf("[PPI] %s: missing full-data Jan-Mar mean dvi_soil baseline for %d rows across %d locations (example locs: %s)",
                 context, sum(need_idx), length(bad_locs), paste(head(bad_locs, 10), collapse = ", ")))
  }
  mapped
}

# Instead of fitting seasonal models, use July-September medians per user request
if (!"date" %in% names(df_raw)) stop("Phenology data must include a 'date' column")
if (!lubridate::is.Date(df_raw$date)) df_raw$date <- as.Date(df_raw$date)
cat(sprintf("[INDEX DETREND] Skipping seasonal detrending; using months %s (July-Sep) median per location instead\n", paste(SUMMER_DETREND_MONTHS, collapse = ",")))

# Keep original index to assign back robustly (create before filtering)
df_raw$.orig_row <- seq_len(nrow(df_raw))

# Prepare masks and containers
summer_mask <- lubridate::month(df_raw$date) %in% SUMMER_DETREND_MONTHS
indices_available <- intersect(INDICES_OF_INTEREST, names(df_raw))
if (length(indices_available) == 0) stop("[INDEX DETREND] No indices available to compute Jul-Sep medians; ensure indices are computed on raw data.")

INDEX_SUMMER_MEDIANS <- list()

# Create "_median" columns but populate them with raw index values only for Jul-Sep rows
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

cat(sprintf("[INDEX DETREND] Created %d '_median' columns populated with Jul-Sep values\n", length(indices_available)))

# Report location_id diagnostics if provided in CSV
if ("location_id" %in% names(df_raw)) {
  unique_loc <- length(unique(df_raw$location_id))
  cat(sprintf("DEBUG: Unique location_id values in CSV: %d\n", unique_loc))
  first5 <- paste(head(unique(df_raw$location_id), 5), collapse = ", ")
  cat(sprintf("DEBUG: First 5 location_id values in CSV: %s \n", first5))
  cat("Using reconstructed 'location_id' from coordinates.\n")
}
# Legacy 'no.soil'/'no_soil' columns are ignored by most of the script,
# but a "no_soil" indicator **is** honoured when computing separability
# (vegetation endmember will be drawn only from rows with no_soil == 1).
# This keeps the earlier behaviour for everything else while satisfying
# the new requirement.

# Require that Veg be present to extract endmembers
if (!("Veg" %in% names(df_raw))) {
  stop("Phenology data must include a 'Veg' column. GeoJSON support has been removed — include this field in landsat_lower.csv.")
}

cat("[NOTICE] Pure endmember delineation remains simplified; separability diagnostics are enabled.\n")

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
    # Use available soil baseline if present.  When the dataframe includes a
    # "fraction_veg" column we treat rows with fraction_veg==0 as pure soil
    # endmember and derive dvi_soil from them.  Otherwise no PPI is computed.
    dvi_soil_val <- NA_real_
    if ("fraction_veg" %in% names(df)) {
      dvi_soil_val <- df$DVI[df$fraction_veg == 0][1]
    }
    if (is.finite(dvi_soil_val)) {
      # Calculate zenith angle (typical mid‑latitude mid‑summer value)
      zenith_rad <- calculate_solar_zenith(lat = 40, doy = 180, hour = 10.5)
      M_val <- 0.7
      df$PPI <- ppi(dvi = df$DVI, zenith.angle = zenith_rad, M = M_val, dvi.soil = dvi_soil_val)
      cat(sprintf("Calculated PPI using soil baseline (dvi_soil=%.6f, zenith=%.4f rad)\n", dvi_soil_val, zenith_rad))
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
# 3. Calculate Separability for Pure Populus vs Soil (No synthetic mixing)
# ==============================================================================

# -----------------------------------------------------------------------------
# Compute PPI first (summer medians used as representative values) so we have one stable index value per loc-year
# -----------------------------------------------------------------------------
# Compute indices (DVI etc.) on raw data if not present
if (!"DVI" %in% names(df_raw) || !all(c("nir","red") %in% names(df_raw))) {
  df_raw <- compute_indices_from_bands(df_raw)
  df_raw <- ensure_evi(df_raw)
  cat("[INDEX SETUP] ensured EVI in df_raw after second compute_indices_from_bands call\n")
}

# PPI must already exist from the strict per-location Jan-Mar mean pipeline above.
if (!"PPI" %in% names(df_raw) || all(!is.finite(df_raw$PPI))) {
  stop("[INDEX SETUP] PPI missing after per-location Jan-Mar mean dvi_soil computation")
}

# Detrend summer indices (June-Sep) using polynomial fit of DOY for each index
if (!"date" %in% names(df_raw)) stop("Phenology data must include a 'date' column for detrending")
if (!lubridate::is.Date(df_raw$date)) df_raw$date <- as.Date(df_raw$date)

# Reset .orig_row to match the current filtered df_raw row indices
df_raw$.orig_row <- seq_len(nrow(df_raw))

cat(sprintf("[INDEX DETREND] Using months %s for summer median representative (no seasonal detrending)\n", paste(SUMMER_DETREND_MONTHS, collapse = ",")))

# Ensure median columns exist (created earlier). If not, create them now using Jul-Sep medians
if (!all(paste0(INDICES_OF_INTEREST, "_median") %in% names(df_raw))) {
  cat("[INDEX DETREND] Median columns not found — creating Jul-Sep '_median' columns now\n")
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
  cat("[INDEX DETREND] Median columns present — proceeding with Jul-Sep medians\n")
  if (!exists("INDEX_SUMMER_MEDIANS", envir = globalenv())) {
    idx_means <- lapply(intersect(INDICES_OF_INTEREST, names(df_raw)), function(idx) median(df_raw[[idx]][lubridate::month(df_raw$date) %in% SUMMER_DETREND_MONTHS], na.rm = TRUE))
    names(idx_means) <- intersect(INDICES_OF_INTEREST, names(df_raw))
    assign("INDEX_SUMMER_MEDIANS", idx_means, envir = globalenv())
  }
}

cat("\n=== CALCULATING PURE ENDMEMBER SEPARABILITY (POPULUS VS SOIL) ===\n")

# Build Jul-Sep representative rows per location and vegetation class.
# IMPORTANT: separability is computed from non-detrended summer values
# (i.e. the raw index columns), NOT from the *_median helper columns.
sep_index_cols <- intersect(INDICES_OF_INTEREST, names(df_raw))
summer_rows_sep <- df_raw |>
  dplyr::filter(lubridate::month(date) %in% SUMMER_DETREND_MONTHS) |>
  dplyr::filter(dplyr::if_any(dplyr::all_of(sep_index_cols), is.finite))

# Normalize no_soil to a simple 0/1 integer when present.
# (Many CSVs store it as factor/character; comparing those to 1 yields NA.)
if ("no_soil" %in% names(summer_rows_sep)) {
  if (is.logical(summer_rows_sep$no_soil)) {
    summer_rows_sep$no_soil <- as.integer(summer_rows_sep$no_soil)
  } else {
    summer_rows_sep$no_soil <- suppressWarnings(as.integer(as.character(summer_rows_sep$no_soil)))
  }
}

if (nrow(summer_rows_sep) == 0 || length(sep_index_cols) == 0) {
  warning("[SEPARABILITY] No Jul-Sep rows with finite index values available; skipping separability output.")
} else {
  # IMPORTANT: filter to the desired endmember populations *before* aggregating.
  # Otherwise a single pure (no_soil==1) observation could get mixed with
  # non-pure observations when we take a group median.
  soil_rows_raw <- summer_rows_sep |>
    dplyr::filter(tolower(Veg) == "barren")

  has_no_soil <- "no_soil" %in% names(summer_rows_sep)
  if (has_no_soil) {
    cat("[SEPARABILITY] restricting vegetation endmembers to rows with no_soil==1\n")
    veg_rows_raw <- summer_rows_sep |>
      dplyr::filter(tolower(Veg) != "barren" & !is.na(Veg) & no_soil == 1)
  } else {
    veg_rows_raw <- summer_rows_sep |>
      dplyr::filter(tolower(Veg) != "barren" & !is.na(Veg))
  }

  soil_rows <- soil_rows_raw |>
    dplyr::group_by(location_id, Veg) |>
    dplyr::summarise(
      n = dplyr::n(),
      dplyr::across(dplyr::all_of(sep_index_cols), ~ median(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  veg_rows <- veg_rows_raw |>
    dplyr::group_by(location_id, Veg) |>
    dplyr::summarise(
      n = dplyr::n(),
      dplyr::across(dplyr::all_of(sep_index_cols), ~ median(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  if (nrow(soil_rows) == 0 || nrow(veg_rows) == 0) {
    warning(sprintf("[SEPARABILITY] Missing required classes for separability (veg rows=%d, barren rows=%d); skipping output.",
                    nrow(veg_rows), nrow(soil_rows)))
  } else {
    separability_results <- data.frame(
      Index = character(),
      Separability = numeric(),
      MeanVeg = numeric(),
      MeanSoil = numeric(),
      SdPooled = numeric(),
      Nveg = integer(),
      Nsoil = integer(),
      stringsAsFactors = FALSE
    )

    for (idx_col in sep_index_cols) {
      veg_vals <- veg_rows[[idx_col]]
      soil_vals <- soil_rows[[idx_col]]
      veg_vals <- veg_vals[is.finite(veg_vals)]
      soil_vals <- soil_vals[is.finite(soil_vals)]

      if (length(veg_vals) > 0 && length(soil_vals) > 0) {
        n_veg <- length(veg_vals)
        n_soil <- length(soil_vals)
        mean_veg <- mean(veg_vals)
        mean_soil <- mean(soil_vals)
        var_veg <- var(veg_vals)
        var_soil <- var(soil_vals)
        if (is.na(var_veg)) var_veg <- 0
        if (is.na(var_soil)) var_soil <- 0
        denom <- sqrt((var_veg + var_soil) / 2)
        separability <- if (is.finite(denom) && denom > 0) (mean_veg - mean_soil) / denom else NA_real_

        separability_results <- rbind(
          separability_results,
          data.frame(
            Index = idx_col,
            Separability = separability,
            MeanVeg = mean_veg,
            MeanSoil = mean_soil,
            SdPooled = denom,
            Nveg = as.integer(n_veg),
            Nsoil = as.integer(n_soil),
            stringsAsFactors = FALSE
          )
        )
      }
    }

    if (nrow(separability_results) == 0) {
      warning("[SEPARABILITY] No valid separability entries could be computed.")
    } else {
      print(separability_results)
      sep_file <- file.path(output_dir, "separability_populus_vs_soil.csv")
      write.csv(separability_results, sep_file, row.names = FALSE)
      cat(sprintf("Saved separability results to: %s\n", sep_file))
    }
  }
}

# Synthetic mixing removed by request. Keep lightweight compatibility objects for downstream bias/trend code paths.
FVC_MODELS <- list()
assign("FVC_CALIBRATION_MODELS", FVC_MODELS, envir = globalenv())

estimate_fvc_from_index <- function(df, index_name, model_list = FVC_MODELS) {
  if (!(index_name %in% names(df))) return(rep(NA_real_, nrow(df)))
  if (!(index_name %in% names(model_list))) return(rep(NA_real_, nrow(df)))

  model_entry <- model_list[[index_name]]
  pred_data <- data.frame(index_val = df[[index_name]])

  if (!is.list(model_entry) && (inherits(model_entry, "rlm") || inherits(model_entry, "lm"))) {
    res <- tryCatch(predict(model_entry, newdata = pred_data), error = function(e) rep(NA_real_, nrow(df)))
    return(pmax(0, pmin(1, res)))
  }

  if (is.list(model_entry) && !is.null(model_entry$model)) {
    if (isTRUE(model_entry$use_menten) && !is.null(model_entry$transform$a) && !is.null(model_entry$transform$b)) {
      a <- model_entry$transform$a
      b <- model_entry$transform$b
      pred_data$index_val <- (a + pred_data$index_val) / (b + pred_data$index_val)
    }
    res <- tryCatch(predict(model_entry$model, newdata = pred_data), error = function(e) rep(NA_real_, nrow(df)))
    res[!is.finite(res)] <- NA_real_
    return(pmax(0, pmin(1, res)))
  }

  rep(NA_real_, nrow(df))
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

# Define input file path
# Reuse the already selected phenology source when available.
INPUT_CSV <- if (exists("data_file") && !is.null(data_file)) {
  data_file
} else {
  "C:\\Users\\yolan\\Downloads\\Landsat_Harmonized_Bands_1985_2025_train (2).csv"
}

# Infer canonical inference scope from filename/tag (MESMA-style low/mid/kon naming)
infer_inference_scope <- function(x) {
  if (is.null(x) || length(x) == 0) return(rep("unknown", length(x)))
  x_chr <- tolower(as.character(x))
  dplyr::case_when(
    grepl("(^|[^a-z0-9])(low|lower)([^a-z0-9]|$)", x_chr) ~ "low",
    grepl("(^|[^a-z0-9])mid([^a-z0-9]|$)", x_chr) ~ "mid",
    grepl("(^|[^a-z0-9])kon([^a-z0-9]|$)", x_chr) ~ "kon",
    TRUE ~ "unknown"
  )
}

# Auto-detect companion low/mid/kon inference files from the selected input file
infer_mid_low_kon_sequence <- function(path) {
  if (is.null(path) || !nzchar(as.character(path))) return(character(0))
  p <- normalizePath(as.character(path), winslash = "/", mustWork = FALSE)
  src_name <- basename(p)
  src_low <- tolower(src_name)

  make_variant <- function(token) {
    sub("(^|[^A-Za-z0-9])mid([^A-Za-z0-9]|$)", paste0("\\1", token, "\\2"), src_name, ignore.case = TRUE)
  }

  if (grepl("(^|[^a-z0-9])mid([^a-z0-9]|$)", src_low)) {
    cands <- c(
      p,
      file.path(dirname(p), make_variant("low")),
      file.path(dirname(p), make_variant("kon"))
    )
  } else {
    cands <- c(p)
  }

  cands <- unique(cands[file.exists(cands)])
  cands
}

# Load the phenology data
# Prefer explicit inference files when present (kon/mid/low), otherwise auto-detect from INPUT_CSV.
preferred_inference_csvs <- c(
  "C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_kon1.csv",
  "C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_mid (2).csv",
  "C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_low (3).csv"
)
existing_preferred_inference_csvs <- preferred_inference_csvs[file.exists(preferred_inference_csvs)]

if (length(existing_preferred_inference_csvs) > 0) {
  input_csvs <- unique(existing_preferred_inference_csvs)
  cat(sprintf("[INPUT] Using explicit inference CSV set (%d file(s))\n", length(input_csvs)))
} else {
  input_csvs <- infer_mid_low_kon_sequence(INPUT_CSV)
  if (length(input_csvs) == 0) input_csvs <- INPUT_CSV
}

if (length(input_csvs) > 1) {
  cat(sprintf("Loading %d inference sources for trend analysis (mid/low/kon sequence when available)\n", length(input_csvs)))
} else {
  cat("Loading data from:", input_csvs[1], "\n")
}

df_list <- lapply(input_csvs, function(csv_path) {
  cat("Loading data from:", csv_path, "\n")
  d <- readr::read_csv(csv_path, show_col_types = FALSE)
  d$.inference_source_csv <- as.character(csv_path)
  d$.inference_scope <- infer_inference_scope(basename(csv_path))
  d
})

df <- dplyr::bind_rows(df_list)
if (!".inference_scope" %in% names(df)) df$.inference_scope <- infer_inference_scope(basename(INPUT_CSV))

# Treat each inference file as one trend location (not each point location).
if (!".trend_location_id" %in% names(df)) {
  if (".inference_source_csv" %in% names(df)) {
    df$.trend_location_id <- tools::file_path_sans_ext(basename(as.character(df$.inference_source_csv)))
  } else {
    df$.trend_location_id <- as.character(df$.inference_scope)
  }
}
df$.trend_location_id <- as.character(df$.trend_location_id)
df$.trend_location_id[!nzchar(df$.trend_location_id) | is.na(df$.trend_location_id)] <- "unknown_source"

scope_counts <- sort(table(df$.inference_scope), decreasing = TRUE)
cat(sprintf("[INPUT] Inference scope counts: %s\n",
            paste(sprintf("%s=%d", names(scope_counts), as.integer(scope_counts)), collapse = ", ")))
trend_counts <- sort(table(df$.trend_location_id), decreasing = TRUE)
cat(sprintf("[INPUT] Trend units (one per inference file): %s\n",
            paste(sprintf("%s=%d", names(trend_counts), as.integer(trend_counts)), collapse = ", ")))

# Determine if input CSV is effectively training data; if so, skip plotting trends.
# Final value is computed AFTER TRAINING_CSV is defined.
SKIP_TRENDS <- FALSE
input_norm <- unique(tryCatch(normalizePath(input_csvs, winslash="/", mustWork=FALSE), error=function(e) input_csvs))

# optionally drop observations from years 1992-1999 (these are often
# problematic Landsat TM/ETM+ years that we historically excluded).  Users
# can disable this behaviour by setting EXCLUDE_PRE2000 <- FALSE before
# sourcing this script.
if (!exists("EXCLUDE_PRE2000", inherits = TRUE)) EXCLUDE_PRE2000 <- TRUE

if (isTRUE(EXCLUDE_PRE2000)) {
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
} else {
  cat("[DATA FILTER] EXCLUDE_PRE2000 is FALSE; retaining all pre‑2000 years\n")
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
  }
}

# Reconstruct location_id from coordinates (fit_veg_mixture_mesma.R policy).
# Do not use existing location_id values from CSV.
if ("location_id" %in% names(df)) {
  df$location_id_orig <- df$location_id
  df$location_id <- NULL
  cat("[NOTICE] Removed existing 'location_id' column - will reconstruct from lat/lon\n")
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

if (all(c("lon", "lat") %in% names(df))) {
  df$location_id <- make_location_id(df$lon, df$lat)
} else if (all(c("target_lon", "target_lat") %in% names(df))) {
  df$location_id <- make_location_id(df$target_lon, df$target_lat)
  df$lat <- df$target_lat
  df$lon <- df$target_lon
} else {
  stop("Cannot reconstruct location_id: neither lon/lat nor target_lon/target_lat are available")
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
# Always use explicit training CSV; no fallback allowed.
TRAINING_CSV <- "C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_train (2).csv"
if (!file.exists(TRAINING_CSV)) {
  stop(sprintf("[TRAINING] Required TRAINING_CSV not found: %s", TRAINING_CSV))
}
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
training_df <- ensure_evi(training_df)
cat(sprintf("[TRAINING] EVI finite values after calculate_indices: %d/%d\n",
            sum(is.finite(training_df$EVI)), nrow(training_df)))
  
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

  # Enforce per-location Jan-Mar mean PPI for training data (always recompute; no global fallback)
  if (!"location_id" %in% names(training_df)) stop("[TRAINING] Cannot compute PPI: missing location_id")
  if (!"DVI" %in% names(training_df) && all(c("nir", "red") %in% names(training_df))) training_df$DVI <- as.numeric(training_df$nir) - as.numeric(training_df$red)
  if (!"DVI" %in% names(training_df) || !any(is.finite(training_df$DVI))) {
    stop("[TRAINING] Cannot compute per-location Jan-Mar mean PPI: DVI missing or invalid")
  }
  # IMPORTANT: compute Jan-Mar mean baseline on full non-detrended training data.
  dvi_soil_vec <- compute_dvi_soil_per_location(training_df, min_samples = 1L, months = PPI_SOIL_BASELINE_MONTHS)
  insufficient_idx <- is.finite(training_df$DVI) & !is.finite(dvi_soil_vec)
  if (any(insufficient_idx)) {
    insufficient_locs <- unique(as.character(training_df$location_id[insufficient_idx]))
    keep_idx <- !(as.character(training_df$location_id) %in% insufficient_locs)
    cat(sprintf("[TRAINING PPI] Skipping %d locations with insufficient Jan-Mar rows for dvi_soil; removed %d rows (example locs: %s)\n",
                length(insufficient_locs), sum(!keep_idx), paste(head(insufficient_locs, 10), collapse = ", ")))
    training_df <- training_df[keep_idx, , drop = FALSE]
    dvi_soil_vec <- dvi_soil_vec[keep_idx]
  }
  if (!any(is.finite(dvi_soil_vec))) {
    stop("[TRAINING PPI] No locations with sufficient Jan-Mar rows remained after skipping")
  }
  ppi_m_train <- compute_ppi_m_from_populus_q90(training_df, months = SUMMER_DETREND_MONTHS)
  if (!is.finite(ppi_m_train)) {
    .diag_veg  <- tolower(trimws(as.character(training_df$Veg)))
    .diag_ns   <- suppressWarnings(as.integer(as.character(training_df$no_soil)))
    .diag_dvi  <- if ("DVI" %in% names(training_df)) as.numeric(training_df$DVI) else rep(NA_real_, nrow(training_df))
    .diag_date <- if ("date" %in% names(training_df)) suppressWarnings(as.Date(training_df$date)) else as.Date(rep(NA_character_, nrow(training_df)))
    .diag_mon  <- lubridate::month(.diag_date)
    .n_rows     <- nrow(training_df)
    .n_no_soil  <- sum(.diag_ns == 1L, na.rm = TRUE)
    .n_no_soil_sum <- sum(.diag_ns == 1L & .diag_mon %in% SUMMER_DETREND_MONTHS, na.rm = TRUE)
    .n_no_soil_fin <- sum(.diag_ns == 1L & .diag_mon %in% SUMMER_DETREND_MONTHS & is.finite(.diag_dvi), na.rm = TRUE)
    stop(sprintf(
      paste0("[PPI] Could not compute no_soil==1 summer DVI max for M on training_df from months %s (need >=5 finite rows).\n",
             "  total rows                   : %d\n",
             "  rows with no_soil==1         : %d\n",
             "  ... and in summer months (%s): %d\n",
             "  ... with finite DVI           : %d\n",
             "Hard requirement violated; aborting."),
      paste(SUMMER_DETREND_MONTHS, collapse = ","),
      .n_pop, .n_pop_ns,
      paste(SUMMER_DETREND_MONTHS, collapse = ","),
      .n_pop_sum, .n_pop_fin
    ))
  }
  PPI_DVI_MAX_DEFAULT <- ppi_m_train
  training_df <- add_ppi_columns(training_df, dvi_soil = dvi_soil_vec)
  cat(sprintf("[TRAINING] Applied per-location Jan-Mar mean dvi_soil baseline with M=%.6f (no_soil==1 Jul-Sep max DVI); soil median=%.6f\n", ppi_m_train, stats::median(dvi_soil_vec, na.rm=TRUE)))
  # Filter training years to the analysis window
  cat("Filtering training data for years 1985-2025...\n")
  training_df <- training_df |> dplyr::filter(year >= 1985 & year <= 2025)
  cat("Training rows after year filtering:", nrow(training_df), "\n")



# Finalize trend plotting skip decision after training fallback choice is known.
TRAINING_CSV_PATH <- tryCatch(normalizePath(TRAINING_CSV, winslash = "/", mustWork = FALSE), error = function(e) NA_character_)
SKIP_TRENDS <- any(!is.na(TRAINING_CSV_PATH) & input_norm == TRAINING_CSV_PATH)
if (isTRUE(SKIP_TRENDS)) {
  cat("[CONFIG] Detected training data as input; trend plotting will be skipped\n")
} else {
  cat("[CONFIG] Input treated as inference data; trend plotting enabled\n")
}

# Reconstruct location_id from coordinates; do not use existing IDs.
if ("location_id" %in% names(df)) {
  df$location_id_orig <- df$location_id
  df$location_id <- NULL
  cat("[NOTICE] Removed existing 'location_id' column - will reconstruct from lat/lon\n")
}
if (all(c("lon", "lat") %in% names(df))) {
  df$location_id <- make_location_id(df$lon, df$lat)
} else if (all(c("target_lon", "target_lat") %in% names(df))) {
  df$location_id <- make_location_id(df$target_lon, df$target_lat)
  df$lat <- df$target_lat
  df$lon <- df$target_lon
} else {
  stop("Cannot reconstruct location_id: neither lon/lat nor target_lon/target_lat are available")
}

if ("Veg" %in% names(df)) {
  df$Veg <- tolower(df$Veg)
  cat("DEBUG: Unique Veg types in df:\n")
  print(table(df$Veg, useNA = "ifany"))
}

# Ensure date and temporal columns
df$date <- as.Date(df$date)
df$year <- lubridate::year(df$date)
if (!"pheno_year" %in% names(df)) df$pheno_year <- assign_pheno_year(df$date)
df$month <- lubridate::month(df$date)

cat("Filtering data for years 1985-2025...\n")
df <- df |> dplyr::filter(year >= 1985 & year <= 2025)
cat("Data rows after year filtering:", nrow(df), "\n")

# Compute indices and strict per-location Jan-Mar mean PPI
df <- calculate_indices(df)
df <- ensure_evi(df)
cat(sprintf("[DATA] EVI finite values after calculate_indices: %d/%d\n",
            sum(is.finite(df$EVI)), nrow(df)))
if (!"DVI" %in% names(df) && all(c("nir", "red") %in% names(df))) {
  df$DVI <- as.numeric(df$nir) - as.numeric(df$red)
}
if (!"DVI" %in% names(df) || !any(is.finite(df$DVI))) {
  stop("[INDEX SETUP] Cannot compute per-location PPI for df because DVI is missing or invalid")
}
# IMPORTANT: compute Jan-Mar mean baseline on full non-detrended df.
dvi_soil_vec_df <- compute_dvi_soil_per_location(df, months = PPI_SOIL_BASELINE_MONTHS)
insufficient_idx_df <- is.finite(df$DVI) & !is.finite(dvi_soil_vec_df)
if (any(insufficient_idx_df)) {
  insufficient_locs_df <- unique(as.character(df$location_id[insufficient_idx_df]))
  keep_idx_df <- !(as.character(df$location_id) %in% insufficient_locs_df)
  cat(sprintf("[DATA PPI] Skipping %d locations with insufficient Jan-Mar rows for dvi_soil; removed %d rows (example locs: %s)\n",
              length(insufficient_locs_df), sum(!keep_idx_df), paste(head(insufficient_locs_df, 10), collapse = ", ")))
  df <- df[keep_idx_df, , drop = FALSE]
  dvi_soil_vec_df <- dvi_soil_vec_df[keep_idx_df]
}
if (!any(is.finite(dvi_soil_vec_df))) {
  stop("[DATA PPI] No locations with sufficient Jan-Mar rows remained after skipping")
}

# Persist lookup from the non-detrended all-month df used for trend analysis.
DVI_SOIL_JFM_MEAN_BY_LOCATION <- tapply(
  dvi_soil_vec_df,
  as.character(df$location_id),
  function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) NA_real_ else as.numeric(x[1])
  }
)
ppi_m_df <- compute_ppi_m_from_populus_q90(df, months = SUMMER_DETREND_MONTHS)
if (!is.finite(ppi_m_df)) {
  # collect same diagnostics used elsewhere to assist debugging
  .dv  <- tolower(trimws(as.character(df$Veg)))
  .ns  <- suppressWarnings(as.integer(as.character(df$no_soil)))
  .dvi <- if ("DVI" %in% names(df)) as.numeric(df$DVI) else rep(NA_real_, nrow(df))
  .dd  <- suppressWarnings(as.Date(df$date))
  .mon <- lubridate::month(.dd)
  .n_rows      <- nrow(df_raw)
  .n_no_soil   <- sum(.ns == 1L, na.rm = TRUE)
  .n_no_soil_sum <- sum(.ns == 1L & .mon %in% SUMMER_DETREND_MONTHS, na.rm = TRUE)
  .n_no_soil_fin <- sum(.ns == 1L & .mon %in% SUMMER_DETREND_MONTHS & is.finite(.dvi), na.rm = TRUE)
  if (is.finite(PPI_DVI_MAX_DEFAULT)) {
    warning(sprintf(
      paste0("[PPI] Could not compute no_soil==1 summer DVI max for M on df from months %s (need >=5 finite rows). ",
             "Falling back to training-derived M=%.6f.\n",
             "  total rows                   : %d\n",
             "  rows with no_soil==1         : %d\n",
             "  ... and in summer months (%s): %d\n",
             "  ... with finite DVI           : %d"),
      paste(SUMMER_DETREND_MONTHS, collapse = ","),
      PPI_DVI_MAX_DEFAULT,
      .n_rows, .n_no_soil,
      paste(SUMMER_DETREND_MONTHS, collapse = ","),
      .n_no_soil_sum, .n_no_soil_fin
    ))
    ppi_m_df <- PPI_DVI_MAX_DEFAULT
  } else {
    stop(sprintf(
      paste0("[PPI] Could not compute no_soil==1 summer DVI max for M on df from months %s (need >=5 finite rows), ",
             "and no valid training-derived fallback M is available.\n",
             "  total rows                   : %d\n",
             "  rows with no_soil==1         : %d\n",
             "  ... and in summer months (%s): %d\n",
             "  ... with finite DVI           : %d\n",
             "Hard requirement violated; aborting."),
      paste(SUMMER_DETREND_MONTHS, collapse = ","),
      .n_rows, .n_no_soil,
      paste(SUMMER_DETREND_MONTHS, collapse = ","),
      .n_no_soil_sum, .n_no_soil_fin
    ))
  }
}
PPI_DVI_MAX_DEFAULT <- ppi_m_df
df <- add_ppi_columns(df, dvi_soil = dvi_soil_vec_df)
cat(sprintf("[INDEX SETUP] Computed df PPI using per-location Jan-Mar mean DVI baseline from full non-detrended data (M=%.6f; no_soil==1 Jul-Sep max DVI); soil median=%.6f\n", ppi_m_df, stats::median(dvi_soil_vec_df, na.rm=TRUE)))



# Filter for winter (Jan-Mar) and summer (Jun-Sep) data in the training set.
# First make sure month column is integer.
if ("month" %in% names(training_df)) {
  if (is.list(training_df$month)) training_df$month <- unlist(training_df$month)
  training_df$month <- as.integer(training_df$month)
} else if ("date" %in% names(training_df)) {
  training_df$month <- lubridate::month(training_df$date)
}

# attach day-of-year for detrending computation
if (!"doy" %in% names(training_df) && "date" %in% names(training_df)) {
  training_df$doy <- lubridate::yday(training_df$date)
}

# compute seasonal polynomial baseline using summer months only.  Previous
# behaviour fitted a smooth spline to the full year; the current requirement is
# that trend metrics are derived solely from detrended summer values.
for (idx in intersect(INDICES_OF_INTEREST, names(training_df))) {
  trend_col <- paste0(idx, "_trend")
  training_df[[trend_col]] <- NA_real_

  # restrict fitting to summer months
  summer_fit <- training_df |> dplyr::filter(month %in% SUMMER_DETREND_MONTHS)
  finite_fit <- is.finite(summer_fit$doy) & is.finite(summer_fit[[idx]])
  if (sum(finite_fit) < MIN_SEASONAL_SAMPLES) {
    cat(sprintf("[TREND] Skipping polynomial for %s: insufficient summer training samples (%d)\n", idx, sum(finite_fit)))
    next
  }

  x_fit <- as.numeric(summer_fit$doy[finite_fit])
  y_fit <- as.numeric(summer_fit[[idx]][finite_fit])
  n_unique <- length(unique(x_fit))
  if (n_unique < 3) {
    cat(sprintf("[TREND] Skipping polynomial for %s: insufficient unique DOY values (%d)\n", idx, n_unique))
    next
  }

  # fit simple polynomial of specified degree rather than a spline
  f <- as.formula(paste(idx, "~ poly(doy,", DETREND_POLY_DEGREE, ")"))
  seasonal_model <- tryCatch(lm(f, data = summer_fit[finite_fit, ]), error = function(e) NULL)
  if (is.null(seasonal_model)) {
    cat(sprintf("[TREND] Skipping polynomial for %s: lm failed on summer values\n", idx))
    next
  }

  # predict only for summer rows; other months remain NA (unused later)
  pred_indices <- which(training_df$month %in% SUMMER_DETREND_MONTHS & is.finite(training_df$doy))
  if (length(pred_indices) > 0) {
    training_df[[trend_col]][pred_indices] <- predict(seasonal_model,
                                                      newdata = training_df[pred_indices, , drop = FALSE])
  }
}

# subsets
winter_data <- training_df |> dplyr::filter(month %in% c(1L,2L,3L))
# maintain backwards compatibility with earlier bias code
january_data <- winter_data

summer_data <- training_df |> dplyr::filter(month %in% c(6L,7L,8L,9L))
# apply detrending to summer subset using precomputed trend
for (idx in intersect(INDICES_OF_INTEREST, names(summer_data))) {
  trop <- paste0(idx, "_trend")
  if (trop %in% names(summer_data)) {
    global_mean <- mean(training_df[[trop]], na.rm = TRUE)
    summer_data[[paste0(idx, "_norm")]] <-
      summer_data[[idx]] - (summer_data[[trop]] - global_mean)
  }
}

cat("Number of winter (Jan-Mar) observations:", nrow(winter_data), "\n")
cat("Number of summer (Jun-Sep) observations:", nrow(summer_data), "\n")

# Compute dataset-level SNR for indices (per-index SNR across locations)
#
# This helper is mostly vestigial but implements the same variance-based approach
# used elsewhere: signal is the square of the seasonal amplitude (July minus
# January mean) and noise is the variance of the residuals around a fitted
# seasonal curve for each location-year.  The overall SNR is the median across
# location-years.
compute_global_index_snr <- function(df, indices, group_col = "location_id", eps = 1e-8) {
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
        if (length(vals_year) < 3) next
        # fit a simple seasonal spline for this year and compute residual variance
        doy <- if ("doy" %in% names(yr)) yr$doy else lubridate::yday(yr$date)
        fit <- tryCatch(stats::smooth.spline(doy, vals_year, df = min(5, length(unique(doy)) - 1)),
                        error = function(e) NULL)
        if (is.null(fit)) {
          noise_var <- stats::var(vals_year, na.rm = TRUE)
        } else {
          pred <- stats::predict(fit, doy)$y
          noise_var <- stats::var(vals_year - pred, na.rm = TRUE)
        }
        if (!is.finite(amp) || !is.finite(noise_var) || noise_var <= 0) next
        # convert amplitude to variance by squaring
        s_vals <- c(s_vals, (amp^2) / (noise_var + eps))
      }
    }
    if (length(s_vals) == 0) { out[idx] <- NA_real_; next }
    out[idx] <- median(s_vals, na.rm = TRUE)
  }
  out
}

# Calculate SNR for the indices of interest using training data
indices_to_snr <- INDICES_OF_INTEREST

## Compute SNR using estimated FVC fractions (variance-based noise)
## Signal variance is approximated by the square of the Dynamic Range
##   (Dynamic Range = 99th - 1st percentile of all data, soil->veg).
## Noise variance is estimated per-location by fitting a simple seasonal curve
## to `df_noise` and computing the variance of residuals; the median of the
## per-location variances is then used.  Final SNR = σ²_signal / σ²_noise.
compute_fvc_snr <- function(df_signal, df_noise, indices, group_col = "location_id", eps = 1e-8) {
  out <- numeric(length(indices)); names(out) <- indices
  for (idx in indices) {
    if (!idx %in% names(df_signal) || !idx %in% names(df_noise)) { out[idx] <- NA_real_; next }
    
    # 1. Estimate Dynamic Range (Signal) from full dataset (Soil to Veg)
    vals_sig <- df_signal[[idx]]
    vals_sig <- vals_sig[is.finite(vals_sig)]
    if (length(vals_sig) < 10) { out[idx] <- NA_real_; next }
    q <- quantile(vals_sig, probs = c(0.01, 0.99), na.rm = TRUE)
    dynamic_range <- as.numeric(diff(q))
    if (dynamic_range <= eps) dynamic_range <- eps
    signal_var <- dynamic_range^2
    
    # 2. Estimate Noise variance from location-wise residuals
    locs <- unique(na.omit(as.character(df_noise[[group_col]])))
    noise_vars <- numeric(0)
    
    for (loc in locs) {
      sub <- df_noise[df_noise[[group_col]] == loc, , drop = FALSE]
      vals <- sub[[idx]]
      vals <- vals[is.finite(vals)]
      if (length(vals) < 3) next
      
      # fit seasonal spline if date information is available
      if ("date" %in% names(sub)) {
        doy <- lubridate::yday(sub$date)
        fit <- tryCatch(stats::smooth.spline(doy, vals, df = min(5, length(unique(doy)) - 1)),
                        error = function(e) NULL)
        if (!is.null(fit)) {
          resid <- vals - stats::predict(fit, doy)$y
          var_loc <- stats::var(resid, na.rm = TRUE)
        } else {
          var_loc <- stats::var(vals, na.rm = TRUE)
        }
      } else {
        var_loc <- stats::var(vals, na.rm = TRUE)
      }
      if (is.finite(var_loc)) noise_vars <- c(noise_vars, var_loc)
    }
    
    if (length(noise_vars) == 0) { out[idx] <- NA_real_; next }
    median_noise_var <- median(noise_vars, na.rm = TRUE)
    
    # SNR = signal variance / noise variance
    out[idx] <- signal_var / (median_noise_var + eps)
  }
  out
}

training_nonbarren <- training_df
if ("Veg" %in% names(training_df)) training_nonbarren <- training_df[!(tolower(training_df$Veg) == "barren"), , drop = FALSE]

# Use training_df (all classes) for Signal (Dynamic Range)
# Use training_nonbarren for Noise
index_snr <- compute_fvc_snr(training_df, training_nonbarren, indices_to_snr, group_col = "location_id")

cat("Index SNR (Signal variance / noise variance, dynamic‑range squared over residual variance):\n")
print(index_snr)
cat("(above values were computed using estimated FVC fractions logic)
")

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
  
  # only use indices actually present in the data (avoids selection errors)
  metrics <- intersect(INDICES_OF_INTEREST, names(data))
  if (length(metrics) == 0) {
    warning(sprintf("No INDICES_OF_INTEREST columns available for %s, returning empty summary", month_name))
    return(data.frame())
  }
  missing <- setdiff(INDICES_OF_INTEREST, metrics)
  if (length(missing) > 0) {
    cat(sprintf("[NOTICE] Dropping missing indices for global %s averages: %s\n", month_name, paste(missing, collapse=", ")))
  }

  global_boot <- bootstrap_hierarchical_means(data, metrics = metrics, B = B)
  
  avg_stats <- data |> 
    dplyr::summarize(
      n_observations = dplyr::n(),
      n_locations = dplyr::n_distinct(location_id),
      across(all_of(metrics), ~ mean_of_means(.x, location_id), .names = "avg_{.col}")
    )
  
  # Bind CIs
  # global_boot is a named vector. Convert to 1-row DF.
  ci_df <- as.data.frame(t(global_boot))
  dplyr::bind_cols(avg_stats, ci_df)
}

# Calculate Global averages for winter (non-detrended)
winter_global_avg <- process_global_averages(winter_data, "Winter (Jan-Mar)")

# Calculate Global averages for summer (using detrended values)
# summer_data already contains '_norm' columns for each index
summer_global_avg <- process_global_averages(summer_data, "Summer (Jun-Sep detrended)")

# --- Averages by Vegetation Type (using mapping CSV at location-level) ---

# Build per-location summaries from observation-level data (do not join mapping to observations)
location_summary <- function(data, metrics = c("MSAVI", "NDVI", "PPI")) {
  data |> dplyr::group_by(location_id) |> dplyr::summarize(
    n_observations = dplyr::n(),
    across(any_of(metrics), ~ mean(.x, na.rm = TRUE), .names = "avg_{.col}"),
    .groups = "drop"
  )
}

loc_winter <- location_summary(winter_data, metrics = INDICES_OF_INTEREST)
loc_summer <- location_summary(summer_data, metrics = INDICES_OF_INTEREST)

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
  # Primary mapping: external map_loc. Fallback: derive Veg per location from the data itself.
  loc_joined <- NULL
  if (!is.null(map_loc)) {
    loc_joined <- loc_tbl |> dplyr::left_join(map_loc, by = "location_id")
  }

  if (is.null(loc_joined) || !"Veg" %in% names(loc_joined) || sum(!is.na(loc_joined$Veg)) == 0) {
    if ("Veg" %in% names(data_tbl)) {
      veg_from_data <- data_tbl |>
        dplyr::filter(!is.na(Veg) & trimws(as.character(Veg)) != "") |>
        dplyr::mutate(Veg = tolower(as.character(Veg))) |>
        dplyr::group_by(location_id, Veg) |>
        dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
        dplyr::group_by(location_id) |>
        dplyr::slice_max(order_by = n, n = 1, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::select(location_id, Veg)

      loc_joined <- loc_tbl |> dplyr::left_join(veg_from_data, by = "location_id")
      cat(sprintf("[VEG SUMMARY] %s: mapping join empty; using Veg labels from period data by location\n", period_name))
    } else {
      cat(sprintf("[VEG SUMMARY] %s: no Veg labels available (mapping and data both missing)\n", period_name))
      return(data.frame())
    }
  }

  loc_joined <- loc_joined |> dplyr::filter(!is.na(Veg))
  if (nrow(loc_joined) == 0) {
    cat(sprintf("[VEG SUMMARY] %s: no matched locations with Veg labels after join/fallback\n", period_name))
    return(data.frame())
  }

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
    if (nrow(sub) == 0) return(NULL)
    
    # Run bootstrap (returns named vector)
    res <- bootstrap_hierarchical_means(sub, metrics = INDICES_OF_INTEREST, B = B)
    
    # Convert res vector to 1-row DF
    res_df <- as.data.frame(t(res))
    res_df$Veg <- v
    res_df
  })
  
  veg_boot <- Filter(Negate(is.null), veg_boot)
  if (length(veg_boot) > 0) {
    veg_boot_df <- dplyr::bind_rows(veg_boot)
  } else {
    veg_boot_df <- data.frame(Veg = character(0))
  }

  dplyr::left_join(veg_loc_stats, veg_boot_df, by = "Veg")
}

cat("Calculating Vegetation Type winter averages using mapping CSV and location summaries...\n")
veg_type_averages <- veg_stats_from_locs(loc_winter, winter_data, "Winter (Jan-Mar)")

cat("Calculating Vegetation Type summer (detrended) averages using mapping CSV and location summaries...\n")
veg_type_averages_summer <- veg_stats_from_locs(loc_summer, summer_data, "Summer (Jun-Sep)")

# Print results
cat("\n=== GLOBAL WINTER AVERAGES (1985-2025) ===\n")
print(winter_global_avg)

cat("\n=== GLOBAL SUMMER AVERAGES (detrended) (1985-2025) ===\n")
print(summer_global_avg)

cat("\n=== WINTER AVERAGES BY VEGETATION TYPE (1985-2025) ===\n")
print(veg_type_averages)

cat("\n=== SUMMER AVERAGES BY VEGETATION TYPE (1985-2025) ===\n")
print(veg_type_averages_summer)

# --- Create Summary Table as requested ---
# Create global averages excluding 'barren' vegetation (for summary table)
exclude_barren_filter <- function(d) {
  if (is.null(d) || nrow(d) == 0) return(d)
  if (!"Veg" %in% names(d)) return(d)
  idx <- tolower(as.character(d$Veg)) == "barren"
  d[is.na(idx) | !idx, , drop = FALSE]
}

nonbarren_winter <- exclude_barren_filter(winter_data)
nonbarren_summer <- exclude_barren_filter(summer_data)

winter_global_avg_nonbarren <- if (nrow(nonbarren_winter) > 0) process_global_averages(nonbarren_winter, "Winter (Jan-Mar, Non-Barren)") else process_global_averages(nonbarren_winter[0,], "Winter (Empty)")
summer_global_avg_nonbarren <- if (nrow(nonbarren_summer) > 0) process_global_averages(nonbarren_summer, "Summer (Jun-Sep, Non-Barren)") else process_global_averages(nonbarren_summer[0,], "Summer (Empty)")
create_summary_entry <- function(idx_name, bias_info, winter_row, summer_row) {
  # Helper to format mean +/- margin
  fmt_val <- function(mean_val, lower, upper) {
    if (length(mean_val) == 0 || length(lower) == 0 || length(upper) == 0) return("NA")
    mean_val <- mean_val[[1]]
    lower <- lower[[1]]
    upper <- upper[[1]]
    if (!is.finite(mean_val) || !is.finite(lower) || !is.finite(upper)) return("NA")
    margin <- (upper - lower) / 2
    # Format with comma decimal
    m_str <- format(round(mean_val, 3), nsmall=3, decimal.mark=",")
    e_str <- format(round(margin, 3), nsmall=3, decimal.mark=",")
    sprintf("%s (+/-%s)*", m_str, e_str)
  }
  
  winter_str <- fmt_val(winter_row[[paste0("avg_", idx_name)]], winter_row[[paste0(idx_name, "_ci_lower")]], winter_row[[paste0(idx_name, "_ci_upper")]])
  summer_str <- fmt_val(summer_row[[paste0("avg_", idx_name)]], summer_row[[paste0(idx_name, "_ci_lower")]], summer_row[[paste0(idx_name, "_ci_upper")]])
  
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
  
  c(idx_name, bias_str, rmse_str, winter_str, summer_str)
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
  create_summary_entry(idx, bias_val, winter_global_avg_nonbarren, summer_global_avg_nonbarren)
})

summary_table <- do.call(rbind, summary_rows)
colnames(summary_table) <- c("INDEX", "Winter Bias", "Winter RMSE", "Winter Average (All areas, excl. barren)", "Summer Average (All areas, excl. barren)")
summary_table <- as.data.frame(summary_table)

cat("\n=== SUMMARY TABLE (excluding Veg=='barren') ===\n")
print(summary_table)

# Optional: Save results to Excel
output_dir <- "C:/MAP/january_averages_results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

results_list <- list(
  "Winter_Global_Avg" = winter_global_avg,
  "Summer_Global_Avg" = summer_global_avg,
  "Winter_VegType_Avg" = veg_type_averages,
  "Summer_VegType_Avg" = veg_type_averages_summer,
  "Summary_Table" = summary_table
)
# write each component of the results list to its own CSV file for easy inspection
for (nm in names(results_list)) {
  df_res <- results_list[[nm]]
  if (is.data.frame(df_res) && nrow(df_res) > 0) {
    try({
      csv_path <- file.path(output_dir, paste0(nm, ".csv"))
      write.csv(df_res, csv_path, row.names = FALSE)
      cat("Saved", nm, "to:", csv_path, "\n")
    }, silent = TRUE)
  }
}
# Add index Bias results to the output bundle
if (exists("index_bias_list") && length(index_bias_list) > 0) {
  bias_df <- do.call(rbind, lapply(names(index_bias_list), function(n) {
    b <- index_bias_list[[n]]
    data.frame(index = n, bias = b$bias, lower = b$lower, upper = b$upper, rmse = b$rmse, stringsAsFactors = FALSE)
  }))
  results_list$Index_Bias <- bias_df
}

openxlsx::write.xlsx(results_list, file = file.path(output_dir, "phenology_averages.xlsx"))
# also record the raw training data so that users have a copy of the "pre-inference" CSV
if (exists("training_df")) {
  try({
    train_csv <- file.path(output_dir, "pre_inference_training_data.csv")
    readr::write_csv(training_df, train_csv)
    cat("Saved training data (pre-inference) to:", train_csv, "\n")
  }, silent = TRUE)
}
if (exists("bias_df") && nrow(bias_df) > 0) {
  # Save the Bias CSV alongside the Excel file
  try({
    write.csv(bias_df, file = file.path(output_dir, "index_bias.csv"), row.names = FALSE)
    cat("Index Bias saved to:", file.path(output_dir, "index_bias.csv"), "\n")
  }, silent = TRUE)
}

cat("\nResults saved to:", file.path(output_dir, "phenology_averages.xlsx"), "\n")

# -----------------------------------------------------------------------------
# plotting and trend generation (skipped if running on training data only)
# -----------------------------------------------------------------------------
if (!SKIP_TRENDS) {

# Plot timeseries only for a selected subset of indices (user requested)
# INDICES_OF_INTEREST remains unchanged for other calculations.  PPI is
# treated identically to the other indices (seasonally normalized along with
# the rest).
indices_to_plot <- intersect(INDICES_OF_INTEREST,
                              c("NDVI", "MSAVI", "EVI", "PPI"))

plot_data_mean <- list()
plot_data_median <- list()

for (idx in indices_to_plot) {
  # skip if this index isn't present in the main dataframe (e.g. EVI might be missing)
  if (!(idx %in% names(df))) {
    cat(sprintf("Index '%s' not present in data - skipping trend computation\n", idx))
    next
  }

  cat(sprintf("Processing trend data for %s using summer-only polynomial detrending...\n", idx))
  
  # restrict to summer months before fitting
  summer_data <- df |> dplyr::filter(month %in% SUMMER_DETREND_MONTHS & is.finite(.data[[idx]]))
  if (nrow(summer_data) == 0) {
    cat(sprintf("Warning: No summer data for %s available. Skipping.\n", idx))
    next
  }

  # Ensure DOY is present
  if (!"doy" %in% names(summer_data)) {
    summer_data$doy <- lubridate::yday(summer_data$date)
  }

  # --- Normalization Step ---
  # fit polynomial trend on summer observations only
  finite_fit <- is.finite(summer_data$doy) & is.finite(summer_data[[idx]])
  if (sum(finite_fit) < MIN_SEASONAL_SAMPLES) {
    cat(sprintf("Warning: not enough summer samples for %s (have %d) - skipping.\n", idx, sum(finite_fit)))
    next
  }
  n_unique <- length(unique(summer_data$doy[finite_fit]))
  if (n_unique < 3) {
    cat(sprintf("Warning: insufficient unique DOY in summer for %s (%d) - skipping.\n", idx, n_unique))
    next
  }

  f <- as.formula(paste(idx, "~ poly(doy,", DETREND_POLY_DEGREE, ")"))
  seasonal_model <- tryCatch(lm(f, data = summer_data[finite_fit, ]), error = function(e) NULL)
  if (is.null(seasonal_model)) {
    cat(sprintf("Warning: polynomial fit failed for %s - skipping.\n", idx))
    next
  }
  summer_data$seasonal_trend <- predict(seasonal_model, newdata = summer_data)
  global_seasonal_mean <- mean(summer_data$seasonal_trend, na.rm = TRUE)
  summer_data[[paste0(idx, "_norm")]] <- summer_data[[idx]] -
      (summer_data$seasonal_trend - global_seasonal_mean)
  norm_metric <- paste0(idx, "_norm")

  # --- Bootstrapping with normalized metric ---
  summer_yearly_boot <- summer_data |> 
    dplyr::group_by(pheno_year) |> 
    dplyr::do({
      tmp <- .
      if (".trend_location_id" %in% names(tmp)) tmp$location_id <- as.character(tmp$.trend_location_id)
      res <- bootstrap_hierarchical_means(tmp, metrics = c(norm_metric), B = B)
      
      # Compute point estimate as mean of file-level means using normalized metric
      loc_means <- sapply(unique(tmp$location_id), function(id) {
        sub <- tmp[tmp$location_id == id, ]
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
      tmp <- .
      if (".trend_location_id" %in% names(tmp)) tmp$location_id <- as.character(tmp$.trend_location_id)
      res <- bootstrap_hierarchical_medians(tmp, metrics = c(norm_metric), B = B)
      
      # Compute point estimate as median of file-level medians using chosen metric
      loc_medians <- sapply(unique(tmp$location_id), function(id) {
        sub <- tmp[tmp$location_id == id, ]
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
       subtitle = "Seasonal normalization applied per index",
       x = "Year",
       y = "Mean Value") +
  theme_minimal()

print(p)

# Save the plot
ggsave(file.path(output_dir, "all_indices_summer_trend_mean.png"), plot = p, width = 16, height = 12)
cat("Saved mean plot to:", file.path(output_dir, "all_indices_summer_trend_mean.png"), "\n")
# also write the underlying data so users can inspect or reuse it later
try({
  readr::write_csv(plot_data_mean, file.path(output_dir, "all_indices_summer_trend_mean.csv"))
  cat("Saved mean trend data CSV to:", file.path(output_dir, "all_indices_summer_trend_mean.csv"), "\n")
}, silent = TRUE)

# Plot median
p_med <- ggplot(plot_data_median, aes(x = pheno_year, y = median_val)) +
  add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
  geom_line(color = "darkgreen") +
  geom_point(color = "orange") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "darkgreen") +
  facet_wrap(~ index, scales = "free_y") +
  labs(title = "Median Seasonally-Normalized June-Sept Vegetation Indices (1985-2025)",
       subtitle = "Seasonal normalization applied per index",
       x = "Year",
       y = "Median Value") +
  theme_minimal()

print(p_med)

# Save the plot
ggsave(file.path(output_dir, "all_indices_summer_trend_median.png"), plot = p_med, width = 16, height = 12)
cat("Saved median plot to:", file.path(output_dir, "all_indices_summer_trend_median.png"), "\n")
# also export the data used for the median graph
try({
  readr::write_csv(plot_data_median, file.path(output_dir, "all_indices_summer_trend_median.csv"))
  cat("Saved median trend data CSV to:", file.path(output_dir, "all_indices_summer_trend_median.csv"), "\n")
}, silent = TRUE)

# -----------------------------------------------------------------------------

# --- Combined trends for all selected indices on a single axes ---
# Both mean and median versions provide an overview of how the four indices
# evolve together without faceting.
combined_mean <- plot_data_mean |> dplyr::filter(index %in% indices_to_plot)
if (nrow(combined_mean) > 0) {
  p_comb_mean <- ggplot(combined_mean, aes(x = pheno_year, y = mean_val, color = index)) +
    add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    labs(title = "All Selected Indices: Mean Seasonally-Normalized Trends",
         subtitle = "Seasonal normalization applied to all indices",
         x = "Year", y = "Value") +
    theme_minimal()
  ggsave(file.path(output_dir, "all_indices_combined_mean.png"), plot = p_comb_mean, width = 10, height = 6)
  cat("Saved combined indices mean plot to:", file.path(output_dir, "all_indices_combined_mean.png"), "\n")
  # CSV export of combined mean
  try({
    readr::write_csv(combined_mean, file.path(output_dir, "all_indices_combined_mean.csv"))
    cat("Saved combined mean trend data CSV to:", file.path(output_dir, "all_indices_combined_mean.csv"), "\n")
  }, silent = TRUE)
  try({
    readr::write_csv(combined_med, file.path(output_dir, "all_indices_combined_median.csv"))
    cat("Saved combined median trend data CSV to:", file.path(output_dir, "all_indices_combined_median.csv"), "\n")
  }, silent = TRUE)
# inference file (one trend location per source CSV).
loc_index_series <- list()
for (idx in indices_to_plot) {
  if (!(idx %in% names(df))) next
    # treat all indices the same when selecting finite observations
    full_loc <- df |> dplyr::filter(is.finite(.data[[idx]]))
  if (!"doy" %in% names(full_loc)) {
    full_loc$doy <- lubridate::yday(full_loc$date)
  }

    # smooth spline on full location data
    n_unique <- length(unique(full_loc$doy))
    spline_df <- min(OUTLIER_SPLINE_MAX_DF, nrow(full_loc)/2, n_unique - 1)
    seasonal_spline_loc <- tryCatch(
      stats::smooth.spline(full_loc$doy, full_loc[[idx]], df = spline_df),
      error = function(e) NULL
    )
    if (is.null(seasonal_spline_loc)) next

    full_loc$seasonal_trend <- predict(seasonal_spline_loc, full_loc$doy)$y
    global_seasonal_mean_loc <- mean(full_loc$seasonal_trend, na.rm = TRUE)
    norm_metric_loc <- paste0(idx, "_norm")
    full_loc[[norm_metric_loc]] <- full_loc[[idx]] - (full_loc$seasonal_trend - global_seasonal_mean_loc)

    summer_data_loc <- full_loc
  loc_fallback <- summer_data_loc |>
    dplyr::group_by(.trend_location_id) |>
    dplyr::summarise(
      n_loc = sum(is.finite(.data[[norm_metric_loc]])),
      sd_loc = stats::sd(.data[[norm_metric_loc]], na.rm = TRUE),
      .groups = "drop"
    )
  loc_fallback$.trend_location_id <- as.character(loc_fallback$.trend_location_id)
  global_sd <- stats::sd(summer_data_loc[[norm_metric_loc]], na.rm = TRUE)
  loc_fallback$se_fallback <- loc_fallback$sd_loc / sqrt(pmax(loc_fallback$n_loc, 1))
  loc_fallback$se_fallback[!is.finite(loc_fallback$se_fallback) | loc_fallback$se_fallback <= 0] <-
    if (is.finite(global_sd) && global_sd > 0) global_sd / sqrt(4) else NA_real_

  # calculate per-location per-year mean and 95% CI
  loc_year <- summer_data_loc |>
    dplyr::group_by(.trend_location_id, pheno_year) |>
    dplyr::do({
      tmp <- .
      vals <- tmp[[norm_metric_loc]]
      mean_val <- mean(vals, na.rm = TRUE)
      if (sum(is.finite(vals)) < 2) {
        loc_id <- as.character(tmp$.trend_location_id[1])
        se_loc <- loc_fallback$se_fallback[loc_fallback$.trend_location_id == loc_id][1]
        if (!is.finite(se_loc) || se_loc <= 0) {
          ci <- c(mean_val, mean_val)
        } else {
          ci <- c(mean_val - 1.96 * se_loc, mean_val + 1.96 * se_loc)
        }
      } else {
        bs <- replicate(B, mean(sample(vals, replace = TRUE), na.rm = TRUE))
        ci <- quantile(bs, probs = c(0.025, 0.975), na.rm = TRUE)
      }
      data.frame(mean_val = mean_val,
                 ci_lower = ci[1],
                 ci_upper = ci[2])
    }) |>
    dplyr::ungroup()
  loc_year$index <- idx
  loc_index_series[[idx]] <- loc_year
}

if (length(loc_index_series) > 0) {
  loc_index_df <- dplyr::bind_rows(loc_index_series)
  if (nrow(loc_index_df) > 0) {
    loc_levels <- sort(unique(as.character(loc_index_df$.trend_location_id)))
    # Default shapes per location; explicitly keep "low" as circle and make
    # "mid" clearly different for readability.
    loc_shapes <- stats::setNames(rep(c(16, 17, 15, 18, 3, 4, 8), length.out = length(loc_levels)), loc_levels)
    loc_shapes[grepl("low", names(loc_shapes), ignore.case = TRUE)] <- 16
    loc_shapes[grepl("mid", names(loc_shapes), ignore.case = TRUE)] <- 17

    p_loc_combined <- ggplot(
      loc_index_df,
      aes(x = pheno_year, y = mean_val, color = index,
          linetype = as.factor(.trend_location_id),
          shape = as.factor(.trend_location_id),
          group = interaction(index, .trend_location_id))
    ) +
      add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
      geom_line(alpha = 0.9, linewidth = 0.8) +
      geom_point(size = 1.8, alpha = 0.9) +
      scale_shape_manual(values = loc_shapes) +
      labs(title = "All Selected Indices (All Locations Combined)",
           subtitle = "Linetype/point shape = location_id",
           x = "Year", y = "Normalized Value",
           color = "Index", linetype = "Location", shape = "Location") +
      theme_minimal()

    ggsave(file.path(output_dir, "all_indices_combined_all_locations_linetype.png"),
           plot = p_loc_combined, width = 13, height = 7)
    cat("Saved all-locations combined indices plot to:",
        file.path(output_dir, "all_indices_combined_all_locations_linetype.png"), "\n")

    # also generate individual-location plots with CI ribbons
    for (loc in unique(loc_index_df$.trend_location_id)) {
      df_loc <- loc_index_df[loc_index_df$.trend_location_id == loc, ]
      p_loc <- ggplot(df_loc, aes(x = pheno_year, y = mean_val, color = index, group = index)) +
        add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
        geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = index),
                    alpha = 0.2, colour = NA) +
        geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, alpha = 0.8) +
        geom_line(linewidth = 1) +
        geom_point(size = 1.5) +
        facet_wrap(~ index, scales = "free_y") +
        labs(title = sprintf("Location %s: Indices with 95%% CI", loc),
             x = "Year", y = "Value") +
        theme_minimal()
      fn <- file.path(output_dir, sprintf("location_%s_trends_ci.png", loc))
      ggsave(fn, plot = p_loc, width = 10, height = 6)
      cat(sprintf("Saved per-location plot for %s to: %s\n", loc, fn))
    }
  } else {
    cat("No valid location-level series to create all-locations combined linetype plot\n")
  }
} else {
  cat("No location-level index data available for all-locations combined linetype plot\n")
}

# -----------------------------------------------------------------------------
# Generate the same trend plots per inference scope (low/mid/kon), each in its
# own child output folder, mirroring MESMA's per-source folder structure.
# -----------------------------------------------------------------------------
run_scope_trend_plots <- function(df_scope, scope_name, scope_output_dir,
                                  indices_to_plot, B,
                                  summer_months = SUMMER_DETREND_MONTHS) {
  if (is.null(df_scope) || nrow(df_scope) == 0) {
    cat(sprintf("[SCOPE %s] No rows; skipping scope trend plots\n", scope_name))
    return(invisible(NULL))
  }
  if (!"month" %in% names(df_scope) && "date" %in% names(df_scope)) {
    df_scope$month <- lubridate::month(as.Date(df_scope$date))
  }
  if (!"pheno_year" %in% names(df_scope) && "date" %in% names(df_scope)) {
    df_scope$pheno_year <- assign_pheno_year(as.Date(df_scope$date))
  }
  if (!".trend_location_id" %in% names(df_scope)) {
    cat(sprintf("[SCOPE %s] Missing .trend_location_id; skipping scope trend plots\n", scope_name))
    return(invisible(NULL))
  }
  if (!dir.exists(scope_output_dir)) dir.create(scope_output_dir, recursive = TRUE)

  plot_data_mean_scope <- list()
  plot_data_median_scope <- list()

  for (idx in indices_to_plot) {
    if (!(idx %in% names(df_scope))) next
    # restrict to summer months for seasonal baseline
    full_data <- df_scope |> dplyr::filter(month %in% SUMMER_DETREND_MONTHS & is.finite(.data[[idx]]))
    if (nrow(full_data) == 0) next
    if (!"doy" %in% names(full_data)) full_data$doy <- lubridate::yday(full_data$date)

    # polynomial seasonal baseline on summer records
    finite_fit <- is.finite(full_data$doy) & is.finite(full_data[[idx]])
    if (sum(finite_fit) < MIN_SEASONAL_SAMPLES) next
    n_unique <- length(unique(full_data$doy[finite_fit]))
    if (n_unique < 3) next
    f <- as.formula(paste(idx, "~ poly(doy,", DETREND_POLY_DEGREE, ")"))
    seasonal_model <- tryCatch(lm(f, data = full_data[finite_fit, ]), error = function(e) NULL)
    if (is.null(seasonal_model)) next

    full_data$seasonal_trend <- predict(seasonal_model, newdata = full_data)
    global_seasonal_mean <- mean(full_data$seasonal_trend, na.rm = TRUE)
    norm_metric <- paste0(idx, "_norm")
    full_data[[norm_metric]] <- full_data[[idx]] - (full_data$seasonal_trend - global_seasonal_mean)
    
    summer_data <- full_data

    summer_yearly_boot <- summer_data |>
      dplyr::group_by(pheno_year) |>
      dplyr::do({
        tmp <- .
        tmp$location_id <- as.character(tmp$.trend_location_id)
        res <- bootstrap_hierarchical_means(tmp, metrics = c(norm_metric), B = B)
        loc_means <- sapply(unique(tmp$location_id), function(id) {
          sub <- tmp[tmp$location_id == id, ]
          mean(sub[[norm_metric]], na.rm = TRUE)
        })
        data.frame(mean_val = mean(loc_means, na.rm = TRUE), ci_lower = res[1], ci_upper = res[2])
      }) |>
      dplyr::ungroup()
    summer_yearly_boot$index <- idx
    plot_data_mean_scope[[idx]] <- summer_yearly_boot

    summer_yearly_boot_median <- summer_data |>
      dplyr::group_by(pheno_year) |>
      dplyr::do({
        tmp <- .
        tmp$location_id <- as.character(tmp$.trend_location_id)
        res <- bootstrap_hierarchical_medians(tmp, metrics = c(norm_metric), B = B)
        loc_medians <- sapply(unique(tmp$location_id), function(id) {
          sub <- tmp[tmp$location_id == id, ]
          median(sub[[norm_metric]], na.rm = TRUE)
        })
        data.frame(median_val = median(loc_medians, na.rm = TRUE), ci_lower = res[1], ci_upper = res[2])
      }) |>
      dplyr::ungroup()
    summer_yearly_boot_median$index <- idx
    plot_data_median_scope[[idx]] <- summer_yearly_boot_median
  }

  if (length(plot_data_mean_scope) == 0 || length(plot_data_median_scope) == 0) {
    cat(sprintf("[SCOPE %s] No valid trend series after filtering; skipping plot export\n", scope_name))
    return(invisible(NULL))
  }
  # write scope-specific CSVs as well
  try({
    readr::write_csv(plot_data_mean_scope, file.path(scope_output_dir, "trend_mean.csv"))
    cat(sprintf("[SCOPE %s] Saved mean trend data CSV to: %s\n", scope_name, file.path(scope_output_dir, "trend_mean.csv")))
  }, silent = TRUE)
  try({
    readr::write_csv(plot_data_median_scope, file.path(scope_output_dir, "trend_median.csv"))
    cat(sprintf("[SCOPE %s] Saved median trend data CSV to: %s\n", scope_name, file.path(scope_output_dir, "trend_median.csv")))
  }, silent = TRUE)

  plot_data_mean_scope <- do.call(rbind, plot_data_mean_scope)
  plot_data_median_scope <- do.call(rbind, plot_data_median_scope)

  p_mean_scope <- ggplot(plot_data_mean_scope, aes(x = pheno_year, y = mean_val)) +
    add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
    geom_line(color = "blue") +
    geom_point(color = "red") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "blue") +
    facet_wrap(~ index, scales = "free_y") +
    labs(title = sprintf("%s: Mean June-Sept Vegetation Indices", toupper(scope_name)),
         subtitle = "Seasonal normalization applied per index",
         x = "Year", y = "Mean Value") +
    theme_minimal()
  ggsave(file.path(scope_output_dir, "all_indices_summer_trend_mean.png"), plot = p_mean_scope, width = 16, height = 12)

  p_median_scope <- ggplot(plot_data_median_scope, aes(x = pheno_year, y = median_val)) +
    add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
    geom_line(color = "darkgreen") +
    geom_point(color = "orange") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "darkgreen") +
    facet_wrap(~ index, scales = "free_y") +
    labs(title = sprintf("%s: Median June-Sept Vegetation Indices", toupper(scope_name)),
         subtitle = "Seasonal normalization applied per index",
         x = "Year", y = "Median Value") +
    theme_minimal()
  ggsave(file.path(scope_output_dir, "all_indices_summer_trend_median.png"), plot = p_median_scope, width = 16, height = 12)

  comb_mean_scope <- plot_data_mean_scope |> dplyr::filter(index %in% indices_to_plot)
  p_comb_mean_scope <- ggplot(comb_mean_scope, aes(x = pheno_year, y = mean_val, color = index)) +
    add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
    geom_line(size = 1) + geom_point(size = 1.5) +
    labs(title = sprintf("%s: All Selected Indices (Mean)", toupper(scope_name)),
         x = "Year", y = "Normalized Value") +
    theme_minimal()
  ggsave(file.path(scope_output_dir, "all_indices_combined_mean.png"), plot = p_comb_mean_scope, width = 10, height = 6)

  comb_med_scope <- plot_data_median_scope |> dplyr::filter(index %in% indices_to_plot)
  p_comb_med_scope <- ggplot(comb_med_scope, aes(x = pheno_year, y = median_val, color = index)) +
    add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
    geom_line(size = 1) + geom_point(size = 1.5) +
    labs(title = sprintf("%s: All Selected Indices (Median)", toupper(scope_name)),
         x = "Year", y = "Normalized Value") +
    theme_minimal()
  ggsave(file.path(scope_output_dir, "all_indices_combined_median.png"), plot = p_comb_med_scope, width = 10, height = 6)

  cat(sprintf("[SCOPE %s] Saved trend plots to: %s\n", scope_name, scope_output_dir))
  invisible(list(mean = plot_data_mean_scope, median = plot_data_median_scope))
}

scope_levels <- intersect(c("low", "mid", "kon"), unique(as.character(df$.inference_scope)))
if (length(scope_levels) > 0) {
  cat(sprintf("[SCOPE] Generating per-scope trend plots for: %s\n", paste(scope_levels, collapse = ", ")))
  for (scope_name in scope_levels) {
    scope_df <- df |> dplyr::filter(.inference_scope == scope_name)
    scope_outdir <- file.path(output_dir, scope_name)
    run_scope_trend_plots(
      df_scope = scope_df,
      scope_name = scope_name,
      scope_output_dir = scope_outdir,
      indices_to_plot = indices_to_plot,
      B = B,
      summer_months = 6:9
    )
    # also compute and save simple pre-inference summary tables for this scope
    try({
      wavg <- process_global_averages(scope_df |> dplyr::filter(month %in% c(1L,2L,3L)),
                                      paste(scope_name, "Winter"))
      savg <- process_global_averages(scope_df |> dplyr::filter(month %in% SUMMER_DETREND_MONTHS),
                                      paste(scope_name, "Summer"))
      veg_w <- veg_stats_from_locs(
                location_summary(scope_df |> dplyr::filter(month %in% c(1L,2L,3L))),
                scope_df, paste(scope_name, "Winter"))
      veg_s <- veg_stats_from_locs(
                location_summary(scope_df |> dplyr::filter(month %in% SUMMER_DETREND_MONTHS)),
                scope_df, paste(scope_name, "Summer"))
      scope_results <- list(
        Winter_Global = wavg,
        Summer_Global = savg,
        Winter_Veg = veg_w,
        Summer_Veg = veg_s
      )
      openxlsx::write.xlsx(scope_results,
                           file = file.path(scope_outdir, "pre_inference_summary.xlsx"))
      cat(sprintf("[SCOPE %s] Saved scope summary workbook to: %s\n", scope_name,
                  file.path(scope_outdir, "pre_inference_summary.xlsx")))
      # also CSV versions
      for (nm in names(scope_results)) {
        dfres <- scope_results[[nm]]
        if (is.data.frame(dfres) && nrow(dfres)>0) {
          csvfile <- file.path(scope_outdir, paste0(nm, ".csv"))
          write.csv(dfres, csvfile, row.names = FALSE)
          cat(sprintf("[SCOPE %s] Saved %s CSV to: %s\n", scope_name, nm, csvfile))
        }
      }
    }, silent = TRUE)
  }
} else {
  cat("[SCOPE] No low/mid/kon scope tags detected in input; skipping per-scope folder outputs\n")
}

# end of SKIP_TRENDS block
}

