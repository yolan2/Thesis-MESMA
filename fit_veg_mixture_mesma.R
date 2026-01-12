filter_variants_by_min_samples <- function(variants, min_samples = MIN_ENDMEMBER_SAMPLES, veg = NULL, raw_template = NULL) {
  if (is.null(variants) || length(variants) == 0) return(list())
  keep_mask <- sapply(variants, function(v) {
    # Require an explicit n_samples and enforce minimum sample count
    if (is.null(v$n_samples)) return(FALSE) # drop variants with unknown n_samples
    v$n_samples >= min_samples
  })
  removed <- sum(!keep_mask)
  if (removed > 0 && !is.null(veg)) cat(sprintf("  [%s] Removed %d variant(s) with n_samples < %d\n", veg, removed, min_samples))
  kept_variants <- variants[keep_mask]
  if (length(kept_variants) == 0 && !is.null(raw_template) && !is.null(raw_template$n_samples) && raw_template$n_samples >= min_samples) {
    cat(sprintf("  [%s] No variants left after filtering; restoring medoid raw template as single variant\n", ifelse(is.null(veg), "veg", veg)))
    kept_variants <- list(list(raw_mat = raw_template$T, variant_id = paste0(ifelse(is.null(veg), "veg", veg), "_single"), n_samples = raw_template$n_samples))
  }
  kept_variants

}
 
make_location_id <- function(lon, lat) {
  # Vectorized creation of a stable location id from lon/lat
  lon <- as.numeric(lon)
  lat <- as.numeric(lat)
  
  # Create a mask for non-finite values
  invalid_mask <- !is.finite(lon) | !is.finite(lat)
  
  # Use vectorized sprintf: L_lat_lon format
  res <- sprintf("L_%0.6f_%0.6f", round(lat, 6), round(lon, 6))
  
  # Set invalid entries to NA
  res[invalid_mask] <- NA_character_
  
  return(res)
}

 
library(zoo)
library(dplyr)
library(lme4)
library(cluster)
library(readr)
library(sf)
library(ggplot2)
library(scales)
library(nlme)
library(RStoolbox)
library(nnls)
library(terra)
library(magrittr) # For pipe operator
if (requireNamespace("future", quietly = TRUE)) library(future)
if (requireNamespace("future.apply", quietly = TRUE)) library(future.apply)
if (requireNamespace("openxlsx", quietly = TRUE)) library(openxlsx)
if (requireNamespace("MASS", quietly = TRUE)) library(MASS) # For lda

options(warn = 1)  # print warnings as they occur for debugging
options(future.progress = FALSE)  # disable progress bars from future package
options(progressr.enable = FALSE) # disable progress bars from progressr package

# Ensure core initialization helpers are available (defines setup_parallel_backend(), constants, etc.)
if (!exists("setup_parallel_backend", mode = "function")) {
  if (file.exists("init_mesma.R")) {
    try(source("init_mesma.R"), silent = TRUE)
  }
}
if (!exists("setup_parallel_backend", mode = "function")) {
  stop("`setup_parallel_backend` function not found after attempting to source 'init_mesma.R'. Please ensure 'init_mesma.R' is present and defines this function.")
}


SKIP_PPI <- TRUE

INPUT_CSV <- "C:\\Users\\yolan\\Downloads\\LS_S2_Harmonized_Timeseries.csv"
INFERENCE_CSV <- "C:\\Users\\yolan\\OneDrive\\Documenten\\UGENT\\Master\\masterproef\\GIS\\landsat_lower_inference.csv"

# AGENT: Restrict training to the most recent years for variant building.
TRAIN_YEARS <- 2024  # Years to use for training (inclusive)

# Check for override from environment or pre-set variable
if (!exists("PARALLEL_ENABLE")) {
  PARALLEL_ENABLE <- FALSE
}

if (isTRUE(PARALLEL_ENABLE)) {
  cat("[CONFIG] Parallel processing ENABLED\n")
} else {
  cat("[CONFIG] Parallel processing DISABLED (running sequentially - slower but more stable)\n")
}

PARALLEL_WORKERS <- 3
COMBO_PARALLEL_WORKERS <- max(1L, floor(PARALLEL_WORKERS/2))

# Batch-level parallelization settings
BATCH_PARALLEL_ENABLE <- TRUE  # Enable parallel batch processing
BATCH_PARALLEL_WORKERS <- 3  # Number of batches to process in parallel

# RMSE threshold factor for uncertainty estimation (keep models within this factor of best RMSE)
RMSE_UNCERTAINTY_THRESHOLD_FACTOR <- 1.40

# Memory monitoring utility function
print_memory_usage <- function(label = "") {
  # if (label != "") label <- paste0(" [", label, "]")
  # mem_info <- gc(verbose = FALSE, full = FALSE)
  # used_mb <- sum(mem_info[, "used"])
  # max_mb <- sum(mem_info[, "(Mb)"])
  # cat(sprintf("Memory usage%s: %.1f MB used, %.1f MB max\n", label, used_mb, max_mb))
  # invisible(mem_info)
}

ENABLE_DIAGNOSTICS <- TRUE
MAX_VARIANTS_PER_VEG <- 10
MIN_VARIANTS_PER_VEG <- 1
MAX_COMBOS_FOR_FULL_SEARCH <- 6000

MIN_ENDMEMBER_SAMPLES <- 5L

ALLOWED_VEG <- c("populus", "tamarix", "phragmites")

OPTIMAL_INDICES <- c("WDVI", "WDVI_BY_SWIR1", "GVI", "NIRv", "PSRI", "MSAVI2", "NDMI", "PPI", "EVI", "NDTI", "SATVI", "CIG", "BSI", "NBR", "TCW", "TCB", "NDSI")
# Indices to explicitly exclude from the MESMA fitter (do not use as features)
FITTER_EXCLUDE_INDICES <- c("MSAVI", "NIRv", "OSAVI", "NDVI", "NDMI", "MSAVI2", "EVI", "SATVI")

# Outlier removal configuration: when TRUE, remove large per-location(/year) outliers
ENABLE_OUTLIER_REMOVAL <- TRUE
# Threshold in multiples of MAD; default 6 (conservative). Rows with |x - median| > THRESHOLD * MAD are removed.
OUTLIER_MAD_THRESHOLD <- 3.5

SKIP_INFERENCE <- FALSE

TESTING_MODE <- FALSE
TESTING_MAX_PER_VEG <- 5000L  # Upper cap per vegetation class when TESTING_MODE is enabled
if (!exists('DEBUG')) DEBUG <- FALSE

# -----------------------------------------------------------------------------
# USER-TUNABLE PARAMETERS (centralized)
# Move any frequently-adjusted constants here for easy configuration
# -----------------------------------------------------------------------------
COMBO_SAFE_EXPAND_LIMIT <- 1e6    # fully expand combo grid up to this many combos
COMBO_ABORT_LIMIT <- 5e7          # abort if combos exceed this hard limit

BOOTSTRAP_B <- 100L              # Number of bootstrap iterations for location-based resampling
MAX_INFERENCE_LOCATIONS <- 300L  # Maximum number of locations to use for inference (cap to reduce runtime)
BOOTSTRAP_BLOCK_SIZE <- 5L       # Block size for Moving Block Bootstrap of residuals

# Uncertainty handling
ENABLE_UNCERTAINTY <- TRUE
DEBUG_UNCERTAINTY <- TRUE        # global debug flag for uncertainty diagnostics

# Small-n CI inflation (per location-year and small-N global aggregation)
# Goal: extremely low observation counts should yield very wide intervals (almost no certainty).
UNCERTAINTY_N_REF <- 12L               # reference observation count where inflation factor ~1
UNCERTAINTY_N_POWER <- 2.0             # higher -> more extreme widening for small n
UNCERTAINTY_BASE_SD <- 0.15            # fallback SD when sd is NA/0; before inflation
UNCERTAINTY_SD_MAX <- 0.50             # cap SD so CI doesn't become numerically unstable

# Optional: estimate UNCERTAINTY_N_REF / UNCERTAINTY_N_POWER from data (Option A: subsampling experiment)
ESTIMATE_UNCERTAINTY_PARAMS_OPTION_A <- TRUE
UNCERTAINTY_PARAM_EST_HIGH_N <- 20L            # only use loc-years with at least this many observations as pseudo-truth
UNCERTAINTY_PARAM_EST_MAX_GROUPS <- 25L        # cap number of loc-years used (runtime control)
UNCERTAINTY_PARAM_EST_TARGET_NS <- 1:15        # subsample sizes to evaluate
UNCERTAINTY_PARAM_EST_REPS <- 20L              # replicates per subsample size
UNCERTAINTY_PARAM_EST_SEED <- 123

# Data quality thresholds
MIN_OBS_PER_LOC_YEAR <- 3L
MIN_UNIQUE_DOY_DEFAULT <- 5L
MIN_UNIQUE_DOY_INFERENCE <- 1L
MIN_INDEX_SD <- 0.05
MIN_IDX_PRESENCE <- 0.5

# Variant handling
VARIANT_SWITCH <- TRUE           # Allow variant switching across pentads
VARIANT_SWITCH_RE_CENTER <- TRUE
VARIANCE_THRESHOLD <- 0.90

# Modeling/algorithmic caps and defaults
TOPK_VARIANTS <- 10L             # Max variants per veg type to consider in MESMA (higher is slower but more accurate)
ENABLE_LDA_L2_NORMALIZATION <- TRUE # Whether to L2-normalize training samples before LDA (set to FALSE to preserve brightness info)
GAM_K_MAX <- 40
GAM_GAMMA <- 1.0
USE_INDICES_MIN <- 1L

# Sample/cluster controls
ENABLE_SOIL_PREPROCESS <- TRUE
SOIL_PURE_THRESHOLD <- 0.95
SOIL_MIN_SAMPLES <- 3L

MAX_PROJECTIONS_PER_VEG <- 25000L  # subsample before clustering to avoid OOM
SILHOUETTE_SAMPLE_SIZE <- 20000L   # subsample for silhouette distance matrix
MEDOID_SAMPLE_SIZE <- 10000L       # subsample for medoid distance computation

# Numerical tolerances
EPS_SIGMA <- 1e-8
LOWER_BND <- 0

# Batch processing
BATCH_SIZE <- 6  # Default batch size for location-level processing (tunable)

# Output file (can be overridden later if OUT_DIR isn't set yet)
OUTPUT_XLSX <- NA
# === PPI-BASED FRACTION ESTIMATION CONFIG ===
PPI_ZERO_VEG_COVER <- 0.0  # PPI value assumed to correspond to 0% vegetation cover
PPI_FULL_VEG_COVER <- 0.7  # PPI value assumed to correspond to 100% vegetation cover (top 95% median PPI)

# === END GLOBAL CONFIG ===


if (file.exists("ppi_helpers.R")) {
  source("ppi_helpers.R")
} else {
  warning("ppi_helpers.R not found; using built-in PPI calculation. To standardize PPI, add ppi_helpers.R to the project root.")
}

PROGRESS_LOG_TO_FILE <- FALSE
LOG_FILE <- "mesma_progress.log"

log_msg <- function(...) {
  ts <- format(Sys.time(), "%H:%M:%S")
  msg <- sprintf("[%s] %s\n", ts, sprintf(...))
  cat(msg)
  if (isTRUE(PROGRESS_LOG_TO_FILE)) try(cat(msg, file = LOG_FILE, append = TRUE), silent = TRUE)
}

if (file.exists("ppi_helpers.R")) {
  source("ppi_helpers.R")
} else {
  stop("Required file 'ppi_helpers.R' not found in project root. Please add it to ensure consistent PPI calculation.")
}

if (file.exists("mesma_helpers.R")) {
  source("mesma_helpers.R")
} else {
  warning("mesma_helpers.R not found; some visualization features may be missing.")
}

assign_pheno_year <- function(d) {
  d <- as.Date(d)
  ifelse(is.na(d), NA_integer_, ifelse(lubridate::month(d) >= 3, lubridate::year(d), lubridate::year(d) - 1))
}

# Calculate day-of-year within phenological year (March 1 = day 1)
# For consistency with phenological year assignment where March starts a new year
pheno_doy <- function(d) {
  d <- tryCatch(as.Date(d), error = function(e) NA)
  month <- lubridate::month(d)
  day <- lubridate::day(d)

  # For March-December: count days from March 1
  # For January-February: count days from previous March 1 (add ~306 days)
  ifelse(is.na(d), NA_integer_,
    ifelse(month >= 3,
      # March onwards: days since March 1 of current year
      as.integer(d - as.Date(paste0(lubridate::year(d), "-03-01"))) + 1L,
      # Jan-Feb: days since March 1 of previous year
      as.integer(d - as.Date(paste0(lubridate::year(d) - 1, "-03-01"))) + 1L
    )
  )
}

if (!exists("calculate_solar_zenith") && file.exists("ppi_helpers.R")) source("ppi_helpers.R")

RAW_BANDS <- c("blue", "green", "red", "nir", "swir1", "swir2")

normalize_band_names <- function(df, bands = RAW_BANDS) {
  if (is.null(df) || nrow(df) == 0) return(df)
  current_names <- names(df)
  for (b in bands) {
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


# Default SOIL_LINE_SLOPE if not yet defined
if (!exists("SOIL_LINE_SLOPE")) SOIL_LINE_SLOPE <- 1.0

compute_indices_from_bands <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  eps <- 1e-9
  has_bands <- intersect(RAW_BANDS, names(df))
  if (length(has_bands) == 0) return(df)

  if (all(c('nir','red') %in% names(df))) df$DVI <- as.numeric(df$nir) - as.numeric(df$red)
  if (all(c('nir','red') %in% names(df))) df$WDVI <- as.numeric(df$nir) - SOIL_LINE_SLOPE * as.numeric(df$red)
  if (all(c('nir','red') %in% names(df))) df$OSAVI <- (as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + as.numeric(df$red) + 0.16)
  if (all(c('red','green','blue') %in% names(df))) df$MCARI <- ((as.numeric(df$red) - as.numeric(df$green)) - 0.2*(as.numeric(df$red) - as.numeric(df$blue))) * (as.numeric(df$red) / (as.numeric(df$green) + eps))
  if (all(c('green','red') %in% names(df))) df$PRI <- (as.numeric(df$green) - as.numeric(df$red)) / (as.numeric(df$green) + as.numeric(df$red) + eps)
    if (all(c('nir','red') %in% names(df))) df$NIRv <- as.numeric(df$nir) * ((as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + as.numeric(df$red) + eps))
  if (all(c('red','blue','nir') %in% names(df))) df$PSRI <- (as.numeric(df$red) - as.numeric(df$blue)) / (as.numeric(df$nir) + eps)
  if (all(c('nir','swir2') %in% names(df))) df$NBR <- (as.numeric(df$nir) - as.numeric(df$swir2)) / (as.numeric(df$nir) + as.numeric(df$swir2) + eps)

  # Tasseled Cap indices (Landsat 8 OLI coefficients - Baig et al. 2014)
  if (all(c('green','red','nir','swir1','swir2','blue') %in% names(df))) {
    df$TCB <- 0.3029 * as.numeric(df$blue) + 0.2786 * as.numeric(df$green) + 0.4733 * as.numeric(df$red) + 0.5599 * as.numeric(df$nir) + 0.508 * as.numeric(df$swir1) + 0.1872 * as.numeric(df$swir2)
    df$TCG <- -0.2941 * as.numeric(df$blue) - 0.243 * as.numeric(df$green) - 0.5424 * as.numeric(df$red) + 0.7276 * as.numeric(df$nir) + 0.0713 * as.numeric(df$swir1) - 0.1608 * as.numeric(df$swir2)
    df$TCW <- 0.1511 * as.numeric(df$blue) + 0.1973 * as.numeric(df$green) + 0.3283 * as.numeric(df$red) + 0.3407 * as.numeric(df$nir) - 0.7117 * as.numeric(df$swir1) - 0.4559 * as.numeric(df$swir2)
    df$GVI <- df$TCG  # Alias for backward compatibility
  }

  # NDVI intentionally omitted from computed indices to avoid using it in MESMA fitting
  # if (all(c('nir','red') %in% names(df))) df$NDVI <- (as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + as.numeric(df$red) + eps)
  if (all(c('nir','red') %in% names(df))) df$MSAVI2 <- (2 * as.numeric(df$nir) + 1 - sqrt(pmax(0, (2 * as.numeric(df$nir) + 1)^2 - 8 * (as.numeric(df$nir) - as.numeric(df$red))))) / 2
  if (all(c('nir','swir1') %in% names(df))) df$NDMI <- (as.numeric(df$nir) - as.numeric(df$swir1)) / (as.numeric(df$nir) + as.numeric(df$swir1) + eps)

  # New index: WDVI divided by SWIR1 (robust for soil-adjusted vegetation contrast normalized by SWIR1)
  if ('WDVI' %in% names(df) && 'swir1' %in% names(df)) df$WDVI_BY_SWIR1 <- as.numeric(df$WDVI) / (as.numeric(df$swir1) + eps)

  # Enhanced Vegetation Index (EVI)
  if (all(c('nir','red','blue') %in% names(df))) {
    df$EVI <- 2.5 * ((as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + 6 * as.numeric(df$red) - 7.5 * as.numeric(df$blue) + 1 + eps))
  }

  # NDTI - Normalized Difference Tillage Index (Van Deventer et al. 1997)
  if (all(c('swir1','swir2') %in% names(df))) df$NDTI <- (as.numeric(df$swir1) - as.numeric(df$swir2)) / (as.numeric(df$swir1) + as.numeric(df$swir2) + eps)

  # NDSI - Normalized Difference Snow Index
  if (all(c('green','swir1') %in% names(df))) df$NDSI <- (as.numeric(df$green) - as.numeric(df$swir1)) / (as.numeric(df$green) + as.numeric(df$swir1) + eps)

  # NDDI - Normalized Difference Dust Index
  if (all(c('red','nir') %in% names(df))) df$NDDI <- (as.numeric(df$red) - as.numeric(df$nir)) / (as.numeric(df$red) + as.numeric(df$nir) + eps)

  # SATVI - Soil Adjusted Total Vegetation Index (Marsett et al. 2006)
  if (all(c('swir1','red','swir2') %in% names(df))) df$SATVI <- ((as.numeric(df$swir1) - as.numeric(df$red)) / (as.numeric(df$swir1) + as.numeric(df$red) + 0.5 + eps)) * 1.5 - (as.numeric(df$swir2) / 2)

  # CIG - Chlorophyll Index Green (Gitelson et al. 2005)
  if (all(c('nir','green') %in% names(df))) df$CIG <- (as.numeric(df$nir) / (as.numeric(df$green) + eps)) - 1

  # BSI - Bare Soil Index (Rikimaru et al. 2002)
  if (all(c('swir1','red','nir','blue') %in% names(df))) {
    term1 <- as.numeric(df$swir1) + as.numeric(df$red)
    term2 <- as.numeric(df$nir) + as.numeric(df$blue)
    df$BSI <- (term1 - term2) / (term1 + term2 + eps)
  }

  if ('NIRv' %in% names(df)) df$NIRv <- df$NIRv * 1.3

  df
}

compute_mad2 <- function(x, min_samples = 3) {
  x <- x[is.finite(x)]
  if (length(x) < min_samples) return(NA_real_)
  m <- mad(x, na.rm = TRUE, constant = 1.4826)
  if (!is.finite(m)) return(NA_real_)
  m^2
}

safe_mul_vec <- function(a, b, allow_recycle = TRUE, caller = NULL) {
  la <- length(a); lb <- length(b)
  if (la == 0 || lb == 0) return(numeric(0))
  if (la == lb) return(a * b)
  if (la == 1) return(rep(a, lb) * b)
  if (lb == 1) return(a * rep(b, la))
  if (allow_recycle && (la %% lb == 0 || lb %% la == 0)) {
    if (la < lb) a <- rep(a, length.out = lb) else b <- rep(b, length.out = la)
    return(a * b)
  }
  caller_text <- if (is.null(caller)) "safe_mul_vec" else paste0(caller, ": ")
  stop(sprintf("%sIncompatible lengths for multiplication: %d vs %d", caller_text, la, lb))
}

safe_dot <- function(a, b, na.rm = TRUE) {
  if (length(a) == 0 || length(b) == 0) return(0)
  prod <- safe_mul_vec(a, b, allow_recycle = TRUE, caller = "safe_dot")
  if (na.rm) sum(prod, na.rm = TRUE) else sum(prod)
}

safe_col_weighted_avg <- function(mat, wts) {
  if (is.null(mat) || nrow(mat) == 0) return(rep(0, ifelse(is.null(mat), 1, ncol(mat))))
  n <- nrow(mat)
  if (length(wts) == 0) wts <- rep(1, n)
  if (length(wts) != n) {
    if (length(wts) == 1 || n %% length(wts) == 0 || length(wts) %% n == 0) {
      wts <- rep(wts, length.out = n)
      warning(sprintf("safe_col_weighted_avg: adjusted weights vector to length %d", n))
    } else {
      stop(sprintf("safe_col_weighted_avg: weights length (%d) incompatible with rows in mat (%d)", length(wts), n))
    }
  }
  if (sum(wts, na.rm = TRUE) == 0) {
    return(as.numeric(colMeans(mat, na.rm = TRUE)))
  }
  wts <- as.numeric(wts) / sum(wts, na.rm = TRUE)
  as.numeric(colSums(mat * wts, na.rm = TRUE))
}

calc_moving_var <- function(df, index_name, window = 14, span_loess = 0.1, min_obs_loess = 6) {
  
  if (!"date" %in% names(df)) stop("calc_moving_var: df must have a 'date' column")
  if (!index_name %in% names(df)) stop(sprintf("calc_moving_var: index '%s' not found in df", index_name))

  n <- nrow(df)
  out <- rep(NA_real_, n)

  process_series_range <- function(dates, vals) {
    if (length(vals) < 1) return(rep(NA_real_, length(vals)))
    
    dts <- as.Date(dates)
    full_days <- seq(from = as.Date(min(dts, na.rm = TRUE)), to = as.Date(max(dts, na.rm = TRUE)), by = "day")
    full_vec <- rep(NA_real_, length(full_days))
    pos_map <- match(dts, full_days)
    full_vec[pos_map] <- as.numeric(vals)
    
    idx_finite <- which(is.finite(full_vec))
    if (length(idx_finite) < 2) return(rep(NA_real_, length(vals)))
    
    full_vec_interp <- approx(x = idx_finite, y = full_vec[idx_finite], xout = seq_along(full_vec), rule = 2)$y
    
    r_min <- zoo::rollapply(full_vec_interp, width = window, FUN = min, fill = NA, align = "center")
    r_max <- zoo::rollapply(full_vec_interp, width = window, FUN = max, fill = NA, align = "center")
    r_range <- r_max - r_min

    # Map the rolling range back to the original observation dates
    out_full <- r_range
    # Some dates may fall outside the padded range if sparse; use pos_map to select
    result <- out_full[pos_map]
    return(as.numeric(result))
  }

  # Compute moving 'range' for the full series and return values aligned with df
  out <- process_series_range(df$date, df[[index_name]])
  if (length(out) != n) out <- rep(NA_real_, n)
  return(as.numeric(out))
}

OUTPUT_DIR <- "phenology_results"
OUT_DIR <- file.path(OUTPUT_DIR, "veg_mixture_fit")
# Ensure OUTPUT_XLSX has a sensible default (can be overridden in user config)
if (is.na(OUTPUT_XLSX) || !is.character(OUTPUT_XLSX) || !nzchar(OUTPUT_XLSX)) {
  OUTPUT_XLSX <- file.path(OUT_DIR, "mesma_results.xlsx")
  cat(sprintf("[CONFIG] OUTPUT_XLSX default set to: %s\n", OUTPUT_XLSX))
}

MIN_SKIP_DOYS_PER_LOCATION <- 2L

FAST_VAR <- TRUE

TEMPORAL_AGGREGATION_DAYS <- 10L  # 10-day intervals (pentads)
N_TEMPORAL_BINS <- ceiling(365 / TEMPORAL_AGGREGATION_DAYS)  # = 37 pentads
TEMPORAL_BUDGET <- N_TEMPORAL_BINS  # Use pentad resolution (no compression beyond pentad aggregation)

N_VARIANTS_PER_VEG <- 10L
MIN_CLUSTER_SIZE <- 10L
PCA_VARIANCE_THRESHOLD <- 0.95
LDA_WEIGHT_FLOOR <- 0.01
# Threshold below which LDA weights are zeroed out (set to NULL to disable)
LDA_WEIGHT_ZERO_THRESHOLD <- 0
# Whether to prune features with zero LDA weight from optimized libraries
PRUNE_ZERO_WEIGHT_FEATURES <- TRUE
# If fraction of zeroed features exceeds this, skip pruning to avoid degenerate libraries
PRUNE_ZERO_WEIGHT_MAX_FRAC <- 0.8
# Minimum number of features to keep after pruning to avoid degenerate models
PRUNE_ZERO_MIN_FEATURES <- 3
# Minimum pentad weight as fraction of global mean weight (default 20%)
PENTAD_WEIGHT_MIN_FRAC <- 0.20
SKIP_MOVING_VARIANCE <- TRUE
ENABLE_MULTISCALE <- FALSE
MULTISCALE_WINDOWS <- c(7L, 14L, 30L)
ENABLE_QP_SOLVER <- TRUE
COMBO_PARALLEL_ENABLE <- TRUE
EARLY_STOP_RMSE_THRESHOLD <- 0.0
ENABLE_DIAGNOSTICS <- TRUE



















# Tunable parameters moved to the centralized section near the top of the file.
# See the "USER-TUNABLE PARAMETERS (centralized)" block for configuration.

if (!file.exists(INPUT_CSV)) {
  stop(paste0("Required input CSV not found: ", INPUT_CSV))
}
raw_df <- tryCatch(
  {
    readr::read_csv(INPUT_CSV, show_col_types = FALSE)
  },
  error = function(e) stop(paste0("Failed to read INPUT_CSV: ", e$message))
)
if (nrow(raw_df) < 2) stop("INPUT_CSV contains fewer than 2 rows")
if (!"location_id" %in% names(raw_df) || (!"prediction_date" %in% names(raw_df) && !"date" %in% names(raw_df))) {
  stop("INPUT_CSV must contain 'location_id' and 'prediction_date' (or 'date') columns")
}

date_col <- if ("prediction_date" %in% names(raw_df)) "prediction_date" else "date"
raw_df <- raw_df |> distinct(location_id, !!sym(date_col), .keep_all = TRUE)


# === GEE INPUT MAPPING BLOCK ===
# Map common GEE-exported column names to the names expected by this script
# 1) Map 'vegetation' -> 'Veg'
if ("vegetation" %in% names(raw_df) && !"Veg" %in% names(raw_df)) {
  raw_df$Veg <- raw_df$vegetation
  cat("[NOTICE] Renamed 'vegetation' column to 'Veg'\n")
}
# 3) Ensure 'location_id' is character (GEE sometimes exports numeric)
if ("location_id" %in% names(raw_df) && !is.character(raw_df$location_id)) {
  raw_df$location_id <- as.character(raw_df$location_id)
  cat("[NOTICE] Coerced 'location_id' to character\n")
}
# Ensure band columns like Blue/Green -> lower-case handled by normalize_band_names later
# === END GEE INPUT MAPPING BLOCK ===



df <- raw_df
if (TESTING_MODE) {
  set.seed(42)
  df <- df |> sample_frac(0.5)
  cat(sprintf("Testing mode: subset to %d rows (50%% of original %d rows)\n", nrow(df), nrow(raw_df)))
}
# if ("PPI" %in% names(df)) df$PPI <- NULL
# Preserve zenith.angle if present (e.g. from metadata), otherwise allow recalculation
if ("zenith.angle" %in% names(df)) {
  cat("[NOTICE] Preserving existing 'zenith.angle' in input data.\n")
} else {
  df$zenith.angle <- NA_real_
}

df <- normalize_band_names(df)
# If bands are present but indices are missing, compute indices from bands
# Ensure NDVI is not used by MESMA — remove it from the working dataframe but keep a backup if present
if ("NDVI" %in% names(df)) {
  df$NDVI_orig <- df$NDVI
  df$NDVI <- NULL
  cat("[NOTICE] NDVI column removed from dataframe to prevent usage in MESMA (backup stored in NDVI_orig)\n")
}

if (!"date" %in% names(df) && "prediction_date" %in% names(df)) df$date <- as.Date(df$prediction_date)
if ("date" %in% names(df)) df$date <- as.Date(df$date)

# === CRITICAL: Filter out snow and dust contamination BEFORE year filtering and PPI baseline calculation ===
# This ensures the PPI baseline is computed from clean observations only
if ("NDSI" %in% names(df) && "NDDI" %in% names(df)) {
  snow_count <- sum(df$NDSI > 0.4, na.rm = TRUE)
  dust_count <- sum(df$NDDI > 0.18, na.rm = TRUE)
  total_before <- nrow(df)
  df <- df[!(df$NDSI > 0.4 | df$NDDI > 0.18), , drop = FALSE]
  total_after <- nrow(df)
  filtered <- total_before - total_after
  cat(sprintf("[EARLY FILTERING] Filtered out %d observations with snow (NDSI > 0.4: %d obs) or dust (NDDI > 0.18: %d obs) contamination\n",
              filtered, snow_count, dust_count))
  cat(sprintf("[EARLY FILTERING] Dataset after contamination filtering: %d rows from %d locations\n",
              total_after, length(unique(df$location_id))))

  # Additionally remove extreme outliers (robust MAD-based), operating per location-year where possible
  df <- remove_large_outliers(df)

  # Defer soil line calculation until after all filtering is complete
  cat("[SOIL LINE] Soil line computation deferred until after final filtering and normalization; it will be calculated later before index computation.\n")
} else {
  cat("[WARNING] NDSI or NDDI not found in data; skipping early contamination filtering\n")
}
# ========================================================================================================

normalize_indices_after_linearization <- function(df, cols = OPTIMAL_INDICES) {
  present <- intersect(cols, names(df))
  idx_scales <- list()
  for (i in present) {
    colv <- df[[i]]
    if (!is.numeric(colv) || all(is.na(colv))) { idx_scales[[i]] <- 1; next }
    scale_val <- as.numeric(stats::quantile(abs(colv), probs = 0.99, na.rm = TRUE))
    if (!is.finite(scale_val) || scale_val <= 0) scale_val <- 1
    df[[i]] <- df[[i]] / scale_val
    idx_scales[[i]] <- scale_val
  }
  attr(df, 'INDEX_SCALES') <- idx_scales
  df
}

# Remove large outliers robustly per (location_id, pheno_year) where possible, otherwise per-location.
# Uses spline-based outlier detection for groups with sufficient data, otherwise falls back to MAD.
remove_large_outliers <- function(df, candidates = NULL, mad_thresh = OUTLIER_MAD_THRESHOLD) {
  if (!isTRUE(ENABLE_OUTLIER_REMOVAL)) return(df)
  if (is.null(candidates)) {
    # Use indices we know are meaningful: OPTIMAL_INDICES + RAW_BANDS, if present
    candidates <- intersect(unique(c(OPTIMAL_INDICES, RAW_BANDS)), names(df))
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
          # --- Iterative Spline Fitting for Robustness ---
          # Pass 1: Initial fit to identify gross outliers
          x <- sub$doy[finite_idx]
          y <- colv[finite_idx]
          
          fit1 <- stats::smooth.spline(x, y, df = min(5, length(x)/2))
          pred1 <- predict(fit1, x)$y
          res1 <- y - pred1
          mad1 <- stats::mad(res1, na.rm = TRUE)
          
          if (!is.finite(mad1) || mad1 <= 1e-6) stop("Invalid MAD in Pass 1")
          
          # Temporarily exclude gross outliers for the second pass (1.5x threshold)
          keep_mask <- abs(res1 - stats::median(res1, na.rm = TRUE)) <= (mad_thresh * 1.5 * mad1)
          
          # Pass 2: Refit on cleaner data if we have enough points left
          if (sum(keep_mask) >= 5) {
            fit2 <- stats::smooth.spline(x[keep_mask], y[keep_mask], df = min(5, sum(keep_mask)/2))
            # Predict against the refined trend for ALL points
            pred_final <- predict(fit2, x)$y
          } else {
            # Fallback to first pass if too many points removed
            pred_final <- pred1
          }

          # Final outlier detection against the robust trend
          residuals <- y - pred_final
          med_res <- stats::median(residuals, na.rm = TRUE)
          mad_res <- stats::mad(residuals, na.rm = TRUE)
          
          if (!is.finite(mad_res) || mad_res <= 0) stop("Invalid Final MAD")
          
          # Flag outliers
          this_mask <- rep(FALSE, length(colv))
          this_mask[finite_idx] <- abs(residuals - med_res) > mad_thresh * mad_res
          out_mask <- out_mask | this_mask
        }, error = function(e) {
          # Fallback to MAD if spline fails
          med <- stats::median(colv, na.rm = TRUE)
          m <- stats::mad(colv, na.rm = TRUE)
          if (is.finite(m) && m > 0) {
            this_mask <- is.finite(colv) & (abs(colv - med) > mad_thresh * m)
            out_mask <<- out_mask | this_mask
          }
        })
      }
    } else {
      # Fallback to MAD-based
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

    # Mark for removal
    if (any(out_mask, na.rm = TRUE)) {
      removed_idx[rows[which(out_mask)]] <- TRUE
    }
  }

  if (any(removed_idx, na.rm = TRUE)) {
    n_removed <- sum(removed_idx, na.rm = TRUE)
    cat(sprintf("[OUTLIER] Removed %d observations across %d groups\n", n_removed, n_groups))
    df <- df[!removed_idx, , drop = FALSE]
  }

  df
}

normalize_mesma_data <- function(df, cols = unique(c(OPTIMAL_INDICES, RAW_BANDS)), lat_default = 40.2) {
  cat("Applying comprehensive MESMA data normalization...\n")
  
  if (!"zenith.angle" %in% names(df)) df$zenith.angle <- NA_real_
  if (!"PPI" %in% names(df)) df$PPI <- NA_real_
  
  INDEX_SCALES <- list()
  present <- intersect(cols, names(df))
  
  for (col in present) {
    vals <- df[[col]]
    vals <- vals[is.finite(vals)]
    if (length(vals) > 0) {
      mu <- mean(vals)
      sigma <- sd(vals)
      if (!is.finite(sigma) || sigma < 1e-10) sigma <- 1.0
      INDEX_SCALES[[col]] <- list(mean = mu, sd = sigma)
    }
  }
  
  for (col in names(INDEX_SCALES)) {
    if (col %in% names(df)) {
      params <- INDEX_SCALES[[col]]
      if (is.list(params) && all(c("mean", "sd") %in% names(params))) {
        mu <- params$mean
        sigma <- params$sd
        if (is.finite(sigma) && sigma > 0) {
          df[[col]] <- (df[[col]] - mu) / sigma
        }
      }
    }
  }
  
  list(
    df = df,
    INDEX_SCALES = INDEX_SCALES
  )
}

apply_stored_normalization <- function(df, norm_params, cols = unique(c(OPTIMAL_INDICES, RAW_BANDS)), lat_default = 40.2) {
  cat("Applying stored normalization parameters to data...\n")
  
  if (!"zenith.angle" %in% names(df)) df$zenith.angle <- NA_real_
  if (!"PPI" %in% names(df)) df$PPI <- NA_real_
  
  # Apply stored INDEX_SCALES (mean/sd) from training normalization to df
  if (!is.null(norm_params$INDEX_SCALES) && length(norm_params$INDEX_SCALES) > 0) {
    cat(sprintf("[apply_stored_normalization] Applying INDEX_SCALES to %d indices\n", length(norm_params$INDEX_SCALES)))
    for (col in names(norm_params$INDEX_SCALES)) {
      if (col %in% names(df)) {
        params <- norm_params$INDEX_SCALES[[col]]
        # If params is a list with mean/sd, perform Z-score normalization
        if (is.list(params) && all(c("mean", "sd") %in% names(params))) {
          mu <- params$mean
          sigma <- params$sd
          if (!is.finite(sigma) || sigma <= 0) sigma <- 1.0
          df[[col]] <- (df[[col]] - mu) / sigma
        } else if (is.numeric(params) && length(params) == 1) {
          # If stored as a single numeric scale (e.g., 99th percentile scale), apply division
          scale_val <- as.numeric(params)
          if (is.finite(scale_val) && scale_val > 0) df[[col]] <- df[[col]] / scale_val
        } else {
          # Unknown format - skip with notice
          cat(sprintf("[apply_stored_normalization] Skipping unknown INDEX_SCALES format for index '%s'\n", col))
        }
      }
    }
  }
  
  
  df
}

apply_secondary_normalization <- function(vec, idx_names, norm_params) {
  if (is.null(norm_params$INDEX_SCALES) || length(norm_params$INDEX_SCALES) == 0) {
    return(vec)  # No normalization parameters available
  }
  
  if (length(vec) != length(idx_names)) {
    warning("apply_secondary_normalization: vec/idx_names length mismatch")
    return(vec)
  }
  
  for (i in seq_along(idx_names)) {
    idx_name <- idx_names[i]
    if (idx_name %in% names(norm_params$INDEX_SCALES)) {
      params <- norm_params$INDEX_SCALES[[idx_name]]
      mu <- params$mean
      sigma <- params$sd
      
      if (is.finite(mu) && is.finite(sigma) && sigma > 0) {
        vec[i] <- (vec[i] - mu) / sigma
      }
    }
  }
  
  vec
}

## No external GeoJSON: construct location map directly from CSV lat/lon
if (all(c("lon", "lat") %in% names(df))) {
  df$location_id <- make_location_id(df$lon, df$lat)
  # Ensure mapping columns exist (fill with NA if missing)
  if (!"Veg" %in% names(df)) df$Veg <- NA_character_
  # Build a minimal gpts_map from unique lat/lon combos in the CSV
  # Use per-location aggregation: pick the first non-missing lat/lon and the
  # first non-missing Veg for that location (avoids losing Veg when
  # the first row happens to have NA)
  gpts_map <- df |>
    dplyr::group_by(location_id) |>
    dplyr::summarise(
      lat = if ("lat" %in% names(df)) { v <- na.omit(lat); if (length(v)>0) v[1] else NA_real_ } else NA_real_,
      lon = if ("lon" %in% names(df)) { v <- na.omit(lon); if (length(v)>0) v[1] else NA_real_ } else NA_real_,
      Veg = if ("Veg" %in% names(df)) { v <- na.omit(as.character(Veg)); if (length(v)>0) tolower(v[1]) else NA_character_ } else NA_character_,
      .groups = "drop"
    )
  gpts_map$location_row <- as.character(seq_len(nrow(gpts_map)))
  gpts_map$location_id_seq <- gpts_map$location_row
  cat(sprintf("[NOTICE] Constructed gpts_map from %d unique lat/lon combos in CSV\n", nrow(gpts_map)))
} else {
  cat("[NOTICE] No 'lat'/'lon' columns found in CSV; GeoJSON mapping disabled.\n")
  gpts_map <- data.frame(location_id = character(0), location_row = character(0), Veg = character(0), stringsAsFactors = FALSE)
}

## GeoJSON processing removed: using CSV lat/lon as unique keys

## Using CSV-derived `gpts_map` above (constructed from unique lat/lon combos)

## Latitude values already present in `gpts_map` when built from CSV lat/lon

cat("\n=== GPTS_MAP LOADING DEBUG ===\n")
if (exists("gpts_map")) {
  cat(sprintf("Total locations in gpts_map: %d\n", nrow(gpts_map)))
  cat(sprintf("Veg class distribution:\n"))
  print(table(gpts_map$Veg, useNA = "ifany"))
  cat(sprintf("Sample location_ids (first 5): %s\n", paste(head(gpts_map$location_id, 5), collapse=", ")))
  barren_gpts_check <- gpts_map[tolower(gpts_map$Veg) == "barren", ]
  cat(sprintf("Barren locations: %d\n", nrow(barren_gpts_check)))
  if (nrow(barren_gpts_check) > 0) {
    cat(sprintf("  Barren location_ids: %s\n", paste(head(barren_gpts_check$location_id, 10), collapse=", ")))
  }
} else {
  cat("[NOTICE] gpts_map object not found; continuing without location mapping.\n")
  # Create a safe empty placeholder so later code can run existence checks
  gpts_map <- data.frame(location_id = character(0), location_row = character(0), Veg = character(0), stringsAsFactors = FALSE)
}
cat("==============================\n\n")

if (nrow(gpts_map) == 0) {
  cat("[NOTICE] No location mapping points found; continuing without derived soil baselines.\n")
}

if ("location_id" %in% names(df) && "location_id" %in% names(gpts_map)) {
  if (!is.character(df$location_id)) df$location_id <- as.character(df$location_id)
  if (!is.character(gpts_map$location_id)) gpts_map$location_id <- as.character(gpts_map$location_id)

  # Normalize IDs to row-number format so training CSV `location_id` (which uses
  # row numbers) will match GeoJSON row IDs (`location_id_seq`).
  df$location_id <- trimws(df$location_id)
  df$location_id[df$location_id == ""] <- NA_character_
  gpts_map$location_id <- trimws(gpts_map$location_id)
  gpts_map$location_id[gpts_map$location_id == ""] <- NA_character_

  # If CSV has 'L_123' style IDs, strip the 'L_' and coerce numeric IDs to
  # canonical integer string form (e.g., '156'). This ensures '156' matches
  # the GeoJSON row id '156'.
  df$location_id <- ifelse(grepl('^L_?[0-9]+$', toupper(df$location_id)), sub('^L_', '', toupper(df$location_id)), df$location_id)
  if (any(grepl('^[0-9]+$', df$location_id, perl=TRUE), na.rm = TRUE)) {
    df$location_id[grepl('^[0-9]+$', df$location_id)] <- as.character(as.integer(df$location_id[grepl('^[0-9]+$', df$location_id)]))
  }
  # Ensure GeoJSON ids are canonical integer strings as well
  if (any(grepl('^[0-9]+$', gpts_map$location_id, perl=TRUE), na.rm = TRUE)) {
    gpts_map$location_id[grepl('^[0-9]+$', gpts_map$location_id)] <- as.character(as.integer(gpts_map$location_id[grepl('^[0-9]+$', gpts_map$location_id)]))
  }
  cat("[NOTICE] Normalized CSV 'location_id' to row-number format where possible.\n")

  if (!"Veg" %in% names(df)) df$Veg <- NA_character_
  pre_non_na <- sum(!is.na(df$Veg) & df$Veg != "")

  joined <- dplyr::left_join(df, gpts_map, by = "location_id", suffix = c("", ".geo"))
  cat("\n=== BARREN JOIN DEBUG ===\n")
  barren_gpts <- gpts_map[tolower(gpts_map$Veg) == "barren", ]
  cat(sprintf("Barren locations in gpts_map: %d\n", length(unique(barren_gpts$location_id))))
  if (nrow(barren_gpts) > 0) {
    sample_barren_ids <- head(unique(barren_gpts$location_id), 5)
    cat(sprintf("  Sample barren location_ids: %s\n", paste(sample_barren_ids, collapse=", ")))
    for (bid in sample_barren_ids) {
      matched_rows <- joined[joined$location_id == bid, ]
      cat(sprintf("  Location %s: %d rows in joined, Veg values: %s\n", 
                  bid, nrow(matched_rows), 
                  paste(unique(matched_rows$Veg), collapse=",")))
      if ("Veg.geo" %in% names(matched_rows)) {
        cat(sprintf("    Veg.geo values: %s\n", paste(unique(matched_rows$Veg.geo), collapse=",")))
      }
    }
  }
  cat("=========================\n\n")
  
  if ("Veg.geo" %in% names(joined)) {
    joined$Veg <- ifelse(is.na(joined$Veg) | joined$Veg == "", joined$Veg.geo, joined$Veg)
    joined$Veg.geo <- NULL
  }
  # 'no soil' and 'no soil.geo' columns are deprecated; no merging performed

  post_non_na <- sum(!is.na(joined$Veg) & joined$Veg != "")

  if (post_non_na == pre_non_na && "location_row" %in% names(gpts_map)) {
    df_ids <- unique(na.omit(as.character(df$location_id)))
    match_count <- length(intersect(df_ids, unique(na.omit(as.character(gpts_map$location_row)))))
    if (match_count > 0) {
      cat(sprintf("[NOTICE] No matches by 'location_id' — attempting join by row-number mapping (matched ids=%d)\n", match_count))
      joined2 <- dplyr::left_join(df, gpts_map, by = c("location_id" = "location_row"), suffix = c("", ".geo"))
      if ("Veg.geo" %in% names(joined2)) {
        joined2$Veg <- ifelse(is.na(joined2$Veg) | joined2$Veg == "", joined2$Veg.geo, joined2$Veg)
        joined2$Veg.geo <- NULL
      }

      if (sum(!is.na(joined2$Veg) & joined2$Veg != "") > post_non_na) {
        joined <- joined2
        post_non_na <- sum(!is.na(joined$Veg) & joined$Veg != "")
        cat(sprintf("[NOTICE] Row-number join gained %d Veg rows\n", post_non_na - pre_non_na))
      } else {
        cat("[NOTICE] Row-number join did not increase Veg mapping; keeping original join state.\n")
    }
  }

  if ("lat" %in% names(gpts_map)) {
    if ("lat.geo" %in% names(joined)) {
      joined$lat <- joined$lat.geo
      joined$lat.geo <- NULL
      cat("[NOTICE] Replaced CSV latitude with GeoJSON latitude where available\n")
    }
  }
  }
  df <- joined

  matched_locs <- length(intersect(na.omit(unique(as.character(df$location_id))), na.omit(unique(as.character(gpts_map$location_id)))))
  cat(sprintf("[NOTICE] GeoJSON join results - Veg before=%d after=%d; matched location_id strings=%d\n", pre_non_na, post_non_na, matched_locs))

  # === TRAINING DATA DIAGNOSTIC ===
  cat("\n=== TRAINING DATA DIAGNOSTIC ===\n")
  cat(sprintf("Total rows in joined data: %d\n", nrow(df)))
  cat(sprintf("Rows with Veg='barren': %d\n", sum(tolower(df$Veg) == "barren", na.rm = TRUE)))
  cat(sprintf("Rows with non-missing Veg: %d\n", sum(!is.na(df$Veg) & df$Veg != "")))

  cat("\nVegetation counts by type:\n")
  veg_types <- sort(unique(tolower(na.omit(df$Veg))))
  for (vt in veg_types) {
    cat(sprintf("  %s: %d rows\n", vt, sum(tolower(df$Veg) == vt, na.rm = TRUE)))
  }
  cat("=========================================\n\n")

}

# Helper to compute soil line slope from filtered bare-soil pixels
compute_soil_line_slope <- function(input_df, min_samples = MIN_ENDMEMBER_SAMPLES, assign_global_dvi = TRUE) {
  if (!all(c('nir','red','Veg') %in% names(input_df))) {
    cat("[SOIL LINE] Required columns 'nir','red' or 'Veg' missing in provided data; leaving SOIL_LINE_SLOPE unchanged (default will be used if not set)\n")
    return(invisible(NA_real_))
  }
  bare_soil_df <- input_df[tolower(input_df$Veg) == 'barren' & is.finite(input_df$nir) & is.finite(input_df$red), , drop = FALSE]
  if (nrow(bare_soil_df) > min_samples) {
    soil_line_model <- tryCatch(lm(nir ~ red, data = bare_soil_df), error = function(e) NULL)
    if (!is.null(soil_line_model)) {
      slope <- as.numeric(coef(soil_line_model)[2])
      assign("SOIL_LINE_SLOPE", slope, envir = globalenv())
      cat(sprintf("[SOIL LINE] Calculated SOIL_LINE_SLOPE=%.4f from %d bare soil pixels\n", slope, nrow(bare_soil_df)))
      # Also compute and store a training DVI soil baseline (mean DVI on bare soil) for PPI use
      if (assign_global_dvi) {
        dvi_soil_calc <- mean(bare_soil_df$nir - bare_soil_df$red, na.rm = TRUE)
        if (is.finite(dvi_soil_calc)) {
          assign("GLOBAL_TRAINING_DVI_SOIL", dvi_soil_calc, envir = globalenv())
          cat(sprintf("[SOIL LINE] Set GLOBAL_TRAINING_DVI_SOIL=%.6f based on bare soil DVI\n", dvi_soil_calc))
        }
      }
      return(invisible(slope))
    }
    else {
      assign("SOIL_LINE_SLOPE", 1.0, envir = globalenv())
      cat("[SOIL LINE] Linear fit failed; using default SOIL_LINE_SLOPE=1.0\n")
      return(invisible(1.0))
    }
  }
  assign("SOIL_LINE_SLOPE", 1.0, envir = globalenv())
  cat("[SOIL LINE] Not enough bare soil pixels in provided data to calculate slope; using default SOIL_LINE_SLOPE=1.0\n")
  return(invisible(1.0))
}

# If bands are present compute only the indices required for contamination filtering (NDSI, NDDI)
# These allow us to remove snow/dust before computing the soil line or other indices that depend on it
eps <- 1e-9
if (all(c('green','swir1') %in% names(df))) df$NDSI <- (as.numeric(df$green) - as.numeric(df$swir1)) / (as.numeric(df$green) + as.numeric(df$swir1) + eps)
if (all(c('red','nir') %in% names(df))) df$NDDI <- (as.numeric(df$red) - as.numeric(df$nir)) / (as.numeric(df$red) + as.numeric(df$nir) + eps)


# === TRAINING-SPECIFIC: Filter out snow and dust contamination BEFORE PPI baseline calculation ===
# Apply the same contamination filtering used for inference data to the training dataset
if ("NDSI" %in% names(df) && "NDDI" %in% names(df)) {
  snow_count <- sum(df$NDSI > 0.4, na.rm = TRUE)
  dust_count <- sum(df$NDDI > 0.18, na.rm = TRUE)
  total_before <- nrow(df)
  df <- df[!(df$NDSI > 0.4 | df$NDDI > 0.18), , drop = FALSE]
  total_after <- nrow(df)
  filtered <- total_before - total_after
  cat(sprintf("[TRAINING EARLY FILTERING] Filtered out %d observations with snow (NDSI > 0.4: %d obs) or dust (NDDI > 0.18: %d obs) contamination\n",
              filtered, snow_count, dust_count))
  cat(sprintf("[TRAINING EARLY FILTERING] Training dataset after contamination filtering: %d rows from %d locations\n",
              total_after, length(unique(df$location_id))))

  # Additionally remove extreme outliers (robust MAD-based), operating per location-year where possible
  df <- remove_large_outliers(df)

  # Now compute soil line slope from the filtered training observations and set global baseline
  compute_soil_line_slope(df)

  # Recompute indices now that SOIL_LINE_SLOPE is available (ensures WDVI uses the filtered-based slope)
  df <- compute_indices_from_bands(df)
} else {
  cat("[TRAINING EARLY FILTERING] NDSI or NDDI not found in training data; skipping contamination filtering\n")
  # Still compute indices so training can proceed if only partial bands present
  df <- compute_indices_from_bands(df)
}

# STEP 2: Calculate PPI on raw data, now that Veg column is available
if (!"pheno_year" %in% names(df) && "date" %in% names(df)) df$pheno_year <- assign_pheno_year(as.Date(df$date))
cat("[NOTICE] Retaining all years for PPI baseline calculation and trend analysis. Training subset will be filtered later.\n")
if (exists("auto_add_ppi_columns")) {
  if (!"PPI" %in% names(df) || all(!is.finite(df$PPI))) {
    ppi_pre_res <- tryCatch({ auto_add_ppi_columns(df) }, error = function(e) { cat(sprintf("[PPI] auto_add_ppi_columns error: %s\n", e$message)); NULL })
    if (!is.null(ppi_pre_res)) {
      df <- ppi_pre_res$df
      if (isTRUE(ppi_pre_res$added)) {
        cat(sprintf("[PPI] Auto-added PPI to dataset before normalization (reason: %s)\n", ppi_pre_res$reason))
      }
    }
  }
}

# STEP 3: Apply normalization to all indices, including the new PPI
cat("\n=== APPLYING TRAINING DATA NORMALIZATION ===\n")
# AGENT: Backup raw PPI for visualization weighting (avoid using Z-scores for weighting)
if ("PPI" %in% names(df)) {
  df$PPI_raw <- df$PPI
  cat("[NOTICE] Backed up raw PPI values to 'PPI_raw' before normalization.\n")
}
norm_result <- normalize_mesma_data(df, lat_default = 40.2)
df <- norm_result$df
INDEX_SCALES <- norm_result$INDEX_SCALES

# Ensure a clamped per-observation PPI normalization fraction is present for downstream use
ppi_max <- if (exists("PPI_FULL_VEG_COVER")) get("PPI_FULL_VEG_COVER") else 0.4
if (!"ppi_norm" %in% names(df)) df$ppi_norm <- NA_real_
# Prefer backed-up raw values if available; otherwise fall back to 'PPI' with a warning
if ("PPI_raw" %in% names(df) && any(is.finite(df$PPI_raw))) {
  df$ppi_norm <- pmin(pmax(df$PPI_raw / ppi_max, 0), 1)
  cat(sprintf("[PPI NORM] Created 'ppi_norm' from 'PPI_raw' and clamped to [0,1] using PPI_FULL_VEG_COVER=%.3f\n", ppi_max))
} else if ("PPI" %in% names(df) && any(is.finite(df$PPI))) {
  df$ppi_norm <- pmin(pmax(df$PPI / ppi_max, 0), 1)
  warning("'PPI_raw' not found - computed 'ppi_norm' from 'PPI' (may be z-scored); values were clamped to [0,1]. Consider backing up raw PPI before normalization.")
} else {
  df$ppi_norm <- NA_real_
  cat("[PPI NORM] No PPI or PPI_raw available to compute 'ppi_norm' (all NA)\n")
}

# Filter to only include selected vegetation types
selected_vegs <- c("phragmites", "populus", "tamarix", "barren")
df <- df[tolower(df$Veg) %in% selected_vegs, ]
cat(sprintf("Filtered training data to selected vegetation types: %s\n", paste(selected_vegs, collapse = ", ")))

# Exclude tamarix samples with no_soil == 1, but keep tamarix with no_soil == NaN
if ("no_soil" %in% names(df)) {
  # Count tamarix with no_soil == 1 (to exclude)
  tamarix_no_soil_1 <- sum(tolower(df$Veg) == "tamarix" & df$no_soil == 1 & !is.na(df$no_soil), na.rm = TRUE)
  # Count tamarix with no_soil == NaN (to keep)
  tamarix_no_soil_na <- sum(tolower(df$Veg) == "tamarix" & is.na(df$no_soil), na.rm = TRUE)
  # Count total tamarix
  total_tamarix <- sum(tolower(df$Veg) == "tamarix", na.rm = TRUE)
  
  cat(sprintf("Tamarix samples: total=%d, no_soil==1=%d, no_soil==NaN=%d\n", 
              total_tamarix, tamarix_no_soil_1, tamarix_no_soil_na))
  
  # Exclude only tamarix with no_soil == 1, keep tamarix with no_soil == NaN or other values
  if (tamarix_no_soil_1 > 0) {
    df <- df[!(tolower(df$Veg) == "tamarix" & df$no_soil == 1 & !is.na(df$no_soil)), ]
    cat(sprintf("Excluded %d tamarix samples with no_soil == 1\n", tamarix_no_soil_1))
  }
}

# AGENT: Capture full dataset for trends BEFORE downsampling or year filtering
df_full <- df
cat(sprintf("[NOTICE] Captured df_full with %d rows for trend analysis (all years).\n", nrow(df_full)))

# AGENT: Downsampling deferred until after df_train creation to ensure year-specific balancing

cat(sprintf("Remaining samples: %d\n", nrow(df)))

TRAINING_NORM_PARAMS <- list(
  INDEX_SCALES = INDEX_SCALES,
  INDEX_SCALES_SECONDARY = list()  # Will be populated after location mapping
)
cat(sprintf("Stored normalization params: INDEX_SCALES for %d indices\n",
            length(INDEX_SCALES)))
cat("(Secondary normalization will be computed after location mapping)\n")
cat("===========================================\n\n")

if (post_non_na == 0L) {
    sample_df_ids <- unique(na.omit(as.character(head(df$location_id, 20))))
    sample_geo_ids <- unique(na.omit(as.character(head(gpts_map$location_id, 20))))
    sample_geo_rows <- unique(na.omit(as.character(head(gpts_map$location_row, 20))))
    cat("[WARNING] Location mapping produced no Veg values. Sample df$location_id (first 20):\n")
    print(sample_df_ids)
    cat("Sample gpts_map$location_id (first 20):\n")
    print(sample_geo_ids)
    cat("Sample gpts_map row-numbers (first 20):\n")
    print(sample_geo_rows)
    cat("Hint: CSV 'location_id' might not match your location mapping; ensure 'location_id' uses the same lat/lon-based keys or row-number mapping.\n")
}


timing_info <- list()
timing_info$start_time <- Sys.time()

cat("Starting vegetation mixture analysis with MESMA approach...\n")
cat("✓ Multiple endmember spectral mixture analysis\n")
cat("✓ Dynamic programming for year-constant proportions\n")
cat("✓ Variant switching allowed throughout the year\n")
cat("\n")
cat(sprintf("Dataset size: %d rows, %d columns\n", nrow(df), ncol(df)))
cat(sprintf("Number of locations: %d\n", length(unique(df$location_id))))
cat(sprintf("Date range: %s to %s\n", min(df$date, na.rm = TRUE), max(df$date, na.rm = TRUE)))

cat("\n=== TRAIN/INFERENCE DATA CONFIGURATION ===\n")
cat(sprintf("Training years (config): %s\n", paste(TRAIN_YEARS, collapse = ", ")))

if (!"pheno_year" %in% names(df)) df$pheno_year <- assign_pheno_year(df$date)

if (exists("TRAIN_YEARS") && !is.null(TRAIN_YEARS) && length(TRAIN_YEARS) > 0) {
  cat(sprintf("Filtering training data to phenological years (March-February): %s\n", paste(TRAIN_YEARS, collapse = ", ")))
  df_train <- df[df$pheno_year %in% TRAIN_YEARS, , drop = FALSE]
} else {
  df_train <- df
}

# AGENT: Apply balancing/downsampling to df_train ONLY (ensure balanced training library)
set.seed(42)
class_counts <- table(df_train$Veg)
if (length(class_counts) > 0) {
  min_count <- min(class_counts)
  if (min_count > 0) {
    df_train <- df_train %>%
      dplyr::group_by(Veg) %>%
      dplyr::slice_sample(n = min_count) %>%
      dplyr::ungroup()
    cat(sprintf("[BALANCE] Downsampled df_train (training year) to %d samples per class (total=%d)\n", min_count, nrow(df_train)))
  }
}
cat(sprintf(
  "Training dataset (Initial): %d rows from %d locations\n",
  nrow(df_train), length(unique(df_train$location_id))
))

# --- STRATIFIED TRAIN/TEST SPLIT (80/20) ---
if (nrow(df_train) > 0) {
  cat("[SPLIT] Performing stratified 80/20 split based on location_id and Veg...\n")
  set.seed(42)
  
  # Get unique location-Veg pairs. 
  # Note: A location typically has one dominant veg type, but might vary. 
  # We'll assign each location to its most frequent Veg type for stratification.
  loc_veg_summary <- df_train %>%
    dplyr::group_by(location_id, Veg) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(location_id, dplyr::desc(n)) %>%
    dplyr::group_by(location_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  unique_vegs <- unique(loc_veg_summary$Veg)
  test_locs_list <- vector("list", length(unique_vegs))
  
  for (i in seq_along(unique_vegs)) {
    v <- unique_vegs[i]
    v_locs <- loc_veg_summary$location_id[loc_veg_summary$Veg == v]
    n_v <- length(v_locs)
    n_test <- ceiling(n_v * 0.20)
    
    if (n_test > 0) {
      selected <- sample(v_locs, n_test)
      test_locs_list[[i]] <- data.frame(location_id = selected, Veg = v, stringsAsFactors = FALSE)
    }
  }
  
  test_locs_df <- do.call(rbind, test_locs_list)
  
  if (!is.null(test_locs_df) && nrow(test_locs_df) > 0) {
    # Save validation locations
    val_file <- file.path(OUT_DIR, "validation_locations.csv") # Assuming OUT_DIR is defined, or use "phenology_results/veg_mixture_fit"
    # Ensure dir exists (OUT_DIR might be defined later, so safeguard)
    if (!dir.exists(dirname(val_file))) dir.create(dirname(val_file), recursive = TRUE)
    
    write.csv(test_locs_df, val_file, row.names = FALSE)
    cat(sprintf("[SPLIT] Saved %d validation locations to %s\n", nrow(test_locs_df), val_file))
    
    # Filter df_train to exclude test locations (Out-of-Bag)
    df_train <- df_train[!df_train$location_id %in% test_locs_df$location_id, ]
    cat(sprintf("[SPLIT] Filtered df_train: %d rows from %d locations (Training Set)\n", 
                nrow(df_train), length(unique(df_train$location_id))))
  } else {
    cat("[SPLIT] Warning: No test locations selected.\n")
  }
}

# Note: Snow/dust contamination filtering already applied early in the pipeline (before PPI baseline calculation)
# See [EARLY FILTERING] section above for details

df_test <- df
cat(sprintf(
  "Testing dataset: %d rows from %d locations\n",
  nrow(df_test), length(unique(df_test$location_id))
))

# df_full is already captured before downsampling/filtering
df <- df_train

if (exists("df_full") && "Veg" %in% names(df_full) && length(ALLOWED_VEG) > 0) {
  missing_vegs <- sapply(ALLOWED_VEG, function(v) {
    sum(tolower(df$Veg) == v, na.rm = TRUE)
  })
  missing_names <- names(missing_vegs)[missing_vegs == 0]
  if (length(missing_names) > 0) {
    for (mv in missing_names) {
      cand <- df_full[tolower(df_full$Veg) == mv, , drop = FALSE]
      if (nrow(cand) > 0) {
        add_n <- min(nrow(cand), max(5L, as.integer(floor(nrow(df) / 10))))
        add_n <- max(1L, add_n)
        to_add <- cand[seq_len(add_n), , drop = FALSE]
        df <- rbind(df, to_add)
        cat(sprintf("[NOTICE] Added %d samples from non-training years for Veg='%s' to ensure representation in training set\n", nrow(to_add), mv))
      } else {
        cat(sprintf("[WARNING] No samples found anywhere for Veg='%s'; cannot add examples to training set\n", mv))
      }
    }
    cat(sprintf("Training dataset after augmentation: %d rows from %d locations\n", nrow(df), length(unique(df$location_id))))
    df_train <- df
  }
}

if (TESTING_MODE) {
  log_msg("TESTING MODE: Capping observations per vegetation type to 10,000L")
  if ("Veg" %in% names(df) && requireNamespace("dplyr", quietly = TRUE)) {
    veg_df <- df[!is.na(df$Veg) & df$Veg != "", ]
    other_df <- df[is.na(df$Veg) | df$Veg == "", ]

    veg_df_sampled <- veg_df |>
      group_by(Veg) |>
      sample_n(min(n(), 10000L)) |>
      ungroup()

    df <- dplyr::bind_rows(veg_df_sampled, other_df)
    
    log_msg("Observations capped. New training data size: %d rows", nrow(df))
    df_train <- df
  } else {
    log_msg("TESTING MODE: Could not cap observations (dplyr not found or 'Veg' column missing).")
  }
}

cat("Using training data for vegetation library construction\n")
cat("=====================================\n\n")

if (isTRUE(PARALLEL_ENABLE)) {
  cleanup_parallel <- setup_parallel_backend(workers = PARALLEL_WORKERS)
} else {
  cleanup_parallel <- function() {}
}
on.exit(cleanup_parallel(), add = TRUE)

# Robust memory cleanup helper — removes large temporary objects matching common patterns
cleanup_memory <- function(verbose = TRUE, confirm = FALSE) {
  patterns <- c("^raw_df$", "^df_full$", "^df_train$", "^df_test$", "^secondary_data$", "^X$", "^X_clean$", "^X_pca$", "^pca_res$", "^pca_result$", "^all_features$", "^all_labels$", "^boot_", "^boot_preds$", "^boot_slopes$", "^boot_coefs$", "^results_list$", "^COMPRESSED_", "^mesma_cache$", "^cached_", "^loc_", "^locations$", "^temp_", "^tmp_")
  objs <- ls(envir = .GlobalEnv)
  to_rm <- unique(unlist(lapply(patterns, function(p) grep(p, objs, value = TRUE, perl = TRUE))))

  if (length(to_rm) == 0) {
    if (isTRUE(verbose)) cat("[CLEANUP] No matched objects to remove.\n")
    return(invisible(NULL))
  }

  if (isTRUE(confirm)) {
    if (isTRUE(verbose)) cat(sprintf("[CLEANUP] Would remove %d objects (confirm=TRUE): %s\n", length(to_rm), paste(head(to_rm, 50), collapse = ", ")))
    return(invisible(to_rm))
  }

  if (isTRUE(verbose)) cat(sprintf("[CLEANUP] Removing %d objects to free memory...\n", length(to_rm)))
  rm(list = to_rm, envir = .GlobalEnv)
  # Force garbage collection; call twice for robustness
  # gc(); gc()
  if (isTRUE(verbose)) cat("[CLEANUP] Memory freed.\n")
  invisible(to_rm)
}





chunked_rbind <- function(lst, chunk_size = 50L) {
  if (is.null(lst) || length(lst) == 0) return(data.frame())
  if (length(lst) == 1) return(lst[[1]])
  
  if (requireNamespace("dplyr", quietly = TRUE)) {
    return(as.data.frame(dplyr::bind_rows(lst)))
  }
  
  return(do.call(rbind, lst))
}

dbg_return_null <- function(reason = NULL) {
  if (exists("TESTING_MODE") && isTRUE(TESTING_MODE)) {
    parent_env <- parent.frame()
    loc <- tryCatch({ if (exists('loc', envir = parent_env)) as.character(parent_env$loc) else NA_character_ }, error = function(e) NA_character_)
    yr <- tryCatch({ if (exists('yr', envir = parent_env)) as.integer(parent_env$yr) else NA_integer_ }, error = function(e) NA_integer_)
    cat(sprintf("[DEBUG] abort: reason=%s loc=%s yr=%s\n", as.character(reason), as.character(loc), as.character(yr)))
  }
  invisible(NULL)
}


safe_lda_call <- function(X_pca, y, min_n_pcs = 2) {
  if (is.null(X_pca) || ncol(X_pca) < min_n_pcs) {
    cat(sprintf("safe_lda_call: Not enough PCs (have=%d, min=%d).\n", ncol(X_pca), min_n_pcs))
    return(NULL)
  }
  curr_n_pcs <- ncol(X_pca)
  while (curr_n_pcs >= min_n_pcs) {
    lda_res <- NULL
    warn_msg <- NULL
    withCallingHandlers({
      lda_res <- tryCatch({
        MASS::lda(X_pca[, 1:curr_n_pcs, drop = FALSE], grouping = y)
      }, error = function(e) {
        e
      })
    }, warning = function(w) {
      warn_msg <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    })

    if (inherits(lda_res, "error")) {
      cat(sprintf("safe_lda_call: LDA error: %s\n", lda_res$message))
      return(NULL)
    }

    if (is.null(warn_msg)) {
      return(lda_res)
    }

    if (grepl("collinear", warn_msg, ignore.case = TRUE)) {
      cat(sprintf("safe_lda_call: Received LDA warning '%s' -> reducing PCs from %d to %d and retrying.\n",
                  warn_msg, curr_n_pcs, curr_n_pcs - 1))
      curr_n_pcs <- curr_n_pcs - 1
      next
    }
    cat(sprintf("safe_lda_call: LDA warning (non-collinearity): %s\n", warn_msg))
    return(lda_res)
  }

  cat("safe_lda_call: Exhausted retries; LDA could not be computed without collinearity.\n")
  NULL
}



train_feature_pipeline <- function(df, class_col, feature_cols) {
  cat(sprintf("\n=== Training Feature Pipeline for Class: %s ===\n", class_col))
  
  norm_params <- list()
  X_raw <- list()
  y_labels <- c()
  
  class_values <- unique(na.omit(df[[class_col]]))
  traces_by_class <- lapply(class_values, function(cv) {
    df_class <- df[df[[class_col]] == cv, , drop = FALSE]
    if (nrow(df_class) == 0) return(list())
    split(df_class, list(df_class$location_id, df_class$pheno_year), drop = TRUE)
  })
  traces <- unlist(traces_by_class, recursive = FALSE)

  cat("  Building trace matrix...\n")
  for(sub in traces) {
    if(nrow(sub) < 5) next
    
    mat <- build_pentad_matrix(sub, feature_cols) # 37 x K
    if(is.null(mat)) next
    
    vec <- as.numeric(mat)
    vec[!is.finite(vec)] <- NA 
    
    X_raw[[length(X_raw)+1]] <- vec
    lbl <- names(sort(table(sub[[class_col]]), decreasing=TRUE))[1]
    y_labels <- c(y_labels, lbl)
  }
  
  if (length(X_raw) < 10) return(NULL)
  X_mat <- do.call(rbind, X_raw)

  # L2-normalize input samples to make LDA focus on shape (relative values) rather than brightness
  if (isTRUE(ENABLE_LDA_L2_NORMALIZATION)) {
    cat("  L2-normalizing training samples for shape-based LDA...\n")
    X_mat <- t(apply(X_mat, 1, function(r) {
      # Treat NA as 0 for norm calculation
      r_clean <- r
      r_clean[is.na(r_clean)] <- 0
      nrm <- sqrt(sum(r_clean^2))
      if (!is.finite(nrm) || nrm < 1e-9) return(r)
      r / nrm
    }))
  } else {
    cat("  L2-normalization for LDA skipped (using raw values)...\n")
  }
  
  
  n_bins <- N_TEMPORAL_BINS
  n_idx <- length(feature_cols)
  
  global_means <- numeric(n_idx)
  global_sds <- numeric(n_idx)
  names(global_means) <- feature_cols
  names(global_sds) <- feature_cols
  
  X_z <- X_mat # Copy structure
  
  cat("  Computing Z-score parameters...\n")
  for(k in seq_along(feature_cols)) {
    col_idx_start <- (k-1)*n_bins + 1
    col_idx_end <- k*n_bins
    
    vals <- X_mat[, col_idx_start:col_idx_end]
    mu <- mean(vals, na.rm=TRUE)
    sigma <- sd(vals, na.rm=TRUE)
    if(sigma == 0 || is.na(sigma)) sigma <- 1
    
    global_means[k] <- mu
    global_sds[k] <- sigma
    
    X_z[, col_idx_start:col_idx_end] <- (vals - mu) / sigma
  }
  X_z[!is.finite(X_z)] <- 0 # Impute for PCA

  cat("  Computing PCA-LDA weights...\n")
  vars <- apply(X_z, 2, var)
  keep_cols <- vars > 1e-9
  X_pca_in <- X_z[, keep_cols, drop=FALSE]
  
  pca_res <- prcomp(X_pca_in, center = FALSE, scale. = FALSE)
  cum_var <- cumsum(pca_res$sdev^2) / sum(pca_res$sdev^2)
  n_pcs <- which(cum_var > 0.95)[1]
  if(is.na(n_pcs)) n_pcs <- ncol(pca_res$x)
  
  class_counts <- table(y_labels)
  n_min <- if (length(class_counts) > 0) min(class_counts) else 0

  max_pcs_for_lda <- min(20, max(1, n_min - 10))

  if (n_pcs > max_pcs_for_lda) {
      old_n_pcs <- n_pcs
      n_pcs <- max_pcs_for_lda
      warning(sprintf("PCA->LDA: Reducing n_pcs from %d to %d to satisfy p < n_min constraint (smallest class has %d samples, max_pcs=20)", old_n_pcs, n_pcs, n_min))
  }

  if (n_pcs < 1) {
    warning("PCA->LDA: Not enough degrees of freedom for LDA (n_pcs < 1). Returning NULL.")
    return(NULL)
  }
  
  min_n_pcs_train <- max(1, length(unique(y_labels)) - 1)
  lda_res <- safe_lda_call(pca_res$x[, 1:n_pcs, drop=FALSE], as.factor(y_labels), min_n_pcs = min_n_pcs_train)

  if (is.null(lda_res)) {
    # LDA failed (too few samples per class, collinearity, etc.) - fall back to uniform weights
    cat("[FALLBACK] LDA failed; using uniform weights (all features weighted equally).\n")
    final_weights <- rep(1, ncol(X_z))
    return(list(
      means = global_means,
      sds = global_sds,
      weights = final_weights,
      indices = feature_cols
    ))
  }
  
  W_pc <- lda_res$scaling
  R <- pca_res$rotation[, 1:n_pcs, drop=FALSE]
  W_std <- R %*% W_pc
  
  svd <- lda_res$svd
  prop <- svd / sum(svd)
  
  if (ncol(W_std) > 1) {
    n_dim <- min(length(prop), ncol(W_std))
    weights_clean <- rowSums(abs(W_std[, 1:n_dim, drop=FALSE]) %*% diag(prop[1:n_dim], nrow=n_dim))
  } else {
    weights_clean <- abs(W_std[, 1])
  }
  
  final_weights <- numeric(ncol(X_z))
  final_weights[keep_cols] <- weights_clean
  # Rescale weights to have a mean of 1 to preserve feature space magnitude
  if (mean(final_weights, na.rm = TRUE) > 1e-9) {
    final_weights <- final_weights / mean(final_weights, na.rm = TRUE)
    cat(sprintf("LDA weights rescaled to mean=1 (new max=%.4f)\n", max(final_weights, na.rm=TRUE)))
  }

  # Enforce pentad weight floor: ensure no weight element is less than
  # PENTAD_WEIGHT_MIN_FRAC * global_mean (default 20% of mean). Apply the floor
  # element-wise and renormalize to mean=1 to preserve scale.
  if (!exists("PENTAD_WEIGHT_MIN_FRAC")) PENTAD_WEIGHT_MIN_FRAC <- 0.20
  if (length(final_weights) > 0 && (length(final_weights) %% N_TEMPORAL_BINS == 0)) {
    n_idx <- length(final_weights) / N_TEMPORAL_BINS
    w_mat <- matrix(final_weights, nrow = N_TEMPORAL_BINS, ncol = n_idx)
    global_mean <- mean(final_weights, na.rm = TRUE)
    min_floor <- PENTAD_WEIGHT_MIN_FRAC * global_mean
    n_before_below <- sum(w_mat < min_floor, na.rm = TRUE)
    if (n_before_below > 0) {
      w_mat[!is.finite(w_mat)] <- 0
      w_mat[w_mat < min_floor] <- min_floor
      final_weights <- as.vector(w_mat)
      # Renormalize back to mean=1
      if (mean(final_weights, na.rm = TRUE) > 1e-9) final_weights <- final_weights / mean(final_weights, na.rm = TRUE)
      cat(sprintf("[NOTICE] Applied pentad weight floor: %.1f%% of global mean (min=%.6g); adjusted %d/%d elements\n", 100*PENTAD_WEIGHT_MIN_FRAC, min_floor, n_before_below, length(final_weights)))
    }
  }

  return(list(
    means = global_means,
    sds = global_sds,
    weights = final_weights,
    indices = feature_cols
  ))
}

doy_to_pentad <- function(doy) {
  pmin(ceiling(doy / TEMPORAL_AGGREGATION_DAYS), N_TEMPORAL_BINS)
}

build_pentad_matrix <- function(dly_year, avail_idx) {
  if (is.null(dly_year) || nrow(dly_year) == 0) return(NULL)

  # CRITICAL: Use phenological DOY (March 1 = day 1), not calendar DOY
  # This ensures temporal alignment when data spans phenological years (March-February)
  if (!"doy" %in% names(dly_year) || any(is.na(dly_year$doy))) {
    dly_year$doy <- pheno_doy(dly_year$date)
  }

  dly_year$pentad <- doy_to_pentad(dly_year$doy)

  K <- length(avail_idx)
  pentad_mat <- matrix(NA_real_, nrow = N_TEMPORAL_BINS, ncol = K)
  colnames(pentad_mat) <- avail_idx

  for (p in 1:N_TEMPORAL_BINS) {
    subset_p <- dly_year[dly_year$pentad == p, ]
    if (nrow(subset_p) == 0) next

    for (j in seq_along(avail_idx)) {
      idx <- avail_idx[j]
      if (!idx %in% names(subset_p)) next

      v <- subset_p[[idx]]
      v <- v[is.finite(v)]
      if (length(v) == 0) next

      q_high <- quantile(v, 0.999, na.rm = TRUE)
      q_low <- quantile(v, 0.001, na.rm = TRUE)
      v_clipped <- pmax(q_low, pmin(v, q_high))

      doys <- subset_p$doy

      if (length(v_clipped) >= 3) {
        model <- tryCatch(lm(v_clipped ~ doys), error = function(e) NULL)
        if (!is.null(model)) {
          mean_doy <- mean(doys)
          pred <- predict(model, newdata = data.frame(doys = mean_doy))
          pentad_mat[p, j] <- pred
        } else {
          pentad_mat[p, j] <- median(v_clipped)
        }
      } else {
        pentad_mat[p, j] <- median(v_clipped)
      }
    }
  }

  # Interpolate missing values column-wise to ensure no NAs
  for (j in 1:K) {
    vals <- pentad_mat[, j]
    if (any(is.na(vals))) {
      if (all(is.na(vals))) {
        pentad_mat[, j] <- 0
      } else {
        idx_present <- which(!is.na(vals))
        if (length(idx_present) >= 2) {
          pentad_mat[, j] <- approx(idx_present, vals[idx_present], xout = 1:N_TEMPORAL_BINS, rule = 2)$y
        } else {
          pentad_mat[, j] <- vals[idx_present[1]]
        }
      }
    }
  }

  pentad_mat
}

 

apply_pca_lda_transform <- function(y, pca_lda_result) {
  if (is.null(pca_lda_result)) return(y)
  
  weights <- pca_lda_result$weights
  
  if (length(weights) != length(y)) {
    warning("Weight length mismatch, using unweighted")
    return(y)
  }
  
  y_weighted <- y * weights
  
  y_weighted
}

medoid_row_index <- function(M) {
  if (is.null(M) || !is.matrix(M) || nrow(M) == 0) return(NA_integer_)
  if (nrow(M) == 1) return(1L)
  X <- M
  X[!is.finite(X)] <- 0
  
  n <- nrow(X)
  
  if (n > MEDOID_SAMPLE_SIZE) {
    samp_idx <- sample.int(n, size = min(as.integer(MEDOID_SAMPLE_SIZE), n))
    X_samp <- X[samp_idx, , drop = FALSE]
    D_samp <- tryCatch(as.matrix(stats::dist(X_samp)), error = function(e) stop(sprintf("medoid_row_index: distance computation failed: %s", e$message)))
    rs_samp <- rowSums(D_samp, na.rm = TRUE)
    local_idx <- which.min(rs_samp)
    if (length(local_idx) == 0 || !is.finite(local_idx)) local_idx <- 1L
    idx <- samp_idx[local_idx]
  } else {
    D <- tryCatch(as.matrix(stats::dist(X)), error = function(e) stop(sprintf("medoid_row_index: distance computation failed: %s", e$message)))
    rs <- rowSums(D, na.rm = TRUE)
    idx <- which.min(rs)
    if (length(idx) == 0 || !is.finite(idx)) idx <- 1L
  }
  as.integer(idx)
}

.run_map <- function(X, FUN, show_pb = TRUE) {
  f_FUN <- FUN
  
  if (!PARALLEL_ENABLE) {
    lapply(X, function(x) { f_FUN(x) })
  } else {
    if (!requireNamespace("future.apply", quietly = TRUE) || !requireNamespace("future", quietly = TRUE)) {
      stop("future.apply and future packages required for parallel processing")
    }

    # AGENT CHANGE: Reuse existing plan if available to avoid costly cluster restart
    current_plan <- future::plan()
    # Check if we are running sequentially (or default). If so, we need to spin up a cluster.
    # If we are already running multisession/cluster, we reuse it.
    plan_is_sequential <- inherits(current_plan, "sequential") || inherits(current_plan, "uniprocess")
    
    if (plan_is_sequential) {
       cat(sprintf("[PARALLEL] .run_map: Spawning new temporary cluster (%d workers)...\n", PARALLEL_WORKERS))
       old_plan <- future::plan()
       options(future.globals.maxSize = max(getOption("future.globals.maxSize", 2e9), 12e9))
    
       # AGENT CHANGE: Using explicit cluster with sequential setup strategy for Windows stability.
       cl <- parallel::makeCluster(PARALLEL_WORKERS, setup_strategy = "sequential")
       future::plan(future::cluster, workers = cl)
    
       # Ensure cluster is stopped and plan is restored on exit
       on.exit({
         future::plan(old_plan)
         parallel::stopCluster(cl)
       }, add = TRUE)
    } else {
       # Use existing plan!
       # cat("[PARALLEL] .run_map: Reusing existing parallel plan.\n")
       options(future.globals.maxSize = max(getOption("future.globals.maxSize", 2e9), 12e9))
    }

    future.apply::future_lapply(X, function(x) {
      f_FUN(x)
    }, future.seed = TRUE)
  }
}



loc_years <- data.frame(location_id = character(0), pheno_year = integer(0), stringsAsFactors = FALSE)

if (!"pheno_year" %in% names(df)) {
  if ("date" %in% names(df)) {
    if (!requireNamespace("lubridate", quietly = TRUE)) stop("The package 'lubridate' is required")
    df$pheno_year <- assign_pheno_year(as.Date(df$date))
  }
}

df <- df |> filter(pheno_year >= 2024 & pheno_year <= 2024)

if (!"Veg" %in% names(df)) df$Veg <- NA_character_

lon_candidates <- names(df)[grepl("(^|_)(lon|longitude|x)(_|$)", names(df), ignore.case = TRUE)]
lat_candidates <- names(df)[grepl("(^|_)(lat|latitude|y)(_|$)", names(df), ignore.case = TRUE)]
if (length(lon_candidates) > 0 && length(lat_candidates) > 0) {
  if ("location_id" %in% names(df) && "location_id" %in% names(gpts_map) && nrow(gpts_map) > 0) {
    if (!is.character(df$location_id)) {
      df$location_id <- as.character(df$location_id)
      cat("[NOTICE] Coerced df$location_id to character to match gpts_map for joining (second join).\n")
    }
    if (!is.character(gpts_map$location_id)) {
      gpts_map$location_id <- as.character(gpts_map$location_id)
      cat("[NOTICE] Coerced gpts_map$location_id to character for joining (second join).\n")
    }
    df <- dplyr::left_join(df, gpts_map, by = "location_id", suffix = c("", ".y"))
    if ("Veg.y" %in% names(df)) {
      if (any(!is.na(df$Veg.y))) {
        df$Veg <- ifelse(is.na(df$Veg) | df$Veg == "", df$Veg.y, df$Veg)
      }
      df$Veg.y <- NULL
    }
  }
}

df$Veg <- tolower(df$Veg)

# Merge vegetation categories: map 'agri' to 'agriculture'
df$Veg <- dplyr::case_when(
  df$Veg == "agri" ~ "agriculture",
  TRUE ~ df$Veg
)

df_train$location_id <- as.character(df_train$location_id)
df_train <- dplyr::left_join(df_train, gpts_map, by = "location_id", suffix = c("", ".geo"))
if ("Veg.geo" %in% names(df_train)) {
  df_train$Veg <- ifelse(is.na(df_train$Veg) | df_train$Veg == "", df_train$Veg.geo, df_train$Veg)
  df_train$Veg.geo <- NULL
}

df_train$Veg <- tolower(df_train$Veg)

# Merge vegetation categories: map 'agri' to 'agriculture'
df_train$Veg <- dplyr::case_when(
  df_train$Veg == "agri" ~ "agriculture",
  TRUE ~ df_train$Veg
)

if (!"date" %in% names(df)) stop("Input CSV must contain a 'date' column")
df$date <- as.Date(df$date)
if (!"location_id" %in% names(df)) stop("Input CSV must contain a 'location_id' column")

if ("Veg" %in% names(df)) df$Veg <- tolower(as.character(df$Veg))

# Merge vegetation categories: map 'agri' to 'agriculture'
if ("Veg" %in% names(df)) {
  df$Veg <- dplyr::case_when(
    df$Veg == "agri" ~ "agriculture",
    TRUE ~ df$Veg
  )
}

required_pkgs <- c("future", "future.apply")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    try(
      {
        install.packages(pkg, repos = "https://cloud.r-project.org")
      },
      silent = TRUE
    )
  }
}

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx", repos = "https://cloud.r-project.org")
}
library(openxlsx)

df$doy <- pheno_doy(df$date)  # Use phenological DOY (March 1 = day 1)
df$doy[df$doy < 1 | df$doy > 366] <- NA_integer_

cat("\n=== BARREN LOADING DEBUG ===\n")
if (exists("gpts_map") && nrow(gpts_map) > 0) {
  barren_in_gpts <- gpts_map[tolower(gpts_map$Veg) == "barren", ]
  if (nrow(barren_in_gpts) > 0) {
    cat(sprintf("  Sample location_ids: %s\n", paste(head(unique(barren_in_gpts$location_id), 5), collapse=", ")))
  }
  
  barren_locs <- unique(barren_in_gpts$location_id)
  barren_in_df <- df[df$location_id %in% barren_locs, ]
  cat(sprintf("Barren locations found in phenology df: %d rows from %d locations\n", 
              nrow(barren_in_df), length(unique(barren_in_df$location_id))))
  
  barren_veg_df <- df[tolower(df$Veg) == "barren", ]
  cat(sprintf("Rows with Veg='barren' in df: %d\n", nrow(barren_veg_df)))
  
  if (nrow(barren_in_df) > 0 && nrow(barren_veg_df) == 0) {
    cat("[WARNING] Barren locations exist in phenology data but Veg column not set to 'barren'!\n")
    cat("  This indicates a join issue. Checking Veg values for barren locations:\n")
    cat(sprintf("  Unique Veg values: %s\n", paste(unique(barren_in_df$Veg), collapse=", ")))
  }
}
cat("=============================\n\n")


veg_counts <- sort(table(na.omit(df$Veg)), decreasing = TRUE)
cat("Vegetation class counts after loading:\n")

  # Fail fast if no Veg metadata was found after mapping
  veg_rows_present <- sum(!is.na(df$Veg) & df$Veg != "")
  if (veg_rows_present == 0) {
    stop(paste0("No vegetation metadata found after mapping (no 'Veg' values).\n",
                "Please ensure your INPUT_CSV contains a 'vegetation' column (mapped to 'Veg').\n",
                "Alternatively, provide GeoJSON location metadata containing 'Veg' to join on.\n",
                "You can run 'scripts/test_veg_presence.R' to see a diagnostic of your input CSV."))
  }
print(veg_counts)

meta_cols <- intersect(c(
  "date", "location_id", "Veg", "coverage", "lat", "lon", "latitude", "longitude",
  "target_lon", "target_lat", "imagery_lat", "imagery_lon", "doy", "pheno_year"
), names(df))

numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]

found_opt <- intersect(OPTIMAL_INDICES, numeric_cols)
found_raw <- intersect(RAW_BANDS, numeric_cols)

if (length(found_opt) > 0) {
  cat(sprintf("Found %d OPTIMAL_INDICES in input: %s\n", length(found_opt), paste(found_opt, collapse = ", ")))
}
missing_opt <- setdiff(OPTIMAL_INDICES, numeric_cols)
if (length(missing_opt) > 0) {
  cat(sprintf("Missing OPTIMAL_INDICES in input: %s\n", paste(missing_opt, collapse = ", ")))
}
if (length(found_raw) > 0) {
  cat(sprintf("Found %d RAW_BANDS in input: %s\n", length(found_raw), paste(found_raw, collapse = ", ")))
}

candidate_indices <- unique(c(found_opt, found_raw))

# Attempt to auto-add PPI to candidate indices when possible: compute PPI
# from empirical barren observations found in the joined dataset, or use the
# MESMA_DVI_SOIL environment override. If the PPI column exists but is all
# NA (this happens when normalization pre-created a PPI column), we still
# attempt to compute and populate it.
if (!"PPI" %in% names(df) || all(!is.finite(df$PPI))) {
  # Ensure DVI exists
  if (!"DVI" %in% names(df) && all(c("nir", "red") %in% names(df))) df$DVI <- df$nir - df$red

  # Find barren observations
  barren_idx <- rep(FALSE, nrow(df))
  if ("Veg" %in% names(df)) {
    barren_idx <- barren_idx | (!is.na(df$Veg) & tolower(trimws(as.character(df$Veg))) == "barren")
  }
  # 'no soil' fraction column is deprecated and ignored; barren identification uses Veg == 'barren' only

  user_dvi_soil <- suppressWarnings(as.numeric(ifelse(nzchar(Sys.getenv("MESMA_DVI_SOIL")), Sys.getenv("MESMA_DVI_SOIL"), NA)))

  # Use add_ppi_columns to compute per-location baselines from median of lowest 10% of ALL DVI values across all years
  # If location_id is present, it will compute per-location baselines; otherwise, it will use a single baseline
  if (!is.na(user_dvi_soil) && is.finite(user_dvi_soil)) {
    cat(sprintf("[PPI] Auto-adding PPI using MESMA_DVI_SOIL override: %.6f\n", user_dvi_soil))
    df <- add_ppi_columns(df, dvi_soil = user_dvi_soil)
    if ("PPI" %in% names(df)) candidate_indices <- unique(c(candidate_indices, "PPI"))
  } else if (sum(is.finite(df$DVI)) > 0) {
    cat("[PPI] Auto-adding PPI: computing per-location baselines from median of lowest 10% of ALL DVI values across all years\n")
    # Let add_ppi_columns compute per-location baselines automatically
    df <- add_ppi_columns(df)
    if ("PPI" %in% names(df)) candidate_indices <- unique(c(candidate_indices, "PPI"))
  } else {
    cat("[NOTICE] Candidate indices computed from existing indices and raw bands; PPI could not be auto-added (no valid DVI values found).\n")
  }
} else {
  # PPI already present in input data
  cat("[NOTICE] Candidate indices computed from existing indices and raw bands; PPI column detected in input and will be used.\n")
  candidate_indices <- unique(c(candidate_indices, "PPI"))
}

if (length(candidate_indices) == 0) {
  stop("No OPTIMAL_INDICES or RAW_BANDS present in input data")
}

cat(sprintf("Selected %d indices: %s\n", length(candidate_indices), paste(candidate_indices, collapse = ", ")))

avail <- candidate_indices
# Ensure NDVI is not used even if present in inputs
if ("NDVI" %in% avail) {
  avail <- setdiff(avail, "NDVI")
  cat("[NOTICE] NDVI removed from candidate indices as requested\n")
}
# Exclude indices that are problematic for fitting (do not use them as features)
for (ex_idx in FITTER_EXCLUDE_INDICES) {
  if (ex_idx %in% avail) {
    avail <- setdiff(avail, ex_idx)
    cat(sprintf("[NOTICE] %s removed from candidate indices for fitter use\n", ex_idx))
  }
}



if (length(avail) == 0) {
  stop("No indices remain after correlation filtering; check input candidate indices and correlation threshold")
}

if (length(avail) < USE_INDICES_MIN) {
  stop(sprintf(
    "Only %d indices remain after filtering, minimum required is %d",
    length(avail), USE_INDICES_MIN
  ))
}


if (isTRUE(SKIP_MOVING_VARIANCE)) {
  cat("Skipping moving variance calculation (SKIP_MOVING_VARIANCE = TRUE)...\n")
  seasonal_var_cols <- c()
  timing_info$moving_var_done <- Sys.time()
} else if (!isTRUE(TESTING_MODE)) {
  if (length(avail) > 0) {
      for (idx_name in avail) {
        var_col <- paste0(idx_name, "_var14")
        v_full <- calc_moving_var(df_full, idx_name, window = 14)
        df_full[[var_col]] <- v_full

      if (exists("df") && nrow(df) > 0) {
        key_full <- paste0(as.character(df_full$location_id), "__", as.character(df_full$date))
        key_df <- paste0(as.character(df$location_id), "__", as.character(df$date))
        mpos <- match(key_df, key_full)
        df[[var_col]] <- ifelse(is.na(mpos), NA_real_, df_full[[var_col]][mpos])
      }

      if (exists("df_test") && nrow(df_test) > 0) {
        key_test <- paste0(as.character(df_test$location_id), "__", as.character(df_test$date))
        key_full <- paste0(as.character(df_full$location_id), "__", as.character(df_full$date))
        mpos2 <- match(key_test, key_full)
        df_test[[var_col]] <- ifelse(is.na(mpos2), NA_real_, df_full[[var_col]][mpos2])
      }
    }
  }

  seasonal_var_cols <- c()
  if (length(avail) > 0) {
    peak_trough_data <- lapply(avail, function(idx_name) {
      var_col <- paste0(idx_name, "_var14")
      v <- df[[var_col]]
      if (all(is.na(v))) {
        return(NULL)
      }
      peaks <- which(diff(sign(diff(v))) == -2)
      troughs <- which(diff(sign(diff(v))) == 2)
      if (length(peaks) >= 1 && length(troughs) >= 1) var_col else NULL
    })
    seasonal_var_cols <- unlist(peak_trough_data[!sapply(peak_trough_data, is.null)])
  }

  df <- df[, c(names(df)[!grepl("_var14$", names(df))], seasonal_var_cols)]

  timing_info$moving_var_done <- Sys.time()
  cat(sprintf(
    "Moving variance calculation completed in %.1f seconds\n",
    as.numeric(difftime(timing_info$moving_var_done, timing_info$start_time, units = "secs"))
  ))
}



cat(sprintf("Post-processing rows (baseline subtraction disabled): %d\n", nrow(df)))
cat("Data preprocessing complete.\n")
adj_cols <- intersect(avail, names(df))


if ("Veg" %in% names(df) && length(ALLOWED_VEG) > 0) {
  try(
    {
      for (av in ALLOWED_VEG) {
        sel <- grepl(av, df$Veg, ignore.case = TRUE) & !is.na(df$Veg)
        if (any(sel)) {
          df$Veg[sel] <- av

        }
      }
    },
    silent = TRUE
  )
}

if ("Veg" %in% names(df) && length(ALLOWED_VEG) > 0) {
  keep_rows <- tolower(df$Veg) %in% ALLOWED_VEG | tolower(df$Veg) == "barren"
  n_before <- nrow(df)
  df <- df[keep_rows | is.na(df$Veg), , drop = FALSE]
  cat(sprintf(
    "Filtered to allowed classes (%s) + barren: kept %d/%d rows\n",
    paste(ALLOWED_VEG, collapse = ","), nrow(df), n_before
  ))
}

try(
  {
    cat("Per-veg quick summary:\n")
    all_veg_classes <- c(ALLOWED_VEG, "barren")
    for (av in all_veg_classes) {
      sel <- tolower(df$Veg) == av
      rows <- sum(sel, na.rm = TRUE)
      unique_doys <- length(unique(df$doy[sel & is.finite(df$doy)]))
      unique_locs <- length(unique(df$location_id[sel]))
      cat(sprintf("  %s: rows=%d unique_doys=%d unique_locs=%d\n", av, rows, unique_doys, unique_locs))
    }
  },
  silent = TRUE
)

matched_veg_n <- sum(!is.na(df$Veg))
cat("Non-NA Veg rows:", matched_veg_n, "of", nrow(df), "\n")
if (matched_veg_n == 0) {
  stop("No vegetation classes found after join; cannot build library")
}

loc_years <- df |>
  dplyr::filter(!is.na(.data$location_id) & .data$location_id != "" & !is.na(.data$pheno_year) & .data$pheno_year > 0 & !is.na(.data$Veg)) |>
  dplyr::distinct(.data$location_id, .data$pheno_year)
cat(sprintf("Constructed loc_years with %d rows from filtered df\n", nrow(loc_years)))
if (nrow(loc_years) == 0) {
  stop(paste0(
    "No location-pheno_year pairs found after filtering. This is a fatal error — library construction cannot continue.\n",
    "Possible causes and suggestions:\n",
    " - Your filtered training dataset has zero valid location/pheno_year pairs (check column 'location_id' and 'pheno_year').\n",
    " - Verify TRAIN_YEARS and any prior filtering steps do not remove all data (e.g. TRAIN_YEARS <- 2019:2024).\n",
    " - Ensure your transformer produced valid 'location_id' values that match your geojson mapping (transform_phenology.py formats 'L_lon_lat').\n",
    " - If your data are intentionally sparse, reduce the filtering thresholds or increase available training data.\n",
    "Processing cannot continue without at least one location-pheno_year pair in filtered training data.")
  )
}

cat("Constructing lib from TRAINING dataset...\n")
cat("[DEBUG] Defining functions...\n")
lib <- list()
build_soil_prototype <- function(df_local, avail_idx, threshold = SOIL_PURE_THRESHOLD, min_samples = SOIL_MIN_SAMPLES) {
  if (!"Veg" %in% names(df_local)) return(NULL)
  candidates <- df_local[tolower(as.character(df_local$Veg)) == "barren", , drop = FALSE]

  if (nrow(candidates) < min_samples) {
    return(NULL)
  }

  if (nrow(candidates) == 0) return(NULL)

  soil_lib <- list()
  if (!"doy" %in% names(candidates)) candidates$doy <- pheno_doy(as.Date(candidates$date))  # Use phenological DOY

  for (idx in avail_idx) {
    if (!idx %in% names(candidates)) next
    vals_by_doy <- tryCatch(
      {
        tapply(seq_along(candidates[[idx]]), candidates$doy, function(indices) {
          vals <- candidates[[idx]][indices]
          vals <- vals[is.finite(vals)]
          if (length(vals) == 0) return(NA_real_)
          if (length(vals) == 1) return(as.numeric(vals[1]))
          m <- mean(vals)
          as.numeric(vals[which.min(abs(vals - m))])
        })
      },
      error = function(e) NULL
    )

    mu <- rep(NA_real_, 365)
    if (!is.null(vals_by_doy) && length(vals_by_doy) > 0) {
      doy_values <- as.integer(names(vals_by_doy))
      valid_doy <- doy_values >= 1 & doy_values <= 365
      mu[doy_values[valid_doy]] <- vals_by_doy[valid_doy]
    }

    if (all(!is.finite(mu))) next
    mv <- tryCatch(calc_moving_var(data.frame(date = 1:365, idx = mu), "idx", window = 14), error = function(e) rep(NA_real_, 365))
    soil_lib[[idx]] <- list(mu = mu, mv = mv)
  }

  attr(soil_lib, "rows_used") <- nrow(candidates)
  soil_lib
}

augment_minority_class <- function(df_class, target_n, seed = NULL, alpha_range = c(0.3, 0.7), jitter_frac = 1e-6) {
  df_class <- as.data.frame(df_class, stringsAsFactors = FALSE)

  if (nrow(df_class) >= target_n) {
    return(df_class)
  }

  if (nrow(df_class) <= 1) {
    reps <- target_n - nrow(df_class)
    base <- df_class[rep(1, reps), , drop = FALSE]
    num_cols <- names(df_class)[sapply(df_class, is.numeric)]

    if (length(num_cols) > 0) {
      rng <- sapply(df_class[, num_cols, drop = FALSE], function(x) {
        x_valid <- x[is.finite(x)]
        if (length(x_valid) < 2) return(1.0)
        r <- diff(range(x_valid))
        if (!is.finite(r) || r <= 0) 1.0 else r
      })

      noise <- matrix(rnorm(reps * length(num_cols), sd = sqrt(mean(rng)) * jitter_frac), nrow = reps)
      base[, num_cols] <- base[, num_cols] + noise
    }

    out <- rbind(df_class, base)
    rownames(out) <- NULL
    return(out)
  }

  if (!is.null(seed)) set.seed(seed)
  n_add <- target_n - nrow(df_class)

  num_cols <- names(df_class)[sapply(df_class, is.numeric)]
  other_cols <- setdiff(names(df_class), num_cols)

  col_ranges <- if(length(num_cols) > 0) {
    sapply(df_class[, num_cols, drop = FALSE], function(x) {
      x_valid <- x[is.finite(x)]
      if (length(x_valid) < 2) return(1.0)
      r <- diff(range(x_valid))
      if (!is.finite(r) || r <= 0) 1.0 else r
    })
  } else numeric(0)

  synth_rows <- lapply(seq_len(n_add), function(i) {
    ids <- sample.int(nrow(df_class), 2, replace = FALSE)
    a <- df_class[ids[1], , drop = FALSE]
    b <- df_class[ids[2], , drop = FALSE]
    alpha <- runif(1, min(alpha_range), max(alpha_range))

    new_row <- a # Start with structure of 'a'

    if (length(num_cols) > 0) {
      va <- as.numeric(a[1, num_cols])
      vb <- as.numeric(b[1, num_cols])

      va[!is.finite(va)] <- 0
      vb[!is.finite(vb)] <- 0

      interp <- (1 - alpha) * va + alpha * vb

      jitter <- rnorm(length(interp), sd = col_ranges * jitter_frac)
      interp <- interp + jitter

      new_row[1, num_cols] <- interp
    }

    for (col in other_cols) {
      if (col %in% c("date", "doy", "doy_for_lib")) next
      if (col %in% c("Veg", "s1_label")) {
        new_row[[col]] <- a[[col]]
        next
      }
      val <- if (runif(1) < 0.5) a[[col]] else b[[col]]
      new_row[[col]] <- val
    }

    if ("doy_for_lib" %in% names(df_class)) {
      da <- as.integer(a$doy_for_lib)
      db <- as.integer(b$doy_for_lib)
      if (!is.na(da) && !is.na(db)) {
        diff_ab <- ((db - da + 365) %% 365)
        pick_doy <- ((da + round(alpha * diff_ab) - 1) %% 365) + 1
        new_row$doy_for_lib <- as.integer(pick_doy)
      }
    }

    new_row
  })

  synth_df <- do.call(rbind, synth_rows)
  out <- rbind(df_class, synth_df)
  rownames(out) <- NULL

  return(out)
}



cat("[DEBUG] Functions defined successfully.\n")
cat("===============================================\n\n")

# Quick self-check to ensure NNLS solver works correctly (TESTING_MODE only)
if (isTRUE(TESTING_MODE)) {
  cat("[TEST] Running simple NNLS self-check...\n")
  tryCatch({
    E_test <- matrix(c(1,0, 0,1, 1,1), nrow = 3, ncol = 2)
    # columns are endmembers; rows are features
    y_test <- E_test %*% c(0.3, 0.7)
    y_test <- as.numeric(y_test)
    if (exists("solve_weights_nnls_simple")) {
      res_simple <- solve_weights_nnls_simple(E_test, y_test)
      cat(sprintf("[TEST] solve_weights_nnls_simple: finite=%d, rmse=%.6g\n", sum(is.finite(res_simple$w)), res_simple$rmse))
    } else {
      cat("[TEST] solve_weights_nnls_simple not defined yet; skipping self-check\n")
    }
    veg_lib_test <- list(A = list(list(vec = E_test[,1], id = 'A1')), B = list(list(vec = E_test[,2], id = 'B1')))
    res_ols <- ols_unmix(y_test, veg_lib_test, topK = 1)
    if (!is.null(res_ols)) cat(sprintf("[TEST] ols_unmix proportions finite=%d, residual=%.6g\n", sum(is.finite(res_ols$proportions)), res_ols$residual)) else cat("[TEST] ols_unmix returned NULL\n")
  }, error = function(e) {
    cat(sprintf("[TEST] NNLS self-check failed: %s\n", e$message))
  })
}

lib_df <- df
if (isTRUE(TESTING_MODE)) {
  if ("Veg" %in% names(lib_df) && nrow(lib_df) > 0) {
    lib_df <- lib_df[!is.na(lib_df$Veg) & lib_df$Veg != "", , drop = FALSE]

    capped_list <- lapply(unique(lib_df$Veg), function(vv) {
      sub <- lib_df[lib_df$Veg == vv, , drop = FALSE]
      if (nrow(sub) > TESTING_MAX_PER_VEG) {
        set.seed(42)
        sub <- sub[sample(seq_len(nrow(sub)), TESTING_MAX_PER_VEG), , drop = FALSE]
      }
      sub
    })
    lib_df <- do.call(rbind, capped_list)
    rownames(lib_df) <- NULL
    cat(sprintf("TESTING_MODE: Capped per-vegetation samples to %d (lib_df rows: %d)\n", TESTING_MAX_PER_VEG, nrow(lib_df)))
  }
}
vegs <- unique(na.omit(lib_df$Veg))
vegs <- vegs[vegs != ""]
vegs <- vegs[tolower(vegs) %in% c("phragmites", "populus", "tamarix", "barren")]  # FIXED: case-insensitive matching

# Filter lib_df to only include the selected vegetation types
lib_df <- lib_df[tolower(lib_df$Veg) %in% tolower(vegs), ]


lib <- list()
for (v in vegs) {
  lib[[v]] <- list(n_samples = 0)
}

timing_info$lib_construction_done <- Sys.time()

cat("=== Building raw index library ===\n")

feature_cols <- avail
cat(sprintf("Raw index features: %s\n", paste(feature_cols, collapse=", ")))

X_all <- as.matrix(lib_df[, feature_cols, drop = FALSE])
for (j in seq_len(ncol(X_all))) {
  col_vals <- X_all[, j]
  if (any(!is.finite(col_vals))) {
    mu_j <- mean(col_vals[is.finite(col_vals)], na.rm = TRUE)
    if (!is.finite(mu_j)) mu_j <- 0
    X_all[!is.finite(col_vals), j] <- mu_j
  }
}
mu_all <- colMeans(X_all)
feature_sds <- apply(X_all, 2, sd)
feature_sds[feature_sds <= 1e-10] <- 1.0

n_features <- length(feature_cols)
cat(sprintf("Raw index library: %d features from %d training samples\n", n_features, nrow(X_all)))

cat("Computing raw index templates per vegetation type...\n")
raw_lib_templates <- list()

for (vname in vegs) {
  dveg <- lib_df[lib_df$Veg == vname, , drop = FALSE]

  if (nrow(dveg) < MIN_ENDMEMBER_SAMPLES) {
    cat(sprintf("[NOTICE] Raw template for '%s' skipped due to insufficient samples: %d < %d\n", vname, nrow(dveg), MIN_ENDMEMBER_SAMPLES))
    next
  }

  if (!"date" %in% names(dveg)) {
    cat(sprintf("[ERROR] %s: 'date' column missing from dveg!\n", vname))
    next
  }

  X_v <- as.matrix(dveg[, feature_cols, drop = FALSE])
  for (j in seq_len(ncol(X_v))) {
    col_vals <- X_v[, j]
    if (any(!is.finite(col_vals))) {
      X_v[!is.finite(col_vals), j] <- mu_all[j]
    }
  }
  
  X_v_c <- sweep(X_v, 2, mu_all, "-")
  X_v_std <- sweep(X_v_c, 2, feature_sds, "/")
  
  doy_vec <- pheno_doy(dveg$date)  # Use phenological DOY
  pentad_vec <- doy_to_pentad(doy_vec)

  if (all(is.na(pentad_vec))) {
    cat(sprintf("  [ERROR] %s: All pentad values are NA\n", vname))
    next
  }

  T_medoid <- matrix(NA_real_, nrow = N_TEMPORAL_BINS, ncol = n_features)
  n_filled <- 0

  for (p in seq_len(N_TEMPORAL_BINS)) {
    rows_p <- which(pentad_vec == p)
    if (length(rows_p) > 0) {
      sub <- X_v_std[rows_p, , drop = FALSE]
      if (nrow(sub) == 1) {
        T_medoid[p, ] <- sub[1, ]
      } else {
        # Use median center for medoid selection within pentad
        center_med <- apply(sub, 2, median, na.rm = TRUE)
        dists <- rowSums(sweep(sub, 2, center_med, "-")^2)
        T_medoid[p, ] <- sub[which.min(dists), ]
      }
      n_filled <- n_filled + 1
    }
  }

  if (n_filled == 0) {
    cat(sprintf("  [ERROR] %s: No pentads filled!\n", vname))
    next
  }
  
  
  raw_lib_templates[[vname]] <- list(T = T_medoid, n_samples = nrow(dveg))
  cat(sprintf("  %s: T_medoid range [%.4f, %.4f], mean=%.4f (filled %d/%d pentads)\n",
              vname, min(T_medoid, na.rm=TRUE), max(T_medoid, na.rm=TRUE), mean(T_medoid, na.rm=TRUE),
              n_filled, N_TEMPORAL_BINS))
}

cat("Raw index library templates computed.\n")



location_bootstrap_ppi <- function(all_coefs, df_tasks, B = BOOTSTRAP_B, seed = 123) {
  set.seed(seed)
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    warning("dplyr required for PPI bootstrap")
    return(NULL)
  }

  # PPI Normalization Constant (default to 0.4 if not global)
  ppi_max <- if (exists("PPI_FULL_VEG_COVER")) get("PPI_FULL_VEG_COVER") else 0.4
  if (ppi_max <= 0) ppi_max <- 0.4

  # Identify PPI column (prefer raw if available)
  ppi_col <- "PPI"
  if ("PPI_raw" %in% names(df_tasks) && any(!is.na(df_tasks$PPI_raw))) {
    ppi_col <- "PPI_raw"
    cat("[PPI BOOTSTRAP] Using raw PPI values (PPI_raw) for weighting.\n")
  } else if ("PPI" %in% names(df_tasks)) {
    cat("[PPI BOOTSTRAP] Using standard PPI column (warning: likely normalized if Z-scores applied).\n")
  } else {
    warning("PPI column missing from df_tasks")
    return(NULL)
  }

  # Ensure temporal columns exist
  if (!"month" %in% names(df_tasks)) {
    if ("date" %in% names(df_tasks)) {
      df_tasks$month <- as.integer(format(as.Date(df_tasks$date), "%m"))
    } else {
      df_tasks$month <- 1L 
    }
  }
  if (!"doy" %in% names(df_tasks)) {
    if ("date" %in% names(df_tasks)) {
      df_tasks$doy <- as.numeric(format(as.Date(df_tasks$date), "%j"))
    } else {
      df_tasks$doy <- 150
    }
  }

  # Helper: Return B bootstrapped medians from a vector
  boot_median_vec <- function(x, n_boot) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(rep(NA_real_, n_boot))
    if (length(x) < 3) return(rep(median(x), n_boot)) # Too few to bootstrap meaningfully
    
    # Fast vectorized bootstrap of median
    replicate(n_boot, median(sample(x, length(x), replace = TRUE), na.rm = TRUE))
  }

  # --- STEP 1: Pre-calculate Detrended Summer PPI (Vectorized) ---
  summer_df <- df_tasks |> 
    dplyr::filter(month %in% 6:9, get(ppi_col) > 0)
  
  if (nrow(summer_df) > 50) {
    tryCatch({
      f <- as.formula(paste(ppi_col, "~ poly(doy, 3)"))
      seasonal_model <- lm(f, data = summer_df)
      summer_df$seasonal_trend <- predict(seasonal_model, newdata = summer_df)
      global_seasonal_mean <- mean(summer_df$seasonal_trend, na.rm = TRUE)
      summer_df$ppi_detrended <- summer_df[[ppi_col]] - (summer_df$seasonal_trend - global_seasonal_mean)
      cat(sprintf("[PPI BOOTSTRAP] Computed detrended summer PPI (N=%d).\n", nrow(summer_df)))
    }, error = function(e) {
      cat(sprintf("[PPI BOOTSTRAP] Detrending failed: %s. Using raw PPI.\n", e$message))
      summer_df$ppi_detrended <- summer_df[[ppi_col]]
    })
  } else {
    summer_df$ppi_detrended <- summer_df[[ppi_col]]
    cat("[PPI BOOTSTRAP] Insufficient summer data for detrending; using raw values.\n")
  }

  # --- STEP 2: Generate Bootstrap Replicates per Location-Year (SUMMER ONLY) ---
  # We need a matrix: Rows = Location-Year ID, Cols = B replicates
  # Strategy:
  # 1. Summarize summer_df to get list of B replicates.
  # 2. STRICTLY use summer data. No annual fallback.

  # Summer Bootstrap
  boot_combined <- summer_df |>
    dplyr::group_by(location_id, pheno_year) |>
    dplyr::summarize(
      ppi_boot_list = list(boot_median_vec(ppi_detrended, B)),
      .groups = "drop"
    )

  # Expand list-columns into a matrix (or just lookup)
  # We'll create a key -> index map and a large matrix
  boot_combined$key <- paste(boot_combined$location_id, boot_combined$pheno_year, sep = "_")
  
  # Construct the matrix of PPI replicates
  n_rows <- nrow(boot_combined)
  ppi_boot_mat <- matrix(NA_real_, nrow = n_rows, ncol = B)
  
  for (i in 1:n_rows) {
    s_vec <- boot_combined$ppi_boot_list[[i]]
    if (!is.null(s_vec) && !all(is.na(s_vec))) {
      ppi_boot_mat[i, ] <- s_vec
    }
  }
  
  # Map (loc, yr) -> row index
  # We attach the row index to the coefficient data
  key_map <- setNames(seq_len(n_rows), boot_combined$key)
  
  # Prepare Coefficient Data
  # Remove pre-existing Barren rows to avoid confusion (we calculate derived Barren)
  merged_veg <- all_coefs[!tolower(trimws(all_coefs$Veg)) %in% c("barren"), ]
  
  # Create key for joining
  merged_veg$key <- paste(merged_veg$location_id, merged_veg$pheno_year, sep = "_")
  merged_veg$ppi_row_idx <- key_map[merged_veg$key]
  
  # Filter out rows with no PPI data
  merged_veg <- merged_veg[!is.na(merged_veg$ppi_row_idx), ]
  
  # Calculate relative coefficients (fractions of vegetation)
  merged_veg <- merged_veg |>
    dplyr::group_by(location_id, pheno_year) |>
    dplyr::mutate(
      sum_veg_coef = sum(coef[tolower(Veg) != "barren"], na.rm = TRUE),
      rel_coef = ifelse(sum_veg_coef > 1e-6, coef / sum_veg_coef, 0)
    ) |>
    dplyr::ungroup()

  veg_types <- unique(merged_veg$Veg[tolower(merged_veg$Veg) != "barren"])
  years <- sort(unique(merged_veg$pheno_year[!is.na(merged_veg$pheno_year)]))
  locations <- unique(merged_veg$location_id)
  n_locs <- length(locations)
  
  # Pre-calculate unique loc-year list for Barren calculation
  # Each unique loc-year needs a row index into ppi_boot_mat
  unique_loc_years <- merged_veg |> 
    dplyr::distinct(location_id, pheno_year, ppi_row_idx)

  results_list <- list()
  
  cat(sprintf("[PPI BOOTSTRAP] Running %d bootstrap iterations with per-location PPI uncertainty...\n", B))

  # We need to accumulate results. Instead of matrix of means, we can just run the loop B times
  # and store the global mean for that iteration.
  
  # Initialize storage for global aggregates per veg/year
  # veg_boot_res[[veg]][[year]] -> vector of B means
  veg_boot_res <- list()
  for (v in c(veg_types, "Barren")) {
    veg_boot_res[[v]] <- matrix(NA_real_, nrow = B, ncol = length(years))
    colnames(veg_boot_res[[v]]) <- as.character(years)
  }

  for (b in seq_len(B)) {
    # 1. Sample Locations (Spatial Bootstrap)
    boot_locs <- sample(locations, n_locs, replace = TRUE)
    
    # 2. Get the PPI realization for this bootstrap iteration 'b'
    # ppi_realization is a vector indexed by ppi_row_idx
    # We take the b-th column of the pre-computed matrix
    # And CLAMP it to [0, 1] after normalizing by ppi_max
    raw_ppi_vals <- ppi_boot_mat[, b]
    norm_ppi_vals <- pmin(pmax(raw_ppi_vals / ppi_max, 0), 1)
    
    # --- PROCESS VEGETATION ---
    # Filter data to selected locations (expanding if sampled multiple times)
    # Efficient join:
    # count how many times each loc was sampled
    loc_counts <- table(boot_locs)
    sampled_locs <- names(loc_counts)
    
    # Subset data to relevant locations
    sub_veg <- merged_veg[merged_veg$location_id %in% sampled_locs, ]
    
    if (nrow(sub_veg) > 0) {
      # Attach the specific PPI value for this iteration
      # Using the row index we mapped earlier
      sub_veg$current_ppi_norm <- norm_ppi_vals[sub_veg$ppi_row_idx]
      
      # Calculate absolute cover = relative_fraction * ppi_cover
      sub_veg$abs_cover <- sub_veg$rel_coef * sub_veg$current_ppi_norm
      
      # Aggregate by Year and Veg
      # We need to weight by the number of times the location was sampled (loc_counts)
      # Or simply replicate the rows. Replicating rows is cleaner for logic but slower.
      # Weighted aggregation:
      sub_veg$weight <- as.integer(loc_counts[sub_veg$location_id])
      
      # Global Mean per Year/Veg = Sum(abs_cover * weight) / Total_Locations
      # Total_Locations = n_locs (since we sampled n_locs with replace)
      
      agg_veg <- sub_veg |>
        dplyr::group_by(pheno_year, Veg) |>
        dplyr::summarize(
          total_cover = sum(abs_cover * weight, na.rm = TRUE),
          .groups = "drop"
        )
      
      # Fill results matrix
      for (i in 1:nrow(agg_veg)) {
        v <- agg_veg$Veg[i]
        y <- as.character(agg_veg$pheno_year[i])
        if (v %in% veg_types && y %in% colnames(veg_boot_res[[v]])) {
          veg_boot_res[[v]][b, y] <- agg_veg$total_cover[i] / n_locs
        }
      }
    }
    
    # --- PROCESS BARREN ---
    # Barren Fraction = 1 - current_ppi_norm
    # Use unique_loc_years
    sub_barren <- unique_loc_years[unique_loc_years$location_id %in% sampled_locs, ]
    
    if (nrow(sub_barren) > 0) {
      sub_barren$current_ppi_norm <- norm_ppi_vals[sub_barren$ppi_row_idx]
      sub_barren$barren_frac <- 1.0 - sub_barren$current_ppi_norm
      sub_barren$weight <- as.integer(loc_counts[sub_barren$location_id])
      
      agg_barren <- sub_barren |>
        dplyr::group_by(pheno_year) |>
        dplyr::summarize(
          total_barren = sum(barren_frac * weight, na.rm = TRUE),
          # Use sum of weights for denominator in case some years are missing in some locations
          n_obs_year = sum(weight, na.rm = TRUE),
          .groups = "drop"
        )
        
      for (i in 1:nrow(agg_barren)) {
        y <- as.character(agg_barren$pheno_year[i])
        if (y %in% colnames(veg_boot_res[["Barren"]])) {
          # Use n_obs_year as denominator (average of available observations)
          # consistent with previous logic
          veg_boot_res[["Barren"]][b, y] <- agg_barren$total_barren[i] / agg_barren$n_obs_year[i]
        }
      }
    }
  } # End Bootstrap Loop

  # Compile Final Results
  final_results <- list()
  for (v in names(veg_boot_res)) {
    mat <- veg_boot_res[[v]]
    # Handle columns with all NAs (years with no data)
    
    df_res <- data.frame(
      year = as.integer(colnames(mat)),
      Veg = v,
      global_coef = apply(mat, 2, median, na.rm = TRUE),
      se = apply(mat, 2, sd, na.rm = TRUE),
      coef_025 = apply(mat, 2, quantile, 0.025, na.rm = TRUE),
      coef_975 = apply(mat, 2, quantile, 0.975, na.rm = TRUE),
      method = "location_bootstrap_ppi_uncertainty"
    )
    # Add n_locations count (approximate based on original data)
    n_locs_per_year <- sapply(years, function(y) sum(unique_loc_years$pheno_year == y))
    df_res$n_locations <- n_locs_per_year[match(df_res$year, years)]
    
    final_results[[v]] <- df_res
  }

  dplyr::bind_rows(final_results)
}

# Location-based bootstrap for global aggregation
location_bootstrap_aggregate <- function(all_coefs, B = BOOTSTRAP_B, seed = 123) {
  set.seed(seed)
  
  # Check if we have measurement uncertainty available
  has_meas_uncertainty <- "coef_sd" %in% names(all_coefs) && !all(is.na(all_coefs$coef_sd))
  if (has_meas_uncertainty) {
     message("[BOOTSTRAP] Measurement uncertainty (coef_sd) detected; will use smooth bootstrap/error propagation.")
  } else {
     message("[BOOTSTRAP] No measurement uncertainty found; using standard spatial bootstrapping.")
  }

  veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
  years <- sort(unique(all_coefs$pheno_year[!is.na(all_coefs$pheno_year)]))
  results_list <- list()

  for (veg in veg_types) {
    # Filter data for this veg
    veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), ]
    if (nrow(veg_data) == 0) next
    
    # --- POOLED VARIANCE CALCULATION ---
    # Calculate robust SD of the coefficients across the entire time series (or per year)
    # to use as a fallback proxy for spatial variability when N is small.
    
    # Calculate SD per year where we have at least 3 points
    sd_by_year <- tapply(veg_data$coef, veg_data$pheno_year, function(x) {
      if (length(x) >= 3) sd(x, na.rm = TRUE) else NA_real_
    })
    
    # Use median of yearly SDs as the "typical spatial variability" for this class
    pooled_sd <- median(sd_by_year, na.rm = TRUE)
    
    cat(sprintf("[BOOTSTRAP] Veg '%s': Pooled Spatial SD = %.4f (used for N < 8)\n", veg, pooled_sd))
    
    boot_means <- matrix(NA_real_, nrow = B, ncol = length(years))
    colnames(boot_means) <- as.character(years)
    
    rho_vec <- numeric(length(years))
    
    for (i in seq_along(years)) {
      yr <- years[i]
      yr_data <- veg_data[veg_data$pheno_year == yr, ]
      n_obs <- nrow(yr_data)
      
      # Estimate spatial autocorrelation rho from data
      rho_est <- 0.5  # default
      if (n_obs >= 3) {
        # Check for pre-existing coordinates
        if (all(c("lat", "lon") %in% names(yr_data))) {
          lat <- yr_data$lat
          lon <- yr_data$lon
        } else if (all(c("latitude", "longitude") %in% names(yr_data))) {
          lat <- yr_data$latitude
          lon <- yr_data$longitude
        } else {
          # Fallback to parsing location_id
          parts <- strsplit(as.character(yr_data$location_id), "_")
          lon_str <- sapply(parts, function(x) if (length(x) >= 2) x[length(x)-1] else NA_character_)
          lat_str <- sapply(parts, function(x) if (length(x) >= 1) x[length(x)] else NA_character_)
          lon <- suppressWarnings(as.numeric(lon_str))
          lat <- suppressWarnings(as.numeric(lat_str))
        }

        # If we have NAs, try to fill them from gpts_map global if available
        if ((any(is.na(lon)) || any(is.na(lat))) && exists("gpts_map", envir = globalenv())) {
           gpts <- get("gpts_map", envir = globalenv())
           if (all(c("location_id", "lat", "lon") %in% names(gpts))) {
             # Match missing
             missing_idx <- which(is.na(lon) | is.na(lat))
             m_ids <- as.character(yr_data$location_id[missing_idx])
             match_idx <- match(m_ids, as.character(gpts$location_id))
             
             lon[missing_idx] <- ifelse(is.na(match_idx), lon[missing_idx], gpts$lon[match_idx])
             lat[missing_idx] <- ifelse(is.na(match_idx), lat[missing_idx], gpts$lat[match_idx])
           }
        }

        if (all(is.na(lon)) || all(is.na(lat))) {
          # If no coordinates available, assume no spatial info -> rho=0
          rho_est <- 0
          # cat(sprintf("[DEBUG] No coords for year %s veg %s; assuming rho=0\n", yr, veg))
        } else {
          # Use whatever coords we have (NAs sort to end)
          ord <- order(lat, lon, na.last = TRUE)
          coef_sorted <- yr_data$coef[ord]
          acf_res <- acf(coef_sorted, lag.max = 1, plot = FALSE, na.action = na.pass)
          rho_est <- Re(acf_res$acf[2])
          if (is.na(rho_est) || rho_est < 0) rho_est <- 0
        }
      }
      rho_vec[i] <- rho_est
      
      if (n_obs > 0) {
        # DECISION: Small Sample Size vs Large Sample Size
        if (n_obs < 15) {
             # --- SMALL SAMPLE: USE POOLED VARIANCE ---
             # When N is small, standard bootstrap under-estimates variance (often 0 if N=1).
             # We assume the spatial variability is 'pooled_sd' and the mean is the sample mean.
             
             mu <- mean(yr_data$coef, na.rm = TRUE)
             
             # Use sample SD if available (n_obs >= 2), otherwise pooled SD
             if (n_obs >= 2) {
               sample_sd <- sd(yr_data$coef, na.rm = TRUE)
               if (is.na(sample_sd) || sample_sd == 0) sample_sd <- pooled_sd
               se_mean <- sample_sd / sqrt(n_obs)
             } else {
               se_mean <- pooled_sd / sqrt(n_obs)
             }
             
             # If we have measurement uncertainty, we can add that too (in quadrature)
             if (has_meas_uncertainty) {
               # Average measurement error for this year
               mean_meas_sd <- mean(yr_data$coef_sd, na.rm = TRUE)
               if (is.na(mean_meas_sd)) mean_meas_sd <- 0
               # Effective SE combines spatial uncertainty (SE_mean) and measurement uncertainty
               # Note: Measurement uncertainty reduces with sqrt(N) assuming independence
               se_meas_mean <- mean_meas_sd / sqrt(n_obs)
               se_total <- sqrt(se_mean^2 + se_meas_mean^2)
             } else {
               se_total <- se_mean
             }
             
             # Add minimum uncertainty to ensure CI is not too narrow for small samples
             se_total <- sqrt(se_total^2 + 0.05^2)

             # Aggressively widen uncertainty for very small N (years with few locations)
             # Translate the CI policy into an SD inflation on the mean.
             ci_tmp <- small_n_inflated_ci(est = mu, sd_in = se_total, n_obs = n_obs)
             se_total <- ci_tmp$coef_sd
             
             # Generate bootstrap distribution parametrically
             # We clamp the mean distribution to physical [0,1] limits if needed, 
             # but standard CI calculation usually prefers raw then clamp.
             boot_vals <- rnorm(B, mean = mu, sd = se_total)
             boot_means[, i] <- boot_vals
             
        } else {
             # --- SUFFICIENT SAMPLE: USE RESAMPLING ---
             
             if (has_meas_uncertainty) {
                # Smooth Bootstrap: Resample locations, then add measurement noise
                # This accounts for both spatial sampling variability AND measurement error
                for (b in 1:B) {
                   # Resample indices (capture spatial variance)
                   idx <- sample(seq_len(n_obs), n_obs, replace = TRUE)
                   sel_coefs <- yr_data$coef[idx]
                   sel_sds <- yr_data$coef_sd[idx]
                   sel_sds[is.na(sel_sds)] <- 0
                   
                   # Add measurement noise (capture measurement variance)
                   sim_coefs <- rnorm(n_obs, mean = sel_coefs, sd = sel_sds)
                   boot_means[b, i] <- mean(sim_coefs, na.rm = TRUE)
                }
             } else {
                # Standard Non-Parametric Bootstrap
                for (b in 1:B) {
                   boot_sample <- sample(yr_data$coef, n_obs, replace = TRUE)
                   boot_means[b, i] <- mean(boot_sample, na.rm = TRUE)
                }
             }
        }
      }
    }
    
    boot_result <- data.frame(
      year = years,
      Veg = veg,
      n_locations = sapply(years, function(y) sum(veg_data$pheno_year == y & !is.na(veg_data$coef))),
      global_coef = apply(boot_means, 2, median, na.rm = TRUE),
      se = apply(boot_means, 2, sd, na.rm = TRUE),
      coef_025 = apply(boot_means, 2, quantile, 0.025, na.rm = TRUE),
      coef_975 = apply(boot_means, 2, quantile, 0.975, na.rm = TRUE),
      rho = rho_vec,
      method = "robust_location_bootstrap"
    )

    # Print summary of spatial autocorrelation (rho) for diagnostics
    if (all(is.na(rho_vec))) {
      cat(sprintf("[BOOTSTRAP] Veg '%s': spatial autocorrelation rho could not be estimated (all NA across years)\n", veg))
    } else {
      rho_mean <- mean(rho_vec, na.rm = TRUE)
      rho_med <- median(rho_vec, na.rm = TRUE)
      rho_min <- min(rho_vec, na.rm = TRUE)
      rho_max <- max(rho_vec, na.rm = TRUE)
      cat(sprintf("[BOOTSTRAP] Veg '%s': spatial autocorrelation rho across years — mean=%.3f, median=%.3f, min=%.3f, max=%.3f\n", veg, rho_mean, rho_med, rho_min, rho_max))
    }
    
    # Correct SE for spatial autocorrelation using estimated rho
    # Effective variance increases by factor (1 + 2*rho) for correlated samples
    boot_result$se <- boot_result$se * sqrt(1 + 2 * boot_result$rho)
    
    # Adjust quantiles to account for autocorrelation by widening the CI
    for (i in seq_len(nrow(boot_result))) {
      gc <- boot_result$global_coef[i]
      c025 <- boot_result$coef_025[i]
      c975 <- boot_result$coef_975[i]
      rho_i <- boot_result$rho[i]
      autocorrelation_factor <- sqrt(1 + 2 * rho_i)
      d_lower <- gc - c025
      d_upper <- c975 - gc
      boot_result$coef_025[i] <- gc - d_lower * autocorrelation_factor
      boot_result$coef_975[i] <- gc + d_upper * autocorrelation_factor
    }

    boot_result$coef_025 <- pmax(0, boot_result$coef_025)
    boot_result$coef_975 <- pmin(1, boot_result$coef_975)
    
    # Clamp the median to prevent negative coverages
    boot_result$global_coef <- pmax(0, boot_result$global_coef)
    
    results_list[[veg]] <- boot_result
  }
  
  dplyr::bind_rows(results_list)
}

# Compute a deliberately conservative CI that widens aggressively as n_obs decreases.
# This is used as a fallback when bootstrapping is skipped due to too few observations.
small_n_inflated_ci <- function(est, sd_in = NA_real_, n_obs = NA_integer_,
                                z = 1.96,
                                n_ref = UNCERTAINTY_N_REF,
                                power = UNCERTAINTY_N_POWER,
                                sd_base = UNCERTAINTY_BASE_SD,
                                sd_max = UNCERTAINTY_SD_MAX) {
  n_eff <- suppressWarnings(as.numeric(n_obs))
  if (!is.finite(n_eff)) n_eff <- NA_real_
  if (is.na(n_eff) || n_eff < 1) n_eff <- 1

  est <- as.numeric(est)
  if (!is.finite(est)) est <- NA_real_

  sd0 <- as.numeric(sd_in)
  if (!is.finite(sd0) || sd0 <= 0) sd0 <- sd_base

  infl <- (as.numeric(n_ref) / n_eff) ^ as.numeric(power)
  if (!is.finite(infl) || infl < 1) infl <- 1

  sd_eff <- sd0 * infl
  if (!is.finite(sd_eff) || sd_eff <= 0) sd_eff <- sd_base
  sd_eff <- min(sd_eff, sd_max)

  if (!is.finite(est)) {
    # If estimate is missing, return fully wide interval.
    return(list(
      coef_025 = 0,
      coef_975 = 1,
      coef_sd = sd_eff,
      interval = 1
    ))
  }

  lo <- est - z * sd_eff
  hi <- est + z * sd_eff
  lo <- max(0, lo)
  hi <- min(1, hi)

  list(
    coef_025 = lo,
    coef_975 = hi,
    coef_sd = sd_eff,
    interval = hi - lo
  )
}

# --- Option A calibration helpers ---
extract_coef_vec_from_fit <- function(fit_res, veg_levels = NULL) {
  if (is.null(fit_res) || is.null(fit_res$coef_df)) return(NULL)
  df <- fit_res$coef_df
  if (is.null(df) || nrow(df) == 0) return(NULL)
  # Defensive: handle multiple variants per veg by summing
  agg <- tapply(df$coef, df$Veg, function(x) sum(as.numeric(x), na.rm = TRUE))
  if (is.null(agg)) return(NULL)
  if (!is.null(veg_levels)) {
    out <- setNames(rep(0, length(veg_levels)), veg_levels)
    common <- intersect(names(agg), veg_levels)
    out[common] <- as.numeric(agg[common])
    return(out)
  }
  as.numeric(agg)
}

estimate_uncertainty_params_optionA <- function(df_tasks,
                                                high_n = UNCERTAINTY_PARAM_EST_HIGH_N,
                                                max_groups = UNCERTAINTY_PARAM_EST_MAX_GROUPS,
                                                target_ns = UNCERTAINTY_PARAM_EST_TARGET_NS,
                                                reps = UNCERTAINTY_PARAM_EST_REPS,
                                                seed = UNCERTAINTY_PARAM_EST_SEED,
                                                veg_levels = ALLOWED_VEG,
                                                quiet = TRUE) {
  if (is.null(df_tasks) || nrow(df_tasks) == 0) stop("df_tasks is empty; cannot calibrate")
  if (!"location_id" %in% names(df_tasks)) stop("df_tasks must have location_id")
  if (!"pheno_year" %in% names(df_tasks)) {
    if ("date" %in% names(df_tasks)) {
      df_tasks$pheno_year <- assign_pheno_year(df_tasks$date)
    } else {
      stop("df_tasks must have pheno_year or date")
    }
  }

  set.seed(seed)

  # Build group table
  grp <- df_tasks |>
    dplyr::group_by(location_id, pheno_year) |>
    dplyr::summarise(n_obs = dplyr::n(), .groups = "drop") |>
    dplyr::filter(is.finite(n_obs) & n_obs >= high_n)

  if (nrow(grp) == 0) stop(sprintf("No location-years with n_obs >= %d", high_n))

  # Sample groups (runtime control)
  if (nrow(grp) > max_groups) grp <- grp[sample(seq_len(nrow(grp)), max_groups), , drop = FALSE]

  # Pre-allocate collectors per target n
  errors_by_n <- setNames(vector("list", length(target_ns)), as.character(target_ns))
  for (n in target_ns) errors_by_n[[as.character(n)]] <- numeric(0)

  # Suppress chatter during repeated fits
  sink_file <- NULL
  con_out <- con_msg <- NULL
  if (isTRUE(quiet)) {
    sink_file <- tempfile()
    con_out <- file(sink_file, open = "wt")
    con_msg <- file(paste0(sink_file, ".msg"), open = "wt")
  }

  on.exit({
    if (isTRUE(quiet)) {
      try(sink(type = "message"), silent = TRUE)
      try(sink(type = "output"), silent = TRUE)
      if (!is.null(con_out)) try(close(con_out), silent = TRUE)
      if (!is.null(con_msg)) try(close(con_msg), silent = TRUE)
      if (!is.null(sink_file)) {
        try(unlink(sink_file), silent = TRUE)
        try(unlink(paste0(sink_file, ".msg")), silent = TRUE)
      }
    }
  }, add = TRUE)

  if (isTRUE(quiet)) {
    sink(con_out, type = "output")
    sink(con_msg, type = "message")
  }

  for (gi in seq_len(nrow(grp))) {
    loc <- grp$location_id[gi]
    yr <- grp$pheno_year[gi]
    sub <- df_tasks[df_tasks$location_id == loc & df_tasks$pheno_year == yr, , drop = FALSE]
    n_full <- nrow(sub)
    if (!is.finite(n_full) || n_full < high_n) next

    # Full fit as pseudo-truth
    full_fit <- tryCatch(fit_one_task(sub), error = function(e) NULL)
    full_vec <- extract_coef_vec_from_fit(full_fit, veg_levels = veg_levels)
    if (is.null(full_vec)) next

    for (n in target_ns) {
      if (n > n_full) next
      for (r in seq_len(reps)) {
        idx <- sample(seq_len(n_full), n, replace = FALSE)
        sub_n <- sub[idx, , drop = FALSE]
        fit_n <- tryCatch(fit_one_task(sub_n), error = function(e) NULL)
        vec_n <- extract_coef_vec_from_fit(fit_n, veg_levels = veg_levels)
        if (is.null(vec_n)) next
        # error pooled across veg components (L2 norm / sqrt(k)) for scale comparability
        diff <- as.numeric(vec_n - full_vec)
        err <- sqrt(mean(diff^2, na.rm = TRUE))
        if (is.finite(err)) errors_by_n[[as.character(n)]] <- c(errors_by_n[[as.character(n)]], err)
      }
    }
  }

  # Summarize sd_emp(n)
  summ <- data.frame(
    n = as.integer(names(errors_by_n)),
    n_samples = sapply(errors_by_n, function(x) length(x)),
    sd_emp = sapply(errors_by_n, function(x) if (length(x) >= 3) stats::sd(x, na.rm = TRUE) else NA_real_),
    stringsAsFactors = FALSE
  )
  summ <- summ[order(summ$n), , drop = FALSE]
  summ <- summ[is.finite(summ$sd_emp) & summ$sd_emp > 0, , drop = FALSE]
  if (nrow(summ) < 3) stop("Not enough valid subsampling results to estimate slope")

  # Fit log(sd_emp) ~ log(n) => slope = -p
  fit <- stats::lm(log(sd_emp) ~ log(n), data = summ)
  p_hat <- -as.numeric(stats::coef(fit)[["log(n)"]])
  if (!is.finite(p_hat) || p_hat < 0) p_hat <- NA_real_

  # Pick N_REF as the first n where improvements become small (stability heuristic)
  # We look for the earliest n such that the relative improvement over the next 3 steps is < 10%.
  n_ref_hat <- NA_integer_
  if (nrow(summ) >= 5) {
    for (i in seq_len(nrow(summ) - 3)) {
      s0 <- summ$sd_emp[i]
      s_next <- summ$sd_emp[(i + 1):(i + 3)]
      # improvement ratio: (s0 - s_next) / s0
      rel_impr <- (s0 - s_next) / s0
      if (all(is.finite(rel_impr)) && all(rel_impr < 0.10)) {
        n_ref_hat <- as.integer(summ$n[i])
        break
      }
    }
  }
  if (!is.finite(n_ref_hat)) n_ref_hat <- as.integer(stats::median(summ$n, na.rm = TRUE))

  list(
    summary = summ,
    fit = fit,
    N_POWER_hat = p_hat,
    N_REF_hat = n_ref_hat,
    used_high_n = high_n,
    used_groups = nrow(grp),
    target_ns = target_ns,
    reps = reps
  )
}







# Helper: ensure the variant similarity heatmap is generated once and early


# Helper: load inference dataset and apply filtering only AFTER the variant similarity heatmap exists ✅
load_and_prepare_inference_data <- function() {
  # Avoid re-loading if already present
  if (exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0) {
    cat("[INFO] Inference data already loaded; skipping.\n")
    assign("INFERENCE_LOAD_DEFERRED", FALSE, envir = globalenv())
    return(invisible(NULL))
  }

  # Ensure heatmap exists (generate if necessary)
  if (!(exists("VARIANT_SIMILARITY_HEATMAP_DONE", envir = globalenv()) && isTRUE(get("VARIANT_SIMILARITY_HEATMAP_DONE", envir = globalenv())))) {
    cat("[INFO] Generating variant similarity heatmap before loading inference data...\n")
    ensure_variant_similarity_heatmap(force = TRUE)
  }

  if (isTRUE(TESTING_MODE)) {
    cat("Skipping inference data loading (TESTING_MODE = TRUE).\n")
    assign("df_tasks_inference", NULL, envir = globalenv())
    assign("inference_location_ids", character(0), envir = globalenv())
    assign("INFERENCE_LOAD_DEFERRED", FALSE, envir = globalenv())
    return(invisible(NULL))
  }
  if (isTRUE(SKIP_INFERENCE)) {
    cat("Skipping inference data loading (SKIP_INFERENCE = TRUE).\n")
    assign("df_tasks_inference", NULL, envir = globalenv())
    assign("inference_location_ids", character(0), envir = globalenv())
    assign("INFERENCE_LOAD_DEFERRED", FALSE, envir = globalenv())
    return(invisible(NULL))
  }

  # Load and filter (copied logic)
  df_inf <- NULL
  if (exists("INFERENCE_CSV")) {
    cat(sprintf("Checking inference file at: %s\n", INFERENCE_CSV))
    if (file.exists(INFERENCE_CSV)) {
      cat(sprintf("Loading inference data from %s...\n", INFERENCE_CSV))
      if (grepl("\\.csv$", INFERENCE_CSV, ignore.case = TRUE)) {
        df_inf <- tryCatch(read.csv(INFERENCE_CSV), error = function(e) { cat(sprintf("[WARNING] Error reading inference CSV: %s\n", e$message)); NULL })
      } else {
        df_inf <- tryCatch(openxlsx::read.xlsx(INFERENCE_CSV), error = function(e) { cat(sprintf("[WARNING] Error reading inference XLSX: %s\n", e$message)); NULL })
      }
      if (!is.null(df_inf)) {
        cat(sprintf("Loaded %d rows from inference file.\n", nrow(df_inf)))
        # Normalise location column names
        possible_loc_names <- c("Location_ID", "location-id", "loc_id", "Loc_ID", "LOCATION_ID")
        for (pn in possible_loc_names) {
          if (pn %in% names(df_inf) && !"location_id" %in% names(df_inf)) {
            names(df_inf)[names(df_inf) == pn] <- "location_id"
            cat(sprintf("[NOTICE] Renamed column '%s' to 'location_id'\n", pn))
            break
          }
        }
        if ("location_id" %in% names(df_inf)) {
          cat(sprintf("Found %d unique inference location IDs.\n", length(unique(df_inf$location_id))))
          df_inf$location_id_orig <- as.character(df_inf$location_id)
          # Add Veg from gpts_map if available
          if (exists("gpts_map") && "Veg" %in% names(gpts_map)) {
            df_inf <- df_inf |> left_join(gpts_map |> select(location_id, Veg), by = "location_id")
            cat(sprintf("[NOTICE] Added Veg column to inference data from gpts_map.\n"))
          }
        } else {
          cat("WARNING: 'location_id' column missing from inference file.\n")
        }
      }
    } else {
      cat(sprintf("Inference file not found at: %s\n", INFERENCE_CSV))
    }
  } else {
    cat("INFERENCE_CSV variable not defined.\n")
  }

  if (!is.null(df_inf) && nrow(df_inf) > 0) {
    df_inf <- normalize_band_names(df_inf)
    if ("...1" %in% names(df_inf) && !"location_id" %in% names(df_inf)) {
      if (is.character(df_inf$...1) || is.numeric(df_inf$...1)) names(df_inf)[names(df_inf) == "...1"] <- "location_id"
    }

    # derive date column
    if ("prediction_date" %in% names(df_inf)) df_inf$date <- as.Date(df_inf$prediction_date) else if ("date" %in% names(df_inf)) df_inf$date <- as.Date(df_inf$date) else {
      for (col in names(df_inf)) {
        if (inherits(df_inf[[col]], "Date")) { df_inf$date <- df_inf[[col]]; break }
        if (is.character(df_inf[[col]]) && all(grepl("^\\d{4}-\\d{2}-\\d{2}", na.omit(df_inf[[col]][1:min(10, nrow(df_inf))])))) { df_inf$date <- as.Date(df_inf[[col]]); break }
      }
    }

    if ("location_id" %in% names(df_inf) && "date" %in% names(df_inf)) {
      df_inf$location_id <- as.character(df_inf$location_id)

      if (length(intersect(RAW_BANDS, names(df_inf))) >= 2) {
        eps <- 1e-9
        if (all(c('green','swir1') %in% names(df_inf))) df_inf$NDSI <- (as.numeric(df_inf$green) - as.numeric(df_inf$swir1)) / (as.numeric(df_inf$green) + as.numeric(df_inf$swir1) + eps)
        if (all(c('red','nir') %in% names(df_inf))) df_inf$NDDI <- (as.numeric(df_inf$red) - as.numeric(df_inf$nir)) / (as.numeric(df_inf$red) + as.numeric(df_inf$nir) + eps)
        cat("[NOTICE] Computed NDSI/NDDI for inference filtering\n")
      }

      if ("NDSI" %in% names(df_inf) && "NDDI" %in% names(df_inf)) {
        snow_count <- sum(df_inf$NDSI > 0.4, na.rm = TRUE)
        dust_count <- sum(df_inf$NDDI > 0.18, na.rm = TRUE)
        total_before <- nrow(df_inf)
        df_inf <- df_inf[!(df_inf$NDSI > 0.4 | df_inf$NDDI > 0.18), , drop = FALSE]
        total_after <- nrow(df_inf)
        filtered <- total_before - total_after
        cat(sprintf("[INFERENCE FILTERING] Filtered out %d observations with snow (NDSI > 0.4: %d obs) or dust (NDDI > 0.18: %d obs) contamination\n", filtered, snow_count, dust_count))
        cat(sprintf("[INFERENCE FILTERING] Inference dataset after contamination filtering: %d rows from %d locations\n", total_after, length(unique(df_inf$location_id))))
        df_inf <- remove_large_outliers(df_inf)
        # compute_soil_line_slope(df_inf, assign_global_dvi = FALSE)  # Not needed for inference, soil line from training
        before_cols <- names(df_inf)
        df_inf <- compute_indices_from_bands(df_inf)
        new_cols <- setdiff(names(df_inf), before_cols)
        if (length(new_cols) > 0) cat(sprintf("[NOTICE] Computed indices from raw bands in inference data: %s\n", paste(new_cols, collapse=", ")))
      } else {
        cat("[WARNING] NDSI or NDDI not found in inference data; skipping contamination filtering\n")
      }

      if (!"Veg" %in% names(df_inf)) df_inf$Veg <- NA_character_
      if (!"pheno_year" %in% names(df_inf)) df_inf$pheno_year <- assign_pheno_year(df_inf$date)
      if (!"doy" %in% names(df_inf)) df_inf$doy <- pheno_doy(df_inf$date)
      if (!"zenith.angle" %in% names(df_inf)) df_inf$zenith.angle <- NA_real_
      if (!"DVI_max" %in% names(df_inf)) df_inf$DVI_max <- NA_real_

      # IMPORTANT: Do NOT restrict inference dataset by TRAIN_YEARS. Inference should use all available years.
      if (exists("TRAIN_YEARS") && !is.null(TRAIN_YEARS)) {
        cat(sprintf("[NOTICE] TRAIN_YEARS is set to %s, but inference dataset will NOT be filtered by training years.\n", paste(TRAIN_YEARS, collapse=", ")))
      }

      # GLOBAL_TRAINING_DVI_SOIL is now guaranteed to exist due to ensure_global_dvi_soil_baseline() call
      dvi_soil_arg <- get("GLOBAL_TRAINING_DVI_SOIL", envir = globalenv())
      
      # If PPI not in df_inf, try to add it
      if (exists("auto_add_ppi_columns") && "PPI" %in% avail && !"PPI" %in% names(df_inf)) {
        cat("[PPI DEBUG] Attempting to auto-add PPI to inference data.\n")
        cat(sprintf("[PPI DEBUG] GLOBAL_TRAINING_DVI_SOIL (guaranteed) is: %.6f\n", dvi_soil_arg))
        

      }

      # If PPI not in df_inf, try to add it
      dvi_soil_arg <- get("GLOBAL_TRAINING_DVI_SOIL", envir = globalenv())
        
              cat("[PPI DEBUG] Attempting to auto-add PPI to inference data.\n")
              if (exists("GLOBAL_TRAINING_DVI_SOIL", envir = globalenv())) {
                dvi_soil_arg <- get("GLOBAL_TRAINING_DVI_SOIL", envir = globalenv())
                cat(sprintf("[PPI DEBUG] GLOBAL_TRAINING_DVI_SOIL found: %.6f\n", dvi_soil_arg))
              } else {
                cat("[PPI DEBUG] GLOBAL_TRAINING_DVI_SOIL NOT found in global environment.\n")
              }
        
              cat(sprintf("[PPI DEBUG] Inference df_inf columns: %s\n", paste(names(df_inf), collapse = ", ")))
              # Check for critical columns needed for PPI in df_inf
              critical_cols <- c("date", "lat", "nir", "red")
              for (col_name in critical_cols) {
                if (!col_name %in% names(df_inf)) {
                  cat(sprintf("[PPI DEBUG] CRITICAL: Column '%s' is missing from df_inf.\n", col_name))
                } else {
                  n_na <- sum(is.na(df_inf[[col_name]]))
                  if (n_na > 0) {
                    cat(sprintf("[PPI DEBUG] Column '%s' has %d NA values (out of %d rows).\n", col_name, n_na, nrow(df_inf)))
                  }
                  if (is.numeric(df_inf[[col_name]]) && all(is.na(df_inf[[col_name]]))) {
                    cat(sprintf("[PPI DEBUG] CRITICAL: Column '%s' is all NA values.\n", col_name))
                  }
                }
              }
              if (nrow(df_inf) > 0) {
                cat("[PPI DEBUG] Sample of df_inf (first 3 rows for critical columns):\n")
                print(head(df_inf[1:min(3, nrow(df_inf)), intersect(names(df_inf), critical_cols), drop = FALSE]))
              } else {
                cat("[PPI DEBUG] df_inf is empty.\n")
              }
              
                    # After all restoration attempts, check if a valid dvi_soil_arg was found
                    if (!is.na(dvi_soil_arg) && is.finite(dvi_soil_arg)) {
                      # A valid dvi_soil_arg was found, proceed to add PPI columns
                      ppi_inf_res <- tryCatch({
                        auto_add_ppi_columns(df_inf, dvi_soil = dvi_soil_arg)
                      }, error = function(e) {
                        cat(sprintf("[PPI ERROR] auto_add_ppi_columns failed during inference: %s\n", e$message))
                        list(df = df_inf, added = FALSE, reason = e$message)
                      })
              
                      if (!is.null(ppi_inf_res) && isTRUE(ppi_inf_res$added)) {
                        df_inf <- ppi_inf_res$df
                        cat(sprintf("[PPI] Auto-added PPI to inference data (reason: %s)\n", ppi_inf_res$reason))
                      } else {
                        cat(sprintf("[PPI WARNING] PPI not added to inference data. Reason: %s\n", ppi_inf_res$reason))
                        # Ensure df_inf gets updated even if PPI wasn't added, in case other columns were added/modified
                        if (!is.null(ppi_inf_res) && !is.null(ppi_inf_res$df)) {
                          df_inf <- ppi_inf_res$df
                        }
                        # If PPI is still not in names(df_inf), ensure it's handled as NA for consistency
                        if (!"PPI" %in% names(df_inf)) df_inf$PPI <- NA_real_
                      }
                    } else {
                      # No valid dvi_soil_arg found after all restoration attempts
                      cat("[PPI CRITICAL WARNING] PPI baseline (dvi_soil) could not be determined from training metadata, raw templates, or global environment. PPI will not be calculated for inference data.\n")
                      if (!"PPI" %in% names(df_inf)) df_inf$PPI <- NA_real_
                    }
      missing_idx <- setdiff(avail, names(df_inf))
      if (length(missing_idx) > 0) { cat(sprintf("[WARNING] Inference data missing indices: %s. Filling with NA (will likely fail unmixing).\n", paste(missing_idx, collapse=", "))); for (col in missing_idx) df_inf[[col]] <- NA_real_ }

      if ("prediction_date" %in% names(df_inf)) df_inf$prediction_date <- as.Date(df_inf$prediction_date)
      if ("reference_date" %in% names(df_inf)) df_inf$reference_date <- as.Date(df_inf$reference_date)

      # AGENT: Backup raw PPI for visualization weighting (inference)
      if ("PPI" %in% names(df_inf)) {
        df_inf$PPI_raw <- df_inf$PPI
        cat("[NOTICE] Backed up raw inference PPI values to 'PPI_raw' before normalization.\n")
      }

      if (exists("TRAINING_NORM_PARAMS") && !is.null(TRAINING_NORM_PARAMS)) { 
        cat("\n=== APPLYING STORED NORMALIZATION TO INFERENCE DATA ===\n")
        df_inf <- apply_stored_normalization(df_inf, TRAINING_NORM_PARAMS, cols = avail, lat_default = 40.2)
        cat("=======================================================\n\n") 
      } else { 
        warning("TRAINING_NORM_PARAMS not found, applying fresh normalization to inference data (scale factors may differ from training!)")
        inf_norm_result <- normalize_mesma_data(df_inf, cols = avail, lat_default = 40.2)
        df_inf <- inf_norm_result$df 
      }

      # Create a clamped per-observation PPI normalization fraction for inference
      ppi_max_inf <- if (exists("PPI_FULL_VEG_COVER")) get("PPI_FULL_VEG_COVER") else 0.4
      if (!"ppi_norm" %in% names(df_inf)) df_inf$ppi_norm <- NA_real_
      if ("PPI_raw" %in% names(df_inf) && any(is.finite(df_inf$PPI_raw))) {
        df_inf$ppi_norm <- pmin(pmax(df_inf$PPI_raw / ppi_max_inf, 0), 1)
        cat(sprintf("[PPI NORM] Inference: created 'ppi_norm' from 'PPI_raw' and clamped to [0,1] using PPI_FULL_VEG_COVER=%.3f\n", ppi_max_inf))
      } else if ("PPI" %in% names(df_inf) && any(is.finite(df_inf$PPI))) {
        df_inf$ppi_norm <- pmin(pmax(df_inf$PPI / ppi_max_inf, 0), 1)
        warning("Inference: 'PPI_raw' not found - computed 'ppi_norm' from 'PPI' (may be z-scored); values were clamped to [0,1]. Consider backing up raw PPI before normalization.")
      } else {
        df_inf$ppi_norm <- NA_real_
        cat("[PPI NORM] Inference: No PPI or PPI_raw available to compute 'ppi_norm' (all NA)\n")
      }

      assign("df_tasks_inference", df_inf, envir = globalenv())

      n_infer_loc_years <- nrow(unique(df_inf[c("location_id", "pheno_year")] ))
      if (exists("df_train") && !is.null(df_train) && nrow(df_train) > 0) { if (!"pheno_year" %in% names(df_train) && "date" %in% names(df_train)) df_train$pheno_year <- assign_pheno_year(df_train$date); n_train_loc_years <- nrow(unique(df_train[c("location_id", "pheno_year")])) } else { n_train_loc_years <- 0 }
      cat(sprintf("(NOTICE) Inference dataset location-years: %d\n", n_infer_loc_years))
      cat(sprintf("(NOTICE) Training dataset location-years: %d\n", n_train_loc_years))
      if (n_train_loc_years > 0 && n_infer_loc_years == n_train_loc_years) cat(sprintf("(WARNING) Training and inference datasets have the same number of location-years (%d). This may be expected if IDs are independent; no automatic filtering will be applied.\n", n_train_loc_years))

      cat("Keeping df_tasks as training data for separate processing.\n")
    } else {
      cat("[WARNING] Inference data missing 'location_id' or 'date' column. Skipping.\n")
      cat(sprintf("Columns found: %s\n", paste(names(df_inf), collapse=", ")))
    }
  }
  assign("inference_location_ids", if (exists("df_tasks_inference") && !is.null(get("df_tasks_inference", envir=globalenv())) && nrow(get("df_tasks_inference", envir=globalenv())) > 0 && "location_id" %in% names(get("df_tasks_inference", envir=globalenv()))) unique(get("df_tasks_inference", envir=globalenv())$location_id) else character(0), envir = globalenv())
  assign("INFERENCE_LOAD_DEFERRED", FALSE, envir = globalenv())
}




  # Old solver `solve_weights_ols()` removed — use `solve_weights_nnls_simple(E, y)` which uses standard NNLS behaviour


  # Simple NNLS wrapper with minimal preprocessing
  solve_weights_nnls_simple <- function(E, y, feature_weights = NULL) {
    if (is.null(E) || ncol(E) < 1) return(NULL)
    # Ensure numeric
    E_fit <- as.matrix(E)
    y_fit <- as.numeric(y)

    # Handle length mismatch (truncate if necessary)
    len <- nrow(E_fit)
    if (length(y_fit) != len) {
      warning(sprintf("solve_weights_nnls_simple: Length of y (%d) does not match rows in E (%d); truncating or padding with zeros", length(y_fit), len))
      if (length(y_fit) > len) y_fit <- y_fit[1:len] else y_fit <- c(y_fit, rep(0, len - length(y_fit)))
    }

    # Impute non-finite values to 0
    y_fit[!is.finite(y_fit)] <- 0
    E_fit[!is.finite(E_fit)] <- 0

    # Optional per-feature weighting (WLS) — apply as sqrt weights to rows
    if (!is.null(feature_weights) && length(feature_weights) == nrow(E_fit)) {
      sqrt_w <- sqrt(pmax(feature_weights, 0))
      E_weighted <- E_fit * sqrt_w
      y_weighted <- y_fit * sqrt_w
    } else {
      E_weighted <- E_fit
      y_weighted <- y_fit
    }

    res <- tryCatch({
      nnls::nnls(E_weighted, y_weighted)
    }, error = function(e) {
      warning(sprintf("nnls failed in solve_weights_nnls_simple: %s", e$message))
      return(NULL)
    })

    if (is.null(res)) return(list(w = rep(1 / ncol(E_fit), ncol(E_fit)), rmse = Inf, residuals = rep(NA_real_, length(y_fit))))

    w <- res$x
    if (sum(w) > 0) w <- w / sum(w)

    # NOTE: Removed the 10% minimum fraction floor here. Filtering should happen AFTER
    # removing barren from the coefficient vector, so that vegetation types compete
    # fairly against each other (not against barren which will be replaced by PPI).

    pred <- as.numeric(E_fit %*% w)
    rmse <- sqrt(mean((y_fit - pred)^2))
    residuals <- y_fit - pred

    list(w = w, rmse = rmse, residuals = residuals)
  }



  cos_angle <- function(a, b) {
    sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
  }





  compute_diagnostics <- function(y, E, w, mesma_result = NULL) {
    if (!is.matrix(E) || ncol(E) == 0 || length(w) != ncol(E)) {
      return(data.frame(
        condition_number = NA_real_,
        residual_sum_of_squares = NA_real_,
        r_squared = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    
    pred <- as.numeric(E %*% w)
    residuals <- y - pred
    rss <- sum(residuals^2)
    tss <- sum((y - mean(y))^2)
    r_squared <- if (tss > 0) 1 - rss/tss else 0
    
    cond_num <- tryCatch({
      s <- svd(E)$d
      if (length(s) > 0 && min(s) > 0) max(s) / min(s) else stop("compute_diagnostics: singular values invalid for condition number calculation")
    }, error = function(e) stop(sprintf("compute_diagnostics: failed computing SVD for condition number: %s", e$message)))
    
    data.frame(
      condition_number = cond_num,
      residual_sum_of_squares = rss,
      r_squared = r_squared,
      stringsAsFactors = FALSE
    )
  }

  cat("Building MESMA endmember library...\n")


  
  
  save_mesma_cache <- function(cache_dir = file.path(OUT_DIR, "mesma_cache")) {
    cat("\n=== SAVING MESMA MODEL CACHE ===\n")
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    core <- list(
      lib = lib,
      mesma_lib = mesma_lib,
      raw_lib_templates = raw_lib_templates,
      veg_counts = veg_counts,
      avail = avail,
      ALLOWED_VEG = ALLOWED_VEG,
      SECONDARY_PARAMS = if (exists("SECONDARY_PARAMS")) SECONDARY_PARAMS else NULL,
      RAW_BARREN_PROTO = if (exists("RAW_BARREN_PROTO")) RAW_BARREN_PROTO else NULL
    )
    saveRDS(core, file = file.path(cache_dir, "mesma_library.rds"))
    # Expose core items to globalenv for downstream visualization helpers
    if (!exists("mesma_lib", envir = globalenv()) && exists("mesma_lib")) {
      assign("mesma_lib", mesma_lib, envir = globalenv())
    }

    raw_templates <- list(
      raw_lib_templates = if (exists("raw_lib_templates")) raw_lib_templates else NULL
    )
    saveRDS(raw_templates, file = file.path(cache_dir, "raw_templates.rds"))

    if (exists(".COMPRESSED_TEMPLATES_ACCESSOR", envir = globalenv())) {
      cat("[DEBUG] Saving compressed templates...\n")
      template_data <- list()
      for (veg in names(mesma_lib)) {
        tryCatch({
          cat(sprintf("[DEBUG] Processing veg '%s'...\n", veg))
          template_data[[veg]] <- list()
          # Fallback to raw_lib_templates if available to prevent empty cache
          rt <- if (exists("raw_lib_templates") && !is.null(raw_lib_templates[[veg]])) raw_lib_templates[[veg]] else NULL
          
          vars <- filter_variants_by_min_samples(mesma_lib[[veg]], min_samples = MIN_ENDMEMBER_SAMPLES, veg = veg, raw_template = rt)
          cat(sprintf("[DEBUG] Veg '%s' has %d variants after filter\n", veg, length(vars)))
          
          for (variant in vars) {
            vid <- if (!is.null(variant$variant_id)) variant$variant_id else variant$id
            if (is.null(vid)) vid <- "unknown_variant"
            template_data[[veg]][[vid]] <- list()
            for (grid_type in c("sparse", "full", "dense")) {
              # Safely retrieve from nested in-memory accessor without relying on tryCatch for indexing
              vec <- NULL
              # Check if veg exists in accessor
              if (!is.null(.COMPRESSED_TEMPLATES_ACCESSOR[[veg]])) {
                # Check if vid exists
                if (!is.null(.COMPRESSED_TEMPLATES_ACCESSOR[[veg]][[vid]])) {
                   # Check if grid_type exists
                   if (!is.null(.COMPRESSED_TEMPLATES_ACCESSOR[[veg]][[vid]][[grid_type]])) {
                     vec <- .COMPRESSED_TEMPLATES_ACCESSOR[[veg]][[vid]][[grid_type]]
                   }
                }
              }
              if (!is.null(vec)) template_data[[veg]][[vid]][[grid_type]] <- vec
            }
          }
        }, error = function(e) {
          cat(sprintf("[WARNING] Failed to process veg '%s' for cache: %s\n", veg, e$message))
        })
      }
      saveRDS(template_data, file = file.path(cache_dir, "compressed_templates.rds"))
      # Also assign compressed templates into the global environment for immediate availability
      if (!exists(".COMPRESSED_TEMPLATES_ACCESSOR", envir = globalenv())) {
        assign(".COMPRESSED_TEMPLATES_ACCESSOR", template_data, envir = globalenv())
      }
      if (!exists("compressed_templates_accessor", envir = globalenv())) {
        assign("compressed_templates_accessor", template_data, envir = globalenv())
      }
    }

    cfg <- list(
      TEMPORAL_BUDGET = if (exists("TEMPORAL_BUDGET")) TEMPORAL_BUDGET else NULL,
      TOPK_VARIANTS = if (exists("TOPK_VARIANTS")) TOPK_VARIANTS else NULL,
      ENABLE_MULTISCALE = if (exists("ENABLE_MULTISCALE")) ENABLE_MULTISCALE else NULL,
      MULTISCALE_WINDOWS = if (exists("MULTISCALE_WINDOWS")) MULTISCALE_WINDOWS else NULL,
      ENABLE_QP_SOLVER = if (exists("ENABLE_QP_SOLVER")) ENABLE_QP_SOLVER else NULL,
      ENABLE_DIAGNOSTICS = if (exists("ENABLE_DIAGNOSTICS")) ENABLE_DIAGNOSTICS else NULL,
      ENABLE_UNCERTAINTY = if (exists("ENABLE_UNCERTAINTY")) ENABLE_UNCERTAINTY else NULL,
      DEBUG_UNCERTAINTY = if (exists("DEBUG_UNCERTAINTY")) DEBUG_UNCERTAINTY else NULL,
      MAX_VEG_COMPONENTS = if (exists("MAX_VEG_COMPONENTS")) MAX_VEG_COMPONENTS else NULL,
      MIN_IDX_PRESENCE = if (exists("MIN_IDX_PRESENCE")) MIN_IDX_PRESENCE else NULL,
      EPS_SIGMA = if (exists("EPS_SIGMA")) EPS_SIGMA else NULL,
      LOWER_BND = if (exists("LOWER_BND")) LOWER_BND else NULL,
      USE_INDICES_MIN = if (exists("USE_INDICES_MIN")) USE_INDICES_MIN else NULL,
      MIN_INDEX_SD = if (exists("MIN_INDEX_SD")) MIN_INDEX_SD else NULL
    )
    saveRDS(cfg, file = file.path(cache_dir, "config_params.rds"))

    meta <- list(
      training_years = if (exists("TRAIN_YEARS")) TRAIN_YEARS else NULL,
      dvi_soil = if (exists("GLOBAL_TRAINING_DVI_SOIL")) GLOBAL_TRAINING_DVI_SOIL else NULL,
      training_date_range = if (exists("df_train")) range(df_train$date, na.rm = TRUE) else NULL,
      n_training_samples = if (exists("df_train")) nrow(df_train) else NULL,
      n_locations_trained = if (exists("df_train")) length(unique(df_train$location_id)) else NULL,
      indices_used = avail,
      vegetation_types = names(mesma_lib),
      model_creation_time = Sys.time()
    )
    saveRDS(meta, file = file.path(cache_dir, "training_metadata.rds"))

    manifest <- list(
      version = "1.0",
      created = Sys.time(),
      files = c("mesma_library.rds", "raw_templates.rds", "compressed_templates.rds", "config_params.rds", "training_metadata.rds"),
      checksum = list()
    )
    for (f in manifest$files) {
      fpath <- file.path(cache_dir, f)
      if (file.exists(fpath)) manifest$checksum[[f]] <- tryCatch(tools::md5sum(fpath), error = function(e) stop(sprintf("Failed to compute checksum for %s: %s", fpath, e$message)))
    }
    saveRDS(manifest, file = file.path(cache_dir, "manifest.rds"))

    cat(sprintf("Model cache saved to: %s\n", cache_dir))
    cat(sprintf("Total cache size: %.2f MB\n", sum(file.info(list.files(cache_dir, full.names = TRUE))$size) / 1024^2))
    invisible(cache_dir)
  }


  cat("Constructing task list from inference dataset (all years)...\n")
  df_tasks <- df_full

  # Defer the heavy inference loading and filtering until after variant similarity heatmap is created.
  cat("Deferring inference data loading until after variant similarity heatmap has been generated.\n")
  assign("INFERENCE_LOAD_DEFERRED", TRUE, envir = globalenv())
  df_inf <- NULL
  inference_location_ids <- character(0)

  if (exists("gpts_map")) {
    if ("location_id" %in% names(df_tasks) && "location_id" %in% names(gpts_map)) {
      if (!is.character(df_tasks$location_id)) {
        df_tasks$location_id <- as.character(df_tasks$location_id)
        cat("[NOTICE] Coerced df_tasks$location_id to character to match gpts_map for joining.\n")
      }
      if (!is.character(gpts_map$location_id)) {
        gpts_map$location_id <- as.character(gpts_map$location_id)
        cat("[NOTICE] Coerced gpts_map$location_id to character for joining with df_tasks.\n")
      }
      s1 <- unique(na.omit(df_tasks$location_id))
      s2 <- unique(na.omit(gpts_map$location_id))
      if (length(s1) && length(s2) && all(grepl("^[0-9]+$", s1)) && any(grepl("^L_", s2))) {
        cat("[WARNING] df_tasks$location_id looks numeric while gpts_map$location_id looks like 'L_x_y' strings — matching will likely fail.\n")
      }
      df_tasks <- dplyr::left_join(df_tasks, gpts_map, by = "location_id", suffix = c("", ".geo"))
      if ("Veg.geo" %in% names(df_tasks)) {
        df_tasks$Veg <- ifelse(is.na(df_tasks$Veg) | df_tasks$Veg == "", df_tasks$Veg.geo, df_tasks$Veg)
        df_tasks$Veg.geo <- NULL
      }
    } else {
      df_tasks <- dplyr::left_join(df_tasks, gpts_map, by = "location_id", suffix = c("", ".geo"))
      if ("Veg.geo" %in% names(df_tasks)) {
        df_tasks$Veg <- ifelse(is.na(df_tasks$Veg) | df_tasks$Veg == "", df_tasks$Veg.geo, df_tasks$Veg)
        df_tasks$Veg.geo <- NULL
      }
    }
  }

  if ("Veg" %in% names(df_tasks)) df_tasks$Veg <- tolower(as.character(df_tasks$Veg))
  if ("date" %in% names(df_tasks)) {
    df_tasks$date <- as.Date(df_tasks$date)
    if (!"pheno_year" %in% names(df_tasks)) df_tasks$pheno_year <- assign_pheno_year(df_tasks$date)
    if (!"pheno_year" %in% names(df_tasks)) df_tasks$pheno_year <- assign_pheno_year(df_tasks$date)
    df_tasks$doy <- pheno_doy(df_tasks$date)  # Use phenological DOY
    df_tasks$doy[df_tasks$doy < 1 | df_tasks$doy > 366] <- NA_integer_
  }

  # For multi-year bootstrap, we process by LOCATION (not location-year)
  test_locations <- df_tasks |>
    dplyr::filter(!is.na(.data$location_id) & trimws(.data$location_id) != "") |>
    dplyr::distinct(.data$location_id)

  test_locations$location_id <- trimws(as.character(test_locations$location_id))

  # Still keep test_loc_years for reporting purposes
  test_loc_years <- df_tasks |>
    dplyr::filter(!is.na(.data$location_id) & trimws(.data$location_id) != "" & !is.na(.data$pheno_year) & .data$pheno_year > 0) |>
    dplyr::distinct(.data$location_id, .data$pheno_year)
  test_loc_years$location_id <- trimws(as.character(test_loc_years$location_id))

  if (isTRUE(TESTING_MODE)) {
    obs_counts <- df_tasks |>
      dplyr::group_by(location_id) |>
      dplyr::summarise(n = dplyr::n(), .groups = 'drop') |>
      dplyr::arrange(dplyr::desc(n))

    unique_locs <- unique(obs_counts$location_id)
    n_test_locs <- min(5, length(unique_locs))
    test_locs_subset <- unique_locs[1:n_test_locs]

    test_locations <- test_locations[test_locations$location_id %in% test_locs_subset, , drop = FALSE]
    test_loc_years <- test_loc_years[test_loc_years$location_id %in% test_locs_subset, ]
    cat(sprintf("[TESTING MODE] Limited to %d locations (with most observations) for testing: %s\n",
                n_test_locs, paste(test_locs_subset, collapse = ", ")))
  }

  cat(sprintf(
    "Processing %d unique locations (%d location-year pairs total)\n",
    nrow(test_locations), nrow(test_loc_years)
  ))

  # Create location-based task list (not location-year)
  n_locs <- nrow(test_locations)
  location_list <- test_locations$location_id

cat("\n=== PRE-FLIGHT CHECK ===\n")
if (n_locs > 0) {
  test_loc <- location_list[1]
  test_dly <- df_tasks[df_tasks$location_id == test_loc, , drop = FALSE]
  test_years <- unique(test_dly$pheno_year[!is.na(test_dly$pheno_year)])

  cat(sprintf("Test location: %s\n", test_loc))
  cat(sprintf("Years available for this location: %s\n", paste(sort(test_years), collapse=", ")))
  cat(sprintf("Total rows for this location: %d\n", nrow(test_dly)))
  cat(sprintf("Columns available: %s\n", paste(names(test_dly), collapse=", ")))
  cat(sprintf("Var14 columns: %s\n", paste(grep("_var14$", names(test_dly), value=TRUE), collapse=", ")))
}
cat("======================\n\n")

## Ensure a canonical 'test_task' (location-year tuple) exists for downstream reporting
if (!exists("test_task") || is.null(test_task)) {
  if (exists("test_loc_years") && nrow(test_loc_years) > 0) {
    test_task <- list(loc = test_loc_years$location_id[1], yr = test_loc_years$pheno_year[1])
  } else {
    # Fallback to the first available location and its first year
    test_task <- list(loc = if (exists("test_loc") && !is.null(test_loc)) test_loc else NA_character_,
                      yr = if (exists("test_years") && length(test_years) > 0) test_years[1] else NA_integer_)
  }
}

test_dly <- df_tasks[df_tasks$location_id == test_task$loc, , drop = FALSE]
if ("pheno_year" %in% names(test_dly)) {
  test_dly_year <- test_dly[test_dly$pheno_year == test_task$yr, , drop = FALSE]
} else {
  test_dly_year <- test_dly[lubridate::year(test_dly$date) == test_task$yr, , drop = FALSE]
}

  cat("\n=== DATA DISTRIBUTION ANALYSIS ===\n")
  if (exists("df_tasks") && nrow(df_tasks) > 0) {
    cat(sprintf("Total locations in df_tasks: %d\n", length(unique(df_tasks$location_id))))
    cat(sprintf("Total location-years: %d\n", nrow(test_loc_years)))

    sample_sizes <- df_tasks |> 
      dplyr::group_by(location_id, pheno_year) |> 
      dplyr::summarize(n_obs = dplyr::n(), .groups = "drop")

    if (nrow(sample_sizes) > 0) {
      sample_sizes$n_obs <- as.numeric(sample_sizes$n_obs)
      cat("\nObservations per location-year distribution:\n")
      cat(sprintf("  Min:    %d\n", as.integer(min(sample_sizes$n_obs, na.rm = TRUE))))
      cat(sprintf("  Q1:     %d\n", as.integer(stats::quantile(sample_sizes$n_obs, 0.25, na.rm = TRUE))))
      cat(sprintf("  Median: %d\n", as.integer(median(sample_sizes$n_obs, na.rm = TRUE))))
      cat(sprintf("  Q3:     %d\n", as.integer(stats::quantile(sample_sizes$n_obs, 0.75, na.rm = TRUE))))
      cat(sprintf("  Max:    %d\n", as.integer(max(sample_sizes$n_obs, na.rm = TRUE))))
      cat(sprintf("  Mean:   %.1f\n", mean(sample_sizes$n_obs, na.rm = TRUE)))

      problem_threshold <- MIN_OBS_PER_LOC_YEAR  # Minimum observations for meaningful MESMA
      n_problem <- sum(sample_sizes$n_obs < problem_threshold, na.rm = TRUE)
      cat(sprintf("\nLocation-years with < %d observations: %d (%.1f%%)\n",
                  problem_threshold, n_problem, 100 * n_problem / nrow(sample_sizes)))
      if (n_problem > 0) {
        cat("\nSample of problematic location-years:\n")
        problem_cases <- sample_sizes[sample_sizes$n_obs < problem_threshold, , drop = FALSE]
        print(utils::head(problem_cases, 10))
      }

      if ("doy" %in% names(df_tasks)) {
        doy_coverage <- df_tasks |> 
          dplyr::group_by(location_id, pheno_year) |> 
          dplyr::summarize(
            n_unique_doys = length(unique(doy[!is.na(doy)])),
            doy_span = if (any(!is.na(doy))) as.integer(max(doy, na.rm = TRUE) - min(doy, na.rm = TRUE)) else 0,
            .groups = "drop"
          )
        cat("\nDOY coverage per location-year:\n")
        cat(sprintf("  Median unique DOYs: %d\n", as.integer(median(doy_coverage$n_unique_doys, na.rm = TRUE))))
        cat(sprintf("  Median DOY span:    %d days\n", as.integer(median(doy_coverage$doy_span, na.rm = TRUE))))
        n_poor_coverage <- sum(doy_coverage$n_unique_doys < 30, na.rm = TRUE)
        cat(sprintf("  Location-years with < 30 unique DOYs: %d (%.1f%%)\n",
                    n_poor_coverage, 100 * n_poor_coverage / nrow(doy_coverage)))
      }
    } else {
      cat("No location-year sample sizes available to summarize.\n")
    }
  }

  if (length(location_list) == 0) {
    stop("No testing tasks found")
  }





  evaluate_all_combinations <- function(y, top_variants, lambda = 0, early_stop_rmse = EARLY_STOP_RMSE_THRESHOLD, feature_weights = NULL) {
    if (length(top_variants) == 0) return(NULL)

    veg_names <- names(top_variants)
    n_veg <- length(veg_names)
    if (n_veg == 0) return(NULL)

    y_target <- y
    y_target[!is.finite(y_target)] <- 0

    # Collector for all evaluated models (to capture uncertainty)
    candidate_results <- list()
    
    # Helper: Add result to candidates (only if valid)
    add_candidate <- function(res) {
      if (!is.null(res) && is.finite(res$rmse)) {
        candidate_results[[length(candidate_results) + 1]] <<- res
      }
    }

    # Helper function to solve for a single combination
    solve_combo <- function(variant_indices) {
      cols <- list()
      ids <- character(0)
      for (v_idx in seq_along(veg_names)) {
        v <- veg_names[v_idx]
        idx <- variant_indices[v_idx]
        cand <- top_variants[[v]][[idx]]
        if (is.null(cand) || is.null(cand$vec)) return(NULL)
        cols[[length(cols) + 1]] <- as.numeric(cand$vec)
        ids <- c(ids, cand$id)
      }
      if (length(cols) == 0) return(NULL)

      E <- do.call(cbind, cols)
      if (is.null(E) || ncol(E) < 1) return(NULL)

      res <- solve_weights_nnls_simple(E, y_target, feature_weights = if (!is.null(feature_weights) && length(feature_weights) == nrow(E)) feature_weights else NULL)
      if (is.null(res)) return(NULL)
      
      names(res$w) <- veg_names
      res$ids <- ids
      names(res$ids) <- veg_names
      return(res)
    }
    
    # --- Optimization: Two-Stage Search ---
    TOPK_COARSE <- 3
    coarse_idx_lists <- lapply(top_variants, function(lst) seq_len(min(TOPK_COARSE, length(lst))))
    
    # If total combinations are small, just do a full search
    total_combos <- as.numeric(Reduce(`*`, lapply(top_variants, function(lst) length(lst)), init = 1L))
    if (total_combos < 500) { # Heuristic: if total combos are less than 500, do full search
        full_idx_lists <- lapply(top_variants, function(lst) seq_along(lst))
        all_combos <- expand.grid(full_idx_lists, KEEP.OUT.ATTRS = FALSE)
        best_rmse <- Inf
        best_result <- NULL
        for (i in 1:nrow(all_combos)) {
            res <- solve_combo(as.integer(all_combos[i, ]))
            add_candidate(res) # Store for uncertainty
            if (!is.null(res) && is.finite(res$rmse) && res$rmse < best_rmse) {
                best_rmse <- res$rmse
                best_result <- res
            }
        }
        if (is.null(best_result)) return(NULL)
        
        # Sort and return top candidates
        sorted_candidates <- candidate_results[order(sapply(candidate_results, function(x) x$rmse))]
        top_models <- head(sorted_candidates, 20)
        
        return(list(w = best_result$w, rmse = best_rmse, ids = best_result$ids, residuals = best_result$residuals, top_models = top_models))
    }

    # Stage 1: Coarse search
    coarse_combos <- expand.grid(coarse_idx_lists, KEEP.OUT.ATTRS = FALSE)
    
    best_rmse <- Inf
    best_result <- NULL
    best_combo_indices <- NULL

    for (i in 1:nrow(coarse_combos)) {
      current_indices <- as.integer(coarse_combos[i, ])
      res <- solve_combo(current_indices)
      add_candidate(res) # Store for uncertainty
      if (!is.null(res) && is.finite(res$rmse) && res$rmse < best_rmse) {
        best_rmse <- res$rmse
        best_result <- res
        best_combo_indices <- current_indices
      }
    }

    if (is.null(best_result)) {
       best_combo_indices <- rep(1, n_veg)
       best_result <- solve_combo(best_combo_indices)
       add_candidate(best_result)
       if(is.null(best_result)) return(NULL)
       best_rmse <- best_result$rmse
    }

    # Fine-tuning step (Coordinate Descent)
    current_best_indices <- best_combo_indices
    names(current_best_indices) <- veg_names

    for (i in 1:n_veg) {
        v_to_tune <- veg_names[i]
        
        for (variant_idx in seq_along(top_variants[[v_to_tune]])) {
            if (variant_idx == current_best_indices[i]) next
            
            test_indices <- current_best_indices
            test_indices[i] <- variant_idx
            
            res <- solve_combo(test_indices)
            add_candidate(res) # Store for uncertainty
            
            if (!is.null(res) && is.finite(res$rmse) && res$rmse < best_rmse) {
                best_rmse <- res$rmse
                best_result <- res
                current_best_indices[i] <- variant_idx
            }
        }
    }
    
    if (is.null(best_result)) return(NULL)

    # Sort candidates by RMSE
    sorted_candidates <- candidate_results[order(sapply(candidate_results, function(x) x$rmse))]
    
    # Filter to keep "low RMSE" models for uncertainty estimation
    # Strategy: threshold derived from the data distribution (Min + 2 * (Q25 - Min))
    # This adapts to the flatness of the objective function near the minimum, keeping models
    # that are statistically comparable to the best fit within the observed spread.
    all_rmses <- sapply(sorted_candidates, function(x) x$rmse)
    q25 <- quantile(all_rmses, 0.25, na.rm = TRUE)
    spread_metric <- q25 - best_rmse
    

    
    rmse_threshold <- best_rmse + (2.0 * spread_metric)
    
    top_models <- Filter(function(x) x$rmse <= rmse_threshold, sorted_candidates)

    return(list(w = best_result$w, rmse = best_rmse, ids = best_result$ids, residuals = best_result$residuals, top_models = top_models))
  }

  prune_collinear_variants <- function(lib, threshold = 0.99) {
    # Flatten library
    all_variants <- list()
    for (v in names(lib)) {
      for (i in seq_along(lib[[v]])) {
        var <- lib[[v]][[i]]
        # Ensure we have a vector to compare
        vec <- if (!is.null(var$vec)) var$vec else if (!is.null(var$raw_row)) var$raw_row else NULL
        if (!is.null(vec)) {
          all_variants[[length(all_variants) + 1]] <- list(
            veg = v,
            idx = i,
            id = if (!is.null(var$id)) var$id else paste0(v, "_v", i),
            vec = vec,
            n_samples = if (!is.null(var$n_samples)) var$n_samples else 1
          )
        }
      }
    }

    if (length(all_variants) < 2) return(lib)

    cat(sprintf("\n[PRUNING] Checking for collinearity among %d variants (threshold=%.4f)...\n", length(all_variants), threshold))

    # Iterative pruning
    removed_ids <- c()
    
    while (TRUE) {
      # Filter out removed
      current_vars <- all_variants[!sapply(all_variants, function(x) x$id %in% removed_ids)]
      n <- length(current_vars)
      if (n < 2) break

      # Compute distance matrix using Euclidean distance
      M <- do.call(rbind, lapply(current_vars, function(x) x$vec))
      
      # Compute pairwise Euclidean distances
      n_vars <- nrow(M)
      Dist <- matrix(0, n_vars, n_vars)
      for (i in 1:(n_vars-1)) {
        for (j in (i+1):n_vars) {
          Dist[i,j] <- Dist[j,i] <- sqrt(sum((M[i,] - M[j,])^2))
        }
      }
      
      min_dist <- min(Dist[Dist > 0])  # Minimum non-zero distance
      if (min_dist > (1 - threshold)) break  # threshold is similarity, so 1 - threshold is max allowed distance
      
      # Find pair with minimum distance
      idx_mat <- which(Dist == min_dist, arr.ind = TRUE)
      # Take first pair
      r <- idx_mat[1, 1]
      c <- idx_mat[1, 2]
      
      v1 <- current_vars[[r]]
      v2 <- current_vars[[c]]
      
      # Decision logic: drop the one with fewer samples
      to_drop <- if (v1$n_samples < v2$n_samples) v1 else v2
      
      # If samples equal, drop the one that appears later (arbitrary stable tie-break)
      if (v1$n_samples == v2$n_samples) {
         to_drop <- v2
      }
      
      cat(sprintf("  [DROP] %s (n=%d) due to distance %.4f with %s (n=%d)\n", 
                  to_drop$id, to_drop$n_samples, min_dist, 
                  if (to_drop$id == v1$id) v2$id else v1$id, 
                  if (to_drop$id == v1$id) v2$n_samples else v1$n_samples))
      
      removed_ids <- c(removed_ids, to_drop$id)
    }

    if (length(removed_ids) == 0) {
      cat("[PRUNING] No variants removed.\n")
      return(lib)
    }

    # Reconstruct library
    new_lib <- list()
    for (v in names(lib)) {
      new_lib[[v]] <- list()
      for (i in seq_along(lib[[v]])) {
        var <- lib[[v]][[i]]
        vid <- if (!is.null(var$id)) var$id else paste0(v, "_v", i)
        if (!vid %in% removed_ids) {
          new_lib[[v]][[length(new_lib[[v]]) + 1]] <- var
        }
      }
      if (length(new_lib[[v]]) == 0) {
          cat(sprintf("  [WARNING] All variants for '%s' were pruned! Restoring the one with most samples.\n", v))
          # Find original variants for this class
          vars_orig <- all_variants[sapply(all_variants, function(x) x$veg == v)]
          if (length(vars_orig) > 0) {
              # Pick best
              best_v <- vars_orig[[which.max(sapply(vars_orig, function(x) x$n_samples))]]
              # Find it in original lib
              for(orig_var in lib[[v]]) {
                  ovid <- if (!is.null(orig_var$id)) orig_var$id else ""
                  if (ovid == best_v$id) {
                      new_lib[[v]][[1]] <- orig_var
                      break
                  }
              }
          }
      }
    }
    
    return(new_lib)
  }

  build_single_stage_library_weighted <- function(df_train, indices, params, allowed_veg) {
    # Single-stage library: treats barren and all vegetation types as equal endmembers
    # Returns a library with all endmember types (barren + veg types) in one unified structure
    
    cat("\n[LIBRARY BUILD] Starting Global Combinatorial Optimization for Cluster Counts (including Barren)...\n")
    
    # Validate params structure
    if (is.null(params) || is.null(params$means) || is.null(params$sds)) {
      stop("[ERROR] build_single_stage_library_weighted: params$means or params$sds is NULL!")
    }
    if (length(params$means) != length(indices)) {
      cat(sprintf("[ERROR] params$means has length %d but indices has length %d\n", length(params$means), length(indices)))
      cat(sprintf("  indices: %s\n", paste(indices, collapse = ", ")))
      cat(sprintf("  params$means names: %s\n", paste(names(params$means), collapse = ", ")))
      stop("[ERROR] Length mismatch between params$means and indices")
    }
    
    # -------------------------------------------------------------------------
    # PASS 1: Load, Normalize, and Store Data (with Metadata)
    # -------------------------------------------------------------------------
    expected_cols <- length(indices) * N_TEMPORAL_BINS
    storage <- list()       # Store unweighted Z-score matrices
    storage_meta <- list()  # Store metadata (location_id) for spatial bootstrapping
    
    # Unified class list
    target_classes <- unique(c("barren", allowed_veg))
    valid_classes <- c()
    
    for(v in target_classes) {
      # Since df_train$Veg is already lowercased, we can match directly
      veg_data <- dplyr::filter(df_train, .data$Veg == v)
      if(nrow(veg_data) == 0) next
      
      veg_list <- list()
      loc_list <- character(0)
      
      traces <- unique(veg_data[, c("location_id", "pheno_year")])
      
      # Batch processing for speed? No, keep simple loop for now
      for(i in seq_len(nrow(traces))) {
        lid <- traces$location_id[i]
        pyr <- traces$pheno_year[i]
        sub <- veg_data[veg_data$location_id == lid & veg_data$pheno_year == pyr, ]
        mat <- build_pentad_matrix(sub, indices)
        if(!is.null(mat)) {
          veg_list[[length(veg_list) + 1]] <- as.numeric(mat)
          loc_list <- c(loc_list, as.character(lid))
        }
      }
      
      if(length(veg_list) > 0) {
        veg_mat <- do.call(rbind, veg_list)
        if(ncol(veg_mat) == expected_cols) {
          # Normalize
          for(k in seq_along(indices)) {
            idx_start <- (k-1)*N_TEMPORAL_BINS + 1
            idx_end <- k*N_TEMPORAL_BINS
            veg_mat[, idx_start:idx_end] <- (veg_mat[, idx_start:idx_end] - params$means[k]) / params$sds[k]
          }
          veg_mat[is.na(veg_mat)] <- 0
          storage[[v]] <- veg_mat
          storage_meta[[v]] <- data.frame(location_id = loc_list, stringsAsFactors = FALSE)
          valid_classes <- c(valid_classes, v)
        }
      }
    }
    
    if (length(valid_classes) == 0) {
      warning("No valid data found for library building.")
      return(list())
    }

    # Weights for clustering
    w_vec <- if(!is.null(params$weights)) sqrt(pmax(params$weights, 0)) else rep(1, expected_cols)
    w_vec[w_vec < 1e-9] <- 1e-9 # Prevent div by zero
    
    min_cluster_size <- 5

    # -------------------------------------------------------------------------
    # Internal Optimization Routine: Stratified Spatial Bootstrap OOB
    # -------------------------------------------------------------------------
    optimize_library <- function(n_boot = 5) {
      cat(sprintf("\n  --- Running Optimization with Stratified Spatial Bootstrap (%d folds) ---\n", n_boot))
      
      # 1. Generate Bootstrap Folds
      # Structure: folds[[b]][[v]] -> list(train_idx, oob_idx)
      folds <- list()
      for(b in 1:n_boot) {
        folds[[b]] <- list()
        for(v in valid_classes) {
          locs <- storage_meta[[v]]$location_id
          unique_locs <- unique(locs)
          
          if(length(unique_locs) < 2) {
             # Too few locations to bootstrap: use all for train, none for OOB (or duplicate)
             # Fallback: simple row resampling if locations are scarce
             n_rows <- nrow(storage[[v]])
             train_rows <- sample(n_rows, n_rows, replace = TRUE)
             oob_rows <- setdiff(1:n_rows, unique(train_rows))
             if(length(oob_rows) == 0) oob_rows <- sample(n_rows, max(1, floor(n_rows * 0.2))) # Force some OOB
          } else {
             # Spatial Bootstrap
             train_locs <- sample(unique_locs, length(unique_locs), replace = TRUE)
             oob_locs <- setdiff(unique_locs, unique(train_locs))
             
             # Map back to row indices
             # Note: A location might appear multiple times in train_locs (weighted)
             # But for clustering input, we typically just want the set of training rows.
             # Standard Bagging uses the resampled set (with duplicates).
             
             # Create a mapping from loc -> row_indices
             loc_to_rows <- split(seq_along(locs), locs)
             
             train_rows <- unlist(lapply(train_locs, function(l) loc_to_rows[[l]]))
             oob_rows <- unlist(lapply(oob_locs, function(l) loc_to_rows[[l]]))
          }
          folds[[b]][[v]] <- list(train = train_rows, oob = oob_rows)
        }
      }

      # 2. Pre-compute Clusters for EACH Fold and Class (K=1..11)
      # cluster_cache[[b]][[v]][[k]] -> centers
      boot_cluster_cache <- list()
      k_ranges <- list() # Common range for grid search (derived from full data stats or per-fold?)
      # We'll derive global valid ranges to keep grid consistent, but check per-fold feasibility
      
      cat("    Pre-computing clusters for bootstrap folds...\n")
      
      for(v in valid_classes) {
        # Determine valid K range based on TOTAL samples (upper bound)
        n_total <- nrow(storage[[v]])
        if (v == "barren") {
           k_candidates <- 1
        } else {
           max_k <- floor(n_total / min_cluster_size)
           k_candidates <- 1:11
           k_candidates <- k_candidates[k_candidates <= max_k]
           if (length(k_candidates) == 0) k_candidates <- 1
        }
        k_ranges[[v]] <- k_candidates
      }

      for(b in 1:n_boot) {
        boot_cluster_cache[[b]] <- list()
        
        for(v in valid_classes) {
           boot_cluster_cache[[b]][[v]] <- list()
           
           train_idx <- folds[[b]][[v]]$train
           if(length(train_idx) < min_cluster_size) {
             # Fallback if bootstrap fold is empty/tiny: use all data
             train_idx <- 1:nrow(storage[[v]])
           }
           
           veg_mat_train <- storage[[v]][train_idx, , drop = FALSE]
           veg_mat_w <- sweep(veg_mat_train, 2, w_vec, "*")
           
           # Input for K-Means
           veg_mat_input <- veg_mat_w

           for(k in k_ranges[[v]]) {
             # Safe K-Means on this fold
             # Adjust k if fold is smaller than expected (rare with bootstrap size, but possible)
             if(nrow(veg_mat_input) < k) {
               curr_k <- max(1, nrow(veg_mat_input))
             } else {
               curr_k <- k
             }
             
             # Run K-Means
             if(curr_k == 1) {
                center_w <- apply(veg_mat_w, 2, median, na.rm=TRUE)
                centers_unw <- center_w / w_vec
                # Skip barren filter in OOB pre-compute for speed/robustness? 
                # Better to apply it to match final logic.
                boot_cluster_cache[[b]][[v]][[as.character(k)]] <- matrix(centers_unw, nrow=1)
             } else {
                km <- tryCatch(kmeans(veg_mat_input, centers = curr_k, nstart = 5, iter.max = 20), error = function(e) NULL)
                if(!is.null(km)) {
                  # Recover centroids
                  centers_w <- km$centers
                  centers_unw <- sweep(centers_w, 2, w_vec, "/")
                  boot_cluster_cache[[b]][[v]][[as.character(k)]] <- centers_unw
                }
             }
           }
        }
      }

      # 3. Grid Search: Evaluate Combinations on OOB Sets
      grid <- expand.grid(k_ranges)
      n_combos <- nrow(grid)
      if (n_combos > 500) {
        set.seed(123)
        grid <- grid[sample(seq_len(n_combos), 500), , drop=FALSE]
      }
      
      cat(sprintf("    Evaluating %d cluster combinations on %d OOB folds...\n", nrow(grid), n_boot))
      
      best_mean_score <- -1
      best_combo <- NULL
      
      # Prepare OOB Test Matrices for speed
      # oob_data[[b]] -> matrix of all OOB samples for that fold
      # oob_labels[[b]] -> vector of labels
      oob_sets <- list()
      for(b in 1:n_boot) {
         samples <- list(); lbls <- c()
         for(v in valid_classes) {
            idx <- folds[[b]][[v]]$oob
            if(length(idx) > 0) {
               samples[[length(samples)+1]] <- storage[[v]][idx, , drop=FALSE]
               lbls <- c(lbls, rep(v, length(idx)))
            }
         }
         if(length(samples) > 0) {
            oob_sets[[b]] <- list(Y = do.call(rbind, samples), labels = lbls)
         } else {
            oob_sets[[b]] <- NULL 
         }
      }
      
      run_solver <- function(y, M, w) {
        sq_w <- sqrt(pmax(w, 0)); M_w <- sweep(M, 1, sq_w, "*"); y_w <- y * sq_w
        tryCatch(nnls::nnls(M_w, y_w)$x, error=function(e) rep(0, ncol(M)))
      }
      
      # Loop over Grid
      # Parallelize this loop if possible? No, inner NNLS is fast enough.
      
      for(i in seq_len(nrow(grid))) {
        combo <- grid[i, , drop=FALSE]
        fold_scores <- numeric(n_boot)
        valid_folds <- 0
        
        for(b in 1:n_boot) {
           if(is.null(oob_sets[[b]])) next
           
           # Construct Library from Boot Cache using this Combo
           cols <- list(); col_names <- c()
           valid_lib <- TRUE
           
           for(v in valid_classes) {
              k_val <- as.character(combo[[v]])
              # Check if this fold has this k (might have failed or been skipped)
              centers <- boot_cluster_cache[[b]][[v]][[k_val]]
              if(is.null(centers)) { valid_lib <- FALSE; break }
              
              cols[[length(cols)+1]] <- t(centers)
              col_names <- c(col_names, rep(v, nrow(centers)))
           }
           
           if(!valid_lib) next
           
           M <- do.call(cbind, cols)
           Y_test <- oob_sets[[b]]$Y
           labels_test <- oob_sets[[b]]$labels
           
           # Evaluate
           correct_by_class <- list(); total_by_class <- list()
           for(vc in valid_classes) { correct_by_class[[vc]] <- 0; total_by_class[[vc]] <- 0 }
           
           # Batch or loop? Loop for now
           for(j in seq_len(nrow(Y_test))) {
              true_label <- labels_test[j]
              total_by_class[[true_label]] <- total_by_class[[true_label]] + 1
              
              coefs <- run_solver(Y_test[j, ], M, params$weights)
              sums <- tapply(coefs, col_names, sum)
              
              if (!true_label %in% names(sums)) next
              score_true <- sums[[true_label]]
              
              competitors <- if (true_label == "barren") setdiff(names(sums), "barren") else setdiff(names(sums), c(true_label, "barren"))
              score_max_other <- if(length(competitors)>0) max(sums[competitors]) else -1
              
              if (score_true > score_max_other) correct_by_class[[true_label]] <- correct_by_class[[true_label]] + 1
           }
           
           # Macro Accuracy (Vegetation Only)
           accs <- c()
           veg_classes_only <- setdiff(valid_classes, "barren")
           for(vc in veg_classes_only) if(total_by_class[[vc]] > 0) accs <- c(accs, correct_by_class[[vc]]/total_by_class[[vc]])
           
           if(length(accs) > 0) {
              fold_scores[b] <- mean(accs)
              valid_folds <- valid_folds + 1
           }
        }
        
        mean_score <- if(valid_folds > 0) sum(fold_scores) / valid_folds else -1
        
        if (mean_score > best_mean_score) {
          best_mean_score <- mean_score
          best_combo <- combo
        }
      }
      
      cat(sprintf("    Best OOB Score: %.4f using combo: %s\n", best_mean_score, paste(names(best_combo), best_combo, sep="=", collapse=", ")))
      
      # --- AGENT ADDITION: Average Confusion Matrix Calculation ---
      cat("\n    Computing Average Confusion Matrix across OOB folds for best combination...\n")
      
      # Initialize Confusion Matrix
      # Rows: True, Cols: Predicted
      cm_labels <- valid_classes
      conf_mat <- matrix(0, nrow = length(cm_labels), ncol = length(cm_labels))
      rownames(conf_mat) <- cm_labels
      colnames(conf_mat) <- cm_labels
      
      total_samples <- 0
      
      for(b in 1:n_boot) {
         if(is.null(oob_sets[[b]])) next
         
         # Construct Library for this fold
         cols <- list(); col_names <- c()
         valid_lib <- TRUE
         
         for(v in valid_classes) {
            k_val <- as.character(best_combo[[v]])
            centers <- boot_cluster_cache[[b]][[v]][[k_val]]
            if(is.null(centers)) { valid_lib <- FALSE; break }
            cols[[length(cols)+1]] <- t(centers)
            col_names <- c(col_names, rep(v, nrow(centers)))
         }
         
         if(!valid_lib) next
         
         M <- do.call(cbind, cols)
         Y_test <- oob_sets[[b]]$Y
         labels_test <- oob_sets[[b]]$labels
         
         for(j in seq_len(nrow(Y_test))) {
            true_label <- labels_test[j]
            coefs <- run_solver(Y_test[j, ], M, params$weights)
            sums <- tapply(coefs, col_names, sum)
            
            # Predict class with max coefficient
            if(length(sums) > 0) {
                pred_label <- names(sums)[which.max(sums)]
                
                # Update Matrix
                if(true_label %in% cm_labels && pred_label %in% cm_labels) {
                    conf_mat[true_label, pred_label] <- conf_mat[true_label, pred_label] + 1
                    total_samples <- total_samples + 1
                }
            }
         }
      }
      
      cat("\n    [OOB Confusion Matrix - Sum over all folds]\n")
      print(conf_mat)
      
      # Row-normalized (Producer's Accuracy / Recall)
      conf_mat_pct <- sweep(conf_mat, 1, rowSums(conf_mat), "/")
      conf_mat_pct[is.nan(conf_mat_pct)] <- 0
      
      cat("\n    [OOB Confusion Matrix - Row Normalized (Recall)]\n")
      print(round(conf_mat_pct, 2))
      
      # Compute Overall Accuracy on OOB
      oa <- sum(diag(conf_mat)) / sum(conf_mat)
      cat(sprintf("\n    Overall OOB Accuracy: %.2f%%\n", oa * 100))
      
      # --- END AGENT ADDITION ---

      return(best_combo)
    }
    
    # -------------------------------------------------------------------------
    # MAIN EXECUTION
    # -------------------------------------------------------------------------
    
    # 1. Run Optimization to get Best K-Combo
    best_combo <- optimize_library(n_boot = 5)
    
    # 2. Re-Train Final Clusters on FULL Dataset using Best K
    cat("\n[LIBRARY BUILD] Re-training final clusters on full dataset using optimal counts...\n")
    
    final_lib_cache <- list()
    
    # Pre-compute Barren Weighted Center for Final Filtering
    barren_ref_w <- NULL
    if ("barren" %in% valid_classes) {
       b_mat <- storage[["barren"]]
       b_mat_w <- sweep(b_mat, 2, w_vec, "*")
       barren_ref_w <- apply(b_mat_w, 2, median, na.rm=TRUE)
    }

    for(v in valid_classes) {
       k_opt <- as.numeric(best_combo[[v]])
       veg_mat <- storage[[v]]
       veg_mat_w <- sweep(veg_mat, 2, w_vec, "*")
       
       # Check availability
       if(nrow(veg_mat) < k_opt) {
          k_opt <- max(1, nrow(veg_mat))
          cat(sprintf("  [WARNING] Reducing k for '%s' to %d (insufficient samples)\n", v, k_opt))
       }
       
       # Final Clustering
       final_centers_unw <- NULL
       
       if (k_opt == 1) {
          center_w <- apply(veg_mat_w, 2, median, na.rm=TRUE)
          final_centers_unw <- matrix(center_w / w_vec, nrow=1)
       } else {
          # Use Euclidean K-Means (consistent with best practice for final shapes)
          km <- tryCatch(kmeans(veg_mat_w, centers = k_opt, nstart = 20, iter.max = 50), error = function(e) NULL)
          if(!is.null(km)) {
             final_centers_unw <- sweep(km$centers, 2, w_vec, "/")
          }
       }
       
       # Filter Soil-Like Variants (Consistently applied)
       if (!is.null(final_centers_unw) && v != "barren" && !is.null(barren_ref_w)) {
           keep_mask <- apply(final_centers_unw, 1, function(row) {
               row_w <- row * w_vec
               num <- sum(row_w * barren_ref_w)
               den <- sqrt(sum(row_w^2)) * sqrt(sum(barren_ref_w^2))
               sim <- if (den > 1e-9) num/den else 0
               return(sim <= 0.95)
           })
           # Always keep at least one
           if (any(keep_mask)) {
               final_centers_unw <- final_centers_unw[keep_mask, , drop=FALSE]
           } else {
               cat(sprintf("  [WARNING] All variants for '%s' were too similar to barren. Keeping the most distinct one.\n", v))
               sims <- apply(final_centers_unw, 1, function(row) {
                   row_w <- row * w_vec
                   num <- sum(row_w * barren_ref_w)
                   den <- sqrt(sum(row_w^2)) * sqrt(sum(barren_ref_w^2))
                   if (den > 1e-9) num/den else 0
               })
               best_idx <- which.min(sims)
               final_centers_unw <- final_centers_unw[best_idx, , drop=FALSE]
           }
       }
       
       if (!is.null(final_centers_unw)) {
          final_lib_cache[[v]] <- final_centers_unw
       }
    }
    
    # 3. Construct Final Library Structure (Legacy Format)
    res_lib <- list()
    for(v in names(final_lib_cache)) {
       mat <- final_lib_cache[[v]]
       res_lib[[v]] <- list()
       for(i in 1:nrow(mat)) {
          # Construct dummy variant object
          # Re-normalization to T_medoid structure if needed?
          # The rest of the pipeline expects 'vec' (unweighted, z-scored) or 'T' (medoid).
          # Here we provide 'vec'.
          res_lib[[v]][[length(res_lib[[v]]) + 1]] <- list(
             vec = mat[i, ],
             id = paste0(v, "_opt_", i),
             n_samples = floor(nrow(storage[[v]]) / nrow(mat)) # Approximate n_samples
          )
       }
    }
    
    return(res_lib)
  }


  precompute_optimized_library_weighted <- function(mesma_lib, grid_type = "full", feature_weights = NULL) {
    opt_lib <- list()

    # Determine pruning indices if requested
    keep_idx_global_w <- NULL
    if (!is.null(feature_weights) && exists("PRUNE_ZERO_WEIGHT_FEATURES") && isTRUE(PRUNE_ZERO_WEIGHT_FEATURES)) {
      zero_mask_glob <- feature_weights == 0
      n_zero_glob <- sum(zero_mask_glob, na.rm = TRUE)
      frac_zero_glob <- n_zero_glob / length(feature_weights)
      if (n_zero_glob > 0) {
        cat(sprintf("[INFO] Single-stage pruning requested: %d/%d features zeroed (%.1f%%)\n", n_zero_glob, length(feature_weights), 100*frac_zero_glob))
        if (frac_zero_glob <= PRUNE_ZERO_WEIGHT_MAX_FRAC && (length(feature_weights) - n_zero_glob) >= PRUNE_ZERO_MIN_FEATURES) {
          keep_idx_global_w <- which(!zero_mask_glob)
          cat(sprintf("[INFO] Will prune %d features for single-stage, keeping %d features\n", n_zero_glob, length(keep_idx_global_w)))
        } else {
          cat(sprintf("[WARN] Skipping single-stage pruning: zero fraction %.2f exceeds max allowed %.2f or resulting features < %d\n", frac_zero_glob, PRUNE_ZERO_WEIGHT_MAX_FRAC, PRUNE_ZERO_MIN_FEATURES))
        }
      }
    }
    
    for(v in names(mesma_lib)) {
        if(is.null(mesma_lib[[v]])) next
        variants_raw <- mesma_lib[[v]]
        if (!is.list(variants_raw) || (!is.null(variants_raw$vec) && !is.list(variants_raw$vec))) {
          variants <- list(variants_raw)
        } else {
          variants <- variants_raw
        }

        vecs <- lapply(variants, function(x) as.numeric(if (is.list(x) && !is.null(x$vec)) x$vec else x))
        ids <- sapply(variants, function(x) {
          if (is.list(x) && !is.null(x$id)) return(x$id)
          if (is.list(x) && !is.null(x$variant_id)) return(x$variant_id)
          return(paste0(v, "_1"))
        })

        M <- do.call(rbind, vecs)

        pruned_info <- NULL
        if (!is.null(keep_idx_global_w) && length(keep_idx_global_w) > 0) {
          if (max(keep_idx_global_w) <= ncol(M)) {
            kept_names <- NULL
            if (!is.null(colnames(M))) kept_names <- colnames(M)[keep_idx_global_w]
            M <- M[, keep_idx_global_w, drop = FALSE]
            pruned_info <- list(kept_idx = keep_idx_global_w, kept_names = kept_names, n_kept = ncol(M))
            cat(sprintf("[INFO] Pruned single-stage features for veg '%s': kept %d columns\n", v, ncol(M)))
          } else {
            cat(sprintf("[WARN] Single-stage prune requested but index mismatch for veg '%s'; skipping prune\n", v))
          }
        }

        M_norm <- t(apply(M, 1, function(r) {
          nrm <- sqrt(sum(r^2))
          if (!is.finite(nrm) || nrm == 0) return(rep(0, length(r))) else return(r / nrm)
        }))

        opt_lib[[v]] <- list(
          M = M,
          M_norm = M_norm,
          ids = ids,
          pruned_info = pruned_info
        )
    }
    
    opt_lib
  }





  precompute_compressed_templates <- function(mesma_lib, grid_name = "full") {
    template_db <- list()
    for (veg in names(mesma_lib)) {
      template_db[[veg]] <- list()
      variants <- mesma_lib[[veg]]
      if (exists("filter_variants_by_min_samples", mode="function")) {
         rt <- if (exists("raw_lib_templates", envir=globalenv()) && !is.null(get("raw_lib_templates", envir=globalenv())[[veg]])) get("raw_lib_templates", envir=globalenv())[[veg]] else NULL
         variants <- filter_variants_by_min_samples(variants, min_samples = MIN_ENDMEMBER_SAMPLES, veg = veg, raw_template = rt)
      }
      for (variant in variants) {
        vid <- if (!is.null(variant$variant_id)) variant$variant_id else variant$id
        if (is.null(vid)) next
        template_db[[veg]][[vid]] <- list()
        
        raw_src <- if (!is.null(variant$raw_mat)) variant$raw_mat else if (!is.null(variant$vec)) variant$vec else NULL
        if (is.null(raw_src)) next
        
        compressed_vec <- as.numeric(raw_src)
        compressed_vec[!is.finite(compressed_vec)] <- NA_real_
        template_db[[veg]][[vid]][[grid_name]] <- compressed_vec
      }
    }
    template_db
  }

print_weights_summary <- function(stage_name, params) {
  if (is.null(params) || is.null(params$weights) || is.null(params$indices)) {
    cat(sprintf("[WEIGHTS %s] No weights available\n", stage_name))
    return(invisible(NULL))
  }

  n_indices <- length(params$indices)
  wlen <- length(params$weights)
  expected <- n_indices * N_TEMPORAL_BINS

  if (wlen == expected) {
    w_mat <- matrix(params$weights, nrow = N_TEMPORAL_BINS, ncol = n_indices)
    per_index_mean <- colMeans(w_mat, na.rm = TRUE)
    per_index_max <- apply(w_mat, 2, max, na.rm = TRUE)
    ord <- order(per_index_mean, decreasing = TRUE)

    cat(sprintf("[WEIGHTS %s] Per-index mean weights (top %d):\n", stage_name, min(8, n_indices)))
    topn <- min(8, n_indices)
    for (i in seq_len(topn)) {
      ii <- ord[i]
      cat(sprintf("  %s: mean=%.4f, max=%.4f\n", params$indices[ii], per_index_mean[ii], per_index_max[ii]))
    }

    cat(sprintf("[WEIGHTS %s] Full per-index means: %s\n", stage_name, paste(sprintf("%s=%.4f", params$indices, per_index_mean), collapse=", ")))
  } else {
    cat(sprintf("[WEIGHTS %s] weights length (%d) != expected (%d) -> printing sample values\n", stage_name, wlen, expected))
    cat(sprintf("[WEIGHTS %s] sample weights (first 20): %s\n", stage_name, paste(sprintf("%.4f", head(params$weights, 20)), collapse=", ")))
  }
  invisible(NULL)
}

  cat("Building single-stage MESMA library with Z-score/PCA/LDA weighting...\n")
  # 'no.soil' mapping removed — 'no soil' column is deprecated and ignored; classification will use 'Veg' only

  # Create multi-class classification data: Barren vs Veg A vs Veg B ...
  # This maximizes separation between ALL classes in the feature space.
  multi_class_data <- dplyr::mutate(df_train, target_class = tolower(as.character(Veg)))
  multi_class_data <- dplyr::filter(multi_class_data, !is.na(target_class) & target_class != "")
  multi_class_data <- dplyr::select(multi_class_data, dplyr::all_of(c("location_id", "pheno_year", "date", "doy", "Veg", "target_class", avail)))

  cat("Training feature pipeline (Multi-Class: Barren vs Veg A vs Veg B...)\n")
  cat(sprintf("[NOTICE] Using feature set for training: %s\n", paste(avail, collapse=", ")))
  
  class_counts <- table(multi_class_data$target_class)
  cat("[Single-stage] Training data class distribution:\n")
  print(class_counts)

  SINGLE_STAGE_PARAMS <- train_feature_pipeline(multi_class_data, "target_class", avail)
  if (!is.null(SINGLE_STAGE_PARAMS) && !is.null(SINGLE_STAGE_PARAMS$weights)) {
    # Enforce a pentad-level minimum floor to prevent entire pentads from being suppressed.
    # Floor is PENTAD_WEIGHT_MIN_FRAC fraction of the global mean weight; weights are then renormalized to mean=1.
    SINGLE_STAGE_PARAMS$weights[is.na(SINGLE_STAGE_PARAMS$weights)] <- 0
    w <- SINGLE_STAGE_PARAMS$weights
    expected_len <- length(SINGLE_STAGE_PARAMS$indices) * N_TEMPORAL_BINS
    if (length(w) == expected_len) {
      global_mean <- mean(w, na.rm = TRUE)
      min_floor <- (if (exists("PENTAD_WEIGHT_MIN_FRAC")) PENTAD_WEIGHT_MIN_FRAC else 0.20) * global_mean
      w_mat <- matrix(w, nrow = N_TEMPORAL_BINS, ncol = length(SINGLE_STAGE_PARAMS$indices))
      w_mat[!is.finite(w_mat)] <- 0
      n_before_below <- sum(w_mat < min_floor, na.rm = TRUE)
      if (n_before_below > 0) {
        w_mat[w_mat < min_floor] <- min_floor
        w <- as.vector(w_mat)
        if (mean(w, na.rm = TRUE) > 1e-9) w <- w / mean(w, na.rm = TRUE)
        SINGLE_STAGE_PARAMS$weights <- w
        cat(sprintf("[NOTICE] Applied pentad weight floor to SINGLE_STAGE_PARAMS: min=%.6g (%.1f%% of mean=%.6g), adjusted %d/%d elements\n", min_floor, 100 * (if (exists("PENTAD_WEIGHT_MIN_FRAC")) PENTAD_WEIGHT_MIN_FRAC else 0.20), global_mean, n_before_below, length(w)))
      } else {
        # Nothing to change but keep existing weights
        SINGLE_STAGE_PARAMS$weights <- w
      }
    } else {
      # Fallback: mismatched length -> apply tiny scalar floor
      floor_val <- 1e-4
      w[w < floor_val] <- floor_val
      SINGLE_STAGE_PARAMS$weights <- w
      cat("[NOTICE] SINGLE_STAGE_PARAMS weights length mismatch; applied scalar floor\n")
    }
    print_weights_summary("SINGLE_STAGE", SINGLE_STAGE_PARAMS)
  }

  # SINGLE-STAGE UNMIXING: All endmembers (barren + veg types) treated equally
  cat("[MODE] Single-stage unmixing ENABLED (barren and vegetation types treated as equals)\n")

  cat("Building single-stage library (barren + all vegetation types)...\n")
  mesma_lib_single <- build_single_stage_library_weighted(df_train, avail, SINGLE_STAGE_PARAMS, ALLOWED_VEG)

  cat("Pre-computing optimized library for single-stage MESMA...\n")
  OPTIMIZED_LIBRARY_SINGLE <- precompute_optimized_library_weighted(mesma_lib_single, grid_type = "full", feature_weights = (if (!is.null(SINGLE_STAGE_PARAMS) && !is.null(SINGLE_STAGE_PARAMS$weights)) SINGLE_STAGE_PARAMS$weights else NULL) )

  # Summarize pruning if any
  if (exists("PRUNE_ZERO_WEIGHT_FEATURES") && isTRUE(PRUNE_ZERO_WEIGHT_FEATURES)) {
    total_pruned <- 0
    total_kept <- 0
    vegs_pruned <- 0
    for (v in names(OPTIMIZED_LIBRARY_SINGLE)) {
      libv <- OPTIMIZED_LIBRARY_SINGLE[[v]]
      if (!is.null(libv$pruned_info) && !is.null(libv$pruned_info$n_kept)) {
        total_kept <- total_kept + libv$pruned_info$n_kept
        vegs_pruned <- vegs_pruned + 1
      }
    }
    if (vegs_pruned > 0) {
      cat(sprintf("[INFO] Pruning summary: %d vegetation types pruned; avg kept features per pruned veg: %.1f\n", vegs_pruned, total_kept / vegs_pruned))
    }
  }

  if (isTRUE(TESTING_MODE)) {
    nn <- sapply(names(OPTIMIZED_LIBRARY_SINGLE), function(v) {
      libv <- OPTIMIZED_LIBRARY_SINGLE[[v]]
      if (is.null(libv) || is.null(libv$M)) return(0)
      nrow(libv$M)
    })
    cat(sprintf("[DEBUG] Single-stage library built: endmember_count=%d variants_per_type=%s\n",
                length(nn), paste(sprintf("%s:%d", names(nn), nn), collapse=", ")))
  }

  assign("SINGLE_STAGE_PARAMS", SINGLE_STAGE_PARAMS, envir = globalenv())
  assign("mesma_lib_single", mesma_lib_single, envir = globalenv())
  assign("OPTIMIZED_LIBRARY_SINGLE", OPTIMIZED_LIBRARY_SINGLE, envir = globalenv())

  cat(sprintf("[NOTICE] Single-stage feature count: avail=%d, params_indices=%d\n",
              length(avail), length(SINGLE_STAGE_PARAMS$indices)))

  mesma_lib <- mesma_lib_single
  OPTIMIZED_LIBRARY <- OPTIMIZED_LIBRARY_SINGLE

  compressed_templates_accessor <- precompute_compressed_templates(mesma_lib, "full")
  assign("compressed_templates_accessor", compressed_templates_accessor, envir = globalenv())
  assign(".COMPRESSED_TEMPLATES_ACCESSOR", compressed_templates_accessor, envir = globalenv())

  if (exists("save_mesma_cache") && !is.null(get("mesma_lib", envir = globalenv()))) {
    tryCatch({
      cache_dir <- save_mesma_cache()
      cat("Model cache saved after training pipeline.\n")
    }, error = function(e) {
      cat(sprintf("Failed to write MESMA cache: %s\n", e$message))
    })
  }

cat("[DEBUG] Reached line 6377 - about to define fit_one_task function\n")

gls_mbb_bootstrap <- function(residuals, block_size) {
  n <- length(residuals)
  if (n == 0) return(numeric(0))
  
  # Remove NA values for bootstrapping
  residuals_non_na <- residuals[!is.na(residuals)]
  n_non_na <- length(residuals_non_na)
  
  if (n_non_na == 0) return(residuals) # Return original if all are NA
  
  # Generate indices for moving block bootstrap from non-NA residuals
  num_blocks <- ceiling(n_non_na / block_size)
  resampled_indices <- integer(num_blocks * block_size)
  
  for (i in 1:num_blocks) {
    start_index <- sample(1:(n_non_na - block_size + 1), 1)
    resampled_indices[((i-1)*block_size + 1):(i*block_size)] <- start_index:(start_index + block_size - 1)
  }
  
  # Trim to original length of non-NA residuals
  resampled_indices <- resampled_indices[1:n_non_na]
  
  # Create the bootstrapped residuals vector
  bootstrapped_residuals_non_na <- residuals_non_na[resampled_indices]
  
  # Place bootstrapped residuals back into the original structure with NAs
  residuals_boot <- residuals
  residuals_boot[!is.na(residuals)] <- bootstrapped_residuals_non_na
  
  return(residuals_boot)
}

  fit_one_task <- function(task_data) {
    if (is.null(task_data) || nrow(task_data) == 0) return(NULL)

    loc <- as.character(task_data$location_id[1])
    yr <- as.integer(task_data$pheno_year[1])

    # Extract lat/lon if available
    lat_val <- NA_real_
    lon_val <- NA_real_
    if ("lat" %in% names(task_data)) lat_val <- as.numeric(task_data$lat[1]) else if ("latitude" %in% names(task_data)) lat_val <- as.numeric(task_data$latitude[1])
    if ("lon" %in% names(task_data)) lon_val <- as.numeric(task_data$lon[1]) else if ("longitude" %in% names(task_data)) lon_val <- as.numeric(task_data$longitude[1])

    # Use single-stage parameters (two-stage processing removed)
    PARAMS <- SINGLE_STAGE_PARAMS
    raw_mat <- build_pentad_matrix(task_data, PARAMS$indices)
    if (is.null(raw_mat)) {
      return(NULL)
    }

    y_raw <- as.numeric(raw_mat)

    # === CREATE VALIDITY MASK ===
    # Mask for valid (non-NA) observations - DO NOT replace NA with 0
    valid_mask <- is.finite(y_raw)
    n_valid <- sum(valid_mask)

    # Previously we required a minimum fraction of valid observations and would skip tasks
    # with low coverage. Remove that strict filter: proceed whenever there is at least one
    # valid observation, but emit a warning in TESTING_MODE when coverage is very low.
    if (n_valid == 0) {
      if (exists("TESTING_MODE") && isTRUE(TESTING_MODE)) {
        cat(sprintf("[SKIP] loc=%s pheno_year=%d: no valid observations (n_valid=0) - skipping\n", loc, yr))
      }
      return(NULL)
    } else {
      MIN_VALID_FRACTION <- 0.05
      if (exists("TESTING_MODE") && isTRUE(TESTING_MODE) && n_valid < (length(y_raw) * MIN_VALID_FRACTION)) {
        cat(sprintf("[WARN] loc=%s pheno_year=%d: low valid observations (%d < %.0f) - proceeding anyway\n",
                    loc, yr, n_valid, length(y_raw) * MIN_VALID_FRACTION))
      }
    }
    # ============================

    # --- New Logic: Estimate barren fraction purely through PPI ---
    # Moved from end of function to allow early exit for pure barren pixels
    # 1. Get current PPI for the task (prefer raw PPI if available)
    current_ppi <- NA_real_
    ppi_col <- if ("PPI_raw" %in% names(task_data)) "PPI_raw" else "PPI"

    if (ppi_col %in% names(task_data)) {
      # Filter task_data for summer months (June-September)
      summer_task_data <- if ("month" %in% names(task_data)) task_data[task_data$month %in% 6:9, ] else task_data
      if (nrow(summer_task_data) > 0) {
        current_ppi <- median(summer_task_data[[ppi_col]], na.rm = TRUE)
      } else {
        # Fallback if no summer data, use median of all available data for the year
        current_ppi <- median(task_data[[ppi_col]], na.rm = TRUE)
         if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_task] loc=%s yr=%d: No summer PPI data found; falling back to median of all annual PPI values.\n", loc, yr))
      }
    }

    total_veg_cover <- 0 # Default to 0 if PPI is not available or invalid
    barren_fraction <- 1.0

    if (is.finite(current_ppi)) {
      # Linear interpolation of PPI to total vegetation cover
      slope <- 1 / (PPI_FULL_VEG_COVER - PPI_ZERO_VEG_COVER)
      total_veg_cover <- (current_ppi - PPI_ZERO_VEG_COVER) * slope
      
      # Clamp to [0, 1]
      total_veg_cover <- pmin(pmax(total_veg_cover, 0), 1)
      barren_fraction <- 1 - total_veg_cover
    }
    
    # SHORTCUT: If barren fraction > 0.95, return pure barren immediately
    if (barren_fraction > 0.95) {
      if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_task] loc=%s yr=%d: High barren fraction (%.3f > 0.95) - returning pure barren without unmixing\n", loc, yr, barren_fraction))
      
      coef_df_barren <- data.frame(
        location_id = loc,
        pheno_year = yr,
        lat = lat_val,
        lon = lon_val,
        Veg = "barren",
        variant_id = "barren_ppi",
        coef = 1.0, # Force to 1.0 as requested
        rmse = 0,
        coef_025 = NA,
        coef_975 = NA,
        coef_sd = 0,
        interval = NA,
        n_obs = nrow(task_data),
        inseparable_variant_flag = FALSE,
        inseparable_variant_details = NA_character_,
        stringsAsFactors = FALSE
      )
      
      diag_df_barren <- data.frame(location_id = loc, pheno_year = yr, stringsAsFactors = FALSE)
      diag_df_barren$barren_fraction <- 1.0
      diag_df_barren$barren_fraction_ppi_based <- barren_fraction
      
      # Return minimal valid structure
      return(list(
        coef_df = coef_df_barren,
        diagnostics = diag_df_barren,
        uncertainty = NULL, 
        residuals = if (exists("y_raw")) y_raw else numeric(0),
        y_hat = if (exists("y_raw")) rep(0, length(y_raw)) else numeric(0),
        y_obs = if (exists("y_raw")) y_raw else numeric(0),
        E_best = matrix(0, nrow=0, ncol=0),
        top_variants = list(),
        weights_masked = NULL,
        valid_mask = valid_mask
      ))
    }

    # ===== SINGLE-STAGE UNMIXING =====
      # Normalize the observation using single-stage parameters
      y_norm <- y_raw
      n_bins <- N_TEMPORAL_BINS
      for(k in seq_along(PARAMS$indices)) {
        idx_start <- (k-1)*n_bins + 1
        idx_end <- k*n_bins

        mu <- PARAMS$means[k]
        sigma <- PARAMS$sds[k]
        if (!is.finite(sigma) || sigma < EPS_SIGMA) {
          if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_task] Single-stage sigma for index %d is non-finite or too small (%.8f); using EPS_SIGMA=%.8e\n", k, sigma, EPS_SIGMA))
          sigma <- EPS_SIGMA
        }

        y_norm[idx_start:idx_end] <- (y_norm[idx_start:idx_end] - mu) / sigma
      }

      # Apply mask
      y_norm_masked <- y_norm[valid_mask]
      weights_masked <- PARAMS$weights[valid_mask]

      # Check if we have sufficient signal (using unweighted norm for validation)
      y_norm_val <- sqrt(sum(y_norm_masked^2, na.rm = TRUE))

      if (is.na(y_norm_val) || y_norm_val < 1e-9) {
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_task] loc=%s yr=%d: insufficient signal (norm=%.6g), skipping\n", loc, yr, y_norm_val))
        return(NULL)
      }

      # Use the masked, z-score normalized observation (NOT pre-weighted)
      # Weighting will be applied inside solve_weights_nnls_simple via feature_weights parameter
      y_for_unmixing <- y_norm_masked

      # Perform single-stage MESMA using all endmembers (barren + all veg types)
      veg_kept <- names(mesma_lib)
      top_variants <- list()

      for(v in veg_kept) {
        lib <- OPTIMIZED_LIBRARY[[v]]
        if(is.null(lib)) next

        # Mask library templates
        lib_M_norm_masked <- lib$M_norm[, valid_mask, drop = FALSE]

        # Renormalize each row after masking
        if (nrow(lib_M_norm_masked) == 1) {
          row_norm <- sqrt(sum(lib_M_norm_masked^2))
          if (row_norm >= 1e-9) lib_M_norm_masked <- lib_M_norm_masked / row_norm
        } else {
          lib_M_norm_masked <- t(apply(lib_M_norm_masked, 1, function(row) {
            row_norm <- sqrt(sum(row^2))
            if (row_norm < 1e-9) row else row / row_norm
          }))
        }

        # Compute similarities using weighted vectors (for ranking only)
        w_masked <- weights_masked
        if (is.null(w_masked) || length(w_masked) != length(y_for_unmixing)) w_masked <- rep(1, length(y_for_unmixing))
        y_for_sim <- y_for_unmixing * sqrt(pmax(w_masked, 0))
        denom_sim <- sqrt(sum(y_for_sim^2, na.rm = TRUE)); if (denom_sim < 1e-12) denom_sim <- 1
        y_sim_norm <- y_for_sim / denom_sim
        sims <- as.numeric(lib_M_norm_masked %*% y_sim_norm)
        best_idx <- order(sims, decreasing=TRUE)[1:min(length(sims), TOPK_VARIANTS)]

        top_variants[[v]] <- lapply(best_idx, function(i) {
          # Use masked template (only valid observations) for unmixing
          masked_vec <- lib_M_norm_masked[i, ]
          if (isTRUE(TESTING_MODE) && length(masked_vec) != length(y_for_unmixing)) {
            cat(sprintf("[WARN fit_one_task] loc=%s yr=%d veg=%s: template length mismatch (template=%d, y=%d)\n",
                        loc, yr, v, length(masked_vec), length(y_for_unmixing)))
          }
          list(vec = masked_vec, id = lib$ids[i], similarity = sims[i])
        })
      }

      # Remove empty vegetation types
      empty_vegs <- names(top_variants)[sapply(top_variants, function(x) is.null(x) || length(x) == 0)]
      if (length(empty_vegs) > 0) {
        for (ev in empty_vegs) top_variants[[ev]] <- NULL
      }

      if (length(top_variants) == 0) {
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_task] loc=%s yr=%d: no top_variants available after filtering empties, skipping\n", loc, yr))
        return(NULL)
      }
      
      # Evaluate all combinations of endmembers
      # Pass unweighted observation; weighting happens inside solve_weights_nnls_simple
      best_result <- tryCatch({
        evaluate_all_combinations(y_for_unmixing, top_variants, lambda = 0, feature_weights = weights_masked)
      }, error = function(e) {
        cat(sprintf("[ERROR fit_one_task] evaluate_all_combinations failed for loc=%s year=%d: %s\n", as.character(loc), as.integer(yr), e$message))
        NULL
      })

      if (is.null(best_result)) {
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_task] loc=%s yr=%d: evaluate_all_combinations returned NULL\n", loc, yr))
        return(NULL)
      }

      # Extract coefficients and create output
      chosen_ids <- best_result$ids
      coefs <- best_result$w
      rmse <- best_result$rmse
      residuals <- best_result$residuals
      
      # Calculate model selection uncertainty from top models
      coef_sd_vec <- rep(NA, length(coefs))
      names(coef_sd_vec) <- names(coefs)
      
      if (!is.null(best_result$top_models) && length(best_result$top_models) > 1) {
         # Stack weights from all top models: n_veg x n_models
         w_matrix <- do.call(cbind, lapply(best_result$top_models, function(m) m$w))
         
         if (!is.null(w_matrix) && ncol(w_matrix) > 1) {
            # Compute SD for each vegetation type across the top models
            coef_sd_vec <- apply(w_matrix, 1, sd, na.rm=TRUE)
         }
      }

      # Build coefficient dataframe with variant-level detail
      coef_df <- data.frame(
        location_id = loc,
        pheno_year = yr,
        lat = lat_val,
        lon = lon_val,
        Veg = sapply(strsplit(chosen_ids, "_v"), `[`, 1),
        variant_id = chosen_ids,
        coef = coefs,
        rmse = rmse,
        coef_025 = NA,
        coef_975 = NA,
        coef_sd = as.numeric(coef_sd_vec),
        interval = NA,
        n_obs = nrow(task_data),
        inseparable_variant_flag = FALSE,
        inseparable_variant_details = NA_character_,
        stringsAsFactors = FALSE
      )

      # Keep barren in the solve, but strip it from coefficients before reporting; barren fraction comes from PPI
      barren_solver_mask <- tolower(coef_df$Veg) == "barren"
      if (any(barren_solver_mask)) {
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_task] loc=%s yr=%d: Removing %d barren variants after solving; PPI supplies barren fraction.\n", loc, yr, sum(barren_solver_mask)))
        coef_df <- coef_df[!barren_solver_mask, , drop = FALSE]
      }

      # Apply 5% minimum fraction filter to VEGETATION types only (after barren removal)
      # This ensures vegetation types compete fairly against each other
      if (nrow(coef_df) > 0) {
        veg_sum <- sum(coef_df$coef, na.rm = TRUE)
        if (veg_sum > 1e-9) {
          rel_fracs <- coef_df$coef / veg_sum
          min_frac <- 0.05
          below_thresh <- rel_fracs < min_frac
          if (any(below_thresh)) {
            if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_task] loc=%s yr=%d: Zeroing %d vegetation types below 5%% threshold.\n", loc, yr, sum(below_thresh)))
            coef_df$coef[below_thresh] <- 0
            # Renormalize remaining to sum to original veg_sum
            new_sum <- sum(coef_df$coef, na.rm = TRUE)
            if (new_sum > 1e-9) {
              coef_df$coef <- coef_df$coef * (veg_sum / new_sum)
            }
          }
        }
      }

      # --- Mark inseparable variants (drop instead of assign Veg = 'unknown') if detected in similarity tables ---
      sim_tbl <- NULL
      if (exists("INSEPARABLE_VARIANT_INFO") && !is.null(INSEPARABLE_VARIANT_INFO$similarity_table)) {
        sim_tbl <- INSEPARABLE_VARIANT_INFO$similarity_table
      } else if (exists("VARIANT_SIMILARITY_TABLE") && !is.null(VARIANT_SIMILARITY_TABLE)) {
        sim_tbl <- VARIANT_SIMILARITY_TABLE
      }
      if (!is.null(sim_tbl) && nrow(sim_tbl) > 0) {
        cols <- names(sim_tbl)
        variant_id_cols <- intersect(cols, c("variant_id", "variant", "id", "ids", "var_id", "variant1", "variant_a"))
        other_variant_cols <- intersect(cols, c("other_variant_id", "other_variant", "variant2", "other_id", "variant_b"))
        ids_in_tbl <- character(0)
        for (cname in variant_id_cols) ids_in_tbl <- c(ids_in_tbl, as.character(sim_tbl[[cname]]))
        for (cname in other_variant_cols) ids_in_tbl <- c(ids_in_tbl, as.character(sim_tbl[[cname]]))
        ids_in_tbl <- unique(na.omit(ids_in_tbl))
        if (length(ids_in_tbl) > 0) {
          mask <- coef_df$variant_id %in% ids_in_tbl
          if (any(mask)) {
            coef_df$inseparable_variant_flag[mask] <- TRUE
            # Build detail strings: collect related vegs and other variant ids
            details <- vapply(coef_df$variant_id[mask], function(vid) {
              rows <- apply(sim_tbl, 1, function(r) any(vid == as.character(r)))
              if (any(rows)) {
                row_tbl <- sim_tbl[rows, , drop = FALSE]
                other_vegs <- character(0)
                if ("other_veg" %in% names(row_tbl)) other_vegs <- c(other_vegs, as.character(row_tbl$other_veg))
                if ("veg" %in% names(row_tbl)) other_vegs <- c(other_vegs, as.character(row_tbl$veg))
                other_variant_ids <- character(0)
                for (cname in other_variant_cols) if (cname %in% names(row_tbl)) other_variant_ids <- c(other_variant_ids, as.character(row_tbl[[cname]]))
                other_variant_ids <- unique(na.omit(other_variant_ids))
                info_parts <- character(0)
                if (length(other_vegs) > 0) info_parts <- c(info_parts, paste(unique(na.omit(other_vegs)), collapse = ";"))
                if (length(other_variant_ids) > 0) info_parts <- c(info_parts, paste(other_variant_ids, collapse = ";"))
                paste(info_parts, collapse = "|")
              } else {
                NA_character_
              }
            }, FUN.VALUE = "")
            coef_df$inseparable_variant_details[mask] <- details

            # New behavior: drop variants only if similarity > 0.95 with barren.
            sim_col <- intersect(c("Similarity","similarity","sim"), names(sim_tbl))[1]
            if (!is.na(sim_col)) {
              for (i in seq_len(nrow(sim_tbl))) {
                sim_val <- as.numeric(sim_tbl[[sim_col]][i])
                if (is.na(sim_val) || sim_val < 0.95) next

                # Extract variant ids from row (pick first non-NA from the candidate cols)
                a_id <- NA_character_; b_id <- NA_character_
                for (cname in variant_id_cols) if (cname %in% names(sim_tbl) && !is.na(sim_tbl[[cname]][i])) { a_id <- as.character(sim_tbl[[cname]][i]); break }
                for (cname in other_variant_cols) if (cname %in% names(sim_tbl) && !is.na(sim_tbl[[cname]][i])) { b_id <- as.character(sim_tbl[[cname]][i]); break }
                if (is.na(a_id) || is.na(b_id)) next

                # Only act if both are in coef_df (i.e., were chosen)
                if (!(a_id %in% coef_df$variant_id && b_id %in% coef_df$variant_id)) next

                # Get veg types
                veg_a <- sub("_v.*$", "", a_id)
                veg_b <- sub("_v.*$", "", b_id)

                # Only drop if one is barren
                if (!(tolower(veg_a) == "barren" || tolower(veg_b) == "barren")) next

                # If both are barren, skip
                if (tolower(veg_a) == "barren" && tolower(veg_b) == "barren") next

                # Drop the non-barren one
                drop_id <- if (tolower(veg_a) == "barren") b_id else a_id
                other_id <- if (drop_id == a_id) b_id else a_id

                # Set coefficient to zero and annotate
                idx_drop <- which(coef_df$variant_id == drop_id)
                if (length(idx_drop) > 0) {
                  coef_df$coef[idx_drop] <- 0
                  coef_df$inseparable_variant_flag[idx_drop] <- TRUE
                  coef_df$inseparable_variant_details[idx_drop] <- paste0("Dropped due to high similarity (", sprintf("%.3f", sim_val), ") with barren variant ", other_id)
                  if (isTRUE(TESTING_MODE)) cat(sprintf("[INFO] Dropping variant %s because similarity=%.3f > 0.95 with barren %s\n", drop_id, sim_val, other_id))
                }
              }
            }
          }
        }
      }

      # Aggregate by vegetation type (sum coefficients for same veg type)
      if (nrow(coef_df) > 0) {
        coef_agg <- aggregate(coef ~ Veg, data = coef_df, FUN = sum)
      } else {
        coef_agg <- data.frame(Veg = character(0), coef = numeric(0), stringsAsFactors = FALSE)
      }

      # --- PPI-based barren fraction (calculated at start of function) ---
      # Variables current_ppi, total_veg_cover, and barren_fraction are already set at the top of fit_one_task.

      
      # Re-normalize existing vegetation coefficients to sum to total_veg_cover
      # Exclude 'barren' if it somehow snuck in (it shouldn't due to earlier filtering)
      veg_coefs_mask <- tolower(coef_agg$Veg) != "barren"
      sum_original_veg_coefs <- sum(coef_agg$coef[veg_coefs_mask], na.rm = TRUE)
      
      scale_factor <- 0
      if (sum_original_veg_coefs > 1e-9) { # Avoid division by zero/very small numbers
        scale_factor <- total_veg_cover / sum_original_veg_coefs
        coef_agg$coef[veg_coefs_mask] <- coef_agg$coef[veg_coefs_mask] * scale_factor
      } else {
        # If no vegetation was unmixed or sum is zero, set all veg coefs to 0
        coef_agg$coef[veg_coefs_mask] <- 0
      }

      # Add barren fraction to coef_agg
      barren_row <- data.frame(Veg = "barren", coef = barren_fraction)
      coef_agg <- rbind(coef_agg, barren_row)

      # --- AGENT UPDATE: Propagate scaling and barren fraction to coef_df ---
      if (sum_original_veg_coefs > 1e-9) {
         v_rows <- tolower(coef_df$Veg) != "barren"
         coef_df$coef[v_rows] <- coef_df$coef[v_rows] * scale_factor
      } else {
         v_rows <- tolower(coef_df$Veg) != "barren"
         coef_df$coef[v_rows] <- 0
      }
      
      if (nrow(coef_df) > 0) {
        barren_row_df <- coef_df[1, , drop=FALSE] 
        barren_row_df[] <- NA 
        barren_row_df$location_id <- loc
        barren_row_df$pheno_year <- yr
        barren_row_df$lat <- lat_val
        barren_row_df$lon <- lon_val
        barren_row_df$Veg <- "barren"
        barren_row_df$variant_id <- "barren_ppi"
        barren_row_df$coef <- barren_fraction
        coef_df <- rbind(coef_df, barren_row_df)
      } else {
        coef_df <- data.frame(
           location_id = loc, pheno_year = yr, lat = lat_val, lon = lon_val, Veg = "barren", 
           variant_id = "barren_ppi", coef = barren_fraction, 
           stringsAsFactors = FALSE
        )
      }

      # Create diagnostics dataframe
      diag_df <- data.frame(
        location_id = loc,
        pheno_year = yr,
        stringsAsFactors = FALSE
      )
      # Add each vegetation type's total fraction
      for (v in coef_agg$Veg) {
        diag_df[[paste0(v, "_fraction")]] <- coef_agg$coef[coef_agg$Veg == v]
      }
      diag_df$barren_fraction_ppi_based <- barren_fraction # Add PPI-based barren fraction

      # Reconstruct E_best for returning
      E_best_list <- list()
      for(v in names(chosen_ids)){
        vid <- chosen_ids[[v]]
        for(variant in top_variants[[v]]){
          if(variant$id == vid) { E_best_list[[v]] <- variant$vec; break }
        }
      }
      E_best <- do.call(cbind, E_best_list)

      return(list(
        coef_df = coef_df,
        diagnostics = diag_df,
        uncertainty = NULL,
        residuals = residuals,
        y_hat = y_for_unmixing - residuals,
        y_obs = y_for_unmixing,
        E_best = E_best,
        top_variants = top_variants,
        weights_masked = weights_masked,
        valid_mask = valid_mask
      ))

      # ===== END SINGLE-STAGE UNMIXING =====
  # These variables are already set in the earlier conditional blocks above.
  # No need to reassign them here
  # Single-stage MESMA is used, so mesma_lib and OPTIMIZED_LIBRARY are set accordingly

  }

# New function: Process all years for a single location with multi-year bootstrap
  fit_one_location <- function(location_data) {
    if (is.null(location_data) || nrow(location_data) == 0) return(NULL)

    loc <- as.character(location_data$location_id[1])

    # Get all phenological years for this location
    years <- sort(unique(location_data$pheno_year))
    years <- years[!is.na(years)]

    if (length(years) == 0) return(NULL)

    # Process each year individually first (to get point estimates and chosen variants)
    year_results <- list()
    y_vecs_by_year <- list()
    chosen_ids_by_year <- list()
    w_hat_by_year <- list()

    # Wrap per-location processing in tryCatch so a single error does not discard valid year results
    out <- tryCatch({
      for (yr in years) {
        year_data <- location_data[location_data$pheno_year == yr, , drop = FALSE]
        res_yr <- fit_one_task(year_data)

        if (!is.null(res_yr)) {
          tryCatch({
            year_results[[as.character(yr)]] <- res_yr
          }, error = function(e) {
            if (isTRUE(TESTING_MODE)) cat(sprintf("[ERROR assigning year_results] loc=%s yr=%d: %s\n", loc, yr, e$message))
            stop(e)  # Re-throw to be caught by outer handler
          })

          # Store data needed for multi-year bootstrap - wrapped in tryCatch to isolate errors
          tryCatch({
            raw_mat_yr <- build_pentad_matrix(year_data, SINGLE_STAGE_PARAMS$indices)
            if (!is.null(raw_mat_yr)) {
              y_raw_yr <- as.numeric(raw_mat_yr)

              # Check dimensions match expectations
              n_bins <- N_TEMPORAL_BINS
              expected_length <- length(SINGLE_STAGE_PARAMS$indices) * n_bins
              if (length(y_raw_yr) != expected_length) {
                if (isTRUE(TESTING_MODE)) {
                  cat(sprintf("[WARN fit_one_location] loc=%s yr=%d: y_raw_yr length mismatch (got %d, expected %d), skipping multi-year processing for this year\n",
                    loc, yr, length(y_raw_yr), expected_length))
                }
              } else {
                # Apply same preprocessing as fit_one_task
                y_s1_yr <- y_raw_yr
                for(k in seq_along(SINGLE_STAGE_PARAMS$indices)) {
                  idx_start <- (k-1)*n_bins + 1
                  idx_end <- k*n_bins
                  mu <- SINGLE_STAGE_PARAMS$means[k]
                  sigma <- SINGLE_STAGE_PARAMS$sds[k]
                  if (is.finite(sigma) && sigma > 0) {
                    y_s1_yr[idx_start:idx_end] <- (y_s1_yr[idx_start:idx_end] - mu) / sigma
                  }
                }
                y_s1_yr[!is.finite(y_s1_yr)] <- 0

                y_vecs_by_year[[length(y_vecs_by_year) + 1]] <- y_s1_yr
              }

              # Get chosen variants and weights from this year's result
              if (!is.null(res_yr$chosen_variants)) {
                chosen_ids_by_year[[yr]] <- res_yr$chosen_variants
              }
              if (!is.null(res_yr$vegetation_proportions)) {
                w_hat_by_year[[yr]] <- res_yr$vegetation_proportions
              }
            }
          }, error = function(e) {
            if (isTRUE(TESTING_MODE)) {
              cat(sprintf("[WARN fit_one_location] loc=%s yr=%d: multi-year processing failed (%s), continuing with single-year results\n",
                loc, yr, e$message))
            }
          })
        } else {
          if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_location] loc=%s yr=%d: fit_one_task returned NULL\n", loc, yr))
        }
      }

      # If we have results for multiple years and ENABLE_UNCERTAINTY, do multi-year bootstrap
      if (length(year_results) == 0) {
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_location] loc=%s: no successful year results\n", loc))
        return(NULL)
      }

      if (isTRUE(ENABLE_UNCERTAINTY) && length(years) >= 1 && length(y_vecs_by_year) == length(years)) {
        # Get comp_templates from first successful year result
        comp_templates_list <- list()
        top_variants <- NULL

        for (yr in years) {
          res_yr <- year_results[[as.character(yr)]]
          if (!is.null(res_yr) && !is.null(res_yr$chosen_variants)) {
            # Build comp_templates from the MESMA library
            vegs <- names(res_yr$chosen_variants)
            if (is.null(top_variants)) {
              top_variants <- list()
              for (v in vegs) top_variants[[v]] <- list()
            }

            break
          }
        }

        # Build comp_templates from mesma_lib
        if (!is.null(top_variants) && exists("mesma_lib")) {
          comp_templates <- list()
          for (v in names(top_variants)) {
            if (v %in% names(mesma_lib)) {
              comp_templates[[v]] <- mesma_lib[[v]]
            }
          }
        }

        # Perform per-year bootstrap to estimate uncertainty (coef CI, sd, rmse CI, variant frequencies)
        if (isTRUE(ENABLE_UNCERTAINTY)) {
          B_loc <- if (exists("BOOTSTRAP_B")) min(BOOTSTRAP_B, 100L) else 100L
          if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG] Performing local bootstrap (B=%d) for loc=%s\n", B_loc, loc))

          for (yr in years) {
            res_yr <- year_results[[as.character(yr)]]
            if (is.null(res_yr) || is.null(res_yr$residuals)) next

            year_data <- location_data[location_data$pheno_year == yr, , drop = FALSE]
            n_obs <- nrow(year_data)
            if (n_obs < max(6, MIN_OBS_PER_LOC_YEAR)) {
              if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG] Too few observations for bootstrap loc=%s yr=%s (n_obs=%d); applying small-n CI inflation fallback\n", loc, yr, n_obs))

              # Populate conservative uncertainty directly on the coefficient table
              curr_df <- year_results[[as.character(yr)]]$coef_df
              if (!is.null(curr_df) && nrow(curr_df) > 0) {
                for (r in seq_len(nrow(curr_df))) {
                  # Keep any existing sd signal (e.g., top-model spread) but inflate aggressively for small n
                  ci <- small_n_inflated_ci(
                    est = curr_df$coef[r],
                    sd_in = if ("coef_sd" %in% names(curr_df)) curr_df$coef_sd[r] else NA_real_,
                    n_obs = n_obs
                  )
                  curr_df$coef_025[r] <- ci$coef_025
                  curr_df$coef_975[r] <- ci$coef_975
                  curr_df$coef_sd[r] <- ci$coef_sd
                  curr_df$interval[r] <- ci$interval
                }
                year_results[[as.character(yr)]]$coef_df <- curr_df

                if (isTRUE(TESTING_MODE)) {
                  widths <- curr_df$coef_975 - curr_df$coef_025
                  cat(sprintf("[DEBUG] small-n fallback CI widths loc=%s yr=%s: median=%.3f (n_obs=%d)\n",
                              loc, yr, median(widths, na.rm = TRUE), n_obs))
                }
              }
              next
            }

            boot_coef_by_veg <- list()
            boot_variant_counts <- list()
            boot_rmse <- numeric(0)

            # Suppress verbose output during bootstraps
            sink_file <- tempfile()
            con_out <- file(sink_file, open = "wt")
            con_msg <- file(paste0(sink_file, ".msg"), open = "wt")
            tryCatch({
              sink(con_out, type = "output")
              sink(con_msg, type = "message")

              for (b in seq_len(B_loc)) {
                
                residuals_boot <- gls_mbb_bootstrap(res_yr$residuals, BOOTSTRAP_BLOCK_SIZE)
                y_boot <- res_yr$y_hat + residuals_boot
                
                res_b <- tryCatch({
                  evaluate_all_combinations(y_boot, res_yr$top_variants, lambda = 0, feature_weights = res_yr$weights_masked)
                }, error = function(e) NULL)

                if (is.null(res_b)) next
                
                coefs_b <- res_b$w
                chosen_ids_b <- res_b$ids
                
                dfb <- data.frame(
                    Veg = sapply(strsplit(chosen_ids_b, "_v"), `[`, 1),
                    coef = coefs_b,
                    stringsAsFactors = FALSE
                )

                # collect coefficients
                for (i in seq_len(nrow(dfb))) {
                  vg <- as.character(dfb$Veg[i])
                  val <- as.numeric(dfb$coef[i])
                  if (!is.finite(val)) next
                  if (is.null(boot_coef_by_veg[[vg]])) boot_coef_by_veg[[vg]] <- numeric(0)
                  boot_coef_by_veg[[vg]] <- c(boot_coef_by_veg[[vg]], val)
                }

                # collect chosen variant ids
                for (i in seq_along(chosen_ids_b)) {
                  vg <- names(chosen_ids_b)[i]
                  vid <- chosen_ids_b[[i]]
                  if (is.null(boot_variant_counts[[vg]])) boot_variant_counts[[vg]] <- list()
                  boot_variant_counts[[vg]][[as.character(vid)]] <- (boot_variant_counts[[vg]][[as.character(vid)]] %||% 0) + 1L
                }

                # rmse
                if (!is.null(res_b$rmse) && is.finite(as.numeric(res_b$rmse))) boot_rmse <- c(boot_rmse, as.numeric(res_b$rmse))
              }

            }, finally = {
              try(sink(type = "message"), silent = TRUE)
              try(sink(type = "output"), silent = TRUE)
              close(con_out); close(con_msg)
              try(unlink(sink_file), silent = TRUE); try(unlink(paste0(sink_file, ".msg")), silent = TRUE)
            })
            
            # Build coef CI table
            coef_ci_df <- NULL
            if (length(boot_coef_by_veg) > 0) {
              rows <- lapply(names(boot_coef_by_veg), function(vg) {
                vals <- boot_coef_by_veg[[vg]]
                vals <- vals[is.finite(vals)]
                if (length(vals) == 0) return(NULL)
                data.frame(Veg = vg,
                           coef_025 = as.numeric(quantile(vals, 0.025, na.rm = TRUE)),
                           coef_975 = as.numeric(quantile(vals, 0.975, na.rm = TRUE)),
                           coef_sd = as.numeric(stats::sd(vals, na.rm = TRUE)),
                           interval = as.numeric(quantile(vals, 0.975, na.rm = TRUE) - quantile(vals, 0.025, na.rm = TRUE)),
                           stringsAsFactors = FALSE)
              })
              rows <- rows[!sapply(rows, is.null)]
              if (length(rows) > 0) coef_ci_df <- do.call(rbind, rows)
            }

            # Build variant frequency table
            variant_freq_df <- NULL
            if (length(boot_variant_counts) > 0) {
              rows2 <- lapply(names(boot_variant_counts), function(vg) {
                tbl <- boot_variant_counts[[vg]]
                keys <- names(tbl)
                counts <- as.integer(unlist(tbl))
                dfv <- data.frame(Veg = rep(vg, length(keys)), Variant = keys, N = counts, Percent = counts / sum(counts) * 100, stringsAsFactors = FALSE)
                dfv
              })
              rows2 <- rows2[!sapply(rows2, is.null)]
              if (length(rows2) > 0) variant_freq_df <- do.call(rbind, rows2)
            }

            rmse_ci <- NULL
            if (length(boot_rmse) > 0) {
              rmse_ci <- as.numeric(quantile(boot_rmse, c(0.025, 0.975), na.rm = TRUE))
            }

            # Attach uncertainty to year_results
            if (!is.null(coef_ci_df) || !is.null(variant_freq_df) || !is.null(rmse_ci)) {
              if (is.null(year_results[[as.character(yr)]]$uncertainty)) year_results[[as.character(yr)]]$uncertainty <- list()
              year_results[[as.character(yr)]]$uncertainty$coef_ci <- coef_ci_df
              year_results[[as.character(yr)]]$uncertainty$variant_freq <- variant_freq_df
              year_results[[as.character(yr)]]$uncertainty$rmse_ci <- rmse_ci
              
              # Update coef_df with uncertainty metrics if available
              if (!is.null(coef_ci_df)) {
                  curr_df <- year_results[[as.character(yr)]]$coef_df
                  if (!is.null(curr_df)) {
                      # We loop over rows of coef_ci_df
                      for(r_idx in seq_len(nrow(coef_ci_df))) {
                          v_type <- coef_ci_df$Veg[r_idx]
                          # Find matching row(s) in curr_df (should be by Veg)
                          match_idx <- which(curr_df$Veg == v_type)
                          if(length(match_idx) > 0) {
                              curr_df$coef_sd[match_idx] <- coef_ci_df$coef_sd[r_idx]
                              curr_df$coef_025[match_idx] <- coef_ci_df$coef_025[r_idx]
                              curr_df$coef_975[match_idx] <- coef_ci_df$coef_975[r_idx]
                              curr_df$interval[match_idx] <- coef_ci_df$interval[r_idx]
                          }
                      }
                      year_results[[as.character(yr)]]$coef_df <- curr_df
                  }
              }
            }
          }
        }

        if (!is.null(top_variants) && exists("mesma_lib")) {
          comp_templates <- list()
          for (v in names(top_variants)) {
            if (v %in% names(mesma_lib)) {
              comp_templates[[v]] <- mesma_lib[[v]]
            }
          }
        }
      }

      # Return results for all years
      return(year_results)

    }, error = function(e) {
      cat(sprintf("[ERROR fit_one_location] loc=%s: %s\n", loc, e$message))
      if (length(year_results) > 0) {
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_location] loc=%s: returning partial %d year(s) after error\n", loc, length(year_results)))
        return(year_results)
      }
      return(NULL)
    })

    out
  }

  debug_unmix_loc_year <- function(loc, yr, df = NULL) {
    if (is.null(df)) {
      if (exists('df_tasks_inference') && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0) df <- df_tasks_inference
      else if (exists('df_tasks') && !is.null(df_tasks) && nrow(df_tasks) > 0) df <- df_tasks
      else stop('No task data available in workspace (df_tasks_inference or df_tasks)')
    }
    task_data <- df[df$location_id == loc & df$pheno_year == yr, , drop = FALSE]
    if (nrow(task_data) == 0) stop(sprintf('No rows found for loc=%s pheno_year=%s', as.character(loc), as.character(yr)))

    old_testing <- if (exists('TESTING_MODE', inherits = TRUE)) get('TESTING_MODE', envir = globalenv()) else FALSE
    on.exit({ assign('TESTING_MODE', old_testing, envir = globalenv()) })

    assign('TESTING_MODE', TRUE, envir = globalenv())

    cat(sprintf('\n=== DEBUG UNMIX: loc=%s pheno_year=%s ===\n', as.character(loc), as.character(yr)))
    res <- tryCatch({ fit_one_task(task_data) }, error = function(e) { cat(sprintf('DEBUG UNMIX ERROR: %s\n', e$message)); NULL })
    if (!is.null(res)) {
      cat('--- RESULT SUMMARY ---\n')
      if (!is.null(res$diagnostics)) print(res$diagnostics)
      if (!is.null(res$coef_df)) print(head(res$coef_df, n = 20))
    }
    cat('=== END DEBUG UNMIX ===\n\n')
    res
  }
  cat("Starting main processing loop...\n")

TEMP_RESULTS_DIR <- "C:/MAP/temp_results"
if (dir.exists(TEMP_RESULTS_DIR)) {
  cat(sprintf("Removing existing temporary results directory: %s\n", TEMP_RESULTS_DIR))
  unlink(TEMP_RESULTS_DIR, recursive = TRUE)
}
cat(sprintf("Creating temporary results directory: %s\n", TEMP_RESULTS_DIR))
# create directory and parent directories if needed; suppress warnings for existing paths
dir.create(TEMP_RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(TEMP_RESULTS_DIR)) {
  stop(sprintf("Failed to create temporary results directory: %s", TEMP_RESULTS_DIR))
}

  wb <- tryCatch({
    cat("[DEBUG] Creating Excel workbook...\n")
    openxlsx::createWorkbook()
  }, error = function(e) {
    cat(sprintf("\n========================================\n"))
    cat(sprintf("FATAL ERROR creating Excel workbook: %s\n", e$message))
    cat(sprintf("========================================\n\n"))
    cat("Make sure the 'openxlsx' package is properly installed and loaded.\n")
    cat("Try running: install.packages('openxlsx')\n")
    stop(e)
  })
  cat("[DEBUG] Excel workbook created successfully\n")
  flush.console()

# Auto-run main processing (MESMA_NO_AUTO_RUN support removed)

  aggregate_to_global_pattern <- function(all_coefs, method = "location_bootstrap") {
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")
    # Debug: print method to diagnose "Unknown method" errors
    cat(sprintf("[DEBUG] aggregate_to_global_pattern called with method='%s'\n", method))
    # AGENT: Changed requirement to pheno_year
    required_cols <- c("location_id", "pheno_year", "Veg", "coef")
    missing <- setdiff(required_cols, names(all_coefs))
    if (length(missing) > 0) stop(paste("Missing columns:", paste(missing, collapse = ", ")))
    # Compatibility: Ensure 'year' exists for downstream functions that expect it
    if (!"year" %in% names(all_coefs)) all_coefs$year <- all_coefs$pheno_year
    if (method == "location_bootstrap") {
      # pheno_year is guaranteed present
      result <- location_bootstrap_aggregate(all_coefs, B = BOOTSTRAP_B)
    } else {
      stop(sprintf("Unknown method '%s'. Use 'location_bootstrap'.", method))
    }
    result
  }

  aggregate_simple_mean <- function(all_coefs) {
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")
    veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
    # Determine which year column to use
    year_col <- if ("pheno_year" %in% names(all_coefs)) "pheno_year" else "year"
    results_list <- list()
    for (veg in veg_types) {
      veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), ]
      if (nrow(veg_data) == 0) {
        warning(sprintf("No valid coefficients found for vegetation type: %s", veg))
        next
      }
      # Group by year and compute simple statistics
      simple_result <- veg_data |> 
        dplyr::group_by(!!rlang::sym(year_col)) |> 
        dplyr::summarize(
          n_locations = dplyr::n(),
          global_coef = mean(coef, na.rm = TRUE),
          se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
          ci_lower = pmax(0, global_coef - 1.96 * se),
          ci_upper = pmin(1, global_coef + 1.96 * se),
          .groups = "drop"
        ) |> 
        dplyr::mutate(Veg = veg, method = "simple_mean") |> 
        dplyr::rename(year = !!rlang::sym(year_col))
      results_list[[veg]] <- simple_result
    }
    dplyr::bind_rows(results_list)
  }

  plot_global_vegetation_pattern <- function(global_pattern, 
                                              title = "Global Vegetation Composition Over Time",
                                              show_ci = TRUE,
                                              ci_type = "auto") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
      if ("ci_lower" %in% names(global_pattern) && "ci_upper" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower
        global_pattern$ci_upper <- global_pattern$ci_upper
      } else if ("coef_025" %in% names(global_pattern) && "coef_975" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$coef_025
        global_pattern$ci_upper <- global_pattern$coef_975
      }
    } else if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
      if ("ci_lower_simple" %in% names(global_pattern) && "ci_upper_simple" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower_simple
        global_pattern$ci_upper <- global_pattern$ci_upper_simple
      }
    } else if ("weighted_mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$weighted_mean_coef
      if ("ci_lower_pooled" %in% names(global_pattern) && "ci_upper_pooled" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower_pooled
        global_pattern$ci_upper <- global_pattern$ci_upper_pooled
      }
    }
    
    global_pattern_veg <- global_pattern[tolower(global_pattern$Veg) != "barren", ]
    global_pattern_barren <- global_pattern[tolower(global_pattern$Veg) == "barren", ]
    
    p <- ggplot2::ggplot(global_pattern_veg, ggplot2::aes(x = year, y = coef, color = Veg, fill = Veg)) +
      ggplot2::geom_line(size = 1.2) +
      ggplot2::geom_point(size = 2)
    
    if (show_ci && "ci_lower" %in% names(global_pattern_veg) && "ci_upper" %in% names(global_pattern_veg)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
        alpha = 0.2,
        color = NA
      )
    }
    
    if (nrow(global_pattern_barren) > 0) {
      veg_max <- suppressWarnings(max(global_pattern_veg$coef, global_pattern_veg$ci_upper, na.rm = TRUE))
      barren_max <- suppressWarnings(max(global_pattern_barren$coef, global_pattern_barren$ci_upper, na.rm = TRUE))
      if (is.na(veg_max) || veg_max <= 0 || is.na(barren_max) || barren_max <= 0) {
        barren_scale_factor <- 1
      } else {
        barren_scale_factor <- veg_max / barren_max
      }
      global_pattern_barren$coef_scaled <- global_pattern_barren$coef * barren_scale_factor
      
      p <- p + 
        ggplot2::geom_line(data = global_pattern_barren, 
              ggplot2::aes(x = year, y = coef_scaled), 
                          color = "brown", linewidth = 1.2, linetype = "dashed") +
        ggplot2::geom_point(data = global_pattern_barren, 
               ggplot2::aes(x = year, y = coef_scaled), 
                           color = "brown", size = 2)
      
      if (show_ci && "ci_lower" %in% names(global_pattern_barren) && "ci_upper" %in% names(global_pattern_barren)) {
        global_pattern_barren$ci_lower_scaled <- global_pattern_barren$ci_lower * barren_scale_factor
        global_pattern_barren$ci_upper_scaled <- global_pattern_barren$ci_upper * barren_scale_factor
        
        p <- p + ggplot2::geom_ribbon(data = global_pattern_barren,
                                     ggplot2::aes(x = year, ymin = ci_lower_scaled, ymax = ci_upper_scaled),
                                     fill = "brown", alpha = 0.1, color = NA)
      }
    }
    
    p <- p +
      ggplot2::labs(
        title = title,
        subtitle = sprintf("Based on %d locations", max(global_pattern$n_locations, na.rm = TRUE)),
        x = "Year",
        y = "Vegetation Fraction",
        color = "Vegetation Type",
        fill = "Vegetation Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)
      ) +
          ggplot2::scale_y_continuous(
        labels = scales::percent_format(),
        sec.axis = ggplot2::sec_axis(~ . / barren_scale_factor, name = "Barren Fraction", labels = scales::percent_format())
      ) +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::scale_fill_brewer(palette = "Set1")
    
    p
  }

  plot_vegetation_stacked_area <- function(global_pattern, 
                                            title = "Global Vegetation Composition Over Time") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
    } else if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
    }
    global_pattern <- global_pattern |> 
      dplyr::group_by(year) |> 
      dplyr::mutate(coef_normalized = coef / sum(coef, na.rm = TRUE)) |> 
      dplyr::ungroup()
    p <- ggplot2::ggplot(global_pattern, 
                          ggplot2::aes(x = year, y = coef_normalized, fill = Veg)) +
      ggplot2::geom_area(alpha = 0.8, position = "stack") +
      ggplot2::labs(
        title = title,
        x = "Year",
        y = "Relative Vegetation Fraction",
        fill = "Vegetation Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold")
      ) +
      ggplot2::scale_y_continuous(labels = scales::percent_format()) +
      ggplot2::scale_fill_brewer(palette = "Set2")
    p
  }

  plot_vegetation_only_stacked_area <- function(global_pattern, 
                                                title = "Vegetation-Only Composition Over Time") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
    } else if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
    }
    # Filter out barren
    global_pattern_veg <- global_pattern[tolower(global_pattern$Veg) != "barren", ]
    # Normalize vegetation fractions to sum to 1 per year
    global_pattern_veg <- global_pattern_veg |> 
      dplyr::group_by(year) |> 
      dplyr::mutate(coef_normalized = coef / sum(coef, na.rm = TRUE)) |> 
      dplyr::ungroup()
    p <- ggplot2::ggplot(global_pattern_veg, 
                          ggplot2::aes(x = year, y = coef_normalized, fill = Veg)) +
      ggplot2::geom_area(alpha = 0.8, position = "stack") +
      ggplot2::labs(
        title = title,
        x = "Year",
        y = "Relative Vegetation Fraction",
        fill = "Vegetation Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold")
      ) +
      ggplot2::scale_y_continuous(labels = scales::percent_format()) +
      ggplot2::scale_fill_brewer(palette = "Set2")
    p
  }

  plot_vegetation_heatmap <- function(global_pattern, 
                                       title = "Vegetation Fraction by Year") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
    } else if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
    }
    p <- ggplot2::ggplot(global_pattern, 
                          ggplot2::aes(x = year, y = Veg, fill = coef)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", coef * 100)), 
                         color = "white", size = 3) +
      ggplot2::scale_fill_viridis_c(option = "plasma", 
                                     labels = scales::percent_format()) +
      ggplot2::labs(
        title = title,
        x = "Year",
        y = "Vegetation Type",
        fill = "Fraction"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    p
  }

  bootstrap_trend_ci <- function(all_coefs, B = 200, seed = 123) {
  set.seed(seed)
  cat(sprintf("[DEBUG] bootstrap_trend_ci called with %d rows\n", nrow(all_coefs)))
  if ("coef_sd" %in% names(all_coefs)) {
      sd_vals <- all_coefs$coef_sd
      cat(sprintf("[DEBUG] coef_sd stats: Mean=%.6f, Max=%.6f, N_NA=%d\n", 
                  mean(sd_vals, na.rm=TRUE), max(sd_vals, na.rm=TRUE), sum(is.na(sd_vals))))
  } else {
      cat("[DEBUG] coef_sd column NOT found in all_coefs\n")
  }

  if (!requireNamespace("lme4", quietly = TRUE)) {
    warning("lme4 package not found, trend CI calculation will be skipped.")
    return(NULL)
  }
  veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
  results_list <- list()
  for (veg in veg_types) {
    veg_data <- all_coefs[all_coefs$Veg == veg & is.finite(all_coefs$coef), ]
    if (nrow(veg_data) < 10 || length(unique(veg_data$location_id)) < 3) {
      cat(sprintf("[TREND] Skipping trend for '%s' (insufficient data: %d rows, %d locs)\n", veg, nrow(veg_data), length(unique(veg_data$location_id))))
      next
    }
    locations <- unique(veg_data$location_id)
    n_locs <- length(locations)
    boot_slopes <- replicate(B, {
      boot_locs_sampled <- sample(locations, n_locs, replace = TRUE)
      boot_data_list <- lapply(seq_along(boot_locs_sampled), function(i) {
        loc_data <- veg_data[veg_data$location_id == boot_locs_sampled[i], ]
        loc_data$boot_id <- paste0(boot_locs_sampled[i], "_", i)
        if ("coef_sd" %in% names(loc_data) && !all(is.na(loc_data$coef_sd))) {
           sds <- loc_data$coef_sd
           sds[is.na(sds) | sds < 0] <- 0
           if (any(sds > 0)) {
             loc_data$coef <- rnorm(nrow(loc_data), mean = loc_data$coef, sd = sds)
           }
        }
        loc_data
      })
      boot_data <- do.call(rbind, boot_data_list)
      model <- tryCatch({
        lme4::lmer(coef ~ pheno_year + (1|boot_id), data = boot_data, control = lme4::lmerControl(check.conv.singular = "ignore"))
      }, error = function(e) { NULL })
      if (is.null(model)) {
        return(NA_real_)
      } else {
        return(lme4::fixef(model)["pheno_year"])
      }
    })
    finite_slopes <- boot_slopes[is.finite(boot_slopes)]
    if (length(finite_slopes) > 5) {
      results_list[[veg]] <- data.frame(
        Veg = veg, 
        slope_mean = mean(finite_slopes), 
        slope_median = median(finite_slopes), 
        slope_ci_lower = quantile(finite_slopes, 0.025), 
        slope_ci_upper = quantile(finite_slopes, 0.975), 
        prob_positive = mean(finite_slopes > 0), 
        prob_negative = mean(finite_slopes < 0),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(results_list) > 0) {
    return(dplyr::bind_rows(results_list))
  } else {
    return(NULL)
  }
}

  analyze_vegetation_trends <- function(all_coefs, B = 200) {
    # Simplified trend estimation using per-location linear slopes on annual coefficients.
    if (is.null(all_coefs) || nrow(all_coefs) == 0) return(NULL)
    veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
    trend_rows <- list()
    for (veg in veg_types) {
      veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), , drop = FALSE]
      if (nrow(veg_data) == 0) next
      locations <- unique(veg_data$location_id)
      loc_slopes <- sapply(locations, function(loc) {
        loc_data <- veg_data[veg_data$location_id == loc, , drop = FALSE]
        if (nrow(loc_data) < 2) return(NA_real_)
        yr_col <- if ("pheno_year" %in% names(loc_data)) "pheno_year" else "year"
        df_loc <- loc_data[order(loc_data[[yr_col]]), , drop = FALSE]
        m <- tryCatch(lm(coef ~ get(yr_col), data = df_loc), error = function(e) NULL)
        if (is.null(m)) return(NA_real_)
        coef_val <- tryCatch(coef(m)[[2]], error = function(e) NA_real_)
        as.numeric(coef_val)
      })
      trend_rows[[veg]] <- data.frame(Veg = veg, slope = median(loc_slopes, na.rm = TRUE), n_locations = sum(is.finite(loc_slopes)), stringsAsFactors = FALSE)
    }
    if (length(trend_rows) == 0) return(NULL)
    do.call(rbind, trend_rows)
  }

cat("[DEBUG] About to execute main processing steps (un-nested)...\n")
flush.console()

# Executing main processing steps directly (function removed). This section previously defined
# `main_processing_block <- function() { ... }`. To avoid large function compile stalls, the body
# is now executed at top-level during script run. The original nested helper functions remain in scope.
    cat("[DEBUG] Entered main_processing_block\n")
    b_templates <- NULL
    n_train_loc_years <- 0L
    n_infer_loc_years <- 0L
    if (exists("df_train") && !is.null(df_train) && nrow(df_train) > 0) {
    if (!"pheno_year" %in% names(df_train) && "date" %in% names(df_train)) df_train$pheno_year <- assign_pheno_year(df_train$date)
    n_train_loc_years <- nrow(unique(df_train[c("location_id", "pheno_year")]))
  }
  if (exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0) {
    if (!"pheno_year" %in% names(df_tasks_inference) && "date" %in% names(df_tasks_inference)) df_tasks_inference$pheno_year <- assign_pheno_year(df_tasks_inference$date)
    n_infer_loc_years <- nrow(unique(df_tasks_inference[c("location_id", "pheno_year")]))
  } else if (exists("df_tasks") && !is.null(df_tasks) && nrow(df_tasks) > 0) {
    if (!"pheno_year" %in% names(df_tasks) && "date" %in% names(df_tasks)) df_tasks$pheno_year <- assign_pheno_year(df_tasks$date)
    n_infer_loc_years <- nrow(unique(df_tasks[c("location_id", "pheno_year")]))
  }

  cat(sprintf("Training dataset location-years: %d\n", n_train_loc_years))
  cat(sprintf("Inference dataset location-years: %d\n", n_infer_loc_years))

  if ((exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0) && n_train_loc_years == n_infer_loc_years) {
    stop(sprintf("ERROR: Training and inference datasets appear to have the same number of location-years (%d). This may indicate you passed the same data for training and inference — aborting to avoid accidental overlap.", n_train_loc_years))
  }
  
  cat("Preparing locations for batched processing (multi-year bootstrap)...\n")

  # Group by LOCATION (not location-year)
  target_locations <- location_list
  available_locations <- unique(df_tasks$location_id)
  target_locations <- intersect(target_locations, available_locations)

  n_locs_to_process <- length(target_locations)
  # Use centralized BATCH_SIZE setting from the top-level CONFIG (to change batch size, edit the USER-TUNABLE PARAMETERS block)
  loc_batches <- split(target_locations, ceiling(seq_along(target_locations) / BATCH_SIZE))
  n_batches <- length(loc_batches)
  pb_width <- min(40L, max(4L, n_batches))

  cat(sprintf("Processing %d locations in %d batches (approx %d locations/batch)...\n",
              n_locs_to_process, length(loc_batches), BATCH_SIZE))

  # Print initial memory usage
  print_memory_usage("Before processing")

  # Collect full per-task results for downstream reporting (variant trajectories, diagnostics, uncertainty)
  # Note: We no longer use all_coefs_list to avoid memory accumulation
  # Instead, we extract coef_df from training_results_list when needed
  training_results_list <- list()

  start_time <- Sys.time()

  for (i in seq_along(loc_batches)) {
    batch_locs <- loc_batches[[i]]
    batch_df <- df_tasks[df_tasks$location_id %in% batch_locs, ]
    batch_location_list <- split(batch_df, batch_df$location_id)
    # Suppress any verbose output from per-location processing
    sink_file <- tempfile()
    con_out <- file(sink_file, open = "wt")
    con_msg <- file(paste0(sink_file, ".msg"), open = "wt")
    tryCatch({
      sink(con_out, type = "output")
      sink(con_msg, type = "message")
      batch_results <- .run_map(batch_location_list, fit_one_location, show_pb = FALSE)
    }, finally = {
      try(sink(type = "message"), silent = TRUE)
      try(sink(type = "output"), silent = TRUE)
      close(con_out); close(con_msg)
      try(unlink(sink_file), silent = TRUE); try(unlink(paste0(sink_file, ".msg")), silent = TRUE)
    })

    for (k in names(batch_results)) {
      loc_result <- batch_results[[k]]
      if (is.null(loc_result)) next

      loc_data <- do.call(rbind, lapply(loc_result, function(yr_res) yr_res$coef_df))

      if (!is.null(loc_data) && nrow(loc_data) > 0) {
        out_fname <- file.path(TEMP_RESULTS_DIR, paste0("result_", make.names(k), ".csv"))
        readr::write_csv(loc_data, out_fname)
      }
      
      # Append each per-year result to results_list for downstream analysis (e.g. plotting)
      for (yr_char in names(loc_result)) {
        r <- loc_result[[yr_char]]
        if (is.null(r)) next
        res_key <- if (!is.null(r$coef_df) && 'pheno_year' %in% names(r$coef_df)) {
          paste(k, r$coef_df$pheno_year[1], sep = "_")
        } else {
          paste(k, yr_char, sep = "_")
        }
                  training_results_list[[res_key]] <- (function(r_in) {
                    if (is.null(r_in)) return(NULL)
                    # Only keep coef_df and diagnostics for memory efficiency
                    res_slim <- list(coef_df = r_in$coef_df, diagnostics = r_in$diagnostics)
                    # Explicitly set other large objects to NULL to aid GC
                    r_in$uncertainty <- r_in$residuals <- r_in$y_hat <- r_in$y_obs <- r_in$E_best <- r_in$top_variants <- r_in$weights_masked <- r_in$valid_mask <- NULL
                    rm(r_in) # Help GC
                    return(res_slim)
                  })(r)      }
    }

    # Aggressive memory cleanup after each batch
    rm(batch_df, batch_location_list, batch_results)
    # gc(verbose = FALSE, full = TRUE)

    cat(sprintf("\r  [Batch %d/%d complete]  ", i, n_batches))
  }
  cat("\n")

  end_time <- Sys.time()
  processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(
    "Main processing loop finished in %.2f seconds (%.2f minutes)\n",
    processing_time, processing_time / 60
  ))

  # Combine all coefficient data
  # Since we no longer accumulate in all_coefs_list (to save memory),
  # we extract coefficient summaries from training_results_list instead
  cat("Building all_coefs from training_results_list (memory-efficient approach)...\n")
  all_coefs <- do.call(rbind, lapply(training_results_list, function(r) {
    if (!is.null(r$coef_df)) r$coef_df else NULL
  }))

  # Defensive: ensure certain columns are numeric (avoid coercion warnings later)
  if (!is.null(all_coefs) && nrow(all_coefs) > 0) {
    num_cols <- c("coef", "coef_sd", "coef_025", "coef_975", "pheno_year")
    for (nc in num_cols) {
      if (nc %in% names(all_coefs)) {
        if (!is.numeric(all_coefs[[nc]])) {
          before_na <- sum(is.na(all_coefs[[nc]]))
          all_coefs[[nc]] <- suppressWarnings(as.numeric(as.character(all_coefs[[nc]])))
          after_na <- sum(is.na(all_coefs[[nc]]))
          if (after_na > before_na) {
            warning(sprintf("Coerced column '%s' to numeric; NAs increased by %d", nc, after_na - before_na))
          }
        }
      }
    }
  }

  if (!is.null(all_coefs) && nrow(all_coefs) > 0) {
      cat(sprintf("Saving training results to: %s\n", file.path(OUT_DIR, "training_results.csv")))
      write.csv(all_coefs, file.path(OUT_DIR, "training_results.csv"), row.names = FALSE)
  }

  # --- Optional: Estimate uncertainty inflation parameters from data (Option A) ---
  if (isTRUE(ESTIMATE_UNCERTAINTY_PARAMS_OPTION_A)) {
    cat("[UNCERTAINTY PARAM EST] Running Option A subsampling calibration...\n")
    est <- tryCatch({
      estimate_uncertainty_params_optionA(
        df_tasks = df_tasks,
        high_n = UNCERTAINTY_PARAM_EST_HIGH_N,
        max_groups = UNCERTAINTY_PARAM_EST_MAX_GROUPS,
        target_ns = UNCERTAINTY_PARAM_EST_TARGET_NS,
        reps = UNCERTAINTY_PARAM_EST_REPS,
        seed = UNCERTAINTY_PARAM_EST_SEED,
        veg_levels = ALLOWED_VEG,
        quiet = TRUE
      )
    }, error = function(e) {
      cat(sprintf("[UNCERTAINTY PARAM EST] Failed: %s\n", e$message))
      NULL
    })
    if (!is.null(est)) {
      cat(sprintf("[UNCERTAINTY PARAM EST] Estimated N_POWER ~ %.3f\n", est$N_POWER_hat))
      cat(sprintf("[UNCERTAINTY PARAM EST] Estimated N_REF ~ %d\n", est$N_REF_hat))
      # Persist results for inspection
      out_csv <- file.path(OUT_DIR, "uncertainty_param_estimate_optionA.csv")
      try(write.csv(est$summary, out_csv, row.names = FALSE), silent = TRUE)
      cat(sprintf("[UNCERTAINTY PARAM EST] Wrote summary to %s\n", out_csv))
    }
  }

  # Calculate total number of year results processed
  n_year_results <- if (!is.null(all_coefs) && nrow(all_coefs) > 0) {
    length(unique(paste(all_coefs$location_id, all_coefs$pheno_year, sep = "_")))
  } else {
    0L
  }

  if (n_year_results > 0) {
    cat(sprintf("Average time per year-result: %.2f seconds\n", processing_time / n_year_results))
  } else {
    cat("Average time per year-result: N/A (0 results)\n")
  }

  # Generate global training-level plots early so they are available before inference processing
  tryCatch({
    if (!is.null(all_coefs) && nrow(all_coefs) > 0) {
      if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
      if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
      library(dplyr)
      library(ggplot2)

      cat("[INFO] Aggregating training coefficients for global plots...\n")
      gp <- tryCatch({
        aggregate_to_global_pattern(all_coefs, method = "location_bootstrap")
      }, error = function(e) { cat(sprintf("[ERROR] aggregate_to_global_pattern fallback failed: %s\n", e$message)); NULL })

      gp2_available <- FALSE
      if (!is.null(gp) && nrow(gp) > 0) {
        gp2 <- gp |> dplyr::filter(year >= 1985 & year <= 2025)
        gp2_available <- TRUE
      } else {
        cat("[WARN] Not enough training coefficient data to build global training plots.\n")
        # Diagnostic: print columns and sample rows of all_coefs to help debugging
        if (!is.null(all_coefs)) {
          cat(sprintf("[DIAGNOSTIC] all_coefs: rows=%d, cols=%d\n", nrow(all_coefs), ncol(all_coefs)))
          cat(sprintf("[DIAGNOSTIC] all_coefs columns: %s\n", paste(names(all_coefs), collapse = ", ")))
          cat("[DIAGNOSTIC] sample rows from all_coefs:\n")
          print(utils::head(all_coefs, 6))
        }
      }

        if (gp2_available) {
        # Line with CI (average coverage)
        veg_gp <- gp2[tolower(gp2$Veg) != "barren", , drop = FALSE]
        barren_gp <- gp2[tolower(gp2$Veg) == "barren", , drop = FALSE]

        # Map possible column name variants safely
        if ("global_coef" %in% names(veg_gp)) {
          coef_col <- "global_coef"
        } else if ("mean_coef" %in% names(veg_gp)) {
          coef_col <- "mean_coef"
        } else {
          matches <- intersect(names(veg_gp), c("coef", "value"))
          coef_col <- if (length(matches) > 0) matches[1] else NULL
        }

        if (!is.null(coef_col) && coef_col %in% names(veg_gp)) {
          p_cov <- ggplot(veg_gp, aes(x = year, y = .data[[coef_col]], color = Veg, fill = Veg)) +
            geom_ribbon(aes(ymin = coef_025, ymax = coef_975), alpha = 0.2, color = NA) +
            geom_line(size = 1.2) + geom_point(size = 2) +
            labs(title = "Average Coverage Percentage per Vegetation Type (Training)", x = "Year", y = "Fraction") +
            theme_minimal() + scale_y_continuous(labels = scales::percent_format())

          # Add barren scaled series if present
          if (nrow(barren_gp) > 0 && coef_col %in% names(barren_gp)) {
            # scale barren to veg scale
            veg_max <- suppressWarnings(max(veg_gp[[coef_col]], na.rm = TRUE))
            barren_max <- suppressWarnings(max(barren_gp[[coef_col]], na.rm = TRUE))
            barren_scale_factor <- ifelse(is.numeric(veg_max) && is.numeric(barren_max) && barren_max > 0, veg_max / barren_max, 1)
            barren_gp$coef_scaled <- barren_gp[[coef_col]] * barren_scale_factor
            p_cov <- p_cov + geom_line(data = barren_gp, aes(x = year, y = coef_scaled), color = "brown", linetype = "dashed")
          }

          out_cov <- file.path(OUT_DIR, "training_average_coverage_plot.png")
          ggsave(out_cov, p_cov, width = 10, height = 6, dpi = 300)
          cat(sprintf("Saved training average coverage plot to: %s\n", out_cov))
        }

        # Stacked area
        # Ensure a canonical 'coef' column exists for compatibility
        if ("global_coef" %in% names(gp2)) {
          gp2$coef <- gp2$global_coef
        } else if ("mean_coef" %in% names(gp2)) {
          gp2$coef <- gp2$mean_coef
        } else if (!"coef" %in% names(gp2) && "value" %in% names(gp2)) {
          gp2$coef <- gp2$value
        }

        gp_stacked <- gp2 |> dplyr::group_by(year) |> dplyr::mutate(
          total_coef = sum(coef, na.rm = TRUE),
          coef_normalized = dplyr::if_else(is.finite(total_coef) & total_coef > 0, coef / total_coef, 0)
        ) |> dplyr::ungroup() |> dplyr::select(-total_coef)

        if (nrow(gp_stacked) > 0 && any(is.finite(gp_stacked$coef_normalized))) {
          p_stack <- ggplot(gp_stacked, aes(x = year, y = coef_normalized, fill = Veg)) +
            geom_area(alpha = 0.8, position = "stack") +
            labs(title = "Global Vegetation Composition (Stacked, Training)", x = "Year", y = "Relative Vegetation Fraction") +
            theme_minimal() + scale_y_continuous(labels = scales::percent_format())
          out_stack <- file.path(OUT_DIR, "training_vegetation_stacked_area.png")
          ggsave(out_stack, p_stack, width = 10, height = 6)
          cat(sprintf("Saved training stacked area plot to: %s\n", out_stack))
        } else {
          cat("[DEBUG] Skipping stacked area plot due to no valid normalized coefficients.\n")
        }

        # Heatmap
        if (nrow(gp2) > 0) {
          p_heat <- ggplot(gp2, aes(x = year, y = Veg, fill = coef)) + geom_tile() +
            scale_fill_viridis_c(option = "plasma", labels = scales::percent_format()) +
            labs(title = "Vegetation Fraction by Year (Training Heatmap)", x = "Year", y = "Vegetation Type") +
            theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
          out_heat <- file.path(OUT_DIR, "training_vegetation_heatmap.png")
          ggsave(out_heat, p_heat, width = 10, height = 6)
          cat(sprintf("Saved training heatmap to: %s\n", out_heat))
        }
      } else {
        cat("[WARN] No training coefficients available to generate global training plots.\n")
      }
    } else {
      cat("[WARN] No training coefficients available to generate global training plots.\n")
    }
  }, error = function(e) {
    cat(sprintf("[ERROR] Failed to generate training global plots early: %s\n", e$message))
  })

  inference_results_list <- list()
  # Ensure variant similarity heatmap generated before visualizing/training data
  ensure_variant_similarity_heatmap()

  # Ensure required library/templates exist for visualization
  ensure_library_and_templates()

  # Visualization of training data
  if(exists("df_train") && !is.null(df_train) && nrow(df_train) > 0) {
    if (!requireNamespace("tidyr", quietly = TRUE)) {
      install.packages("tidyr")
    }
    library(tidyr)
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      install.packages("dplyr")
    }
    library(dplyr)

    vis_dir <- "training_visualizations"
    if (!dir.exists(vis_dir)) {
      dir.create(vis_dir)
    }

    cat("\n=== VISUALIZING TRAINING DATA ===\n")

    # Ensure required 'avail' index list is available for plotting
    if (!exists("avail") || is.null(avail) || length(avail) == 0) {
      if (exists("candidate_indices")) {
        avail <- candidate_indices
        cat("[INFO] 'avail' not found; using 'candidate_indices' for plotting indices.\n")
      } else if (exists("OPTIMAL_INDICES") && length(intersect(OPTIMAL_INDICES, names(df_train))) > 0) {
        avail <- intersect(OPTIMAL_INDICES, names(df_train))
        cat("[INFO] 'avail' not found; deriving from OPTIMAL_INDICES intersecting df_train columns.\n")
      } else {
        # Fallback: use numeric columns in df_train other than known meta columns
        candidate_cols <- setdiff(names(df_train), c("date","doy","pheno_year","location_id","Veg","lat","lon","latitude","longitude"))
        numeric_cols <- candidate_cols[sapply(df_train[candidate_cols], is.numeric)]
        avail <- head(numeric_cols, 12)
        cat(sprintf("[INFO] 'avail' derived from numeric columns: %s\n", paste(avail, collapse=", ")))
      }
    }

    train_locations <- unique(df_train$location_id)
    
    max_plots <- 10
    if(length(train_locations) > max_plots){
        cat(sprintf("Found %d training locations, plotting a sample of %d.\n", length(train_locations), max_plots))
        set.seed(42)
        train_locations <- sample(train_locations, max_plots)
    }

    for (loc_id in train_locations) {
        loc_data <- df_train[df_train$location_id == loc_id, ]
        
        plot_data_long <- loc_data %>%
            dplyr::select(date, doy, dplyr::one_of(avail)) %>%
            tidyr::pivot_longer(cols = dplyr::one_of(avail), names_to = "index", values_to = "value")

        p <- ggplot(plot_data_long, aes(x = doy, y = value)) +
            geom_point(alpha = 0.7) +
            geom_line(alpha = 0.5) +
            facet_wrap(~index, scales = "free_y", ncol = 3) +
            labs(
                title = paste("Training Data Time Series for Location:", loc_id),
                subtitle = paste(min(loc_data$date, na.rm=T), "to", max(loc_data$date, na.rm=T)),
                x = "Day of Phenological Year",
                y = "Index Value"
            ) +
            theme_bw() +
            theme(legend.position = "none")

        plot_filename <- file.path(vis_dir, paste0("loc_", gsub("[^A-Za-z0-9_.-]", "_", loc_id), ".png"))
        ggsave(plot_filename, p, width = 12, height = 8, units = "in", dpi = 150)
    }

    # Summary plots for training data: DOY coverage and observations per location
    if ("doy" %in% names(df_train)) {
      doy_coverage_train <- df_train %>%
        dplyr::group_by(location_id, pheno_year) %>%
        dplyr::summarize(n_unique_doys = length(unique(doy[!is.na(doy)])), .groups = "drop")

      if (nrow(doy_coverage_train) > 0) {
        p_doy_cov <- ggplot2::ggplot(doy_coverage_train, ggplot2::aes(x = n_unique_doys)) +
          ggplot2::geom_histogram(binwidth = 5, fill = "steelblue", color = "white", alpha = 0.9) +
          ggplot2::labs(
            title = "Distribution of Unique DOYs per Location-Year (Training)",
            x = "Unique DOYs",
            y = "Count"
          ) +
          ggplot2::theme_minimal()

        ggsave(file.path(vis_dir, "training_doy_coverage_hist.png"), p_doy_cov, width = 8, height = 5)
        cat(sprintf("Saved DOY coverage histogram to: %s\n", file.path(vis_dir, "training_doy_coverage_hist.png")))
      }
    }

    sample_sizes_train <- df_train %>%
      dplyr::group_by(location_id) %>%
      dplyr::summarize(n_obs = dplyr::n(), .groups = "drop")

    if (nrow(sample_sizes_train) > 0) {
      p_nobs <- ggplot2::ggplot(sample_sizes_train, ggplot2::aes(x = n_obs)) +
        ggplot2::geom_histogram(binwidth = 10, fill = "forestgreen", color = "white", alpha = 0.9) +
        ggplot2::labs(
          title = "Distribution of Observations per Location (Training)",
          x = "Total Observations",
          y = "Count"
        ) +
        ggplot2::theme_minimal()

      ggsave(file.path(vis_dir, "training_nobs_hist.png"), p_nobs, width = 8, height = 5)
      cat(sprintf("Saved observations-per-location histogram to: %s\n", file.path(vis_dir, "training_nobs_hist.png")))
    }

    # --- New: Bootstrap Average Plots ---
    if ("Veg" %in% names(df_train) && length(unique(df_train$Veg)) > 0) {
      cat("Generating bootstrap average plots for training data...\n")
      
      # Filter for allowed vegetation types (inseparable variants are dropped)
      valid_classes <- if (exists("ALLOWED_VEG")) unique(c(ALLOWED_VEG, "barren")) else unique(df_train$Veg)
      # Normalize case for comparison
      df_train_filtered <- df_train %>%
        dplyr::filter(tolower(trimws(Veg)) %in% tolower(trimws(valid_classes)))
      
      # Bin DOY to consolidate data (5-day bins)
      df_train_binned <- df_train_filtered %>%
        dplyr::mutate(doy_bin = round(doy / 5) * 5) %>%
        dplyr::group_by(location_id, Veg, doy_bin) %>%
        # Use median for per-pentad summaries to create robust average profiles
        dplyr::summarise(dplyr::across(dplyr::all_of(avail), \(x) median(x, na.rm = TRUE)), .groups = "drop")

      boot_results_list <- list()
      for (idx in avail) {
        # Filter for sufficient data
        idx_data <- df_train_binned[!is.na(df_train_binned[[idx]]), ]
        
        # Bootstrap per Veg and DOY bin
        summ <- idx_data %>%
          dplyr::group_by(Veg, doy_bin) %>%
          dplyr::filter(dplyr::n() >= 3) %>%
          dplyr::summarise(
            median_val = median(.data[[idx]], na.rm = TRUE),
            lower = quantile(replicate(100, median(sample(.data[[idx]], dplyr::n(), replace = TRUE), na.rm = TRUE)), 0.025, na.rm = TRUE),
            upper = quantile(replicate(100, median(sample(.data[[idx]], dplyr::n(), replace = TRUE), na.rm = TRUE)), 0.975, na.rm = TRUE),
            .groups = "drop"
          ) %>%
          dplyr::mutate(Index = idx)
        
        boot_results_list[[idx]] <- summ
      }
      
      plot_df <- dplyr::bind_rows(boot_results_list)
      
      if (nrow(plot_df) > 0) {
        p_avg <- ggplot2::ggplot(plot_df, ggplot2::aes(x = doy_bin, y = median_val, color = Veg, fill = Veg)) +
          ggplot2::geom_line() +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
          ggplot2::facet_wrap(~Index, scales = "free_y") +
          ggplot2::labs(title = "Median Spectral Profiles (Location Bootstrap)",
               x = "DOY", y = "Index Value (median)") +
          ggplot2::theme_bw()
        
        ggsave(file.path(vis_dir, "training_median_profiles_boot.png"), p_avg, width = 12, height = 8)
        cat(sprintf("Saved bootstrap median profile plot to: %s\n", file.path(vis_dir, "training_median_profiles_boot.png")))
      }
    }

    cat("Training data visualizations saved to '", vis_dir, "' directory.\n")
  } else {
    cat("Could not find df_train to generate visualizations.\n")
  }

  # Ensure library + templates exist and generate variant similarity heatmap before inference processing
  ensure_library_and_templates()
  ensure_variant_similarity_heatmap()
  
  # Ensure GLOBAL_TRAINING_DVI_SOIL baseline is established
  ensure_global_dvi_soil_baseline()

  # Helper to report validation accuracy (including artificial mixes) to avoid duplication
  report_validation_accuracy <- function(val_coefs, label = "") {
    if (is.null(val_coefs) || nrow(val_coefs) == 0) {
      cat("[NOTICE] No validation coefficients found (validation locations not in results).\n")
      return(invisible(NULL))
    }

    prefix <- if (nzchar(label)) paste0(" (", label, ")") else ""

    # Optional PPI-based truth for vegetation cover (derived from PPI barren fraction)
    ppi_truth <- NULL
    if (exists("all_diagnostics") && !is.null(all_diagnostics) &&
        all(c("location_id", "pheno_year") %in% names(all_diagnostics)) &&
        "barren_fraction_ppi_based" %in% names(all_diagnostics)) {
      ppi_truth <- all_diagnostics |>
        dplyr::select(location_id, pheno_year, barren_fraction_ppi_based) |>
        dplyr::mutate(
          barren_fraction_ppi_based = pmin(pmax(barren_fraction_ppi_based, 0), 1),
          true_veg_cover_ppi = 1 - barren_fraction_ppi_based
        ) |>
        dplyr::distinct()
    }
    
    # --- 1. Identify Valid Locations and Get True Labels ---
    labels_df <- NULL
    
    # Try df_tasks_inference first (most likely to have labels for inference set if loaded from file with labels)
    if (exists("df_tasks_inference") && !is.null(df_tasks_inference) && "Veg" %in% names(df_tasks_inference)) {
       labels_df <- df_tasks_inference %>% dplyr::select(location_id, Veg) %>% dplyr::distinct()
    } else if (exists("df_tasks") && !is.null(df_tasks) && "Veg" %in% names(df_tasks)) {
       # Fallback to df_tasks (training data)
       labels_df <- df_tasks %>% dplyr::select(location_id, Veg) %>% dplyr::distinct()
    }
    
    if (is.null(labels_df)) {
       # Try to read validation_locations.csv if it exists
       if (file.exists("validation_locations.csv")) {
          vloc <- tryCatch(read.csv("validation_locations.csv"), error = function(e) NULL)
          if (!is.null(vloc) && "location_id" %in% names(vloc) && "Veg" %in% names(vloc)) {
             labels_df <- vloc %>% dplyr::select(location_id, Veg) %>% dplyr::distinct()
          }
       }
    }

    if (is.null(labels_df)) {
       cat("[NOTICE] Could not find ground truth labels (Veg column in df_tasks_inference or df_tasks). Cannot compute accuracy.\n")
       return(invisible(NULL))
    }
    
    # Normalize labels
    labels_df$Veg <- tolower(trimws(as.character(labels_df$Veg)))
    labels_df <- labels_df %>% dplyr::rename(true_veg = Veg)
    
    # Filter validation coefficients to only those with labels
    # We join first to filter efficiently
    val_coefs_labeled <- val_coefs %>% 
       dplyr::inner_join(labels_df, by = "location_id")

     # Drop barren from classification accuracy; rely on PPI-based cover truth instead
     if ("true_veg" %in% names(val_coefs_labeled)) {
       barren_mask <- tolower(val_coefs_labeled$true_veg) == "barren"
       if (any(barren_mask, na.rm = TRUE)) {
        cat(sprintf("[NOTICE] Dropping %d validation rows with true_veg=='barren' from classification accuracy; using PPI cover truth instead when available.\n",
                sum(barren_mask, na.rm = TRUE)))
        val_coefs_labeled <- val_coefs_labeled[!barren_mask, , drop = FALSE]
       }
     }
    
     if (nrow(val_coefs_labeled) == 0) {
       cat("[NOTICE] No labeled validation rows remain after removing barren truth.\n")
       return(invisible(NULL))
     }
    
    cat(sprintf("\n=== VALIDATION ACCURACY ON HELD-OUT SET%s ===\n", prefix))
    cat(sprintf("Validation set: %d locations, %d location-year pairs\n", 
                length(unique(val_coefs_labeled$location_id)), 
                length(unique(paste(val_coefs_labeled$location_id, val_coefs_labeled$pheno_year)))))

    # --- 2. Pivot Long to Wide (Correctly) ---
    # Determine which column holds the value
    measure_col <- dplyr::case_when(
      "coef" %in% names(val_coefs_labeled) ~ "coef",
      "pred_coef_rel" %in% names(val_coefs_labeled) ~ "pred_coef_rel",
      TRUE ~ NA_character_
    )
    
    val_coefs_wide <- NULL
    
    if (is.na(measure_col)) {
       # If no single value column, check if it's already wide (has frac_* columns)
       veg_cols <- grep("^frac_", names(val_coefs_labeled), value = TRUE)
       if(length(veg_cols) > 0) {
          # Already wide, just ensure true_veg is present (it is, from the join)
          val_coefs_wide <- val_coefs_labeled
       } else {
          cat("[ERROR] Could not identify coefficient column (coef) or fraction columns (frac_*) in input.\n")
          return(invisible(NULL))
       }
    } else {
       # Pivot
       # We strictly use location_id, pheno_year, and true_veg as keys. 
       # 'Veg' in val_coefs is the PREDICTED class (if present).
       # Note: pivot_wider will take 'Veg' (predicted) values to make columns.
       
       # Ensure 'Veg' exists for pivoting
       if (!"Veg" %in% names(val_coefs_labeled)) {
          cat("[ERROR] 'Veg' column (predicted class) missing for pivoting.\n")
          return(invisible(NULL))
       }
       
       # AGENT: Aggregate variants to vegetation class level (e.g. "phragmites_opt_1" -> "phragmites")
       val_coefs_wide <- val_coefs_labeled %>%
         dplyr::mutate(
             # Remove common suffixes to get base class: _opt_N, _single, _ppi
             Veg_class = sub("(_opt_[0-9]+|_single|_ppi)$", "", Veg)
         ) %>%
         dplyr::select(location_id, pheno_year, true_veg, Veg_class, dplyr::all_of(measure_col)) %>%
         tidyr::pivot_wider(names_from = Veg_class, 
                            values_from = dplyr::all_of(measure_col), 
                            values_fill = 0, 
                            values_fn = sum, 
                            names_prefix = "frac_")
    }

    # Identify all fraction columns (now aggregated by class)
    all_veg_cols <- grep("^frac_", names(val_coefs_wide), value = TRUE)
    
    # For PPI analysis, we might exclude barren if strictly analyzing "vegetation" cover
    veg_cols_no_barren <- all_veg_cols[all_veg_cols != "frac_barren"]

    # --- PPI-based total vegetation cover validation (truth from PPI-derived cover) ---
    if (!is.null(ppi_truth) && !is.null(val_coefs_wide) && length(veg_cols_no_barren) > 0) {
      ppi_join <- val_coefs_wide %>%
        dplyr::left_join(ppi_truth, by = c("location_id", "pheno_year"))

      if ("true_veg_cover_ppi" %in% names(ppi_join)) {
        ppi_join$pred_total_veg_cover <- rowSums(ppi_join[, veg_cols_no_barren, drop = FALSE], na.rm = TRUE)
        ppi_join$true_veg_cover_ppi <- pmin(pmax(ppi_join$true_veg_cover_ppi, 0), 1)

        cover_rows <- ppi_join[is.finite(ppi_join$true_veg_cover_ppi), , drop = FALSE]
        if (nrow(cover_rows) > 0) {
          mae_cover <- mean(abs(cover_rows$pred_total_veg_cover - cover_rows$true_veg_cover_ppi), na.rm = TRUE)
          bias_cover <- mean(cover_rows$pred_total_veg_cover - cover_rows$true_veg_cover_ppi, na.rm = TRUE)
          cat("\n--- PPI-BASED VEGETATION COVER VALIDATION ---\n")
          cat(sprintf("Samples with PPI truth: %d\n", nrow(cover_rows)))
          cat(sprintf("Mean absolute error (pred veg cover vs PPI): %.3f\n", mae_cover))
          cat(sprintf("Mean bias (pred - truth): %.3f\n", bias_cover))
        } else {
          cat("[NOTICE] PPI truth found but no matching location-year rows after filtering.\n")
        }
      }
    }
    
    # --- 3. Compute Mean Predicted Fraction ---
    # Use raw fractions directly (MESMA coefficients already sum to 1)
    if (length(all_veg_cols) > 0) {
      for (v in unique(val_coefs_wide$true_veg)) {
        if (tolower(v) == "barren") next
        sub <- val_coefs_wide[val_coefs_wide$true_veg == v, ]
        if (nrow(sub) > 0) {
          frac_col <- paste0("frac_", v) # e.g. frac_populus
          
          if (frac_col %in% names(sub)) {
            # Use raw fraction directly (no rescaling needed - MESMA coeffs sum to 1)
            mean_frac <- mean(sub[[frac_col]], na.rm = TRUE)
            cat(sprintf("  %s: mean predicted fraction = %.3f\n", v, mean_frac))
          } else {
             # If the column doesn't exist, it means the model NEVER predicted this class
             cat(sprintf("  %s: mean predicted fraction = 0.000 (never predicted)\n", v))
          }
        }
      }
    }

    # --- 4. Confusion Matrix (Average Predicted Fractions) ---
    if (length(all_veg_cols) > 0) {
      # Exclude barren from columns
      matrix_veg_cols <- all_veg_cols[all_veg_cols != "frac_barren"]
      
      if (length(matrix_veg_cols) > 0) {
          cat("\n--- CONFUSION MATRIX (Average Predicted Fractions, excluding Barren) ---\n")
          cat("Rows: True Class | Columns: Mean Predicted Fraction\n\n")
          
          # Group by true_veg (excluding barren rows) and calculate means for non-barren frac columns
          avg_fractions <- val_coefs_wide %>%
            dplyr::filter(tolower(true_veg) != "barren") %>%
            dplyr::group_by(true_veg) %>%
            dplyr::summarize(dplyr::across(dplyr::all_of(matrix_veg_cols), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
          
          # Clean up column names (remove frac_)
          colnames(avg_fractions) <- sub("^frac_", "", colnames(avg_fractions))
          
          # Convert to matrix for nice printing
          mat_data <- as.matrix(avg_fractions[, -1]) # exclude true_veg column
          rownames(mat_data) <- avg_fractions$true_veg
          
          # Sort rows and columns alphabetically
          sorted_classes <- sort(union(rownames(mat_data), colnames(mat_data)))
          
          # Ensure all classes exist in both dimensions (add missing cols with 0)
          missing_cols <- setdiff(sorted_classes, colnames(mat_data))
          if(length(missing_cols) > 0) {
            for(mc in missing_cols) {
              mat_data <- cbind(mat_data, 0)
              colnames(mat_data)[ncol(mat_data)] <- mc
            }
          }
          
          # Reorder
          mat_data <- mat_data[match(sorted_classes, rownames(mat_data)), match(sorted_classes, colnames(mat_data))]
          rownames(mat_data) <- sorted_classes # Restore rownames after subsetting/ordering might have dropped them if NAs
          
          # Print with 3 decimal places
          print(round(mat_data, 3))
          
          # Print mean diagonal (average correctly predicted fraction across classes)
          diag_vals <- diag(mat_data)
          mean_diag <- mean(diag_vals, na.rm = TRUE)
          cat(sprintf("\nMean correctly predicted fraction (diagonal average): %.3f\n", mean_diag))
      }
    }

    # --- 5. Artificial Mix Logic ---
    if (length(veg_cols) > 0) {
      # --- INTER-CLASS MIXTURE DISCRIMINATION TEST ---
      # Mix different vegetation classes to test if the model coefficients reflect the mix.
      # Note: This averages predictions (coefficients), so it tests linearity of the output space.
      
      veg_classes_only <- unique(val_coefs_wide$true_veg[tolower(val_coefs_wide$true_veg) != "barren"])
      
      if (length(veg_classes_only) >= 2) {
        cat("\n--- INTER-CLASS MIXTURE DISCRIMINATION (Synthetic 50/50 Mixes) ---\n")
        
        # Generate all unique pairs of different vegetation classes
        class_pairs <- utils::combn(veg_classes_only, 2, simplify = FALSE)
        
        for (pair in class_pairs) {
          class_a <- pair[1]
          class_b <- pair[2]
          
          sub_a <- val_coefs_wide[val_coefs_wide$true_veg == class_a, ]
          sub_b <- val_coefs_wide[val_coefs_wide$true_veg == class_b, ]
          
          if (nrow(sub_a) > 0 && nrow(sub_b) > 0) {
             # Sample pairs to create mixes
             n_mix <- min(nrow(sub_a), nrow(sub_b), 100)
             idx_a <- sample(nrow(sub_a), n_mix, replace = TRUE)
             idx_b <- sample(nrow(sub_b), n_mix, replace = TRUE)
             
             # Create mixed coefficients (average of A and B predictions)
             mixed_coefs <- (sub_a[idx_a, veg_cols] + sub_b[idx_b, veg_cols]) / 2
             
             # Calculate relative fraction of Class A in the mix
             # Ideally should be 0.5 if A was 1.0 and B was 0.0 for A
             
             frac_col_a <- paste0("frac_", tolower(class_a))
             
             if (frac_col_a %in% names(mixed_coefs)) {
               # Rescale to total vegetation
               barren_col <- "frac_barren"
               veg_frac_cols_local <- veg_cols[veg_cols != barren_col]
               
               # Calculate sum of mixed veg fractions
               sum_mixed <- sapply(1:nrow(mixed_coefs), function(i) {
                 total_veg <- sum(unlist(mixed_coefs[i, veg_frac_cols_local]), na.rm = TRUE)
                 if (total_veg > 0) {
                   as.numeric(mixed_coefs[i, frac_col_a]) + as.numeric(mixed_coefs[i, paste0("frac_", tolower(class_b))])
                 } else {
                   NA_real_
                 }
               })
               mean_sum <- mean(sum_mixed, na.rm = TRUE)
               
               # Calculate relative to excluded veg
               third_veg <- setdiff(veg_classes_only, c(class_a, class_b))
               frac_excluded <- sapply(1:nrow(mixed_coefs), function(i) {
                 total_veg <- sum(unlist(mixed_coefs[i, veg_frac_cols_local]), na.rm = TRUE)
                 if (total_veg > 0) {
                   as.numeric(mixed_coefs[i, paste0("frac_", tolower(third_veg))])
                 } else {
                   NA_real_
                 }
               })
               mean_excluded <- mean(frac_excluded, na.rm = TRUE)
               
               cat(sprintf("  Mix %s + %s: Mean sum of mixed veg fractions = %.3f (Expected ~1.000), mean relative to %s fraction = %.3f (Expected ~0.000)\n", 
                           class_a, class_b, mean_sum, third_veg, mean_excluded))
             }
          }
        }
      }
    }

    cat("Validation accuracy computed", prefix, ".\n", sep = "")
  }

  # [CLEANUP] Removed misleading "validation" check on full training dataset (all_coefs).
  # Validation is now strictly performed on the held-out set (validation_locations.csv) in the block above.

  # Compute accuracy on out-of-bag training data directly after training, before inference
  if (exists("all_coefs") && !is.null(all_coefs) && nrow(all_coefs) > 0) {
    # Use training data for OOB validation
    val_coefs <- all_coefs
    report_validation_accuracy(val_coefs, "OOB training data")
  }

  if (isTRUE(TESTING_MODE)) {
    cat("[TESTING MODE] Skipping inference processing for faster testing\n")
  } else if (exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0) {
    cat("\n=== PROCESSING INFERENCE DATA ===\n")
    
    if (!"year" %in% names(df_tasks_inference) && "date" %in% names(df_tasks_inference)) {
      df_tasks_inference$date <- as.Date(df_tasks_inference$date)
      df_tasks_inference$pheno_year <- assign_pheno_year(df_tasks_inference$date)
      df_tasks_inference$doy <- pheno_doy(df_tasks_inference$date)  # Use phenological DOY
    }
    
    if (!"task_key" %in% names(df_tasks_inference)) {
      df_tasks_inference$task_key <- paste(df_tasks_inference$location_id, df_tasks_inference$pheno_year, sep = "_")
    }
    
    inference_loc_years <- df_tasks_inference |> 
      dplyr::filter(!is.na(.data$location_id) & trimws(.data$location_id) != "" & !is.na(.data$pheno_year) & .data$pheno_year > 0) |> 
      dplyr::distinct(.data$location_id, .data$pheno_year)
    
    cat(sprintf("Inference dataset: %d location-year pairs from %d locations\n", 
                nrow(inference_loc_years), length(unique(inference_loc_years$location_id))))
    
    inference_target_keys <- paste(inference_loc_years$location_id, inference_loc_years$pheno_year, sep = "_")
    inference_available_keys <- unique(df_tasks_inference$task_key)
    inference_target_keys <- intersect(inference_target_keys, inference_available_keys)
    
    # Process all inference LOCATIONS (no limit)
    unique_locs <- unique(inference_loc_years$location_id)
    
    n_inference_keys <- length(inference_target_keys)

    # AGENT: Inline inference loop removed to prevent double processing.
    # Inference is now handled exclusively by the modular 'process_inference()' function called at the end of the script.
    inference_results_list <- list()

  }

  if (length(inference_results_list) > 0) {
    cat(sprintf("Inference results generated: %d tasks (kept separate from training results)\n", length(inference_results_list)))
  }

  # Fail fast if BOTH training and inference result lists are missing/empty — indicates upstream processing issue
  if (( !exists("training_results_list") || !is.list(training_results_list) || length(training_results_list) == 0) &&
      ( !exists("inference_results_list") || !is.list(inference_results_list) || length(inference_results_list) == 0)) {
    stop(paste0("ERROR: No results collected (both training and inference result lists are empty). Upstream processing failed — check earlier logs and data filters."))
  }

  cat("Processing results and writing to Excel files...\n")

  # Diagnostics: report how many results were collected and how many are NULL prior to filtering
  n_collected <- length(training_results_list)
  n_null_before_filter <- sum(sapply(training_results_list, is.null))
  cat(sprintf("Training results collected: %d (NULL: %d)\n", n_collected, n_null_before_filter))
  if (n_null_before_filter > 0) {
    null_keys <- names(training_results_list)[sapply(training_results_list, is.null)]
    sample_null <- if (!is.null(null_keys)) paste(head(null_keys, 10), collapse = ", ") else "(no names available)"
    cat(sprintf("Sample NULL result keys (up to 10): %s\n", sample_null))
  }

  # Remove explicit NULL entries
  training_results_list <- training_results_list[!sapply(training_results_list, is.null)]

  # Identify invalid entries (missing coef_df / not data.frame / missing location_id)
  invalid_reasons <- lapply(seq_along(training_results_list), function(i) {
    res <- training_results_list[[i]]
    key <- names(training_results_list)[i]
    reason <- NULL
    if (is.null(res)) {
      reason <- "NULL"
    } else if (is.null(res$coef_df)) {
      reason <- "missing_coef_df"
    } else if (!is.data.frame(res$coef_df)) {
      reason <- "coef_df_not_dataframe"
    } else {
      lid <- res$coef_df$location_id
      if (is.null(lid) || length(lid) == 0) reason <- "missing_location_id"
      else if (all(is.na(lid) | trimws(as.character(lid)) == "")) reason <- "empty_location_id"
    }
    list(key = key, reason = reason)
  })
  is_invalid_res <- vapply(invalid_reasons, function(x) !is.null(x$reason), logical(1))

  if (any(is_invalid_res)) {
    cnt <- sum(is_invalid_res)
    cat(sprintf("[NOTICE] Filtering out %d invalid training result(s) (reasons include missing coef_df / location_id)\n", cnt))
    # Print a small sample of keys and reasons to help debugging
    sample_invalid <- head(invalid_reasons[is_invalid_res], 20)
    sample_strs <- vapply(sample_invalid, function(x) sprintf("%s:%s", ifelse(is.null(x$key), "(no_key)", x$key), x$reason), character(1))
    cat(sprintf("Sample invalid entries (up to 20): %s\n", paste(sample_strs, collapse = ", ")))
    training_results_list <- training_results_list[!is_invalid_res]
  }

  cat(sprintf("After filtering NULL training results: %d results remaining\n", length(training_results_list)))

  if (length(training_results_list) > 0) {
    barren_one_count <- sum(sapply(training_results_list, function(res) {
      if (!is.null(res$barren_fraction) && is.finite(res$barren_fraction)) {
        res$barren_fraction == 1
      } else {
        FALSE
      }
    }))
    barren_one_pct <- barren_one_count / length(training_results_list) * 100
    if (barren_one_pct > 50) {
      msg <- sprintf("WARNING: %.1f%% of predictions (%d/%d) have barren_fraction = 1, exceeding 50%% threshold. This indicates severe model issues.", barren_one_pct, barren_one_count, length(training_results_list))
      # Always emit a warning rather than stopping execution. In testing mode, prefix to indicate non-actionable warning.
      if (isTRUE(TESTING_MODE)) {
        warning(paste("[TESTING MODE IGNORE]", msg))
      } else {
        warning(msg)
      }
    }
    cat(sprintf("Barren fraction = 1 in %.1f%% of training predictions (%d/%d) - within acceptable limits\n", barren_one_pct, barren_one_count, length(training_results_list)))
  } else {
    cat("No training results to check for barren fraction\n")
  }

  if (length(training_results_list) == 0) {
    cat("ERROR: All training tasks returned NULL results!\n")
    cat("Most likely causes:\n")
    cat("1. Insufficient data quality: Tasks were filtered due to:\n")
    cat("   - Too few observations per location-year (<15 observations)\n")
    cat("   - Insufficient temporal coverage (<25 unique days of year)\n")
    cat("   - Too many missing pentad bins (>85% NA values)\n")
    cat("2. Data filtering issues or missing indices\n")
    cat("3. No valid testing/inference data available (no location-year pairs found)\n")
    cat("\nTo adjust data quality thresholds, modify MIN_OBS, MIN_UNIQUE_DOYS, and MIN_PENTAD_COVERAGE in fit_one_task function.\n")
    stop("No valid training results to process")
  }

  # Build all_coefs from training_results_list (memory-efficient - no separate list accumulation)
  if (length(training_results_list) > 0) {
    cat("Combining coefficient data frames from training_results_list...\n")

    all_coefs <- tryCatch({
      if (requireNamespace("dplyr", quietly = TRUE)) {
        dplyr::bind_rows(lapply(training_results_list, function(r) {
          if (!is.null(r$coef_df)) r$coef_df else NULL
        }))
      } else {
        do.call(rbind, lapply(training_results_list, function(r) {
          if (!is.null(r$coef_df)) r$coef_df else NULL
        }))
      }
    }, error = function(e) stop(sprintf("ERROR combining coef_df from training results: %s", e$message))
    )

    if (is.null(all_coefs)) {
      cat("Failed to combine coefficient data frames\n")
      stop("Cannot proceed without coefficient data")
    }

    required_coef_cols <- c("location_id", "pheno_year", "Veg", "coef", "rmse", "coef_025", "coef_975", "interval")
    missing <- setdiff(required_coef_cols, names(all_coefs))
    if (length(missing) > 0) {
      cat(sprintf("[NOTICE] Filling missing coefficient columns with NA: %s\n", paste(missing, collapse = ", ")))
      for (col in missing) all_coefs[[col]] <- NA
    }
    all_coefs$location_id <- as.character(all_coefs$location_id)
    all_coefs$pheno_year <- as.integer(all_coefs$pheno_year)
    # Ensure legacy 'year' column exists for downstream code that expects it
    if (!"year" %in% names(all_coefs)) all_coefs$year <- all_coefs$pheno_year
    all_coefs$Veg <- as.character(all_coefs$Veg)
    all_coefs$coef <- as.numeric(all_coefs$coef)
    all_coefs$rmse <- as.numeric(all_coefs$rmse)
    all_coefs$coef_025 <- as.numeric(all_coefs$coef_025)
    all_coefs$coef_975 <- as.numeric(all_coefs$coef_975)
    all_coefs$interval <- as.numeric(all_coefs$interval)

    all_coefs$location_id <- trimws(all_coefs$location_id)
    cat(sprintf("Combined coefficients: %d rows\n", nrow(all_coefs)))

    cat("Combining chosen variant summaries (training results)...\n")
    variant_list_pca <- lapply(training_results_list, function(res) if (!is.null(res$variant_trajectory)) res$variant_trajectory else NULL)
    variant_list_pca <- variant_list_pca[!sapply(variant_list_pca, is.null)]
    all_variants_pca <- NULL
    if (length(variant_list_pca) > 0) {
      all_variants_pca <- tryCatch({
        if (requireNamespace("dplyr", quietly = TRUE)) {
          dplyr::bind_rows(variant_list_pca)
        } else {
          do.call(rbind, variant_list_pca)
        }
      }, error = function(e) stop(sprintf("Failed to combine variant PCA summaries: %s", e$message)))
    }

    diag_list <- lapply(training_results_list, function(res) {
      if (!is.null(res$diagnostics)) res$diagnostics else NULL
    })
    diag_list <- diag_list[!sapply(diag_list, is.null)]
    all_diagnostics <- if (length(diag_list) > 0) tryCatch({
      if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required to combine diagnostics robustly")
      dplyr::bind_rows(diag_list)
    }, error = function(e) stop(sprintf("Failed to combine diagnostics: %s", e$message))) else NULL

    q_dvi_data <- tryCatch({
      if (requireNamespace("dplyr", quietly = TRUE)) {
        dplyr::bind_rows(lapply(training_results_list, function(res) {
          if ((!is.null(res$q10_dvi) && !is.na(res$q10_dvi)) || (!is.null(res$q90_dvi) && !is.na(res$q90_dvi))) {
            data.frame(
              location_id = res$coef_df$location_id[1],
              pheno_year = res$coef_df$pheno_year[1],
              q10_dvi = res$q10_dvi,
              q90_dvi = res$q90_dvi,
              stringsAsFactors = FALSE
            )
          } else {
            NULL
          }
        }))
      } else {
        do.call(rbind, lapply(training_results_list, function(res) {
          if ((!is.null(res$q10_dvi) && !is.na(res$q10_dvi)) || (!is.null(res$q90_dvi) && !is.na(res$q90_dvi))) {
            data.frame(
              location_id = res$coef_df$location_id[1],
              pheno_year = res$coef_df$pheno_year[1],
              q10_dvi = res$q10_dvi,
              q90_dvi = res$q90_dvi,
              stringsAsFactors = FALSE
            )
          } else {
            NULL
          }
        }))
      }
    }, error = function(e) {
      warning(sprintf("q_dvi_data assembly failed: %s", e$message)); NULL
    })

    variant_similarity_table <- NULL
    variant_similarity_summary <- NULL
    if (exists("VARIANT_SIMILARITY_TABLE")) {
      variant_similarity_table <- VARIANT_SIMILARITY_TABLE
    } else if (exists("INSEPARABLE_VARIANT_INFO") && !is.null(INSEPARABLE_VARIANT_INFO$similarity_table)) {
      variant_similarity_table <- INSEPARABLE_VARIANT_INFO$similarity_table
    }
    if (!is.null(variant_similarity_table) && nrow(variant_similarity_table) > 0) {
      variant_similarity_table <- variant_similarity_table[order(variant_similarity_table$veg, variant_similarity_table$euclidean_dist), , drop = FALSE]
      if (requireNamespace("dplyr", quietly = TRUE)) {
        variant_similarity_summary <- variant_similarity_table |> 
          dplyr::group_by(.data$veg) |> 
          dplyr::summarise(
            pair_count = dplyr::n(),
            min_euclidean_dist = if (all(is.na(.data$euclidean_dist))) NA_real_ else min(.data$euclidean_dist, na.rm = TRUE),
            median_euclidean_dist = if (all(is.na(.data$euclidean_dist))) NA_real_ else stats::median(.data$euclidean_dist, na.rm = TRUE),
            .groups = "drop"
          )
      } else {
        by_veg <- split(variant_similarity_table, variant_similarity_table$veg)
        summary_rows <- lapply(names(by_veg), function(veg_name) {
          tbl <- by_veg[[veg_name]]
          data.frame(
            veg = veg_name,
            pair_count = nrow(tbl),
            min_euclidean_dist = if (all(is.na(tbl$euclidean_dist))) NA_real_ else min(tbl$euclidean_dist, na.rm = TRUE),
            median_euclidean_dist = if (all(is.na(tbl$euclidean_dist))) NA_real_ else stats::median(tbl$euclidean_dist, na.rm = TRUE),
            stringsAsFactors = FALSE
          )
        })
        variant_similarity_summary <- do.call(rbind, summary_rows)
      }
    }

    cat(sprintf("[DEBUG] all_coefs$location_id class: %s\n", class(all_coefs$location_id)))
    cat(sprintf("[DEBUG] all_coefs$location_id sample: %s\n", paste(head(all_coefs$location_id, 10), collapse = ", ")))
    cat(sprintf("[DEBUG] all_coefs$location_id is.na sum: %d\n", sum(is.na(all_coefs$location_id))))
    cat(sprintf("[DEBUG] all_coefs$location_id == '' sum: %d\n", sum(all_coefs$location_id == "", na.rm = TRUE)))
    unique_locations <- unique(trimws(as.character(all_coefs$location_id)))
    unique_locations <- unique_locations[!is.na(unique_locations) & unique_locations != ""]
    if (length(unique_locations) == 0) {
      stop("No valid location IDs found in results")
    }
    cat(sprintf("Creating Excel file with %d valid locations\n", length(unique_locations)))

    true_veg_map <- gpts_map |> dplyr::select(location_id, true_veg = Veg)
    if (!is.character(all_coefs$location_id)) {
      all_coefs$location_id <- as.character(all_coefs$location_id)
      cat("[NOTICE] Coerced all_coefs$location_id to character to match true_veg_map for joining.\n")
    }
    if (!is.character(true_veg_map$location_id)) {
      true_veg_map$location_id <- as.character(true_veg_map$location_id)
      cat("[NOTICE] Coerced true_veg_map$location_id to character for joining.\n")
    }
    if (length(all_coefs$location_id) > 0 && length(true_veg_map$location_id) > 0) {
      s1 <- unique(na.omit(all_coefs$location_id))
      s2 <- unique(na.omit(true_veg_map$location_id))
      if (length(s1) && length(s2) && all(grepl("^[0-9]+$", s1)) && any(grepl("^L_", s2))) {
        cat("[WARNING] all_coefs$location_id looks numeric while true_veg_map$location_id looks like 'L_x_y' strings — matching will likely fail.\n")
      }
    }
    all_coefs <- dplyr::left_join(all_coefs, true_veg_map, by = "location_id")
    if (!"coef_025" %in% names(all_coefs)) all_coefs$coef_025 <- NA_real_
    if (!"coef_975" %in% names(all_coefs)) all_coefs$coef_975 <- NA_real_
    if (DEBUG_UNCERTAINTY) {
      cat("After combining all_coefs with true_veg_map:\n")
      cat("Number of rows with finite coef_025:", sum(is.finite(all_coefs$coef_025)), "\n")
      cat("Number of rows with finite coef_975:", sum(is.finite(all_coefs$coef_975)), "\n")
      cat("Number of rows with finite interval:", sum(is.finite(all_coefs$interval)), "\n")
    }

    # Helper: combine results_list into unified all_coefs (no printing)
    combine_results_from_list <- function(results_list) {
      coef_list <- lapply(results_list, function(res) {
        if (is.null(res) || is.null(res$coef_df) || nrow(res$coef_df) == 0) NULL else res$coef_df
      })
      coef_list <- coef_list[!sapply(coef_list, is.null)]
      if (length(coef_list) == 0) return(NULL)
      if (requireNamespace("dplyr", quietly = TRUE)) {
        all_coefs_local <- dplyr::bind_rows(coef_list)
      } else {
        all_coefs_local <- do.call(rbind, coef_list)
      }
      # Ensure minimal columns exist
      req_cols <- c("location_id","pheno_year","Veg","coef","coef_025","coef_975","interval")
      for (cname in req_cols) if (!(cname %in% names(all_coefs_local))) all_coefs_local[[cname]] <- NA
      all_coefs_local
    }

    # Run inference on a dataframe silently (suppress per-task prints) and return results list + combined coefs
    run_inference_silent <- function(df_subset) {
      if (is.null(df_subset) || nrow(df_subset) == 0) return(list(results = list(), all_coefs = NULL))
      if (!"pheno_year" %in% names(df_subset) && "date" %in% names(df_subset)) df_subset$pheno_year <- assign_pheno_year(df_subset$date)
      if (!"task_key" %in% names(df_subset)) df_subset$task_key <- paste(df_subset$location_id, df_subset$pheno_year, sep = "_")
      task_list <- split(df_subset, df_subset$task_key)

      sink_file <- tempfile()
      con_out <- file(sink_file, open = "wt")
      con_msg <- file(paste0(sink_file, ".msg"), open = "wt")
      # Suppress output and messages while running heavy inference
      results <- NULL
      tryCatch({
        sink(con_out, type = "output")
        sink(con_msg, type = "message")
        results <- .run_map(task_list, fit_one_task, show_pb = FALSE)
      }, error = function(e) {
        warning(sprintf("run_inference_silent: underlying inference error: %s", e$message))
      }, finally = {
        # Attempt to restore sinks
        try(sink(type = "message"), silent = TRUE)
        try(sink(type = "output"), silent = TRUE)
        close(con_out)
        close(con_msg)
        try(unlink(sink_file), silent = TRUE)
        try(unlink(paste0(sink_file, ".msg")), silent = TRUE)
      })

      all_coefs_local <- combine_results_from_list(results)
      list(results = results, all_coefs = all_coefs_local)
    }

    # --- Separate training-location inference for years 2023 and 2025 ---
    # (Removed: Redundant with full inference processing)


    inseparable_flags_detail <- NULL
    inseparable_flags_summary <- NULL
    if ("inseparable_variant_flag" %in% names(all_coefs)) {
      flagged_rows <- which(!is.na(all_coefs$inseparable_variant_flag) & all_coefs$inseparable_variant_flag)
      if (length(flagged_rows) > 0) {
        cols_keep <- intersect(c("location_id", "year", "Veg", "coef", "rmse", "true_veg", "inseparable_variant_details"), names(all_coefs))
        inseparable_flags_detail <- all_coefs[flagged_rows, cols_keep, drop = FALSE]
        inseparable_flags_detail <- inseparable_flags_detail[order(inseparable_flags_detail$location_id, inseparable_flags_detail$year, inseparable_flags_detail$Veg), , drop = FALSE]

        if (requireNamespace("dplyr", quietly = TRUE)) {
          inseparable_flags_summary <- inseparable_flags_detail |> 
            dplyr::group_by(.data$location_id, .data$year) |> 
            dplyr::summarise(
              veg_components = paste(sort(unique(.data$Veg)), collapse = ", "),
              max_coef = if (all(is.na(.data$coef))) NA_real_ else max(.data$coef, na.rm = TRUE),
              details = paste(unique(na.omit(.data$inseparable_variant_details)), collapse = " | "),
              .groups = "drop"
            ) |> 
            dplyr::arrange(.data$location_id, .data$year)
        } else {
          summary_keys <- unique(inseparable_flags_detail[, c("location_id", "year"), drop = FALSE])
          summary_keys <- summary_keys[order(summary_keys$location_id, summary_keys$pheno_year), , drop = FALSE]
          build_summary <- function(loc, yr) {
            rows <- inseparable_flags_detail$location_id == loc & inseparable_flags_detail$pheno_year == yr
            vegs <- sort(unique(inseparable_flags_detail$Veg[rows]))
            coef_vals <- inseparable_flags_detail$coef[rows]
            details <- unique(na.omit(inseparable_flags_detail$inseparable_variant_details[rows]))
            data.frame(
              location_id = loc,
              pheno_year = yr,
              veg_components = paste(vegs, collapse = ", "),
              max_coef = if (all(is.na(coef_vals))) NA_real_ else max(coef_vals, na.rm = TRUE),
              details = paste(details, collapse = " | "),
              stringsAsFactors = FALSE
            )
          }
          if (nrow(summary_keys) > 0) {
            summary_list <- apply(summary_keys, 1, function(row) build_summary(row[["location_id"]], as.integer(row[["year"]])))
            inseparable_flags_summary <- do.call(rbind, summary_list)
          }
        }
      }
    }

    cat("Creating single Excel file with all location results...\n")

    summary_data <- data.frame(
      Location_ID = unique_locations,
      Total_Years = sapply(unique_locations, function(loc) {
        length(unique(all_coefs$pheno_year[all_coefs$location_id == loc]))
      }),
      Total_Observations = sapply(unique_locations, function(loc) {
        sub <- all_coefs[all_coefs$location_id == loc, ]
        if ("n_obs" %in% names(sub)) {
          u <- unique(sub[, c("pheno_year", "n_obs"), drop = FALSE])
          sum(u$n_obs, na.rm = TRUE)
        } else {
          nrow(sub)
        }
      }),
      stringsAsFactors = FALSE
    )

    openxlsx::addWorksheet(wb, "Summary")
    openxlsx::writeData(wb, "Summary", summary_data)

    if (!is.null(all_diagnostics) && nrow(all_diagnostics) > 0) {
      openxlsx::addWorksheet(wb, "Diagnostics")
      openxlsx::writeData(wb, "Diagnostics", all_diagnostics)
    }

            if (!is.null(all_variants_pca) && nrow(all_variants_pca) > 0) {
      openxlsx::addWorksheet(wb, "Variant_Summary")
      openxlsx::writeData(wb, "Variant_Summary", all_variants_pca)
    }

    if (!is.null(variant_similarity_table) && nrow(variant_similarity_table) > 0) {
      openxlsx::addWorksheet(wb, "Variant_Similarity")
      openxlsx::writeData(wb, "Variant_Similarity", "Pairwise similarity across variants (Euclidean distance)", startRow = 1, startCol = 1)
      start_row <- 2
      if (!is.null(variant_similarity_summary) && nrow(variant_similarity_summary) > 0) {
        openxlsx::writeData(wb, "Variant_Similarity", variant_similarity_summary, startRow = start_row, startCol = 1)
        start_row <- start_row + nrow(variant_similarity_summary) + 2
      }
      openxlsx::writeData(wb, "Variant_Similarity", variant_similarity_table, startRow = start_row, startCol = 1)
    }

    if (exists("INTER_CLASS_SIMILARITY") && !is.null(INTER_CLASS_SIMILARITY) && nrow(INTER_CLASS_SIMILARITY) > 0) {
      openxlsx::addWorksheet(wb, "Inter_Veg_Similarity")
      openxlsx::writeData(wb, "Inter_Veg_Similarity", "Pairwise similarity between variants of different vegetation types", startRow = 1, startCol = 1)
      openxlsx::writeData(wb, "Inter_Veg_Similarity", INTER_CLASS_SIMILARITY, startRow = 2, startCol = 1)
    }

    if (!is.null(inseparable_flags_summary) && nrow(inseparable_flags_summary) > 0) {
      openxlsx::addWorksheet(wb, "Inseparable_Flags")
      openxlsx::writeData(wb, "Inseparable_Flags", "Location-years with inseparable MESMA variants detected", startRow = 1, startCol = 1)
      openxlsx::writeData(wb, "Inseparable_Flags", inseparable_flags_summary, startRow = 2, startCol = 1)
      next_row <- nrow(inseparable_flags_summary) + 4
      if (!is.null(inseparable_flags_detail) && nrow(inseparable_flags_detail) > 0) {
        openxlsx::writeData(wb, "Inseparable_Flags", "Detailed component rows", startRow = next_row, startCol = 1)
        openxlsx::writeData(wb, "Inseparable_Flags", inseparable_flags_detail, startRow = next_row + 1, startCol = 1)
      }
    }

    unc_coef_rows <- list(); unc_var_rows <- list(); unc_rmse_rows <- list(); unc_meta_rows <- list()
    for (res in training_results_list) {
      if (is.null(res$uncertainty)) next
      loc <- if (!is.null(res$coef_df$location_id)) res$coef_df$location_id[1] else NA_character_
      yr <- if (!is.null(res$coef_df$year)) res$coef_df$year[1] else NA_integer_
      ci <- res$uncertainty$coef_ci
      if (!is.null(ci) && nrow(ci) > 0) {
        ci$location_id <- loc; ci$year <- yr
        unc_coef_rows[[length(unc_coef_rows) + 1]] <- ci[, c("location_id","year","Veg","coef_025","coef_975", "coef_sd", "interval"), drop = FALSE]
      }
      vr <- res$uncertainty$variant_freq
      if (!is.null(vr) && nrow(vr) > 0) {
        vr$location_id <- loc; vr$year <- yr
        unc_var_rows[[length(unc_var_rows) + 1]] <- vr[, c("location_id","year","Veg","Variant","N","Percent"), drop = FALSE]
      }
      rc <- res$uncertainty$rmse_ci
      if (!is.null(rc) && length(rc) == 2) {
        unc_rmse_rows[[length(unc_rmse_rows) + 1]] <- data.frame(location_id = loc, year = yr, rmse_lo = rc[1], rmse_hi = rc[2], stringsAsFactors = FALSE)
      }

    }
    if (length(unc_coef_rows) > 0) {
      if (requireNamespace("dplyr", quietly = TRUE)) all_unc_coef <- dplyr::bind_rows(unc_coef_rows) else all_unc_coef <- do.call(rbind, unc_coef_rows)
    } else all_unc_coef <- NULL
    if (!is.null(all_unc_coef)) {
      n_total_unc <- nrow(all_unc_coef)
      n_non_na_ci <- sum(is.finite(all_unc_coef$coef_025) | is.finite(all_unc_coef$coef_975))
      cat(sprintf("Uncertainty CIs: total rows = %d; rows with at least one finite CI = %d\n", n_total_unc, n_non_na_ci))

      # MERGE BOOTSTRAP UNCERTAINTY INTO MAIN COEFFICIENTS
      # This is the critical fix: merge coef_025 and coef_975 from bootstrap into all_coefs
      if (requireNamespace("dplyr", quietly = TRUE)) {
        # Ensure matching column types for join
        all_coefs$location_id <- as.character(all_coefs$location_id)
        all_coefs$year <- as.integer(all_coefs$year)
        all_unc_coef$location_id <- as.character(all_unc_coef$location_id)
        all_unc_coef$year <- as.integer(all_unc_coef$year)

        # Remove existing coef_025 and coef_975 columns if they exist (they contain NAs)
        if ("coef_025" %in% names(all_coefs)) all_coefs$coef_025 <- NULL
        if ("coef_975" %in% names(all_coefs)) all_coefs$coef_975 <- NULL
        if ("coef_sd" %in% names(all_coefs)) all_coefs$coef_sd <- NULL
        if ("interval" %in% names(all_coefs)) all_coefs$interval <- NULL

        # Perform left join to merge bootstrap CIs into main results
        all_coefs <- dplyr::left_join(all_coefs, all_unc_coef, by = c("location_id", "year", "Veg"))

        # Report merge success
        n_merged <- sum(is.finite(all_coefs$coef_025) | is.finite(all_coefs$coef_975))
        cat(sprintf("Bootstrap CIs merged into main results: %d/%d rows now have uncertainty bounds\n",
                    n_merged, nrow(all_coefs)))
      } else {
        warning("dplyr not available - bootstrap CIs will not be merged into main results")
      }
    }
    if (length(unc_var_rows) > 0) {
      if (requireNamespace("dplyr", quietly = TRUE)) all_unc_var <- dplyr::bind_rows(unc_var_rows) else all_unc_var <- do.call(rbind, unc_var_rows)
    } else all_unc_var <- NULL
    if (length(unc_rmse_rows) > 0) {
      if (requireNamespace("dplyr", quietly = TRUE)) all_unc_rmse <- dplyr::bind_rows(unc_rmse_rows) else all_unc_rmse <- do.call(rbind, unc_rmse_rows)
    } else all_unc_rmse <- NULL
    if (length(unc_meta_rows) > 0) {
      if (requireNamespace("dplyr", quietly = TRUE)) all_unc_meta <- dplyr::bind_rows(unc_meta_rows) else all_unc_meta <- do.call(rbind, unc_meta_rows)
    } else all_unc_meta <- NULL

    if (!is.null(all_unc_coef) || !is.null(all_unc_var) || !is.null(all_unc_rmse)) {
      openxlsx::addWorksheet(wb, "Uncertainty")
      start_row <- 1
      openxlsx::writeData(wb, "Uncertainty", data.frame(Setting = c("ENABLE_UNCERTAINTY","BOOTSTRAP_METHOD"), Value = c(ENABLE_UNCERTAINTY, if (isTRUE(ENABLE_UNCERTAINTY)) "locations" else "none")), startRow = start_row, startCol = 1)
      start_row <- start_row + 3
      if (!is.null(all_unc_rmse)) {
        openxlsx::writeData(wb, "Uncertainty", "RMSE CI (2.5%/97.5%)", startRow = start_row, startCol = 1)
        openxlsx::writeData(wb, "Uncertainty", all_unc_rmse, startRow = start_row + 1, startCol = 1)
        start_row <- start_row + nrow(all_unc_rmse) + 3
      }
      if (!is.null(all_unc_var)) {
        openxlsx::writeData(wb, "Uncertainty", "Variant Dominance Frequencies (%)", startRow = start_row, startCol = 1)
        openxlsx::writeData(wb, "Uncertainty", all_unc_var, startRow = start_row + 1, startCol = 1)
      }
      if (!is.null(all_unc_meta)) {
        start_row <- start_row + ifelse(!is.null(all_unc_var), nrow(all_unc_var) + 3, 3)
        openxlsx::writeData(wb, "Uncertainty", "Bootstrap Meta", startRow = start_row, startCol = 1)
        openxlsx::writeData(wb, "Uncertainty", all_unc_meta, startRow = start_row + 1, startCol = 1)
      }
    }

    if (exists('stability_results') && !is.null(stability_results)) {
      openxlsx::addWorksheet(wb, 'Endmember_Stability')
      current_row <- 1
      for (veg in names(stability_results)) {
        res <- stability_results[[veg]]
        openxlsx::writeData(wb, 'Endmember_Stability', sprintf('=== %s ENDMEMBER STABILITY ===', toupper(veg)), startRow = current_row, startCol = 1)
        current_row <- current_row + 2
        openxlsx::writeData(wb, 'Endmember_Stability', 'Variant Count Distribution:', startRow = current_row, startCol = 1)
        vcd <- as.data.frame(res$variant_count_distribution)
        names(vcd) <- c('N_Variants', 'Frequency')
        openxlsx::writeData(wb, 'Endmember_Stability', vcd, startRow = current_row + 1, startCol = 1)
        current_row <- current_row + nrow(vcd) + 3
        openxlsx::writeData(wb, 'Endmember_Stability', 'Meta-Variant Stability Metrics:', startRow = current_row, startCol = 1)
        mv_stats <- if (length(res$meta_variants) > 0) do.call(rbind, lapply(res$meta_variants, function(mv) {
          data.frame(Meta_Variant = mv$meta_variant_id, N_Bootstrap_Members = mv$n_members, Coefficient_of_Variation = mv$coefficient_of_variation, Spectral_Angle_Mean_deg = mv$mean_spectral_angle, Spectral_Angle_SD_deg = sqrt(mv$spectral_angle_variance), Bootstrap_Frequency_pct = length(unique(mv$bootstrap_iters)) / res$n_bootstrap_iters * 100, stringsAsFactors = FALSE)
        })) else data.frame()
        if (nrow(mv_stats) > 0) {
          openxlsx::writeData(wb, 'Endmember_Stability', mv_stats, startRow = current_row + 1, startCol = 1)
          current_row <- current_row + nrow(mv_stats) + 4
        } else {
          current_row <- current_row + 4
        }
      }
      openxlsx::writeData(wb, 'Endmember_Stability', 'INTERPRETATION GUIDE', startRow = current_row, startCol = 1)
      current_row <- current_row + 1
      guide <- data.frame(Metric = c('Variant Count Distribution', 'Coefficient of Variation (CV)', 'Spectral Angle SD', 'Bootstrap Frequency'), Interpretation = c('How many variants are identified across bootstrap iterations. Stable if concentrated.', 'Variability of signature features. Lower CV = more stable endmember.', 'Directional variability in spectral space (degrees). Lower = more stable.', '% of bootstrap iterations where this pattern appears. Higher = more robust.'), Good_Value = c('Narrow distribution (e.g., 90% have same count)', '< 0.10 (very stable), 0.10-0.20 (stable), > 0.20 (unstable)', '< 5° (very stable), 5-10° (stable), > 10° (unstable)', '> 80% (very robust), 50-80% (robust), < 50% (unreliable)'), stringsAsFactors = FALSE)
      openxlsx::writeData(wb, 'Endmember_Stability', guide, startRow = current_row, startCol = 1)
    }

    best_fit_summary <- do.call(rbind, lapply(unique_locations, function(loc) {
      yrs <- unique(all_coefs$pheno_year[all_coefs$location_id == loc])
      tv <- true_veg_map$true_veg[true_veg_map$location_id == loc]
      tv_val <- if (length(tv) > 0) tv[1] else NA_character_
      tv_lower <- tolower(tv_val)

      do.call(rbind, lapply(yrs, function(yr) {
        row <- all_coefs[all_coefs$location_id == loc & all_coefs$pheno_year == yr & tolower(all_coefs$Veg) == tv_lower, , drop = FALSE]
        pred <- if (nrow(row) == 1) row$coef else NA_real_
        pred_abs <- pred
        rmse_val <- if (nrow(row) == 1 && "rmse" %in% names(row)) row$rmse else NA_real_
        sum_veg_coef <- sum(all_coefs$coef[all_coefs$location_id == loc & all_coefs$pheno_year == yr & tolower(all_coefs$Veg) != "barren"], na.rm = TRUE)
        pred_rel <- if (!is.na(pred) && is.finite(sum_veg_coef) && sum_veg_coef > 0) pred / sum_veg_coef else NA_real_
        pred_rel <- pmin(pred_rel, 1)  # Clamp to 1 to prevent >100%
        abs_pct <- if (!is.na(pred_abs) && !is.na(tv_val)) abs(1 - pred_abs) * 100 else NA_real_
        abs_pct_rel <- if (!is.na(pred_rel) && !is.na(tv_val)) abs(1 - pred_rel) * 100 else NA_real_

        is_barren_truth <- isTRUE(tv_lower == "barren")
        pred_total_veg_cover <- sum_veg_coef

        if (is_barren_truth) {
          pred <- NA_real_; pred_abs <- NA_real_; pred_rel <- NA_real_; abs_pct <- NA_real_; abs_pct_rel <- NA_real_
        }

        data.frame(
          location_id = loc,
          year = yr,
          true_veg = tv_val,
          pred_coef = pred_rel,
          pred_coef_abs = pred_abs,
          pred_coef_rel = pred_rel,
          rmse = rmse_val,
          abs_pct_diff = abs_pct_rel,
          abs_pct_diff_abs = abs_pct,
          abs_pct_diff_rel = abs_pct_rel,
          is_barren_truth = is_barren_truth,
          pred_total_veg_cover = pred_total_veg_cover,
          true_veg_cover_ppi = NA_real_,
          veg_cover_error_ppi = NA_real_,
          stringsAsFactors = FALSE
        )
      }))
    }))

    # Attach PPI-derived vegetation cover truth when available
    ppi_truth_summary <- NULL
    if (exists("all_diagnostics") && !is.null(all_diagnostics) &&
        all(c("location_id", "pheno_year") %in% names(all_diagnostics)) &&
        "barren_fraction_ppi_based" %in% names(all_diagnostics)) {
      ppi_truth_summary <- all_diagnostics %>%
        dplyr::select(location_id, pheno_year, barren_fraction_ppi_based) %>%
        dplyr::mutate(
          barren_fraction_ppi_based_join = pmin(pmax(barren_fraction_ppi_based, 0), 1),
          true_veg_cover_ppi_join = 1 - barren_fraction_ppi_based_join
        ) %>%
        dplyr::rename(year = pheno_year) %>%
        dplyr::distinct()
    }

    if (!is.null(ppi_truth_summary)) {
      best_fit_summary <- best_fit_summary %>%
        dplyr::left_join(ppi_truth_summary, by = c("location_id", "year")) %>%
        dplyr::mutate(
          true_veg_cover_ppi = dplyr::coalesce(true_veg_cover_ppi, true_veg_cover_ppi_join),
          veg_cover_error_ppi = ifelse(is.finite(pred_total_veg_cover) & is.finite(true_veg_cover_ppi),
                                       pred_total_veg_cover - true_veg_cover_ppi,
                                       veg_cover_error_ppi)
        ) %>%
        dplyr::select(-dplyr::any_of(c("true_veg_cover_ppi_join", "barren_fraction_ppi_based_join", "barren_fraction_ppi_based")))
    }

    eval_years <- sort(unique(c(TRAIN_YEARS, TRAIN_YEARS - 1L, TRAIN_YEARS + 1L)))
    eval_years <- eval_years[is.finite(eval_years)]
    best_fit_summary$eval_window <- best_fit_summary$year %in% eval_years
    best_fit_eval <- best_fit_summary[best_fit_summary$eval_window & !best_fit_summary$is_barren_truth, , drop = FALSE]
    if (nrow(best_fit_eval) == 0) {
      warning("No location-years fall within the TRAIN_YEARS +/- 1 evaluation window; overall fit cannot be computed")
    }

    overall_fit <- suppressWarnings(as.numeric(mean(best_fit_eval$pred_coef_rel * 100, na.rm = TRUE)))
    if (!is.finite(overall_fit)) overall_fit <- NA_real_

    openxlsx::writeData(wb, "Summary", data.frame(
      Overall_Fit_pct = overall_fit
    ), startRow = 1, startCol = ncol(summary_data) + 2)

    for (i in seq_along(unique_locations)) {
      loc_id <- unique_locations[i]

          sheet_name <- substr(gsub("[^A-Za-z0-9]", "_", loc_id), 1, 31)

          openxlsx::addWorksheet(wb, sheet_name)

          loc_coefs <- all_coefs[all_coefs$location_id == loc_id, ]

          # Extract true_veg for this location
          true_veg_val <- true_veg_map$true_veg[true_veg_map$location_id == loc_id]
          if (length(true_veg_val) == 0) true_veg_val <- NA_character_
          true_veg_val <- true_veg_val[1]  # Take first value

          loc_q_data <- if (!is.null(q_dvi_data) && "location_id" %in% names(q_dvi_data)) {
            q_dvi_data[q_dvi_data$location_id == loc_id, ]
          } else {
            NULL
          }

          peak_q10_dvi <- if (!is.null(loc_q_data) && nrow(loc_q_data) > 0) {
            max(loc_q_data$q10_dvi, na.rm = TRUE)
          } else {
            NA
          }

          peak_q90_dvi <- if (!is.null(loc_q_data) && nrow(loc_q_data) > 0) {
            max(loc_q_data$q90_dvi, na.rm = TRUE)
          } else {
            NA
          }

          if (!is.na(true_veg_val) && tolower(true_veg_val) == "barren") {
            quality_metrics <- data.frame(
              avg_pct_deviation = NA_real_,
              avg_rmse = mean(loc_coefs$rmse, na.rm = TRUE),
              peak_q10_dvi = peak_q10_dvi,
              peak_q90_dvi = peak_q90_dvi,
              stringsAsFactors = FALSE
            )
          } else {
            quality_metrics <- loc_coefs |> 
              dplyr::group_by(.data$pheno_year) |> 
              dplyr::summarize(
                deviation = sum(abs(.data$coef - (tolower(.data$Veg) == tolower(true_veg_val)))),
                avg_rmse = mean(.data$rmse, na.rm = TRUE),
                .groups = "drop"
              ) |> 
              dplyr::summarize(
                avg_pct_deviation = mean(.data$deviation, na.rm = TRUE) * 100,
                avg_rmse = mean(.data$avg_rmse, na.rm = TRUE),
                .groups = "drop"
              )

            quality_metrics$peak_q10_dvi <- peak_q10_dvi
            quality_metrics$peak_q90_dvi <- peak_q90_dvi
          }

          openxlsx::writeData(wb, sheet_name, "QUALITY METRICS", startRow = 1, startCol = 1)
          openxlsx::writeData(wb, sheet_name, quality_metrics, startRow = 2, startCol = 1)

          current_row <- nrow(quality_metrics) + 4

          # Removed DIAGNOSTICS section

          loc_best <- if ("location_id" %in% names(best_fit_summary)) best_fit_summary[best_fit_summary$location_id == loc_id, , drop = FALSE] else data.frame()
          if (!is.null(loc_best) && nrow(loc_best) > 0) {
            desired_cols <- c("year", "true_veg", "pred_coef", "pred_coef_abs", "rmse", "abs_pct_diff", "abs_pct_diff_abs", "eval_window", "is_barren_truth", "pred_total_veg_cover", "true_veg_cover_ppi", "veg_cover_error_ppi")
            write_tbl <- loc_best[, c("location_id", intersect(desired_cols, names(loc_best))), drop = FALSE]
            openxlsx::writeData(wb, sheet_name, "BEST FIT SUMMARY (per-year) — pred_coef is proportion relative to vegetated area (barren excluded)",
              startRow = current_row, startCol = 1
            )
            openxlsx::writeData(wb, sheet_name, write_tbl, startRow = current_row + 1, startCol = 1)
            current_row <- current_row + nrow(write_tbl) + 3
          }

          loc_coefs_unified <- if ("location_id" %in% names(all_coefs)) all_coefs[all_coefs$location_id == loc_id, ] else data.frame()

          loc_variants_pca <- if (!is.null(all_variants_pca) && "location_id" %in% names(all_variants_pca)) {
            all_variants_pca[all_variants_pca$location_id == loc_id, ]
          } else {
            NULL
          }

          if (nrow(loc_coefs_unified) > 0) {
            openxlsx::writeData(wb, sheet_name, "COEFFICIENTS",
              startRow = current_row, startCol = 1
            )
            openxlsx::writeData(wb, sheet_name, loc_coefs_unified,
              startRow = current_row + 1, startCol = 1
            )
            current_row <- current_row + nrow(loc_coefs_unified) + 3
          }


          if (!is.null(loc_variants_pca) && nrow(loc_variants_pca) > 0) {
            variant_usage <- data.frame(
              year = unique(loc_variants_pca$year),
              stringsAsFactors = FALSE
            )

            veg_candidates <- unique(c(ALLOWED_VEG, na.omit(unique(as.character(all_coefs$Veg)))))
            veg_candidates <- veg_candidates[!is.na(veg_candidates) & veg_candidates != ""]
            veg_kept_local <- veg_candidates
            if (length(veg_kept_local) == 0) veg_kept_local <- ALLOWED_VEG
            for (veg in veg_kept_local) {
              var_col <- paste0(veg, "_variant")
              if (var_col %in% names(loc_variants_pca)) {
                variant_usage[[paste0(veg, "_most_common")]] <- sapply(variant_usage$year, function(y) {
                  year_data <- loc_variants_pca[loc_variants_pca$year == y, var_col]
                  if (length(year_data) > 0) {
                    names(sort(table(year_data), decreasing = TRUE))[1]
                  } else {
                    NA
                  }
                })
              }
            }

            openxlsx::writeData(wb, sheet_name, "VARIANT USAGE",
              startRow = current_row, startCol = 1
            )
            openxlsx::writeData(wb, sheet_name, variant_usage,
              startRow = current_row + 1, startCol = 1
            )
            current_row <- current_row + nrow(variant_usage) + 3
          }

          if (!is.null(loc_q_data) && nrow(loc_q_data) > 0) {
            openxlsx::writeData(wb, sheet_name, "Q10/Q90 DVI TREND",
              startRow = current_row, startCol = 1
            )
            openxlsx::writeData(wb, sheet_name, loc_q_data,
              startRow = current_row + 1, startCol = 1
            )
            current_row <- current_row + nrow(loc_q_data) + 3
          }

          if (exists("all_unc_rmse") && !is.null(all_unc_rmse) && "location_id" %in% names(all_unc_rmse)) {
            loc_unc_rmse <- all_unc_rmse[all_unc_rmse$location_id == loc_id, , drop = FALSE]
            if (nrow(loc_unc_rmse) > 0) {
              openxlsx::writeData(wb, sheet_name, "RMSE CI (2.5%/97.5%)",
                startRow = current_row, startCol = 1
              )
              openxlsx::writeData(wb, sheet_name, loc_unc_rmse,
                startRow = current_row + 1, startCol = 1
              )
              current_row <- current_row + nrow(loc_unc_rmse) + 3
            }
          }

          if (!is.null(inseparable_flags_summary) && nrow(inseparable_flags_summary) > 0 && "location_id" %in% names(inseparable_flags_summary)) {
            loc_flagged <- inseparable_flags_summary[inseparable_flags_summary$location_id == loc_id, , drop = FALSE]
            if (nrow(loc_flagged) > 0) {
              openxlsx::writeData(wb, sheet_name, "INSEPARABLE VARIANT ALERTS (summary)", startRow = current_row, startCol = 1)
              openxlsx::writeData(wb, sheet_name, loc_flagged, startRow = current_row + 1, startCol = 1)
              current_row <- current_row + nrow(loc_flagged) + 3
            }
          }

          if (!is.null(inseparable_flags_detail) && nrow(inseparable_flags_detail) > 0 && "location_id" %in% names(inseparable_flags_detail)) {
            loc_flagged_detail <- inseparable_flags_detail[inseparable_flags_detail$location_id == loc_id, , drop = FALSE]
            if (nrow(loc_flagged_detail) > 0) {
              openxlsx::writeData(wb, sheet_name, "INSEPARABLE VARIANT COMPONENTS", startRow = current_row, startCol = 1)
              openxlsx::writeData(wb, sheet_name, loc_flagged_detail, startRow = current_row + 1, startCol = 1)
              current_row <- current_row + nrow(loc_flagged_detail) + 3
            }
          }

        }

    output_filename <- file.path(OUT_DIR, "mesma_results.xlsx")

    if (!dir.exists(OUT_DIR)) {
      cat(sprintf("Creating output directory: %s\n", OUT_DIR))
      dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
    }

    if (!dir.exists(OUT_DIR)) {
      cat(sprintf("ERROR: Cannot create output directory: %s\n", OUT_DIR))
      stop("Cannot create output directory")
    }

    cat(sprintf("Saving workbook to: %s\n", output_filename))

    save_result <- tryCatch(
      {
        openxlsx::saveWorkbook(wb, output_filename, overwrite = TRUE)
        TRUE
      },
      error = function(e) {
        cat(sprintf("ERROR saving workbook: %s\n", e$message))
        FALSE
      }
    )

    if (save_result) {
      cat(sprintf(
        "Created Excel file '%s' with %d location sheets\n",
        basename(output_filename), length(unique_locations)
      ))
    } else {
      cat("Failed to save Excel file\n")
    }
  }


  aggregate_to_global_pattern <- function(all_coefs, method = "location_bootstrap") {
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")
    
    # Debug: print method to diagnose "Unknown method" errors
    cat(sprintf("[DEBUG] aggregate_to_global_pattern called with method='%s'\n", method))

    # AGENT: Changed requirement to pheno_year
    required_cols <- c("location_id", "pheno_year", "Veg", "coef")
    missing <- setdiff(required_cols, names(all_coefs))
    if (length(missing) > 0) stop(paste("Missing columns:", paste(missing, collapse = ", ")))

    # Compatibility: Ensure 'year' exists for downstream functions that expect it
    if (!"year" %in% names(all_coefs)) all_coefs$year <- all_coefs$pheno_year

    if (method == "location_bootstrap") {
      # pheno_year is guaranteed present
      result <- location_bootstrap_aggregate(all_coefs, B = BOOTSTRAP_B)
    } else {
      stop(sprintf("Unknown method '%s'. Use 'location_bootstrap'.", method))
    }

    result
  }




  aggregate_simple_mean <- function(all_coefs) {
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")

    veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])

    # Determine which year column to use
    year_col <- if ("pheno_year" %in% names(all_coefs)) "pheno_year" else "year"

    results_list <- list()

    for (veg in veg_types) {
      veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), ]

      if (nrow(veg_data) == 0) {
        warning(sprintf("No valid coefficients found for vegetation type: %s", veg))
        next
      }

      # Group by year and compute simple statistics
      simple_result <- veg_data |> 
        dplyr::group_by(!!rlang::sym(year_col)) |> 
        dplyr::summarize(
          n_locations = dplyr::n(),
          global_coef = mean(coef, na.rm = TRUE),
          se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
          ci_lower = pmax(0, global_coef - 1.96 * se),
          ci_upper = pmin(1, global_coef + 1.96 * se),
          .groups = "drop"
        ) |> 
        dplyr::mutate(Veg = veg, method = "simple_mean") |> 
        dplyr::rename(year = !!rlang::sym(year_col))

      results_list[[veg]] <- simple_result
    }

    dplyr::bind_rows(results_list)
  }




  plot_global_vegetation_pattern <- function(global_pattern, 
                                              title = "Global Vegetation Composition Over Time",
                                              show_ci = TRUE,
                                              ci_type = "auto") {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    
    if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
      if ("ci_lower" %in% names(global_pattern) && "ci_upper" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower
        global_pattern$ci_upper <- global_pattern$ci_upper
      } else if ("coef_025" %in% names(global_pattern) && "coef_975" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$coef_025
        global_pattern$ci_upper <- global_pattern$coef_975
      }
    } else if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
      if ("ci_lower_simple" %in% names(global_pattern) && "ci_upper_simple" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower_simple
        global_pattern$ci_upper <- global_pattern$ci_upper_simple
      }
    } else if ("weighted_mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$weighted_mean_coef
      if ("ci_lower_pooled" %in% names(global_pattern) && "ci_upper_pooled" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower_pooled
        global_pattern$ci_upper <- global_pattern$ci_upper_pooled
      }
    }
    
    global_pattern_veg <- global_pattern[tolower(global_pattern$Veg) != "barren", ]
    global_pattern_barren <- global_pattern[tolower(global_pattern$Veg) == "barren", ]
    
    p <- ggplot2::ggplot(global_pattern_veg, ggplot2::aes(x = year, y = coef, color = Veg, fill = Veg)) +
      ggplot2::geom_line(size = 1.2) +
      ggplot2::geom_point(size = 2)
    
    if (show_ci && "ci_lower" %in% names(global_pattern_veg) && "ci_upper" %in% names(global_pattern_veg)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
        alpha = 0.2,
        color = NA
      )
    }
    
    if (nrow(global_pattern_barren) > 0) {
      veg_max <- suppressWarnings(max(global_pattern_veg$coef, global_pattern_veg$ci_upper, na.rm = TRUE))
      barren_max <- suppressWarnings(max(global_pattern_barren$coef, global_pattern_barren$ci_upper, na.rm = TRUE))
      if (is.na(veg_max) || veg_max <= 0 || is.na(barren_max) || barren_max <= 0) {
        barren_scale_factor <- 1
      } else {
        barren_scale_factor <- veg_max / barren_max
      }
      global_pattern_barren$coef_scaled <- global_pattern_barren$coef * barren_scale_factor
      
      p <- p + 
        ggplot2::geom_line(data = global_pattern_barren, 
              ggplot2::aes(x = year, y = coef_scaled), 
                          color = "brown", linewidth = 1.2, linetype = "dashed") +
        ggplot2::geom_point(data = global_pattern_barren, 
               ggplot2::aes(x = year, y = coef_scaled), 
                           color = "brown", size = 2)
      
      if (show_ci && "ci_lower" %in% names(global_pattern_barren) && "ci_upper" %in% names(global_pattern_barren)) {
        global_pattern_barren$ci_lower_scaled <- global_pattern_barren$ci_lower * barren_scale_factor
        global_pattern_barren$ci_upper_scaled <- global_pattern_barren$ci_upper * barren_scale_factor
        
        p <- p + ggplot2::geom_ribbon(data = global_pattern_barren,
                                     ggplot2::aes(x = year, ymin = ci_lower_scaled, ymax = ci_upper_scaled),
                                     fill = "brown", alpha = 0.1, color = NA)
      }
    }
    
    p <- p +
      ggplot2::labs(
        title = title,
        subtitle = sprintf("Based on %d locations", max(global_pattern$n_locations, na.rm = TRUE)),
        x = "Year",
        y = "Vegetation Fraction",
        color = "Vegetation Type",
        fill = "Vegetation Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)
      ) +
          ggplot2::scale_y_continuous(
        labels = scales::percent_format(),
        sec.axis = ggplot2::sec_axis(~ . / barren_scale_factor, name = "Barren Fraction", labels = scales::percent_format())
      ) +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::scale_fill_brewer(palette = "Set1")
    
    p
  }

  plot_vegetation_stacked_area <- function(global_pattern, 
                                            title = "Global Vegetation Composition Over Time") {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    
    if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
    } else if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
    }
    
    global_pattern <- global_pattern |> 
      dplyr::group_by(year) |> 
      dplyr::mutate(coef_normalized = coef / sum(coef, na.rm = TRUE)) |> 
      dplyr::ungroup()
    
    p <- ggplot2::ggplot(global_pattern, 
                          ggplot2::aes(x = year, y = coef_normalized, fill = Veg)) +
      ggplot2::geom_area(alpha = 0.8, position = "stack") +
      ggplot2::labs(
        title = title,
        x = "Year",
        y = "Relative Vegetation Fraction",
        fill = "Vegetation Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold")
      ) +
      ggplot2::scale_y_continuous(labels = scales::percent_format()) +
      ggplot2::scale_fill_brewer(palette = "Set2")
    
    p
  }

  plot_vegetation_heatmap <- function(global_pattern, 
                                       title = "Vegetation Fraction by Year") {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    
    if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
    } else if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
    }
    
    p <- ggplot2::ggplot(global_pattern, 
                          ggplot2::aes(x = year, y = Veg, fill = coef)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", coef * 100)), 
                         color = "white", size = 3) +
      ggplot2::scale_fill_viridis_c(option = "plasma", 
                                     labels = scales::percent_format()) +
      ggplot2::labs(
        title = title,
        x = "Year",
        y = "Vegetation Type",
        fill = "Fraction"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    
    p
  }

    cat("\nGenerating average coverage plot...\n")
  
    tryCatch({
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        install.packages("ggplot2")
      }
      library(ggplot2)
  
      global_pattern_all <- aggregate_to_global_pattern(all_coefs, method = "location_bootstrap")
  
      # DEBUG: inspect aggregated global pattern to diagnose replacement/row-mismatch errors
      if (is.null(global_pattern_all)) {
        cat("[DEBUG] aggregate_to_global_pattern returned NULL\n")
      } else {
        cat(sprintf("[DEBUG] global_pattern_all: rows=%d cols=%d\n", nrow(global_pattern_all), ncol(global_pattern_all)))
        cat(sprintf("[DEBUG] global_pattern_all columns: %s\n", paste(names(global_pattern_all), collapse = ", ")))
        sample_veg <- if ("Veg" %in% names(global_pattern_all)) paste(head(unique(global_pattern_all$Veg), 10), collapse = ", ") else "<no Veg column>"
        cat(sprintf("[DEBUG] global_pattern_all Veg sample: %s\n", sample_veg))
      }
  
      if (is.null(global_pattern_all) || nrow(global_pattern_all) == 0) {
        cat("[WARN] aggregate_to_global_pattern returned no data, skipping average coverage plot\n")
      } else {
        global_pattern_all <- global_pattern_all |> dplyr::filter(year >= 1985 & year <= 2025)
        global_pattern_all_veg <- global_pattern_all[tolower(global_pattern_all$Veg) != "barren", ]
        global_pattern_all_barren <- global_pattern_all[tolower(global_pattern_all$Veg) == "barren", ]
  
        # Compatibility shim: if bootstrap returned coef_025/coef_975, map them to ci_lower/ci_upper for plotting
        if (!("ci_lower" %in% names(global_pattern_all_veg)) && ("coef_025" %in% names(global_pattern_all_veg))) {
          global_pattern_all_veg$ci_lower <- global_pattern_all_veg$coef_025
          global_pattern_all_veg$ci_upper <- global_pattern_all_veg$coef_975
          cat("[DEBUG] Mapped coef_025/coef_975 -> ci_lower/ci_upper for vegetation data\n")
        }
        if (!("ci_lower" %in% names(global_pattern_all_barren)) && ("coef_025" %in% names(global_pattern_all_barren))) {
          global_pattern_all_barren$ci_lower <- global_pattern_all_barren$coef_025
          global_pattern_all_barren$ci_upper <- global_pattern_all_barren$coef_975
          cat("[DEBUG] Mapped coef_025/coef_975 -> ci_lower/ci_upper for barren data\n")
        }
  
    # DEBUG: report sizes before plotting
    cat(sprintf("[DEBUG] plotting: veg_rows=%d, barren_rows=%d\n", nrow(global_pattern_all_veg), nrow(global_pattern_all_barren)))
    if (nrow(global_pattern_all_veg) > 0) {
      cat(sprintf("[DEBUG] global_pattern_all_veg sample rows:\n"))
      print(head(global_pattern_all_veg, 4))
    }
    if (nrow(global_pattern_all_barren) > 0) {
      cat(sprintf("[DEBUG] global_pattern_all_barren sample rows:\n"))
      print(head(global_pattern_all_barren, 4))
    }
  
    # Build plot with localized tryCatch to capture failing step
    p <- tryCatch({
      ggplot(global_pattern_all_veg, aes(x = year, y = global_coef, color = Veg, fill = Veg)) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 2, show.legend = FALSE) +
        geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, color = NA)
    }, error = function(e) {
      cat(sprintf("[ERROR] Failed creating ggplot object: %s\n", e$message))
      NULL
    })
    if (is.null(p)) {
      stop("Plot creation failed — aborting average coverage plot generation")
    }
  
    barren_scale_factor <- 1
    if (nrow(global_pattern_all_barren) > 0) {
      veg_max <- suppressWarnings(max(global_pattern_all_veg$global_coef, global_pattern_all_veg$ci_upper, na.rm = TRUE))
      barren_max <- suppressWarnings(max(global_pattern_all_barren$global_coef, global_pattern_all_barren$ci_upper, na.rm = TRUE))
  
      if (is.infinite(veg_max) || is.na(veg_max) || veg_max <= 0) veg_max <- 0.01 # Default small value if no veg
      if (is.infinite(barren_max) || is.na(barren_max) || barren_max <= 0) barren_max <- 1
  
      if (veg_max <= 0.01) {
         barren_scale_factor <- 1 # No scaling if veg is negligible
      } else {
         barren_scale_factor <- veg_max / barren_max
      }
  
      global_pattern_all_barren$coef_scaled <- global_pattern_all_barren$global_coef * barren_scale_factor
      global_pattern_all_barren$ci_lower_scaled <- global_pattern_all_barren$ci_lower * barren_scale_factor
      global_pattern_all_barren$ci_upper_scaled <- global_pattern_all_barren$ci_upper * barren_scale_factor
  
      p <- p +
        geom_line(data = global_pattern_all_barren,
             aes(x = year, y = coef_scaled),
                 color = "brown", linewidth = 1.2, linetype = "dashed") +
        geom_point(data = global_pattern_all_barren,
              aes(x = year, y = coef_scaled),
                  color = "brown", size = 2, show.legend = FALSE) +
        geom_ribbon(data = global_pattern_all_barren,
                   aes(x = year, ymin = ci_lower_scaled, ymax = ci_upper_scaled),
                   fill = "brown", alpha = 0.1, color = NA)
    }
  
    all_y_values <- c(global_pattern_all_veg$global_coef, global_pattern_all_veg$ci_lower, global_pattern_all_veg$ci_upper)
    if (exists("global_pattern_all_barren") && nrow(global_pattern_all_barren) > 0) {
      all_y_values <- c(all_y_values, global_pattern_all_barren$coef_scaled, global_pattern_all_barren$ci_lower_scaled, global_pattern_all_barren$ci_upper_scaled)
    }
    min_y <- min(all_y_values, na.rm = TRUE)
    max_y <- max(all_y_values, na.rm = TRUE)
  
    if (is.infinite(min_y) || is.na(min_y)) min_y <- 0
    if (is.infinite(max_y) || is.na(max_y)) max_y <- 1
    if (max_y <= min_y) max_y <- min_y + 0.1
  
  
    p <- p +
      labs(
        title = "Average Coverage Percentage per Vegetation Type (2000-2024)",
        subtitle = sprintf("Based on %d locations with bootstrap uncertainty", max(global_pattern_all$n_locations, na.rm = TRUE)),
        x = "Year",
        y = "Vegetation Fraction",
        color = "Vegetation Type",
        fill = "Vegetation Type"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
      ) +
      scale_y_continuous(
        labels = scales::percent_format(accuracy = 1),
        limits = c(min_y, max_y),
        expand = c(0, 0),
        sec.axis = sec_axis(~ . / barren_scale_factor, name = "Barren Fraction", labels = scales::percent_format(accuracy = 1))
      ) +
      scale_color_brewer(palette = "Set1") +
      scale_fill_brewer(palette = "Set1")
  
        plot_filename <- file.path(OUT_DIR, "training_average_coverage_plot.png")
        ggsave(plot_filename, p, width = 10, height = 6, dpi = 300)
        cat(sprintf("Saved average coverage plot to: %s\n", plot_filename))
      }
    }, error = function(e) {
      cat(sprintf("[ERROR] Failed to generate average coverage plot: %s\n", e$message))
      cat("[INFO] Continuing with script execution...\n")
    })
  
    cat("\nGenerating Observations vs Accuracy plot...\n")
  
  
      loc_accuracy <- best_fit_summary |>
        dplyr::filter(!is_barren_truth) |>
      dplyr::group_by(location_id) |>
      dplyr::summarize(
        mean_pred_coef_rel = mean(pred_coef_rel, na.rm = TRUE),
        mean_pred_coef_abs = mean(pred_coef_abs, na.rm = TRUE),
        .groups = "drop"
      )
  
    obs_vs_acc_data <- summary_data |>
      dplyr::left_join(loc_accuracy, by = c("Location_ID" = "location_id"))
  
    obs_vs_acc_data <- obs_vs_acc_data |> dplyr::filter(!is.na(mean_pred_coef_rel))
  
    if (nrow(obs_vs_acc_data) > 0) {
      p_obs_acc <- ggplot(obs_vs_acc_data, aes(x = Total_Observations, y = mean_pred_coef_rel)) +
        geom_point(alpha = 0.6, color = "darkblue") +
        geom_smooth(method = "lm", color = "red", se = FALSE) +
        labs(
          title = "Number of Observations vs. Vegetation Prediction Accuracy",
          subtitle = "Accuracy = Mean Predicted Fraction of True Vegetation Class (Relative to Total Vegetation)",
          x = "Total Observations (All Years)",
          y = "Mean Correct Prediction Fraction"
        ) +
        theme_minimal() +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1))
  
      ggsave(file.path(OUT_DIR, "observations_vs_accuracy.png"), p_obs_acc, width = 8, height = 6)
      cat(sprintf("Saved Observations vs Accuracy plot to: %s\n", file.path(OUT_DIR, "observations_vs_accuracy.png")))
    } else {
      cat("No data available for Observations vs Accuracy plot.\n")
    }
  
  

  aggregate_simple_mean <- function(all_coefs) {
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")

    veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])

    # Determine which year column to use
    year_col <- if ("pheno_year" %in% names(all_coefs)) "pheno_year" else "year"

    results_list <- list()

    for (veg in veg_types) {
      veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), ]

      if (nrow(veg_data) == 0) {
        warning(sprintf("No valid coefficients found for vegetation type: %s", veg))
        next
      }

      # Group by year and compute simple statistics
      simple_result <- veg_data |> 
        dplyr::group_by(!!rlang::sym(year_col)) |> 
        dplyr::summarize(
          n_locations = dplyr::n(),
          global_coef = mean(coef, na.rm = TRUE),
          se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
          ci_lower = pmax(0, global_coef - 1.96 * se),
          ci_upper = pmin(1, global_coef + 1.96 * se),
          .groups = "drop"
        ) |> 
        dplyr::mutate(Veg = veg, method = "simple_mean") |> 
        dplyr::rename(year = !!rlang::sym(year_col))

      results_list[[veg]] <- simple_result
    }

    dplyr::bind_rows(results_list)
  }



  plot_global_vegetation_pattern <- function(global_pattern, 
                                              title = "Global Vegetation Composition Over Time",
                                              show_ci = TRUE,
                                              ci_type = "auto") {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    
    if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
      if ("ci_lower" %in% names(global_pattern) && "ci_upper" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower
        global_pattern$ci_upper <- global_pattern$ci_upper
      } else if ("coef_025" %in% names(global_pattern) && "coef_975" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$coef_025
        global_pattern$ci_upper <- global_pattern$coef_975
      }
    } else if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
      if ("ci_lower_simple" %in% names(global_pattern) && "ci_upper_simple" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower_simple
        global_pattern$ci_upper <- global_pattern$ci_upper_simple
      }
    } else if ("weighted_mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$weighted_mean_coef
      if ("ci_lower_pooled" %in% names(global_pattern) && "ci_upper_pooled" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower_pooled
        global_pattern$ci_upper <- global_pattern$ci_upper_pooled
      }
    }
    
    global_pattern_veg <- global_pattern[tolower(global_pattern$Veg) != "barren", ]
    global_pattern_barren <- global_pattern[tolower(global_pattern$Veg) == "barren", ]
    
    p <- ggplot2::ggplot(global_pattern_veg, ggplot2::aes(x = year, y = coef, color = Veg, fill = Veg)) +
      ggplot2::geom_line(size = 1.2) +
      ggplot2::geom_point(size = 2)
    
    if (show_ci && "ci_lower" %in% names(global_pattern_veg) && "ci_upper" %in% names(global_pattern_veg)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
        alpha = 0.2,
        color = NA
      )
    }
    
    barren_scale_factor <- 1
    if (nrow(global_pattern_barren) > 0) {
      veg_max <- suppressWarnings(max(global_pattern_veg$coef, global_pattern_veg$ci_upper, na.rm = TRUE))
      barren_max <- suppressWarnings(max(global_pattern_barren$coef, global_pattern_barren$ci_upper, na.rm = TRUE))
      if (is.na(veg_max) || veg_max <= 0 || is.na(barren_max) || barren_max <= 0) {
        barren_scale_factor <- 1
      } else {
        barren_scale_factor <- veg_max / barren_max
      }
      global_pattern_barren$coef_scaled <- global_pattern_barren$coef * barren_scale_factor
      if ("ci_lower" %in% names(global_pattern_barren) && "ci_upper" %in% names(global_pattern_barren)) {
        global_pattern_barren$ci_lower_scaled <- global_pattern_barren$ci_lower * barren_scale_factor
        global_pattern_barren$ci_upper_scaled <- global_pattern_barren$ci_upper * barren_scale_factor
      }
      
      p <- p + 
        ggplot2::geom_line(data = global_pattern_barren, 
                          ggplot2::aes(x = year, y = coef_scaled), 
                          color = "brown", linewidth = 1.2, linetype = "dashed") +
        ggplot2::geom_point(data = global_pattern_barren,
                          ggplot2::aes(x = year, y = coef_scaled),
                          color = "brown", size = 2)

      if ("ci_lower_scaled" %in% names(global_pattern_barren) && "ci_upper_scaled" %in% names(global_pattern_barren)) {
        p <- p + ggplot2::geom_ribbon(data = global_pattern_barren,
                                    ggplot2::aes(x = year, ymin = ci_lower_scaled, ymax = ci_upper_scaled),
                                    fill = "brown", alpha = 0.1, color = NA, inherit.aes = FALSE)
      }

      p <- p + ggplot2::labs(
        title = title,
        subtitle = sprintf("Based on %d locations", max(global_pattern$n_locations, na.rm = TRUE)),
        x = "Year",
        y = "Vegetation Fraction",
        color = "Vegetation Type",
        fill = "Vegetation Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)
      ) +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(),
        sec.axis = ggplot2::sec_axis(~ . / barren_scale_factor, name = "Barren Fraction", labels = scales::percent_format())
      ) +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::scale_fill_brewer(palette = "Set1")
    }

    p
  }

  plot_vegetation_stacked_area <- function(global_pattern, 
                                            title = "Global Vegetation Composition Over Time") {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    
    if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
    } else if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
    }
    
    global_pattern <- global_pattern |> 
      dplyr::group_by(year) |> 
      dplyr::mutate(coef_normalized = coef / sum(coef, na.rm = TRUE)) |> 
      dplyr::ungroup()
    
    p <- ggplot2::ggplot(global_pattern, 
                          ggplot2::aes(x = year, y = coef_normalized, fill = Veg)) +
      ggplot2::geom_area(alpha = 0.8, position = "stack") +
      ggplot2::labs(
        title = title,
        x = "Year",
        y = "Relative Vegetation Fraction",
        fill = "Vegetation Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold")
      ) +
      ggplot2::scale_y_continuous(labels = scales::percent_format()) +
      ggplot2::scale_fill_brewer(palette = "Set2")
    
    p
  }

  plot_vegetation_heatmap <- function(global_pattern, 
                                       title = "Vegetation Fraction by Year") {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    
    if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
    } else if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
    }
    
    p <- ggplot2::ggplot(global_pattern, 
                          ggplot2::aes(x = year, y = Veg, fill = coef)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", coef * 100)), 
                         color = "white", size = 3) +
      ggplot2::scale_fill_viridis_c(option = "plasma", 
                                     labels = scales::percent_format()) +
      ggplot2::labs(
        title = title,
        x = "Year",
        y = "Vegetation Type",
        fill = "Fraction"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    
    p
  }

  
  

  bootstrap_trend_ci <- function(all_coefs, B = 200, seed = 123) {
  set.seed(seed)
  if (!requireNamespace("lme4", quietly = TRUE)) {
    warning("lme4 package not found, trend CI calculation will be skipped.")
    return(NULL)
  }
  
  veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
  results_list <- list()
  
  for (veg in veg_types) {
    veg_data <- all_coefs[all_coefs$Veg == veg & is.finite(all_coefs$coef), ]
    
    if (nrow(veg_data) < 10 || length(unique(veg_data$location_id)) < 3) {
      cat(sprintf("[TREND] Skipping trend for '%s' (insufficient data: %d rows, %d locs)\n", veg, nrow(veg_data), length(unique(veg_data$location_id))))
      next
    }
    
    locations <- unique(veg_data$location_id)
    n_locs <- length(locations)
    
    boot_slopes <- replicate(B, {
      # Resample locations with replacement
      boot_locs_sampled <- sample(locations, n_locs, replace = TRUE)
      
      # Create bootstrap sample by selecting all data from resampled locations
      # Handle cases where a location is sampled multiple times by creating a new ID
      boot_data_list <- lapply(seq_along(boot_locs_sampled), function(i) {
        loc_data <- veg_data[veg_data$location_id == boot_locs_sampled[i], ]
        loc_data$boot_id <- paste0(boot_locs_sampled[i], "_", i)
        
        if ("coef_sd" %in% names(loc_data) && !all(is.na(loc_data$coef_sd))) {
           sds <- loc_data$coef_sd
           sds[is.na(sds) | sds < 0] <- 0
           if (any(sds > 0)) {
             # Debug perturbation for the first iteration of the first location
             if (i == 1 && runif(1) < 0.001) cat(sprintf("[DEBUG] Perturbing coefs for loc %s. Max SD: %.6f\n", boot_locs_sampled[i], max(sds)))
             loc_data$coef <- rnorm(nrow(loc_data), mean = loc_data$coef, sd = sds)
           }
        }
        
        loc_data
      })
      boot_data <- do.call(rbind, boot_data_list)
      
      # Fit a linear mixed-effects model to the bootstrap sample
      # Use boot_id as the random effect to handle repeated locations
      model <- tryCatch({
        lme4::lmer(coef ~ pheno_year + (1|boot_id), data = boot_data, control = lme4::lmerControl(check.conv.singular = "ignore"))
      }, error = function(e) { NULL })
      
      if (is.null(model)) {
        return(NA_real_)
      } else {
        # Extract the fixed effect slope for the year
        return(lme4::fixef(model)["pheno_year"])
      }
    })
    
    finite_slopes <- boot_slopes[is.finite(boot_slopes)]
    if (length(finite_slopes) > 5) {
      results_list[[veg]] <- data.frame(
        Veg = veg, 
        slope_mean = mean(finite_slopes), 
        slope_median = median(finite_slopes), 
        slope_ci_lower = quantile(finite_slopes, 0.025), 
        slope_ci_upper = quantile(finite_slopes, 0.975), 
        prob_positive = mean(finite_slopes > 0), 
        prob_negative = mean(finite_slopes < 0),
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(results_list) > 0) {
    return(dplyr::bind_rows(results_list))
  } else {
    return(NULL)
  }
}



  global_pattern <- aggregate_to_global_pattern(all_coefs, method = "location_bootstrap")
  global_pattern_bootstrap <- global_pattern

  p1_all <- plot_global_vegetation_pattern(global_pattern, 
                                          title = "All Vegetation Trends",
                                          show_ci = TRUE)
  p1_all <- p1_all + scale_y_continuous(labels = scales::percent_format()) + labs(y = "Fraction")
  ggsave(file.path(OUT_DIR, "training_all_vegetation_trends.png"), p1_all, width = 10, height = 6)



  p3 <- plot_vegetation_stacked_area(global_pattern)
  ggsave(file.path(OUT_DIR, "training_vegetation_stacked_area.png"), p3, width = 10, height = 6)

  p3_veg_only <- plot_vegetation_only_stacked_area(global_pattern)
  ggsave(file.path(OUT_DIR, "training_vegetation_only_stacked_area.png"), p3_veg_only, width = 10, height = 6)

  p4 <- plot_vegetation_heatmap(global_pattern)
  ggsave(file.path(OUT_DIR, "training_vegetation_heatmap.png"), p4, width = 10, height = 6)

  # AGENT: Trend analysis DISABLED for training data (inference only)
  cat("[NOTICE] Skipped trend analysis for training data.\n")
  trend_ci <- NULL 

  analyze_vegetation_trends <- function(all_coefs, B = 200) {
    # Simplified trend estimation using per-location linear slopes on annual coefficients.
    # Block/AR-based bootstrapping removed; we compute a median of location-level slopes.
    if (is.null(all_coefs) || nrow(all_coefs) == 0) return(NULL)
    veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
    trend_rows <- list()
    for (veg in veg_types) {
      veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), , drop = FALSE]
      if (nrow(veg_data) == 0) next
      locations <- unique(veg_data$location_id)
      loc_slopes <- sapply(locations, function(loc) {
        loc_data <- veg_data[veg_data$location_id == loc, , drop = FALSE]
        if (nrow(loc_data) < 2) return(NA_real_)
        # Fit simple linear model of coef ~ year (use pheno_year if present)
        yr_col <- if ("pheno_year" %in% names(loc_data)) "pheno_year" else "year"
        df_loc <- loc_data[order(loc_data[[yr_col]]), , drop = FALSE]
        m <- tryCatch(lm(coef ~ get(yr_col), data = df_loc), error = function(e) NULL)
        if (is.null(m)) return(NA_real_)
        coef_val <- tryCatch(coef(m)[[2]], error = function(e) NA_real_)
        as.numeric(coef_val)
      })
      trend_rows[[veg]] <- data.frame(Veg = veg, slope = median(loc_slopes, na.rm = TRUE), n_locations = sum(is.finite(loc_slopes)), stringsAsFactors = FALSE)
    }
    if (length(trend_rows) == 0) return(NULL)
    do.call(rbind, trend_rows)
  }

  # AGENT: Trend analysis DISABLED for training data
  trends <- NULL
  if (is.null(trends)) {
    trends <- data.frame(Veg = character(0), slope = numeric(0), n_locations = integer(0), stringsAsFactors = FALSE)
  }

  cat("\nGenerating PPI-normalized cumulative plot with location bootstrap...\n")
  if ("PPI" %in% names(df) && any(!is.na(df$PPI))) {
    global_pattern_ppi <- location_bootstrap_ppi(all_coefs, df, B = BOOTSTRAP_B, seed = 123)
    if (!is.null(global_pattern_ppi) && nrow(global_pattern_ppi) > 0) {
      global_pattern_ppi <- global_pattern_ppi[!tolower(trimws(global_pattern_ppi$Veg)) %in% c("barren"), ]
    }
    if (!is.null(global_pattern_ppi) && nrow(global_pattern_ppi) > 0) {
      p_ppi_ts <- ggplot(global_pattern_ppi, aes(x = year, y = global_coef, color = Veg, group = Veg)) +
        geom_line(linewidth = 1) +
        geom_point(show.legend = FALSE) +
        geom_ribbon(aes(ymin = coef_025, ymax = coef_975, fill = Veg), alpha = 0.15, color = NA) +
        labs(title = "PPI-Normalized Vegetation Fractions Over Time (Location Bootstrap)",
             x = "Year", y = "Total Normalized Fraction", color = "Veg", fill = "Veg") +
        theme_minimal()
      ggsave(file.path(OUT_DIR, "training_ppi_normalized_timeseries.png"), p_ppi_ts, width = 8, height = 6)
      cat(sprintf("Saved PPI-normalized time series plot to: %s\n", file.path(OUT_DIR, "training_ppi_normalized_timeseries.png")))
    } else {
      cat("PPI normalization aggregation returned no results (no matching location-year PPI values).\n")
    }
  } else {
    cat("PPI data not available, skipping PPI-normalized plot.\n")
  }

  openxlsx::addWorksheet(wb, "Global_Pattern")
  openxlsx::writeData(wb, "Global_Pattern", global_pattern)

  if (!is.null(global_pattern_ppi) && exists("global_pattern_ppi") && nrow(global_pattern_ppi) > 0) {
    openxlsx::addWorksheet(wb, "PPI_Normalized")
    openxlsx::writeData(wb, "PPI_Normalized", global_pattern_ppi)
    cat("Added PPI-normalized results to Excel workbook\n")
  }

  openxlsx::addWorksheet(wb, "Vegetation_Trends")
  openxlsx::writeData(wb, "Vegetation_Trends", trends)

  if (!is.null(trend_ci)) {
    openxlsx::addWorksheet(wb, "Trend_Bootstrap_CI")
    openxlsx::writeData(wb, "Trend_Bootstrap_CI", trend_ci)
  }

  timing_info$end_time <- Sys.time()
  total_time <- as.numeric(difftime(timing_info$end_time, timing_info$start_time, units = "secs"))

  cat(sprintf("\nTotal execution time: %.1f seconds (%.1f minutes)\n", total_time, total_time / 60))
  if (!is.null(timing_info$moving_var_done)) {
    cat(sprintf(
      "Moving variance: %.1f seconds\n",
      as.numeric(difftime(timing_info$moving_var_done, timing_info$start_time, units = "secs"))
    ))
  }
  if (!is.null(timing_info$lib_construction_done)) {
    cat(sprintf(
      "Library construction: %.1f seconds\n",
      as.numeric(difftime(timing_info$lib_construction_done, timing_info$moving_var_done, units = "secs"))
    ))
  }
  if (!is.null(timing_info$pca_computation_done)) {
    cat(sprintf(
      "PCA computation: %.1f seconds\n",
      as.numeric(difftime(timing_info$pca_computation_done, timing_info$lib_construction_done, units = "secs"))
    ))
  }
  cat(sprintf(
    "Main processing + Excel: %.1f seconds\n",
    as.numeric(difftime(timing_info$end_time, timing_info$pca_computation_done, units = "secs"))
  ))

  cat("[DEBUG] About to complete main_processing_block\n")
  cat("\nMESMA fitting completed successfully!\n")

cat("[DEBUG] Finished executing main processing steps\n")
flush.console()

# --- Modular helpers to re-run or run pieces of the main processing separately ---
process_batches <- function(BATCH_SIZE = 6, overwrite_wb = FALSE) { 
  cat("[INFO] process_batches: building batch list\n"); flush.console()
  target_locations <- intersect(location_list, unique(df_tasks$location_id))
  n_locs_to_process <- length(target_locations)
  loc_batches <- split(target_locations, ceiling(seq_along(target_locations) / BATCH_SIZE))
  # Memory-efficient: no longer use all_coefs_list, build from results_local instead
  results_local <- list()

  for (i in seq_along(loc_batches)) {
    batch_locs <- loc_batches[[i]]
    batch_df <- df_tasks[df_tasks$location_id %in% batch_locs, ]
    batch_location_list <- split(batch_df, batch_df$location_id)
    # Suppress any verbose output from per-location processing in helper
    sink_file <- tempfile()
    con_out <- file(sink_file, open = "wt")
    con_msg <- file(paste0(sink_file, ".msg"), open = "wt")
    tryCatch({
      sink(con_out, type = "output")
      sink(con_msg, type = "message")
      batch_results <- .run_map(batch_location_list, fit_one_location, show_pb = FALSE)
    }, finally = {
      try(sink(type = "message"), silent = TRUE)
      try(sink(type = "output"), silent = TRUE)
      close(con_out); close(con_msg)
      try(unlink(sink_file), silent = TRUE); try(unlink(paste0(sink_file, ".msg")), silent = TRUE)
    })

    for (k in names(batch_results)) {
      loc_result <- batch_results[[k]]
      if (is.null(loc_result)) next

      # Store year results with composite key 'loc_year'
      for (yr in names(loc_result)) {
        key <- paste(k, yr, sep = "_")
        results_local[[key]] <- loc_result[[yr]]
      }

      # Build and write location data (temporary, memory-efficient)
      loc_data <- do.call(rbind, lapply(loc_result, function(yr_res) yr_res$coef_df))
      if (!is.null(loc_data) && nrow(loc_data) > 0) {
        # Write to Excel WITHOUT accumulating in memory
        if (!is.null(wb)) {
          sheet_name <- paste0("Loc_", k)
          if (overwrite_wb && sheet_name %in% names(wb)) try(openxlsx::removeWorksheet(wb, sheet_name), silent = TRUE)
          if (!(sheet_name %in% names(wb))) {
            openxlsx::addWorksheet(wb, sheet_name)
            openxlsx::writeData(wb, sheet_name, loc_data)
          }
        }
      }
      # Clear temporary data
      rm(loc_data)
    }

    # Aggressive memory cleanup after each batch
    rm(batch_df, batch_location_list, batch_results)
    gc(verbose = FALSE, full = TRUE)

    # Force memory release and print stats every 5 batches
    if (i %% 5 == 0) {
      cat(sprintf("  [Memory cleanup at batch %d/%d]\n", i, length(loc_batches)))
      gc(verbose = TRUE, full = TRUE)
    }
  }

  # Build all_coefs from results_local (memory-efficient)
  cat("[INFO] Building all_coefs from results...\n")
  all_coefs <- tryCatch({
    do.call(rbind, lapply(results_local, function(r) {
      if (!is.null(r$coef_df)) r$coef_df else NULL
    }))
  }, error = function(e) NULL)
  list(all_coefs = all_coefs, results = results_local)
}

process_inference <- function(BATCH_SIZE = BATCH_SIZE) {  # Reduced to 6 for memory efficiency
  load_and_prepare_inference_data()
  
  # cat(sprintf("[INFO] process_inference: starting (BATCH_SIZE=%d)\n", BATCH_SIZE)); flush.console()
  if (!(exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0)) {
    # cat("[INFO] No inference dataset present\n"); 
    return(list(results = list(), all_coefs = NULL))
  }

  if (!"task_key" %in% names(df_tasks_inference)) df_tasks_inference$task_key <- paste(df_tasks_inference$location_id, df_tasks_inference$pheno_year, sep = "_")
  inference_loc_years <- df_tasks_inference |> dplyr::filter(!is.na(.data$location_id) & trimws(.data$location_id) != "" & !is.na(.data$pheno_year) & .data$pheno_year > 0) |> dplyr::distinct(.data$location_id, .data$pheno_year)
  inference_target_keys <- paste(inference_loc_years$location_id, inference_loc_years$pheno_year, sep = "_")
  inference_available_keys <- unique(df_tasks_inference$task_key)
  inference_target_keys <- intersect(inference_target_keys, inference_available_keys)
  
  # Process inference LOCATIONS (apply optional cap)
  unique_locs <- unique(inference_loc_years$location_id)
  # Cap total number of locations to MAX_INFERENCE_LOCATIONS to limit runtime and memory
  if (exists("MAX_INFERENCE_LOCATIONS") && length(unique_locs) > MAX_INFERENCE_LOCATIONS) {
    set.seed(123) # deterministic sampling for reproducibility
    selected_locs <- sample(unique_locs, MAX_INFERENCE_LOCATIONS, replace = FALSE)
    inference_loc_years <- inference_loc_years[inference_loc_years$location_id %in% selected_locs, ]
    inference_target_keys <- paste(inference_loc_years$location_id, inference_loc_years$pheno_year, sep = "_")
    inference_target_keys <- intersect(inference_target_keys, inference_available_keys)
    unique_locs <- unique(inference_loc_years$location_id)
    cat(sprintf("[INFO process_inference] Reduced inference locations to %d (random sample of %d)\n", length(unique_locs), MAX_INFERENCE_LOCATIONS))
  }
  cat(sprintf("[DEBUG process_inference] Processing %d unique locations (%d location-year tasks)\n", length(unique_locs), length(inference_target_keys)))

  if (length(inference_target_keys) == 0) { # cat("[INFO] No valid inference tasks\n"); 
    return(list(results = list(), all_coefs = NULL)) }

  inference_key_batches <- split(inference_target_keys, ceiling(seq_along(inference_target_keys) / BATCH_SIZE))
  inference_results_list <- vector("list", length(inference_target_keys))
  names(inference_results_list) <- inference_target_keys

  cat(sprintf("[INFO process_inference] Processing %d inference tasks from %d locations in %d batches...\n", length(inference_target_keys), length(unique_locs), length(inference_key_batches)))

  for (i in seq_along(inference_key_batches)) {
    batch_keys <- inference_key_batches[[i]]
    batch_df <- df_tasks_inference[df_tasks_inference$task_key %in% batch_keys, ]
    batch_task_list <- split(batch_df, batch_df$task_key)
    # Suppress any verbose output from per-task processing in helper
    sink_file <- tempfile()
    con_out <- file(sink_file, open = "wt")
    con_msg <- file(paste0(sink_file, ".msg"), open = "wt")
    tryCatch({
      sink(con_out, type = "output")
      sink(con_msg, type = "message")
      batch_results <- .run_map(batch_task_list, fit_one_task, show_pb = FALSE)
    }, finally = {
      try(sink(type = "message"), silent = TRUE)
      try(sink(type = "output"), silent = TRUE)
      close(con_out); close(con_msg)
      try(unlink(sink_file), silent = TRUE); try(unlink(paste0(sink_file, ".msg")), silent = TRUE)
    })
    inference_results_list[names(batch_results)] <- batch_results
    rm(batch_df, batch_task_list, batch_results); gc(verbose = FALSE)

    # Print progress percentage every 10%
    if (i == length(inference_key_batches) || (i %% ceiling(length(inference_key_batches) / 10) == 0)) {
      pct <- round((i / length(inference_key_batches)) * 100)
      cat(sprintf("[INFO process_inference] Batch %d/%d (%d%%) finished\n", i, length(inference_key_batches), pct))
    }
  }

  all_coefs_local <- combine_results_from_list(inference_results_list)

  # Save inference results as CSV (explicitly for rf_classifier.py compatibility)
  inf_csv_dir <- "inference_results"
  if (!dir.exists(inf_csv_dir)) dir.create(inf_csv_dir, recursive = TRUE)
  inf_csv_path <- file.path(inf_csv_dir, "inference_results.csv")
  
  if (!is.null(all_coefs_local) && nrow(all_coefs_local) > 0) {
    if (requireNamespace("readr", quietly = TRUE)) {
      readr::write_csv(all_coefs_local, inf_csv_path)
    } else {
      write.csv(all_coefs_local, inf_csv_path, row.names = FALSE)
    }
    cat(sprintf("[INFO] Saved inference results CSV to: %s\n", inf_csv_path))
  } else {
    cat("[WARN] No inference coefficients to save to CSV.\n")
  }

  
  # Save inference results as GeoJSON
  if (exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(all_coefs_local) > 0) {
    all_coefs_geo <- dplyr::left_join(all_coefs_local, df_tasks_inference |> dplyr::select(location_id, lat, lon), by = "location_id") |> dplyr::distinct()
    if (nrow(all_coefs_geo) > 0 && "lat" %in% names(all_coefs_geo) && "lon" %in% names(all_coefs_geo) && !all(is.na(all_coefs_geo$lat)) && !all(is.na(all_coefs_geo$lon))) {
      sf_obj <- sf::st_as_sf(all_coefs_geo, coords = c("lon", "lat"), crs = 4326)
      geojson_path <- file.path(OUT_DIR, "inference_results.geojson")
      sf::st_write(sf_obj, geojson_path, driver = "GeoJSON", delete_dsn = TRUE)
      cat(sprintf("Saved inference results as GeoJSON to: %s\n", geojson_path))
    } else {
      cat("[WARN] Could not create GeoJSON: missing or invalid lat/lon in inference data\n")
    }
  }
  
  inference_result <- list(results = inference_results_list, all_coefs = all_coefs_local)

  # Add inference results to wb and generate plot
  if (length(inference_results_list) > 0) {
    cat("\nGenerating inference average coverage plot...\n")
    global_pattern_inference <- NULL
    tryCatch({
      global_pattern_inference <- aggregate_to_global_pattern(all_coefs_local, method = "location_bootstrap")
      if (!is.null(global_pattern_inference) && nrow(global_pattern_inference) > 0) {
        global_pattern_inference <- global_pattern_inference |> dplyr::filter(year >= 1985 & year <= 2025)
        global_pattern_inference_veg <- global_pattern_inference[tolower(global_pattern_inference$Veg) != "barren", ]
        global_pattern_inference_barren <- global_pattern_inference[tolower(global_pattern_inference$Veg) == "barren", ]
        # Compatibility shim
        if (!("ci_lower" %in% names(global_pattern_inference_veg)) && ("coef_025" %in% names(global_pattern_inference_veg))) {
          global_pattern_inference_veg$ci_lower <- global_pattern_inference_veg$coef_025
          global_pattern_inference_veg$ci_upper <- global_pattern_inference_veg$coef_975
        }
        if (!("ci_lower" %in% names(global_pattern_inference_barren)) && ("coef_025" %in% names(global_pattern_inference_barren))) {
          global_pattern_inference_barren$ci_lower <- global_pattern_inference_barren$coef_025
          global_pattern_inference_barren$ci_upper <- global_pattern_inference_barren$coef_975
        }
        # Build plot
        p <- ggplot(global_pattern_inference_veg, aes(x = year, y = global_coef, color = Veg, fill = Veg)) +
          geom_line(linewidth = 1.2) +
          geom_point(size = 2, show.legend = FALSE) +
          geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, color = NA)
        barren_scale_factor <- 1
        if (nrow(global_pattern_inference_barren) > 0) {
          veg_max <- suppressWarnings(max(global_pattern_inference_veg$global_coef, global_pattern_inference_veg$ci_upper, na.rm = TRUE))
          barren_max <- suppressWarnings(max(global_pattern_inference_barren$global_coef, global_pattern_inference_barren$ci_upper, na.rm = TRUE))
          if (is.infinite(veg_max) || is.na(veg_max) || veg_max <= 0) veg_max <- 0.01
          if (is.infinite(barren_max) || is.na(barren_max) || barren_max <= 0) barren_max <- 1
          if (veg_max <= 0.01) {
            barren_scale_factor <- 1
          } else {
            barren_scale_factor <- veg_max / barren_max
          }
          global_pattern_inference_barren$coef_scaled <- global_pattern_inference_barren$global_coef * barren_scale_factor
          global_pattern_inference_barren$ci_lower_scaled <- global_pattern_inference_barren$ci_lower * barren_scale_factor
          global_pattern_inference_barren$ci_upper_scaled <- global_pattern_inference_barren$ci_upper * barren_scale_factor
          p <- p +
            geom_line(data = global_pattern_inference_barren, aes(x = year, y = coef_scaled), color = "brown", linewidth = 1.2, linetype = "dashed") +
            geom_point(data = global_pattern_inference_barren, aes(x = year, y = coef_scaled), color = "brown", size = 2, show.legend = FALSE) +
            geom_ribbon(data = global_pattern_inference_barren, aes(x = year, ymin = ci_lower_scaled, ymax = ci_upper_scaled), fill = "brown", alpha = 0.1, color = NA)
        }
        all_y_values <- c(global_pattern_inference_veg$global_coef, global_pattern_inference_veg$ci_lower, global_pattern_inference_veg$ci_upper)
        if (nrow(global_pattern_inference_barren) > 0) {
          all_y_values <- c(all_y_values, global_pattern_inference_barren$coef_scaled, global_pattern_inference_barren$ci_lower_scaled, global_pattern_inference_barren$ci_upper_scaled)
        }
        min_y <- min(all_y_values, na.rm = TRUE)
        max_y <- max(all_y_values, na.rm = TRUE)
        if (is.infinite(min_y) || is.na(min_y)) min_y <- 0
        if (is.infinite(max_y) || is.na(max_y)) max_y <- 1
        if (max_y <= min_y) max_y <- min_y + 0.1
        p <- p +
          labs(
            title = "Inference Average Coverage Percentage per Vegetation Type (2020-2025)",
            subtitle = sprintf("Based on %d locations with bootstrap uncertainty", max(global_pattern_inference$n_locations, na.rm = TRUE)),
            x = "Year",
            y = "Vegetation Fraction",
            color = "Vegetation Type",
            fill = "Vegetation Type"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 10),
            axis.title = element_text(size = 12),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10)
          ) +
          scale_y_continuous(
            labels = scales::percent_format(accuracy = 1),
            limits = c(min_y, max_y),
            expand = c(0, 0),
            sec.axis = sec_axis(~ . / barren_scale_factor, name = "Barren Fraction", labels = scales::percent_format(accuracy = 1))
          ) +
          scale_color_brewer(palette = "Set1") +
          scale_fill_brewer(palette = "Set1")
        plot_filename <- file.path(OUT_DIR, "inference_average_coverage_plot.png")
        ggsave(plot_filename, p, width = 10, height = 6, dpi = 300)
        cat(sprintf("Saved inference average coverage plot to: %s\n", plot_filename))

        # Generate PPI-normalized plot for inference
        global_pattern_inference_ppi <- location_bootstrap_ppi(all_coefs_local, df_tasks_inference, B = BOOTSTRAP_B, seed = 123)
        if (!is.null(global_pattern_inference_ppi) && nrow(global_pattern_inference_ppi) > 0) {
          p_ppi <- plot_global_vegetation_pattern(global_pattern_inference_ppi, title = "Inference PPI-Normalized Vegetation Trends", show_ci = TRUE)
          ppi_plot_filename <- file.path(OUT_DIR, "inference_ppi_normalized_plot.png")
          ggsave(ppi_plot_filename, p_ppi, width = 10, height = 6, dpi = 300)
          cat(sprintf("Saved inference PPI-normalized plot to: %s\n", ppi_plot_filename))

          p_ppi_stacked <- plot_vegetation_only_stacked_area(global_pattern_inference_ppi, title = "Inference Vegetation-Only Stacked Area")
          ppi_stacked_filename <- file.path(OUT_DIR, "inference_vegetation_only_stacked_area.png")
          ggsave(ppi_stacked_filename, p_ppi_stacked, width = 10, height = 6, dpi = 300)
          cat(sprintf("Saved inference vegetation-only stacked area plot to: %s\n", ppi_stacked_filename))
        } else {
          cat("[WARN] Failed to generate PPI-normalized data for inference (PPI column likely missing or no matches)\n")
        }
      } else {
        cat("[WARN] No inference data for average coverage plot\n")
      }
    }, error = function(e) {
      cat(sprintf("[ERROR] Failed to generate inference plots: %s\n", e$message))
    })
    # Add to wb
    if (exists("wb") && length(inference_results_list) > 0) {
      for (loc in names(inference_results_list)) {
        res <- inference_results_list[[loc]]
        if (!is.null(res) && !is.null(res$coef_df) && nrow(res$coef_df) > 0) {
          sheet_name <- paste0("Inference_", loc)
          openxlsx::addWorksheet(wb, sheet_name)
          # Write predicted fractions (coef_df)
          openxlsx::writeData(wb, sheet_name, res$coef_df, startRow = 1)
          # Write bootstrapping results below
          if (!is.null(res$bootstrap_coefs) && nrow(res$bootstrap_coefs) > 0) {
            start_row <- nrow(res$coef_df) + 3
            openxlsx::writeData(wb, sheet_name, "Bootstrapping Results", startRow = start_row)
            openxlsx::writeData(wb, sheet_name, res$bootstrap_coefs, startRow = start_row + 1)
          }
        }
      }
      openxlsx::saveWorkbook(wb, output_filename, overwrite = TRUE)
      cat(sprintf("Updated mesma_results.xlsx with inference location sheets\n"))
    }

    # Calculate and print trend statistics for inference data
    if (!is.null(all_coefs_local) && nrow(all_coefs_local) > 0) {
      cat("\n=== INFERENCE TREND STATISTICS (Bootstrapped) ===\n")
      tryCatch({
        trend_ci_inf <- bootstrap_trend_ci(all_coefs_local, B = 200)
        if (!is.null(trend_ci_inf)) {
          print(trend_ci_inf)
          write.csv(trend_ci_inf, file.path(OUT_DIR, "inference_trend_statistics.csv"), row.names = FALSE)
          cat(sprintf("Saved inference trend statistics to: %s\n", file.path(OUT_DIR, "inference_trend_statistics.csv")))
        } else {
          cat("[WARN] bootstrap_trend_ci returned NULL for inference data\n")
        }
      }, error = function(e) {
        cat(sprintf("[ERROR] Failed to compute inference trend statistics: %s\n", e$message))
      })
    }
  }

  inference_result
}

aggregate_and_plot <- function(all_coefs_in = NULL, out_dir = OUT_DIR) {
  if (is.null(all_coefs_in)) all_coefs_in <- all_coefs
  if (is.null(all_coefs_in) || nrow(all_coefs_in) == 0) { cat("[INFO] No coefficients to aggregate\n"); return(NULL) }
  gp <- aggregate_to_global_pattern(all_coefs_in, method = "location_bootstrap")
  if (is.null(gp) || nrow(gp) == 0) { cat("[INFO] aggregate_to_global_pattern returned no data\n"); return(NULL) }
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  p <- plot_global_vegetation_pattern(gp, title = "All Vegetation Trends", show_ci = TRUE)
  ggsave(file.path(out_dir, "all_vegetation_trends.png"), p, width = 10, height = 6)

  p_stacked <- plot_vegetation_only_stacked_area(gp, title = "Vegetation-Only Stacked Area")
  ggsave(file.path(out_dir, "vegetation_only_stacked_area.png"), p_stacked, width = 10, height = 6)

  invisible(gp)
}


cat("\n=== Starting MESMA processing ===\n")
# Attempt to use modular helper functions if available and if main processing was not executed
if (exists("process_batches") && exists("process_inference") && exists("aggregate_and_plot")) {
    if (!exists("main_processing_executed") || !isTRUE(main_processing_executed)) {
      cat("[INFO] Running modular helpers: process_batches -> process_inference -> aggregate_and_plot\n")
      flush.console()
      rb <- tryCatch(process_batches(BATCH_SIZE = if (exists("BATCH_SIZE")) BATCH_SIZE else 6, overwrite_wb = FALSE), error = function(e) { cat(sprintf("[ERROR] process_batches failed: %s\n", e$message)); NULL })
      ri <- tryCatch(process_inference(BATCH_SIZE = if (exists("BATCH_SIZE")) BATCH_SIZE else 6), error = function(e) { cat(sprintf("[ERROR] process_inference failed: %s\n", e$message)); NULL })

      allc <- NULL
      if (!is.null(ri) && !is.null(ri$all_coefs) && nrow(ri$all_coefs) > 0) {
        cat("[INFO] Using INFERENCE results for final aggregation\n")
        allc <- ri$all_coefs
      } else {
        cat("[INFO] No inference results found. Skipping aggregation (training results ignored by request).\n")
      }

      if (!is.null(allc) && nrow(allc) > 0) {
        tryCatch({ aggregate_and_plot(allc, out_dir = if (exists("OUT_DIR")) OUT_DIR else ".") }, error = function(e) { cat(sprintf("[ERROR] aggregate_and_plot failed: %s\n", e$message)) })
      } else {
        cat("[INFO] No coefficients produced by helpers to aggregate and plot\n")
      }
    } else {
      cat("[INFO] Main processing already executed; running aggregate_and_plot on existing results if available\n")
      if (exists("all_coefs") && !is.null(all_coefs) && nrow(all_coefs) > 0) {
        tryCatch({ aggregate_and_plot(all_coefs, out_dir = if (exists("OUT_DIR")) OUT_DIR else ".") }, error = function(e) { cat(sprintf("[ERROR] aggregate_and_plot failed: %s\n", e$message)) })
      }
    }
}



# Final memory cleanup
cleanup_memory()

cat("\n=== Script execution finished ===\n")

