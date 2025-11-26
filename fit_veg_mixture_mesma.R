 
library(zoo)
library(dplyr)
library(cluster)

fit_cost_mkl <- function(obs, weight, t_row) {
  # fit_cost_mkl: vectorized cost calculation helper
  if (length(obs) != length(t_row)) t_row <- t_row[seq_len(length(obs))]
  if (any(!is.finite(c(obs, weight, t_row)))) stop("fit_cost_mkl: inputs contain non-finite values; clean data before calling")
  
  # Use vectorized operations
  pred <- weight * t_row
  residuals <- (obs * weight) - pred
  sum(residuals^2)
}

# Set up progressr handlers for ETA display (use single-line text handler)
if (requireNamespace("progressr", quietly = TRUE)) {
  # Use the txt handler to avoid multi-line progress chatter.
  # Width adapts to console or uses 60 characters as default for neat display.
  width <- tryCatch(getOption("width", 60), error = function(e) 60)
  progressr::handlers(progressr::handler_txt(width = min(max(as.integer(width), 40), 120)))
}

# Canonical optimal indices (required - no fallback alternatives)
OPTIMAL_INDICES <- c(
  #"DVI", # Difference Vegetation Index
  #"OSAVI", # Optimized Soil-Adjusted Vegetation Index
  "MCARI", # Modified Chlorophyll Absorption in Reflectance Index
  "CRI", # Carotenoid Reflectance Index
  "PRI", # Photochemical Reflectance Index
  "NIRv", # NIR*NDVI vigor
  "PSRI", # Plant Senescence Reflectance Index
  "NBR", # Burn ratio (NIR-SWIR2)/(NIR+SWIR2)
  "TCW", # Total Column Water
  "TCG"#, # Triangular Chlorophyll Index (green-based)
  #"MNDWI" # Modified Normalized Difference Water Index
)

# Raw spectral bands (optional - included if present)
RAW_BANDS <- c("blue", "green", "red", "nir", "swir1", "swir2")

# Function to calculate robust variance using STL decomposition + MAD
# Helper: compute MAD^2 with minimal sample requirement (top-level helper, no nested defs)
compute_mad2 <- function(x, min_samples = 3) {
  x <- x[is.finite(x)]
  if (length(x) < min_samples) return(NA_real_)
  m <- mad(x, na.rm = TRUE, constant = 1.4826)
  if (!is.finite(m)) return(NA_real_)
  m^2
}

calc_moving_var <- function(df, index_name, window = 14, span_loess = 0.1, min_obs_loess = 6) {
  # STL-only implementation: expects enough data to estimate seasonal components.
  # If a location/time-series does not have sufficient length for STL, returns NA for those rows.
  if (!"date" %in% names(df)) stop("calc_moving_var: df must have a 'date' column")
  if (!index_name %in% names(df)) stop(sprintf("calc_moving_var: index '%s' not found in df", index_name))

  n <- nrow(df)
  out <- rep(NA_real_, n)

  process_series_stl <- function(dates, vals) {
    # dates: Date or integer-like; vals: numeric vector aligned to dates
    if (length(vals) < 1) return(rep(NA_real_, length(vals)))
    finite_idx <- which(is.finite(vals) & !is.na(dates))
    if (length(finite_idx) == 0) return(rep(NA_real_, length(vals)))

    # Require a minimum length to run STL (two full seasonal periods recommended)
    # We use frequency = 365 (daily seasonality). STL requires length >= 2*frequency.
    freq <- 365
    min_len_for_stl <- 2 * freq

    if (length(vals) < min_len_for_stl) {
      # Not enough data to run STL robustly: return NA for this series
      return(rep(NA_real_, length(vals)))
    }

    # Build a regularly spaced series over the dates (fill missing days with NA)
    dts <- as.Date(dates)
    doy <- as.integer(lubridate::yday(dts))
    full_days <- seq(from = as.Date(min(dts, na.rm = TRUE)), to = as.Date(max(dts, na.rm = TRUE)), by = "day")
    if (length(full_days) < min_len_for_stl) return(rep(NA_real_, length(vals)))
    full_vec <- rep(NA_real_, length(full_days))
    pos_map <- match(dts, full_days)
    full_vec[pos_map] <- as.numeric(vals)

    # Require a minimum number of finite observations to interpolate
    if (sum(is.finite(full_vec)) < max(10, floor(0.1 * length(full_vec)))) return(rep(NA_real_, length(vals)))

    # Simple linear interpolation to fill gaps (STL cannot handle NA). This is a minimal pre-step.
    idx_finite <- which(is.finite(full_vec))
    full_vec_interp <- full_vec
    if (any(!is.finite(full_vec_interp))) {
      interp_vals <- tryCatch(
        {
          approx(x = idx_finite, y = full_vec[idx_finite], xout = seq_along(full_vec), rule = 2)$y
        },
        error = function(e) stop(sprintf("calc_moving_var: interpolation failed: %s", e$message))
      )
      if (is.null(interp_vals) || any(!is.finite(interp_vals))) return(rep(NA_real_, length(vals)))
      full_vec_interp <- interp_vals
    }

    # Run STL; fail loudly on error
    stl_out <- tryCatch(
      {
        stats::stl(ts(full_vec_interp, frequency = freq), s.window = "periodic", robust = TRUE)
      },
      error = function(e) stop(sprintf("calc_moving_var: STL decomposition failed: %s", e$message))
    )

    resid <- as.numeric(stl_out$time.series[, "remainder"])
    # Compute rolling MAD^2 using compute_mad2
    if (length(resid) < window) return(rep(NA_real_, length(vals)))
    rv <- zoo::rollapply(resid, width = window, FUN = function(x) compute_mad2(x, min_samples = max(3, floor(window/2))), fill = NA, align = "center")
    # Map back to original dates
    res_for_rows <- rv[pos_map]
    as.numeric(res_for_rows)
  }

  if ("location_id" %in% names(df)) {
    locs <- unique(df$location_id)
    for (loc in locs) {
      rows <- which(df$location_id == loc)
      if (length(rows) < 1) next
      sub_dates <- df$date[rows]
      ord <- order(sub_dates)
      rrows <- rows[ord]
      vals <- df[[index_name]][rrows]
      out_vals <- process_series_stl(sub_dates, vals)
      out[rrows] <- out_vals
    }
  } else {
    out <- process_series_stl(df$date, df[[index_name]])
  }

  out[!is.finite(out)] <- NA_real_
  out
}

# Helper: normalize variants of the 'no soil' column name to a single
# canonical name 'no soil' (space). This accepts 'no_soil', 'no.soil',
# '.__no soil__' etc and copies the first matching column into 'no soil'
# if present. This makes downstream access with backticks consistent.
normalize_no_soil_col <- function(tbl) {
  if (is.null(tbl) || !is.data.frame(tbl)) return(tbl)
  nm <- names(tbl)
  candidates <- c("no soil", "no_soil", "no.soil", "__no soil__", "__no_soil__", ".__no soil__", ".__no_soil__")
  # Keep the user's exact 'no soil' when present; otherwise pick the first matching candidate
  # If present, coerce the column to safe numeric form (handles factors/strings)
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


# Helper: safely coerce values to numeric. Handles factors, character "TRUE"/"FALSE",
# trims whitespace and avoids direct as.numeric(factor) pitfalls.
safe_as_numeric <- function(x) {
  if (is.null(x)) return(x)
  # Factors -> characters first
  if (is.factor(x)) x <- as.character(x)
  # Characters: trim, handle logical strings
  if (is.character(x)) {
    s <- trimws(x)
    lower <- tolower(s)
    # Map boolean-like strings to 1/0
    lower[lower %in% c("true", "t")] <- "1"
    lower[lower %in% c("false", "f")] <- "0"
    # Attempt numeric conversion
    suppressWarnings(num <- as.numeric(lower))
    return(num)
  }
  # Already numeric (or integer) -> coerce to numeric
  if (is.numeric(x)) return(as.numeric(x))
  # Otherwise try character roundtrip
  suppressWarnings(num <- as.numeric(as.character(x)))
  num
}

# Persistent parallel backend (set up once)
setup_parallel_backend <- function() {
  if (isTRUE(PARALLEL_ENABLE) && requireNamespace("future", quietly = TRUE)) {
    old_plan <- future::plan()
    options(future.globals.maxSize = max(getOption("future.globals.maxSize", 2e9), 4e9))
    future::plan(future::multisession, workers = PARALLEL_WORKERS, gc = TRUE, earlySignal = TRUE)
    return(function() {
      try(future::plan(old_plan), silent = TRUE)
    })
  }
  function() {}
}

# User-configurable defaults
OUTPUT_DIR <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/phenology_results"
INPUT_CSV <- file.path(OUTPUT_DIR, "hls_phenology_data.csv")
OUT_DIR <- file.path(OUTPUT_DIR, "veg_mixture_fit")

# Training configuration defaults
# The script uses all available years for training by default.
TRAIN_YEARS <- 2019:2024

# Allow skipping a minimal number of DOYs per-location when computing sufficiency
# This is useful for small gaps — set to 0 to preserve the strict 50% requirement
MIN_SKIP_DOYS_PER_LOCATION <- 2L

# PCA / factor projection tuning
MAX_PCA_COMPONENTS <- 20L
MAX_FACTORS_CAP <- 20L

# General algorithm toggles
FAST_VAR <- TRUE
factor_mode <- FALSE
SHAPE_NORMALIZATION_ENABLE <- FALSE
TEMPORAL_BUDGET <- 10L
TOPK_VARIANTS <- 2L
ENABLE_PHASE_ALIGNMENT <- FALSE
REFERENCE_PHASE_MARKERS <- c(1, 90, 180, 270, 365)
ENABLE_MULTISCALE <- FALSE
MULTISCALE_WINDOWS <- c(7L, 14L, 30L)
ENABLE_QP_SOLVER <- TRUE
COMBO_PARALLEL_ENABLE <- FALSE
EARLY_STOP_RMSE_THRESHOLD <- 0.0
ENABLE_DIAGNOSTICS <- TRUE

# Combination expansion safety thresholds
COMBO_SAFE_EXPAND_LIMIT <- 1e6    # fully expand grid up to this many combos
COMBO_ABORT_LIMIT <- 5e7          # abort if combos exceed this hard limit

# Bootstrap settings 
BOOTSTRAP_B <- 200L
BOOT_MIN_REPS_PER_VEG <- 10L
MIN_OBS_FOR_BOOT <- 5L
ENABLE_UNCERTAINTY <- TRUE
GLSBB_MIN_BLOCK <- 5L
GLSBB_MAX_BLOCK <- 60L

# Optimization/solver defaults
VARIANCE_THRESHOLD <- 0.90
MAX_VEG_COMPONENTS <- 8
GAM_K_MAX <- 40
GAM_GAMMA <- 1.0

# Index selection and prefiltering
USE_INDICES_MIN <- 1L
MIN_INDEX_SD <- 0.05

# Sample-balancing and augmentation
ENABLE_SAMPLE_BALANCING <- TRUE

# Soil-preprocessing: run a preliminary soil/no-soil MESMA and subtract soil
# contribution from all spectral time series before building vegetation libs.
# This helps avoid soil contamination when clustering vegetation time-series.
ENABLE_SOIL_PREPROCESS <- TRUE
# Threshold (0-1) to treat a geo-labelled sample as pure soil when building the soil prototype
SOIL_PURE_THRESHOLD <- 0.95
# Minimum number of pure-soil rows required to build a robust soil prototype
SOIL_MIN_SAMPLES <- 3L
# When insufficient pure soil points are present, the soil prototype will not
# be built (returns NULL) unless a minimum number of candidate rows meet
# the SOIL_PURE_THRESHOLD.

# Memory-safe clustering limits (sampling parameters to avoid RAM overload)
MAX_PROJECTIONS_PER_VEG <- 25000L  # subsample before clustering to avoid OOM
SILHOUETTE_SAMPLE_SIZE <- 20000L   # subsample for silhouette distance matrix
MEDOID_SAMPLE_SIZE <- 10000L       # subsample for medoid distance computation

# Vegetation whitelist
ALLOWED_VEG <- c("populus", "tamarix", "phragmites")

# Numeric safety constants
EPS_SIGMA <- 1e-8
LOWER_BND <- 0

# DOY / index presence tuning
MIN_IDX_PRESENCE <- 0.5

# Parallel / progress settings
PROGRESS_EVERY_TASK <- 25
PROGRESS_LOG_TO_FILE <- TRUE
PROGRESS_BAR <- TRUE
PARALLEL_ENABLE <- TRUE
PARALLEL_WORKERS <- tryCatch(
  {
    if (requireNamespace("parallel", quietly = TRUE)) max(1L, parallel::detectCores(logical = TRUE) - 1L) else 1L
  },
  error = function(...) 1L
)
COMBO_PARALLEL_WORKERS <- max(1L, floor(PARALLEL_WORKERS/2))
PERSISTENT_PARALLEL_BACKEND <- TRUE

# Strict validation: enforce extractor/predictor output format only
if (!file.exists(INPUT_CSV)) stop(paste0("Required input CSV not found: ", INPUT_CSV))
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

# Ensure all OPTIMAL_INDICES are present in the CSV
missing_idx <- setdiff(OPTIMAL_INDICES, names(raw_df))
if (length(missing_idx) > 0) {
  stop(paste0("INPUT_CSV missing required indices: ", paste(missing_idx, collapse = ", ")))
}

df <- raw_df
# Normalize any variants of the no-soil column name for consistent downstream access
df <- normalize_no_soil_col(df)

# Add timing for major operations
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

# Train/Test Split Implementation
cat("\n=== TRAIN/INFERENCE DATA CONFIGURATION ===\n")
cat(sprintf("Training years (config): %s\n", paste(TRAIN_YEARS, collapse = ", ")))
cat("Note: the script will always use all available years for data; TEST_YEARS is not supported.\n")

df$year <- lubridate::year(df$date)

## Use all available data for training by default. Keep a df_train object to
## preserve a clear training/inference distinction for downstream code.
df_train <- df
cat(sprintf(
  "Training dataset: %d rows from %d locations\n",
  nrow(df_train), length(unique(df_train$location_id))
))

## Inference/testing uses the full dataset (all years).
df_test <- df
cat(sprintf(
  "Testing dataset: %d rows from %d locations\n",
  nrow(df_test), length(unique(df_test$location_id))
))

df_full <- df
df_full <- normalize_no_soil_col(df_full)
df <- df_train

# Ensure training subset has at least some samples for each allowed vegetation class
# If a class from ALLOWED_VEG has zero samples in the available training set,
# add a small number of samples for that class from the full dataset.
if (exists("df_full") && "Veg" %in% names(df_full) && length(ALLOWED_VEG) > 0) {
  missing_vegs <- sapply(ALLOWED_VEG, function(v) {
    sum(tolower(df$Veg) == v, na.rm = TRUE)
  })
  missing_names <- names(missing_vegs)[missing_vegs == 0]
  if (length(missing_names) > 0) {
    for (mv in missing_names) {
      # find candidates anywhere in the full dataset
      cand <- df_full[tolower(df_full$Veg) == mv, , drop = FALSE]
      if (nrow(cand) > 0) {
        # Add a modest number of examples so the training set has representation
        add_n <- min(nrow(cand), max(5L, as.integer(floor(nrow(df) / 10))))
        add_n <- max(1L, add_n)
        to_add <- cand[seq_len(add_n), , drop = FALSE]
        df <- rbind(df, to_add)
        cat(sprintf("[NOTICE] Added %d samples from non-training years for Veg='%s' to ensure representation in training set\n", nrow(to_add), mv))
      } else {
        cat(sprintf("[WARNING] No samples found anywhere for Veg='%s'; cannot add examples to training set\n", mv))
      }
    }
    # re-report training dataset size after augmentation
    cat(sprintf("Training dataset after augmentation: %d rows from %d locations\n", nrow(df), length(unique(df$location_id))))
    # keep df_train consistent with df (augmented)
    df_train <- df
  }
}

cat("Using training data for vegetation library construction\n")
cat("=====================================\n\n")

# Initialize persistent parallel backend if enabled
cleanup_parallel <- setup_parallel_backend()
on.exit(cleanup_parallel(), add = TRUE)

# Helper: project a vector onto the probability simplex
project_to_simplex <- function(v) {
  if (!is.numeric(v)) stop("project_to_simplex: input must be numeric")
  n <- length(v)
  if (n == 0) {
    v
  }
  u <- sort(v, decreasing = TRUE)
  cssv <- cumsum(u)
  rho <- max(which(u - (cssv - 1) / seq_along(u) > 0))
  theta <- (cssv[rho] - 1) / rho
  w <- pmax(v - theta, 0)
  if (sum(w) <= 0) {
    rep(1 / n, n)
  }
  w / sum(w)
}

# Progress logging helpers
LOG_FILE <- tryCatch(
  {
    if (exists("OUT_DIR") && is.character(OUT_DIR) && nchar(OUT_DIR) > 0) {
      file.path(OUT_DIR, "fit_veg_mesma.log")
    } else {
      file.path(getwd(), "fit_veg_mesma.log")
    }
  },
  error = function(e) stop(sprintf("Failed to compute LOG_FILE path: %s", e$message))
)

log_msg <- function(...) {
  ts <- format(Sys.time(), "%H:%M:%S")
  msg <- sprintf("[%s] %s\n", ts, sprintf(...))
  cat(msg)
  if (isTRUE(PROGRESS_LOG_TO_FILE)) try(cat(msg, file = LOG_FILE, append = TRUE), silent = TRUE)
}

# Helper: canonical location_id creation from lon/lat (centralize formatting)
make_location_id <- function(lon, lat) {
  if (length(lon) == 1 && length(lat) == 1) {
    lon_num <- as.numeric(lon); lat_num <- as.numeric(lat)
    if (!is.finite(lon_num) || !is.finite(lat_num)) return(NA_character_)
    sprintf("L_%0.4f_%0.4f", round(lon_num, 4), round(lat_num, 4))
  } else {
    # vectorized form: return character vector
    sapply(seq_along(lon), function(i) {
      lon_i <- as.numeric(lon[i]); lat_i <- as.numeric(lat[i])
      if (!is.finite(lon_i) || !is.finite(lat_i)) NA_character_ else sprintf("L_%0.4f_%0.4f", round(lon_i,4), round(lat_i,4))
    }, USE.NAMES = FALSE)
  }
}

# Memory helpers: chunked rbind to avoid building huge intermediate objects
chunked_rbind <- function(lst, chunk_size = 50L) {
  if (length(lst) == 0) return(matrix(nrow = 0, ncol = 0))
  if (length(lst) == 1) return(lst[[1]])
  n <- length(lst)
  # Iteratively bind in chunks
  i <- 1L
  acc <- NULL
  while (i <= n) {
    end <- min(n, i + chunk_size - 1L)
    block <- do.call(rbind, lst[i:end])
    if (is.null(acc)) acc <- block else acc <- rbind(acc, block)
    # free memory for block
    rm(block); gc()
    i <- end + 1L
  }
  acc
}

# Debug helper: produce a nil return object for the main processing loop
dbg_return_null <- function(reason = NULL) {
  if (!is.null(reason)) cat(sprintf("[DEBUG] abort: %s\n", as.character(reason)))
  coef_df <- data.frame(location_id = NA_character_, year = NA_integer_, Veg = NA_character_, coef = NA_real_, rmse = NA_real_, stringsAsFactors = FALSE)
  coef_df$coef_lo <- NA_real_
  coef_df$coef_hi <- NA_real_
  diag_df <- data.frame(location_id = NA_character_, year = NA_integer_, vegetated_fraction = NA_real_, barren_fraction = NA_real_, stringsAsFactors = FALSE)
  list(
    coef_df = coef_df,
    variant_trajectory = NULL,
    diagnostics = diag_df,
    uncertainty = NULL,
    q10_dvi = NA_real_,
    q90_dvi = NA_real_,
    vegetated_fraction = NA_real_,
    barren_fraction = NA_real_
  )
}

# Clamp for parallel workers memory usage: limit workers according to available RAM
recommended_workers <- function(max_workers_hint = NULL) {
  # crude heuristic: check memory.limit on Windows or try to estimate via gc
  total_mem_gb <- tryCatch({
    if (.Platform$OS.type == "windows") {
      mem <- memory.limit() / 1024
      mem
    } else {
      NA_real_
    }
  }, error = function(e) stop(sprintf("recommended_workers: memory inspection failed: %s", e$message)))
  # conservative default
  if (!is.na(total_mem_gb) && total_mem_gb > 0) {
    # one worker per ~2 GB
    max(1L, floor(total_mem_gb / 2))
  } else {
    if (!is.null(max_workers_hint)) max(1L, as.integer(max_workers_hint)) else PARALLEL_WORKERS
  }
}

# Temporal compression and template helpers (promoted to top-level so they're available everywhere)
compute_information_content <- function(pca_mat) {
  n_days <- nrow(pca_mat); n_pcs <- ncol(pca_mat)
  info <- rep(0, n_days)
  if (n_days <= 1 || n_pcs == 0) return(info)
  w <- 7L; half <- floor(w/2)
  for (d in seq_len(n_days)) {
    s <- max(1L, d - half); e <- min(n_days, d + half)
    if (e - s + 1 >= 2) {
      local_var <- sum(apply(pca_mat[s:e, , drop = FALSE], 2, stats::var, na.rm = TRUE), na.rm = TRUE)
    } else local_var <- 0
    grad <- if (d > 1 && d < n_days) sum(abs(pca_mat[d+1, ] - pca_mat[d-1, ]), na.rm = TRUE) else 0
    spec_cx <- stats::var(as.numeric(pca_mat[d, ]), na.rm = TRUE)
    info[d] <- sum(c(local_var, grad, ifelse(is.finite(spec_cx), spec_cx, 0)), na.rm = TRUE)
  }
  info <- zoo::rollmean(info, k = 3, fill = "extend", align = "center")
  rng <- range(info, na.rm = TRUE)
  if (is.finite(rng[1]) && is.finite(rng[2]) && rng[2] > rng[1]) info <- (info - rng[1])/(rng[2]-rng[1]) else info[] <- 0
  info[!is.finite(info)] <- 0
  info
}

create_adaptive_grid <- function(information_content, budget) {
  n_days <- length(information_content)
  if (n_days < 1) return(integer(0))
  cum <- cumsum(information_content)
  if (!is.finite(cum[length(cum)]) || cum[length(cum)] <= 0) {
    pts <- unique(round(seq(1, n_days, length.out = max(2, budget))))
  } else {
    cum <- cum / cum[length(cum)]
    targets <- if (budget <= 1) 0 else seq(0, 1, length.out = budget)
    pts <- sapply(targets, function(tg) which.min(abs(cum - tg)))
    pts[1] <- 1L; pts[length(pts)] <- n_days
    pts <- sort(unique(as.integer(pts)))
  }
  as.integer(pts)
}

extract_grid_features <- function(pca_mat, temporal_grid, information_content) {
  n_days <- nrow(pca_mat); k <- ncol(pca_mat)
  feats <- c()
  if (length(temporal_grid) == 0) return(rep(0, k))
  for (i in seq_along(temporal_grid)) {
    if (i == 1) {
      next_mid <- if (length(temporal_grid) > 1) floor((temporal_grid[i] + temporal_grid[i + 1]) / 2) else n_days
      s <- 1L; e <- max(1L, min(n_days, next_mid))
    } else if (i == length(temporal_grid)) {
      prev_mid <- floor((temporal_grid[i - 1] + temporal_grid[i]) / 2)
      s <- max(1L, min(n_days, prev_mid)); e <- n_days
    } else {
      prev_mid <- floor((temporal_grid[i - 1] + temporal_grid[i]) / 2)
      next_mid <- floor((temporal_grid[i] + temporal_grid[i + 1]) / 2)
      s <- max(1L, min(n_days, prev_mid)); e <- max(1L, min(n_days, next_mid))
    }
    if (e < s) { tmp <- s; s <- e; e <- tmp }
    window <- pca_mat[s:e, , drop = FALSE]
    wts <- information_content[s:e]
    if (nrow(window) > 0) {
      if (sum(wts, na.rm = TRUE) > 0) {
        wts <- wts / sum(wts, na.rm = TRUE)
        wavg <- colSums(window * wts)
      } else {
        wavg <- colMeans(window, na.rm = TRUE)
      }
    } else {
      wavg <- rep(0, k)
    }
    feats <- c(feats, as.numeric(wavg))
  }
  feats[!is.finite(feats)] <- 0
  feats
}

extract_multiscale_features <- function(mat, temporal_grid, information_content, windows = MULTISCALE_WINDOWS) {
  out <- list()
  out[["fine"]] <- extract_grid_features(mat, temporal_grid, information_content)
  if (length(windows) > 0) {
    for (w in windows) {
      w <- as.integer(w)
      if (!is.finite(w) || w < 2) next
      smoothed <- apply(mat, 2, function(x) zoo::rollmean(x, k = w, fill = "extend", align = "center"))
      out[[paste0("scale_", w)]] <- extract_grid_features(smoothed, temporal_grid, information_content)
    }
  }
  as.numeric(unlist(out))
}
precompute_compressed_templates <- function(mesma_lib, budget = TEMPORAL_BUDGET) {
  if (length(mesma_lib) == 0) return(function(...) NULL)
  uniform_info <- rep(0.5, 365)
  standard_grids <- list(
    sparse = create_adaptive_grid(uniform_info, ceiling(budget * 0.7)),
    medium = create_adaptive_grid(uniform_info, budget),
    dense = create_adaptive_grid(uniform_info, ceiling(budget * 1.3))
  )
  template_db <- new.env(hash = TRUE, parent = emptyenv())
  for (veg in names(mesma_lib)) {
    variants <- mesma_lib[[veg]]
    for (variant in variants) {
      tmat <- variant$T_pca
      if (is.null(tmat) || !is.matrix(tmat) || nrow(tmat) == 0) next
      tcur <- tmat
      if (isTRUE(ENABLE_PHASE_ALIGNMENT)) {
        tcur <- align_to_phenological_phase(tcur, reference_phase = REFERENCE_PHASE_MARKERS)
      }
      info <- compute_information_content(tcur)
      for (grid_name in names(standard_grids)) {
        grid <- standard_grids[[grid_name]]
        feats <- if (isTRUE(ENABLE_MULTISCALE)) {
          extract_multiscale_features(tcur, grid, info, windows = MULTISCALE_WINDOWS)
        } else {
          extract_grid_features(tcur, grid, info)
        }
        
        key <- paste(veg, variant$variant_id, grid_name, sep = "|")
        assign(key, feats, envir = template_db)
      }
    }
  }
  
  function(veg, variant_id, grid_type = "medium") {
      key <- paste(veg, variant_id, grid_type, sep = "|")
      if (exists(key, envir = template_db, inherits = FALSE)) get(key, envir = template_db, inherits = FALSE) else NULL
  }
}

# Whitening helper: moved to top-level so other stages can reuse it
whiten_matrix <- function(X, epsilon = 1e-6) {
  if (is.null(X) || nrow(X) < 2 || ncol(X) < 1) return(list(Xw = X, W = diag(ncol(X)), mu = rep(0, ncol(X))))
  mu <- colMeans(X)
  Xc <- sweep(X, 2, mu, "-")
  Sigma <- cov(Xc)
  eig <- eigen(Sigma)
  vals <- eig$values
  vals[vals < epsilon] <- epsilon
  D_inv_sqrt <- diag(1 / sqrt(vals))
  W <- eig$vectors %*% D_inv_sqrt %*% t(eig$vectors)
  Xw <- Xc %*% W
  list(Xw = Xw, W = W, mu = mu)
}

# Helper: index of medoid row (row with minimal total distance)
# For large matrices, sample to avoid O(N^2) memory allocation
medoid_row_index <- function(M) {
  if (is.null(M) || !is.matrix(M) || nrow(M) == 0) return(NA_integer_)
  if (nrow(M) == 1) return(1L)
  X <- M
  X[!is.finite(X)] <- 0
  
  # For large matrices, compute medoid on a random subsample
  n <- nrow(X)
  
  if (n > MEDOID_SAMPLE_SIZE) {
    # Sample rows, find medoid among sample, return index in original matrix
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

# Cosine similarity helper (promoted to top-level)
cos_sim <- function(a, b) {
  da <- sqrt(sum(a * a)); db <- sqrt(sum(b * b))
  if (da == 0 || db == 0) return(0)
  sum(a * b) / (da * db)
}

# Helper: parallel map
.run_map <- function(X, FUN) {
  f_FUN <- FUN
  # Use txt handler explicitly for the local run to avoid multi-line spamming
  if (requireNamespace("progressr", quietly = TRUE)) progressr::handlers("txt")

  if (!PARALLEL_ENABLE) {
    if (!requireNamespace("progressr", quietly = TRUE)) {
      stop("progressr package required for serial processing")
    }
    progressr::with_progress({
      p <- progressr::progressor(steps = length(X))
      lapply(X, function(x) {
        res <- f_FUN(x)
        p()
        res
      })
    })
  } else {
    if (!requireNamespace("future.apply", quietly = TRUE) || !requireNamespace("future", quietly = TRUE)) {
      stop("future.apply and future packages required for parallel processing")
    }
    if (!requireNamespace("progressr", quietly = TRUE)) {
      stop("progressr package required for progress tracking")
    }

    if (!isTRUE(PERSISTENT_PARALLEL_BACKEND)) {
      old_plan <- future::plan()
      options(future.globals.maxSize = max(getOption("future.globals.maxSize", 2e9), 4e9))
      future::plan(future::multisession, workers = PARALLEL_WORKERS)
      on.exit(future::plan(old_plan))
    }

    progressr::with_progress({
      p <- progressr::progressor(steps = length(X))
      future.apply::future_lapply(X, function(x) {
        res <- f_FUN(x)
        p()
        res
      }, future.seed = TRUE)
    })
  }
}

# Read GeoJSON for vegetation ground truth
if (!"date" %in% names(df) && "Date" %in% names(df)) df$date <- df$Date
if ("date" %in% names(df)) {
  df$date <- as.Date(df$date)
  if (!"year" %in% names(df)) df$year <- as.integer(lubridate::year(df$date))
}

geojson_path <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/updated_zuizer_zonder_foto_UTM.geojson"
if (!file.exists(geojson_path)) stop(paste0("GeoJSON points not found at ", geojson_path))
gpts_raw <- sf::st_read(geojson_path, quiet = TRUE)

geojson_names <- names(gpts_raw)
normalized_names <- gsub("[^a-z0-9]+", "_", tolower(geojson_names))
no_soil_col <- geojson_names[normalized_names == "no_soil"]
if (length(no_soil_col) == 0) stop("GeoJSON point data requires a 'no soil' column")
no_soil_raw <- gpts_raw[[no_soil_col[1]]]
no_soil_vals <- safe_as_numeric(no_soil_raw)
if (is.null(no_soil_vals) || (!is.numeric(no_soil_vals) && !is.logical(no_soil_vals))) {
  stop("'no soil' column must be numeric, logical, or character coercible to numeric")
}
  if (any(is.na(no_soil_vals))) {
  # Leave missing 'no soil' values as NA instead of coercing to 0.
  # Coercing to 0 earlier caused many locations to be treated as fully barren/mixed
  # when GeoJSON contained missing entries. Keeping NA helps us detect missing
  # metadata and prevents silently forcing binary fractions.
  cat("[NOTICE] Found missing/non-numeric 'no soil' values in GeoJSON; keeping as NA (do not coerce to 0).\n")
}
## Validate numeric range only for non-missing values; keep NA values allowed
non_na_no_soil <- no_soil_vals[!is.na(no_soil_vals)]
if (length(non_na_no_soil) > 0 && any(non_na_no_soil < 0 | non_na_no_soil > 1)) stop("'no soil' values must lie within [0,1]")
gpts_raw$`.__no soil__` <- no_soil_vals

matched_cols <- names(gpts_raw)[tolower(names(gpts_raw)) %in% c("vegetation", "veg", "class")]
if (length(matched_cols) > 0) {
  veg_col_orig <- matched_cols[1]
  gpts_raw$.__veg__ <- as.character(gpts_raw[[veg_col_orig]])
} else {
  gpts_raw$.__veg__ <- NA_character_
}
coords <- sf::st_coordinates(gpts_raw)
gpts_raw$.__lon__ <- coords[, 1]
gpts_raw$.__lat__ <- coords[, 2]
# Preserve the coordinate-based id under a separate name so we can use row numbers as the canonical id
gpts_raw$location_id_geo <- make_location_id(gpts_raw$.__lon__, gpts_raw$.__lat__)
# Use geojson row number as the canonical location_id for matching with CSV indices
gpts_raw$location_id <- as.character(seq_len(nrow(gpts_raw)))

gpts_map <- sf::st_drop_geometry(gpts_raw) %>%
  dplyr::select(location_id, Veg = .__veg__, `no soil` = `.__no soil__`) %>%
  dplyr::mutate(location_row = as.character(seq_len(dplyr::n()))) %>%
  dplyr::distinct(location_id, .keep_all = TRUE)

if (nrow(gpts_map) == 0) stop("GeoJSON mapping produced no valid points")

if ("location_id" %in% names(df) && "location_id" %in% names(gpts_map)) {
  # Robust left-join: normalize and coerce types so mismatched formats (whitespace/case)
  # do not prevent the mapping.
  if (!is.character(df$location_id)) df$location_id <- as.character(df$location_id)
  if (!is.character(gpts_map$location_id)) gpts_map$location_id <- as.character(gpts_map$location_id)

  # Trim whitespace and coerce empty->NA for stable joining
  df$location_id <- trimws(df$location_id)
  df$location_id[df$location_id == ""] <- NA_character_
  gpts_map$location_id <- trimws(gpts_map$location_id)
  gpts_map$location_id[gpts_map$location_id == ""] <- NA_character_

  # For safety, normalize simple differences (e.g. lowercase/uppercase 'L_' prefix)
  df$location_id <- ifelse(grepl('^l_', df$location_id, ignore.case = TRUE), sub('^l_', 'L_', df$location_id, ignore.case = TRUE), df$location_id)
  gpts_map$location_id <- ifelse(grepl('^l_', gpts_map$location_id, ignore.case = TRUE), sub('^l_', 'L_', gpts_map$location_id, ignore.case = TRUE), gpts_map$location_id)

  # Try joining in two ways:
  #  1) by 'location_id' (preferable), 2) by row-number (if CSV location_id are numeric row indices)
  # Ensure `Veg` exists so subsequent checks don't raise 'unknown column' warnings
  if (!"Veg" %in% names(df)) df$Veg <- NA_character_
  pre_non_na <- sum(!is.na(df$Veg) & df$Veg != "")

  joined <- dplyr::left_join(df, gpts_map, by = "location_id", suffix = c("", ".geo"))
  # Prefer existing Veg values in df; fall back to GeoJSON mapping
  if ("Veg.geo" %in% names(joined)) {
    joined$Veg <- ifelse(is.na(joined$Veg) | joined$Veg == "", joined$Veg.geo, joined$Veg)
    joined$Veg.geo <- NULL
  }
  # Merge 'no soil' from gpts_map into df when df has missing values
  if ("no soil.geo" %in% names(joined)) {
    # If df already has a 'no soil' column, prefer it when non-NA; otherwise use geo value
    if ("no soil" %in% names(joined)) {
      joined$`no soil` <- ifelse(is.na(joined$`no soil`), joined$`no soil.geo`, joined$`no soil`)
    } else {
      joined$`no soil` <- joined$`no soil.geo`
    }
    joined$`no soil.geo` <- NULL
  }

  post_non_na <- sum(!is.na(joined$Veg) & joined$Veg != "")

  # If no Veg found, try matching by geo row-number values (CSV may contain numeric row ids)
  if (post_non_na == pre_non_na && "location_row" %in% names(gpts_map)) {
    # check for intersection between CSV values and location_row
    df_ids <- unique(na.omit(as.character(df$location_id)))
    match_count <- length(intersect(df_ids, unique(na.omit(as.character(gpts_map$location_row)))))
    if (match_count > 0) {
      cat(sprintf("[NOTICE] No matches by 'location_id' — attempting join by geojson row-number mapping (matched ids=%d)\n", match_count))
      joined2 <- dplyr::left_join(df, gpts_map, by = c("location_id" = "location_row"), suffix = c("", ".geo"))
      if ("Veg.geo" %in% names(joined2)) {
        joined2$Veg <- ifelse(is.na(joined2$Veg) | joined2$Veg == "", joined2$Veg.geo, joined2$Veg)
        joined2$Veg.geo <- NULL
      }
      if ("no soil.geo" %in% names(joined2)) {
        if ("no soil" %in% names(joined2)) {
          joined2$`no soil` <- ifelse(is.na(joined2$`no soil`), joined2$`no soil.geo`, joined2$`no soil`)
        } else {
          joined2$`no soil` <- joined2$`no soil.geo`
        }
        joined2$`no soil.geo` <- NULL
      }
      if (sum(!is.na(joined2$Veg) & joined2$Veg != "") > post_non_na) {
        joined <- normalize_no_soil_col(joined2)
        post_non_na <- sum(!is.na(joined$Veg) & joined$Veg != "")
        cat(sprintf("[NOTICE] GeoJSON join (row-number) gained %d Veg rows\n", post_non_na - pre_non_na))
      } else {
        cat("[NOTICE] Row-number join did not increase Veg mapping; keeping original join state.\n")
      }
    }
  }

  df <- normalize_no_soil_col(joined)

  matched_locs <- length(intersect(na.omit(unique(as.character(df$location_id))), na.omit(unique(as.character(gpts_map$location_id)))))
  cat(sprintf("[NOTICE] GeoJSON join results - Veg before=%d after=%d; matched location_id strings=%d\n", pre_non_na, post_non_na, matched_locs))

  if (post_non_na == 0L) {
    sample_df_ids <- unique(na.omit(as.character(head(df$location_id, 20))))
    sample_geo_ids <- unique(na.omit(as.character(head(gpts_map$location_id, 20))))
    sample_geo_rows <- unique(na.omit(as.character(head(gpts_map$location_row, 20))))
    cat("[WARNING] GeoJSON join produced no Veg values. Sample df$location_id (first 20):\n")
    print(sample_df_ids)
    cat("Sample gpts_map$location_id (first 20):\n")
    print(sample_geo_ids)
    cat("Sample gpts_map row-numbers (first 20):\n")
    print(sample_geo_rows)
    cat("Hint: CSV 'location_id' might be numeric row indices; try setting location IDs in the CSV to match the GeoJSON row order or use the row-number mapping.\n")
  }
}

loc_years <- data.frame(location_id = character(0), year = integer(0), stringsAsFactors = FALSE)

if (!"year" %in% names(df)) {
  if ("date" %in% names(df)) {
    if (!requireNamespace("lubridate", quietly = TRUE)) stop("The package 'lubridate' is required")
    df$year <- as.integer(lubridate::year(as.Date(df$date)))
  }
}

if (!"Veg" %in% names(df)) df$Veg <- NA_character_

lon_candidates <- names(df)[grepl("(^|_)(lon|longitude|x)(_|$)", names(df), ignore.case = TRUE)]
lat_candidates <- names(df)[grepl("(^|_)(lat|latitude|y)(_|$)", names(df), ignore.case = TRUE)]
if (length(lon_candidates) > 0 && length(lat_candidates) > 0) {
  if ("location_id" %in% names(df) && "location_id" %in% names(gpts_map) && nrow(gpts_map) > 0) {
    # coerce consistent type before join
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
if (!"date" %in% names(df)) stop("Input CSV must contain a 'date' column")
df$date <- as.Date(df$date)
if (!"location_id" %in% names(df)) stop("Input CSV must contain a 'location_id' column")

if ("Veg" %in% names(df)) df$Veg <- tolower(as.character(df$Veg))

required_pkgs <- c("future", "future.apply", "progressr")
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

df$doy <- lubridate::yday(df$date)
df$doy[df$doy < 1 | df$doy > 366] <- NA_integer_

veg_counts <- sort(table(na.omit(df$Veg)), decreasing = TRUE)

# Add synthetic barren rows if barren locations are missing from phenology CSV
if (exists("gpts_map") && nrow(gpts_map) > 0) {
  barren_locs <- unique(gpts_map$location_id[tolower(gpts_map$Veg) == "barren" | gpts_map$`no soil` == 0])
  existing_locs <- unique(df$location_id)
  missing_barren_locs <- setdiff(barren_locs, existing_locs)
  if (length(missing_barren_locs) > 0) {
    cat(sprintf("[NOTICE] Adding %d synthetic barren rows for missing barren locations\n", length(missing_barren_locs)))
    synthetic_rows <- list()
    for (loc in missing_barren_locs) {
      # Create a base row
      new_row <- df[1, , drop = FALSE]  # copy structure
      new_row[1, ] <- NA  # clear
      new_row$location_id <- loc
      new_row$Veg <- "barren"
      new_row$`no soil` <- 0
      new_row$date <- as.Date("2000-01-01")
      new_row$year <- 2000
      new_row$doy <- 1
      # Set spectral columns to 0
      numeric_cols <- names(df)[sapply(df, is.numeric)]
      spectral_cols <- setdiff(numeric_cols, c("doy", "year"))
      for (col in spectral_cols) {
        if (col %in% names(new_row)) {
          new_row[[col]] <- 0
        }
      }
      synthetic_rows[[length(synthetic_rows) + 1]] <- new_row
    }
    if (length(synthetic_rows) > 0) {
      synthetic_df <- do.call(rbind, synthetic_rows)
      df <- rbind(df, synthetic_df)
      cat(sprintf("[NOTICE] Added %d synthetic barren rows\n", nrow(synthetic_df)))
    }
  }
}

# Define indices to use
meta_cols <- intersect(c(
  "date", "location_id", "Veg", "coverage", "lat", "lon", "latitude", "longitude",
  "target_lon", "target_lat", "imagery_lat", "imagery_lon", "doy", "year"
), names(df))
numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
candidate_indices <- intersect(OPTIMAL_INDICES, numeric_cols)
found_opt <- intersect(OPTIMAL_INDICES, numeric_cols)
missing_opt <- setdiff(OPTIMAL_INDICES, numeric_cols)
found_raw <- intersect(RAW_BANDS, numeric_cols)
missing_raw <- setdiff(RAW_BANDS, numeric_cols)

if (length(found_opt) > 0) {
  cat(sprintf("Found %d OPTIMAL_INDICES in input: %s\n", length(found_opt), paste(found_opt, collapse = ", ")))
}
if (length(missing_opt) > 0) {
  cat(sprintf("Missing OPTIMAL_INDICES in input: %s\n", paste(missing_opt, collapse = ", ")))
}
if (length(found_raw) > 0) {
  cat(sprintf("Found %d RAW_BANDS in input: %s\n", length(found_raw), paste(found_raw, collapse = ", ")))
}
if (length(missing_raw) > 0) {
  cat(sprintf("Missing RAW_BANDS in input: %s\n", paste(missing_raw, collapse = ", ")))
}

candidate_indices <- c(found_opt, found_raw)
if (length(candidate_indices) == 0) {
  stop("No OPTIMAL_INDICES present in input CSV")
}

# Global amplitude-based index selection
cat("Calculating global index amplitudes...\n")
index_amplitudes <- vapply(candidate_indices, function(idx) {
  vals <- df[[idx]][is.finite(df[[idx]])]
  if (length(vals) < 10) {
    return(0)
  }
  amp <- diff(range(vals, na.rm = TRUE))
  if (!is.finite(amp)) 0 else amp
}, numeric(1))

amplitude_threshold <- 0.05
high_amplitude_indices <- names(index_amplitudes)[index_amplitudes > amplitude_threshold]
low_amplitude_indices <- names(index_amplitudes)[index_amplitudes <= amplitude_threshold]

if (length(low_amplitude_indices) > 0) {
  cat(sprintf(
    "Filtered out %d low-amplitude indices: %s\n",
    length(low_amplitude_indices), paste(low_amplitude_indices, collapse = ", ")
  ))
}

avail <- high_amplitude_indices

# Correlation-based filtering
if (length(avail) > 1) {
  cat("Applying correlation-based filtering...\n")

  corr_data <- df[, avail, drop = FALSE]
  corr_matrix <- cor(corr_data, use = "pairwise.complete.obs")
  corr_matrix[is.na(corr_matrix)] <- 0

  high_corr_threshold <- 0.90
  to_remove <- c()
  if (ncol(corr_matrix) >= 2) {
    for (i in seq_len(ncol(corr_matrix) - 1)) {
      for (j in seq((i + 1), ncol(corr_matrix))) {
        if (abs(corr_matrix[i, j]) > high_corr_threshold) {
          idx1 <- colnames(corr_matrix)[i]
          idx2 <- colnames(corr_matrix)[j]
          amp1 <- index_amplitudes[idx1]
          amp2 <- index_amplitudes[idx2]

          if (abs(amp1 - amp2) / max(amp1, amp2) < 0.01) {
            to_remove <- c(to_remove, idx2)
          } else {
            if (amp1 < amp2) {
              to_remove <- c(to_remove, idx1)
            } else {
              to_remove <- c(to_remove, idx2)
            }
          }
        }
      }
    }
  }

  to_remove <- unique(to_remove)
  if (length(to_remove) > 0) {
    avail <- setdiff(avail, to_remove)
    cat(sprintf(
      "Removed %d highly correlated indices: %s\n",
      length(to_remove), paste(to_remove, collapse = ", ")
    ))
  }
}

if (length(avail) == 0) {
  stop("No indices found with amplitude > 0.05 after correlation filtering")
}

if (length(avail) < USE_INDICES_MIN) {
  stop(sprintf(
    "Only %d indices remain after filtering, minimum required is %d",
    length(avail), USE_INDICES_MIN
  ))
}

cat(sprintf("Selected %d indices: %s\n", length(avail), paste(avail, collapse = ", ")))

  # Data sufficiency check DISABLED per user request - all locations are kept regardless of DOY coverage
  
  if (FALSE && (is.null(TRAIN_YEARS) || length(TRAIN_YEARS) == 0)) stop("TRAIN_YEARS must be defined and non-empty for data sufficiency checks; no fallback permitted.")
  if (FALSE && "location_id" %in% names(df_full)) {
    # For each location and each selected index, count finite observations and error if below threshold
    offenders <- list()
    for (loc in unique(df_full$location_id)) {
      sub <- df_full[df_full$location_id == loc, , drop = FALSE]
      # Determine number of training years present for this location
      loc_years <- unique(sub$year[is.finite(sub$year) & !is.na(sub$year)])
      train_loc_years <- intersect(as.integer(loc_years), as.integer(TRAIN_YEARS))
      n_years <- length(train_loc_years)
      if (n_years == 0) {
        stop(sprintf("Data sufficiency check failed: location '%s' has no observations in TRAIN_YEARS (%s) — cannot compute required DOY threshold; aborting (no fallback permitted).", loc, paste(TRAIN_YEARS, collapse = ", ")))
      }

      # required days = 50% of all possible DOYs across the location's training years
      # Allow skipping a minimal number of DOYs per-location (controlled by MIN_SKIP_DOYS_PER_LOCATION)
      base_required <- as.integer(ceiling(n_years * 365 / 2))
      required_days_for_loc <- as.integer(max(1L, base_required - as.integer(MIN_SKIP_DOYS_PER_LOCATION)))

      # Count unique DOYs in the training years where at least one of the selected indices is finite
      if (nrow(sub) == 0) {
        n_doys_present <- 0L
      } else {
        sub_tr <- sub[sub$year %in% train_loc_years, , drop = FALSE]
        if (nrow(sub_tr) == 0) {
          n_doys_present <- 0L
        } else {
          # Mark rows where any of the selected indices are finite
          has_any_index <- apply(sub_tr[, intersect(avail, names(sub_tr)), drop = FALSE], 1, function(r) any(is.finite(r)))
          # Count unique dates (or doy) present with any index
          if ('date' %in% names(sub_tr)) {
            n_doys_present <- length(unique(as.character(sub_tr$date[which(has_any_index)])))
          } else if ('doy' %in% names(sub_tr)) {
            n_doys_present <- length(unique(sub_tr$doy[which(has_any_index)]))
          } else {
            # fallback to counting rows with any index finite
            n_doys_present <- sum(has_any_index)
          }
        }
      }

      if (n_doys_present < required_days_for_loc) {
        offenders[[length(offenders) + 1]] <- list(location = loc, n = n_doys_present, required = required_days_for_loc, years = paste(train_loc_years, collapse = ","))
      }
    }
    if (length(offenders) > 0) {
      # Build detailed messages per offending location (aggregate across indices)
      msg_lines <- vapply(offenders, function(x) sprintf("loc=%s n=%d required=%d train_years=%s", x$location, x$n, x$required, x$years), character(1))
      bad_locs <- unique(vapply(offenders, function(x) as.character(x$location), character(1)))

      warning(sprintf(
        "Data sufficiency: skipping %d locations because they do not meet the per-location DOY threshold (50%% of DOYs across training years) after allowing up to %d missing DOYs per-location.\nSkipped locations: %s\nDetails (loc,obs,required,train_years):\n%s",
        length(bad_locs), as.integer(MIN_SKIP_DOYS_PER_LOCATION), paste(head(bad_locs, 100), collapse = ", "), paste(msg_lines, collapse = "\n")
      ))

      
      # we flag them instead so downstream code can decide how to handle them.
      flag_col <- ".insufficient_data_by_DOY_threshold"
      if (!flag_col %in% names(df_full)) df_full[[flag_col]] <- FALSE
      df_full[[flag_col]][df_full$location_id %in% bad_locs] <- TRUE
      if (exists("df")) {
        if (!flag_col %in% names(df)) df[[flag_col]] <- FALSE
        df[[flag_col]][df$location_id %in% bad_locs] <- TRUE
      }
      if (exists("df_train")) {
        if (!flag_col %in% names(df_train)) df_train[[flag_col]] <- FALSE
        df_train[[flag_col]][df_train$location_id %in% bad_locs] <- TRUE
      }
      if (exists("df_test")) {
        if (!flag_col %in% names(df_test)) df_test[[flag_col]] <- FALSE
        df_test[[flag_col]][df_test$location_id %in% bad_locs] <- TRUE
      }

      # Keep all locations but notify if this leaves nothing usable
      if (nrow(df_full) == 0 || length(unique(na.omit(df_full$location_id))) == 0) {
        stop("Data sufficiency check left no usable data — aborting.")
      }
    }
  } else {
    # Global time series: require TRAIN_YEARS and enforce proportional check globally
    if (is.null(TRAIN_YEARS) || length(TRAIN_YEARS) == 0) stop("TRAIN_YEARS must be defined for global sufficiency checks; no fallback permitted.")
    base_req_global <- as.integer(ceiling(length(TRAIN_YEARS) * 365 / 2))
    required_days_global <- as.integer(max(1L, base_req_global - as.integer(MIN_SKIP_DOYS_PER_LOCATION)))
    if (nrow(df_full) < required_days_global) stop(sprintf("Data sufficiency check failed: overall data contains fewer than %d rows (50%% of days across TRAIN_YEARS after allowing up to %d missing DOYs)", required_days_global, MIN_SKIP_DOYS_PER_LOCATION))
  }

#
if (length(avail) > 0) {
  pb <- utils::txtProgressBar(min = 0, max = length(avail), style = 3)
  i <- 0
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

    i <- i + 1
    try(utils::setTxtProgressBar(pb, i), silent = TRUE)
  }
  try(close(pb), silent = TRUE)
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

# Post-variance cleaning
var_cols <- names(df)[grepl("(_var14$|_mv$)", names(df))]
if (length(var_cols) > 0) {
  cat(sprintf("Post-variance cleaning: %d variance columns detected\n", length(var_cols)))

  ac1 <- numeric(length(var_cols))
  names(ac1) <- var_cols
  cv <- numeric(length(var_cols))
  names(cv) <- var_cols
  finite_frac <- numeric(length(var_cols))
  names(finite_frac) <- var_cols

  for (i in seq_along(var_cols)) {
    v <- df[[var_cols[i]]]
    v_f <- v[is.finite(v)]
    finite_frac[i] <- if (length(v) == 0) 0 else sum(is.finite(v)) / length(v)
    if (length(v_f) > 3) {
      ac1[i] <- tryCatch(
        {
          cor(v_f[-1], v_f[-length(v_f)], use = "complete.obs")
        },
        error = function(e) stop(sprintf("Post-variance cleaning: failed to compute lag-1 autocorrelation for %s: %s", var_cols[i], e$message))
      )
      mm <- median(v_f, na.rm = TRUE)
      msd <- stats::sd(v_f, na.rm = TRUE)
      cv[i] <- if (is.finite(mm) && mm != 0) msd / abs(mm) else Inf
    } else {
      ac1[i] <- NA_real_
      cv[i] <- NA_real_
    }
  }

  keep_mask <- (finite_frac >= 0.30) & (is.finite(ac1) & ac1 > 0.10 | is.finite(cv) & cv > 0.05)
  remove_rand <- var_cols[!keep_mask]
  if (length(remove_rand) > 0) {
    cat(sprintf("Removing %d random-looking variance cols: %s\n", length(remove_rand), paste(remove_rand, collapse = ", ")))
    df[remove_rand] <- NULL
    var_cols <- setdiff(var_cols, remove_rand)
  }

if (length(var_cols) > 1) {
    vmat <- as.data.frame(lapply(var_cols, function(nm) df[[nm]]))
    names(vmat) <- var_cols
    
    # Check for non-finite values (store results with print)
    print(paste("NAs:", sum(is.na(vmat))))           # Count NAs
    print(paste("NaNs:", sum(is.nan(as.matrix(vmat)))))          # Count NaNs (convert to matrix first)
    print(paste("Infs:", sum(is.infinite(as.matrix(vmat)))))     # Count Inf/-Inf

    # Check proportion of missing data per column
    na_summary <- data.frame(
      column = colnames(vmat),
      na_count = colSums(is.na(vmat)),
      na_prop = colMeans(is.na(vmat)))
    print(na_summary[order(-na_summary$na_prop), ])

    # Check for infinite values
    inf_summary <- colSums(is.infinite(as.matrix(vmat)))
    print(inf_summary[inf_summary > 0])
    
    # Clean the data
    # First, replace infinite values with NA
    vmat_clean <- vmat
    vmat_clean[sapply(vmat_clean, is.infinite)] <- NA
    
    # Remove columns with too many NAs (e.g., >50%)
    na_prop <- colMeans(is.na(vmat_clean))
    vmat_clean <- vmat_clean[, na_prop < 0.5, drop = FALSE]
    
    # Remove rows with any remaining NAs
    vmat_clean <- vmat_clean[complete.cases(vmat_clean), , drop = FALSE]
    
    # Check if we have enough data left
    if (nrow(vmat_clean) < 2 || ncol(vmat_clean) < 2) {
      stop("Not enough clean data remaining after filtering")
    }
    
    # Convert to matrix before calling cor
    vmat_matrix <- as.matrix(vmat_clean)
    
    # Then use the cleaned matrix
    cm <- cor(vmat_matrix, use = "pairwise.complete.obs")
    cm[is.na(cm)] <- 0
    
    high_corr_thresh <- 0.95
    to_remove <- c()
    for (i in seq_len(ncol(cm) - 1)) {
      for (j in (i + 1):ncol(cm)) {
        if (abs(cm[i, j]) > high_corr_thresh) {
          c1 <- var_cols[i]
          c2 <- var_cols[j]
          f1 <- sum(is.finite(df[[c1]]))
          f2 <- sum(is.finite(df[[c2]]))
          if (f1 >= f2) to_remove <- c(to_remove, c2) else to_remove <- c(to_remove, c1)
        }
      }
    }

    to_remove <- unique(to_remove)
    if (length(to_remove) > 0) {
      cat(sprintf("Removing %d highly correlated variance cols: %s\n", length(to_remove), paste(to_remove, collapse = ", ")))
      df[to_remove] <- NULL
    }
  }
}

cat(sprintf("Post baseline normalization rows=%d\n", nrow(df)))
cat("Data preprocessing complete.\n")
adj_cols <- intersect(avail, names(df))

# Simple feature pruning
if (length(avail) > 0) {
  idx_sd <- vapply(avail, function(nm) {
    x <- df[[nm]]
    x <- x[is.finite(x)]
    if (length(x) < 5) {
      return(0)
    }
    stats::sd(x, na.rm = TRUE)
  }, numeric(1))
  keep <- idx_sd >= MIN_INDEX_SD
  dropped <- avail[!keep]
  if (length(dropped) > 0) cat("Pruned low-variance indices:", paste(dropped, collapse = ", "), "\n")
  avail <- avail[keep]
}

# Restrict dataset to allowed vegetation classes
if ("Veg" %in% names(df) && length(ALLOWED_VEG) > 0) {
  try(
    {
      for (av in ALLOWED_VEG) {
        sel <- grepl(av, df$Veg, ignore.case = TRUE) & !is.na(df$Veg)
        if (any(sel)) {
          df$Veg[sel] <- av
          cat(sprintf("Normalized %d rows to canonical Veg '%s'\n", sum(sel, na.rm = TRUE), av))
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

# Per-vegetation quick summary
try(
  {
    cat("Per-veg quick summary:\n")
    for (av in ALLOWED_VEG) {
      sel <- tolower(df$Veg) == av
      rows <- sum(sel, na.rm = TRUE)
      unique_doys <- length(unique(df$doy[sel & is.finite(df$doy)]))
      cat(sprintf("  %s: rows=%d unique_doys=%d\n", av, rows, unique_doys))
    }
  },
  silent = TRUE
)

matched_veg_n <- sum(!is.na(df$Veg))
cat("Non-NA Veg rows:", matched_veg_n, "of", nrow(df), "\n")
if (matched_veg_n == 0) {
  stop("No vegetation classes found after join; cannot build library")
}

loc_years <- df %>%
  dplyr::filter(!is.na(.data$location_id) & .data$location_id != "" & !is.na(.data$year) & .data$year > 0) %>%
  dplyr::distinct(.data$location_id, .data$year)
cat(sprintf("Constructed loc_years with %d rows from filtered df\n", nrow(loc_years)))
if (nrow(loc_years) == 0) {
  # Hard fail: do not fall back silently to the full dataset. Training data
  # must contain at least one location-year pair after filtering.
  stop(paste0(
    "No location-year pairs found after filtering. This is a fatal error — library construction cannot continue.\n",
    "Possible causes and suggestions:\n",
    " - Your filtered training dataset has zero valid location/year pairs (check column 'location_id' and 'year').\n",
    " - Verify TRAIN_YEARS and any prior filtering steps do not remove all data (e.g. TRAIN_YEARS <- 2019:2024).\n",
    " - Ensure your transformer produced valid 'location_id' values that match your geojson mapping (transform_phenology.py formats 'L_lon_lat').\n",
    " - If your data are intentionally sparse, reduce the filtering thresholds or increase available training data.\n",
    "Processing cannot continue without at least one location-year pair in filtered training data.")
  )
}

# Construct vegetation library
cat("Constructing lib from TRAINING dataset...\n")
lib <- list()
## Soil-first MESMA preprocessing: build a soil prototype from geojson-labelled
## 'no soil' points and subtract estimated soil fraction from every observation
## in the working data. This runs before vegetation library construction.
build_soil_prototype <- function(df_local, avail_idx, threshold = SOIL_PURE_THRESHOLD, min_samples = SOIL_MIN_SAMPLES) {
  # Prefer using explicit 'barren' Veg-labeled rows as soil (bare ground) references.
  # If no 'barren' rows are present, do not attempt to create a soil prototype.
  if (!"Veg" %in% names(df_local)) return(NULL)
  candidates <- df_local[tolower(as.character(df_local$Veg)) == "barren", , drop = FALSE]

  ## Strict requirement: do not fallback to selecting an arbitrary top-N set.
  ## If there are fewer than min_samples candidates that meet the threshold,
  ## return NULL and let the caller handle the absence of a soil prototype.
  if (nrow(candidates) < min_samples) {
    return(NULL)
  }

  if (nrow(candidates) == 0) return(NULL)

  soil_lib <- list()
  # ensure 'doy' is present
  if (!"doy" %in% names(candidates)) candidates$doy <- lubridate::yday(as.Date(candidates$date))

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

## Build a two-endmember library for stage 1 MESMA: barren (no soil==0) vs pure vegetation (no soil==1)
## This library is used to unmix the vegetated fraction for all pixels before decomposing vegetation types
build_barren_veg_library <- function(df_local, avail_idx, min_samples = 5) {
  if (!"no soil" %in% names(df_local)) return(NULL)
  if (!"doy" %in% names(df_local)) df_local$doy <- lubridate::yday(as.Date(df_local$date))
  
  # Debug: Check what values are in the 'no soil' column
  no_soil_vals <- df_local$`no soil`
  cat(sprintf("[Stage1] 'no soil' column summary:\n"))
  cat(sprintf("  - Class: %s\n", class(no_soil_vals)))
  cat(sprintf("  - Unique values (first 10): %s\n", paste(head(unique(no_soil_vals), 10), collapse=", ")))
  cat(sprintf("  - Range: %s to %s\n", min(no_soil_vals, na.rm=TRUE), max(no_soil_vals, na.rm=TRUE)))
  cat(sprintf("  - NAs: %d\n", sum(is.na(no_soil_vals))))
  
  # Extract barren endmember (no soil ≈ 0)
  # Be more robust: handle strings, very small values, and different formats
  barren_rows <- df_local[!is.na(df_local$`no soil`) & {
    val <- df_local$`no soil`
    # Convert to numeric if it's a string
    if (is.character(val)) val <- as.numeric(val)
    # Check if close to 0 (within 0.01) or exactly 0
    abs(val - 0) < 0.01
  }, , drop = FALSE]
  
  # Extract pure vegetation endmember (no soil ≈ 1)
  veg_rows <- df_local[!is.na(df_local$`no soil`) & {
    val <- df_local$`no soil`
    # Convert to numeric if it's a string
    if (is.character(val)) val <- as.numeric(val)
    # Check if close to 1 (within 0.01) or exactly 1
    abs(val - 1) < 0.01
  }, , drop = FALSE]
  
  cat(sprintf("[Stage1] After filtering: barren=%d, veg=%d rows\n", 
              nrow(barren_rows), nrow(veg_rows)))
  
  if (nrow(barren_rows) < min_samples || nrow(veg_rows) < min_samples) {
    cat(sprintf("[Stage1] Insufficient training data: barren=%d, veg=%d (need >=%d each)\n", 
                nrow(barren_rows), nrow(veg_rows), min_samples))
    return(NULL)
  }
  
  cat(sprintf("[Stage1] Building barren-veg library from barren=%d, veg=%d rows\n", 
              nrow(barren_rows), nrow(veg_rows)))
  
  stage1_lib <- list()
  
  # Build barren endmember
  barren_lib <- list()
  for (idx in avail_idx) {
    if (!idx %in% names(barren_rows)) next
    vals_by_doy <- tapply(seq_along(barren_rows[[idx]]), barren_rows$doy, function(indices) {
      vals <- barren_rows[[idx]][indices]
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) return(NA_real_)
      median(vals)
    })
    mu <- rep(NA_real_, 365)
    if (length(vals_by_doy) > 0) {
      doy_values <- as.integer(names(vals_by_doy))
      valid_doy <- doy_values >= 1 & doy_values <= 365
      mu[doy_values[valid_doy]] <- vals_by_doy[valid_doy]
    }
    if (all(!is.finite(mu))) next
    mv <- tryCatch(calc_moving_var(data.frame(date = 1:365, idx = mu), "idx", window = 14), 
                   error = function(e) rep(NA_real_, 365))
    barren_lib[[idx]] <- list(mu = mu, mv = mv)
  }
  
  # Build vegetation endmember
  veg_lib <- list()
  for (idx in avail_idx) {
    if (!idx %in% names(veg_rows)) next
    vals_by_doy <- tapply(seq_along(veg_rows[[idx]]), veg_rows$doy, function(indices) {
      vals <- veg_rows[[idx]][indices]
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) return(NA_real_)
      median(vals)
    })
    mu <- rep(NA_real_, 365)
    if (length(vals_by_doy) > 0) {
      doy_values <- as.integer(names(vals_by_doy))
      valid_doy <- doy_values >= 1 & doy_values <= 365
      mu[doy_values[valid_doy]] <- vals_by_doy[valid_doy]
    }
    # Interpolate missing DOYs for complete temporal coverage
    if (sum(is.finite(mu)) > 2) {
      mu <- approx(x = which(is.finite(mu)), y = mu[is.finite(mu)], 
                   xout = 1:365, rule = 2)$y
    }
    if (all(!is.finite(mu))) next
    mv <- tryCatch(calc_moving_var(data.frame(date = 1:365, idx = mu), "idx", window = 14), 
                   error = function(e) rep(NA_real_, 365))
    veg_lib[[idx]] <- list(mu = mu, mv = mv)
  }
  
  if (length(barren_lib) == 0 || length(veg_lib) == 0) return(NULL)
  
  stage1_lib$barren <- barren_lib
  stage1_lib$vegetation <- veg_lib
  stage1_lib
}

## Unmix vegetated fraction using stage 1 MESMA (barren vs vegetation endmembers)
## Returns the estimated vegetated fraction for this location-year
unmix_vegetated_fraction <- function(dly_local, stage1_lib, avail_idx) {
  if (is.null(stage1_lib) || length(stage1_lib) == 0) return(NA_real_)
  if (!"doy" %in% names(dly_local)) dly_local$doy <- lubridate::yday(as.Date(dly_local$date))
  
  # For each row, perform constrained least squares: pixel = alpha*barren + (1-alpha)*veg
  # We solve: min ||X - alpha*B - (1-alpha)*V||^2 subject to alpha in [0,1]
  # This simplifies to: min ||X - V - alpha*(B - V)||^2 --> alpha = (X-V)·(B-V) / ||B-V||^2, clipped to [0,1]
  
  alphas <- numeric(nrow(dly_local))
  
  for (r in seq_len(nrow(dly_local))) {
    row <- dly_local[r, , drop = FALSE]
    doy <- as.integer(row$doy)
    if (is.na(doy) || doy < 1 || doy > 365) {
      alphas[r] <- NA_real_
      next
    }
    
    # Build vectors for this DOY
    X <- numeric(length(avail_idx))
    B <- numeric(length(avail_idx))
    V <- numeric(length(avail_idx))
    valid <- logical(length(avail_idx))
    
    for (i in seq_along(avail_idx)) {
      idx <- avail_idx[i]
      if (!idx %in% names(row)) next
      X[i] <- as.numeric(row[[idx]])
      if (!is.null(stage1_lib$barren[[idx]])) B[i] <- stage1_lib$barren[[idx]]$mu[doy]
      if (!is.null(stage1_lib$vegetation[[idx]])) V[i] <- stage1_lib$vegetation[[idx]]$mu[doy]
      valid[i] <- is.finite(X[i]) && is.finite(B[i]) && is.finite(V[i])
    }
    
    if (sum(valid) < 2) {
      alphas[r] <- NA_real_
      next
    }
    
    X <- X[valid]
    B <- B[valid]
    V <- V[valid]
    
    # Solve for alpha (barren fraction)
    diff_BV <- B - V
    diff_XV <- X - V
    denom <- sum(diff_BV^2)
    
    if (denom < 1e-10) {
      # Endmembers are too similar, default to mid-range
      alphas[r] <- 0.5
    } else {
      alpha <- sum(diff_XV * diff_BV) / denom
      alphas[r] <- max(0, min(1, alpha))
    }
  }
  
  # Return median vegetated fraction across all valid observations
  valid_alphas <- alphas[is.finite(alphas)]
  if (length(valid_alphas) == 0) return(NA_real_)
  
  barren_frac <- median(valid_alphas)
  veg_frac <- 1 - barren_frac
  
  cat(sprintf("[Stage1] Unmixed vegetated fraction: %.3f (barren: %.3f, n=%d obs)\n", 
              veg_frac, barren_frac, length(valid_alphas)))
  
  veg_frac
}

 


## Build a provisional global PCA from raw data (no averaging)
## This ensures the soil subtraction logic operates in a space that captures full dataset variance.
build_prelim_global_pca <- function(df_local, avail_idx) {
  # REVISED: Build PCA from raw data, not averaged location-year medoids.
  
  # Filter valid data
  # We need columns: avail_idx
  missing_cols <- setdiff(avail_idx, names(df_local))
  if (length(missing_cols) > 0) return(NULL)
  
  # Extract matrix
  X_all <- as.matrix(df_local[, avail_idx, drop = FALSE])
  
  # Handle NAs
  # Simple mean imputation for PCA training
  for (j in seq_len(ncol(X_all))) {
    col_vals <- X_all[, j]
    if (any(!is.finite(col_vals))) {
      mu_j <- mean(col_vals[is.finite(col_vals)], na.rm = TRUE)
      if (!is.finite(mu_j)) mu_j <- 0
      X_all[!is.finite(col_vals), j] <- mu_j
    }
  }
  
  # PCA
  mu_all <- colMeans(X_all)
  Xc <- sweep(X_all, 2, mu_all, "-")
  feature_sds <- apply(Xc, 2, sd)
  feature_sds[feature_sds <= 1e-10] <- 1.0
  
  Xs_std <- sweep(Xc, 2, feature_sds, "/")
  
  sv <- try(svd(Xs_std, nu = 0, nv = min(ncol(Xs_std), nrow(Xs_std))), silent = TRUE)
  if (inherits(sv, "try-error") || length(sv$d) == 0) return(NULL)
  
  eigenvalues <- sv$d^2 / (nrow(Xs_std) - 1)
  keep_idx <- which(eigenvalues > 1)
  pca_rank <- if (length(keep_idx) == 0) 1 else min(length(keep_idx), MAX_FACTORS_CAP, ncol(sv$v))
  
  V_pca <- sv$v[, seq_len(pca_rank), drop = FALSE]
  if (ncol(V_pca) > 1) V_pca <- qr.Q(qr(V_pca))[, seq_len(ncol(V_pca)), drop = FALSE]
  
  list(
    global = list(V = V_pca, col_means = mu_all, feature_sds = feature_sds, rank = ncol(V_pca)),
    idx_order = avail_idx
  )
}


## Using a provisional gpca, compute a single alpha per location-year using
## the factor projections for all rows and the soil prototype projected into the
## same factor space (s_t per DOY). Returns df updated with soil_frac per-row
## (constant across location-year) and subtracts alpha * soil_mu per-row.
subtract_soil_by_location_year <- function(df_in, soil_lib, avail_idx, loc_year_pairs, gpca) {
  if (is.null(soil_lib) || length(soil_lib) == 0) return(df_in)
  if (is.null(gpca)) return(NULL)

  # Precompute soil factors for each DOY
  soil_factors <- vector("list", 366)
  names(soil_factors) <- as.character(0:365)
  for (doy in 1:365) {
    raw_vec <- rep(NA_real_, length(gpca$idx_order))
    names(raw_vec) <- gpca$idx_order
    for (idx in avail_idx) {
      kpos <- match(idx, gpca$idx_order)
      if (!is.na(kpos) && !is.null(soil_lib[[idx]])) raw_vec[kpos] <- soil_lib[[idx]]$mu[doy]
    }
    if (all(!is.finite(raw_vec))) next
    # center and scale
    col_means <- gpca$global$col_means
    feature_sds <- gpca$global$feature_sds
    nas <- which(!is.finite(raw_vec))
    if (length(nas) > 0) raw_vec[nas] <- col_means[nas]
    centered <- (raw_vec - col_means) / feature_sds
    s_factor <- as.numeric(centered %*% gpca$global$V)
    soil_factors[[as.character(doy)]] <- s_factor
  }

  df_local <- df_in
  if (!"doy" %in% names(df_local)) df_local$doy <- lubridate::yday(df_local$date)
  # map for alphas
  alpha_map <- list()

  for (i in seq_len(nrow(loc_year_pairs))) {
    loc <- as.character(loc_year_pairs$location_id[i])
    yr <- as.integer(loc_year_pairs$year[i])
    sel_rows <- which(df_local$location_id == loc & df_local$year == yr)
    if (length(sel_rows) == 0) next

    Zs <- list(); Ss <- list()
    for (r in sel_rows) {
      row <- df_local[r, , drop = FALSE]
      # project observation into factor space using gpca
      z <- tryCatch(project_row_to_factors(row, gpca, avail_idx, veg_type = row$Veg), error = function(e) NULL)
      if (is.null(z)) next
      doy <- as.integer(row$doy)
      if (is.na(doy) || doy < 1 || doy > 365) next
      s_t <- soil_factors[[as.character(doy)]]
      if (is.null(s_t) || !is.finite(sum(s_t))) next
      Zs[[length(Zs) + 1]] <- z
      Ss[[length(Ss) + 1]] <- s_t
    }

    if (length(Zs) < 2) next
    Zm <- do.call(rbind, Zs)
    Sm <- do.call(rbind, Ss)
    # Solve scalar alpha minimizing ||Zm - alpha * Sm||_F^2 --> alpha = sum(vec(Zm) dot vec(Sm)) / sum(vec(Sm)^2)
    numer <- sum(Zm * Sm, na.rm = TRUE)
    denom <- sum(Sm * Sm, na.rm = TRUE)
    if (!is.finite(denom) || denom <= 0) next
    alpha <- numer / denom
    if (!is.finite(alpha)) next
    alpha <- max(0, min(1, alpha))
    key <- paste0(loc, "__", yr)
    alpha_map[[key]] <- alpha
  }

  # Apply the computed alphas (per location-year) to subtract from index space
  for (r in seq_len(nrow(df_local))) {
    row <- df_local[r, , drop = FALSE]
    loc <- as.character(row$location_id)
    yr <- as.integer(row$year)
    key <- paste0(loc, "__", yr)
    if (!key %in% names(alpha_map)) next
    alpha <- alpha_map[[key]]
    if (!is.finite(alpha)) next

    doy <- as.integer(row$doy)
    if (is.na(doy) || doy < 1 || doy > 365) next

    for (idx in avail_idx) {
      if (!idx %in% names(df_local)) next
      raw_col <- paste0("raw_", idx)
      if (!raw_col %in% names(df_local)) df_local[[raw_col]] <- NA_real_
      df_local[[raw_col]][r] <- df_local[[idx]][r]
      s_val <- if (!is.null(soil_lib[[idx]])) soil_lib[[idx]]$mu[doy] else NA_real_
      if (is.finite(s_val)) df_local[[idx]][r] <- df_local[[idx]][r] - alpha * s_val
    }
    df_local$soil_frac[r] <- alpha
  }

  df_local
}

# Build stage 1 library for vegetated fraction unmixing (barren vs pure vegetation)
STAGE1_LIB <- NULL
try({
  STAGE1_LIB <- build_barren_veg_library(df, avail, min_samples = 5)
}, silent = TRUE)

if (is.null(STAGE1_LIB)) {
  stop("[Stage1] Could not build barren-veg library (insufficient no soil==0 or no soil==1 rows). Two-stage MESMA requires sufficient stage-1 training data. Cannot proceed with spectral-only unmixing.")
} else {
  cat("[Stage1] Barren-vegetation library built successfully for stage 1 unmixing\n")
}

if (FALSE && isTRUE(ENABLE_SOIL_PREPROCESS)) {
  cat("Running soil-preprocessing: building soil prototype and subtracting soil contributions from signals...\n")
  soil_proto <- NULL
  try({
    soil_proto <- build_soil_prototype(df, avail)
  }, silent = TRUE)

  if (is.null(soil_proto) || length(soil_proto) == 0) {
    cat("[NOTICE] Soil prototype could not be constructed (insufficient/no 'no soil' rows) — skipping soil subtraction\n")
  } else {
    cat(sprintf("Soil prototype created from %d geo rows; estimating constant soil fractions per location-year using provisional factor projection...\n", attr(soil_proto, "rows_used")))

    gpca <- tryCatch({ build_prelim_global_pca(df, avail) }, error = function(e) NULL)
    if (!is.null(gpca)) {
      cat("Provisional global PCA built — computing per-location-year alpha and subtracting soil (constant across time-dimension)\n")
      new_df <- tryCatch({ subtract_soil_by_location_year(df, soil_proto, avail, loc_years, gpca) }, error = function(e) NULL)
      if (!is.null(new_df)) {
        df <- new_df
        df_full <- tryCatch({ subtract_soil_by_location_year(df_full, soil_proto, avail, loc_years, gpca) }, error = function(e) df_full)
        if (exists("df_test") && nrow(df_test) > 0) df_test <- tryCatch({ subtract_soil_by_location_year(df_test, soil_proto, avail, loc_years, gpca) }, error = function(e) df_test)
      } else {
        cat("[NOTICE] Per-location-year subtraction failed. Soil subtraction will be skipped (no fallback)\n")
      }
    } else {
      cat("[NOTICE] Could not build provisional PCA (not enough data). Soil subtraction will be skipped (no fallback)\n")
    }
  }
}

lib_df <- df
vegs <- unique(na.omit(lib_df$Veg))
vegs <- vegs[vegs != ""]

# Sample Size Balancing
cat("Applying sample size balancing to training data...\n")
lib_df$doy_for_lib <- lubridate::yday(lib_df$date)

veg_sample_sizes <- sapply(vegs, function(v) {
  dveg <- lib_df[lib_df$Veg == v & is.finite(lib_df$doy_for_lib), , drop = FALSE]
  nrow(dveg)
})

min_samples <- min(veg_sample_sizes)
max_samples <- max(veg_sample_sizes)
median_samples <- median(veg_sample_sizes)

target_size <- pmax(min_samples, floor(median_samples * 0.8))
target_size <- pmin(target_size, max_samples)

cat(sprintf(
  "Sample sizes - Min: %.0f, Max: %.0f, Median: %.1f\n",
  min_samples, max_samples, median_samples
))
cat(sprintf("Balancing to target size: %.0f\n", target_size))

# SMOTE-like augmentation for minority classes
augment_minority_class <- function(df_class, target_n, seed = NULL, alpha_range = c(0.3, 0.7), jitter_frac = 1e-6) {
  if (!is.data.frame(df_class)) stop("augment_minority_class requires a data.frame")
  if (nrow(df_class) >= target_n) {
    return(df_class)
  }
  if (nrow(df_class) <= 1) {
    reps <- target_n - nrow(df_class)
    base <- df_class[rep(1, reps), , drop = FALSE]
    num_cols <- sapply(base, is.numeric)
    if (any(num_cols)) {
      rng <- apply(df_class[, num_cols, drop = FALSE], 2, function(x) {
        if (all(!is.finite(x))) 0 else diff(range(x, na.rm = TRUE))
      })
      rng[rng == 0] <- 1.0
      base[, num_cols] <- base[, num_cols, drop = FALSE] +
        matrix(rnorm(reps * sum(num_cols), sd = sqrt(mean(rng)) * jitter_frac), nrow = reps)
    }
    out <- rbind(df_class, base)
    rownames(out) <- NULL
    return(out)
  }

  if (!is.null(seed)) set.seed(seed)
  n_add <- target_n - nrow(df_class)
  synth_rows <- vector("list", n_add)

  num_cols <- names(df_class)[sapply(df_class, is.numeric)]
  other_cols <- setdiff(names(df_class), num_cols)

  for (i in seq_len(n_add)) {
    ids <- sample.int(nrow(df_class), 2, replace = FALSE)
    a <- df_class[ids[1], , drop = FALSE]
    b <- df_class[ids[2], , drop = FALSE]
    alpha <- runif(1, min(alpha_range), max(alpha_range))
    new_row <- a

    if (length(num_cols) > 0) {
      va <- as.numeric(a[1, num_cols, drop = FALSE])
      vb <- as.numeric(b[1, num_cols, drop = FALSE])
      interp <- (1 - alpha) * va + alpha * vb
      rng <- apply(df_class[, num_cols, drop = FALSE], 2, function(x) {
        if (all(!is.finite(x))) 0 else diff(range(x, na.rm = TRUE))
      })
      jitter <- rnorm(length(interp), sd = pmax(abs(rng), 1e-8) * jitter_frac)
      interp <- interp + jitter
      # Assign numeric interpolated values element-wise to avoid tibble/vctrs
      # assignment errors when new_row is a tibble-like object
      if (length(num_cols) == length(interp)) {
        for (k in seq_along(num_cols)) {
          colname <- num_cols[k]
          new_row[[colname]] <- interp[k]
        }
      } else {
        # Do not silently recycle/interpolate when shapes don't match: this
        # hides subtle bugs in numeric column handling. Fail loudly so calling
        # code can be fixed instead of producing silently corrupted samples.
        stop(sprintf("augment_minority_class: numeric assignment length mismatch (num_cols=%d, interp_len=%d)", length(num_cols), length(interp)))
      }
    }

    for (col in other_cols) {
      if (col %in% c("date", "doy", "doy_for_lib")) next
      # Preserve Veg column (do not interpolate label)
      if (col == "Veg") {
        new_row[[col]] <- a[[col]]
        next
      }
      val <- if (runif(1) < 0.5) a[[col]] else b[[col]]
      new_row[[col]] <- val
    }

    if ("doy_for_lib" %in% names(df_class)) {
      da <- as.integer(a$doy_for_lib)
      db <- as.integer(b$doy_for_lib)
      if (is.finite(da) && is.finite(db)) {
        diff_ab <- ((db - da + 365) %% 365)
        pick_doy <- ((da + round(alpha * diff_ab) - 1) %% 365) + 1
        new_row$doy_for_lib <- as.integer(pick_doy)
      } else {
        new_row$doy_for_lib <- if (!is.na(da) && is.finite(da)) da else if (!is.na(db) && is.finite(db)) db else NA_integer_
      }
    }

    if ("date" %in% names(df_class)) {
      if (!is.null(new_row$doy_for_lib) && is.finite(new_row$doy_for_lib)) {
        new_row$date <- as.Date(sprintf("1970-01-01")) + (as.integer(new_row$doy_for_lib) - 1)
      } else {
        new_row$date <- NA
      }
    }

    synth_rows[[i]] <- new_row[1, , drop = FALSE]
  }

  synth_df <- chunked_rbind(synth_rows, chunk_size = 50L)
  gc()
  out <- rbind(df_class, synth_df)
  for (col in other_cols) {
    if (is.factor(df_class[[col]])) out[[col]] <- factor(out[[col]], levels = levels(df_class[[col]]))
  }
  rownames(out) <- NULL
  out
}

balanced_dfs <- list()
for (v in vegs) {
  dveg <- lib_df[lib_df$Veg == v & is.finite(lib_df$doy_for_lib), , drop = FALSE]
  if (nrow(dveg) == 0) next
  if (nrow(dveg) < target_size) {
    augmented <- tryCatch(
      {
        augment_minority_class(dveg, target_size)
      },
      error = function(e) stop(sprintf("augment_minority_class failed for %s: %s", v, e$message))
    )
    balanced_dfs[[v]] <- augmented
    cat(sprintf("Augmented %s from %d to %d samples\n", v, nrow(dveg), nrow(augmented)))
  } else {
    balanced_dfs[[v]] <- dveg
    cat(sprintf("Kept %s with %d samples\n", v, nrow(dveg)))
  }
}

if (length(balanced_dfs) > 0) {
  lib_df <- chunked_rbind(balanced_dfs, chunk_size = 50L)
  gc()
  # DEBUG: Check Veg column
  if (!"Veg" %in% names(lib_df)) {
    cat("[ERROR] Veg column missing after balancing!\n")
    cat("Available columns:", paste(names(lib_df), collapse=", "), "\n")
  } else {
    cat(sprintf("[DEBUG] After balancing: %d rows with Veg column\n", nrow(lib_df)))
    cat(sprintf("[DEBUG] Veg values: %s\n", paste(unique(lib_df$Veg), collapse=",")))
  }
  cat(sprintf(
    "Training data balanced: %d total samples across %d vegetation types\n",
    nrow(lib_df), length(balanced_dfs)
  ))
}

# Build basic library
# REVISED: We skip the per-DOY averaging step here.
# The 'lib' object will be constructed later from the variant centroids.
# We initialize an empty list for now to satisfy downstream checks until we rebuild it.
lib <- list()
# Placeholder for legacy code compatibility
for (v in vegs) {
  lib[[v]] <- list(n_samples = 0)
}

cat("Skipping per-DOY averaging (will build variants from raw traces)...\n")

timing_info$lib_construction_done <- Sys.time()
cat(sprintf(
  "Library construction (skipped) completed in %.1f seconds\n",
  as.numeric(difftime(timing_info$lib_construction_done, timing_info$moving_var_done, units = "secs"))
))

# Build PCA projections
# REVISED: Build Global PCA directly from raw training data (lib_df)
# instead of using the averaged M_list.

cat("Building Global PCA from raw training data...\n")

# 1. Construct X_all from lib_df
# We need to stack all valid rows for the selected indices
# lib_df has columns: avail (indices) and potentially _var14 columns

# Identify columns to include in PCA
pca_cols <- avail
# Check if we should include variance columns
# The original code included _mv columns if they existed in lib.
# Here we check if _var14 columns exist in lib_df and are valid.
var_cols_present <- paste0(avail, "_var14")
var_cols_present <- intersect(var_cols_present, names(lib_df))

# We will build a matrix X_all where columns are [avail, var_cols_present]
# But we need to be careful about column ordering to match 'avail_aug' logic later.
# The original code used 'avail' then 'avail_mv'. Let's stick to that pattern if possible,
# or just use whatever columns we have.
# Let's define the feature set explicitly.

feature_cols <- c(avail, var_cols_present)
cat(sprintf("PCA features: %s\n", paste(feature_cols, collapse=", ")))

# Extract data matrix
X_all <- as.matrix(lib_df[, feature_cols, drop = FALSE])

# Handle NAs/Infs in X_all
# Replace non-finite with column means
for (j in seq_len(ncol(X_all))) {
  col_vals <- X_all[, j]
  if (any(!is.finite(col_vals))) {
    mu_j <- mean(col_vals[is.finite(col_vals)], na.rm = TRUE)
    if (!is.finite(mu_j)) mu_j <- 0
    X_all[!is.finite(col_vals), j] <- mu_j
  }
}

# 2. Compute PCA
mu_all <- colMeans(X_all)
Xc <- sweep(X_all, 2, mu_all, "-")
feature_sds <- apply(Xc, 2, sd)
feature_sds[feature_sds <= 1e-10] <- 1.0

Xs_std <- sweep(Xc, 2, feature_sds, "/")

sv <- try(svd(Xs_std, nu = 0, nv = min(ncol(Xs_std), nrow(Xs_std))), silent = TRUE)
if (inherits(sv, "try-error") || length(sv$d) == 0) {
  stop("PCA failed: SVD did not produce usable singular values")
}

eigenvalues <- sv$d^2 / (nrow(Xs_std) - 1)
keep_idx <- which(eigenvalues > 1)

if (length(keep_idx) == 0) {
  warning("Kaiser criterion: No components with eigenvalue > 1, using first component")
  pca_rank <- 1
} else {
  pca_rank <- min(length(keep_idx), MAX_FACTORS_CAP, ncol(sv$v))
}

V_pca <- sv$v[, seq_len(pca_rank), drop = FALSE]
# Ensure orthogonality (SVD V is already orthogonal, but qr.Q is safe)
if (ncol(V_pca) > 1) {
  V_pca <- qr.Q(qr(V_pca))[, seq_len(ncol(V_pca)), drop = FALSE]
}

variance_explained <- sv$d^2 / sum(sv$d^2)
cumulative_variance <- cumsum(variance_explained)

cat(sprintf(
  "PCA trained on raw data (%d rows); retained %d components (cum var: %.1f%%)\n",
  nrow(X_all), ncol(V_pca), 100 * cumulative_variance[pca_rank]
))

# 3. Construct GLOBAL_PCA object
# We need to map the column names of X_all to the expected 'idx_order'
# The original code expected 'idx' and 'idx_mv'.
# Our feature_cols are 'idx' and 'idx_var14'.
# We should rename 'idx_var14' to 'idx_mv' in the idx_order to match downstream expectations
# if downstream code expects '_mv'.
# Looking at 'project_row_to_factors', it constructs names as paste0(idx, "_mv").
# So we should ensure our idx_order uses "_mv" suffix for variance columns.

idx_order <- feature_cols
# Rename _var14 to _mv in idx_order
idx_order <- gsub("_var14$", "_mv", idx_order)

GLOBAL_PCA <- list(
  global = list(
    V = V_pca,
    col_means = mu_all,
    global_avg_scales = rep(1.0, length(mu_all)),
    feature_sds = feature_sds,
    rank = ncol(V_pca)
  ),
  V = V_pca,
  col_means = mu_all,
  global_avg_scales = rep(1.0, length(mu_all)),
  feature_sds = feature_sds,
  per_bin = NULL,
  rank = ncol(V_pca),
  idx_order = idx_order
)

# Update lib_factor_pca (Legacy support / Fallback)
# We can project the raw data to get T scores, but lib_factor_pca expects
# a single T matrix (365 x k) per vegetation type (representing the "average" or "medoid").
# Since we skipped building the averaged 'lib', we don't have 'M_list'.
# We can construct a simple 'lib_factor_pca' by taking the medoid of the projected raw traces.

cat("Computing legacy lib_factor_pca from raw projections...\n")
lib_factor_pca <- list()
gpca_order <- GLOBAL_PCA$idx_order
gpca_means <- GLOBAL_PCA$col_means
vpca <- GLOBAL_PCA$global$V
pca_rank <- GLOBAL_PCA$global$rank

for (vname in vegs) {
  # Get raw data for this veg
  dveg <- lib_df[lib_df$Veg == vname, , drop = FALSE]
  if (nrow(dveg) < 10) {
     lib_factor_pca[[vname]] <- list(T = matrix(0, nrow = 365, ncol = pca_rank), n_samples = 0)
     next
  }
  
  # Project all rows
  # We can do this efficiently using matrix multiplication since we have X_all logic
  # But dveg is a subset.
  
  # Extract features for dveg
  X_v <- as.matrix(dveg[, feature_cols, drop = FALSE])
  # Handle NAs
  for (j in seq_len(ncol(X_v))) {
    col_vals <- X_v[, j]
    if (any(!is.finite(col_vals))) {
      X_v[!is.finite(col_vals), j] <- mu_all[j]
    }
  }
  
  # Standardize and Project
  X_v_c <- sweep(X_v, 2, mu_all, "-")
  X_v_std <- sweep(X_v_c, 2, feature_sds, "/")
  T_scores <- X_v_std %*% vpca
  
  # Now we need to form a single 365-day trace (medoid)
  # Group by DOY
  doy_vec <- lubridate::yday(dveg$date)
  T_medoid <- matrix(NA_real_, nrow = 365, ncol = pca_rank)
  
  for (d in 1:365) {
    rows_d <- which(doy_vec == d)
    if (length(rows_d) > 0) {
      sub <- T_scores[rows_d, , drop = FALSE]
      # Find medoid
      if (nrow(sub) == 1) {
        T_medoid[d, ] <- sub[1, ]
      } else {
        # Geometric median or medoid
        # Simple medoid: row with min sum of squared distances
        center <- colMeans(sub)
        dists <- rowSums(sweep(sub, 2, center, "-")^2)
        T_medoid[d, ] <- sub[which.min(dists), ]
      }
    }
  }
  
  # Interpolate missing days
  for (j in 1:pca_rank) {
    y <- T_medoid[, j]
    if (any(is.na(y))) {
      x <- which(!is.na(y))
      if (length(x) >= 2) {
        T_medoid[is.na(y), j] <- approx(x, y[x], xout = which(is.na(y)), rule = 2)$y
      } else if (length(x) == 1) {
        T_medoid[, j] <- y[x[1]]
      } else {
        T_medoid[, j] <- 0
      }
    }
  }
  
  lib_factor_pca[[vname]] <- list(T = T_medoid, n_samples = nrow(dveg))
}

cat("Global projections computed.\n")

timing_info$pca_computation_done <- Sys.time()
        cat(sprintf(
          "PCA computation completed in %.1f seconds\n",
          as.numeric(difftime(timing_info$pca_computation_done,
            timing_info$lib_construction_done, units = "secs"))
        ))

factor_mode <- exists("GLOBAL_PCA") && !is.null(GLOBAL_PCA) && is.list(GLOBAL_PCA) &&
  !is.null(GLOBAL_PCA$rank) && GLOBAL_PCA$rank >= 1L

if (!factor_mode) {
  stop("Between-veg factor projection unavailable")
}

# ================== MESMA-SPECIFIC FUNCTIONS ==================

# Optional: seasonal phase alignment before temporal compression
safe_diff <- function(x) {
  if (length(x) <= 1) return(numeric(0))
  dx <- diff(x)
  dx[!is.finite(dx)] <- 0
  dx
}

detect_phase_markers <- function(dvi_proxy) {
  n <- length(dvi_proxy)
  if (n < 5) return(c(1L, round(n/4), round(n/2), round(3*n/4), n))
  dvi <- as.numeric(dvi_proxy)
  dvi[!is.finite(dvi)] <- 0
  d1 <- safe_diff(dvi)
  d2 <- safe_diff(d1)
  # Windows for green-up and senescence
  gup_win <- 1:max(2L, floor(n * 0.5))
  sen_win <- (max(1L, floor(n * 0.5))):n
  green_up <- if (length(gup_win) > 1 && length(d1) >= max(gup_win)) which.max(d1[gup_win]) else floor(n * 0.25)
  green_up <- max(2L, min(green_up, n - 3L))
  peak <- suppressWarnings(which.max(dvi)); if (!is.finite(peak) || length(peak) == 0) peak <- floor(n * 0.5)
  senescence <- if (length(sen_win) > 1 && length(d1) >= max(sen_win) && min(sen_win) <= length(d1)) {
    which.max(-d1[sen_win]) + (min(sen_win) - 1)
  } else floor(n * 0.75)
  senescence <- max(green_up + 1L, min(senescence, n - 1L))
  sort(unique(as.integer(c(1L, green_up, peak, senescence, n))))
}

align_to_phenological_phase <- function(pca_mat, reference_phase = REFERENCE_PHASE_MARKERS) {
  # pca_mat: 365 x k matrix of factor time series
  if (is.null(pca_mat) || !is.matrix(pca_mat) || nrow(pca_mat) < 5) return(pca_mat)
  k <- ncol(pca_mat)
  # DVI proxy from first up to 3 components
  use_cols <- seq_len(min(3L, k))
  dvi_proxy <- rowSums(abs(pca_mat[, use_cols, drop = FALSE]))
  dvi_proxy[!is.finite(dvi_proxy)] <- 0
  phase_markers <- detect_phase_markers(dvi_proxy)
  ref <- reference_phase
  if (is.null(ref) || length(ref) < 2) ref <- c(1, 90, 180, 270, 365)
  ref <- sort(unique(as.integer(pmax(1, pmin(nrow(pca_mat), ref)))))
  # Build piecewise-linear mapping from phase_markers -> ref
  aligned_mat <- matrix(NA_real_, nrow = nrow(pca_mat), ncol = k)
  for (j in seq_len(k)) {
    y <- pca_mat[, j]
    y[!is.finite(y)] <- 0
    # Sample values at phase markers
    pm <- phase_markers
    rm <- ref
    # Interpolate segment-wise by mapping original day indices to reference markers
    # Construct a monotone mapping function f: [1..n] -> [1..n]
    x_src <- pm
    x_dst <- rm[seq_len(min(length(rm), length(pm)))]
    # Ensure strictly increasing
    x_src <- unique(sort(x_src))
    x_dst <- unique(sort(x_dst))
    len_map <- min(length(x_src), length(x_dst))
    x_src <- x_src[seq_len(len_map)]
    x_dst <- x_dst[seq_len(len_map)]
    if (length(x_src) < 2) {
      aligned_mat[, j] <- y
    } else {
      # Map destination grid 1..n back to source via inverse approx
      fwd <- stats::approx(x = x_dst, y = x_src, xout = seq_len(nrow(pca_mat)), rule = 2)$y
      aligned_mat[, j] <- stats::approx(x = seq_len(nrow(pca_mat)), y = y, xout = fwd, rule = 2)$y
    }
  }
  aligned_mat[!is.finite(aligned_mat)] <- 0
  aligned_mat
}

# ===== Uncertainty: GLS + Block Bootstrap Helpers =====
estimate_block_size <- function(resid_series, min_block = GLSBB_MIN_BLOCK, max_block = GLSBB_MAX_BLOCK) {
  x <- as.numeric(resid_series)
  x <- x[is.finite(x)]
  if (length(x) < 10) return(max(min_block, 5L))
  acf_vals <- tryCatch(as.numeric(stats::acf(x, plot = FALSE, na.action = na.pass)$acf)[-1], error = function(e) stop(sprintf("estimate_block_size: acf computation failed: %s", e$message)))
  if (length(acf_vals) == 0) stop("estimate_block_size: acf returned no values; cannot estimate block size")
  thr <- 0.1
  k <- which(abs(acf_vals) < thr)[1]
  if (is.na(k)) k <- length(acf_vals)
  b <- max(min_block, min(max_block, 2L * k + 1L))
  as.integer(b)
}

solve_weights_gls <- function(E, y, lambda = 1e-6, sum_to_one = TRUE, non_negative = TRUE) {
  # Estimate HAC covariance of residuals under OLS fit, then do GLS via whitening
  n <- length(y); p <- ncol(E)
  if (!is.matrix(E) || nrow(E) != n || p < 1) return(solve_weights_simplex(E, y, ridge = lambda))
  # Initial weights via ridge simplex
  w0 <- tryCatch(solve_weights_simplex(E, y, ridge = lambda), error = function(e) stop(sprintf("solve_weights_gls: initial simplex solve failed: %s", e$message)))
  r0 <- as.numeric(y - E %*% w0)
  # Newey-West HAC covariance (Toeplitz) approximation from acf
  ac <- tryCatch(as.numeric(stats::acf(r0, plot = FALSE, na.action = na.pass)$acf), error = function(e) stop(sprintf("solve_weights_gls: acf computation failed: %s", e$message)))
  if (length(ac) < 2) stop("solve_weights_gls: insufficient autocorrelation values to estimate HAC covariance")
  ac <- ac[-1]
  L <- min(length(ac), max(1L, floor(sqrt(n))))
  gamma <- c(1, ac[seq_len(L)])
  # Build Toeplitz covariance Sigma with taper (Bartlett weights)
  w_bart <- c(1, 1 - (seq_len(L) / (L + 1)))
  gamma_t <- gamma * w_bart
  Sigma <- tryCatch(stats::toeplitz(gamma_t), error = function(e) stop(sprintf("solve_weights_gls: Toeplitz covariance construction failed: %s", e$message)))
  if (any(!is.finite(Sigma))) stop("solve_weights_gls: constructed covariance contains non-finite values")
  # Regularize Sigma
  Sigma <- Sigma + diag(EPS_SIGMA, nrow(Sigma))
  # Whitening
  cinv <- tryCatch(chol(Sigma), error = function(e) stop(sprintf("solve_weights_gls: Cholesky failed on Sigma: %s", e$message)))
  A <- backsolve(cinv, cbind(y, E), transpose = TRUE)
  y_w <- A[,1]; E_w <- A[,-1, drop = FALSE]
  # Constrained solve in whitened space
  if (isTRUE(ENABLE_QP_SOLVER) && requireNamespace("quadprog", quietly = TRUE)) {
    return(solve_weights_constrained(E_w, y_w, lambda = lambda, sum_to_one = sum_to_one, non_negative = non_negative))
  }
  solve_weights_simplex(E_w, y_w, ridge = lambda)
}

gls_block_bootstrap <- function(y_vec, comp_templates, top_variants, chosen_ids, w_hat, B = BOOTSTRAP_B, lambda_star = 1e-6, seed = 123) {
  # Returns list with coef_ci (data.frame), rmse_ci (numeric vector), variant_freq (data.frame)
  if (length(y_vec) < 5) return(NULL)
  set.seed(seed)
  # Rebuild E_best with columns per veg in names(top_variants) in fixed order
  vegs <- names(top_variants)
  cols <- list(); id_map <- list()
  for (v in vegs) {
    cid <- chosen_ids[[v]]
    vec <- NULL
    for (cand in comp_templates[[v]]) { if (!is.null(cand$id) && cand$id == cid) { vec <- cand$vec; break } }
    if (is.null(vec)) return(NULL)
    cols[[length(cols) + 1]] <- as.numeric(vec)
    id_map[[v]] <- sapply(comp_templates[[v]], function(c) c$id)
  }
  E_best <- do.call(cbind, cols)
  for (j in seq_len(ncol(E_best))) { nj <- sqrt(sum(E_best[, j]^2)); if (nj > 0) E_best[, j] <- E_best[, j] / nj }
  # Residuals from current fit to size block
  r <- as.numeric(y_vec - E_best %*% as.numeric(w_hat))
  bsz <- estimate_block_size(r)
  n <- length(y_vec)
  # Precompute blocks (circular to preserve ends)
  starts <- seq_len(n)
  res_mat <- matrix(NA_real_, nrow = n, ncol = B)
  coef_mat <- matrix(NA_real_, nrow = length(vegs), ncol = B, dimnames = list(vegs, NULL))
  rmse_vec <- rep(NA_real_, B)
  dom_counts <- lapply(vegs, function(.) list())
  names(dom_counts) <- vegs

  for (b in seq_len(B)) {
    # Circular block bootstrap of residuals
    idx <- integer(0)
    while (length(idx) < n) {
      s <- sample(starts, 1)
      block <- ((s - 1) + seq_len(bsz) - 1) %% n + 1
      idx <- c(idx, block)
    }
    idx <- idx[seq_len(n)]
    r_b <- r[idx]
    # Perturb y: y* = E_best w_hat + r_b
    yb <- as.numeric(E_best %*% as.numeric(w_hat)) + r_b
    # Variant reselection: pick top per-veg by cosine, then solve GLS
    chosen_b <- list(); cols_b <- list()
    for (v in vegs) {
      cand <- comp_templates[[v]]
      sims <- sapply(cand, function(x) cos_sim(yb, x$vec))
      i_best <- which.max(sims)
      chosen_b[[v]] <- cand[[i_best]]$id
      cols_b[[length(cols_b) + 1]] <- cand[[i_best]]$vec
    }
    E_b <- do.call(cbind, cols_b)
    for (j in seq_len(ncol(E_b))) { nj <- sqrt(sum(E_b[, j]^2)); if (nj > 0) E_b[, j] <- E_b[, j] / nj }
    w_b <- tryCatch(solve_weights_gls(E_b, yb, lambda = lambda_star, sum_to_one = TRUE, non_negative = TRUE), error = function(e) stop(sprintf("gls_block_bootstrap: GLS solve failed during bootstrap: %s", e$message)))
    coef_mat[, b] <- as.numeric(w_b)
    pred_b <- as.numeric(E_b %*% w_b)
    rmse_vec[b] <- sqrt(mean((yb - pred_b)^2))
    # Track variant dominance
    for (v in vegs) {
      key <- chosen_b[[v]]
      prev <- dom_counts[[v]][[key]]
      if (is.null(prev)) prev <- 0
      dom_counts[[v]][[key]] <- prev + 1
    }
  }
  # Compute CIs (2.5/97.5) for coefficients and RMSE
  ci_low <- ci_high <- rep(NA_real_, length(vegs)); names(ci_low) <- vegs; names(ci_high) <- vegs
  for (i in seq_along(vegs)) {
    x <- coef_mat[i, ]; x <- x[is.finite(x)]
    if (length(x) >= 10) { ci <- stats::quantile(x, c(0.025, 0.975), na.rm = TRUE); ci_low[i] <- ci[1]; ci_high[i] <- ci[2] }
  }
  rmse_ci <- if (sum(is.finite(rmse_vec)) >= 10) as.numeric(stats::quantile(rmse_vec, c(0.025, 0.975), na.rm = TRUE)) else c(NA_real_, NA_real_)
  coef_ci_df <- data.frame(Veg = vegs, coef_lo = ci_low, coef_hi = ci_high, stringsAsFactors = FALSE)
  # Variant frequencies
  var_rows <- chunked_rbind(lapply(names(dom_counts), function(v) {
    freqs <- dom_counts[[v]]
    if (length(freqs) == 0) return(NULL)
    data.frame(Veg = v, Variant = names(freqs), Freq = as.numeric(unlist(freqs)), stringsAsFactors = FALSE)
  }))
  if (!is.null(var_rows) && nrow(var_rows) > 0) {
    var_rows <- var_rows %>% dplyr::group_by(.data$Veg, .data$Variant) %>% dplyr::summarize(N = sum(.data$Freq), .groups = "drop")
    var_rows$Percent <- 100 * var_rows$N / sum(var_rows$N)
  }
  list(coef_ci = coef_ci_df, rmse_ci = rmse_ci, variant_freq = var_rows, block_size = bsz)
}

# Helper: project a single row's raw indices to factor space (PCA)
project_row_to_factors <- function(row, gproj, avail_idx, veg_type = NULL) {
  if (is.null(gproj)) return(NULL)
  idx_order <- gproj$idx_order
  col_means <- gproj$col_means
  
  # Use global average scaling instead of vegetation-specific
  # No per-band scaling in use
  
  # Get standardization scales
  feature_sds <- if (!is.null(gproj$feature_sds)) gproj$feature_sds else rep(1, length(idx_order))
  
  if (is.null(idx_order) || is.null(col_means)) return(NULL)
  
  raw_vec <- rep(NA_real_, length(idx_order))
  names(raw_vec) <- idx_order
  
  # Process indices
  for (idx in avail_idx) {
    kpos <- match(idx, idx_order)
    if (is.na(kpos)) next
    val <- row[[idx]]
    if (!is.finite(val)) val <- col_means[kpos]
    raw_vec[kpos] <- as.numeric(val)
  }
  
  # Process variance features
  for (idx in avail_idx) {
    mv_name <- paste0(idx, "_mv")
    kpos_mv <- match(mv_name, idx_order)
    if (is.na(kpos_mv)) next
    src_col <- paste0(idx, "_var14")
    mv_val <- if (src_col %in% names(row)) row[[src_col]] else NA_real_
    if (!is.finite(mv_val)) mv_val <- col_means[kpos_mv]
    raw_vec[kpos_mv] <- as.numeric(mv_val)
  }
  
  nas <- which(!is.finite(raw_vec))
  if (length(nas) > 0) raw_vec[nas] <- col_means[nas]
  
  # Apply standardization using global feature standard deviations only
  centered <- (raw_vec - col_means) / feature_sds
  centered[!is.finite(centered)] <- 0
  
  V <- gproj$global$V
  if (is.null(V) || ncol(V) == 0) return(NULL)
  as.numeric(centered %*% V)
}

# Perform time-dimension reduction on all traces BEFORE variant construction
# This function projects raw data to PCA space and reduces each trace to a fixed-grid feature vector
reduce_all_traces <- function(lib_df, veg_types, global_pca, avail_idx, fixed_grid_size = TEMPORAL_BUDGET, 
                               enable_phase_alignment = FALSE, reference_phase = NULL,
                               enable_multiscale = FALSE, multiscale_windows = NULL) {
  cat("Performing time-dimension reduction on all traces...\n")
  
  # Fixed grid for all traces (ensures comparability)
  fixed_grid <- unique(round(seq(1, 365, length.out = fixed_grid_size)))
  
  reduced_data <- list()
  
  # Resolve which column in lib_df holds vegetation labels
  veg_col_name <- NULL
  if ("Veg" %in% names(lib_df)) veg_col_name <- "Veg"
  else if ("veg" %in% names(lib_df)) veg_col_name <- "veg"
  else {
    # Try a case-insensitive search for a likely vegetation column
    lcnames <- tolower(names(lib_df))
    idx <- which(lcnames %in% c("veg", "vegetation", "vegetation_type", "type"))
    if (length(idx) > 0) veg_col_name <- names(lib_df)[idx[1]]
  }

  if (is.null(veg_col_name)) {
    cat("[WARN] reduce_all_traces: training data has no vegetation label column ('Veg' or similar). No per-veg reduction will be performed.\n")
    return(list())
  }

  for (veg in veg_types) {
    # subset by resolved vegetation column
    veg_data <- lib_df[tolower(as.character(lib_df[[veg_col_name]])) == tolower(as.character(veg)), , drop = FALSE]
    
    if (!"location_id" %in% names(veg_data) || !"year" %in% names(veg_data)) {
      cat(sprintf("  [%s] Missing location_id or year, skipping\n", veg))
      next
    }
    
    traces <- unique(veg_data[, c("location_id", "year")])
    cat(sprintf("  [%s] Reducing %d traces to %d time points...\n", veg, nrow(traces), length(fixed_grid)))
    
    feature_list <- list()
    Z_list <- list()
    trace_info <- list()
    
    for (i in seq_len(nrow(traces))) {
      loc <- traces$location_id[i]
      yr <- traces$year[i]
      
      dly_year <- veg_data[veg_data$location_id == loc & veg_data$year == yr, , drop = FALSE]
      if (nrow(dly_year) < 5) next
      
      # 1. Build a full 365 x K raw-index matrix for this trace
      # For each index in avail_idx we choose the medoid per DOY and interpolate missing days
      idxs <- avail_idx
      K <- length(idxs)
      raw_mat <- matrix(NA_real_, nrow = 365, ncol = K)
      colnames(raw_mat) <- idxs

      for (j in seq_along(idxs)) {
        idn <- idxs[j]
        if (!idn %in% names(dly_year)) next
        vals_by_doy <- tryCatch(
          {
            tapply(seq_along(dly_year[[idn]]), dly_year$doy, function(indices) {
              vals <- dly_year[[idn]][indices]
              vals <- vals[is.finite(vals)]
              if (length(vals) == 0) return(NA_real_)
              if (length(vals) == 1) return(as.numeric(vals[1]))
              m <- mean(vals)
              as.numeric(vals[which.min(abs(vals - m))])
            })
          }, error = function(e) NULL)

        if (!is.null(vals_by_doy) && length(vals_by_doy) > 0) {
          doy_values <- as.integer(names(vals_by_doy))
          valid_doy <- doy_values >= 1 & doy_values <= 365
          raw_mat[doy_values[valid_doy], j] <- vals_by_doy[valid_doy]
        }

        # Interpolate missing days per index column
        colv <- raw_mat[, j]
        finite_idx <- which(is.finite(colv))
        if (length(finite_idx) == 0) {
          raw_mat[, j] <- 0
        } else if (length(finite_idx) == 1) {
          raw_mat[, j] <- colv[finite_idx]
        } else {
          x <- c(finite_idx - 365, finite_idx, finite_idx + 365)
          y <- rep(colv[finite_idx], 3)
          interp <- tryCatch(stats::approx(x = x, y = y, xout = seq_len(365), rule = 2)$y, error = function(e) rep(colv[finite_idx[1]], 365))
          interp <- as.numeric(interp)
          interp[!is.finite(interp)] <- median(colv[finite_idx], na.rm = TRUE)
          raw_mat[, j] <- interp
        }
      }

      # 2. PCA projection (keep for downstream variant T_pca matrices)
      date_list <- tryCatch(prepare_factor_data(dly_year, global_pca, avail_idx, veg), error = function(e) NULL)
      if (length(date_list) == 0) next
      Z <- tryCatch(build_Z365(date_list, global_pca$rank), error = function(e) NULL)
      if (is.null(Z)) next
      if (isTRUE(enable_phase_alignment) && !is.null(reference_phase)) Z <- align_to_phenological_phase(Z, reference_phase = reference_phase)

      # 3. Time-dimension reduction performed ON RAW INDEX MATRIX (ensures features = fixed_grid × indices)
      info <- compute_information_content(raw_mat)
      feat <- if (isTRUE(enable_multiscale) && !is.null(multiscale_windows)) {
        extract_multiscale_features(raw_mat, fixed_grid, info, windows = multiscale_windows)
      } else {
        extract_grid_features(raw_mat, fixed_grid, info)
      }
      
      if (any(!is.finite(feat))) next
      
      # Ensure feature vector has consistent names: index_t1..index_tN
      if (is.numeric(feat) && length(feat) == length(idxs) * length(fixed_grid)) {
        # Create column names in the order: idx1_t1..idx1_tN, idx2_t1..idx2_tN, ...
        grid_count <- length(fixed_grid)
        nm <- unlist(lapply(idxs, function(x) paste0(x, "_t", seq_len(grid_count))))
        names(feat) <- nm
      }
      feature_list[[length(feature_list) + 1]] <- feat
      Z_list[[length(Z_list) + 1]] <- Z
      trace_info[[length(trace_info) + 1]] <- list(
        location_id = loc,
        year = yr,
        trace_index = i
      )
    }
    
    if (length(feature_list) > 0) {
      reduced_data[[veg]] <- list(
        features = do.call(rbind, feature_list),  # Matrix: n_traces x n_features
        Z_matrices = Z_list,                       # List of full 365 x k matrices
        trace_info = trace_info,                   # Metadata for each trace
        n_samples = nrow(veg_data)                 # Total samples for this veg
      )
      cat(sprintf("  [%s] Reduced %d traces successfully\n", veg, length(feature_list)))
      # Sanity checks: expected feature width = TEMPORAL_BUDGET * n_indices
      feat_mat <- reduced_data[[veg]]$features
      grid_count <- length(fixed_grid)
      expected_cols <- grid_count * length(idxs)
      if (!is.null(dim(feat_mat)) && ncol(feat_mat) != expected_cols) {
        warning(sprintf("reduce_all_traces: feature column count mismatch for '%s' (got %d, expected %d).", veg, ncol(feat_mat), expected_cols))
      }
      keys <- sapply(reduced_data[[veg]]$trace_info, function(x) paste0(as.character(x$location_id), "__", as.character(x$year)))
      if (any(duplicated(keys))) {
        warning(sprintf("reduce_all_traces: duplicate location-year traces found for '%s' (%d duplicates) — expected exactly one row per loc/year.", veg, sum(duplicated(keys))))
      }
    }
  }
  
  cat(sprintf("Time-dimension reduction complete: %d vegetation types processed\n", length(reduced_data)))
  return(reduced_data)
}

# Build MESMA library with endmember variants
# REVISED: Accepts pre-reduced features (PCA projection and time-reduction already applied)
build_mesma_variants <- function(reduced_data, lib_factor_pca, min_cluster_size = 10) {
  mesma_lib <- list()
  
  # Whitening is performed by the top-level whiten_matrix() function.

  # [NEW] Identify barren prototype from reduced data if available
  barren_proto_raw <- NULL
  if ("barren" %in% names(reduced_data)) {
    b_feat <- reduced_data[["barren"]]$features
    if (!is.null(b_feat) && nrow(b_feat) > 0) {
      # Use median of barren features as prototype
      barren_proto_raw <- apply(b_feat, 2, median, na.rm = TRUE)
      cat(sprintf("[Soil Correction] Identified barren prototype from %d traces\n", nrow(b_feat)))
    }
  }

  for (veg in names(reduced_data)) {
    veg_info <- reduced_data[[veg]]
    
    # Extract pre-reduced features and full Z matrices
    X_feat <- veg_info$features        # n_traces x n_features (time-reduced)
    Z_list <- veg_info$Z_matrices      # List of 365 x k matrices (full resolution)
    n_samples <- veg_info$n_samples
    
    if (nrow(X_feat) < min_cluster_size) {
      # Fallback to single variant
      mesma_lib[[veg]] <- list(
        list(
          T_pca = lib_factor_pca[[veg]]$T,
          variant_id = paste0(veg, "_single"),
          n_samples = n_samples
        )
      )
      cat(sprintf("  [%s] Insufficient traces (%d), using single variant\n", veg, nrow(X_feat)))
      next
    }
    
    cat(sprintf("  [%s] Building variants from %d reduced traces...\n", veg, nrow(X_feat)))
    
    # Subsample if too many traces
    if (nrow(X_feat) > MAX_PROJECTIONS_PER_VEG) {
      set.seed(123)
      idx <- sample(nrow(X_feat), MAX_PROJECTIONS_PER_VEG)
      X_feat <- X_feat[idx, , drop = FALSE]
      Z_list <- Z_list[idx]
      cat(sprintf("  [%s] Subsampled to %d traces\n", veg, nrow(X_feat)))
    }
    # Whitening on the reduced features (if Xw precomputed in the reduced_data, use that)
    if (!is.null(veg_info$Xw)) {
      X_w <- as.matrix(veg_info$Xw)
      whitened <- list(Xw = X_w, W = if (!is.null(veg_info$W)) veg_info$W else diag(ncol(X_w)), mu = if (!is.null(veg_info$mu)) veg_info$mu else rep(0, ncol(X_w)))
    } else {
      whitened <- whiten_matrix(as.matrix(X_feat))
      X_w <- whitened$Xw
    }
    
    # Soil subtraction is only performed here when data are NOT pre-whitened and we have a raw barren prototype
    if (is.null(veg_info$Xw) && !is.null(barren_proto_raw) && veg != "barren") {
      if (length(barren_proto_raw) == ncol(X_w)) {
        b_centered <- barren_proto_raw - whitened$mu
        b_proj <- as.numeric(b_centered %*% whitened$W)
        b_norm2 <- sum(b_proj^2)
        if (b_norm2 > 1e-9) {
          alphas <- (X_w %*% b_proj) / b_norm2
          # Orthogonalization: remove component parallel to soil
          X_w <- X_w - (alphas %*% t(b_proj))
          cat(sprintf("  [%s] Applied soil subtraction (orthogonalization) in whitened space\n", veg))
        }
      }
    }
    
    # Endmember Construction: Cluster the whitened, time-reduced features
    
    max_k <- min(5, floor(nrow(X_w) / 5))
    if (max_k < 2) max_k <- 2
    
    best_k <- 2
    best_sil <- -Inf
    
    # Silhouette selection
    for (k in 2:max_k) {
      km <- kmeans(X_w, centers = k, nstart = 10, iter.max = 100)
      
      # Silhouette on subsample if needed
      nproj <- nrow(X_w)
      if (nproj > SILHOUETTE_SAMPLE_SIZE) {
        samp_idx <- sample.int(nproj, size = min(as.integer(SILHOUETTE_SAMPLE_SIZE), nproj))
        dsub <- dist(X_w[samp_idx, , drop = FALSE])
        labels_sub <- km$cluster[samp_idx]
        sil_obj <- tryCatch(cluster::silhouette(labels_sub, dsub), error = function(e) NULL)
        sil <- if (!is.null(sil_obj) && is.matrix(sil_obj)) mean(sil_obj[, 3], na.rm = TRUE) else NA_real_
      } else {
        sil_obj <- tryCatch(cluster::silhouette(km$cluster, dist(X_w)), error = function(e) NULL)
        sil <- if (!is.null(sil_obj) && is.matrix(sil_obj)) mean(sil_obj[, 3], na.rm = TRUE) else NA_real_
      }
      
      if (is.finite(sil) && sil > best_sil) {
        best_sil <- sil
        best_k <- k
      }
    }
    
    km_final <- kmeans(X_w, centers = best_k, nstart = 25, iter.max = 100)
    
    # Create variants
    variants <- list()
    for (clust in seq_len(best_k)) {
      clust_members <- which(km_final$cluster == clust)
      if (length(clust_members) < 1) next
      
      # Find medoid in the WHITENED space
      # The medoid is the sample closest to the cluster center (or geometric median)
      # km_final$centers[clust, ] is the centroid in whitened space.
      # Find sample closest to centroid
      
      dists <- rowSums(sweep(X_w[clust_members, , drop=FALSE], 2, km_final$centers[clust, ], "-")^2)
      best_idx_local <- which.min(dists)
      best_idx_global <- clust_members[best_idx_local]
      
      # The variant is defined by the full Z matrix of the medoid trace
      # This allows it to be re-compressed dynamically during inference
      medoid_Z <- Z_list[[best_idx_global]]
      
      variants[[length(variants) + 1]] <- list(
        T_pca = medoid_Z,
        variant_id = paste0(veg, "_v", clust),
        n_samples = length(clust_members),
        cluster_center = km_final$centers[clust, ] # Store whitened center if needed
        , whitening_W = whitened$W
        , whitening_mu = whitened$mu
      )
    }
    
    mesma_lib[[veg]] <- variants
  }
  
  return(mesma_lib)
}

  

  solve_weights_simplex <- function(E, y, ridge = 1e-6) {
    EtE <- crossprod(E)
    Ety <- as.numeric(crossprod(E, matrix(y, ncol = 1)))
    p <- ncol(E)
    EtE_reg <- EtE + ridge * diag(p)
    w <- tryCatch({ as.numeric(solve(EtE_reg, Ety)) }, error = function(e) stop(sprintf("solve_weights_constrained: direct solve failed: %s", e$message)))
    if (any(!is.finite(w))) w[!is.finite(w)] <- 0
    project_to_simplex(w)
  }

  # Quadratic programming solver with constraints: sum(w)=1, w>=0
  solve_weights_constrained <- function(E, y, lambda = 1e-6, sum_to_one = TRUE, non_negative = TRUE) {
    p <- ncol(E)
    if (p < 1) return(numeric(0))
    if (!requireNamespace("quadprog", quietly = TRUE)) {
      return(solve_weights_simplex(E, y, ridge = lambda))
    }
    D <- crossprod(E)
    if (!is.matrix(D)) D <- as.matrix(D)
    # Ensure symmetry and positive definiteness
    D <- 0.5 * (D + t(D)) + (lambda + 1e-8) * diag(p)
    d <- as.numeric(crossprod(E, matrix(y, ncol = 1)))

    # Build constraints Amat (columns are constraints), bvec, and meq count
    Amat <- matrix(, nrow = p, ncol = 0)
    bvec <- numeric(0)
    meq <- 0L
    if (isTRUE(sum_to_one)) {
      Amat <- cbind(Amat, rep(1, p))
      bvec <- c(bvec, 1)
      meq <- meq + 1L
    }
    if (isTRUE(non_negative)) {
      Amat <- cbind(Amat, diag(p))
      bvec <- c(bvec, rep(0, p))
    }
    # Solve QP: min 1/2 w^T D w - d^T w s.t. t(Amat) %*% w >= bvec, with first meq equalities
    sol <- tryCatch(
      {
        quadprog::solve.QP(Dmat = D, dvec = d, Amat = Amat, bvec = bvec, meq = meq)
      },
      error = function(e) stop(sprintf("extract_grid_features / template accessor failed: %s", e$message))
    )
    if (is.null(sol) || is.null(sol$solution)) {
      return(solve_weights_simplex(E, y, ridge = lambda))
    }
    w <- as.numeric(sol$solution)
    # Numeric repair: clip negatives and renormalize if required
    if (isTRUE(non_negative)) w[w < 0] <- 0
    if (isTRUE(sum_to_one)) {
      s <- sum(w)
      if (s <= 0) {
        w <- rep(1 / p, p)
      } else {
        w <- w / s
      }
    }
    w
  }

  # Adaptive ridge selection using K-fold CV on pooled template columns
  estimate_optimal_ridge <- function(E, y, cv_folds = 5L, ridge_candidates = 10^seq(-8, -2, length.out = 20)) {
    n <- length(y)
    if (!is.matrix(E) || nrow(E) != n || n <= 3 || ncol(E) < 1) return(1e-6)
    cv_folds <- max(2L, min(cv_folds, n))
    idx <- seq_len(n)
    # Randomized, approximately equal folds
    idx <- sample(idx, length(idx))
    splits <- split(idx, cut(seq_along(idx), breaks = cv_folds, labels = FALSE))
    cv_err <- rep(Inf, length(ridge_candidates))
    for (ci in seq_along(ridge_candidates)) {
      lam <- ridge_candidates[ci]
      errs <- c()
      for (fold in seq_along(splits)) {
        test_idx <- splits[[fold]]
        train_idx <- setdiff(seq_len(n), test_idx)
        if (length(train_idx) < 2 || length(test_idx) < 1) next
        E_train <- E[train_idx, , drop = FALSE]
        y_train <- y[train_idx]
        E_test <- E[test_idx, , drop = FALSE]
        y_test <- y[test_idx]
        w <- tryCatch(solve_weights_simplex(E_train, y_train, ridge = lam), error = function(e) stop(sprintf("estimate_optimal_ridge: solve_weights_simplex failed: %s", e$message)))
        pred <- as.numeric(E_test %*% w)
        errs <- c(errs, mean((y_test - pred)^2, na.rm = TRUE))
      }
      cv_err[ci] <- if (length(errs) > 0 && all(is.finite(errs))) mean(errs) else Inf
    }
    best <- which.min(cv_err)
    lam_star <- if (length(best) == 1L && is.finite(cv_err[best])) ridge_candidates[best] else 1e-6
    lam_star
  }

  # Retrieve precomputed templates into the structure used by downstream code
  get_variant_templates <- function(veg_types, observation_grid_type, compressed_templates, mesma_lib) {
    templates <- list()
    for (veg in veg_types) {
      veg_templates <- list()
      for (variant in mesma_lib[[veg]]) {
        vec <- compressed_templates(veg, variant$variant_id, grid_type = observation_grid_type)
        if (!is.null(vec)) {
          veg_templates[[length(veg_templates) + 1]] <- list(vec = vec, id = variant$variant_id)
        }
      }
      templates[[veg]] <- veg_templates
    }
    templates
  }

  # Vectorized evaluation of combinations with chunking
  evaluate_all_combinations <- function(y, templates, lambda, early_stop_rmse = EARLY_STOP_RMSE_THRESHOLD) {
    veg_names <- names(templates)
    n_veg <- length(veg_names)
    if (n_veg == 0) return(NULL)
    # all indices per veg
    idx_lists <- lapply(templates, function(lst) seq_along(lst))
    if (any(lengths(idx_lists) == 0)) return(NULL)

    # Avoid materializing the full Cartesian product in memory (expand.grid)
    # which can easily explode (OOM) when many veg types * candidates per veg.
    radices <- as.integer(lengths(idx_lists))
    n_veg <- length(radices)
    n_combos <- as.numeric(Reduce(`*`, as.list(radices), init = 1L))
    chunk_size <- 100L

    # Safety thresholds
    safe_expand_limit <- COMBO_SAFE_EXPAND_LIMIT    # small enough to expand safely in memory
    abort_limit <- COMBO_ABORT_LIMIT                # too many combinations — abort rather than OOM

    use_expanded_grid <- !is.na(n_combos) && n_combos <= safe_expand_limit
    if (use_expanded_grid) combos <- do.call(expand.grid, idx_lists)
    if (!is.na(n_combos) && n_combos > abort_limit) {
      stop(sprintf("Combination space is too large (%d combinations). Reduce TOPK_VARIANTS, MAX_VEG_COMPONENTS, or limit vegetation types to avoid huge memory allocations.", as.integer(n_combos)))
    }
    best_rmse <- Inf
    best_result <- NULL
    for (start in seq(1, n_combos, by = chunk_size)) {
      end <- min(start + chunk_size - 1L, n_combos)
      if (use_expanded_grid) {
        # small search space — use the materialized combos directly
        chunk <- combos[start:end, , drop = FALSE]
      } else {
        # lazily generate the mixed-radix indices for this chunk to avoid full expand.grid
        nrows <- end - start + 1L
        chunk <- matrix(NA_integer_, nrow = nrows, ncol = n_veg)
        colnames(chunk) <- names(idx_lists)
        # 0-based start for easier mixed-radix math
        base_start <- as.numeric(start - 1L)
        for (r in seq_len(nrows)) {
          x <- base_start + (r - 1L)
          # compute indices in each digit position
          rem <- x
          for (j in seq_len(n_veg)) {
            base <- radices[j]
            # mixed-radix digit (1-based index into candidate list)
            digit <- (rem %% base) + 1L
            chunk[r, j] <- as.integer(digit)
            rem <- rem %/% base
          }
        }
        # convert to data.frame-like access used later
        chunk <- as.data.frame(chunk, stringsAsFactors = FALSE, check.names = FALSE)
      }
      # build E matrices for the chunk
      results <- vector("list", nrow(chunk))
      for (i in seq_len(nrow(chunk))) {
        cols <- list(); ids <- character(0)
        for (v in veg_names) {
          idx <- as.integer(chunk[i, v])
          cols[[length(cols) + 1]] <- templates[[v]][[idx]]$vec
          ids <- c(ids, templates[[v]][[idx]]$id)
        }
        E <- do.call(cbind, cols)
        for (j in seq_len(ncol(E))) { nj <- sqrt(sum(E[, j]^2)); if (nj > 0) E[, j] <- E[, j] / nj }
        w <- if (isTRUE(ENABLE_QP_SOLVER) && requireNamespace("quadprog", quietly = TRUE)) {
          solve_weights_constrained(E, y, lambda = lambda, sum_to_one = TRUE, non_negative = TRUE)
        } else {
          solve_weights_simplex(E, y, ridge = lambda)
        }
        pred <- as.numeric(E %*% w)
        rmse <- sqrt(mean((y - pred)^2))
        results[[i]] <- list(w = w, rmse = rmse, ids = ids)
      }
      # pick best in chunk
      rmses <- vapply(results, function(r) r$rmse, numeric(1))
      i_best <- which.min(rmses)
      if (is.finite(rmses[i_best]) && rmses[i_best] < best_rmse) {
        best_rmse <- rmses[i_best]
        best_result <- results[[i_best]]
        names(best_result$w) <- veg_names
        names(best_result$ids) <- veg_names
      }
      if (is.finite(early_stop_rmse) && early_stop_rmse > 0 && best_rmse <= early_stop_rmse) break
    }
    best_result
  }

# New helper: Compress a single trace to a feature vector
compress_trace <- function(dly_year, global_pca, avail_idx, budget = TEMPORAL_BUDGET) {
  # Resolve vegetation type (not used for scaling, but kept for API consistency)
  veg_ref <- if ("Veg" %in% names(dly_year)) unique(na.omit(as.character(dly_year$Veg)))[1] else NA_character_
  if (is.na(veg_ref) || !nzchar(as.character(veg_ref))) veg_ref <- NA_character_

  date_list <- tryCatch({
    prepare_factor_data(dly_year, global_pca, avail_idx, veg_ref)
  }, error = function(e) NULL)
  
  if (length(date_list) == 0) return(NULL)
  
  k <- ncol(global_pca$global$V)
  Z <- tryCatch({ build_Z365(date_list, k) }, error = function(e) NULL)
  if (is.null(Z)) return(NULL)
  
  if (isTRUE(ENABLE_PHASE_ALIGNMENT)) {
    Z <- align_to_phenological_phase(Z, reference_phase = REFERENCE_PHASE_MARKERS)
  }
  
  info <- compute_information_content(Z)
  obs_grid <- create_adaptive_grid(info, budget)
  
  y <- if (isTRUE(ENABLE_MULTISCALE)) {
    extract_multiscale_features(Z, obs_grid, info, windows = MULTISCALE_WINDOWS)
  } else {
    extract_grid_features(Z, obs_grid, info)
  }

  gl <- length(obs_grid)
  grid_type <- if (gl <= ceiling(budget * 0.85)) "sparse" else if (gl >= ceiling(budget * 1.15)) "dense" else "medium"
  
  list(y = y, grid_type = grid_type, info = info, obs_grid = obs_grid)
}

# New helper: Compress Stage 1 Library (Barren/Veg prototypes)
compress_stage1_lib <- function(stage1_lib, global_pca, avail_idx, budget = TEMPORAL_BUDGET) {
  if (is.null(stage1_lib)) return(NULL)
  
  compress_proto <- function(proto_list) {
    # proto_list is list of indices, each with $mu (365 vector)
    # Construct a mock dly dataframe
    doy <- 1:365
    df_mock <- data.frame(doy = doy, date = as.Date("1970-01-01") + (doy - 1))
    
    valid_indices <- 0
    for (idx in avail_idx) {
      if (!is.null(proto_list[[idx]])) {
        df_mock[[idx]] <- proto_list[[idx]]$mu
        valid_indices <- valid_indices + 1
      } else {
        df_mock[[idx]] <- NA_real_
      }
      # Add dummy variance columns if needed by prepare_factor_data
      df_mock[[paste0(idx, "_var14")]] <- 0 
    }
    
    if (valid_indices == 0) return(NULL)
    
    res <- compress_trace(df_mock, global_pca, avail_idx, budget)
    if (!is.null(res)) res$y else NULL
  }
  
  y_barren <- compress_proto(stage1_lib$barren)
  y_veg <- compress_proto(stage1_lib$vegetation)
  
  if (is.null(y_barren) || is.null(y_veg)) return(NULL)
  
  list(barren = y_barren, vegetation = y_veg)
}

# New helper: Stage 1 Unmixing on Compressed Data
unmix_stage1_compressed <- function(y, compressed_stage1_lib) {
  if (is.null(compressed_stage1_lib)) return(NA_real_)
  
  B <- compressed_stage1_lib$barren
  V <- compressed_stage1_lib$vegetation
  X <- y
  
  # Ensure lengths match
  len <- min(length(X), length(B), length(V))
  X <- X[1:len]; B <- B[1:len]; V <- V[1:len]
  
  # Solve min ||X - alpha*B - (1-alpha)*V||^2
  # min ||(X - V) - alpha*(B - V)||^2
  
  diff_BV <- B - V
  diff_XV <- X - V
  denom <- sum(diff_BV^2)
  
  if (denom < 1e-10) return(0.5)
  
  alpha <- sum(diff_XV * diff_BV) / denom
  alpha <- max(0, min(1, alpha))
  
  return(1 - alpha) # Return vegetated fraction
}

# New helper: Stage 2 Unmixing on Compressed Data
unmix_stage2_compressed <- function(y, grid_type, mesma_lib, topK = TOPK_VARIANTS) {
  veg_types <- names(mesma_lib)
  if (length(veg_types) == 0) return(NULL)
  
  # Access compressed templates
  if (!exists(".COMPRESSED_TEMPLATES_ACCESSOR", envir = globalenv())) return(NULL)
  compressed_templates <- get(".COMPRESSED_TEMPLATES_ACCESSOR", envir = globalenv())
  
  comp_templates <- tryCatch({
    get_variant_templates(veg_types, grid_type, compressed_templates, mesma_lib)
  }, error = function(e) NULL)
  
  if (is.null(comp_templates)) return(NULL)
  
  top_variants <- list()
  for (v in veg_types) {
    cand <- comp_templates[[v]]
    if (length(cand) == 0) next
    sims <- sapply(cand, function(x) cos_sim(y, x$vec))
    ord <- order(sims, decreasing = TRUE)
    keep <- ord[seq_len(min(topK, length(ord)))]
    top_variants[[v]] <- cand[keep]
  }
  
  pool_cols <- list()
  for (v in names(top_variants)) {
    for (i in seq_along(top_variants[[v]])) pool_cols[[length(pool_cols) + 1]] <- top_variants[[v]][[i]]$vec
  }
  
  ny <- sqrt(sum(y^2))
  y_norm <- if (ny > 0) y / ny else y
  lambda_star <- 1e-6
  
  if (length(pool_cols) > 0) {
    E_pool <- do.call(cbind, lapply(pool_cols, function(col) {
      col <- as.numeric(col)
      nn <- sqrt(sum(col^2))
      if (nn > 0) col / nn else col
    }))
    lambda_star <- estimate_optimal_ridge(E_pool, y_norm)
    if (!is.finite(lambda_star) || lambda_star <= 0) lambda_star <- 1e-6
  }

  best <- list(rmse = Inf, w = NULL, chosen = NULL)
  res <- tryCatch({
    evaluate_all_combinations(y_norm, top_variants, lambda = lambda_star, early_stop_rmse = EARLY_STOP_RMSE_THRESHOLD)
  }, error = function(e) NULL)
  
  if (!is.null(res)) best <- list(rmse = res$rmse, w = res$w, chosen = res$ids)
  
  diagnostics <- NULL
  uncertainty <- NULL
  
  if (isTRUE(ENABLE_DIAGNOSTICS) && !is.null(best$w) && !is.null(best$chosen)) {
    E_best_cols <- list()
    for (v in names(best$chosen)) {
      for (variant in top_variants[[v]]) {
        if (variant$id == best$chosen[[v]]) {
          E_best_cols[[v]] <- variant$vec
          break
        }
      }
    }
    E_best <- do.call(cbind, E_best_cols)
    for (j in seq_len(ncol(E_best))) {
      nj <- sqrt(sum(E_best[, j]^2))
      if (nj > 0) E_best[, j] <- E_best[, j] / nj
    }
    diagnostics <- tryCatch({
      compute_diagnostics(y_norm, E_best, best$w, mesma_result = NULL)
    }, error = function(e) NULL)
  }

  if (isTRUE(ENABLE_UNCERTAINTY) && !is.null(best$w) && !is.null(best$chosen)) {
    uncertainty <- tryCatch({
      gls_block_bootstrap(
        y_vec = y_norm,
        comp_templates = comp_templates,
        top_variants = top_variants,
        chosen_ids = best$chosen,
        w_hat = best$w,
        B = BOOTSTRAP_B,
        lambda_star = lambda_star,
        seed = 123
      )
    }, error = function(e) NULL)
  }
  
  list(vegetation_proportions = best$w, chosen_variants = best$chosen, rmse = best$rmse, diagnostics = diagnostics, uncertainty = uncertainty)
}

compress_and_unmix_year <- function(dly_year, mesma_lib, budget = TEMPORAL_BUDGET, topK = TOPK_VARIANTS) {
  # Legacy wrapper if needed, but we will use the split functions in fit_one_task
  res <- compress_trace(dly_year, GLOBAL_PCA, avail, budget)
  if (is.null(res)) return(NULL)
  unmix_stage2_compressed(res$y, res$grid_type, mesma_lib, topK)
}

  # End of mesma_dynamic_programming inner helpers
  # Helper to prepare factor data
  prepare_factor_data <- function(dly, gpca, avail_idx, veg_type) {
    date_list <- list()
    dly$doy <- lubridate::yday(dly$date)
    dts <- sort(unique(dly$date))

    for (i_dt in seq_along(dts)) {
      dt <- dts[i_dt]
      sub <- dly[dly$date == dt, , drop = FALSE]
      sub <- sub[is.finite(sub$doy), , drop = FALSE]
      if (nrow(sub) == 0) next

      idx_present_orig <- intersect(avail_idx, names(sub))
      min_required <- ceiling(length(avail_idx) * MIN_IDX_PRESENCE)
      if (length(idx_present_orig) < min_required) {
        next
      }

      vals_aug <- c()
      vrows_list <- list()

      for (idx in idx_present_orig) {
        vals_idx <- sub[[idx]]
        vals_idx <- vals_idx[is.finite(vals_idx)]
        if (length(vals_idx) == 0) {
          vals_aug <- NULL
          break
        }

        # Medoid selection for duplicate observations (closest to mean)
        yy <- if (length(vals_idx) == 1) {
          as.numeric(vals_idx[1])
        } else {
          m <- mean(vals_idx, na.rm = TRUE)
          as.numeric(vals_idx[which.min(abs(vals_idx - m))])
        }
        if (!is.finite(yy)) next

        kpos <- match(idx, gpca$idx_order)
        if (is.na(kpos)) next

        mi <- gpca$col_means[kpos]
        # Do not apply per-band scaling — assume identity scaling (1.0)
        valid_veg <- !is.null(veg_type) && !is.na(veg_type) && nzchar(as.character(veg_type))
        veg_type_chr <- if (valid_veg) as.character(veg_type) else NA_character_
        veg_type_lc <- if (valid_veg) tolower(veg_type_chr) else NA_character_

        si <- 1.0
        vals_aug <- c(vals_aug, as.numeric((yy - mi) / si))
        vrows_list[[length(vrows_list) + 1]] <- gpca$V[kpos, , drop = FALSE]
      }
      if (is.null(vals_aug)) next

      for (idx in idx_present_orig) {
        mv_col <- paste0(idx, "_var14")
        if (!mv_col %in% names(sub)) next

        vals_mv <- sub[[mv_col]]
        vals_mv <- vals_mv[is.finite(vals_mv)]
        if (length(vals_mv) == 0) {
          vals_aug <- NULL
          break
        }

        # Medoid selection for duplicate moving-variance observations
        yy_mv <- if (length(vals_mv) == 1) {
          as.numeric(vals_mv[1])
        } else {
          m_mv <- mean(vals_mv, na.rm = TRUE)
          as.numeric(vals_mv[which.min(abs(vals_mv - m_mv))])
        }
        if (!is.finite(yy_mv)) next

        mv_name <- paste0(idx, "_mv")
        kpos_mv <- match(mv_name, gpca$idx_order)
        if (is.na(kpos_mv)) next

        mi_mv <- gpca$col_means[kpos_mv]
        # No per-band scaling: assume scale = 1
        si_mv <- 1.0
        vals_aug <- c(vals_aug, as.numeric((yy_mv - mi_mv) / si_mv))
        vrows_list[[length(vrows_list) + 1]] <- gpca$V[kpos_mv, , drop = FALSE]
      }
      if (is.null(vals_aug)) next

      if (length(vrows_list) == 0) next
      vsub <- do.call(rbind, vrows_list)
      if (ncol(vsub) < 1) next

      xtx <- crossprod(vsub)
      xty <- as.numeric(crossprod(vsub, matrix(vals_aug, ncol = 1)))

      ridge_small <- 0.01
      xtx_reg <- xtx + ridge_small * diag(ncol(xtx))
      z <- tryCatch(
        {
          solve(xtx_reg, xty)
        },
        error = function(e) {
          qr.solve(xtx_reg, xty)
        }
      )

      if (any(!is.finite(z))) {
        if (!requireNamespace("MASS", quietly = TRUE)) {
          stop("Matrix inversion failed and MASS package not available")
        }
        z <- as.numeric(MASS::ginv(xtx_reg) %*% xty)
      }

      if (any(!is.finite(z))) {
        stop("All matrix inversion methods failed")
      }

      date_list[[length(date_list) + 1]] <- list(doy = as.integer(sub$doy[1]), z = as.numeric(z))
    }
    date_list
  }

  # Missing utility functions for MESMA
  build_Z365 <- function(date_list, k) {
    # Build a 365 x k matrix from the date_list factor projections
    Z <- matrix(NA_real_, nrow = 365, ncol = k)
    for (entry in date_list) {
      if (!is.null(entry$doy) && !is.null(entry$z) && 
          entry$doy >= 1 && entry$doy <= 365 && 
          length(entry$z) == k) {
        Z[entry$doy, ] <- entry$z
      }
    }
    # Interpolate missing DOYs per component
    for (col_idx in seq_len(k)) {
      finite_mask <- is.finite(Z[, col_idx])
      if (sum(finite_mask) >= 2) {
        Z[!finite_mask, col_idx] <- approx(
          x = which(finite_mask),
          y = Z[finite_mask, col_idx],
          xout = which(!finite_mask),
          rule = 2
        )$y
      } else if (sum(finite_mask) == 1) {
        Z[, col_idx] <- Z[finite_mask, col_idx][1]
      } else {
        Z[, col_idx] <- 0
      }
    }
    Z[!is.finite(Z)] <- 0
    Z
  }

  precompute_whitening <- function(mesma_lib, compressed_templates, combinations) {
    # Placeholder for whitening precomputation - returns a simple accessor
    function(combo_key) {
      # Simple identity whitening for now
      list(whiten = function(x) x, unwhiten = function(x) x)
    }
  }

  compute_diagnostics <- function(y, E, w, mesma_result = NULL) {
    # Compute basic diagnostic statistics
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
    
    # Condition number of E
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

  # Build MESMA library
  # NOTE (2025-11-23): Fix to veg lookup when preparing factor data.
  # If a per-row `Veg` value is missing, we attempt to resolve the vegetation
  # type for a location from `gpts_map` (location -> Veg mapping created from
  # the source geojson), and finally fall back to the first available
  # vegetation in the trained library to reduce repeated "veg_ref=NA" cases.
  # calls which used to trigger 'subscript out of bounds' inside
  # `prepare_factor_data`.
  cat("Building MESMA endmember library...\n")
  
  # STEP 1: Time-dimension reduction on all traces (PCA projection + fixed-grid reduction)
  # Debug: print available columns to verify Veg column preserved
  cat(sprintf("[DEBUG] Columns in lib_df before reduction: %s\n", paste(names(lib_df), collapse = ", ")))
  # This happens BEFORE variant construction to ensure all traces are comparable
  reduced_traces <- reduce_all_traces(
    lib_df = lib_df,
    veg_types = names(lib),
    global_pca = GLOBAL_PCA,
    avail_idx = avail,
    fixed_grid_size = TEMPORAL_BUDGET,
    enable_phase_alignment = ENABLE_PHASE_ALIGNMENT,
    reference_phase = REFERENCE_PHASE_MARKERS,
    enable_multiscale = ENABLE_MULTISCALE,
    multiscale_windows = MULTISCALE_WINDOWS
  )
  
  # [NEW] Print rows per class after time compression
  cat("\n=== Rows per class after time compression ===\n")
  for (v in names(reduced_traces)) {
    cat(sprintf("  %s: %d rows\n", v, nrow(reduced_traces[[v]]$features)))
  }
  cat("=============================================\n\n")

  # STEP: Whitening (per-vegetation) and Soil Correction in whitened space
  cat("\n=== Performing per-vegetation whitening and soil correction ===\n")
  whitened_traces <- list()
  for (veg in names(reduced_traces)) {
    veg_info <- reduced_traces[[veg]]
    X_feat <- veg_info$features
    if (is.null(X_feat) || nrow(X_feat) == 0) next

    if (nrow(X_feat) < 10) {
      cat(sprintf("  [%s] Skipping whitening (insufficient data)\n", veg))
      whitened_traces[[veg]] <- list(
        Xw = X_feat,
        W = diag(ncol(X_feat)),
        mu = rep(0, ncol(X_feat)),
        Z_matrices = veg_info$Z_matrices,
        trace_info = veg_info$trace_info,
        n_samples = veg_info$n_samples
      )
      next
    }

    whitened <- whiten_matrix(as.matrix(X_feat))
    whitened_traces[[veg]] <- list(
      Xw = whitened$Xw,
      W = whitened$W,
      mu = whitened$mu,
      Z_matrices = veg_info$Z_matrices,
      trace_info = veg_info$trace_info,
      n_samples = veg_info$n_samples
    )
    cat(sprintf("  [%s] Whitened %d traces (features: %d -> %d dims)\n", veg, nrow(whitened$Xw), ncol(X_feat), ncol(whitened$Xw)))
  }

  # Identify barren prototype in whitened space for soil orthogonalization
  barren_proto_whitened <- NULL
  if ("barren" %in% names(whitened_traces)) {
    b_xw <- whitened_traces[["barren"]]$Xw
    if (!is.null(b_xw) && nrow(b_xw) > 0) {
      barren_proto_whitened <- apply(b_xw, 2, median, na.rm = TRUE)
      cat(sprintf("[Soil Correction] Identified barren prototype from %d whitened traces\n", nrow(b_xw)))
    }
  }

  # Apply orthogonalization (soil subtraction) in whitened space
  soil_corrected_traces <- list()
  for (veg in names(whitened_traces)) {
    veg_data <- whitened_traces[[veg]]
    X_w <- veg_data$Xw
    if (veg == "barren") {
      soil_corrected_traces[[veg]] <- veg_data
      cat(sprintf("  [%s] Kept without soil correction\n", veg))
      next
    }

    if (is.null(barren_proto_whitened)) {
      cat(sprintf("  [%s] No barren prototype available, skipping soil correction\n", veg))
      soil_corrected_traces[[veg]] <- veg_data
      next
    }

    if (length(barren_proto_whitened) == ncol(X_w)) {
      b_proj <- barren_proto_whitened
      b_norm2 <- sum(b_proj^2)
      if (b_norm2 > 1e-9) {
        alphas <- (X_w %*% b_proj) / b_norm2
        X_w_corrected <- X_w - (alphas %*% t(b_proj))
        veg_data$Xw <- X_w_corrected
        soil_corrected_traces[[veg]] <- veg_data
        cat(sprintf("  [%s] Applied soil orthogonalization (mean alpha: %.3f)\n", veg, mean(alphas)))
      } else {
        soil_corrected_traces[[veg]] <- veg_data
        cat(sprintf("  [%s] Barren norm too small, skipped correction\n", veg))
      }
    } else {
      soil_corrected_traces[[veg]] <- veg_data
      cat(sprintf("  [%s] Dimension mismatch, skipped correction\n", veg))
    }
  }

  # Build reduced data structure for building variants
  reduced_data_corrected <- list()
  for (veg in names(soil_corrected_traces)) {
    veg_data <- soil_corrected_traces[[veg]]
    reduced_data_corrected[[veg]] <- list(
      features = veg_data$Xw,  # whitened & soil-corrected features
      Z_matrices = veg_data$Z_matrices,
      trace_info = veg_data$trace_info,
      n_samples = veg_data$n_samples
    )
  }

  cat("Whitening and soil correction complete\n")

  # Safety fallback: if no corrected data produced, fall back to raw reduced_traces
  if (length(reduced_data_corrected) == 0) {
    cat("[NOTICE] No soil-corrected reduced data generated; falling back to raw reduced traces for variant construction\n")
    reduced_data_corrected <- reduced_traces
  }
  if (length(reduced_data_corrected) == 0) {
    stop("No reduced traces available for variant construction — cannot build MESMA library")
  }

  # [NEW] Compress Stage 1 Library for Inference
  cat("Compressing Stage 1 Library for inference...\n")
  COMPRESSED_STAGE1_LIB <- NULL
  if (exists("STAGE1_LIB") && !is.null(STAGE1_LIB)) {
    COMPRESSED_STAGE1_LIB <- compress_stage1_lib(STAGE1_LIB, GLOBAL_PCA, avail, TEMPORAL_BUDGET)
  }
  
  # STEP 2: Build variants from the soil-corrected whitened reduced traces (clustering in whitened space)
  mesma_lib <- build_mesma_variants(reduced_data_corrected, lib_factor_pca, min_cluster_size = 10)

  # Precompute compressed templates for all veg/variant/projection once
  .COMPRESSED_TEMPLATES_ACCESSOR <- precompute_compressed_templates(mesma_lib, TEMPORAL_BUDGET)
  assign(".COMPRESSED_TEMPLATES_ACCESSOR", .COMPRESSED_TEMPLATES_ACCESSOR, envir = globalenv())
  # Optionally precompute whitening for common combinations (3-way if available)
  if (length(mesma_lib) >= 3) {
    combos <- utils::combn(names(mesma_lib), min(3, length(mesma_lib)), simplify = FALSE)
    .WHITENING_DB_ACCESSOR <- precompute_whitening(mesma_lib, .COMPRESSED_TEMPLATES_ACCESSOR, combos)
    assign(".WHITENING_DB_ACCESSOR", .WHITENING_DB_ACCESSOR, envir = globalenv())
  }

  # Report library statistics
  for (veg in names(mesma_lib)) {
    n_variants <- length(mesma_lib[[veg]])
    cat(sprintf("%s: %d variants\n", veg, n_variants))
    for (var in mesma_lib[[veg]]) {
      cat(sprintf("  %s: %d samples\n", var$variant_id, var$n_samples))
    }
  }

  # Save a MESMA model cache so the trained model (library + projections + templates)
  # can be re-used later for testing/inference without re-training.
  save_mesma_cache <- function(cache_dir = file.path(OUT_DIR, "mesma_cache")) {
    cat("\n=== SAVING MESMA MODEL CACHE ===\n")
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    # 1. Core library components
    core <- list(
      lib = lib,
      mesma_lib = mesma_lib,
      lib_factor_pca = lib_factor_pca,
      veg_counts = veg_counts,
      avail = avail,
      ALLOWED_VEG = ALLOWED_VEG,
      BAND_SCALE = NULL,
      COMPRESSED_STAGE1_LIB = if (exists("COMPRESSED_STAGE1_LIB")) COMPRESSED_STAGE1_LIB else NULL
    )
    saveRDS(core, file = file.path(cache_dir, "mesma_library.rds"))

    # 2. Projection matrices and PCA state
    proj <- list(
      GLOBAL_PCA = if (exists("GLOBAL_PCA")) GLOBAL_PCA else NULL,
      pca_rank = if (exists("pca_rank")) pca_rank else NULL
    )
    saveRDS(proj, file = file.path(cache_dir, "projection_matrices.rds"))

    # 3. Compressed templates (precomputed feature vectors) accessed via accessor
    if (exists(".COMPRESSED_TEMPLATES_ACCESSOR", envir = globalenv())) {
      template_data <- list()
      for (veg in names(mesma_lib)) {
        template_data[[veg]] <- list()
        for (variant in mesma_lib[[veg]]) {
          for (grid_type in c("sparse", "medium", "dense")) {
            key <- paste(veg, variant$variant_id, grid_type, sep = "|")
            vec <- tryCatch(.COMPRESSED_TEMPLATES_ACCESSOR(veg, variant$variant_id, grid_type), error = function(e) stop(sprintf("Compressed templates accessor failed for %s|%s|%s: %s", veg, variant$variant_id, grid_type, e$message)))
            if (!is.null(vec)) template_data[[veg]][[key]] <- vec
          }
        }
      }
      saveRDS(template_data, file = file.path(cache_dir, "compressed_templates.rds"))
    }

    # 4. Configuration parameters
    cfg <- list(
      TEMPORAL_BUDGET = if (exists("TEMPORAL_BUDGET")) TEMPORAL_BUDGET else NULL,
      TOPK_VARIANTS = if (exists("TOPK_VARIANTS")) TOPK_VARIANTS else NULL,
      ENABLE_PHASE_ALIGNMENT = if (exists("ENABLE_PHASE_ALIGNMENT")) ENABLE_PHASE_ALIGNMENT else NULL,
      REFERENCE_PHASE_MARKERS = if (exists("REFERENCE_PHASE_MARKERS")) REFERENCE_PHASE_MARKERS else NULL,
      ENABLE_MULTISCALE = if (exists("ENABLE_MULTISCALE")) ENABLE_MULTISCALE else NULL,
      MULTISCALE_WINDOWS = if (exists("MULTISCALE_WINDOWS")) MULTISCALE_WINDOWS else NULL,
      ENABLE_QP_SOLVER = if (exists("ENABLE_QP_SOLVER")) ENABLE_QP_SOLVER else NULL,
      ENABLE_DIAGNOSTICS = if (exists("ENABLE_DIAGNOSTICS")) ENABLE_DIAGNOSTICS else NULL,
      ENABLE_UNCERTAINTY = if (exists("ENABLE_UNCERTAINTY")) ENABLE_UNCERTAINTY else NULL,
      BOOTSTRAP_B = if (exists("BOOTSTRAP_B")) BOOTSTRAP_B else NULL,
      GLSBB_MIN_BLOCK = if (exists("GLSBB_MIN_BLOCK")) GLSBB_MIN_BLOCK else NULL,
      GLSBB_MAX_BLOCK = if (exists("GLSBB_MAX_BLOCK")) GLSBB_MAX_BLOCK else NULL,
      MAX_VEG_COMPONENTS = if (exists("MAX_VEG_COMPONENTS")) MAX_VEG_COMPONENTS else NULL,
      MIN_IDX_PRESENCE = if (exists("MIN_IDX_PRESENCE")) MIN_IDX_PRESENCE else NULL,
      EPS_SIGMA = if (exists("EPS_SIGMA")) EPS_SIGMA else NULL,
      LOWER_BND = if (exists("LOWER_BND")) LOWER_BND else NULL,
      USE_INDICES_MIN = if (exists("USE_INDICES_MIN")) USE_INDICES_MIN else NULL,
      MIN_INDEX_SD = if (exists("MIN_INDEX_SD")) MIN_INDEX_SD else NULL,
      MAX_PCA_COMPONENTS = if (exists("MAX_PCA_COMPONENTS")) MAX_PCA_COMPONENTS else NULL,
      MAX_FACTORS_CAP = if (exists("MAX_FACTORS_CAP")) MAX_FACTORS_CAP else NULL
    )
    saveRDS(cfg, file = file.path(cache_dir, "config_params.rds"))

    # 5. Training metadata
    meta <- list(
      training_years = if (exists("TRAIN_YEARS")) TRAIN_YEARS else NULL,
      training_date_range = if (exists("df_train")) range(df_train$date, na.rm = TRUE) else NULL,
      n_training_samples = if (exists("df_train")) nrow(df_train) else NULL,
      n_locations_trained = if (exists("df_train")) length(unique(df_train$location_id)) else NULL,
      indices_used = avail,
      vegetation_types = names(mesma_lib),
      model_creation_time = Sys.time()
    )
    saveRDS(meta, file = file.path(cache_dir, "training_metadata.rds"))

    # 6. Manifest with checksums
    manifest <- list(
      version = "1.0",
      created = Sys.time(),
      files = c("mesma_library.rds", "projection_matrices.rds", "compressed_templates.rds", "config_params.rds", "training_metadata.rds"),
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

  # Persist cache now so the trained model can be used later without re-training
  tryCatch({
    cache_dir <- save_mesma_cache()
  }, error = function(e) {
    cat(sprintf("Failed to write MESMA cache: %s\n", e$message))
  })

  # Prepare inference/test tasks — always use the full dataset (all years)
  cat("Constructing task list from inference dataset (all years)...\n")
  df_tasks <- df_full

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
      # warn if formats look mismatched
      s1 <- unique(na.omit(df_tasks$location_id))
      s2 <- unique(na.omit(gpts_map$location_id))
      if (length(s1) && length(s2) && all(grepl("^[0-9]+$", s1)) && any(grepl("^L_", s2))) {
        cat("[WARNING] df_tasks$location_id looks numeric while gpts_map$location_id looks like 'L_x_y' strings — matching will likely fail.\n")
      }
      # Use suffix to avoid overwriting existing columns and prefer existing task values
      df_tasks <- dplyr::left_join(df_tasks, gpts_map, by = "location_id", suffix = c("", ".geo"))
      if ("Veg.geo" %in% names(df_tasks)) {
        df_tasks$Veg <- ifelse(is.na(df_tasks$Veg) | df_tasks$Veg == "", df_tasks$Veg.geo, df_tasks$Veg)
        df_tasks$Veg.geo <- NULL
      }
      if ("no soil.geo" %in% names(df_tasks)) {
        if ("no soil" %in% names(df_tasks)) {
          df_tasks$`no soil` <- ifelse(is.na(df_tasks$`no soil`), df_tasks$`no soil.geo`, df_tasks$`no soil`)
        } else {
          df_tasks$`no soil` <- df_tasks$`no soil.geo`
        }
        df_tasks$`no soil.geo` <- NULL
      }
      # Normalize possible alternate column names (e.g. no.soil, no_soil)
      df_tasks <- normalize_no_soil_col(df_tasks)
    } else {
      # perform a straight join attempt (this will error if keys missing)
      df_tasks <- dplyr::left_join(df_tasks, gpts_map, by = "location_id", suffix = c("", ".geo"))
      if ("Veg.geo" %in% names(df_tasks)) {
        df_tasks$Veg <- ifelse(is.na(df_tasks$Veg) | df_tasks$Veg == "", df_tasks$Veg.geo, df_tasks$Veg)
        df_tasks$Veg.geo <- NULL
      }
      if ("no soil.geo" %in% names(df_tasks)) {
        if ("no soil" %in% names(df_tasks)) {
          df_tasks$`no soil` <- ifelse(is.na(df_tasks$`no soil`), df_tasks$`no soil.geo`, df_tasks$`no soil`)
        } else {
          df_tasks$`no soil` <- df_tasks$`no soil.geo`
        }
        df_tasks$`no soil.geo` <- NULL
      }
      df_tasks <- normalize_no_soil_col(df_tasks)
    }
  }

  if ("Veg" %in% names(df_tasks)) df_tasks$Veg <- tolower(as.character(df_tasks$Veg))
  if ("date" %in% names(df_tasks)) {
    df_tasks$date <- as.Date(df_tasks$date)
    df_tasks$year <- as.integer(lubridate::year(df_tasks$date))
    df_tasks$doy <- lubridate::yday(df_tasks$date)
    df_tasks$doy[df_tasks$doy < 1 | df_tasks$doy > 366] <- NA_integer_
  }

  test_loc_years <- df_tasks %>%
    dplyr::filter(!is.na(.data$location_id) & .data$location_id != "" & !is.na(.data$year) & .data$year > 0) %>%
    dplyr::distinct(.data$location_id, .data$year)

  cat(sprintf(
    "Final test_loc_years: %d rows from %d unique locations\n",
    nrow(test_loc_years), length(unique(test_loc_years$location_id))
  ))

  n_loc_years <- nrow(test_loc_years)
  task_list <- lapply(seq_len(n_loc_years), function(i) {
    list(loc = test_loc_years$location_id[i], yr = test_loc_years$year[i])
  })

  # Before task processing, validate one location
cat("\n=== PRE-FLIGHT CHECK ===\n")
test_task <- task_list[[1]]
test_dly <- df_tasks[df_tasks$location_id == test_task$loc, , drop = FALSE]
test_dly_year <- test_dly[lubridate::year(test_dly$date) == test_task$yr, , drop = FALSE]

cat(sprintf("Test location: %s, year: %d\n", test_task$loc, test_task$yr))
cat(sprintf("Rows in test_dly_year: %d\n", nrow(test_dly_year)))
cat(sprintf("Columns available: %s\n", paste(names(test_dly_year), collapse=", ")))
cat(sprintf("Var14 columns: %s\n", paste(grep("_var14$", names(test_dly_year), value=TRUE), collapse=", ")))
cat("======================\n\n")

  if (length(task_list) == 0) {
    stop("No testing tasks found")
  }

  # Define fit_one_task function for MESMA
  fit_one_task <- function(task) {
    loc <- task$loc
    yr <- task$yr

    res_safe <- tryCatch(
      {
        dly <- df_tasks[df_tasks$location_id == loc, , drop = FALSE]
        dly_year <- dly[lubridate::year(dly$date) == yr, , drop = FALSE]
        if (nrow(dly_year) == 0) return(NULL)
        if (nrow(dly) == 0) return(NULL)

        # Calculate Q10 and Q90 DVI
        dly_train <- df_tasks[df_tasks$location_id == loc & lubridate::year(df_tasks$date) == yr, , drop = FALSE]
        dvi_vals <- dly_train$DVI[is.finite(dly_train$DVI)]
        q10_dvi <- if (length(dvi_vals) > 0) stats::quantile(dvi_vals, 0.10, na.rm = TRUE) else NA
        q90_dvi <- if (length(dvi_vals) > 0) stats::quantile(dvi_vals, 0.90, na.rm = TRUE) else NA

        if (!factor_mode) return(NULL)

        # 1. Time Compression (BEFORE MESMA)
        comp_res <- compress_trace(dly_year, GLOBAL_PCA, avail, TEMPORAL_BUDGET)
        if (is.null(comp_res)) return(NULL)
        y <- comp_res$y
        
        # 2. STAGE 1 MESMA: Unmix vegetated fraction using compressed data
        veg_fraction_total <- NA_real_
        if (exists("COMPRESSED_STAGE1_LIB") && !is.null(COMPRESSED_STAGE1_LIB)) {
           veg_fraction_total <- unmix_stage1_compressed(y, COMPRESSED_STAGE1_LIB)
        } else {
           # Fallback to raw unmixing if compressed lib not available
           if (!is.null(STAGE1_LIB)) {
             veg_fraction_total <- unmix_vegetated_fraction(dly_year, STAGE1_LIB, avail)
           }
        }

        if (!is.finite(veg_fraction_total)) return(NULL)
        if (veg_fraction_total < 0 || veg_fraction_total > 1) return(NULL)
        
        barren_fraction_total <- max(0, min(1, 1 - veg_fraction_total))

        if (veg_fraction_total <= 0) {
          coef_df <- data.frame(
            location_id = loc,
            year = yr,
            Veg = "barren",
            coef = 1,
            rmse = NA_real_,
            stringsAsFactors = FALSE
          )
          coef_df$coef_lo <- NA_real_
          coef_df$coef_hi <- NA_real_
          diag_df <- data.frame(
            location_id = loc,
            year = yr,
            vegetated_fraction = veg_fraction_total,
            barren_fraction = barren_fraction_total,
            stringsAsFactors = FALSE
          )
          return(list(
            coef_df = coef_df,
            variant_trajectory = NULL,
            diagnostics = diag_df,
            uncertainty = NULL,
            q10_dvi = q10_dvi,
            q90_dvi = q90_dvi,
            vegetated_fraction = veg_fraction_total,
            barren_fraction = barren_fraction_total
          ))
        }

        veg_names <- names(lib)
        veg_kept <- intersect(ALLOWED_VEG, veg_names)

        if (length(veg_kept) < MAX_VEG_COMPONENTS && length(veg_counts) > 0) {
          global_order <- names(veg_counts)
          global_order <- intersect(global_order, veg_names)
          add <- setdiff(global_order, veg_kept)
          add <- add[tolower(add) != "barren"]
          if (length(add) > 0) {
            veg_kept <- c(veg_kept, add[seq_len(min(length(add), MAX_VEG_COMPONENTS - length(veg_kept)))])
          }
        }

        if (length(veg_kept) > MAX_VEG_COMPONENTS) {
          veg_kept <- veg_kept[seq_len(MAX_VEG_COMPONENTS)]
        }

        if (length(veg_kept) == 0) {
          return(dbg_return_null("no_veg_kept"))
        }

        # 3. STAGE 2 MESMA: Unmix vegetation using compressed data
        mesma_result <- unmix_stage2_compressed(y, comp_res$grid_type, mesma_lib[veg_kept], topK = TOPK_VARIANTS)
        if (is.null(mesma_result)) return(dbg_return_null("mesma_pca_failed"))

        coef_stage2 <- mesma_result$vegetation_proportions
        coef_df <- data.frame(
          location_id = loc,
          year = yr,
          Veg = names(coef_stage2),
          coef = as.numeric(coef_stage2) * veg_fraction_total,
          rmse = mesma_result$rmse,
          stringsAsFactors = FALSE
        )
        coef_df$coef_lo <- NA_real_
        coef_df$coef_hi <- NA_real_
        if (!is.null(mesma_result$uncertainty) && !is.null(mesma_result$uncertainty$coef_ci)) {
          ci_tbl <- mesma_result$uncertainty$coef_ci
          for (i in seq_len(nrow(coef_df))) {
            vname <- coef_df$Veg[i]
            row_ci <- ci_tbl[ci_tbl$Veg == vname, , drop = FALSE]
            if (nrow(row_ci) == 1) {
              coef_df$coef_lo[i] <- veg_fraction_total * row_ci$coef_lo[1]
              coef_df$coef_hi[i] <- veg_fraction_total * row_ci$coef_hi[1]
            }
          }
        }
        variant_info_pca <- if (!is.null(mesma_result$chosen_variants)) {
          vi <- as.list(mesma_result$chosen_variants)
          names(vi) <- paste0(names(vi), "_variant")
          data.frame(location_id = loc, year = yr, vi, stringsAsFactors = FALSE, check.names = FALSE)
        } else NULL
        diag_df <- if (!is.null(mesma_result$diagnostics)) {
          dd <- as.data.frame(mesma_result$diagnostics, stringsAsFactors = FALSE)
          dd$location_id <- loc; dd$year <- yr
          dd <- dd[, c("location_id", "year", setdiff(names(dd), c("location_id", "year"))) , drop = FALSE]
          dd
        } else {
          data.frame(location_id = loc, year = yr, stringsAsFactors = FALSE)
        }
        diag_df$vegetated_fraction <- veg_fraction_total
        diag_df$barren_fraction <- barren_fraction_total

        if (barren_fraction_total > 0) {
          barren_row <- data.frame(
            location_id = loc,
            year = yr,
            Veg = "barren",
            coef = barren_fraction_total,
            rmse = mesma_result$rmse,
            stringsAsFactors = FALSE
          )
          barren_row$coef_lo <- NA_real_
          barren_row$coef_hi <- NA_real_
          coef_df <- rbind(coef_df, barren_row)
        }

        unc <- mesma_result$uncertainty
        if (!is.null(unc) && !is.null(unc$coef_ci)) {
          unc$coef_ci$coef_lo <- veg_fraction_total * unc$coef_ci$coef_lo
          unc$coef_ci$coef_hi <- veg_fraction_total * unc$coef_ci$coef_hi
        }

        return(list(
          coef_df = coef_df,
          variant_trajectory = variant_info_pca,
          diagnostics = diag_df,
          uncertainty = unc,
          q10_dvi = q10_dvi,
          q90_dvi = q90_dvi,
          vegetated_fraction = veg_fraction_total,
          barren_fraction = barren_fraction_total
        ))
      },
      error = function(e) {
        dbg_return_null(paste0("error:", as.character(e$message)))
      }
    )
    res_safe
  }

  # Prepare task function environment
  required_globals <- list(
    df = df,
    lib = lib,
    GLOBAL_PCA = GLOBAL_PCA,
    lib_factor_pca = lib_factor_pca,
    mesma_lib = mesma_lib,
    avail = avail,
    LOWER_BND = LOWER_BND,
    EPS_SIGMA = EPS_SIGMA,
    MIN_OBS_FOR_BOOT = MIN_OBS_FOR_BOOT,
    MAX_VEG_COMPONENTS = MAX_VEG_COMPONENTS,
    veg_counts = veg_counts,
    OUT_DIR = OUT_DIR,
    factor_mode = factor_mode,
    df_tasks = df_tasks,
    prepare_factor_data = prepare_factor_data,
    compress_trace = compress_trace,
    unmix_stage1_compressed = unmix_stage1_compressed,
    unmix_stage2_compressed = unmix_stage2_compressed,
    COMPRESSED_STAGE1_LIB = if (exists("COMPRESSED_STAGE1_LIB")) COMPRESSED_STAGE1_LIB else NULL,
    STAGE1_LIB = if (exists("STAGE1_LIB")) STAGE1_LIB else NULL,
    # GLS / bootstrap helpers (exposed for worker env)
    gls_block_bootstrap = gls_block_bootstrap,
    solve_weights_gls = solve_weights_gls,
    estimate_block_size = estimate_block_size,
    solve_weights_constrained = solve_weights_constrained,
    solve_weights_simplex = solve_weights_simplex,
    cos_sim = cos_sim,
    BOOTSTRAP_B = BOOTSTRAP_B,
    ENABLE_UNCERTAINTY = ENABLE_UNCERTAINTY,
    GLSBB_MIN_BLOCK = GLSBB_MIN_BLOCK,
    GLSBB_MAX_BLOCK = GLSBB_MAX_BLOCK,
    ENABLE_QP_SOLVER = ENABLE_QP_SOLVER,
    compute_diagnostics = compute_diagnostics,
    project_to_simplex = project_to_simplex,
    dbg_return_null = dbg_return_null
  )

  env_task <- list2env(required_globals, parent = globalenv())
  environment(fit_one_task) <- env_task

  # Main processing loop
  cat("Starting main processing loop...\n")
  start_time <- Sys.time()
  results_list <- .run_map(task_list, fit_one_task)
  end_time <- Sys.time()
  processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(
    "Main processing loop finished in %.2f seconds (%.2f minutes)\n",
    processing_time, processing_time / 60
  ))
  cat(sprintf("Average time per task: %.2f seconds\n", processing_time / length(task_list)))

  # Process results and write to Excel
  cat("Processing results and writing to Excel files...\n")

  # Filter out null results
  results_list <- results_list[!sapply(results_list, is.null)]
  cat(sprintf("After filtering NULL results: %d results remaining\n", length(results_list)))

  if (length(results_list) == 0) {
    cat("ERROR: All tasks returned NULL results!\n")
    cat("Most likely causes:\n")
    cat("1. No valid testing/inference data available (no location-year pairs found)\n")
    cat("2. Data filtering issues or missing indices\n")
    cat("3. All locations have insufficient data for fitting\n")
    stop("No valid results to process")
  }

  if (length(results_list) > 0) {
    # Combine all coefficient data frames
    cat("Combining coefficient data frames...\n")
    coef_list <- lapply(results_list, function(res) {
      if (is.null(res$coef_df) || nrow(res$coef_df) == 0) {
        NULL
      } else {
        res$coef_df
      }
    })
    coef_list <- coef_list[!sapply(coef_list, is.null)]

    if (length(coef_list) == 0) {
      cat("ERROR: No valid coefficient data frames found!\n")
      stop("No coefficient data to process")
    }

    all_coefs <- tryCatch(
      {
        do.call(rbind, coef_list)
      },
      error = function(e) stop(sprintf("ERROR combining coef_df: %s", e$message))
    )

    if (is.null(all_coefs)) {
      cat("Failed to combine coefficient data frames\n")
      stop("Cannot proceed without coefficient data")
    }

    cat(sprintf("Combined coefficients: %d rows\n", nrow(all_coefs)))

  # Combine chosen variant summaries (no per-DOY trajectories in new approach)
    cat("Combining chosen variant summaries...\n")
    variant_list_pca <- lapply(results_list, function(res) res$variant_trajectory$pca)
    variant_list_pca <- variant_list_pca[!sapply(variant_list_pca, is.null)]
    all_variants_pca <- if (length(variant_list_pca) > 0) tryCatch({ do.call(rbind, variant_list_pca) }, error = function(e) stop(sprintf("Failed to combine variant PCA summaries: %s", e$message))) else NULL

    # Collect diagnostics
    diag_list <- lapply(results_list, function(res) {
      if (!is.null(res$diagnostics)) res$diagnostics else NULL
    })
    diag_list <- diag_list[!sapply(diag_list, is.null)]
    # Use dplyr::bind_rows so differing diagnostic column sets are safely merged (missing columns filled with NA)
    all_diagnostics <- if (length(diag_list) > 0) tryCatch({
      if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required to combine diagnostics robustly")
      dplyr::bind_rows(diag_list)
    }, error = function(e) stop(sprintf("Failed to combine diagnostics: %s", e$message))) else NULL

    # Collect Q10 and Q90 DVI values
    q_dvi_data <- do.call(rbind, lapply(results_list, function(res) {
      if ((!is.null(res$q10_dvi) && !is.na(res$q10_dvi)) ||
        (!is.null(res$q90_dvi) && !is.na(res$q90_dvi))) {
        data.frame(
          location_id = res$coef_df$location_id[1],
          year = res$coef_df$year[1],
          q10_dvi = res$q10_dvi,
          q90_dvi = res$q90_dvi,
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    }))

    # Get unique locations
    unique_locations <- unique(all_coefs$location_id)

    # Join with ground truth vegetation data
    true_veg_map <- gpts_map %>% dplyr::select(location_id, true_veg = Veg)
    # Coerce types to avoid incompatible-type join errors
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

    # Create Excel file with all locations as separate sheets
    cat("Creating single Excel file with all location results...\n")

    wb <- tryCatch(
      {
        openxlsx::createWorkbook()
      },
      error = function(e) {
        cat(sprintf("ERROR creating workbook: %s\n", e$message))
        return(NULL)
      }
    )

    if (is.null(wb)) {
      cat("Failed to create workbook\n")
      stop("Cannot create Excel output")
    }

    # Add summary sheet
    summary_data <- data.frame(
      Location_ID = unique_locations,
      Total_Years = sapply(unique_locations, function(loc) {
        length(unique(all_coefs$year[all_coefs$location_id == loc]))
      }),
      Total_Observations = sapply(unique_locations, function(loc) {
        nrow(all_coefs[all_coefs$location_id == loc, ])
      }),
      # Projection method column not included: a single method is used per run
      stringsAsFactors = FALSE
    )

    openxlsx::addWorksheet(wb, "Summary")
    openxlsx::writeData(wb, "Summary", summary_data)

    # Diagnostics sheet
    if (!is.null(all_diagnostics) && nrow(all_diagnostics) > 0) {
      openxlsx::addWorksheet(wb, "Diagnostics")
      openxlsx::writeData(wb, "Diagnostics", all_diagnostics)
    }

    # Add MESMA Variant Summary sheet if variant data exists (chosen variants only)
    if (!is.null(all_variants_pca) && nrow(all_variants_pca) > 0) {
      openxlsx::addWorksheet(wb, "Variant_Summary")
      openxlsx::writeData(wb, "Variant_Summary", all_variants_pca)
    }

    # Collect uncertainty across locations if available
    unc_coef_rows <- list(); unc_var_rows <- list(); unc_rmse_rows <- list(); unc_meta_rows <- list()
    for (res in results_list) {
      if (is.null(res$uncertainty)) next
      loc <- if (!is.null(res$coef_df$location_id)) res$coef_df$location_id[1] else NA_character_
      yr <- if (!is.null(res$coef_df$year)) res$coef_df$year[1] else NA_integer_
      ci <- res$uncertainty$coef_ci
      if (!is.null(ci) && nrow(ci) > 0) {
        ci$location_id <- loc; ci$year <- yr
        unc_coef_rows[[length(unc_coef_rows) + 1]] <- ci[, c("location_id","year","Veg","coef_lo","coef_hi"), drop = FALSE]
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
      if (!is.null(res$uncertainty$block_size)) {
        unc_meta_rows[[length(unc_meta_rows) + 1]] <- data.frame(location_id = loc, year = yr, block_size = res$uncertainty$block_size, B = BOOTSTRAP_B, stringsAsFactors = FALSE)
      }
    }
    all_unc_coef <- if (length(unc_coef_rows) > 0) do.call(rbind, unc_coef_rows) else NULL
    all_unc_var <- if (length(unc_var_rows) > 0) do.call(rbind, unc_var_rows) else NULL
    all_unc_rmse <- if (length(unc_rmse_rows) > 0) do.call(rbind, unc_rmse_rows) else NULL
    all_unc_meta <- if (length(unc_meta_rows) > 0) do.call(rbind, unc_meta_rows) else NULL

    if (!is.null(all_unc_coef) || !is.null(all_unc_var) || !is.null(all_unc_rmse)) {
      openxlsx::addWorksheet(wb, "Uncertainty")
      start_row <- 1
      openxlsx::writeData(wb, "Uncertainty", data.frame(Setting = c("ENABLE_UNCERTAINTY","BOOTSTRAP_B"), Value = c(ENABLE_UNCERTAINTY, BOOTSTRAP_B)), startRow = start_row, startCol = 1)
      start_row <- start_row + 3
      if (!is.null(all_unc_coef)) {
        openxlsx::writeData(wb, "Uncertainty", "Coefficient CIs (2.5%/97.5%)", startRow = start_row, startCol = 1)
        openxlsx::writeData(wb, "Uncertainty", all_unc_coef, startRow = start_row + 1, startCol = 1)
        start_row <- start_row + nrow(all_unc_coef) + 3
      }
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

    # Build best fit summary per location-year (single projection per run)
    best_fit_summary <- do.call(rbind, lapply(unique_locations, function(loc) {
      yrs <- unique(all_coefs$year[all_coefs$location_id == loc])
      tv <- true_veg_map$true_veg[true_veg_map$location_id == loc]
      if (length(tv) == 0) tv <- NA_character_

      do.call(rbind, lapply(yrs, function(yr) {
        row <- all_coefs[all_coefs$location_id == loc & all_coefs$year == yr & tolower(all_coefs$Veg) == tolower(tv), , drop = FALSE]
        pred <- if (nrow(row) == 1) row$coef else NA_real_
        rmse_val <- if (nrow(row) == 1 && "rmse" %in% names(row)) row$rmse else NA_real_
        abs_pct <- if (!is.na(pred) && !is.na(tv)) abs(1 - pred) * 100 else NA_real_

        data.frame(
          location_id = loc,
          year = yr,
          true_veg = tv,
          pred_coef = pred,
          rmse = rmse_val,
          abs_pct_diff = abs_pct,
          stringsAsFactors = FALSE
        )
      }))
    }))

    # Compute overall fit score (single projection)
    overall_fit <- suppressWarnings(as.numeric(mean(best_fit_summary$abs_pct_diff, na.rm = TRUE)))
    if (!is.finite(overall_fit)) overall_fit <- NA_real_

    # Write overall fit score at top of Summary sheet
    openxlsx::writeData(wb, "Summary", data.frame(
      Overall_Fit_pct = overall_fit
    ), startRow = 1, startCol = ncol(summary_data) + 2)

    # Add location sheets
    if (requireNamespace("progressr", quietly = TRUE)) {
      progressr::with_progress({
        p_excel <- progressr::progressor(
          steps = length(unique_locations),
          message = "Adding location sheets"
        )

        for (i in seq_along(unique_locations)) {
          loc_id <- unique_locations[i]

          # Create sheet name (Excel limits sheet names to 31 characters)
          sheet_name <- substr(gsub("[^A-Za-z0-9]", "_", loc_id), 1, 31)

          # Add worksheet for this location
          openxlsx::addWorksheet(wb, sheet_name)

          # Quality Metrics Calculation
          loc_coefs <- all_coefs[all_coefs$location_id == loc_id, ]

          quality_metrics <- loc_coefs %>%
            dplyr::group_by(.data$year) %>%
            dplyr::summarize(
              deviation = sum(abs(.data$coef - (tolower(.data$Veg) == tolower(true_veg[1])))),
              avg_rmse = mean(.data$rmse, na.rm = TRUE),
              .groups = "drop"
            ) %>%
            dplyr::summarize(
              avg_pct_deviation = mean(.data$deviation, na.rm = TRUE) * 100,
              avg_rmse = mean(.data$avg_rmse, na.rm = TRUE),
              .groups = "drop"
            )

          # Calculate peak Q10 and Q90 DVI for this location
          loc_q_data <- if (!is.null(q_dvi_data)) {
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

          quality_metrics$peak_q10_dvi <- peak_q10_dvi
          quality_metrics$peak_q90_dvi <- peak_q90_dvi

          # Write quality metrics to the sheet
          openxlsx::writeData(wb, sheet_name, "QUALITY METRICS", startRow = 1, startCol = 1)
          openxlsx::writeData(wb, sheet_name, quality_metrics, startRow = 2, startCol = 1)

          # Data Sections
          current_row <- nrow(quality_metrics) + 4

          # Optional: per-location diagnostics
          loc_diag <- if (exists("all_diagnostics") && !is.null(all_diagnostics)) all_diagnostics[all_diagnostics$location_id == loc_id, , drop = FALSE] else NULL
          if (!is.null(loc_diag) && nrow(loc_diag) > 0) {
            openxlsx::writeData(wb, sheet_name, "DIAGNOSTICS", startRow = current_row, startCol = 1)
            openxlsx::writeData(wb, sheet_name, loc_diag, startRow = current_row + 1, startCol = 1)
            current_row <- current_row + nrow(loc_diag) + 3
          }

          # BEST FIT SUMMARY for this location
          loc_best <- best_fit_summary[best_fit_summary$location_id == loc_id, , drop = FALSE]
          if (!is.null(loc_best) && nrow(loc_best) > 0) {
            desired_cols <- c("year", "true_veg", "pred_coef", "rmse", "abs_pct_diff")
            write_tbl <- loc_best[, c("location_id", intersect(desired_cols, names(loc_best))), drop = FALSE]
            openxlsx::writeData(wb, sheet_name, "BEST FIT SUMMARY (per-year)",
              startRow = current_row, startCol = 1
            )
            openxlsx::writeData(wb, sheet_name, write_tbl, startRow = current_row + 1, startCol = 1)
            current_row <- current_row + nrow(write_tbl) + 3
          }

          # Filter data for the current location
          # Unified coefficients for this location (single projection per run)
          loc_coefs_unified <- all_coefs[all_coefs$location_id == loc_id, ]

          loc_variants_pca <- if (!is.null(all_variants_pca)) {
            all_variants_pca[all_variants_pca$location_id == loc_id, ]
          } else {
            NULL
          }

          # Add MESMA PCA coefficients
          if (nrow(loc_coefs_unified) > 0) {
            openxlsx::writeData(wb, sheet_name, "COEFFICIENTS",
              startRow = current_row, startCol = 1
            )
            openxlsx::writeData(wb, sheet_name, loc_coefs_unified,
              startRow = current_row + 1, startCol = 1
            )
            current_row <- current_row + nrow(loc_coefs_unified) + 3
          }

          # If uncertainty coefficients exist, write them next
          if (exists("all_unc_coef") && !is.null(all_unc_coef)) {
            loc_unc_coef <- all_unc_coef[all_unc_coef$location_id == loc_id, , drop = FALSE]
            if (nrow(loc_unc_coef) > 0) {
              openxlsx::writeData(wb, sheet_name, "COEFFICIENT CIs (2.5%/97.5%)",
                startRow = current_row, startCol = 1
              )
              openxlsx::writeData(wb, sheet_name, loc_unc_coef,
                startRow = current_row + 1, startCol = 1
              )
              current_row <- current_row + nrow(loc_unc_coef) + 3
            }
          }

          # Add MESMA PCA variant trajectory summary
          if (!is.null(loc_variants_pca) && nrow(loc_variants_pca) > 0) {
            # Create summary of variant usage
            variant_usage <- data.frame(
              year = unique(loc_variants_pca$year),
              stringsAsFactors = FALSE
            )

            for (veg in veg_kept) {
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

          # If uncertainty variant dominance exists, write it
          if (exists("all_unc_var") && !is.null(all_unc_var)) {
            loc_unc_var <- all_unc_var[all_unc_var$location_id == loc_id, , drop = FALSE]
            if (nrow(loc_unc_var) > 0) {
              openxlsx::writeData(wb, sheet_name, "VARIANT DOMINANCE (%)",
                startRow = current_row, startCol = 1
              )
              openxlsx::writeData(wb, sheet_name, loc_unc_var,
                startRow = current_row + 1, startCol = 1
              )
              current_row <- current_row + nrow(loc_unc_var) + 3
            }
          }

          # Add Q10/Q90 DVI data if available
          if (!is.null(loc_q_data) && nrow(loc_q_data) > 0) {
            openxlsx::writeData(wb, sheet_name, "Q10/Q90 DVI TREND",
              startRow = current_row, startCol = 1
            )
            openxlsx::writeData(wb, sheet_name, loc_q_data,
              startRow = current_row + 1, startCol = 1
            )
            current_row <- current_row + nrow(loc_q_data) + 3
          }

          # If RMSE CI exists for this location
          if (exists("all_unc_rmse") && !is.null(all_unc_rmse)) {
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

          p_excel(sprintf("Added sheet for %s (%d/%d)", loc_id, i, length(unique_locations)))
        }
      })
    } else {
      stop("progressr package required for Excel generation")
    }

    # Save the workbook
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

  # Final timing summary
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

  cat("\nMESMA fitting completed successfully!\n")

print("Script fit_veg_mixture_mesma.R execution finished.")