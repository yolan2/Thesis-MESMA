library(zoo)
library(dplyr)
library(cluster)

suppressPackageStartupMessages({
  suppressWarnings({
    library(zoo)
    library(dplyr)
    library(cluster)
  })
})


# Declare known globals for R CMD check / static linters
utils::globalVariables(c(
  "GLOBAL_PCA", "PARALLEL_ENABLE", "PARALLEL_WORKERS",
    "PROGRESS_LOG_TO_FILE", "PROGRESS_EVERY_TASK", "OUT_DIR", "ALLOWED_VEG",
    "MAX_VEG_COMPONENTS", "MIN_IDX_PRESENCE", "lib_factor_pca", "veg_counts", "mesma_lib", "FACTOR_MODE",
    "veg_kept", "switches", "all_variants_pca", "variant_usage",
    "unique_locations", "results_list", "best_fit_summary",
    "variant_list_pca", "loc_variants_pca",
    "FitMethod", "true_veg", "rmse", "deviation", "avg_rmse", "q_dvi_data",
    "loc_q_data", "loc_coefs", "loc_best", "summary_data",
    "EPS_SIGMA", "LOWER_BND", "MIN_OBS_FOR_BOOT", "BOOTSTRAP_B",
    "TEST_YEARS", "TRAIN_YEARS", "PROGRESS_BAR", "LOG_FILE",
    "V_names", "pca_rank",
    "t_pca_med", "variant_t_pca", "all_coefs",
    "location_id", "year", "Veg", "coef",
    "MAX_PCA_COMPONENTS", "MAX_FACTORS_CAP", "FAST_VAR",
    "BOOT_MIN_REPS_PER_VEG", "VARIANCE_THRESHOLD",
    "GAM_K_MAX", "GAM_GAMMA", "USE_INDICES_MIN", "MIN_INDEX_SD",
    "ENABLE_SAMPLE_BALANCING", "RAW_BANDS",
    "ENABLE_PHASE_ALIGNMENT", "REFERENCE_PHASE_MARKERS",
    "ENABLE_MULTISCALE", "MULTISCALE_WINDOWS",
    "ENABLE_DIAGNOSTICS",
    "PERSISTENT_PARALLEL_BACKEND",
    ".COMPRESSED_TEMPLATES_ACCESSOR", ".WHITENING_DB_ACCESSOR",
    "ENABLE_UNCERTAINTY", "BOOTSTRAP_B", "GLSBB_MIN_BLOCK", "GLSBB_MAX_BLOCK",
    "COMBO_PARALLEL_ENABLE", "COMBO_PARALLEL_WORKERS", "EARLY_STOP_RMSE_THRESHOLD",
    "ENABLE_QP_SOLVER", "TEMPORAL_BUDGET", "TOPK_VARIANTS"
  ))



svd_mkl <- function(mat, nu = 0, nv = NULL) {
  # High-performance SVD using R's optimized BLAS/LAPACK (Intel MKL)
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if (is.null(nv)) nv <- min(ncol(mat), nrow(mat))
  if (any(!is.finite(mat))) stop("svd_mkl: input matrix contains non-finite values; clean data before calling")
  
  # Use R's built-in SVD which will leverage Intel MKL automatically
  svd_res <- svd(mat, nu = nu, nv = nv)
  list(d = svd_res$d, u = svd_res$u, v = svd_res$v)
}

crossprod_mkl <- function(A, B = NULL) {
  # High-performance crossprod using R's optimized BLAS (Intel MKL)
  if (!is.matrix(A)) A <- as.matrix(A)
  if (!is.null(B) && !is.matrix(B)) B <- as.matrix(B)
  if (any(!is.finite(A)) || (!is.null(B) && any(!is.finite(B)))) stop("crossprod_mkl: inputs contain non-finite values; clean data before calling")
  
  # Use R's built-in crossprod which will leverage Intel MKL automatically
  if (is.null(B)) {
    crossprod(A)
  } else {
    crossprod(A, B)
  }
}

cor_mkl <- function(df, use = "pairwise.complete.obs") {
  # High-performance correlation using R's optimized functions (Intel MKL)
  if (!is.data.frame(df) && !is.matrix(df)) df <- as.data.frame(df)
  X <- as.matrix(df)
  if (any(!is.finite(X))) stop("cor_mkl: input contains non-finite values; clean or impute before calling")
  
  # Use R's built-in cor function which will leverage Intel MKL for matrix operations
  cor(X, use = use)
}

kmeans_mkl <- function(x, centers, nstart = 1, max_iter = 100, tol = 1e-6) {
  # High-performance k-means using R's built-in functions (Intel MKL optimized)
  if (!is.matrix(x)) x <- as.matrix(x)
  if (any(!is.finite(x))) stop("kmeans_mkl: input contains non-finite values; clean data before calling")
  
  # Use R's built-in kmeans which will leverage Intel MKL for matrix operations
  if (is.numeric(centers) && length(centers) == 1) {
    # centers is the number of clusters
    km_res <- kmeans(x, centers = centers, nstart = nstart, iter.max = max_iter)
  } else if (is.matrix(centers)) {
    # centers is initial cluster centers
    km_res <- kmeans(x, centers = centers, nstart = nstart, iter.max = max_iter)
  } else {
    stop("kmeans_mkl: invalid centers")
  }
  
  return(km_res)
}

fit_cost_mkl <- function(obs, weight, t_row) {
  # High-performance cost calculation using vectorized operations (Intel MKL optimized)
  if (length(obs) != length(t_row)) t_row <- t_row[seq_len(length(obs))]
  if (any(!is.finite(c(obs, weight, t_row)))) stop("fit_cost_mkl: inputs contain non-finite values; clean data before calling")
  
  # Use vectorized operations which will leverage Intel MKL
  pred <- weight * t_row
  residuals <- (obs * weight) - pred
  sum(residuals^2)
}

# Set up progressr handlers for ETA display
if (requireNamespace("progressr", quietly = TRUE)) {
  progressr::handlers(global = TRUE)
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

# Training/Test configuration defaults
if (!exists("TEST_YEARS")) TEST_YEARS <- NULL
if (!exists("TRAIN_YEARS")) TRAIN_YEARS <- 2019:2024

# Allow skipping a minimal number of DOYs per-location when computing sufficiency
# This is useful for small gaps — set to 0 to preserve the strict 50% requirement
if (!exists("MIN_SKIP_DOYS_PER_LOCATION")) MIN_SKIP_DOYS_PER_LOCATION <- 2L

# PCA / factor projection tuning
if (!exists("MAX_PCA_COMPONENTS")) MAX_PCA_COMPONENTS <- 20L
if (!exists("MAX_FACTORS_CAP")) MAX_FACTORS_CAP <- 20L

# General algorithm toggles
if (!exists("FAST_VAR")) FAST_VAR <- TRUE
if (!exists("factor_mode")) factor_mode <- FALSE
if (!exists("SHAPE_NORMALIZATION_ENABLE")) SHAPE_NORMALIZATION_ENABLE <- FALSE
if (!exists("TEMPORAL_BUDGET")) TEMPORAL_BUDGET <- 10L
if (!exists("TOPK_VARIANTS")) TOPK_VARIANTS <- 2L
if (!exists("ENABLE_PHASE_ALIGNMENT")) ENABLE_PHASE_ALIGNMENT <- FALSE
if (!exists("REFERENCE_PHASE_MARKERS")) REFERENCE_PHASE_MARKERS <- c(1, 90, 180, 270, 365)
if (!exists("ENABLE_MULTISCALE")) ENABLE_MULTISCALE <- FALSE
if (!exists("MULTISCALE_WINDOWS")) MULTISCALE_WINDOWS <- c(7L, 14L, 30L)
if (!exists("ENABLE_QP_SOLVER")) ENABLE_QP_SOLVER <- TRUE
if (!exists("COMBO_PARALLEL_ENABLE")) COMBO_PARALLEL_ENABLE <- FALSE
if (!exists("COMBO_PARALLEL_WORKERS")) COMBO_PARALLEL_WORKERS <- max(1L, if (exists("PARALLEL_WORKERS")) floor(PARALLEL_WORKERS/2) else 1L)
if (!exists("EARLY_STOP_RMSE_THRESHOLD")) EARLY_STOP_RMSE_THRESHOLD <- 0.0
if (!exists("ENABLE_DIAGNOSTICS")) ENABLE_DIAGNOSTICS <- TRUE

# Bootstrap settings 
if (!exists("BOOTSTRAP_B")) BOOTSTRAP_B <- 200L
if (!exists("BOOT_MIN_REPS_PER_VEG")) BOOT_MIN_REPS_PER_VEG <- 10L
if (!exists("MIN_OBS_FOR_BOOT")) MIN_OBS_FOR_BOOT <- 5L
if (!exists("ENABLE_UNCERTAINTY")) ENABLE_UNCERTAINTY <- TRUE
if (!exists("GLSBB_MIN_BLOCK")) GLSBB_MIN_BLOCK <- 5L
if (!exists("GLSBB_MAX_BLOCK")) GLSBB_MAX_BLOCK <- 60L

# Optimization/solver defaults
if (!exists("VARIANCE_THRESHOLD")) VARIANCE_THRESHOLD <- 0.90
if (!exists("MAX_VEG_COMPONENTS")) MAX_VEG_COMPONENTS <- 8
if (!exists("GAM_K_MAX")) GAM_K_MAX <- 40
if (!exists("GAM_GAMMA")) GAM_GAMMA <- 1.0

# Index selection and prefiltering
if (!exists("USE_INDICES_MIN")) USE_INDICES_MIN <- 1L
if (!exists("MIN_INDEX_SD")) MIN_INDEX_SD <- 0.05

# Sample-balancing and augmentation
if (!exists("ENABLE_SAMPLE_BALANCING")) ENABLE_SAMPLE_BALANCING <- TRUE

# Vegetation whitelist
if (!exists("ALLOWED_VEG")) ALLOWED_VEG <- c("populus", "tamarix", "phragmites")

# Numeric safety constants
if (!exists("EPS_SIGMA")) EPS_SIGMA <- 1e-8
if (!exists("LOWER_BND")) LOWER_BND <- 0

# DOY / index presence tuning
if (!exists("MIN_IDX_PRESENCE")) MIN_IDX_PRESENCE <- 0.5

# Parallel / progress settings
if (!exists("PROGRESS_EVERY_TASK")) PROGRESS_EVERY_TASK <- 25
if (!exists("PROGRESS_LOG_TO_FILE")) PROGRESS_LOG_TO_FILE <- TRUE
if (!exists("PROGRESS_BAR")) PROGRESS_BAR <- TRUE
if (!exists("PARALLEL_ENABLE")) PARALLEL_ENABLE <- TRUE
if (!exists("PARALLEL_WORKERS")) {
  PARALLEL_WORKERS <- tryCatch(
    {
      if (requireNamespace("parallel", quietly = TRUE)) max(1L, parallel::detectCores(logical = TRUE) - 1L) else 1L
    },
    error = function(...) 1L
  )
}
if (!exists("PERSISTENT_PARALLEL_BACKEND")) PERSISTENT_PARALLEL_BACKEND <- TRUE

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
cat("\n=== TRAIN/TEST SPLIT CONFIGURATION ===\n")
cat(sprintf("Training years: %s\n", paste(TRAIN_YEARS, collapse = ", ")))
if (!is.null(TEST_YEARS)) {
  cat(sprintf("Testing years: %s\n", paste(TEST_YEARS, collapse = ", ")))
} else {
  cat("Testing years: All available years\n")
}

df$year <- lubridate::year(df$date)

df_train <- df[df$year %in% TRAIN_YEARS, , drop = FALSE]
cat(sprintf(
  "Training dataset: %d rows from %d locations\n",
  nrow(df_train), length(unique(df_train$location_id))
))

if (!is.null(TEST_YEARS)) {
  df_test <- df[df$year %in% TEST_YEARS, , drop = FALSE]
} else {
  df_test <- df
}
cat(sprintf(
  "Testing dataset: %d rows from %d locations\n",
  nrow(df_test), length(unique(df_test$location_id))
))

df_full <- df
df <- df_train

# Ensure training subset has at least some samples for each allowed vegetation class
# If a class from ALLOWED_VEG has zero samples in the selected TRAIN_YEARS, try to
# add a small number of samples for that class from the full dataset (non-TRAIN years)
if (exists("df_full") && "Veg" %in% names(df_full) && length(ALLOWED_VEG) > 0) {
  missing_vegs <- sapply(ALLOWED_VEG, function(v) {
    sum(tolower(df$Veg) == v, na.rm = TRUE)
  })
  missing_names <- names(missing_vegs)[missing_vegs == 0]
  if (length(missing_names) > 0) {
    for (mv in missing_names) {
      # find candidates in the full dataset (outside TRAIN_YEARS)
      cand <- df_full[tolower(df_full$Veg) == mv & !(df_full$year %in% TRAIN_YEARS), , drop = FALSE]
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
if (!exists("PROGRESS_EVERY_TASK")) PROGRESS_EVERY_TASK <- 25
if (!exists("PROGRESS_LOG_TO_FILE")) PROGRESS_LOG_TO_FILE <- TRUE
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
        # Fixed: Remove 'projection' from key since we're using single method
        key <- paste(veg, variant$variant_id, grid_name, sep = "|")
        assign(key, feats, envir = template_db)
      }
    }
  }
  # Fixed: Remove 'projection' parameter from accessor function
  function(veg, variant_id, grid_type = "medium") {
      key <- paste(veg, variant_id, grid_type, sep = "|")
      if (exists(key, envir = template_db, inherits = FALSE)) get(key, envir = template_db, inherits = FALSE) else NULL
  }
}

# Helper: index of medoid row (row with minimal total distance)
medoid_row_index <- function(M) {
  if (is.null(M) || !is.matrix(M) || nrow(M) == 0) return(NA_integer_)
  if (nrow(M) == 1) return(1L)
  X <- M
  X[!is.finite(X)] <- 0
  D <- tryCatch(as.matrix(stats::dist(X)), error = function(e) stop(sprintf("medoid_row_index: distance computation failed: %s", e$message)))
  # dist returns lower triangular; as.matrix rebuilds full symmetric
  rs <- rowSums(D, na.rm = TRUE)
  idx <- which.min(rs)
  if (length(idx) == 0 || !is.finite(idx)) idx <- 1L
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

geojson_path <- file.path("C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/fotos", "zuiver_zonder_foto.geojson")
if (!file.exists(geojson_path)) stop(paste0("GeoJSON points not found at ", geojson_path))
gpts_raw <- sf::st_read(geojson_path, quiet = TRUE)

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
gpts_raw$location_id <- make_location_id(gpts_raw$.__lon__, gpts_raw$.__lat__)

gpts_map <- sf::st_drop_geometry(gpts_raw) %>%
  dplyr::select(location_id, Veg = .__veg__) %>%
  dplyr::distinct(location_id, .keep_all = TRUE)

if (nrow(gpts_map) == 0) stop("GeoJSON mapping produced no valid points")

if ("location_id" %in% names(df) && "location_id" %in% names(gpts_map)) {
  # Ensure both sides are character to avoid incompatible-type join errors
  # If df has numeric IDs but gpts_map uses 'L_x_y' format, prefer building location_id from lon/lat when available
  if (!is.character(df$location_id)) {
  sample_df <- unique(na.omit(as.character(df$location_id)))
  sample_g  <- unique(na.omit(as.character(gpts_map$location_id)))
    if (length(sample_df) && length(sample_g) && all(grepl("^\\s*[0-9]+\\s*$", sample_df)) && any(grepl("^L_", sample_g))) {
      # Attempt to reconstruct df$location_id from lon/lat columns if available
      lon_candidates_local <- names(df)[grepl("(^|_)(lon|longitude|x)(_|$)", names(df), ignore.case = TRUE)]
      lat_candidates_local <- names(df)[grepl("(^|_)(lat|latitude|y)(_|$)", names(df), ignore.case = TRUE)]
      if (length(lon_candidates_local) > 0 && length(lat_candidates_local) > 0) {
        lon_col_local <- lon_candidates_local[1]
        lat_col_local <- lat_candidates_local[1]
        df$location_id <- make_location_id(df[[lon_col_local]], df[[lat_col_local]])
        cat("[NOTICE] Recreated df$location_id from lon/lat to match gpts_map format before joining.\n")
      } else {
        df$location_id <- as.character(df$location_id)
        cat("[NOTICE] Coerced df$location_id to character to match gpts_map for joining.\n")
        cat("[WARNING] df$location_id looks numeric while gpts_map$location_id looks like formatted 'L_x_y' strings — join will not match without remapping (e.g., use lon/lat or consistent IDs).\n")
      }
    } else {
      df$location_id <- as.character(df$location_id)
      cat("[NOTICE] Coerced df$location_id to character to match gpts_map for joining.\n")
    }
  }
  if (!is.character(gpts_map$location_id)) {
    gpts_map$location_id <- as.character(gpts_map$location_id)
    cat("[NOTICE] Coerced gpts_map$location_id to character for joining.\n")
  }
  # Warn if the values are likely incompatible (numeric IDs versus L_x_y strings)
  sample_df <- unique(na.omit(df$location_id))
  sample_g <- unique(na.omit(gpts_map$location_id))
  if (length(sample_df) && length(sample_g) && all(grepl("^[0-9]+$", sample_df)) && any(grepl("^L_", sample_g))) {
    cat("[WARNING] df$location_id looks numeric while gpts_map$location_id looks like formatted 'L_x_y' strings — join will not match these values.\n")
  }
  df <- dplyr::left_join(df, gpts_map, by = "location_id")
}

loc_years <- data.frame(location_id = character(0), year = integer(0), stringsAsFactors = FALSE)

lon_candidates <- names(df)[grepl("(^|_)(lon|longitude|x)(_|$)", names(df), ignore.case = TRUE)]
lat_candidates <- names(df)[grepl("(^|_)(lat|latitude|y)(_|$)", names(df), ignore.case = TRUE)]
if (length(lon_candidates) > 0 && length(lat_candidates) > 0) {
  lon_col <- lon_candidates[1]
  lat_col <- lat_candidates[1]
  df$location_id <- make_location_id(df[[lon_col]], df[[lat_col]])
}

if (!"year" %in% names(df)) {
  if ("date" %in% names(df)) {
    if (!requireNamespace("lubridate", quietly = TRUE)) stop("The package 'lubridate' is required")
    df$year <- as.integer(lubridate::year(as.Date(df$date)))
  }
}

if (!"Veg" %in% names(df)) df$Veg <- NA_character_

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
  corr_matrix <- cor_mkl(corr_data, use = "pairwise.complete.obs")
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
  # Previously this enforced 50% DOY threshold and removed/flagged insufficient locations
  if (FALSE && !exists("TRAIN_YEARS") || is.null(TRAIN_YEARS) || length(TRAIN_YEARS) == 0) stop("TRAIN_YEARS must be defined and non-empty for data sufficiency checks; no fallback permitted.")
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

      # Previously offending locations were removed here. To preserve all locations
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
    if (!exists("TRAIN_YEARS") || is.null(TRAIN_YEARS) || length(TRAIN_YEARS) == 0) stop("TRAIN_YEARS must be defined for global sufficiency checks; no fallback permitted.")
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
    
    # Convert to matrix before calling cor_mkl
    vmat_matrix <- as.matrix(vmat_clean)
    
    # Then use the cleaned matrix
    cm <- cor_mkl(vmat_matrix, use = "pairwise.complete.obs")
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
if (!exists("MIN_INDEX_SD")) MIN_INDEX_SD <- 0.05
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
  keep_rows <- tolower(df$Veg) %in% ALLOWED_VEG
  n_before <- nrow(df)
  df <- df[keep_rows | is.na(df$Veg), , drop = FALSE]
  cat(sprintf(
    "Filtered to allowed classes (%s): kept %d/%d rows\n",
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

# Recompute standardized location_id
lon_candidates <- names(df)[grepl("(^|_)(lon|longitude|x)(_|$)", names(df), ignore.case = TRUE)]
lat_candidates <- names(df)[grepl("(^|_)(lat|latitude|y)(_|$)", names(df), ignore.case = TRUE)]
if (length(lon_candidates) > 0 && length(lat_candidates) > 0) {
  lon_col <- lon_candidates[1]
  lat_col <- lat_candidates[1]
  df$location_id <- make_location_id(df[[lon_col]], df[[lat_col]])
}

loc_years <- df %>%
  dplyr::filter(!is.na(.data$location_id) & .data$location_id != "" & !is.na(.data$year) & .data$year > 0) %>%
  dplyr::distinct(.data$location_id, .data$year)
cat(sprintf("Constructed loc_years with %d rows from filtered df\n", nrow(loc_years)))
if (nrow(loc_years) == 0) {
  # Try a last-ditch fallback: use the full (unfiltered) dataframe to construct location-year pairs
  alt_loc_years <- df_full %>%
    dplyr::filter(!is.na(.data$location_id) & .data$location_id != "" & !is.na(.data$year) & .data$year > 0) %>%
    dplyr::distinct(.data$location_id, .data$year)
  if (nrow(alt_loc_years) > 0) {
    cat(sprintf("[NOTICE] No location-year pairs found in filtered training data; falling back to full dataset and using %d loc-year pairs for library construction.\n", nrow(alt_loc_years)))
    loc_years <- alt_loc_years
  } else {
    stop(paste0(
      "No location-year pairs found after filtering, and none were available in the full dataset.\n",
      "Suggestions:\n",
      " - Verify TRAIN_YEARS is correct and includes the years you expect (e.g. TRAIN_YEARS <- 2019:2024).\n",
      " - If your sampling cadence is sparse, consider adjusting training years or the fractional threshold (e.g. reduce MIN_DAYS_FRACTION from 0.5).\n",
      " - Re-run the transformer to ensure 'location_id' values align with your geojson mapping (transform_phenology.py ensures 'L_lon_lat' formatting).\n",
      "Processing cannot continue without at least one location-year pair.")
    )
  }
}

# Construct vegetation library
cat("Constructing lib from TRAINING dataset...\n")
lib <- list()
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
      new_row[1, num_cols] <- interp
    }

    for (col in other_cols) {
      if (col %in% c("date", "doy", "doy_for_lib")) next
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
  cat(sprintf(
    "Training data balanced: %d total samples across %d vegetation types\n",
    nrow(lib_df), length(balanced_dfs)
  ))
}

# Build basic library
if (length(vegs) > 0) {
  if (requireNamespace("progressr", quietly = TRUE)) {
    progressr::with_progress({
      p_lib <- progressr::progressor(steps = length(vegs), message = "Building vegetation library")

      lib_df$doy_for_lib <- lubridate::yday(lib_df$date)

      for (v in vegs) {
        dveg <- lib_df[lib_df$Veg == v & is.finite(lib_df$doy_for_lib), , drop = FALSE]
        if (nrow(dveg) == 0) {
          next
        }

        lib[[v]] <- list()
        p_lib()

        kept_idx <- intersect(avail, names(dveg))
        for (idx in kept_idx) {
          # Use medoid per DOY (closest-to-mean observed value) to preserve real signatures
          vals_by_doy <- tapply(seq_along(dveg[[idx]]), dveg$doy_for_lib, function(indices) {
            vals <- dveg[[idx]][indices]
            vals <- vals[is.finite(vals)]
            if (length(vals) == 0) return(NA_real_)
            if (length(vals) == 1) return(as.numeric(vals[1]))
            m <- mean(vals)
            as.numeric(vals[which.min(abs(vals - m))])
          })

          mu <- rep(NA_real_, 365)
          doy_values <- as.integer(names(vals_by_doy))
          valid_doy <- doy_values >= 1 & doy_values <= 365
          mu[doy_values[valid_doy]] <- vals_by_doy[valid_doy]

          if (all(!is.finite(mu))) next
          mv <- calc_moving_var(data.frame(date = 1:365, idx = mu), "idx", window = 14)
          lib[[v]][[idx]] <- list(mu = mu, mv = mv)
        }
        lib[[v]]$n_samples <- nrow(dveg)
        cat(sprintf("Processed prototype for %s (n=%d)\n", v, nrow(dveg)))
      }

      lib_df$doy_for_lib <- NULL
    })
  } else {
    stop("progressr package required for library construction")
  }
}

if (length(lib) == 0) stop("lib could not be constructed from training data")
cat(sprintf("Constructed simple lib from training dataset: vegetations=%d\n", length(lib)))

timing_info$lib_construction_done <- Sys.time()
cat(sprintf(
  "Library construction completed in %.1f seconds\n",
  as.numeric(difftime(timing_info$lib_construction_done, timing_info$moving_var_done, units = "secs"))
))

# Build PCA projections
if (length(avail) >= 2) {
  M_list <- list()
  for (v in names(lib)) {
    Mv <- matrix(NA_real_, nrow = 365, ncol = length(avail))
    colnames(Mv) <- avail
    for (j in seq_along(avail)) {
      idx <- avail[j]
      if (!is.null(lib[[v]][[idx]])) Mv[, j] <- lib[[v]][[idx]]$mu
    }

    for (j in seq_along(avail)) {
      idx <- avail[j]
      if (!is.null(lib[[v]][[idx]]) && !is.null(lib[[v]][[idx]]$mv)) {
        mv_col <- paste0(idx, "_mv")
        Mv <- cbind(Mv, lib[[v]][[idx]]$mv)
        colnames(Mv)[ncol(Mv)] <- mv_col
      }
    }

    days <- seq_len(365)
    for (i in seq_len(ncol(Mv))) {
      col <- as.numeric(Mv[, i])
      if (any(!is.finite(col))) {
        finite_idx <- which(is.finite(col))
        if (length(finite_idx) == 0) {
          col[] <- 0
        } else if (length(finite_idx) == 1) {
          col[] <- col[finite_idx]
        } else {
          x <- c(finite_idx - 365, finite_idx, finite_idx + 365)
          y <- rep(col[finite_idx], 3)
          interp <- tryCatch(
            {
              stats::approx(x = x, y = y, xout = days, rule = 2)$y
            },
            error = function(e) stop(sprintf("interpolation failed while building vegetation matrix: %s", e$message))
          )
          interp <- as.numeric(interp)
          if (!all(is.finite(interp))) interp[!is.finite(interp)] <- median(col[finite_idx], na.rm = TRUE)
          col <- interp
        }
        Mv[, i] <- col
      }
    }

    M_list[[v]] <- Mv
  }

  # Shape-based normalization (optional)
  BAND_SCALE <- list()
  if (isTRUE(SHAPE_NORMALIZATION_ENABLE)) {
    for (vn in names(lib)) {
      BAND_SCALE[[vn]] <- list()
      for (idx in avail) {
        if (!is.null(lib[[vn]][[idx]]) && !is.null(lib[[vn]][[idx]]$mu) && length(lib[[vn]][[idx]]$mu) >= 80) {
          ts_data <- lib[[vn]][[idx]]$mu
          n_days <- length(ts_data)

          if (n_days >= 50) {
            window_medians_low <- sapply(1:(n_days - 49), function(i) {
              median(ts_data[i:(i + 49)], na.rm = TRUE)
            })
            lowest_50_median <- min(window_medians_low, na.rm = TRUE)
          } else {
            lowest_50_median <- median(ts_data, na.rm = TRUE)
          }

          if (n_days >= 30) {
            window_medians_high <- sapply(1:(n_days - 29), function(i) {
              median(ts_data[i:(i + 29)], na.rm = TRUE)
            })
            top_30_median <- max(window_medians_high, na.rm = TRUE)
          } else {
            top_30_median <- median(ts_data, na.rm = TRUE)
          }

          scale_value <- median(c(lowest_50_median, top_30_median), na.rm = TRUE)
          if (!is.finite(scale_value) || scale_value <= 1e-6) scale_value <- 1.0

          BAND_SCALE[[vn]][[idx]] <- scale_value
        } else {
          BAND_SCALE[[vn]][[idx]] <- 1.0
        }
      }
    }
  } else {
    # When disabled, use identity scaling (no shape normalization)
    for (vn in names(lib)) {
      BAND_SCALE[[vn]] <- list()
      for (idx in avail) {
        BAND_SCALE[[vn]][[idx]] <- 1.0
      }
    }
  }

  # Global shape projection PCA
# Global shape projection PCA
V_names <- names(M_list)
if (length(V_names) == 0) {
  stop("Shape projection requires at least one vegetation class with data")
}

K_idx <- length(avail)
if (K_idx < 1) stop("No indices available for shape projection")

X_blocks <- list()
expected_cols <- NULL

for (vn in names(M_list)) {
  Mv <- M_list[[vn]]
  if (!is.null(Mv) && is.matrix(Mv)) {
    expected_cols <- ncol(Mv)
    break
  }
}

if (is.null(expected_cols)) {
  stop("Shape projection requires valid vegetation matrices")
}

for (vn in names(M_list)) {
  Mv <- M_list[[vn]]
  if (is.null(Mv) || !is.matrix(Mv)) next
  if (ncol(Mv) != expected_cols) {
    stop(sprintf("Column count mismatch for vegetation '%s'", vn))
  }
  X_blocks[[length(X_blocks) + 1]] <- Mv
}

if (length(X_blocks) == 0) {
  stop("Shape projection requires per-vegetation matrices")
}

X_all <- chunked_rbind(X_blocks, chunk_size = 25L)
gc()

if (!is.matrix(X_all) || nrow(X_all) < 2) {
  stop("Failed to assemble stacked X_all matrix")
}
if (is.null(colnames(X_all))) {
  stop("X_all matrix has no column names")
}

avail_aug <- colnames(X_all)
mu_all <- colMeans(X_all, na.rm = TRUE)
Xc <- sweep(X_all, 2, mu_all, "-")

band_scale_matrix <- matrix(1, nrow = length(V_names), ncol = length(avail_aug))
rownames(band_scale_matrix) <- V_names
colnames(band_scale_matrix) <- avail_aug

n_days <- 365

for (v_idx in seq_along(V_names)) {
  vn <- V_names[v_idx]
  if (!vn %in% names(BAND_SCALE)) next
  for (idx in avail) {
    col_name <- idx
    if (col_name %in% colnames(band_scale_matrix)) {
      scale_val <- BAND_SCALE[[vn]][[idx]]
      if (is.finite(scale_val) && scale_val > 0) {
        band_scale_matrix[v_idx, col_name] <- scale_val
      }
    }
    col_name_mv <- paste0(idx, "_mv")
    if (col_name_mv %in% colnames(band_scale_matrix)) {
      scale_val <- BAND_SCALE[[vn]][[idx]]
      if (is.finite(scale_val) && scale_val > 0) {
        band_scale_matrix[v_idx, col_name_mv] <- scale_val
      }
    }
  }
}

Xs <- Xc
for (v_idx in seq_along(V_names)) {
  start_row <- (v_idx - 1) * n_days + 1
  end_row <- v_idx * n_days
  if (end_row <= nrow(Xs)) {
    veg_scales <- band_scale_matrix[v_idx, ]
    Xs[start_row:end_row, ] <- sweep(Xc[start_row:end_row, , drop = FALSE], 2, veg_scales, "/")
  }
}
Xs[!is.finite(Xs)] <- 0

feature_sds <- apply(Xs, 2, sd, na.rm = TRUE)
feature_sds[feature_sds <= 1e-10] <- 1.0

cat(sprintf(
  "Feature std devs before final standardization (first 5): %s\n",
  paste(round(feature_sds[1:min(5, length(feature_sds))], 3), collapse = ", ")
))

Xs_std <- sweep(Xs, 2, feature_sds, "/")
Xs_std[!is.finite(Xs_std)] <- 0

sv <- try(svd_mkl(Xs_std, nu = 0, nv = min(ncol(Xs_std), nrow(Xs_std))), silent = TRUE)
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

if (pca_rank == 0) stop("PCA selection resulted in 0 factors")

V_pca <- sv$v[, seq_len(pca_rank), drop = FALSE]
if (ncol(V_pca) > 1) {
  V_pca <- qr.Q(qr(V_pca))[, seq_len(ncol(V_pca)), drop = FALSE]
}

variance_explained <- sv$d^2 / sum(sv$d^2)
cumulative_variance <- cumsum(variance_explained)

cat(sprintf(
  "PCA (Kaiser criterion) trained on standardized data across %d vegetation classes; retained %d components (eigenvalues > 1, cum var: %.1f%%)\n",
  length(V_names), ncol(V_pca), 100 * cumulative_variance[pca_rank]
))
cat(sprintf(
  "  First %d eigenvalues: %s\n",
  min(5, length(eigenvalues)),
  paste(round(eigenvalues[1:min(5, length(eigenvalues))], 2), collapse = ", ")
))

global_avg_scales <- colMeans(band_scale_matrix, na.rm = TRUE)
global_avg_scales[!is.finite(global_avg_scales)] <- 1.0

cat(sprintf(
  "Global average scales (first 5): %s\n",
  paste(round(global_avg_scales[1:min(5, length(global_avg_scales))], 3), collapse = ", ")
))

GLOBAL_PCA <- list(
  global = list(
    V = if (is.matrix(V_pca)) V_pca else matrix(0, nrow = length(avail_aug), ncol = 0),
    col_means = mu_all,
    col_scales = band_scale_matrix,
    global_avg_scales = global_avg_scales,
    feature_sds = feature_sds,
    rank = if (is.matrix(V_pca)) ncol(V_pca) else 0
  ),
  V = if (is.matrix(V_pca)) V_pca else matrix(0, nrow = length(avail_aug), ncol = 0),
  col_means = mu_all,
  col_scales = band_scale_matrix,
  global_avg_scales = global_avg_scales,
  feature_sds = feature_sds,
  per_bin = NULL,
  rank = if (is.matrix(V_pca)) ncol(V_pca) else 0,
  idx_order = avail_aug
)

lib_factor_pca <- list()
gpca_order <- GLOBAL_PCA$idx_order
gpca_means <- GLOBAL_PCA$col_means
gpca_scales <- GLOBAL_PCA$col_scales
vpca <- GLOBAL_PCA$global$V
pca_rank <- GLOBAL_PCA$global$rank

for (vname in names(M_list)) {
  dveg <- lib_df[lib_df$Veg == vname & !is.na(lib_df$date), , drop = FALSE]
  n_samples_v <- if (!is.null(dveg)) nrow(dveg) else 0

  if (n_samples_v <= 0) {
    lib_factor_pca[[vname]] <- list(T = matrix(0, nrow = 365, ncol = max(0, pca_rank)), n_samples = 0)
    next
  }

  proj_list <- list()
  proj_doy <- integer(0)

  for (i_row in seq_len(nrow(dveg))) {
    row <- dveg[i_row, , drop = FALSE]
    doy_row <- as.integer(lubridate::yday(row$date))
    if (is.na(doy_row) || doy_row < 1 || doy_row > 365) next

    raw_vec <- rep(NA_real_, length(gpca_order))
    names(raw_vec) <- gpca_order

    for (idx in avail) {
      kpos <- match(idx, gpca_order)
      if (is.na(kpos)) next
      val <- row[[idx]]
      if (!is.finite(val)) val <- gpca_means[kpos]
      raw_vec[kpos] <- as.numeric(val)
    }

    for (idx in avail) {
      mv_name <- paste0(idx, "_mv")
      kpos_mv <- match(mv_name, gpca_order)
      if (is.na(kpos_mv)) next
      mv_col <- paste0(idx, "_var14")
      mv_val <- if (mv_col %in% names(dveg)) row[[mv_col]] else NA_real_
      if (!is.finite(mv_val)) mv_val <- gpca_means[kpos_mv]
      raw_vec[kpos_mv] <- as.numeric(mv_val)
    }

    nas <- which(!is.finite(raw_vec))
    if (length(nas) > 0) raw_vec[nas] <- gpca_means[nas]

    if (pca_rank > 0 && ncol(vpca) >= pca_rank) {
      centered <- (raw_vec - gpca_means) / gpca_scales[vname, ]
      centered[!is.finite(centered)] <- 0
      t_scores <- as.numeric(centered %*% vpca)
    } else {
      t_scores <- numeric(0)
    }

    proj_list[[length(proj_list) + 1]] <- t_scores
    proj_doy <- c(proj_doy, doy_row)
  }

  if (length(proj_list) > 0 && pca_rank > 0) {
    obs_mat <- chunked_rbind(proj_list, chunk_size = 50L)
    gc()
    T_pca <- matrix(NA_real_, nrow = 365, ncol = ncol(obs_mat))
    for (d in seq_len(365)) {
      idx <- which(proj_doy == d)
      if (length(idx) > 0) {
        sub <- obs_mat[idx, , drop = FALSE]
        if (is.vector(sub)) sub <- matrix(sub, nrow = 1)
        mi <- medoid_row_index(sub)
        T_pca[d, ] <- sub[mi, ]
      }
    }
    for (col_idx in seq_len(ncol(T_pca))) {
      finite_mask <- is.finite(T_pca[, col_idx])
      if (sum(finite_mask) >= 2) {
        T_pca[!finite_mask, col_idx] <- approx(
          x = which(finite_mask),
          y = T_pca[finite_mask, col_idx],
          xout = which(!finite_mask),
          rule = 2
        )$y
      } else if (sum(finite_mask) == 1) {
        T_pca[, col_idx] <- T_pca[finite_mask, col_idx][1]
      } else {
        T_pca[, col_idx] <- 0
      }
    }
    lib_factor_pca[[vname]] <- list(T = T_pca, n_samples = n_samples_v)
  } else {
    lib_factor_pca[[vname]] <- list(T = matrix(0, nrow = 365, ncol = pca_rank), n_samples = n_samples_v)
  }
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
  if (!is.null(gproj$global_avg_scales)) {
    col_scales <- gproj$global_avg_scales
  } else {
    # Fallback: compute average if not pre-computed
    if (is.matrix(gproj$col_scales)) {
      col_scales <- colMeans(gproj$col_scales, na.rm = TRUE)
    } else {
      col_scales <- rep(1, length(idx_order))
    }
  }
  
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
  
  # Apply global scaling + standardization
  centered <- (raw_vec - col_means) / col_scales / feature_sds
  centered[!is.finite(centered)] <- 0
  
  V <- gproj$global$V
  if (is.null(V) || ncol(V) == 0) return(NULL)
  as.numeric(centered %*% V)
}
# Build MESMA library with endmember variants
build_mesma_variants <- function(lib_df, lib_factor_pca, veg_types, min_cluster_size = 10) {
  mesma_lib <- list()

  for (veg in veg_types) {
    veg_data <- lib_df[tolower(lib_df$Veg) == tolower(veg), ]
    if (nrow(veg_data) < min_cluster_size) {
      # Single variant for small samples
      mesma_lib[[veg]] <- list(
        list(
          T_pca = lib_factor_pca[[veg]]$T,
          variant_id = paste0(veg, "_single"),
          n_samples = nrow(veg_data)
        )
      )
      next
    }

    # Project vegetation data to PCA space for clustering
    veg_projections <- list()
    for (i in seq_len(nrow(veg_data))) {
      date_data <- prepare_factor_data(
        veg_data[i, , drop = FALSE],
        GLOBAL_PCA,
        avail,
        veg
      )
      if (length(date_data) > 0) {
        veg_projections[[length(veg_projections) + 1]] <- date_data[[1]]$z
      }
    }

    if (length(veg_projections) < 2) {
      mesma_lib[[veg]] <- list(
        list(
          T_pca = lib_factor_pca[[veg]]$T,
          variant_id = paste0(veg, "_single"),
          n_samples = nrow(veg_data)
        )
      )
      next
    }

    # Cluster to find variants
  proj_matrix <- chunked_rbind(veg_projections, chunk_size = 25L)
  gc()
    max_k <- min(5, floor(nrow(proj_matrix) / 5))
    if (max_k < 2) max_k <- 2

    # Find optimal number of clusters using silhouette
    best_k <- 2
    best_sil <- -1
    for (k in 2:max_k) {
      km <- kmeans_mkl(proj_matrix, centers = k, nstart = 10)
      sil <- mean(silhouette(km$cluster, dist(proj_matrix))[, 3])
      if (sil > best_sil) {
        best_sil <- sil
        best_k <- k
      }
    }

    km_final <- kmeans_mkl(proj_matrix, centers = best_k, nstart = 25)

    # Create variants based on clusters
    variants <- list()
    if (best_k >= 1L) {
      for (clust in seq_len(best_k)) {
        clust_members <- which(km_final$cluster == clust)
        if (length(clust_members) < 3) next

        # Build variant-specific temporal signature using medoids per DOY
        variant_t_pca <- matrix(0, 365, ncol(lib_factor_pca[[veg]]$T))

        for (doy in 1:365) {
          # rows in this cluster and DOY
          doy_rows <- which(lubridate::yday(veg_data$date) == doy)
          doy_rows <- intersect(doy_rows, clust_members)
          if (length(doy_rows) == 0) next
          # Project all those rows to factor spaces and pick medoid row
          z_p_list <- list(); z_l_list <- list()
          for (rid in doy_rows) {
            row <- veg_data[rid, , drop = FALSE]
            zp <- project_row_to_factors(row, GLOBAL_PCA, avail, veg)
            if (!is.null(zp)) z_p_list[[length(z_p_list) + 1]] <- zp
          }
          if (length(z_p_list) > 0) {
            ZP <- do.call(rbind, z_p_list); if (is.vector(ZP)) ZP <- matrix(ZP, nrow = 1)
            mi <- medoid_row_index(ZP)
            variant_t_pca[doy, ] <- ZP[mi, ]
          }
        }
        # Interpolate missing values per column
        if (ncol(variant_t_pca) > 0) {
          for (col_idx in seq_len(ncol(variant_t_pca))) {
            finite_mask <- is.finite(variant_t_pca[, col_idx])
            if (sum(finite_mask) >= 2) {
              variant_t_pca[!finite_mask, col_idx] <- approx(
                x = which(finite_mask),
                y = variant_t_pca[finite_mask, col_idx],
                xout = which(!finite_mask),
                rule = 2
              )$y
            } else if (sum(finite_mask) == 1) {
              variant_t_pca[, col_idx] <- variant_t_pca[finite_mask, col_idx][1]
            }
          }
        }


        variants[[length(variants) + 1]] <- list(
          T_pca = variant_t_pca,
          variant_id = paste0(veg, "_v", clust),
          n_samples = length(clust_members),
          cluster_center = km_final$centers[clust, ]
        )
      }

      if (length(variants) == 0) {
        variants <- list(list(
          T_pca = lib_factor_pca[[veg]]$T,
          variant_id = paste0(veg, "_single"),
          n_samples = nrow(veg_data)
        ))
      }

      mesma_lib[[veg]] <- variants
    }

    # finished building variants for this veg type
  }
  # finished building variants for all veg types - return the assembled library
  return(mesma_lib)
}

  

  solve_weights_simplex <- function(E, y, ridge = 1e-6) {
    EtE <- crossprod_mkl(E)
    Ety <- as.numeric(crossprod_mkl(E, matrix(y, ncol = 1)))
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
    D <- crossprod_mkl(E)
    if (!is.matrix(D)) D <- as.matrix(D)
    # Ensure symmetry and positive definiteness
    D <- 0.5 * (D + t(D)) + (lambda + 1e-8) * diag(p)
    d <- as.numeric(crossprod_mkl(E, matrix(y, ncol = 1)))

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
    combos <- do.call(expand.grid, idx_lists)
    n_combos <- nrow(combos)
    chunk_size <- 100L
    best_rmse <- Inf
    best_result <- NULL
    for (start in seq(1, n_combos, by = chunk_size)) {
      end <- min(start + chunk_size - 1L, n_combos)
      chunk <- combos[start:end, , drop = FALSE]
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

compress_and_unmix_year <- function(dly_year, mesma_lib, budget = TEMPORAL_BUDGET, topK = TOPK_VARIANTS) {
  veg_types <- names(mesma_lib)
  if (length(veg_types) == 0) return(NULL)
  
  g_proj <- GLOBAL_PCA 
  if (is.null(g_proj)) return(NULL)
  
  veg_ref <- unique(dly_year$Veg)[1]
  
  cat(sprintf("[3a] Calling prepare_factor_data with veg_ref=%s\n", veg_ref))
  date_list <- tryCatch({
    prepare_factor_data(dly_year, g_proj, avail, veg_ref)
  }, error = function(e) {
    cat(sprintf("[3a ERROR] prepare_factor_data failed: %s\n", e$message))
    return(NULL)
  })
  
  if (length(date_list) == 0) {
    cat("[3a] prepare_factor_data returned empty list\n")
    return(NULL)
  }
  
  k <- ncol(g_proj$global$V)
  if (is.null(k) || k == 0) return(NULL)
  
  cat(sprintf("[3b] Building Z365 matrix with k=%d\n", k))
  Z <- tryCatch({
    build_Z365(date_list, k)
  }, error = function(e) {
    cat(sprintf("[3b ERROR] build_Z365 failed: %s\n", e$message))
    return(NULL)
  })
  
  if (is.null(Z)) return(NULL)
  
  if (isTRUE(ENABLE_PHASE_ALIGNMENT)) {
    cat("[3c] Applying phase alignment\n")
    Z <- align_to_phenological_phase(Z, reference_phase = REFERENCE_PHASE_MARKERS)
  }
  
  cat("[3d] Computing information content\n")
  info <- compute_information_content(Z)
  
  cat(sprintf("[3e] Creating adaptive grid with budget=%d\n", budget))
  obs_grid <- create_adaptive_grid(info, budget)
  
  cat(sprintf("[3f] Extracting features (ENABLE_MULTISCALE=%s)\n", ENABLE_MULTISCALE))
  y <- if (isTRUE(ENABLE_MULTISCALE)) {
    extract_multiscale_features(Z, obs_grid, info, windows = MULTISCALE_WINDOWS)
  } else {
    extract_grid_features(Z, obs_grid, info)
  }

  # Map observation grid to template grid type
  gl <- length(obs_grid)
  grid_type <- if (gl <= ceiling(budget * 0.85)) "sparse" else if (gl >= ceiling(budget * 1.15)) "dense" else "medium"
  cat(sprintf("[3g] Grid type selected: %s (grid length=%d)\n", grid_type, gl))

  # Access compressed templates
  if (!exists(".COMPRESSED_TEMPLATES_ACCESSOR", envir = globalenv())) {
    cat("[3h ERROR] .COMPRESSED_TEMPLATES_ACCESSOR not found\n")
    return(NULL)
  }
  compressed_templates <- get(".COMPRESSED_TEMPLATES_ACCESSOR", envir = globalenv())
  
  cat("[3i] Getting variant templates\n")
  comp_templates <- tryCatch({
    get_variant_templates(veg_types, grid_type, compressed_templates, mesma_lib)
  }, error = function(e) {
    cat(sprintf("[3i ERROR] get_variant_templates failed: %s\n", e$message))
    return(NULL)
  })
  
  if (is.null(comp_templates)) return(NULL)
  
  allowed_veg <- veg_types
  
  cat("[3j] Computing top variants\n")
  top_variants <- list()
  for (v in veg_types) {
    if (!(v %in% allowed_veg)) next
    cand <- comp_templates[[v]]
    if (length(cand) == 0) next
    sims <- sapply(cand, function(x) cos_sim(y, x$vec))
    ord <- order(sims, decreasing = TRUE)
    keep <- ord[seq_len(min(topK, length(ord)))]
    top_variants[[v]] <- cand[keep]
  }
  
  cat("[3k] Building pooled matrix for ridge estimation\n")
  pool_cols <- list()
  for (v in names(top_variants)) {
    for (i in seq_along(top_variants[[v]])) pool_cols[[length(pool_cols) + 1]] <- top_variants[[v]][[i]]$vec
  }
  
  ny <- sqrt(sum(y^2))
  y_norm <- if (ny > 0) y / ny else y
  lambda_star <- 1e-6
  
  if (length(pool_cols) > 0) {
    cat(sprintf("[3l] Estimating optimal ridge from %d pooled columns\n", length(pool_cols)))
    E_pool <- do.call(cbind, lapply(pool_cols, function(col) {
      col <- as.numeric(col)
      nn <- sqrt(sum(col^2))
      if (nn > 0) col / nn else col
    }))
    lambda_star <- estimate_optimal_ridge(E_pool, y_norm)
    if (!is.finite(lambda_star) || lambda_star <= 0) lambda_star <- 1e-6
    cat(sprintf("Selected lambda=%.2e\n", lambda_star))
  }

  cat("[3m] Evaluating all combinations\n")
  best <- list(rmse = Inf, w = NULL, chosen = NULL)
  res <- tryCatch({
    evaluate_all_combinations(y_norm, top_variants, lambda = lambda_star, early_stop_rmse = EARLY_STOP_RMSE_THRESHOLD)
  }, error = function(e) {
    cat(sprintf("[3m ERROR] evaluate_all_combinations failed: %s\n", e$message))
    return(NULL)
  })
  
  if (!is.null(res)) best <- list(rmse = res$rmse, w = res$w, chosen = res$ids)
  
  diagnostics <- NULL
  uncertainty <- NULL
  
  if (isTRUE(ENABLE_DIAGNOSTICS) && !is.null(best$w) && !is.null(best$chosen)) {
    cat("[3n] Computing diagnostics\n")
    # ... rest of diagnostics code
  }
  
  cat("[3o] Returning results\n")
  list(vegetation_proportions = best$w, chosen_variants = best$chosen, rmse = best$rmse, diagnostics = diagnostics, uncertainty = uncertainty)
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
      if (!exists("min_idx_presence")) min_idx_presence <- 0.5
      min_required <- ceiling(length(avail_idx) * min_idx_presence)
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
        si <- if (is.null(gpca$col_scales)) {
          1
        } else if (is.list(gpca$col_scales)) {
          if (veg_type %in% names(gpca$col_scales) && kpos <= length(gpca$col_scales[[veg_type]])) {
            gpca$col_scales[[veg_type]][kpos]
          } else {
            1
          }
        } else if (is.matrix(gpca$col_scales)) {
          if (veg_type %in% rownames(gpca$col_scales) && kpos <= ncol(gpca$col_scales)) {
            gpca$col_scales[veg_type, kpos]
          } else {
            1
          }
        } else {
          1
        }
        if (!is.finite(si) || si <= 0) si <- 1

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
        si_mv <- if (is.list(gpca$col_scales)) {
          gpca$col_scales[[veg_type]][kpos_mv]
        } else {
          gpca$col_scales[veg_type, kpos_mv]
        }
        if (!is.finite(si_mv) || si_mv <= 0) si_mv <- 1

        vals_aug <- c(vals_aug, as.numeric((yy_mv - mi_mv) / si_mv))
        vrows_list[[length(vrows_list) + 1]] <- gpca$V[kpos_mv, , drop = FALSE]
      }
      if (is.null(vals_aug)) next

      if (length(vrows_list) == 0) next
      vsub <- do.call(rbind, vrows_list)
      if (ncol(vsub) < 1) next

      xtx <- crossprod_mkl(vsub)
      xty <- as.numeric(crossprod_mkl(vsub, matrix(vals_aug, ncol = 1)))

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

  rcompute_diagnostics <- function(y, E, w, mesma_result = NULL) {
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
      if (length(s) > 0 && min(s) > 0) max(s) / min(s) else stop("rcompute_diagnostics: singular values invalid for condition number calculation")
    }, error = function(e) stop(sprintf("rcompute_diagnostics: failed computing SVD for condition number: %s", e$message)))
    
    data.frame(
      condition_number = cond_num,
      residual_sum_of_squares = rss,
      r_squared = r_squared,
      stringsAsFactors = FALSE
    )
  }

  # Build MESMA library
  cat("Building MESMA endmember library...\n")
  mesma_lib <- build_mesma_variants(lib_df, lib_factor_pca, names(lib), min_cluster_size = 10)

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
      BAND_SCALE = if (exists("BAND_SCALE")) BAND_SCALE else NULL
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

  # Prepare test tasks
  cat("Constructing task list from TEST dataset...\n")

  df_tasks <- if (!is.null(TEST_YEARS)) {
    df_full[df_full$year %in% TEST_YEARS, , drop = FALSE]
  } else {
    df_full
  }

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
      df_tasks <- dplyr::left_join(df_tasks, gpts_map, by = "location_id")
    } else {
      # perform a straight join attempt (this will error if keys missing)
      df_tasks <- dplyr::left_join(df_tasks, gpts_map, by = "location_id")
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

    dbg_return_null <- function(reason) {
      # Normalize reason string
      reason_str <- as.character(reason)
      if (!grepl("^(ERROR:|error:)", reason_str, ignore.case = TRUE)) {
        reason_str <- paste0("ERROR:", reason_str)
      }

      # Prepare output CSV row
      row <- data.frame(
        location_id = as.character(loc),
        year = as.integer(yr),
        reason = reason_str,
        ts = as.character(Sys.time()),
        stringsAsFactors = FALSE
      )

      outf <- NULL
      success <- FALSE
      try({
        if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
        outf <<- file.path(OUT_DIR, "fit_fail_reasons.csv")
        if (!file.exists(outf)) {
          write.table(row, outf, sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
        } else {
          write.table(row, outf, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        }
        # Basic check: file exists and has non-zero size
        success <<- file.exists(outf) && file.info(outf)$size > 0
      }, silent = TRUE)

      return(NULL)
    }

    res_safe <- tryCatch(
      {
        dly <- df_tasks[df_tasks$location_id == loc, , drop = FALSE]
        dly_year <- dly[lubridate::year(dly$date) == yr, , drop = FALSE]
        if (nrow(dly_year) == 0) {
          return(dbg_return_null("no_rows_year"))
        }
        if (nrow(dly) == 0) {
          return(dbg_return_null("no_rows"))
        }

        # Calculate Q10 and Q90 DVI
        dly_train <- df_tasks[df_tasks$location_id == loc & lubridate::year(df_tasks$date) == yr, , drop = FALSE]
        dvi_vals <- dly_train$DVI[is.finite(dly_train$DVI)]
        q10_dvi <- if (length(dvi_vals) > 0) {
          stats::quantile(dvi_vals, 0.10, na.rm = TRUE)
        } else {
          NA
        }
        q90_dvi <- if (length(dvi_vals) > 0) {
          stats::quantile(dvi_vals, 0.90, na.rm = TRUE)
        } else {
          NA
        }

        if (!factor_mode) {
          return(dbg_return_null("not_factor_mode"))
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

        coef_df <- NULL
        variant_info_pca <- NULL

        mesma_result <- compress_and_unmix_year(dly_year, mesma_lib[veg_kept], budget = TEMPORAL_BUDGET, topK = TOPK_VARIANTS)
        if (is.null(mesma_result)) return(dbg_return_null("mesma_pca_failed"))
        coef <- mesma_result$vegetation_proportions
        coef_df <- data.frame(
          location_id = loc,
          year = yr,
          Veg = names(coef),
          coef = as.numeric(coef),
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
              coef_df$coef_lo[i] <- row_ci$coef_lo[1]
              coef_df$coef_hi[i] <- row_ci$coef_hi[1]
            }
          }
        }
        variant_info_pca <- if (!is.null(mesma_result$chosen_variants)) {
          # Create a one-row data frame with columns like '<veg>_variant'
          vi <- as.list(mesma_result$chosen_variants)
          names(vi) <- paste0(names(vi), "_variant")
          data.frame(location_id = loc, year = yr, vi, stringsAsFactors = FALSE, check.names = FALSE)
        } else NULL
        diag_df <- if (!is.null(mesma_result$diagnostics)) {
          dd <- as.data.frame(mesma_result$diagnostics, stringsAsFactors = FALSE)
          dd$location_id <- loc; dd$year <- yr
          dd <- dd[, c("location_id", "year", setdiff(names(dd), c("location_id", "year"))) , drop = FALSE]
          dd
        } else NULL
        unc <- mesma_result$uncertainty
      

        return(list(
          coef_df = coef_df,
          variant_trajectory = variant_info_pca,
          diagnostics = diag_df,
          uncertainty = unc,
          q10_dvi = q10_dvi,
          q90_dvi = q90_dvi
        ))
      },
      error = function(e) {
        dbg_return_null(paste0("error:", as.character(e$message)))
      }
    )
    dbg_return_null <- function(reason) {
    cat(sprintf("[FAIL] %s (year=%s): %s\n", loc, yr, reason), file = stderr())
    suppressMessages(suppressWarnings(tryCatch({
      # ... CSV write code
    }, error = function(e) invisible(NULL))))
    return(NULL)
  }
  
  res_safe <- tryCatch({
    cat(sprintf("[1] Getting data for %s/%s\n", loc, yr))
    dly <- df_tasks[df_tasks$location_id == loc, , drop = FALSE]
    cat(sprintf("    Found %d total rows\n", nrow(dly)))
    
    dly_year <- dly[lubridate::year(dly$date) == yr, , drop = FALSE]
    cat(sprintf("    Found %d rows for year %s\n", nrow(dly_year), yr))
    
    if (nrow(dly_year) == 0) return(dbg_return_null("no_rows_year"))
    
    cat(sprintf("[2] Calling compress_and_unmix_year\n"))
    mesma_result <- compress_and_unmix_year(dly_year, mesma_lib[veg_kept], 
                                            budget = TEMPORAL_BUDGET, topK = TOPK_VARIANTS)
    
    if (is.null(mesma_result)) {
      cat(sprintf("[3] compress_and_unmix_year returned NULL\n"))
      return(dbg_return_null("mesma_pca_failed"))
    }
    
    cat(sprintf("[4] Success for %s/%s\n", loc, yr))
    # ... rest of processing
  }, error = function(e) {
    cat(sprintf("[ERROR] %s/%s: %s\n", loc, yr, e$message))
    dbg_return_null(paste0("error:", e$message))
  })

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
    compress_and_unmix_year = compress_and_unmix_year,
    project_to_simplex = project_to_simplex
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
    cat("ERROR: All tasks returned NULL results! Check fit_fail_reasons.csv for details.\n")
    cat("Most likely causes:\n")
    cat("1. TEST_YEARS not set or empty\n")
    cat("2. No data available for the specified test years\n")
    cat("3. Data filtering issues\n")
    cat("4. All locations have insufficient data for fitting\n")
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
    all_diagnostics <- if (length(diag_list) > 0) tryCatch({ do.call(rbind, diag_list) }, error = function(e) stop(sprintf("Failed to combine diagnostics: %s", e$message))) else NULL

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
      # Projection method column removed per requirement; single method per run
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

}

