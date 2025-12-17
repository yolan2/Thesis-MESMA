 
filter_variants_by_min_samples <- function(variants, min_samples = MIN_ENDMEMBER_SAMPLES, veg = NULL, raw_template = NULL) {
  if (is.null(variants) || length(variants) == 0) return(list())
  keep_mask <- sapply(variants, function(v) {
    if (is.null(v$n_samples)) return(TRUE) # keep if no count provided
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
  if (length(lon) == 1 && length(lat) == 1) {
    if (!is.finite(lon) || !is.finite(lat)) return(NA_character_)
    sprintf("L_%0.6f_%0.6f", round(lon, 6), round(lat, 6))
  } else {
    res <- mapply(function(x, y) {
      if (!is.finite(x) || !is.finite(y)) return(NA_character_)
      sprintf("L_%0.6f_%0.6f", round(as.numeric(x), 6), round(as.numeric(y), 6))
    }, lon, lat, USE.NAMES = FALSE)
    as.character(res)
  }
}

 
library(zoo)
library(dplyr)

# Ensure legacy pipe `%>%` is available for any remaining files that haven't yet migrated to `|>`
if (!exists("%>%", mode = "function")) {
  if (requireNamespace("magrittr", quietly = TRUE)) {
    try(suppressPackageStartupMessages(library(magrittr)), silent = TRUE)
    if (exists("%>%", mode = "function")) cat("[NOTICE] Loaded magrittr to provide legacy `%>%` operator\n")
  } else {
    warning("magrittr not installed and legacy `%>%` operator not found. Consider migrating all pipelines to `|>` or installing 'magrittr'.")
  }
}

library(cluster)
if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}
library(readr)
if (!requireNamespace("sf", quietly = TRUE)) {
  install.packages("sf")
}
library(sf)
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)
if (!requireNamespace("scales", quietly = TRUE)) {
  install.packages("scales")
}
library(scales)
if (!requireNamespace("nlme", quietly = TRUE)) {
  install.packages("nlme")
}
library(nlme)
if (!requireNamespace("RStoolbox", quietly = TRUE)) {
  install.packages("RStoolbox")
}
library(RStoolbox)
if (!requireNamespace("nnls", quietly = TRUE)) {
  install.packages("nnls")
}
library(nnls)
if (!requireNamespace("terra", quietly = TRUE)) {
  install.packages("terra")
}
library(terra)
options(warn = 1)  # print warnings as they occur for debugging
options(future.progress = FALSE)  # disable progress bars from future package

# DISABLE global error handler - let tryCatch in main block handle errors
options(error = NULL)


# Ensure core initialization helpers are available (defines setup_parallel_backend(), constants, etc.)
if (!exists("setup_parallel_backend", mode = "function")) {
  if (file.exists("init_mesma.R")) {
    try(source("init_mesma.R"), silent = TRUE)
  }
}
if (!exists("setup_parallel_backend", mode = "function")) {
  warning("setup_parallel_backend() not found after attempting to source 'init_mesma.R'; defining a harmless fallback.")
  setup_parallel_backend <- function() {
    function() {}
  }
}


SKIP_PPI <- TRUE

INPUT_CSV <- "C:\\Users\\yolan\\Downloads\\landsat_timeseries_vegetation_filtered (4)_fixed.csv"
INFERENCE_CSV <- "C:\\Users\\yolan\\OneDrive\\Documenten\\UGENT\\Master\\masterproef\\GIS\\landsat_lower_inference.csv"

# AGENT: Expanded training years to provide more data for the weighting model,
# which was failing due to an insufficient number of unique time-series traces.
TRAIN_YEARS <- c(2022, 2023, 2024)  # Years to use for training (note: script uses all data)

# Check for override from environment or pre-set variable
if (!exists("PAENABLERALLEL_")) {
  PARALLEL_ENABLE <- if (nzchar(Sys.getenv("MESMA_PARALLEL"))) {
    as.logical(Sys.getenv("MESMA_PARALLEL"))
  } else {
    # CHANGE THIS TO FALSE IF YOU EXPERIENCE CRASHES:
    TRUE  # Default: enabled (set to FALSE to disable parallel processing)
  }
}

if (isTRUE(PARALLEL_ENABLE)) {
  cat("[CONFIG] Parallel processing ENABLED\n")
} else {
  cat("[CONFIG] Parallel processing DISABLED (running sequentially - slower but more stable)\n")
}

PARALLEL_WORKERS <- min(4, parallel::detectCores() - 1)  # Limit to max 4 workers for stability
COMBO_PARALLEL_WORKERS <- max(1L, floor(PARALLEL_WORKERS/2))

ENABLE_DIAGNOSTICS <- TRUE
TOP_K_CANDIDATES <- 10L
MAX_VARIANTS_PER_VEG <- 10
MIN_VARIANTS_PER_VEG <- 1
MAX_COMBOS_FOR_FULL_SEARCH <- 5000
STAGE1_WEIGHTING_ENABLED <- TRUE

MIN_ENDMEMBER_SAMPLES <- 5L

ALLOWED_VEG <- c("populus", "tamarix", "phragmites")

OPTIMAL_INDICES <- c("DVI", "GVI", "NIRv", "MSAVI2", "PSRI", "TCW", "PPI")

SKIP_INFERENCE <- FALSE

TESTING_MODE <- FALSE  # Enable debug messages to diagnose Stage 2 issue
TESTING_MAX_PER_VEG <- 5000L  # Upper cap per vegetation class when TESTING_MODE is enabled
if (!exists('DEBUG_STAGE1_VERBOSE')) DEBUG_STAGE1_VERBOSE <- FALSE
if (!exists('DEBUG')) DEBUG <- FALSE

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

# Removed legacy helper: fit_cost_mkl() — no longer used in the codebase.

if (file.exists("ppi_helpers.R")) {
  source("ppi_helpers.R")
} else {
  stop("Required file 'ppi_helpers.R' not found in project root. Please add it to ensure consistent PPI calculation.")
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

linearize_indices <- function(df) {
  cat("Applying linearization transformations to indices...\n")
  
  if ("NIRv" %in% names(df)) {
    cat("  Linearizing NIRv (2x - x^2)...\n")
    df$NIRv <- 2 * df$NIRv - df$NIRv^2
  }
  
  if ("MSAVI2" %in% names(df)) {
    cat("  Linearizing MSAVI2 (log(x+1))...\n")
    df$MSAVI2 <- log(df$MSAVI2 + 1)
  }
  
  if ("PSRI" %in% names(df)) {
    cat("  Linearizing PSRI (signed sqrt)...\n")
    df$PSRI <- sign(df$PSRI) * sqrt(abs(df$PSRI))
  }
  
  if ("TCW" %in% names(df)) {
    cat("  Linearizing TCW (log(x+1))...\n")
    df$TCW <- log(pmax(df$TCW + 1, 1e-6))
  }
  
  df
}

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

compute_indices_from_bands <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  eps <- 1e-9
  has_bands <- intersect(RAW_BANDS, names(df))
  if (length(has_bands) == 0) return(df)

  if (all(c('nir','red') %in% names(df))) df$DVI <- as.numeric(df$nir) - as.numeric(df$red)
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

  if (all(c('nir','red') %in% names(df))) df$NDVI <- (as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + as.numeric(df$red) + eps)
  if (all(c('nir','red') %in% names(df))) df$MSAVI2 <- (2 * as.numeric(df$nir) + 1 - sqrt(pmax(0, (2 * as.numeric(df$nir) + 1)^2 - 8 * (as.numeric(df$nir) - as.numeric(df$red))))) / 2
  if (all(c('nir','swir1') %in% names(df))) df$NDMI <- (as.numeric(df$nir) - as.numeric(df$swir1)) / (as.numeric(df$nir) + as.numeric(df$swir1) + eps)

  # NDTI - Normalized Difference Tillage Index (Van Deventer et al. 1997)
  if (all(c('swir1','swir2') %in% names(df))) df$NDTI <- (as.numeric(df$swir1) - as.numeric(df$swir2)) / (as.numeric(df$swir1) + as.numeric(df$swir2) + eps)

  # NDSI - Normalized Difference Snow Index
  if (all(c('green','swir1') %in% names(df))) df$NDSI <- (as.numeric(df$green) - as.numeric(df$swir1)) / (as.numeric(df$green) + as.numeric(df$swir1) + eps)

  # NDDI - Normalized Difference Dust Index
  if (all(c('red','nir') %in% names(df))) df$NDDI <- (as.numeric(df$red) - as.numeric(df$nir)) / (as.numeric(df$red) + as.numeric(df$nir) + eps)

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

MIN_SKIP_DOYS_PER_LOCATION <- 2L

FAST_VAR <- TRUE

TEMPORAL_AGGREGATION_DAYS <- 5L  # 5-day intervals (pentads)
N_TEMPORAL_BINS <- ceiling(365 / TEMPORAL_AGGREGATION_DAYS)  # = 73 pentads
TEMPORAL_BUDGET <- N_TEMPORAL_BINS  # Use pentad resolution (no compression beyond pentad aggregation)

STAGE1_WEIGHT_MAX_RATIO <- 6.0
N_STAGE1_VEG_ENDMEMBERS <- 3


TOPK_VARIANTS <- 5L
N_VARIANTS_PER_VEG <- 10L
MIN_CLUSTER_SIZE <- 10L
PCA_VARIANCE_THRESHOLD <- 0.95
LDA_WEIGHT_FLOOR <- 0.01
STAGE1_INDEX_WEIGHT_THRESHOLD <- 0.02
STAGE1_MIN_VEG_FRACTION <- 0.10
STAGE1_RESIDUAL_NORM_THRESHOLD <- 0.2  # fraction of y norm; if residual is larger, consider OLS fallback
SKIP_MOVING_VARIANCE <- TRUE
ENABLE_MULTISCALE <- FALSE
MULTISCALE_WINDOWS <- c(7L, 14L, 30L)
ENABLE_QP_SOLVER <- TRUE
COMBO_PARALLEL_ENABLE <- TRUE
EARLY_STOP_RMSE_THRESHOLD <- 0.0
ENABLE_DIAGNOSTICS <- TRUE














project_to_simplex <- function(v) {
  n <- length(v)
  if (n == 0) return(numeric(0))
  if (n == 1) return(1)

  if (any(!is.finite(v))) {
    warning("[project_to_simplex] Non-finite values in input, replacing with 0")
    v[!is.finite(v)] <- 0
  }

  u <- sort(v, decreasing = TRUE)

  cssv <- cumsum(u)
  cond_vec <- u + (1 - cssv) / seq_along(u) > 0

  if (!any(cond_vec)) {
    warning("[project_to_simplex] No valid rho found, returning equal weights. This often indicates that unconstrained fractions are all very negative or non-finite — please check endmember normalization and template quality (see FIX_SIMPLEX_PROJECTION.md).")
    return(rep(1 / n, n))
  }

  rho <- max(which(cond_vec))

  theta <- (cssv[rho] - 1) / rho

  w <- pmax(v - theta, 0)

  w_sum <- sum(w)
  if (w_sum > 0) {
    w <- w / w_sum
  } else {
    warning("[project_to_simplex] Projection resulted in zero vector; returning equal weights. Check normalization/template matching (see FIX_SIMPLEX_PROJECTION.md)")
    w <- rep(1 / n, n)
  }

  return(w)
}


ols_unmix_n_endmembers <- function(y, M, weights = NULL) {
  y <- as.numeric(y)
  N <- ncol(M)  # Number of endmembers
  P <- nrow(M)  # Number of bands/features
  
  # When weights are provided, compute a weighted residual
  effective_residual <- function(y_true, y_pred) {
    if (!is.null(weights) && length(weights) == P) {
      sqrt(sum(weights * (y_true - y_pred)^2, na.rm = TRUE))
    } else {
      sqrt(sum((y_true - y_pred)^2, na.rm = TRUE))
    }
  }

  if (N == 1) {
    y_proj <- M[, 1]
    residual <- effective_residual(y, y_proj)
    return(list(f = 1, residual = residual, y_proj = y_proj))
  }

  M_fit <- M
  y_fit <- y

  if (!is.null(weights) && length(weights) == P) {
    sqrt_w <- sqrt(pmax(weights, 0)) # Ensure non-negative weights
    M_fit <- M_fit * sqrt_w
    y_fit <- y_fit * sqrt_w
  }

  res <- nnls::nnls(M_fit, y_fit)
  f <- res$x
  
  # sum-to-one constraint post-NNLS
  if (sum(f) > 0) {
    f <- f / sum(f)
  }
  
  y_proj <- as.numeric(M %*% f)
  residual <- effective_residual(y, y_proj)
  
  return(list(f = f, residual = residual, y_proj = y_proj))
}


compute_lda_weights <- function(mesma_lib, compressed_templates_accessor, grid_type = "full") {
  cat("\n=== COMPUTING LDA FEATURE WEIGHTS ===\n")
  
  X_list <- list()
  y_list <- list()
  
  for (veg in names(mesma_lib)) {
    variants <- mesma_lib[[veg]]
    for (variant in variants) {
      vid <- variant$variant_id
      vec <- tryCatch(compressed_templates_accessor[[veg]][[vid]][[grid_type]], error = function(e) NULL)
      
      if (!is.null(vec)) {
        X_list[[length(X_list) + 1]] <- vec
        y_list[[length(y_list) + 1]] <- veg
      }
    }
  }
  
  if (length(X_list) == 0) return(NULL)
  
  X <- do.call(rbind, X_list)
  y <- unlist(y_list)
  
  counts <- table(y)
  if (any(counts < 2) || length(unique(y)) < 2) {
    cat("LDA skipped: insufficient classes or samples per class.\n")
    return(NULL)
  }
  
  # Compute per-column standard deviation (ignore NA). If all values are NA,
  # sd(..., na.rm=TRUE) returns NA and we set that to 0 (will be dropped).
  sds <- apply(X, 2, function(col) sd(col, na.rm = TRUE))
  sds[is.na(sds)] <- 0
  keep_cols <- sds > 1e-6
  # Diagnostic for low-variance features (helpful for PPI debugging)
  if ("PPI" %in% colnames(X)) {
    ppi_col <- X[, "PPI"]
    ppi_count <- sum(is.finite(ppi_col))
    ppi_sd <- sd(ppi_col, na.rm = TRUE)
    ppi_min <- min(ppi_col, na.rm = TRUE)
    ppi_max <- max(ppi_col, na.rm = TRUE)
    cat(sprintf("[PPI DIAG] PPI finite count=%d, sd=%.6g, range=[%.6g, %.6g]\n", ppi_count, ppi_sd, ppi_min, ppi_max))
  }
  if (sum(keep_cols) < 2) {
    cat("LDA skipped: too few variable features.\n")
    return(NULL)
  }
  
  X_clean <- X[, keep_cols, drop = FALSE]
  # Impute remaining NA values with column means (na.rm=TRUE) to allow PCA
  # and LDA to operate; columns with all-NA were removed above.
  if (any(!is.finite(as.matrix(X_clean)))) {
    col_means <- apply(X_clean, 2, function(col) mean(col, na.rm = TRUE))
    for (j in seq_len(ncol(X_clean))) {
      missing_idx <- !is.finite(X_clean[, j])
      if (any(missing_idx)) X_clean[missing_idx, j] <- col_means[j]
    }
  }
  
  pca_res <- tryCatch({
    prcomp(X_clean, center = TRUE, scale. = TRUE)
  }, error = function(e) {
    cat(sprintf("PCA failed: %s\n", e$message))
    NULL
  })
  
  if (is.null(pca_res)) return(NULL)
  
  cum_var <- cumsum(pca_res$sdev^2) / sum(pca_res$sdev^2)
  n_pcs <- which(cum_var >= 0.99)[1]
  if (is.na(n_pcs)) n_pcs <- length(cum_var)
  
  n_pcs <- max(min(n_pcs, ncol(X_clean)), 2)
  if (n_pcs > ncol(X_clean)) n_pcs <- ncol(X_clean)
  
  n_total <- nrow(X_clean)
  n_classes <- length(unique(y))
  max_pcs_allowed <- max(1, n_total - n_classes)
  if (n_pcs > max_pcs_allowed) {
    cat(sprintf("LDA pre-filter: reducing n_pcs from %d to %d to ensure p <= N - K (N=%d,K=%d)\n", n_pcs, max_pcs_allowed, n_total, n_classes))
    n_pcs <- max_pcs_allowed
  }

  if (n_pcs < 1) {
    cat("LDA pre-filter: not enough effective samples to compute LDA after PCA (n_pcs < 1). Skipping LDA weights.\n")
    return(NULL)
  }

  X_pca <- pca_res$x[, 1:n_pcs, drop = FALSE]
  cat(sprintf("LDA Pre-processing: Reduced %d features to %d PCs (%.1f%% variance)\n", 
              ncol(X_clean), n_pcs, cum_var[n_pcs]*100))
  
  rank_pca <- qr(X_pca)$rank
  n_total_after <- nrow(X_pca)
  n_classes_after <- length(unique(y))
  cat(sprintf("Pre-LDA diagnostic: N=%d, K=%d, n_pcs=%d, rank(X_pca)=%d\n", n_total_after, n_classes_after, ncol(X_pca), rank_pca))
  if (rank_pca < ncol(X_pca)) {
    cat(sprintf("X_pca rank-deficient: reducing PCs from %d to %d (rank) before LDA\n", ncol(X_pca), rank_pca))
    n_pcs <- rank_pca
    if (n_pcs < 1) {
      cat("After rank reduction, no PCs remain for LDA; skipping compute_lda_weights\n")
      return(NULL)
    }
    X_pca <- X_pca[, 1:n_pcs, drop = FALSE]
  }
  max_pcs_allowed_after <- max(1, n_total_after - n_classes_after)
  if (ncol(X_pca) > max_pcs_allowed_after) {
    cat(sprintf("Reducing PCs from %d to %d to satisfy p <= N - K (N=%d,K=%d)\n", ncol(X_pca), max_pcs_allowed_after, n_total_after, n_classes_after))
    n_pcs <- max_pcs_allowed_after
    X_pca <- X_pca[, 1:n_pcs, drop = FALSE]
  }
  min_n_pcs <- max(2, length(unique(y)) - 1)
  lda_res <- tryCatch({
    safe_lda_call(X_pca, y, min_n_pcs = min_n_pcs)
  }, error = function(e) {
    cat(sprintf("LDA failed: %s\n", e$message))
    NULL
  })
  
  if (is.null(lda_res)) return(NULL)
  
  
  W_pc <- lda_res$scaling
  R <- pca_res$rotation[, 1:n_pcs, drop = FALSE]
  
  W_std <- R %*% W_pc
  
  svd <- lda_res$svd
  prop <- svd / sum(svd)
  
  if (ncol(W_std) > 1) {
    n_dim <- min(length(prop), ncol(W_std))
    weights_clean <- rowSums(abs(W_std[, 1:n_dim, drop=FALSE]) %*% diag(prop[1:n_dim], nrow=n_dim))
  } else {
    weights_clean <- abs(W_std[, 1])
  }
  
  weights <- rep(0, ncol(X))
  weights[keep_cols] <- weights_clean

  # Diagnostic means before dampening
  weights_mean_before_threshold <- mean(weights, na.rm = TRUE)
  weights_mean_after_threshold <- mean(weights, na.rm = TRUE)

  # === AGENT FIX: DAMPEN WEIGHTS TO PREVENT VECTOR COLLAPSE ===
  # 1) Take sqrt(abs(W_LDA)) to dampen extreme ratios
  w_sqrt <- sqrt(pmax(abs(weights), 0))

  # 2) Apply a floor so that no feature is effectively zeroed (floor = 0.2 * max)
  floor_val <- 0.2 * max(w_sqrt, na.rm = TRUE)
  if (!is.finite(floor_val) || floor_val <= 0) floor_val <- 0
  w_floor <- pmax(w_sqrt, floor_val)

  # 3) Mean-normalize so mean(w) == 1.0 (preserves energy)
  avg_w <- mean(w_floor[w_floor > 0], na.rm = TRUE)
  if (is.finite(avg_w) && avg_w > 0) {
    w_final <- w_floor / avg_w
  } else {
    w_final <- w_floor
  }

  # 4) Hard cap to avoid dominance (safety)
  cap_val <- 3.0
  w_final[w_final > cap_val] <- cap_val

  # Attach diagnostics
  attr(w_final, "dampened") <- TRUE
  attr(w_final, "mean_after_dampening") <- mean(w_final, na.rm = TRUE)
  attr(w_final, "floor_value") <- floor_val
  attr(w_final, "cap_value") <- cap_val

  cat("LDA weights computed (dampened+normalized).\n")
  cat(sprintf("[DIAGNOSTIC] Weight means: before_threshold=%.6f, after_threshold=%.6f, after_dampening=%.6f, floor=%.6f, cap=%.6f\n", weights_mean_before_threshold, weights_mean_after_threshold, attr(w_final, "mean_after_dampening"), floor_val, cap_val))

  weights <- w_final
  return(weights)
}


enforce_stage1_weight_clip <- function(weights) {
  clip_factor <- if (exists("STAGE1_WEIGHT_MAX_RATIO", inherits = TRUE)) STAGE1_WEIGHT_MAX_RATIO else NA
  if (is.na(clip_factor) || !is.finite(clip_factor) || clip_factor <= 0) {
    return(list(weights = weights, clipped = 0, threshold = NA, old_max = max(weights, na.rm=TRUE), new_max = max(weights, na.rm=TRUE)))
  }

  clip_threshold <- clip_factor
  old_max <- max(weights, na.rm = TRUE)
  n_clipped <- sum(weights > clip_threshold, na.rm = TRUE)
  if (n_clipped > 0) {
    new_weights <- pmin(weights, clip_threshold)
    # No normalization - keep weights at raw scale after clipping
    return(list(weights = new_weights, clipped = n_clipped, threshold = clip_threshold, old_max = old_max, new_max = max(new_weights, na.rm = TRUE)))
  }

  list(weights = weights, clipped = 0, threshold = clip_threshold, old_max = old_max, new_max = old_max)
}




COMBO_SAFE_EXPAND_LIMIT <- 1e6    # fully expand grid up to this many combos
COMBO_ABORT_LIMIT <- 5e7          # abort if combos exceed this hard limit

BOOTSTRAP_B <- 200L  # Number of bootstrap iterations for location-based resampling
ENABLE_UNCERTAINTY <- TRUE
DEBUG_UNCERTAINTY <- FALSE  # global debug flag for uncertainty diagnostics
MIN_OBS_PER_LOC_YEAR <- 1L
MIN_UNIQUE_DOY_DEFAULT <- 5L
MIN_UNIQUE_DOY_INFERENCE <- 1L
VARIANT_SWITCH <- TRUE  # Always use variant switching in bootstraps
VARIANT_SWITCH_RE_CENTER <- TRUE


VARIANCE_THRESHOLD <- 0.90
MAX_VEG_COMPONENTS <- 8
GAM_K_MAX <- 40
GAM_GAMMA <- 1.0


USE_INDICES_MIN <- 1L
MIN_INDEX_SD <- 0.05

ENABLE_SAMPLE_BALANCING <- TRUE

ENABLE_SOIL_PREPROCESS <- TRUE
SOIL_PURE_THRESHOLD <- 0.95
SOIL_MIN_SAMPLES <- 3L

MAX_PROJECTIONS_PER_VEG <- 25000L  # subsample before clustering to avoid OOM
SILHOUETTE_SAMPLE_SIZE <- 20000L   # subsample for silhouette distance matrix
MEDOID_SAMPLE_SIZE <- 10000L       # subsample for medoid distance computation

EPS_SIGMA <- 1e-8
LOWER_BND <- 0

MIN_IDX_PRESENCE <- 0.5

if (!file.exists(INPUT_CSV)) {
  # Try plausible alternative filenames (previous versions)
  alt_candidates <- c(
    sub("\\(4\\)", "(3)", INPUT_CSV, fixed = TRUE),
    sub("\\(4\\)", "(2)", INPUT_CSV, fixed = TRUE),
    sub("landsat_timeseries_vegetation_filtered (4).csv", "landsat_timeseries_vegetation_filtered.csv", INPUT_CSV, fixed = TRUE)
  )
  alt_found <- alt_candidates[file.exists(alt_candidates)]
  if (length(alt_found) > 0) {
    INPUT_CSV <- alt_found[1]
    cat(sprintf("[NOTICE] INPUT_CSV not found; using fallback: %s\n", INPUT_CSV))
  } else {
    stop(paste0("Required input CSV not found: ", INPUT_CSV, " (and no alternatives found)"))
  }
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
cat(sprintf("After deduplication: %d rows remaining from original %d rows.\n", nrow(raw_df), nrow(raw_df)))

# === GEE INPUT MAPPING BLOCK ===
# Map common GEE-exported column names to the names expected by this script
# 1) Map 'vegetation' -> 'Veg'
if ("vegetation" %in% names(raw_df) && !"Veg" %in% names(raw_df)) {
  raw_df$Veg <- raw_df$vegetation
  cat("[NOTICE] Renamed 'vegetation' column to 'Veg'\n")
}
# 2) Disallow legacy 'no.soil' / 'no_soil' / 'no-soil' variants — require 'no soil' (space)
bad_nosoil_variants <- intersect(c("no.soil", "no_soil", "no-soil"), names(raw_df))
if (length(bad_nosoil_variants) > 0) {
  stop(sprintf("Unsupported column name(s) present: %s. Please rename to 'no soil' (with a space) if you intend to provide soil fraction metadata.", paste(bad_nosoil_variants, collapse = ", ")))
}
# Keep "no soil" as provided by the CSV (no auto-mapping)

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
if ("zenith.angle" %in% names(df)) df$zenith.angle <- NULL

df <- normalize_band_names(df)
# If bands are present but indices are missing, compute indices from bands
df <- compute_indices_from_bands(df)

if (!"date" %in% names(df) && "prediction_date" %in% names(df)) df$date <- as.Date(df$prediction_date)
if ("date" %in% names(df)) df$date <- as.Date(df$date)


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

normalize_mesma_data <- function(df, cols = OPTIMAL_INDICES, lat_default = 40.2) {
  cat("Applying comprehensive MESMA data normalization...\n")
  
  if (!"zenith.angle" %in% names(df)) df$zenith.angle <- NA_real_
  if (!"PPI" %in% names(df)) df$PPI <- NA_real_
  
  df <- linearize_indices(df)
  
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

apply_stored_normalization <- function(df, norm_params, cols = OPTIMAL_INDICES, lat_default = 40.2) {
  cat("Applying stored normalization parameters to data...\n")
  
  if (!"zenith.angle" %in% names(df)) df$zenith.angle <- NA_real_
  if (!"PPI" %in% names(df)) df$PPI <- NA_real_
  
  df <- linearize_indices(df)
  
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

apply_stage2_normalization <- function(vec, idx_names, norm_params) {
  if (is.null(norm_params$INDEX_SCALES) || length(norm_params$INDEX_SCALES) == 0) {
    return(vec)  # No normalization parameters available
  }
  
  if (length(vec) != length(idx_names)) {
    warning("apply_stage2_normalization: vec/idx_names length mismatch")
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

linearity_file <- "index_linearity_scores.csv"
if (file.exists(linearity_file)) {
  linearity_scores <- read.csv(linearity_file)
  good_indices <- linearity_scores$Index[linearity_scores$Normalized_Max_Dev <= 0.2]
  cat(sprintf("Filtering indices with linearity > 0.2 after linearization. Retained indices: %s\n", paste(good_indices, collapse = ", ")))
  
  # Ensure that Veg/'no soil' are preserved as part of metadata selection so
  # vegetation mapping survives the index filtering step
  meta_cols <- c("location_id", "lat", "lon", "imagery_lat", "imagery_lon", "target_lat", "target_lon", "date", "pheno_year", "doy", "prediction_date", "DUSTI", "PPI", "zenith.angle", "DVI_max", "Veg", "no soil")
  raw_bands_present <- intersect(RAW_BANDS, names(df))
  index_cols <- setdiff(names(df), c(meta_cols, raw_bands_present))
  keep_cols <- c(meta_cols, raw_bands_present, intersect(index_cols, good_indices))
  keep_cols <- intersect(keep_cols, names(df))  # Filter to only existing columns
  df <- df[, keep_cols, drop = FALSE]
  cat(sprintf("Preserved %d raw bands: %s\n", length(raw_bands_present), paste(raw_bands_present, collapse=", ")))
  
  OPTIMAL_INDICES <- good_indices
  if ("PPI" %in% names(df)) OPTIMAL_INDICES <- c(OPTIMAL_INDICES, "PPI")
} else {
  warning("Linearity scores file not found, proceeding without filtering")
}

missing_idx <- setdiff(OPTIMAL_INDICES, names(df))
if (length(missing_idx) > 0) {
  cat(sprintf("[NOTICE] Some OPTIMAL_INDICES are missing: %s. Attempting to compute from raw bands...\n", paste(missing_idx, collapse = ", ")))
  df <- compute_indices_from_bands(df)
  missing_idx <- setdiff(OPTIMAL_INDICES, names(df))
  if (length(missing_idx) > 0) {
    # Reduce OPTIMAL_INDICES to available indices and continue with a notice
    OPTIMAL_INDICES <- intersect(OPTIMAL_INDICES, names(df))
    if (length(OPTIMAL_INDICES) == 0) {
      stop(paste0("INPUT_CSV missing all OPTIMAL_INDICES and could not compute them from bands: ", paste(missing_idx, collapse = ", ")))
    } else {
      cat(sprintf("[NOTICE] Reduced OPTIMAL_INDICES to available set: %s\n", paste(OPTIMAL_INDICES, collapse = ", ")))
    }
  }
}

cat("\n=== APPLYING TRAINING DATA NORMALIZATION ===\n")
# Attempt to auto-add PPI *before* normalization. normalize_mesma_data will
# create a PPI column filled with NA if missing, which creates a catch-22
# where later checks for the presence of a PPI column incorrectly assume
# it should not be computed. Try to populate PPI from joined barren
# observations (or MESMA_DVI_SOIL override) now so the normalization step
# preserves the computed values.
## --- FIX: Filter dataset to years used in January script BEFORE computing soil baseline ---
if (!"pheno_year" %in% names(df) && "date" %in% names(df)) df$pheno_year <- assign_pheno_year(as.Date(df$date))
cat("Filtering training data to 1985-2025 (phenological years: March-February) before PPI baseline calculation...\n")
  df <- df |> dplyr::filter(pheno_year >= 1985 & pheno_year <= 2025)
if (nrow(df) == 0) cat("[NOTICE] Year filter removed all rows; PPI baseline calculation will proceed on empty data and likely be skipped.\n")
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
norm_result <- normalize_mesma_data(df, cols = OPTIMAL_INDICES, lat_default = 40.2)
df <- norm_result$df
INDEX_SCALES <- norm_result$INDEX_SCALES

TRAINING_NORM_PARAMS <- list(
  INDEX_SCALES = INDEX_SCALES,
  INDEX_SCALES_STAGE2 = list()  # Will be populated after location mapping
)
cat(sprintf("Stored normalization params: INDEX_SCALES for %d indices\n",
            length(INDEX_SCALES)))
cat("(Stage 2 normalization will be computed after location mapping)\n")
cat("===========================================\n\n")


if ("location_id" %in% names(df)) {
  df$location_id <- as.character(df$location_id)
  bad_idx <- which(is.na(df$location_id) | df$location_id == "L_NA_NA")
  if (length(bad_idx) > 0L) {
    df$location_id[bad_idx] <- as.character(bad_idx)
    cat(sprintf("[NOTICE] Replaced %d invalid location_id entries with row numbers\n", length(bad_idx)))
  }
}

if (!"date" %in% names(df) && "Date" %in% names(df)) df$date <- df$Date
if ("date" %in% names(df)) {
  df$date <- as.Date(df$date)
  if (!"pheno_year" %in% names(df)) df$pheno_year <- assign_pheno_year(df$date)
}

## No external GeoJSON: construct location map directly from CSV lat/lon
if (all(c("lon", "lat") %in% names(df))) {
  df$location_id <- make_location_id(df$lon, df$lat)
  # Ensure mapping columns exist (fill with NA if missing)
  if (!"Veg" %in% names(df)) df$Veg <- NA_character_
  if (!"no soil" %in% names(df)) df$`no soil` <- NA_real_
  # Build a minimal gpts_map from unique lat/lon combos in the CSV
  # Use per-location aggregation: pick the first non-missing lat/lon and the
  # first non-missing Veg / no soil for that location (avoids losing Veg when
  # the first row happens to have NA)
  gpts_map <- df |>
    dplyr::group_by(location_id) |>
    dplyr::summarise(
      lat = if ("lat" %in% names(df)) { v <- na.omit(lat); if (length(v)>0) v[1] else NA_real_ } else NA_real_,
      lon = if ("lon" %in% names(df)) { v <- na.omit(lon); if (length(v)>0) v[1] else NA_real_ } else NA_real_,
      Veg = if ("Veg" %in% names(df)) { v <- na.omit(as.character(Veg)); if (length(v)>0) tolower(v[1]) else NA_character_ } else NA_character_,
      `no soil` = if ("no soil" %in% names(df)) { v <- na.omit(as.numeric(`no soil`)); if (length(v)>0) v[1] else NA_real_ } else NA_real_,
      .groups = "drop"
    )
  gpts_map$location_row <- as.character(seq_len(nrow(gpts_map)))
  gpts_map$location_id_seq <- gpts_map$location_row
  cat(sprintf("[NOTICE] Constructed gpts_map from %d unique lat/lon combos in CSV\n", nrow(gpts_map)))
} else {
  cat("[NOTICE] No 'lat'/'lon' columns found in CSV; GeoJSON mapping disabled.\n")
  gpts_map <- data.frame(location_id = character(0), location_row = character(0), Veg = character(0), `no soil` = numeric(0), stringsAsFactors = FALSE)
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
  gpts_map <- data.frame(location_id = character(0), location_row = character(0), Veg = character(0), `no soil` = numeric(0), stringsAsFactors = FALSE)
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
  if ("no soil.geo" %in% names(joined)) {
    if ("no soil" %in% names(joined)) {
      joined$`no soil` <- ifelse(is.na(joined$`no soil`), joined$`no soil.geo`, joined$`no soil`)
    } else {
      joined$`no soil` <- joined$`no soil.geo`
    }
    joined$`no soil.geo` <- NULL
  }

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
      if ("no soil.geo" %in% names(joined2)) {
        if ("no soil" %in% names(joined2)) {
          joined2$`no soil` <- ifelse(is.na(joined2$`no soil`), joined2$`no soil.geo`, joined2$`no soil`)
        } else {
          joined2$`no soil` <- joined2$`no soil.geo`
        }
        joined2$`no soil.geo` <- NULL
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

  df <- joined

  matched_locs <- length(intersect(na.omit(unique(as.character(df$location_id))), na.omit(unique(as.character(gpts_map$location_id)))))
  cat(sprintf("[NOTICE] GeoJSON join results - Veg before=%d after=%d; matched location_id strings=%d\n", pre_non_na, post_non_na, matched_locs))

  # === STAGE 1 TRAINING DATA DIAGNOSTIC ===
  cat("\n=== STAGE 1 TRAINING DATA DIAGNOSTIC ===\n")
  if ("no soil" %in% names(df)) {
    cat(sprintf("Total rows in joined data: %d\n", nrow(df)))
    cat(sprintf("Rows with Veg='barren': %d\n", sum(tolower(df$Veg) == "barren", na.rm = TRUE)))
    cat(sprintf("Rows with non-NA no soil: %d\n", sum(!is.na(df$`no soil`))))
    cat(sprintf("Rows with no soil ≈ 1 (pure veg, threshold 0.01): %d\n",
                sum(!is.na(df$`no soil`) & abs(df$`no soil` - 1) < 0.01, na.rm = TRUE)))
    cat(sprintf("Rows with no soil > 0.95: %d\n",
                sum(!is.na(df$`no soil`) & df$`no soil` > 0.95, na.rm = TRUE)))
    cat(sprintf("Rows with no soil ≈ 0 (barren, threshold 0.01): %d\n",
                sum(!is.na(df$`no soil`) & abs(df$`no soil`) < 0.01, na.rm = TRUE)))

    # Show distribution of no soil values
    no_soil_summary <- summary(df$`no soil`)
    cat("no soil value distribution:\n")
    print(no_soil_summary)

    # Check for specific vegetation types with high no soil
    veg_types <- unique(df$Veg[!is.na(df$Veg) & df$Veg != ""])
    cat("\nPure vegetation (no soil ≈ 1) by type:\n")
    for (veg_type in veg_types) {
      if (tolower(veg_type) != "barren") {
        pure_veg_count <- sum(!is.na(df$`no soil`) &
                               df$Veg == veg_type &
                               abs(df$`no soil` - 1) < 0.01, na.rm = TRUE)
        cat(sprintf("  %s: %d rows\n", veg_type, pure_veg_count))
      }
    }
  } else {
    cat("[WARNING] 'no soil' column not found in joined data!\n")
  }
  cat("=========================================\n\n")
  }

  # Replace the "PPI calculation disabled" block with:
  if (exists("MESMA_NO_TOPLEVEL") && isTRUE(MESMA_NO_TOPLEVEL)) {
    cat("MESMA_NO_TOPLEVEL set: skipping PPI calculation and related top-level processing\n")
  } else if (exists("add_ppi_columns")) {
    cat("Calculating PPI...\n")
    # Delegate the auto-add logic to the helper in ppi_helpers.R
    if (exists("auto_add_ppi_columns")) {
      ppi_res <- auto_add_ppi_columns(df)
      df <- ppi_res$df
      if (isTRUE(ppi_res$added)) {
        cat(sprintf("[PPI] Auto-added PPI to dataset (reason: %s)\n", ppi_res$reason))
      } else {
        cat("[NOTICE] No barren observations found after location mapping; PPI not auto-added. To force PPI computation set MESMA_DVI_SOIL or ensure input provides Veg/'no soil' values.\n")
      }
    } else {
      cat("[NOTICE] PPI helper not available; skipping auto-add of PPI.\n")
    }
  }
  if (!"zenith.angle" %in% names(df)) df$zenith.angle <- NA_real_

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
  
  
  cat("\n=== STAGE 2 NORMALIZATION SKIPPED (Using single Z-score normalization) ===\n")
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
cat(sprintf(
  "Training dataset: %d rows from %d locations\n",
  nrow(df_train), length(unique(df_train$location_id))
))

# Filter out observations with critical snow or dust contamination
if ("NDSI" %in% names(df_train) && "NDDI" %in% names(df_train)) {
  snow_count <- sum(df_train$NDSI > 0.4, na.rm = TRUE)
  dust_count <- sum(df_train$NDDI > 0.18, na.rm = TRUE)
  total_before <- nrow(df_train)
  df_train <- df_train[!(df_train$NDSI > 0.4 | df_train$NDDI > 0.18), , drop = FALSE]
  total_after <- nrow(df_train)
  filtered <- total_before - total_after
  cat(sprintf("Filtered out %d observations with snow (NDSI > 0.4) or dust (NDDI > 0.18) contamination\n", filtered))
  cat(sprintf("Training dataset after contamination filtering: %d rows from %d locations\n", total_after, length(unique(df_train$location_id))))
} else {
  cat("[WARNING] NDSI or NDDI not found in training data; skipping contamination filtering\n")
}

df_test <- df
cat(sprintf(
  "Testing dataset: %d rows from %d locations\n",
  nrow(df_test), length(unique(df_test$location_id))
))

df_full <- df
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

cleanup_parallel <- setup_parallel_backend()
on.exit(cleanup_parallel(), add = TRUE)

# Robust memory cleanup helper — removes large temporary objects matching common patterns
cleanup_memory <- function(verbose = TRUE, confirm = FALSE) {
  patterns <- c("^raw_df$", "^df_full$", "^df_train$", "^df_test$", "^stage1_data$", "^stage2_data$", "^X$", "^X_clean$", "^X_pca$", "^pca_res$", "^pca_result$", "^all_features$", "^all_labels$", "^boot_", "^boot_preds$", "^boot_slopes$", "^boot_coefs$", "^results_list$", "^inference_", "^COMPRESSED_", "^mesma_cache$", "^cached_", "^training_", "^loc_", "^locations$", "^temp_", "^tmp_")
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
  gc(); gc()
  if (isTRUE(verbose)) cat("[CLEANUP] Memory freed (gc called).\n")
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

# Removed legacy functions: `recommended_workers()` and `compute_inter_class_similarity_table()` — no longer used.


# (Removed leftover portion of `compute_inter_class_similarity_table`.)





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

compute_pca_lda_weights <- function(lib_df, avail_idx, pca_variance_threshold = PCA_VARIANCE_THRESHOLD, lda_weight_floor = LDA_WEIGHT_FLOOR) {
  
  cat("\n=== COMPUTING PCA → LDA FEATURE WEIGHTS ===\n")
  
  veg_types <- unique(na.omit(lib_df$Veg))
  veg_types <- veg_types[veg_types != "" & tolower(veg_types) != "barren"]
  
  if (length(veg_types) < 2) {
    warning("Need at least 2 vegetation classes for LDA")
    return(NULL)
  }
  
  all_features <- list()
  all_labels <- c()
  
  for (veg in veg_types) {
    veg_data <- lib_df[tolower(lib_df$Veg) == tolower(veg), ]
    traces <- unique(veg_data[, c("location_id", "pheno_year")])

    for (i in seq_len(nrow(traces))) {
      loc <- traces$location_id[i]
      yr <- traces$pheno_year[i]

      dly_year <- veg_data[veg_data$location_id == loc & veg_data$pheno_year == yr, ]
      n_unique_doys <- length(unique(dly_year$doy))
      if (n_unique_doys < 5) next
      
      raw_mat <- build_pentad_matrix(dly_year, avail_idx)
      if (is.null(raw_mat)) next
      
      compressed <- as.numeric(raw_mat)  # N_TEMPORAL_BINS * K vector

      # Removed valid_frac check - process all available data

      compressed[!is.finite(compressed)] <- 0
      
      all_features[[length(all_features) + 1]] <- compressed
      all_labels <- c(all_labels, veg)
    }
  }
  
  if (length(all_features) < 20) {
    warning("Insufficient traces for PCA-LDA")
    return(NULL)
  }
  
  X <- do.call(rbind, all_features)
  y <- factor(all_labels)
  
  cat(sprintf("Feature matrix: %d samples × %d features (pentad resolution)\n", nrow(X), ncol(X)))
  cat(sprintf("Class distribution: %s\n", 
              paste(sprintf("%s=%d", levels(y), table(y)), collapse=", ")))
  
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  X_centered[!is.finite(X_centered)] <- 0
  
  col_vars <- apply(X_centered, 2, var, na.rm = TRUE)
  keep_cols <- col_vars > 1e-10
  if (sum(keep_cols) < 10) {
    warning("Too few variable features after filtering")
    return(NULL)
  }
  X_filtered <- X_centered[, keep_cols, drop = FALSE]
  
  cat(sprintf("After variance filtering: %d features\n", ncol(X_filtered)))
  
  pca_result <- prcomp(X_filtered, center = FALSE, scale. = FALSE)
  
  cum_var <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
  n_pcs <- which(cum_var >= pca_variance_threshold)[1]
  if (is.na(n_pcs)) n_pcs <- length(cum_var)
  n_pcs <- max(n_pcs, length(veg_types))  # At least as many PCs as classes
  n_pcs <- min(n_pcs, ncol(X_filtered), nrow(X_filtered) - 1)
  
  cat(sprintf("PCA: Retaining %d PCs (%.1f%% variance)\n", n_pcs, 100 * cum_var[n_pcs]))
  
  X_pca <- pca_result$x[, 1:n_pcs, drop = FALSE]
  
  class_counts <- table(y)
  if (any(class_counts < 3)) {
    warning("Some classes have fewer than 3 samples")
    keep_classes <- names(class_counts)[class_counts >= 3]
    keep_rows <- y %in% keep_classes
    X_pca <- X_pca[keep_rows, , drop = FALSE]
    y <- factor(y[keep_rows])
  }
  
  if (length(levels(y)) < 2) {
    warning("Fewer than 2 classes after filtering")
    return(NULL)
  }
  
  n_total_after <- nrow(X_pca)
  n_classes_after <- length(unique(y))
  max_pcs_allowed_after <- max(1, n_total_after - n_classes_after)
  if (n_pcs > max_pcs_allowed_after) {
    warning(sprintf("Reducing number of PCs from %d to %d to satisfy p <= N - K after class filtering (N=%d,K=%d)", n_pcs, max_pcs_allowed_after, n_total_after, n_classes_after))
    n_pcs <- max_pcs_allowed_after
    if (n_pcs < 1) {
      warning("Not enough samples to perform LDA after class filtering; skipping PCA-LDA weights")
      return(NULL)
    }
    X_pca <- X_pca[, 1:n_pcs, drop = FALSE]
  }

  rank_pca <- qr(X_pca)$rank
  n_total_after <- nrow(X_pca)
  n_classes_after <- length(unique(y))
  cat(sprintf("Pre-LDA diagnostic (pca_lda): N=%d, K=%d, n_pcs=%d, rank(X_pca)=%d\n", n_total_after, n_classes_after, ncol(X_pca), rank_pca))
  if (rank_pca < ncol(X_pca)) {
    cat(sprintf("X_pca rank-deficient: reducing n_pcs from %d to %d (rank)\n", ncol(X_pca), rank_pca))
    n_pcs <- rank_pca
    if (n_pcs < 1) {
      warning("Not enough PCs for LDA after rank reduction; skipping PCA-LDA weight computation")
      return(NULL)
    }
    X_pca <- X_pca[, 1:n_pcs, drop = FALSE]
  }
  max_pcs_allowed_after <- max(1, n_total_after - n_classes_after)
  if (ncol(X_pca) > max_pcs_allowed_after) {
    warning(sprintf("Reducing PCs from %d to %d to satisfy p <= N - K (N=%d,K=%d)", ncol(X_pca), max_pcs_allowed_after, n_total_after, n_classes_after))
    n_pcs <- max_pcs_allowed_after
    X_pca <- X_pca[, 1:n_pcs, drop = FALSE]
  }

  min_n_pcs_pca <- max(2, length(unique(y)) - 1)
  lda_result <- tryCatch({
    safe_lda_call(X_pca, y, min_n_pcs = min_n_pcs_pca)
  }, error = function(e) {
    warning(sprintf("LDA failed: %s", e$message))
    NULL
  })
  
  if (is.null(lda_result)) return(NULL)
  
  
  W_lda <- lda_result$scaling  # n_pcs × n_discriminants
  R_pca <- pca_result$rotation[, 1:n_pcs, drop = FALSE]  # n_filtered × n_pcs
  
  W_combined <- R_pca %*% W_lda  # n_filtered × n_discriminants
  
  svd_vals <- lda_result$svd
  prop_var <- svd_vals / sum(svd_vals)
  
  if (ncol(W_combined) > 1) {
    feature_importance <- rowSums(abs(W_combined) %*% diag(prop_var))
  } else {
    feature_importance <- abs(W_combined[, 1])
  }
  
  full_weights <- rep(0, ncol(X_centered))
  full_weights[keep_cols] <- feature_importance

  # No normalization or thresholding - keep weights at raw LDA-derived scale

  cat(sprintf("LDA weights computed: %d non-zero features\n", sum(full_weights > 0)))
  
  list(
    weights = full_weights,
    pca = pca_result,
    lda = lda_result,
    n_pcs = n_pcs,
    keep_cols = keep_cols,
    feature_names = colnames(X)
  )
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
    
    mat <- build_pentad_matrix(sub, feature_cols) # 73 x K
    if(is.null(mat)) next
    
    vec <- as.numeric(mat)
    vec[!is.finite(vec)] <- NA 
    
    X_raw[[length(X_raw)+1]] <- vec
    lbl <- names(sort(table(sub[[class_col]]), decreasing=TRUE))[1]
    y_labels <- c(y_labels, lbl)
  }
  
  if(length(X_raw) < 10) return(NULL)
  X_mat <- do.call(rbind, X_raw)
  
  
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

  # Debugging: report PPI coverage/variance if requested
  if (exists('DEBUG_STAGE1_VERBOSE', envir = globalenv()) && isTRUE(get('DEBUG_STAGE1_VERBOSE', envir = globalenv()))) {
    n_bins <- N_TEMPORAL_BINS
    ppi_idx <- which(feature_cols == 'PPI')
    if (length(ppi_idx) == 1) {
      col_start <- (ppi_idx - 1) * n_bins + 1
      col_end <- ppi_idx * n_bins
      ppi_block <- X_mat[, col_start:col_end, drop = FALSE]
      finite_counts <- apply(is.finite(ppi_block), 2, sum)
      ppi_sds <- apply(ppi_block, 2, sd, na.rm = TRUE)
      cat(sprintf("[DEBUG_STAGE1_PPI] PPI pentad finite counts (min,max) = (%d,%d); SDs (min,max,mean) = (%.6g, %.6g, %.6g)\n",
                  min(finite_counts), max(finite_counts), min(ppi_sds, na.rm = TRUE), max(ppi_sds, na.rm = TRUE), mean(ppi_sds, na.rm = TRUE)))
      # Show unique values summary across PPI block
      ppi_vals <- as.numeric(ppi_block)
      ppi_vals <- ppi_vals[is.finite(ppi_vals)]
      if (length(ppi_vals) > 0) {
        cat(sprintf("[DEBUG_STAGE1_PPI] PPI values: unique_count=%d, min=%.6g, max=%.6g, median=%.6g\n",
                    length(unique(ppi_vals)), min(ppi_vals), max(ppi_vals), median(ppi_vals)))
      } else {
        cat("[DEBUG_STAGE1_PPI] No finite PPI values in X_mat\n")
      }
    }
  }
  
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
    warning("LDA failed within train_feature_pipeline, returning NULL.")
    return(NULL)
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

  # If PPI ended up with a very small weight (and thus may be zeroed later),
  # print diagnostics when DEBUG_STAGE1_VERBOSE is enabled so we can trace
  # why PPI is not contributing to weights.
  if (exists('DEBUG_STAGE1_VERBOSE', envir = globalenv()) && isTRUE(get('DEBUG_STAGE1_VERBOSE', envir = globalenv()))) {
    idx_ppi <- which(feature_cols == 'PPI')
    if (length(idx_ppi) == 1) {
      # compute mean weight for the PPI index across its pentads
      ppi_mean_w <- mean(final_weights[((idx_ppi-1)*N_TEMPORAL_BINS + 1):(idx_ppi*N_TEMPORAL_BINS)], na.rm = TRUE)
      cat(sprintf("[DEBUG_STAGE1_PPI] PPI mean weight after adjustment = %.6g\n", ppi_mean_w))
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

  for (j in seq_along(avail_idx)) {
    idx <- avail_idx[j]
    if (!idx %in% names(dly_year)) next

    vals_by_pentad <- tapply(dly_year[[idx]], dly_year$pentad, function(v) {
      v <- v[is.finite(v)]
      if (length(v) == 0) return(NA_real_)

      q_high <- quantile(v, 0.999, na.rm = TRUE)
      q_low <- quantile(v, 0.001, na.rm = TRUE)
      v_clipped <- pmax(q_low, pmin(v, q_high))

      median(v_clipped)
    })

    pentad_values <- as.integer(names(vals_by_pentad))
    valid <- pentad_values >= 1 & pentad_values <= N_TEMPORAL_BINS
    pentad_mat[pentad_values[valid], j] <- vals_by_pentad[valid]
    
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

weighted_cosine_similarity <- function(a, b, weights = NULL) {
  if (is.null(weights) || length(weights) != length(a)) {
    return(sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2))))
  }
  
  a_w <- a * sqrt(weights)
  b_w <- b * sqrt(weights)
  
  sum(a_w * b_w) / (sqrt(sum(a_w^2)) * sqrt(sum(b_w^2)))
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

cos_sim <- function(a, b) {
  da <- sqrt(safe_dot(a, a)); db <- sqrt(safe_dot(b, b))
  if (da == 0 || db == 0) return(0)
  safe_dot(a, b) / (da * db)
}

.run_map <- function(X, FUN, show_pb = TRUE) {
  f_FUN <- FUN
  
  if (!PARALLEL_ENABLE) {
    lapply(X, function(x) { f_FUN(x) })
  } else {
    if (!requireNamespace("future.apply", quietly = TRUE) || !requireNamespace("future", quietly = TRUE)) {
      stop("future.apply and future packages required for parallel processing")
    }

    old_plan <- future::plan()
    options(future.globals.maxSize = max(getOption("future.globals.maxSize", 2e9), 12e9))
    future::plan(future::multisession, workers = PARALLEL_WORKERS)
    on.exit(future::plan(old_plan))

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

df <- df |> filter(pheno_year >= 1985 & pheno_year <= 2025)

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

df_train$location_id <- as.character(df_train$location_id)
df_train <- dplyr::left_join(df_train, gpts_map, by = "location_id", suffix = c("", ".geo"))
if ("Veg.geo" %in% names(df_train)) {
  df_train$Veg <- ifelse(is.na(df_train$Veg) | df_train$Veg == "", df_train$Veg.geo, df_train$Veg)
  df_train$Veg.geo <- NULL
}
if ("no soil.geo" %in% names(df_train)) {
  if ("no soil" %in% names(df_train)) {
    df_train$`no soil` <- ifelse(is.na(df_train$`no soil`), df_train$`no soil.geo`, df_train$`no soil`)
  } else {
    df_train$`no soil` <- df_train$`no soil.geo`
  }
  df_train$`no soil.geo` <- NULL
}
df_train$Veg <- tolower(df_train$Veg)

if (!"date" %in% names(df)) stop("Input CSV must contain a 'date' column")
df$date <- as.Date(df$date)
if (!"location_id" %in% names(df)) stop("Input CSV must contain a 'location_id' column")

if ("Veg" %in% names(df)) df$Veg <- tolower(as.character(df$Veg))

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

df$doy <- pheno_doy(df$date)  # Use phenological DOY (March 1 = day 1)
df$doy[df$doy < 1 | df$doy > 366] <- NA_integer_

cat("\n=== BARREN LOADING DEBUG ===\n")
if (exists("gpts_map") && nrow(gpts_map) > 0) {
  barren_in_gpts <- gpts_map[tolower(gpts_map$Veg) == "barren" | (!is.na(gpts_map$`no soil`) & gpts_map$`no soil` == 0), ]
  cat(sprintf("Barren locations in gpts_map: %d\n", nrow(barren_in_gpts)))
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

  # Fail fast if the input & joined metadata do not contain Vegetation or no soil info
  veg_rows_present <- sum(!is.na(df$Veg) & df$Veg != "")
  no_soil_present <- sum(!is.na(df$`no soil`))
  if (veg_rows_present == 0 && no_soil_present == 0) {
    stop(paste0("No vegetation metadata found after mapping (no 'Veg' values and no 'no soil' values).\n",
                "Please ensure your INPUT_CSV contains a 'vegetation' column (mapped to 'Veg') or a 'no soil' column,\n",
                "or provide GeoJSON location metadata containing 'Veg' or 'no soil' to join on.\n",
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
  # Identify based on "no soil" fraction if available
  ns_cols <- c("no soil", "no.soil", "no_soil")
  ns_col <- intersect(ns_cols, names(df))[1]
  if (!is.na(ns_col)) {
    vals <- as.numeric(as.character(df[[ns_col]]))
    barren_idx <- barren_idx | (is.finite(vals) & vals > 0.5)
  }

  user_dvi_soil <- suppressWarnings(as.numeric(ifelse(nzchar(Sys.getenv("MESMA_DVI_SOIL")), Sys.getenv("MESMA_DVI_SOIL"), NA)))

  if (any(barren_idx, na.rm = TRUE)) {
    dvi_soil_from_join <- mean(df$DVI[barren_idx], na.rm = TRUE)
    cat(sprintf("[PPI] Auto-adding PPI: using soil baseline from joined data (dvi_soil=%.6f, n=%d)\n", dvi_soil_from_join, sum(barren_idx, na.rm = TRUE)))
    ## --- NEW: Save this value globally for later use ---
    assign("GLOBAL_TRAINING_DVI_SOIL", dvi_soil_from_join, envir = globalenv())
    cat(sprintf("[PPI] Saved dynamic baseline (%.6f) for inference step.\n", dvi_soil_from_join))
    ## --------------------------------------------------
    df <- add_ppi_columns(df, dvi_soil = dvi_soil_from_join)
    if ("PPI" %in% names(df)) candidate_indices <- unique(c(candidate_indices, "PPI"))
  } else if (!is.na(user_dvi_soil) && is.finite(user_dvi_soil)) {
    cat(sprintf("[PPI] Auto-adding PPI using MESMA_DVI_SOIL override: %.6f\n", user_dvi_soil))
    df <- add_ppi_columns(df, dvi_soil = user_dvi_soil)
    if ("PPI" %in% names(df)) candidate_indices <- unique(c(candidate_indices, "PPI"))
  } else {
    cat("[NOTICE] Candidate indices computed from existing indices and raw bands; PPI could not be auto-added (no barren observations found after join).\n")
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

if (!isTRUE(SKIP_MOVING_VARIANCE)) {
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


  
  if (length(var_cols) > 0) {
    cat(sprintf("Adding %d moving variance features to avail: %s\n", length(var_cols), paste(var_cols, collapse=", ")))
    avail <- unique(c(avail, var_cols))
  }
}



} # End of if (!isTRUE(SKIP_MOVING_VARIANCE))

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

build_barren_veg_library <- function(df_local, avail_idx, min_samples = 5) {
  if (!"doy" %in% names(df_local)) df_local$doy <- pheno_doy(as.Date(df_local$date))  # Use phenological DOY
  
  # Ensure alternate no-soil column spellings are recognized locally. The
  # top-level script normally maps these earlier, but be defensive in case
  # this function is called in isolation.
  if ("no_soil" %in% names(df_local) && !"no soil" %in% names(df_local)) {
    df_local$`no soil` <- df_local$no_soil
    cat("[Stage1] Mapped 'no_soil' -> 'no soil' inside build_barren_veg_library\n")
  }
  if ("no-soil" %in% names(df_local) && !"no soil" %in% names(df_local)) {
    df_local$`no soil` <- df_local$`no-soil`
    cat("[Stage1] Mapped 'no-soil' -> 'no soil' inside build_barren_veg_library\n")
  }
  if ("no soil" %in% names(df_local) && !"no.soil" %in% names(df_local)) {
    df_local$`no.soil` <- df_local$`no soil`
    cat("[Stage1] Mapped 'no soil' -> 'no.soil' inside build_barren_veg_library\n")
  }
  
  barren_by_veg <- if ("Veg" %in% names(df_local)) {
    df_local[!is.na(df_local$Veg) & tolower(df_local$Veg) == "barren", , drop = FALSE]
  } else {
    df_local[FALSE, , drop = FALSE]  # empty
  }
  
  barren_rows <- barren_by_veg
  
  cat(sprintf("[Stage1] Barren rows: %d from Veg=='barren'\n",
              nrow(barren_rows)))
  
  if (nrow(barren_rows) < min_samples) {
    cat(sprintf("[Stage1] Insufficient barren training data. Found %d barren observations, but require at least %d. Returning NULL.\n", nrow(barren_rows), min_samples))
    return(NULL)
  }
  
  veg_rows <- if ("no soil" %in% names(df_local)) {
    pure_veg <- df_local[!is.na(df_local$`no soil`) & {
      val <- as.numeric(as.character(df_local$`no soil`))
      abs(val - 1) < 0.01
    }, , drop = FALSE]
    
    if (nrow(pure_veg) < min_samples) {
      stop(sprintf("[Stage1] CRITICAL ERROR: Insufficient pure vegetation training data. Found %d pure vegetation observations (no soil ≈ 1), but require at least %d. Cannot proceed with MESMA analysis.", nrow(pure_veg), min_samples))
    }
    
    pure_veg
  } else {
    # Attempt to infer pure vegetation rows when 'no soil' is missing by
    # selecting the top observations per vegetation class by a robust
    # vegetation-sensitive index (prefer MSAVI2, then NDVI, then DVI/NIRv).
    candidate_idx <- intersect(c("MSAVI2", "NDVI", "DVI", "NIRv"), avail_idx)
    if (length(candidate_idx) == 0) {
      stop("[Stage1] CRITICAL ERROR: 'no soil' column is missing and no suitable index (MSAVI2/NDVI/DVI/NIRv) is available to infer pure vegetation samples.")
    }

    cat(sprintf("[Stage1] 'no soil' missing; attempting to infer pure vegetation using index '%s' (candidates: %s)\n",
                candidate_idx[1], paste(candidate_idx, collapse = ", ")))

    inferred_pure <- list()
    for (veg_type in unique(df_local$Veg[!is.na(df_local$Veg) & tolower(df_local$Veg) != "barren"])) {
      sub <- df_local[!is.na(df_local$Veg) & tolower(df_local$Veg) == tolower(veg_type), , drop = FALSE]
      if (nrow(sub) == 0) next
      idx_col <- candidate_idx[1]
      if (!(idx_col %in% names(sub))) next

      # Try iteratively less strict quantiles until we have enough samples
      probs <- c(0.98, 0.95, 0.90, 0.85, 0.80, 0.75)
      sel <- data.frame()
      for (p in probs) {
        thr <- suppressWarnings(stats::quantile(sub[[idx_col]], probs = p, na.rm = TRUE))
        if (!is.finite(thr)) next
        sel <- sub[which(!is.na(sub[[idx_col]]) & sub[[idx_col]] >= thr), , drop = FALSE]
        if (nrow(sel) >= min_samples) break
      }

      # As a last resort, take top-N values by idx_col
      if (nrow(sel) < min_samples) {
        ord <- order(-as.numeric(sub[[idx_col]]), na.last = NA)
        if (length(ord) >= min_samples) sel <- sub[ord[seq_len(min_samples)], , drop = FALSE]
      }

      if (nrow(sel) >= min_samples) {
        inferred_pure[[veg_type]] <- sel[seq_len(min(nrow(sel), min_samples)), , drop = FALSE]
        cat(sprintf("[Stage1] Inferred %d pure-veg samples for veg='%s' using index '%s'\n", nrow(inferred_pure[[veg_type]]), veg_type, idx_col))
      } else {
        cat(sprintf("[Stage1] Could not infer enough pure-veg for veg='%s' (found %d)\n", veg_type, nrow(sel)))
      }
    }

    if (length(inferred_pure) == 0) {
      stop("[Stage1] CRITICAL ERROR: Could not infer any pure vegetation samples and 'no soil' is not present; cannot proceed.")
    }

    veg_rows <- do.call(rbind, inferred_pure)
    if (nrow(veg_rows) < min_samples) {
      stop(sprintf("[Stage1] CRITICAL ERROR: Inferred only %d pure vegetation rows (minimum required %d). Cannot proceed.", nrow(veg_rows), min_samples))
    }
    cat(sprintf("[Stage1] Using %d inferred pure vegetation rows across veg classes\n", nrow(veg_rows)))
  }
  
  cat(sprintf("[Stage1] After filtering: barren=%d, veg=%d rows\n", 
              nrow(barren_rows), nrow(veg_rows)))

  n_barren_orig <- nrow(barren_rows)
  n_veg_orig <- nrow(veg_rows)
  stage1_target <- NULL
  if (exists("STAGE1_TARGET_SIZE")) {
    stage1_target <- as.integer(STAGE1_TARGET_SIZE)
  } else {
    stage1_target <- 2 * min(n_barren_orig, n_veg_orig)
    if (exists('TESTING_MAX_PER_VEG') && isTRUE(TESTING_MODE)) {
      stage1_target <- min(as.integer(stage1_target), as.integer(TESTING_MAX_PER_VEG))
    }
  }
  if (is.null(stage1_target) || stage1_target < 1) {
    stage1_target <- max(n_barren_orig, n_veg_orig)
  }

  cat(sprintf("[Stage1] Original sample sizes - Barren: %d, VegPure: %d\n", n_barren_orig, n_veg_orig))

  if (n_barren_orig != stage1_target || n_veg_orig != stage1_target) {
    cat(sprintf("[Stage1] Stage1 balanced sample sizes - Barren: %d, VegPure: %d (target=%d)\n", stage1_target, stage1_target, stage1_target))

    prepare_and_balance <- function(df_class, target_n) {
      if (nrow(df_class) == 0) return(df_class)
      
      df_clean <- as.data.frame(df_class, stringsAsFactors = FALSE)
      
      cols_to_keep <- sapply(df_clean, function(x) !all(is.na(x)))
      df_clean <- df_clean[, cols_to_keep, drop = FALSE]
      
      if (nrow(df_clean) >= target_n) {
        set.seed(123)
        if (nrow(df_clean) > target_n) df_clean <- df_clean[sample.int(nrow(df_clean), target_n), , drop = FALSE]
        return(df_clean)
      } else {
        return(augment_minority_class(df_clean, target_n))
      }
    }

    barren_rows <- prepare_and_balance(barren_rows, stage1_target)
    veg_rows <- prepare_and_balance(veg_rows, stage1_target)
    
    cat(sprintf("[Stage1] Balanced sample sizes - Barren: %d, VegPure: %d (target=%d)\n", 
                nrow(barren_rows), nrow(veg_rows), stage1_target))

    cat("[Stage1 Balancing] Rebuilt STAGE1_LIB from balanced stage-1 data\n")
  }
  
  if (nrow(barren_rows) < min_samples) {
     cat(sprintf("[Stage1] Insufficient training data: barren=%d (need >=%d)\n", nrow(barren_rows), min_samples))
     return(NULL)
  }

  if (nrow(veg_rows) < min_samples) {
    cat(sprintf("[Stage1] Insufficient training data: veg=%d (need >=%d)\n", 
                nrow(veg_rows), min_samples))
    return(NULL)
  }
  
  cat(sprintf("[Stage1] Building barren-veg library from barren=%d, veg=%d rows\n", 
              nrow(barren_rows), nrow(veg_rows)))
  
  stage1_lib <- list()
  
  barren_rows$pentad <- doy_to_pentad(barren_rows$doy)
  veg_rows$pentad <- doy_to_pentad(veg_rows$doy)
  
  barren_lib <- list()
  dropped_barren <- character(0)
  for (idx in avail_idx) {
    if (!idx %in% names(barren_rows)) {
      dropped_barren <- c(dropped_barren, paste0(idx, " (not in data)"))
      next
    }
    vals_by_pentad <- tapply(seq_along(barren_rows[[idx]]), barren_rows$pentad, function(indices) {
      vals <- barren_rows[[idx]][indices]
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) return(NA_real_)
      median(vals)
    })
    mu <- rep(NA_real_, N_TEMPORAL_BINS)
    if (length(vals_by_pentad) > 0) {
      pentad_values <- as.integer(names(vals_by_pentad))
      valid_pentad <- pentad_values >= 1 & pentad_values <= N_TEMPORAL_BINS
      mu[pentad_values[valid_pentad]] <- vals_by_pentad[valid_pentad]
    }
    if (all(!is.finite(mu))) {
      n_finite <- sum(is.finite(barren_rows[[idx]]))
      n_total <- length(barren_rows[[idx]])
      dropped_barren <- c(dropped_barren, sprintf("%s (all pentads NA; %d/%d finite values in raw data)", idx, n_finite, n_total))
      next
    }
    barren_lib[[idx]] <- list(mu = mu, mv = rep(NA_real_, N_TEMPORAL_BINS))
  }
  if (length(dropped_barren) > 0) {
    cat(sprintf("[Stage1] Dropped %d indices from barren library:\n", length(dropped_barren)))
    for (msg in dropped_barren) cat(sprintf("  - %s\n", msg))
  }
  
  veg_lib <- list()
  dropped_veg <- character(0)
  for (idx in avail_idx) {
    if (!idx %in% names(veg_rows)) {
      dropped_veg <- c(dropped_veg, paste0(idx, " (not in data)"))
      next
    }
    vals_by_pentad <- tapply(seq_along(veg_rows[[idx]]), veg_rows$pentad, function(indices) {
      vals <- veg_rows[[idx]][indices]
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) return(NA_real_)
      median(vals)
    })
    mu <- rep(NA_real_, N_TEMPORAL_BINS)
    if (length(vals_by_pentad) > 0) {
      pentad_values <- as.integer(names(vals_by_pentad))
      valid_pentad <- pentad_values >= 1 & pentad_values <= N_TEMPORAL_BINS
      mu[pentad_values[valid_pentad]] <- vals_by_pentad[valid_pentad]
    }
    if (all(!is.finite(mu))) {
      n_finite <- sum(is.finite(veg_rows[[idx]]))
      n_total <- length(veg_rows[[idx]])
      dropped_veg <- c(dropped_veg, sprintf("%s (all pentads NA; %d/%d finite values in raw data)", idx, n_finite, n_total))
      next
    }
    veg_lib[[idx]] <- list(mu = mu, mv = rep(NA_real_, N_TEMPORAL_BINS))
  }
  if (length(dropped_veg) > 0) {
    cat(sprintf("[Stage1] Dropped %d indices from vegetation library:\n", length(dropped_veg)))
    for (msg in dropped_veg) cat(sprintf("  - %s\n", msg))
  }
  
  if (length(barren_lib) == 0 || length(veg_lib) == 0) return(NULL)
  
  stage1_lib$barren <- barren_lib
  stage1_lib$vegetation <- veg_lib
  stage1_lib
}

unmix_vegetated_fraction <- function(dly_local, stage1_lib, avail_idx) {
  if (is.null(stage1_lib) || length(stage1_lib) == 0) return(NA_real_)
  if (!"doy" %in% names(dly_local)) dly_local$doy <- pheno_doy(as.Date(dly_local$date))  # Use phenological DOY
  
  veg_fractions <- numeric(nrow(dly_local))
  valid_count <- 0
  
  for (r in seq_len(nrow(dly_local))) {
    row <- dly_local[r, , drop = FALSE]
    doy <- as.integer(row$doy)
    if (is.na(doy) || doy < 1 || doy > 366) {
      veg_fractions[r] <- NA_real_
      next
    }
    
    pentad <- doy_to_pentad(doy)
    
    y <- numeric(0)  # Observation
    B <- numeric(0)  # Barren endmember
    V <- numeric(0)  # Vegetation endmember
    
    used_idx <- character(0)
    for (idx in avail_idx) {
      if (!idx %in% names(row)) next
      
      y_val <- as.numeric(row[[idx]])
      
      b_val <- if (!is.null(stage1_lib$barren[[idx]]) && 
                   !is.null(stage1_lib$barren[[idx]]$mu)) {
        stage1_lib$barren[[idx]]$mu[pentad]
      } else NA_real_
      v_val <- if (!is.null(stage1_lib$vegetation[[idx]]) && 
                   !is.null(stage1_lib$vegetation[[idx]]$mu)) {
        stage1_lib$vegetation[[idx]]$mu[pentad]
      } else NA_real_
      
      if (is.finite(y_val) && is.finite(b_val) && is.finite(v_val)) {
        y <- c(y, y_val)
        B <- c(B, b_val)
        V <- c(V, v_val)
        used_idx <- c(used_idx, idx)
      }
    }
    
    if (length(y) < 2) {
      veg_fractions[r] <- NA_real_
      next
    }
    
    separation <- sqrt(sum((B - V)^2))
    if (length(used_idx) > 0 && (exists('DEBUG_STAGE1_VERBOSE', envir = globalenv()) && isTRUE(get('DEBUG_STAGE1_VERBOSE', envir = globalenv())) || separation < 0.01)) {
      ratios <- sapply(seq_along(used_idx), function(i){
        b_val <- B[i]; v_val <- V[i]
        if (abs(b_val) < .Machine$double.eps && abs(v_val) < .Machine$double.eps) {
          ratio <- 1
        } else {
          ratio <- max(abs(v_val), abs(b_val)) / (min(abs(v_val), abs(b_val)) + .Machine$double.eps)
        }
        ratio
      })
      dist <- ratios - 1
      cat(sprintf("[Stage1 DOY=%d] Index ratios for %d indices (idx:ratio,dist):\n", pentad, length(used_idx)))
      for (i in seq_along(used_idx)) {
        cat(sprintf("    %s: ratio=%.3f, dist=%.3f\n", used_idx[i], ratios[i], dist[i]))
      }
    }
    if (separation < 0.01) {
      veg_fractions[r] <- NA_real_
      next
    }
    
    E <- cbind(B, V)
    sol <- solve(t(E) %*% E, t(E) %*% y)
    if (exists("project_to_simplex")) {
      w <- project_to_simplex(as.numeric(sol))
    } else {
      w <- pmax(as.numeric(sol), 0); if (sum(w) > 0) w <- w / sum(w) else w <- rep(0.5, 2)
    }
    result <- list(f2 = w[2])
    
    veg_fractions[r] <- max(result$f2, STAGE1_MIN_VEG_FRACTION)
    valid_count <- valid_count + 1
  }
  
  valid_fracs <- veg_fractions[is.finite(veg_fractions)]
  
  if (length(valid_fracs) == 0) {
    cat(sprintf("[Stage1 OLS RAW] No valid unmixing results (checked %d rows)\n", nrow(dly_local)))
    return(NA_real_)
  }
  
  veg_frac <- median(valid_fracs)
  
  cat(sprintf("[Stage1 OLS RAW] Veg fraction: %.3f (from %d/%d valid obs)\n", 
              veg_frac, length(valid_fracs), nrow(dly_local)))
  
  veg_frac
}

 


select_stage1_indices <- function(barren_lib, veg_lib, avail_idx, min_separation = 0.1) {
  separations <- sapply(avail_idx, function(idx) {
    if (!idx %in% names(barren_lib) || !idx %in% names(veg_lib)) return(NA_real_)
    b_mean <- mean(barren_lib[[idx]]$mu, na.rm = TRUE)
    v_mean <- mean(veg_lib[[idx]]$mu, na.rm = TRUE)
    if (!is.finite(b_mean) || !is.finite(v_mean)) return(NA_real_)
    if (abs(b_mean) < .Machine$double.eps && abs(v_mean) < .Machine$double.eps) {
      ratio <- 1
    } else {
      ratio <- max(abs(v_mean), abs(b_mean)) / (min(abs(v_mean), abs(b_mean)) + .Machine$double.eps)
    }
    dist <- ratio - 1
    cat(sprintf("[Stage1] Index %s: barren=%.4g, veg=%.4g, ratio=%.3f (dist=%.3f)\n", idx, b_mean, v_mean, ratio, dist))
    dist
  })
  keep <- avail_idx[!is.na(separations) & separations >= min_separation]
  if (length(keep) < 3) keep <- avail_idx[order(separations, decreasing = TRUE)[1:min(5, length(avail_idx))]]
  keep
}


build_stage1_lib <- function(STAGE1_LIB, grid_name = "full") {
  if (is.null(STAGE1_LIB)) return(NULL)
  
  barren_raw <- STAGE1_LIB$barren
  vegetation_raw <- STAGE1_LIB$vegetation
  
  avail_idx <- intersect(names(barren_raw), names(vegetation_raw))
  if (length(avail_idx) == 0) return(NULL)
  
  
  build_raw_trace <- function(endmember_lib, avail_idx) {
    raw_mat <- matrix(NA_real_, nrow = N_TEMPORAL_BINS, ncol = length(avail_idx))
    colnames(raw_mat) <- avail_idx
    
    for (j in seq_along(avail_idx)) {
      idx <- avail_idx[j]
      if (!is.null(endmember_lib[[idx]]) && !is.null(endmember_lib[[idx]]$mu)) {
        mu_vec <- endmember_lib[[idx]]$mu
        if (length(mu_vec) == 365) {
          for (p in seq_len(N_TEMPORAL_BINS)) {
            doy_start <- (p - 1) * TEMPORAL_AGGREGATION_DAYS + 1
            doy_end <- min(p * TEMPORAL_AGGREGATION_DAYS, 365)
            vals <- mu_vec[doy_start:doy_end]
            vals <- vals[is.finite(vals)]
            raw_mat[p, j] <- if (length(vals) > 0) median(vals) else NA_real_
          }
        } else if (length(mu_vec) == N_TEMPORAL_BINS) {
          raw_mat[, j] <- mu_vec
        }
      }
    }
    
    raw_mat
  }
  
  raw_barren <- build_raw_trace(barren_raw, avail_idx)
  raw_veg <- build_raw_trace(vegetation_raw, avail_idx)
  
  avail_idx_all <- intersect(names(barren_raw), names(vegetation_raw))
  raw_barren_full <- build_raw_trace(barren_raw, avail_idx_all)
  
  sim <- cos_sim(as.numeric(raw_barren), as.numeric(raw_veg))
  cat(sprintf("[build_stage1_lib] Full resolution: %d days x %d indices. Cosine Similarity(Barren, Veg) = %.4f\n", nrow(raw_barren), ncol(raw_barren), sim))

  if (is.finite(sim) && sim > 0.95) {
    cat(sprintf("[build_stage1_lib] [NOTICE] Barren and Vegetation endmembers are very similar (cos_sim=%.4f). Stage 1 unmixing may be unreliable.\n", sim))
  }

  euclidean_dist <- sqrt(sum((as.numeric(raw_barren) - as.numeric(raw_veg))^2))
  cat(sprintf("[build_stage1_lib] Euclidean distance(Barren, Veg) = %.4f\n", euclidean_dist))
  if (euclidean_dist < 0.1) {
    cat(sprintf("[build_stage1_lib] [NOTICE] Barren and Vegetation endmembers have very small Euclidean separation: %.4f\n", euclidean_dist))
  }
  norm_b <- sqrt(sum(as.numeric(raw_barren)^2, na.rm = TRUE))
  norm_v <- sqrt(sum(as.numeric(raw_veg)^2, na.rm = TRUE))
  if (is.finite(norm_b) && is.finite(norm_v) && (norm_b > 0) && (norm_v > 0)) {
    ratio <- max(norm_b / norm_v, norm_v / norm_b)
    if (ratio > 10) {
      cat(sprintf("[build_stage1_lib] [NOTICE] Barren and Vegetation full-trace norms differ substantially (barren=%.4g, veg=%.4g, ratio=%.4g). Consider normalizing stage1 endmembers to avoid projection bias.\n", norm_b, norm_v, ratio))
    }

  }
  
  cat(sprintf("[build_stage1_lib] Stage 1 indices: %s\n", paste(avail_idx, collapse = ", ")))
  
  list(barren = raw_barren, vegetation = raw_veg, raw_barren = raw_barren, raw_vegetation = raw_veg, indices = avail_idx, norm_b = norm_b, norm_v = norm_v, barren_full = raw_barren_full, indices_full = avail_idx_all)
}

if (exists('DEBUG_STAGE1_VERBOSE', envir = globalenv()) && isTRUE(get('DEBUG_STAGE1_VERBOSE', envir = globalenv())) && exists('COMPRESSED_STAGE1_LIB') && !is.null(COMPRESSED_STAGE1_LIB)) {
  cat("\n=== STAGE 1 DEBUG - Sample Unmix Diagnostics ===\n")
  if (exists('df_train') && nrow(df_train) > 0) {
    sample_rows <- df_train[df_train$Veg %in% ALLOWED_VEG, ]
    sample_rows <- head(sample_rows, 20)
    if (nrow(sample_rows) > 0) {
      for (i in seq_len(nrow(sample_rows))) {
        row <- sample_rows[i, , drop = FALSE]
        dly_year <- df_train[df_train$location_id == row$location_id & df_train$pheno_year == row$pheno_year, , drop = FALSE]
        res_compress <- compress_trace(dly_year, COMPRESSED_STAGE1_LIB$indices, budget = TEMPORAL_BUDGET)
        if (!is.null(res_compress)) {
          y_vec_raw <- res_compress$y
          # Apply same preprocessing as in fit_one_task
          y_vec <- y_vec_raw
          n_bins <- N_TEMPORAL_BINS
          for (k in seq_along(COMPRESSED_STAGE1_LIB$indices)) {
            idx_start <- (k-1)*n_bins + 1
            idx_end <- k*n_bins
            mu <- if (exists("STAGE1_PARAMS")) STAGE1_PARAMS$means[k] else 0
            sigma <- if (exists("STAGE1_PARAMS") && is.finite(STAGE1_PARAMS$sds[k]) && STAGE1_PARAMS$sds[k] > 1e-10) STAGE1_PARAMS$sds[k] else 1
            y_vec[idx_start:idx_end] <- (y_vec[idx_start:idx_end] - mu) / sigma
          }
          y_vec[!is.finite(y_vec)] <- 0
          if (exists("STAGE1_PARAMS") && !is.null(STAGE1_PARAMS$weights)) {
            # Apply sqrt(weights) to match how templates are weighted (WLS consistency)
            y_vec <- y_vec * sqrt(pmax(STAGE1_PARAMS$weights, 0))
          }
          res_stage1 <- ols_stage1_unmix(y_vec, COMPRESSED_STAGE1_LIB$barren, COMPRESSED_STAGE1_LIB$vegetation)
          cat(sprintf("Sample %d: loc=%s pheno_year=%d t=%.6f veg_frac=%.6f resid=%.6f y_norm_raw=%.4f barren_norm_raw=%.4f veg_norm_raw=%.4f\n",
                      i, as.character(row$location_id), as.integer(row$pheno_year), res_stage1$projection[1], res_stage1$veg_frac, res_stage1$residual, sqrt(sum(y_vec_raw^2)), COMPRESSED_STAGE1_LIB$norm_b, COMPRESSED_STAGE1_LIB$norm_v))
        }
      }
    }
  }
  cat("=== STAGE 1 DEBUG END ===\n\n")
}

precompute_optimized_library <- function(mesma_lib, compressed_templates_accessor, grid_type, feature_weights = NULL, stage2_norm_params = NULL) {
  opt_lib <- list()
  
  template_indices <- NULL
  for (veg in names(mesma_lib)) {
    if (length(mesma_lib[[veg]]) > 0 && !is.null(mesma_lib[[veg]][[1]]$raw_mat)) {
      template_indices <- colnames(mesma_lib[[veg]][[1]]$raw_mat)
      break
    }
  }
  
  scaling_means <- NULL
  scaling_sds <- NULL
  
  if (!is.null(stage2_norm_params) && !is.null(stage2_norm_params$INDEX_SCALES) && !is.null(template_indices)) {
    n_bins <- N_TEMPORAL_BINS
    scaling_means <- numeric(0)
    scaling_sds <- numeric(0)
    
    for (idx_name in template_indices) {
      mu <- 0
      sigma <- 1
      if (idx_name %in% names(stage2_norm_params$INDEX_SCALES)) {
        params <- stage2_norm_params$INDEX_SCALES[[idx_name]]
        mu <- params$mean
        sigma <- params$sd
      }
      scaling_means <- c(scaling_means, rep(mu, n_bins))
      scaling_sds <- c(scaling_sds, rep(sigma, n_bins))
    }
    cat(sprintf("[precompute_optimized_library] Applying Z-score scaling to templates using %d indices\n", length(template_indices)))
  }
  
  for (veg in names(mesma_lib)) {
    variants <- mesma_lib[[veg]]
    n_vars <- length(variants)
    if (n_vars == 0) next
    
    vec_list <- list()
    ids <- character(n_vars)
    valid_idx <- integer(0)
    
    for (i in seq_along(variants)) {
      vid <- variants[[i]]$variant_id
      vec <- NULL
      if (!is.null(compressed_templates_accessor[[veg]][[vid]][[grid_type]])) {
        vec <- compressed_templates_accessor[[veg]][[vid]][[grid_type]]
      } else if (!is.null(compressed_templates_accessor[[veg]][[vid]][["full"]])) {
        vec <- compressed_templates_accessor[[veg]][[vid]][["full"]]
      }
      
      if (!is.null(vec)) {
        if (!is.null(scaling_means) && length(vec) == length(scaling_means)) {
          vec <- (vec - scaling_means) / scaling_sds
          vec[!is.finite(vec)] <- 0
        }
        
        vec_list[[length(vec_list) + 1]] <- vec
        ids[length(vec_list)] <- vid
        valid_idx <- c(valid_idx, i)
      }
    }
    
    if (length(vec_list) == 0) next
    
    M <- do.call(rbind, vec_list)
    ids <- ids[1:length(vec_list)]
    
    # For similarity we apply sqrt(weights) to each feature, but keep the original
    # unweighted vectors in `M` for NNLS so the solver can apply the same sqrt weights
    # to both `E` and `y` (WLS consistency).
    M_weighted_for_norm <- M
    if (!is.null(feature_weights)) {
      if (length(feature_weights) == ncol(M)) {
        sqrt_w <- sqrt(pmax(feature_weights, 0))
        M_weighted_for_norm <- t(t(M) * sqrt_w)
      }
    }

    row_norms <- sqrt(rowSums(M_weighted_for_norm^2))
    row_norms[row_norms < 1e-9] <- 1
    M_norm <- M_weighted_for_norm / row_norms
    
    opt_lib[[veg]] <- list(
      M = M,              # Original vectors (for unmixing)
      M_norm = M_norm,    # Weighted & Normalized vectors (for similarity search)
      ids = ids
    )
  }
  return(opt_lib)
}

DEBUG_STAGE1_VERBOSE <- if (exists('DEBUG_STAGE1_VERBOSE', envir = globalenv())) {
  isTRUE(get('DEBUG_STAGE1_VERBOSE', envir = globalenv()))
} else {
  FALSE
}


unmix_stage2_compressed <- function(veg_kept, veg_frac, y, grid_type, compressed_templates_accessor, mesma_lib, topK = TOPK_VARIANTS, feature_weights = NULL, optimized_library = NULL) {
  compressed_stage1_lib <- if (exists("COMPRESSED_STAGE1_LIB")) COMPRESSED_STAGE1_LIB else NULL
  if (is.null(optimized_library)) {
    stop("[unmix_stage2_compressed] OPTIMIZED_LIBRARY is required. Ensure precompute_optimized_library() was called.")
  }
  

  comp_templates <- NULL

  if (!exists(".GEOM_DEBUG_COUNTER", envir = globalenv())) assign(".GEOM_DEBUG_COUNTER", 0L, envir = globalenv())
  debug_counter <- get(".GEOM_DEBUG_COUNTER", envir = globalenv())
  verbose_this_call <- (debug_counter < 5)
  if (verbose_this_call) assign(".GEOM_DEBUG_COUNTER", debug_counter + 1L, envir = globalenv())
  
  veg_types <- veg_kept
  if (length(veg_types) == 0) return(NULL)
  
  veg_types <- intersect(veg_types, names(mesma_lib))
  if (length(veg_types) == 0) return(NULL)
  
  top_variants <- list()
  
  y_vec <- as.numeric(y)
  # Do not apply extra weighting or L2-normalization; use raw observation vector
  y_vec[!is.finite(y_vec)] <- 0
    
    for (v in veg_types) {
      lib_v <- optimized_library[[v]]
      if (verbose_this_call || isTRUE(TESTING_MODE)) {
        if (is.null(lib_v)) {
          cat(sprintf("[DEBUG Stage2 optimized] veg=%s: optimized_library missing entry\n", v))
        } else {
          cat(sprintf("[DEBUG Stage2 optimized] veg=%s: %d variants in optimized lib\n", v, length(lib_v$ids)))
        }
      }
      if (is.null(lib_v)) next
      
      if (isTRUE(TESTING_MODE)) {
         cat(sprintf("[DEBUG Stage2] M dim: %s, y_vec len: %d\n", paste(dim(lib_v$M), collapse="x"), length(y_vec)))
      }
      sims <- tryCatch({
        # Use weighted & normalized templates for similarity (precomputed with sqrt(feature_weights))
        # Weight the observation vector consistently by sqrt(feature_weights) when provided
        y_vec_for_sim <- y_vec
        if (!is.null(feature_weights) && length(feature_weights) == length(y_vec)) {
          y_vec_for_sim <- y_vec_for_sim * sqrt(pmax(feature_weights, 0))
        }
        as.numeric(lib_v$M_norm %*% y_vec_for_sim)
      }, error = function(e) {
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG Stage2] Matrix mult failed: %s\n", e$message))
        return(NULL)
      })
      if (is.null(sims)) next
      
      ord <- order(sims, decreasing = TRUE)
      keep_idx <- ord[seq_len(min(topK, length(ord)))]
      
      selected_vecs <- lib_v$M[keep_idx, , drop = FALSE]
      selected_ids <- lib_v$ids[keep_idx]
      
      variant_list <- list()
      for (k in seq_along(keep_idx)) {
        variant_list[[k]] <- list(vec = selected_vecs[k, ], id = selected_ids[k])
      }
      top_variants[[v]] <- variant_list
    }
  
  if (length(top_variants) == 0) {
    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG Stage2] top_variants empty for veg_types=%s ; returning NULL\n", paste(veg_types, collapse=",")))
    return(NULL)
  }


  # Use raw observation vector y_vec for fitting (no preprocessing)
  if (verbose_this_call || isTRUE(TESTING_MODE)) {
    cat(sprintf("[DEBUG Stage2] y_vec length=%d, finite=%d, nan=%d\n", length(y_vec), sum(is.finite(y_vec)), sum(is.nan(y_vec))))
  }

  veg_libraries <- list()
  for (v in names(top_variants)) {
    veg_libraries[[v]] <- top_variants[[v]]
  }

  if (isTRUE(TESTING_MODE)) cat("[DEBUG Stage2] Calling stage2_ols_unmix...\n")
  geom_result <- stage2_ols_unmix(y_vec, veg_libraries, topK = topK, feature_weights = feature_weights)
  if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG Stage2] stage2_ols_unmix returned %s\n", if(is.null(geom_result)) "NULL" else "List"))
  
  if (is.null(geom_result)) {
    return(NULL)
  }
  
  w <- geom_result$proportions
  chosen <- geom_result$chosen_variants
  rmse <- geom_result$residual
  
  MIN_VEG_FRACTION <- 0.10
  w[w < MIN_VEG_FRACTION] <- 0
  if (sum(w) > 0) {
    w <- w / sum(w)
  }
  
  if (verbose_this_call) {
  
  diagnostics <- NULL
  uncertainty <- NULL
  
  if (isTRUE(ENABLE_DIAGNOSTICS) && !is.null(w) && !is.null(chosen)) {
    E_best_cols <- list()
    for (v in names(chosen)) {
      for (variant in top_variants[[v]]) {
        if (variant$id == chosen[[v]]) {
          E_best_cols[[v]] <- variant$vec
          break
        }
      }
    }
    if (length(E_best_cols) > 0) {
      E_best <- do.call(cbind, E_best_cols)
      for (j in seq_len(ncol(E_best))) {
        nj <- sqrt(sum(E_best[, j]^2))
        if (nj > 0) E_best[, j] <- E_best[, j] / nj
      }
      mesma_result <- list(w = w, chosen = chosen, rmse = rmse)
      diagnostics <- tryCatch({
        compute_diagnostics(y, E_best, w, mesma_result = mesma_result)
      }, error = function(e) NULL)
    }
  }

  }  # Close if (verbose_this_call)

  list(vegetation_proportions = w, proportions = w, fractions = w, chosen_variants = chosen, rmse = rmse, diagnostics = diagnostics, uncertainty = uncertainty)
}

cat("[DEBUG] Functions defined successfully.\n")
cat("===============================================\n\n")

# Quick self-check to ensure Stage 2 uses same NNLS behaviour as Stage 1 (TESTING_MODE only)
if (isTRUE(TESTING_MODE)) {
  cat("[TEST] Running simple Stage2 NNLS self-check...\n")
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
    res_stage2 <- stage2_ols_unmix(y_test, veg_lib_test, topK = 1)
    if (!is.null(res_stage2)) cat(sprintf("[TEST] stage2_ols_unmix proportions finite=%d, residual=%.6g\n", sum(is.finite(res_stage2$proportions)), res_stage2$residual)) else cat("[TEST] stage2_ols_unmix returned NULL\n")
  }, error = function(e) {
    cat(sprintf("[TEST] Stage2 NNLS self-check failed: %s\n", e$message))
  })
}


STAGE1_LIB <- NULL
cat("[DEBUG] About to call build_barren_veg_library...\n")
cat(sprintf("[DEBUG] df has %d rows, 'no soil' column exists: %s\n",
            nrow(df), "no.soil" %in% names(df)))
if ("no.soil" %in% names(df)) {
  n_pure_veg <- sum(!is.na(df$no.soil) & abs(as.numeric(as.character(df$no.soil)) - 1) < 0.01, na.rm = TRUE)
  cat(sprintf("[DEBUG] Found %d rows with no soil ≈ 1 (pure vegetation)\n", n_pure_veg))
}
try({
  # Attempt to ensure PPI is available for Stage-1 library construction.
  # Two cases to handle:
  # 1) PPI column does not exist at all -> try to auto-add based on joined barren observations or MESMA_DVI_SOIL
  # 2) PPI column exists but is entirely NA within barren or veg training subsets -> try to compute missing values
  if (exists("auto_add_ppi_columns")) {
    if (!"PPI" %in% names(df)) {
      ppi_try <- auto_add_ppi_columns(df)
      df <- ppi_try$df
      if (isTRUE(ppi_try$added)) {
        cat(sprintf("[PPI] Auto-added PPI prior to Stage-1 library (reason: %s)\n", ppi_try$reason))
        avail <- unique(c(avail, "PPI"))
      } else {
        cat("[PPI] PPI not available for Stage-1 (no barren observations in current training data).\n")
      }
    } else {
      # PPI exists in the dataframe: ensure at least some finite PPI values are present
      # in the barren and vegetation training subsets. If not, compute PPI using
      # an empirical soil baseline from the joined data (preferentially) or the
      # MESMA_DVI_SOIL override.
      barren_global_idx <- !is.na(df$Veg) & tolower(trimws(as.character(df$Veg))) == "barren"
      veg_global_idx <- !is.na(df$`no.soil`) & is.finite(as.numeric(as.character(df$`no.soil`))) & abs(as.numeric(as.character(df$`no.soil`)) - 1) < 0.01
      need_fill <- FALSE
      if (any(barren_global_idx) && sum(is.finite(df$PPI[barren_global_idx])) == 0) need_fill <- TRUE
      if (any(veg_global_idx) && sum(is.finite(df$PPI[veg_global_idx])) == 0) need_fill <- TRUE
      if (need_fill) {
        # Try to compute a baseline from joined barren observations (across df)
        dvi_soil_from_join <- NA_real_
        if (any(barren_global_idx)) dvi_soil_from_join <- mean(df$DVI[barren_global_idx], na.rm = TRUE)
        user_dvi <- suppressWarnings(as.numeric(ifelse(nzchar(Sys.getenv("MESMA_DVI_SOIL")), Sys.getenv("MESMA_DVI_SOIL"), NA)))
        if (is.finite(dvi_soil_from_join) && !is.na(dvi_soil_from_join)) {
            cat(sprintf("[PPI] Filling missing PPI values using joined barren baseline: dvi_soil=%.6f\n", dvi_soil_from_join))
            ## --- NEW: Save this value globally for inference use ---
            assign("GLOBAL_TRAINING_DVI_SOIL", dvi_soil_from_join, envir = globalenv())
            cat(sprintf("[PPI] Saved dynamic baseline (%.6f) for inference step.\n", dvi_soil_from_join))
            ## -----------------------------------------------------
          filled <- add_ppi_columns(df, dvi_soil = dvi_soil_from_join)
        } else if (!is.na(user_dvi) && is.finite(user_dvi)) {
          cat(sprintf("[PPI] Filling missing PPI values using MESMA_DVI_SOIL override: %.6f\n", user_dvi))
          filled <- add_ppi_columns(df, dvi_soil = user_dvi)
        } else {
          cat("[PPI] Could not fill missing PPI values: no barren baseline available and no MESMA_DVI_SOIL set.\n")
          filled <- NULL
        }
        if (!is.null(filled)) {
          # Only fill NA PPI entries to avoid overwriting existing valid values
          fill_idx <- is.na(df$PPI) & is.finite(filled$PPI)
          if (any(fill_idx)) {
            df$PPI[fill_idx] <- filled$PPI[fill_idx]
            cat(sprintf("[PPI] Filled %d PPI NA values from computed baseline\n", sum(fill_idx)))
          }
          if ("PPI" %in% names(df)) avail <- unique(c(avail, "PPI"))
        }
      }
    }
  }
  STAGE1_LIB <- build_barren_veg_library(df, avail, min_samples = 5)
}, silent = FALSE)  # Changed to FALSE to see the actual error

if (is.null(STAGE1_LIB)) {
  stop("[Stage1] Could not build barren-veg library (insufficient 'no soil' == 0 or 'no soil' == 1 rows). Two-stage MESMA requires sufficient stage-1 training data. Cannot proceed with spectral-only unmixing.")
} else {
  cat("[Stage1] Barren-vegetation library built successfully for stage 1 unmixing\n")
  
  cat("\n=== STAGE 1 LIBRARY CHECK ===\n")
  
  barren_lib <- STAGE1_LIB$barren
  barren_vals <- sapply(barren_lib, function(x) mean(x$mu, na.rm = TRUE))
  cat(sprintf("Barren: mean across indices = %.4f (range: [%.4f, %.4f])\n",
              mean(barren_vals), min(barren_vals), max(barren_vals)))
  
  veg_lib <- STAGE1_LIB$vegetation
  veg_vals <- sapply(veg_lib, function(x) mean(x$mu, na.rm = TRUE))
  cat(sprintf("Vegetation: mean across indices = %.4f (range: [%.4f, %.4f])\n",
              mean(veg_vals), min(veg_vals), max(veg_vals)))
  
  cat("\nPer-index separation (ratio-based; dist = ratio - 1):\n")
  for (idx in names(barren_lib)) {
    if (idx %in% names(veg_lib)) {
      b_mean <- mean(barren_lib[[idx]]$mu, na.rm = TRUE)
      v_mean <- mean(veg_lib[[idx]]$mu, na.rm = TRUE)
      if (!is.finite(b_mean) || !is.finite(v_mean)) {
        cat(sprintf("  %s: missing values (skipping)\n", idx))
        next
      }
      if (abs(b_mean) < .Machine$double.eps && abs(v_mean) < .Machine$double.eps) {
        ratio <- 1
      } else {
        ratio <- max(abs(v_mean), abs(b_mean)) / (min(abs(v_mean), abs(b_mean)) + .Machine$double.eps)
      }
      dist <- ratio - 1
      cat(sprintf("  %s: barren=%.4f, veg=%.4f, ratio=%.3f, dist=%.3f\n", idx, b_mean, v_mean, ratio, dist))
      if (dist < 0.1) {
        cat(sprintf("    [NOTICE] Index %s has low relative barren-veg separation (dist=%.3f < 0.1)\n", idx, dist))
      }
    }
  }
  cat("=============================\n\n")
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
vegs <- vegs[tolower(vegs) %in% c("phragmites", "populus", "tamarix")]  # FIXED: case-insensitive matching

if (ENABLE_SAMPLE_BALANCING) {
  lib_df <- lib_df[!is.na(lib_df$Veg), ]
  cat("Applying aggressive sample size balancing to training data...\n")
  lib_df$doy_for_lib <- pheno_doy(lib_df$date)  # Use phenological DOY

  veg_sample_sizes <- sapply(vegs, function(v) {
    dveg <- lib_df[lib_df$Veg == v & is.finite(lib_df$doy_for_lib), , drop = FALSE]
    nrow(dveg)
  })

  min_samples <- min(veg_sample_sizes)
  max_samples <- max(veg_sample_sizes)
  median_samples <- median(veg_sample_sizes)

  cat(sprintf(
    "Original sample sizes - Min: %.0f, Max: %.0f, Median: %.1f\n",
    min_samples, max_samples, median_samples
  ))

  final_balancing_size <- 2 * min(veg_sample_sizes)
  if (exists('TESTING_MAX_PER_VEG') && isTRUE(TESTING_MODE)) {
    final_balancing_size <- min(as.integer(final_balancing_size), as.integer(TESTING_MAX_PER_VEG))
  }

  cat(sprintf("Final balancing size (2x minimum original sample size): %d\n", final_balancing_size))


balanced_dfs <- list()
for (v in vegs) {
  dveg <- lib_df[lib_df$Veg == v & is.finite(lib_df$doy_for_lib), , drop = FALSE]
  if (nrow(dveg) == 0) next
  if (isTRUE(TESTING_MODE)) {
    balanced_dfs[[v]] <- dveg
  } else {
    orig_n <- nrow(dveg)
    target_for_class <- final_balancing_size

    if (target_for_class > orig_n) {
      augmented <- tryCatch(
        {
          augment_minority_class(dveg, target_for_class)
        },
        error = function(e) stop(sprintf("augment_minority_class failed for %s: %s", v, e$message))
      )
      balanced_dfs[[v]] <- augmented
      cat(sprintf("Augmented %s from %d to %d samples\n", v, orig_n, nrow(augmented)))
    } else {
      if (orig_n > target_for_class) {
        set.seed(123)
        dveg <- dveg[sample(1:nrow(dveg), target_for_class), , drop = FALSE]
        cat(sprintf("Downsampled %s from %d to %d samples\n", v, orig_n, target_for_class))
      }
      balanced_dfs[[v]] <- dveg
      cat(sprintf("Using %s with %d samples\n", v, nrow(dveg)))
    }
  }
}

if (length(balanced_dfs) > 0) {
  lib_df <- chunked_rbind(balanced_dfs, chunk_size = 50L)
  if (is.matrix(lib_df)) {
    lib_df <- as.data.frame(lib_df, stringsAsFactors = FALSE)
  }
  gc()
  cat(sprintf(
    "Training data balanced: %d total samples across %d vegetation types (exactly equal)\n",
    nrow(lib_df), length(balanced_dfs)
  ))
}

}  # end ENABLE_SAMPLE_BALANCING


lib <- list()
for (v in vegs) {
  lib[[v]] <- list(n_samples = 0)
}

cat("Skipping per-DOY averaging (will build variants from raw traces)...\n")

timing_info$lib_construction_done <- Sys.time()


cat("Skipping first (old) training pipeline - using second pipeline instead...\n")

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
        center <- colMeans(sub)
        dists <- rowSums(sweep(sub, 2, center, "-")^2)
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



COMPRESSED_STAGE1_LIB <- NULL


# Removed legacy helper: safe_diff() — not used anywhere in the codebase.


# Location-based bootstrap for nested two-stage MESMA with variant switching
# Resamples location-years to compute uncertainty CIs
location_bootstrap_nested_mesma <- function(all_results, comp_templates, compressed_stage1_lib, B = 100, seed = 123) {
  set.seed(seed)

  if (length(all_results) < 2) return(NULL)

  years <- sapply(all_results, function(x) x$year)
  unique_years <- unique(years)

  if (length(unique_years) < 2) return(NULL)

  # Bootstrap resamples of years
  year_results_list <- lapply(1:B, function(b) {
    sampled_years <- sample(unique_years, replace = TRUE)
    sampled_results <- all_results[match(sampled_years, years)]

    # For each sampled result, use the coefs to set mean_med
    res_list <- lapply(sampled_results, function(res) {
      mean_med <- res$coefs
      list(mean_med = mean_med)
    })

    # Average the mean_med across the resample
    veg_names <- unique(unlist(lapply(res_list, function(x) names(x$mean_med))))
    if (length(veg_names) == 0) return(NULL)

    mean_med_avg <- sapply(veg_names, function(v) mean(sapply(res_list, function(x) x$mean_med[[v]]), na.rm = TRUE))
    list(mean_med = setNames(mean_med_avg, veg_names))
  })

  # Filter out NULL
  year_results_list <- year_results_list[!sapply(year_results_list, is.null)]

  if (length(year_results_list) == 0) return(NULL)

  # Now compute CIs
  veg_names <- unique(unlist(lapply(year_results_list, function(x) names(x$mean_med))))
  if (length(veg_names) == 0) return(NULL)

  mean_med <- sapply(veg_names, function(v) mean(sapply(year_results_list, function(x) x$mean_med[[v]]), na.rm = TRUE))
  lower <- sapply(veg_names, function(v) quantile(sapply(year_results_list, function(x) x$mean_med[[v]]), 0.025, na.rm = TRUE))
  upper <- sapply(veg_names, function(v) quantile(sapply(year_results_list, function(x) x$mean_med[[v]]), 0.975, na.rm = TRUE))

  ci_results <- data.frame(Veg = veg_names, coef_025 = lower, coef_975 = upper)
  ci_results
}

# Location-based bootstrap for PPI-normalized aggregation
location_bootstrap_ppi <- function(all_coefs, df_tasks, B = BOOTSTRAP_B, seed = 123) {
  set.seed(seed)
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    warning("dplyr required for PPI bootstrap")
    return(NULL)
  }

  if (!"PPI" %in% names(df_tasks)) {
    warning("PPI column missing from df_tasks")
    return(NULL)
  }

  # Compute median PPI per location-year
  ppi_per_loc_year <- df_tasks |>
    dplyr::group_by(location_id, pheno_year) |>
    dplyr::summarize(median_ppi = median(PPI, na.rm = TRUE), .groups = "drop")

  merged <- dplyr::left_join(all_coefs, ppi_per_loc_year, by = c("location_id", "pheno_year"))
  if (nrow(merged) == 0) return(NULL)

  merged <- merged[!is.na(merged$median_ppi), ]
  veg_types <- unique(merged$Veg[!is.na(merged$Veg)])
  veg_types <- veg_types[!tolower(trimws(veg_types)) %in% c("barren")]

  if (length(veg_types) == 0) return(NULL)

  years <- sort(unique(merged$pheno_year[!is.na(merged$pheno_year)]))
  locations <- unique(merged$location_id)
  n_locs <- length(locations)

  results_list <- list()

  for (veg in veg_types) {
    veg_data <- merged[merged$Veg == veg & !is.na(merged$coef) & !is.na(merged$median_ppi), ]
    if (nrow(veg_data) == 0) next

    boot_sums <- matrix(NA_real_, nrow = B, ncol = length(years))
    colnames(boot_sums) <- as.character(years)

    for (b in seq_len(B)) {
      # Resample locations with replacement
      boot_locs <- sample(locations, n_locs, replace = TRUE)
      boot_data <- veg_data[veg_data$location_id %in% boot_locs, ]

      for (i in seq_along(years)) {
        yr <- years[i]
        yr_data <- boot_data[boot_data$pheno_year == yr, ]
        if (nrow(yr_data) > 0) {
          vals <- yr_data$coef * yr_data$median_ppi
          boot_sums[b, i] <- sum(vals, na.rm = TRUE)
        }
      }
    }

    boot_result <- data.frame(
      year = years,
      Veg = veg,
      n_locations = sapply(years, function(y) sum(veg_data$pheno_year == y & !is.na(veg_data$coef))),
      global_coef = apply(boot_sums, 2, median, na.rm = TRUE),
      se = apply(boot_sums, 2, sd, na.rm = TRUE),
      coef_025 = apply(boot_sums, 2, quantile, 0.025, na.rm = TRUE),
      coef_975 = apply(boot_sums, 2, quantile, 0.975, na.rm = TRUE),
      method = "location_bootstrap_ppi"
    )

    results_list[[veg]] <- boot_result
  }

  dplyr::bind_rows(results_list)
}

# Location-based bootstrap for global aggregation
location_bootstrap_aggregate <- function(all_coefs, B = BOOTSTRAP_B, seed = 123) {
  set.seed(seed)

  veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
  years <- sort(unique(all_coefs$pheno_year[!is.na(all_coefs$pheno_year)]))

  if (!"location_id" %in% names(all_coefs)) {
    warning("location_id column not found in all_coefs")
    return(NULL)
  }

  locations <- unique(all_coefs$location_id)
  n_locs <- length(locations)

  if (n_locs == 0 || length(years) == 0) {
    warning("No valid locations or years found for bootstrapping")
    return(NULL)
  }

  results_list <- list()

  for (veg in veg_types) {
    veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), ]

    boot_means <- matrix(NA_real_, nrow = B, ncol = length(years))
    colnames(boot_means) <- as.character(years)

    for (b in seq_len(B)) {
      # Resample locations with replacement
      boot_locs <- sample(locations, n_locs, replace = TRUE)
      boot_data <- veg_data[veg_data$location_id %in% boot_locs, ]

      for (i in seq_along(years)) {
        yr <- years[i]
        yr_data <- boot_data[boot_data$pheno_year == yr, ]
        if (nrow(yr_data) > 0) {
          boot_means[b, i] <- mean(yr_data$coef, na.rm = TRUE)
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
      method = "location_bootstrap"
    )

    boot_result$coef_025 <- pmax(0, boot_result$coef_025)
    boot_result$coef_975 <- pmin(1, boot_result$coef_975)

    results_list[[veg]] <- boot_result
  }

  dplyr::bind_rows(results_list)
}


reduce_all_traces <- function(lib_df, veg_types, avail_idx, fixed_grid_size = TEMPORAL_BUDGET, 
                               enable_multiscale = FALSE, multiscale_windows = NULL) {
  cat("Performing full-resolution processing on all traces using raw indices...\n")
  
  fixed_grid <- seq_len(N_TEMPORAL_BINS)  # 1:73 for 5-day pentads
  
  reduced_data <- list()
  
  if (is.matrix(lib_df)) lib_df <- as.data.frame(lib_df, stringsAsFactors = FALSE)
  veg_col_name <- NULL
  if ("Veg" %in% names(lib_df)) veg_col_name <- "Veg"
  else if ("veg" %in% names(lib_df)) veg_col_name <- "veg"

  if (is.null(veg_col_name)) {
    cat("[WARN] reduce_all_traces: training data has no vegetation label column ('Veg' or similar). No per-veg reduction will be performed.\n")
    if (is.data.frame(lib_df) || is.matrix(lib_df)) {
    }
    return(list())
  }

  for (veg in veg_types) {
    veg_data <- lib_df[tolower(as.character(lib_df[[veg_col_name]])) == tolower(as.character(veg)), , drop = FALSE]
    
    if (!"location_id" %in% names(veg_data) || !"pheno_year" %in% names(veg_data)) {
      cat(sprintf("  [%s] Missing location_id or pheno_year, skipping\n", veg))
      next
    }

    traces <- unique(veg_data[, c("location_id", "pheno_year")])
    cat(sprintf("  [%s] Reducing %d traces to %d time points...\n", veg, nrow(traces), length(fixed_grid)))

    feature_list <- list()
    Z_list <- list()
    trace_info <- list()

    for (i in seq_len(nrow(traces))) {
      loc <- traces$location_id[i]
      yr <- traces$pheno_year[i]

      dly_year <- veg_data[veg_data$location_id == loc & veg_data$pheno_year == yr, , drop = FALSE]
      n_unique_doys <- length(unique(dly_year$doy))
      if (n_unique_doys < 5) next
      
      idxs <- avail_idx
      K <- length(idxs)
      raw_mat <- matrix(NA_real_, nrow = N_TEMPORAL_BINS, ncol = K)
      colnames(raw_mat) <- idxs
      
      dly_year$pentad <- doy_to_pentad(dly_year$doy)

      for (j in seq_along(idxs)) {
        idn <- idxs[j]
        if (!idn %in% names(dly_year)) next
        
        vals_by_pentad <- tapply(dly_year[[idn]], dly_year$pentad, function(v) {
          v <- v[is.finite(v)]
          if (length(v) == 0) NA_real_ else median(v)
        })
        
        if (!is.null(vals_by_pentad) && length(vals_by_pentad) > 0) {
          pentad_values <- as.integer(names(vals_by_pentad))
          valid_pentad <- pentad_values >= 1 & pentad_values <= N_TEMPORAL_BINS
          raw_mat[pentad_values[valid_pentad], j] <- vals_by_pentad[valid_pentad]
        }
      }

      feat <- as.numeric(raw_mat)

      feat[!is.finite(feat)] <- 0
      if (sum(as.numeric(raw_mat) == 0 | !is.finite(raw_mat)) > 0.5 * length(raw_mat)) next
      
      if (is.numeric(feat) && length(feat) == length(idxs) * length(fixed_grid)) {
        grid_count <- length(fixed_grid)
        nm <- unlist(lapply(idxs, function(x) paste0(x, "_t", seq_len(grid_count))))
        names(feat) <- nm
      }
      feature_list[[length(feature_list) + 1]] <- feat
      Z_list[[length(Z_list) + 1]] <- raw_mat  # Store raw_mat for discriminative templates
      trace_info[[length(trace_info) + 1]] <- list(
        location_id = loc,
        pheno_year = yr,
        trace_index = i
      )
    }
    
    if (length(feature_list) > 0) {
      reduced_data[[veg]] <- list(
        features = do.call(rbind, feature_list),  # Matrix: n_traces x n_features
        Z_matrices = Z_list,                       # List of N_TEMPORAL_BINS x K RAW INDEX matrices
        trace_info = trace_info,                   # Metadata for each trace
        n_samples = nrow(veg_data)                 # Total samples for this veg
      )
      cat(sprintf("  [%s] Reduced %d traces successfully\n", veg, length(feature_list)))
      feat_mat <- reduced_data[[veg]]$features
      grid_count <- length(fixed_grid)
      expected_cols <- grid_count * length(idxs)
      if (!is.null(dim(feat_mat)) && ncol(feat_mat) != expected_cols) {
        warning(sprintf("reduce_all_traces: feature column count mismatch for '%s' (got %d, expected %d).", veg, ncol(feat_mat), expected_cols))
      }
      keys <- sapply(reduced_data[[veg]]$trace_info, function(x) paste0(as.character(x$location_id), "__", as.character(x$pheno_year)))
      if (any(duplicated(keys))) {
        warning(sprintf("reduce_all_traces: duplicate location-pheno_year traces found for '%s' (%d duplicates) — expected exactly one row per loc/pheno_year.", veg, sum(duplicated(keys))))
      }
    }
  }
  
  cat(sprintf("Full-resolution processing complete: %d vegetation types processed\n", length(reduced_data)))
  small_vegs <- names(reduced_data)[sapply(reduced_data, function(x) is.null(x$n_samples) || x$n_samples < MIN_ENDMEMBER_SAMPLES)]
  if (length(small_vegs) > 0) {
    cat(sprintf("[NOTICE] Removing %d vegetation types from reduced_data due to insufficient samples: %s\n", length(small_vegs), paste(small_vegs, collapse=", ")))
    reduced_data[small_vegs] <- NULL
  }
  return(reduced_data)
}

analyze_library_similarity <- function(mesma_lib, compressed_templates_accessor, grid_type = "full") {
  cat("\n=== INTER-CLASS VARIANT SIMILARITY ANALYSIS ===\n")
  
  all_variants <- list()
  variant_names <- c()
  variant_vegs <- c()
  
  for (veg in names(mesma_lib)) {
    for (variant in mesma_lib[[veg]]) {
      vid <- variant$variant_id
      vec <- compressed_templates_accessor[[veg]][[vid]][["full"]]
      if (!is.null(vec)) {
        all_variants[[length(all_variants) + 1]] <- list(vec = vec, variant = variant, veg = veg)
        variant_names <- c(variant_names, vid)
        variant_vegs <- c(variant_vegs, veg)
      }
    }
  }
  
  if (length(all_variants) < 2) {
    cat("Not enough variants for similarity analysis.\n")
    return(NULL)
  }
  
  cat("\n=== SMALL VARIANTS (potential outliers) ===\n")
  for (i in seq_along(all_variants)) {
    var <- all_variants[[i]]
    if (var$variant$n_samples < 10) {
      cat(sprintf("Variant %s (%s): %d samples\n", var$variant$variant_id, var$veg, var$variant$n_samples))
      if (!is.null(var$variant$sample_ids)) {
        cat(sprintf("  Samples: %s\n", paste(var$variant$sample_ids, collapse = ", ")))
      }
    }
  }
  cat("=======================================\n\n")
  
  n_vars <- length(all_variants)
  
  high_sim_pairs <- data.frame(
    Variant1 = character(),
    Veg1 = character(),
    Variant2 = character(),
    Veg2 = character(),
    Similarity = numeric(),
    threshold_tier = character(),
    euclidean_dist = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:(n_vars - 1)) {
    for (j in (i + 1):n_vars) {
      sim <- cos_sim(all_variants[[i]]$vec, all_variants[[j]]$vec)
      
      if (sim == 0) {
        cat(sprintf("DEBUG: Zero similarity between %s (%s) and %s (%s)\n", 
                    variant_names[i], variant_vegs[i], variant_names[j], variant_vegs[j]))
        
        vec1 <- all_variants[[i]]$vec
        vec2 <- all_variants[[j]]$vec
        vec1_norm <- sqrt(sum(vec1^2))
        vec2_norm <- sqrt(sum(vec2^2))
        vec1_zero <- all(vec1 == 0)
        vec2_zero <- all(vec2 == 0)
        
        cat(sprintf("  Vec1 (%s): norm=%.6f, all_zero=%s, range=[%.6f, %.6f]\n", 
                    variant_names[i], vec1_norm, vec1_zero, min(vec1), max(vec1)))
        cat(sprintf("  Vec2 (%s): norm=%.6f, all_zero=%s, range=[%.6f, %.6f]\n", 
                    variant_names[j], vec2_norm, vec2_zero, min(vec2), max(vec2)))
      }
      
      dist <- sqrt(sum((all_variants[[i]]$vec - all_variants[[j]]$vec)^2))
      
      if (variant_vegs[i] != variant_vegs[j]) {
        tier <- NA_character_
        if (sim > 0.98) tier <- "CRITICAL (>0.98)"
        else if (sim > 0.95) tier <- "HIGH (>0.95)"
        else if (sim > 0.90) tier <- "WARN (>0.90)"
        
        if (!is.na(tier)) {
           high_sim_pairs[nrow(high_sim_pairs) + 1, ] <- list(
             variant_names[i], variant_vegs[i],
             variant_names[j], variant_vegs[j],
             sim, tier, dist
           )
        }
      }
    }
  }
  
  if (nrow(high_sim_pairs) > 0) {
    high_sim_pairs <- high_sim_pairs[order(-high_sim_pairs$Similarity), ]
    cat(sprintf("Found %d pairs of variants from DIFFERENT vegetation types with similarity > 0.90:\n", nrow(high_sim_pairs)))
    print(head(high_sim_pairs, 20))
    if (nrow(high_sim_pairs) > 20) cat(sprintf("... and %d more.\n", nrow(high_sim_pairs) - 20))
    
    assign("INSEPARABLE_VARIANT_INFO", list(pairs = high_sim_pairs), envir = globalenv())
    
    inseparable_vars <- list()
    for (i in 1:nrow(high_sim_pairs)) {
       if (high_sim_pairs$Similarity[i] > 0.95) {
         v1 <- high_sim_pairs$Veg1[i]
         v2 <- high_sim_pairs$Veg2[i]
         vid1 <- high_sim_pairs$Variant1[i]
         vid2 <- high_sim_pairs$Variant2[i]
         
         if (is.null(inseparable_vars[[v1]])) inseparable_vars[[v1]] <- c()
         inseparable_vars[[v1]] <- unique(c(inseparable_vars[[v1]], vid1))
         
         if (is.null(inseparable_vars[[v2]])) inseparable_vars[[v2]] <- c()
         inseparable_vars[[v2]] <- unique(c(inseparable_vars[[v2]], vid2))
       }
    }
    assign("INSEPARABLE_VARIANTS", inseparable_vars, envir = globalenv())
    
  } else {
    cat("No high similarity (> 0.90) pairs found between different vegetation types.\n")
    assign("INSEPARABLE_VARIANTS", NULL, envir = globalenv())
  }

  cat("Generating full variant similarity heatmap...\n")

  all_vecs <- list()
  all_ids <- c()
  all_vegs <- c()

  for (veg in names(mesma_lib)) {
    for (variant in mesma_lib[[veg]]) {
      vid <- variant$variant_id
      vec <- NULL
      if (!is.null(compressed_templates_accessor[[veg]]) &&
          !is.null(compressed_templates_accessor[[veg]][[vid]]) &&
          !is.null(compressed_templates_accessor[[veg]][[vid]][["full"]])) {
         vec <- compressed_templates_accessor[[veg]][[vid]][["full"]]
      }

      if (!is.null(vec)) {
         all_vecs[[length(all_vecs)+1]] <- vec
         all_ids <- c(all_ids, vid)
         all_vegs <- c(all_vegs, veg)
      }
    }
  }

  if (length(all_vecs) > 1) {
     n_v <- length(all_vecs)
     sim_mat <- matrix(NA, n_v, n_v)
     rownames(sim_mat) <- all_ids
     colnames(sim_mat) <- all_ids

     for (i in 1:n_v) {
       for (j in 1:n_v) {
          sim_mat[i,j] <- pmax(pmin(cos_sim(all_vecs[[i]], all_vecs[[j]]), 1.01), -0.01)
       }
     }

     sim_df <- as.data.frame(as.table(sim_mat))
     colnames(sim_df) <- c("Var1", "Var2", "Similarity")

     ord_idx <- order(all_vegs, all_ids)
     ordered_ids <- all_ids[ord_idx]
     ordered_vegs <- all_vegs[ord_idx]

     sim_df$Var1 <- factor(sim_df$Var1, levels = ordered_ids)
     sim_df$Var2 <- factor(sim_df$Var2, levels = rev(ordered_ids)) # Reverse for y-axis to match matrix layout

     veg_changes <- which(diff(as.numeric(factor(ordered_vegs))) != 0)
     if (length(veg_changes) > 0) {
       vline_positions <- veg_changes + 0.5
       hline_positions <- length(ordered_ids) - veg_changes + 0.5
     } else {
       vline_positions <- numeric(0)
       hline_positions <- numeric(0)
     }

     p_heat <- ggplot(sim_df, aes(Var1, Var2, fill = Similarity)) +
       geom_tile() +
       scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.9, limits = c(-0.01, 1.01), na.value = "white") +
       theme_minimal() +
       theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
             axis.text.y = element_text(size = 6)) +
       labs(title = "Pairwise Cosine Similarity of All Variants", x = NULL, y = NULL)

     if (length(vline_positions) > 0) {
       p_heat <- p_heat +
         geom_vline(xintercept = vline_positions, color = "black", linewidth = 0.5) +
         geom_hline(yintercept = hline_positions, color = "black", linewidth = 0.5)
     }

     ggsave(file.path(OUT_DIR, "variant_similarity_heatmap.png"), p_heat, width = 12, height = 10)
     cat(sprintf("Saved similarity heatmap to: %s\n", file.path(OUT_DIR, "variant_similarity_heatmap.png")))

     write.csv(sim_mat, file.path(OUT_DIR, "variant_similarity_matrix.csv"))
  }

  cat("==============================================\n\n")
}

reduce_all_traces_simple <- function(lib_df, veg_types, avail_idx) {
  
  reduced_data <- list()
  
  for (veg in veg_types) {
    veg_data <- lib_df[tolower(lib_df$Veg) == tolower(veg), ]
    if (nrow(veg_data) == 0) next
    
    traces <- unique(veg_data[, c("location_id", "pheno_year")])

    feature_list <- list()
    Z_list <- list()
    trace_info <- list()

    for (i in seq_len(nrow(traces))) {
      loc <- traces$location_id[i]
      yr <- traces$pheno_year[i]

      dly_year <- veg_data[veg_data$location_id == loc & veg_data$pheno_year == yr, ]
      n_unique_doys <- length(unique(dly_year$doy))
      if (n_unique_doys < 5) next
      
      raw_mat <- build_pentad_matrix(dly_year, avail_idx)
      if (is.null(raw_mat)) next
      
      n_unique_doys <- length(unique(dly_year$doy))
      if (n_unique_doys < 5) next
      
      full_features <- as.numeric(raw_mat)  # 365 * K vector
      full_features[!is.finite(full_features)] <- 0
      if (sum(as.numeric(raw_mat) == 0 | !is.finite(raw_mat)) > 0.5 * length(raw_mat)) next
      
      feature_list[[length(feature_list) + 1]] <- full_features
      Z_list[[length(Z_list) + 1]] <- raw_mat
      trace_info[[length(trace_info) + 1]] <- list(
        location_id = loc, pheno_year = yr, trace_index = i
      )
    }
    
    if (length(feature_list) > 0) {
      reduced_data[[veg]] <- list(
        features = do.call(rbind, feature_list),
        Z_matrices = Z_list,
        trace_info = trace_info,
        n_samples = nrow(veg_data)
      )
      cat(sprintf("  [%s] Processed %d traces\n", veg, length(feature_list)))
    }
  }
  
  reduced_data
}


build_mesma_variants_weighted <- function(reduced_data, raw_lib_templates, 
                                           pca_lda_weights = NULL,
                                           n_variants = 10,
                                           min_cluster_size = 10,
                                           avail_idx = NULL,
                                           norm_params = NULL) {
  mesma_lib <- list()
  
  for (veg in names(reduced_data)) {
    veg_info <- reduced_data[[veg]]
    if (!is.null(veg_info$n_samples) && veg_info$n_samples < MIN_ENDMEMBER_SAMPLES) {
      cat(sprintf("[NOTICE] Skipping veg '%s' in variant construction due to insufficient samples: %d < %d\n", veg, veg_info$n_samples, MIN_ENDMEMBER_SAMPLES))
      next
    }
    if (!is.null(veg_info$n_samples) && veg_info$n_samples < MIN_ENDMEMBER_SAMPLES) {
      cat(sprintf("[NOTICE] Skipping veg '%s' in variant construction due to insufficient samples: %d < %d\n", veg, veg_info$n_samples, MIN_ENDMEMBER_SAMPLES))
      next
    }
    X_feat <- veg_info$features
    Z_list <- veg_info$Z_matrices
    
    if (nrow(X_feat) < min_cluster_size) {
      mesma_lib[[veg]] <- list(list(
        raw_mat = raw_lib_templates[[veg]]$T,
        variant_id = paste0(veg, "_single"),
        n_samples = veg_info$n_samples,
        sample_ids = sapply(veg_info$trace_info, function(x) paste0(as.character(x$location_id), "__", as.character(x$pheno_year)))
      ))
      next
    }
    
    if (!is.null(Z_list) && length(Z_list) > 0) {
      n_timepoints <- nrow(Z_list[[1]])
      K <- ncol(Z_list[[1]])
    } else {
      n_timepoints <- 365
      K <- as.integer(ncol(X_feat) / n_timepoints)
    }
    
    X_std <- matrix(NA_real_, nrow = nrow(X_feat), ncol = ncol(X_feat))
    
    idx_names <- if (!is.null(avail_idx)) avail_idx else paste0("idx", 1:K)
    
    for (k_idx in seq_len(K)) {
      col_start <- (k_idx - 1) * n_timepoints + 1
      col_end <- k_idx * n_timepoints
      cols <- seq(col_start, col_end)
      if (max(cols) > ncol(X_feat)) cols <- cols[cols <= ncol(X_feat)]
      
      index_data <- X_feat[, cols, drop = FALSE]
      
      idx_name <- idx_names[k_idx]
      mu <- 0
      sigma <- 1
      
      if (!is.null(norm_params) && !is.null(norm_params$INDEX_SCALES) && idx_name %in% names(norm_params$INDEX_SCALES)) {
        params <- norm_params$INDEX_SCALES[[idx_name]]
        mu <- params$mean
        sigma <- params$sd
      } else {
        mu <- mean(index_data, na.rm = TRUE)
        sigma <- sd(as.vector(index_data), na.rm = TRUE)
        if (!is.finite(sigma) || sigma < 1e-10) sigma <- 1.0
      }
      
      X_std[, cols] <- (index_data - mu) / sigma
    }
    X_std[!is.finite(X_std)] <- 0
    
    X_weighted <- X_std
    if (!is.null(pca_lda_weights) && length(pca_lda_weights) == ncol(X_feat)) {
      # Weights are already dampened & normalized in compute_lda_weights; apply directly
      X_weighted <- sweep(X_std, 2, pca_lda_weights, "*")
    }
    
    col_vars <- apply(X_weighted, 2, var, na.rm = TRUE)
    keep_cols <- col_vars > 1e-10
    if (sum(keep_cols) < 10) {
      warning(sprintf("[%s] Too few variable features, using single variant", veg))
      mesma_lib[[veg]] <- list(list(
        raw_mat = raw_lib_templates[[veg]]$T,
        variant_id = paste0(veg, "_single"),
        n_samples = veg_info$n_samples
      ))
      next
    }
    X_weighted <- X_weighted[, keep_cols, drop = FALSE]
    
    X_std <- X_weighted
    
    if (ncol(X_std) > 50) {
      pca_clust <- prcomp(X_std, center = FALSE, scale. = FALSE)
      cum_var <- cumsum(pca_clust$sdev^2) / sum(pca_clust$sdev^2)
      n_pcs <- min(50, which(cum_var >= 0.95)[1])
      if (is.na(n_pcs)) n_pcs <- min(50, ncol(pca_clust$x))
      X_clust <- pca_clust$x[, 1:n_pcs, drop = FALSE]
      cat(sprintf("  [%s] Clustering in %d-PC space\n", veg, n_pcs))
    } else {
      X_clust <- X_std
    }
    
    k <- min(n_variants, nrow(X_clust))
    k <- max(k, 1)

    km <- kmeans(X_clust, centers = k, nstart = 25, iter.max = 100)
    
    variants <- list()
    for (clust in seq_len(k)) {
      members <- which(km$cluster == clust)
      if (length(members) == 0) next
      
      X_members <- X_feat[members, , drop = FALSE]
      center <- colMeans(X_members)
      dists <- rowSums(sweep(X_members, 2, center, "-")^2)
      medoid_idx <- members[which.min(dists)]
      
      raw_mat <- Z_list[[medoid_idx]]
      
      raw_vec <- as.numeric(raw_mat)
      raw_vec[is.na(raw_vec)] <- 0
      raw_norm <- sqrt(sum(raw_vec^2))
      if (raw_norm < 1e-10) {
        warning(sprintf("Skipping zero norm variant %s_v%d (norm=%.2e)", veg, clust, raw_norm))
        next
      }
      raw_sd <- sd(raw_vec, na.rm = TRUE)
      raw_range <- diff(range(raw_vec, na.rm = TRUE))
      raw_unique <- length(unique(raw_vec[is.finite(raw_vec)]))
      
      if (!is.finite(raw_sd) || raw_sd < 1e-10) {
        warning(sprintf("Skipping degenerate variant %s_v%d (SD=%.2e)", veg, clust, raw_sd))
        next
      }
      
      if (!is.finite(raw_range) || raw_range < 1e-10) {
        warning(sprintf("Skipping degenerate variant %s_v%d (range=%.2e)", veg, clust, raw_range))
        next
      }
      
      if (raw_unique < 3) {
        warning(sprintf("Skipping degenerate variant %s_v%d (only %d unique values)", veg, clust, raw_unique))
        next
      }
      
      variants[[length(variants) + 1]] <- list(
        raw_mat = raw_mat,
        variant_id = paste0(veg, "_v", clust),
        n_samples = length(members),
        cluster_center = km$centers[clust, ],
        sample_ids = sapply(veg_info$trace_info[members], function(x) paste0(as.character(x$location_id), "__", as.character(x$pheno_year)))
      )
    }

    if (length(variants) >= 3) {
      var_vecs <- do.call(rbind, lapply(variants, function(v) {
        vec <- as.numeric(v$raw_mat)
        vec[is.na(vec)] <- 0
        vec
      }))
      n_v <- nrow(var_vecs)
      
      sim_mat <- matrix(1, n_v, n_v)
      for (i in 1:(n_v-1)) {
        for (j in (i+1):n_v) {
          dot_prod <- sum(var_vecs[i,] * var_vecs[j,])
          norm_i <- sqrt(sum(var_vecs[i,]^2))
          norm_j <- sqrt(sum(var_vecs[j,]^2))
          
          if (norm_i > 0 && norm_j > 0) {
            s <- dot_prod / (norm_i * norm_j)
          } else {
            s <- 0
          }
          sim_mat[i,j] <- s
          sim_mat[j,i] <- s
        }
      }
      
      avg_sim <- rowSums(sim_mat - diag(n_v)) / (n_v - 1)
      
      outlier_mask <- rep(FALSE, length(avg_sim))
      
      if (any(outlier_mask)) {
         cat(sprintf("  [%s] Removed %d outlier variants (avg pairwise sim < 0.5)\n", veg, sum(outlier_mask)))
         variants <- variants[!outlier_mask]
      }
    }
    
    removed_variants <- which(sapply(variants, function(v) !is.null(v$n_samples) && v$n_samples < MIN_ENDMEMBER_SAMPLES))
    if (length(removed_variants) > 0) {
      cat(sprintf("  [%s] Removing %d variant(s) with samples < %d\n", veg, length(removed_variants), MIN_ENDMEMBER_SAMPLES))
      variants <- variants[-removed_variants]
    }
    if (length(variants) == 0) {
      if (!is.null(raw_lib_templates[[veg]]) && !is.null(raw_lib_templates[[veg]]$n_samples) && raw_lib_templates[[veg]]$n_samples >= MIN_ENDMEMBER_SAMPLES) {
        cat(sprintf("  [%s] No valid cluster variants remain; using single medoid fallback\n", veg))
        variants <- list(list(raw_mat = raw_lib_templates[[veg]]$T, variant_id = paste0(veg, "_single"), n_samples = raw_lib_templates[[veg]]$n_samples))
      } else {
        cat(sprintf("  [%s] No valid variants and no suitable raw template: skipping veg type\n", veg))
        mesma_lib[[veg]] <- list()  # Insert empty list to indicate skipped
        next
      }
    }

    desired_k <- n_variants
    if (length(variants) > desired_k) {
      ord <- order(sapply(variants, function(v) if (!is.null(v$n_samples)) v$n_samples else 0), decreasing = TRUE)
      variants <- variants[ord[1:desired_k]]
      cat(sprintf("  [%s] Trimmed variants to %d (kept top by n_samples)\n", veg, desired_k))
    } else if (length(variants) < desired_k) {
      if (length(variants) == 0) {
        cat(sprintf("  [%s] No valid variants to augment; skipping augmentation\n", veg))
      } else {
        base_ord <- order(sapply(variants, function(v) if (!is.null(v$n_samples)) v$n_samples else 0), decreasing = TRUE)
        i_dup <- 1
        while (length(variants) < desired_k) {
          src <- variants[[base_ord[((i_dup - 1) %% length(base_ord)) + 1]]]
          dup_id <- paste0(src$variant_id, "_dup", i_dup)
          dup <- list(raw_mat = src$raw_mat, variant_id = dup_id, n_samples = 0, sample_ids = src$sample_ids)
          variants[[length(variants) + 1]] <- dup
          i_dup <- i_dup + 1
        }
        cat(sprintf("  [%s] Augmented variants to %d by duplicating existing prototypes\n", veg, desired_k))
      }
    }

    mesma_lib[[veg]] <- variants
    cat(sprintf("  [%s] Created %d variants\n", veg, length(variants)))
  }
  
  mesma_lib
}

build_spectral_library <- function(lib_df, avail_idx) {
  
  lib <- list()
  veg_types <- ALLOWED_VEG
  
  for (veg in veg_types) {
    veg_data <- lib_df[tolower(lib_df$Veg) == tolower(veg), ]
    if (nrow(veg_data) == 0) {
      cat(sprintf("WARNING: No training data for veg class '%s', skipping.\n", veg))
      next
    }
    
    lib[[veg]] <- list(
      raw_traces = veg_data,
      n_samples = nrow(veg_data)
    )
  }
  
  cat("\n=== COMPUTING GLOBAL PCA-LDA WEIGHTS ===\n")
  
  PCA_LDA_RESULT <- compute_pca_lda_weights(
    lib_df = lib_df,
    avail_idx = avail_idx,
    pca_variance_threshold = 0.95,
    lda_weight_floor = 0.01
  )
  
  if (!is.null(PCA_LDA_RESULT)) {
    LDA_FEATURE_WEIGHTS <- PCA_LDA_RESULT$weights
    cat(sprintf("PCA-LDA weights computed: %d features, %d PCs retained\n", 
                length(LDA_FEATURE_WEIGHTS), PCA_LDA_RESULT$n_pcs))
  } else {
    LDA_FEATURE_WEIGHTS <- NULL
    warning("PCA-LDA weight computation failed, using unweighted features")
  }
  
  assign("LDA_FEATURE_WEIGHTS", LDA_FEATURE_WEIGHTS, envir = globalenv())
  assign("PCA_LDA_RESULT", PCA_LDA_RESULT, envir = globalenv())
  
  cat("\n=== PROCESSING TRACES AT FULL YEAR RESOLUTION ===\n")
  
  reduced_traces <- reduce_all_traces_simple(
    lib_df = lib_df,
    veg_types = names(lib),
    avail_idx = avail_idx
  )
  
  cat("\n=== Rows per class after full-resolution processing ===\n")
  for (v in names(reduced_traces)) {
    cat(sprintf("  %s: %d rows\n", v, nrow(reduced_traces[[v]]$features)))
  }
  cat("=============================================\n\n")
  
  mesma_lib <- build_mesma_variants_weighted(
    reduced_data = reduced_traces, 
    raw_lib_templates = lib,
    pca_lda_weights = LDA_FEATURE_WEIGHTS,
    n_variants = N_VARIANTS_PER_VEG,
    min_cluster_size = MIN_CLUSTER_SIZE,
    avail_idx = avail_idx,
    norm_params = TRAINING_NORM_PARAMS
  )

  precompute_compressed_templates <- function(mesma_lib, grid_name = "full") {
    template_db <- list()
    for (veg in names(mesma_lib)) {
      template_db[[veg]] <- list()
      for (variant in mesma_lib[[veg]]) {
        template_db[[veg]][[variant$variant_id]] <- list()
        compressed_vec <- as.numeric(variant$raw_mat)
        # Preserve NA for missing values so downstream weight computation
        # can compute statistics using na.rm=TRUE instead of treating NA as 0.
        compressed_vec[!is.finite(compressed_vec)] <- NA_real_
        template_db[[veg]][[variant$variant_id]][["full"]] <- compressed_vec
      }
    }
    template_db
  }
  compressed_templates_accessor <- precompute_compressed_templates(mesma_lib, "full")
  assign(".COMPRESSED_TEMPLATES_ACCESSOR", compressed_templates_accessor, envir = globalenv())

  analyze_library_similarity(mesma_lib, compressed_templates_accessor, grid_type = "full")

  LDA_FEATURE_WEIGHTS <- compute_lda_weights(mesma_lib, compressed_templates_accessor, grid_type = "full")
  assign("LDA_FEATURE_WEIGHTS", LDA_FEATURE_WEIGHTS, envir = globalenv())

  mesma_lib
}

build_mesma_variants <- function(reduced_data, raw_lib_templates, min_cluster_size = 10) {
  mesma_lib <- list()

  barren_proto_raw <- NULL
  if ("barren" %in% names(reduced_data)) {
    b_feat <- reduced_data[["barren"]]$features
    if (!is.null(b_feat) && nrow(b_feat) > 0) {
      barren_proto_raw <- apply(b_feat, 2, median, na.rm = TRUE)
      cat(sprintf("[Soil Correction] Identified barren prototype from %d traces\n", nrow(b_feat)))
    }
  }

  for (veg in names(reduced_data)) {
    veg_info <- reduced_data[[veg]]
    
    X_feat <- veg_info$features        # n_traces x n_features (time-reduced)
    Z_list <- veg_info$Z_matrices      # List of 365 x n_indices raw matrices
    n_samples <- veg_info$n_samples
    
    if (nrow(X_feat) < min_cluster_size) {
      mesma_lib[[veg]] <- list(
        list(
          raw_mat = raw_lib_templates[[veg]]$T,
          variant_id = paste0(veg, "_single"),
          n_samples = n_samples
        )
      )
      cat(sprintf("  [%s] Insufficient traces (%d), using single variant\n", veg, nrow(X_feat)))
      next
    }
    
    cat(sprintf("  [%s] Building variants from %d reduced traces...\n", veg, nrow(X_feat)))
    
    if (nrow(X_feat) > MAX_PROJECTIONS_PER_VEG) {
      set.seed(123)
      idx <- sample(nrow(X_feat), MAX_PROJECTIONS_PER_VEG)
      X_feat <- X_feat[idx, , drop = FALSE]
      Z_list <- Z_list[idx]
      cat(sprintf("  [%s] Subsampled to %d traces\n", veg, nrow(X_feat)))
    }
    if (!is.null(veg_info$Xw)) {
      X_w <- as.matrix(veg_info$Xw)
      whitened <- list(Xw = X_w, W = if (!is.null(veg_info$W)) veg_info$W else diag(ncol(X_w)), mu = if (!is.null(veg_info$mu)) veg_info$mu else rep(0, ncol(X_w)))
    } else {
      cat(sprintf("\n[CLUSTERING PREP in build_mesma_variants] %s: %d traces × %d features\n", 
                  veg, nrow(X_feat), ncol(X_feat)))
      
      if (!is.null(Z_list) && length(Z_list) > 0) {
        n_timepoints <- nrow(Z_list[[1]])
        K <- ncol(Z_list[[1]])
      } else {
        n_timepoints <- 365
        K <- as.integer(ncol(X_feat) / n_timepoints)
      }
      X_std <- matrix(NA_real_, nrow = nrow(X_feat), ncol = ncol(X_feat))
      for (k_idx in seq_len(K)) {
        col_start <- (k_idx - 1) * n_timepoints + 1
        col_end <- k_idx * n_timepoints
        cols <- seq(col_start, col_end)
        if (max(cols) > ncol(X_feat)) {
          cols <- cols[cols <= ncol(X_feat)]
        }
        index_data <- X_feat[, cols, drop = FALSE]
        global_mean <- mean(index_data, na.rm = TRUE)
        global_sd <- sd(as.vector(index_data), na.rm = TRUE)
        if (!is.finite(global_sd) || global_sd < 1e-10) global_sd <- 1.0
        X_std[, cols] <- (index_data - global_mean) / global_sd
      }
      X_std[!is.finite(X_std)] <- 0
      
      if (ncol(X_std) > 50 && nrow(X_std) < ncol(X_std) * 3) {
        pca_result <- prcomp(X_std, center = FALSE, scale. = FALSE)
        var_explained <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
        n_pcs <- min(50, which(var_explained > 0.95)[1])
        if (is.na(n_pcs)) n_pcs <- min(50, ncol(pca_result$x))
        X_w <- pca_result$x[, 1:n_pcs, drop = FALSE]
        cat(sprintf("  Reduced to %d PCs (%.1f%% variance) for clustering\n", 
                    n_pcs, 100 * var_explained[n_pcs]))
      } else {
        X_w <- X_std
        cat(sprintf("  Using all %d standardized features for clustering\n", ncol(X_std)))
      }
      
      X_means <- colMeans(X_feat, na.rm = TRUE)
      whitened <- list(Xw = X_w, W = diag(ncol(X_feat)), mu = X_means)
    }
    
    if (is.null(veg_info$Xw) && !is.null(barren_proto_raw) && veg != "barren") {
      if (length(barren_proto_raw) == ncol(X_w)) {
        b_centered <- barren_proto_raw - whitened$mu
        b_proj <- as.numeric(b_centered %*% whitened$W)
        b_norm2 <- sum(b_proj^2)
        if (b_norm2 > 1e-9) {
          alphas <- (X_w %*% b_proj) / b_norm2
          X_w <- X_w - (alphas %*% t(b_proj))
          cat(sprintf("  [%s] Applied soil subtraction (orthogonalization) in whitened space\n", veg))
        }
      }
    }
    
    if (veg == "barren") {
      clust_members <- seq_len(nrow(X_w))
      median_center <- apply(X_w, 2, median, na.rm = TRUE)
      dists <- rowSums(sweep(X_w[clust_members, , drop=FALSE], 2, median_center, "-")^2)
      best_idx_local <- which.min(dists)
      best_idx_global <- clust_members[best_idx_local]
      medoid_Z <- Z_list[[best_idx_global]]
      
      raw_vec <- as.numeric(medoid_Z)
      raw_sd <- sd(raw_vec, na.rm = TRUE)
      raw_range <- diff(range(raw_vec, na.rm = TRUE))
      raw_unique <- length(unique(raw_vec[is.finite(raw_vec)]))
      
      if (!is.finite(raw_sd) || raw_sd < 1e-10 ||
          !is.finite(raw_range) || raw_range < 1e-10 ||
          raw_unique < 3) {
        warning(sprintf("Barren prototype is degenerate, using fallback"))
        medoid_Z <- raw_lib_templates[["barren"]]$T
      }
      
      variants <- list(list(
        raw_mat = medoid_Z,
        variant_id = paste0(veg, "_single"),
        n_samples = length(clust_members),
        cluster_center = median_center
      ))
      mesma_lib[[veg]] <- variants
      next
    }
    
    
    best_k <- 7
    
    if (best_k > nrow(X_w)) {
       best_k <- nrow(X_w)
       cat(sprintf("  [%s] Warning: Requested 10 variants but only %d samples available. Using %d variants.\n", veg, nrow(X_w), best_k))
    }
    
    best_sil <- NA_real_ # Silhouette not calculated
    
    km_final <- kmeans(X_w, centers = best_k, nstart = 25, iter.max = 100)
    
    cat(sprintf("\n[CLUSTERING DEBUG] %s: best_k=%d (silhouette=%.4f)\n", 
                veg, best_k, best_sil))
    
    cluster_sizes <- table(km_final$cluster)
    cat(sprintf("  Cluster sizes: %s\n", 
                paste(sprintf("C%d=%d", 1:best_k, cluster_sizes), collapse=", ")))
    
    bss_tss <- km_final$betweenss / km_final$totss
    cat(sprintf("  Between-cluster variance: %.1f%%\n", 100 * bss_tss))
    
    if (bss_tss < 0.3) {
      warning(sprintf("%s: Clusters are poorly separated (%.1f%% variance explained)", 
                      veg, 100 * bss_tss))
    }
    
    variants <- list()
    for (clust in seq_len(best_k)) {
      clust_members <- which(km_final$cluster == clust)
      if (length(clust_members) < 1) next
      
      dists <- rowSums(sweep(X_w[clust_members, , drop=FALSE], 2, km_final$centers[clust, ], "-")^2)
      best_idx_local <- which.min(dists)
      best_idx_global <- clust_members[best_idx_local]
      
      medoid_Z <- Z_list[[best_idx_global]]
      
      raw_vec <- as.numeric(medoid_Z)
      raw_sd <- sd(raw_vec, na.rm = TRUE)
      raw_range <- diff(range(raw_vec, na.rm = TRUE))
      raw_unique <- length(unique(raw_vec[is.finite(raw_vec)]))
      
      if (!is.finite(raw_sd) || raw_sd < 1e-10) {
        warning(sprintf("Skipping degenerate variant %s_v%d (SD=%.2e)", veg, clust, raw_sd))
        next
      }
      
      if (!is.finite(raw_range) || raw_range < 1e-10) {
        warning(sprintf("Skipping degenerate variant %s_v%d (range=%.2e)", veg, clust, raw_range))
        next
      }
      
      if (raw_unique < 3) {
        warning(sprintf("Skipping degenerate variant %s_v%d (only %d unique values)", veg, clust, raw_unique))
        next
      }
      
      variants[[length(variants) + 1]] <- list(
        raw_mat = medoid_Z,
        variant_id = paste0(veg, "_v", clust),
        n_samples = length(clust_members),
        cluster_center = km_final$centers[clust, ]
      )
    }
    
    mesma_lib[[veg]] <- variants
  }
  
  return(mesma_lib)
}

  # Old solver `solve_weights_ols()` removed — use `solve_weights_nnls_simple(E, y)` which mirrors Stage 1 NNLS behaviour


  # Simple NNLS wrapper that mirrors ols_stage1_unmix behaviour (no extra preprocessing)
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

    # Impute non-finite values to 0 (same as stage1)
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

    if (is.null(res)) return(list(w = rep(1 / ncol(E_fit), ncol(E_fit)), rmse = Inf))

    w <- res$x
    if (sum(w) > 0) w <- w / sum(w)

    pred <- as.numeric(E_fit %*% w)
    rmse <- sqrt(mean((y_fit - pred)^2))

    list(w = w, rmse = rmse)
  }

  spectral_angle <- function(a, b) {
    cs <- sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
    acos(pmax(-1, pmin(1, cs)))
  }

  cos_angle <- function(a, b) {
    sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
  }

  solve_fclsu <- function(E, y, lambda = 1e-6) {
    # Use stage-1 style NNLS solver for FCLSU behaviour (no extra preprocessing)
    result <- solve_weights_nnls_simple(E, y)
    list(w = result$w, residual = result$rmse * sqrt(nrow(E)))
  }

  stage2_ols_unmix <- function(y, vegetation_libraries, topK = 2, feature_weights = NULL) {

    veg_types <- names(vegetation_libraries)
    if (length(veg_types) == 0) return(NULL)
    X <- as.numeric(y)
    # Do NOT re-normalize here - y is z-score normalized from fit_one_task; we'll apply
    # sqrt(feature_weights) to both observations and templates for consistent WLS similarity.

    sqrt_w <- NULL
    if (!is.null(feature_weights) && length(feature_weights) == length(X)) {
      sqrt_w <- sqrt(pmax(feature_weights, 0))
    }

    top_variants <- list()
    for (v in veg_types) {
      cand <- vegetation_libraries[[v]]; if (length(cand) == 0) next
      # Compute similarities using dot product in the weighted space
      sims <- sapply(cand, function(x) {
        vec <- if (is.list(x) && !is.null(x$vec)) x$vec else x
        if (!is.null(sqrt_w) && length(vec) == length(sqrt_w)) {
          sum((vec * sqrt_w) * (X * sqrt_w))
        } else {
          sum(vec * X)
        }
      })
      ord <- order(sims, decreasing = TRUE);
      keep <- ord[seq_len(min(topK, length(ord)))]
      top_variants[[v]] <- lapply(keep, function(i) { x <- cand[[i]]; if (is.list(x) && !is.null(x$vec)) x else list(vec = x, id = paste0(v, "_", i)) })
    }
    if (length(top_variants) == 0) return(NULL)

    result <- evaluate_all_combinations(X, top_variants, lambda = 0, feature_weights = feature_weights)

    if (is.null(result)) return(NULL)
    list(proportions = result$w, chosen_variants = result$ids, residual = result$rmse)
  }







  angle_based_mesma <- function(y, library_list, max_components = NULL) {
    veg_names <- names(library_list)
    n_libs <- length(veg_names)
    if (n_libs == 0) return(NULL)
    if (n_libs == 1) {
      lib <- library_list[[1]]
      best_dist <- Inf
      best_idx <- 1
      for (i in seq_along(lib)) {
        dist <- sqrt(sum((y - lib[[i]]$vec)^2))
        if (dist < best_dist) {
          best_dist <- dist; best_idx <- i
        }
      }
      return(list(fractions = setNames(1.0, veg_names[1]), chosen = setNames(lib[[best_idx]]$id, veg_names[1]), residual = best_dist))
    }
    if (n_libs == 2) {
      if (is.null(res$m1) || is.null(res$m2)) return(NULL)
      fractions <- c(unmix_res$f1, unmix_res$f2)
      names(fractions) <- veg_names
      chosen <- c(res$m1$id, res$m2$id); names(chosen) <- veg_names
      return(list(fractions = fractions, chosen = chosen, residual = unmix_res$residual))
    }
    lib_sizes <- sapply(library_list, length)
    lib_order <- order(lib_sizes)
    v1 <- veg_names[lib_order[1]]; v2 <- veg_names[lib_order[2]]
    if (is.null(pair_result$m1)) return(NULL)
    selected_endmembers <- list(); selected_endmembers[[v1]] <- pair_result$m1; selected_endmembers[[v2]] <- pair_result$m2
    for (k in seq(3, n_libs)) {
      vk <- veg_names[lib_order[k]]
      best_dist <- Inf; best_em <- NULL
      for (candidate in library_list[[vk]]) {
        M_cols <- lapply(selected_endmembers, function(em) em$vec)
        M_cols[[length(M_cols) + 1]] <- candidate$vec
        M <- do.call(cbind, M_cols)
        if (unmix_result$residual < best_dist) {
          best_dist <- unmix_result$residual; best_em <- candidate
        }
      }
      if (!is.null(best_em)) selected_endmembers[[vk]] <- best_em
    }
    M_final <- do.call(cbind, lapply(selected_endmembers, function(em) em$vec))
    colnames(M_final) <- names(selected_endmembers)
    fractions <- final_unmix$f; names(fractions) <- names(selected_endmembers)
    chosen <- sapply(selected_endmembers, function(em) em$id)
    list(fractions = fractions, chosen = chosen, residual = final_unmix$residual)
  }

  stage1_debug_env <- new.env()
  stage1_debug_env$count <- 0
  
  ols_stage1_unmix <- function(y, barren_endmember, vegetation_endmember) {
    # Stage-1 unmix using weighted Euclidean distance (no L2 normalization)
    if (is.null(barren_endmember) || is.null(vegetation_endmember)) return(list(veg_frac = NA_real_, barren_frac = NA_real_, residual = NA_real_, y_proj = NA))
    
    B <- as.numeric(barren_endmember)
    X <- as.numeric(y)

    # Handle single or multiple vegetation endmembers
    V <- if (is.vector(vegetation_endmember)) {
        matrix(as.numeric(vegetation_endmember), ncol = 1)
    } else if (is.matrix(vegetation_endmember)) {
        vegetation_endmember
    } else {
        stop("vegetation_endmember must be a vector or a matrix")
    }

    E <- cbind(B, V)
    n_endmembers <- ncol(E)
    
    len <- nrow(E)
    if(length(X) != len) {
        warning(sprintf("Length of observation y (%d) does not match number of features in endmembers (%d), truncating", length(X), len))
        X <- X[1:len]
    }
    X[!is.finite(X)] <- 0
    E[!is.finite(E)] <- 0

    # Check for degenerate inputs
    if (sqrt(sum(X^2)) < 1e-9 || any(sqrt(colSums(E^2)) < 1e-9)) {
      barren_frac <- 1.0
      veg_frac <- 0.0
      w <- c(barren_frac, rep(veg_frac, n_endmembers - 1))
      return(list(veg_frac = veg_frac, barren_frac = barren_frac, residual = Inf, y_proj = X, w = w))
    }

    res <- nnls::nnls(E, X)
    w <- res$x
    
    # sum-to-one constraint post-NNLS
    if (sum(w) > 0) {
      w <- w / sum(w)
    }

    y_proj <- as.numeric(E %*% w)
    residual <- sqrt(sum((X - y_proj)^2))

    barren_frac <- as.numeric(w[1])
    veg_frac <- sum(w[-1])
    
    return(list(veg_frac = veg_frac, barren_frac = barren_frac, residual = residual, y_proj = y_proj, w = w))
  }

  ols_stage2_unmix <- function(y, vegetation_libraries, topK = 2) {
    res <- stage2_ols_unmix(y, vegetation_libraries, topK = topK)
    if (is.null(res)) return(NULL)
    list(proportions = res$proportions, chosen_variants = res$chosen_variants, residual = res$residual)
  }

  fit_one_task <- function(task_data) {

    X <- as.numeric(y)
    if (is.list(barren_endmember) || is.list(vegetation_endmember)) {
      B_vars <- if (is.list(barren_endmember)) barren_endmember else list(list(vec = as.numeric(barren_endmember), id = 'B1'))
      V_vars <- if (is.list(vegetation_endmember)) vegetation_endmember else list(list(vec = as.numeric(vegetation_endmember), id = 'V1'))
      b_rep <- if (!is.null(B_vars[[1]]$vec)) B_vars[[1]]$vec else as.numeric(B_vars[[1]])
      v_rep <- if (!is.null(V_vars[[1]]$vec)) V_vars[[1]]$vec else as.numeric(V_vars[[1]])
      s1_res <- ols_stage1_unmix(X, b_rep, v_rep)
      veg_fraction_total <- s1_res$veg_frac; barren_fraction_total <- s1_res$barren_frac
    } else {
      stage1_unmix <- ols_stage1_unmix(X, barren_endmember, vegetation_endmember)
      veg_fraction_total <- stage1_unmix$veg_frac; barren_fraction_total <- stage1_unmix$barren_frac
    }
    if (!is.finite(veg_fraction_total)) {
      if (isTRUE(compute_stage2_diagnostics)) {
        return(list(location_id = loc, pheno_year = yr, n_raw_bins = length(y_raw), n_valid_mask = sum(valid_mask), y_s2_norm_val = NA_real_, y_s2_masked_mean = NA_real_, weights_s2_fallback_used = isTRUE(weights_s2_fallback_used)))
      }
      return(list(fractions = NULL, stage1 = stage1_result, stage2 = NULL, error = 'Stage 1 failed'))
    }
    if (isTRUE(veg_fraction_total <= 0.01)) {
      if (isTRUE(compute_stage2_diagnostics)) {
        return(list(location_id = loc, pheno_year = yr, n_raw_bins = length(y_raw), n_valid_mask = sum(valid_mask), y_s2_norm_val = NA_real_, y_s2_masked_mean = NA_real_, weights_s2_fallback_used = isTRUE(weights_s2_fallback_used)))
      }
      return(list(fractions = c(barren = 1.0), stage1 = stage1_result, stage2 = NULL, veg_fraction = 0, barren_fraction = 1))
    }
    stage2_result <- stage2_ols_unmix(X, vegetation_libraries, topK = topK)
    if (is.null(stage2_result)) {
      if (isTRUE(compute_stage2_diagnostics)) {
        return(list(location_id = loc, pheno_year = yr, n_raw_bins = length(y_raw), n_valid_mask = sum(valid_mask), y_s2_norm_val = NA_real_, y_s2_masked_mean = NA_real_, weights_s2_fallback_used = isTRUE(weights_s2_fallback_used)))
      }
      return(list(fractions = c(barren = barren_fraction_total, vegetation = veg_fraction_total), stage1 = stage1_result, stage2 = NULL, veg_fraction = veg_fraction_total, barren_fraction = barren_fraction_total))
    }
    stage2_proportions <- stage2_result$proportions
    final_fractions <- stage2_proportions * veg_fraction_total
    final_fractions <- c(final_fractions, barren = barren_fraction_total)
    total <- sum(final_fractions); if (abs(total - 1) > 1e-6 && total > 0) final_fractions <- final_fractions / total
    list(fractions = final_fractions, stage1 = stage1_result, stage2 = stage2_result, veg_fraction = veg_fraction_total, barren_fraction = barren_fraction_total, chosen_variants = stage2_result$chosen_variants, residual = list(stage1 = if (exists('unmix_stage1')) unmix_stage1$residual else NA_real_, stage2 = stage2_result$residual))
  }

  get_variant_templates <- function(veg_types, grid_type, compressed_templates, mesma_lib) {
    templates <- list()
    for (veg in veg_types) {
      veg_templates <- list()
      for (variant in mesma_lib[[veg]]) {
        vec <- compressed_templates[[veg]][[variant$variant_id]][[grid_type]]
        if (!is.null(vec) && length(vec) > 0) {
          veg_templates[[length(veg_templates) + 1]] <- list(vec = vec, id = variant$variant_id)
        }
      }
      templates[[veg]] <- veg_templates
    }
    templates
  }

  evaluate_all_combinations <- function(y, templates, lambda, early_stop_rmse = EARLY_STOP_RMSE_THRESHOLD, feature_weights = NULL) {
    veg_names <- names(templates)
    n_veg <- length(veg_names)
    if (n_veg == 0) return(NULL)
    idx_lists <- lapply(templates, function(lst) seq_along(lst))
    if (any(lengths(idx_lists) == 0)) return(NULL)

    y_fit <- y
    # No extra preprocessing of y — use raw observations (non-finite -> 0)
    y[!is.finite(y)] <- 0
    y_fit <- y

    radices <- as.integer(lengths(idx_lists))
    n_combos <- as.numeric(Reduce(`*`, as.list(radices), init = 1L))
    chunk_size <- 100L

    safe_expand_limit <- COMBO_SAFE_EXPAND_LIMIT
    abort_limit <- COMBO_ABORT_LIMIT

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
        chunk <- combos[start:end, , drop = FALSE]
      } else {
        nrows <- end - start + 1L
        chunk <- matrix(NA_integer_, nrow = nrows, ncol = n_veg)
        colnames(chunk) <- names(idx_lists)
        base_start <- as.numeric(start - 1L)
        for (r in seq_len(nrows)) {
          x <- base_start + (r - 1L)
          rem <- x
          for (j in seq_len(n_veg)) {
            base <- radices[j]
            digit <- (rem %% base) + 1L
            chunk[r, j] <- as.integer(digit)
            rem <- rem %/% base
          }
        }
        chunk <- as.data.frame(chunk, stringsAsFactors = FALSE, check.names = FALSE)
      }
      results <- vector("list", nrow(chunk))
      for (i in seq_len(nrow(chunk))) {
        cols <- list(); ids <- character(0)
        for (v in veg_names) {
          idx <- as.integer(chunk[i, v])
          cols[[length(cols) + 1]] <- templates[[v]][[idx]]$vec
          ids <- c(ids, templates[[v]][[idx]]$id)
        }
        E <- do.call(cbind, cols)
        
        E_fit <- E

        # Use simple NNLS solver consistent with Stage 1 (no extra preprocessing).
        # Pass `feature_weights` so the solver applies sqrt(weights) to rows of E and to y
        # (WLS). Only pass them when they match the number of features (rows of E).
        w_res <- solve_weights_nnls_simple(E_fit, y, feature_weights = if (!is.null(feature_weights) && length(feature_weights) == nrow(E)) feature_weights else NULL)
        w <- if (is.list(w_res) && !is.null(w_res$w)) as.numeric(w_res$w) else as.numeric(w_res)
        rmse <- if (is.list(w_res) && !is.null(w_res$rmse)) w_res$rmse else NA_real_
        results[[i]] <- list(w = w, rmse = rmse, ids = ids)
      }
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

compress_trace <- function(dly_year, avail_idx, budget = N_TEMPORAL_BINS) {
  raw_mat <- build_pentad_matrix(dly_year, avail_idx)
  if (is.null(raw_mat)) return(NULL)
  
  y <- as.numeric(raw_mat)  # Flatten to vector
  
  list(y = y, grid_type = "full", raw_mat = raw_mat)
}



compress_and_unmix_year <- function(dly_year, mesma_lib, budget = TEMPORAL_BUDGET, topK = TOPK_VARIANTS) {
  res <- compress_trace(dly_year, avail, budget)
  if (is.null(res)) return(NULL)
  unmix_stage2_compressed(res$y, res$grid_type, mesma_lib, topK, feature_weights = LDA_FEATURE_WEIGHTS)
}
if (exists('DEBUG_STAGE1_VERBOSE', envir = globalenv()) && isTRUE(get('DEBUG_STAGE1_VERBOSE', envir = globalenv())) && exists('COMPRESSED_STAGE1_LIB') && !is.null(COMPRESSED_STAGE1_LIB)) {
  cat('\n=== STAGE1 RUNTIME DEBUG BLOCK ===\n')
  sample_rows <- tryCatch({ df_train |> dplyr::filter(.data$Veg %in% ALLOWED_VEG) |> dplyr::distinct(location_id, pheno_year) |> head(10) }, error = function(e) NULL)
  if (!is.null(sample_rows) && nrow(sample_rows) > 0) {
    for (i in seq_len(nrow(sample_rows))) {
      loc <- sample_rows$location_id[i]; yr <- sample_rows$pheno_year[i]
      dly_year <- df_train |> dplyr::filter(.data$location_id == loc & .data$pheno_year == yr)
      raw_mat <- build_pentad_matrix(dly_year, COMPRESSED_STAGE1_LIB$indices)
      if (!is.null(raw_mat)) {
        y_vec_raw <- as.numeric(raw_mat)
        # Apply same preprocessing as in fit_one_task
        y_vec <- y_vec_raw
        n_bins <- N_TEMPORAL_BINS
        for (k in seq_along(COMPRESSED_STAGE1_LIB$indices)) {
          idx_start <- (k-1)*n_bins + 1
          idx_end <- k*n_bins
          mu <- if (exists("STAGE1_PARAMS")) STAGE1_PARAMS$means[k] else 0
          sigma <- if (exists("STAGE1_PARAMS") && is.finite(STAGE1_PARAMS$sds[k]) && STAGE1_PARAMS$sds[k] > 1e-10) STAGE1_PARAMS$sds[k] else 1
          y_vec[idx_start:idx_end] <- (y_vec[idx_start:idx_end] - mu) / sigma
        }
        y_vec[!is.finite(y_vec)] <- 0
        if (exists("STAGE1_PARAMS") && !is.null(STAGE1_PARAMS$weights)) {
          # Apply sqrt(weights) to match how templates are weighted (WLS consistency)
          y_vec <- y_vec * sqrt(pmax(STAGE1_PARAMS$weights, 0))
        }
        stage1_res <- ols_stage1_unmix(y_vec, COMPRESSED_STAGE1_LIB$barren, COMPRESSED_STAGE1_LIB$vegetation)
        cat(sprintf("sample %d: loc=%s pheno_year=%d y_norm_raw=%.4f t=%.6f f2=%.6f resid=%.6f barren_norm_raw=%.4f veg_norm_raw=%.4f\n", i, as.character(loc), as.integer(yr), sqrt(sum(y_vec_raw^2)), ifelse(is.null(stage1_res$projection), NA_real_, stage1_res$projection[1]), stage1_res$veg_frac, stage1_res$residual, COMPRESSED_STAGE1_LIB$norm_b, COMPRESSED_STAGE1_LIB$norm_v))
      }
    }
  }
  cat('=== END STAGE1 RUNTIME DEBUG BLOCK ===\n\n')
}

prepare_factor_data <- function(dly, gpca, avail_idx, veg_type) {
    date_list <- list()
    dly$doy <- pheno_doy(dly$date)  # Use phenological DOY
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
      z <- solve(xtx_reg, xty)

      if (any(!is.finite(z))) {
        stop("All matrix inversion methods failed")
      }

      date_list[[length(date_list) + 1]] <- list(doy = as.integer(sub$doy[1]), z = as.numeric(z))
    }
    date_list
  }

  build_Z365 <- function(date_list, k) {
    Z <- matrix(NA_real_, nrow = N_TEMPORAL_BINS, ncol = k)
    
    pentad_data <- list()
    for (entry in date_list) {
      if (!is.null(entry$doy) && !is.null(entry$z) && 
          entry$doy >= 1 && entry$doy <= 366 && 
          length(entry$z) == k) {
        pentad <- doy_to_pentad(entry$doy)
        if (is.null(pentad_data[[as.character(pentad)]])) {
          pentad_data[[as.character(pentad)]] <- list()
        }
        pentad_data[[as.character(pentad)]][[length(pentad_data[[as.character(pentad)]]) + 1]] <- entry$z
      }
    }
    
    for (pentad_str in names(pentad_data)) {
      pentad <- as.integer(pentad_str)
      if (pentad >= 1 && pentad <= N_TEMPORAL_BINS) {
        z_vals <- do.call(rbind, pentad_data[[pentad_str]])
        if (!is.null(z_vals) && nrow(z_vals) > 0) {
          Z[pentad, ] <- apply(z_vals, 2, median, na.rm = TRUE)
        }
      }
    }
    
    Z[!is.finite(Z)] <- NA_real_  # Keep NA instead of 0
    Z
  }

  precompute_whitening <- function(mesma_lib, compressed_templates, combinations) {
    function(combo_key) {
      list(whiten = function(x) x, unwhiten = function(x) x)
    }
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
      COMPRESSED_STAGE1_LIB = if (exists("COMPRESSED_STAGE1_LIB")) COMPRESSED_STAGE1_LIB else NULL,
      STAGE1_PARAMS = if (exists("STAGE1_PARAMS")) STAGE1_PARAMS else NULL,
      STAGE2_PARAMS = if (exists("STAGE2_PARAMS")) STAGE2_PARAMS else NULL,
      RAW_BARREN_PROTO = if (exists("RAW_BARREN_PROTO")) RAW_BARREN_PROTO else NULL
    )
    saveRDS(core, file = file.path(cache_dir, "mesma_library.rds"))

    raw_templates <- list(
      raw_lib_templates = if (exists("raw_lib_templates")) raw_lib_templates else NULL
    )
    saveRDS(raw_templates, file = file.path(cache_dir, "raw_templates.rds"))

    if (exists(".COMPRESSED_TEMPLATES_ACCESSOR", envir = globalenv())) {
      template_data <- list()
      for (veg in names(mesma_lib)) {
        template_data[[veg]] <- list()
        for (variant in mesma_lib[[veg]]) {
          for (grid_type in c("sparse", "full", "dense")) {
            key <- paste(veg, variant$variant_id, grid_type, sep = "|")
            vec <- tryCatch(.COMPRESSED_TEMPLATES_ACCESSOR[[veg]][[variant$variant_id]][[grid_type]], error = function(e) stop(sprintf("Compressed templates accessor failed for %s|%s|%s: %s", veg, variant$variant_id, grid_type, e$message)))
            if (!is.null(vec)) template_data[[veg]][[key]] <- vec
          }
        }
      }
      saveRDS(template_data, file = file.path(cache_dir, "compressed_templates.rds"))
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

if (isTRUE(TESTING_MODE)) {
  cat("Skipping inference data loading (TESTING_MODE = TRUE).\n")
  df_inf <- NULL
  inference_location_ids <- character(0)
} else if (!isTRUE(SKIP_INFERENCE)) {
  df_inf <- NULL
  if (exists("INFERENCE_CSV")) {
    cat(sprintf("Checking inference file at: %s\n", INFERENCE_CSV))
    if (file.exists(INFERENCE_CSV)) {
      cat(sprintf("Loading inference data from %s...\n", INFERENCE_CSV))
      
      if (grepl("\\.csv$", INFERENCE_CSV, ignore.case = TRUE)) {
         df_inf <- tryCatch(read.csv(INFERENCE_CSV), error = function(e) {
            cat(sprintf("[WARNING] Error reading inference CSV: %s\n", e$message))
            NULL
         })
      } else {
         df_inf <- tryCatch(openxlsx::read.xlsx(INFERENCE_CSV), error = function(e) {
            cat(sprintf("[WARNING] Error reading inference XLSX: %s\n", e$message))
            NULL
         })
      }
      
      if (!is.null(df_inf)) {
         cat(sprintf("Loaded %d rows from inference file.\n", nrow(df_inf)))
         
         # Normalize location_id column name
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
          df_inf$location_id <- paste0(df_inf$location_id_orig, "a")
          cat(sprintf("Appended 'a' to inference location_IDs for uniqueness. New unique count: %d\n", length(unique(df_inf$location_id))))
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
       if (is.character(df_inf$...1) || is.numeric(df_inf$...1)) {
         names(df_inf)[names(df_inf) == "...1"] <- "location_id"
       }
    }
    
    if ("prediction_date" %in% names(df_inf)) {
      df_inf$date <- as.Date(df_inf$prediction_date)
    } else if ("date" %in% names(df_inf)) {
      df_inf$date <- as.Date(df_inf$date)
    } else {
      for (col in names(df_inf)) {
         if (inherits(df_inf[[col]], "Date")) {
           df_inf$date <- df_inf[[col]]
           break
         }
         if (is.character(df_inf[[col]]) && all(grepl("^\\d{4}-\\d{2}-\\d{2}", na.omit(df_inf[[col]][1:min(10, nrow(df_inf))])))) {
           df_inf$date <- as.Date(df_inf[[col]])
           break
         }
      }
    }
    
    if ("location_id" %in% names(df_inf) && "date" %in% names(df_inf)) {
      df_inf$location_id <- as.character(df_inf$location_id)
      
      original_rows <- nrow(df_inf)
      df_inf <- df_inf |> distinct(location_id, date, .keep_all = TRUE)
      if (nrow(df_inf) < original_rows) {
        cat(sprintf("Inference data deduplicated: %d rows remaining from %d original rows.\n", nrow(df_inf), original_rows))
      }
      
      if (length(intersect(RAW_BANDS, names(df_inf))) >= 2) {
        before_cols <- names(df_inf)
        df_inf <- compute_indices_from_bands(df_inf)
        new_cols <- setdiff(names(df_inf), before_cols)
        if (length(new_cols) > 0) cat(sprintf("[NOTICE] Computed indices from raw bands in inference data: %s\n", paste(new_cols, collapse=", ")))
      }

      # Filter out observations with critical snow or dust contamination
      if ("NDSI" %in% names(df_inf) && "NDDI" %in% names(df_inf)) {
        snow_count <- sum(df_inf$NDSI > 0.4, na.rm = TRUE)
        dust_count <- sum(df_inf$NDDI > 0.18, na.rm = TRUE)
        total_before <- nrow(df_inf)
        df_inf <- df_inf[!(df_inf$NDSI > 0.4 | df_inf$NDDI > 0.18), , drop = FALSE]
        total_after <- nrow(df_inf)
        filtered <- total_before - total_after
        cat(sprintf("Filtered out %d observations with snow (NDSI > 0.4) or dust (NDDI > 0.18) contamination from inference data\n", filtered))
        cat(sprintf("Inference dataset after contamination filtering: %d rows from %d locations\n", total_after, length(unique(df_inf$location_id))))
      } else {
        cat("[WARNING] NDSI or NDDI not found in inference data; skipping contamination filtering\n")
      }

      # Add required columns before attempting PPI calculation
      if (!"Veg" %in% names(df_inf)) df_inf$Veg <- NA_character_
      if (!"pheno_year" %in% names(df_inf)) df_inf$pheno_year <- assign_pheno_year(df_inf$date)
      if (!"doy" %in% names(df_inf)) df_inf$doy <- pheno_doy(df_inf$date)  # Use phenological DOY
      if (!"zenith.angle" %in% names(df_inf)) df_inf$zenith.angle <- NA_real_
      if (!"DVI_max" %in% names(df_inf)) df_inf$DVI_max <- NA_real_

      # Attempt to auto-add PPI to inference data using training baseline BEFORE checking for missing indices
      if (exists("auto_add_ppi_columns") && "PPI" %in% avail && !"PPI" %in% names(df_inf)) {
        ## --- Check for the saved training baseline ---
        dvi_soil_arg <- NA_real_
        if (exists("GLOBAL_TRAINING_DVI_SOIL", envir = globalenv())) {
           dvi_soil_arg <- get("GLOBAL_TRAINING_DVI_SOIL", envir = globalenv())
           cat(sprintf("[PPI] Inference: Using dynamic baseline from training data (%.6f)\n", dvi_soil_arg))
        }
        ## ------------------------------------------------
        # Pass the retrieved baseline to the function
        ppi_inf_res <- tryCatch({
            if (!is.na(dvi_soil_arg) && is.finite(dvi_soil_arg)) {
                auto_add_ppi_columns(df_inf, dvi_soil = dvi_soil_arg)
            } else {
                auto_add_ppi_columns(df_inf) # Fallback to default behavior
            }
        }, error = function(e) {
            cat(sprintf("[PPI] auto_add_ppi_columns failed: %s\n", e$message))
            NULL
        })

        if (!is.null(ppi_inf_res) && isTRUE(ppi_inf_res$added)) {
          df_inf <- ppi_inf_res$df
          cat(sprintf("[PPI] Auto-added PPI to inference data (reason: %s)\n", ppi_inf_res$reason))
        } else if (!is.null(ppi_inf_res) && identical(ppi_inf_res$reason, "no_baseline")) {
          cat("[PPI] Inference PPI not computed: no barren baseline available in inference data and no training baseline found.\n")
        }
      }

      # Now check for missing indices AFTER attempting to compute PPI
      missing_idx <- setdiff(avail, names(df_inf))
      if (length(missing_idx) > 0) {
        cat(sprintf("[WARNING] Inference data missing indices: %s. Filling with NA (will likely fail unmixing).\n", paste(missing_idx, collapse=", ")))
        for (col in missing_idx) df_inf[[col]] <- NA_real_
      }
      
      if ("prediction_date" %in% names(df_inf)) {
        df_inf$prediction_date <- as.Date(df_inf$prediction_date)
      }
      if ("reference_date" %in% names(df_inf)) {
        df_inf$reference_date <- as.Date(df_inf$reference_date)
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
      
      df_inf <- df_inf |> filter(pheno_year >= 1985 & pheno_year <= 2025)

      for (col in avail) {
        if (col %in% names(df_inf)) df_inf[[col]] <- as.numeric(df_inf[[col]])
      }

      df_tasks_inference <- df_inf
      cat(sprintf("Created separate inference task list with %d rows from %d locations.\n", nrow(df_tasks_inference), length(unique(df_tasks_inference$location_id))))
      
      before_dedup <- nrow(df_tasks_inference)
      df_tasks_inference <- df_tasks_inference |> dplyr::distinct(location_id, date, .keep_all = TRUE)
      after_dedup <- nrow(df_tasks_inference)
      if (before_dedup > after_dedup) {
        cat(sprintf("Removed %d duplicate rows from inference data.\n", before_dedup - after_dedup))
      }

      if (!"pheno_year" %in% names(df_tasks_inference) && "date" %in% names(df_tasks_inference)) df_tasks_inference$pheno_year <- assign_pheno_year(df_tasks_inference$date)
      n_infer_loc_years <- nrow(unique(df_tasks_inference[c("location_id", "pheno_year")]))
      if (exists("df_train") && !is.null(df_train) && nrow(df_train) > 0) {
        if (!"pheno_year" %in% names(df_train) && "date" %in% names(df_train)) df_train$pheno_year <- assign_pheno_year(df_train$date)
        n_train_loc_years <- nrow(unique(df_train[c("location_id", "pheno_year")]))
      } else {
        n_train_loc_years <- 0
      }
      cat(sprintf("(NOTICE) Inference dataset location-years: %d\n", n_infer_loc_years))
      cat(sprintf("(NOTICE) Training dataset location-years: %d\n", n_train_loc_years))
      if (n_train_loc_years > 0 && n_infer_loc_years == n_train_loc_years) {
        cat(sprintf("(WARNING) Training and inference datasets have the same number of location-years (%d). This may be expected if IDs are independent; no automatic filtering will be applied.\n", n_train_loc_years))
      }

      cat("Keeping df_tasks as training data for separate processing.\n")
    } else {
      cat("[WARNING] Inference data missing 'location_id' or 'date' column. Skipping.\n")
      cat(sprintf("Columns found: %s\n", paste(names(df_inf), collapse=", ")))
    }
  }
  inference_location_ids <- if (exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0 && "location_id" %in% names(df_tasks_inference)) unique(df_tasks_inference$location_id) else character(0)
} else {
  cat("Skipping inference data loading (SKIP_INFERENCE = TRUE).\n")
  df_inf <- NULL
  inference_location_ids <- character(0)
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
      if ("no.soil.geo" %in% names(df_tasks)) {
        if ("no.soil" %in% names(df_tasks)) {
          df_tasks$`no.soil` <- ifelse(is.na(df_tasks$`no.soil`), df_tasks$`no.soil.geo`, df_tasks$`no.soil`)
        } else {
          df_tasks$`no.soil` <- df_tasks$`no.soil.geo`
        }
        df_tasks$`no.soil.geo` <- NULL
      }
    } else {
      df_tasks <- dplyr::left_join(df_tasks, gpts_map, by = "location_id", suffix = c("", ".geo"))
      if ("Veg.geo" %in% names(df_tasks)) {
        df_tasks$Veg <- ifelse(is.na(df_tasks$Veg) | df_tasks$Veg == "", df_tasks$Veg.geo, df_tasks$Veg)
        df_tasks$Veg.geo <- NULL
      }
      if ("no.soil.geo" %in% names(df_tasks)) {
        if ("no.soil" %in% names(df_tasks)) {
          df_tasks$`no.soil` <- ifelse(is.na(df_tasks$`no.soil`), df_tasks$`no.soil.geo`, df_tasks$`no.soil`)
        } else {
          df_tasks$`no.soil` <- df_tasks$`no.soil.geo`
        }
        df_tasks$`no.soil.geo` <- NULL
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

#     cat(sprintf("\n=== SPECIFIC TEST LOCATION ===\n"))
#     cat(sprintf("Location %s, Year %d:\n", test_task$loc, test_task$yr))
#     cat(sprintf("  Total observations: %d\n", nrow(test_dly_year)))
#     if (nrow(test_dly_year) > 0) {
#       if ("date" %in% names(test_dly_year)) {
#         cat(sprintf("  Date range: %s to %s\n", min(test_dly_year$date, na.rm = TRUE), max(test_dly_year$date, na.rm = TRUE)))
#       }
#       cat(sprintf("  Unique DOYs: %d\n", length(unique(test_dly_year$doy))))
#       if (any(!is.na(test_dly_year$doy))) {
#         cat(sprintf("  DOY range: %d to %d\n", min(test_dly_year$doy, na.rm = TRUE), max(test_dly_year$doy, na.rm = TRUE)))
#       }
#     }
#     cat("======================\n\n")
#   } else {
#     cat("df_tasks not present or empty — skipping data distribution analysis.\n")
#   }
# 
# 

  if (length(location_list) == 0) {
    stop("No testing tasks found")
  }


  build_temporal_matrix <- function(dly_year, avail) {
    if (nrow(dly_year) == 0 || length(avail) == 0) return(NULL)
    
    if (!"doy" %in% names(dly_year)) {
      if ("date" %in% names(dly_year)) {
        dly_year$doy <- pheno_doy(dly_year$date)  # Use phenological DOY
      } else {
        return(NULL)
      }
    }
    
    dly_year$pentad <- doy_to_pentad(dly_year$doy)

    temporal_matrix <- matrix(NA_real_, nrow = N_TEMPORAL_BINS, ncol = length(avail))
    colnames(temporal_matrix) <- avail

    avail_present <- intersect(avail, names(dly_year))
    if (length(avail_present) == 0) return(NULL)

    for (idx in avail_present) {
      vals_by_pentad <- tapply(dly_year[[idx]], dly_year$pentad, function(v) {
        v <- v[is.finite(v)]
        if (length(v) == 0) NA_real_ else median(v)
      })
      pentad_values <- as.integer(names(vals_by_pentad))
      valid <- pentad_values >= 1 & pentad_values <= N_TEMPORAL_BINS
      col_idx <- match(idx, avail)
      if (!is.na(col_idx)) {
        temporal_matrix[pentad_values[valid], col_idx] <- vals_by_pentad[valid]
      }
    }

    temporal_matrix
  }

  build_raw_index_matrix <- function(temporal_matrix, avail) {
    if (is.null(temporal_matrix) || length(avail) == 0) return(NULL)
    
    if (!is.null(colnames(temporal_matrix))) {
      available_indices <- intersect(avail, colnames(temporal_matrix))
      if (length(available_indices) > 0) {
        raw_matrix <- temporal_matrix[, available_indices, drop = FALSE]
      } else {
        if (ncol(temporal_matrix) >= length(avail)) {
          raw_matrix <- temporal_matrix[, seq_along(avail), drop = FALSE]
          colnames(raw_matrix) <- avail
        } else {
          raw_matrix <- temporal_matrix
        }
      }
    } else {
      raw_matrix <- temporal_matrix
      if (ncol(raw_matrix) == length(avail)) {
        colnames(raw_matrix) <- avail
      }
    }
    
    raw_matrix[!is.finite(raw_matrix)] <- NA
    
    raw_matrix
  }

compress_temporal_matrix <- function(data, temporal_budget, avail) {
  list(
    compressed_matrix = data,
    grid = seq_len(N_TEMPORAL_BINS),  # 1:73 for 5-day pentads
    grid_type = "full"
  )
}
  compress_stage1_lib_unified <- function(stage1_lib, avail, temporal_budget) {
    if (is.null(stage1_lib)) return(NULL)

    compressed_lib <- list()

    for (endmember_name in names(stage1_lib)) {
      endmember_data <- stage1_lib[[endmember_name]]

      compressed <- compress_temporal_matrix(endmember_data, temporal_budget, avail)

      if (!is.null(compressed)) {
        compressed_lib[[endmember_name]] <- list(
          compressed = compressed$compressed_matrix,
          grid = compressed$grid,
          grid_type = compressed$grid_type,
          veg_frac = compressed$veg_frac
        )
      }
    }

    compressed_lib
  }

  precompute_compressed_templates_unified <- function(mesma_lib, avail, temporal_budget) {
    if (is.null(mesma_lib)) return(NULL)

    compressed_templates <- list()

    for (veg_type in names(mesma_lib)) {
      veg_variants <- mesma_lib[[veg_type]]

      compressed_variants <- list()

      for (variant_name in names(veg_variants)) {
        variant_data <- veg_variants[[variant_name]]

        compressed <- compress_temporal_matrix(variant_data, temporal_budget, avail)

        if (!is.null(compressed)) {
          compressed_variants[[variant_name]] <- list(
            compressed = compressed$compressed_matrix,
            grid = compressed$grid,
            grid_type = compressed$grid_type,
            veg_frac = compressed$veg_frac
          )
        }
      }

      compressed_templates[[veg_type]] <- compressed_variants
    }

    compressed_templates
  }

  reduce_all_traces_unified <- function(traces, avail, temporal_budget) {
    if (is.null(traces) || length(traces) == 0) return(list())

    compressed_traces <- list()

    for (trace_name in names(traces)) {
      trace_data <- traces[[trace_name]]

      compressed <- compress_trace(trace_data, avail, temporal_budget)

      if (!is.null(compressed)) {
        compressed_traces[[trace_name]] <- compressed
      }
    }

    compressed_traces
  }

  build_mesma_variants_unified <- function(lib, avail, temporal_budget) {
    if (is.null(lib)) return(NULL)

    mesma_variants <- list()

    for (veg_type in names(lib)) {
        veg_data <- lib[[veg_type]]
        n_samp <- NA_integer_
        if (!is.null(veg_data$n_samples)) n_samp <- veg_data$n_samples
        if (is.na(n_samp) && !is.null(veg_data$features) && !is.null(dim(veg_data$features))) {
          n_samp <- nrow(veg_data$features)
        }
        if (!is.na(n_samp) && n_samp < MIN_ENDMEMBER_SAMPLES) {
          cat(sprintf("[NOTICE] Skipping unified variant for '%s' due to insufficient samples: %d < %d\n", veg_type, n_samp, MIN_ENDMEMBER_SAMPLES))
          next
        }

      compressed <- compress_temporal_matrix(veg_data, temporal_budget, avail)

      if (!is.null(compressed)) {
        mesma_variants[[veg_type]] <- list(
          compressed = compressed$compressed_matrix,
          grid = compressed$grid,
          grid_type = compressed$grid_type,
          veg_frac = compressed$veg_frac
        )
      }
    }

    mesma_variants
  }

  precompute_templates_from_unified_lib <- function(unified_lib, veg_types) {
    if (is.null(unified_lib)) return(NULL)

    function(veg_type, variant_id = NULL) {
      if (!veg_type %in% names(unified_lib)) return(NULL)

      veg_templates <- unified_lib[[veg_type]]

      if (is.null(variant_id)) {
        return(veg_templates)
      } else {
        return(veg_templates[[variant_id]])
      }
    }
  }


  evaluate_all_combinations <- function(y, top_variants, lambda = 0, early_stop_rmse = EARLY_STOP_RMSE_THRESHOLD, feature_weights = NULL) {
    if(length(top_variants) == 0) return(NULL)
    
    veg_names <- names(top_variants)
    n_veg <- length(veg_names)
    
    # Use raw observation vector (no L2-normalization or preprocessing)
    y_target <- y
    y_target[!is.finite(y_target)] <- 0

    idx_lists <- lapply(top_variants, function(lst) seq_along(lst))
    radices <- as.integer(lengths(idx_lists))
    n_veg <- length(radices)
    n_combos <- as.numeric(Reduce(`*`, as.list(radices), init = 1L))

    safe_expand_limit <- COMBO_SAFE_EXPAND_LIMIT
    abort_limit <- COMBO_ABORT_LIMIT
    if (!is.na(n_combos) && n_combos > abort_limit) {
      best_variants <- lapply(top_variants, function(x) x[[1]])
      if (any(sapply(best_variants, function(b) is.null(b) || is.null(b$vec)))) {
        w_unif <- rep(1/n_veg, n_veg); names(w_unif) <- veg_names; return(list(w = w_unif, rmse = Inf, ids = sapply(best_variants, function(b) if (is.null(b)) NA else b$id)))
      }
      E <- do.call(cbind, lapply(best_variants, function(b) as.numeric(b$vec)))
      # Use stage-1 style NNLS solver on raw E/y (no L2-normalization)
      res <- solve_weights_nnls_simple(E, y_target)
      names(res$w) <- veg_names
      return(list(w = res$w, rmse = res$rmse, ids = sapply(best_variants, function(b) b$id)))
    }

    chunk_size <- 100L
    best_rmse <- Inf; best_w <- NULL; best_ids <- NULL
    n_veg <- length(idx_lists)
    if (n_veg == 0) return(NULL)
    for (start in seq(1, n_combos, by = chunk_size)) {
      end <- min(start + chunk_size - 1L, n_combos)
      nrows <- end - start + 1L
      chunk <- matrix(NA_integer_, nrow = nrows, ncol = n_veg)
      colnames(chunk) <- names(idx_lists)
      base_start <- as.numeric(start - 1L)
      for (r in seq_len(nrows)) {
        x <- base_start + (r - 1L)
        rem <- x
        for (j in seq_len(n_veg)) {
          base <- radices[j]
          digit <- (rem %% base) + 1L
          chunk[r, j] <- as.integer(digit)
          rem <- rem %/% base
        }
      }

        for (i in seq_len(nrow(chunk))) {
          cols <- list(); ids <- character(0); skip_combo <- FALSE
          for (v in names(idx_lists)) {
            idx <- as.integer(chunk[i, v])
            cand <- top_variants[[v]][[idx]]
            if (is.null(cand) || is.null(cand$vec)) {
              skip_combo <- TRUE
              break
            }
            cols[[length(cols) + 1]] <- as.numeric(cand$vec)
            ids <- c(ids, cand$id)
          }
          if (isTRUE(skip_combo) || length(cols) == 0) next
          E <- tryCatch({ do.call(cbind, cols) }, error = function(e) {
            return(NULL)
          })
          if (is.null(E) || ncol(E) < 1) next
          
          # Solve using stage-1 style NNLS solver on raw E/y (no L2-normalization)
          # Pass feature_weights so the solver applies sqrt(weights) to rows of E and to y (WLS)
          res <- solve_weights_nnls_simple(E, y_target, feature_weights = if (!is.null(feature_weights) && length(feature_weights) == nrow(E)) feature_weights else NULL)
          
          if (is.list(res) && !is.null(res$rmse) && res$rmse < best_rmse) {
            best_rmse <- res$rmse
            best_w <- res$w
            best_ids <- ids
          }
        }
      if (is.finite(best_rmse) && !is.null(early_stop_rmse) && is.finite(early_stop_rmse) && best_rmse <= early_stop_rmse) break
    }

    if (is.null(best_w)) return(NULL)
    names(best_w) <- veg_names
    list(w = best_w, rmse = best_rmse, ids = best_ids)
  }


  build_weighted_stage1_lib <- function(df_train, indices, params) {
    class_col <- if ("stage1_class" %in% names(df_train)) "stage1_class" else "Veg"
    barren_data <- df_train |> dplyr::filter(.data[[class_col]] == "barren")

    barren_list <- list()
    traces <- unique(barren_data[, c("location_id", "pheno_year")])
    for(i in seq_len(nrow(traces))) {
      sub <- barren_data[barren_data$location_id == traces$location_id[i] &
                         barren_data$pheno_year == traces$pheno_year[i], ]

      mat <- build_pentad_matrix(sub, indices)
      if(!is.null(mat)) {
        vec <- as.numeric(mat)  # Flattens to length(indices) * N_TEMPORAL_BINS
        barren_list[[length(barren_list) + 1]] <- vec
      }
    }

    if(length(barren_list) == 0) {
      stop("No valid barren data found for Stage 1 library building")
    }
    barren_mat <- do.call(rbind, barren_list)

    expected_cols <- length(indices) * N_TEMPORAL_BINS
    if(ncol(barren_mat) != expected_cols) {
      stop(sprintf("barren_mat has %d columns, expected %d (%d indices × %d bins)",
                   ncol(barren_mat), expected_cols, length(indices), N_TEMPORAL_BINS))
    }

    for(k in seq_along(indices)) {
      idx_start <- (k-1)*N_TEMPORAL_BINS + 1
      idx_end <- k*N_TEMPORAL_BINS
      barren_mat[, idx_start:idx_end] <- (barren_mat[, idx_start:idx_end] - params$means[k]) / params$sds[k]
    }
    barren_mat[is.na(barren_mat)] <- 0

      veg_data <- df_train |> dplyr::filter(.data[[class_col]] == "vegetation")
    veg_list <- list()
    traces <- unique(veg_data[, c("location_id", "pheno_year")])
    for(i in seq_len(nrow(traces))) {
      sub <- veg_data[veg_data$location_id == traces$location_id[i] &
                      veg_data$pheno_year == traces$pheno_year[i], ]

      mat <- build_pentad_matrix(sub, indices)
      if(!is.null(mat)) {
        vec <- as.numeric(mat)
        veg_list[[length(veg_list) + 1]] <- vec
      }
    }

    if(length(veg_list) == 0) {
      stop("No valid vegetation data found for Stage 1 library building")
    }
    veg_mat <- do.call(rbind, veg_list)

    if(ncol(veg_mat) != expected_cols) {
      stop(sprintf("veg_mat has %d columns, expected %d (%d indices × %d bins)",
                   ncol(veg_mat), expected_cols, length(indices), N_TEMPORAL_BINS))
    }

    for(k in seq_along(indices)) {
      idx_start <- (k-1)*N_TEMPORAL_BINS + 1
      idx_end <- k*N_TEMPORAL_BINS
      veg_mat[, idx_start:idx_end] <- (veg_mat[, idx_start:idx_end] - params$means[k]) / params$sds[k]
    }
    veg_mat[is.na(veg_mat)] <- 0

    barren_z <- apply(barren_mat, 2, mean)
    
    n_veg_endmembers <- N_STAGE1_VEG_ENDMEMBERS
    if (nrow(veg_mat) > n_veg_endmembers) {
        set.seed(123)
        km <- kmeans(veg_mat, centers = n_veg_endmembers, nstart = 10)
        veg_protos_unweighted <- km$centers
    } else {
        veg_protos_unweighted <- veg_mat
    }

    barren_raw <- colMeans(do.call(rbind, barren_list))
    veg_raw <- colMeans(do.call(rbind, veg_list))

    cos_sim_before <- sum(barren_z * apply(veg_mat, 2, mean)) / (sqrt(sum(barren_z^2)) * sqrt(sum(apply(veg_mat, 2, mean)^2)))
    eucl_dist_before <- sqrt(sum((barren_z - apply(veg_mat, 2, mean))^2))
    cat(sprintf("[DIAGNOSTIC] BEFORE weighting (using mean veg) - Cosine similarity: %.4f, Euclidean distance: %.4f\n", cos_sim_before, eucl_dist_before))

    mean_before_attr <- attr(params$weights, "mean_before_threshold")
    if (!is.null(mean_before_attr)) {
      cat(sprintf("[DIAGNOSTIC] Weight statistics: min=%.6f, max=%.6f, mean_after=%.6f, sd=%.6f, mean_before_threshold=%.6f\n",
                  min(params$weights), max(params$weights), mean(params$weights), sd(params$weights), mean_before_attr))
    } else {
      cat(sprintf("[DIAGNOSTIC] Weight statistics: min=%.6f, max=%.6f, mean=%.6f, sd=%.6f\n",
                  min(params$weights), max(params$weights), mean(params$weights), sd(params$weights)))
    }
    cat(sprintf("[DIAGNOSTIC] Number of weights: %d, Expected: %d\n", length(params$weights), length(barren_z)))


    clip_res <- enforce_stage1_weight_clip(params$weights)
    if (!is.null(clip_res) && is.list(clip_res) && clip_res$clipped > 0) {
      params$weights <- clip_res$weights
      cat(sprintf("[DIAGNOSTIC] Clipped %d weights at threshold=%.3f (old_max=%.3f -> new_max=%.3f), re-normalized to mean=%.6f\n",
                  clip_res$clipped, clip_res$threshold, clip_res$old_max, clip_res$new_max, mean(params$weights, na.rm=TRUE)))
    }

    # Apply sqrt of weights for weighted least squares
    sqrt_w <- sqrt(pmax(params$weights, 0))
    barren_proto <- barren_z * sqrt_w
    veg_protos <- sweep(veg_protos_unweighted, 2, sqrt_w, "*")


    cos_sim_after <- sum(barren_proto * apply(veg_protos, 2, mean)) / (sqrt(sum(barren_proto^2)) * sqrt(sum(apply(veg_protos, 2, mean)^2)))
    eucl_dist_after <- sqrt(sum((barren_proto - apply(veg_protos, 2, mean))^2))
    cat(sprintf("[DIAGNOSTIC] AFTER weighting (using mean veg) - Cosine similarity: %.4f, Euclidean distance: %.4f\n", cos_sim_after, eucl_dist_after))

    try({
      if (!is.na(cos_sim_before) && !is.na(cos_sim_after) && cos_sim_after < 0.7 * cos_sim_before) {
        cat(sprintf("[WARNING] Stage 1 weighting reduced cosine similarity from %.4f -> %.4f (%.1f%%). Consider inspecting weights or disabling Stage1 weighting for diagnostics.\n",
                    cos_sim_before, cos_sim_after, 100 * (cos_sim_after / cos_sim_before)))
      }
      if (!is.na(eucl_dist_before) && !is.na(eucl_dist_after) && eucl_dist_after < 0.5 * eucl_dist_before) {
        cat(sprintf("[WARNING] Stage 1 weighting reduced Euclidean distance from %.4f -> %.4f (%.1f%%). Consider inspecting weights or disabling Stage1 weighting for diagnostics.\n",
                    eucl_dist_before, eucl_dist_after, 100 * (eucl_dist_after / eucl_dist_before)))
      }
    }, silent = TRUE)

    cat(sprintf("[DIAGNOSTIC] Barren mean before weighting: %.6f, after weighting: %.6f\n", mean(barren_z), mean(barren_proto)))
    cat(sprintf("[DIAGNOSTIC] Veg mean before weighting: %.6f, after weighting: %.6f\n", mean(apply(veg_mat, 2, mean)), mean(veg_protos)))

    cat(sprintf("[build_weighted_stage1_lib] Barren samples: %d traces, Vegetation samples: %d traces\n", nrow(barren_mat), nrow(veg_mat)))

    list(
      barren = barren_proto,
      vegetation = t(veg_protos),
      barren_raw = barren_raw,
      veg_raw = veg_raw,
      indices = indices
    )
  }

  build_mesma_library_weighted <- function(df_train, indices, params, allowed_veg) {
    lib <- list()
    expected_cols <- length(indices) * N_TEMPORAL_BINS

    for(v in allowed_veg) {
      veg_data <- dplyr::filter(df_train, .data$Veg == v)
      if(nrow(veg_data) == 0) next

      veg_list <- list()
      traces <- unique(veg_data[, c("location_id", "pheno_year")])
      for(i in seq_len(nrow(traces))) {
        sub <- veg_data[veg_data$location_id == traces$location_id[i] &
                        veg_data$pheno_year == traces$pheno_year[i], ]

        mat <- build_pentad_matrix(sub, indices)
        if(!is.null(mat)) {
          vec <- as.numeric(mat)
          veg_list[[length(veg_list) + 1]] <- vec
        }
      }

      if(length(veg_list) == 0) next

      veg_mat <- do.call(rbind, veg_list)

      if(ncol(veg_mat) != expected_cols) {
        warning(sprintf("Skipping %s: matrix has %d columns, expected %d", v, ncol(veg_mat), expected_cols))
        next
      }

      for(k in seq_along(indices)) {
        idx_start <- (k-1)*N_TEMPORAL_BINS + 1
        idx_end <- k*N_TEMPORAL_BINS
        veg_mat[, idx_start:idx_end] <- (veg_mat[, idx_start:idx_end] - params$means[k]) / params$sds[k]
      }
      veg_mat[is.na(veg_mat)] <- 0

      max_variants <- if (exists("MAX_VARIANTS_PER_VEG", inherits = TRUE)) MAX_VARIANTS_PER_VEG else 7
      n_variants <- min(max_variants, nrow(veg_mat))
      variants_out <- list()
      if (n_variants <= 1) {
        # Apply sqrt of weights for weighted least squares
        vec <- apply(veg_mat, 2, mean) * sqrt(pmax(params$weights, 0))
        variants_out[[1]] <- list(vec = vec, id = paste0(v, "_v1"), raw_mat = veg_mat)
      } else {
        set.seed(1)
        km <- tryCatch(kmeans(veg_mat, centers = n_variants, nstart = 5), error = function(e) NULL)
        if (is.null(km)) {
          centers_idx <- seq_len(n_variants)
        } else {
          centers_idx <- seq_len(n_variants)
        }
        for (ci in seq_len(n_variants)) {
          if (!is.null(km)) {
            members <- which(km$cluster == ci)
            if (length(members) == 0) next
            center <- km$centers[ci,]
            dists <- rowSums((veg_mat[members, , drop = FALSE] - matrix(center, nrow = length(members), ncol = ncol(veg_mat), byrow = TRUE))^2)
            med_idx <- members[which.min(dists)]
          } else {
            med_idx <- ci
          }
          # Apply sqrt of weights for weighted least squares
          vec <- veg_mat[med_idx, ] * sqrt(pmax(params$weights, 0))
          variants_out[[length(variants_out) + 1]] <- list(vec = vec, id = paste0(v, "_v", length(variants_out) + 1), raw_row = veg_mat[med_idx, , drop = TRUE])
        }
      }

      lib[[v]] <- variants_out
    }

    lib
  }

  build_single_stage_library_weighted <- function(df_train, indices, params, allowed_veg) {
    # Single-stage library: treats barren and all vegetation types as equal endmembers
    # Returns a library with all endmember types (barren + veg types) in one unified structure
    lib <- list()
    expected_cols <- length(indices) * N_TEMPORAL_BINS

    # Build barren endmember library
    barren_data <- df_train |> dplyr::filter(tolower(trimws(as.character(.data$Veg))) == "barren")
    if (nrow(barren_data) > 0) {
      barren_list <- list()
      traces <- unique(barren_data[, c("location_id", "pheno_year")])
      for(i in seq_len(nrow(traces))) {
        sub <- barren_data[barren_data$location_id == traces$location_id[i] &
                           barren_data$pheno_year == traces$pheno_year[i], ]
        mat <- build_pentad_matrix(sub, indices)
        if(!is.null(mat)) {
          vec <- as.numeric(mat)
          barren_list[[length(barren_list) + 1]] <- vec
        }
      }

      if(length(barren_list) > 0) {
        barren_mat <- do.call(rbind, barren_list)
        if(ncol(barren_mat) == expected_cols) {
          # Normalize barren data
          for(k in seq_along(indices)) {
            idx_start <- (k-1)*N_TEMPORAL_BINS + 1
            idx_end <- k*N_TEMPORAL_BINS
            barren_mat[, idx_start:idx_end] <- (barren_mat[, idx_start:idx_end] - params$means[k]) / params$sds[k]
          }
          barren_mat[is.na(barren_mat)] <- 0

          # Create barren variants (same logic as vegetation)
          max_variants <- if (exists("MAX_VARIANTS_PER_VEG", inherits = TRUE)) MAX_VARIANTS_PER_VEG else 7
          n_variants <- min(max_variants, nrow(barren_mat))
          variants_out <- list()

          if (n_variants <= 1) {
            vec <- apply(barren_mat, 2, mean) * sqrt(pmax(params$weights, 0))
            variants_out[[1]] <- list(vec = vec, id = "barren_v1", raw_mat = barren_mat)
          } else {
            set.seed(1)
            km <- tryCatch(kmeans(barren_mat, centers = n_variants, nstart = 5), error = function(e) NULL)
            for (ci in seq_len(n_variants)) {
              if (!is.null(km)) {
                members <- which(km$cluster == ci)
                if (length(members) == 0) next
                center <- km$centers[ci,]
                dists <- rowSums((barren_mat[members, , drop = FALSE] - matrix(center, nrow = length(members), ncol = ncol(barren_mat), byrow = TRUE))^2)
                med_idx <- members[which.min(dists)]
              } else {
                med_idx <- ci
              }
              vec <- barren_mat[med_idx, ] * sqrt(pmax(params$weights, 0))
              variants_out[[length(variants_out) + 1]] <- list(vec = vec, id = paste0("barren_v", length(variants_out) + 1), raw_row = barren_mat[med_idx, , drop = TRUE])
            }
          }
          lib[["barren"]] <- variants_out
        }
      }
    }

    # Build vegetation endmember libraries (same as two-stage approach)
    for(v in allowed_veg) {
      veg_data <- dplyr::filter(df_train, .data$Veg == v)
      if(nrow(veg_data) == 0) next

      veg_list <- list()
      traces <- unique(veg_data[, c("location_id", "pheno_year")])
      for(i in seq_len(nrow(traces))) {
        sub <- veg_data[veg_data$location_id == traces$location_id[i] &
                        veg_data$pheno_year == traces$pheno_year[i], ]

        mat <- build_pentad_matrix(sub, indices)
        if(!is.null(mat)) {
          vec <- as.numeric(mat)
          veg_list[[length(veg_list) + 1]] <- vec
        }
      }

      if(length(veg_list) == 0) next

      veg_mat <- do.call(rbind, veg_list)

      if(ncol(veg_mat) != expected_cols) {
        warning(sprintf("Skipping %s: matrix has %d columns, expected %d", v, ncol(veg_mat), expected_cols))
        next
      }

      for(k in seq_along(indices)) {
        idx_start <- (k-1)*N_TEMPORAL_BINS + 1
        idx_end <- k*N_TEMPORAL_BINS
        veg_mat[, idx_start:idx_end] <- (veg_mat[, idx_start:idx_end] - params$means[k]) / params$sds[k]
      }
      veg_mat[is.na(veg_mat)] <- 0

      max_variants <- if (exists("MAX_VARIANTS_PER_VEG", inherits = TRUE)) MAX_VARIANTS_PER_VEG else 7
      n_variants <- min(max_variants, nrow(veg_mat))
      variants_out <- list()
      if (n_variants <= 1) {
        vec <- apply(veg_mat, 2, mean) * sqrt(pmax(params$weights, 0))
        variants_out[[1]] <- list(vec = vec, id = paste0(v, "_v1"), raw_mat = veg_mat)
      } else {
        set.seed(1)
        km <- tryCatch(kmeans(veg_mat, centers = n_variants, nstart = 5), error = function(e) NULL)
        if (is.null(km)) {
          centers_idx <- seq_len(n_variants)
        } else {
          centers_idx <- seq_len(n_variants)
        }
        for (ci in seq_len(n_variants)) {
          if (!is.null(km)) {
            members <- which(km$cluster == ci)
            if (length(members) == 0) next
            center <- km$centers[ci,]
            dists <- rowSums((veg_mat[members, , drop = FALSE] - matrix(center, nrow = length(members), ncol = ncol(veg_mat), byrow = TRUE))^2)
            med_idx <- members[which.min(dists)]
          } else {
            med_idx <- ci
          }
          vec <- veg_mat[med_idx, ] * sqrt(pmax(params$weights, 0))
          variants_out[[length(variants_out) + 1]] <- list(vec = vec, id = paste0(v, "_v", length(variants_out) + 1), raw_row = veg_mat[med_idx, , drop = TRUE])
        }
      }

      lib[[v]] <- variants_out
    }

    lib
  }

  precompute_optimized_library_weighted <- function(mesma_lib, grid_type = "full") {
    opt_lib <- list()
    
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
        M_norm <- t(apply(M, 1, function(r) {
          nrm <- sqrt(sum(r^2))
          if (!is.finite(nrm) || nrm == 0) return(rep(0, length(r))) else return(r / nrm)
        }))

        opt_lib[[v]] <- list(
          M = M,
          M_norm = M_norm,
          ids = ids
        )
    }
    
    opt_lib
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

  cat("Building separated Stage 1 and Stage 2 libraries with Z-score/PCA/LDA weighting...\n")
  # Ensure df_train has a consistent 'no.soil' column for Stage1 processing
  if (!"no.soil" %in% names(df_train)) {
    if ("no soil" %in% names(df_train)) {
      df_train$`no.soil` <- df_train$`no soil`
      cat("[NOTICE] Mapped 'no soil' -> 'no.soil' in df_train for Stage1 processing\n")
    } else {
      df_train$`no.soil` <- NA_real_
    }
  }

  stage1_data <- dplyr::mutate(df_train, stage1_class = dplyr::if_else(Veg == "barren", "barren",
                                                 dplyr::if_else(!is.na(no.soil) & abs(as.numeric(as.character(no.soil)) - 1) < 0.01, "vegetation", NA_character_)))
  stage1_data <- dplyr::filter(stage1_data, !is.na(stage1_class))
  stage1_data <- dplyr::select(stage1_data, dplyr::all_of(c("location_id", "pheno_year", "date", "doy", "Veg", "stage1_class", avail)))

  stage2_data <- dplyr::filter(df_train, Veg %in% ALLOWED_VEG)
  stage2_data <- dplyr::select(stage2_data, dplyr::all_of(c("location_id", "pheno_year", "date", "doy", "Veg", avail)))

  cat("Training Stage 1 feature pipeline (barren vs vegetation)...\n")
  cat(sprintf("[NOTICE] Using feature set for Stage1 training: %s\n", paste(avail, collapse=", ")))
  cat(sprintf("[Stage1 Weighted] Stage1 data: %d rows with barren=%d, vegetation=%d\n",
              nrow(stage1_data),
              sum(stage1_data$stage1_class == "barren", na.rm=TRUE),
              sum(stage1_data$stage1_class == "vegetation", na.rm=TRUE)))

  STAGE1_PARAMS <- train_feature_pipeline(stage1_data, "stage1_class", avail)
  if (!is.null(STAGE1_PARAMS) && !is.null(STAGE1_PARAMS$weights)) {
    # Ensure tiny non-zero floor to prevent complete zeroing of indices like PPI
    floor_val <- 1e-4
    STAGE1_PARAMS$weights[is.na(STAGE1_PARAMS$weights)] <- 0
    STAGE1_PARAMS$weights[STAGE1_PARAMS$weights < floor_val] <- floor_val
    # No normalization - keep weights at raw LDA-derived scale
    cat(sprintf("[NOTICE] STAGE1_PARAMS weights (unnormalized): mean=%.6f, max=%.6f\n", mean(STAGE1_PARAMS$weights, na.rm = TRUE), max(STAGE1_PARAMS$weights, na.rm = TRUE)))
    print_weights_summary("STAGE1", STAGE1_PARAMS)
  }

  # SINGLE-STAGE UNMIXING: All endmembers (barren + veg types) treated equally
  cat("[MODE] Single-stage unmixing ENABLED (barren and vegetation types treated as equals)\n")

  # Use the same normalization parameters for all endmembers
  SINGLE_STAGE_PARAMS <- STAGE1_PARAMS  # Use Stage1 params for single-stage normalization

  cat("Building single-stage library (barren + all vegetation types)...\n")
  mesma_lib_single <- build_single_stage_library_weighted(df_train, avail, SINGLE_STAGE_PARAMS, ALLOWED_VEG)

  cat("Pre-computing optimized library for single-stage MESMA...\n")
  OPTIMIZED_LIBRARY_SINGLE <- precompute_optimized_library_weighted(mesma_lib_single, grid_type = "full")

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
  assign("TOPK_VARIANTS", 10, envir = globalenv())

  cat(sprintf("[NOTICE] Single-stage feature count: avail=%d, params_indices=%d\n",
              length(avail), length(SINGLE_STAGE_PARAMS$indices)))

  mesma_lib <- mesma_lib_single
  OPTIMIZED_LIBRARY <- OPTIMIZED_LIBRARY_SINGLE

  if (exists("save_mesma_cache")) {
    tryCatch({
      cache_dir <- save_mesma_cache()
      cat("Model cache saved after training pipeline.\n")
    }, error = function(e) {
      cat(sprintf("Failed to write MESMA cache: %s\n", e$message))
    })
  }

  fit_one_task <- function(task_data, compute_stage2_diagnostics = FALSE) {
    if (is.null(task_data) || nrow(task_data) == 0) return(NULL)

    loc <- as.character(task_data$location_id[1])
    yr <- as.integer(task_data$pheno_year[1])

    # Use single-stage parameters
    PARAMS <- SINGLE_STAGE_PARAMS
    raw_mat <- build_pentad_matrix(task_data, PARAMS$indices)
    if (is.null(raw_mat)) {
      if (isTRUE(compute_stage2_diagnostics)) {
        return(list(location_id = loc, pheno_year = yr, n_raw_bins = NA_integer_, n_valid_mask = NA_integer_, y_s2_norm_val = NA_real_, y_s2_masked_mean = NA_real_, weights_s2_fallback_used = NA))
      }
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
      if (isTRUE(compute_stage2_diagnostics)) {
        return(list(location_id = loc, pheno_year = yr, n_raw_bins = length(y_raw), n_valid_mask = 0L, y_s2_norm_val = NA_real_, y_s2_masked_mean = NA_real_, weights_s2_fallback_used = NA))
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

    # Branch based on unmixing mode
    if (!isTRUE(enable_two_stage)) {
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
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_task] loc=%s yr=%d: no top_variants available, skipping\n", loc, yr))
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

      # Build coefficient dataframe with variant-level detail
      coef_df <- data.frame(
        location_id = loc,
        pheno_year = yr,
        Veg = sapply(strsplit(chosen_ids, "_v"), `[`, 1),
        variant_id = chosen_ids,
        coef = coefs,
        rmse = rmse,
        coef_025 = NA,
        coef_975 = NA,
        coef_sd = NA,
        interval = NA,
        inseparable_variant_flag = FALSE,
        inseparable_variant_details = NA_character_,
        stringsAsFactors = FALSE
      )

      # Aggregate by vegetation type (sum coefficients for same veg type)
      coef_agg <- aggregate(coef ~ Veg, data = coef_df, FUN = sum)

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

      return(list(
        coef_df = coef_df,
        diagnostics = diag_df,
        uncertainty = NULL
      ))

      # ===== END SINGLE-STAGE UNMIXING =====
    }
  # These variables are already set in the conditional blocks above (lines 6363-6364 or 6397-6398)
  # No need to reassign them here
  # Single-stage MESMA is used, so mesma_lib and OPTIMIZED_LIBRARY are set accordingly

  required_globals <- list(
    lib = lib,
    raw_lib_templates = raw_lib_templates,
    mesma_lib = mesma_lib,
    avail = avail,
    ALLOWED_VEG = ALLOWED_VEG,
    LOWER_BND = LOWER_BND,
    EPS_SIGMA = EPS_SIGMA,
    MAX_VEG_COMPONENTS = MAX_VEG_COMPONENTS,
    veg_counts = veg_counts,
    OUT_DIR = OUT_DIR,
    prepare_factor_data = prepare_factor_data,
    compress_temporal_matrix = compress_temporal_matrix,
    build_temporal_matrix = build_temporal_matrix,
    build_raw_index_matrix = build_raw_index_matrix,
    compress_stage1_lib_unified = compress_stage1_lib_unified,
    precompute_compressed_templates_unified = precompute_compressed_templates_unified,
    reduce_all_traces_unified = reduce_all_traces_unified,
    build_mesma_variants_unified = build_mesma_variants_unified,
    precompute_templates_from_unified_lib = precompute_templates_from_unified_lib,
    unmix_stage2_compressed = unmix_stage2_compressed,
    unmix_stage1_ols = if (exists("ols_stage1_unmix")) ols_stage1_unmix else NULL,
    stage2_ols_unmix = if (exists("stage2_ols_unmix")) stage2_ols_unmix else NULL,
    evaluate_all_combinations = if (exists("evaluate_all_combinations")) evaluate_all_combinations else NULL,
    COMPRESSED_STAGE1_LIB = if (exists("COMPRESSED_STAGE1_LIB")) COMPRESSED_STAGE1_LIB else NULL,
    COMPRESSED_STAGE1_LIB_WEIGHTED = if (exists("COMPRESSED_STAGE1_LIB_WEIGHTED")) COMPRESSED_STAGE1_LIB_WEIGHTED else NULL,
    STAGE1_LIB = if (exists("STAGE1_LIB")) STAGE1_LIB else NULL,
    cos_sim = cos_sim,
    ENABLE_UNCERTAINTY = ENABLE_UNCERTAINTY,
    ENABLE_QP_SOLVER = ENABLE_QP_SOLVER,
    compute_diagnostics = compute_diagnostics,
    project_to_simplex = project_to_simplex,
    dbg_return_null = dbg_return_null,
    inference_location_ids = if (exists("inference_location_ids")) inference_location_ids else character(0),
    MIN_UNIQUE_DOY_DEFAULT = MIN_UNIQUE_DOY_DEFAULT,
    MIN_UNIQUE_DOY_INFERENCE = MIN_UNIQUE_DOY_INFERENCE,
    build_pentad_matrix = if (exists("build_pentad_matrix")) build_pentad_matrix else NULL,
    ols_stage1_unmix = if (exists("ols_stage1_unmix")) ols_stage1_unmix else NULL,
    solve_weights_nnls_simple = if (exists("solve_weights_nnls_simple")) solve_weights_nnls_simple else NULL,
    safe_mul_vec = if (exists("safe_mul_vec")) safe_mul_vec else NULL,
    build_Z365 = if (exists("build_Z365")) build_Z365 else NULL,
    DEBUG_UNCERTAINTY = DEBUG_UNCERTAINTY,
    INSEPARABLE_VARIANT_INFO = if (exists("INSEPARABLE_VARIANT_INFO")) INSEPARABLE_VARIANT_INFO else NULL,
    OPTIMIZED_LIBRARY = if (exists("OPTIMIZED_LIBRARY")) OPTIMIZED_LIBRARY else NULL,
    TRAINING_NORM_PARAMS = if (exists("TRAINING_NORM_PARAMS")) TRAINING_NORM_PARAMS else NULL,
    apply_stage2_normalization = apply_stage2_normalization,
    INSEPARABLE_VARIANTS = if (exists("INSEPARABLE_VARIANTS")) INSEPARABLE_VARIANTS else NULL,
    VARIANT_SIMILARITY_TABLE = if (exists("VARIANT_SIMILARITY_TABLE")) VARIANT_SIMILARITY_TABLE else NULL,
    STAGE1_PARAMS = if (exists("STAGE1_PARAMS")) STAGE1_PARAMS else NULL,
    TESTING_MODE = if (exists("TESTING_MODE")) TESTING_MODE else FALSE,
    DEBUG = if (exists("DEBUG")) DEBUG else FALSE
  )

  env_task <- list2env(required_globals, parent = globalenv())
  environment(fit_one_task) <- env_task

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
            raw_mat_yr <- build_pentad_matrix(year_data, STAGE1_PARAMS$indices)
            if (!is.null(raw_mat_yr)) {
              y_raw_yr <- as.numeric(raw_mat_yr)

              # Check dimensions match expectations
              n_bins <- N_TEMPORAL_BINS
              expected_length <- length(STAGE1_PARAMS$indices) * n_bins
              if (length(y_raw_yr) != expected_length) {
                if (isTRUE(TESTING_MODE)) {
                  cat(sprintf("[WARN fit_one_location] loc=%s yr=%d: y_raw_yr length mismatch (got %d, expected %d), skipping multi-year processing for this year\n",
                    loc, yr, length(y_raw_yr), expected_length))
                }
              } else {
                # Apply same preprocessing as fit_one_task
                y_s1_yr <- y_raw_yr
                for(k in seq_along(STAGE1_PARAMS$indices)) {
                  idx_start <- (k-1)*n_bins + 1
                  idx_end <- k*n_bins
                  mu <- STAGE1_PARAMS$means[k]
                  sigma <- STAGE1_PARAMS$sds[k]
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

          # REMOVED: nested_two_stage_bootstrap_multiyear is no longer available
          compressed_stage1_lib <- if (exists("COMPRESSED_STAGE1_LIB")) COMPRESSED_STAGE1_LIB else NULL

          if (!is.null(compressed_stage1_lib) && length(comp_templates) > 0 && length(y_vecs_by_year) >= 2) {
            # Build all_results for bootstrap
            all_results <- list()
            for (yr in years) {
              yr_char <- as.character(yr)
              res_yr <- year_results[[yr_char]]
              if (!is.null(res_yr) && !is.null(res_yr$coef_df)) {
                coefs <- setNames(res_yr$coef_df$coef, res_yr$coef_df$Veg)
                y_vec <- y_vecs_by_year[[yr_char]]
                if (!is.null(y_vec) && length(coefs) > 0) {
                  all_results <- c(all_results, list(list(
                    location_id = loc,
                    year = yr,
                    coefs = coefs,
                    y_vec = y_vec
                  )))
                }
              }
            }

            if (length(all_results) >= 2) {
              ci_results <- tryCatch({
                location_bootstrap_nested_mesma(all_results, comp_templates, compressed_stage1_lib, B = BOOTSTRAP_B, seed = 123)
              }, error = function(e) {
                cat(sprintf("[ERROR bootstrap] loc=%s: %s\n", as.character(loc), as.character(e$message)))
                NULL
              })
              if (!is.null(ci_results)) {
                for (yr in years) {
                  yr_char <- as.character(yr)
                  if (!is.null(year_results[[yr_char]])) {
                    year_results[[yr_char]]$uncertainty <- list(coef_ci = ci_results)
                  }
                }
              }
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

  debug_unmix_loc_year <- function(loc, yr, df = NULL, disable_stage1_weighting = FALSE) {
    if (is.null(df)) {
      if (exists('df_tasks_inference') && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0) df <- df_tasks_inference
      else if (exists('df_tasks') && !is.null(df_tasks) && nrow(df_tasks) > 0) df <- df_tasks
      else stop('No task data available in workspace (df_tasks_inference or df_tasks)')
    }
    task_data <- df[df$location_id == loc & df$pheno_year == yr, , drop = FALSE]
    if (nrow(task_data) == 0) stop(sprintf('No rows found for loc=%s pheno_year=%s', as.character(loc), as.character(yr)))

    old_testing <- if (exists('TESTING_MODE', inherits = TRUE)) get('TESTING_MODE', envir = globalenv()) else FALSE
    old_stage1 <- if (exists('STAGE1_WEIGHTING_ENABLED', inherits = TRUE)) get('STAGE1_WEIGHTING_ENABLED', envir = globalenv()) else TRUE
    on.exit({ assign('TESTING_MODE', old_testing, envir = globalenv()); assign('STAGE1_WEIGHTING_ENABLED', old_stage1, envir = globalenv()) })

    assign('TESTING_MODE', TRUE, envir = globalenv())
    if (isTRUE(disable_stage1_weighting)) assign('STAGE1_WEIGHTING_ENABLED', FALSE, envir = globalenv())

    cat(sprintf('\n=== DEBUG UNMIX: loc=%s pheno_year=%s (disable_stage1_weighting=%s) ===\n', as.character(loc), as.character(yr), as.character(disable_stage1_weighting)))
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

  main_processing_block <- function() {
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
  BATCH_SIZE <- 50 # Smaller batches for location-level processing (each has multiple years)

  loc_batches <- split(target_locations, ceiling(seq_along(target_locations) / BATCH_SIZE))
  n_batches <- length(loc_batches)
  pb_width <- min(40L, max(4L, n_batches))

  cat(sprintf("Processing %d locations in %d batches (approx %d locations/batch)...\n",
              n_locs_to_process, length(loc_batches), BATCH_SIZE))

  # Results will be keyed by location, but each contains results for multiple years
  results_by_location <- vector("list", n_locs_to_process)
  names(results_by_location) <- target_locations

  start_time <- Sys.time()

  for (i in seq_along(loc_batches)) {
    batch_locs <- loc_batches[[i]]

    batch_df <- df_tasks[df_tasks$location_id %in% batch_locs, ]

    # Split by location (not location-year)
    batch_location_list <- split(batch_df, batch_df$location_id)

    # Use fit_one_location instead of fit_one_task
    batch_results <- .run_map(batch_location_list, fit_one_location, show_pb = FALSE)

    results_by_location[names(batch_results)] <- batch_results

    if (isTRUE(TESTING_MODE)) {
      for (k in names(batch_results)) {
        loc_result <- batch_results[[k]]
        if (is.null(loc_result)) {
          cat(sprintf("[DEBUG batch_result] location %s returned NULL\n", k))
        } else if (is.list(loc_result)) {
          n_years <- length(loc_result)
          cat(sprintf("[DEBUG batch_result] location %s returned %d year(s)\n", k, n_years))
          for (yr_char in names(loc_result)) {
            r <- loc_result[[yr_char]]
            if (!is.null(r)) {
              ca <- as.numeric(r$vegetated_fraction)
              cb <- as.numeric(r$barren_fraction)
              coef_n <- if (!is.null(r$coef_df) && is.data.frame(r$coef_df)) nrow(r$coef_df) else 0
              has_unc <- !is.null(r$uncertainty)
              cat(sprintf("  Year %s: veg_frac=%.4f barren_frac=%.4f coef_rows=%d has_uncertainty=%s\n",
                          yr_char, ca, cb, coef_n, has_unc))
            }
          }
        }
      }
    }
    
    rm(batch_df, batch_location_list, batch_results)
    gc(verbose = FALSE)
  }
  
  end_time <- Sys.time()
  processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(
    "Main processing loop finished in %.2f seconds (%.2f minutes)\n",
    processing_time, processing_time / 60
  ))

  # Flatten location-based results into year-based results
  cat("Flattening location-based results into year-based results...\n")
  results_list <- list()
  for (loc in names(results_by_location)) {
    loc_results <- results_by_location[[loc]]
    if (!is.null(loc_results) && is.list(loc_results)) {
      for (yr_char in names(loc_results)) {
        task_key <- paste(loc, yr_char, sep = "_")
        results_list[[task_key]] <- loc_results[[yr_char]]
      }
    }
  }

  n_year_results <- length(results_list)
  cat(sprintf("Flattened %d location results into %d year results\n", length(results_by_location), n_year_results))

  if (n_year_results > 0) {
    cat(sprintf("Average time per year-result: %.2f seconds\n", processing_time / n_year_results))
  } else {
    cat("Average time per year-result: N/A (0 results)\n")
  }

  inference_results_list <- list()
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
    
    n_inference_keys <- length(inference_target_keys)

    if (n_inference_keys > 0) {
      cat(sprintf("[MASK-BASED PROCESSING] Using validity masks - processing all tasks with any valid observations (>5%% of bins)\n"))

      inference_key_batches <- split(inference_target_keys, ceiling(seq_along(inference_target_keys) / BATCH_SIZE))
      inference_n_batches <- length(inference_key_batches)
      inference_pb_width <- min(40L, max(4L, inference_n_batches))

      cat(sprintf("Processing %d inference tasks in %d batches...\n", n_inference_keys, length(inference_key_batches)))
      
      inference_results_list <- vector("list", n_inference_keys)
      names(inference_results_list) <- inference_target_keys
      
      inference_start_time <- Sys.time()
      
      for (i in seq_along(inference_key_batches)) {
        batch_keys <- inference_key_batches[[i]]
        
        batch_df <- df_tasks_inference[df_tasks_inference$task_key %in% batch_keys, ]
        
        batch_task_list <- split(batch_df, batch_df$task_key)
        
        batch_results <- .run_map(batch_task_list, fit_one_task, show_pb = FALSE)
        
        inference_results_list[names(batch_results)] <- batch_results
        if (isTRUE(TESTING_MODE)) {
          for (k in names(batch_results)) {
            r <- batch_results[[k]]
            if (is.null(r)) {
              cat(sprintf("[DEBUG batch_result] inference task %s returned NULL\n", k))
            } else {
              ca <- as.numeric(r$vegetated_fraction); cb <- as.numeric(r$barren_fraction)
              coef_n <- if (!is.null(r$coef_df) && is.data.frame(r$coef_df)) nrow(r$coef_df) else 0
              veg_list <- if (!is.null(r$coef_df)) paste(unique(r$coef_df$Veg), collapse = ",") else NA_character_
              diag_vf <- if (!is.null(r$diagnostics) && 'vegetated_fraction' %in% names(r$diagnostics)) r$diagnostics$vegetated_fraction else NA_real_
              cat(sprintf("[DEBUG batch_result] inference task %s returned non-NULL: veg_frac=%.4f barren_frac=%.4f coef_rows=%d vegs=%s diag_vf=%s\n", k, ca, cb, coef_n, veg_list, as.character(diag_vf)))
            }
          }
        }
        
        rm(batch_df, batch_task_list, batch_results)
        gc(verbose = FALSE)
      }
      
      inference_end_time <- Sys.time()
      inference_processing_time <- as.numeric(difftime(inference_end_time, inference_start_time, units = "secs"))
      cat(sprintf("Inference processing finished in %.2f seconds (%.2f minutes)\n",
                  inference_processing_time, inference_processing_time / 60))

      n_null_before_filter <- sum(sapply(inference_results_list, is.null))
      inference_results_list <- inference_results_list[!sapply(inference_results_list, is.null)]
      n_valid <- length(inference_results_list)

      cat(sprintf("Valid inference results: %d out of %d tasks (%.1f%%)\n",
                  n_valid, n_inference_keys, 100*n_valid/n_inference_keys))
      cat(sprintf("Skipped/filtered tasks: %d (%.1f%%)\n",
                  n_null_before_filter, 100*n_null_before_filter/n_inference_keys))
    } else {
      cat("No valid inference location-year pairs found.\n")
    }
  }

  if (length(inference_results_list) > 0) {
    cat(sprintf("Adding %d inference results to results list\n", length(inference_results_list)))
    results_list <- c(results_list, inference_results_list)
  }

  cat("Processing results and writing to Excel files...\n")

  results_list <- results_list[!sapply(results_list, is.null)]
  is_invalid_res <- sapply(results_list, function(res) {
    if (is.null(res)) return(TRUE)
    if (is.null(res$coef_df) || !is.data.frame(res$coef_df)) return(TRUE)
    lid <- res$coef_df$location_id
    if (is.null(lid) || length(lid) == 0) return(TRUE)
    all(is.na(lid) | trimws(as.character(lid)) == "")
  })
  if (any(is_invalid_res)) {
    cat(sprintf("[NOTICE] Filtering out %d invalid result(s) with missing location_id\n", sum(is_invalid_res)))
    results_list <- results_list[!is_invalid_res]
  }
  cat(sprintf("After filtering NULL results: %d results remaining\n", length(results_list)))

  if (length(results_list) > 0) {
    barren_one_count <- sum(sapply(results_list, function(res) {
      if (!is.null(res$barren_fraction) && is.finite(res$barren_fraction)) {
        res$barren_fraction == 1
      } else {
        FALSE
      }
    }))
    barren_one_pct <- barren_one_count / length(results_list) * 100
    if (barren_one_pct > 50) {
      msg <- sprintf("WARNING: %.1f%% of predictions (%d/%d) have barren_fraction = 1, exceeding 50%% threshold. This indicates severe model issues.", barren_one_pct, barren_one_count, length(results_list))
      # Always emit a warning rather than stopping execution. In testing mode, prefix to indicate non-actionable warning.
      if (isTRUE(TESTING_MODE)) {
        warning(paste("[TESTING MODE IGNORE]", msg))
      } else {
        warning(msg)
      }
    }
    cat(sprintf("Barren fraction = 1 in %.1f%% of predictions (%d/%d) - within acceptable limits\n", barren_one_pct, barren_one_count, length(results_list)))
  } else {
    cat("No results to check for barren fraction\n")
  }

  if (length(results_list) == 0) {
    cat("ERROR: All tasks returned NULL results!\n")
    cat("Most likely causes:\n")
    cat("1. Insufficient data quality: Tasks were filtered due to:\n")
    cat("   - Too few observations per location-year (<15 observations)\n")
    cat("   - Insufficient temporal coverage (<25 unique days of year)\n")
    cat("   - Too many missing pentad bins (>85% NA values)\n")
    cat("2. Data filtering issues or missing indices\n")
    cat("3. No valid testing/inference data available (no location-year pairs found)\n")
    cat("\nTo adjust data quality thresholds, modify MIN_OBS, MIN_UNIQUE_DOYS, and MIN_PENTAD_COVERAGE in fit_one_task function.\n")
    stop("No valid results to process")
  }

  if (length(results_list) > 0) {
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

    all_coefs <- tryCatch({
      if (requireNamespace("dplyr", quietly = TRUE)) {
        if (exists("DEBUG_UNCERTAINTY") && isTRUE(DEBUG_UNCERTAINTY)) {
          col_sets <- lapply(seq_along(coef_list), function(i) { list(i = i, cols = names(coef_list[[i]]), ncol = ncol(coef_list[[i]])) })
          uniq <- unique(sapply(col_sets, function(x) paste(sort(x$cols), collapse = ",")))
          cat(sprintf("[DEBUG] combining %d coefficient frames; unique column sets: %d\n", length(col_sets), length(uniq)))
          if (length(uniq) > 1) {
            cat("[DEBUG] Column sets per frame (first 10 shown):\n")
            for (i in seq_len(min(length(col_sets), 10))) {
              cat(sprintf("  frame %d: %d cols: %s\n", col_sets[[i]]$i, col_sets[[i]]$ncol, paste(col_sets[[i]]$cols, collapse = ",")))
            }
          }
        }
        dplyr::bind_rows(coef_list)
      } else {
        do.call(rbind, coef_list)
      }
    }, error = function(e) stop(sprintf("ERROR combining coef_df: %s", e$message))
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

    cat("Combining chosen variant summaries...\n")
    variant_list_pca <- lapply(results_list, function(res) if (!is.null(res$variant_trajectory)) res$variant_trajectory else NULL)
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

    diag_list <- lapply(results_list, function(res) {
      if (!is.null(res$diagnostics)) res$diagnostics else NULL
    })
    diag_list <- diag_list[!sapply(diag_list, is.null)]
    all_diagnostics <- if (length(diag_list) > 0) tryCatch({
      if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required to combine diagnostics robustly")
      dplyr::bind_rows(diag_list)
    }, error = function(e) stop(sprintf("Failed to combine diagnostics: %s", e$message))) else NULL

    q_dvi_data <- tryCatch({
      if (requireNamespace("dplyr", quietly = TRUE)) {
        dplyr::bind_rows(lapply(results_list, function(res) {
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
        do.call(rbind, lapply(results_list, function(res) {
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
      variant_similarity_table <- variant_similarity_table[order(variant_similarity_table$veg, -variant_similarity_table$cos_sim, variant_similarity_table$euclidean_dist), , drop = FALSE]
      if (requireNamespace("dplyr", quietly = TRUE)) {
        variant_similarity_summary <- variant_similarity_table |> 
          dplyr::group_by(.data$veg) |> 
          dplyr::summarise(
            pair_count = dplyr::n(),
            max_cos_sim = if (all(is.na(.data$cos_sim))) NA_real_ else max(.data$cos_sim, na.rm = TRUE),
            min_cos_sim = if (all(is.na(.data$cos_sim))) NA_real_ else min(.data$cos_sim, na.rm = TRUE),
            median_cos_sim = if (all(is.na(.data$cos_sim))) NA_real_ else stats::median(.data$cos_sim, na.rm = TRUE),
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
            max_cos_sim = if (all(is.na(tbl$cos_sim))) NA_real_ else max(tbl$cos_sim, na.rm = TRUE),
            min_cos_sim = if (all(is.na(tbl$cos_sim))) NA_real_ else min(tbl$cos_sim, na.rm = TRUE),
            median_cos_sim = if (all(is.na(tbl$cos_sim))) NA_real_ else stats::median(tbl$cos_sim, na.rm = TRUE),
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
    if (exists("df_train") && !is.null(df_train) && nrow(df_train) > 0) {
      training_locs <- unique(as.character(df_train$location_id))
      df_train_infer <- df_full |> dplyr::filter(location_id %in% training_locs & pheno_year %in% c(2023, 2025))
      if (nrow(df_train_infer) > 0) {
        cat(sprintf("[NOTICE] Running separate (silent) inference on %d training-location rows (years 2023 & 2025)...\n", nrow(df_train_infer)))
        res_train_infer <- run_inference_silent(df_train_infer)
        all_coefs_train_infer <- res_train_infer$all_coefs
        if (is.null(all_coefs_train_infer) || nrow(all_coefs_train_infer) == 0) {
          cat("[WARNING] Training-location inference produced no coefficients; skipping accuracy summary.\n")
        } else {
          # True veg mapping
          true_map <- if (exists("true_veg_map")) true_veg_map else if (exists("gpts_map")) gpts_map |> dplyr::select(location_id, true_veg = Veg) else NULL
          if (is.null(true_map)) {
            cat("[WARNING] No true_veg mapping available; cannot compute veg-type accuracy.\n")
          } else {
            # Merge all predictions with true vegetation labels
            merged <- all_coefs_train_infer |> dplyr::left_join(true_map, by = "location_id")

            # Mark correct predictions (predicted Veg matches true_veg)
            merged$correct <- tolower(merged$Veg) == tolower(merged$true_veg)

            # Calculate fraction-weighted accuracy per location-year
            location_year_accuracy <- merged |>
              dplyr::filter(tolower(true_veg) != "barren") |>  # Only vegetation samples
              dplyr::group_by(location_id, pheno_year, true_veg) |>
              dplyr::summarise(
                total_veg_fraction = sum(coef[tolower(Veg) != "barren"], na.rm = TRUE),  # Total predicted vegetation
                correct_veg_fraction = sum(coef[correct & tolower(Veg) != "barren"], na.rm = TRUE),  # Correctly predicted vegetation
                accuracy = ifelse(total_veg_fraction > 0, 100 * correct_veg_fraction / total_veg_fraction, NA_real_),
                .groups = "drop"
              )

            # Overall vegetation accuracy by year
            summary_by_year_veg_only <- location_year_accuracy |>
              dplyr::group_by(pheno_year) |>
              dplyr::summarise(
                n = dplyr::n(),
                mean_accuracy = mean(accuracy, na.rm = TRUE),
                total_correct_frac = sum(correct_veg_fraction, na.rm = TRUE),
                total_veg_frac = sum(total_veg_fraction, na.rm = TRUE),
                weighted_accuracy = 100 * total_correct_frac / pmax(0.001, total_veg_frac),
                .groups = "drop"
              )

            # Per-vegetation breakdown
            summary_by_veg <- location_year_accuracy |>
              dplyr::group_by(pheno_year, true_veg) |>
              dplyr::summarise(
                n = dplyr::n(),
                mean_accuracy = mean(accuracy, na.rm = TRUE),
                total_correct_frac = sum(correct_veg_fraction, na.rm = TRUE),
                total_veg_frac = sum(total_veg_fraction, na.rm = TRUE),
                weighted_accuracy = 100 * total_correct_frac / pmax(0.001, total_veg_frac),
                .groups = "drop"
              )

            cat("\n=== Training-locations inference accuracy (fraction-weighted) ===\n")
            cat("Accuracy = (correctly predicted veg fraction) / (total predicted veg fraction)\n\n")
            for (r in seq_len(nrow(summary_by_year_veg_only))) {
              cat(sprintf("Year %d: VEGETATION classification accuracy = %.1f%% [mean across %d locations]\n",
                summary_by_year_veg_only$pheno_year[r],
                summary_by_year_veg_only$mean_accuracy[r],
                summary_by_year_veg_only$n[r]))
              cat(sprintf("         Weighted accuracy = %.1f%% (%.3f correct / %.3f total veg)\n",
                summary_by_year_veg_only$weighted_accuracy[r],
                summary_by_year_veg_only$total_correct_frac[r],
                summary_by_year_veg_only$total_veg_frac[r]))
            }
            cat("\nPer-vegetation accuracy breakdown:\n")
            print(summary_by_veg)
            # Save concise summary
            pf_file <- file.path(OUT_DIR, "training_locations_inference_accuracy.csv")
            dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
            write.csv(summary_by_veg, pf_file, row.names = FALSE)
            cat(sprintf("Saved training-locations inference accuracy summary to: %s\n", pf_file))
          }
        }
      } else {
        cat("[NOTICE] No rows found for training-location inference on years 2023/2025.\n")
      }
    }

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

    summary_data <- data.frame(
      Location_ID = unique_locations,
      Total_Years = sapply(unique_locations, function(loc) {
        length(unique(all_coefs$pheno_year[all_coefs$location_id == loc]))
      }),
      Total_Observations = sapply(unique_locations, function(loc) {
        nrow(all_coefs[all_coefs$location_id == loc, ])
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
      openxlsx::writeData(wb, "Variant_Similarity", "Pairwise similarity across variants (cosine similarity and Euclidean distance)", startRow = 1, startCol = 1)
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
    for (res in results_list) {
      if (is.null(res$uncertainty)) next
      loc <- if (!is.null(res$coef_df$location_id)) res$coef_df$location_id[1] else NA_character_
      yr <- if (!is.null(res$coef_df$year)) res$coef_df$year[1] else NA_integer_
      ci <- res$uncertainty$coef_ci
      if (!is.null(ci) && nrow(ci) > 0) {
        ci$location_id <- loc; ci$year <- yr
        unc_coef_rows[[length(unc_coef_rows) + 1]] <- ci[, c("location_id","year","Veg","coef_025","coef_975"), drop = FALSE]
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
      if (length(tv) == 0) tv <- NA_character_

        do.call(rbind, lapply(yrs, function(yr) {
        row <- all_coefs[all_coefs$location_id == loc & all_coefs$pheno_year == yr & tolower(all_coefs$Veg) == tolower(tv), , drop = FALSE]
        pred <- if (nrow(row) == 1) row$coef else NA_real_
        pred_abs <- pred
        rmse_val <- if (nrow(row) == 1 && "rmse" %in% names(row)) row$rmse else NA_real_
          sum_veg_coef <- sum(all_coefs$coef[all_coefs$location_id == loc & all_coefs$pheno_year == yr & tolower(all_coefs$Veg) != "barren"], na.rm = TRUE)
          pred_rel <- if (!is.na(pred) && is.finite(sum_veg_coef) && sum_veg_coef > 0) pred / sum_veg_coef else NA_real_
          pred_rel <- pmin(pred_rel, 1)  # Clamp to 1 to prevent >100%
          abs_pct <- if (!is.na(pred_abs) && !is.na(tv)) abs(1 - pred_abs) * 100 else NA_real_
          abs_pct_rel <- if (!is.na(pred_rel) && !is.na(tv)) abs(1 - pred_rel) * 100 else NA_real_

        data.frame(
          location_id = loc,
          year = yr,
          true_veg = tv,
          pred_coef = pred_rel,
          pred_coef_abs = pred_abs,
          pred_coef_rel = pred_rel,
          rmse = rmse_val,
          abs_pct_diff = abs_pct_rel,
          abs_pct_diff_abs = abs_pct,
          abs_pct_diff_rel = abs_pct_rel,
          stringsAsFactors = FALSE
        )
      }))
    }))

    eval_years <- sort(unique(c(TRAIN_YEARS, TRAIN_YEARS - 1L, TRAIN_YEARS + 1L)))
    eval_years <- eval_years[is.finite(eval_years)]
    best_fit_summary$eval_window <- best_fit_summary$year %in% eval_years
    best_fit_eval <- best_fit_summary[best_fit_summary$eval_window, , drop = FALSE]
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

          quality_metrics$peak_q10_dvi <- peak_q10_dvi
          quality_metrics$peak_q90_dvi <- peak_q90_dvi

          openxlsx::writeData(wb, sheet_name, "QUALITY METRICS", startRow = 1, startCol = 1)
          openxlsx::writeData(wb, sheet_name, quality_metrics, startRow = 2, startCol = 1)

          current_row <- nrow(quality_metrics) + 4

          loc_diag <- if (exists("all_diagnostics") && !is.null(all_diagnostics) && "location_id" %in% names(all_diagnostics)) all_diagnostics[all_diagnostics$location_id == loc_id, , drop = FALSE] else NULL
          if (!is.null(loc_diag) && nrow(loc_diag) > 0) {
            openxlsx::writeData(wb, sheet_name, "DIAGNOSTICS", startRow = current_row, startCol = 1)
            openxlsx::writeData(wb, sheet_name, loc_diag, startRow = current_row + 1, startCol = 1)
            current_row <- current_row + nrow(loc_diag) + 3
          }

          loc_best <- if ("location_id" %in% names(best_fit_summary)) best_fit_summary[best_fit_summary$location_id == loc_id, , drop = FALSE] else data.frame()
          if (!is.null(loc_best) && nrow(loc_best) > 0) {
            desired_cols <- c("year", "true_veg", "pred_coef", "pred_coef_abs", "rmse", "abs_pct_diff", "abs_pct_diff_abs", "eval_window")
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

          if (exists("all_unc_var") && !is.null(all_unc_var) && "location_id" %in% names(all_unc_var)) {
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

    required_cols <- c("location_id", "year", "Veg", "coef")
    missing <- setdiff(required_cols, names(all_coefs))
    if (length(missing) > 0) stop(paste("Missing columns:", paste(missing, collapse = ", ")))

    if (method == "location_bootstrap") {
      # Rename year to pheno_year for location_bootstrap_aggregate
      if ("year" %in% names(all_coefs) && !"pheno_year" %in% names(all_coefs)) {
        all_coefs$pheno_year <- all_coefs$year
      }
      result <- location_bootstrap_aggregate(all_coefs, B = BOOTSTRAP_B)
    } else if (method == "hierarchical") {
      if (!"coef_025" %in% names(all_coefs)) all_coefs$coef_025 <- NA_real_
      if (!"coef_975" %in% names(all_coefs)) all_coefs$coef_975 <- NA_real_

      all_coefs$interval <- all_coefs$coef_975 - all_coefs$coef_025
      all_coefs$interval[!is.finite(all_coefs$interval)] <- NA_real_

      all_coefs$se_proxy <- all_coefs$interval / 3.92
      all_coefs$se_proxy[all_coefs$se_proxy <= 0 | !is.finite(all_coefs$se_proxy)] <- NA_real_
      all_coefs$weight <- 1 / (all_coefs$se_proxy^2)
      all_coefs$weight[!is.finite(all_coefs$weight)] <- 1  # Default weight if no CI

      result <- aggregate_hierarchical(all_coefs)
    } else {
      stop("Unknown method. Use 'location_bootstrap' or 'hierarchical'.")
    }

    result
  }


  aggregate_hierarchical <- function(all_coefs) {

    if (!requireNamespace("lme4", quietly = TRUE)) {
      warning("lme4 not available, falling back to simple mean aggregation")
      # Simple aggregation fallback
      result <- all_coefs |>
        dplyr::group_by(Veg, year) |>
        dplyr::summarize(
          n_locations = dplyr::n(),
          global_coef = mean(coef, na.rm = TRUE),
          se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
          ci_lower = global_coef - 1.96 * se,
          ci_upper = global_coef + 1.96 * se,
          .groups = "drop"
        ) |>
        dplyr::mutate(method = "simple_mean")
      return(result)
    }
    
    veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
    results_list <- list()
    
    for (veg in veg_types) {
      veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), ]
      
      if (nrow(veg_data) == 0) {
        warning(sprintf("No valid coefficients found for vegetation type: %s", veg))
        next
      }
      
      if (nrow(veg_data) < 10) {
        simple_result <- veg_data |> 
          dplyr::group_by(year) |> 
          dplyr::summarize(
            n_locations = dplyr::n(),
            global_coef = mean(coef, na.rm = TRUE),
            se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
            ci_lower = global_coef - 1.96 * se,
            ci_upper = global_coef + 1.96 * se,
            .groups = "drop"
          ) |> 
          dplyr::mutate(Veg = veg, method = "simple_mean")
        
        results_list[[veg]] <- simple_result
        next
      }
      
      tryCatch({
        veg_data$year_factor <- as.factor(veg_data$year)
        
        model <- suppressWarnings(lme4::lmer(coef ~ year_factor + (1|location_id), data = veg_data))
        
        fe <- lme4::fixef(model)
        vcov_fe <- as.matrix(vcov(model))
        
        years <- sort(unique(veg_data$year))
        pred_data <- data.frame(year_factor = as.factor(years))
        
        preds <- predict(model, newdata = pred_data, re.form = NA)
        
        boot_preds <- lme4::bootMer(model, function(m) {
          predict(m, newdata = pred_data, re.form = NA)
        }, nsim = 100, type = "parametric")
        
        ci_lower <- apply(boot_preds$t, 2, quantile, 0.025)
        ci_upper <- apply(boot_preds$t, 2, quantile, 0.975)
        
        hier_result <- data.frame(
          year = years,
          Veg = veg,
          n_locations = sapply(years, function(y) sum(veg_data$year == y)),
          global_coef = preds,
          se = apply(boot_preds$t, 2, sd),
          ci_lower = pmax(0, ci_lower),
          ci_upper = pmin(1, ci_upper),
          method = "hierarchical"
        )
        
        results_list[[veg]] <- hier_result
        
      }, error = function(e) {
        warning(sprintf("Hierarchical model failed for %s: %s. Using simple mean.", veg, e$message))
        simple_result <- veg_data |> 
          dplyr::group_by(year) |> 
          dplyr::summarize(
            n_locations = dplyr::n(),
            global_coef = mean(coef, na.rm = TRUE),
            se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
            ci_lower = pmax(0, global_coef - 1.96 * se),
            ci_upper = pmin(1, global_coef + 1.96 * se),
            .groups = "drop"
          ) |> 
          dplyr::mutate(Veg = veg, method = "simple_mean")
        
        results_list[[veg]] <- simple_result
      })
    }
    
    dplyr::bind_rows(results_list)
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

  if (exists("inference_results_list") && length(inference_results_list) > 0) {
    cat("\nProcessing inference results for separate output...\n")
    
    inference_coef_list <- lapply(inference_results_list, function(res) {
      if (is.null(res$coef_df) || nrow(res$coef_df) == 0) {
        NULL
      } else {
        res$coef_df
      }
    })
    inference_coef_list <- inference_coef_list[!sapply(inference_coef_list, is.null)]
    
    if (length(inference_coef_list) > 0) {
      inference_coefs <- tryCatch({
        dplyr::bind_rows(inference_coef_list)
      }, error = function(e) {
        do.call(rbind, inference_coef_list)
      })
      
      cat(sprintf("Combined %d inference coefficient rows from %d location-years\n", 
                  nrow(inference_coefs), length(inference_coef_list)))
      wb_inference <- openxlsx::createWorkbook()
      
      unique_inference_locations <- unique(inference_coefs$location_id)
      for (loc_id in unique_inference_locations) {
        loc_data <- inference_coefs[inference_coefs$location_id == loc_id, ]
        sheet_name <- paste0("Loc_", loc_id)
        openxlsx::addWorksheet(wb_inference, sheet_name)
        openxlsx::writeData(wb_inference, sheet_name, loc_data)
      }
      
      inference_output_filename <- file.path(OUT_DIR, "inference_results.xlsx")
      openxlsx::saveWorkbook(wb_inference, inference_output_filename, overwrite = TRUE)
      cat(sprintf("Saved inference Excel file to: %s\n", inference_output_filename))
      
      inference_global_pattern <- aggregate_to_global_pattern(inference_coefs, method = "location_bootstrap")
      
      if (nrow(inference_global_pattern) > 0) {
        if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
        if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")
        
        inference_global_pattern_veg <- inference_global_pattern[tolower(inference_global_pattern$Veg) != "barren", ]
        inference_global_pattern_barren <- inference_global_pattern[tolower(inference_global_pattern$Veg) == "barren", ]
        
        veg_max <- suppressWarnings(max(inference_global_pattern_veg$global_coef, inference_global_pattern_veg$ci_upper, na.rm = TRUE))
        barren_max <- suppressWarnings(max(inference_global_pattern_barren$global_coef, inference_global_pattern_barren$ci_upper, na.rm = TRUE))
        if (is.na(veg_max) || veg_max <= 0 || is.na(barren_max) || barren_max <= 0) {
          inference_barren_scale <- 1
        } else {
          inference_barren_scale <- veg_max / barren_max
        }
        if (nrow(inference_global_pattern_barren) > 0) {
          inference_global_pattern_barren$global_coef_scaled <- inference_global_pattern_barren$global_coef * inference_barren_scale
          inference_global_pattern_barren$ci_lower_scaled <- inference_global_pattern_barren$ci_lower * inference_barren_scale
          inference_global_pattern_barren$ci_upper_scaled <- inference_global_pattern_barren$ci_upper * inference_barren_scale
        }
        
        p_inference <- ggplot(inference_global_pattern_veg, aes(x = year, y = global_coef, color = Veg, fill = Veg)) +
          geom_line(linewidth = 1.2) +
          geom_point(size = 2, show.legend = FALSE) +
          geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, color = NA)
        
        if (nrow(inference_global_pattern_barren) > 0) {
          p_inference <- p_inference + 
            geom_line(data = inference_global_pattern_barren, 
                     aes(x = year, y = global_coef_scaled), 
                     color = "brown", linewidth = 1.2, linetype = "dashed") +
            geom_point(data = inference_global_pattern_barren, 
                      aes(x = year, y = global_coef_scaled), 
                      color = "brown", size = 2, show.legend = FALSE) +
            geom_ribbon(data = inference_global_pattern_barren,
                       aes(x = year, ymin = ci_lower_scaled, ymax = ci_upper_scaled),
                       fill = "brown", alpha = 0.1, color = NA)
        }
        
        all_y_values <- c(inference_global_pattern_veg$global_coef, inference_global_pattern_veg$ci_lower, inference_global_pattern_veg$ci_upper)
        if (nrow(inference_global_pattern_barren) > 0) {
          all_y_values <- c(all_y_values, inference_global_pattern_barren$global_coef_scaled, inference_global_pattern_barren$ci_lower_scaled, inference_global_pattern_barren$ci_upper_scaled)
        }
        min_y <- min(all_y_values, na.rm = TRUE)
        max_y <- max(all_y_values, na.rm = TRUE)
        
        p_inference <- p_inference +
          labs(
            title = "Inference Locations: Average Coverage Percentage per Vegetation Type (2020-2025)",
            subtitle = sprintf("Based on %d locations with bootstrap uncertainty", max(inference_global_pattern$n_locations, na.rm = TRUE)),
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
          scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                             limits = c(min_y, max_y),
                             expand = c(0, 0),
                             sec.axis = ggplot2::sec_axis(~ . / inference_barren_scale, name = "Barren Fraction", labels = scales::percent_format(accuracy = 1))) +
          scale_x_continuous(limits = c(2020, 2025)) +
          scale_color_brewer(palette = "Set1") +
          scale_fill_brewer(palette = "Set1")
        
        inference_plot_filename <- file.path(OUT_DIR, "inference_average_coverage_plot.png")
        ggsave(inference_plot_filename, p_inference, width = 10, height = 6, dpi = 300)
        cat(sprintf("Saved inference plot to: %s\n", inference_plot_filename))
        
        if (exists("bootstrap_trend_ci")) {
          cat("Computing trend CI via location bootstrap for inference results...\n")
          inference_trend_ci <- tryCatch({
            bootstrap_trend_ci(inference_coefs, B = 200)
          }, error = function(e) {
            warning(sprintf("bootstrap_trend_ci failed for inference: %s", e$message))
            NULL
          })
          if (!is.null(inference_trend_ci) && nrow(inference_trend_ci) > 0) {
            openxlsx::addWorksheet(wb_inference, "Trend_CI")
            openxlsx::writeData(wb_inference, "Trend_CI", inference_trend_ci)
            cat("Saved inference trend CI to workbook sheet 'Trend_CI'\n")
          }
        }
        
        if (exists("inference_coefs") && nrow(inference_coefs) > 0 && "PPI" %in% names(df_tasks_inference) && any(!is.na(df_tasks_inference$PPI))) {
          inference_pattern_ppi <- location_bootstrap_ppi(inference_coefs, df_tasks_inference, B = BOOTSTRAP_B, seed = 123)
          if (!is.null(inference_pattern_ppi) && nrow(inference_pattern_ppi) > 0) {
            inference_pattern_ppi <- inference_pattern_ppi[!tolower(trimws(inference_pattern_ppi$Veg)) %in% c("barren"), ]
          }
          if (!is.null(inference_pattern_ppi) && nrow(inference_pattern_ppi) > 0) {
            p_inf_ppi_ts <- ggplot(inference_pattern_ppi, aes(x = year, y = global_coef, color = Veg, group = Veg)) +
              geom_line(linewidth = 1) +
              geom_point(show.legend = FALSE) +
              geom_ribbon(aes(ymin = coef_025, ymax = coef_975, fill = Veg), alpha = 0.15, color = NA) +
              labs(title = "Inference: PPI-Normalized Vegetation Fractions Over Time (Location Bootstrap)",
                   x = "Year", y = "Total Normalized Fraction") +
              theme_minimal()
            inf_ppi_plot_filename <- file.path(OUT_DIR, "inference_ppi_normalized_timeseries.png")
            ggsave(inf_ppi_plot_filename, p_inf_ppi_ts, width = 8, height = 6)
            cat(sprintf("Saved inference PPI-normalized time series plot to: %s\n", inf_ppi_plot_filename))
          } else {
            cat("Inference PPI normalization aggregation returned no results.\n")
          }
        }
      }
    } else {
      cat("No inference coefficient data extracted.\n")
    }
  } else {
    cat("No inference results to process.\n")
  }

  cat("\nGenerating average coverage plot...\n")

  tryCatch({
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      install.packages("ggplot2")
    }
    library(ggplot2)

    global_pattern_all <- aggregate_to_global_pattern(all_coefs, method = "location_bootstrap")

    if (is.null(global_pattern_all) || nrow(global_pattern_all) == 0) {
      cat("[WARN] aggregate_to_global_pattern returned no data, skipping average coverage plot\n")
    } else {
      global_pattern_all <- global_pattern_all |> dplyr::filter(year >= 1985 & year <= 2025)
      global_pattern_all_veg <- global_pattern_all[tolower(global_pattern_all$Veg) != "barren", ]
      global_pattern_all_barren <- global_pattern_all[tolower(global_pattern_all$Veg) == "barren", ]

  p <- ggplot(global_pattern_all_veg, aes(x = year, y = global_coef, color = Veg, fill = Veg)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2, show.legend = FALSE) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, color = NA)
  
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

      plot_filename <- file.path(OUT_DIR, "average_coverage_plot.png")
      ggsave(plot_filename, p, width = 10, height = 6, dpi = 300)
      cat(sprintf("Saved average coverage plot to: %s\n", plot_filename))
    }
  }, error = function(e) {
    cat(sprintf("[ERROR] Failed to generate average coverage plot: %s\n", e$message))
    cat("[INFO] Continuing with script execution...\n")
  })

  cat("\nGenerating Observations vs Accuracy plot...\n")
  
  
  loc_accuracy <- best_fit_summary |> 
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




  aggregate_to_global_pattern <- function(all_coefs, method = "location_bootstrap") {

    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")

    required_cols <- c("location_id", "year", "Veg", "coef")
    missing <- setdiff(required_cols, names(all_coefs))
    if (length(missing) > 0) stop(paste("Missing columns:", paste(missing, collapse = ", ")))

    if (method == "location_bootstrap") {
      # Rename year to pheno_year for location_bootstrap_aggregate
      if ("year" %in% names(all_coefs) && !"pheno_year" %in% names(all_coefs)) {
        all_coefs$pheno_year <- all_coefs$year
      }
      result <- location_bootstrap_aggregate(all_coefs, B = BOOTSTRAP_B)
    } else if (method == "hierarchical") {
      if (!"coef_025" %in% names(all_coefs)) all_coefs$coef_025 <- NA_real_
      if (!"coef_975" %in% names(all_coefs)) all_coefs$coef_975 <- NA_real_

      all_coefs$interval <- all_coefs$coef_975 - all_coefs$coef_025
      all_coefs$interval[!is.finite(all_coefs$interval)] <- NA_real_

      all_coefs$se_proxy <- all_coefs$interval / 3.92
      all_coefs$se_proxy[all_coefs$se_proxy <= 0 | !is.finite(all_coefs$se_proxy)] <- NA_real_
      all_coefs$weight <- 1 / (all_coefs$se_proxy^2)
      all_coefs$weight[!is.finite(all_coefs$weight)] <- 1  # Default weight if no CI

      result <- aggregate_hierarchical(all_coefs)
    } else {
      stop("Unknown method. Use 'location_bootstrap' or 'hierarchical'.")
    }

    result
  }


  aggregate_hierarchical <- function(all_coefs) {

    if (!requireNamespace("lme4", quietly = TRUE)) {
      warning("lme4 not available, falling back to simple mean aggregation")
      # Simple aggregation fallback
      result <- all_coefs |>
        dplyr::group_by(Veg, year) |>
        dplyr::summarize(
          n_locations = dplyr::n(),
          global_coef = mean(coef, na.rm = TRUE),
          se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
          ci_lower = global_coef - 1.96 * se,
          ci_upper = global_coef + 1.96 * se,
          .groups = "drop"
        ) |>
        dplyr::mutate(method = "simple_mean")
      return(result)
    }
    
    veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
    results_list <- list()
    
    for (veg in veg_types) {
      veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), ]
      
      if (nrow(veg_data) == 0) {
        warning(sprintf("No valid coefficients found for vegetation type: %s", veg))
        next
      }
      
      if (nrow(veg_data) < 10) {
        simple_result <- veg_data |> 
          dplyr::group_by(year) |> 
          dplyr::summarize(
            n_locations = dplyr::n(),
            global_coef = mean(coef, na.rm = TRUE),
            se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
            ci_lower = global_coef - 1.96 * se,
            ci_upper = global_coef + 1.96 * se,
            .groups = "drop"
          ) |> 
          dplyr::mutate(Veg = veg, method = "simple_mean")
        
        results_list[[veg]] <- simple_result
        next
      }
      
      tryCatch({
        veg_data$year_factor <- as.factor(veg_data$year)
        
        model <- suppressWarnings(lme4::lmer(coef ~ year_factor + (1|location_id), data = veg_data))
        
        fe <- lme4::fixef(model)
        vcov_fe <- as.matrix(vcov(model))
        
        years <- sort(unique(veg_data$year))
        pred_data <- data.frame(year_factor = as.factor(years))
        
        preds <- predict(model, newdata = pred_data, re.form = NA)
        
        boot_preds <- lme4::bootMer(model, function(m) {
          predict(m, newdata = pred_data, re.form = NA)
        }, nsim = 100, type = "parametric")
        
        ci_lower <- apply(boot_preds$t, 2, quantile, 0.025)
        ci_upper <- apply(boot_preds$t, 2, quantile, 0.975)
        
        hier_result <- data.frame(
          year = years,
          Veg = veg,
          n_locations = sapply(years, function(y) sum(veg_data$year == y)),
          global_coef = preds,
          se = apply(boot_preds$t, 2, sd),
          ci_lower = pmax(0, ci_lower),
          ci_upper = pmin(1, ci_upper),
          method = "hierarchical"
        )
        
        results_list[[veg]] <- hier_result
        
      }, error = function(e) {
        warning(sprintf("Hierarchical model failed for %s: %s. Using simple mean.", veg, e$message))
        simple_result <- veg_data |> 
          dplyr::group_by(year) |> 
          dplyr::summarize(
            n_locations = dplyr::n(),
            global_coef = mean(coef, na.rm = TRUE),
            se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
            ci_lower = pmax(0, global_coef - 1.96 * se),
            ci_upper = pmin(1, global_coef + 1.96 * se),
            .groups = "drop"
          ) |> 
          dplyr::mutate(Veg = veg, method = "simple_mean")
        
        results_list[[veg]] <- simple_result
      })
    }
    
    dplyr::bind_rows(results_list)
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
      
      p <- p + 
        ggplot2::geom_line(data = global_pattern_barren, 
                          ggplot2::aes(x = year, y = coef_scaled), 
                          color = "brown", linewidth = 1.2, linetype = "dashed") +
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
        loc_data
      })
      boot_data <- do.call(rbind, boot_data_list)
      
      # Fit a linear mixed-effects model to the bootstrap sample
      # Use boot_id as the random effect to handle repeated locations
      model <- tryCatch({
        lme4::lmer(coef ~ pheno_year + (1|boot_id), data = boot_data)
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
  ggsave(file.path(OUT_DIR, "all_vegetation_trends.png"), p1_all, width = 10, height = 6)



  p3 <- plot_vegetation_stacked_area(global_pattern)
  ggsave(file.path(OUT_DIR, "vegetation_stacked_area.png"), p3, width = 10, height = 6)

  p4 <- plot_vegetation_heatmap(global_pattern)
  ggsave(file.path(OUT_DIR, "vegetation_heatmap.png"), p4, width = 10, height = 6)

  # Use a smaller bootstrap B for interactive/debug runs to avoid long hangs;
  # keep it moderate for production use (can be tuned by setting TREND_BOOT_B)
  TREND_BOOT_B <- if (exists("TREND_BOOT_B")) TREND_BOOT_B else 200
  trend_ci <- bootstrap_trend_ci(all_coefs, B = TREND_BOOT_B)
  # Ensure all ALLOWED_VEG are present in trend_ci (fill missing with NA rows)
  if (exists("ALLOWED_VEG")) {
    missing_vegs <- setdiff(ALLOWED_VEG, trend_ci$Veg)
    if (length(missing_vegs) > 0) {
      cat(sprintf("[NOTICE] No trend data for: %s. Adding placeholder NA rows for reporting.\n", paste(missing_vegs, collapse = ", ")))
      for (mv in missing_vegs) {
        trend_ci <- rbind(trend_ci, data.frame(Veg = mv, slope_mean = NA_real_, slope_median = NA_real_, slope_ci_lower = NA_real_, slope_ci_upper = NA_real_, prob_positive = NA_real_, prob_negative = NA_real_, stringsAsFactors = FALSE))
      }
    }
  }
  print(trend_ci)

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

  trends <- tryCatch(analyze_vegetation_trends(all_coefs, B = 200), error = function(e) { warning(sprintf("analyze_vegetation_trends failed: %s", e$message)); NULL })
  if (is.null(trends)) {
    cat("[NOTICE] No trends computed (insufficient data). Writing empty placeholder to Excel sheet.\n")
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
      ggsave(file.path(OUT_DIR, "ppi_normalized_timeseries.png"), p_ppi_ts, width = 8, height = 6)
      cat(sprintf("Saved PPI-normalized time series plot to: %s\n", file.path(OUT_DIR, "ppi_normalized_timeseries.png")))
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

  openxlsx::addWorksheet(wb, "Trend_Bootstrap_CI")
  openxlsx::writeData(wb, "Trend_Bootstrap_CI", trend_ci)

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

  # Only auto-run the main processing when the script is executed normally.
  # Tests can set MESMA_NO_AUTO_RUN <- TRUE to source this file without executing.
  if (!exists("MESMA_NO_AUTO_RUN") || !isTRUE(MESMA_NO_AUTO_RUN)) {
    cat("\n=== Starting MESMA processing ===\n")

    result <- tryCatch(
      {
        main_processing_block()
        list(success = TRUE)
      },
      error = function(e) {
        cat(sprintf("\n========================================\n"))
        cat(sprintf("FATAL ERROR: %s\n", e$message))
        cat(sprintf("========================================\n\n"))


        cat("Stack trace:\n")

        try({
          dump.frames("mesma_error_dump", to.file = TRUE)
          cat("\nDiagnostic frames dumped to 'mesma_error_dump.rda'\n")
          cat("To inspect: load('mesma_error_dump.rda'); debugger(mesma_error_dump)\n")
          flush.console()
        }, silent = TRUE)

        cat("\nScript terminated due to error.\n")
        flush.console()

        list(success = FALSE, error = e$message)
      },
      warning = function(w) {
        cat(sprintf("\n[WARNING] %s\n", w$message))
        flush.console()
        # Treat warnings as non-fatal here; do not attempt to invoke a muffle restart
        invisible(NULL)
      }
    )

    if (!is.null(result) && !isTRUE(result$success)) {
      cat("\n========================================\n")
      cat("[SCRIPT FAILED] Check error messages above\n")
      cat("========================================\n")

    } else {
      cat("\n========================================\n")
      cat("[SCRIPT COMPLETED SUCCESSFULLY]\n")
      cat("========================================\n")

    }
  }

cat("[CLEANUP] Finalizing outputs and freeing temporary memory...\n")
# Final cleanup of large temporaries before exiting
tryCatch({
  cleanup_memory(verbose = TRUE)
}, error = function(e) {
  cat(sprintf("[CLEANUP] Final cleanup failed: %s\n", e$message))
})
}
}

cat("\n=== Script execution finished ===\n")
