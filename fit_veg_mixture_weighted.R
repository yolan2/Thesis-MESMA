library(zoo)
library(dplyr)

suppressPackageStartupMessages({
  suppressWarnings({
    library(zoo)
    library(dplyr)
  })
})
# Set up progressr handlers for ETA display
if(requireNamespace("progressr", quietly = TRUE)) {
  progressr::handlers(global = TRUE)
}

# Canonical optimal indices (required - no fallback alternatives)
OPTIMAL_INDICES <- c(
  "DVI",        # Difference Vegetation Index
  "OSAVI",      # Optimized Soil-Adjusted Vegetation Index
  "MCARI",      # Modified Chlorophyll Absorption in Reflectance Index
  "CRI",        # Carotenoid Reflectance Index
  "PRI",        # Photochemical Reflectance Index
  "NIRv",       # NIR*NDVI vigor
  "PSRI",       # Plant Senescence Reflectance Index
  "NBR",        # Burn ratio (NIR-SWIR2)/(NIR+SWIR2)
  "TCW",        # Total Column Water
  "TCG",        # Triangular Chlorophyll Index (green-based)
  "MNDWI"       # Modified Normalized Difference Water Index
)

# Function to calculate robust variance using STL decomposition + MAD
## Helper: compute MAD^2 with minimal sample requirement (top-level)
compute_mad2 <- function(x, min_samples = 3) {
  x <- x[is.finite(x)]
  if (length(x) < min_samples) return(NA_real_)
  m <- mad(x, na.rm = TRUE, constant = 1.4826)
  if (!is.finite(m)) return(NA_real_)
  m^2
}

## STL-only calc_moving_var: removes loess fallback and nested safe_mad2
calc_moving_var <- function(df, index_name, window = 14, span_loess = 0.1, min_obs_loess = 6) {
  if(!"date" %in% names(df)) stop("calc_moving_var: df must have a 'date' column")
  if(!index_name %in% names(df)) stop(sprintf("calc_moving_var: index '%s' not found in df", index_name))

  n <- nrow(df)
  out <- rep(NA_real_, n)

  process_series_stl <- function(dates, vals) {
    if (length(vals) < 1) return(rep(NA_real_, length(vals)))
    finite_idx <- which(is.finite(vals) & !is.na(dates))
    if (length(finite_idx) == 0) return(rep(NA_real_, length(vals)))

    freq <- 365
    min_len_for_stl <- 2 * freq
    if (length(vals) < min_len_for_stl) return(rep(NA_real_, length(vals)))

    dts <- as.Date(dates)
    full_days <- seq(from = min(dts, na.rm = TRUE), to = max(dts, na.rm = TRUE), by = "day")
    if (length(full_days) < min_len_for_stl) return(rep(NA_real_, length(vals)))
    full_vec <- rep(NA_real_, length(full_days))
    pos_map <- match(dts, full_days)
    full_vec[pos_map] <- as.numeric(vals)

    if (sum(is.finite(full_vec)) < max(10, floor(0.1 * length(full_vec)))) return(rep(NA_real_, length(vals)))

    idx_finite <- which(is.finite(full_vec))
    full_vec_interp <- full_vec
    if (any(!is.finite(full_vec_interp))) {
      interp_vals <- tryCatch({ approx(x = idx_finite, y = full_vec[idx_finite], xout = seq_along(full_vec), rule = 2)$y }, error = function(e) stop(sprintf("calc_moving_var: interpolation failed: %s", e$message)))
      if (is.null(interp_vals) || any(!is.finite(interp_vals))) return(rep(NA_real_, length(vals)))
      full_vec_interp <- interp_vals
    }

    stl_out <- tryCatch({ stats::stl(ts(full_vec_interp, frequency = freq), s.window = "periodic", robust = TRUE) }, error = function(e) stop(sprintf("calc_moving_var: STL decomposition failed: %s", e$message)))
    resid <- as.numeric(stl_out$time.series[, "remainder"])
    if (length(resid) < window) return(rep(NA_real_, length(vals)))
    rv <- zoo::rollapply(resid, width = window, FUN = function(x) compute_mad2(x, min_samples = max(3, floor(window/2))), fill = NA, align = "center")
    pos_map_res <- pos_map
    res_for_rows <- rv[pos_map_res]
    as.numeric(res_for_rows)
  }

  if("location_id" %in% names(df)) {
    locs <- unique(df$location_id)
    for(loc in locs) {
      rows <- which(df$location_id == loc)
      if(length(rows) < 1) next
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

# -----------------------------------------------------------------
# User-configurable defaults (placed early so you can override by
# setting variables in your environment before sourcing this file)
# -----------------------------------------------------------------
OUTPUT_DIR <- "phenology_results"
INPUT_CSV <- file.path(OUTPUT_DIR, "hls_phenology_data.csv")



# --- Strict validation: enforce extractor/predictor output format only ---
if(!file.exists(INPUT_CSV)) stop(paste0("Required input CSV not found: ", INPUT_CSV, ". Run the predictor and then R_extract_hls.R to generate this file."))
raw_df <- tryCatch({ readr::read_csv(INPUT_CSV, show_col_types = FALSE) }, error = function(e) stop(paste0("Failed to read INPUT_CSV: ", e$message)))
if(nrow(raw_df) < 2) stop("INPUT_CSV contains fewer than 2 rows; the fitter requires multi-row per-location time series produced by the predictor/extractor.")
if(!'location_id' %in% names(raw_df) || !'prediction_date' %in% names(raw_df) && !'date' %in% names(raw_df)) stop("INPUT_CSV must contain 'location_id' and 'prediction_date' (or 'date') columns produced by the extractor.")

# Ensure all OPTIMAL_INDICES are present in the CSV (fail hard)
missing_idx <- setdiff(OPTIMAL_INDICES, names(raw_df))
if(length(missing_idx) > 0) stop(paste0('INPUT_CSV missing required indices: ', paste(missing_idx, collapse=", "), '. Run R_extract_hls.R to compute these indices from predictor outputs.'))

# Replace df with the validated raw_df for the rest of the script
df <- raw_df

# Performance optimization: Reduce PCA complexity if needed
if(!exists("PCA_VARIANCE_THRESHOLD")) PCA_VARIANCE_THRESHOLD <- 0.70  # Default 70%
## MAX_PCA_COMPONENTS default moved to consolidated tunables at top

# Fast variance mode: use rolling MAD instead of per-row STL (much faster, approximate)
## FAST_VAR default moved to consolidated tunables at top

# Add timing for major operations
timing_info <- list()
timing_info$start_time <- Sys.time()

# Helper: canonical location_id creation from lon/lat (centralize formatting)
make_location_id <- function(lon, lat) {
  if (length(lon) == 1 && length(lat) == 1) {
    lon_num <- as.numeric(lon); lat_num <- as.numeric(lat)
    if (!is.finite(lon_num) || !is.finite(lat_num)) return(NA_character_)
    sprintf("L_%0.4f_%0.4f", round(lon_num, 4), round(lat_num, 4))
  } else {
    sapply(seq_along(lon), function(i) {
      lon_i <- as.numeric(lon[i]); lat_i <- as.numeric(lat[i])
      if (!is.finite(lon_i) || !is.finite(lat_i)) NA_character_ else sprintf("L_%0.4f_%0.4f", round(lon_i,4), round(lat_i,4))
    }, USE.NAMES = FALSE)
  }
}

cat("Starting vegetation mixture analysis with ROBUST IMPROVEMENTS...\n")
cat("✓ Robust outlier detection using MAD (Median Absolute Deviation)\n")
cat("✓ Correlation-based index filtering (>95% correlation removed)\n")
cat("✓ MAD-based variance estimation for improved robustness\n")
cat("\n")
cat(sprintf("Dataset size: %d rows, %d columns\n", nrow(df), ncol(df)))
cat(sprintf("Number of locations: %d\n", length(unique(df$location_id))))
cat(sprintf("Date range: %s to %s\n", min(df$date, na.rm=TRUE), max(df$date, na.rm=TRUE)))

# -----------------------------------------------------------------
# Consolidated tunable parameters (edit here before sourcing the script)
# These defaults are set with if(!exists(...)) so you can override them
# by setting environment variables or defining them prior to sourcing.
# -----------------------------------------------------------------

# I/O and dataset selection
if(!exists("OUTPUT_DIR")) OUTPUT_DIR <- "phenology_results"
if(!exists("INPUT_CSV")) INPUT_CSV <- file.path(OUTPUT_DIR, "hls_phenology_data.csv")
if(!exists("OUT_DIR")) OUT_DIR <- file.path(OUTPUT_DIR, "veg_mixture_fit")

# Train/test defaults
if(!exists("TEST_YEARS")) TEST_YEARS <- NULL
if(!exists("TRAIN_YEARS")) TRAIN_YEARS <- 2019:2024

# PCA / factor projection tuning
if(!exists("PCA_VARIANCE_THRESHOLD")) PCA_VARIANCE_THRESHOLD <- 0.70
if(!exists("MAX_PCA_COMPONENTS")) MAX_PCA_COMPONENTS <- 20L
if(!exists("MAX_FACTORS_CAP")) MAX_FACTORS_CAP <- 20L

# General algorithm toggles
if(!exists("FAST_VAR")) FAST_VAR <- TRUE
if(!exists("FACTOR_MODE")) FACTOR_MODE <- FALSE
if(!exists("SHAPE_LDA_ENABLE")) SHAPE_LDA_ENABLE <- TRUE
if(!exists("ENABLE_KERNEL_PCA")) ENABLE_KERNEL_PCA <- FALSE

# Ridge / regularization
## RIDGE parameters (dynamic ridge computation requires these to be set externally)
## Provide RIDGE_K_MAIN, RIDGE_MIN_MAIN, RIDGE_MAX_MAIN, RIDGE_K_BOOT, RIDGE_MIN_BOOT, RIDGE_MAX_BOOT
## before sourcing this script. compute_dynamic_ridge(...) will use them to compute per-problem ridge penalties.

# Bootstrap settings
if(!exists("BOOTSTRAP_B")) BOOTSTRAP_B <- 1L
if(!exists("BOOT_MIN_REPS_PER_VEG")) BOOT_MIN_REPS_PER_VEG <- 10L
if(!exists("MIN_OBS_FOR_BOOT")) MIN_OBS_FOR_BOOT <- 5L

# Optimization/solver defaults
if(!exists("VARIANCE_THRESHOLD")) VARIANCE_THRESHOLD <- 0.90
if(!exists("MAX_VEG_COMPONENTS")) MAX_VEG_COMPONENTS <- 8
if(!exists("GAM_K_MAX")) GAM_K_MAX <- 40
if(!exists("GAM_GAMMA")) GAM_GAMMA <- 1.0

# Index selection and prefiltering
if(!exists("USE_INDICES_MIN")) USE_INDICES_MIN <- 1L
if(!exists("MIN_INDEX_SD")) MIN_INDEX_SD <- 0.05

# Sample-balancing and augmentation
if(!exists("ENABLE_SAMPLE_BALANCING")) ENABLE_SAMPLE_BALANCING <- TRUE

# Vegetation whitelist and string handling
if(!exists("ALLOWED_VEG")) ALLOWED_VEG <- c("populus","tamarix","phragmites")

# Numeric safety constants
if(!exists("EPS_SIGMA")) EPS_SIGMA <- 1e-8
if(!exists("LOWER_BND")) LOWER_BND <- 0

# DOY / index presence tuning
if(!exists("MIN_IDX_PRESENCE")) MIN_IDX_PRESENCE <- 0.5

# Parallel / progress settings
if(!exists("PROGRESS_EVERY_TASK")) PROGRESS_EVERY_TASK <- 25
if(!exists("PROGRESS_LOG_TO_FILE")) PROGRESS_LOG_TO_FILE <- TRUE
if(!exists("PROGRESS_BAR")) PROGRESS_BAR <- TRUE
if(!exists("PARALLEL_ENABLE")) PARALLEL_ENABLE <- TRUE
if(!exists("PARALLEL_WORKERS")) PARALLEL_WORKERS <- tryCatch({ if(requireNamespace("parallel", quietly = TRUE)) max(1L, parallel::detectCores(logical = TRUE) - 1L) else 1L }, error = function(...) 1L)

# Misc defaults (kept for backwards compatibility)
if(!exists("BOOTSTRAP_B")) BOOTSTRAP_B <- 1L

# ---- Train/Test Split Implementation ----
cat("\n=== TRAIN/TEST SPLIT CONFIGURATION ===\n")
cat(sprintf("Training years: %s\n", paste(TRAIN_YEARS, collapse=", ")))
if(!is.null(TEST_YEARS)) {
  cat(sprintf("Testing years: %s\n", paste(TEST_YEARS, collapse=", ")))
} else {
  cat("Testing years: All available years\n")
}

# Create training and testing datasets
df$year <- lubridate::year(df$date)

# Training data: specified years
df_train <- df[df$year %in% TRAIN_YEARS, , drop = FALSE]
cat(sprintf("Training dataset: %d rows from %d locations\n",
            nrow(df_train), length(unique(df_train$location_id))))

# Testing data: either specified years or all data
if(!is.null(TEST_YEARS)) {
  df_test <- df[df$year %in% TEST_YEARS, , drop = FALSE]
} else {
  df_test <- df  # Use all data for testing
}
cat(sprintf("Testing dataset: %d rows from %d locations\n",
            nrow(df_test), length(unique(df_test$location_id))))

# Store original df for later use, use training data for library construction
df_full <- df
df <- df_train

cat("Using training data for vegetation library construction\n")
cat("=====================================\n\n")



# Gradient descent solver with sum-to-one constraint
# Helper: project a vector onto the probability simplex (non-negative, sum to 1)
project_to_simplex <- function(v) {
  # Implementation of the algorithm from: Wang & Carreira-Perpinan (2013)
  if(!is.numeric(v)) stop("project_to_simplex: input must be numeric")
  n <- length(v)
  if(n == 0) return(v)
  u <- sort(v, decreasing = TRUE)
  cssv <- cumsum(u)
  rho <- max(which(u - (cssv - 1) / seq_along(u) > 0))
  theta <- (cssv[rho] - 1) / rho
  w <- pmax(v - theta, 0)
  if(sum(w) <= 0) return(rep(1/n, n))
  w / sum(w)
}

# Smart weight initialization using spectral angle-like similarity
initialize_weights_smart <- function(observation, library_signatures) {
  # observation: numeric vector (observed factor vector)
  # library_signatures: list or matrix of candidate signature vectors
  if(is.null(observation) || length(observation) == 0) stop("initialize_weights_smart: empty observation")
  # Accept either list of vectors or a matrix (cols = signatures)
  if(is.matrix(library_signatures)) {
    sigs <- lapply(seq_len(ncol(library_signatures)), function(i) library_signatures[, i])
  } else if(is.list(library_signatures)) {
    sigs <- library_signatures
  } else {
    stop("library_signatures must be a matrix or list of numeric vectors")
  }

  similarities <- sapply(sigs, function(sig) {
    if(!is.numeric(sig) || length(sig) != length(observation)) return(0)
    denom <- sqrt(sum(observation^2, na.rm = TRUE)) * sqrt(sum(sig^2, na.rm = TRUE))
    if(!is.finite(denom) || denom == 0) return(0)
    cos_sim <- sum(observation * sig, na.rm = TRUE) / denom
    max(0, cos_sim)
  })

  if(all(similarities <= 0)) {
    n <- length(similarities)
    return(rep(1/n, n))
  }
  weights <- similarities^2
  weights / sum(weights)
}

# Adaptive gradient-descent solver with momentum, adaptive learning rate,
# early stopping (patience), and optional smart initialization.
rsolve_mixture_adaptive <- function(B_w, y_w, init = NULL, library_signatures = NULL,
                                    max_iter = 5000, tol = 1e-8, init_lr = 0.1,
                                    beta = 0.9, max_patience = 50, lambda2 = 0) {
  n_veg <- ncol(B_w)
  if(is.null(colnames(B_w))) colnames(B_w) <- paste0("V", seq_len(n_veg))

  # Initialize weights
  if(!is.null(init) && length(init) == n_veg) {
    weights <- as.numeric(init)
    weights <- weights / sum(weights)
  } else if(!is.null(library_signatures) && !is.null(y_w) && length(y_w) > 0) {
    # If provided, compute smart initialization in factor-space using y_w as observation
    # y_w here should correspond to the observed factor vector; when not available
    # fall back to uniform init
    try({
      weights <- initialize_weights_smart(as.numeric(y_w), library_signatures)
    }, silent = TRUE)
    if(!exists("weights") || length(weights) != n_veg || sum(weights) <= 0) {
      weights <- rep(1/n_veg, n_veg)
    }
  } else {
    weights <- rep(1/n_veg, n_veg)
  }

  learning_rate <- init_lr
  momentum <- rep(0, n_veg)
  best_weights <- weights
  best_loss <- Inf
  patience <- 0
  grad_prev <- rep(0, n_veg)

  # Optional ridge augmentation to B_w/y_w to implement lambda2
  if(is.finite(lambda2) && lambda2 > 0) {
    B_aug <- rbind(B_w, sqrt(lambda2) * diag(n_veg))
    y_aug <- c(y_w, rep(0, n_veg))
  } else {
    B_aug <- B_w; y_aug <- y_w
  }

  for(iter in seq_len(max_iter)) {
    pred <- as.numeric(B_aug %*% weights)
    residuals <- as.numeric(y_aug - pred)
    loss <- sum(residuals^2, na.rm = TRUE)

    # Compute gradient of SSE
    grad <- -2 * as.numeric(t(B_aug) %*% residuals)

    # Momentum update (exponential moving average of gradients)
    momentum <- beta * momentum + (1 - beta) * grad

    # Detect oscillation: gradient dot previous gradient negative indicates sign change
    if(iter > 1 && sum(grad * grad_prev, na.rm = TRUE) < 0) {
      learning_rate <- learning_rate * 0.9
    }

    # Parameter update and projection to simplex
    weights_new <- weights - learning_rate * momentum
    weights_new <- project_to_simplex(weights_new)

    # Evaluate improvement
    if(is.finite(loss) && loss < best_loss - tol) {
      best_loss <- loss
      best_weights <- weights_new
      patience <- 0
    } else {
      patience <- patience + 1
      if(patience > max_patience) break
    }

    # Convergence by parameter change
    if(max(abs(weights_new - weights), na.rm = TRUE) < tol) {
      weights <- weights_new
      break
    }

    weights <- weights_new
    grad_prev <- grad
  }

  # Return in the same shape as original solver
  if(any(!is.finite(best_weights)) || sum(best_weights) <= 0) {
    best_weights <- rep(1/n_veg, n_veg)
  }
  names(best_weights) <- colnames(B_w)
  return(list(X = as.numeric(best_weights)))
}

# Alias for backward compatibility: make rsolve_mixture_adaptive a drop-in replacement
# for the original solve_mixture_gd name used elsewhere in the codebase.
solve_mixture_gd <- function(B_w, y_w, base_mask = rep(TRUE, ncol(B_w)), lambda2 = 0,
                             learning_rate = 0.01, max_iter = 1000, tol = 1e-6) {
  # Map previous arguments to the adaptive solver signature with reasonable defaults.
  # base_mask is ignored by the adaptive solver (kept for compatibility).
  res <- rsolve_mixture_adaptive(B_w = B_w, y_w = y_w, init = NULL, library_signatures = NULL,
                                 max_iter = max_iter, tol = tol, init_lr = learning_rate,
                                 beta = 0.9, max_patience = 50, lambda2 = lambda2)
  return(res)
}

# GPU-backed solver prototype using torch (optional). This is a low-risk, additive
# function: it tries to use torch on GPU if available, otherwise falls back to CPU
# (so it won't break environments without CUDA). It parameterizes the simplex via
# logits -> softmax and uses Adam to minimize SSE. Returns a list(X = weights)
torch_solve_mixture_gpu <- function(B_w, y_w, init = NULL, library_signatures = NULL,
                                    max_iter = 5000, tol = 1e-8, init_lr = 0.1,
                                    beta = 0.9, max_patience = 50, lambda2 = 0,
                                    device = NULL, verbose = FALSE) {
  # Mirror the CPU rsolve_mixture_adaptive algorithm as closely as possible
  # but perform matrix multiplies on torch tensors to allow GPU acceleration.
  if(!requireNamespace("torch", quietly = TRUE)) {
    if(verbose) message("torch not installed - falling back to rsolve_mixture_adaptive")
    return(rsolve_mixture_adaptive(B_w = B_w, y_w = y_w, init = init, library_signatures = library_signatures,
                                   max_iter = max_iter, tol = tol, init_lr = init_lr,
                                   beta = beta, max_patience = max_patience, lambda2 = lambda2))
  }
  library(torch)

  # Prepare numeric inputs
  B_w_num <- as.matrix(B_w)
  y_w_num <- as.numeric(y_w)
  n_veg <- ncol(B_w_num)
  if(n_veg == 0) return(list(X = numeric(0)))

  # Initialize weights exactly like CPU solver
  if(!is.null(init) && length(init) == n_veg) {
    weights <- as.numeric(init); weights <- weights / sum(weights)
  } else if(!is.null(library_signatures) && !is.null(y_w) && length(y_w) > 0) {
    try({ weights <- initialize_weights_smart(as.numeric(y_w), library_signatures) }, silent = TRUE)
    if(!exists("weights") || length(weights) != n_veg || sum(weights) <= 0) weights <- rep(1/n_veg, n_veg)
  } else {
    weights <- rep(1/n_veg, n_veg)
  }

  learning_rate <- init_lr
  momentum <- rep(0, n_veg)
  best_weights <- weights
  best_loss <- Inf
  patience <- 0
  grad_prev <- rep(0, n_veg)

  # Ridge augmentation exactly like CPU solver
  if(is.finite(lambda2) && lambda2 > 0) {
    B_aug_num <- rbind(B_w_num, sqrt(lambda2) * diag(n_veg))
    y_aug_num <- c(y_w_num, rep(0, n_veg))
  } else {
    B_aug_num <- B_w_num; y_aug_num <- y_w_num
  }

  # Prepare torch tensors for matrix multiplies (on chosen device)
  if(is.null(device)) device <- if(torch::cuda_is_available()) torch_device("cuda") else torch_device("cpu")
  B_aug_t <- torch_tensor(B_aug_num, dtype = torch_float(), device = device)
  # y_aug will be used as numeric R vector for SSE computations (to match CPU path)

  for(iter in seq_len(max_iter)) {
    # Compute pred and residuals using torch matmul but convert to numeric for identical CPU arithmetic where needed
    # pred = B_aug %*% weights
    w_t <- torch_tensor(weights, dtype = torch_float(), device = device)
    pred_t <- B_aug_t$matmul(w_t$unsqueeze(2))$squeeze()
    pred_num <- as.numeric(pred_t$to(device = torch_device("cpu"))$cpu()$numpy())
    residuals <- as.numeric(y_aug_num - pred_num)
    loss <- sum(residuals^2, na.rm = TRUE)

    # Compute gradient exactly as CPU: grad = -2 * t(B_aug) %*% residuals
    grad <- -2 * as.numeric(t(B_aug_num) %*% residuals)

    # Momentum update
    momentum <- beta * momentum + (1 - beta) * grad

    # Oscillation detection: gradient dot previous gradient negative
    if(iter > 1 && sum(grad * grad_prev, na.rm = TRUE) < 0) {
      learning_rate <- learning_rate * 0.9
    }

    # Parameter update and projection to simplex (identical to CPU)
    weights_new <- weights - learning_rate * momentum
    weights_new <- project_to_simplex(weights_new)

    # Evaluate improvement and patience
    if(is.finite(loss) && loss < best_loss - tol) {
      best_loss <- loss
      best_weights <- weights_new
      patience <- 0
    } else {
      patience <- patience + 1
      if(patience > max_patience) break
    }

    # Convergence by parameter change
    if(max(abs(weights_new - weights), na.rm = TRUE) < tol) {
      weights <- weights_new
      break
    }

    weights <- weights_new
    grad_prev <- grad
  }

  if(any(!is.finite(best_weights)) || sum(best_weights) <= 0) best_weights <- rep(1/n_veg, n_veg)
  names(best_weights) <- colnames(B_w_num)
  return(list(X = as.numeric(best_weights)))
}

# Optimize weights using unconstrained parameters mapped to the simplex via softmax
# Uses Nelder-Mead on the unconstrained params to minimize weighted SSE: sum(w_all * (target - F %*% weights)^2)
optimize_weights_simplex <- function(F_matrix, w_all, target = rep(1, nrow(F_matrix)), init = NULL, maxit = 2000) {
  nveg <- ncol(F_matrix)
  if(is.null(init)) {
    # Initialize with inverse-variance weighting across vegetation types
    # to start from a more balanced position than equal weights
    col_vars <- apply(F_matrix, 2, function(col) var(col, na.rm = TRUE))
    col_vars[!is.finite(col_vars) | col_vars <= 0] <- 1
    init_weights <- 1 / col_vars
    init <- init_weights / sum(init_weights)
  }
  # map params -> weights via softmax
  softmax <- function(p) {
    ex <- exp(p - max(p))
    ex / sum(ex)
  }

  # map softmax output into bounded-simplex with per-coordinate lower bound LOWER_BND
  # IMPORTANT: do not force a per-coordinate LOWER_BND here.
  # Treat LOWER_BND as a post-process threshold (entries below are zeroed and remaining renormalized).
  bounded_from_softmax <- function(pvec) {
    s <- softmax(pvec)
    # return softmax directly (standard simplex). Iterative thresholding will handle the LOWER_BND behavior.
    s / sum(s)
  }

  # objective: weighted SSE against target
  obj <- function(p) {
    w <- softmax(p)
    pred <- as.numeric(F_matrix %*% w)
    res <- target - pred
    sse <- sum(w_all * (res^2), na.rm = TRUE)
    if(!is.finite(sse)) sse <- .Machine$double.xmax
    sse
  }

  # initialize in parameter space (map init back to unconstrained params)
  # Robust initialization: add epsilon, re-normalize, then log
  # This avoids log(0) = -Inf which stalls the optimizer.
  eps <- 1e-9
  init_safe <- init + eps
  init_safe <- init_safe / sum(init_safe)
  p0 <- log(init_safe)
  res <- tryCatch({
    stats::optim(par = p0, fn = obj, method = "Nelder-Mead", control = list(maxit = maxit, reltol = 1e-8))
  }, error = function(e) stop(sprintf("optimize_weights_simplex: Nelder-Mead optimization failed: %s", e$message)))

  if(is.null(res) || !is.list(res) || is.null(res$par)) {
    return(init)
  }
  w_opt <- bounded_from_softmax(res$par)
  # Ensure non-negative and sums to one (numeric safety)
  w_opt[w_opt < 0] <- 0
  if(!is.finite(sum(w_opt)) || sum(w_opt) <= 0) w_opt <- rep(1/length(w_opt), length(w_opt))
  w_opt <- w_opt / sum(w_opt)
  w_opt
}


# Enhanced IVW solver with gradient descent and sum-to-one constraint
ivw_solve_mixture_gd <- function(y_all, g_proj, date_list, lib_factor, veg_kept,
                               init_weights = NULL, learning_rate = 0.01, max_iter = 1000) {
  # Build B_all (rows = dates, cols = veg_kept) projected through factors
  if(is.null(date_list) || length(date_list) == 0) stop("date_list cannot be null or empty")

  # Determine rank from the projection object (PCA or LDA)
  rnk <- if (!is.null(g_proj$rank)) g_proj$rank else ncol(g_proj$fit$scaling)

  F_matrix <- matrix(0, nrow = length(date_list), ncol = length(veg_kept))
  for(i in seq_along(date_list)) {
    doy <- date_list[[i]]$doy
    for(j in seq_along(veg_kept)) {
      vkey <- veg_kept[j]
      proto <- lib_factor[[vkey]]
      if(!is.null(proto) && !is.null(proto$T) && doy > 0 && doy <= nrow(proto$T) && ncol(proto$T) >= rnk) {
        F_matrix[i, j] <- sum(proto$T[doy, seq_len(rnk)] * date_list[[i]]$z[seq_len(rnk)])
      }
    }
  }

  # Use inverse-variance weights if available from proto$T_var, else uniform
  w_all <- rep(1, nrow(F_matrix))
  if(length(veg_kept) > 0) {
    date_var <- rep(NA_real_, nrow(F_matrix))
    for(i in seq_len(nrow(F_matrix))) {
      vars <- numeric(0)
      for(j in seq_along(veg_kept)) {
        vkey <- veg_kept[j]
        proto <- lib_factor[[vkey]]
        if(!is.null(proto) && !is.null(proto$T_var) && length(proto$T_var) >= date_list[[i]]$doy) {
          d <- date_list[[i]]$doy
          if(d >=1 && d <= length(proto$T_var) && is.finite(proto$T_var[d])) {
            # Use MAD-based variance for robustness
            mad_val <- mad(proto$T_var[max(1, d-7):min(length(proto$T_var), d+7)],
                          na.rm = TRUE, constant = 1.4826)
            var_val <- if(is.finite(mad_val)) mad_val^2 else proto$T_var[d]
            # Use actual variance without sample-size correction to avoid bias toward common classes
            # The sample-size correction was creating bias by giving lower variance to classes with more samples
            vars <- c(vars, var_val)
          }
        }
      }
      if(length(vars) > 0) date_var[i] <- mean(vars, na.rm = TRUE)
    }
    date_var[!is.finite(date_var) | date_var <= 0] <- NA_real_
    if(all(is.na(date_var))) {
      w_all <- rep(1, nrow(F_matrix))
    } else {
      # Use KDE densities directly as weights
      w_all <- date_var
      w_all[!is.finite(w_all) | w_all <= 0] <- 1
    }
  }

  # DIAGNOSTIC: Check F_matrix for potential issues and apply critical normalization
  if(nrow(F_matrix) > 0 && ncol(F_matrix) > 0) {
    f_colmeans <- colMeans(F_matrix, na.rm = TRUE)
    f_colsds <- apply(F_matrix, 2, sd, na.rm = TRUE)
    names(f_colmeans) <- veg_kept
    names(f_colsds) <- veg_kept
    cat(sprintf("[F_MATRIX_DIAG] F_matrix %dx%d, veg_means: %s\n", nrow(F_matrix), ncol(F_matrix),
                paste(sprintf("%s=%.4f", names(f_colmeans), f_colmeans), collapse=", ")))
    
    # Check if one vegetation type dominates - this indicates corrupted spectral library
    if(length(f_colmeans) >= 2) {
      max_idx <- which.max(abs(f_colmeans))
      other_mean <- mean(abs(f_colmeans[-max_idx]))
      dominance_ratio <- abs(f_colmeans[max_idx]) / max(other_mean, 1e-6)
      if(dominance_ratio > 5) {
        # Collect normalization data instead of printing individual warnings
        spectral_normalizations <<- spectral_normalizations + 1
        normalized_vegetations <<- c(normalized_vegetations, names(f_colmeans)[max_idx])
        
        # CRITICAL FIX: Normalize spectral signatures to prevent bias
        # Scale each vegetation column to have similar magnitude (unit variance)
        for(j in 1:ncol(F_matrix)) {
          col_vals <- F_matrix[, j]
          col_mean <- mean(col_vals, na.rm = TRUE)
          col_sd <- sd(col_vals, na.rm = TRUE)
          if(is.finite(col_sd) && col_sd > 0) {
            F_matrix[, j] <- (col_vals - col_mean) / col_sd
          } else {
            F_matrix[, j] <- col_vals - col_mean  # just center if no variance
          }
        }
        
        # Verify normalization worked
        norm_means <- colMeans(F_matrix, na.rm = TRUE)
        norm_sds <- apply(F_matrix, 2, sd, na.rm = TRUE)
        names(norm_means) <- veg_kept
        names(norm_sds) <- veg_kept
        cat(sprintf("    FIXED: Normalized means: %s\n", 
                   paste(sprintf("%s=%.4f", names(norm_means), norm_means), collapse=", ")))
        cat(sprintf("    FIXED: Normalized sds: %s\n", 
                   paste(sprintf("%s=%.4f", names(norm_sds), norm_sds), collapse=", ")))
      }
    }
  }

  sw <- sqrt(w_all)
  # The response is now implicitly 1, as we are solving for weights that sum to 1
  y_w <- rep(1, length(sw)) * sw
  B_w <- sweep(F_matrix, 1, sw, "*")

  col_scale <- sqrt(colSums(B_w^2)); col_scale[!is.finite(col_scale) | col_scale <= 0] <- 1
  B_ws <- sweep(B_w, 2, col_scale, "/")
  base_mask <- rep(TRUE, ncol(B_ws))
  # Compute condition number from weighted, scaled Gram matrix to choose ridge
  cond_val <- tryCatch({
    k <- stats::kappa(crossprod(B_ws))
    if(!is.finite(k) || k <= 0) NA_real_ else as.numeric(k)
  }, error = function(e) {
    svd_vals <- tryCatch(svd(crossprod(B_ws), nu = 0, nv = 0)$d, error = function(e) numeric(0))
    if(length(svd_vals) >= 2 && all(is.finite(svd_vals)) && svd_vals[1] > 0 && tail(svd_vals,1) > 0) svd_vals[1] / tail(svd_vals,1) else NA_real_
  })
  ridge_val <- compute_dynamic_ridge(cond_val, 1e10, 0.01, 1e2)
  if(!is.finite(ridge_val) || ridge_val <= 0) ridge_val <- compute_dynamic_ridge(1e6, 1e10, 0.01, 1e2)

  # Use gradient descent solver
  sol <- tryCatch({
    solve_mixture_gd(B_ws, y_w, base_mask = base_mask, lambda2 = ridge_val,
                    learning_rate = learning_rate, max_iter = max_iter)
  }, error = function(e) stop("solve_mixture_gd failed: ", e$message))

  if(is.null(sol) || is.null(sol$X)) stop("solve_mixture_gd returned invalid solution")
  coef_gd <- as.numeric(sol$X) / col_scale
  names(coef_gd) <- veg_kept
  # Gradient descent already ensures simplex constraint

  # Try Nelder-Mead simplex on original F_matrix with weights w_all
  coef_init <- if(!is.null(init_weights)) init_weights else coef_gd
  target_obs <- sapply(date_list, function(d) sqrt(sum(d$z^2)))
  coef_simplex <- optimize_weights_simplex(F_matrix, w_all, target = target_obs, init = coef_init)
  # optimize_weights_simplex already ensures simplex constraint

  # Compute weighted SSE using vector reconstruction objective (block matrix design)
  sse_calc <- function(w) {
    # Build block design matrix for vector reconstruction
    n_dates <- nrow(F_matrix)
    n_veg <- ncol(F_matrix)
    n_comp <- rnk

    # X_block: (n_dates * n_comp) x n_veg
    X_block <- matrix(0, nrow = n_dates * n_comp, ncol = n_veg)
    # y_block: stacked observed factor vectors
    y_block <- numeric(n_dates * n_comp)
    # w_block: expanded weights for each component
    w_block <- numeric(n_dates * n_comp)

    for(i in seq_len(n_dates)) {
      for(j in seq_len(n_veg)) {
        vkey <- veg_kept[j]
        proto <- lib_factor[[vkey]]
        if(!is.null(proto) && !is.null(proto$T) && date_list[[i]]$doy > 0 &&
           date_list[[i]]$doy <= nrow(proto$T) && ncol(proto$T) >= n_comp) {
          # Fill X_block: each block of n_comp rows for date i, veg j
          row_start <- (i-1) * n_comp + 1
          row_end <- i * n_comp
          X_block[row_start:row_end, j] <- proto$T[date_list[[i]]$doy, seq_len(n_comp)]
        }
      }
      # Fill y_block: observed factor vector for date i
      y_block[((i-1)*n_comp + 1):(i*n_comp)] <- date_list[[i]]$z[seq_len(n_comp)]
      # Fill w_block: repeat weight for each component
      w_block[((i-1)*n_comp + 1):(i*n_comp)] <- w_all[i]
    }

    # Vector reconstruction objective: minimize sum(w_block * (X_block %*% w - y_block)^2)
    pred_block <- X_block %*% w
    res_block <- pred_block - y_block
    s <- sum(w_block * (res_block^2), na.rm = TRUE)
    if(!is.finite(s)) s <- .Machine$double.xmax
    s
  }
  sse_gd <- tryCatch(sse_calc(coef_gd), error = function(e) stop(sprintf("sse_calc failed for gradient-descent solution: %s", e$message)))
  sse_simplex <- tryCatch(sse_calc(coef_simplex), error = function(e) stop(sprintf("sse_calc failed for simplex solution: %s", e$message)))

  if(is.finite(sse_simplex) && sse_simplex < sse_gd) {
    coef_best <- coef_simplex
  } else {
    coef_best <- coef_gd
  }

  # Optimization methods already ensure simplex constraint - no additional projection needed
  coef_final <- coef_best
  
  # Ensure returned vector follows original veg_kept ordering and is numeric
  coef_final <- as.numeric(coef_final)
  names(coef_final) <- veg_kept
  
  # Print spectral normalization summary
  if(spectral_normalizations > 0) {
    cat(sprintf("SPECTRAL_NORMALIZATION_SUMMARY: Applied spectral normalization to %d vegetation types (%s)\n", 
                spectral_normalizations, paste(unique(normalized_vegetations), collapse=", ")))
  }
  
  return(list(X = coef_final, fmat = NULL))
}

  # date_list: list of date entries used to build B_all; used to map obs -> doy
  # lib_factor_pca: list of per-veg prototypes with T_var (variances) for IVW weighting
  # veg_kept: character vector of vegetation names corresponding to columns

# Progress logging helpers (defaults centralized at top; guard here to avoid override)
if(!exists("PROGRESS_EVERY_TASK")) PROGRESS_EVERY_TASK <- 25
# Also write logs to a shared file so messages from parallel workers are captured on Windows
if(!exists("PROGRESS_LOG_TO_FILE")) PROGRESS_LOG_TO_FILE <- TRUE
LOG_FILE <- tryCatch({
  if(exists("OUT_DIR") && is.character(OUT_DIR) && nchar(OUT_DIR) > 0) {
    file.path(OUT_DIR, "fit_veg_bootstrap.log")
  } else {
    file.path(getwd(), "fit_veg_bootstrap.log")
  }
}, error = function(e) file.path(getwd(), "fit_veg_bootstrap.log"))
log_msg <- function(...) {
  ts <- format(Sys.time(), "%H:%M:%S")
  msg <- sprintf("[%s] %s\n", ts, sprintf(...))
  cat(msg)
  if(isTRUE(PROGRESS_LOG_TO_FILE)) try(cat(msg, file = LOG_FILE, append = TRUE), silent = TRUE)
}


# Parallelization controls (Windows-friendly via multisession)
if(!exists("PARALLEL_ENABLE")) PARALLEL_ENABLE <- TRUE
if(!exists("PROGRESS_BAR")) PROGRESS_BAR <- TRUE
if(!exists("PARALLEL_WORKERS")) PARALLEL_WORKERS <- tryCatch({
  if(requireNamespace("parallel", quietly = TRUE)) max(1L, parallel::detectCores(logical = TRUE) - 1L) else 1L
}, error = function(...) 1L)

# Helper: parallel map - fail hard if requirements not met
.run_map <- function(X, FUN) {
  # Bind local copy so futures capture the actual function object
  f_FUN <- FUN
  
  if(!PARALLEL_ENABLE) {
    # Serial path - progressr required
    if(!requireNamespace("progressr", quietly = TRUE)) {
      stop("progressr package required for serial processing - no fallback allowed")
    }
    progressr::with_progress({
      p <- progressr::progressor(steps = length(X))
      lapply(X, function(x) { res <- f_FUN(x); p(); res })
    })
  } else {
    # Parallel path - future packages required
    if(!requireNamespace("future.apply", quietly = TRUE) || !requireNamespace("future", quietly = TRUE)) {
      stop("future.apply and future packages required for parallel processing - no fallback allowed")
    }
    if(!requireNamespace("progressr", quietly = TRUE)) {
      stop("progressr package required for progress tracking - no fallback allowed")
    }
    
    old_plan <- future::plan()
    options(future.globals.maxSize = max(getOption("future.globals.maxSize", 2e9), 4e9))
    future::plan(future::multisession, workers = PARALLEL_WORKERS)
    on.exit(future::plan(old_plan))

    progressr::with_progress({
      p <- progressr::progressor(steps = length(X))
      future.apply::future_lapply(X, function(x) { res <- f_FUN(x); p(); res }, future.seed = TRUE)
    })
  }
}

# Ensure date/year columns exist for downstream filtering and grouping
if(!"date" %in% names(df) && "Date" %in% names(df)) df$date <- df$Date
if("date" %in% names(df)) {
  df$date <- as.Date(df$date)
  if(!"year" %in% names(df)) df$year <- as.integer(lubridate::year(df$date))
geojson_path <- file.path("fotos","zuiver_zonder_foto.geojson")
if(!file.exists(geojson_path)) stop(paste0("GeoJSON points not found at ", geojson_path))
gpts_raw <- sf::st_read(geojson_path, quiet = TRUE)
# Identify vegetation attribute (case-insensitive among vegetation, veg, class)
# Use original column name matching so we can index the data.frame safely.
matched_cols <- names(gpts_raw)[tolower(names(gpts_raw)) %in% c("vegetation","veg","class")]
if(length(matched_cols) > 0) {
  veg_col_orig <- matched_cols[1]
  gpts_raw$.__veg__ <- as.character(gpts_raw[[veg_col_orig]])
} else {
  gpts_raw$.__veg__ <- NA_character_
}
coords <- sf::st_coordinates(gpts_raw)
gpts_raw$.__lon__ <- coords[,1]; gpts_raw$.__lat__ <- coords[,2]
gpts_raw$location_id <- make_location_id(gpts_raw$.__lon__, gpts_raw$.__lat__)
# Create a lightweight mapping table (location_id -> Veg) for joins below.
gpts_map <- sf::st_drop_geometry(gpts_raw) %>%
  dplyr::select(location_id, Veg = .__veg__) %>%
  dplyr::distinct(location_id, .keep_all = TRUE)

if(nrow(gpts_map) == 0) stop(paste0("GeoJSON mapping produced no valid points (gpts_map empty). Check ", geojson_path))

# Join df with gpts_map using the standardized location_id
if ("location_id" %in% names(df) && "location_id" %in% names(gpts_map)) {
  df <- dplyr::left_join(df, gpts_map, by = "location_id")
}

## NOTE: loc_years will be constructed after all filtering and canonicalization
## to avoid producing task entries for location-years that are removed by
## flood/index/veg filters. A placeholder is defined here so references
## earlier in the file do not error; the real loc_years is built later.
loc_years <- data.frame(location_id = character(0), year = integer(0), stringsAsFactors = FALSE)

# Always (re)derive a standardized location_id in df from any lon/lat style columns to maximize matches.
lon_candidates <- names(df)[grepl("(^|_)(lon|longitude|x)(_|$)", names(df), ignore.case=TRUE)]
lat_candidates <- names(df)[grepl("(^|_)(lat|latitude|y)(_|$)", names(df), ignore.case=TRUE)]
if(length(lon_candidates)>0 && length(lat_candidates)>0) {
  # Use the first available pair
  lon_col <- lon_candidates[1]
  lat_col <- lat_candidates[1]
  df$location_id <- make_location_id(df[[lon_col]], df[[lat_col]])
}

# Ensure a `year` column exists (derive from `date` if missing) to prevent filter failures
if(!"year" %in% names(df)) {
  if("date" %in% names(df)) {
    # Use lubridate::year safely without attaching the package
    if(!requireNamespace("lubridate", quietly = TRUE)) stop("The package 'lubridate' is required but not installed")
    df$year <- as.integer(lubridate::year(as.Date(df$date)))
  } else {
    # attempt to find a date-like column
    possible_date <- grep("date|time|obs", names(df), ignore.case = TRUE, value = TRUE)
    if(length(possible_date) > 0) {
      if(!requireNamespace("lubridate", quietly = TRUE)) stop("The package 'lubridate' is required but not installed")
      df$year <- as.integer(lubridate::year(as.Date(df[[possible_date[1]]])))
    }
  }
}

# Simplified spatial join: associate each point in df with the nearest feature in gpts_raw.
# This replaces the multi-step logic of direct merge and appending missing locations.
if(!"Veg" %in% names(df)) df$Veg <- NA_character_

# PRIORITIZING GeoJSON OVER CSV: Vegetation ground truth data is now loaded from GeoJSON file instead of CSV
# The GeoJSON reading code is located earlier in the script and provides more accurate vegetation classification
# CSV reading has been disabled to ensure consistent use of GeoJSON data source

# Ensure df has coordinates for joining.
lon_cols <- names(df)[grepl("(^|_)(lon|longitude|x)(_|$)", names(df), ignore.case=TRUE)]
lat_cols <- names(df)[grepl("(^|_)(lat|latitude|y)(_|$)", names(df), ignore.case=TRUE)]
if (length(lon_cols) > 0 && length(lat_cols) > 0) {
  # Join again to fill in any missing Veg values using the explicit 2-col map
  if ("location_id" %in% names(df) && "location_id" %in% names(gpts_map) && nrow(gpts_map) > 0) {
    df <- dplyr::left_join(df, gpts_map, by = "location_id", suffix = c("", ".y"))
    # If existing Veg is missing, replace with mapped Veg
    if ("Veg.y" %in% names(df)) {
      # Only attempt to merge if there's something to merge
      if(any(!is.na(df$Veg.y))) {
          df$Veg <- ifelse(is.na(df$Veg) | df$Veg == "", df$Veg.y, df$Veg)
      }
      df$Veg.y <- NULL
    }
  }
}

df$Veg <- tolower(df$Veg)
if(!"date" %in% names(df)) stop("Input CSV must contain a 'date' column")
df$date <- as.Date(df$date)
if(!"location_id" %in% names(df)) stop("Input CSV must contain a 'location_id' column")
if(!"Veg" %in% names(df)) df$Veg <- NA_character_

# Standardize vegetation names to lower-case early
if("Veg" %in% names(df)) df$Veg <- tolower(as.character(df$Veg))

# Ensure required packages are available: install if missing (best-effort)
required_pkgs <- c("future", "future.apply", "progressr")
for(pkg in required_pkgs) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    try({
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }, silent = TRUE)
  }
}

# Ensure openxlsx is available for Excel output
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx", repos = "https://cloud.r-project.org")
}

# Diagnostic: show unique vegetation values to help debug collapsed classes



# Helper: compute DOY
df$doy <- lubridate::yday(df$date)
df$doy[df$doy < 1 | df$doy > 366] <- NA_integer_

# Global veg frequency excluding excluded classes
veg_counts <- sort(table(na.omit(df$Veg)), decreasing = TRUE)

# ---- Define indices to use: restrict to OPTIMAL_INDICES by default ----
# Strict: only OPTIMAL_INDICES are accepted - no fallback to raw bands
meta_cols <- intersect(c(
  "date","location_id","Veg","coverage","lat","lon","latitude","longitude",
  "target_lon","target_lat","imagery_lat","imagery_lon","doy","year"
), names(df))
numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
candidate_indices <- intersect(OPTIMAL_INDICES, numeric_cols)
found_opt <- intersect(OPTIMAL_INDICES, numeric_cols)
missing_opt <- setdiff(OPTIMAL_INDICES, numeric_cols)
if(length(found_opt) > 0) {
  cat(sprintf("Found %d OPTIMAL_INDICES in input: %s\n", length(found_opt), paste(found_opt, collapse = ", ")))
}
if(length(missing_opt) > 0) {
  cat(sprintf("Missing OPTIMAL_INDICES in input: %s\n", paste(missing_opt, collapse = ", ")))
}
candidate_indices <- found_opt
if(length(candidate_indices) == 0) {
  stop("No OPTIMAL_INDICES present in input CSV. This fitter requires the canonical OPTIMAL_INDICES to be present in the input CSV. Please run the predictor/extractor to generate them.")
}
# ---- Global amplitude-based index selection ----
# Calculate amplitude for each candidate index across the entire dataset
# and filter to keep only indices with amplitude > 0.05
cat("Calculating global index amplitudes...\n")
index_amplitudes <- vapply(candidate_indices, function(idx) {
  vals <- df[[idx]][is.finite(df[[idx]])]
  if(length(vals) < 10) return(0)  # Need minimum observations
  # Use range (max - min) as amplitude measure
  amp <- diff(range(vals, na.rm = TRUE))
  if(!is.finite(amp)) 0 else amp
}, numeric(1))

# Filter indices with amplitude > 0.05
amplitude_threshold <- 0.05
high_amplitude_indices <- names(index_amplitudes)[index_amplitudes > amplitude_threshold]
low_amplitude_indices <- names(index_amplitudes)[index_amplitudes <= amplitude_threshold]

if(length(low_amplitude_indices) > 0) {
  cat(sprintf("Filtered out %d low-amplitude indices (amplitude <= %.3f): %s\n",
              length(low_amplitude_indices), amplitude_threshold,
              paste(low_amplitude_indices, collapse = ", ")))
}

# Use high-amplitude indices as available indices
avail <- high_amplitude_indices

# ---- Correlation-based filtering: remove highly correlated indices ----
# This reduces redundancy and improves computational efficiency
# while maintaining the most informative spectral indices
if(length(avail) > 1) {
  cat("Applying correlation-based filtering...\n")

  # Calculate correlation matrix for available indices
  corr_data <- df[, avail, drop = FALSE]
  corr_matrix <- cor(corr_data, use = "pairwise.complete.obs")
  corr_matrix[is.na(corr_matrix)] <- 0

  # Remove indices with correlation > 0.90
  high_corr_threshold <- 0.90
  to_remove <- c()

  for(i in 1:(ncol(corr_matrix)-1)) {
    for(j in (i+1):ncol(corr_matrix)) {
      if(abs(corr_matrix[i,j]) > high_corr_threshold) {
        # Keep the index with higher amplitude (more informative)
        idx1 <- colnames(corr_matrix)[i]
        idx2 <- colnames(corr_matrix)[j]

        amp1 <- index_amplitudes[idx1]
        amp2 <- index_amplitudes[idx2]

        # If amplitudes are very close (within 1%), keep the first one
        if(abs(amp1 - amp2) / max(amp1, amp2) < 0.01) {
          # Keep idx1 (arbitrary choice when amplitudes are equal)
          to_remove <- c(to_remove, idx2)
          cat(sprintf("Keeping %s (equal amplitude to %s, correlation=%.3f)\n",
                      idx1, idx2, corr_matrix[i,j]))
        } else {
          # Remove the one with lower amplitude
          if(amp1 < amp2) {
            to_remove <- c(to_remove, idx1)
          } else {
            to_remove <- c(to_remove, idx2)
          }
        }
      }
    }
  }

  to_remove <- unique(to_remove)
  if(length(to_remove) > 0) {
    avail <- setdiff(avail, to_remove)
    cat(sprintf("Removed %d highly correlated indices: %s\n",
                length(to_remove), paste(to_remove, collapse = ", ")))
  }

  cat(sprintf("Final selection: %d indices after correlation filtering\n", length(avail)))
}

if(length(avail) == 0) {
  stop("No indices found with amplitude > 0.05 after correlation filtering; cannot proceed with fitting.")
}

if(length(avail) < USE_INDICES_MIN) {
  stop(sprintf("Only %d indices remain after correlation filtering, but minimum required is %d.",
               length(avail), USE_INDICES_MIN))
}

cat(sprintf("Selected %d indices with amplitude > %.3f: %s\n",
            length(avail), amplitude_threshold, paste(avail, collapse = ", ")))

# Apply to each available index (after amplitude filtering) and add as new columns
cat("Computing moving variances for indices (computed on full dataset and propagated to train/test)...\n")
# Strict data sufficiency check: require at least two seasonal cycles (default 2*365 days)
min_required_days <- 2 * 365
if ("location_id" %in% names(df_full)) {
  offenders <- list()
  for (loc in unique(df_full$location_id)) {
    sub <- df_full[df_full$location_id == loc, , drop = FALSE]
    for (idx_name in avail) {
      if (!idx_name %in% names(sub)) next
      nobs <- sum(is.finite(sub[[idx_name]]))
      if (nobs < min_required_days) {
        offenders[[length(offenders) + 1]] <- list(location = loc, index = idx_name, n = nobs)
      }
    }
  }
  if (length(offenders) > 0) {
    msg_lines <- vapply(offenders, function(x) sprintf("loc=%s index=%s n=%d", x$location, x$index, x$n), character(1))
    stop(sprintf("Data sufficiency check failed: some location/index pairs have fewer than %d observations:\n%s", min_required_days, paste(msg_lines, collapse = "\n")))
  }
} else {
  if (nrow(df_full) < min_required_days) stop(sprintf("Data sufficiency check failed: overall data contains fewer than %d rows", min_required_days))
}

# Compute moving variance on the full dataset so both training and test rows receive the _var14 columns
if(length(avail) > 0) {
  pb <- utils::txtProgressBar(min = 0, max = length(avail), style = 3)
  i <- 0
  for(idx_name in avail) {
    var_col <- paste0(idx_name, "_var14")
    # compute on df_full (original full dataframe)
    v_full <- calc_moving_var(df_full, idx_name, window = 14)
    df_full[[var_col]] <- v_full

    # propagate to current df (training subset) by matching location_id+date
    if(exists("df") && nrow(df) > 0) {
      key_full <- paste0(as.character(df_full$location_id), "__", as.character(df_full$date))
      key_df   <- paste0(as.character(df$location_id), "__", as.character(df$date))
      mpos <- match(key_df, key_full)
      df[[var_col]] <- ifelse(is.na(mpos), NA_real_, df_full[[var_col]][mpos])
    }

    # propagate to df_test as well (if present)
    if(exists("df_test") && nrow(df_test) > 0) {
      key_test <- paste0(as.character(df_test$location_id), "__", as.character(df_test$date))
      key_full <- paste0(as.character(df_full$location_id), "__", as.character(df_full$date))
      mpos2 <- match(key_test, key_full)
      df_test[[var_col]] <- ifelse(is.na(mpos2), NA_real_, df_full[[var_col]][mpos2])
    }

    i <- i + 1; try(utils::setTxtProgressBar(pb, i), silent = TRUE)
  }
  try(close(pb), silent = TRUE)
} else {
  cat("No indices available to compute moving variances.\n")
}

# Optimize seasonal variance selection
seasonal_var_cols <- c()
if(length(avail) > 0) {
  # Pre-compute peaks and troughs for all available indices at once
  peak_trough_data <- lapply(avail, function(idx_name) {
    var_col <- paste0(idx_name, "_var14")
    v <- df[[var_col]]
    if(all(is.na(v))) return(NULL)
    peaks <- which(diff(sign(diff(v))) == -2)
    troughs <- which(diff(sign(diff(v))) == 2)
    if(length(peaks) >= 1 && length(troughs) >= 1) var_col else NULL
  })
  seasonal_var_cols <- unlist(peak_trough_data[!sapply(peak_trough_data, is.null)])
}

df <- df[, c(names(df)[!grepl("_var14$", names(df))], seasonal_var_cols)]

timing_info$moving_var_done <- Sys.time()
cat(sprintf("Moving variance calculation completed in %.1f seconds\n", 
           as.numeric(difftime(timing_info$moving_var_done, timing_info$start_time, units="secs"))))

# -------------------------
# Post-variance cleaning
# - remove "random-looking" variance columns (low autocorrelation / near-constant)
# - remove highly correlated variance columns (redundant)
# -------------------------
{
  # Identify moving-variance columns (either _var14 or _mv from earlier steps)
  var_cols <- names(df)[grepl("(_var14$|_mv$)", names(df))]
  if(length(var_cols) > 0) {
    cat(sprintf("Post-variance cleaning: %d variance columns detected\n", length(var_cols)))

    # Compute heuristics: lag-1 autocorrelation and coefficient of variation for each var col
    ac1 <- numeric(length(var_cols)); names(ac1) <- var_cols
    cv  <- numeric(length(var_cols)); names(cv) <- var_cols
    finite_frac <- numeric(length(var_cols)); names(finite_frac) <- var_cols

    for(i in seq_along(var_cols)) {
      v <- df[[var_cols[i]]]
      v_f <- v[is.finite(v)]
      finite_frac[i] <- if(length(v) == 0) 0 else sum(is.finite(v)) / length(v)
      if(length(v_f) > 3) {
        # lag-1 autocorr
        ac1[i] <- tryCatch({ cor(v_f[-1], v_f[-length(v_f)], use = "complete.obs") }, error = function(e) stop(sprintf("autocorrelation computation failed for %s: %s", var_cols[i], e$message)))
        # coefficient of variation (sd/mean) robustified by median
        mm <- median(v_f, na.rm = TRUE)
        msd <- stats::sd(v_f, na.rm = TRUE)
        cv[i] <- if(is.finite(mm) && mm != 0) msd / abs(mm) else Inf
      } else {
        ac1[i] <- NA_real_;
        cv[i] <- NA_real_;
      }
    }

    # Heuristic thresholds: require at least 30% finite coverage, ac1 > 0.1 and cv > 0.05
    keep_mask <- (finite_frac >= 0.30) & (is.finite(ac1) & ac1 > 0.10 | is.finite(cv) & cv > 0.05)
    remove_rand <- var_cols[!keep_mask]
    if(length(remove_rand) > 0) {
      cat(sprintf("Removing %d random-looking variance cols: %s\n", length(remove_rand), paste(remove_rand, collapse=", ")))
      df[remove_rand] <- NULL
      var_cols <- setdiff(var_cols, remove_rand)
    }

    # Recompute if any remain, then remove highly correlated variance columns (>0.90)
    if(length(var_cols) > 1) {
      vmat <- as.data.frame(lapply(var_cols, function(nm) df[[nm]]))
      names(vmat) <- var_cols
      # Use pairwise complete obs correlation
      cm <- stats::cor(vmat, use = "pairwise.complete.obs")
      cm[is.na(cm)] <- 0
      high_corr_thresh <- 0.95
      to_remove <- c()
      for(i in seq_len(ncol(cm)-1)) {
        for(j in (i+1):ncol(cm)) {
          if(abs(cm[i,j]) > high_corr_thresh) {
            # prefer to keep the column with higher finite fraction
            c1 <- var_cols[i]; c2 <- var_cols[j]
            f1 <- sum(is.finite(df[[c1]])); f2 <- sum(is.finite(df[[c2]]))
            if(f1 >= f2) to_remove <- c(to_remove, c2) else to_remove <- c(to_remove, c1)
          }
        }
      }
      to_remove <- unique(to_remove)
      if(length(to_remove) > 0) {
        cat(sprintf("Removing %d highly correlated variance cols: %s\n", length(to_remove), paste(to_remove, collapse=", ")))
        df[to_remove] <- NULL
      }
    }
  }
}


cat(sprintf("Post baseline normalization rows=%d\n", nrow(df)))
cat("Data preprocessing complete.\n")
adj_cols <- intersect(avail, names(df))

# --- Simple feature pruning: drop indices (including *_var30) with very low variability or sparse coverage ---
if(!exists("MIN_INDEX_SD")) MIN_INDEX_SD <- 0.05  # more lenient
if(length(avail) > 0) {
  idx_sd <- vapply(avail, function(nm) {
    x <- df[[nm]]; x <- x[is.finite(x)]
    if(length(x) < 5) return(0)
    stats::sd(x, na.rm=TRUE)
  }, numeric(1))
  keep <- idx_sd >= MIN_INDEX_SD
  dropped <- avail[!keep]
  if(length(dropped) > 0) cat("Pruned low-variance indices:", paste(dropped, collapse=", "), "\n")
  avail <- avail[keep]
}

# Restrict dataset to allowed vegetation classes (if Veg present)
# Normalize variant names: map any Veg that contains an allowed substring to that canonical name
if("Veg" %in% names(df) && length(ALLOWED_VEG) > 0) {
  try({
    for(av in ALLOWED_VEG) {
      sel <- grepl(av, df$Veg, ignore.case = TRUE) & !is.na(df$Veg)
      if(any(sel)) {
        df$Veg[sel] <- av
        cat(sprintf("Normalized %d rows to canonical Veg '%s'\n", sum(sel, na.rm = TRUE), av))
      }
    }
  }, silent = TRUE)
}
if("Veg" %in% names(df) && length(ALLOWED_VEG) > 0) {
  keep_rows <- tolower(df$Veg) %in% ALLOWED_VEG
  n_before <- nrow(df)
  df <- df[keep_rows | is.na(df$Veg), , drop = FALSE]
  cat(sprintf("Filtered to allowed classes (%s): kept %d/%d rows (NAs retained)\n", paste(ALLOWED_VEG, collapse=","), nrow(df), n_before))
}

# Per-vegetation quick summary to diagnose why some vegs may be dropped later
try({
  cat("Per-veg quick summary (rows, unique_doys, sample index finite counts):\n")
  for(av in ALLOWED_VEG) {
    sel <- tolower(df$Veg) == av
    rows <- sum(sel, na.rm = TRUE)
    unique_doys <- length(unique(df$doy[sel & is.finite(df$doy)]))
    idx_counts <- integer(0)
    if(length(avail) > 0) idx_counts <- sapply(avail, function(ix) sum(is.finite(df[[ix]][sel])), USE.NAMES = TRUE)
    samp_idx <- if(length(idx_counts) > 0) paste(names(idx_counts)[which(idx_counts>0)[seq_len(min(5, sum(idx_counts>0)))]], collapse = ",") else "none"
    cat(sprintf("SUMMARY: veg=%s rows=%d unique_doys=%d idx_present_sample=%s\n", av, rows, unique_doys, samp_idx))
  }
}, silent = TRUE)

# Quick check: do we have any vegetation labels?
matched_veg_n <- sum(!is.na(df$Veg))
cat("Non-NA Veg rows:", matched_veg_n, "of", nrow(df), "\n")
if(matched_veg_n == 0) {
  stop("No vegetation classes found after join; cannot build library. Ensure the extractor wrote Veg, or fix the GeoJSON (rounded 4 decimals).")
}

# Recompute standardized location_id from lon/lat if present (ensure same format used by GeoJSON mapping)
lon_candidates <- names(df)[grepl("(^|_)(lon|longitude|x)(_|$)", names(df), ignore.case=TRUE)]
lat_candidates <- names(df)[grepl("(^|_)(lat|latitude|y)(_|$)", names(df), ignore.case=TRUE)]
if(length(lon_candidates)>0 && length(lat_candidates)>0) {
  lon_col <- lon_candidates[1]; lat_col <- lat_candidates[1]
  df$location_id <- make_location_id(df[[lon_col]], df[[lat_col]])
}

# Now construct loc_years from the filtered and canonicalized df so tasks only cover existing rows
loc_years <- df %>% dplyr::filter(!is.na(location_id) & location_id != "" & !is.na(year) & year > 0) %>% dplyr::distinct(location_id, year)
cat(sprintf("Constructed loc_years with %d rows from filtered df\n", nrow(loc_years)))
if(nrow(loc_years) == 0) stop("No location-year pairs found after filtering and canonicalization; aborting. Check INPUT_CSV and Veg/geolocation joins.")


## Always construct the vegetation library `lib` from the TRAINING dataset `df`.
# Build for each vegetation class a per-index 365-length median "mu" vector.
cat("Constructing lib from TRAINING dataset (optimized)...\n")
lib <- list()
# Use training data for library construction
lib_df <- df
vegs <- unique(na.omit(lib_df$Veg))
vegs <- vegs[vegs != ""]



# ---- Sample Size Balancing (Always Enabled) ----
# Balance training data to prevent bias toward vegetation types with more samples
cat("Applying sample size balancing to training data...\n")

# Pre-compute DOY for sample balancing
lib_df$doy_for_lib <- lubridate::yday(lib_df$date)

# Calculate sample sizes per vegetation type
veg_sample_sizes <- sapply(vegs, function(v) {
  dveg <- lib_df[lib_df$Veg == v & is.finite(lib_df$doy_for_lib), , drop = FALSE]
  nrow(dveg)
})

# Calculate min and max sample sizes
min_samples <- min(veg_sample_sizes)
max_samples <- max(veg_sample_sizes)
median_samples <- median(veg_sample_sizes)

# Use a balanced target size that doesn't overly favor the smallest class
# Target size should be sufficient for robust estimation but balanced across classes
target_size <- pmax(min_samples, floor(median_samples * 0.8))  # 80% of median, but at least min
target_size <- pmin(target_size, max_samples)  # Don't exceed what any class has

cat(sprintf("Sample sizes - Min: %.0f, Max: %.0f, Median: %.1f\n", min_samples, max_samples, median_samples))
cat(sprintf("Balancing to target size: %.0f (80%% of median, clamped to [min, max])\n", target_size))

# Show sample size distribution
cat("Sample sizes by vegetation:\n")
for(v in names(veg_sample_sizes)) {
  cat(sprintf("  %s: %.0f samples\n", v, veg_sample_sizes[v]))
}

# Balance each vegetation type
## SMOTE-like augmentation for minority classes
augment_minority_class <- function(df_class, target_n, seed = NULL, alpha_range = c(0.3, 0.7), jitter_frac = 1e-6) {
  # df_class: data.frame of existing samples for a single class
  # target_n: desired number of rows after augmentation
  if(!is.data.frame(df_class)) stop("augment_minority_class requires a data.frame")
  if(nrow(df_class) >= target_n) return(df_class)
  if(nrow(df_class) <= 1) {
    # Not enough samples to interpolate: replicate with small jitter
    reps <- target_n - nrow(df_class)
    base <- df_class[rep(1, reps), , drop = FALSE]
    # jitter numeric columns
    num_cols <- sapply(base, is.numeric)
    if(any(num_cols)) {
      rng <- apply(df_class[, num_cols, drop = FALSE], 2, function(x) if(all(!is.finite(x))) 0 else diff(range(x, na.rm = TRUE)))
      rng[rng == 0] <- 1.0
      base[, num_cols] <- base[, num_cols, drop = FALSE] + matrix(rnorm(reps * sum(num_cols), sd = sqrt(mean(rng)) * jitter_frac), nrow = reps)
    }
    out <- rbind(df_class, base)
    rownames(out) <- NULL
    return(out)
  }

  if(!is.null(seed)) set.seed(seed)
  n_add <- target_n - nrow(df_class)
  synth_rows <- vector("list", n_add)

  # Precompute numeric columns to interpolate and keep others as most-common or sampled
  num_cols <- names(df_class)[sapply(df_class, is.numeric)]
  other_cols <- setdiff(names(df_class), num_cols)

  for(i in seq_len(n_add)) {
    # pick two distinct random indices
    ids <- sample.int(nrow(df_class), 2, replace = FALSE)
    a <- df_class[ids[1], , drop = FALSE]
    b <- df_class[ids[2], , drop = FALSE]
    alpha <- runif(1, min(alpha_range), max(alpha_range))
    new_row <- a
    # linear interpolation for numeric columns
    if(length(num_cols) > 0) {
      va <- as.numeric(a[1, num_cols, drop = FALSE])
      vb <- as.numeric(b[1, num_cols, drop = FALSE])
      interp <- (1 - alpha) * va + alpha * vb
      # add tiny jitter proportional to column range to avoid exact duplicates
      rng <- apply(df_class[, num_cols, drop = FALSE], 2, function(x) if(all(!is.finite(x))) 0 else diff(range(x, na.rm = TRUE)))
      jitter <- rnorm(length(interp), sd = pmax(abs(rng), 1e-8) * jitter_frac)
      interp <- interp + jitter
      new_row[1, num_cols] <- interp
    }

    # For non-numeric columns, try to pick one of the parents (a or b) randomly
    for(col in other_cols) {
      if(col %in% c("date", "doy", "doy_for_lib")) next
      val <- if(runif(1) < 0.5) a[[col]] else b[[col]]
      new_row[[col]] <- val
    }

    # Handle date/doy: interpolate dates by choosing one parent's date or sampling within year using DOY
    if("doy_for_lib" %in% names(df_class)) {
      da <- as.integer(a$doy_for_lib)
      db <- as.integer(b$doy_for_lib)
      if(is.finite(da) && is.finite(db)) {
        # circular interpolation for DOY
        diff_ab <- ((db - da + 365) %% 365)
        pick_doy <- ((da + round(alpha * diff_ab) - 1) %% 365) + 1
        new_row$doy_for_lib <- as.integer(pick_doy)
      } else {
        new_row$doy_for_lib <- if(!is.na(da) && is.finite(da)) da else if(!is.na(db) && is.finite(db)) db else NA_integer_
      }
    }

    # If 'date' column present, set to NA (we don't have a real date) or reconstruct from doy if possible
    if("date" %in% names(df_class)) {
      if(!is.null(new_row$doy_for_lib) && is.finite(new_row$doy_for_lib)) {
        # create a placeholder date within arbitrary non-leap year (1970)
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
  # Preserve factor types for non-numeric columns where applicable
  for(col in other_cols) {
    if(is.factor(df_class[[col]])) out[[col]] <- factor(out[[col]], levels = levels(df_class[[col]]))
  }
  rownames(out) <- NULL
  out
}


balanced_dfs <- list()
for(v in vegs) {
  dveg <- lib_df[lib_df$Veg == v & is.finite(lib_df$doy_for_lib), , drop = FALSE]
  if(nrow(dveg) == 0) next
  if(nrow(dveg) < target_size) {
    # Augment minority class up to target_size
    augmented <- tryCatch({
      augment_minority_class(dveg, target_size)
    }, error = function(e) {
      cat(sprintf("Augmentation failed for %s: %s - falling back to replication\n", v, e$message))
      # fallback: simple replication
      dveg[rep(seq_len(nrow(dveg)), length.out = target_size), , drop = FALSE]
    })
    balanced_dfs[[v]] <- augmented
    cat(sprintf("Augmented %s from %d to %d samples\n", v, nrow(dveg), nrow(augmented)))
  } else {
    # Keep classes that have >= target_size as-is (no downsampling)
    balanced_dfs[[v]] <- dveg
    cat(sprintf("Kept %s with %d samples (>= target)\n", v, nrow(dveg)))
  }
}


# Combine balanced data
if(length(balanced_dfs) > 0) {
  lib_df <- chunked_rbind(balanced_dfs, chunk_size = 50L)
  gc()
  cat(sprintf("Training data balanced: %d total samples across %d vegetation types\n",
              nrow(lib_df), length(balanced_dfs)))
}

if(length(vegs) > 0) {
  if(requireNamespace("progressr", quietly = TRUE)) {
    # Wrap progressor creation in with_progress() to avoid global environment issues
    progressr::with_progress({
      p_lib <- progressr::progressor(steps = length(vegs), message = "Building vegetation library")
      
      # Pre-compute DOY for all data to avoid repeated calculations
      lib_df$doy_for_lib <- lubridate::yday(lib_df$date)

      for(v in vegs) {
        dveg <- lib_df[lib_df$Veg == v & is.finite(lib_df$doy_for_lib), , drop = FALSE]
        if(nrow(dveg) == 0) {
          cat(sprintf("[DIAG] veg=%s rows=0 -- skipping\n", v))
          next
        }

        # Forced console diagnostic: rows, unique locations, per-index finite counts
        idx_list <- intersect(avail, names(dveg))
        idx_counts <- if(length(idx_list) > 0) sapply(idx_list, function(ix) sum(is.finite(dveg[[ix]])), USE.NAMES = TRUE) else integer(0)
        nlocs <- length(unique(na.omit(dveg$location_id)))
        cat(sprintf("[DIAG] veg=%s rows=%d unique_locs=%d idx_present=%d\n", v, nrow(dveg), nlocs, length(idx_counts)))
        if(length(idx_counts) > 0) cat(sprintf("[DIAG] veg=%s idx_counts=%s\n", v, paste(sprintf("%s:%d", names(idx_counts), idx_counts), collapse=",")))

        lib[[v]] <- list()
        p_lib()

        # Optimized: Use data.table or matrix operations for median calculation
        kept_idx <- intersect(avail, names(dveg))
        for(idx in kept_idx) {
          # Use tapply for faster grouped median calculation
          vals_by_doy <- tapply(dveg[[idx]], dveg$doy_for_lib, function(x) {
            x <- x[is.finite(x)]
            if(length(x) > 0) stats::median(x, na.rm = TRUE) else NA_real_
          })

          # Fill the 365-length vector
          mu <- rep(NA_real_, 365)
          doy_values <- as.integer(names(vals_by_doy))
          valid_doy <- doy_values >= 1 & doy_values <= 365
          mu[doy_values[valid_doy]] <- vals_by_doy[valid_doy]

          # If the index has no finite values for all days, skip it
          if(all(!is.finite(mu))) next
          mv <- calc_moving_var(data.frame(date = 1:365, idx = mu), "idx", window=14)
          lib[[v]][[idx]] <- list(mu = mu, mv = mv)
        }
        # Store sample size for variance normalization
        lib[[v]]$n_samples <- nrow(dveg)
        cat(sprintf("Processed prototype for %s (n=%d)\n", v, nrow(dveg)))
      }

      # Clean up temporary column
      lib_df$doy_for_lib <- NULL
    })
  } else {
    stop("progressr package required for library construction - no fallback allowed")
  }
}
if(length(lib) == 0) stop("lib could not be constructed from training data: no Veg labels or indices available.")
cat(sprintf("Constructed simple lib from training dataset: vegetations=%d\n", length(lib)))

# CRITICAL DIAGNOSTIC: Check if libraries are potentially mislabeled by comparing spectral signatures
if(length(lib) >= 2) {
  veg_names <- names(lib)
  cat("[SPECTRAL_DIAG] Checking for potential library mislabeling:\n")
  for(i in 1:(length(veg_names)-1)) {
    for(j in (i+1):length(veg_names)) {
      v1 <- veg_names[i]; v2 <- veg_names[j]
      has_t1 <- !is.null(lib[[v1]]$T) && is.matrix(lib[[v1]]$T) && nrow(lib[[v1]]$T) > 0 && ncol(lib[[v1]]$T) >= 1
      has_t2 <- !is.null(lib[[v2]]$T) && is.matrix(lib[[v2]]$T) && nrow(lib[[v2]]$T) > 0 && ncol(lib[[v2]]$T) >= 1
      if(has_t1 && has_t2) {
        # Compare first component across peak growing season (DOY 120-240)
        doys <- 120:240
        valid_doys <- doys[doys <= nrow(lib[[v1]]$T) & doys <= nrow(lib[[v2]]$T)]
        if(length(valid_doys) > 10) {
          sig1 <- lib[[v1]]$T[valid_doys, 1]
          sig2 <- lib[[v2]]$T[valid_doys, 1]
          correlation <- cor(sig1, sig2, use="complete.obs")
          mean1 <- mean(sig1, na.rm=TRUE)
          mean2 <- mean(sig2, na.rm=TRUE)
          if(is.finite(correlation)) {
            cat(sprintf("  %s vs %s: correlation=%.3f, mean1=%.4f, mean2=%.4f\n", v1, v2, correlation, mean1, mean2))
            if(correlation > 0.9) {
              cat(sprintf("    WARNING: Very high correlation (%.3f) between %s and %s - possible mislabeling!\n", correlation, v1, v2))
            }
            if(abs(correlation) < 0.1) {
              cat(sprintf("    INFO: Very low correlation (%.3f) between %s and %s - good spectral separation\n", correlation, v1, v2))
            }
          } else {
            cat(sprintf("  %s vs %s: correlation=NA (insufficient valid data)\n", v1, v2))
          }
        }
      } else {
        cat(sprintf("  %s vs %s: skipped (missing T matrix)\n", v1, v2))
      }
    }
  }
}

cat("Vegetation library built.\n")

timing_info$lib_construction_done <- Sys.time()
cat(sprintf("Library construction completed in %.1f seconds\n", 
           as.numeric(difftime(timing_info$lib_construction_done, timing_info$moving_var_done, units="secs"))))

# Diagnostic: report lib contents and per-veg index coverage
try({
  lib_names <- names(lib)
  cat(sprintf("LIB DIAG: vegs_in_lib=%d\n", length(lib_names)))
  if(length(lib_names) > 0) cat(sprintf("LIB DIAG: names=%s\n", paste(lib_names, collapse = ",")))
  for(vn in lib_names) {
    idxs <- names(lib[[vn]])
    counts <- sapply(idxs, function(ix) sum(is.finite(lib[[vn]][[ix]]$mu)), USE.NAMES = TRUE)
    cat(sprintf("LIB DIAG: veg=%s indices=%d (%s)\n", vn, length(idxs), paste(sprintf("%s:%d", idxs, counts), collapse=",")))
  }
}, silent = TRUE)

if(length(avail) >= 2) {
  # Assemble per-vegetation mu matrices with common column order = avail
  M_list <- list()
  for(v in names(lib)) {
    Mv <- matrix(NA_real_, nrow = 365, ncol = length(avail))
    colnames(Mv) <- avail
    for(j in seq_along(avail)) {
      idx <- avail[j]
      if(!is.null(lib[[v]][[idx]])) Mv[, j] <- lib[[v]][[idx]]$mu
    }
    # NA handling deferred: we'll perform DOY-based interpolation after
    # moving-variance columns have been appended so interpolation covers
    # both spectral indices and moving-variance features.
    # Add moving variances
    for(j in seq_along(avail)) {
      idx <- avail[j]
      if(!is.null(lib[[v]][[idx]]) && !is.null(lib[[v]][[idx]]$mv)) {
        mv_col <- paste0(idx, "_mv")
        Mv <- cbind(Mv, lib[[v]][[idx]]$mv)
        colnames(Mv)[ncol(Mv)] <- mv_col
      }
    }
    # Perform time interpolation for any remaining non-finite values
    # across the 365-day DOY axis. This uses linear interpolation with
    # wrap-around (periodic extension) to avoid edge artifacts near
    # the start/end of year. Fallbacks:
    #  - 0 finite values -> fill with 0 (robust default)
    #  - 1 finite value  -> replicate that value for all days
    #  - >=2 finite values -> interpolate with circular extension
    days <- seq_len(365)
    for(i in seq_len(ncol(Mv))) {
      col <- as.numeric(Mv[, i])
      if(any(!is.finite(col))) {
        finite_idx <- which(is.finite(col))
        if(length(finite_idx) == 0) {
          # No data: fallback to 0 (matches previous median fallback default)
          col[] <- 0
        } else if(length(finite_idx) == 1) {
          # Single observation: replicate across year
          col[] <- col[finite_idx]
        } else {
          # Use periodic extension to interpolate across year boundary
          x <- c(finite_idx - 365, finite_idx, finite_idx + 365)
          y <- rep(col[finite_idx], 3)
          interp <- tryCatch({
            stats::approx(x = x, y = y, xout = days, rule = 2)$y
          }, error = function(e) rep(median(col[finite_idx], na.rm = TRUE), length(days)))
          # Ensure numeric and finite
          interp <- as.numeric(interp)
          if(!all(is.finite(interp))) interp[!is.finite(interp)] <- median(col[finite_idx], na.rm = TRUE)
          col <- interp
        }
        Mv[, i] <- col
      }
    }

    M_list[[v]] <- Mv
  }

  # ---- Shape-based normalization scale per vegetation type and index ----
  # Compute normalization scales using lowest 50-day and top 30-day periods for each index
  BAND_SCALE <- list()
  for(vn in names(lib)) {
    # Compute a single DVI-based scale per vegetation and apply to all indices
    BAND_SCALE[[vn]] <- list()
    # Default fallback
    dvi_scale_value <- 1.0
    if(!is.null(lib[[vn]][["DVI"]]) && !is.null(lib[[vn]][["DVI"]]$mu) && length(lib[[vn]][["DVI"]]$mu) >= 1) {
      ts_dvi <- lib[[vn]][["DVI"]]$mu
      n_days_dvi <- length(ts_dvi)
      # Find lowest 50-day period (minimum median of 50-day windows)
      if(n_days_dvi >= 50) {
        window_medians_low <- sapply(1:(n_days_dvi - 49), function(i) median(ts_dvi[i:(i+49)], na.rm = TRUE))
        lowest_50_median <- min(window_medians_low, na.rm = TRUE)
      } else {
        lowest_50_median <- median(ts_dvi, na.rm = TRUE)
      }
      # Find top 30-day period (maximum median of 30-day windows)
      if(n_days_dvi >= 30) {
        window_medians_high <- sapply(1:(n_days_dvi - 29), function(i) median(ts_dvi[i:(i+29)], na.rm = TRUE))
        top_30_median <- max(window_medians_high, na.rm = TRUE)
      } else {
        top_30_median <- median(ts_dvi, na.rm = TRUE)
      }
      dvi_scale_value <- median(c(lowest_50_median, top_30_median), na.rm = TRUE)
      if(!is.finite(dvi_scale_value) || dvi_scale_value <= 1e-6) dvi_scale_value <- 1.0
    }

    # Apply the DVI-based scale to all available indices (including _mv variants)
    for(idx in avail) {
      BAND_SCALE[[vn]][[idx]] <- dvi_scale_value
    }
  }
  # Log the scales
  for(vn in names(BAND_SCALE)) {
    scale_str <- paste(sapply(avail, function(idx) {
      sprintf("%s=%.4f", idx, BAND_SCALE[[vn]][[idx]])
    }), collapse = ", ")
    cat("Shape-based normalization scales for", vn, ":", scale_str, "\n", sep = "")
  }

  # DEBUG: Check for potential normalization bias
  if(length(names(BAND_SCALE)) > 1) {
    cat("DEBUG_NORMALIZATION: Checking for scale differences between vegetation types\n")
    for(idx in avail) {
      scales_for_idx <- sapply(names(BAND_SCALE), function(vn) BAND_SCALE[[vn]][[idx]])
      scale_range <- max(scales_for_idx, na.rm=TRUE) / min(scales_for_idx, na.rm=TRUE)
      if(is.finite(scale_range) && scale_range > 2) {
        cat(sprintf("DEBUG_NORMALIZATION: Large scale difference for %s: range=%.2f (min=%.4f, max=%.4f)\n",
                   idx, scale_range, min(scales_for_idx, na.rm=TRUE), max(scales_for_idx, na.rm=TRUE)))
      }
    }
  }

  # ---- Global shape projection: compute PCA and LDA once on all vegetation prototypes ----
  V_names <- names(M_list)
  if(length(V_names) < 2) {
    cat("Shape projection skipped: need >= 2 vegetation classes with data.\n")
  } else {
    K_idx <- length(avail)
    if(K_idx < 1) stop("No indices available for shape projection")

    # Build stacked data matrix X_all (rows = veg * 365, cols = all_features) and labels
    X_blocks <- list(); y_lbls <- character(0)
    expected_cols <- NULL
    
    # First pass: determine the expected number of columns from the first valid matrix
    for(vn in names(M_list)) {
      Mv <- M_list[[vn]]
      if(!is.null(Mv) && is.matrix(Mv)) {
        expected_cols <- ncol(Mv)
        break
      }
    }
    
    if(is.null(expected_cols)) {
      cat("Shape projection skipped: no valid vegetation matrices found.\n")
    } else {
      # Second pass: collect all matrices, ensuring they have the same number of columns
      for(vn in names(M_list)) {
        Mv <- M_list[[vn]]
        if(is.null(Mv) || !is.matrix(Mv)) next
        
        # All vegetation matrices should have the same structure (same column names)
        if(ncol(Mv) != expected_cols) {
          stop(sprintf("Column count mismatch for vegetation '%s': expected %d, got %d", 
                      vn, expected_cols, ncol(Mv)))
        }
        
        X_blocks[[length(X_blocks)+1]] <- Mv
        y_lbls <- c(y_lbls, rep(vn, nrow(Mv)))
      }
    }
      if(length(X_blocks) == 0) {
        stop("Shape projection requires per-vegetation matrices but none were assembled")
      } else {
        # Strictly assemble stacked matrix from per-vegetation blocks; fail hard on any error
  X_all <- chunked_rbind(X_blocks, chunk_size = 25L)
  gc()
        if(!is.matrix(X_all) || nrow(X_all) < 2) {
          stop("Failed to assemble stacked X_all matrix from X_blocks; check dimensions of M_list entries")
        }

        # Ensure X_all has proper column names
        if(is.null(colnames(X_all))) {
          stop("X_all matrix has no column names - check M_list matrix construction")
        }

  # Quadratic terms removed
  avail_aug <- colnames(X_all)

        # Vegetation-specific standardization: center globally but scale per vegetation using shape-based normalization
        mu_all <- colMeans(X_all, na.rm = TRUE)
        Xc <- sweep(X_all, 2, mu_all, "-")
        
        # Use shape-based normalization scales (BAND_SCALE) instead of standard deviation
        # Convert BAND_SCALE from list of lists to matrix format for easier application
        band_scale_matrix <- matrix(1, nrow = length(V_names), ncol = length(avail_aug))
        rownames(band_scale_matrix) <- V_names
        colnames(band_scale_matrix) <- avail_aug
        
        for(v_idx in seq_along(V_names)) {
          vn <- V_names[v_idx]
          if(vn %in% names(BAND_SCALE)) {
            for(idx in avail) {
              # Find the column index for this index in avail_aug
              col_name <- idx
              if(col_name %in% colnames(band_scale_matrix)) {
                scale_val <- BAND_SCALE[[vn]][[idx]]
                if(is.finite(scale_val) && scale_val > 0) {
                  band_scale_matrix[v_idx, col_name] <- scale_val
                }
              }
              
              # Handle moving variance terms
              col_name_mv <- paste0(idx, "_mv")
              if(col_name_mv %in% colnames(band_scale_matrix)) {
                scale_val <- BAND_SCALE[[vn]][[idx]]
                if(is.finite(scale_val) && scale_val > 0) {
                  band_scale_matrix[v_idx, col_name_mv] <- scale_val
                }
              }
            }
          }
        }
        
        # Apply vegetation-specific scaling using BAND_SCALE
        Xs <- Xc
        for(v_idx in seq_along(V_names)) {
          vn <- V_names[v_idx]
          start_row <- (v_idx - 1) * n_days + 1
          end_row <- v_idx * n_days
          if(end_row <= nrow(Xs)) {
            veg_scales <- band_scale_matrix[v_idx, ]
            Xs[start_row:end_row, ] <- sweep(Xc[start_row:end_row, , drop = FALSE], 2, veg_scales, "/")
          }
        }
        Xs[!is.finite(Xs)] <- 0
        
        sv <- try(svd(Xs, nu = 0, nv = min(ncol(Xs), nrow(Xs))), silent = TRUE)
        if(inherits(sv, "try-error") || length(sv$d) == 0) {
          stop("PCA failed: SVD did not produce usable singular values; cannot build GLOBAL_PCA.")
        } else {
          variance_explained <- sv$d^2 / sum(sv$d^2)
          cumulative_variance <- cumsum(variance_explained)
          keep_idx <- which(cumulative_variance >= 0.70)
          
          if(length(keep_idx) == 0) {
            stop("Cumulative variance criterion (70%) not met; insufficient data for PCA")
          } else {
            pca_rank <- min(keep_idx[1], MAX_FACTORS_CAP, ncol(sv$v))  
          }
          
          if(pca_rank == 0) stop("PCA selection (70% variance) resulted in 0 factors; cannot build projection.")

          V_pca <- sv$v[, seq_len(pca_rank), drop = FALSE]
          if(ncol(V_pca) > 1) V_pca <- qr.Q(qr(V_pca))[, seq_len(ncol(V_pca)), drop = FALSE]
        }
        
        if(isTRUE(SHAPE_LDA_ENABLE)) {
          # Fail hard if there are not enough classes to perform LDA
          if(length(unique(y_lbls)) < 2) {
            stop("CRITICAL: LDA fit skipped due to insufficient classes (need >= 2). Check input data and Veg labels.")
          }
          # Fail hard if MASS is not available
          if(!requireNamespace("MASS", quietly = TRUE)) {
            stop("CRITICAL: LDA fit skipped because MASS package is not installed.")
          }

          # Use per-vegetation scaling like PCA instead of shape-based normalization
          # This prevents constant variables within groups that cause LDA failures
          lda_dat <- data.frame(Xc, .cls = factor(y_lbls))

          # Apply per-vegetation scaling to prevent constant variables within groups
          for(v_idx in seq_along(V_names)) {
            vname <- V_names[v_idx]
            veg_rows <- y_lbls == vname
            if(any(veg_rows)) {
              # Use the same per-vegetation scaling as PCA
              veg_scales <- band_scale_matrix[v_idx, ]
              lda_dat[veg_rows, seq_len(ncol(Xc))] <- sweep(Xc[veg_rows, , drop = FALSE], 2, veg_scales, "/")
            }
          }
          # Fix non-finite values in spectral columns using matrix conversion
          lda_matrix_part <- as.matrix(lda_dat[, seq_len(ncol(Xc))])
          lda_matrix_part[!is.finite(lda_matrix_part)] <- 0
          lda_dat[, seq_len(ncol(Xc))] <- lda_matrix_part

          # Remove any remaining constant variables (should be minimal now)
          var_info <- sapply(colnames(lda_dat)[-ncol(lda_dat)], function(var_name) {
            if(var_name == ".cls") return(TRUE)
            var_data <- lda_dat[[var_name]]
            if(!is.numeric(var_data)) return(TRUE)

            # Check if variable is constant within any group
            group_vars <- tapply(var_data, lda_dat$.cls, function(x) length(unique(x)))
            any(group_vars <= 1)  # TRUE if constant in any group
          })

          constant_vars <- names(var_info)[!var_info]
          if(length(constant_vars) > 0) {
            cat(sprintf("DEBUG_LDA: Removing %d constant variables within groups: %s\n",
                       length(constant_vars), paste(constant_vars, collapse = ", ")))
            lda_dat <- lda_dat[, !(colnames(lda_dat) %in% constant_vars), drop = FALSE]
          }

          # Additional check: remove variables with very low variance across all data
          if(ncol(lda_dat) > 2) {  # Keep at least .cls column + 1 variable
            numeric_cols <- sapply(lda_dat, is.numeric)
            if(any(numeric_cols)) {
              variances <- sapply(lda_dat[, numeric_cols, drop = FALSE], var, na.rm = TRUE)
              low_var_threshold <- 1e-10  # Very small variance threshold
              low_var_vars <- names(variances)[variances < low_var_threshold]
              if(length(low_var_vars) > 0) {
                cat(sprintf("DEBUG_LDA: Removing %d variables with very low variance (< %g): %s\n",
                           length(low_var_vars), low_var_threshold, paste(low_var_vars, collapse = ", ")))
                lda_dat <- lda_dat[, !(colnames(lda_dat) %in% low_var_vars), drop = FALSE]
              }
            }
          }

          # Final check: ensure we have enough variables and they're not all constant
          if(ncol(lda_dat) < 3) {  # Need .cls + at least 1 predictor  
            stop("LDA fit failed: insufficient non-constant variables after preprocessing")
          } else {
            lda_fit <- tryCatch({
              MASS::lda(.cls ~ ., data = lda_dat)
            }, error = function(e) {
              # Enhanced error handling for LDA failures
              error_msg <- as.character(e$message)
              if(grepl("constant within groups", error_msg)) {
                cat(sprintf("DEBUG_LDA: LDA failed due to constant variables. Error: %s\n", error_msg))
                cat("DEBUG_LDA: This should not happen with per-vegetation scaling - check scaling computation\n")
              } else {
                if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
                diagf <- file.path(OUT_DIR, sprintf('lda_fit_error_%s.txt', format(Sys.time(), "%Y%m%d_%H%M%S")))
                cat(sprintf('LDA_FIT_ERROR: msg=%s\n', error_msg), file = diagf)
              }
              return(NULL)  # Return NULL instead of stopping
            })

            if(is.null(lda_fit)) {
              stop("LDA fit failed - cannot proceed without discriminant analysis")
            } else {
              # Apply 70% variance rule to LDA
              prop_trace <- lda_fit$svd^2 / sum(lda_fit$svd^2)
              cum_trace <- cumsum(prop_trace)
              n_lda <- min(which(cum_trace >= 0.70)[1], length(lda_fit$svd))
              if (is.na(n_lda)) n_lda <- length(lda_fit$svd)

              lda_fit$scaling <- lda_fit$scaling[, 1:n_lda, drop = FALSE]

              GLOBAL_LDA <- list(fit = lda_fit, col_means = mu_all, col_scales = band_scale_matrix, idx_order = avail_aug)
              cat(sprintf("LDA fit successful: %d components retained from %d variables\n", n_lda, ncol(lda_dat) - 1))
            }
          }
        }
        }
        # ---- Alternative Dimension Reduction: Kernel PCA ----
        # Use when spectral relationships are non-linear (e.g., vegetation stress responses)
        if(exists("ENABLE_KERNEL_PCA") && isTRUE(ENABLE_KERNEL_PCA)) {
          if(!requireNamespace("kernlab", quietly = TRUE)) {
            stop("kernlab package required for Kernel PCA; install with: install.packages('kernlab')")
          } else {
            cat("Computing Kernel PCA for non-linear spectral relationships...\n")

            # First, reduce spectral dimension with kernel PCA
            spectral_features <- min(MAX_FACTORS_CAP, ncol(Xs))
            kpca_spectral <- tryCatch({
              kernlab::kpca(Xs, kernel = "rbfdot", kpar = list(sigma = 0.1), features = spectral_features)
            }, error = function(e) {
              stop(sprintf("Kernel PCA spectral reduction failed: %s", e$message))
            })

            if(!is.null(kpca_spectral)) {
              # Extract spectral kernel principal components
              KPC_spectral <- kernlab::rotated(kpca_spectral)

              if(is.matrix(KPC_spectral) && ncol(KPC_spectral) > 0) {
                # The KPC_spectral matrix has dimensions (n_veg * n_days, n_spectral)
                # Data is structured as: [veg1_day1:365, veg2_day1:365, ...]
                n_veg <- length(V_names)
                n_days <- 365
                n_spectral <- ncol(KPC_spectral)
                
                # Verify the input dimensions match what we expect
                expected_rows <- n_veg * n_days
                actual_rows <- nrow(KPC_spectral)
                
                if(actual_rows != expected_rows) {
                  stop(sprintf("KPCA dimension mismatch: KPC_spectral has %d rows, expected %d (%d veg * %d days). Check vegetation library construction.", 
                              actual_rows, expected_rows, n_veg, n_days))
                }
                
                # Reshape: KPC_spectral is structured as blocks [veg1_days, veg2_days, ...]
                # Convert from (n_veg * n_days, n_spectral) to (n_veg, n_days, n_spectral)
                KPC_reshaped <- array(dim = c(n_veg, n_days, n_spectral))
                
                # Fill the array by extracting each vegetation block
                for(v_idx in 1:n_veg) {
                  start_row <- (v_idx - 1) * n_days + 1
                  end_row <- v_idx * n_days
                  KPC_reshaped[v_idx, , ] <- KPC_spectral[start_row:end_row, ]
                }
                
                # Set dimnames for the reshaped array
                dimnames(KPC_reshaped) <- list(
                  vegetation = V_names,
                  day_of_year = as.character(1:n_days),
                  spectral_component = paste0("PC", 1:n_spectral)
                )

                # Apply temporal kernel PCA to reduce time dimension
                # Flatten spectral components for each vegetation type and apply temporal KPCA
                temporal_reduction <- 50  # Reduce 365 days to ~50 temporal components (good for phenology)

                KPC_temporal_list <- list()
                for(v_idx in 1:n_veg) {
                  # Extract time series for this vegetation type: (365, n_spectral)
                  veg_timeseries <- KPC_reshaped[v_idx, , , drop = FALSE]
                  veg_timeseries <- matrix(veg_timeseries, nrow = n_days, ncol = n_spectral)

                  # Apply kernel PCA on temporal dimension - fail if unsuccessful
                  kpca_temp <- kernlab::kpca(veg_timeseries, kernel = "rbfdot", kpar = list(sigma = 0.05),
                                           features = min(temporal_reduction, n_days))
                  
                  KPC_temporal_list[[v_idx]] <- kernlab::rotated(kpca_temp)
                }

                # Combine temporal components across vegetation types
                if(length(KPC_temporal_list) > 0) {
                  # Stack all vegetation temporal components
                  KPC_combined <- do.call(rbind, KPC_temporal_list)

                  GLOBAL_KPCA <- list(
                    fit = kpca_spectral,  # Keep spectral KPCA for reference
                    rotated = KPC_combined,
                    col_means = mu_all,
                    col_scales = band_scale_matrix,
                    idx_order = avail_aug,
                    rank = ncol(KPC_combined),
                    temporal_components = temporal_reduction,
                    spectral_components = n_spectral
                  )
                  cat(sprintf("Kernel PCA computed: %d spectral + %d temporal components (%d total)\n",
                              n_spectral, temporal_reduction, ncol(KPC_combined)))
                } else {
                  GLOBAL_KPCA <- NULL
                }
              } else {
                GLOBAL_KPCA <- NULL
              }
            } else {
              GLOBAL_KPCA <- NULL
            }
          }
        } else {
          GLOBAL_KPCA <- NULL
        }
        
        GLOBAL_PCA <- list(
          global = list(
            V = if(is.matrix(V_pca)) V_pca else matrix(0, nrow = length(avail_aug), ncol = 0),
            col_means = mu_all,
            col_scales = band_scale_matrix,  # Now using shape-based normalization
            rank = if(is.matrix(V_pca)) ncol(V_pca) else 0
          ),
          V = if(is.matrix(V_pca)) V_pca else matrix(0, nrow = length(avail_aug), ncol = 0),
          col_means = mu_all,
          col_scales = band_scale_matrix,  # Now using shape-based normalization
          per_bin = NULL,
          rank = if(is.matrix(V_pca)) ncol(V_pca) else 0,
          idx_order = avail_aug
        )

        # Build per-vegetation daily factor time series for PCA and LDA
        # New approach: project original training rows for each vegetation into
        # the global PCA/LDA spaces, then compute per-DOY medians of the
        # component scores. Also compute a per-DOY robust variance scalar
        # (mean of component-wise MAD^2) to be used as T_var.
        lib_factor_pca <- list(); lib_factor_lda <- list()
        gpca_order <- GLOBAL_PCA$idx_order
        gpca_means <- GLOBAL_PCA$col_means
        gpca_scales <- GLOBAL_PCA$col_scales
        Vpca <- GLOBAL_PCA$global$V
        pca_rank <- GLOBAL_PCA$global$rank

        has_lda <- !is.null(GLOBAL_LDA) && is.list(GLOBAL_LDA) && !is.null(GLOBAL_LDA$fit)
        if(has_lda) {
          glda_order <- GLOBAL_LDA$idx_order
          glda_means <- GLOBAL_LDA$col_means
          glda_scales <- GLOBAL_LDA$col_scales
          lda_scaling <- GLOBAL_LDA$fit$scaling
        }

        for(vname in names(M_list)) {
          # Gather training rows for this vegetation
          dveg <- lib_df[lib_df$Veg == vname & !is.na(lib_df$date), , drop = FALSE]
          n_samples_v <- if(!is.null(dveg)) nrow(dveg) else 0
          if(n_samples_v <= 0) {
            # Stop execution - no fallback allowed
            lib_factor_pca[[vname]] <- list(T = matrix(0, nrow = 365, ncol = max(0, pca_rank)), T_var = rep(EPS_SIGMA, 365))
            lib_factor_pca[[vname]]$n_samples <- 0
            if(has_lda) {
              lib_factor_lda[[vname]] <- list(T = matrix(0, nrow = 365, ncol = ncol(lda_scaling)), T_var = rep(EPS_SIGMA, 365))
              lib_factor_lda[[vname]]$n_samples <- 0
            }
            next
          }

          # Prepare storage for projected scores and DOY
          # PCA
          proj_list <- list()
          proj_doy <- integer(0)

          # Loop through training rows and build augmented feature vectors
          for(i_row in seq_len(nrow(dveg))) {
            row <- dveg[i_row, , drop = FALSE]
            doy_row <- as.integer(lubridate::yday(row$date))
            if(is.na(doy_row) || doy_row < 1 || doy_row > 365) next

            # Build augmented raw vector following gpca_order
            raw_vec <- rep(NA_real_, length(gpca_order))
            names(raw_vec) <- gpca_order

            # Linear terms
            for(idx in avail) {
              kpos <- match(idx, gpca_order)
              if(is.na(kpos)) next
              val <- row[[idx]]
              if(!is.finite(val)) val <- gpca_means[kpos]
              raw_vec[kpos] <- as.numeric(val)
            }

            # Quadratic terms removed

            # Moving-variance terms (use existing _var14 columns if present)
            for(idx in avail) {
              mv_name <- paste0(idx, "_mv")
              kpos_mv <- match(mv_name, gpca_order)
              if(is.na(kpos_mv)) next
              mv_col <- paste0(idx, "_var14")
              mv_val <- if(mv_col %in% names(dveg)) row[[mv_col]] else NA_real_
              if(!is.finite(mv_val)) mv_val <- gpca_means[kpos_mv]
              raw_vec[kpos_mv] <- as.numeric(mv_val)
            }

            # Use time interpolation instead of global means for remaining NAs
            nas <- which(!is.finite(raw_vec))
            if(length(nas) > 0) {
              # For interpolation, collect all available data for this vegetation type
              valid_rows <- dveg[!is.na(dveg$date), , drop = FALSE]
              if(nrow(valid_rows) > 1) {
                valid_doys <- as.integer(lubridate::yday(valid_rows$date))
                
                for(na_idx in nas) {
                  var_name <- gpca_order[na_idx]
                  # Find corresponding column in original data
                  if(var_name %in% names(valid_rows)) {
                    valid_vals <- valid_rows[[var_name]]
                    finite_mask <- is.finite(valid_vals) & is.finite(valid_doys)
                    
                    if(sum(finite_mask) >= 2) {
                      # Use linear interpolation based on DOY
                      interp_val <- tryCatch({
                        stats::approx(x = valid_doys[finite_mask], 
                                     y = valid_vals[finite_mask], 
                                     xout = doy_row, 
                                     method = "linear", 
                                     rule = 2)$y
                      }, error = function(e) gpca_means[na_idx])
                      
                      if(is.finite(interp_val)) {
                        raw_vec[na_idx] <- as.numeric(interp_val)
                      } else {
                        raw_vec[na_idx] <- gpca_means[na_idx]
                      }
                    } else {
                      # Fall back to global mean if insufficient data
                      raw_vec[na_idx] <- gpca_means[na_idx]
                    }
                  } else {
                    # Fall back to global mean if variable not found
                    raw_vec[na_idx] <- gpca_means[na_idx]
                  }
                }
              } else {
                # Fall back to global means if insufficient temporal data
                raw_vec[nas] <- gpca_means[nas]
              }
            }

            # Standardize and project if PCA rank > 0
            if(pca_rank > 0 && ncol(Vpca) >= pca_rank) {
              centered <- (raw_vec - gpca_means) / gpca_scales[vname, ]
              centered[!is.finite(centered)] <- 0
              t_scores <- as.numeric(centered %*% Vpca)
            } else {
              t_scores <- numeric(0)
            }

            proj_list[[length(proj_list) + 1]] <- t_scores
            proj_doy <- c(proj_doy, doy_row)
          }

          # Assemble KDE-based probabilistic prototypes for PCA
          if(length(proj_list) > 0 && pca_rank > 0) {
            obs_mat <- chunked_rbind(proj_list, chunk_size = 50L)
            gc()
            # obs_mat rows = observations, cols = components
            T_pca_med <- matrix(NA_real_, nrow = 365, ncol = ncol(obs_mat))
            T_pca_kde <- rep(NA_real_, 365)
            doy_grid <- seq_len(365)
            
            # Use median projected values per DOY (no GAM smoothing)
            for(comp in seq_len(ncol(obs_mat))) {
              for(d in seq_len(365)) {
                idx <- which(proj_doy == d)
                if(length(idx) > 0) {
                  T_pca_med[d, comp] <- median(obs_mat[idx, comp], na.rm = TRUE)
                } else {
                  T_pca_med[d, comp] <- 0
                }
              }
            }
            
            # Fit KDE to DOY distribution for this vegetation type
            if(length(unique(proj_doy)) < 3) {
              stop("Insufficient data for KDE computation: need at least 3 unique DOY values, got ", 
                   length(unique(proj_doy)), " for vegetation type ", vname)
            }
            
            kde_fit <- density(proj_doy, bw = "SJ", from = 1, to = 365, n = 365)
            T_pca_kde <- kde_fit$y
            # Normalize KDE densities
            T_pca_kde <- T_pca_kde / sum(T_pca_kde)
            
            lib_factor_pca[[vname]] <- list(T = T_pca_med, T_var = T_pca_kde)
            lib_factor_pca[[vname]]$n_samples <- n_samples_v
          } else {
            lib_factor_pca[[vname]] <- list(T = matrix(0, nrow = 365, ncol = pca_rank), T_var = rep(1/365, 365))
            lib_factor_pca[[vname]]$n_samples <- n_samples_v
          }

          # LDA prototypes: similar flow but using GLOBAL_LDA projection
            if(has_lda) {
            # Build projections using same raw vectors but standardize with per-vegetation scaling like PCA
            proj_list_lda <- list(); proj_doy_lda <- integer(0)
            for(i_row in seq_len(nrow(dveg))) {
              row <- dveg[i_row, , drop = FALSE]
              doy_row <- as.integer(lubridate::yday(row$date))
              if(is.na(doy_row) || doy_row < 1 || doy_row > 365) next

              raw_vec <- rep(NA_real_, length(glda_order)); names(raw_vec) <- glda_order
              for(idx in avail) {
                kpos <- match(idx, glda_order)
                if(is.na(kpos)) next
                val <- row[[idx]]
                if(!is.finite(val)) val <- glda_means[kpos]
                raw_vec[kpos] <- as.numeric(val)
              }
              for(idx in avail) {
                # Quadratic terms removed
              }
              for(idx in avail) {
                mv_name <- paste0(idx, "_mv")
                kpos_mv <- match(mv_name, glda_order)
                if(is.na(kpos_mv)) next
                mv_col <- paste0(idx, "_var14")
                mv_val <- if(mv_col %in% names(dveg)) row[[mv_col]] else NA_real_
                if(!is.finite(mv_val)) mv_val <- glda_means[kpos_mv]
                raw_vec[kpos_mv] <- as.numeric(mv_val)
              }
              nas <- which(!is.finite(raw_vec))
              if(length(nas) > 0) raw_vec[nas] <- glda_means[nas]

              # Use per-vegetation scaling like PCA (not global scaling)
              centered <- (raw_vec - glda_means) / glda_scales[vname, ]
              centered[!is.finite(centered)] <- 0
              t_scores_lda <- as.numeric(centered %*% lda_scaling)

              proj_list_lda[[length(proj_list_lda) + 1]] <- t_scores_lda
              proj_doy_lda <- c(proj_doy_lda, doy_row)
            }
          }

          # LDA prototypes: create timeseries prototypes like PCA
          if(has_lda && length(proj_list_lda) > 0) {
            obs_lda <- do.call(rbind, proj_list_lda)
            T_lda_med <- matrix(NA_real_, nrow = 365, ncol = ncol(obs_lda))
            T_lda_kde <- rep(NA_real_, 365)
            
            # KDE-based density estimation for LDA prototypes
            if(length(unique(proj_doy_lda)) < 3) {
              stop("Insufficient LDA data for KDE computation: need at least 3 unique DOY values, got ", 
                   length(unique(proj_doy_lda)), " for vegetation type ", vname)
            }
            
            kde_lda <- density(proj_doy_lda, bw = "SJ", from = 1, to = 365, n = 365)
            T_lda_kde <- kde_lda$y
            T_lda_kde <- T_lda_kde / sum(T_lda_kde)  # Normalize to sum to 1
            
            # Use median per DOY for T matrix (no GAM smoothing)
            for(comp in seq_len(ncol(obs_lda))) {
              for(d in 1:365) {
                doy_vals <- obs_lda[proj_doy_lda == d, comp]
                T_lda_med[d, comp] <- if(length(doy_vals) > 0) median(doy_vals, na.rm = TRUE) else 0
              }
            }
            
            lib_factor_lda[[vname]] <- list(T = T_lda_med, T_var = T_lda_kde)
            lib_factor_lda[[vname]]$n_samples <- n_samples_v
          } else if(has_lda) {
            stop("LDA projections failed for vegetation type: ", vname, " - no valid projections computed")
          } # end has_lda
        }
      }
    }
  }
cat("Global projections computed.\n")

timing_info$pca_computation_done <- Sys.time()
cat(sprintf("PCA/LDA computation completed in %.1f seconds\n", 
           as.numeric(difftime(timing_info$pca_computation_done, timing_info$lib_construction_done, units="secs"))))

FACTOR_MODE <- exists("GLOBAL_PCA") && !is.null(GLOBAL_PCA) && is.list(GLOBAL_PCA) && !is.null(GLOBAL_PCA$rank) && GLOBAL_PCA$rank >= 1L

# Enforce factor-mode only; abort with a clear message if projection couldn't be built
if(!FACTOR_MODE) {
  stop("Between-veg factor projection unavailable; ensure >=2 vegetation classes with data to build the discriminative projection (GLOBAL_PCA missing or rank<1).")
}

# Helper to compute dynamic ridge penalty based on condition number
compute_dynamic_ridge <- function(cond, k, min_ridge, max_ridge) {
  if(!is.finite(cond) || is.na(cond)) return(min_ridge)
  ridge <- k / cond
  pmin(max_ridge, pmax(min_ridge, ridge))
}
## GCV-based ridge selection: choose ridge (lambda) that minimises generalized cross-validation
select_ridge_gcv <- function(X, y, lambdas = exp(seq(log(1e-6), log(1e6), length.out = 80)),
                             min_ridge = 0.01, max_ridge = 1e2) {
  # X: design matrix (n x p), y: response vector (length n)
  X <- tryCatch(as.matrix(X), error = function(e) stop(sprintf("select_ridge_gcv: failed to coerce X to matrix: %s", e$message)))
  y <- tryCatch(as.numeric(y), error = function(e) stop(sprintf("select_ridge_gcv: failed to coerce y to numeric: %s", e$message)))
  if(length(y) != nrow(X) || nrow(X) < 1) stop("select_ridge_gcv: invalid input shapes - length(y) must match nrow(X) and be >= 1")

  svdX <- tryCatch(svd(X, nu = min(nrow(X), ncol(X)), nv = 0), error = function(e) stop(sprintf("select_ridge_gcv: SVD failed on X: %s", e$message)))
  d <- svdX$d
  U <- svdX$u
  # projected responses onto left singular vectors
  Uty <- as.numeric(t(U) %*% y)
  n <- nrow(X)

  gcv_vals <- rep(NA_real_, length(lambdas))
  for(i in seq_along(lambdas)) {
    lam <- lambdas[i]
    # trace(S) = sum(d^2/(d^2 + lam))
    trS <- sum((d^2) / (d^2 + lam))
    # residual squared norm: ||(I - S) y||^2
    # contribution from rank components
    res_rank <- sum( ((lam / (d^2 + lam))^2) * (Uty^2) )
    # contribution from orthogonal complement (if any)
    if(n > length(d)) {
      total_y_sq <- sum(y^2)
      res_sq <- res_rank + (total_y_sq - sum(Uty^2))
    } else {
      res_sq <- res_rank
    }
    denom <- (1 - trS / n)^2
    gcv_vals[i] <- ifelse(denom <= 0, Inf, res_sq / denom)
  }

  best <- which.min(gcv_vals)
  if(length(best) == 0 || !is.finite(gcv_vals[best])) return(NA_real_)
  lambda_best <- lambdas[best]
  # enforce bounds
  lambda_best <- pmin(max_ridge, pmax(min_ridge, lambda_best))
  return(as.numeric(lambda_best))
}
# Helper function to prepare factor data for a location-year
prepare_factor_data <- function(dly, gpca, avail_idx, veg_type) {
    date_list <- list()
    dly$doy <- lubridate::yday(dly$date)
    dts <- sort(unique(dly$date))
    for(i_dt in seq_along(dts)) {
        dt <- dts[i_dt]
        sub <- dly[dly$date == dt, , drop = FALSE]
        sub <- sub[is.finite(sub$doy), , drop = FALSE]
        if(nrow(sub) == 0) next
        
    # Get original indices present
    idx_present_orig <- intersect(avail_idx, names(sub))
    # Allow partial index coverage: require a minimum fraction of indices present
    if(!exists("MIN_IDX_PRESENCE")) MIN_IDX_PRESENCE <- 0.5
    min_required <- ceiling(length(avail_idx) * MIN_IDX_PRESENCE)
    if(length(idx_present_orig) < min_required) {
      # Not enough indices available for this date; skip but record for diagnostics
      next
    }

        # Create augmented feature vector with quadratic terms
        vals_aug <- c()
        vrows_list <- list()
        
    # Original terms (use only indices present for this date)
    for(idx in idx_present_orig) {
            vals_idx <- sub[[idx]]
            vals_idx <- vals_idx[is.finite(vals_idx)]
            if(length(vals_idx) == 0) { vals_aug <- NULL; break; }
            
            yy <- suppressWarnings(stats::median(vals_idx, na.rm = TRUE))
            if(!is.finite(yy)) continue
            
            kpos <- match(idx, gpca$idx_order)
            if(is.na(kpos)) continue

            mi <- gpca$col_means[kpos]
            si <- if(is.list(gpca$col_scales)) {
              gpca$col_scales[[veg_type]][kpos]
            } else {
              gpca$col_scales[veg_type, kpos]
            }
            if(!is.finite(si) || si <= 0) si <- 1
            
            vals_aug <- c(vals_aug, as.numeric((yy - mi) / si))
            vrows_list[[length(vrows_list) + 1]] <- gpca$V[kpos, , drop = FALSE]
        }
        if(is.null(vals_aug)) next

    # Quadratic terms removed

    # Moving variance terms (only for indices present)
    for(idx in idx_present_orig) {
            mv_col <- paste0(idx, "_var14")
            if(!mv_col %in% names(sub)) next

            vals_mv <- sub[[mv_col]]
            vals_mv <- vals_mv[is.finite(vals_mv)]
            if(length(vals_mv) == 0) { vals_aug <- NULL; break; }

            yy_mv <- suppressWarnings(stats::median(vals_mv, na.rm = TRUE))
            if(!is.finite(yy_mv)) continue

            mv_name <- paste0(idx, "_mv")
            kpos_mv <- match(mv_name, gpca$idx_order)
            if(is.na(kpos_mv)) continue

            mi_mv <- gpca$col_means[kpos_mv]
            si_mv <- if(is.list(gpca$col_scales)) {
              gpca$col_scales[[veg_type]][kpos_mv]
            } else {
              gpca$col_scales[veg_type, kpos_mv]
            }
            if(!is.finite(si_mv) || si_mv <= 0) si_mv <- 1

            vals_aug <- c(vals_aug, as.numeric((yy_mv - mi_mv) / si_mv))
            vrows_list[[length(vrows_list) + 1]] <- gpca$V[kpos_mv, , drop = FALSE]
        }
        if(is.null(vals_aug)) next

        if(length(vrows_list) == 0) next
        Vsub <- do.call(rbind, vrows_list)
        if(ncol(Vsub) < 1) next
        
        XtX <- crossprod(Vsub)
        Xty <- crossprod(Vsub, vals_aug)
        # Compute condition number from XtX to drive dynamic ridge selection
        cond_val <- tryCatch({
          k <- stats::kappa(XtX)
          if(!is.finite(k) || k <= 0) stop("condition number computation returned non-finite or non-positive kappa")
          as.numeric(k)
        }, error = function(e) {
          svd_vals <- tryCatch(svd(XtX, nu = 0, nv = 0)$d, error = function(e2) stop(sprintf("Failed to compute condition number: %s; SVD failed: %s", e$message, e2$message)))
          if(length(svd_vals) >= 2 && all(is.finite(svd_vals)) && svd_vals[1] > 0 && tail(svd_vals,1) > 0) as.numeric(svd_vals[1] / tail(svd_vals,1)) else stop("SVD produced invalid singular values for condition number calculation")
        })
        # Prefer GCV-based selection using the original Vsub and response vector (vals_aug).
        # Fall back to the previous heuristic if GCV fails.
        ridge_small <- tryCatch({
          rl <- select_ridge_gcv(Vsub, vals_aug, min_ridge = 0.01, max_ridge = 1e2)
          if(!is.finite(rl) || rl <= 0) stop("select_ridge_gcv returned an invalid ridge value")
          rl
        }, error = function(e) stop(sprintf("select_ridge_gcv failed: %s", e$message)))
        XtX_reg <- XtX + ridge_small * diag(ncol(XtX))
        z <- tryCatch({
          solve(XtX_reg, Xty)
        }, silent = TRUE)
        if(inherits(z, "try-error") || any(!is.finite(z))) {
            z <- qr.solve(XtX_reg, Xty)
        }
        if(any(!is.finite(z))) {
            if(!requireNamespace("MASS", quietly = TRUE)) {
                stop("Matrix inversion failed and MASS package not available for generalized inverse")
            }
            z <- as.numeric(MASS::ginv(XtX_reg) %*% Xty)
        }
        if(any(!is.finite(z))) {
            stop("All matrix inversion methods failed - matrix is severely ill-conditioned")
        }
        date_list[[length(date_list)+1]] <- list(doy = as.integer(sub$doy[1]), z = as.numeric(z))
    }
    return(date_list)
}

## KDE helpers removed. Replaced with a simple inverse-variance-weighted (IVW)
## constrained solver that performs weighted least-squares without assuming
## normal residual distributions. This produces non-negative weights and
## uses the existing solve_mixture helper for constraints/ridge handling.
ivw_solve_mixture <- function(y_all, g_proj, date_list, lib_factor, veg_kept, init_weights = NULL) {
  # Initialize counters for spectral normalization summary
  spectral_normalizations <- 0
  normalized_vegetations <- character(0)
  
  # Build B_all (rows = dates, cols = veg_kept) projected through factors
  if(is.null(date_list) || length(date_list) == 0) return(NULL)
  
  # Determine rank from the projection object (PCA or LDA)
  rnk <- if (!is.null(g_proj$rank)) g_proj$rank else ncol(g_proj$fit$scaling)

  F_matrix <- matrix(0, nrow = length(date_list), ncol = length(veg_kept))
  for(i in seq_along(date_list)) {
    doy <- date_list[[i]]$doy
    for(j in seq_along(veg_kept)) {
      vkey <- veg_kept[j]
      proto <- lib_factor[[vkey]]
      if(!is.null(proto) && !is.null(proto$T) && doy > 0 && doy <= nrow(proto$T) && ncol(proto$T) >= rnk) {
        F_matrix[i, j] <- sum(proto$T[doy, seq_len(rnk)] * date_list[[i]]$z[seq_len(rnk)])
      }
    }
  }

  # Use KDE density weights from proto$T_var (normalized densities), else uniform
  # NOTE: Using KDE-based probabilistic weighting instead of inverse variance
  w_all <- rep(1, nrow(F_matrix))
  if(length(veg_kept) > 0) {
    date_var <- rep(NA_real_, nrow(F_matrix))
    for(i in seq_len(nrow(F_matrix))) {
      vars <- numeric(0)
      for(j in seq_along(veg_kept)) {
        vkey <- veg_kept[j]
        proto <- lib_factor[[vkey]]
        if(!is.null(proto) && !is.null(proto$T_var) && length(proto$T_var) >= date_list[[i]]$doy) {
          d <- date_list[[i]]$doy
          if(d >=1 && d <= length(proto$T_var) && is.finite(proto$T_var[d])) {
            # For KDE, T_var contains normalized density values, use directly
            density_val <- proto$T_var[d]
            vars <- c(vars, density_val)
          }
        }
      }
      if(length(vars) > 0) date_var[i] <- mean(vars, na.rm = TRUE)  # mean KDE density across veg types
    }
    date_var[!is.finite(date_var) | date_var <= 0] <- NA_real_
    if(all(is.na(date_var))) {
      w_all <- rep(1, nrow(F_matrix))
    } else {
      # Use KDE densities directly as weights (no inversion needed)
      w_all <- date_var
      w_all[!is.finite(w_all) | w_all <= 0] <- 1

      # ROBUSTNESS: Transform KDE densities to prevent extreme weights
      # Use square root to soften the distribution and prevent binary-like results
      w_all <- sqrt(w_all)

  # ROBUSTNESS: Cap maximum weight to prevent dominance by any single date
  # This prevents bias when some dates have extremely high KDE density
  q95_w <- stats::quantile(w_all, 0.95, na.rm = TRUE)
  max_weight_cap <- as.numeric(q95_w) * 2  # Allow up to 2x the 95th percentile
  w_all <- pmin(w_all, max_weight_cap)

      # Ensure weights sum to a reasonable value and renormalize
      w_all <- w_all / mean(w_all, na.rm = TRUE)
    }
  }

  # DEBUG: Log weighting information to understand bias
  if(length(veg_kept) > 0 && any(w_all != 1) && FALSE) {  # Disabled debug output
    cat(sprintf("DEBUG_WEIGHTING: date_count=%d, mean_weight=%.4f, weight_range=[%.4f, %.4f]\n",
               length(w_all), mean(w_all, na.rm=TRUE), min(w_all, na.rm=TRUE), max(w_all, na.rm=TRUE)))
    # Log which dates have highest/lowest weights
    if(length(w_all) > 1) {
      high_idx <- which.max(w_all)
      low_idx <- which.min(w_all)
      cat(sprintf("DEBUG_WEIGHTING: highest_weight_date=%d (doy=%d, weight=%.4f), lowest_weight_date=%d (doy=%d, weight=%.4f)\n",
                 high_idx, date_list[[high_idx]]$doy, w_all[high_idx],
                 low_idx, date_list[[low_idx]]$doy, w_all[low_idx]))
    }
  }

  sw <- sqrt(w_all)
  # The response is now implicitly 1, as we are solving for weights
  y_w <- rep(1, length(sw)) * sw 
  B_w <- sweep(F_matrix, 1, sw, "*")

  col_scale <- sqrt(colSums(B_w^2)); col_scale[!is.finite(col_scale) | col_scale <= 0] <- 1
  B_ws <- sweep(B_w, 2, col_scale, "/")
  base_mask <- rep(TRUE, ncol(B_ws))
  # Compute condition number from the scaled Gram matrix to drive dynamic ridge
  cond_val <- tryCatch({
    k <- stats::kappa(crossprod(B_ws))
    if(!is.finite(k) || k <= 0) NA_real_ else as.numeric(k)
  }, error = function(e) {
    svd_vals <- tryCatch(svd(crossprod(B_ws), nu = 0, nv = 0)$d, error = function(e) stop(sprintf("Failed to compute SVD for condition number fallback: %s", e$message)))
    if(length(svd_vals) >= 2 && all(is.finite(svd_vals)) && svd_vals[1] > 0 && tail(svd_vals,1) > 0) svd_vals[1] / tail(svd_vals,1) else NA_real_
  })
  # Use GCV to select ridge for the (weighted) design matrix B_ws and response y_w.
  ridge_val <- tryCatch({
    rl <- select_ridge_gcv(B_ws, y_w, min_ridge = 0.01, max_ridge = 1e2)
    if(!is.finite(rl) || rl <= 0) rl <- compute_dynamic_ridge(cond_val, 1e10, 0.01, 1e2)
    rl
  }, error = function(e) compute_dynamic_ridge(cond_val, 1e10, 0.01, 1e2))
  
  # Use the new gradient-descent-based solver for ridge-constrained solution (no legacy solve_mixture)
  sol <- solve_mixture_gd(B_ws, y_w, base_mask = base_mask, lambda2 = ridge_val)
  coef_ridge <- NULL
  if(!is.null(sol) && !is.null(sol$X)) {
    coef_ridge_raw <- as.numeric(sol$X) / col_scale
    names(coef_ridge_raw) <- veg_kept
    # enforce simplex via optimize mapping (no post-fit normalization)
    target_obs <- sapply(date_list, function(d) sqrt(sum(d$z^2)))
    coef_ridge <- optimize_weights_simplex(F_matrix, w_all, target = target_obs, init = coef_ridge_raw)
    # optimize_weights_simplex already ensures simplex constraint
  }

  # Try gradient-descent solver as alternative, then map to simplex
  coef_gd <- NULL
  try({
    sol_gd <- solve_mixture_gd(B_ws, y_w, base_mask = base_mask, lambda2 = ridge_val)
    if(!is.null(sol_gd) && !is.null(sol_gd$X)) {
      coef_gd_raw <- as.numeric(sol_gd$X) / col_scale
      names(coef_gd_raw) <- veg_kept
      target_obs <- sapply(date_list, function(d) sqrt(sum(d$z^2)))
      coef_gd <- optimize_weights_simplex(F_matrix, w_all, target = target_obs, init = coef_gd_raw)
      # optimize_weights_simplex already ensures simplex constraint
    }
  }, silent = TRUE)

  # Score function (weighted SSE)
  sse_calc <- function(w) {
    # Build block design matrix for vector reconstruction
    n_dates <- nrow(F_matrix)
    n_veg <- ncol(F_matrix)
    n_comp <- rnk

    # X_block: (n_dates * n_comp) x n_veg
    X_block <- matrix(0, nrow = n_dates * n_comp, ncol = n_veg)
    # y_block: stacked observed factor vectors
    y_block <- numeric(n_dates * n_comp)
    # w_block: expanded weights for each component
    w_block <- numeric(n_dates * n_comp)

    for(i in seq_len(n_dates)) {
      for(j in seq_len(n_veg)) {
        vkey <- veg_kept[j]
        proto <- lib_factor[[vkey]]
        if(!is.null(proto) && !is.null(proto$T) && date_list[[i]]$doy > 0 &&
           date_list[[i]]$doy <= nrow(proto$T) && ncol(proto$T) >= n_comp) {
          # Fill X_block: each block of n_comp rows for date i, veg j
          row_start <- (i-1) * n_comp + 1
          row_end <- i * n_comp
          X_block[row_start:row_end, j] <- proto$T[date_list[[i]]$doy, seq_len(n_comp)]
        }
      }
      # Fill y_block: observed factor vector for date i
      y_block[((i-1)*n_comp + 1):(i*n_comp)] <- date_list[[i]]$z[seq_len(n_comp)]
      # Fill w_block: repeat weight for each component
      w_block[((i-1)*n_comp + 1):(i*n_comp)] <- w_all[i]
    }

    # Vector reconstruction objective: minimize sum(w_block * (X_block %*% w - y_block)^2)
    pred_block <- X_block %*% w
    res_block <- pred_block - y_block
    s <- sum(w_block * (res_block^2), na.rm = TRUE)
    if(!is.finite(s)) s <- .Machine$double.xmax
    s
  }

  candidates <- list()
  if(!is.null(coef_ridge)) candidates[["ridge"]] <- coef_ridge
  if(!is.null(coef_gd)) candidates[["gd"]] <- coef_gd

  best_name <- NULL; best_score <- Inf; best_coef <- NULL
  for(nm in names(candidates)) {
    w <- candidates[[nm]]
    if(is.null(w) || any(!is.finite(w))) next
    # ensure candidate is a proper simplex (non-negative, sums ~1)
    # optimization methods already ensure simplex constraint
    w <- as.numeric(w)
    if(any(!is.finite(w))) next
    if(any(w < -1e-8)) next
    if(abs(sum(w) - 1) > 1e-6) next

    s <- tryCatch(sse_calc(w), error = function(e) stop(sprintf("sse_calc failed while scoring candidate '%s': %s", nm, e$message)))

    if(is.finite(s) && s < best_score) {
      best_score <- s; best_name <- nm; best_coef <- w
    }
  }

  if(is.null(best_coef)) return(NULL)
  # optimization methods already ensure simplex constraint - no additional processing needed
  coef_final <- best_coef
  
  names(coef_final) <- veg_kept
  return(list(X = coef_final, fmat = NULL))
}

# Function to generate predictions for all years using fitted coefficients
generate_all_year_predictions <- function(coef, g_proj, lib_factor, veg_kept, df_location, avail_idx, veg_type) {
  # Get all unique years for this location
  all_years <- sort(unique(lubridate::year(df_location$date)))
  
  predictions <- list()
  
  for(yr in all_years) {
    # Subset data for this year
    df_year <- df_location[lubridate::year(df_location$date) == yr, , drop = FALSE]
    if(nrow(df_year) == 0) next
    
    # Prepare factor data for this year
    date_list <- prepare_factor_data(df_year, g_proj, avail_idx, veg_type)
    if(length(date_list) == 0) next
    
    # Generate predictions using fitted coefficients
    pred_matrix <- matrix(0, nrow = length(date_list), ncol = length(veg_kept))
    for(i in seq_along(date_list)) {
      doy <- date_list[[i]]$doy
      for(j in seq_along(veg_kept)) {
        vkey <- veg_kept[j]
        proto <- lib_factor[[vkey]]
        if(!is.null(proto) && !is.null(proto$T) && doy > 0 && doy <= nrow(proto$T) && ncol(proto$T) >= length(date_list[[i]]$z)) {
          pred_matrix[i, j] <- sum(proto$T[doy, seq_along(date_list[[i]]$z)] * date_list[[i]]$z)
        }
      }
    }
    
    # Apply fitted coefficients to get final prediction
    final_pred <- pred_matrix %*% coef
    
    predictions[[as.character(yr)]] <- data.frame(
      date = sapply(date_list, function(x) df_year$date[df_year$doy == x$doy][1]),
      doy = sapply(date_list, `[[`, "doy"),
      prediction = as.numeric(final_pred),
      year = yr,
      stringsAsFactors = FALSE
    )
  }
  
  return(predictions)
}

glsbb <- function(B, y, X, w, groups = NULL, progressor = NULL) {
  if (!requireNamespace("MASS", quietly = TRUE)) stop("MASS package is required for ginv")

  # Fit the initial GLS model to get residuals (weights already provided per-observation)
  w_sqrt <- sqrt(w)
  y_w <- y * w_sqrt
  X_w <- X * w_sqrt

  # Use pseudo-inverse for robustness
  XwX_inv <- MASS::ginv(crossprod(X_w))
  beta_hat <- XwX_inv %*% crossprod(X_w, y_w)

  # Calculate residuals on the original scale and ensure numeric vector
  residuals <- as.numeric(y - X %*% beta_hat)
  residuals[!is.finite(residuals)] <- 0

  n <- length(residuals)

  # Helper: estimate block size from autocorrelation structure of residuals
  estimate_blocksize_autocorr <- function(resid) {
    nres <- length(resid)
    if(nres < 4) return(1)
    max_lag <- min(50, floor(nres / 4))
    if(max_lag < 1) return(1)
    acf_obj <- tryCatch({ stats::acf(resid, plot = FALSE, lag.max = max_lag, na.action = na.pass) }, error = function(e) stop(sprintf("estimate_blocksize_autocorr: acf computation failed: %s", e$message)))
    if(is.null(acf_obj) || is.null(acf_obj$acf)) stop("estimate_blocksize_autocorr: acf result missing 'acf' component")
    # acf includes lag-0 at index 1
    acf_vals <- as.numeric(acf_obj$acf)
    acf_lags <- acf_vals[-1]
    # Use the sum of absolute autocorrelations as a simple proxy for dependence
    s <- sum(abs(acf_lags), na.rm = TRUE)
    # Heuristic block length: scale sum by factor, clamp to reasonable bounds
    l1 <- ceiling(2 * s)
    l2 <- ceiling(nres^(1/3))
    blk <- max(1, min(l1, l2, floor(nres / 2)))
    return(as.integer(blk))
  }

  # Decide whether to use moving-block bootstrap (MBB) or group/block resampling
  use_mbb <- TRUE
  if(!is.null(groups)) {
    # If groups are provided but represent true grouping (few unique groups << n), prefer group-based resampling
    u <- unique(groups)
    if(length(u) < max(3, n / 10)) {
      use_mbb <- FALSE
    }
  }

  boot_list <- vector("list", B)

  if(use_mbb) {
    blocksize <- estimate_blocksize_autocorr(residuals)
    starts <- seq_len(max(1, n - blocksize + 1))

    for(bi in seq_len(B)) {
      if(blocksize <= 1 || length(starts) <= 1) {
        # IID residual bootstrap fallback
        boot_residuals <- sample(residuals, n, replace = TRUE)
      } else {
        # Sample starting indices for blocks, concatenate, truncate to length n
        nblocks <- ceiling(n / blocksize)
        sampled_starts <- sample(starts, nblocks, replace = TRUE)
        idxs <- unlist(lapply(sampled_starts, function(s) seq.int(s, length.out = blocksize)))
        idxs <- idxs[seq_len(n)]
        boot_residuals <- residuals[idxs]
      }

      # Create bootstrap response variable
      y_boot <- as.numeric(X %*% beta_hat) + boot_residuals
      y_boot_w <- y_boot * w_sqrt

      # Re-fit the model with the bootstrap response
      boot_beta_hat <- XwX_inv %*% crossprod(X_w, y_boot_w)
      boot_list[[bi]] <- boot_beta_hat

      # Progress callback
      try({ if(!is.null(progressor) && is.function(progressor)) progressor() }, silent = TRUE)
    }
  } else {
    # Group-based block resampling (existing behaviour)
    res_grouped <- split(residuals, groups)
    group_names <- names(res_grouped)
    for(bi in seq_len(B)) {
      boot_groups <- sample(group_names, length(group_names), replace = TRUE)
      boot_res_list <- res_grouped[boot_groups]
      boot_residuals <- unsplit(boot_res_list, rep(names(boot_res_list), lengths(boot_res_list)))

      y_boot <- as.numeric(X %*% beta_hat) + boot_residuals
      y_boot_w <- y_boot * w_sqrt
      boot_beta_hat <- XwX_inv %*% crossprod(X_w, y_boot_w)
      boot_list[[bi]] <- boot_beta_hat

      try({ if(!is.null(progressor) && is.function(progressor)) progressor() }, silent = TRUE)
    }
  }

  boot_betas <- do.call(cbind, boot_list)
  return(t(boot_betas))
}

# Helper to perform bootstrap for IVW results
perform_ivw_bootstrap <- function(y_all, g_proj, date_list, lib_factor, veg_kept, coef, B, loc, yr, progress_every) {
  if (B <= 1) return(NULL)

  # Construct the design matrix for the bootstrap
  B_all <- do.call(rbind, lapply(date_list, function(dl) {
    sapply(veg_kept, function(vv) {
      proto <- lib_factor[[vv]]
      if (dl$doy > 0 && dl$doy <= nrow(proto$T)) {
        sum(proto$T[dl$doy, ])
      } else {
        0
      }
    })
  }))

  # Adaptive bootstrap size based on available data
  n_obs <- nrow(B_all)
  if (n_obs < MIN_OBS_FOR_BOOT) return(NULL)

  # Use smaller bootstrap size for sparse data, but at least 10 reps
  adaptive_B <- max(10, min(B, n_obs * 2))
  if (adaptive_B != B && FALSE) {  # Disabled verbose note
    cat(sprintf("NOTE: Reduced bootstrap size from %d to %d due to limited observations (%d)\n", B, adaptive_B, n_obs))
  }

  # Define groups for block bootstrap (e.g., by day of year)
  groups <- sapply(date_list, `[[`, "doy")

  # Get weights from KDE densities (normalized), no inversion needed
  w_all <- sapply(date_list, function(dl) {
    vars <- sapply(veg_kept, function(v) {
      if(!is.null(lib_factor[[v]]$T_var) && length(lib_factor[[v]]$T_var) >= dl$doy && dl$doy >= 1) {
        # Use KDE density directly as weight
        density_val <- lib_factor[[v]]$T_var[dl$doy]
        if(is.finite(density_val) && density_val > 0) density_val else NA
      } else {
        NA
      }
    })
    valid <- is.finite(vars)
    if(any(valid)) {
      mean(vars[valid], na.rm = TRUE)
    } else {
      NA
    }
  })

  # Apply the same weight transformation as in ivw_solve_mixture
  w_all <- ifelse(is.finite(w_all) & w_all > 0, w_all, 1)
  w_all <- sqrt(w_all)  # Soften extreme KDE densities
  w_all <- w_all / mean(w_all, na.rm = TRUE)  # Renormalize

  # Use progressr properly: create and use the progressor inside with_progress
  boot_coefs <- NULL
  if(requireNamespace("progressr", quietly = TRUE) && is.numeric(progress_every) && progress_every > 0 && adaptive_B >= progress_every) {
    try({
      boot_coefs <- progressr::with_progress({
        p_local <- progressr::progressor(steps = adaptive_B, message = sprintf("Bootstrap %s_%s", loc, yr))
        # pass the progressor object directly to glsbb; calls to p_local() will be listened to
        glsbb(adaptive_B, y_all, B_all, w_all, groups, progressor = p_local)
      })
    }, silent = TRUE)
  } else {
    # no progressr or not requested -> regular call
    boot_coefs <- glsbb(adaptive_B, y_all, B_all, w_all, groups, progressor = NULL)
  }

  boot_rows <- as.data.frame(boot_coefs)
  names(boot_rows) <- veg_kept
  boot_rows$location_id <- loc
  boot_rows$year <- yr
  boot_rows$boot_rep <- 1:adaptive_B

  return(boot_rows)
}

# Declare commonly used global-like variables to silence static-analysis "no visible binding" warnings
  utils::globalVariables(c(
  "PROGRESS_LOG_TO_FILE", "PARALLEL_ENABLE", "PROGRESS_BAR", "PARALLEL_WORKERS",
  "BOOTSTRAP_B", "BOOT_MIN_REPS_PER_VEG", "MIN_INDEX_SD", "GLOBAL_PCA", "lib_factor_pca",
  "avail", "OUT_DIR", "LOWER_BND", "EPS_SIGMA", "MIN_OBS_FOR_BOOT", "MAX_VEG_COMPONENTS", "veg_counts",
compute_dynamic_ridge
))
  # Duplicate processing block removed (it was accidentally inserted and referenced an undefined 'v').)
  n_loc_years <- nrow(loc_years)

## Clustering data structures removed; worker code will use lib_factor_pca directly.

# Use TEST years for processing: tasks should cover test-year location-year pairs
cat("Constructing task list from TEST dataset...\n")

cat("DEBUG: Checking TEST_YEARS and data...\n")
if(exists("TEST_YEARS")) {
  cat(sprintf("TEST_YEARS exists: %s\n", paste(TEST_YEARS, collapse = ", ")))
} else {
  cat("TEST_YEARS does not exist\n")
}

# Check original data size
cat(sprintf("Original df_full has %d rows\n", nrow(df_full)))
cat(sprintf("Unique locations in df_full: %d\n", length(unique(df_full$location_id))))
cat(sprintf("Unique years in df_full: %s\n", paste(sort(unique(df_full$year)), collapse = ", ")))

# df_full was saved earlier (original full dataset). Use it to derive test tasks.
df_tasks <- if(!is.null(TEST_YEARS)) {
  filtered <- df_full[df_full$year %in% TEST_YEARS, , drop = FALSE]
  cat(sprintf("After TEST_YEARS filtering: %d rows\n", nrow(filtered)))
  filtered
} else {
  cat("No TEST_YEARS filtering applied\n")
  df_full
}

# Ensure Veg mapping from gpts_map is available for test data (join if present)
if(exists("gpts_map")) {
  cat(sprintf("gpts_map exists with %d rows\n", nrow(gpts_map)))
  df_tasks <- dplyr::left_join(df_tasks, gpts_map, by = "location_id")
} else {
  cat("gpts_map does not exist\n")
}

# Canonicalize minimal fields used downstream
if("Veg" %in% names(df_tasks)) df_tasks$Veg <- tolower(as.character(df_tasks$Veg))
if("date" %in% names(df_tasks)) {
  df_tasks$date <- as.Date(df_tasks$date)
  df_tasks$year <- as.integer(lubridate::year(df_tasks$date))
  df_tasks$doy <- lubridate::yday(df_tasks$date)
  df_tasks$doy[df_tasks$doy < 1 | df_tasks$doy > 366] <- NA_integer_
}

test_loc_years <- df_tasks %>% dplyr::filter(!is.na(location_id) & location_id != "" & !is.na(year) & year > 0) %>% dplyr::distinct(location_id, year)
cat(sprintf("Final test_loc_years: %d rows from %d unique locations\n", nrow(test_loc_years), length(unique(test_loc_years$location_id))))
if(nrow(test_loc_years) > 0) {
  print(head(test_loc_years, min(10, nrow(test_loc_years))))
}

n_loc_years <- nrow(test_loc_years)
task_list <- lapply(seq_len(n_loc_years), function(i) {
  list(loc = test_loc_years$location_id[i], yr = test_loc_years$year[i])
})

cat(sprintf("Final task_list length: %d\n", length(task_list)))

# Fail hard if there are no test tasks to process
if(length(task_list) == 0) {
  stop("CRITICAL: No testing tasks found. Ensure TEST_YEARS and input data produce at least one location-year pair for testing.")
}


  # Ensure loc_years looks sane
  if(nrow(loc_years) > 0) {
    cat(sprintf("loc_years sample: %d rows, first: %s_%s\n", nrow(loc_years), as.character(loc_years$location_id[1]), as.character(loc_years$year[1])))
    flush.console()
  }

  # Validate loc_years: remove any rows missing location_id or year
  bad_idx <- which(is.na(loc_years$location_id) | !is.finite(loc_years$year))
  if(length(bad_idx) > 0) {
    cat(sprintf("loc_years: dropping %d invalid rows (missing location_id or year)\n", length(bad_idx)))
    loc_years <- loc_years[-bad_idx, , drop = FALSE]
  }

  # Pre-bind required globals into a dedicated environment so the closure has reliable access
  required_globals <- list(
    df = df, lib = lib, GLOBAL_PCA = GLOBAL_PCA, lib_factor_pca = lib_factor_pca,
    GLOBAL_LDA = if(exists("GLOBAL_LDA")) GLOBAL_LDA else NULL,
    lib_factor_lda = if(exists("lib_factor_lda")) lib_factor_lda else list(),
    avail = avail, LOWER_BND = LOWER_BND, EPS_SIGMA = EPS_SIGMA,
    MIN_OBS_FOR_BOOT = MIN_OBS_FOR_BOOT, BOOTSTRAP_B = BOOTSTRAP_B,
    BOOT_MIN_REPS_PER_VEG = BOOT_MIN_REPS_PER_VEG,
    MAX_VEG_COMPONENTS = MAX_VEG_COMPONENTS, veg_counts = veg_counts, OUT_DIR = OUT_DIR,
    FACTOR_MODE = FACTOR_MODE, PROGRESS_EVERY_TASK = PROGRESS_EVERY_TASK,
    df_tasks = df_tasks,
    # Functions
  prepare_factor_data = prepare_factor_data,
  ivw_solve_mixture = ivw_solve_mixture,
  ivw_solve_mixture_gd = ivw_solve_mixture_gd,
  solve_mixture_gd = solve_mixture_gd,
  generate_all_year_predictions = generate_all_year_predictions,
  perform_ivw_bootstrap = perform_ivw_bootstrap,
  glsbb = glsbb,
  compute_dynamic_ridge = compute_dynamic_ridge
  )
  # Helper: compute weighted fit score for a set of coefficients on a date_list
  compute_fit_score <- function(date_list, lib_factor, veg_kept, coef, g_proj) {
    if(is.null(coef) || length(coef) == 0) return(Inf)
    # Build F_matrix as in ivw solvers
    rnk <- if (!is.null(g_proj$rank)) g_proj$rank else if(!is.null(g_proj$fit$scaling)) ncol(g_proj$fit$scaling) else 0
    if(rnk == 0) return(Inf)
    F_matrix <- matrix(0, nrow = length(date_list), ncol = length(veg_kept))
    for(i in seq_along(date_list)) {
      doy <- date_list[[i]]$doy
      for(j in seq_along(veg_kept)) {
        vkey <- veg_kept[j]
        proto <- lib_factor[[vkey]]
        if(!is.null(proto) && !is.null(proto$T) && doy > 0 && doy <= nrow(proto$T) && ncol(proto$T) >= rnk) {
          F_matrix[i, j] <- sum(proto$T[doy, seq_len(rnk)] * date_list[[i]]$z[seq_len(rnk)])
        }
      }
    }

    # compute KDE densities as weights (no inversion needed)
    date_var <- rep(NA_real_, nrow(F_matrix))
    for(i in seq_len(nrow(F_matrix))) {
      vars <- numeric(0)
      for(j in seq_along(veg_kept)) {
        vkey <- veg_kept[j]
        proto <- lib_factor[[vkey]]
        if(!is.null(proto) && !is.null(proto$T_var) && length(proto$T_var) >= date_list[[i]]$doy) {
          d <- date_list[[i]]$doy
          if(d >=1 && d <= length(proto$T_var) && is.finite(proto$T_var[d])) {
            # Use KDE density directly as weight
            density_val <- proto$T_var[d]
            vars <- c(vars, density_val)
          }
        }
      }
      if(length(vars) > 0) date_var[i] <- mean(vars, na.rm = TRUE)
    }
    date_var[!is.finite(date_var) | date_var <= 0] <- NA_real_
    if(all(is.na(date_var))) {
      w_all <- rep(1, nrow(F_matrix))
    } else {
      # Use KDE densities directly as weights
      w_all <- date_var
      w_all[!is.finite(w_all) | w_all <= 0] <- 1
      # Apply the same transformation as in ivw_solve_mixture
      w_all <- sqrt(w_all)  # Soften extreme KDE densities
    }

    pred <- as.numeric(F_matrix %*% coef)
    residuals <- 1 - pred
    # weighted SSE (ignore NA rows)
    mask <- is.finite(residuals) & is.finite(w_all)
    if(!any(mask)) return(Inf)
    sse <- sum((residuals[mask]^2) * w_all[mask], na.rm = TRUE)
    # normalize by sum weights so scores comparable across date counts
    sse / sum(w_all[mask])
  }
  env_task <- list2env(required_globals, parent = globalenv())
  # Expose helper functions and small utilities into env_task for safe worker evaluation
  # (prepares functions used by fit_one_task that may be defined later in the file)
  # This block is intentionally left empty for now. Helpers are defined later
  
  fit_one_task <- function(task) {
    # Unpack task variables
    loc <- task$loc
    yr <- task$yr
    
    # Helper to log why a task returned NULL so we can debug failed fits.
    dbg_return_null <- function(reason) {
      reason_str <- as.character(reason)
      if (!grepl("^(ERROR:|error:)", reason_str, ignore.case = TRUE)) {
        reason_str <- paste0("ERROR:", reason_str)
      }

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
        if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
        outf <<- file.path(OUT_DIR, "fit_fail_reasons.csv")
        if(!file.exists(outf)) write.table(row, outf, sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
        else write.table(row, outf, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        success <<- file.exists(outf) && file.info(outf)$size > 0
      }, silent = TRUE)

      if (!isTRUE(success)) {
        msg <- sprintf("[dbg_return_null] Failed to write failure row for loc=%s year=%s reason=%s", as.character(loc), as.character(yr), reason_str)
        try(cat(msg, "\n"), silent = TRUE)
        try(log_msg(msg), silent = TRUE)
      }

      return(NULL)
    }
        
    res_safe <- tryCatch({
  # Subset from testing dataset - include ALL test years for this location for full predictions
  dly <- df_tasks[df_tasks$location_id == loc, , drop = FALSE]
  # Create a year-specific subset for fitting so each year is fitted independently (use test years)
  dly_year <- dly[lubridate::year(dly$date) == yr, , drop = FALSE]
  if(nrow(dly_year) == 0) return(dbg_return_null("no_rows_year"))
      if(nrow(dly) == 0) return(dbg_return_null("no_rows"))

  # Calculate Q10 and Q90 DVI for this location-year (testing year)
  dly_train <- df_tasks[df_tasks$location_id == loc & lubridate::year(df_tasks$date) == yr, , drop = FALSE]
      dvi_vals <- dly_train$DVI[is.finite(dly_train$DVI)]
      q10_dvi <- if(length(dvi_vals) > 0) {
        stats::quantile(dvi_vals, 0.10, na.rm = TRUE)
      } else {
        NA
      }
      q90_dvi <- if(length(dvi_vals) > 0) {
        stats::quantile(dvi_vals, 0.90, na.rm = TRUE)
      } else {
        NA
      }

      if(!FACTOR_MODE) return(dbg_return_null("not_factor_mode"))

      veg_names <- names(lib)
      veg_task <- unique(na.omit(dly$Veg))
      veg_kept <- intersect(ALLOWED_VEG, veg_names)
      veg_kept <- veg_kept[tolower(veg_kept) != "barren"]

      if(length(veg_kept) < MAX_VEG_COMPONENTS && length(veg_counts) > 0) {
        global_order <- names(veg_counts)
        global_order <- intersect(global_order, veg_names)
        add <- setdiff(global_order, veg_kept)
        add <- add[tolower(add) != "barren"]
        if(length(add) > 0) veg_kept <- c(veg_kept, add[seq_len(min(length(add), MAX_VEG_COMPONENTS - length(veg_kept)))])
      }
      if(length(veg_kept) > MAX_VEG_COMPONENTS) veg_kept <- veg_kept[seq_len(MAX_VEG_COMPONENTS)]
      if(length(veg_kept) == 0) return(dbg_return_null("no_veg_kept"))

      # --- PCA Path with Gradient Descent ---
  # Prepare factor data using only the target year so fits are independent per year
  date_list_pca <- prepare_factor_data(dly_year, GLOBAL_PCA, avail, veg_task[1])
      if(length(date_list_pca) == 0) return(dbg_return_null("no_date_list_pca"))
      
      y_all_pca <- unlist(lapply(date_list_pca, function(x) x$z))
      if(length(y_all_pca) == 0) return(dbg_return_null("no_y_all_pca"))

      # Try both gradient-descent IVW and the constrained solver; pick best by weighted SSE
      sol_gd <- ivw_solve_mixture_gd(y_all_pca, GLOBAL_PCA, date_list_pca, lib_factor_pca, veg_kept)
      sol_alt <- ivw_solve_mixture(y_all_pca, GLOBAL_PCA, date_list_pca, lib_factor_pca, veg_kept)

      # Extract candidate coefficients and compute scores
      cand <- list()
      if(!is.null(sol_gd$X)) cand[['PCA_GD']] <- list(coef = sol_gd$X, sol = sol_gd)
      if(!is.null(sol_alt$X)) cand[['PCA_ALT']] <- list(coef = sol_alt$X, sol = sol_alt)

      if(length(cand) == 0) {
        stop("Both PCA solvers failed to produce valid coefficients for ", loc, "_", yr)
      }
      
      scores <- sapply(names(cand), function(nm) {
        compute_fit_score(date_list_pca, lib_factor_pca, veg_kept, cand[[nm]]$coef, GLOBAL_PCA)
      })
      best_nm <- names(which.min(scores))
      best <- cand[[best_nm]]
      coef_pca <- best$coef
  # Ensure names (solver returns simplex weights already)
  if(is.null(names(coef_pca))) names(coef_pca) <- veg_kept
        coef_df_pca <- data.frame(location_id = loc, year = yr, Veg = names(coef_pca), coef = as.numeric(coef_pca), FitMethod = best_nm, stringsAsFactors = FALSE)

  # Bootstrap using chosen solver - BOOTSTRAP_B must be properly configured
  if(!is.null(coef_pca) && length(coef_pca) > 0) {
    if(!is.numeric(BOOTSTRAP_B) || BOOTSTRAP_B <= 1) {
      stop("BOOTSTRAP_B must be configured as a positive integer > 1, got: ", as.character(BOOTSTRAP_B))
    }
    B_used_pca <- BOOTSTRAP_B
    boot_rows_pca <- perform_ivw_bootstrap(y_all_pca, GLOBAL_PCA, date_list_pca, lib_factor_pca, veg_kept, coef_pca, B_used_pca, loc, yr, PROGRESS_EVERY_TASK)
        if(is.null(boot_rows_pca) || nrow(boot_rows_pca) == 0) {
          boot_rows_pca <- as.data.frame(matrix(NA_real_, nrow = 1, ncol = length(veg_kept)))
          names(boot_rows_pca) <- veg_kept
          boot_rows_pca$location_id <- loc; boot_rows_pca$year <- yr; boot_rows_pca$boot_rep <- NA_integer_; boot_rows_pca$FitMethod <- best_nm
        } else {
          boot_rows_pca$FitMethod <- best_nm
        }
    # Diagnostic log: show bootstrap counts per veg for PCA
    try({
      if(!is.null(boot_rows_pca) && nrow(boot_rows_pca) > 0) {
        n_per_veg <- sapply(veg_kept, function(v) sum(!is.na(boot_rows_pca[[v]])))
        cat(sprintf("BOOT_DIAG PCA %s_%s: n_per_veg=%s\n", loc, yr, paste(names(n_per_veg), n_per_veg, sep='=', collapse=',')))
        flush.console()
      }
    }, silent = TRUE)
  }

  # --- LDA Path with Gradient Descent ---
  coef_df_lda <- NULL
  boot_rows_lda <- NULL
  coef_lda <- NULL
  if(!is.null(GLOBAL_LDA) && is.list(GLOBAL_LDA) && !is.null(GLOBAL_LDA$fit)) {
        # Use the same year-specific subset for LDA fitting
        date_list_lda <- prepare_factor_data(dly_year, GLOBAL_LDA, avail, veg_task[1])
        if(length(date_list_lda) > 0) {
          y_all_lda <- unlist(lapply(date_list_lda, function(x) x$z))
          if(length(y_all_lda) > 0) {
            # Basic LDA preprocessing diagnostics
            try({
              if(FALSE) {  # Disabled debug output
                cat(sprintf("DBG_LDA: date_list_lda entries=%d, y_all_lda length=%d\n", length(date_list_lda), length(y_all_lda)))
                flush.console()
              }
            }, silent = TRUE)

            sol_gd_lda <- ivw_solve_mixture_gd(y_all_lda, GLOBAL_LDA, date_list_lda, lib_factor_lda, veg_kept)
            sol_alt_lda <- ivw_solve_mixture(y_all_lda, GLOBAL_LDA, date_list_lda, lib_factor_lda, veg_kept)
            cand_lda <- list()
            if(!is.null(sol_gd_lda$X)) cand_lda[['LDA_GD']] <- list(coef = sol_gd_lda$X, sol = sol_gd_lda)
            if(!is.null(sol_alt_lda$X)) cand_lda[['LDA_ALT']] <- list(coef = sol_alt_lda$X, sol = sol_alt_lda)

            if(length(cand_lda) == 0) {
              stop("Both LDA solvers failed to produce valid coefficients for ", loc, "_", yr)
              # Build detailed diagnostics to help trace the failure
              diag_parts <- list()
              has_glda <- !is.null(GLOBAL_LDA) && is.list(GLOBAL_LDA) && !is.null(GLOBAL_LDA$fit)
              diag_parts[['GLOBAL_LDA_present']] <- as.character(has_glda)
              if(has_glda) {
                try({
                  diag_parts[['LDA_scaling_dim']] <- paste(dim(GLOBAL_LDA$fit$scaling), collapse = "x")
                }, silent = TRUE)
              }
              diag_parts[['date_list_lda_len']] <- as.character(length(date_list_lda))
              diag_parts[['y_all_lda_len']] <- as.character(length(y_all_lda))
              # summary of y_all_lda
              try({
                diag_parts[['y_all_lda_summary']] <- paste0("min=", signif(min(y_all_lda, na.rm=TRUE),6), ", mean=", signif(mean(y_all_lda, na.rm=TRUE),6), ", max=", signif(max(y_all_lda, na.rm=TRUE),6))
              }, silent = TRUE)
              # show first few z values from date_list_lda
              try({
                zvals <- unlist(lapply(date_list_lda[seq_len(min(5, length(date_list_lda)))], function(x) x$z))
                diag_parts[['date_list_lda_sample_z']] <- paste(signif(zvals,6), collapse = ",")
              }, silent = TRUE)
              # prototype availability per veg
              try({
                proto_summary <- sapply(veg_kept, function(v) {
                  if(!exists("lib_factor_lda") || is.null(lib_factor_lda[[v]])) return(NA_integer_)
                  ns <- lib_factor_lda[[v]]$n_samples
                  ns <- if(is.null(ns)) NA_integer_ else ns
                  # check if T is finite and non-zero
                  ok <- FALSE
                  Tmat <- lib_factor_lda[[v]]$T
                  if(!is.null(Tmat) && is.matrix(Tmat)) {
                    ok <- any(is.finite(Tmat)) && any(abs(Tmat) > 0)
                  }
                  paste0("n_samples=", ns, ",has_valid_T=", ok)
                })
                diag_parts[['proto_per_veg']] <- paste(names(proto_summary), proto_summary, sep='=', collapse=';')
              }, silent = TRUE)

              # Log diagnostics instead of stopping
              diag_msg <- paste(unlist(diag_parts), collapse = "; ")
              cat(sprintf("WARNING: LDA fit failed (no candidate solutions) for %s_%s. Diagnostics: %s\n", loc, yr, diag_msg))
              flush.console()

              # Create NA coefficients instead of failing
              coef_lda <- rep(NA_real_, length(veg_kept)); names(coef_lda) <- veg_kept
              coef_df_lda <- data.frame(location_id = loc, year = yr, Veg = names(coef_lda), coef = as.numeric(coef_lda), FitMethod = "LDA_GD", stringsAsFactors = FALSE)
              boot_rows_lda <- as.data.frame(matrix(NA_real_, nrow = 1, ncol = length(veg_kept)))
              names(boot_rows_lda) <- veg_kept
              boot_rows_lda$location_id <- loc; boot_rows_lda$year <- yr; boot_rows_lda$boot_rep <- NA_integer_; boot_rows_lda$FitMethod <- "LDA_GD"
            } else {
              scores_lda <- sapply(names(cand_lda), function(nm) {
                compute_fit_score(date_list_lda, lib_factor_lda, veg_kept, cand_lda[[nm]]$coef, GLOBAL_LDA)
              })
              best_nm_lda <- names(which.min(scores_lda))
              best_lda <- cand_lda[[best_nm_lda]]
              coef_lda <- best_lda$coef
              if(is.null(names(coef_lda))) names(coef_lda) <- veg_kept
              coef_df_lda <- data.frame(location_id = loc, year = yr, Veg = names(coef_lda), coef = as.numeric(coef_lda), FitMethod = best_nm_lda, stringsAsFactors = FALSE)

              # Use the same bootstrap configuration as PCA (already validated)
              B_used_lda <- BOOTSTRAP_B
              boot_rows_lda <- perform_ivw_bootstrap(y_all_lda, GLOBAL_LDA, date_list_lda, lib_factor_lda, veg_kept, coef_lda, B_used_lda, loc, yr, PROGRESS_EVERY_TASK)
              if(is.null(boot_rows_lda) || nrow(boot_rows_lda) == 0) {
                boot_rows_lda <- as.data.frame(matrix(NA_real_, nrow = 1, ncol = length(veg_kept)))
                names(boot_rows_lda) <- veg_kept
                boot_rows_lda$location_id <- loc; boot_rows_lda$year <- yr; boot_rows_lda$boot_rep <- NA_integer_; boot_rows_lda$FitMethod <- best_nm_lda
              } else {
                boot_rows_lda$FitMethod <- best_nm_lda
              }
              # Diagnostic log: show bootstrap counts per veg for LDA
              try({
                if(!is.null(boot_rows_lda) && nrow(boot_rows_lda) > 0) {
                  n_per_veg_lda <- sapply(veg_kept, function(v) sum(!is.na(boot_rows_lda[[v]])))
                  cat(sprintf("BOOT_DIAG LDA %s_%s: n_per_veg=%s\n", loc, yr, paste(names(n_per_veg_lda), n_per_veg_lda, sep='=', collapse=',')))
                  flush.console()
                }
              }, silent = TRUE)
            }
          }
        }
      }

      # Generate predictions for ALL years using fitted coefficients
      predictions_pca <- generate_all_year_predictions(coef_pca, GLOBAL_PCA, lib_factor_pca, veg_kept, dly, avail, veg_task[1])
      predictions_lda <- if(!is.null(coef_lda)) {
        generate_all_year_predictions(coef_lda, GLOBAL_LDA, lib_factor_lda, veg_kept, dly, avail, veg_task[1])
      } else {
        NULL
      }

      # Observed-based 20-day rolling median maxima (single value per location-year)
      compute_max20_observed <- function(dly_df) {
        # Strict observed-based: compute 20-day centered rolling median; if no finite
        # rolling-median values exist, return NA (do NOT fallback to raw max).
        if(is.null(dly_df) || nrow(dly_df) == 0) return(NA_real_)
        if(!"DVI" %in% names(dly_df)) return(NA_real_)
        series <- as.numeric(dly_df$DVI)
        if(all(!is.finite(series))) return(NA_real_)
        ma_obs <- NULL
        try({
          ma_obs <- zoo::rollmedian(series, k = 20, fill = NA, align = "center")
        }, silent = TRUE)
        if(!is.null(ma_obs) && any(is.finite(ma_obs))) {
          return(max(ma_obs[is.finite(ma_obs)], na.rm = TRUE))
        }
        # If rollmedian produced no finite values, return NA (no fallback)
        return(NA_real_)
      }

      # PCA/LDA-derived max20 helper removed; use observed-only
      max20_observed <- compute_max20_observed(dly_year)

      # Combine results (both PCA and LDA must succeed)
      if(is.null(coef_df_pca)) {
        stop("PCA coefficient computation failed for ", loc, "_", yr)
      }
      if(is.null(coef_df_lda)) {
        stop("LDA coefficient computation failed for ", loc, "_", yr)
      }
      coef_df <- rbind(coef_df_pca, coef_df_lda)
      boot_rows <- rbind(boot_rows_pca, boot_rows_lda)
      
      # Return a list of data frames including Q10 and Q90 DVI and all-year predictions
  return(list(coef_df = coef_df, boot_rows = boot_rows, q10_dvi = q10_dvi, q90_dvi = q90_dvi,
    predictions_pca = predictions_pca, predictions_lda = predictions_lda,
    max20_observed = max20_observed))

    }, error = function(e) {
      dbg_return_null(paste0('error:', as.character(e$message)))
    })
    
    # FAIL HARD if LDA produced no coefficients: provide diagnostics and stop
    try({
      if(!is.null(res_safe) && is.data.frame(res_safe$coef_df)) {
        lda_rows <- res_safe$coef_df[grepl("^LDA", res_safe$coef_df$FitMethod), , drop = FALSE]
        if(nrow(lda_rows) > 0 && all(is.na(lda_rows$coef))) {
          # Build diagnostics
          has_glda <- !is.null(GLOBAL_LDA) && is.list(GLOBAL_LDA) && !is.null(GLOBAL_LDA$fit)
          # Count test-year rows for this location-year
          n_test_rows_locyr <- 0
          try({
            n_test_rows_locyr <- nrow(df_test[df_test$location_id == loc & lubridate::year(df_test$date) == yr, , drop = FALSE])
          }, silent = TRUE)
          # Sum of prototype sample counts available in lib_factor_lda
          total_proto_samples <- NA_integer_
          try({
            if(exists("lib_factor_lda") && is.list(lib_factor_lda)) {
              total_proto_samples <- sum(unlist(lapply(lib_factor_lda, function(x) if(!is.null(x$n_samples)) x$n_samples else 0)), na.rm = TRUE)
            }
          }, silent = TRUE)
          msg <- paste0("LDA fit failure at ", loc, "_", yr, ": no LDA coefficients produced. ",
                        "GLOBAL_LDA_present=", has_glda, ", test_rows_locyr=", n_test_rows_locyr, ", total_proto_samples=", total_proto_samples,
                        ". Check that TEST_YEARS contain observations for this loc-year and that GLOBAL_LDA and lib_factor_lda were built correctly.")
          cat(paste0("WARNING: ", msg, "\n")); flush.console()
        }
      }
    }, silent = TRUE)

    cat(sprintf("Completed fit for %s_%s\n", loc, yr))  # Disabled - too verbose
    return(res_safe)
  }
  # Ensure the task function can see the prebound globals
  environment(fit_one_task) <- env_task

  # Main processing loop
cat("Starting main processing loop...\n")
start_time <- Sys.time()
results_list <- .run_map(task_list, fit_one_task)
end_time <- Sys.time()
processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat(sprintf("Main processing loop finished in %.2f seconds (%.2f minutes)\n", processing_time, processing_time/60))
cat(sprintf("Average time per task: %.2f seconds\n", processing_time / length(task_list)))

# Process results and write to Excel
cat("Processing results and writing to Excel files...\n")
  
# Filter out null results
results_list <- results_list[!sapply(results_list, is.null)]
cat(sprintf("After filtering NULL results: %d results remaining\n", length(results_list)))

# Debug: Analyze what happened to the tasks
if(exists("task_list")) {
  null_count <- length(task_list) - length(results_list)
  cat(sprintf("Task summary: %d successful, %d failed (NULL results)\n", length(results_list), null_count))
}

# Comprehensive diagnostics helper: RMSE, MAE, R2, AIC/BIC, leverage, Cook's D, VIF, condition number
compute_extended_diagnostics <- function(y, E, w = NULL, templates = NULL) {
  # y: response vector (n)
  # E: design matrix (n x p)
  # w: optional observation weights (length n) - currently unused except in residual weighting
  # templates: optional list of additional info (kept for future extension)
  if(is.null(E) || is.null(y)) stop("compute_extended_diagnostics: E and y must be provided and non-null")
  y <- as.numeric(y)
  E <- tryCatch(as.matrix(E), error = function(e) stop(sprintf("compute_extended_diagnostics: failed to coerce E to matrix: %s", e$message)))
  if(length(y) != nrow(E) || nrow(E) < 1) stop("compute_extended_diagnostics: invalid input shapes - length(y) must match nrow(E) and be >= 1")

  # Fit via least squares (unweighted) to compute diagnostics
  fit <- tryCatch({
    qr_res <- qr(E)
    coef_hat <- tryCatch(qr.solve(E, y), error = function(e) stop(sprintf("compute_extended_diagnostics: linear solve failed: %s", e$message)))
    pred <- as.numeric(E %*% coef_hat)
    resid <- y - pred
    list(coef = coef_hat, pred = pred, resid = resid)
  }, error = function(e) stop(sprintf("compute_extended_diagnostics: least-squares fit failed: %s", e$message)))

  resid <- fit$resid
  n <- length(y)
  p <- ncol(E)

  # Basic metrics
  rmse <- sqrt(mean(resid^2, na.rm = TRUE))
  mae <- mean(abs(resid), na.rm = TRUE)
  ss_tot <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  ss_res <- sum(resid^2, na.rm = TRUE)
  r2 <- ifelse(ss_tot <= 0, NA_real_, 1 - ss_res / ss_tot)

  # Information criteria (use residual variance)
  sigma2 <- ss_res / max(1, (n - p))
  aic <- n * log(ss_res / n) + 2 * p
  bic <- n * log(ss_res / n) + log(n) * p

  # Leverage and Cook's distance
  H <- tryCatch({
    E %*% solve(crossprod(E)) %*% t(E)
  }, error = function(e) stop(sprintf("compute_extended_diagnostics: failed to compute hat matrix H: %s", e$message)))
  leverage <- tryCatch({ if(is.matrix(H)) diag(H) else as.numeric(H) }, error = function(e) stop(sprintf("compute_extended_diagnostics: failed to extract leverage from H: %s", e$message)))
  leverage[!is.finite(leverage)] <- NA_real_
  cooks_d <- rep(NA_real_, n)
  if(all(is.finite(resid)) && all(is.finite(leverage))) {
    cooks_d <- (resid^2) * leverage / (p * (1 - leverage)^2)
  }

  # Multicollinearity: VIF from correlation matrix as a proxy
  vif <- tryCatch({
    cr <- cor(E, use = "pairwise.complete.obs")
    diag(solve(cr))
  }, error = function(e) stop(sprintf("compute_extended_diagnostics: VIF computation failed: %s", e$message)))

  cond_num <- tryCatch(kappa(E), error = function(e) stop(sprintf("compute_extended_diagnostics: condition number computation failed: %s", e$message)))

  list(
    rmse = as.numeric(rmse), mae = as.numeric(mae), r2 = as.numeric(r2),
    aic = as.numeric(aic), bic = as.numeric(bic),
    max_leverage = as.numeric(max(leverage, na.rm = TRUE)),
    max_cooks_d = as.numeric(max(cooks_d, na.rm = TRUE)),
    mean_vif = as.numeric(mean(vif, na.rm = TRUE)),
    condition_number = as.numeric(cond_num)
  )
}

if(length(results_list) == 0) {
  cat("ERROR: All tasks returned NULL results! Check fit_fail_reasons.csv for details.\n")
  cat("Most likely causes:\n")
  cat("1. TEST_YEARS not set or empty\n")
  cat("2. No data available for the specified test years\n")
  cat("3. Data filtering issues (missing location_id, year, or Veg)\n")
  cat("4. All locations have insufficient data for fitting\n")
  return(NULL)
}

if (length(results_list) > 0) {
    # Debug: Check structure of first few results
    cat("Checking result structures...\n")
    for(i in seq_len(min(3, length(results_list)))) {
      res <- results_list[[i]]
      cat(sprintf("Result %d: coef_df has %d rows, boot_rows is %s\n", 
                 i, 
                 if(!is.null(res$coef_df)) nrow(res$coef_df) else 0,
                 if(!is.null(res$boot_rows)) paste(nrow(res$boot_rows), "rows") else "NULL"))
    }
    
    # Combine all coefficient and bootstrap data frames with error handling
    cat("Combining coefficient data frames...\n")
    coef_list <- lapply(results_list, function(res) {
      if(is.null(res$coef_df) || nrow(res$coef_df) == 0) return(NULL)
      return(res$coef_df)
    })
    coef_list <- coef_list[!sapply(coef_list, is.null)]
    
    if(length(coef_list) == 0) {
      cat("ERROR: No valid coefficient data frames found!\n")
      return(NULL)
    }
    
    all_coefs <- tryCatch({
      do.call(rbind, coef_list)
    }, error = function(e) {
      cat(sprintf("ERROR combining coef_df: %s\n", e$message))
      cat("First coef_df structure:\n")
      print(str(coef_list[[1]]))
      if(length(coef_list) > 1) {
        cat("Second coef_df structure:\n")
        print(str(coef_list[[2]]))
      }
      return(NULL)
    })
    
    if(is.null(all_coefs)) {
      cat("Failed to combine coefficient data frames\n")
      return(NULL)
    }
    
    cat(sprintf("Combined coefficients: %d rows\n", nrow(all_coefs)))
    
    # Combine bootstrap data with error handling
    cat("Combining bootstrap data frames...\n")
    boot_list <- lapply(results_list, function(res) {
      if(is.null(res$boot_rows) || nrow(res$boot_rows) == 0) return(NULL)
      return(res$boot_rows)
    })
    boot_list <- boot_list[!sapply(boot_list, is.null)]
    
    all_boots <- if(length(boot_list) > 0) {
      tryCatch({
        do.call(rbind, boot_list)
      }, error = function(e) {
        cat(sprintf("ERROR combining boot_rows: %s\n", e$message))
        cat("First boot_rows structure:\n")
        print(str(boot_list[[1]]))
        return(NULL)
      })
    } else {
      NULL
    }
    
    cat(sprintf("Combined bootstrap: %s rows\n", if(!is.null(all_boots)) nrow(all_boots) else "NULL"))
    
    # Collect Q10 and Q90 DVI values
    q_dvi_data <- do.call(rbind, lapply(results_list, function(res) {
      if((!is.null(res$q10_dvi) && !is.na(res$q10_dvi)) || 
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

  # Collect observed max20 per-veg across all tasks to compute reference maxima per veg.
  # Only the observed 20-day centered rolling-median maxima (res$max20_observed)
  # are used to build reference maxima for FVC. PCA/LDA-derived maxima are
  # intentionally NOT used for reference FVC calculations.
    max20_records <- do.call(rbind, lapply(results_list, function(res) {
      if(is.null(res) || is.null(res$coef_df)) return(NULL)
      loc <- res$coef_df$location_id[1]
      yr <- res$coef_df$year[1]
      out <- list()
      # Observed series maps to the TRUE vegetation for this location (if available in gpts_map)
      try({
        tv <- NA_character_
        if(exists("gpts_map")) {
          tvs <- gpts_map$Veg[gpts_map$location_id == loc]
          if(length(tvs) >= 1) tv <- as.character(tvs[1])
        }
        if(!is.na(tv) && !is.null(res$max20_observed) && is.finite(res$max20_observed)) {
          dfo <- data.frame(location_id = loc, year = yr, Veg = tv, max20_obs = as.numeric(res$max20_observed), stringsAsFactors = FALSE)
          out[[length(out)+1]] <- dfo
        }
      }, silent = TRUE)
      if(length(out) == 0) return(NULL)
      Reduce(function(a,b) merge(a,b, by = c("location_id","year","Veg"), all = TRUE), out)
    }))

    # Compute per-veg reference maxima: for each Veg, take the maximum observed max20 across all location-years
    ref_max_by_veg <- NULL
    if(!is.null(max20_records) && nrow(max20_records) > 0) {
      max20_long <- max20_records %>% dplyr::mutate(
        max20_obs = ifelse(is.finite(max20_obs), max20_obs, NA_real_)
      )
      ref_max_by_veg <- max20_long %>% dplyr::group_by(Veg) %>% dplyr::summarize(ref_max = max(max20_obs, na.rm = TRUE), .groups = 'drop')
      # If ref_max is -Inf due to all NA, set to NA
      ref_max_by_veg$ref_max[!is.finite(ref_max_by_veg$ref_max)] <- NA_real_
    }

    # Compute FVC per location-year and attach to results_list entries for later reporting
    if(!is.null(ref_max_by_veg) && nrow(ref_max_by_veg) > 0) {
      results_list <- lapply(results_list, function(res) {
        if(is.null(res) || is.null(res$coef_df)) return(res)
        loc <- res$coef_df$location_id[1]
        yr <- res$coef_df$year[1]
        # Prepare per-veg FVC estimates: observed-based FVC assigned to the TRUE veg only
        fvc_vec <- rep(NA_real_, length(veg_kept))
        names(fvc_vec) <- veg_kept
        # Determine true veg for this location (if available)
        tv <- NA_character_
        try({
          if(exists("gpts_map")) {
            tvs <- gpts_map$Veg[gpts_map$location_id == loc]
            if(length(tvs) >= 1) tv <- as.character(tvs[1])
          }
        }, silent = TRUE)
        if(!is.na(tv) && !is.null(res$max20_observed)) {
          ref_row <- ref_max_by_veg$ref_max[ref_max_by_veg$Veg == tv]
          ref_val <- if(length(ref_row) == 1 && is.finite(ref_row)) ref_row else NA_real_
          if(is.finite(res$max20_observed) && is.finite(ref_val) && ref_val != 0) {
            if(tv %in% names(fvc_vec)) fvc_vec[[tv]] <- res$max20_observed / ref_val
          }
        }
        res$fvc_by_veg <- fvc_vec
        res
      })
    }
    
    # Get unique locations
    unique_locations <- unique(all_coefs$location_id)
    
    # Join with ground truth vegetation data
    true_veg_map <- gpts_map %>% dplyr::select(location_id, true_veg = Veg)
    all_coefs <- dplyr::left_join(all_coefs, true_veg_map, by = "location_id")

    # Create a single Excel file with all locations as separate sheets
    cat("Creating single Excel file with all location results...\n")

    # Create a single workbook for all locations
    wb <- tryCatch({
      openxlsx::createWorkbook()
    }, error = function(e) {
      cat(sprintf("ERROR creating workbook: %s\n", e$message))
      return(NULL)
    })
    
    if(is.null(wb)) {
      cat("Failed to create workbook\n")
      return(NULL)
    }

    # Add a summary sheet first
    summary_data <- data.frame(
      Location_ID = unique_locations,
      Total_Years = sapply(unique_locations, function(loc) {
        length(unique(all_coefs$year[all_coefs$location_id == loc]))
      }),
      Total_Observations = sapply(unique_locations, function(loc) {
        nrow(all_coefs[all_coefs$location_id == loc, ])
      }),
      Fit_Methods = sapply(unique_locations, function(loc) {
        methods <- unique(all_coefs$FitMethod[all_coefs$location_id == loc])
        paste(methods, collapse = ", ")
      }),
      stringsAsFactors = FALSE
    )

    openxlsx::addWorksheet(wb, "Summary")
    openxlsx::writeData(wb, "Summary", summary_data)

    # If bootstrap rows exist, compute a high-level summary sheet
    if(!is.null(all_boots) && nrow(all_boots) > 0) {
      # Melt bootstrap rows into long form and compute mean/sd per veg per loc-year
      boot_long <- NULL
      try({
        boot_long <- stats::reshape(as.data.frame(all_boots), direction = "long",
                                    idvar = c("location_id","year","boot_rep","FitMethod"),
                                    varying = setdiff(names(all_boots), c("location_id","year","boot_rep","FitMethod")),
                                    timevar = "Veg",
                                    times = setdiff(names(all_boots), c("location_id","year","boot_rep","FitMethod")),
                                    v.names = "coef")
      })
      
      # Ensure reshape succeeded
      if(is.null(boot_long) || !"coef" %in% names(boot_long)) {
        stop("Bootstrap data reshaping failed - could not create proper long format data structure")
      }

      boot_summary <- boot_long %>% dplyr::group_by(location_id, year, FitMethod, Veg) %>% dplyr::summarize(
        n = sum(!is.na(coef)),
        mean_coef = ifelse(n > 0, mean(coef, na.rm = TRUE), NA_real_),
        sd_coef = ifelse(n > 1, stats::sd(coef, na.rm = TRUE), NA_real_),
        q025_coef = ifelse(n > 0, stats::quantile(coef, 0.025, na.rm = TRUE), NA_real_),
        q975_coef = ifelse(n > 0, stats::quantile(coef, 0.975, na.rm = TRUE), NA_real_),
        .groups = 'drop')
      openxlsx::addWorksheet(wb, "Bootstrap_Summary")
      openxlsx::writeData(wb, "Bootstrap_Summary", boot_summary)
    }

    # Compute best-fit summary per location-year including BOTH PCA_GD and LDA_GD results
    true_veg_map <- gpts_map %>% dplyr::select(location_id, true_veg = Veg)
    # Build best_fit_summary with two rows per location-year: PCA and LDA fits.
    # Do NOT include FVC results in this sheet (FVC is reported elsewhere).
    best_fit_summary <- do.call(rbind, lapply(unique_locations, function(loc) {
      yrs <- unique(all_coefs$year[all_coefs$location_id == loc])
      tv <- true_veg_map$true_veg[true_veg_map$location_id == loc]
      if(length(tv) == 0) tv <- NA_character_
      do.call(rbind, lapply(yrs, function(yr) {
        # PCA prediction for true veg (if present)
        pca_row <- all_coefs[all_coefs$location_id == loc & all_coefs$year == yr & all_coefs$FitMethod == "PCA_GD" & tolower(all_coefs$Veg) == tolower(tv), , drop = FALSE]
        pca_pred <- if(nrow(pca_row) == 1) pca_row$coef else NA_real_
        pca_abs_pct <- if(!is.na(pca_pred) && !is.na(tv)) abs(1 - pca_pred) * 100 else NA_real_

        # LDA prediction for true veg (if present)
        lda_row <- all_coefs[all_coefs$location_id == loc & all_coefs$year == yr & all_coefs$FitMethod == "LDA_GD" & tolower(all_coefs$Veg) == tolower(tv), , drop = FALSE]
        lda_pred <- if(nrow(lda_row) == 1) lda_row$coef else NA_real_
        lda_abs_pct <- if(!is.na(lda_pred) && !is.na(tv)) abs(1 - lda_pred) * 100 else NA_real_

        df_pca <- data.frame(
          location_id = loc,
          year = yr,
          FitMethod = "PCA_GD",
          true_veg = tv,
          pred_coef = pca_pred,
          abs_pct_diff = pca_abs_pct,
          stringsAsFactors = FALSE
        )

        df_lda <- data.frame(
          location_id = loc,
          year = yr,
          FitMethod = "LDA_GD",
          true_veg = tv,
          pred_coef = lda_pred,
          abs_pct_diff = lda_abs_pct,
          stringsAsFactors = FALSE
        )

        rbind(df_pca, df_lda)
      }))
    }))

  # Compute overall fit scores for PCA, LDA, and combined (mean of available values)
  overall_fit_pca <- suppressWarnings(as.numeric(mean(best_fit_summary$abs_pct_diff[best_fit_summary$FitMethod == "PCA_GD"], na.rm = TRUE)))
  if(!is.finite(overall_fit_pca)) overall_fit_pca <- NA_real_
  overall_fit_lda <- suppressWarnings(as.numeric(mean(best_fit_summary$abs_pct_diff[best_fit_summary$FitMethod == "LDA_GD"], na.rm = TRUE)))
  if(!is.finite(overall_fit_lda)) overall_fit_lda <- NA_real_
  # Combined: mean across both methods, ignoring NAs
  combined_abs <- best_fit_summary$abs_pct_diff
  overall_fit_combined <- suppressWarnings(as.numeric(mean(combined_abs, na.rm = TRUE)))
  if(!is.finite(overall_fit_combined)) overall_fit_combined <- NA_real_

    # Write overall fit scores at top of Summary sheet
    openxlsx::writeData(wb, "Summary", data.frame(
      Overall_Fit_PCA_pct = overall_fit_pca,
      Overall_Fit_LDA_pct = overall_fit_lda,
      Overall_Fit_Combined_pct = overall_fit_combined
    ), startRow = 1, startCol = ncol(summary_data) + 2)

    if(requireNamespace("progressr", quietly = TRUE)) {
      # Wrap progressor creation in with_progress() to avoid global environment issues
      progressr::with_progress({
        p_excel <- progressr::progressor(steps = length(unique_locations), message = "Adding location sheets")
        for (i in seq_along(unique_locations)) {
          loc_id <- unique_locations[i]

          # Create sheet name (Excel limits sheet names to 31 characters)
          sheet_name <- substr(gsub("[^A-Za-z0-9]", "_", loc_id), 1, 31)

          # Add worksheet for this location
          openxlsx::addWorksheet(wb, sheet_name)

          # --- Quality Metrics Calculation ---
          loc_coefs <- all_coefs[all_coefs$location_id == loc_id, ]

          quality_metrics <- loc_coefs %>%
            dplyr::group_by(year, FitMethod) %>%
            dplyr::summarize(
              deviation = sum(abs(coef - (tolower(Veg) == tolower(true_veg[1])))),
              .groups = 'drop'
            ) %>%
            dplyr::group_by(FitMethod) %>%
            dplyr::summarize(
              avg_pct_deviation = mean(deviation, na.rm = TRUE) * 100,
              .groups = 'drop'
            )

          # Calculate peak Q10 and Q90 DVI for this location
          loc_q_data <- if(!is.null(q_dvi_data)) q_dvi_data[q_dvi_data$location_id == loc_id, ] else NULL
          peak_q10_dvi <- if(!is.null(loc_q_data) && nrow(loc_q_data) > 0) {
            max(loc_q_data$q10_dvi, na.rm = TRUE)
          } else {
            NA
          }
          peak_q90_dvi <- if(!is.null(loc_q_data) && nrow(loc_q_data) > 0) {
            max(loc_q_data$q90_dvi, na.rm = TRUE)
          } else {
            NA
          }

          # Add peak Q10 and Q90 DVI to quality metrics
          quality_metrics$peak_q10_dvi <- peak_q10_dvi
          quality_metrics$peak_q90_dvi <- peak_q90_dvi

          # Add summary FVC metrics for this location (mean FVC across veg types where available)
          loc_fvc_vals <- unlist(lapply(results_list, function(rr) {
            if(is.null(rr) || is.null(rr$coef_df)) return(NULL)
            if(rr$coef_df$location_id[1] != loc_id) return(NULL)
            if(!is.null(rr$fvc_by_veg)) return(rr$fvc_by_veg) else return(NULL)
          }))
          # loc_fvc_vals may be a named vector repeated; coerce to numeric and compute mean across available entries per veg
          mean_fvc <- NA_real_
          try({
            # flatten to numeric values
            vals <- as.numeric(loc_fvc_vals)
            if(length(vals) > 0) mean_fvc <- mean(vals, na.rm = TRUE)
          }, silent = TRUE)
          quality_metrics$mean_fvc_by_veg <- mean_fvc

          # Write quality metrics to the sheet
          openxlsx::writeData(wb, sheet_name, "QUALITY METRICS", startRow = 1, startCol = 1)
          openxlsx::writeData(wb, sheet_name, quality_metrics, startRow = 2, startCol = 1)

          # --- Data Sections ---
          current_row <- nrow(quality_metrics) + 4

          # BEST FIT SUMMARY for this location (year-level, two rows per year: PCA_GD and LDA_GD)
          loc_best <- best_fit_summary[best_fit_summary$location_id == loc_id, , drop = FALSE]
          if(!is.null(loc_best) && nrow(loc_best) > 0) {
            # Ensure columns exist in expected order for readability; do NOT export any FVC_* columns here
            desired_cols <- c("year", "FitMethod", "true_veg", "pred_coef", "abs_pct_diff")
            # Keep location_id as first column
            write_tbl <- loc_best[, c("location_id", intersect(desired_cols, names(loc_best))), drop = FALSE]
            openxlsx::writeData(wb, sheet_name, "BEST FIT SUMMARY (per-year)", startRow = current_row, startCol = 1)
            openxlsx::writeData(wb, sheet_name, write_tbl, startRow = current_row + 1, startCol = 1)
            current_row <- current_row + nrow(write_tbl) + 3
          }

          # Filter data for the current location
          loc_coefs_pca <- all_coefs[all_coefs$location_id == loc_id & all_coefs$FitMethod == "PCA_GD", ]
          loc_coefs_lda <- all_coefs[all_coefs$location_id == loc_id & all_coefs$FitMethod == "LDA_GD", ]

          loc_boots_pca <- if (!is.null(all_boots)) all_boots[all_boots$location_id == loc_id & all_boots$FitMethod == "PCA_GD", ] else NULL
          loc_boots_lda <- if (!is.null(all_boots)) all_boots[all_boots$location_id == loc_id & all_boots$FitMethod == "LDA_GD", ] else NULL

          # Add PCA coefficients
          if (nrow(loc_coefs_pca) > 0) {
            openxlsx::writeData(wb, sheet_name, "PCA COEFFICIENTS", startRow = current_row, startCol = 1)
            openxlsx::writeData(wb, sheet_name, loc_coefs_pca, startRow = current_row + 1, startCol = 1)
            current_row <- current_row + nrow(loc_coefs_pca) + 3
          }

          # PCA bootstrap results removed to reduce file size - see Bootstrap_Summary sheet for confidence intervals

          # Add LDA coefficients
          if (nrow(loc_coefs_lda) > 0) {
            openxlsx::writeData(wb, sheet_name, "LDA COEFFICIENTS", startRow = current_row, startCol = 1)
            openxlsx::writeData(wb, sheet_name, loc_coefs_lda, startRow = current_row + 1, startCol = 1)
            current_row <- current_row + nrow(loc_coefs_lda) + 3
          }

          # LDA bootstrap results removed to reduce file size - see Bootstrap_Summary sheet for confidence intervals

          # Add Q10/Q90 DVI data if available
          if(!is.null(loc_q_data) && nrow(loc_q_data) > 0) {
            openxlsx::writeData(wb, sheet_name, "Q10/Q90 DVI TREND", startRow = current_row, startCol = 1)
            openxlsx::writeData(wb, sheet_name, loc_q_data, startRow = current_row + 1, startCol = 1)
            current_row <- current_row + nrow(loc_q_data) + 3
          }

          # DOY-level PCA/LDA predictions intentionally omitted (per request)
          if(FALSE) {
            # Placeholder: PCA/LDA DOY predictions suppressed
          }

          p_excel(sprintf("Added sheet for %s (%d/%d)", loc_id, i, length(unique_locations)))
        }
      })
    } else {
      stop("progressr package required for Excel generation - no fallback allowed")
    }
    
    # Directory and file validation - fail hard if issues found
    if(!dir.exists(OUT_DIR)) {
      stop(paste0("ERROR: Cannot create output directory: ", OUT_DIR))
    }
    
    # Check if we can write to the directory
    test_file <- file.path(OUT_DIR, "test_write.tmp")
    write_test <- tryCatch({
      writeLines("test", test_file)
      file.remove(test_file)
      TRUE
    }, error = function(e) {
      stop(paste0("ERROR: Cannot write to output directory: ", e$message))
    })
    
    cat(sprintf("Saving workbook to: %s\n", output_filename))

        # --- Data Sections ---
        current_row <- nrow(quality_metrics) + 4

        # Filter data for the current location
        loc_coefs_pca <- all_coefs[all_coefs$location_id == loc_id & all_coefs$FitMethod == "PCA_GD", ]
        loc_coefs_lda <- all_coefs[all_coefs$location_id == loc_id & all_coefs$FitMethod == "LDA_GD", ]

        loc_boots_pca <- if (!is.null(all_boots)) all_boots[all_boots$location_id == loc_id & all_boots$FitMethod == "PCA_GD", ] else NULL
        loc_boots_lda <- if (!is.null(all_boots)) all_boots[all_boots$location_id == loc_id & all_boots$FitMethod == "LDA_GD", ] else NULL

        # Add PCA coefficients
        if (nrow(loc_coefs_pca) > 0) {
          openxlsx::writeData(wb, sheet_name, "PCA COEFFICIENTS", startRow = current_row, startCol = 1)
          openxlsx::writeData(wb, sheet_name, loc_coefs_pca, startRow = current_row + 1, startCol = 1)
          current_row <- current_row + nrow(loc_coefs_pca) + 3
        }

        # Add LDA coefficients
        if (nrow(loc_coefs_lda) > 0) {
          openxlsx::writeData(wb, sheet_name, "LDA COEFFICIENTS", startRow = current_row, startCol = 1)
          openxlsx::writeData(wb, sheet_name, loc_coefs_lda, startRow = current_row + 1, startCol = 1)
          current_row <- current_row + nrow(loc_coefs_lda) + 3
        }

        # LDA bootstrap results removed to reduce file size - see Bootstrap_Summary sheet for confidence intervals

        # Add Q10/Q90 DVI data if available
          if(!is.null(loc_q_data) && nrow(loc_q_data) > 0) {
            openxlsx::writeData(wb, sheet_name, "Q10/Q90 DVI TREND", startRow = current_row, startCol = 1)
            openxlsx::writeData(wb, sheet_name, loc_q_data, startRow = current_row + 1, startCol = 1)
            current_row <- current_row + nrow(loc_q_data) + 3
          }

        # Add PCA all-year predictions if available
        loc_predictions_pca <- lapply(results_list, function(res) {
          if(res$coef_df$location_id[1] == loc_id && !is.null(res$predictions_pca)) {
            return(res$predictions_pca)
          }
          return(NULL)
        })
        loc_predictions_pca <- loc_predictions_pca[!sapply(loc_predictions_pca, is.null)]
        
        if(length(loc_predictions_pca) > 0) {
          openxlsx::writeData(wb, sheet_name, "PCA ALL-YEAR PREDICTIONS", startRow = current_row, startCol = 1)
          current_row <- current_row + 1
          
          for(pred_year in names(loc_predictions_pca[[1]])) {
            pred_data <- loc_predictions_pca[[1]][[pred_year]]
            if(!is.null(pred_data) && nrow(pred_data) > 0) {
              openxlsx::writeData(wb, sheet_name, paste("Year", pred_year), startRow = current_row, startCol = 1)
              openxlsx::writeData(wb, sheet_name, pred_data, startRow = current_row + 1, startCol = 1)
              current_row <- current_row + nrow(pred_data) + 3
            }
          }
        }

        # Add LDA all-year predictions if available
        loc_predictions_lda <- lapply(results_list, function(res) {
          if(res$coef_df$location_id[1] == loc_id && !is.null(res$predictions_lda)) {
            return(res$predictions_lda)
          }
          return(NULL)
        })
        loc_predictions_lda <- loc_predictions_lda[!sapply(loc_predictions_lda, is.null)]
        
        if(length(loc_predictions_lda) > 0) {
          openxlsx::writeData(wb, sheet_name, "LDA ALL-YEAR PREDICTIONS", startRow = current_row, startCol = 1)
          current_row <- current_row + 1
          
          for(pred_year in names(loc_predictions_lda[[1]])) {
            pred_data <- loc_predictions_lda[[1]][[pred_year]]
            if(!is.null(pred_data) && nrow(pred_data) > 0) {
              openxlsx::writeData(wb, sheet_name, paste("Year", pred_year), startRow = current_row, startCol = 1)
              openxlsx::writeData(wb, sheet_name, pred_data, startRow = current_row + 1, startCol = 1)
              current_row <- current_row + nrow(pred_data) + 3
            }
          }
        }
      }

    # Save the single workbook with all locations
    output_filename <- file.path(OUT_DIR, "all_locations_results.xlsx")
    
    # Check if output directory exists and is writable
    if(!dir.exists(OUT_DIR)) {
      cat(sprintf("Creating output directory: %s\n", OUT_DIR))
      dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
    }
    
    if(!dir.exists(OUT_DIR)) {
      cat(sprintf("ERROR: Cannot create output directory: %s\n", OUT_DIR))
      return(NULL)
    }
    
    # Check if we can write to the directory
    test_file <- file.path(OUT_DIR, "test_write.tmp")
    write_test <- tryCatch({
      writeLines("test", test_file)
      file.remove(test_file)
      TRUE
    }, error = function(e) {
      cat(sprintf("ERROR: Cannot write to output directory: %s\n", e$message))
      FALSE
    })
    
    if(!write_test) {
      return(NULL)
    }
    
    cat(sprintf("Saving workbook to: %s\n", output_filename))
    
    save_result <- tryCatch({
      openxlsx::saveWorkbook(wb, output_filename, overwrite = TRUE)
      TRUE
    }, error = function(e) {
      cat(sprintf("ERROR saving workbook: %s\n", e$message))
      FALSE
    })
    
    if(save_result) {
      cat(sprintf("Created single Excel file '%s' with %d location sheets\n", basename(output_filename), length(unique_locations)))
    } else {
      cat("Failed to save Excel file\n")
    }


# Final timing summary
timing_info$end_time <- Sys.time()
total_time <- as.numeric(difftime(timing_info$end_time, timing_info$start_time, units="secs"))

cat(sprintf("Total execution time: %.1f seconds (%.1f minutes)\n", total_time, total_time/60))
if(!is.null(timing_info$moving_var_done)) {
  cat(sprintf("Moving variance: %.1f seconds\n", 
             as.numeric(difftime(timing_info$moving_var_done, timing_info$start_time, units="secs"))))
}
if(!is.null(timing_info$lib_construction_done)) {
  cat(sprintf("Library construction: %.1f seconds\n", 
             as.numeric(difftime(timing_info$lib_construction_done, timing_info$moving_var_done, units="secs"))))
}
if(!is.null(timing_info$pca_computation_done)) {
  cat(sprintf("PCA/LDA computation: %.1f seconds\n", 
             as.numeric(difftime(timing_info$pca_computation_done, timing_info$lib_construction_done, units="secs"))))
}
cat(sprintf("Main processing + Excel: %.1f seconds\n",
           as.numeric(difftime(timing_info$end_time, timing_info$pca_computation_done, units="secs"))))
