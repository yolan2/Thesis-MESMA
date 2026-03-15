# =============================================================================
# mesma_config.R — Global configuration constants and infrastructure utilities
# =============================================================================
# This file is sourced early in the pipeline and provides all tunable parameters,
# logging, and parallel backend setup used by downstream modules.


# --- Logging ---
log_msg <- function(...) {
  ts <- format(Sys.time(), "%H:%M:%S")
  msg <- sprintf("[%s] %s\n", ts, sprintf(...))
  cat(msg)
  if (isTRUE(PROGRESS_LOG_TO_FILE)) try(cat(msg, file = LOG_FILE, append = TRUE), silent = TRUE)
}

# --- Parallel backend ---
setup_parallel_backend <- function(workers = NULL) {
  # Capture current plan to allow restoration
  old_plan <- future::plan()

  # Determine the number of workers: respect MESMA_PARALLEL env var; otherwise use detectCores()-1
  if (is.null(workers)) {
    env <- Sys.getenv("MESMA_PARALLEL", "")
    if (nzchar(env)) {
      w <- suppressWarnings(as.integer(env))
      if (is.na(w) || w <= 0) w <- max(1L, parallel::detectCores() - 1L)
      workers <- w
    } else {
      workers <- max(1L, parallel::detectCores() - 1L)
    }
  }

  if (is.null(workers) || workers <= 1L) {
    cat("[PARALLEL] Parallel disabled (workers <= 1). Running sequentially.\n")
    return(function() { invisible(future::plan(old_plan)) })
  }

  future::plan(future::multisession, workers = workers)
  cat(sprintf("[PARALLEL] Set up 'future' multisession with %d workers\n", workers))

  cleanup <- function() {
    future::plan(old_plan)
    cat("[PARALLEL] Restored previous plan\n")
    invisible(NULL)
  }
  cleanup
}

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

INPUT_CSV <- "C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_train (3).csv"
                                 # Path to the input CSV used for training. Replace with your own file path.
INFERENCE_CSV <- "C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_kon (1).csv"
#choose "C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_kon (1).csv"
#or "C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_low (3).csv"
#or "C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_mid (2).csv"
ADDITIONAL_BIAS_CSVS <- "C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_low (3).csv"


# Training year selection
# If you only want to build variants using the most recent data, set TRAIN_YEARS to a single year (e.g., 2024) or a vector of years.
TRAIN_YEARS <- c(2023, 2024)  # Years to use for training (inclusive)
# Quiet mode suppresses verbose informational printing when TRUE
# Set to TRUE to reduce console noise (useful for CI or large runs)
QUIET_MODE <- TRUE


PARALLEL_ENABLE <- TRUE

PARALLEL_WORKERS <- 4            # # workers to spawn when PARALLEL_ENABLE=TRUE for general tasks (keep small for memory reasons)
PERMUTATION_PARALLEL_WORKERS <- 30  # # workers to use specifically for permutation testing (can be larger than PARALLEL_WORKERS)
# Stage-2 (pentad-level) significance threshold - less strict than index-level

# Note: the permutation worker count will be capped to available cores at runtime.
COMBO_PARALLEL_WORKERS <- max(1L, floor(PARALLEL_WORKERS/2))


MIN_ENDMEMBER_SAMPLES <- 5L      # Minimum raw samples required to form an endmember variant. Raise to be stricter.

ALLOWED_VEG <- c("populus", "tamarix", "herbs")
                                 # Restrict vegetation classes to these labels during training; edit to match your dataset.

# Index defaults
OPTIMAL_INDICES <- c(
  "EVI",      # Enhanced Vegetation Index (requires blue)
  "PSRI",     # Plant Senescence Reflectance Index (requires blue)
  "NDMI",     # moisture / water index (NIR-SWIR1)
  "NDTI",     # SWIR1-SWIR2 difference
  "MSI"       # moisture stress (SWIR1/NIR)
)
                                 # Suggested default indices to include in fitting; add or remove indices depending on sensor and quality.

# Outlier removal
ENABLE_OUTLIER_REMOVAL <- TRUE   # Remove strong outliers prior to prototype building (recommended TRUE for noisy sensors).
OUTLIER_MAD_THRESHOLD <- 3.5     # How aggressively to trim outliers (higher -> fewer removals). ~3-4 is typical.
OUTLIER_SPLINE_MAX_DF <- 10L     # Max degrees-of-freedom for spline trend in outlier detection (caps wiggliness)

# Prototype plot options
# When feature_weights are provided during plotting, the background of each
# temporal bin will be shaded darker proportional to the relative weight
# (higher weight = darker).  This highlights the most significant bins.
GENERATE_PROTO_PLOTS <- TRUE     # Turn on to save one plot per prototype (useful for inspection; can generate many files)
GENERATE_PROTO_PLOTS_VARIANTS_ONLY <- TRUE  # Also save a variants-only version (no median overlays)

# Run / testing flags
# Training enabled by default (do not permanently disable training here)



# -----------------------------------------------------------------------------
# USER-TUNABLE PARAMETERS (centralized)
# Move any frequently-adjusted constants here for easy configuration
# NOTE: each parameter below has a short tuning note explaining its purpose and how
#       to change it safely. Prefer conservative changes and test on a small dataset.
# -----------------------------------------------------------------------------

# Combinatorics safeguards
COMBO_SAFE_EXPAND_LIMIT <- 1e6    # Max number of total combinations to expand fully before switching to random sampling.
                                  # Increase if you want exhaustive search and have time + memory; lower to reduce runtime.
COMBO_ABORT_LIMIT <- 5e7          # Hard abort threshold: if combos exceed this, the grid search will stop to avoid OOM.
                                  # Only raise with extreme caution on large-memory hosts.

# Bootstrap / inference sizing
BOOTSTRAP_B <- 100L              # Number of bootstrap iterations (location-level). Higher -> more stable CIs but slower. 100-500 typical.
MAX_INFERENCE_LOCATIONS <- 2000L # Max locations processed per inference CSV file (single-file run cap)

# Spatial dependence handling in bootstrap
# NOTE: The default location bootstrap resamples locations i.i.d., which can
#       underestimate long-distance spatial correlation and therefore understate
#       uncertainty of spatial aggregates.
#       When enabled, we resample *spatial blocks* of locations instead.
#       BloALSEck size is determined from the data via empirical variogram fitting
#       (estimate_autocorrelation_range). BOOTSTRAP_BLOCK_KM serves only as a
#       fallback when the variogram cannot be estimated (too few locations, etc.).
ENABLE_SPATIAL_BLOCK_BOOTSTRAP <- TRUE
BOOTSTRAP_BLOCK_KM <- 30.0            # Fallback block width (km) used only when variogram estimation fails.
BOOTSTRAP_BLOCK_MAX_MISSING_FRAC <- 0.30  # If > this fraction of locations lack coords, fall back to i.i.d. bootstrap.

# Uncertainty handling
ENABLE_UNCERTAINTY <- TRUE       # Turn off to skip expensive uncertainty estimation (bootstrap, MC propagation) for faster runs.
ENABLE_MULTI_YEAR_BOOTSTRAP <- FALSE # If TRUE, run per-location multi-year bootstrap; set FALSE to skip for speed/stability.

# MC propagation
ENABLE_MONTE_CARLO <- TRUE       # Propagate observation noise into coefficient uncertainty via repeated unmixing draws. Disable to save time.
MC_N_DRAWS <- 50L                # Number of MC draws per location-year (>=10 recommended). Lower to speed up, increase to stabilize CI estimates.
MC_NOISE_SCALE <- 1.0            # Multiplier applied to per-task residual RMSE to set MC noise SD. Increase to model larger noise.
# Optional explicit MC noise SD can be set via MC_NOISE_SD (if present, it overrides MC_NOISE_SCALE * rmse).

ENABLE_ENDMEMBER_BUNDLES <- TRUE # If TRUE, represent endmember variability as a mean+cov and sample endmembers during MC. Improves uncertainty realism at cost of CPU.

# OOB fraction prediction uncertainty
# Uses the OOB validation residuals (predicted - true fractions) to add systematic prediction bias/error to MC draws
ENABLE_OOB_FRACTION_UNCERTAINTY <- TRUE  # If TRUE, add OOB-derived fraction errors to MC draws



# --- IWLMM (Iteratively Weighted Linear Mixing Model) ---
# Experimental alternate unmixing solver that perturbs endmembers within
# bounds derived from within-class variance (Li et al. 2021).  Enabling this
# mode may improve robustness to endmember variability but increases
# computational cost.  This mode can now be used together with sparse
# unmixing; when both flags are TRUE the pipeline will run a combined solver
# that alternates perturbation updates with a sparse coefficient solver.
USE_IWLMM <- TRUE        # Set TRUE to enable IWLMM unmixing during MESMA
IWLMM_MAX_ITER <- 5            # Maximum alternating optimization iterations
IWLMM_TOL <- 1e-4               # Convergence tolerance on coefficients
IWLMM_BOUND_SIGMA <- 2.0        # Multiplier on per-feature sigma to bound perturbations
IWLMM_REGULARIZE_DELTA <- 0.0   # Optional L2 regularization on endmember perturbations

# When the iterative solver computes feature weights it can produce extreme
# values that destabilise the QP.  These limits are applied after the normal
# inverse-coefficient update.  Setting either bound to NA disables clamping.
IWLMM_FEAT_W_MIN <- 0.1           # minimum allowed weight (use NA for no lower bound)
IWLMM_FEAT_W_MAX <- 10.0          # maximum allowed weight (use NA for no upper bound)

# Solver selection
DEFAULT_UNMIX_SOLVER <- "fcls"   # default fallback when MESMA_PARAMS$solver unset if IWLMM and PSRASE are not selected

# Sparse mixing controls (subset selection + sparsity-aware model score)
# When performing per-sample unmixing the default constrained QP solver may be
# replaced by an L1-penalized regression.  The flag below controls that
# behaviour; sparse library search/selection used to be governed by
# ENABLE_SPARSE_MIXING, which has been removed because it had no effect.
# Sparse subset behaviour is now driven entirely by the parameters below.
USE_SPARSE_UNMIXING <- TRUE  # Set TRUE to apply L1-penalized solver to every unmixing call.
# (removed) ENABLE_SPARSE_MIXING <- FALSE
SPARSE_MIXING_LAMBDA <- 0.01     # Sparsity penalty added to score as: rmse + lambda * n_active_components.
# Notes: n_active_components counts coefficients exceeding SPARSE_COEF_THRESHOLD.
# Barren endmember is excluded from this count so that sparsity penalties only
# affect vegetation classes.
SPARSE_MIN_COMPONENTS <- 1L      # Minimum number of active classes allowed in sparse subset search.
SPARSE_MAX_COMPONENTS <- 4L      # Maximum number of active classes allowed in sparse subset search.
SPARSE_COEF_THRESHOLD <- 0.02    # Coefficients > threshold are counted as active components for sparsity penalty.
SPARSE_UNMIX_K <- 3L              # Number of nonzero endmembers to keep per-sample when using sparse solver

# Feature selection / pruning
ENABLE_FEATURE_PRUNING <- TRUE  # Automatically drop highly correlated features before training?
FEATURE_PRUNING_THRESHOLD <- 0.95 # Correlation threshold above which a feature is considered redundant and dropped.



# Data quality thresholds (filtering / skipping)
MIN_OBS_PER_LOC_YEAR <- 8L      # Minimum rows required per location-year for processing (raised from 3). Increase to be stricter on noisy cases.
MIN_UNIQUE_DOY_DEFAULT <- 5L    # Minimum unique DOYs per loc-year during training to consider time-series adequate.
MIN_UNIQUE_DOY_INFERENCE <- 3L  # Less strict for inference to allow sparse inputs.
MIN_PENTADS_PER_TRAIN_SAMPLE <- 8L  # Minimum observations required per location-year trace to construct a training sample.

# Modeling/algorithmic caps and defaults
USE_LDA_SPACE_SOLVER <- TRUE        # Set TRUE to project observations and library into LDA space before unmixing (experimental, disabled).

ENABLE_LDA_L2_NORMALIZATION <- TRUE # L2-normalize training samples to focus on temporal shape rather than amplitude.
                                     # Set TRUE to emphasize shape; FALSE to preserve brightness differences.
ENABLE_ZSCORE_AFTER_L2 <- TRUE      # Z-score features after L2-normalization (equalizes variance across indices).
                                     # Set FALSE to skip z-scoring and preserve natural post-L2 variance differences.
                                     # Skipping can be beneficial if post-L2 variances are similar or you want LDA
                                     # to handle scale natively without amplifying low-variance (flat) indices.

GAM_K_MAX <- 10                  # Max basis dimension for GAM fits. Larger -> more flexible curves but risk overfitting.
GAM_GAMMA <- 1.0                 # Regularization / smoothing strength for GAM (higher -> smoother).

# Sample/cluster controls
# Support multiple barren/soil prototypes extracted from training bare-soil observations
RAW_BARREN_N_PROTOTYPES <- 2L  # Maximum number of barren endmembers to consider during optimization (like MAX_K_EAR for vegetation).
                                # The optimizer will search k=1..RAW_BARREN_N_PROTOTYPES for barren class.
                                # Tuning: set >1 if you expect distinct soil types/brightness regimes in the region
                                # (e.g., sandy vs. clay surfaces). Increasing this improves representativeness but
                                # will increase library size and runtime; typical values: 1-7.

# Numerical tolerances
EPS_SIGMA <- 1e-8                # Small epsilon used in numeric ops to avoid division by zero
LOWER_BND <- 0                   # Lower bound used for non-negativity constraints

# Batch processing
BATCH_SIZE <- 6  # Default batch size for location-level processing. Lower for memory-constrained runs; higher for throughput.

# Temporal (pentad) settings
TEMPORAL_AGGREGATION_DAYS <- 10L  # Temporal aggregation window in days (~pentad length). Smaller -> finer but sparser.
TEMPORAL_BUDGET <- ceiling(365 / TEMPORAL_AGGREGATION_DAYS)  # Number of temporal bins (pentads) used in modeling

# Variant / clustering settings
MIN_CLUSTER_SIZE <- 4L           # Minimum cluster size to accept a variant prototype
# Temporal-filling method used during inference/validation (and in any
# functions that accept an `interpolate` argument).  Accepted values are:
#   * "linear" (or TRUE)   - linearly interpolate missing pentads (legacy default)
#   * "whittaker"          - apply Whittaker penalized smoothing (fills gaps and
#                             smooths the resulting pentad series)
#   * "none" (or FALSE)    - leave missing pentads as NA
# The value may also be specified via `MESMA_PARAMS$interpolate_inference`.
INTERPOLATE_INFERENCE <- "linear"
# Smoothing penalty parameter used by Whittaker smoothing when
# INTERPOLATE_INFERENCE == "linear".  Larger values produce smoother curves.
WHITTAKER_LAMBDA <- 1600

PCA_VARIANCE_THRESHOLD <- 0.95    # PCA energy to retain when reducing endmember dimensionality prior to clustering

# Whether to prune library features whose LDA-derived weight is exactly zero.
# When TRUE, features (temporal bins) with zero weight are removed from each
# vegetation matrix during `precompute_optimized_library_weighted()` and the
# same pruning is repeated during inference so that the solver works on the
# reduced feature set. Pruned columns are recorded in `pruned_info` for
# diagnostics and bookkeeping. See PRUNE_ZERO_WEIGHT_MAX_FRAC and
# PRUNE_ZERO_MIN_FEATURES for guards that prevent over‑pruning.
PRUNE_ZERO_WEIGHT_FEATURES <- TRUE
# If fraction of zeroed features exceeds this, skip pruning to avoid degenerate libraries
PRUNE_ZERO_WEIGHT_MAX_FRAC <- 0.8
# Minimum number of features to keep after pruning to avoid degenerate models
PRUNE_ZERO_MIN_FEATURES <- 3

# === PERMUTATION IMPORTANCE FEATURE PRUNING ===
# Measures feature contribution by randomizing each feature and testing performance degradation.
# Features whose permutation doesn't significantly degrade performance are pruned.
OOB_TUNING_FRACTION <- 0.2         # Fraction held out from training data for cluster optimization evaluation
VALIDATION_FRACTION <- 0.20         # Fraction held out for validation (stratified by location/Veg)
PERMUTATION_N_ITER <- 500            # Number of permutations per feature for significance testing
PERMUTATION_MIN_SAMPLES <- 20       # Minimum OOB samples required for testing
PERMUTATION_ALPHA <- 0.1          # Significance level for permutation p-values (features with p >= alpha are pruned)
INDEX_PERM_ALPHA <- 0.1           # Significance level for dropping entirely index groups

# Endmember selection tuning
MAX_K_EAR <- 2L                  # Maximum number of endmembers to consider per vegetation class in EAR selection (conservative increase to explore one more k).
CLUSTER_COMPLEXITY_LAMBDA <- 0.005  # Complexity penalty per total endmember: S = min(R_oob, R_train) - λ * sum(k)
BARREN_SIM_THRESHOLD <- 1     # Pre-filter: drop vegetation training observations whose cosine similarity to barren mean exceeds this threshold

# === END GLOBAL CONFIG ===
PROGRESS_LOG_TO_FILE <- FALSE
LOG_FILE <- "mesma_progress.log"

# -----------------------------------------------------------------------------
# Spectral band configuration
# 
# NEW LOGIC: "blue" is always included in RAW_BANDS for index calculation
# (EVI, PSRI, etc.) and will never be dropped.
RAW_BANDS <- c("blue", "green", "red", "nir", "swir1", "swir2")
# EXCLUDE_BLUE_BAND is deprecated and ignored.

# Sensor-bias correction behavior in preprocessing
# When TRUE, shift LANDSAT_89 (OLI) raw bands to ETM+ scale using affine
# coefficients from satellite_bias_check/bias_stats_features.csv before index
# computation.  No additive 'mean bias' corrections are applied.
# Recommended TRUE even for Collection 2: C2 corrects radiometric calibration
# but residual surface-reflectance biases persist due to differing spectral
# response functions (OLI vs ETM+) and atmospheric correction algorithms
# (LaSRC vs LEDAPS).  These matter for temporal MESMA where both sensor
# families contribute observations to the same phenological time series.
ENABLE_BAND_BIAS_CORRECTION <- TRUE

# Sensor-bias correction behavior in preprocessing
# (deprecated) direct index-level shift for LANDSAT_89 using stats from
# satellite_bias_check/bias_stats_features.csv.  This mechanism is no longer
# supported and the value remains FALSE for compatibility.
ENABLE_DIRECT_INDEX_BIAS_CORRECTION <- FALSE

# Features excluded from direct index correction (kept case-insensitive).
# Keep raw bands excluded to avoid double-correcting them.
DIRECT_INDEX_BIAS_EXCLUDE <- c("blue", "green", "red", "nir", "swir1", "swir2")


OUTPUT_DIR <- "C:/MAP/phenology_results"
OUT_DIR <- file.path(OUTPUT_DIR, "veg_mixture_fit")
