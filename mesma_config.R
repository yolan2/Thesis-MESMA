# =============================================================================
# mesma_config.R — Global configuration constants and infrastructure utilities
# =============================================================================
# This file is sourced early in the pipeline and provides all tunable parameters,
# logging, and parallel backend setup used by downstream modules.

# Debug logging stub (no-op; legacy call sites remain harmless)
write_debug <- function(...) { invisible(NULL) }

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

INPUT_CSV <- "C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_train (2).csv"
                                 # Path to the input CSV used for training. Replace with your own file path.
INFERENCE_CSV <- "C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_mid (1).csv"
                                 # Path for spatial inference input (set to NA to disable inference steps).
                                 #
                                 # When the script is invoked with an environment variable
                                 # MESMA_INFERENCE_CSV, that path will override this value.
                                 # You may also provide a directory path or a comma-separated
                                 # list of CSV files; the pipeline will treat each file as a
                                 # separate inference run and process them sequentially.  An
                                 # alternate environment variable MESMA_INFERENCE_DIR may be
                                 # used to specify a folder containing CSVs.  This makes
                                 # "batch mode" the default when working with multiple inputs.


# Training year selection
# If you only want to build variants using the most recent data, set TRAIN_YEARS to a single year (e.g., 2024) or a vector of years.
TRAIN_YEARS <- 2024  # Years to use for training (inclusive)
# Quiet mode suppresses verbose informational/debug printing when TRUE
# Set to TRUE to reduce console noise (useful for CI or large runs)
QUIET_MODE <- TRUE


PARALLEL_ENABLE <- TRUE

PARALLEL_WORKERS <- 4            # # workers to spawn when PARALLEL_ENABLE=TRUE for general tasks (keep small for memory reasons)
PERMUTATION_PARALLEL_WORKERS <- 30  # # workers to use specifically for permutation testing (can be larger than PARALLEL_WORKERS)
# Stage-2 (pentad-level) significance threshold - less strict than index-level

# Note: the permutation worker count will be capped to available cores at runtime.
COMBO_PARALLEL_WORKERS <- max(1L, floor(PARALLEL_WORKERS/2))

# Variant generation
MAX_VARIANTS_PER_VEG <- 7       # Max variants generated per vegetation type; increase to allow greater intra-class variability.

MIN_ENDMEMBER_SAMPLES <- 5L      # Minimum raw samples required to form an endmember variant. Raise to be stricter.

ALLOWED_VEG <- c("populus", "tamarix", "herbs")
                                 # Restrict vegetation classes to these labels during training; edit to match your dataset.

# Index defaults
OPTIMAL_INDICES <- c(
  "PSRI", "NDMI", "EVI",
  "NDTI", "BSI", "TCW", "TCB",
  # Complementary indices (measure different properties than the core set)
  "MSI",   # moisture stress (SWIR1/NIR)
  "SIPI"  # pigment ratio (structure-insensitive)
)
                                 # Suggested default indices to include in fitting; add or remove indices depending on sensor and quality.

# Outlier removal
ENABLE_OUTLIER_REMOVAL <- TRUE   # Remove strong outliers prior to prototype building (recommended TRUE for noisy sensors).
OUTLIER_MAD_THRESHOLD <- 3.5     # How aggressively to trim outliers (higher -> fewer removals). ~3-4 is typical.
OUTLIER_SPLINE_MAX_DF <- 10L     # Max degrees-of-freedom for spline trend in outlier detection (caps wiggliness)

# Prototype plot options
GENERATE_PROTO_PLOTS <- TRUE     # Turn on to save one plot per prototype (useful for inspection; can generate many files)
GENERATE_PROTO_PLOTS_VARIANTS_ONLY <- TRUE  # Also save a variants-only version (no median overlays)

# Run / testing flags
# Training enabled by default (do not permanently disable training here)
# Testing mode disabled for production runs: TESTING_MODE set to FALSE
TESTING_MODE <- FALSE  # Do not set to TRUE in production; reserved for manual debugging only

# Explicitly disable DEBUG to prevent accidental verbose debugging in production
DEBUG <- FALSE


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
BOOTSTRAP_B <- 200L              # Number of bootstrap iterations (location-level). Higher -> more stable CIs but slower. 100-500 typical.
MAX_INFERENCE_LOCATIONS <- 2000L # Max locations processed per inference CSV file (single-file run cap)

# Spatial dependence handling in bootstrap
# NOTE: The default location bootstrap resamples locations i.i.d., which can
#       underestimate long-distance spatial correlation and therefore understate
#       uncertainty of spatial aggregates.
#       When enabled, we resample *spatial blocks* of locations instead.
#       Block size is determined from the data via empirical variogram fitting
#       (estimate_autocorrelation_range). BOOTSTRAP_BLOCK_KM serves only as a
#       fallback when the variogram cannot be estimated (too few locations, etc.).
ENABLE_SPATIAL_BLOCK_BOOTSTRAP <- TRUE
BOOTSTRAP_BLOCK_KM <- 30.0            # Fallback block width (km) used only when variogram estimation fails.
BOOTSTRAP_BLOCK_MAX_MISSING_FRAC <- 0.20  # If > this fraction of locations lack coords, fall back to i.i.d. bootstrap.

# Uncertainty handling
ENABLE_UNCERTAINTY <- TRUE       # Turn off to skip expensive uncertainty estimation (bootstrap, MC propagation) for faster runs.
ENABLE_MULTI_YEAR_BOOTSTRAP <- FALSE # If TRUE, run per-location multi-year bootstrap; set FALSE to skip for speed/stability.
DEBUG_UNCERTAINTY <- TRUE        # Verbose diagnostics for uncertainty steps; set FALSE for quiet production runs.

# MC propagation
ENABLE_MONTE_CARLO <- TRUE       # Propagate observation noise into coefficient uncertainty via repeated unmixing draws. Disable to save time.
MC_N_DRAWS <- 50L                # Number of MC draws per location-year (>=10 recommended). Lower to speed up, increase to stabilize CI estimates.
MC_NOISE_SCALE <- 1.0            # Multiplier applied to per-task residual RMSE to set MC noise SD. Increase to model larger noise.
# Optional explicit MC noise SD can be set via MC_NOISE_SD (if present, it overrides MC_NOISE_SCALE * rmse).

ENABLE_ENDMEMBER_BUNDLES <- TRUE # If TRUE, represent endmember variability as a mean+cov and sample endmembers during MC. Improves uncertainty realism at cost of CPU.

# OOB fraction prediction uncertainty
# Uses the OOB validation residuals (predicted - true fractions) to add systematic prediction bias/error to MC draws
ENABLE_OOB_FRACTION_UNCERTAINTY <- TRUE  # If TRUE, add OOB-derived fraction errors to MC draws

# Huber loss for FCLS (robust to outliers)
USE_HUBER_LOSS <- TRUE         # If TRUE, use Huber loss instead of RMSE in FCLS solver (more robust to outliers)
HUBER_DELTA <- 1.345             # Huber delta parameter (threshold for switching from L2 to L1).
                                 # 1.345 * sigma gives 95% efficiency for Gaussian data.
                                 # Lower values = more robust but less efficient. Higher = closer to RMSE.
HUBER_MAX_ITER <- 20             # Max iterations for IRLS (iteratively reweighted least squares)
HUBER_TOL <- 1e-4                # Convergence tolerance for IRLS

# Feature selection / pruning
ENABLE_FEATURE_PRUNING <- TRUE  # Automatically drop highly correlated features before training?
FEATURE_PRUNING_THRESHOLD <- 0.95 # Correlation threshold above which a feature is considered redundant and dropped.



# Data quality thresholds (filtering / skipping)
MIN_OBS_PER_LOC_YEAR <- 3L      # Minimum rows required per location-year for processing. Increase to be stricter on noisy cases.
MIN_UNIQUE_DOY_DEFAULT <- 5L    # Minimum unique DOYs per loc-year during training to consider time-series adequate.
MIN_UNIQUE_DOY_INFERENCE <- 3L  # Less strict for inference to allow sparse inputs.
MIN_PENTADS_PER_TRAIN_SAMPLE <- 10L  # Minimum observations required per location-year trace to construct a training sample.

# Modeling/algorithmic caps and defaults
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
INTERPOLATE_INFERENCE <- TRUE   # If TRUE, linearly interpolate missing pentads for inference/validation (default: FALSE — preserves validation integrity)
PCA_VARIANCE_THRESHOLD <- 0.95    # PCA energy to retain when reducing endmember dimensionality prior to clustering

# Whether to prune features with zero LDA weight from optimized libraries
# Disabled for now - diagnostics only, no actual pruning applied
PRUNE_ZERO_WEIGHT_FEATURES <- TRUE
# If fraction of zeroed features exceeds this, skip pruning to avoid degenerate libraries
PRUNE_ZERO_WEIGHT_MAX_FRAC <- 0.8
# Minimum number of features to keep after pruning to avoid degenerate models
PRUNE_ZERO_MIN_FEATURES <- 3

# === PERMUTATION IMPORTANCE FEATURE PRUNING ===
# Measures feature contribution by randomizing each feature and testing performance degradation.
# Features whose permutation doesn't significantly degrade performance are pruned.
OOB_TUNING_FRACTION <- 0.15         # Fraction held out from training data for cluster optimization evaluation
VALIDATION_FRACTION <- 0.20         # Fraction held out for validation (stratified by location/Veg)
PERMUTATION_N_ITER <- 400            # Number of permutations per feature for significance testing
PERMUTATION_MIN_SAMPLES <- 20       # Minimum OOB samples required for testing

# Endmember selection tuning
MAX_K_EAR <- 7L                  # Maximum number of endmembers to consider per vegetation class in EAR selection (conservative increase to explore one more k).
CLUSTER_COMPLEXITY_LAMBDA <- 0.005  # Complexity penalty per total endmember: S = min(R_oob, R_train) - λ * sum(k)
BARREN_SIM_THRESHOLD <- 0.75     # Pre-filter: drop vegetation training observations whose cosine similarity to barren mean exceeds this threshold


# === END GLOBAL CONFIG ===
PROGRESS_LOG_TO_FILE <- FALSE
LOG_FILE <- "mesma_progress.log"

RAW_BANDS <- c("blue", "green", "red", "nir", "swir1", "swir2")

OUTPUT_DIR <- "C:/MAP/phenology_results"
OUT_DIR <- file.path(OUTPUT_DIR, "veg_mixture_fit")
