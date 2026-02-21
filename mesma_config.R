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

ALLOWED_VEG <- c("populus", "tamarix", "herbs", "agriculture")
                                 # Restrict vegetation classes to these labels during training; edit to match your dataset.

# Index defaults
OPTIMAL_INDICES <- c(
  "PSRI", "NDMI", "EVI","TCW",
  "NDTI", "BSI", 
  # Complementary indices (measure different properties than the core set)
  "MSI",   # moisture stress (SWIR1/NIR)
  "SIPI"  # pigment ratio (structure-insensitive)
)
                                 # Suggested default indices to include in fitting; add or remove indices depending on sensor and quality.

# PPI normalization / canopy-maximum settings
# PPI_FULL_VEG_COVER is the index value that corresponds to full vegetation cover.
# Set to 0.7 and DO NOT compute dynamically for now (fixed canopy maximum used by add_ppi_columns()).
PPI_FULL_VEG_COVER <- 0.7


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
MAX_INFERENCE_LOCATIONS <- 2000L    # Restored to full inference scale

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
DEBUG_UNCERTAINTY <- TRUE        # Verbose diagnostics for uncertainty steps; set FALSE for quiet production runs.

# Robust loss configuration
HUBER_DELTA <- 1.0               # Huber transition point (|residual| <= delta uses L2; larger residuals use linear penalty).

# Outlier removal
ENABLE_OUTLIER_REMOVAL <- FALSE  # Enable/disable robust outlier removal in time-series data



# MC propagation
ENABLE_MONTE_CARLO <- TRUE       # Propagate observation noise into coefficient uncertainty via repeated unmixing draws. Disable to save time.
MC_N_DRAWS <- 50L                # Number of MC draws per location-year (>=10 recommended). Lower to speed up, increase to stabilize CI estimates.
MC_NOISE_SCALE <- 1.0            # Multiplier applied to per-task residual RMSE to set MC noise SD. Increase to model larger noise.
# Optional explicit MC noise SD can be set via MC_NOISE_SD (if present, it overrides MC_NOISE_SCALE * rmse).

ENABLE_ENDMEMBER_BUNDLES <- TRUE # If TRUE, represent endmember variability as a mean+cov and sample endmembers during MC. Improves uncertainty realism at cost of CPU.

# OOB fraction prediction uncertainty
# Uses the OOB validation residuals (predicted - true fractions) to add systematic prediction bias/error to MC draws
ENABLE_OOB_FRACTION_UNCERTAINTY <- TRUE  # If TRUE, add OOB-derived fraction errors to MC draws

# IWLMM (Inequality-constrained Weighted Linear Mixture Model, Li et al. 2021)
USE_IWLMM <- TRUE               # If TRUE, use IWLMM instead of standard FCLS. Allows endmember perturbations
                                 # within bounded intervals derived from within-class variance, improving fit
                                 # for pixels where fixed endmembers are suboptimal.
IWLMM_BOUND_SIGMA <- 2.0        # Perturbation bound as multiple of within-class SD per feature per endmember.
                                 # ±IWLMM_BOUND_SIGMA * sigma_j_k. Higher = more flexibility, risk of overfitting.
                                 # Typical range: 1.0 - 3.0. Set to 2.0 for a moderate trade-off.
IWLMM_MAX_ITER <- 15            # Max alternating iterations between solving fractions and perturbations.
IWLMM_TOL <- 1e-4               # Convergence tolerance for alternating optimization (max change in fractions).
IWLMM_REGULARIZE_DELTA <- 0.01  # L2 penalty on perturbation magnitude to prevent overfitting.
                                 # J_total = J_fit + lambda_delta * ||ΔM||^2. Set to 0 to disable.

# Sparse unmixing (L1-penalized solver)
USE_SPARSE_UNMIXING <- TRUE         # If TRUE, use L1-penalized unmixing: min ||Ef-y||² + λ·||f||₁  s.t. f≥0, sum(f)≤1.
                                    # Relaxes sum-to-one to sum-to-at-most-one; λ drives irrelevant endmembers to exactly 0.
SPARSE_LAMBDA <- 0.01              # L1 penalty strength. Higher = more aggressive sparsity.
                                    # Typical range: 0.001 (mild) - 0.1 (aggressive). Set to 0 for standard FCLS behavior.

# Feature selection / pruning
ENABLE_FEATURE_PRUNING <- TRUE  # Automatically drop highly correlated features before training?
FEATURE_PRUNING_THRESHOLD <- 0.95 # Correlation threshold above which a feature is considered redundant and dropped.

# PCA->LDA PC pruning: enable automatic removal of weak PCs by default.
# Set to FALSE to disable PC pruning even when LDA_PC_CUM_CONTRIB < 1.0
ENABLE_LDA_PC_PRUNING <- TRUE



# Data quality thresholds (filtering / skipping)
MIN_OBS_PER_LOC_YEAR <- 3L      # Minimum rows required per location-year for processing. Increase to be stricter on noisy cases.
MIN_UNIQUE_DOY_DEFAULT <- 5L    # Minimum unique DOYs per loc-year during training to consider time-series adequate.
MIN_UNIQUE_DOY_INFERENCE <- 3L  # Less strict for inference to allow sparse inputs.
# Minimum number of non-empty temporal pentads required for a training trace to be
# considered for endmember construction. This prefilters very sparse training
# samples before they enter the prototype/variant pipeline.
# Relaxed from 30 -> 10 to allow OOB holdout traces with sparser coverage to
# participate in cluster/OOB evaluation during library optimization.
MIN_PENTADS_PER_TRAIN_SAMPLE <- 10L

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
# Fraction of LDA discriminative strength to retain when pruning PCA components (0.95 = 95%)
LDA_PC_CUM_CONTRIB <- 0.95

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
VALIDATION_FRACTION <- 0.30         # Fraction held out for validation (stratified by location/Veg)
PERMUTATION_N_ITER <- 400            # Number of permutations per feature for significance testing
PERMUTATION_MIN_SAMPLES <- 20       # Minimum OOB samples required for testing

# Endmember selection tuning
MAX_K_EAR <- 6L                  # Maximum number of endmembers to consider per vegetation class in EAR selection (increased to allow up to 8 clusters per class).
CLUSTER_COMPLEXITY_LAMBDA <- 0.005  # Complexity penalty per total endmember: S = min(R_oob, R_train) - λ * sum(k)
BARREN_SIM_THRESHOLD <- 0.8     # Pre-filter: drop vegetation training observations whose cosine similarity to barren mean exceeds this threshold


# === END GLOBAL CONFIG ===
PROGRESS_LOG_TO_FILE <- FALSE
LOG_FILE <- "mesma_progress.log"

RAW_BANDS <- c("blue", "green", "red", "nir", "swir1", "swir2")

OUTPUT_DIR <- "C:/MAP/phenology_results"
OUT_DIR <- file.path(OUTPUT_DIR, "veg_mixture_fit")
