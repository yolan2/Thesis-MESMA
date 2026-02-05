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
library(quadprog)
library(terra)
library(magrittr) # For pipe operator
if (requireNamespace("future", quietly = TRUE)) library(future)
if (requireNamespace("future.apply", quietly = TRUE)) library(future.apply)
if (requireNamespace("MASS", quietly = TRUE)) library(MASS) # For lda

# Provide a global fallback logger so functions can call write_debug() safely
if (!exists("write_debug", mode = "function")) {
  write_debug <- function(msg) {
    # Basic fallback: timestamped message to stdout for visibility
    cat(paste0(Sys.time(), " - ", as.character(msg), "\n"))
  }
}

options(warn = 1)  # print warnings as they occur for debugging
options(future.progress = FALSE)  # disable progress bars from future package
# Ensure core initialization helpers are available (defines setup_parallel_backend(), constants, etc.)
if (!exists("setup_parallel_backend", mode = "function")) {
  if (file.exists("init_mesma.R")) {
    try(source("init_mesma.R"), silent = TRUE)
  }
}
if (!exists("setup_parallel_backend", mode = "function")) {
  stop("`setup_parallel_backend` function not found after attempting to source 'init_mesma.R'. Please ensure 'init_mesma.R' is present and defines this function.")
}

INPUT_CSV <- "C:\\Users\\yolan\\Downloads\\LS_S2_Harmonized_Timeseries_training.csv"
                                 # Path to the input CSV used for training. Replace with your own file path.
INFERENCE_CSV <- "C:\\Users\\yolan\\Downloads\\LS_S2_Harmonized_Timeseries_kon.csv"
                                 # Path for spatial inference input (set to NA to disable inference steps).


# Training year selection
# If you only want to build variants using the most recent data, set TRAIN_YEARS to a single year (e.g., 2024) or a vector of years.
TRAIN_YEARS <- 2024  # Years to use for training (inclusive)
# Quiet mode suppresses verbose informational/debug printing when TRUE
# Set to TRUE to reduce console noise (useful for CI or large runs)
QUIET_MODE <- TRUE


PARALLEL_ENABLE <- TRUE

PARALLEL_WORKERS <- 3            # # workers to spawn when PARALLEL_ENABLE=TRUE for general tasks (keep small for memory reasons)
PERMUTATION_PARALLEL_WORKERS <- 30  # # workers to use specifically for permutation testing (can be larger than PARALLEL_WORKERS)
# Stage-2 (pentad-level) significance threshold - less strict than index-level

# Note: the permutation worker count will be capped to available cores at runtime.
COMBO_PARALLEL_WORKERS <- max(1L, floor(PARALLEL_WORKERS/2))

# Variant generation
MAX_VARIANTS_PER_VEG <- 7       # Max variants generated per vegetation type; increase to allow greater intra-class variability.

# Library Optimization (Random Search)
SEARCH_ITERATIONS <- 1000        # Number of random cluster combinations to evaluate during library optimization (user-requested increase)

MIN_ENDMEMBER_SAMPLES <- 5L      # Minimum raw samples required to form an endmember variant. Raise to be stricter.

ALLOWED_VEG <- c("populus", "tamarix", "herbs")
                                 # Restrict vegetation classes to these labels during training; edit to match your dataset.

# Index defaults
OPTIMAL_INDICES <- c(
  "WDVI", "GVI", "NIRv", "PSRI", "MSAVI2", "NDMI", "EVI",
  "NDTI", "SATVI", "CIG", "BSI", "TCW", "TCB", "TCG",
  # Complementary indices (measure different properties than the core set)
  "NDWI",  # open water / surface wetness (McFeeters)
  "NDBI",  # built-up / bare soil proxy
  "MSI",   # moisture stress (SWIR1/NIR)
  "SIPI",  # pigment ratio (structure-insensitive)
  "ARVI",  # atmospherically resistant vegetation
  "GNDVI"  # green NDVI (chlorophyll proxy)
)
                                 # Suggested default indices to include in fitting; add or remove indices depending on sensor and quality.
# Indices to explicitly exclude from the MESMA fitter (some indices are unstable or redundant for MESMA)
FITTER_EXCLUDE_INDICES <- c("NIRv", "OSAVI", "NDVI", "MSAVI2", "SATVI")

# Outlier removal
ENABLE_OUTLIER_REMOVAL <- TRUE   # Remove strong outliers prior to prototype building (recommended TRUE for noisy sensors).
OUTLIER_MAD_THRESHOLD <- 3.5     # How aggressively to trim outliers (higher -> fewer removals). ~3-4 is typical.

# Prototype plot options
GENERATE_PROTO_PLOTS <- TRUE     # Turn on to save one plot per prototype (useful for inspection; can generate many files)

# Run / testing flags
SKIP_INFERENCE <- FALSE        # If TRUE, skip the inference/prediction stage (useful to only perform training steps)
# Training enabled by default (do not permanently disable training here)
# Testing mode disabled for production runs: TESTING_MODE set to FALSE
TESTING_MODE <- FALSE  # Do not set to TRUE in production; reserved for manual debugging only

# Explicitly disable DEBUG to prevent accidental verbose debugging in production
DEBUG <- FALSE

# === SPLINE ENDMEMBER MODE ===
# When USE_SPLINE_ENDMEMBERS = TRUE, each vegetation class is represented by a single 
# smooth spline fitted to all training samples for that class across DOY.
# This eliminates endmember variants and uses continuous phenology curves for unmixing.
# Motivation: MESMA results often show best performance with 1 endmember per class;
# splines provide a natural continuous representation of class phenology.
USE_SPLINE_ENDMEMBERS <- FALSE   # Set TRUE to use spline-based unmixing instead of discrete endmembers
SPLINE_SPAR <- 0.7               # Smoothing parameter for splines (0-1). Lower = more wiggly, higher = smoother.
                                 # 0.7 is typical for phenology; reduce for noisier data, increase for smoother curves.
SPLINE_LAMBDA <- 10.0            # Constraint weight for sum-to-one in FCLS solver (higher = stricter constraint)

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
MAX_INFERENCE_LOCATIONS <- 2000L    # TEMP DEBUG: Reduce to e.g. 2 locations to reproduce/inspect issue quickly

# Uncertainty handling
ENABLE_UNCERTAINTY <- TRUE       # Turn off to skip expensive uncertainty estimation (bootstrap, MC propagation) for faster runs.
ENABLE_MULTI_YEAR_BOOTSTRAP <- FALSE # If TRUE, run per-location multi-year bootstrap; set FALSE to skip for speed/stability.
DEBUG_UNCERTAINTY <- TRUE        # Verbose diagnostics for uncertainty steps; set FALSE for quiet production runs.

# Classification uncertainty propagation (Dirichlet perturbation)
ENABLE_CLASSIFICATION_UNCERTAINTY <- TRUE  # If TRUE, apply Dirichlet perturbation based on confusion matrix during bootstrap.
DIRICHLET_CONCENTRATION_SCALE <- 10.0      # Multiplier for validation sample size to set Dirichlet concentration (alpha).
                                           # Higher values = tighter concentration = less classification noise.
                                           # Lower values = more spread = larger classification uncertainty.

# MC propagation
ENABLE_MONTE_CARLO <- TRUE       # Propagate observation noise into coefficient uncertainty via repeated unmixing draws. Disable to save time.
MC_N_DRAWS <- 50L                # Number of MC draws per location-year (>=10 recommended). Lower to speed up, increase to stabilize CI estimates.
MC_NOISE_SCALE <- 1.0            # Multiplier applied to per-task residual RMSE to set MC noise SD. Increase to model larger noise.
# Optional explicit MC noise SD can be set via MC_NOISE_SD (if present, it overrides MC_NOISE_SCALE * rmse).

ENABLE_ENDMEMBER_BUNDLES <- TRUE # If TRUE, represent endmember variability as a mean+cov and sample endmembers during MC. Improves uncertainty realism at cost of CPU.

# OOB fraction prediction uncertainty
# Uses the OOB validation residuals (predicted - true fractions) to add systematic prediction bias/error to MC draws
ENABLE_OOB_FRACTION_UNCERTAINTY <- TRUE  # If TRUE, add OOB-derived fraction errors to MC draws

# Endmember selection tuning
MAX_K_EAR <- 7L                  # Maximum number of endmembers to consider per vegetation class in EAR selection (conservative increase to explore one more k).
MIN_CLUSTER_SIZE <- 1L           # Minimum samples per endmember cluster. Decrease to allow more endmembers from smaller datasets.
BARREN_SIM_THRESHOLD <- 0.7     # Drop vegetation clusters whose cosine similarity to barren exceeds this threshold

# Robust clustering options
ROBUST_CLUSTERING <- TRUE      # Use robust statistics (median, MAD) in endmember selection to reduce outlier influence

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



# Small-N CI inflation (prevent tiny-sample overconfidence)
UNCERTAINTY_N_REF <- 12L         # Reference observation count; above this, no inflation. Increase if your typical loc-year has more obs.
UNCERTAINTY_N_POWER <- 1.0       # Strength of small-n inflation; reduced from 2.0 to avoid over-inflation (spatial autocorrelation now handled via n_eff).
UNCERTAINTY_BASE_SD <- 0.12      # Fallback SD for locations with no empirical SD (before inflation). Reduced from 0.15.
UNCERTAINTY_SD_MAX <- 0.40       # Cap to avoid absurdly large SD values from inflation. Reduced from 0.50.

# Option A: estimate UNCERTAINTY params by subsampling experiment (auto-calibrate)
ESTIMATE_UNCERTAINTY_PARAMS_OPTION_A <- TRUE
UNCERTAINTY_PARAM_EST_HIGH_N <- 20L            # Only use location-years with >= this many obs as 'truth' in the subsampling experiment.
UNCERTAINTY_PARAM_EST_MAX_GROUPS <- 25L        # Cap the number of loc-years used to limit runtime of the estimation.
UNCERTAINTY_PARAM_EST_TARGET_NS <- 1:15        # Subsample sizes to evaluate during estimation (typical small-n range).
UNCERTAINTY_PARAM_EST_REPS <- 20L              # Repeats per subsample size (increase for smoother estimates; slower).
UNCERTAINTY_PARAM_EST_SEED <- 123              # RNG seed for reproducible estimation runs.

# Data quality thresholds (filtering / skipping)
MIN_OBS_PER_LOC_YEAR <- 3L      # Minimum rows required per location-year for processing. Increase to be stricter on noisy cases.
MIN_UNIQUE_DOY_DEFAULT <- 5L    # Minimum unique DOYs per loc-year during training to consider time-series adequate.
MIN_UNIQUE_DOY_INFERENCE <- 3L  # Less strict for inference to allow sparse inputs.

# Modeling/algorithmic caps and defaults
ENABLE_LDA_L2_NORMALIZATION <- TRUE # L2-normalize training samples to focus on temporal shape rather than amplitude.
                                     # Set TRUE to emphasize shape; FALSE to preserve brightness differences.
                                     # Ignored if ENABLE_DUAL_REPRESENTATION is TRUE.

ENABLE_DUAL_REPRESENTATION <- FALSE  # If TRUE, include BOTH raw and L2-normalized features in the model.
                                     # This doubles the feature count but lets the model learn which
                                     # representation works best for each index. Overrides ENABLE_LDA_L2_NORMALIZATION.

GAM_K_MAX <- 10                  # Max basis dimension for GAM fits. Larger -> more flexible curves but risk overfitting.
GAM_GAMMA <- 1.0                 # Regularization / smoothing strength for GAM (higher -> smoother).

# Sample/cluster controls
# Support multiple barren/soil prototypes extracted from training bare-soil observations
RAW_BARREN_N_PROTOTYPES <- 4L  # Maximum number of barren endmembers to consider during optimization (like MAX_K_EAR for vegetation).
                                # The optimizer will search k=1..RAW_BARREN_N_PROTOTYPES for barren class.
                                # Tuning: set >1 if you expect distinct soil types/brightness regimes in the region
                                # (e.g., sandy vs. clay surfaces). Increasing this improves representativeness but
                                # will increase library size and runtime; typical values: 1-7.

MAX_PROJECTIONS_PER_VEG <- 25000L  # Subsample cap per vegetation class before clustering to avoid OOM; reduce if memory is constrained

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
PCA_VARIANCE_THRESHOLD <- 0.95    # PCA energy to retain when reducing endmember dimensionality prior to clustering

# Whether to prune features with zero LDA weight from optimized libraries
# Disabled for now - diagnostics only, no actual pruning applied
PRUNE_ZERO_WEIGHT_FEATURES <- FALSE
# If fraction of zeroed features exceeds this, skip pruning to avoid degenerate libraries
PRUNE_ZERO_WEIGHT_MAX_FRAC <- 0.8
# Minimum number of features to keep after pruning to avoid degenerate models
PRUNE_ZERO_MIN_FEATURES <- 3

# === PERMUTATION IMPORTANCE FEATURE PRUNING ===
# Measures feature contribution by randomizing each feature and testing performance degradation.
# Features whose permutation doesn't significantly degrade performance are pruned.
OOB_TUNING_FRACTION <- 0.10         # Fraction held out for cluster optimization evaluation
PERMUTATION_N_ITER <- 200            # Number of permutations per feature for significance testing
PERMUTATION_ALPHA <- 0.1          # Significance level (one-sided test)
# Higher alpha = easier to keep pentads (more lenient)
PERMUTATION_PENTAD_ALPHA <- 0.40
PERMUTATION_MIN_SAMPLES <- 20       # Minimum OOB samples required for testing


# === END GLOBAL CONFIG ===
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


RAW_BANDS <- c("blue", "green", "red", "nir", "swir1", "swir2")

# NOTE: Keep default args literal (not RAW_BANDS) so this remains worker-safe in parallel futures.
normalize_band_names <- function(df, bands = c("blue", "green", "red", "nir", "swir1", "swir2")) {
  if (is.null(df) || nrow(df) == 0) return(df)
  current_names <- names(df)
  for (b in bands) {
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


# Default SOIL_LINE_SLOPE if not yet defined (main process). Workers may not execute this.
if (!exists("SOIL_LINE_SLOPE")) SOIL_LINE_SLOPE <- 1.0

# Compute all supported indices from raw bands.
# IMPORTANT: This must be worker-safe (parallel futures), so avoid relying on globals.
compute_indices_from_bands <- function(df, raw_bands = c("blue", "green", "red", "nir", "swir1", "swir2")) {
  if (is.null(df) || nrow(df) == 0) return(df)
  eps <- 1e-9
  has_bands <- intersect(raw_bands, names(df))
  if (length(has_bands) == 0) return(df)

  soil_slope <- 1.0
  if (exists("SOIL_LINE_SLOPE", inherits = TRUE)) {
    ss <- suppressWarnings(as.numeric(get("SOIL_LINE_SLOPE", inherits = TRUE)))
    if (is.finite(ss)) soil_slope <- ss
  }

  if (all(c('nir','red') %in% names(df))) df$DVI <- as.numeric(df$nir) - as.numeric(df$red)
  if (all(c('nir','red') %in% names(df))) df$WDVI <- as.numeric(df$nir) - soil_slope * as.numeric(df$red)
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
    df$GVI <- df$TCG

  }

  # NDVI intentionally omitted from computed indices to avoid using it in MESMA fitting
  if (all(c('nir','red') %in% names(df))) df$MSAVI2 <- (2 * as.numeric(df$nir) + 1 - sqrt(pmax(0, (2 * as.numeric(df$nir) + 1)^2 - 8 * (as.numeric(df$nir) - as.numeric(df$red))))) / 2
  if (all(c('nir','swir1') %in% names(df))) df$NDMI <- (as.numeric(df$nir) - as.numeric(df$swir1)) / (as.numeric(df$nir) + as.numeric(df$swir1) + eps)

  # Enhanced Vegetation Index (EVI)
  if (all(c('nir','red','blue') %in% names(df))) {
    df$EVI <- 2.5 * ((as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + 6 * as.numeric(df$red) - 7.5 * as.numeric(df$blue) + 1 + eps))
  }

  # NDTI - Normalized Difference Tillage Index (Van Deventer et al. 1997)
  if (all(c('swir1','swir2') %in% names(df))) df$NDTI <- (as.numeric(df$swir1) - as.numeric(df$swir2)) / (as.numeric(df$swir1) + as.numeric(df$swir2) + eps)

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

  # --- Additional complementary indices ---
  # NDWI (McFeeters 1996): open water / surface wetness (green vs nir)
  if (all(c('green','nir') %in% names(df))) {
    df$NDWI <- (as.numeric(df$green) - as.numeric(df$nir)) / (as.numeric(df$green) + as.numeric(df$nir) + eps)
  }

  # NDBI (Zha et al. 2003): built-up / bare soil proxy (swir1 vs nir)
  if (all(c('swir1','nir') %in% names(df))) {
    df$NDBI <- (as.numeric(df$swir1) - as.numeric(df$nir)) / (as.numeric(df$swir1) + as.numeric(df$nir) + eps)
  }

  # MSI: moisture stress index (SWIR1/NIR)
  if (all(c('swir1','nir') %in% names(df))) {
    df$MSI <- as.numeric(df$swir1) / (as.numeric(df$nir) + eps)
  }

  # VARI: visible atmospherically resistant index (green-red / green+red-blue)
  if (all(c('green','red','blue') %in% names(df))) {
    df$VARI <- (as.numeric(df$green) - as.numeric(df$red)) / (as.numeric(df$green) + as.numeric(df$red) - as.numeric(df$blue) + eps)
  }

  # SIPI: structure insensitive pigment index (nir-blue)/(nir-red)
  if (all(c('nir','blue','red') %in% names(df))) {
    df$SIPI <- (as.numeric(df$nir) - as.numeric(df$blue)) / (as.numeric(df$nir) - as.numeric(df$red) + eps)
  }

  # ARVI: atmospherically resistant vegetation index
  if (all(c('nir','red','blue') %in% names(df))) {
    rb <- (2 * as.numeric(df$red) - as.numeric(df$blue))
    df$ARVI <- (as.numeric(df$nir) - rb) / (as.numeric(df$nir) + rb + eps)
  }

  # GNDVI: green NDVI (chlorophyll proxy)
  if (all(c('nir','green') %in% names(df))) {
    df$GNDVI <- (as.numeric(df$nir) - as.numeric(df$green)) / (as.numeric(df$nir) + as.numeric(df$green) + eps)
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


# L2-normalize a feature vector per-index
# Input: vec with n_indices * n_bins values (organized as [idx1_p1, idx1_p2, ..., idx1_pN, idx2_p1, ...])
# Output: vector of same length with each index L2-normalized independently
# This captures the SHAPE (relative temporal pattern) of each spectral index
l2_normalize_perindex <- function(vec, n_indices, n_bins) {
  if (length(vec) != n_indices * n_bins) {
    warning(sprintf("l2_normalize_perindex: vec length %d != n_indices*n_bins (%d*%d=%d)",
                    length(vec), n_indices, n_bins, n_indices * n_bins))
    return(vec)
  }

  result <- vec

  for (k in seq_len(n_indices)) {
    idx_start <- (k - 1) * n_bins + 1
    idx_end <- k * n_bins
    vals <- vec[idx_start:idx_end]
    vals_clean <- vals
    vals_clean[is.na(vals_clean)] <- 0

    nrm <- sqrt(sum(vals_clean^2))
    if (is.finite(nrm) && nrm >= 1e-9) {
      result[idx_start:idx_end] <- vals / nrm
    }
  }

  return(result)
}

# Compute diagonal covariance for endmember bundle sampling
# Uses per-band variance only (no cross-band covariance) which is more reliable
# with small sample sizes typical in endmember bundles (n=3-7)
compute_bundle_covariance <- function(Mv, verbose = FALSE) {
  # Mv: n x p matrix of endmember variant vectors (rows = variants, cols = features)
  # Returns: list(C = diagonal covariance, A = transformation matrix for sampling, mu = mean)

  Mv <- as.matrix(Mv)
  n <- nrow(Mv)
  p <- ncol(Mv)

  if (n < 2 || p < 1) {
    return(NULL)
  }

  # Compute mean
  mu <- colMeans(Mv, na.rm = TRUE)

  # Compute per-band variance (diagonal covariance)
  vars <- apply(Mv, 2, var, na.rm = TRUE)
  vars[!is.finite(vars) | vars < 1e-12] <- 1e-12

  # Diagonal covariance matrix
  C <- diag(vars, nrow = p)

  # Transformation matrix for sampling: A = diag(sqrt(vars))
  # For sampling: x = mu + A * z where z ~ N(0, I)
  A <- diag(sqrt(vars), nrow = p)

  if (verbose) {
    cat(sprintf("    [BUNDLE] n=%d variants, p=%d bands, var range=[%.4g, %.4g]\n",
                n, p, min(vars), max(vars)))
  }

  return(list(
    C = C,
    A = A,
    mu = mu,
    method = "diagonal"
  ))
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

OUTPUT_DIR <- "C:/MAP/phenology_results"
OUT_DIR <- file.path(OUTPUT_DIR, "veg_mixture_fit")

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

# Treat certain classes as herbs
veg_col <- if ("Veg" %in% names(raw_df)) "Veg" else if ("vegetation" %in% names(raw_df)) "vegetation" else NULL
if (!is.null(veg_col)) {
  # Fix typos first
  raw_df[[veg_col]] <- ifelse(raw_df[[veg_col]] == "tamairx", "tamarix", raw_df[[veg_col]])
  # Then map certain classes to herbs
  raw_df[[veg_col]] <- ifelse(raw_df[[veg_col]] %in% c("herbs", "alhagi", "salicornia", "halocnemum", "phragmites"), "herbs", raw_df[[veg_col]])
} else {
  cat("[WARNING] Neither 'Veg' nor 'vegetation' column found in raw_df. Skipping Veg treatment for training data.\n")
}

# === GEE INPUT MAPPING BLOCK ===
# Map common GEE-exported column names to the names expected by this script
# 1) Map 'vegetation' -> 'Veg'
if ("vegetation" %in% names(raw_df) && !"Veg" %in% names(raw_df)) {
  raw_df$Veg <- raw_df$vegetation
}
# 3) Ensure 'location_id' is character (GEE sometimes exports numeric)
if ("location_id" %in% names(raw_df) && !is.character(raw_df$location_id)) {
  raw_df$location_id <- as.character(raw_df$location_id)
}
# Ensure band columns like Blue/Green -> lower-case handled by normalize_band_names later
# === END GEE INPUT MAPPING BLOCK ===



df <- raw_df

# === FILTER OUT SENTINEL-2 RECORDS ===
# Keep only Landsat observations (LANDSAT_8, LANDSAT_9), exclude SENTINEL_2
if ("satellite" %in% names(df)) {
  total_before <- nrow(df)
  sentinel_count <- sum(grepl("SENTINEL", df$satellite, ignore.case = TRUE), na.rm = TRUE)
  df <- df[!grepl("SENTINEL", df$satellite, ignore.case = TRUE), , drop = FALSE]
  cat(sprintf("[SATELLITE FILTER] Filtered out %d Sentinel-2 observations (keeping Landsat only)\n", sentinel_count))
  cat(sprintf("[SATELLITE FILTER] Dataset after Sentinel filtering: %d rows from %d locations\n",
              nrow(df), length(unique(df$location_id))))
}

# Preserve zenith.angle if present (e.g. from metadata), otherwise allow recalculation
if ("zenith.angle" %in% names(df)) {
  cat("[NOTICE] Preserving existing 'zenith.angle' in input data.\n")
} else {
  df$zenith.angle <- NA_real_
}

df <- normalize_band_names(df)

if (!"date" %in% names(df) && "prediction_date" %in% names(df)) df$date <- as.Date(df$prediction_date)
if ("date" %in% names(df)) df$date <- as.Date(df$date)

# === CRITICAL: Filter out dust contamination BEFORE year filtering and PPI baseline calculation ===
# This ensures the PPI baseline is computed from clean observations only
if ("NDDI" %in% names(df)) {
  dust_count <- sum(df$NDDI > 0.18, na.rm = TRUE)
  total_before <- nrow(df)
  df <- df[!(df$NDDI > 0.18), , drop = FALSE]
  total_after <- nrow(df)
  filtered <- total_before - total_after
  cat(sprintf("[EARLY FILTERING] Filtered out %d observations with dust (NDDI > 0.18) contamination\n", filtered))
  cat(sprintf("[EARLY FILTERING] Dataset after contamination filtering: %d rows from %d locations\n",
              total_after, length(unique(df$location_id))))

  # Additionally remove extreme outliers (robust MAD-based), operating per location-year where possible
  df <- remove_large_outliers(df)

  # Defer soil line calculation until after all filtering is complete
  cat("[SOIL LINE] Soil line computation deferred until after final filtering and normalization; it will be calculated later before index computation.\n")

  # Print average images per location for each year, excluding years with < 10 observations
  if ("date" %in% names(df)) {
    year_stats <- df %>%
      dplyr::mutate(.year = lubridate::year(date)) %>%
      dplyr::group_by(.year) %>%
      dplyr::summarise(total_images = dplyr::n(),
                       n_locations = dplyr::n_distinct(location_id),
                       avg_images_per_location = total_images / n_locations,
                       .groups = "drop") %>%
      dplyr::filter(total_images >= 10)
    if (nrow(year_stats) == 0) {
      cat("[DATA STATS] No years with >= 10 observations to summarize.\n")
    } else {
      cat("[DATA STATS] Average images per location per year (years with >= 10 total obs):\n")
      print(year_stats)
    }
  } else if ("year" %in% names(df)) {
    year_stats <- df %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(total_images = dplyr::n(),
                       n_locations = dplyr::n_distinct(location_id),
                       avg_images_per_location = total_images / n_locations,
                       .groups = "drop") %>%
      dplyr::filter(total_images >= 10)
    if (nrow(year_stats) == 0) {
      cat("[DATA STATS] No years with >= 10 observations to summarize (using 'year' column).\n")
    } else {
      cat("[DATA STATS] Average images per location per year (years with >= 10 total obs, using 'year' column):\n")
      print(year_stats)
    }
  }

} else {
  cat("[WARNING] NDDI not found in data; skipping early contamination filtering\n")
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
    # Remove location-years with fewer than 5 observations entirely
    if (length(rows) < 5) {
      removed_idx[rows] <- TRUE
      next
    }
    out_mask <- rep(FALSE, nrow(sub))

    # Check if we have date for spline - require at least 10 observations
    has_date <- "date" %in% names(sub) && any(!is.na(sub$date))
    if (!has_date || length(rows) < 10) next  # skip outlier removal if spline fitting is not possible

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

        n_unique <- length(unique(x))
        fit1 <- stats::smooth.spline(x, y, df = min(10, length(x)/2, n_unique - 1))
        pred1 <- predict(fit1, x)$y
        res1 <- y - pred1
        mad1 <- stats::mad(res1, na.rm = TRUE)

        if (!is.finite(mad1) || mad1 <= 1e-6) next

        # Temporarily exclude gross outliers for the second pass (1.5x threshold)
        keep_mask <- abs(res1 - stats::median(res1, na.rm = TRUE)) <= (mad_thresh * 1.5 * mad1)

        # Pass 2: Refit on cleaner data if we have enough points left
        if (sum(keep_mask) >= 5) {
          n_unique2 <- length(unique(x[keep_mask]))
          fit2 <- stats::smooth.spline(x[keep_mask], y[keep_mask], df = min(10, sum(keep_mask)/2, n_unique2 - 1))
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

        if (!is.finite(mad_res) || mad_res <= 0) next

        # Flag outliers
        this_mask <- rep(FALSE, length(colv))
        this_mask[finite_idx] <- abs(residuals - med_res) > mad_thresh * mad_res
        out_mask <- out_mask | this_mask
      }, error = function(e) {
        # Skip this column if spline fitting fails
      })
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


# Robustly prune collinear features based on correlation matrix
prune_collinear_features <- function(df, features, threshold = 0.95) {
  if (is.null(features) || length(features) < 2) return(features)
  
  common_feats <- intersect(names(df), features)
  if (length(common_feats) < 2) return(common_feats)
  
  cat(sprintf("[FEATURE PRUNE] Checking %d features for collinearity > %.3f...\n", length(common_feats), threshold))
  
  # Compute correlation matrix
  # Use spearman for robustness to non-linearity? or pearson. Pearson is standard for strict linear redundancy.
  cm <- cor(df[, common_feats], use = "pairwise.complete.obs", method = "pearson")
  
  if (any(is.na(cm))) {
    cm[is.na(cm)] <- 0
  }
  
  diag(cm) <- 0 # Ignore self-correlation
  
  dropped <- character(0)
  
  # Greedy removal: prioritize keeping features earlier in the list (assuming user preference order)
  for (i in seq_along(common_feats)) {
    f1 <- common_feats[i]
    if (f1 %in% dropped) next
    
    for (j in (i + 1):length(common_feats)) {
      if (j > length(common_feats)) break
      f2 <- common_feats[j]
      if (f2 %in% dropped) next
      
      score <- abs(cm[f1, f2])
      if (score > threshold) {
        cat(sprintf("  Dropping '%s' (corr %.3f with '%s')\n", f2, score, f1))
        dropped <- c(dropped, f2)
      }
    }
  }
  
  kept <- setdiff(common_feats, dropped)
  cat(sprintf("[FEATURE PRUNE] Kept %d features, dropped %d.\n", length(kept), length(dropped)))
  return(kept)
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
  # Ensure PPI is present and, if possible, computed from raw bands BEFORE any normalization is applied
  if (!"PPI" %in% names(df) || all(!is.finite(df$PPI))) {
    # If we have raw bands available, derive DVI so auto_add_ppi_columns can compute PPI from raw bands
    if (all(c("nir", "red") %in% names(df))) {
      df$DVI <- as.numeric(df$nir) - as.numeric(df$red)
    }
    if (exists("auto_add_ppi_columns")) {
      ppi_try <- tryCatch({ auto_add_ppi_columns(df) }, error = function(e) { cat(sprintf("[PPI] auto_add_ppi_columns error during stored normalization: %s\n", e$message)); NULL })
      if (!is.null(ppi_try) && isTRUE(ppi_try$added)) {
        df <- ppi_try$df
        cat("[PPI] Auto-added PPI to data before applying stored normalization.\n")
      } else {
        # If auto_add failed, ensure a PPI column exists so normalization logic is consistent
        if (!"PPI" %in% names(df)) df$PPI <- NA_real_
      }
    } else {
      if (!"PPI" %in% names(df)) df$PPI <- NA_real_
    }
  } else {
    # PPI exists and appears finite - back it up as raw copy before normalization
    if (!"PPI_raw" %in% names(df)) df$PPI_raw <- df$PPI
  }
  
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
  
  # If PPI was computed but not yet backed up, ensure a 'PPI_raw' backup exists for downstream use
  if ("PPI" %in% names(df) && !"PPI_raw" %in% names(df)) df$PPI_raw <- df$PPI
  
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
          cat(sprintf("[SOIL LINE] Computed training DVI soil baseline (local only): dvi_soil = %.6f\n", dvi_soil_calc))
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

# If bands are present compute only the indices required for contamination filtering (NDDI)
# This allows us to remove dust before computing the soil line or other indices that depend on it
eps <- 1e-9
if (all(c('red','nir') %in% names(df))) df$NDDI <- (as.numeric(df$red) - as.numeric(df$nir)) / (as.numeric(df$red) + as.numeric(df$nir) + eps)


# === TRAINING-SPECIFIC: Filter out dust contamination BEFORE PPI baseline calculation ===
# Apply the same contamination filtering used for inference data to the training dataset
if ("NDDI" %in% names(df)) {
  dust_count <- sum(df$NDDI > 0.18, na.rm = TRUE)
  total_before <- nrow(df)
  df <- df[!(df$NDDI > 0.18), , drop = FALSE]
  total_after <- nrow(df)
  filtered <- total_before - total_after
  cat(sprintf("[TRAINING EARLY FILTERING] Filtered out %d observations with dust (NDDI > 0.18) contamination\n", filtered))
  cat(sprintf("[TRAINING EARLY FILTERING] Training dataset after contamination filtering: %d rows from %d locations\n",
              total_after, length(unique(df$location_id))))

  # Additionally remove extreme outliers (robust MAD-based), operating per location-year where possible
  df <- remove_large_outliers(df)

  # Now compute soil line slope from the filtered training observations and set global baseline
  compute_soil_line_slope(df)

  # Recompute indices now that SOIL_LINE_SLOPE is available (ensures WDVI uses the filtered-based slope)
  df <- compute_indices_from_bands(df)
} else {
  cat("[TRAINING EARLY FILTERING] NDDI not found in training data; skipping contamination filtering\n")
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

# --- FEATURE PRUNING (Optional) ---
if (exists("ENABLE_FEATURE_PRUNING") && isTRUE(ENABLE_FEATURE_PRUNING)) {
   thresh <- if(exists("FEATURE_PRUNING_THRESHOLD")) FEATURE_PRUNING_THRESHOLD else 0.95
   OPTIMAL_INDICES <- prune_collinear_features(df, OPTIMAL_INDICES, threshold = thresh)
}

# Backup raw PPI for visualization weighting


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
selected_vegs <- c("herbs", "populus", "tamarix", "barren")
df <- df[tolower(df$Veg) %in% selected_vegs, ]
cat(sprintf("Filtered training data to selected vegetation types: %s\n", paste(selected_vegs, collapse = ", ")))

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

cat("\n")
cat(sprintf("Dataset size: %d rows, %d columns\n", nrow(df), ncol(df)))
cat(sprintf("Number of locations: %d\n", length(unique(df$location_id))))
cat(sprintf("Date range: %s to %s\n", min(df$date, na.rm = TRUE), max(df$date, na.rm = TRUE)))
cat(sprintf("Temporal aggregation: %d days (%d bins)\n", TEMPORAL_AGGREGATION_DAYS, TEMPORAL_BUDGET))

cat("\n=== TRAIN/INFERENCE DATA CONFIGURATION ===\n")
cat(sprintf("Training years (config): %s\n", paste(TRAIN_YEARS, collapse = ", ")))
if (isTRUE(ENABLE_DUAL_REPRESENTATION)) {
  cat("Feature representation: DUAL (raw + L2-normalized)\n")
} else if (isTRUE(ENABLE_LDA_L2_NORMALIZATION)) {
  cat("Feature representation: L2-normalized only\n")
} else {
  cat("Feature representation: Raw only\n")
}

if (!"pheno_year" %in% names(df)) df$pheno_year <- assign_pheno_year(df$date)

if (exists("TRAIN_YEARS") && !is.null(TRAIN_YEARS) && length(TRAIN_YEARS) > 0) {
  cat(sprintf("Filtering training data to phenological years (March-February): %s\n", paste(TRAIN_YEARS, collapse = ", ")))
  df_train <- df[df$pheno_year %in% TRAIN_YEARS, , drop = FALSE]
} else {
  df_train <- df
}
# Apply balancing/downsampling to df_train for balanced training library

# AGENT: Apply balancing/downsampling to df_train ONLY (ensure balanced training library)
set.seed(4)
class_counts <- table(df_train$Veg)
if (length(class_counts) > 0) {
  # Downsample vegetation classes to the minimum count among NON-barren classes.
  # This prevents barren from driving the downsampling target.
  non_barren_counts <- class_counts[names(class_counts) != "barren"]
  if (length(non_barren_counts) > 0) {
    min_count_non_barren <- min(non_barren_counts)
    if (is.finite(min_count_non_barren) && min_count_non_barren > 0) {
      df_train_non_barren <- df_train %>%
        dplyr::filter(.data$Veg != "barren") %>%
        dplyr::group_by(.data$Veg) %>%
        dplyr::slice_sample(n = min_count_non_barren) %>%
        dplyr::ungroup()

      df_train_barren <- df_train %>% dplyr::filter(.data$Veg == "barren")
      df_train <- dplyr::bind_rows(df_train_non_barren, df_train_barren)

      cat(sprintf(
        "[BALANCE] Downsampled NON-barren classes to %d samples/class; kept barren as-is (total=%d)\n",
        min_count_non_barren, nrow(df_train)
      ))
    }
  }
}
cat(sprintf(
  "Training dataset (Initial): %d rows from %d locations\n",
  nrow(df_train), length(unique(df_train$location_id))
))

# --- STRATIFIED TRAIN/TEST SPLIT (80/20) ---
if (nrow(df_train) > 0) {
  cat("[SPLIT] Performing stratified 80/20 split based on location_id and Veg...\n")
  set.seed(123)
  
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
    # Log validation set class distribution
    val_class_dist <- table(test_locs_df$Veg)
    cat(sprintf("[SPLIT] Validation set class distribution: %s\n",
                paste(names(val_class_dist), "=", val_class_dist, collapse=", ")))
    
    # Verify all classes are represented
    missing_classes <- setdiff(unique_vegs, names(val_class_dist))
    if (length(missing_classes) > 0) {
      cat(sprintf("[SPLIT] WARNING: Missing classes in validation set: %s\n", 
                  paste(missing_classes, collapse=", ")))
    }
    
    # Save validation locations
    val_file <- file.path(OUT_DIR, "validation_locations.csv") # Assuming OUT_DIR is defined, or use "phenology_results/veg_mixture_fit"
    # Ensure dir exists (OUT_DIR might be defined later, so safeguard)
    if (!dir.exists(dirname(val_file))) dir.create(dirname(val_file), recursive = TRUE)
    
    # Ensure location_id is character for consistency
    test_locs_df$location_id <- as.character(test_locs_df$location_id)
    df_train$location_id <- as.character(df_train$location_id)
    
    write.csv(test_locs_df, val_file, row.names = FALSE)
    if (!isTRUE(QUIET_MODE)) cat(sprintf("[SPLIT] Saved %d validation locations to %s\n", nrow(test_locs_df), val_file))
    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG SPLIT] Sample validation location_ids: %s\n", paste(head(test_locs_df$location_id, 3), collapse=", ")))
    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG SPLIT] Sample df_train location_ids before filter: %s\n", paste(head(unique(df_train$location_id), 3), collapse=", ")))
    
    # Filter df_train to exclude test locations (Out-of-Bag)
    before_filter <- nrow(df_train)
    df_train <- df_train[!df_train$location_id %in% test_locs_df$location_id, ]
    after_filter <- nrow(df_train)
    cat(sprintf("[SPLIT] Filtered df_train: %d rows from %d locations (Training Set) - removed %d rows\n", 
                nrow(df_train), length(unique(df_train$location_id)), before_filter - after_filter))
    
    # ==========================================================================
    # CREATE SEPARATE VALIDATION AND INFERENCE DATASETS DIRECTLY
    # We create these here, immediately after validation locations are known
    # ==========================================================================
    validation_location_ids <- as.character(test_locs_df$location_id)
    
    # df_validation: validation locations, ALL years (from the original df before year filtering)
    df_validation <- df[df$location_id %in% validation_location_ids, , drop = FALSE]
    if (!isTRUE(QUIET_MODE)) cat(sprintf("[VALIDATION] Created df_validation directly: %d rows from %d locations (all years)\n",
                nrow(df_validation), length(unique(df_validation$location_id))))

    # Store validation location IDs globally for later use
    assign("validation_location_ids", validation_location_ids, envir = globalenv())
    # No validation split performed (this is the expected branch when stratified split found test locations)

    # ==========================================================================
    # CREATE OOB HOLDOUT FOR THRESHOLD TUNING (from remaining training data)
    # This OOB set is used to find optimal PCA-LDA threshold and cluster sizes
    # It is NOT the validation set - validation remains truly independent
    # ==========================================================================
    cat("\n[OOB SPLIT] Creating OOB holdout from training data for threshold/cluster tuning...\n")
    set.seed(43)  # Different seed from validation split

    # Get unique training location-Veg pairs for stratified OOB split
    train_loc_veg <- df_train %>%
      dplyr::group_by(location_id, Veg) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::arrange(location_id, dplyr::desc(n)) %>%
      dplyr::group_by(location_id) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()

    oob_locs_list <- vector("list", length(unique(train_loc_veg$Veg)))
    unique_train_vegs <- unique(train_loc_veg$Veg)

    for (i in seq_along(unique_train_vegs)) {
      v <- unique_train_vegs[i]
      v_locs <- train_loc_veg$location_id[train_loc_veg$Veg == v]
      n_v <- length(v_locs)
      n_oob <- ceiling(n_v * OOB_TUNING_FRACTION)

      if (n_oob > 0 && n_v > 1) {
        # Ensure at least 1 location remains for model training
        n_oob <- min(n_oob, n_v - 1)
        if (n_oob > 0) {
          selected <- sample(v_locs, n_oob)
          oob_locs_list[[i]] <- data.frame(location_id = selected, Veg = v, stringsAsFactors = FALSE)
        }
      }
    }

    oob_locs_df <- do.call(rbind, oob_locs_list)

    if (!is.null(oob_locs_df) && nrow(oob_locs_df) > 0) {
      oob_location_ids <- as.character(oob_locs_df$location_id)

      # Create df_train_oob (for threshold tuning and cluster optimization)
      df_train_oob <- df_train[df_train$location_id %in% oob_location_ids, , drop = FALSE]

      # Create df_train_model (for initial PCA-LDA fitting)
      df_train_model <- df_train[!df_train$location_id %in% oob_location_ids, , drop = FALSE]

      cat(sprintf("[OOB SPLIT] OOB tuning set: %d rows from %d locations\n",
                  nrow(df_train_oob), length(unique(df_train_oob$location_id))))
      cat(sprintf("[OOB SPLIT] Model training set: %d rows from %d locations\n",
                  nrow(df_train_model), length(unique(df_train_model$location_id))))

      # Log OOB class distribution
      oob_class_dist <- table(oob_locs_df$Veg)
      cat(sprintf("[OOB SPLIT] OOB class distribution: %s\n",
                  paste(names(oob_class_dist), "=", oob_class_dist, collapse=", ")))

      # Store OOB data globally
      assign("df_train_oob", df_train_oob, envir = globalenv())
      assign("df_train_model", df_train_model, envir = globalenv())
      assign("oob_location_ids", oob_location_ids, envir = globalenv())
    } else {
      cat("[OOB SPLIT] Warning: No OOB locations could be selected, using full training set\n")
      df_train_model <- df_train
      df_train_oob <- data.frame()
      oob_location_ids <- character(0)
      assign("df_train_oob", df_train_oob, envir = globalenv())
      assign("df_train_model", df_train_model, envir = globalenv())
      assign("oob_location_ids", oob_location_ids, envir = globalenv())
    }

  } else {
    cat("[SPLIT] Warning: No test locations selected.\n")
    # Fallback: empty validation, inference will be loaded from INFERENCE_CSV
    df_validation <- data.frame()
    df_inference <- NULL
    validation_location_ids <- character(0)

    # No OOB split either when no validation split
    df_train_model <- df_train
    df_train_oob <- data.frame()
    oob_location_ids <- character(0)
    assign("df_train_model", df_train_model, envir = globalenv())
    assign("df_train_oob", df_train_oob, envir = globalenv())
    assign("oob_location_ids", oob_location_ids, envir = globalenv())
  }
}

# Note: Snow/dust contamination filtering already applied early in the pipeline (before PPI baseline calculation)
# See [EARLY FILTERING] section above for details

df_test <- df
cat(sprintf(
  "Testing dataset: %d rows from %d locations\n",
  nrow(df_test), length(unique(df_test$location_id))
))

# Look for missing vegetation in original df (before train/test split) for augmentation
if ("Veg" %in% names(df) && length(ALLOWED_VEG) > 0) {
  # Check which classes are missing from df_train and add from df (original full dataset)
  original_df <- df  # df still has all data before it was assigned to df_train
  missing_vegs <- sapply(ALLOWED_VEG, function(v) {
    sum(tolower(df$Veg) == v, na.rm = TRUE)
  })
  missing_names <- names(missing_vegs)[missing_vegs == 0]
  if (length(missing_names) > 0) {
    for (mv in missing_names) {
      # Look in df_validation and df_inference combined for missing class samples
      cand <- rbind(
        if (exists("df_validation") && nrow(df_validation) > 0) df_validation[tolower(df_validation$Veg) == mv, , drop = FALSE] else data.frame(),
        if (exists("df_inference") && nrow(df_inference) > 0) df_inference[tolower(df_inference$Veg) == mv, , drop = FALSE] else data.frame()
      )
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


cat("Using training data for vegetation library construction\n")
cat("=====================================\n\n")

if (isTRUE(PARALLEL_ENABLE)) {
  cleanup_parallel <- setup_parallel_backend(workers = PARALLEL_WORKERS)
} else {
  cleanup_parallel <- function() {}
}
on.exit(cleanup_parallel(), add = TRUE)


chunked_rbind <- function(lst, chunk_size = 50L) {
  if (is.null(lst) || length(lst) == 0) return(data.frame())
  if (length(lst) == 1) return(lst[[1]])
  
  if (requireNamespace("dplyr", quietly = TRUE)) {
    return(as.data.frame(dplyr::bind_rows(lst)))
  }
  
  return(do.call(rbind, lst))
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
  X_mat_raw <- do.call(rbind, X_raw)

  n_bins_local <- TEMPORAL_BUDGET
  n_idx_local <- length(feature_cols)

  # Determine representation mode
  dual_mode <- isTRUE(ENABLE_DUAL_REPRESENTATION)
  l2_only_mode <- !dual_mode && isTRUE(ENABLE_LDA_L2_NORMALIZATION)

  if (dual_mode) {
    # DUAL REPRESENTATION: include both raw and L2-normalized features
    cat(sprintf("  DUAL REPRESENTATION: creating both raw and L2-normalized features for %d indices...\n", n_idx_local))

    # Create L2-normalized version
    X_mat_l2 <- t(apply(X_mat_raw, 1, function(r) {
      l2_normalize_perindex(r, n_idx_local, n_bins_local)
    }))

    # Concatenate: [raw_idx1, raw_idx2, ..., l2_idx1, l2_idx2, ...]
    X_mat <- cbind(X_mat_raw, X_mat_l2)

    # Create combined feature names
    l2_feature_cols <- paste0("L2norm_", feature_cols)
    all_feature_cols <- c(feature_cols, l2_feature_cols)

    cat(sprintf("  Created %d total features: %d raw + %d L2-normalized\n",
                length(all_feature_cols), length(feature_cols), length(l2_feature_cols)))
  } else if (l2_only_mode) {
    # L2 ONLY: replace raw with L2-normalized
    cat(sprintf("  L2-normalizing training samples for %d indices (%d pentads each)...\n",
                n_idx_local, n_bins_local))
    X_mat <- t(apply(X_mat_raw, 1, function(r) {
      l2_normalize_perindex(r, n_idx_local, n_bins_local)
    }))
    all_feature_cols <- feature_cols
    cat(sprintf("  L2 normalization ENABLED: using %d indices.\n", n_idx_local))
  } else {
    # RAW ONLY: use raw features as-is
    X_mat <- X_mat_raw
    all_feature_cols <- feature_cols
    cat(sprintf("  L2 normalization DISABLED: using %d raw indices.\n", n_idx_local))
  }

  n_bins <- TEMPORAL_BUDGET

  # Compute Z-score parameters for all indices
  n_total_indices <- length(all_feature_cols)
  global_means <- numeric(n_total_indices)
  global_sds <- numeric(n_total_indices)
  names(global_means) <- all_feature_cols
  names(global_sds) <- all_feature_cols

  X_z <- X_mat  # Copy structure

  cat(sprintf("  Computing Z-score parameters for %d indices...\n", n_total_indices))

  # Z-score all indices
  for(k in seq_along(all_feature_cols)) {
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
  n_classes <- length(unique(y_labels))

  # LDA requires at least (n_classes - 1) PCs for discriminant functions
  min_pcs_for_lda <- max(1, n_classes - 1)

  # Allow up to n_min - 2 PCs (need some samples beyond features for LDA stability)
  # But ensure we have at least min_pcs_for_lda
  max_pcs_for_lda <- min(20, max(min_pcs_for_lda, n_min - 2))

  if (n_pcs > max_pcs_for_lda) {
      old_n_pcs <- n_pcs
      n_pcs <- max_pcs_for_lda
      warning(sprintf("PCA->LDA: Reducing n_pcs from %d to %d to satisfy p < n_min constraint (smallest class has %d samples, max_pcs=20)", old_n_pcs, n_pcs, n_min))
  }

  if (n_pcs < min_pcs_for_lda) {
    warning(sprintf("PCA->LDA: Not enough degrees of freedom for LDA (n_pcs=%d < min_required=%d, n_min=%d). Need more samples in smallest class.", n_pcs, min_pcs_for_lda, n_min))
    return(NULL)
  }

  lda_res <- safe_lda_call(pca_res$x[, 1:n_pcs, drop=FALSE], as.factor(y_labels), min_n_pcs = min_pcs_for_lda)

  if (is.null(lda_res)) {
    # LDA failed (too few samples per class, collinearity, etc.) - fall back to uniform weights
    cat("[FALLBACK] LDA failed; using uniform weights (all features weighted equally).\n")
    final_weights <- rep(1, ncol(X_z))
    return(list(
      means = global_means,
      sds = global_sds,
      weights = final_weights,
      indices = all_feature_cols,
      base_indices = feature_cols,
      dual_mode = dual_mode,
      l2_normalize = l2_only_mode
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
  cat(sprintf("LDA weights (no normalization): min=%.4f, max=%.4f, mean=%.4f\n",
      min(final_weights[final_weights > 0], na.rm=TRUE), max(final_weights, na.rm=TRUE), mean(final_weights, na.rm=TRUE)))

  return(list(
    means = global_means,
    sds = global_sds,
    weights = final_weights,
    indices = all_feature_cols,
    base_indices = feature_cols,
    dual_mode = dual_mode,
    l2_normalize = l2_only_mode
  ))
}

doy_to_pentad <- function(doy) {
  pmin(ceiling(doy / TEMPORAL_AGGREGATION_DAYS), TEMPORAL_BUDGET)
}

build_pentad_matrix <- function(dly_year, avail_idx, interpolate = TRUE) {
  if (is.null(dly_year) || nrow(dly_year) == 0) return(NULL)

  # CRITICAL: Use phenological DOY (March 1 = day 1), not calendar DOY
  # This ensures temporal alignment when data spans phenological years (March-February)
  if (!"doy" %in% names(dly_year) || any(is.na(dly_year$doy))) {
    # Local implementation (keeps this function self-contained for future workers)
    local_pheno_doy <- function(d) {
      d <- as.Date(d)
      yr <- as.integer(format(d, "%Y"))
      march1 <- as.Date(paste0(yr, "-03-01"))
      pheno_start <- ifelse(d >= march1, march1, as.Date(paste0(yr - 1L, "-03-01")))
      as.integer(d - as.Date(pheno_start) + 1L)
    }
    dly_year$doy <- local_pheno_doy(dly_year$date)
  }

  dly_year$pentad <- doy_to_pentad(dly_year$doy)

  K <- length(avail_idx)
  pentad_mat <- matrix(NA_real_, nrow = TEMPORAL_BUDGET, ncol = K)
  colnames(pentad_mat) <- avail_idx

  for (p in 1:TEMPORAL_BUDGET) {
    subset_p <- dly_year[dly_year$pentad == p, ]
    if (nrow(subset_p) == 0) next

    # Pentad center in phenological DOY space (March 1 = 1), using the nominal bin boundaries.
    # Use a linear intra-pentad trend for the representative pentad value:
    # yi,b = beta0 + beta1*(ti - t_center) + eps, and take beta0 as the pentad center value.
    t_start <- (p - 1) * TEMPORAL_AGGREGATION_DAYS + 1
    t_end <- min(p * TEMPORAL_AGGREGATION_DAYS, TEMPORAL_BUDGET * TEMPORAL_AGGREGATION_DAYS)
    t_center <- (t_start + t_end) / 2

    for (j in seq_along(avail_idx)) {
      idx <- avail_idx[j]
      if (!idx %in% names(subset_p)) next

      v <- subset_p[[idx]]
      v <- v[is.finite(v)]
      if (length(v) == 0) next

      # Always use Case B (linear intra-pentad trend) and remove quantile clipping.
      # When time variation is insufficient (e.g., all doys identical), this reduces to an intercept-only fit.
      doys <- subset_p$doy
      doys <- doys[is.finite(doys)]
      if (length(doys) != length(v)) {
        n_use <- min(length(doys), length(v))
        doys <- doys[seq_len(n_use)]
        v <- v[seq_len(n_use)]
      }

      if (length(v) == 1 || length(unique(doys)) < 2) {
        # Degenerate case: no slope information; representative is the intercept-only estimate.
        pentad_mat[p, j] <- mean(v, na.rm = TRUE)
      } else {
        x <- doys - t_center
        model <- tryCatch(stats::lm(v ~ x), error = function(e) NULL)
        if (!is.null(model)) {
          b0 <- tryCatch(stats::coef(model)[[1]], error = function(e) NA_real_)
          if (is.finite(b0)) {
            pentad_mat[p, j] <- as.numeric(b0)
          } else {
            pentad_mat[p, j] <- mean(v, na.rm = TRUE)
          }
        } else {
          pentad_mat[p, j] <- mean(v, na.rm = TRUE)
        }
      }
    }
  }

  # Only interpolate missing values if interpolate=TRUE (for training endmembers)
  # For inference/validation data, keep NAs to use only actual observations
  if (interpolate) {
    for (j in 1:K) {
      vals <- pentad_mat[, j]
      if (any(is.na(vals))) {
        if (all(is.na(vals))) {
          pentad_mat[, j] <- 0
        } else {
          idx_present <- which(!is.na(vals))
          if (length(idx_present) >= 2) {
            pentad_mat[, j] <- approx(idx_present, vals[idx_present], xout = 1:TEMPORAL_BUDGET, rule = 2)$y
          } else {
            pentad_mat[, j] <- vals[idx_present[1]]
          }
        }
      }
    }
  }

  pentad_mat
}

 

apply_pca_lda_transform <- function(y, pca_lda_result) {
  # NO-OP: PCA-LDA weights are applied ONLY in the solver (via feature_weights)
  # This function now just returns the z-scored observation unchanged
  # The solver will apply weights directly to both endmembers and observations
  return(y)
}


.run_map <- function(X, FUN, show_pb = TRUE) {
  f_FUN <- FUN
  
  if (!PARALLEL_ENABLE) {
    lapply(X, function(x) { f_FUN(x) })
  } else {
    # Reuse existing plan if available
    if (requireNamespace("future.apply", quietly = TRUE)) {
      stop("future.apply and future packages required for parallel processing")
    }

    # AGENT CHANGE: Reuse existing plan if available to avoid costly cluster restart
    current_plan <- future::plan()
    # Check if we are running sequentially (or default). If so, we need to spin up a cluster.
    # If we are already running multisession/cluster, we reuse it.
    plan_is_sequential <- inherits(current_plan, "sequential") || inherits(current_plan, "uniprocess")
    
    # Using explicit cluster with sequential setup for Windows stability
    if (plan_is_sequential) {
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
    }
    if (!is.character(gpts_map$location_id)) {
      gpts_map$location_id <- as.character(gpts_map$location_id)
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

# Ensure derived indices are available from raw bands (do not rely on precomputed columns)
df <- normalize_band_names(df)
df <- compute_indices_from_bands(df)

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




  user_dvi_soil <- suppressWarnings(as.numeric(ifelse(nzchar(Sys.getenv("MESMA_DVI_SOIL")), Sys.getenv("MESMA_DVI_SOIL"), NA)))

  # Use add_ppi_columns to compute per-location baselines from median of lowest 10% of ALL DVI values across all years
  # If location_id is present, it will compute per-location baselines; otherwise, it will use a single baseline
  if (!is.na(user_dvi_soil) && is.finite(user_dvi_soil)) {
    cat(sprintf("[PPI] Auto-adding PPI using MESMA_DVI_SOIL override: %.6f\n", user_dvi_soil))
    df <- add_ppi_columns(df, dvi_soil = user_dvi_soil)
  } else if (sum(is.finite(df$DVI)) > 0) {
    cat("[PPI] Auto-adding PPI: computing per-location baselines from median of lowest 10% of ALL DVI values across all years\n")
    # Let add_ppi_columns compute per-location baselines automatically
    df <- add_ppi_columns(df)
  } else {
    cat("[NOTICE] Candidate indices computed from existing indices and raw bands; PPI could not be auto-added (no valid DVI values found).\n")
  }
} else {
  # PPI already present in input data
  cat("[NOTICE] Candidate indices computed from existing indices and raw bands; PPI column detected in input and will be used.\n")
}

if (length(candidate_indices) == 0) {
  stop("No OPTIMAL_INDICES or RAW_BANDS present in input data")
}

cat(sprintf("Selected %d indices: %s\n", length(candidate_indices), paste(candidate_indices, collapse = ", ")))

avail <- candidate_indices
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

timing_info$moving_var_done <- Sys.time()



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
if (FALSE) cat("[DEBUG] Defining functions...\n")
lib <- list()




if (FALSE) cat("[DEBUG] Functions defined successfully.\n")
cat("===============================================\n\n")

lib_df <- df
# No artificial per-vegetation sampling applied; use full lib_df for variant construction

vegs <- unique(na.omit(lib_df$Veg))
vegs <- vegs[vegs != ""]
vegs <- vegs[tolower(vegs) %in% c("herbs", "populus", "tamarix", "barren")]  # FIXED: case-insensitive matching

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

  T_medoid <- matrix(NA_real_, nrow = TEMPORAL_BUDGET, ncol = n_features)
  n_filled <- 0

  for (p in seq_len(TEMPORAL_BUDGET)) {
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
              n_filled, TEMPORAL_BUDGET))
}

cat("Raw index library templates computed.\n")

# --- GLOBAL HELPER: Bootstrap Median Vector ---
# Defined at global scope so all bootstrap functions can use it
boot_median_vec <- function(x, n_boot) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(rep(NA_real_, n_boot))
  if (length(x) < 3) return(rep(median(x), n_boot)) # Too few to bootstrap meaningfully

  # Fast vectorized bootstrap of median
  replicate(n_boot, median(sample(x, length(x), replace = TRUE), na.rm = TRUE))
}

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

  # Defensive: ensure ppi_detrended column exists to avoid data.table warnings when referenced
  if (!"ppi_detrended" %in% names(df_tasks)) df_tasks$ppi_detrended <- NA_real_

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
      # Ensure numeric and finite handling
      summer_df$ppi_detrended <- as.numeric(summer_df$ppi_detrended)
      cat(sprintf("[PPI BOOTSTRAP] Computed detrended summer PPI (N=%d).\n", nrow(summer_df)))
    }, error = function(e) {
      cat(sprintf("[PPI BOOTSTRAP] Detrending failed: %s. Using raw PPI.\n", e$message))
      summer_df$ppi_detrended <- as.numeric(summer_df[[ppi_col]])
    })
  } else {
    summer_df$ppi_detrended <- as.numeric(summer_df[[ppi_col]])
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
  # Remove pre-existing Barren rows and NA/empty Veg values to avoid confusion
  veg_char <- as.character(all_coefs$Veg)
  valid_veg_mask <- !is.na(veg_char) &
                    nchar(trimws(veg_char)) > 0 &
                    trimws(veg_char) != "NA" &
                    !tolower(trimws(veg_char)) %in% c("barren")
  merged_veg <- all_coefs[valid_veg_mask, ]

  # Debug: show veg types before filtering
  veg_before <- unique(merged_veg$Veg)
  cat(sprintf("[PPI BOOTSTRAP DEBUG] Veg types before PPI join: %s (n=%d rows)\n", paste(veg_before, collapse=", "), nrow(merged_veg)))

  # Create key for joining
  merged_veg$key <- paste(merged_veg$location_id, merged_veg$pheno_year, sep = "_")
  merged_veg$ppi_row_idx <- key_map[merged_veg$key]

  # Debug: show how many rows match PPI data
  n_matched <- sum(!is.na(merged_veg$ppi_row_idx))
  n_total <- nrow(merged_veg)
  cat(sprintf("[PPI BOOTSTRAP DEBUG] PPI key match: %d/%d rows (%.1f%%)\n", n_matched, n_total, 100*n_matched/max(1,n_total)))

  # Filter out rows with no PPI data
  merged_veg <- merged_veg[!is.na(merged_veg$ppi_row_idx), ]

  # Debug: show veg types after filtering
  veg_after <- unique(merged_veg$Veg)
  cat(sprintf("[PPI BOOTSTRAP DEBUG] Veg types after PPI join: %s\n", paste(veg_after, collapse=", ")))
  
  # Calculate relative coefficients (fractions of vegetation)
  # IMPORTANT: Skip normalization for 100% barren cases to preserve relative proportions
  merged_veg <- merged_veg |>
    dplyr::group_by(location_id, pheno_year) |>
    dplyr::mutate(
      sum_veg_coef = sum(coef[tolower(Veg) != "barren"], na.rm = TRUE),
      # Only normalize if there's actual vegetation (sum > threshold)
      # If 100% barren (sum_veg_coef ~ 0), keep original coefficients as rel_coef (preserves relative proportions)
      rel_coef = ifelse(sum_veg_coef > 1e-6, coef / sum_veg_coef, coef)
    ) |>
    dplyr::ungroup() |>
    # Renormalize to sum-to-one ONLY if the sum deviates from 1 AND there's actual vegetation
    # This preserves relative proportions - we scale all values by the same factor
    dplyr::group_by(location_id, pheno_year) |>
    dplyr::mutate(
      total_rel = sum(rel_coef, na.rm = TRUE),
      # Only renormalize if: (1) sum is not already ~1, (2) sum is > threshold (not 100% barren)
      rel_coef = ifelse(total_rel > 1e-6 & abs(total_rel - 1) > 1e-9,
                        rel_coef / total_rel,
                        rel_coef)
    ) |>
    dplyr::select(-total_rel) |>
    dplyr::ungroup()

  # Diagnostic: report any location-years that required correction
  rels_check <- merged_veg |>
    dplyr::group_by(location_id, pheno_year) |>
    dplyr::summarize(total_rel = sum(rel_coef, na.rm = TRUE), .groups = "drop")
  bad_rel <- sum(abs(rels_check$total_rel - 1) > 1e-6 & rels_check$total_rel > 1e-6, na.rm = TRUE)
  if (bad_rel > 0) cat(sprintf("[NORMALIZE COEF] Corrected rel_coef for %d location-year(s) to sum-to-one\n", bad_rel))

  # Report 100% barren cases that were skipped
  barren_100_count <- sum(rels_check$total_rel < 1e-6, na.rm = TRUE)
  if (barren_100_count > 0) cat(sprintf("[PPI BOOTSTRAP] Skipped normalization for %d location-year(s) with 100%% barren\n", barren_100_count))

  veg_types <- unique(merged_veg$Veg[tolower(merged_veg$Veg) != "barren"])
  veg_types <- veg_types[!is.na(veg_types) & nchar(trimws(veg_types)) > 0]  # Remove NA and empty strings
  years <- sort(unique(merged_veg$pheno_year[!is.na(merged_veg$pheno_year)]))
  locations <- unique(merged_veg$location_id)
  n_locs <- length(locations)

  # Guard against empty data
  if (length(years) == 0 || n_locs == 0 || nrow(merged_veg) == 0) {
    warning("[PPI BOOTSTRAP] No valid data after filtering - skipping bootstrap")
    return(all_coefs)
  }

  # Additional guard: ensure B is positive
  if (!is.numeric(B) || length(B) != 1 || B < 1) {
    warning("[PPI BOOTSTRAP] Invalid B parameter - skipping bootstrap")
    return(all_coefs)
  }
  B <- as.integer(B)

  # Pre-calculate unique loc-year list for Barren calculation
  # Each unique loc-year needs a row index into ppi_boot_mat
  unique_loc_years <- merged_veg |>
    dplyr::distinct(location_id, pheno_year, ppi_row_idx)

  cat(sprintf("[PPI BOOTSTRAP] Running %d bootstrap iterations with per-location PPI uncertainty...\n", B))

  # We need to accumulate results. Instead of matrix of means, we can just run the loop B times
  # and store the global mean for that iteration.

  # Initialize storage for global aggregates per veg/year
  # veg_boot_res[[veg]][[year]] -> vector of B means
  n_years <- length(years)
  cat(sprintf("[PPI BOOTSTRAP DEBUG] n_years=%d, B=%d, n_veg_types=%d\n", n_years, B, length(veg_types)))

  if (n_years == 0) {
    warning("[PPI BOOTSTRAP] No years found after filtering - cannot create result matrices")
    return(all_coefs)
  }

  veg_boot_res <- list()
  all_veg_plus_barren <- c(veg_types, "Barren")
  for (v in all_veg_plus_barren) {
    # Create matrix and explicitly ensure it's 2D
    mat <- matrix(NA_real_, nrow = as.integer(B), ncol = as.integer(n_years))
    if (is.matrix(mat) && ncol(mat) > 0) {
      colnames(mat) <- as.character(years)
    }
    veg_boot_res[[v]] <- mat
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
  # Skip if matrix is empty

  # --- Compute "woody" (populus + tamarix + woody_unknown) aggregation with proper CI propagation ---
  # Since most uncertainty is between populus and tamarix, summing their bootstrap
  # replicates BEFORE computing quantiles will reduce CI width appropriately.
  # woody_unknown is included as it represents indistinguishable populus/tamarix variants
  woody_types <- c("populus", "tamarix", "woody_unknown")
  woody_mats <- veg_boot_res[tolower(names(veg_boot_res)) %in% woody_types]
  if (length(woody_mats) >= 1) {
    # Sum the bootstrap matrices element-wise
    woody_mat <- Reduce(`+`, lapply(woody_mats, function(m) {
      m[is.na(m)] <- 0  # Treat NA as 0 for summing
      m
    }))
    # Restore NAs where ALL constituent types were NA
    all_na_mask <- Reduce(`&`, lapply(woody_mats, is.na))
    woody_mat[all_na_mask] <- NA_real_
    colnames(woody_mat) <- colnames(woody_mats[[1]])
    veg_boot_res[["woody"]] <- woody_mat
    cat("[PPI BOOTSTRAP] Added 'woody' category (populus + tamarix + woody_unknown) with combined bootstrap CIs\n")
  }

  # Compile Final Results
  final_results <- list()
  for (v in names(veg_boot_res)) {
    mat <- veg_boot_res[[v]]
    if (is.null(mat) || ncol(mat) == 0 || is.null(colnames(mat)) || length(colnames(mat)) == 0) next

    df_res <- data.frame(
      year = as.integer(colnames(mat)),
      Veg = v,
      global_coef = apply(mat, 2, mean, na.rm = TRUE),
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

  result <- dplyr::bind_rows(final_results)
  # Replace "phragmites" with "herbs" in results
  if (!is.null(result) && "Veg" %in% names(result)) {
    result$Veg <- ifelse(tolower(trimws(result$Veg)) == "phragmites", "herbs", result$Veg)
  }
  result
}

# --- NEW: NDVI-BASED LOCATION BOOTSTRAP (Analogous to PPI) ---
# Uses per-location summer NDVI medians and bootstraps medians like PPI.
location_bootstrap_ndvi <- function(all_coefs, df_tasks, B = BOOTSTRAP_B, seed = 123, ndvi_max = 0.6) {
  set.seed(seed)
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    warning("dplyr required for NDVI bootstrap")
    return(NULL)
  }

  ndvi_col <- "NDVI"
  if ("NDVI_raw" %in% names(df_tasks) && any(!is.na(df_tasks$NDVI_raw))) {
    ndvi_col <- "NDVI_raw"
    cat("[NDVI BOOTSTRAP] Using raw NDVI values (NDVI_raw) for weighting.\n")
  } else if ("NDVI" %in% names(df_tasks)) {
    cat("[NDVI BOOTSTRAP] Using standard NDVI column.\n")
  } else {
    warning("NDVI column missing from df_tasks")
    return(NULL)
  }

  # Prepare summer detrended NDVI similar to PPI pipeline
  if (!"month" %in% names(df_tasks)) {
    if ("date" %in% names(df_tasks)) df_tasks$month <- as.integer(format(as.Date(df_tasks$date), "%m")) else df_tasks$month <- 1L
  }
  if (!"doy" %in% names(df_tasks)) {
    if ("date" %in% names(df_tasks)) df_tasks$doy <- as.numeric(format(as.Date(df_tasks$date), "%j")) else df_tasks$doy <- 150
  }

  summer_df <- df_tasks |> dplyr::filter(month %in% 6:9, get(ndvi_col) > -Inf)
  if (nrow(summer_df) > 50) {
    tryCatch({
      f <- as.formula(paste(ndvi_col, "~ poly(doy, 3)"))
      seasonal_model <- lm(f, data = summer_df)
      summer_df$seasonal_trend <- predict(seasonal_model, newdata = summer_df)
      global_seasonal_mean <- mean(summer_df$seasonal_trend, na.rm = TRUE)
      summer_df$ndvi_detrended <- summer_df[[ndvi_col]] - (summer_df$seasonal_trend - global_seasonal_mean)
      summer_df$ndvi_detrended <- as.numeric(summer_df$ndvi_detrended)
      cat(sprintf("[NDVI BOOTSTRAP] Computed detrended summer NDVI (N=%d).\n", nrow(summer_df)))
    }, error = function(e) {
      cat(sprintf("[NDVI BOOTSTRAP] Detrending failed: %s. Using raw NDVI.\n", e$message))
      summer_df$ndvi_detrended <- as.numeric(summer_df[[ndvi_col]])
    })
  } else {
    summer_df$ndvi_detrended <- as.numeric(summer_df[[ndvi_col]])
    cat("[NDVI BOOTSTRAP] Insufficient summer data for detrending; using raw values.\n")
  }

  # Safe summarise: ensure we only boot groups with finite detrended values and provide diagnostics on failure
  # Build per-group list via split+lapply to avoid dplyr evaluation scoping issues
  if (nrow(summer_df) == 0) {
    boot_combined <- data.frame(location_id = character(0), pheno_year = integer(0), ndvi_boot_list = I(list()))
  } else {
    # Use interaction with a unique separator to avoid issues with dots in location_id
    summer_df$grp_key <- paste(summer_df$location_id, summer_df$pheno_year, sep = "|||")
    grp <- split(summer_df, summer_df$grp_key, drop = TRUE)
    grp_names <- names(grp)
    ndvi_lists <- vector("list", length(grp))
    for (gi in seq_along(grp)) {
      g <- grp[[gi]]
      vec <- as.numeric(g$ndvi_detrended)
      if (sum(is.finite(vec)) > 0) {
        ndvi_lists[[gi]] <- boot_median_vec(vec, B)
      } else {
        ndvi_lists[[gi]] <- rep(NA_real_, as.integer(B))
      }
    }
    # Parse names back to location_id / pheno_year using the unique separator
    parsed <- do.call(rbind, strsplit(grp_names, "\\|\\|\\|"))
    boot_combined <- data.frame(location_id = parsed[,1], pheno_year = as.integer(parsed[,2]), ndvi_boot_list = I(ndvi_lists), stringsAsFactors = FALSE)
  }
  boot_combined$key <- paste(boot_combined$location_id, boot_combined$pheno_year, sep = "_")
  n_rows <- nrow(boot_combined)
  ndvi_boot_mat <- matrix(NA_real_, nrow = n_rows, ncol = as.integer(B))
  for (i in 1:n_rows) {
    s_vec <- boot_combined$ndvi_boot_list[[i]]
    if (!is.null(s_vec) && !all(is.na(s_vec))) ndvi_boot_mat[i, ] <- s_vec
  }
  key_map <- setNames(seq_len(n_rows), boot_combined$key)

  # Prepare merged veg data like PPI (exclude barren rows)
  veg_char <- as.character(all_coefs$Veg)
  valid_veg_mask <- !is.na(veg_char) & nchar(trimws(veg_char)) > 0 & trimws(veg_char) != "NA" & !tolower(trimws(veg_char)) %in% c("barren")
  merged_veg <- all_coefs[valid_veg_mask, ]
  merged_veg$key <- paste(merged_veg$location_id, merged_veg$pheno_year, sep = "_")
  merged_veg$ndvi_row_idx <- key_map[merged_veg$key]
  merged_veg <- merged_veg[!is.na(merged_veg$ndvi_row_idx), ]

  # Calculate relative coefficients (fractions of vegetation) - same as PPI bootstrap
  # IMPORTANT: Skip normalization for 100% barren cases to preserve relative proportions
  merged_veg <- merged_veg |>
    dplyr::group_by(location_id, pheno_year) |>
    dplyr::mutate(
      sum_veg_coef = sum(coef[tolower(Veg) != "barren"], na.rm = TRUE),
      # Only normalize if there's actual vegetation (sum > threshold)
      # If 100% barren (sum_veg_coef ~ 0), keep original coefficients as rel_coef (preserves relative proportions)
      rel_coef = ifelse(sum_veg_coef > 1e-6, coef / sum_veg_coef, coef)
    ) |>
    dplyr::ungroup() |>
    # Renormalize to sum-to-one ONLY if the sum deviates from 1 AND there's actual vegetation
    # This preserves relative proportions - we scale all values by the same factor
    dplyr::group_by(location_id, pheno_year) |>
    dplyr::mutate(
      total_rel = sum(rel_coef, na.rm = TRUE),
      # Only renormalize if: (1) sum is not already ~1, (2) sum is > threshold (not 100% barren)
      rel_coef = ifelse(total_rel > 1e-6 & abs(total_rel - 1) > 1e-9,
                        rel_coef / total_rel,
                        rel_coef)
    ) |>
    dplyr::select(-total_rel) |>
    dplyr::ungroup()

  # Diagnostic: report any location-years that required correction
  rels_check <- merged_veg |>
    dplyr::group_by(location_id, pheno_year) |>
    dplyr::summarize(total_rel = sum(rel_coef, na.rm = TRUE), .groups = "drop")
  bad_rel <- sum(abs(rels_check$total_rel - 1) > 1e-6 & rels_check$total_rel > 1e-6, na.rm = TRUE)
  if (bad_rel > 0) cat(sprintf("[NORMALIZE COEF] Corrected rel_coef for %d location-year(s) to sum-to-one\n", bad_rel))

  # Report 100% barren cases that were skipped
  barren_100_count <- sum(rels_check$total_rel < 1e-6, na.rm = TRUE)
  if (barren_100_count > 0) cat(sprintf("[NDVI BOOTSTRAP] Skipped normalization for %d location-year(s) with 100%% barren\n", barren_100_count))

  veg_types <- unique(merged_veg$Veg[!tolower(merged_veg$Veg) %in% c("barren")])
  years <- sort(unique(merged_veg$pheno_year[!is.na(merged_veg$pheno_year)]))
  locations <- unique(merged_veg$location_id)
  n_locs <- length(locations)
  if (length(years) == 0 || n_locs == 0 || nrow(merged_veg) == 0) {
    warning("[NDVI BOOTSTRAP] No valid data after filtering - skipping bootstrap")
    return(NULL)  # Return NULL so caller knows aggregation failed
  }
  if (!is.numeric(B) || length(B) != 1 || B < 1) { warning("[NDVI BOOTSTRAP] Invalid B parameter - skipping bootstrap"); return(NULL) }
  B <- as.integer(B)

  unique_loc_years <- merged_veg |> dplyr::distinct(location_id, pheno_year, ndvi_row_idx)

  veg_boot_res <- list()
  all_veg_plus_barren <- c(veg_types, "Barren")
  for (v in all_veg_plus_barren) {
    mat <- matrix(NA_real_, nrow = as.integer(B), ncol = as.integer(length(years)))
    if (is.matrix(mat) && ncol(mat) > 0) colnames(mat) <- as.character(years)
    veg_boot_res[[v]] <- mat
  }

  for (b in seq_len(B)) {
    boot_locs <- sample(locations, n_locs, replace = TRUE)
    raw_ndvi_vals <- ndvi_boot_mat[, b]
    norm_ndvi_vals <- pmin(pmax(raw_ndvi_vals / ndvi_max, 0), 1)
    loc_counts <- table(boot_locs)
    sampled_locs <- names(loc_counts)
    sub_veg <- merged_veg[merged_veg$location_id %in% sampled_locs, ]
    if (nrow(sub_veg) > 0) {
      sub_veg$current_ndvi_norm <- norm_ndvi_vals[sub_veg$ndvi_row_idx]
      sub_veg$abs_cover <- sub_veg$rel_coef * sub_veg$current_ndvi_norm
      sub_veg$weight <- as.integer(loc_counts[sub_veg$location_id])
      agg_veg <- sub_veg |> dplyr::group_by(pheno_year, Veg) |> dplyr::summarize(total_cover = sum(abs_cover * weight, na.rm = TRUE), .groups = "drop")
      for (i in 1:nrow(agg_veg)) {
        v <- agg_veg$Veg[i]; y <- as.character(agg_veg$pheno_year[i])
        if (v %in% veg_types && y %in% colnames(veg_boot_res[[v]])) veg_boot_res[[v]][b, y] <- agg_veg$total_cover[i] / n_locs
      }
    }
    sub_barren <- unique_loc_years[unique_loc_years$location_id %in% sampled_locs, ]
    if (nrow(sub_barren) > 0) {
      sub_barren$current_ndvi_norm <- norm_ndvi_vals[sub_barren$ndvi_row_idx]
      sub_barren$barren_frac <- 1.0 - sub_barren$current_ndvi_norm
      sub_barren$weight <- as.integer(loc_counts[sub_barren$location_id])
      agg_barren <- sub_barren |> dplyr::group_by(pheno_year) |> dplyr::summarize(total_barren = sum(barren_frac * weight, na.rm = TRUE), n_obs_year = sum(weight, na.rm = TRUE), .groups = "drop")
      for (i in 1:nrow(agg_barren)) {
        y <- as.character(agg_barren$pheno_year[i]); if (y %in% colnames(veg_boot_res[["Barren"]])) veg_boot_res[["Barren"]][b, y] <- agg_barren$total_barren[i] / agg_barren$n_obs_year[i]
      }
    }
  }

  # --- Compute "woody" (populus + tamarix + woody_unknown) aggregation with proper CI propagation ---
  woody_types <- c("populus", "tamarix", "woody_unknown")
  woody_mats <- veg_boot_res[tolower(names(veg_boot_res)) %in% woody_types]
  if (length(woody_mats) >= 1) {
    woody_mat <- Reduce(`+`, lapply(woody_mats, function(m) { m[is.na(m)] <- 0; m }))
    all_na_mask <- Reduce(`&`, lapply(woody_mats, is.na))
    woody_mat[all_na_mask] <- NA_real_
    colnames(woody_mat) <- colnames(woody_mats[[1]])
    veg_boot_res[["woody"]] <- woody_mat
    cat("[NDVI BOOTSTRAP] Added 'woody' category (populus + tamarix + woody_unknown) with combined bootstrap CIs\n")
  }

  final_results <- list()
  for (v in names(veg_boot_res)) {
    mat <- veg_boot_res[[v]]
    if (is.null(mat) || ncol(mat) == 0 || is.null(colnames(mat)) || length(colnames(mat)) == 0) next
    df_res <- data.frame(year = as.integer(colnames(mat)), Veg = v, global_coef = apply(mat, 2, mean, na.rm = TRUE), se = apply(mat, 2, sd, na.rm = TRUE), coef_025 = apply(mat, 2, quantile, 0.025, na.rm = TRUE), coef_975 = apply(mat, 2, quantile, 0.975, na.rm = TRUE), method = "location_bootstrap_ndvi")
    n_locs_per_year <- sapply(years, function(y) sum(unique_loc_years$pheno_year == y))
    df_res$n_locations <- n_locs_per_year[match(df_res$year, years)]
    final_results[[v]] <- df_res
  }
  dplyr::bind_rows(final_results)
}

# --- NEW: MSAVI-BASED LOCATION BOOTSTRAP (Analogous to PPI) ---
# Uses per-location summer MSAVI medians and bootstraps medians like PPI.
location_bootstrap_msavi <- function(all_coefs, df_tasks, B = BOOTSTRAP_B, seed = 123, msavi_max = 0.6) {
  set.seed(seed)
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    warning("dplyr required for MSAVI bootstrap")
    return(NULL)
  }

  msavi_col <- "MSAVI"
  if ("MSAVI_raw" %in% names(df_tasks) && any(!is.na(df_tasks$MSAVI_raw))) {
    msavi_col <- "MSAVI_raw"
    cat("[MSAVI BOOTSTRAP] Using raw MSAVI values (MSAVI_raw) for weighting.\n")
  } else if ("MSAVI" %in% names(df_tasks)) {
    cat("[MSAVI BOOTSTRAP] Using standard MSAVI column.\n")
  } else {
    warning("MSAVI column missing from df_tasks")
    return(NULL)
  }

  # Prepare summer detrended MSAVI similar to PPI pipeline
  if (!"month" %in% names(df_tasks)) {
    if ("date" %in% names(df_tasks)) df_tasks$month <- as.integer(format(as.Date(df_tasks$date), "%m")) else df_tasks$month <- 1L
  }
  if (!"doy" %in% names(df_tasks)) {
    if ("date" %in% names(df_tasks)) df_tasks$doy <- as.numeric(format(as.Date(df_tasks$date), "%j")) else df_tasks$doy <- 150
  }

  summer_df <- df_tasks |> dplyr::filter(month %in% 6:9, get(msavi_col) > -Inf)
  if (nrow(summer_df) > 50) {
    tryCatch({
      f <- as.formula(paste(msavi_col, "~ poly(doy, 3)"))
      seasonal_model <- lm(f, data = summer_df)
      summer_df$seasonal_trend <- predict(seasonal_model, newdata = summer_df)
      global_seasonal_mean <- mean(summer_df$seasonal_trend, na.rm = TRUE)
      summer_df$msavi_detrended <- summer_df[[msavi_col]] - (summer_df$seasonal_trend - global_seasonal_mean)
      summer_df$msavi_detrended <- as.numeric(summer_df$msavi_detrended)
      cat(sprintf("[MSAVI BOOTSTRAP] Computed detrended summer MSAVI (N=%d).\n", nrow(summer_df)))
    }, error = function(e) {
      cat(sprintf("[MSAVI BOOTSTRAP] Detrending failed: %s. Using raw MSAVI.\n", e$message))
      summer_df$msavi_detrended <- as.numeric(summer_df[[msavi_col]])
    })
  } else {
    summer_df$msavi_detrended <- as.numeric(summer_df[[msavi_col]])
    cat("[MSAVI BOOTSTRAP] Insufficient summer data for detrending; using raw values.\n")
  }

  # Safe summarise: ensure we only boot groups with finite detrended values and provide diagnostics on failure
  # Build per-group list via split+lapply to avoid dplyr evaluation scoping issues
  if (nrow(summer_df) == 0) {
    boot_combined <- data.frame(location_id = character(0), pheno_year = integer(0), msavi_boot_list = I(list()))
  } else {
    # Use interaction with a unique separator to avoid issues with dots in location_id
    summer_df$grp_key <- paste(summer_df$location_id, summer_df$pheno_year, sep = "|||")
    grp <- split(summer_df, summer_df$grp_key, drop = TRUE)
    grp_names <- names(grp)
    msavi_lists <- vector("list", length(grp))
    for (gi in seq_along(grp)) {
      g <- grp[[gi]]
      vec <- as.numeric(g$msavi_detrended)
      if (sum(is.finite(vec)) > 0) {
        msavi_lists[[gi]] <- boot_median_vec(vec, B)
      } else {
        msavi_lists[[gi]] <- rep(NA_real_, as.integer(B))
      }
    }
    # Parse names back to location_id / pheno_year using the unique separator
    parsed <- do.call(rbind, strsplit(grp_names, "\\|\\|\\|"))
    boot_combined <- data.frame(location_id = parsed[,1], pheno_year = as.integer(parsed[,2]), msavi_boot_list = I(msavi_lists), stringsAsFactors = FALSE)
  }
  boot_combined$key <- paste(boot_combined$location_id, boot_combined$pheno_year, sep = "_")
  n_rows <- nrow(boot_combined)
  msavi_boot_mat <- matrix(NA_real_, nrow = n_rows, ncol = as.integer(B))
  for (i in 1:n_rows) {
    s_vec <- boot_combined$msavi_boot_list[[i]]
    if (!is.null(s_vec) && !all(is.na(s_vec))) msavi_boot_mat[i, ] <- s_vec
  }
  key_map <- setNames(seq_len(n_rows), boot_combined$key)

  # Prepare merged veg data like PPI (exclude barren rows)
  veg_char <- as.character(all_coefs$Veg)
  valid_veg_mask <- !is.na(veg_char) & nchar(trimws(veg_char)) > 0 & trimws(veg_char) != "NA" & !tolower(trimws(veg_char)) %in% c("barren")
  merged_veg <- all_coefs[valid_veg_mask, ]
  merged_veg$key <- paste(merged_veg$location_id, merged_veg$pheno_year, sep = "_")
  merged_veg$msavi_row_idx <- key_map[merged_veg$key]

  # Debug: show key matching statistics
  n_coef_keys <- length(unique(merged_veg$key))
  n_msavi_keys <- length(names(key_map))
  n_matched <- sum(!is.na(merged_veg$msavi_row_idx))
  cat(sprintf("[MSAVI BOOTSTRAP DEBUG] Coef keys: %d unique, MSAVI keys: %d, Matched: %d/%d rows\n",
              n_coef_keys, n_msavi_keys, n_matched, nrow(merged_veg)))
  if (n_matched == 0 && n_coef_keys > 0 && n_msavi_keys > 0) {
    cat(sprintf("[MSAVI BOOTSTRAP DEBUG] Sample coef keys: %s\n", paste(head(unique(merged_veg$key), 3), collapse=", ")))
    cat(sprintf("[MSAVI BOOTSTRAP DEBUG] Sample MSAVI keys: %s\n", paste(head(names(key_map), 3), collapse=", ")))
  }

  merged_veg <- merged_veg[!is.na(merged_veg$msavi_row_idx), ]

  # Calculate relative coefficients (fractions of vegetation) - same as PPI bootstrap
  # IMPORTANT: Skip normalization for 100% barren cases to preserve relative proportions
  merged_veg <- merged_veg |>
    dplyr::group_by(location_id, pheno_year) |>
    dplyr::mutate(
      sum_veg_coef = sum(coef[tolower(Veg) != "barren"], na.rm = TRUE),
      # Only normalize if there's actual vegetation (sum > threshold)
      # If 100% barren (sum_veg_coef ~ 0), keep original coefficients as rel_coef (preserves relative proportions)
      rel_coef = ifelse(sum_veg_coef > 1e-6, coef / sum_veg_coef, coef)
    ) |>
    dplyr::ungroup() |>
    # Renormalize to sum-to-one ONLY if the sum deviates from 1 AND there's actual vegetation
    # This preserves relative proportions - we scale all values by the same factor
    dplyr::group_by(location_id, pheno_year) |>
    dplyr::mutate(
      total_rel = sum(rel_coef, na.rm = TRUE),
      # Only renormalize if: (1) sum is not already ~1, (2) sum is > threshold (not 100% barren)
      rel_coef = ifelse(total_rel > 1e-6 & abs(total_rel - 1) > 1e-9,
                        rel_coef / total_rel,
                        rel_coef)
    ) |>
    dplyr::select(-total_rel) |>
    dplyr::ungroup()

  # Diagnostic: report any location-years that required correction
  rels_check <- merged_veg |>
    dplyr::group_by(location_id, pheno_year) |>
    dplyr::summarize(total_rel = sum(rel_coef, na.rm = TRUE), .groups = "drop")
  bad_rel <- sum(abs(rels_check$total_rel - 1) > 1e-6 & rels_check$total_rel > 1e-6, na.rm = TRUE)
  if (bad_rel > 0) cat(sprintf("[NORMALIZE COEF] Corrected rel_coef for %d location-year(s) to sum-to-one\n", bad_rel))

  # Report 100% barren cases that were skipped
  barren_100_count <- sum(rels_check$total_rel < 1e-6, na.rm = TRUE)
  if (barren_100_count > 0) cat(sprintf("[MSAVI BOOTSTRAP] Skipped normalization for %d location-year(s) with 100%% barren\n", barren_100_count))

  veg_types <- unique(merged_veg$Veg[!tolower(merged_veg$Veg) %in% c("barren")])
  years <- sort(unique(merged_veg$pheno_year[!is.na(merged_veg$pheno_year)]))
  locations <- unique(merged_veg$location_id)
  n_locs <- length(locations)
  if (length(years) == 0 || n_locs == 0 || nrow(merged_veg) == 0) {
    warning("[MSAVI BOOTSTRAP] No valid data after filtering - skipping bootstrap")
    return(NULL)  # Return NULL so caller knows aggregation failed
  }
  if (!is.numeric(B) || length(B) != 1 || B < 1) { warning("[MSAVI BOOTSTRAP] Invalid B parameter - skipping bootstrap"); return(NULL) }
  B <- as.integer(B)

  unique_loc_years <- merged_veg |> dplyr::distinct(location_id, pheno_year, msavi_row_idx)

  veg_boot_res <- list()
  all_veg_plus_barren <- c(veg_types, "Barren")
  for (v in all_veg_plus_barren) {
    mat <- matrix(NA_real_, nrow = as.integer(B), ncol = as.integer(length(years)))
    if (is.matrix(mat) && ncol(mat) > 0) colnames(mat) <- as.character(years)
    veg_boot_res[[v]] <- mat
  }

  for (b in seq_len(B)) {
    boot_locs <- sample(locations, n_locs, replace = TRUE)
    raw_msavi_vals <- msavi_boot_mat[, b]
    norm_msavi_vals <- pmin(pmax(raw_msavi_vals / msavi_max, 0), 1)
    loc_counts <- table(boot_locs)
    sampled_locs <- names(loc_counts)
    sub_veg <- merged_veg[merged_veg$location_id %in% sampled_locs, ]
    if (nrow(sub_veg) > 0) {
      sub_veg$current_msavi_norm <- norm_msavi_vals[sub_veg$msavi_row_idx]
      sub_veg$abs_cover <- sub_veg$rel_coef * sub_veg$current_msavi_norm
      sub_veg$weight <- as.integer(loc_counts[sub_veg$location_id])
      agg_veg <- sub_veg |> dplyr::group_by(pheno_year, Veg) |> dplyr::summarize(total_cover = sum(abs_cover * weight, na.rm = TRUE), .groups = "drop")
      for (i in 1:nrow(agg_veg)) {
        v <- agg_veg$Veg[i]; y <- as.character(agg_veg$pheno_year[i])
        if (v %in% veg_types && y %in% colnames(veg_boot_res[[v]])) veg_boot_res[[v]][b, y] <- agg_veg$total_cover[i] / n_locs
      }
    }
    sub_barren <- unique_loc_years[unique_loc_years$location_id %in% sampled_locs, ]
    if (nrow(sub_barren) > 0) {
      sub_barren$current_msavi_norm <- norm_msavi_vals[sub_barren$msavi_row_idx]
      sub_barren$barren_frac <- 1.0 - sub_barren$current_msavi_norm
      sub_barren$weight <- as.integer(loc_counts[sub_barren$location_id])
      agg_barren <- sub_barren |> dplyr::group_by(pheno_year) |> dplyr::summarize(total_barren = sum(barren_frac * weight, na.rm = TRUE), n_obs_year = sum(weight, na.rm = TRUE), .groups = "drop")
      for (i in 1:nrow(agg_barren)) {
        y <- as.character(agg_barren$pheno_year[i]); if (y %in% colnames(veg_boot_res[["Barren"]])) veg_boot_res[["Barren"]][b, y] <- agg_barren$total_barren[i] / agg_barren$n_obs_year[i]
      }
    }
  }

  # --- Compute "woody" (populus + tamarix + woody_unknown) aggregation with proper CI propagation ---
  woody_types <- c("populus", "tamarix", "woody_unknown")
  woody_mats <- veg_boot_res[tolower(names(veg_boot_res)) %in% woody_types]
  if (length(woody_mats) >= 1) {
    woody_mat <- Reduce(`+`, lapply(woody_mats, function(m) { m[is.na(m)] <- 0; m }))
    all_na_mask <- Reduce(`&`, lapply(woody_mats, is.na))
    woody_mat[all_na_mask] <- NA_real_
    colnames(woody_mat) <- colnames(woody_mats[[1]])
    veg_boot_res[["woody"]] <- woody_mat
    cat("[MSAVI BOOTSTRAP] Added 'woody' category (populus + tamarix + woody_unknown) with combined bootstrap CIs\n")
  }

  final_results <- list()
  for (v in names(veg_boot_res)) {
    mat <- veg_boot_res[[v]]
    if (is.null(mat) || ncol(mat) == 0 || is.null(colnames(mat)) || length(colnames(mat)) == 0) next
    df_res <- data.frame(year = as.integer(colnames(mat)), Veg = v, global_coef = apply(mat, 2, mean, na.rm = TRUE), se = apply(mat, 2, sd, na.rm = TRUE), coef_025 = apply(mat, 2, quantile, 0.025, na.rm = TRUE), coef_975 = apply(mat, 2, quantile, 0.975, na.rm = TRUE), method = "location_bootstrap_msavi")
    n_locs_per_year <- sapply(years, function(y) sum(unique_loc_years$pheno_year == y))
    df_res$n_locations <- n_locs_per_year[match(df_res$year, years)]
    final_results[[v]] <- df_res
  }
  dplyr::bind_rows(final_results)
}

# --- NEW: No-Index (unscaled) location bootstrap ---
# Uses raw relative coefficients as absolute cover (no scaling by an index)
location_bootstrap_noindex <- function(all_coefs, df_tasks = NULL, B = BOOTSTRAP_B, seed = 123) {
  set.seed(seed)
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    warning("dplyr required for NoIndex bootstrap")
    return(NULL)
  }

  veg_char <- as.character(all_coefs$Veg)
  valid_veg_mask <- !is.na(veg_char) & nchar(trimws(veg_char)) > 0 & trimws(veg_char) != "NA" & !tolower(trimws(veg_char)) %in% c("barren")
  merged_veg <- all_coefs[valid_veg_mask, ]
  merged_veg$rel_coef <- merged_veg$coef

  veg_types <- unique(merged_veg$Veg[!tolower(merged_veg$Veg) %in% c("barren")])
  years <- sort(unique(merged_veg$pheno_year[!is.na(merged_veg$pheno_year)]))
  locations <- unique(merged_veg$location_id)
  n_locs <- length(locations)
  if (length(years) == 0 || n_locs == 0 || nrow(merged_veg) == 0) { warning("[NOINDEX BOOTSTRAP] No valid data after filtering - skipping bootstrap"); return(all_coefs) }
  if (!is.numeric(B) || length(B) != 1 || B < 1) { warning("[NOINDEX BOOTSTRAP] Invalid B parameter - skipping bootstrap"); return(all_coefs) }
  B <- as.integer(B)

  unique_loc_years <- merged_veg |> dplyr::distinct(location_id, pheno_year)

  veg_boot_res <- list()
  all_veg_plus_barren <- c(veg_types, "barren")
  for (v in all_veg_plus_barren) {
    mat <- matrix(NA_real_, nrow = as.integer(B), ncol = as.integer(length(years)))
    if (is.matrix(mat) && ncol(mat) > 0) colnames(mat) <- as.character(years)
    veg_boot_res[[v]] <- mat
  }

  for (b in seq_len(B)) {
    boot_locs <- sample(locations, n_locs, replace = TRUE)
    loc_counts <- table(boot_locs)
    sampled_locs <- names(loc_counts)
    sub_veg <- merged_veg[merged_veg$location_id %in% sampled_locs, ]
    if (nrow(sub_veg) > 0) {
      sub_veg$abs_cover <- sub_veg$rel_coef
      sub_veg$weight <- as.integer(loc_counts[sub_veg$location_id])
      agg_veg <- sub_veg |> dplyr::group_by(pheno_year, Veg) |> dplyr::summarize(total_cover = sum(abs_cover * weight, na.rm = TRUE), .groups = "drop")
      for (i in 1:nrow(agg_veg)) {
        v <- agg_veg$Veg[i]; y <- as.character(agg_veg$pheno_year[i])
        if (v %in% veg_types && y %in% colnames(veg_boot_res[[v]])) veg_boot_res[[v]][b, y] <- agg_veg$total_cover[i] / n_locs
      }
    }

    # Barren = 1 - sum(veg fractions) per year, averaged across locations
    # First compute total vegetation cover per location-year, then average across locations per year
    if (nrow(sub_veg) > 0) {
      # Sum vegetation coefficients per location-year
      loc_year_sums <- sub_veg |>
        dplyr::group_by(location_id, pheno_year) |>
        dplyr::summarize(total_veg = sum(rel_coef, na.rm = TRUE), .groups = "drop")
      loc_year_sums$weight <- as.integer(loc_counts[as.character(loc_year_sums$location_id)])

      # Compute barren = 1 - total_veg, clamped to [0, 1]
      loc_year_sums$barren <- pmax(0, pmin(1, 1.0 - loc_year_sums$total_veg))

      # Average barren across locations per year (weighted by bootstrap counts)
      barren_by_year <- loc_year_sums |>
        dplyr::group_by(pheno_year) |>
        dplyr::summarize(mean_barren = sum(barren * weight, na.rm = TRUE) / sum(weight, na.rm = TRUE), .groups = "drop")

      for (i in 1:nrow(barren_by_year)) {
        y <- as.character(barren_by_year$pheno_year[i])
        if (y %in% colnames(veg_boot_res[["barren"]])) {
          veg_boot_res[["barren"]][b, y] <- barren_by_year$mean_barren[i]
        }
      }
    }
  }

  # --- Compute "woody" (populus + tamarix + woody_unknown) aggregation with proper CI propagation ---
  woody_types <- c("populus", "tamarix", "woody_unknown")
  woody_mats <- veg_boot_res[tolower(names(veg_boot_res)) %in% woody_types]
  if (length(woody_mats) >= 1) {
    woody_mat <- Reduce(`+`, lapply(woody_mats, function(m) { m[is.na(m)] <- 0; m }))
    all_na_mask <- Reduce(`&`, lapply(woody_mats, is.na))
    woody_mat[all_na_mask] <- NA_real_
    colnames(woody_mat) <- colnames(woody_mats[[1]])
    veg_boot_res[["woody"]] <- woody_mat
    cat("[NOINDEX BOOTSTRAP] Added 'woody' category (populus + tamarix + woody_unknown) with combined bootstrap CIs\n")
  }

  final_results <- list()
  for (v in names(veg_boot_res)) {
    mat <- veg_boot_res[[v]]
    if (is.null(mat) || ncol(mat) == 0 || is.null(colnames(mat)) || length(colnames(mat)) == 0) next
    df_res <- data.frame(year = as.integer(colnames(mat)), Veg = v, global_coef = apply(mat, 2, mean, na.rm = TRUE), se = apply(mat, 2, sd, na.rm = TRUE), coef_025 = apply(mat, 2, quantile, 0.025, na.rm = TRUE), coef_975 = apply(mat, 2, quantile, 0.975, na.rm = TRUE), method = "location_bootstrap_noindex")
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

  # =========================================================================
  # CLASSIFICATION UNCERTAINTY: Apply Dirichlet perturbation if enabled
 # This perturbs the coefficient vectors based on the confusion matrix
  # to propagate classification/misclassification uncertainty
  # =========================================================================
  has_classification_uncertainty <- isTRUE(ENABLE_CLASSIFICATION_UNCERTAINTY) &&
                                    exists(".CONFUSION_MATRIX", envir = globalenv())

  if (has_classification_uncertainty) {
    conf_matrix <- get(".CONFUSION_MATRIX", envir = globalenv())
    sample_sizes <- if (exists(".VALIDATION_SAMPLE_SIZES", envir = globalenv())) {
      get(".VALIDATION_SAMPLE_SIZES", envir = globalenv())
    } else NULL

    message("[BOOTSTRAP] Classification uncertainty enabled; applying Dirichlet perturbation.")

    # Get vegetation classes from confusion matrix
    conf_classes <- tolower(rownames(conf_matrix))

    # For each unique location-year, perturb the coefficient vector
    # We need to work with wide-format data (one row per location-year, columns for each veg)
    loc_yr_keys <- unique(paste(all_coefs$location_id, all_coefs$pheno_year, sep = "___"))

    perturbed_coefs <- all_coefs

    for (key in loc_yr_keys) {
      parts <- strsplit(key, "___")[[1]]
      loc_id <- parts[1]
      p_year <- as.numeric(parts[2])

      # Get rows for this location-year
      mask <- all_coefs$location_id == loc_id & all_coefs$pheno_year == p_year
      loc_yr_data <- all_coefs[mask, ]

      if (nrow(loc_yr_data) == 0) next

      # Build the coefficient vector for this location-year
      veg_classes <- as.character(loc_yr_data$Veg)
      coef_vals <- loc_yr_data$coef

      # Create named vector
      frac_vec <- setNames(coef_vals, tolower(veg_classes))

      # Find the dominant class (highest coefficient) as proxy for "true class"
      # In practice, training locations have known true class, but for inference
      # locations we use the dominant predicted class
      dominant_class <- names(which.max(frac_vec))

      # Check if dominant class is in confusion matrix
      if (!is.null(dominant_class) && dominant_class %in% conf_classes) {
        # Apply Dirichlet perturbation
        perturbed <- apply_dirichlet_perturbation(
          fractions = frac_vec,
          true_class = dominant_class,
          conf_matrix = conf_matrix,
          sample_sizes = sample_sizes,
          concentration_scale = DIRICHLET_CONCENTRATION_SCALE
        )

        # Map perturbed values back to original rows
        for (i in seq_len(nrow(loc_yr_data))) {
          veg_lower <- tolower(as.character(loc_yr_data$Veg[i]))
          if (veg_lower %in% names(perturbed)) {
            row_idx <- which(mask)[i]
            perturbed_coefs$coef[row_idx] <- perturbed[veg_lower]
          }
        }
      }
    }

    all_coefs <- perturbed_coefs
    message(sprintf("[BOOTSTRAP] Dirichlet perturbation applied to %d location-years", length(loc_yr_keys)))
  } else if (isTRUE(ENABLE_CLASSIFICATION_UNCERTAINTY)) {
    message("[BOOTSTRAP] Classification uncertainty enabled but no confusion matrix available; skipping Dirichlet perturbation.")
  }

  veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
  years <- sort(unique(all_coefs$pheno_year[!is.na(all_coefs$pheno_year)]))
  results_list <- list()

  # Guard against empty data
  if (length(years) == 0 || length(veg_types) == 0) {
    warning("[BOOTSTRAP] No valid years or veg types - returning empty results")
    return(data.frame())
  }

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
    if (length(years) > 0) {
      colnames(boot_means) <- as.character(years)
    }
    
    n_eff_vec <- numeric(length(years))

    for (i in seq_along(years)) {
      yr <- years[i]
      yr_data <- veg_data[veg_data$pheno_year == yr, ]
      n_obs <- nrow(yr_data)

      # Estimate effective sample size using pairwise spatial correlation
      # Default: n_eff = n_obs (assume independence if no spatial info)
      n_eff_est <- n_obs
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
             missing_idx <- which(is.na(lon) | is.na(lat))
             m_ids <- as.character(yr_data$location_id[missing_idx])
             match_idx <- match(m_ids, as.character(gpts$location_id))
             lon[missing_idx] <- ifelse(is.na(match_idx), lon[missing_idx], gpts$lon[match_idx])
             lat[missing_idx] <- ifelse(is.na(match_idx), lat[missing_idx], gpts$lat[match_idx])
           }
        }

        if (all(is.na(lon)) || all(is.na(lat))) {
          n_eff_est <- n_obs  # No spatial info: assume independence
        } else {
          # Compute pairwise great-circle distances (km) using Haversine
          valid <- which(!is.na(lon) & !is.na(lat))
          if (length(valid) >= 3) {
            coords <- cbind(lon[valid], lat[valid])
            coefs_valid <- yr_data$coef[valid]

            # Distance matrix in km (Haversine approximation)
            rad <- pi / 180
            n_v <- length(valid)
            dist_mat <- matrix(0, n_v, n_v)
            for (ii in 1:(n_v - 1)) {
              for (jj in (ii + 1):n_v) {
                dlat <- (coords[jj, 2] - coords[ii, 2]) * rad
                dlon <- (coords[jj, 1] - coords[ii, 1]) * rad
                a <- sin(dlat / 2)^2 + cos(coords[ii, 2] * rad) * cos(coords[jj, 2] * rad) * sin(dlon / 2)^2
                dist_mat[ii, jj] <- 2 * 6371 * asin(sqrt(a))
                dist_mat[jj, ii] <- dist_mat[ii, jj]
              }
            }

            # Estimate spatial correlation range via empirical variogram
            # Fit exponential model: C(d) = exp(-d / range)
            # Use Moran-style approach: mean pairwise correlation as function of distance
            dists <- dist_mat[upper.tri(dist_mat)]
            coef_diffs_sq <- outer(coefs_valid, coefs_valid, function(a, b) (a - b)^2)
            gamma_vals <- coef_diffs_sq[upper.tri(coef_diffs_sq)] / 2  # semivariance

            if (length(dists) > 0 && var(coefs_valid, na.rm = TRUE) > 0) {
              total_var <- var(coefs_valid, na.rm = TRUE)
              # Estimate range by finding distance at which semivariance reaches ~63% of sill
              # Simple robust estimate: fit exponential variogram gamma(d) = sill * (1 - exp(-d/range))
              # Use median-based bins for robustness
              n_bins <- min(10, max(3, length(dists) %/% 5))
              bin_breaks <- quantile(dists, probs = seq(0, 1, length.out = n_bins + 1))
              bin_breaks <- unique(bin_breaks)
              if (length(bin_breaks) >= 3) {
                bin_mid <- (bin_breaks[-length(bin_breaks)] + bin_breaks[-1]) / 2
                bin_gamma <- numeric(length(bin_mid))
                for (bb in seq_along(bin_mid)) {
                  in_bin <- dists >= bin_breaks[bb] & dists < bin_breaks[bb + 1]
                  if (sum(in_bin) > 0) bin_gamma[bb] <- median(gamma_vals[in_bin])
                  else bin_gamma[bb] <- NA
                }
                valid_bins <- !is.na(bin_gamma)
                if (sum(valid_bins) >= 2) {
                  # Fit exponential variogram: gamma(d) = sill * (1 - exp(-d/range))
                  # Using NLS with reasonable starting values
                  fit_ok <- tryCatch({
                    nls_fit <- nls(g ~ s * (1 - exp(-d / r)),
                                   data = data.frame(d = bin_mid[valid_bins], g = bin_gamma[valid_bins]),
                                   start = list(s = total_var, r = median(dists)),
                                   lower = list(s = total_var * 0.1, r = max(dists) * 0.01),
                                   upper = list(s = total_var * 3, r = max(dists) * 2),
                                   algorithm = "port",
                                   control = list(maxiter = 50, warnOnly = TRUE))
                    range_est <- coef(nls_fit)["r"]
                    TRUE
                  }, error = function(e) FALSE)

                  if (!fit_ok) {
                    # Fallback: use distance at which semivariance first exceeds 0.5 * total_var
                    thresh_idx <- which(bin_gamma[valid_bins] >= 0.5 * total_var)
                    if (length(thresh_idx) > 0) {
                      range_est <- bin_mid[valid_bins][thresh_idx[1]]
                    } else {
                      range_est <- max(dists)  # All correlated -> conservative
                    }
                  }

                  # Compute mean pairwise correlation: C(d) = exp(-d / range)
                  # Effective n = n / (1 + (n-1) * mean_corr)  [Kish formula]
                  mean_corr <- mean(exp(-dists / range_est))
                  n_eff_est <- max(1, n_obs / (1 + (n_obs - 1) * mean_corr))
                } else {
                  n_eff_est <- n_obs
                }
              } else {
                n_eff_est <- n_obs
              }
            } else {
              n_eff_est <- n_obs
            }
          } else {
            n_eff_est <- n_obs
          }
        }
      }
      n_eff_vec[i] <- n_eff_est

      if (n_obs > 0) {
        # DECISION: Small Sample Size vs Large Sample Size
        if (n_obs < 15) {
             # --- SMALL SAMPLE: USE POOLED VARIANCE ---
             # When N is small, standard bootstrap under-estimates variance (often 0 if N=1).
             # We assume the spatial variability is 'pooled_sd' and the mean is the sample mean.
             
             mu <- mean(yr_data$coef, na.rm = TRUE)
             
             # Use effective sample size (accounts for spatial autocorrelation)
             n_eff_i <- n_eff_vec[i]

             # Use sample SD if available (n_obs >= 2), otherwise pooled SD
             if (n_obs >= 2) {
               sample_sd <- sd(yr_data$coef, na.rm = TRUE)
               if (is.na(sample_sd) || sample_sd == 0) sample_sd <- pooled_sd
               se_mean <- sample_sd / sqrt(n_eff_i)
             } else {
               se_mean <- pooled_sd / sqrt(n_eff_i)
             }

             # If we have measurement uncertainty, we can add that too (in quadrature)
             if (has_meas_uncertainty) {
               # Average measurement error for this year
               mean_meas_sd <- mean(yr_data$coef_sd, na.rm = TRUE)
               if (is.na(mean_meas_sd)) mean_meas_sd <- 0
               # Effective SE combines spatial uncertainty (SE_mean) and measurement uncertainty
               se_meas_mean <- mean_meas_sd / sqrt(n_eff_i)
               se_total <- sqrt(se_mean^2 + se_meas_mean^2)
             } else {
               se_total <- se_mean
             }

             # Moderate inflation for very small effective N
             ci_tmp <- small_n_inflated_ci(est = mu, sd_in = se_total, n_obs = n_eff_i)
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
                na_warning_logged <- FALSE  # Flag to log NA warning only once per year/veg
                rnorm_na_warning_logged <- FALSE  # Flag for rnorm NA warning
                for (b in 1:B) {
                   # Resample indices (capture spatial variance)
                   idx <- sample(seq_len(n_obs), n_obs, replace = TRUE)
                   sel_coefs <- yr_data$coef[idx]
                   sel_sds <- yr_data$coef_sd[idx]

                   # Defensive cleaning: replace NA/Inf values with mean of finite values
                   if (all(is.na(sel_coefs))) {
                     next  # Skip this iteration instead of setting to 0
                   }
                   sel_coefs[is.na(sel_coefs) | !is.finite(sel_coefs)] <- mean(sel_coefs, na.rm = TRUE)

                   sel_sds[is.na(sel_sds) | !is.finite(sel_sds) | sel_sds < 0] <- 0

                   # Add measurement noise (capture measurement variance)
                   sim_coefs <- rnorm(n_obs, mean = sel_coefs, sd = sel_sds)

                   if (any(is.na(sim_coefs))) {
                     if (!rnorm_na_warning_logged) {
                       write_debug(sprintf("Year %d veg=%s: rnorm produced %d NAs (n_obs=%d); sel_coefs finite=%d sel_sds finite=%d (suppressing further warnings)", years[i], veg, sum(is.na(sim_coefs)), n_obs, sum(is.finite(sel_coefs)), sum(is.finite(sel_sds))))
                       rnorm_na_warning_logged <- TRUE
                     }
                     # Replace NAs with mean to avoid propagating NA further
                     sim_coefs[is.na(sim_coefs)] <- mean(sim_coefs, na.rm = TRUE)
                   }

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
      global_coef = apply(boot_means, 2, mean, na.rm = TRUE),
      se = apply(boot_means, 2, sd, na.rm = TRUE),
      coef_025 = apply(boot_means, 2, quantile, 0.025, na.rm = TRUE),
      coef_975 = apply(boot_means, 2, quantile, 0.975, na.rm = TRUE),
      n_eff = n_eff_vec,
      method = "robust_location_bootstrap"
    )

    # Print spatial effective sample size diagnostics
    cat(sprintf("[BOOTSTRAP] Veg '%s': effective sample size (n_eff) across years — mean=%.1f, median=%.1f, min=%.1f, max=%.1f (of n_obs mean=%.1f)\n",
        veg,
        mean(n_eff_vec, na.rm = TRUE), median(n_eff_vec, na.rm = TRUE),
        min(n_eff_vec, na.rm = TRUE), max(n_eff_vec, na.rm = TRUE),
        mean(sapply(years, function(y) sum(veg_data$pheno_year == y & !is.na(veg_data$coef)))))
    )

    # NOTE: No post-hoc autocorrelation inflation applied.
    # For the large-N path (n>=15), the bootstrap already captures between-location variance.
    # For the small-N path (n<15), spatial autocorrelation is accounted for via n_eff
    # in the SE calculation above (dividing by sqrt(n_eff) instead of sqrt(n_obs)).

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
    con_out <- tryCatch(file(sink_file, open = "wt"), error = function(e) NULL)
    con_msg <- tryCatch(file(paste0(sink_file, ".msg"), open = "wt"), error = function(e) NULL)
  }

  # Monkey-patch close to be safe on NULL connections in case other code calls close(NULL)
  if (!exists("_close_patched", envir = globalenv())) {
    assign("_close_patched", TRUE, envir = globalenv())
    assign("_orig_close", base::close, envir = globalenv())
    close <- function(conn, ...) {
      if (is.null(conn)) return(invisible(NULL))
      try(.GlobalEnv$`_orig_close`(conn, ...), silent = TRUE)
    }
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

  if (isTRUE(quiet) && !is.null(con_out) && !is.null(con_msg)) {
    tryCatch({
      if (inherits(con_out, "connection") && isOpen(con_out)) {
        sink(con_out, type = "output")
      }
      if (inherits(con_msg, "connection") && isOpen(con_msg)) {
        sink(con_msg, type = "message")
      }
    }, error = function(e) {
      warning(sprintf("Failed to redirect output: %s", e$message))
    })
  }

  # Helper function to handle both spline and standard modes for uncertainty calibration
  fit_task_for_calib <- function(task_data) {
    if (is.null(task_data) || nrow(task_data) == 0) return(NULL)
    loc <- as.character(task_data$location_id[1])
    yr <- as.integer(task_data$pheno_year[1])
    
    if (isTRUE(USE_SPLINE_ENDMEMBERS) && exists("SPLINE_LIBRARY") && exists("SPLINE_PARAMS")) {
      res <- tryCatch({
        fit_one_task_spline(task_data, SPLINE_LIBRARY, SPLINE_PARAMS$indices, SPLINE_PARAMS, loc, yr)
      }, error = function(e) NULL)
    } else {
      res <- fit_one_task(task_data)
    }
    return(res)
  }

  for (gi in seq_len(nrow(grp))) {
    loc <- grp$location_id[gi]
    yr <- grp$pheno_year[gi]
    sub <- df_tasks[df_tasks$location_id == loc & df_tasks$pheno_year == yr, , drop = FALSE]
    n_full <- nrow(sub)
    if (!is.finite(n_full) || n_full < high_n) next

    # Full fit as pseudo-truth
    full_fit <- tryCatch(fit_task_for_calib(sub), error = function(e) NULL)
    full_vec <- extract_coef_vec_from_fit(full_fit, veg_levels = veg_levels)
    if (is.null(full_vec)) next

    for (n in target_ns) {
      if (n > n_full) next
      for (r in seq_len(reps)) {
        idx <- sample(seq_len(n_full), n, replace = FALSE)
        sub_n <- sub[idx, , drop = FALSE]
        fit_n <- tryCatch(fit_task_for_calib(sub_n), error = function(e) NULL)
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
  # Always reload inference data from CSV (do not use cached version)
  if (exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0) {
    cat("[INFO] Clearing cached inference data to force reload from CSV.\n")
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
    if (file.exists(INFERENCE_CSV)) {
      cat(sprintf("Loading inference data from %s...\n", INFERENCE_CSV))
      df_inf <- tryCatch(read.csv(INFERENCE_CSV), error = function(e) { cat(sprintf("[WARNING] Error reading inference CSV: %s\n", e$message)); NULL })
      if (!is.null(df_inf)) {
        cat(sprintf("Loaded %d rows from inference file.\n", nrow(df_inf)))
        
        # Treat certain classes as herbs
        veg_col_inf <- if ("Veg" %in% names(df_inf)) "Veg" else if ("vegetation" %in% names(df_inf)) "vegetation" else NULL
        if (!is.null(veg_col_inf)) {
          # Fix typos first
          df_inf[[veg_col_inf]] <- ifelse(df_inf[[veg_col_inf]] == "tamairx", "tamarix", df_inf[[veg_col_inf]])
          # Then map certain classes to herbs
          df_inf[[veg_col_inf]] <- ifelse(df_inf[[veg_col_inf]] %in% c("herbs", "alhagi", "salicornia", "halocnemum", "phragmites"), "herbs", df_inf[[veg_col_inf]])
        } else {
          cat("[WARNING] Neither 'Veg' nor 'vegetation' column found in df_inf. Skipping Veg treatment for inference data.\n")
        }
        
        # IMPORTANT: Reconstruct location_id from lat/lon coordinates
        # Do NOT use the existing location_id column (if present) as it may be incorrect
        if ("location_id" %in% names(df_inf)) {
          df_inf$location_id_orig <- df_inf$location_id  # Keep original for debugging
          df_inf$location_id <- NULL  # Remove to force reconstruction
          cat("[NOTICE] Removed existing 'location_id' column - will reconstruct from lat/lon\n")
        }
        
        # Normalize coordinate column names
        lon_candidates <- names(df_inf)[grepl("(^|_)(lon|longitude|x)(_|$)", names(df_inf), ignore.case = TRUE)]
        lat_candidates <- names(df_inf)[grepl("(^|_)(lat|latitude|y)(_|$)", names(df_inf), ignore.case = TRUE)]
        
        if (length(lon_candidates) > 0 && length(lat_candidates) > 0) {
          lon_col <- lon_candidates[1]
          lat_col <- lat_candidates[1]
          
          # Ensure numeric
          df_inf[[lon_col]] <- as.numeric(df_inf[[lon_col]])
          df_inf[[lat_col]] <- as.numeric(df_inf[[lat_col]])
          
          # Standardize to 'lon' and 'lat' column names
          if (lon_col != "lon") {
            df_inf$lon <- df_inf[[lon_col]]
            cat(sprintf("[NOTICE] Using '%s' as longitude\n", lon_col))
          }
          if (lat_col != "lat") {
            df_inf$lat <- df_inf[[lat_col]]
            cat(sprintf("[NOTICE] Using '%s' as latitude\n", lat_col))
          }
          
          # Reconstruct location_id from coordinates
          df_inf$location_id <- make_location_id(df_inf$lon, df_inf$lat)
          cat(sprintf("[NOTICE] Reconstructed location_id from lat/lon: %d unique locations\n",
                      length(unique(df_inf$location_id[!is.na(df_inf$location_id)]))))

          # Apply MAX_INFERENCE_LOCATIONS limit immediately after loading
          if (exists("MAX_INFERENCE_LOCATIONS") && is.numeric(MAX_INFERENCE_LOCATIONS)) {
            unique_locations <- unique(df_inf$location_id[!is.na(df_inf$location_id)])
            n_unique <- length(unique_locations)

            if (n_unique > MAX_INFERENCE_LOCATIONS) {
              set.seed(123)  # Deterministic sampling for reproducibility
              sampled_locations <- sample(unique_locations, MAX_INFERENCE_LOCATIONS, replace = FALSE)

              # Filter the dataframe to only include sampled locations
              df_inf <- df_inf[df_inf$location_id %in% sampled_locations, ]

              cat(sprintf("[INFERENCE LOADING] Reduced from %d to %d locations (MAX_INFERENCE_LOCATIONS=%d)\n",
                          n_unique, MAX_INFERENCE_LOCATIONS, MAX_INFERENCE_LOCATIONS))
              cat(sprintf("[INFERENCE LOADING] Filtered dataset: %d rows from %d locations\n",
                          nrow(df_inf), length(unique(df_inf$location_id[!is.na(df_inf$location_id)]))))
            } else {
              cat(sprintf("[INFERENCE LOADING] Using all %d locations (MAX_INFERENCE_LOCATIONS=%d)\n",
                          n_unique, MAX_INFERENCE_LOCATIONS))
            }
          }

          # Add Veg from gpts_map if available
          if (exists("gpts_map") && "Veg" %in% names(gpts_map)) {
            gpts_map$location_id <- as.character(gpts_map$location_id)
            df_inf$location_id <- as.character(df_inf$location_id)
            df_inf <- df_inf |>
              dplyr::left_join(
                gpts_map |> dplyr::select(location_id, Veg),
                by = "location_id"
              )
            cat(sprintf("[NOTICE] Added Veg column to inference data from gpts_map.\n"))
          }
        } else {
          stop("ERROR: Inference CSV must contain latitude and longitude columns for location_id reconstruction")
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
        if (all(c('red','nir') %in% names(df_inf))) df_inf$NDDI <- (as.numeric(df_inf$red) - as.numeric(df_inf$nir)) / (as.numeric(df_inf$red) + as.numeric(df_inf$nir) + eps)
        cat("[NOTICE] Computed NDDI for inference filtering\n")
      }

      if ("NDDI" %in% names(df_inf)) {
        dust_count <- sum(df_inf$NDDI > 0.18, na.rm = TRUE)
        total_before <- nrow(df_inf)
        df_inf <- df_inf[!(df_inf$NDDI > 0.18), , drop = FALSE]
        total_after <- nrow(df_inf)
        filtered <- total_before - total_after
        cat(sprintf("[INFERENCE FILTERING] Filtered out %d observations with dust (NDDI > 0.18) contamination\n", filtered))
        cat(sprintf("[INFERENCE FILTERING] Inference dataset after contamination filtering: %d rows from %d locations\n", total_after, length(unique(df_inf$location_id))))
        df_inf <- remove_large_outliers(df_inf)
        # compute_soil_line_slope(df_inf, assign_global_dvi = FALSE)  # Not needed for inference, soil line from training
        before_cols <- names(df_inf)
        df_inf <- compute_indices_from_bands(df_inf)
        new_cols <- setdiff(names(df_inf), before_cols)
        if (length(new_cols) > 0) cat(sprintf("[NOTICE] Computed indices from raw bands in inference data: %s\n", paste(new_cols, collapse=", ")))
      } else {
        cat("[WARNING] NDDI not found in inference data; skipping contamination filtering\n")
        # Ensure indices are computed even if contamination filtering is skipped
        before_cols <- names(df_inf)
        df_inf <- compute_indices_from_bands(df_inf)
        new_cols <- setdiff(names(df_inf), before_cols)
        if (length(new_cols) > 0) cat(sprintf("[NOTICE] Computed indices from raw bands in inference data (fallback): %s\n", paste(new_cols, collapse=", ")))
      }

      if (!"Veg" %in% names(df_inf)) df_inf$Veg <- NA_character_
      if (!"pheno_year" %in% names(df_inf)) df_inf$pheno_year <- assign_pheno_year(df_inf$date)
      if (!"doy" %in% names(df_inf)) df_inf$doy <- pheno_doy(df_inf$date)
      if (!"zenith.angle" %in% names(df_inf)) df_inf$zenith.angle <- NA_real_
      if (!"DVI_max" %in% names(df_inf)) df_inf$DVI_max <- NA_real_

      # IMPORTANT: Use ALL years for inference - inference should cover the full temporal range
      if (exists("TRAIN_YEARS") && !is.null(TRAIN_YEARS)) {
        if (!isTRUE(QUIET_MODE)) cat(sprintf("[NOTICE] TRAIN_YEARS is set to %s for training, but inference uses ALL available years.\n", paste(TRAIN_YEARS, collapse=", ")))
        cat(sprintf("[NOTICE] Inference dataset has %d rows from %d locations (all years)\n", nrow(df_inf), length(unique(df_inf$location_id))))
        cat("[NOTICE] Trend computations will use this inference dataset (INFERENCE_CSV) and will NOT use training data.\n")
      }

      # -----------------------------------------------------------------------------
      # TREND (INFERENCE DATA ONLY): compute per-pheno_year means for key indices
      # Do not use training data for trend computation — use the INFERENCE_CSV file provided
      # -----------------------------------------------------------------------------
      tryCatch({
        trend_indices <- c("MSAVI", "NDVI", "PPI")
        trend_indices <- intersect(trend_indices, names(df_inf))
        # Detrend PPI using saved seasonal model if available
        if ("PPI" %in% trend_indices && exists("INDEX_SEASONAL_MODELS") && exists("INDEX_SEASONAL_MEANS")) {
          if ("PPI" %in% names(INDEX_SEASONAL_MODELS)) {
            seasonal_model <- INDEX_SEASONAL_MODELS[["PPI"]]
            global_seasonal_mean <- INDEX_SEASONAL_MEANS[["PPI"]]
            summer_inf <- df_inf |> dplyr::filter(lubridate::month(date) %in% c(6,7,8,9))
            if (nrow(summer_inf) > 0) {
              summer_inf$.orig_row <- seq_len(nrow(df_inf))[which(lubridate::month(df_inf$date) %in% c(6,7,8,9))]
              summer_inf$doy <- lubridate::yday(summer_inf$date)
              summer_inf$seasonal_trend <- predict(seasonal_model, newdata = summer_inf)
              summer_inf$ppi_detrended <- summer_inf$PPI - (summer_inf$seasonal_trend - global_seasonal_mean)
              # assign back
              df_inf$ppi_detrended <- NA_real_
              df_inf$ppi_detrended[summer_inf$.orig_row] <- summer_inf$ppi_detrended
              if (!"ppi_detrended" %in% trend_indices) trend_indices <- c(trend_indices, "ppi_detrended")
            } else {
              df_inf$ppi_detrended <- df_inf$PPI
              if (!"ppi_detrended" %in% trend_indices) trend_indices <- c(trend_indices, "ppi_detrended")
            }
          } else {
            df_inf$ppi_detrended <- df_inf$PPI
            if (!"ppi_detrended" %in% trend_indices) trend_indices <- c(trend_indices, "ppi_detrended")
          }
        }

        # Aggregate summer (Jun-Sep) per pheno_year
        summer_inf_all <- df_inf |> dplyr::filter(lubridate::month(date) %in% c(6,7,8,9))
        if (nrow(summer_inf_all) > 0 && length(trend_indices) > 0) {
          sum_by_year <- summer_inf_all |> dplyr::group_by(pheno_year) |>
            dplyr::summarise(across(dplyr::all_of(trend_indices), ~ mean(.x, na.rm = TRUE)), n = dplyr::n(), .groups = "drop")

          if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
          # Use fixed filenames (no date/time) so outputs are deterministic and do not embed system timestamps
          csvfile <- file.path(OUT_DIR, "inference_trend_summary.csv")
          tryCatch({ utils::write.csv(sum_by_year, file = csvfile, row.names = FALSE); cat(sprintf("[TREND] Saved inference trend summary CSV: %s\n", csvfile)) }, error = function(e) cat(sprintf("[TREND] Failed to save inference trend CSV: %s\n", e$message)))

          # Plot trends
          library(ggplot2)
          plot_dt <- tidyr::pivot_longer(sum_by_year, cols = dplyr::starts_with("ppi") | dplyr::starts_with("MSAVI") | dplyr::starts_with("NDVI"), names_to = "index", values_to = "mean_val")
          if (nrow(plot_dt) > 0) {
            p_trend <- ggplot(plot_dt, aes(x = pheno_year, y = mean_val)) +
              add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
              geom_line() + geom_point() + facet_wrap(~index, scales = "free_y") + theme_minimal() +
              coord_cartesian(ylim = c(0, NA))
            # Use fixed filename (no date/time)
            pfile <- file.path(OUT_DIR, "inference_trends.png")
            tryCatch({ ggplot2::ggsave(pfile, plot = p_trend, width = 10, height = 6); cat(sprintf("[TREND] Saved inference trend plot: %s\n", pfile)) }, error = function(e) cat(sprintf("[TREND] Failed to save trend plot: %s\n", e$message)))
          }
        } else {
          cat("[TREND] No summer observations in inference data or no indices available; skipping trend summary\n")
        }
      }, error = function(e) {
        cat(sprintf("[TREND] Error during inference trend computation: %s\n", e$message))
      })

      # Determine available indices from training parameters
      if (exists("TRAINING_NORM_PARAMS") && !is.null(TRAINING_NORM_PARAMS) && !is.null(TRAINING_NORM_PARAMS$INDEX_SCALES)) {
        avail <- names(TRAINING_NORM_PARAMS$INDEX_SCALES)
      } else {
        avail <- OPTIMAL_INDICES
      }

      # Ensure avail includes all indices required by the trained model (MESMA_PARAMS)
      if (exists("MESMA_PARAMS") && !is.null(MESMA_PARAMS) && !is.null(MESMA_PARAMS$indices)) {
        avail <- unique(c(avail, MESMA_PARAMS$indices))
      }

      # Calculate missing indices
      missing_idx <- setdiff(avail, names(df_inf))

      dvi_soil_arg <- NA_real_
      
      # If PPI not in df_inf, try to add it
      if (exists("auto_add_ppi_columns") && "PPI" %in% avail && !"PPI" %in% names(df_inf)) {
                    # Attempt to auto-add PPI to inference data using per-location baselines only.
                    dvi_soil_arg <- NA_real_
                    # Check for critical columns needed for PPI in df_inf
                    critical_cols <- c("date", "lat", "nir", "red")
                    missing_cols <- setdiff(critical_cols, names(df_inf))
                    if (length(missing_cols) > 0) {
                      ppi_inf_res <- list(df = df_inf, added = FALSE, reason = sprintf("Missing critical columns: %s", paste(missing_cols, collapse = ", ")))
                    } else {
                      ppi_inf_res <- tryCatch({
                        auto_add_ppi_columns(df_inf, dvi_soil = dvi_soil_arg)
                      }, error = function(e) {
                        list(df = df_inf, added = FALSE, reason = e$message)
                      })
                    }
                  }
          
                    if (!is.null(ppi_inf_res) && isTRUE(ppi_inf_res$added)) {
                      df_inf <- ppi_inf_res$df
                      cat(sprintf("[PPI] Auto-added PPI to inference data (reason: %s)\n", ppi_inf_res$reason))
                    } else {
                      if (!is.null(ppi_inf_res) && !is.null(ppi_inf_res$df)) {
                        df_inf <- ppi_inf_res$df
                      }
                      if (!"PPI" %in% names(df_inf)) df_inf$PPI <- NA_real_
                    }
      }

      # Re-calculate missing indices in case PPI was added
      missing_idx <- setdiff(avail, names(df_inf))

      # Backup raw PPI for inference visualization

      if (length(missing_idx) > 0) { cat(sprintf("[WARNING] Inference data missing indices: %s. Filling with NA (will likely fail unmixing).\n", paste(missing_idx, collapse=", "))); for (col in missing_idx) df_inf[[col]] <- NA_real_ }

      if ("prediction_date" %in% names(df_inf)) df_inf$prediction_date <- as.Date(df_inf$prediction_date)
      if ("reference_date" %in% names(df_inf)) df_inf$reference_date <- as.Date(df_inf$reference_date)

      # AGENT: Backup raw PPI for visualization weighting (inference)
      if ("PPI" %in% names(df_inf)) {
        df_inf$PPI_raw <- df_inf$PPI
        cat("[NOTICE] Backed up raw inference PPI values to 'PPI_raw' before normalization.\n")
      }

      # CRITICAL: Check which indices will be used for inference
      if (exists("MESMA_PARAMS") && !is.null(MESMA_PARAMS) && !is.null(MESMA_PARAMS$indices)) {
        required_indices <- MESMA_PARAMS$indices
        available_indices <- intersect(required_indices, names(df_inf))
        missing_indices <- setdiff(required_indices, names(df_inf))
        
        cat(sprintf("\n[INFERENCE FEATURE SPACE CHECK]\n"))
        cat(sprintf("  Required indices (from MESMA_PARAMS): %d indices\n", length(required_indices)))
        cat(sprintf("    %s\n", paste(required_indices, collapse=", ")))
        cat(sprintf("  Available in inference data: %d/%d\n", length(available_indices), length(required_indices)))
        if (length(missing_indices) > 0) {
          cat(sprintf("  [ERROR] MISSING indices: %s\n", paste(missing_indices, collapse=", ")))
          stop(sprintf("Inference data is missing %d required indices. Cannot proceed.", length(missing_indices)))
        }
        cat(sprintf("  [OK] All required indices present in inference data\n\n"))
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

      # Compute average value of first pentad of first index for inference data
      if (exists("MESMA_PARAMS") && !is.null(MESMA_PARAMS$indices) && length(MESMA_PARAMS$indices) > 0) {
        first_index <- MESMA_PARAMS$indices[1]
        col_name <- paste0(first_index, "_1")
        if (col_name %in% names(df_inf)) {
          avg_inference <- mean(df_inf[[col_name]], na.rm = TRUE)
          cat(sprintf("[INFO] Average value of first pentad of first index (inference, after normalization): %.6f\n", avg_inference))
        } else {
          cat("[WARNING] Column for first pentad of first index not found in inference data\n")
        }
      }

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





  # =============================================================================
  # FCLS: Fully Constrained Least Squares Solver
  # Enforces BOTH non-negativity AND sum-to-one constraints simultaneously
  # Uses iterative active-set method (Heinz & Chang, 2001)
  # =============================================================================
  
  solve_weights_fcls <- function(E, y, feature_weights = NULL, max_iter = 500, tol = 1e-8) {
    if (is.null(E) || ncol(E) < 1) return(NULL)

    # Ensure numeric
    E_fit <- as.matrix(E)
    y_fit <- as.numeric(y)
    n_endmembers <- ncol(E_fit)
    n_bands <- nrow(E_fit)

    # Handle length mismatch
    if (length(y_fit) != n_bands) {
      if (length(y_fit) > n_bands) y_fit <- y_fit[1:n_bands] else y_fit <- c(y_fit, rep(0, n_bands - length(y_fit)))
    }

    # Impute non-finite values
    y_fit[!is.finite(y_fit)] <- 0
    E_fit[!is.finite(E_fit)] <- 0

    # Base feature weights (from PCA-LDA or uniform)
    if (!is.null(feature_weights) && length(feature_weights) == n_bands) {
      feature_weights <- as.numeric(feature_weights)
      feature_weights[!is.finite(feature_weights)] <- 0
      base_weights <- pmax(feature_weights, 0)
    } else {
      base_weights <- rep(1, n_bands)
    }

    # Check if Huber loss is enabled
    use_huber <- exists("USE_HUBER_LOSS") && isTRUE(USE_HUBER_LOSS)

    if (use_huber) {
      # === HUBER LOSS via IRLS (Iteratively Reweighted Least Squares) ===
      # Huber loss: L(r) = 0.5*r^2 if |r| <= delta, else delta*(|r| - 0.5*delta)
      # IRLS weights: w_i = 1 if |r_i| <= delta, else delta/|r_i|

      huber_delta <- if (exists("HUBER_DELTA")) HUBER_DELTA else 1.345
      huber_max_iter <- if (exists("HUBER_MAX_ITER")) HUBER_MAX_ITER else 20
      huber_tol <- if (exists("HUBER_TOL")) HUBER_TOL else 1e-4

      # Initialize with uniform weights for IRLS
      irls_weights <- rep(1, n_bands)
      w_qp_prev <- rep(1/n_endmembers, n_endmembers)

      for (irls_iter in 1:huber_max_iter) {
        # Combine base weights with IRLS weights
        combined_weights <- base_weights * irls_weights
        E_w <- E_fit * combined_weights
        y_w <- y_fit * combined_weights

        # Solve weighted QP
        Dmat <- 2 * crossprod(E_w)
        dvec <- 2 * crossprod(E_w, y_w)

        qp_scale <- mean(diag(Dmat))
        if(is.na(qp_scale) || qp_scale < 1e-12) qp_scale <- 1.0

        Dmat <- Dmat / qp_scale
        dvec <- dvec / qp_scale

        ridge <- 1e-6
        Dmat <- Dmat + diag(n_endmembers) * ridge
        Dmat <- (Dmat + t(Dmat)) / 2

        Amat <- cbind(rep(1, n_endmembers), diag(n_endmembers))
        bvec <- c(1, rep(0, n_endmembers))

        res_qp <- tryCatch({
          quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
        }, error = function(e) NULL)

        if (is.null(res_qp)) {
          # QP failed, use previous solution
          break
        }

        w_qp <- res_qp$solution
        w_qp[!is.finite(w_qp)] <- 0
        w_qp[w_qp < 0] <- 0
        w_sum <- sum(w_qp, na.rm = TRUE)
        if(w_sum > 0) w_qp <- w_qp / w_sum else w_qp <- rep(1/n_endmembers, n_endmembers)

        # Calculate residuals
        pred <- as.numeric(E_fit %*% w_qp)
        resid <- y_fit - pred

        # Check convergence
        if (max(abs(w_qp - w_qp_prev)) < huber_tol) {
          break
        }
        w_qp_prev <- w_qp

        # Update IRLS weights based on Huber loss
        # Scale residuals by MAD for robust delta scaling
        mad_resid <- mad(resid, na.rm = TRUE)
        if (!is.finite(mad_resid) || mad_resid < 1e-10) mad_resid <- 1
        scaled_resid <- abs(resid) / mad_resid

        # Huber weights: 1 for small residuals, delta/|r| for large ones
        irls_weights <- ifelse(scaled_resid <= huber_delta,
                               1,
                               huber_delta / pmax(scaled_resid, 1e-10))
        irls_weights[!is.finite(irls_weights)] <- 1
      }

      # Final solution
      pred <- as.numeric(E_fit %*% w_qp)
      resid <- y_fit - pred

      # Compute Huber loss instead of RMSE
      mad_resid <- mad(resid, na.rm = TRUE)
      if (!is.finite(mad_resid) || mad_resid < 1e-10) mad_resid <- 1
      scaled_resid <- abs(resid) / mad_resid
      huber_losses <- ifelse(scaled_resid <= huber_delta,
                             0.5 * (resid/mad_resid)^2,
                             huber_delta * (scaled_resid - 0.5 * huber_delta))
      loss <- sqrt(mean(huber_losses))  # Return sqrt for consistency with RMSE scale

      return(list(w = w_qp, rmse = loss, residuals = resid, loss_type = "huber"))

    } else {
      # === STANDARD RMSE (Original behavior) ===
      E_w <- E_fit * base_weights
      y_w <- y_fit * base_weights

      # --- DEDICATED QP SOLVER (ALWAYS used) ---
      # Problem: min ||Ex - y||^2  s.t. sum(x)=1, x>=0
      # ||Ex - y||^2 = x'E'Ex - 2y'Ex + y'y
      # quadprog solves: min 1/2 b^T Dmat b - dvec^T b
      # Dmat = 2 * E'E, dvec = 2 * E'y

      Dmat <- 2 * crossprod(E_w)
      dvec <- 2 * crossprod(E_w, y_w)

      # Scale QP to prevent inconsistent constraints
      qp_scale <- mean(diag(Dmat))
      if(is.na(qp_scale) || qp_scale < 1e-12) qp_scale <- 1.0

      Dmat <- Dmat / qp_scale
      dvec <- dvec / qp_scale

      # Add stronger ridge to diagonal for numerical stability (PD requirement)
      ridge <- 1e-6
      Dmat <- Dmat + diag(n_endmembers) * ridge

      # Ensure Dmat is symmetric
      Dmat <- (Dmat + t(Dmat)) / 2

      # Constraints: A^T b >= b_0
      # 1. sum(x) = 1  => use meq=1. Row 1 of A^T is [1, 1, ..., 1].
      # 2. x >= 0      => I * x >= 0. Rows 2..N+1 of A^T are I.

      Amat <- cbind(rep(1, n_endmembers), diag(n_endmembers))
      bvec <- c(1, rep(0, n_endmembers))

      res_qp <- tryCatch({
        quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
      }, error = function(e) {
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG] Quadprog failed in solve_weights_fcls: %s\n", e$message))
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG] Dmat dim: %s, dvec length: %d, Amat dim: %s, bvec length: %d\n",
                    paste(dim(Dmat), collapse="x"), length(dvec), paste(dim(Amat), collapse="x"), length(bvec)))
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG] Dmat eigenvalues: %s\n", paste(eigen(Dmat)$values[1:5], collapse=", ")))
        NULL
      })

      if (!is.null(res_qp)) {
        w_qp <- res_qp$solution
        w_qp[!is.finite(w_qp)] <- 0
        w_qp[w_qp < 0] <- 0
        w_sum <- sum(w_qp, na.rm = TRUE)
        if(w_sum > 0) w_qp <- w_qp / w_sum else w_qp <- rep(1/n_endmembers, n_endmembers)

        # Calculate RMSE for consistency
        pred <- as.numeric(E_fit %*% w_qp)
        resid <- y_fit - pred
        rmse <- sqrt(mean(resid^2))

        return(list(w = w_qp, rmse = rmse, residuals = resid, loss_type = "rmse"))
      }

      # Fallback to uniform weights if QP fails (should affect very few pixels)
      w_qp <- rep(1 / n_endmembers, n_endmembers)
      pred <- as.numeric(E_fit %*% w_qp)
      resid <- y_fit - pred
      rmse <- sqrt(mean(resid^2))

      return(list(w = w_qp, rmse = rmse, residuals = resid, loss_type = "rmse"))
    }
  }

  # Batch FCLS Solver using Quadprog (Optimization for GA/Grid Search)
  # E: Features x Endmembers matrix
  # Y: Samples x Features matrix
  # Returns: Samples x Endmembers weight matrix
  solve_batch_fcls <- function(E, Y, feature_weights = NULL) {
    if (is.null(E) || ncol(E) < 1) return(NULL)
    n_endmembers <- ncol(E)
    n_samples <- nrow(Y)
    if (n_samples == 0) return(matrix(0, 0, n_endmembers))

    # Convert Y to Features x Samples for matrix math
    Y_t <- t(Y)

    # Determine E and Y used for fitting (weighted or raw)
    if (!is.null(feature_weights) && length(feature_weights) == nrow(E)) {
      w <- pmax(feature_weights, 0)
      w[!is.finite(w)] <- 0

      # Apply weights directly: E_w = W * E, Y_w = W * Y
      # Since W is diagonal, just multiply rows
      E_fit <- E * w
      Y_fit <- Y_t * w
    } else {
      E_fit <- E
      Y_fit <- Y_t
    }

    # Precompute constant QP matrices
    # Dmat = 2 * E'E
    Dmat <- 2 * crossprod(E_fit)

    # Scale QP problem to improve numerical stability for quadprog
    qp_scale <- mean(diag(Dmat))
    if(is.na(qp_scale) || qp_scale < 1e-12) qp_scale <- 1.0

    Dmat <- Dmat / qp_scale

    # Add stronger ridge
    ridge <- 1e-6 # Relative to scaled matrix
    Dmat <- Dmat + diag(n_endmembers) * ridge
    # Ensure symmetric
    Dmat <- (Dmat + t(Dmat)) / 2

    # Amat: [1s; I]
    Amat <- cbind(rep(1, n_endmembers), diag(n_endmembers))
    bvec <- c(1, rep(0, n_endmembers))

    # Precompute all linear terms: dvec = 2 * E'y
    # Dvecs: Endmembers x Samples
    # Must also be scaled by same factor
    Dvecs <- (2 * crossprod(E_fit, Y_fit)) / qp_scale

    w_out <- matrix(0, nrow=n_samples, ncol=n_endmembers)

    # Loop over samples - solve.QP is fast when Dmat is precomputed
    for(i in 1:n_samples) {
      dvec <- Dvecs[, i]

      res <- tryCatch({
        quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq=1)
      }, error = function(e) {
        NULL
      })

      if(!is.null(res)) {
        w <- res$solution
        w[!is.finite(w)] <- 0
        w[w < 0] <- 0
        s <- sum(w)
        if(s > 0) w <- w / s else w <- rep(1/n_endmembers, n_endmembers)
        w_out[i, ] <- w
      } else {
        w_out[i, ] <- rep(1/n_endmembers, n_endmembers)
      }
    }
    return(w_out)
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


  

  # ==========================================================================
  # Validation and inference datasets already created during stratified split
  # Validation pipeline reconstructed in the dedicated section below

  cat("\n=== DATA DISTRIBUTION ANALYSIS ===\n")
  if (exists("df_tasks") && nrow(df_tasks) > 0) {
    cat(sprintf("Total locations in df_tasks: %d\n", length(unique(df_tasks$location_id))))

    sample_sizes <- df_tasks |> 
      dplyr::group_by(location_id, pheno_year) |> 
      dplyr::summarize(n_obs = dplyr::n(), .groups = "drop")

    # Also compute distribution excluding barren (so min sample size is not driven by barren rows)
    sample_sizes_no_barren <- NULL
    if ("Veg" %in% names(df_tasks)) {
      sample_sizes_no_barren <- df_tasks |>
        dplyr::filter(!is.na(.data$Veg) & tolower(trimws(as.character(.data$Veg))) != "barren") |>
        dplyr::group_by(location_id, pheno_year) |>
        dplyr::summarize(n_obs = dplyr::n(), .groups = "drop")
    }

    if (nrow(sample_sizes) > 0) {
      sample_sizes$n_obs <- as.numeric(sample_sizes$n_obs)
      cat("\nObservations per location-year distribution:\n")
      cat(sprintf("  Min:    %d\n", as.integer(min(sample_sizes$n_obs, na.rm = TRUE))))
      cat(sprintf("  Q1:     %d\n", as.integer(stats::quantile(sample_sizes$n_obs, 0.25, na.rm = TRUE))))
      cat(sprintf("  Median: %d\n", as.integer(median(sample_sizes$n_obs, na.rm = TRUE))))
      cat(sprintf("  Q3:     %d\n", as.integer(stats::quantile(sample_sizes$n_obs, 0.75, na.rm = TRUE))))
      cat(sprintf("  Max:    %d\n", as.integer(max(sample_sizes$n_obs, na.rm = TRUE))))
      cat(sprintf("  Mean:   %.1f\n", mean(sample_sizes$n_obs, na.rm = TRUE)))

      if (!is.null(sample_sizes_no_barren) && nrow(sample_sizes_no_barren) > 0) {
        sample_sizes_no_barren$n_obs <- as.numeric(sample_sizes_no_barren$n_obs)
        cat("\nObservations per location-year distribution (excluding barren):\n")
        cat(sprintf("  Min (no barren): %d\n", as.integer(min(sample_sizes_no_barren$n_obs, na.rm = TRUE))))
      }

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

  # NOTE: location_list check moved to after validation split creates it (around line 5705)





  evaluate_all_combinations <- function(
    y,
    top_variants,
    lambda = 0,
    early_stop_rmse = 0,
    feature_weights = NULL
  ) {
    if (length(top_variants) == 0) return(NULL)

    full_veg_names <- names(top_variants)
    n_veg_full <- length(full_veg_names)
    if (n_veg_full == 0) return(NULL)

    y_target <- y
    y_target[!is.finite(y_target)] <- 0

    score_from_solution <- function(res) {
      if (is.null(res) || is.null(res$residuals)) return(Inf)
      as.numeric(res$rmse)
    }

    # Solve a combination of variant indices for a SUBSET of vegetation types
    # veg_subset: character vector of vegetation type names to include
    # variant_indices: integer vector of variant indices (one per veg type in subset)
    solve_combo_subset <- function(veg_subset, variant_indices) {
      cols <- list()
      ids <- character(length(veg_subset))
      names(ids) <- veg_subset

      for (v_idx in seq_along(veg_subset)) {
        v <- veg_subset[v_idx]
        idx <- variant_indices[v_idx]
        cand <- top_variants[[v]][[idx]]
        if (is.null(cand) || is.null(cand$vec)) return(NULL)
        cols[[length(cols) + 1]] <- as.numeric(cand$vec)
        ids[v] <- cand$id
      }
      if (length(cols) == 0) return(NULL)

      E <- do.call(cbind, cols)
      if (is.null(E) || ncol(E) < 1) return(NULL)

      res <- solve_weights_fcls(E, y_target, feature_weights = if (!is.null(feature_weights) && length(feature_weights) == nrow(E)) feature_weights else NULL)
      if (is.null(res)) return(NULL)

      # Expand weights to full veg list (zeros for excluded types)
      w_full <- rep(0, n_veg_full)
      names(w_full) <- full_veg_names
      w_full[veg_subset] <- res$w

      ids_full <- rep(NA_character_, n_veg_full)
      names(ids_full) <- full_veg_names
      ids_full[veg_subset] <- ids

      res$w <- w_full
      res$ids <- ids_full
      res$E <- E
      res$veg_subset <- veg_subset
      res$score <- score_from_solution(res)
      return(res)
    }

    # --- Search over all subsets of 1, 2, or 3 vegetation types (MAX_VEG_TYPES) ---
    MAX_VEG_TYPES <- 3L
    n_variants_per_veg <- sapply(top_variants, length)

    global_best_res <- NULL
    global_best_score <- Inf

    # Generate all subsets of vegetation types with size 1 to MAX_VEG_TYPES
    all_subsets <- list()
    for (k in 1:min(MAX_VEG_TYPES, n_veg_full)) {
      subsets_k <- combn(full_veg_names, k, simplify = FALSE)
      all_subsets <- c(all_subsets, subsets_k)
    }

    # For each subset, search variant combinations
    for (veg_subset in all_subsets) {
      n_in_subset <- length(veg_subset)
      n_variants_subset <- n_variants_per_veg[veg_subset]
      total_combos_subset <- as.numeric(prod(n_variants_subset))

      if (total_combos_subset <= 500) {
        # Full exhaustive search for this subset
        combos <- expand.grid(lapply(n_variants_subset, seq_len), KEEP.OUT.ATTRS = FALSE)
        for (i in seq_len(nrow(combos))) {
          r <- solve_combo_subset(veg_subset, as.integer(combos[i, ]))
          if (!is.null(r) && r$score < global_best_score) {
            global_best_score <- r$score
            global_best_res <- r
          }
        }
      } else {
        # Coarse search + coordinate descent for large search spaces
        TOPK_COARSE <- 3
        coarse_idx <- lapply(n_variants_subset, function(n) seq_len(min(TOPK_COARSE, n)))
        combos <- expand.grid(coarse_idx, KEEP.OUT.ATTRS = FALSE)
        best_idx <- rep(1, n_in_subset)

        for (i in seq_len(nrow(combos))) {
          r <- solve_combo_subset(veg_subset, as.integer(combos[i, ]))
          if (!is.null(r) && r$score < global_best_score) {
            global_best_score <- r$score
            global_best_res <- r
            best_idx <- as.integer(combos[i, ])
          }
        }

        # Coordinate descent refinement within this subset
        if (!is.null(global_best_res) && identical(global_best_res$veg_subset, veg_subset)) {
          curr_idx <- best_idx
          improved <- TRUE
          while (improved) {
            improved <- FALSE
            for (k in seq_len(n_in_subset)) {
              for (v_opt in seq_len(n_variants_subset[k])) {
                if (v_opt == curr_idx[k]) next
                t_idx <- curr_idx
                t_idx[k] <- v_opt
                r <- solve_combo_subset(veg_subset, t_idx)
                if (!is.null(r) && r$score < global_best_score) {
                  global_best_score <- r$score
                  global_best_res <- r
                  curr_idx <- t_idx
                  improved <- TRUE
                }
              }
            }
          }
        }
      }
    }
    
    if (is.null(global_best_res)) return(NULL)

    # [MODIFIED] Return only the single best model as top_models per request
    # Previous logic sorted accumulated candidates and kept near-optimal ones.
    top_models <- list(global_best_res)

    return(list(
      w = global_best_res$w,
      rmse = global_best_res$rmse,
      score = global_best_score,
      ids = global_best_res$ids,
      residuals = global_best_res$residuals,
      E_best = global_best_res$E,
      top_models = top_models
    ))
  }

  
  build_mesma_library_weighted <- function(df_train, indices, params, allowed_veg, precomputed_clusters = NULL) {
    # MESMA library: treats barren and all vegetation types as equal endmembers
    # Returns a library with all endmember types (barren + veg types) in one unified structure
    # If precomputed_clusters is provided (list with class -> k mapping), skip cluster optimization

    if (!is.null(precomputed_clusters)) {
      cat("\n[LIBRARY BUILD] Using precomputed cluster counts (skipping optimization)...\n")
      cat(sprintf("[LIBRARY BUILD] Precomputed clusters: %s\n", paste(names(precomputed_clusters), precomputed_clusters, sep="=", collapse=", ")))
    } else {
      cat("\n[LIBRARY BUILD] Starting Global Combinatorial Optimization for Cluster Counts (including Barren)...\n")
    }

    # Validate params structure
    if (is.null(params) || is.null(params$means) || is.null(params$sds)) {
      stop("[ERROR] build_mesma_library_weighted: params$means or params$sds is NULL!")
    }
    if (length(params$means) != length(indices)) {
      cat(sprintf("[ERROR] params$means has length %d but indices has length %d\n", length(params$means), length(indices)))
      cat(sprintf("  indices: %s\n", paste(indices, collapse = ", ")))
      cat(sprintf("  params$means names: %s\n", paste(names(params$means), collapse = ", ")))
      stop("[ERROR] Length mismatch between params$means and indices")
    }

    # -------------------------------------------------------------------------
    # PASS 1: Load, Normalize, and Store Data (with Metadata)
    # Separate storage for model training and OOB evaluation
    # -------------------------------------------------------------------------
    expected_cols <- length(indices) * TEMPORAL_BUDGET
    storage <- list()       # Store unweighted Z-score matrices (for model training)
    storage_meta <- list()  # Store metadata (location_id) for spatial bootstrapping
    storage_oob <- list()   # Store OOB data for cluster optimization evaluation
    storage_oob_meta <- list()

    # Verify OOB data is available (required, no fallback)
    if (!exists("df_train_oob", envir = globalenv()) ||
        is.null(get("df_train_oob", envir = globalenv())) ||
        nrow(get("df_train_oob", envir = globalenv())) == 0) {
      stop("[ERROR] OOB holdout data (df_train_oob) is required but not found. Ensure OOB split was performed during training data preparation.")
    }

    df_oob <- get("df_train_oob", envir = globalenv())
    oob_locs <- get("oob_location_ids", envir = globalenv())
    has_oob_data <- TRUE  # OOB data verified above
    cat(sprintf("[LIBRARY BUILD] OOB data available: %d rows from %d locations\n",
                nrow(df_oob), length(oob_locs)))
    cat("[LIBRARY BUILD] Will use OOB data for cluster optimization evaluation\n")

    # Unified class list
    target_classes <- unique(c("barren", allowed_veg))
    valid_classes <- c()

    # Helper function to build storage from a dataframe
    build_storage_from_df <- function(df_source, storage_list, meta_list, source_name = "train") {
      dual_mode <- isTRUE(params$dual_mode)
      l2_normalize <- isTRUE(params$l2_normalize)
      base_indices <- if (!is.null(params$base_indices)) params$base_indices else indices
      n_base_idx <- length(base_indices)
      expected_base_cols <- n_base_idx * TEMPORAL_BUDGET

      for(v in target_classes) {
        veg_data <- dplyr::filter(df_source, .data$Veg == v)
        if(nrow(veg_data) == 0) next

        veg_list <- list()
        loc_list <- character(0)

        traces <- unique(veg_data[, c("location_id", "pheno_year")])

        for(i in seq_len(nrow(traces))) {
          lid <- traces$location_id[i]
          pyr <- traces$pheno_year[i]
          sub <- veg_data[veg_data$location_id == lid & veg_data$pheno_year == pyr, ]
          # Always build from base indices
          mat <- build_pentad_matrix(sub, base_indices)
          if(!is.null(mat)) {
            veg_list[[length(veg_list) + 1]] <- as.numeric(mat)
            loc_list <- c(loc_list, as.character(lid))
          }
        }

        if(length(veg_list) > 0) {
          veg_mat_raw <- do.call(rbind, veg_list)

          if(ncol(veg_mat_raw) == expected_base_cols) {
            if (dual_mode) {
              # DUAL MODE: create both raw and L2-normalized, then concatenate
              veg_mat_l2 <- t(apply(veg_mat_raw, 1, function(r) {
                l2_normalize_perindex(r, n_base_idx, TEMPORAL_BUDGET)
              }))
              veg_mat <- cbind(veg_mat_raw, veg_mat_l2)
            } else if (l2_normalize) {
              # L2 ONLY MODE
              veg_mat <- t(apply(veg_mat_raw, 1, function(r) {
                l2_normalize_perindex(r, n_base_idx, TEMPORAL_BUDGET)
              }))
            } else {
              # RAW ONLY MODE
              veg_mat <- veg_mat_raw
            }

            # Z-score all indices (raw, L2, or both)
            for(k in seq_along(indices)) {
              idx_start <- (k-1)*TEMPORAL_BUDGET + 1
              idx_end <- k*TEMPORAL_BUDGET
              idx_name <- indices[k]
              param_idx <- which(names(params$means) == idx_name)
              if (length(param_idx) > 0) {
                veg_mat[, idx_start:idx_end] <- (veg_mat[, idx_start:idx_end] - params$means[param_idx]) / params$sds[param_idx]
              }
            }
            veg_mat[!is.finite(veg_mat)] <- 0

            storage_list[[v]] <- veg_mat
            meta_list[[v]] <- data.frame(location_id = loc_list, stringsAsFactors = FALSE)
          }
        }
      }
      return(list(storage = storage_list, meta = meta_list))
    }

    # Build storage from training data (excluding OOB if available)
    if (has_oob_data) {
      # Use df_train_model for building endmembers
      df_train_model_local <- get("df_train_model", envir = globalenv())
      cat(sprintf("[LIBRARY BUILD] Building endmember storage from model data: %d rows\n", nrow(df_train_model_local)))
      result <- build_storage_from_df(df_train_model_local, storage, storage_meta, "model")
      storage <- result$storage
      storage_meta <- result$meta

      # Build OOB storage for evaluation
      cat(sprintf("[LIBRARY BUILD] Building OOB storage for evaluation: %d rows\n", nrow(df_oob)))
      result_oob <- build_storage_from_df(df_oob, storage_oob, storage_oob_meta, "oob")
      storage_oob <- result_oob$storage
      storage_oob_meta <- result_oob$meta

      # Log OOB storage sizes
      for (v in names(storage_oob)) {
        cat(sprintf("[LIBRARY BUILD] OOB storage for %s: %d samples\n", v, nrow(storage_oob[[v]])))
      }
    }

    valid_classes <- names(storage)

    if (length(valid_classes) == 0) {
      warning("No valid data found for library building.")
      return(list())
    }

    # ------------------------------------------------------------------------- 
    # Pre-filter vegetation spectra for barren-like signatures BEFORE optimization
    # -------------------------------------------------------------------------
    if ("barren" %in% valid_classes && !is.null(storage[["barren"]]) && nrow(storage[["barren"]]) > 0) {
      barren_ref <- colMeans(storage[["barren"]], na.rm = TRUE)
      barren_norm <- barren_ref / sqrt(sum(barren_ref^2))
      
      cosine_threshold <- BARREN_SIM_THRESHOLD
      
      for (v in valid_classes) {
        if (v == "barren" || is.null(storage[[v]]) || nrow(storage[[v]]) == 0) next
        
        veg_mat <- storage[[v]]
        veg_norm <- t(apply(veg_mat, 1, function(r) {
          nrm <- sqrt(sum(r^2))
          if (!is.finite(nrm) || nrm == 0) return(rep(0, length(r))) else return(r / nrm)
        }))
        
        sims <- veg_norm %*% barren_norm
        keep_mask <- sims <= cosine_threshold
        
        n_before <- nrow(veg_mat)
        if (any(!keep_mask)) {
          storage[[v]] <- veg_mat[keep_mask, , drop = FALSE]
          storage_meta[[v]] <- storage_meta[[v]][keep_mask, , drop = FALSE]
          n_after <- nrow(storage[[v]])
          cat(sprintf("[BARREN FILTER] Pre-filtered %s: removed %d/%d spectra (similarity > %.2f)\n", 
                      v, n_before - n_after, n_before, cosine_threshold))
        }
      }
    }

    # Compute average value of first pentad of first index for training data
    if (length(valid_classes) > 0) {
      all_vals <- c()
      for (v in valid_classes) {
        if (!is.null(storage[[v]]) && nrow(storage[[v]]) > 0) {
          all_vals <- c(all_vals, storage[[v]][, 1])
        }
      }
      if (length(all_vals) > 0) {
        avg_train <- mean(all_vals, na.rm = TRUE)
        cat(sprintf("[INFO] Average value of first pentad of first index (training, after normalization): %.6f\n", avg_train))
      }
    }

    # Weights for clustering (PCA-LDA weights for distance calculations)
    # Storage matrices are built to match `indices` exactly, so they always have `expected_cols` columns.
    n_storage_cols <- expected_cols
    w_vec <- if(!is.null(params$weights) && length(params$weights) == n_storage_cols) {
      pmax(params$weights, 0)
    } else {
      rep(1, n_storage_cols)
    }
    w_vec[w_vec < 1e-9] <- 1e-9 # Prevent div by zero
    
    min_cluster_size <- MIN_CLUSTER_SIZE

    # Greedy EAR-based endmember extraction (minimizes reconstruction error)
    ear_extract_all_levels_greedy <- function(data_mat, max_k, w_vec = NULL) {
      n <- nrow(data_mat)
      if (n == 0) return(NULL)

      
      # Apply weighting
      if (!is.null(w_vec) && length(w_vec) == ncol(data_mat)) {
        data_w <- sweep(data_mat, 2, w_vec, "*")
      } else {
        data_w <- data_mat
      }
      
      # Determine max_k effectively
      max_k <- min(max_k, n)
      
      # Precompute Distance Matrix (O(N^2))
      # This is heavy but less than doing it k times
      if (n > 4000) {
         # If N is huge, random subsample to avoid memory OOM
         n_sub <- 4000
         set.seed(42)
         idx <- sort(sample(n, n_sub))
         data_sub <- data_w[idx, , drop = FALSE]
         dist_mat <- as.matrix(dist(data_sub))^2
         candidates <- idx
         mapping <- setNames(seq_along(idx), idx) # map original idx to subsample idx
         is_subsampled <- TRUE
         orig_to_sub <- function(id) mapping[[as.character(id)]]
      } else {
         dist_mat <- as.matrix(dist(data_w))^2
         candidates <- 1:n
         is_subsampled <- FALSE
      }

      use_robust <- exists("ROBUST_CLUSTERING") && isTRUE(ROBUST_CLUSTERING)
      
      selected <- integer(0)
      min_sq_dists <- rep(Inf, if (is_subsampled) nrow(data_sub) else n)
      
      results <- list()
      
      for (i in 1:max_k) {
        if (i == 1) {
          # Step 1: Find Medoid
          # For squared Euclidean, minimizing sum of squared dists = Sample closest to Mean
          # Robust = Sample closest to Geometric Median
          
          if (is_subsampled) d_matrix <- data_sub else d_matrix <- data_w
          
          if (use_robust) {
            robust_centroid <- apply(d_matrix, 2, median, na.rm = TRUE)
          } else {
            robust_centroid <- colMeans(d_matrix, na.rm = TRUE)
          }
          
          dists_to_centroid <- rowSums(sweep(d_matrix, 2, robust_centroid, "-")^2)
          best_local_idx <- which.min(dists_to_centroid)
          best_idx <- candidates[best_local_idx]
          
        } else {
          remaining_local_idx <- seq_along(candidates) 
          remaining_local_idx <- remaining_local_idx[!candidates[remaining_local_idx] %in% selected]
          
          # We want to maximize total gain = sum(pmax(0, current_err - new_err))
          # new_err = dist(pt, cand)
          # gain = pmax(0, min_sq_dists - D[pt, cand])
          
          # Vectorized over all candidates
          # Calculate gains for all remaining candidates in one matrix op if poss,
          # or loop over candidates but use the precomputed dist_mat
          
          best_gain <- -1
          best_local_idx <- -1
          
          # Fast scan using precomputed distance matrix
          # dist_mat is (N_sub x N_sub)
          # min_sq_dists is (N_sub)
          
          # Total Error with current set = sum(min_sq_dists)
          # If we add candidate C: new_error_sum = sum(min(min_sq_dists[j], dist_mat[j, C]))
          
          # We can compute this for all C in 'remaining' quickly
          current_err_sum <- sum(min_sq_dists)
          
          # This loop is effectively vectorized compared to re-calculating distances
          for(cand_loc in remaining_local_idx) {
             d_col <- dist_mat[, cand_loc]
             # Total error if we add this candidate
             new_total <- sum(pmin(min_sq_dists, d_col))
             gain <- current_err_sum - new_total
             
             if(gain > best_gain) {
               best_gain <- gain
               best_local_idx <- cand_loc
             }
          }
          best_idx <- candidates[best_local_idx]
        }
        
        selected <- c(selected, best_idx)
        
        # Update distances
        # dist_mat is [local_idx, local_idx]
        # best_local_idx corresponds to the selected candidate
        if (i == 1) {
           # Initial distances logic
           if (is_subsampled) best_local <- which(candidates == best_idx) else best_local <- best_idx
           min_sq_dists <- dist_mat[, best_local]
        } else {
           if (is_subsampled) best_local <- which(candidates == best_idx) else best_local <- best_idx
           min_sq_dists <- pmin(min_sq_dists, dist_mat[, best_local])
        }
        
        # Save snapshot
        results[[as.character(i)]] <- data_mat[selected, , drop = FALSE]
      }
      return(results)
    }

    # Greedy EAR-based endmember extraction wrapper
    ear_extract_all_levels <- function(data_mat, max_k, w_vec = NULL) {
      cat(sprintf("[ENDMEMBER EXTRACTION] Using greedy EAR (n=%d, max_k=%d)\n", nrow(data_mat), max_k))
      return(ear_extract_all_levels_greedy(data_mat, max_k, w_vec))
    }

    # Endmember extraction wrapper: returns endmembers for a specific k
    ear_extract_endmembers <- function(data_mat, k, w_vec = NULL) {
      if (k < 1) return(NULL)
      res <- ear_extract_all_levels(data_mat, k, w_vec)
      return(res[[as.character(k)]])
    }
    
    optimize_library <- function(n_boot = 5) {
      cat(sprintf("\n  --- Running EAR-Based Endmember Selection with Cross-Validation (%d folds) ---\n", n_boot))
      
      # Pre-filter storage data to remove barren-similar variants before optimization
      cat("    Pre-filtering vegetation variants to remove barren-similar spectra...\n")
      
      # Prepare barren reference for filtering
      barren_refs_filter <- NULL
      if ("barren" %in% valid_classes) {
        b_mat <- storage[["barren"]]
        if (!is.null(b_mat) && nrow(b_mat) > 0) {
          barren_refs_filter <- matrix(colMeans(b_mat, na.rm = TRUE), nrow = 1)
        }
      }
      
      if (!is.null(barren_refs_filter)) {
        # Normalize barren reference
        barren_refs_filter_norm <- barren_refs_filter / sqrt(sum(barren_refs_filter^2))
        
        for (v in valid_classes) {
          if (v == "barren") next  # Don't filter barren itself
          
          veg_mat <- storage[[v]]
          if (is.null(veg_mat) || nrow(veg_mat) == 0) next
          
          n_before <- nrow(veg_mat)
          
          # Normalize vegetation variants for cosine similarity
          veg_mat_norm <- t(apply(veg_mat, 1, function(r) {
            nrm <- sqrt(sum(r^2))
            if (!is.finite(nrm) || nrm == 0) return(rep(0, length(r))) else return(r / nrm)
          }))
          
          # Compute cosine similarities to barren
          sim_mat <- veg_mat_norm %*% t(barren_refs_filter_norm)
          max_sims <- apply(sim_mat, 1, function(x) suppressWarnings(max(x, na.rm = TRUE)))
          keep_mask <- !(max_sims > BARREN_SIM_THRESHOLD)
          
          # Always keep at least one variant
          if (any(keep_mask)) {
            storage[[v]] <- veg_mat[keep_mask, , drop=FALSE]
            storage_meta[[v]] <- storage_meta[[v]][keep_mask, , drop=FALSE]
          } else {
            # Keep the least barren-like variant
            min_idx <- which.min(max_sims)
            storage[[v]] <- veg_mat[min_idx, , drop=FALSE]
            storage_meta[[v]] <- storage_meta[[v]][min_idx, , drop=FALSE]
          }
          
          n_after <- nrow(storage[[v]])
          if (n_after < n_before) {
            cat(sprintf("      %s: removed %d/%d variants (barren similarity > %.2f)\n",
                        v, n_before - n_after, n_before, BARREN_SIM_THRESHOLD))
          }
        }
      }
      
      # Pre-compute EAR-based endmembers on full model data
      boot_endmember_cache <- list()
      k_ranges <- list()

      cat("    Pre-computing EAR-based endmembers on full model data...\n")

      for(v in valid_classes) {
        n_total <- nrow(storage[[v]])
        max_k <- floor(n_total / min_cluster_size)
        if (v == "barren") {
           # Allow multiple barren endmembers controlled by RAW_BARREN_N_PROTOTYPES
           max_k_barren <- if (exists("RAW_BARREN_N_PROTOTYPES") && is.finite(RAW_BARREN_N_PROTOTYPES) && RAW_BARREN_N_PROTOTYPES >= 1) as.integer(RAW_BARREN_N_PROTOTYPES) else MAX_K_EAR
           k_candidates <- 1:max_k_barren
           k_candidates <- k_candidates[k_candidates <= max_k]
           if (length(k_candidates) == 0) k_candidates <- 1
        } else {
           k_candidates <- 1:MAX_K_EAR
           k_candidates <- k_candidates[k_candidates <= max_k]
           if (length(k_candidates) == 0) k_candidates <- 1
        }
        k_ranges[[v]] <- k_candidates
      }

      # Build endmembers from full storage (model data)
      boot_endmember_cache[[1]] <- list()

      for(v in valid_classes) {
         boot_endmember_cache[[1]][[v]] <- list()

         veg_mat_train <- storage[[v]]

         # Determine max K needed for this class
         max_k_needed <- max(k_ranges[[v]])

         cat(sprintf("[DEBUG] Extracting endmembers for '%s': n_samples=%d, max_k=%d\n",
                     v, nrow(veg_mat_train), max_k_needed))

         # Compute ALL levels at once
         all_levels <- ear_extract_all_levels(veg_mat_train, max_k_needed, w_vec)

         for(k in k_ranges[[v]]) {
           k_char <- as.character(k)
           if(!is.null(all_levels[[k_char]])) {
             boot_endmember_cache[[1]][[v]][[k_char]] <- all_levels[[k_char]]
           }
         }
      }

      # --- Random Search with OOB Evaluation ---
      cat(sprintf("    Optimizing library structure using OOB holdout data for evaluation...\n"))
      cat(sprintf("    (OOB classes available: %s)\n", paste(names(storage_oob), collapse = ", ")))

      # 3a. Prepare Test Matrices using OOB data
      train_sets <- list()
      oob_sets <- list()

      # Build OOB evaluation set from storage_oob
      oob_samples <- list(); oob_lbls <- c()
      for (v in valid_classes) {
        if (v %in% names(storage_oob) && !is.null(storage_oob[[v]]) && nrow(storage_oob[[v]]) > 0) {
          oob_samples[[length(oob_samples)+1]] <- storage_oob[[v]]
          oob_lbls <- c(oob_lbls, rep(v, nrow(storage_oob[[v]])))
        }
      }

      if (length(oob_samples) == 0) {
        stop("[ERROR] No OOB samples found in storage_oob. Cannot proceed with cluster optimization.")
      }

      oob_sets[[1]] <- list(Y = do.call(rbind, oob_samples), labels = oob_lbls)
      cat(sprintf("    OOB evaluation set: %d samples across %d classes\n",
                  nrow(oob_sets[[1]]$Y), length(unique(oob_lbls))))

      # Use FCLS solver for evaluation (wrapper)
      run_solver_fcls <- function(y, M, w) {
        res <- solve_weights_fcls(M, y, feature_weights = w)
        if (is.null(res)) return(rep(0, ncol(M)))
        return(res$w)
      }

      # Pre-compute Barren References for filtering during evaluation
      # Use representative barren endmembers (via EAR) to support multiple barren prototypes
      barren_refs_eval_norm <- NULL
      if ("barren" %in% valid_classes) {
         b_mat <- storage[["barren"]]
         if (!is.null(b_mat) && nrow(b_mat) > 0) {
           # Get representative barren endmembers using EAR (up to RAW_BARREN_N_PROTOTYPES)
           n_bproto <- if (exists("RAW_BARREN_N_PROTOTYPES") && is.finite(RAW_BARREN_N_PROTOTYPES) && RAW_BARREN_N_PROTOTYPES >= 1) as.integer(RAW_BARREN_N_PROTOTYPES) else 1L
           if (n_bproto > 1 && nrow(b_mat) >= n_bproto) {
             barren_refs <- ear_extract_endmembers(b_mat, min(n_bproto, nrow(b_mat)), w_vec)
             if (is.null(barren_refs) || nrow(barren_refs) == 0) {
               barren_refs <- matrix(colMeans(b_mat, na.rm = TRUE), nrow = 1)
             }
           } else {
             barren_refs <- matrix(colMeans(b_mat, na.rm = TRUE), nrow = 1)
           }
           # Normalize each barren reference
           barren_refs_eval_norm <- t(apply(barren_refs, 1, function(r) {
             nrm <- sqrt(sum(r^2))
             if (!is.finite(nrm) || nrm == 0) return(rep(0, length(r))) else return(r / nrm)
           }))
           cat(sprintf("[BARREN FILTER SETUP] %d barren reference(s) prepared for filtering (dim: %dx%d)\n",
                       nrow(barren_refs_eval_norm), nrow(barren_refs_eval_norm), ncol(barren_refs_eval_norm)))
         } else {
           cat("[WARNING] Barren class exists but storage is NULL or empty - barren filtering disabled\n")
         }
      } else {
        cat("[NOTICE] Barren class not in valid_classes - barren filtering disabled\n")
      }

      # 3b. Fitness Function - Evaluate on OOB holdout data
      # Optimizes VEGETATION ACCURACY (populus, tamarix, herbs) - excludes barren
      evaluate_config <- function(combo_list) {
          # Single evaluation on OOB holdout (no bootstrap CV)
          if(is.null(oob_sets[[1]])) return(-1.0)
          target_set <- oob_sets[[1]]

          cols <- list(); col_names <- c()
          valid_lib <- TRUE

          for(v in valid_classes) {
             # combo_list can be list or named vector
             k_val <- as.character(combo_list[[v]])
             # Use endmembers from the single cache entry (index 1)
             endmembers <- boot_endmember_cache[[1]][[v]][[k_val]]
             if(is.null(endmembers)) { valid_lib <- FALSE; break }

             cols[[length(cols)+1]] <- t(endmembers)
             col_names <- c(col_names, rep(v, nrow(endmembers)))
          }

          if(!valid_lib) return(-1.0)

          M <- do.call(cbind, cols)

          # Use OOB data for fitness evaluation
          Y_test <- target_set$Y
          labels_test <- target_set$labels

          # Vegetation classes only (exclude barren)
          veg_classes <- setdiff(valid_classes, "barren")

          n_test <- nrow(Y_test)

          # Batch process using solve_batch_fcls for speed
          all_coefs <- solve_batch_fcls(M, Y_test, params$weights)

          # Compute average correctly predicted fraction for vegetation classes
          # Row-normalized after barren subtraction: veg_frac_true / sum(veg_fracs)
          veg_norm_frac_sums <- setNames(rep(0, length(veg_classes)), veg_classes)
          veg_counts <- setNames(rep(0L, length(veg_classes)), veg_classes)

          for(j in 1:n_test) {
             true_label <- labels_test[j]
             coefs <- all_coefs[j, ]
             sums <- tapply(coefs, col_names, sum)

             for (vc in valid_classes) if (!(vc %in% names(sums))) sums[[vc]] <- 0

             # Only count vegetation samples (exclude barren)
             if (true_label %in% veg_classes) {
               # Sum of vegetation fractions (excluding barren)
               veg_total <- sum(sapply(veg_classes, function(vc) sums[[vc]]), na.rm = TRUE)

               # Row-normalized fraction for true class (after barren subtraction)
               if (veg_total > 1e-10) {
                 norm_frac <- sums[[true_label]] / veg_total
               } else {
                 norm_frac <- 0
               }

               veg_norm_frac_sums[true_label] <- veg_norm_frac_sums[true_label] + norm_frac
               veg_counts[true_label] <- veg_counts[true_label] + 1L
             }
          }

          # Return avg row-normalized veg fraction correctly predicted across veg classes
          veg_diag_fracs <- ifelse(veg_counts > 0, veg_norm_frac_sums / veg_counts, NA)
          mean_veg_frac <- mean(veg_diag_fracs, na.rm = TRUE)
          if (is.finite(mean_veg_frac)) {
            return(mean_veg_frac)
          } else {
            return(0)
          }
      }

      # 3c. Random Search Loop (User requested ~50 random tries, no generations)
      # Skip if precomputed clusters are provided
      if (!is.null(precomputed_clusters)) {
        cat("[CLUSTER] Using precomputed cluster counts - skipping random search optimization\n")
        best_combo <- precomputed_clusters
        best_mean_score <- NA  # Not computed when using precomputed
        cat(sprintf("    Using precomputed combo: %s\n", paste(names(best_combo), best_combo, sep="=", collapse=", ")))
      } else {
        search_iter <- 50
        if(exists("SEARCH_ITERATIONS")) {
           search_iter <- SEARCH_ITERATIONS
        }

        best_mean_score <- -1
        best_combo <- NULL
        n_valid_combos_evaluated <- 0

        cat(sprintf("      [Random Search] Evaluating %d random cluster combinations...\n", search_iter))

        for(i in 1:search_iter) {
           # Generate Random Combination
           combo <- list()
           for(v in valid_classes) combo[[v]] <- sample(k_ranges[[v]], 1)

           score <- evaluate_config(combo)

           if(score > best_mean_score) {
               best_mean_score <- score
               best_combo <- combo
           }

           # Progress report
           if(i %% 10 == 0 || i == search_iter) {
              cat(sprintf("      [Iter %d/%d] Current Best Score: %.4f\n", i, search_iter, best_mean_score))
           }

           n_valid_combos_evaluated <- n_valid_combos_evaluated + 1
        }

        # If no valid combination was found, return NULL to trigger hard failure
        if (is.null(best_combo)) {
          cat("[ERROR] No valid combinations found during grid search.\n")
          return(NULL)
        }

        cat(sprintf("    Best OOB Score: %.4f using combo: %s\n", best_mean_score, paste(names(best_combo), best_combo, sep="=", collapse=", ")))

        # Store optimal cluster counts globally for reuse (avoids re-optimization after threshold tuning)
        assign("OPTIMAL_CLUSTER_COUNTS", best_combo, envir = globalenv())
        cat(sprintf("[CLUSTER] Stored optimal cluster counts: %s\n", paste(names(best_combo), best_combo, sep="=", collapse=", ")))
      }

      # NOTE: Confusion matrix is computed AFTER threshold optimization (Step 3)
      # to use the final optimized weights. Skipping here.
      cat("\n    [CLUSTER] Confusion matrix will be computed after threshold optimization (Step 3)\n")

      # Return best combination (or NULL if none found)
      return(best_combo)
    }
    
    # -------------------------------------------------------------------------
    # MAIN EXECUTION
    # -------------------------------------------------------------------------
    
    # 1. Run Optimization to get Best K-Combo
    best_combo <- optimize_library(n_boot = 5)
    
    # Fail hard if optimization failed to find a valid combo
    if (is.null(best_combo) || length(best_combo) == 0) {
      stop("[ERROR] optimize_library() failed to find any valid cluster combination. Aborting library build.")
    } else {
      cat(sprintf("[CLUSTER OPT] Best combination found: %s\n", 
                  paste(names(best_combo), "=", sapply(best_combo, as.integer), collapse = ", ")))
    }
    
    # 2. Extract Final Endmembers on FULL Dataset
    # IMPORTANT: Now that we've found optimal cluster sizes using OOB evaluation,
    # we merge OOB data back into the training storage for final endmember extraction
    cat("\n[LIBRARY BUILD] Extracting final endmembers on full dataset...\n")

    # Merge OOB data back into storage for final clustering
    storage_final <- storage
    if (length(storage_oob) > 0 && any(sapply(storage_oob, function(x) !is.null(x) && nrow(x) > 0))) {
      cat("[LIBRARY BUILD] Merging OOB data back into storage for final endmember extraction...\n")
      for (v in names(storage_oob)) {
        if (!is.null(storage_oob[[v]]) && nrow(storage_oob[[v]]) > 0) {
          if (v %in% names(storage_final) && !is.null(storage_final[[v]])) {
            n_before <- nrow(storage_final[[v]])
            storage_final[[v]] <- rbind(storage_final[[v]], storage_oob[[v]])
            n_after <- nrow(storage_final[[v]])
            cat(sprintf("  %s: merged %d OOB samples (total: %d -> %d)\n",
                        v, nrow(storage_oob[[v]]), n_before, n_after))
          } else {
            storage_final[[v]] <- storage_oob[[v]]
            cat(sprintf("  %s: added %d OOB samples (was empty)\n", v, nrow(storage_oob[[v]])))
          }
        }
      }
    } else {
      cat("[LIBRARY BUILD] No OOB data to merge, using storage as-is\n")
    }

    final_lib_cache <- list()

    # Pre-compute Barren Reference(s) for Filtering (support multiple prototypes)
    barren_refs_mat <- NULL
    if ("barren" %in% valid_classes) {
       b_mat <- storage_final[["barren"]]
       if (!is.null(b_mat) && nrow(b_mat) > 0) {
         n_bproto <- if (exists("RAW_BARREN_N_PROTOTYPES") && is.finite(RAW_BARREN_N_PROTOTYPES) && RAW_BARREN_N_PROTOTYPES >= 1) as.integer(RAW_BARREN_N_PROTOTYPES) else 1L
         if (n_bproto > 1 && nrow(b_mat) >= n_bproto) {
           # Use EAR selection to get representative barren prototypes (brightness-invariant)
           pb <- ear_extract_endmembers(b_mat, min(n_bproto, nrow(b_mat)), w_vec)
           if (!is.null(pb) && nrow(pb) > 0) {
             barren_refs_mat <- pb
           } else {
             barren_refs_mat <- matrix(colMeans(b_mat, na.rm = TRUE), nrow = 1)
           }
         } else {
           barren_refs_mat <- matrix(colMeans(b_mat, na.rm = TRUE), nrow = 1)
         }
       }
    }

    for(v in valid_classes) {
       k_opt <- as.numeric(best_combo[[v]])
       veg_mat <- storage_final[[v]]

       # Check availability
       if(nrow(veg_mat) < k_opt) {
          k_opt <- max(1, nrow(veg_mat))
          cat(sprintf("  [WARNING] Reducing k for '%s' to %d (insufficient samples)\n", v, k_opt))
       }

       # EAR-based endmember extraction (brightness-invariant)
       final_endmembers <- ear_extract_endmembers(veg_mat, k_opt, w_vec)

       # Note: Barren filtering now happens before optimization, so final library uses pre-filtered data

       if (!is.null(final_endmembers)) {
          final_lib_cache[[v]] <- final_endmembers
       }
    }
    
    # 3. Construct Final Library Structure (Standard Format)
    res_lib <- list()
    for(v in names(final_lib_cache)) {
       mat <- final_lib_cache[[v]]
       res_lib[[v]] <- list()
       for(i in 1:nrow(mat)) {
          # Construct dummy variant object
          # Re-normalization to T_medoid structure if needed?
          # The rest of the pipeline expects 'vec' (unweighted, z-scored) or 'T' (medoid).
          # Here we provide 'vec'.
          # Use storage_final for n_samples to reflect merged OOB data
          n_total <- if (v %in% names(storage_final)) nrow(storage_final[[v]]) else nrow(storage[[v]])
          res_lib[[v]][[length(res_lib[[v]]) + 1]] <- list(
             vec = mat[i, ],
             id = paste0(v, "_opt_", i),
             n_samples = floor(n_total / nrow(mat)) # Approximate n_samples
          )
       }
    }

    # -------------------------------------------------------------------------
    # Final barren filtering on the library variants
    # -------------------------------------------------------------------------
    if ("barren" %in% names(res_lib) && !is.null(res_lib[["barren"]]) && length(res_lib[["barren"]]) > 0) {
      # Compute barren references from ALL the library's barren variants (support multiple endmembers)
      barren_vecs <- do.call(rbind, lapply(res_lib[["barren"]], function(x) as.numeric(x$vec)))
      # Use all barren endmembers as references (not just the mean)
      barren_refs_mat <- barren_vecs

      cosine_threshold <- BARREN_SIM_THRESHOLD

      for (v in names(res_lib)) {
        if (v == "barren" || is.null(res_lib[[v]]) || length(res_lib[[v]]) == 0) next

        variants <- res_lib[[v]]
        vecs <- do.call(rbind, lapply(variants, function(x) as.numeric(x$vec)))

        # Normalize each row
        vecs_norm <- matrix(0, nrow = nrow(vecs), ncol = ncol(vecs))
        for (i in 1:nrow(vecs)) {
          r <- vecs[i, ]
          nrm <- sqrt(sum(r^2))
          vecs_norm[i, ] <- if (!is.finite(nrm) || nrm == 0) rep(0, length(r)) else r / nrm
        }
        # Normalize all barren references
        barren_refs_norm <- t(apply(barren_refs_mat, 1, function(r) {
          nrm <- sqrt(sum(r^2))
          if (!is.finite(nrm) || nrm == 0) return(rep(0, length(r))) else return(r / nrm)
        }))
        
        sim_mat <- vecs_norm %*% t(barren_refs_norm)
        max_sims <- apply(sim_mat, 1, function(x) suppressWarnings(max(x, na.rm = TRUE)))
        keep_mask <- !(max_sims > cosine_threshold)
        
        # Always keep at least one variant
        if (any(keep_mask)) {
          kept_variants <- variants[keep_mask]
        } else {
          # Keep the least barren-like
          min_idx <- which.min(max_sims)
          kept_variants <- variants[min_idx]
        }
        
        n_before <- length(variants)
        n_after <- length(kept_variants)
        if (n_after < n_before) {
          cat(sprintf("[BARREN FILTER] Final library %s: removed %d/%d variants (similarity > %.2f)\n", 
                      v, n_before - n_after, n_before, cosine_threshold))
        }
        
        res_lib[[v]] <- kept_variants
        best_combo[[v]] <- n_after
      }
      
      # Print updated best combination after final filtering
      cat(sprintf("[CLUSTER OPT] Best combination after barren filtering: %s\n", 
                  paste(names(best_combo), "=", sapply(best_combo, as.integer), collapse = ", ")))
    }

    # Optional: Generate prototype plots (one plot per index/band) showing endmember centers across pentads
    plot_vegetation_prototypes <- function(lib, indices = NULL, out_dir = if (exists("OUT_DIR")) OUT_DIR else ".", prefix = "veg_prototypes", save_png = TRUE, dpi = 150) {
      if (is.null(lib) || length(lib) == 0) return(NULL)
      if (is.null(indices)) {
        if (exists("MESMA_PARAMS") && !is.null(MESMA_PARAMS) && !is.null(MESMA_PARAMS$indices)) indices <- MESMA_PARAMS$indices else indices <- OPTIMAL_INDICES
      }
      if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required for prototype plotting")
      if (!dir.exists(out_dir) && save_png) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      # Build long dataframe for plotting
      rows <- list()
      for (v in names(lib)) {
        variants <- lib[[v]]
        for (var in variants) {
          vec <- as.numeric(var$vec)
          vid <- if (!is.null(var$id)) var$id else if (!is.null(var$variant_id)) var$variant_id else paste0(v, "_unknown")
          n_idx <- length(indices)
          expected_len <- n_idx * TEMPORAL_BUDGET
          if (length(vec) < expected_len) {
            # Skip malformed variants but warn
            warning(sprintf("Skipping variant %s for veg %s: length(vec)=%d != expected=%d", vid, v, length(vec), expected_len))
            next
          }
          for (k in seq_len(n_idx)) {
            idx_name <- indices[k]
            start <- (k-1) * TEMPORAL_BUDGET + 1
            end <- k * TEMPORAL_BUDGET
            vals <- vec[start:end]
            df_tmp <- data.frame(pentad = seq_len(TEMPORAL_BUDGET), value = vals, Veg = v, variant_id = vid, index = idx_name, stringsAsFactors = FALSE)
            rows[[length(rows) + 1]] <- df_tmp
          }
        }
      }

      if (length(rows) == 0) return(NULL)
      proto_df <- do.call(rbind, rows)

      plots <- list()
      for (idx in unique(proto_df$index)) {
        df_idx <- proto_df[proto_df$index == idx, , drop = FALSE]

        # Prepare a veg -> color mapping: prefer user-provided VEG_CALIBRATION_COLORS (case-insensitive),
        # otherwise fall back to RColorBrewer Set1 palette.
        veg_levels <- unique(df_idx$Veg)
        veg_palette <- NULL
        if (exists("VEG_CALIBRATION_COLORS", envir = globalenv())) {
          supplied <- get("VEG_CALIBRATION_COLORS", envir = globalenv())
          # Match by case-insensitive names
          matched <- sapply(veg_levels, function(v) {
            nm <- names(supplied)
            im <- which(tolower(nm) == tolower(v))
            if (length(im) > 0) supplied[im[1]] else NA_character_
          }, USE.NAMES = FALSE)
          if (!all(is.na(matched))) {
            veg_palette <- setNames(matched, veg_levels)
          }
        }

        # Fallback: a simple Set1 palette sized to available veg levels
        if (is.null(veg_palette) || any(is.na(veg_palette))) {
          if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
            stop("RColorBrewer required for default veg palette")
          }
          nveg <- max(3, length(veg_levels))
          brewer_cols <- RColorBrewer::brewer.pal(n = nveg, name = "Set1")
          brewer_cols <- brewer_cols[seq_len(length(veg_levels))]
          names(brewer_cols) <- veg_levels
          if (is.null(veg_palette)) veg_palette <- brewer_cols else {
            na_idx <- which(is.na(veg_palette))
            if (length(na_idx) > 0) veg_palette[na_idx] <- brewer_cols[na_idx]
          }
        }

        # Compute y-axis limits to include all data (including negative values)
        y_min <- min(df_idx$value, na.rm = TRUE)
        y_max <- max(df_idx$value, na.rm = TRUE)
        # Add 5% padding on each side
        y_range <- y_max - y_min
        y_pad <- if (y_range > 0) y_range * 0.05 else 0.1

        p <- ggplot2::ggplot(df_idx, ggplot2::aes(x = pentad, y = value, color = Veg, group = variant_id)) +
             ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.4) +
             ggplot2::geom_line(alpha = 0.9, size = 0.8) +
             ggplot2::labs(title = sprintf("Vegetation prototypes: %s", idx), x = "Pentad", y = sprintf("%s (endmember center)", idx)) +
             ggplot2::theme_minimal() +
             ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
             ggplot2::scale_color_manual(values = veg_palette) +
             ggplot2::coord_cartesian(ylim = c(y_min - y_pad, y_max + y_pad))
        plots[[idx]] <- p
        if (save_png) {
          fn <- file.path(out_dir, sprintf("%s_%s.png", prefix, idx))
          ggplot2::ggsave(filename = fn, plot = p, width = 8, height = 4, dpi = dpi)
        }
      }
      invisible(plots)
    }

    if (exists("GENERATE_PROTO_PLOTS") && isTRUE(GENERATE_PROTO_PLOTS)) {
      tryCatch({
        # Use params$indices (passed to this function) instead of MESMA_PARAMS which may not exist yet
        plot_indices <- if (!is.null(params) && !is.null(params$indices)) params$indices else indices
        plot_vegetation_prototypes(res_lib, indices = plot_indices, out_dir = file.path(if (exists("OUT_DIR")) OUT_DIR else ".", "prototype_plots"))
        cat(sprintf("[NOTICE] Generated prototype plots to %s\n", file.path(if (exists("OUT_DIR")) OUT_DIR else ".", "prototype_plots")))
      }, error = function(e) {
        cat(sprintf("[WARN] Failed to generate prototype plots: %s\n", e$message))
      })
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
        cat(sprintf("[INFO] MESMA pruning requested: %d/%d features zeroed (%.1f%%)\n", n_zero_glob, length(feature_weights), 100*frac_zero_glob))

        # Helpful diagnostics: summarize which indices were effectively ignored (all pentads zero-weighted)
        if (exists("MESMA_PARAMS") && !is.null(MESMA_PARAMS) &&
            !is.null(MESMA_PARAMS$indices) &&
            length(feature_weights) == (TEMPORAL_BUDGET * length(MESMA_PARAMS$indices))) {
          idxs <- MESMA_PARAMS$indices
          bins <- TEMPORAL_BUDGET
          zero_by_idx <- vapply(seq_along(idxs), function(k) {
            rng <- ((k - 1) * bins + 1):(k * bins)
            sum(feature_weights[rng] == 0, na.rm = TRUE)
          }, numeric(1))
          names(zero_by_idx) <- idxs

          all_zero <- names(zero_by_idx)[zero_by_idx >= bins]
          if (length(all_zero) > 0) {
            cat(sprintf("[INFO] Indices with ALL %d pentads zero-weighted: %s\n", bins, paste(all_zero, collapse = ", ")))
          }
        }

        if (frac_zero_glob <= PRUNE_ZERO_WEIGHT_MAX_FRAC && (length(feature_weights) - n_zero_glob) >= PRUNE_ZERO_MIN_FEATURES) {
          keep_idx_global_w <- which(!zero_mask_glob)
          cat(sprintf("[INFO] Will prune %d features for MESMA, keeping %d features\n", n_zero_glob, length(keep_idx_global_w)))
        } else {
          cat(sprintf("[WARN] Skipping MESMA pruning: zero fraction %.2f exceeds max allowed %.2f or resulting features < %d\n", frac_zero_glob, PRUNE_ZERO_WEIGHT_MAX_FRAC, PRUNE_ZERO_MIN_FEATURES))
        }
      }
    }

    # Precompute normalized barren references (after any global pruning) for similarity filtering
    barren_refs_norm <- NULL
    if ("barren" %in% names(mesma_lib) && !is.null(mesma_lib[["barren"]])) {
      variants_barren_raw <- mesma_lib[["barren"]]
      if (!is.list(variants_barren_raw) || (!is.null(variants_barren_raw$vec) && !is.list(variants_barren_raw$vec))) {
        variants_barren <- list(variants_barren_raw)
      } else {
        variants_barren <- variants_barren_raw
      }

      barr_vecs <- lapply(variants_barren, function(x) as.numeric(if (is.list(x) && !is.null(x$vec)) x$vec else x))
      B <- do.call(rbind, barr_vecs)

      if (!is.null(keep_idx_global_w) && length(keep_idx_global_w) > 0 && max(keep_idx_global_w) <= ncol(B)) {
        B <- B[, keep_idx_global_w, drop = FALSE]
      }

      barren_refs_norm <- t(apply(B, 1, function(r) {
        nrm <- sqrt(sum(r^2))
        if (!is.finite(nrm) || nrm == 0) return(rep(0, length(r))) else return(r / nrm)
      }))

      # Drop empty rows
      if (!is.null(barren_refs_norm) && nrow(barren_refs_norm) > 0) {
        keep_barr <- rowSums(abs(barren_refs_norm)) > 0
        barren_refs_norm <- barren_refs_norm[keep_barr, , drop = FALSE]
        if (nrow(barren_refs_norm) == 0) barren_refs_norm <- NULL
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
            cat(sprintf("[INFO] Pruned MESMA features for veg '%s': kept %d columns\n", v, ncol(M)))
          } else {
            cat(sprintf("[WARN] MESMA prune requested but index mismatch for veg '%s'; skipping prune\n", v))
          }
        }

        M_norm <- t(apply(M, 1, function(r) {
          nrm <- sqrt(sum(r^2))
          if (!is.finite(nrm) || nrm == 0) return(rep(0, length(r))) else return(r / nrm)
        }))

        # Note: Barren similarity filtering now happens during cluster evaluation (not here)
        # This ensures accuracy scores reflect performance of filtered endmembers

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
  expected <- n_indices * TEMPORAL_BUDGET

  if (wlen == expected) {
    w_mat <- matrix(params$weights, nrow = TEMPORAL_BUDGET, ncol = n_indices)
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

# =============================================================================
# SPLINE-BASED ENDMEMBER LIBRARY (Alternative Mode)
# =============================================================================
# When USE_SPLINE_ENDMEMBERS = TRUE, this function builds a smooth spline for
# each class/index combination. The spline captures the phenological trajectory
# and allows continuous evaluation at any DOY.
# =============================================================================

build_spline_library <- function(df_train, indices, allowed_veg, spar = SPLINE_SPAR, norm_params = NULL) {
  cat("\n=== Building Spline Library (One Curve Per Class Per Index) ===\n")
  
  # Normalize training data if norm_params provided
  if (!is.null(norm_params) && !is.null(norm_params$means) && !is.null(norm_params$sds)) {
    cat("[SPLINE] Normalizing training data using provided params...\n")
    cat(sprintf("[SPLINE] Available means names: %s\n", paste(head(names(norm_params$means), 5), collapse=", ")))
    cat(sprintf("[SPLINE] Available sds names: %s\n", paste(head(names(norm_params$sds), 5), collapse=", ")))
    n_normalized <- 0
    for (idx in indices) {
      if (!idx %in% names(df_train)) next
      # Access by name, not by numeric index
      mu <- norm_params$means[idx]
      sigma <- norm_params$sds[idx]
      if (is.null(mu) || is.na(mu) || !is.finite(mu)) {
        cat(sprintf("[SPLINE WARNING] No mean for index '%s', using 0\n", idx))
        mu <- 0
      }
      if (is.null(sigma) || is.na(sigma) || !is.finite(sigma) || sigma < 1e-8) sigma <- 1e-8
      df_train[[idx]] <- (df_train[[idx]] - mu) / sigma
      n_normalized <- n_normalized + 1
    }
    cat(sprintf("[SPLINE] Successfully normalized %d indices.\n", n_normalized))
  } else {
    cat("[SPLINE] No normalization params provided - fitting splines on raw data.\n")
  }
  
  classes <- unique(c("barren", allowed_veg))
  library_splines <- list()
  
  # Prepare plotting data for diagnostics
  plot_data_list <- list()
  
  for (veg in classes) {
    cat(sprintf("Fitting splines for class: %s\n", veg))
    sub <- df_train[tolower(df_train$Veg) == veg, ]
    
    if (nrow(sub) < 10) {
      cat(sprintf("  [WARNING] Not enough samples (%d) for %s. Using constant mean.\n", nrow(sub), veg))
      # Create constant functions for all indices
      class_splines <- list()
      for (idx in indices) {
        if (!idx %in% names(sub)) next
        mean_val <- mean(sub[[idx]], na.rm = TRUE)
        if (!is.finite(mean_val)) mean_val <- 0
        class_splines[[idx]] <- list(
          predict = local({
            mv <- mean_val
            function(new_doy) rep(mv, length(new_doy))
          }),
          n_samples = nrow(sub),
          type = "constant"
        )
      }
      library_splines[[veg]] <- class_splines
      next
    }
    
    class_splines <- list()
    n_samples_veg <- nrow(sub)
    
    for (idx in indices) {
      if (!idx %in% names(sub)) next
      
      # Prepare data: DOY and Value
      d <- sub$doy
      v <- sub[[idx]]
      mask <- is.finite(d) & is.finite(v)
      d <- d[mask]
      v <- v[mask]
      
      if (length(d) < 5) {
        mean_val <- mean(v, na.rm = TRUE)
        if (!is.finite(mean_val)) mean_val <- 0
        class_splines[[idx]] <- list(
          predict = local({
            mv <- mean_val
            function(new_doy) rep(mv, length(new_doy))
          }),
          n_samples = length(d),
          type = "constant"
        )
        next
      }
      
      # Handle circularity: Pad data (Year-1, Year, Year+1) to ensure smooth wrapping
      # DOY range: 1 to 365
      d_aug <- c(d - 365, d, d + 365)
      v_aug <- c(v, v, v)
      
      # Fit smoothing spline
      fit <- tryCatch({
        smooth.spline(x = d_aug, y = v_aug, spar = spar)
      }, error = function(e) {
        warning(paste("Spline fit failed for", veg, idx, ":", e$message))
        NULL
      })

      if (is.null(fit)) {
        mean_val <- mean(v, na.rm = TRUE)
        if (!is.finite(mean_val)) mean_val <- 0
        class_splines[[idx]] <- list(
          predict = local({
            mv <- mean_val
            function(new_doy) rep(mv, length(new_doy))
          }),
          n_samples = length(d),
          type = "constant"
        )
        next
      }
      
      # Create a closure function that predicts for new DOYs
      # Store fit object in the closure environment
      spline_predict <- local({
        fit_obj <- fit
        function(new_doy) {
          predict(fit_obj, new_doy)$y
        }
      })
      
      class_splines[[idx]] <- list(
        predict = spline_predict,
        fit = fit,
        n_samples = length(d),
        type = "spline"
      )
      
      # For plotting/debug
      pred_doy <- 1:365
      pred_val <- spline_predict(pred_doy)
      plot_data_list[[paste(veg, idx, sep = "_")]] <- data.frame(
        Veg = veg, Index = idx, DOY = pred_doy, Value = pred_val, stringsAsFactors = FALSE
      )
    }
    
    library_splines[[veg]] <- class_splines
    cat(sprintf("  Fitted %d indices for %s (n_samples=%d)\n", length(class_splines), veg, n_samples_veg))
  }
  
  # Plot the splines
  if (length(plot_data_list) > 0) {
    all_plots <- do.call(rbind, plot_data_list)
    if (!is.null(all_plots) && nrow(all_plots) > 0) {
      tryCatch({
        p <- ggplot(all_plots, aes(x = DOY, y = Value, color = Veg)) +
          geom_line(size = 1) +
          facet_wrap(~Index, scales = "free_y") +
          theme_minimal() +
          labs(title = "Class Phenology Splines", x = "Day of Year", y = "Index Value") +
          coord_cartesian(ylim = c(0, NA))
        out_path <- file.path(OUT_DIR, "spline_library_curves.png")
        if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
        ggsave(out_path, p, width = 14, height = 10, dpi = 150)
        cat(sprintf("[SPLINE] Saved spline library plot to %s\n", out_path))
      }, error = function(e) {
        cat(sprintf("[WARN] Failed to save spline plot: %s\n", e$message))
      })
    }
  }
  
  return(library_splines)
}

# Spline-based unmixing solver: FCLS with sum-to-one constraint
solve_weights_fcls_spline <- function(E, y, lambda = SPLINE_LAMBDA) {
  if (is.null(E) || ncol(E) < 1 || length(y) == 0) return(rep(0, ncol(E)))
  
  # Handle NA values
  valid_mask <- is.finite(y)
  if (sum(valid_mask) < 2) return(rep(1/ncol(E), ncol(E)))
  
  E_valid <- E[valid_mask, , drop = FALSE]
  y_valid <- y[valid_mask]
  
  # --- DEDICATED QP SOLVER (ALWAYS used) ---
  n_endmembers <- ncol(E_valid)
  Dmat <- 2 * crossprod(E_valid)
  dvec <- 2 * crossprod(E_valid, y_valid)
  
  # Scale QP to prevent inconsistent constraints
  qp_scale <- mean(diag(Dmat))
  if(is.na(qp_scale) || qp_scale < 1e-12) qp_scale <- 1.0
  
  Dmat <- Dmat / qp_scale
  dvec <- dvec / qp_scale
  
  # Add stronger ridge
  ridge <- 1e-6
  Dmat <- Dmat + diag(n_endmembers) * ridge
  # Ensure symmetric
  Dmat <- (Dmat + t(Dmat)) / 2
  
  # Constraints: sum=1, x>=0
  Amat <- cbind(rep(1, n_endmembers), diag(n_endmembers))
  bvec <- c(1, rep(0, n_endmembers))

  res <- tryCatch({
    quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  }, error = function(e) {
    NULL
  })
  
  if (is.null(res)) return(rep(1/n_endmembers, n_endmembers))
  
  w <- res$solution
  w[!is.finite(w)] <- 0
  w[w < 0] <- 0
  
  # Normalize to sum to 1 strictly
  if (sum(w) > 0) w <- w / sum(w) else w <- rep(1/n_endmembers, n_endmembers)
  return(w)
}

# Fit a single observation (row) using spline library
fit_observation_spline <- function(obs_row, doy, spline_lib, indices, params = NULL) {
  # obs_row: named vector with index values
  # doy: day of year for this observation
  # spline_lib: output from build_spline_library
  # indices: vector of index names to use
  # params: normalization parameters (means, sds) - if NULL, use raw values
  
  classes <- names(spline_lib)
  n_classes <- length(classes)
  
  if (is.na(doy) || !is.finite(doy)) {
    return(list(weights = rep(NA, n_classes), rmse = NA, classes = classes))
  }
  
  # Build endmember matrix E (rows = indices, cols = classes)
  E <- matrix(0, nrow = length(indices), ncol = n_classes)
  colnames(E) <- classes
  rownames(E) <- indices
  
  y <- numeric(length(indices))
  names(y) <- indices
  valid_idx_mask <- rep(TRUE, length(indices))
  
  for (k in seq_along(indices)) {
    idx <- indices[k]
    
    # Get observed value
    val <- if (idx %in% names(obs_row)) obs_row[[idx]] else NA
    if (!is.finite(val)) {
      valid_idx_mask[k] <- FALSE
      next
    }
    
    # Normalize if params provided - access by index name, not numeric position
    if (!is.null(params) && !is.null(params$means) && !is.null(params$sds)) {
      mu <- params$means[idx]
      sigma <- params$sds[idx]
      if (!is.finite(mu)) mu <- 0
      if (!is.finite(sigma) || sigma < 1e-8) sigma <- 1e-8
      val <- (val - mu) / sigma
    }
    y[k] <- val
    
    # Evaluate each class's spline at this DOY
    for (c_idx in seq_along(classes)) {
      cls <- classes[c_idx]
      spline_obj <- spline_lib[[cls]][[idx]]
      
      if (!is.null(spline_obj) && !is.null(spline_obj$predict)) {
        pred_val <- spline_obj$predict(doy)
        # NOTE: Splines are already fitted on normalized data, so do NOT normalize again
        E[k, c_idx] <- pred_val
      } else {
        valid_idx_mask[k] <- FALSE
      }
    }
  }
  
  # Filter to valid indices
  if (sum(valid_idx_mask) < 2) {
    # DEBUG: Log when we have insufficient valid indices
    if (exists("TESTING_MODE") && isTRUE(TESTING_MODE)) {
      cat(sprintf("[SPLINE OBS DEBUG] DOY=%d: Only %d valid indices (need >=2), skipping\n", doy, sum(valid_idx_mask)))
    }
    return(list(weights = rep(NA, n_classes), rmse = NA, classes = classes))
  }
  
  E_curr <- E[valid_idx_mask, , drop = FALSE]
  y_curr <- y[valid_idx_mask]
  
  # DEBUG: Print sample of E and y to understand the unmixing inputs
  if (exists("SPLINE_DEBUG_COUNTER", envir = .GlobalEnv)) {
    counter <- get("SPLINE_DEBUG_COUNTER", envir = .GlobalEnv)
  } else {
    counter <- 0
  }
  if (counter < 3) {
    cat(sprintf("\n[SPLINE UNMIX DEBUG] DOY=%d, n_valid=%d\n", doy, sum(valid_idx_mask)))
    cat(sprintf("  y (observation): %s\n", paste(sprintf("%.3f", head(y_curr, 5)), collapse=", ")))
    cat("  E (endmember predictions per class):\n")
    for (ci in seq_along(classes)) {
      cat(sprintf("    %s: %s\n", classes[ci], paste(sprintf("%.3f", head(E_curr[, ci], 5)), collapse=", ")))
    }
    assign("SPLINE_DEBUG_COUNTER", counter + 1, envir = .GlobalEnv)
  }
  
  # Solve unmixing
  w <- solve_weights_fcls_spline(E_curr, y_curr)
  
  # DEBUG: Print the resulting weights
  if (counter < 3) {
    cat(sprintf("  Weights: %s\n", paste(sprintf("%s=%.3f", classes, w), collapse=", ")))
  }
  
  # Calculate predicted values and residuals
  pred <- E_curr %*% w
  residuals_valid <- y_curr - as.numeric(pred)
  rmse <- sqrt(mean(residuals_valid^2, na.rm = TRUE))
  
  # Construct full residual vector with NAs for invalid indices
  residuals_full <- rep(NA_real_, length(indices))
  residuals_full[valid_idx_mask] <- residuals_valid
  names(residuals_full) <- indices
  
  return(list(weights = w, rmse = rmse, classes = classes, residuals = residuals_full, 
              E = E_curr, y = y_curr, valid_mask = valid_idx_mask))
}

# Fit all observations for a location-year using spline library
fit_one_task_spline <- function(task_data, spline_lib, indices, params, loc, yr) {
  n_obs <- nrow(task_data)
  classes <- names(spline_lib)
  n_classes <- length(classes)
  
  # DEBUG: Log function entry and write diagnostic file (first 3 calls only)
  if (!exists(".SPLINE_TASK_COUNTER", envir = .GlobalEnv)) {
    assign(".SPLINE_TASK_COUNTER", 0L, envir = .GlobalEnv)
  }
  counter <- get(".SPLINE_TASK_COUNTER", envir = .GlobalEnv)
  if (counter < 3) {
    # Write to a fixed file to bypass sink suppression
    debug_file <- "spline_debug_output.txt"
    con <- file(debug_file, open = if (counter == 0) "wt" else "at")
    writeLines(sprintf("\n=== SPLINE TASK DEBUG: loc=%s yr=%d n_obs=%d ===", loc, yr, n_obs), con)
    writeLines(sprintf("Indices (%d): %s", length(indices), paste(head(indices, 5), collapse=", ")), con)
    writeLines(sprintf("Classes: %s", paste(classes, collapse=", ")), con)
    
    # Test first observation
    if (n_obs > 0) {
      test_doy <- task_data$doy[1]
      test_row <- as.list(task_data[1, ])
      writeLines(sprintf("\nTest observation DOY=%d:", test_doy), con)
      
      # Show raw values for first 3 indices
      for (idx in head(indices, 3)) {
        raw_val <- if (idx %in% names(test_row)) test_row[[idx]] else NA
        mu <- if (!is.null(params$means) && idx %in% names(params$means)) params$means[idx] else NA
        sd <- if (!is.null(params$sds) && idx %in% names(params$sds)) params$sds[idx] else NA
        norm_val <- if (!is.na(raw_val) && !is.na(mu) && !is.na(sd) && sd > 0) (raw_val - mu) / sd else NA
        writeLines(sprintf("  %s: raw=%.4f, mu=%.4f, sd=%.4f, normalized=%.4f", idx, raw_val, mu, sd, norm_val), con)
        
        # Show spline predictions for each class
        for (cls in classes) {
          spline_obj <- spline_lib[[cls]][[idx]]
          if (!is.null(spline_obj) && !is.null(spline_obj$predict)) {
            pred_val <- spline_obj$predict(test_doy)
            writeLines(sprintf("    -> %s spline prediction: %.4f", cls, pred_val), con)
          }
        }
      }
    }
    close(con)
    assign(".SPLINE_TASK_COUNTER", counter + 1L, envir = .GlobalEnv)
  }
  
  # Pre-allocate results
  coefs <- matrix(NA, nrow = n_obs, ncol = n_classes)
  colnames(coefs) <- classes
  rmses <- numeric(n_obs)
  doys <- numeric(n_obs)
  
  # Store per-observation residuals for bootstrap
  obs_residuals <- vector("list", n_obs)
  obs_E <- vector("list", n_obs)  # Store endmember matrices
  obs_y <- vector("list", n_obs)  # Store observed values
  obs_valid_mask <- vector("list", n_obs)  # Store valid index masks
  
  # Get lat/lon from task_data if available
  lat_val <- if ("lat" %in% names(task_data)) task_data$lat[1] else NA
  lon_val <- if ("lon" %in% names(task_data)) task_data$lon[1] else NA
  
  for (i in seq_len(n_obs)) {
    doy <- task_data$doy[i]
    doys[i] <- doy
    obs_row <- as.list(task_data[i, ])
    
    result <- fit_observation_spline(obs_row, doy, spline_lib, indices, params)
    
    if (!is.null(result) && !all(is.na(result$weights))) {
      coefs[i, ] <- result$weights
      rmses[i] <- result$rmse
      obs_residuals[[i]] <- result$residuals
      obs_E[[i]] <- result$E
      obs_y[[i]] <- result$y
      obs_valid_mask[[i]] <- result$valid_mask
    }
  }
  
  # DEBUG: Check if any coefficients were computed
  n_valid_rows <- sum(apply(coefs, 1, function(r) !all(is.na(r))))
  if (n_valid_rows == 0) {
    cat(sprintf("[SPLINE DEBUG] loc=%s yr=%d: NO valid coefficients computed from %d observations!\n", loc, yr, n_obs))
  } else if (isTRUE(TESTING_MODE) || n_obs <= 5) {
    cat(sprintf("[SPLINE DEBUG] loc=%s yr=%d: %d/%d observations produced valid coefficients\n", loc, yr, n_valid_rows, n_obs))
    cat(sprintf("[SPLINE DEBUG] Sample coefs (first valid): %s\n", paste(sprintf("%.3f", coefs[which(!is.na(coefs[,1]))[1], ]), collapse=", ")))
  }
  
  # Aggregate to location-year level: mean fractions across observations
  mean_coefs <- colMeans(coefs, na.rm = TRUE)
  mean_rmse <- mean(rmses, na.rm = TRUE)
  
  # DEBUG: Check aggregated coefficients
  if (all(is.na(mean_coefs)) || sum(mean_coefs, na.rm = TRUE) == 0) {
    cat(sprintf("[SPLINE DEBUG] loc=%s yr=%d: mean_coefs all NA or zero! mean_coefs=%s\n", loc, yr, paste(sprintf("%.4f", mean_coefs), collapse=", ")))
  }
  
  # Normalize mean_coefs to sum to 1
  if (sum(mean_coefs, na.rm = TRUE) > 0) {
    mean_coefs <- mean_coefs / sum(mean_coefs, na.rm = TRUE)
  }
  
  # Build coefficient dataframe (same structure as standard MESMA output)
  coef_rows <- lapply(seq_along(classes), function(ci) {
    data.frame(
      location_id = loc,
      pheno_year = yr,
      lat = lat_val,
      lon = lon_val,
      Veg = classes[ci],
      variant_id = paste0(classes[ci], "_spline"),
      coef = mean_coefs[ci],
      rmse = mean_rmse,
      coef_025 = NA_real_,
      coef_975 = NA_real_,
      coef_sd = sd(coefs[, ci], na.rm = TRUE),
      interval = NA_real_,
      n_obs = n_obs,
      inseparable_variant_flag = FALSE,
      inseparable_variant_details = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  
  coef_df <- do.call(rbind, coef_rows)
  
  # Build per-observation results for detailed output
  obs_results <- cbind(task_data[, c("location_id", "pheno_year", "doy", "date")], 
                       as.data.frame(coefs), 
                       RMSE = rmses)
  
  # Diagnostics
  diag_df <- data.frame(
    location_id = loc,
    pheno_year = yr,
    n_obs = n_obs,
    mean_rmse = mean_rmse,
    stringsAsFactors = FALSE
  )
  for (ci in seq_along(classes)) {
    diag_df[[paste0("frac_", classes[ci])]] <- mean_coefs[ci]
  }
  
  # Concatenate all observation residuals into a single vector for bootstrap
  # Order by DOY to preserve temporal autocorrelation structure
  doy_order <- order(doys, na.last = TRUE)
  residuals_concat <- unlist(obs_residuals[doy_order])
  residuals_concat <- residuals_concat[!is.na(residuals_concat)]
  
  return(list(
    coef_df = coef_df,
    obs_results = obs_results,
    diagnostics = diag_df,
    mean_coefs = mean_coefs,
    mean_rmse = mean_rmse,
    # For bootstrap
    residuals = residuals_concat,
    obs_residuals = obs_residuals,
    obs_E = obs_E,
    obs_y = obs_y,
    obs_valid_mask = obs_valid_mask,
    doys = doys,
    obs_coefs = coefs,
    spline_lib = spline_lib,
    indices = indices,
    params = params
  ))
}

  cat("Building MESMA library with Z-score/PCA/LDA weighting...\n")


  # ==========================================================================
  # OOB-OPTIMIZED PCA-LDA THRESHOLD SELECTION
  # Step 1: Train initial PCA-LDA on df_train_model (excluding OOB tuning set)
  # Step 2: Find optimal threshold for zeroing weights using OOB data
  # Step 3: Remove variables with zero weights and re-run PCA-LDA
  # ==========================================================================

  # Create multi-class classification data: Barren vs Veg A vs Veg B ...
  # This maximizes separation between ALL classes in the feature space.
  # Use df_train_model for initial fitting if OOB split was performed
  train_for_pipeline <- NULL
  if (exists("df_train_model") && !is.null(df_train_model) && nrow(df_train_model) > 0) {
    train_for_pipeline <- df_train_model
    cat("[PCA/LDA] Using df_train_model (OOB holdout excluded) for initial PCA-LDA fitting\n")
  } else if (exists("df_train") && !is.null(df_train) && nrow(df_train) > 0) {
    train_for_pipeline <- df_train
    cat("[PCA/LDA] Using df_train for PCA-LDA fitting (no OOB split)\n")
  } else if (exists("df") && !is.null(df) && nrow(df) > 0) {
    train_for_pipeline <- df
  } else {
    stop("No training data available for PCA/LDA training pipeline")
  }

  # Exclude validation locations
  if (exists("df_validation") && !is.null(df_validation) && nrow(df_validation) > 0) {
    val_locs <- unique(df_validation$location_id)
    before_n <- nrow(train_for_pipeline)
    train_for_pipeline <- train_for_pipeline[!(train_for_pipeline$location_id %in% val_locs), , drop = FALSE]
    after_n <- nrow(train_for_pipeline)
    cat(sprintf("[PCA/LDA] Excluded %d rows belonging to %d validation locations from PCA/LDA training (rows: %d -> %d)\n", before_n - after_n, length(val_locs), before_n, after_n))
  } else {
    cat(sprintf("[PCA/LDA] Training rows used for PCA/LDA: %d\n", nrow(train_for_pipeline)))
  }

  multi_class_data <- dplyr::mutate(train_for_pipeline, target_class = tolower(as.character(Veg)))
  multi_class_data <- dplyr::filter(multi_class_data, !is.na(target_class) & target_class != "")
  multi_class_data <- dplyr::select(multi_class_data, dplyr::all_of(c("location_id", "pheno_year", "date", "doy", "Veg", "target_class", avail)))

  cat("Training feature pipeline (Multi-Class: Barren vs Veg A vs Veg B...)\n")
  cat(sprintf("[NOTICE] Using feature set for training: %s\n", paste(avail, collapse=", ")))

  class_counts <- table(multi_class_data$target_class)
  cat("[MESMA] Training data class distribution:\n")
  print(class_counts)

  # Step 1: Initial PCA-LDA training on model data (excluding OOB)
  cat("\n=== STEP 1: Initial PCA-LDA Training (on model data, excluding OOB) ===\n")
  MESMA_PARAMS_INITIAL <- train_feature_pipeline(multi_class_data, "target_class", avail)

  if (!is.null(MESMA_PARAMS_INITIAL) && !is.null(MESMA_PARAMS_INITIAL$weights)) {
    MESMA_PARAMS_INITIAL$weights[is.na(MESMA_PARAMS_INITIAL$weights)] <- 0
    MESMA_PARAMS_INITIAL$weights[!is.finite(MESMA_PARAMS_INITIAL$weights)] <- 0
    print_weights_summary("INITIAL_PCA_LDA", MESMA_PARAMS_INITIAL)
  }

  # Step 2: Optimal cluster sizing FIRST (before threshold optimization)
  # This ensures cluster optimization uses the full feature set
  cat("\n=== STEP 2: Optimal Cluster Sizing (before threshold optimization) ===\n")
  cat("Building initial MESMA library with full feature weights for cluster optimization...\n")

  # Use indices from params (includes L2norm if enabled) instead of original avail
  indices_for_library <- if (!is.null(MESMA_PARAMS_INITIAL$indices)) MESMA_PARAMS_INITIAL$indices else avail

  # Build library with initial (unthresholded) weights to determine optimal cluster counts
  mesma_lib_initial <- build_mesma_library_weighted(df_train, indices_for_library, MESMA_PARAMS_INITIAL, ALLOWED_VEG)

  # Store the optimal cluster counts discovered during library building
  if (exists("OPTIMAL_CLUSTER_COUNTS", envir = globalenv())) {
    cat("[CLUSTER] Optimal cluster counts determined from full-feature optimization\n")
  }

  # Step 3: Find optimal threshold using OOB data
  optimal_threshold <- 0
  pruned_indices <- avail
  MESMA_PARAMS <- MESMA_PARAMS_INITIAL

  if (exists("df_train_oob") && !is.null(df_train_oob) && nrow(df_train_oob) > 0 &&
      !is.null(MESMA_PARAMS_INITIAL) && !is.null(MESMA_PARAMS_INITIAL$weights)) {

    cat("\n=== STEP 3: Permutation Importance Feature Selection Using OOB Data ===\n")

    # Prepare OOB data for evaluation
    oob_for_eval <- dplyr::mutate(df_train_oob, target_class = tolower(as.character(Veg)))
    oob_for_eval <- dplyr::filter(oob_for_eval, !is.na(target_class) & target_class != "")

    # Debug: verify OOB and training data are separate
    oob_locs <- unique(oob_for_eval$location_id)
    train_locs <- unique(multi_class_data$location_id)
    overlap_locs <- intersect(oob_locs, train_locs)
    cat(sprintf("[PERMUTATION] OOB locations: %d, Training locations: %d, Overlap: %d\n",
                length(oob_locs), length(train_locs), length(overlap_locs)))
    if (length(overlap_locs) > 0) {
      cat(sprintf("[PERMUTATION] WARNING: Overlapping locations found: %s\n",
                  paste(head(overlap_locs, 5), collapse=", ")))
    }
    cat(sprintf("[PERMUTATION] OOB class distribution: %s\n",
                paste(names(table(oob_for_eval$target_class)), "=", table(oob_for_eval$target_class), collapse=", ")))

    # Check minimum sample size
    n_oob_samples <- nrow(unique(oob_for_eval[, c("location_id", "pheno_year")]))
    if (n_oob_samples < PERMUTATION_MIN_SAMPLES) {
      cat(sprintf("[PERMUTATION] WARNING: Only %d OOB samples (< min=%d). Skipping permutation testing.\n",
                  n_oob_samples, PERMUTATION_MIN_SAMPLES))
      optimal_threshold <- 0
      pruned_indices <- avail
      MESMA_PARAMS <- MESMA_PARAMS_INITIAL
    } else {

      # Use feature list from params for permutation testing
      perm_indices <- if (!is.null(MESMA_PARAMS_INITIAL$indices)) MESMA_PARAMS_INITIAL$indices else avail
      dual_mode <- isTRUE(MESMA_PARAMS_INITIAL$dual_mode)
      l2_normalize <- isTRUE(MESMA_PARAMS_INITIAL$l2_normalize)

      # Permutation importance support
      # We precompute the OOB feature matrix once (rows = location_id×pheno_year traces)
      # and then permute feature values across rows. This tests actual feature importance.
      build_endmember_matrix <- function(mesma_library, unique_classes) {
        cols <- list()
        col_names <- character(0)
        for (cls in names(mesma_library)) {
          if (!(cls %in% unique_classes)) next
          cls_variants <- mesma_library[[cls]]
          if (is.null(cls_variants) || length(cls_variants) == 0) next
          for (variant in cls_variants) {
            vec <- variant$vec
            if (!is.null(vec) && length(vec) > 0) {
              cols[[length(cols) + 1]] <- vec
              col_names <- c(col_names, cls)
            }
          }
        }
        if (length(cols) < 2) stop("[PERMUTATION] ERROR: Insufficient endmembers built from library")
        list(E = do.call(cbind, cols), col_names = col_names)
      }

      build_oob_Y <- function(oob_data, params, n_bins) {
        oob_traces <- unique(oob_data[, c("location_id", "pheno_year", "target_class")])
        if (nrow(oob_traces) == 0) stop("[PERMUTATION] ERROR: No OOB traces found")

        indices <- params$indices
        base_indices <- if (!is.null(params$base_indices)) params$base_indices else indices
        n_base_idx <- length(base_indices)
        dual_mode <- isTRUE(params$dual_mode)
        l2_normalize <- isTRUE(params$l2_normalize)

        oob_vecs <- vector("list", nrow(oob_traces))
        oob_labels <- as.character(oob_traces$target_class)

        for (j in seq_len(nrow(oob_traces))) {
          lid <- oob_traces$location_id[j]
          pyr <- oob_traces$pheno_year[j]
          sub <- oob_data[oob_data$location_id == lid & oob_data$pheno_year == pyr, ]
          # Build from base indices
          mat <- build_pentad_matrix(sub, base_indices)
          if (is.null(mat)) {
            oob_vecs[[j]] <- NULL
            next
          }
          vec_raw <- as.numeric(mat)

          # Build feature vector based on representation mode
          if (dual_mode) {
            vec_l2 <- l2_normalize_perindex(vec_raw, n_base_idx, n_bins)
            vec <- c(vec_raw, vec_l2)
          } else if (l2_normalize) {
            vec <- l2_normalize_perindex(vec_raw, n_base_idx, n_bins)
          } else {
            vec <- vec_raw
          }

          # Z-score values (per-index scalar mean/sd)
          vec_zscored <- vec
          for (k in seq_along(indices)) {
            idx_start <- (k - 1) * n_bins + 1
            idx_end <- k * n_bins
            idx_name <- indices[k]
            param_idx <- which(names(params$means) == idx_name)
            if (length(param_idx) > 0) {
              vec_zscored[idx_start:idx_end] <- (vec[idx_start:idx_end] - params$means[param_idx]) / params$sds[param_idx]
            }
          }
          vec_zscored[!is.finite(vec_zscored)] <- 0
          oob_vecs[[j]] <- vec_zscored
        }

        keep <- vapply(oob_vecs, function(x) !is.null(x) && length(x) > 0, logical(1))
        oob_vecs <- oob_vecs[keep]
        oob_labels <- oob_labels[keep]
        if (length(oob_vecs) == 0) stop("[PERMUTATION] ERROR: No OOB vectors built")
        list(Y = do.call(rbind, oob_vecs), labels = oob_labels)
      }

      compute_score_from_Y <- function(E, col_names, Y, labels, unique_classes, weights) {
        all_coefs <- solve_batch_fcls(E, Y, weights)
        if (is.null(all_coefs)) stop("[PERMUTATION] ERROR: FCLS solver returned NULL")

        veg_classes <- setdiff(unique_classes, "barren")
        veg_norm_frac_sums <- setNames(rep(0, length(veg_classes)), veg_classes)
        veg_counts <- setNames(rep(0L, length(veg_classes)), veg_classes)

        total <- nrow(Y)
        for (j in seq_len(total)) {
          true_cls <- labels[j]
          if (!(true_cls %in% veg_classes)) next
          coefs <- all_coefs[j, ]
          class_sums <- tapply(coefs, col_names, sum)
          for (uc in unique_classes) {
            if (!(uc %in% names(class_sums))) class_sums[[uc]] <- 0
          }
          veg_total <- sum(sapply(veg_classes, function(vc) class_sums[[vc]]), na.rm = TRUE)
          norm_frac <- if (veg_total > 1e-10) class_sums[[true_cls]] / veg_total else 0
          veg_norm_frac_sums[true_cls] <- veg_norm_frac_sums[true_cls] + norm_frac
          veg_counts[true_cls] <- veg_counts[true_cls] + 1L
        }

        veg_diag_fracs <- ifelse(veg_counts > 0, veg_norm_frac_sums / veg_counts, NA)
        mean(veg_diag_fracs, na.rm = TRUE)
      }

      .perm_get_block <- function(index_name, all_index_names, n_bins) {
        pos <- which(all_index_names == index_name)
        if (length(pos) != 1) return(NULL)
        s <- (pos - 1) * n_bins + 1
        e <- pos * n_bins
        list(pos = pos, start = s, end = e)
      }

      # Precompute OOB feature matrix + endmember matrix for permutation importance
      unique_classes <- unique(oob_for_eval$target_class)
      if (length(unique_classes) < 2) {
        stop(sprintf("[PERMUTATION] ERROR: Only %d unique class(es) in OOB data", length(unique_classes)))
      }

      endm <- build_endmember_matrix(mesma_lib_initial, unique_classes)
      oob_built <- build_oob_Y(oob_for_eval, MESMA_PARAMS_INITIAL, TEMPORAL_BUDGET)
      Y_base <- oob_built$Y
      oob_labels <- oob_built$labels

      # Compute baseline score with all features
      cat("[PERMUTATION] Computing baseline score with all features...\n")
      baseline_score <- compute_score_from_Y(endm$E, endm$col_names, Y_base, oob_labels, unique_classes, MESMA_PARAMS_INITIAL$weights)
      cat(sprintf("[PERMUTATION] Baseline score (all features): %.4f\n", baseline_score))

      # === STAGE 1: Test indices ===
      mode_str <- if (dual_mode) "DUAL" else if (l2_normalize) "L2" else "RAW"
      cat(sprintf("\n[PERMUTATION STAGE 1] Testing %d indices (%d permutations each) [mode=%s]...\n",
                  length(perm_indices), PERMUTATION_N_ITER, mode_str))

      index_results <- data.frame(
        index = character(),
        baseline_score = numeric(),
        mean_perm_score = numeric(),
        score_drop = numeric(),
        p_value = numeric(),
        stringsAsFactors = FALSE
      )

      test_index_stage1 <- function(idx_i, index_names, base_weights, Y_base, oob_labels,
                                    endm, unique_classes,
                                    TEMPORAL_BUDGET, PERMUTATION_N_ITER) {
        idx_name <- index_names[idx_i]
        blk <- .perm_get_block(idx_name, index_names, TEMPORAL_BUDGET)
        if (is.null(blk)) {
          return(list(idx_name = idx_name, perm_scores = rep(NA_real_, PERMUTATION_N_ITER), n_failed = PERMUTATION_N_ITER, first_error = "missing block"))
        }

        perm_scores <- rep(NA_real_, PERMUTATION_N_ITER)
        n_failed <- 0L
        first_error <- NA_character_

        for (perm_iter in seq_len(PERMUTATION_N_ITER)) {
          perm_scores[perm_iter] <- tryCatch({
            Y_work <- Y_base
            perm_rows <- sample.int(nrow(Y_base))
            Y_work[, blk$start:blk$end] <- Y_base[perm_rows, blk$start:blk$end, drop = FALSE]
            compute_score_from_Y(endm$E, endm$col_names, Y_work, oob_labels, unique_classes, base_weights)
          }, error = function(e) {
            n_failed <<- n_failed + 1L
            if (is.na(first_error)) first_error <<- conditionMessage(e)
            NA_real_
          })
        }

        list(idx_name = idx_name, perm_scores = perm_scores, n_failed = n_failed, first_error = first_error)
      }

      # Run Stage 1 in parallel or sequentially
      if (isTRUE(PARALLEL_ENABLE) && requireNamespace("future.apply", quietly = TRUE)) {
        # Use a temporary future plan with a larger worker pool for permutation testing only.
        max_cores <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) NA_integer_)
        perm_workers_env <- suppressWarnings(as.integer(Sys.getenv("MESMA_PERMUTATION_WORKERS", unset = NA)))
        perm_workers_req <- if (!is.na(perm_workers_env)) perm_workers_env else PERMUTATION_PARALLEL_WORKERS
        perm_workers_eff <- perm_workers_req
        if (is.finite(max_cores) && !is.na(max_cores)) perm_workers_eff <- min(perm_workers_req, max(1L, max_cores - 1L))
        perm_workers_eff <- max(1L, as.integer(perm_workers_eff))

        cat(sprintf("[PERMUTATION STAGE 1] Running in PARALLEL mode with %d workers (requested %d)...\n",
                    perm_workers_eff, perm_workers_req))
        old_plan <- future::plan()
        tryCatch({
          # Exporting MESMA libraries/params to workers can exceed future's default globals size.
          options(future.globals.maxSize = max(getOption("future.globals.maxSize", 2e9), 12e9))
          future::plan(future::multisession, workers = perm_workers_eff)

          # Submit one future per index
          futures_list <- vector("list", length(perm_indices))
          for (idx_i in seq_along(perm_indices)) {
            idx_name_local <- perm_indices[idx_i]
            local_idx_i <- idx_i
            futures_list[[idx_i]] <- future::future({
              test_index_stage1(local_idx_i, perm_indices, MESMA_PARAMS_INITIAL$weights, Y_base, oob_labels,
                                endm, unique_classes,
                                TEMPORAL_BUDGET, PERMUTATION_N_ITER)
            }, seed = TRUE)
            cat(sprintf("  Submitted index job [%d/%d]: %s\n", idx_i, length(perm_indices), idx_name_local)); flush.console()
          }

          stage1_results <- vector("list", length(futures_list))
          resolved_mask <- rep(FALSE, length(futures_list))
          start_time <- Sys.time()
          last_heartbeat <- start_time
          cat("[PERMUTATION STAGE 1] Submitted all index jobs — streaming results as they finish\n"); flush.console()

          while(!all(resolved_mask)) {
            for (i in seq_along(futures_list)) {
              if (resolved_mask[i]) next
              if (future::resolved(futures_list[[i]])) {
                res <- tryCatch(future::value(futures_list[[i]]), error = function(e) e)
                if (inherits(res, "error")) {
                  cat(sprintf("  [%d/%d] %s: ERROR during permutations — keeping by default\n", i, length(futures_list), perm_indices[i]))
                  stage1_results[[i]] <- list(idx_name = perm_indices[i], perm_scores = rep(NA_real_, PERMUTATION_N_ITER), n_failed = PERMUTATION_N_ITER, first_error = conditionMessage(res))
                } else {
                  perm_scores <- res$perm_scores
                  perm_scores <- perm_scores[is.finite(perm_scores)]
                  p_val <- if (length(perm_scores) == 0) 0 else sum(perm_scores >= baseline_score) / length(perm_scores)
                  decision <- if (p_val < PERMUTATION_ALPHA) "KEEP" else "PRUNE"
                  cat(sprintf("  [%d/%d] %s: p=%.4f -> %s\n", i, length(futures_list), perm_indices[i], p_val, decision)); flush.console()
                  stage1_results[[i]] <- res
                }
                resolved_mask[i] <- TRUE
              }
            }

            # periodic heartbeat so user sees progress even when individual jobs are long
            if (as.numeric(difftime(Sys.time(), last_heartbeat, units = "secs")) > 5) {
              done <- sum(resolved_mask)
              running <- which(!resolved_mask)
              running_names <- if (length(running) > 0) paste(head(perm_indices[running], 6), collapse = ", ") else "-"
              cat(sprintf("  [heartbeat] %d/%d finished — running: %s\n", done, length(futures_list), running_names)); flush.console()
              last_heartbeat <- Sys.time()
            }

            Sys.sleep(0.15)
          }
          elapsed <- difftime(Sys.time(), start_time, units = "secs")
          cat(sprintf("[PERMUTATION STAGE 1] All index jobs completed (%.1fs)\n", as.numeric(elapsed)))
        }, finally = {
          # restore user's plan
          try(future::plan(old_plan), silent = TRUE)
        })
      } else {
        cat("[PERMUTATION STAGE 1] Running in SEQUENTIAL mode...\n")
        stage1_results <- lapply(seq_along(perm_indices), function(idx_i) {
          cat(sprintf("  [%d/%d] Testing index: %s ... ", idx_i, length(perm_indices), perm_indices[idx_i]))
          res <- test_index_stage1(idx_i, perm_indices, MESMA_PARAMS_INITIAL$weights, Y_base, oob_labels,
                                   endm, unique_classes,
                                   TEMPORAL_BUDGET, PERMUTATION_N_ITER)
          perm_scores <- res$perm_scores
          perm_scores <- perm_scores[is.finite(perm_scores)]
          p_val <- if (length(perm_scores) == 0) 0 else sum(perm_scores >= baseline_score) / length(perm_scores)
          decision <- if (p_val < PERMUTATION_ALPHA) "KEEP" else "PRUNE"
          cat(sprintf("p=%.4f -> %s\n", p_val, decision))
          res
        })
      }

      # Process Stage 1 results
      for (res in stage1_results) {
        idx_name <- res$idx_name
        if (is.null(idx_name) || is.na(idx_name) || !(idx_name %in% perm_indices)) next
        perm_scores <- res$perm_scores
        perm_scores <- perm_scores[is.finite(perm_scores)]
        p_val <- if (length(perm_scores) == 0) 0 else sum(perm_scores >= baseline_score) / length(perm_scores)
        mean_perm <- if (length(perm_scores) == 0) NA_real_ else mean(perm_scores, na.rm = TRUE)
        drop_val <- if (is.finite(mean_perm)) baseline_score - mean_perm else NA_real_
        index_results <- rbind(index_results, data.frame(
          index = idx_name,
          baseline_score = baseline_score,
          mean_perm_score = mean_perm,
          score_drop = drop_val,
          p_value = p_val,
          stringsAsFactors = FALSE
        ))
      }

      indices_to_keep_stage1 <- index_results$index[is.finite(index_results$p_value) & index_results$p_value < PERMUTATION_ALPHA]
      indices_to_keep_stage1 <- unique(indices_to_keep_stage1)
      indices_to_prune_stage1_full <- setdiff(perm_indices, indices_to_keep_stage1)

      cat(sprintf("\n[PERMUTATION STAGE 1] Results (alpha=%.3f):\n", PERMUTATION_ALPHA))
      cat(sprintf("  Indices to prune (p >= %.3f): %d\n", PERMUTATION_ALPHA, length(indices_to_prune_stage1_full)))
      if (length(indices_to_prune_stage1_full) > 0) cat(sprintf("    %s\n", paste(indices_to_prune_stage1_full, collapse=", ")))
      cat(sprintf("  Indices kept for pentad testing: %d\n", length(indices_to_keep_stage1)))

      # ALWAYS rebuild a Stage-2 baseline AFTER removing Stage-1 indices
      # (this ensures pentad testing measures incremental value on the reduced feature set)
      pruned_weights_stage1 <- MESMA_PARAMS_INITIAL$weights
      if (length(indices_to_prune_stage1_full) > 0) {
        for (idx_name in indices_to_prune_stage1_full) {
          blk <- .perm_get_block(idx_name, perm_indices, TEMPORAL_BUDGET)
          if (!is.null(blk)) pruned_weights_stage1[blk$start:blk$end] <- 0
        }
      }

      baseline_stage0 <- baseline_score
      cat(sprintf("[PERMUTATION] Baseline before Stage-1 pruning: %.4f\n", baseline_stage0))
      baseline_after_stage1 <- compute_score_from_Y(endm$E, endm$col_names, Y_base, oob_labels, unique_classes, pruned_weights_stage1)
      cat(sprintf("[PERMUTATION] Baseline AFTER Stage-1 pruning: %.4f (delta=%.4f)\n",
                  baseline_after_stage1, baseline_after_stage1 - baseline_stage0))

      # No representation selection stage in single-mode operation.
      repr_prune_full <- character(0)
      pruned_weights_stage2 <- pruned_weights_stage1
      baseline_stage3 <- baseline_after_stage1
      cat(sprintf("[PERMUTATION] Baseline used for pentad tests: %.4f\n", baseline_stage3))

      # === STAGE 3: Test pentads of chosen representations ===
      pentad_results <- data.frame(
        index = character(),
        pentad = integer(),
        feature_name = character(),
        baseline_score = numeric(),
        mean_perm_score = numeric(),
        score_drop = numeric(),
        p_value = numeric(),
        stringsAsFactors = FALSE
      )

      # In single-representation mode, indices_to_keep_stage3 is simply the kept indices.
      indices_to_keep_stage3 <- indices_to_keep_stage1

      if (length(indices_to_keep_stage3) > 0) {
        n_features_stage3 <- length(indices_to_keep_stage3) * TEMPORAL_BUDGET
        cat(sprintf("\n[PERMUTATION STAGE 3] Testing %d features (%d chosen representations × %d pentads) with %d permutations each...\n",
                    n_features_stage3, length(indices_to_keep_stage3), TEMPORAL_BUDGET, PERMUTATION_N_ITER))

        # Helper function for Stage 3: test one feature (index + pentad combination)
        test_feature_stage2 <- function(feature_info, avail, base_weights, Y_base, oob_labels,
                                        endm, unique_classes,
                                        TEMPORAL_BUDGET, PERMUTATION_N_ITER) {
          idx_name <- feature_info$idx_name
          idx_i <- feature_info$idx_i
          pentad_j <- feature_info$pentad_j
          feature_pos <- (idx_i - 1) * TEMPORAL_BUDGET + pentad_j
          feature_name <- sprintf("%s_p%d", idx_name, pentad_j)

          perm_scores <- numeric(PERMUTATION_N_ITER)
          n_failed <- 0L
          first_error <- NA_character_

          for (perm_iter in seq_len(PERMUTATION_N_ITER)) {
            perm_score <- tryCatch({
              Y_work <- Y_base
              perm_rows <- sample.int(nrow(Y_base))
              Y_work[, feature_pos] <- Y_base[perm_rows, feature_pos]
              compute_score_from_Y(endm$E, endm$col_names, Y_work, oob_labels, unique_classes, base_weights)
            }, error = function(e) {
              n_failed <<- n_failed + 1L
              if (is.na(first_error)) first_error <<- conditionMessage(e)
              NA_real_
            })
            perm_scores[perm_iter] <- perm_score
          }

          list(
            idx_name = idx_name,
            pentad_j = pentad_j,
            feature_name = feature_name,
            perm_scores = perm_scores,
            n_failed = n_failed,
            first_error = first_error
          )
        }

        # Build list of all features to test
        features_to_test <- list()
        for (idx_name in indices_to_keep_stage3) {
          idx_i <- which(perm_indices == idx_name)
          for (pentad_j in seq_len(TEMPORAL_BUDGET)) {
            features_to_test[[length(features_to_test) + 1]] <- list(
              idx_name = idx_name,
              idx_i = idx_i,
              pentad_j = pentad_j
            )
          }
        }

        # Run Stage 2 in parallel or sequentially
        if (isTRUE(PARALLEL_ENABLE) && requireNamespace("future", quietly = TRUE) && requireNamespace("future.apply", quietly = TRUE)) {
          # Temporary larger worker pool for permutation stage 2 — submit one batch per index
          max_cores <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) NA_integer_)
          perm_workers_env <- suppressWarnings(as.integer(Sys.getenv("MESMA_PERMUTATION_WORKERS", unset = NA)))
          perm_workers_req <- if (!is.na(perm_workers_env)) perm_workers_env else PERMUTATION_PARALLEL_WORKERS
          perm_workers_eff <- perm_workers_req
          if (is.finite(max_cores) && !is.na(max_cores)) perm_workers_eff <- min(perm_workers_req, max(1L, max_cores - 1L))
          perm_workers_eff <- max(1L, as.integer(perm_workers_eff))

          cat(sprintf("[PERMUTATION STAGE 2] Running in PARALLEL mode with %d workers (requested %d)...\n",
                      perm_workers_eff, perm_workers_req))

          old_plan <- future::plan()
          tryCatch({
            # Exporting MESMA libraries/params to workers can exceed future's default globals size.
            options(future.globals.maxSize = max(getOption("future.globals.maxSize", 2e9), 12e9))
            future::plan(future::multisession, workers = perm_workers_eff)

            idx_batches <- indices_to_keep_stage3
            batch_futures <- vector("list", length(idx_batches))
            for (ii in seq_along(idx_batches)) {
              idx_name_local <- idx_batches[ii]
              idx_i_local <- which(perm_indices == idx_name_local)
              local_idx_name <- idx_name_local
              local_idx_i <- idx_i_local
              batch_futures[[ii]] <- future::future({
                res_list <- vector("list", TEMPORAL_BUDGET)
                for (pentad_j in seq_len(TEMPORAL_BUDGET)) {
                  res_list[[pentad_j]] <- test_feature_stage2(list(idx_name = local_idx_name, idx_i = local_idx_i, pentad_j = pentad_j),
                                                              perm_indices, pruned_weights_stage2, Y_base, oob_labels,
                                                              endm, unique_classes,
                                                              TEMPORAL_BUDGET, PERMUTATION_N_ITER)
                }
                res_list
              }, seed = TRUE)
              cat(sprintf("  Submitted batch [%d/%d]: %s\n", ii, length(idx_batches), idx_name_local))
            }

            # Poll and stream results as index-batches complete
            stage2_results <- list()
            resolved_mask2 <- rep(FALSE, length(batch_futures))
            start2 <- Sys.time()
            last_heartbeat2 <- start2
            cat("[PERMUTATION STAGE 2] Submitted index batches — streaming pentad results as batches finish\n"); flush.console()
            while (!all(resolved_mask2)) {
              for (bi in seq_along(batch_futures)) {
                if (resolved_mask2[bi]) next
                if (future::resolved(batch_futures[[bi]])) {
                  val <- tryCatch(future::value(batch_futures[[bi]]), error = function(e) e)
                  idx_name_local <- idx_batches[bi]
                  if (inherits(val, "error")) {
                    cat(sprintf("  Testing index: %s -> ERROR on worker (marking pentads as failed)\n", idx_name_local))
                    for (pj in seq_len(TEMPORAL_BUDGET)) {
                      stage2_results[[length(stage2_results) + 1]] <- list(idx_name = idx_name_local, pentad_j = pj, feature_name = sprintf("%s_p%d", idx_name_local, pj), perm_scores = rep(NA_real_, PERMUTATION_N_ITER))
                    }
                  } else {
                    for (res in val) {
                      perm_scores <- res$perm_scores[is.finite(res$perm_scores)]
                      if (length(perm_scores) == 0) {
                        err_msg <- if (!is.null(res$first_error) && is.character(res$first_error) && !is.na(res$first_error)) res$first_error else "(unknown error)"
                        cat(sprintf("  Testing index: %s -> pentad %2d: p=NA (failed) -> KEEP | first_error: %s\n",
                                    res$idx_name, res$pentad_j, err_msg)); flush.console()
                      } else {
                        p_val <- sum(perm_scores >= baseline_stage3) / length(perm_scores)
                        decision <- if (p_val < PERMUTATION_PENTAD_ALPHA) "KEEP" else "PRUNE"
                        cat(sprintf("  Testing index: %s -> pentad %2d: p=%.4f -> %s\n", res$idx_name, res$pentad_j, p_val, decision)); flush.console()
                      }
                      stage2_results[[length(stage2_results) + 1]] <- res
                    }
                  }
                  resolved_mask2[bi] <- TRUE
                }
              }

              if (as.numeric(difftime(Sys.time(), last_heartbeat2, units = "secs")) > 5) {
                done2 <- sum(resolved_mask2)
                running2 <- which(!resolved_mask2)
                running_names2 <- if (length(running2) > 0) paste(head(idx_batches[running2], 6), collapse = ", ") else "-"
                cat(sprintf("  [heartbeat] %d/%d batches finished — running: %s\n", done2, length(batch_futures), running_names2))
                last_heartbeat2 <- Sys.time()
              }

              Sys.sleep(0.15)
            }
            elapsed2 <- difftime(Sys.time(), start2, units = "secs")
            cat(sprintf("[PERMUTATION STAGE 2] All index batches completed (%.1fs)\n", as.numeric(elapsed2)))

          }, finally = {
            try(future::plan(old_plan), silent = TRUE)
          })
        } else {
          cat("[PERMUTATION STAGE 2] Running in SEQUENTIAL mode...\n")
          stage2_results <- lapply(seq_along(features_to_test), function(i) {
            feature_info <- features_to_test[[i]]
            if (i == 1 || features_to_test[[i-1]]$idx_name != feature_info$idx_name) {
              cat(sprintf("  Testing index: %s\n", feature_info$idx_name))
            }
            cat(sprintf("    [pentad %2d/%d] ... ", feature_info$pentad_j, TEMPORAL_BUDGET))
            res <- test_feature_stage2(feature_info, perm_indices, pruned_weights_stage2, oob_for_eval,
                                Y_base, oob_labels,
                                endm, unique_classes,
                                TEMPORAL_BUDGET, PERMUTATION_N_ITER)
            # Compute and display p-value immediately
            perm_scores <- res$perm_scores[is.finite(res$perm_scores)]
            if (length(perm_scores) == 0) {
              cat("p=NA (failed) -> KEEP\n")
            } else {
              p_val <- sum(perm_scores >= baseline_stage3) / length(perm_scores)
              decision <- if (p_val < PERMUTATION_PENTAD_ALPHA) "KEEP" else "PRUNE"
              cat(sprintf("p=%.4f -> %s\n", p_val, decision))
            }
            res
          })
        }

        # Process Stage 2 results
        for (res in stage2_results) {
          idx_name <- res$idx_name
          pentad_j <- res$pentad_j
          feature_name <- res$feature_name
          perm_scores <- res$perm_scores

          # Remove NAs
          perm_scores <- perm_scores[is.finite(perm_scores)]
          if (length(perm_scores) == 0) {
            pentad_results <- rbind(pentad_results, data.frame(
              index = idx_name,
              pentad = pentad_j,
              feature_name = feature_name,
              baseline_score = baseline_stage3,
              mean_perm_score = NA_real_,
              score_drop = NA_real_,
              p_value = 0,  # p=0 means highly significant -> keep
              stringsAsFactors = FALSE
            ))
            next
          }

          mean_perm_score <- mean(perm_scores, na.rm = TRUE)
          score_drop <- baseline_stage3 - mean_perm_score
          p_value <- sum(perm_scores >= baseline_stage3) / length(perm_scores)

          pentad_results <- rbind(pentad_results, data.frame(
            index = idx_name,
            pentad = pentad_j,
            feature_name = feature_name,
            baseline_score = baseline_stage3,
            mean_perm_score = mean_perm_score,
            score_drop = score_drop,
            p_value = p_value,
            stringsAsFactors = FALSE
          ))
        }

        # Report summary per index
        for (idx_name in indices_to_keep_stage3) {
          idx_pentad_results <- pentad_results[pentad_results$index == idx_name, ]
          n_sig <- sum(idx_pentad_results$p_value < PERMUTATION_PENTAD_ALPHA, na.rm = TRUE)
          cat(sprintf("  %s: %d/%d pentads significant (p<%.2f)\n",
                      idx_name, n_sig, TEMPORAL_BUDGET, PERMUTATION_PENTAD_ALPHA))
        }
        # Compute baseline after pentad pruning
        pruned_weights_after_pentads <- pruned_weights_stage2
        pentads_pruned <- pentad_results[pentad_results$p_value >= PERMUTATION_PENTAD_ALPHA, ]
        if (nrow(pentads_pruned) > 0) {
          for (r in seq_len(nrow(pentads_pruned))) {
            idx_name <- pentads_pruned$index[r]
            pentad_j <- pentads_pruned$pentad[r]
            idx_i <- which(perm_indices == idx_name)
            if (length(idx_i) == 1) {
              feat_pos <- (idx_i - 1) * TEMPORAL_BUDGET + pentad_j
              pruned_weights_after_pentads[feat_pos] <- 0
            }
          }
        }
        baseline_after_pentads <- compute_score_from_Y(endm$E, endm$col_names, Y_base, oob_labels, unique_classes, pruned_weights_after_pentads)
        cat(sprintf("[PERMUTATION] Baseline AFTER pentad pruning (Stage 2): %.4f (delta from Stage 1=%.4f)\n",
                    baseline_after_pentads, baseline_after_pentads - baseline_stage3))
      } else {
        cat("[PERMUTATION STAGE 3] No indices passed Stage 1/2, skipping pentad testing.\n")
      }

      # Combine results for saving
      perm_results <- rbind(
        data.frame(
          index = index_results$index,
          pentad = NA_integer_,
          feature_name = paste0(index_results$index, "_all"),
          baseline_score = index_results$baseline_score,
          mean_perm_score = index_results$mean_perm_score,
          score_drop = index_results$score_drop,
          p_value = index_results$p_value,
          stage = "index",
          stringsAsFactors = FALSE
        ),
        if (nrow(pentad_results) > 0) {
          data.frame(
            pentad_results,
            stage = "pentad",
            stringsAsFactors = FALSE
          )
        } else {
          data.frame(
            index = character(),
            pentad = integer(),
            feature_name = character(),
            baseline_score = numeric(),
            mean_perm_score = numeric(),
            score_drop = numeric(),
            p_value = numeric(),
            stage = character(),
            stringsAsFactors = FALSE
          )
        }
      )

      # === PRUNING LOGIC ===
      # Stage 1 base indices that failed are completely removed (raw + L2norm)
      # Stage 2 representation selection removes either raw or L2norm for each kept index
      # Stage 3 pentads that failed are zeroed out
      
      cat(sprintf("\n[PERMUTATION] Final Results (alpha=%.3f):\n", PERMUTATION_ALPHA))
      
      # Create pruned weight vector
      pruned_weights <- MESMA_PARAMS_INITIAL$weights
      
      # First, zero out all pentads from Stage 1 pruned indices (full names)
      for (idx_name in indices_to_prune_stage1_full) {
        idx_i <- which(perm_indices == idx_name)
        if (length(idx_i) == 1) {
          idx_start <- (idx_i - 1) * TEMPORAL_BUDGET + 1
          idx_end <- idx_i * TEMPORAL_BUDGET
          pruned_weights[idx_start:idx_end] <- 0
        }
      }

      # Then, remove the non-selected representation (Stage 2)
      if (length(repr_prune_full) > 0) {
        for (idx_name in repr_prune_full) {
          idx_i <- which(perm_indices == idx_name)
          if (length(idx_i) == 1) {
            idx_start <- (idx_i - 1) * TEMPORAL_BUDGET + 1
            idx_end <- idx_i * TEMPORAL_BUDGET
            pruned_weights[idx_start:idx_end] <- 0
          }
        }
      }
      
      # Then, zero out specific pentads from Stage 2 results
      pentads_to_prune <- data.frame(
        index = character(),
        pentad = integer(),
        stringsAsFactors = FALSE
      )
      if (nrow(pentad_results) > 0) {
        pentads_to_prune <- pentad_results[pentad_results$p_value >= PERMUTATION_PENTAD_ALPHA, ]
        for (i in seq_len(nrow(pentads_to_prune))) {
          idx_name <- pentads_to_prune$index[i]
          pentad_j <- pentads_to_prune$pentad[i]
          idx_i <- which(perm_indices == idx_name)
          feature_pos <- (idx_i - 1) * TEMPORAL_BUDGET + pentad_j
          pruned_weights[feature_pos] <- 0
        }
      }
      
      # Check which indices have ALL pentads zeroed
      indices_to_remove_all <- c()
      for (k in seq_along(perm_indices)) {
        idx_start <- (k-1)*TEMPORAL_BUDGET + 1
        idx_end <- k*TEMPORAL_BUDGET
        if (all(pruned_weights[idx_start:idx_end] == 0)) {
          indices_to_remove_all <- c(indices_to_remove_all, perm_indices[k])
        }
      }

      # Determine kept indices for reporting / minimum-feature enforcement
      indices_to_keep <- character(0)
      for (idx_name in perm_indices) {
        idx_pos <- which(perm_indices == idx_name)
        if (length(idx_pos) == 1) {
          idx_start <- (idx_pos - 1) * TEMPORAL_BUDGET + 1
          idx_end <- idx_pos * TEMPORAL_BUDGET
          if (any(pruned_weights[idx_start:idx_end] != 0)) {
            indices_to_keep <- c(indices_to_keep, idx_name)
          }
        }
      }
      indices_to_keep <- unique(indices_to_keep)
      
      features_to_keep_count <- sum(pruned_weights > 0)
      
      cat(sprintf("  Total indices removed: %d\n", length(indices_to_remove_all)))
      if (length(indices_to_remove_all) > 0) {
        cat(sprintf("    %s\n", paste(indices_to_remove_all, collapse=", ")))
      }
      cat(sprintf("  Features kept (non-zero weights): %d/%d\n",
                  features_to_keep_count, length(MESMA_PARAMS_INITIAL$weights)))
      cat(sprintf("  Indices kept: %d/%d\n", length(indices_to_keep), length(perm_indices)))

      # Enforce minimum feature count
      if (length(indices_to_keep) < PRUNE_ZERO_MIN_FEATURES) {
        cat(sprintf("[PERMUTATION] WARNING: Only %d indices kept (< min=%d). Keeping all features and ignoring pruning decisions.\n",
                    length(indices_to_keep), PRUNE_ZERO_MIN_FEATURES))
        # IMPORTANT: ensure downstream pruning logic does not zero-out all weights.
        indices_to_prune_stage1_full <- character(0)
        repr_prune_full <- character(0)
        if (exists("pentads_to_prune") && is.data.frame(pentads_to_prune)) {
          pentads_to_prune <- pentads_to_prune[0, , drop = FALSE]
        }
        pruned_indices <- avail
        optimal_threshold <- 0
      } else {
        pruned_indices <- indices_to_keep
        optimal_threshold <- 0  # No threshold-based pruning in permutation approach
      }

      # Save permutation results
      write.csv(perm_results, file = "permutation_importance_results.csv", row.names = FALSE)
      cat("[PERMUTATION] Saved results to permutation_importance_results.csv\n")

      # Helper to apply pruning rules to a given weight vector
      apply_pruning_rules <- function(weights, index_names, n_bins, indices_to_prune_stage1_full, repr_prune_full, pentads_to_prune) {
        if (is.null(weights) || length(weights) == 0) return(weights)

        # Stage 1: prune entire indices
        for (idx_name in indices_to_prune_stage1_full) {
          idx_i <- which(index_names == idx_name)
          if (length(idx_i) == 1) {
            idx_start <- (idx_i - 1) * n_bins + 1
            idx_end <- idx_i * n_bins
            weights[idx_start:idx_end] <- 0
          }
        }

        # Stage 2: prune non-selected representation blocks
        if (!is.null(repr_prune_full) && length(repr_prune_full) > 0) {
          for (idx_name in repr_prune_full) {
            idx_i <- which(index_names == idx_name)
            if (length(idx_i) == 1) {
              idx_start <- (idx_i - 1) * n_bins + 1
              idx_end <- idx_i * n_bins
              weights[idx_start:idx_end] <- 0
            }
          }
        }

        # Stage 2: prune specific pentads
        if (!is.null(pentads_to_prune) && nrow(pentads_to_prune) > 0) {
          for (i in seq_len(nrow(pentads_to_prune))) {
            idx_name <- pentads_to_prune$index[i]
            pentad_j <- pentads_to_prune$pentad[i]
            idx_i <- which(index_names == idx_name)
            if (length(idx_i) == 1 && is.finite(pentad_j)) {
              feature_pos <- (idx_i - 1) * n_bins + pentad_j
              if (feature_pos >= 1 && feature_pos <= length(weights)) {
                weights[feature_pos] <- 0
              }
            }
          }
        }

        weights
      }

      # Re-train PCA-LDA with pruned features if pruning occurred
      # IMPORTANT: train_feature_pipeline expects BASE indices only (not L2norm prefixed)
      # In dual mode, it will add L2norm variants itself. So we need to extract base indices
      # from pruned_indices (which may contain both raw and L2norm names).
      pruned_base_indices <- unique(gsub("^L2norm_", "", pruned_indices))

      if (length(pruned_base_indices) < length(avail)) {
        cat(sprintf("\n[PERMUTATION] Re-training PCA-LDA with %d pruned base features...\n", length(pruned_base_indices)))
        MESMA_PARAMS <- train_feature_pipeline(multi_class_data, "target_class", pruned_base_indices)
        if (is.null(MESMA_PARAMS)) {
          cat("[PERMUTATION] ERROR: Re-training failed. Using original features.\n")
          MESMA_PARAMS <- MESMA_PARAMS_INITIAL
          MESMA_PARAMS$weights <- apply_pruning_rules(MESMA_PARAMS$weights, MESMA_PARAMS$indices, TEMPORAL_BUDGET,
                                                      indices_to_prune_stage1_full, repr_prune_full, pentads_to_prune)
        } else {
          cat(sprintf("[PERMUTATION] Re-training complete. New feature weights computed.\n"))
          MESMA_PARAMS$weights <- apply_pruning_rules(MESMA_PARAMS$weights, MESMA_PARAMS$indices, TEMPORAL_BUDGET,
                                                      indices_to_prune_stage1_full, repr_prune_full, pentads_to_prune)
          # Update avail to use pruned base indices
          avail <- pruned_base_indices
          cat(sprintf("[PERMUTATION] Updated avail to %d base features\n", length(avail)))
        }
      } else {
        cat("[PERMUTATION] No pruning performed. Using all features.\n")
        MESMA_PARAMS <- MESMA_PARAMS_INITIAL
        MESMA_PARAMS$weights <- apply_pruning_rules(MESMA_PARAMS$weights, MESMA_PARAMS$indices, TEMPORAL_BUDGET,
                                                    indices_to_prune_stage1_full, repr_prune_full, pentads_to_prune)
      }
    }

  } else {
    cat("[PERMUTATION] No OOB data available, using initial weights without permutation testing\n")
    MESMA_PARAMS <- MESMA_PARAMS_INITIAL
    pruned_indices <- avail
  }

  if (!is.null(MESMA_PARAMS) && !is.null(MESMA_PARAMS$weights)) {
    print_weights_summary("MESMA_FINAL", MESMA_PARAMS)
  }

  # Store pruning results
  assign("PRUNED_INDICES", pruned_indices, envir = globalenv())
  if (exists("perm_results")) {
    assign("PERMUTATION_RESULTS", perm_results, envir = globalenv())
  }

  # === STEP 4: Compute Confusion Matrix with Final Optimized Weights ===
  # This uses the permutation-pruned weights for accurate evaluation
  if (exists("df_train_oob") && !is.null(df_train_oob) && nrow(df_train_oob) > 0 &&
      exists("OPTIMAL_CLUSTER_COUNTS") && !is.null(OPTIMAL_CLUSTER_COUNTS) &&
      !is.null(MESMA_PARAMS) && !is.null(MESMA_PARAMS$weights)) {

    cat("\n=== STEP 4: Computing Confusion Matrix with Final Weights ===\n")

    tryCatch({
      # Prepare OOB data
      oob_for_cm <- dplyr::mutate(df_train_oob, target_class = tolower(as.character(Veg)))
      oob_for_cm <- dplyr::filter(oob_for_cm, !is.na(target_class) & target_class != "")

      unique_classes <- unique(oob_for_cm$target_class)
      veg_classes <- setdiff(unique_classes, "barren")

      # Build class endmembers from training data
      base_indices_cm <- if (!is.null(MESMA_PARAMS$base_indices)) MESMA_PARAMS$base_indices else avail
      indices_cm <- MESMA_PARAMS$indices
      n_base_idx_cm <- length(base_indices_cm)
      dual_mode_cm <- isTRUE(MESMA_PARAMS$dual_mode)
      l2_normalize_cm <- isTRUE(MESMA_PARAMS$l2_normalize)

      class_endmembers <- list()
      for (cls in unique_classes) {
        cls_data <- multi_class_data[multi_class_data$target_class == cls, ]
        if (nrow(cls_data) == 0) next

        cls_traces <- unique(cls_data[, c("location_id", "pheno_year")])
        cls_vecs <- list()

        for (j in seq_len(nrow(cls_traces))) {
          lid <- cls_traces$location_id[j]
          pyr <- cls_traces$pheno_year[j]
          sub <- cls_data[cls_data$location_id == lid & cls_data$pheno_year == pyr, ]
          # Build from base indices
          mat <- build_pentad_matrix(sub, base_indices_cm)
          if (!is.null(mat)) {
            vec_raw <- as.numeric(mat)

            # Build feature vector based on representation mode
            if (dual_mode_cm) {
              vec_l2 <- l2_normalize_perindex(vec_raw, n_base_idx_cm, TEMPORAL_BUDGET)
              vec <- c(vec_raw, vec_l2)
            } else if (l2_normalize_cm) {
              vec <- l2_normalize_perindex(vec_raw, n_base_idx_cm, TEMPORAL_BUDGET)
            } else {
              vec <- vec_raw
            }

            # Z-score
            for (k in seq_along(indices_cm)) {
              idx_start <- (k-1)*TEMPORAL_BUDGET + 1
              idx_end <- k*TEMPORAL_BUDGET
              idx_name <- indices_cm[k]
              param_idx <- which(names(MESMA_PARAMS$means) == idx_name)
              if (length(param_idx) > 0) {
                vec[idx_start:idx_end] <- (vec[idx_start:idx_end] - MESMA_PARAMS$means[param_idx]) / MESMA_PARAMS$sds[param_idx]
              }
            }
            vec[!is.finite(vec)] <- 0
            cls_vecs[[length(cls_vecs) + 1]] <- vec
          }
        }

        if (length(cls_vecs) > 0) {
          cls_mat <- do.call(rbind, cls_vecs)
          if (nrow(cls_mat) > 1) {
            dist_mat <- as.matrix(dist(cls_mat))
            medoid_idx <- which.min(rowSums(dist_mat))
            class_endmembers[[cls]] <- cls_mat[medoid_idx, ]
          } else {
            class_endmembers[[cls]] <- cls_mat[1, ]
          }
        }
      }

      if (length(class_endmembers) >= 2) {
        E <- do.call(cbind, class_endmembers)
        colnames(E) <- names(class_endmembers)

        # Build OOB samples (reuse base_indices_cm, indices_cm, dual_mode_cm, l2_normalize_cm from above)
        oob_traces <- unique(oob_for_cm[, c("location_id", "pheno_year", "target_class")])
        oob_vecs <- list()
        oob_labels <- c()

        for (j in seq_len(nrow(oob_traces))) {
          lid <- oob_traces$location_id[j]
          pyr <- oob_traces$pheno_year[j]
          true_cls <- oob_traces$target_class[j]

          sub <- oob_for_cm[oob_for_cm$location_id == lid & oob_for_cm$pheno_year == pyr, ]
          # Build from base indices
          mat <- build_pentad_matrix(sub, base_indices_cm)
          if (!is.null(mat)) {
            vec_raw <- as.numeric(mat)

            # Build feature vector based on representation mode
            if (dual_mode_cm) {
              vec_l2 <- l2_normalize_perindex(vec_raw, n_base_idx_cm, TEMPORAL_BUDGET)
              vec <- c(vec_raw, vec_l2)
            } else if (l2_normalize_cm) {
              vec <- l2_normalize_perindex(vec_raw, n_base_idx_cm, TEMPORAL_BUDGET)
            } else {
              vec <- vec_raw
            }

            # Z-score
            for (k in seq_along(indices_cm)) {
              idx_start <- (k-1)*TEMPORAL_BUDGET + 1
              idx_end <- k*TEMPORAL_BUDGET
              idx_name <- indices_cm[k]
              param_idx <- which(names(MESMA_PARAMS$means) == idx_name)
              if (length(param_idx) > 0) {
                vec[idx_start:idx_end] <- (vec[idx_start:idx_end] - MESMA_PARAMS$means[param_idx]) / MESMA_PARAMS$sds[param_idx]
              }
            }
            vec[!is.finite(vec)] <- 0
            oob_vecs[[length(oob_vecs) + 1]] <- vec
            oob_labels <- c(oob_labels, true_cls)
          }
        }

        if (length(oob_vecs) > 0) {
          Y <- do.call(rbind, oob_vecs)

          # Run FCLS with FINAL optimized weights
          all_coefs <- solve_batch_fcls(E, Y, MESMA_PARAMS$weights)

          if (!is.null(all_coefs)) {
            # Compute confusion matrix
            cm_labels <- names(class_endmembers)
            frac_mat <- matrix(0, nrow = length(cm_labels), ncol = length(cm_labels))
            rownames(frac_mat) <- cm_labels
            colnames(frac_mat) <- cm_labels
            class_counts <- setNames(rep(0, length(cm_labels)), cm_labels)

            correct_veg <- 0
            total_veg <- 0

            for (j in seq_len(nrow(Y))) {
              true_cls <- oob_labels[j]
              coefs <- all_coefs[j, ]
              names(coefs) <- colnames(E)

              for (cc in names(coefs)) {
                if (cc %in% cm_labels && true_cls %in% cm_labels) {
                  frac_mat[true_cls, cc] <- frac_mat[true_cls, cc] + coefs[cc]
                }
              }
              if (true_cls %in% cm_labels) {
                class_counts[true_cls] <- class_counts[true_cls] + 1
              }

              pred_cls <- names(which.max(coefs))
              if (true_cls %in% veg_classes) {
                total_veg <- total_veg + 1
                if (pred_cls == true_cls) correct_veg <- correct_veg + 1
              }
            }

            # Normalize and display
            avg_pred_frac <- frac_mat
            for (r in seq_len(nrow(frac_mat))) {
              if (class_counts[rownames(frac_mat)[r]] > 0) {
                avg_pred_frac[r, ] <- frac_mat[r, ] / class_counts[rownames(frac_mat)[r]]
              }
            }

            # Row-normalize
            row_sums <- rowSums(avg_pred_frac, na.rm = TRUE)
            for (ri in seq_len(nrow(avg_pred_frac))) {
              if (row_sums[ri] > 0) avg_pred_frac[ri, ] <- avg_pred_frac[ri, ] / row_sums[ri]
            }

            cat("\n[CONFUSION MATRIX] OOB Predicted Fractions - Row Normalized (All Classes):\n")
            print(round(avg_pred_frac, 3))

            # Vegetation-only matrix (row-normalized after barren subtraction)
            veg_rows <- intersect(veg_classes, rownames(avg_pred_frac))
            veg_cols <- intersect(veg_classes, colnames(avg_pred_frac))
            avg_cor_pred_veg_frac <- NA
            if (length(veg_rows) > 0 && length(veg_cols) > 0) {
              veg_matrix <- avg_pred_frac[veg_rows, veg_cols, drop = FALSE]
              row_sums_veg <- rowSums(veg_matrix, na.rm = TRUE)
              for (ri in seq_len(nrow(veg_matrix))) {
                if (row_sums_veg[ri] > 0) veg_matrix[ri, ] <- veg_matrix[ri, ] / row_sums_veg[ri]
              }
              cat("\n[CONFUSION MATRIX] OOB Predicted Fractions - Vegetation Only (Row Normalized):\n")
              print(round(veg_matrix, 3))

              # Compute avg correctly predicted veg fraction from diagonal of row-normalized veg matrix
              veg_diag_fracs <- sapply(veg_rows, function(cls) {
                if (cls %in% colnames(veg_matrix)) veg_matrix[cls, cls] else NA
              })
              avg_cor_pred_veg_frac <- mean(veg_diag_fracs, na.rm = TRUE)
            }

            # Also report classification accuracy for reference
            veg_class_accuracy <- if (total_veg > 0) correct_veg / total_veg else NA

            cat(sprintf("\n[CONFUSION MATRIX] Avg Row-Normalized Veg Fraction (diagonal): %.1f%%\n", 100 * avg_cor_pred_veg_frac))
            cat(sprintf("[CONFUSION MATRIX] Veg Classification Accuracy: %.1f%% (%d/%d correct)\n",
                        100 * veg_class_accuracy, correct_veg, total_veg))
          }
        }
      }
    }, error = function(e) {
      cat(sprintf("[CONFUSION MATRIX] Error computing confusion matrix: %s\n", e$message))
    })
  }

  # MESMA UNMIXING: All endmembers (barren + veg types) treated equally
  cat("[MODE] MESMA unmixing ENABLED (barren and vegetation types treated as equals)\n")

  # Log loss function choice
  if (exists("USE_HUBER_LOSS") && isTRUE(USE_HUBER_LOSS)) {
    cat(sprintf("[MODE] Using HUBER LOSS for FCLS (delta=%.3f, robust to outliers)\n",
                if (exists("HUBER_DELTA")) HUBER_DELTA else 1.345))
  } else {
    cat("[MODE] Using standard RMSE loss for FCLS\n")
  }

  # === SPLINE MODE BRANCH ===
  # If USE_SPLINE_ENDMEMBERS is TRUE, build spline library instead of discrete endmembers
  if (isTRUE(USE_SPLINE_ENDMEMBERS)) {
    cat("\n=== SPLINE ENDMEMBER MODE ENABLED ===\n")
    cat("Preparing training data for spline library (EXCLUDING validation locations)...\n")

    # Prefer explicit df_train (should already exclude validation locations). If not present, fall back to filtering df
    train_for_spline <- NULL
    if (exists("df_train") && !is.null(df_train) && nrow(df_train) > 0) {
      train_for_spline <- df_train
    } else if (exists("df") && !is.null(df) && nrow(df) > 0) {
      train_for_spline <- df
    } else {
      stop("No training data available for spline library building")
    }

    # If df_validation exists, ensure its locations are excluded
    if (exists("df_validation") && !is.null(df_validation) && nrow(df_validation) > 0) {
      val_locs <- unique(df_validation$location_id)
      before_n <- nrow(train_for_spline)
      train_for_spline <- train_for_spline[!(train_for_spline$location_id %in% val_locs), , drop = FALSE]
      after_n <- nrow(train_for_spline)
      cat(sprintf("  Excluded %d rows belonging to %d validation locations (training rows: %d -> %d)\n", before_n - after_n, length(val_locs), before_n, after_n))
    } else {
      cat(sprintf("  Training rows used for spline building: %d\n", nrow(train_for_spline)))
    }

    # Compute SPLINE-SPECIFIC normalization parameters from raw observations
    # IMPORTANT: MESMA_PARAMS$means/sds are computed from pentad matrices (37 bins × K),
    # but spline mode works with individual observations. We need per-index stats from raw data.
    cat("[SPLINE] Computing observation-level normalization parameters...\n")
    spline_means <- numeric(length(avail))
    spline_sds <- numeric(length(avail))
    names(spline_means) <- avail
    names(spline_sds) <- avail
    
    for (idx in avail) {
      if (idx %in% names(train_for_spline)) {
        vals <- train_for_spline[[idx]]
        vals <- vals[is.finite(vals)]
        spline_means[idx] <- mean(vals, na.rm = TRUE)
        spline_sds[idx] <- sd(vals, na.rm = TRUE)
        if (!is.finite(spline_sds[idx]) || spline_sds[idx] < 1e-8) spline_sds[idx] <- 1
      } else {
        spline_means[idx] <- 0
        spline_sds[idx] <- 1
      }
    }
    cat(sprintf("[SPLINE] Computed normalization for %d indices from %d observations\n", 
                length(avail), nrow(train_for_spline)))
    cat(sprintf("[SPLINE] Sample means: %s\n", 
                paste(sprintf("%s=%.3f", head(avail, 4), spline_means[head(avail, 4)]), collapse=", ")))
    cat(sprintf("[SPLINE] Sample sds: %s\n", 
                paste(sprintf("%s=%.3f", head(avail, 4), spline_sds[head(avail, 4)]), collapse=", ")))

    # Store spline-specific normalization parameters
    SPLINE_PARAMS <- list(
      indices = avail,
      means = spline_means,
      sds = spline_sds,
      weights = NULL  # Not used in spline mode
    )

    # Build spline library using the filtered training data WITH normalization
    SPLINE_LIBRARY <- build_spline_library(train_for_spline, avail, ALLOWED_VEG, spar = SPLINE_SPAR, norm_params = SPLINE_PARAMS)

    assign("SPLINE_LIBRARY", SPLINE_LIBRARY, envir = globalenv())
    assign("SPLINE_PARAMS", SPLINE_PARAMS, envir = globalenv())

    cat(sprintf("[SPLINE] Library built with %d classes: %s\n", 
                length(SPLINE_LIBRARY), paste(names(SPLINE_LIBRARY), collapse = ", ")))
    cat(sprintf("[SPLINE] Using %d indices: %s\n", length(avail), paste(avail, collapse = ", ")))
    
    # DIAGNOSTIC: Check class separability at key DOYs (summer peak DOY ~180)
    cat("\n[SPLINE DIAGNOSTIC] Class separability check (DOY=180, peak growing season):\n")
    test_doy <- 180
    classes <- names(SPLINE_LIBRARY)
    for (idx in head(avail, 5)) {
      cat(sprintf("  %s: ", idx))
      vals <- sapply(classes, function(cls) {
        sp <- SPLINE_LIBRARY[[cls]][[idx]]
        if (!is.null(sp) && !is.null(sp$predict)) sp$predict(test_doy) else NA
      })
      cat(paste(sprintf("%s=%.3f", names(vals), vals), collapse=", "))
      cat(sprintf(" | range=%.3f\n", diff(range(vals, na.rm=TRUE))))
    }
    
    cat("[SPLINE] Skipping endmember extraction (not needed for spline mode)\n")

    # Create minimal placeholder structures for compatibility
    mesma_lib <- list()
    OPTIMIZED_LIBRARY <- list()
    mesma_lib <- list()
    OPTIMIZED_LIBRARY <- list()
    compressed_templates_accessor <- list()

    assign("MESMA_PARAMS", MESMA_PARAMS, envir = globalenv())
    assign("mesma_lib", mesma_lib, envir = globalenv())
    assign("OPTIMIZED_LIBRARY", OPTIMIZED_LIBRARY, envir = globalenv())
    assign("compressed_templates_accessor", compressed_templates_accessor, envir = globalenv())
    assign(".COMPRESSED_TEMPLATES_ACCESSOR", compressed_templates_accessor, envir = globalenv())
    
  } else {
  # === STANDARD ENDMEMBER MODE ===

  cat("Building final MESMA library (barren + all vegetation types)...\n")

  # Use precomputed cluster counts from Step 2 (if available) to skip re-optimization
  precomputed_clusters <- NULL
  if (exists("OPTIMAL_CLUSTER_COUNTS", envir = globalenv())) {
    precomputed_clusters <- get("OPTIMAL_CLUSTER_COUNTS", envir = globalenv())
  }

  # Use indices from MESMA_PARAMS (includes L2norm if enabled)
  indices_for_final_library <- if (!is.null(MESMA_PARAMS$indices)) MESMA_PARAMS$indices else avail

  mesma_lib <- build_mesma_library_weighted(df_train, indices_for_final_library, MESMA_PARAMS, ALLOWED_VEG, precomputed_clusters = precomputed_clusters)

  cat("Pre-computing optimized library for MESMA...\n")
  OPTIMIZED_LIBRARY <- precompute_optimized_library_weighted(mesma_lib, grid_type = "full", feature_weights = (if (!is.null(MESMA_PARAMS) && !is.null(MESMA_PARAMS$weights)) MESMA_PARAMS$weights else NULL) )

  # Summarize pruning if any
  if (exists("PRUNE_ZERO_WEIGHT_FEATURES") && isTRUE(PRUNE_ZERO_WEIGHT_FEATURES)) {
    total_pruned <- 0
    total_kept <- 0
    vegs_pruned <- 0
    for (v in names(OPTIMIZED_LIBRARY)) {
      libv <- OPTIMIZED_LIBRARY[[v]]
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
    nn <- sapply(names(OPTIMIZED_LIBRARY), function(v) {
      libv <- OPTIMIZED_LIBRARY[[v]]
      if (is.null(libv) || is.null(libv$M)) return(0)
      nrow(libv$M)
    })
    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG] MESMA library built: endmember_count=%d variants_per_type=%s\n",
                length(nn), paste(sprintf("%s:%d", names(nn), nn), collapse=", ")))
  }

  assign("MESMA_PARAMS", MESMA_PARAMS, envir = globalenv())
  assign("mesma_lib", mesma_lib, envir = globalenv())
  assign("OPTIMIZED_LIBRARY", OPTIMIZED_LIBRARY, envir = globalenv())

  cat(sprintf("[NOTICE] MESMA feature count: avail=%d, params_indices=%d\n",
              length(avail), length(MESMA_PARAMS$indices)))
  cat(sprintf("[NOTICE] MESMA feature list: %s\n", paste(MESMA_PARAMS$indices, collapse=", ")))

  mesma_lib <- mesma_lib

  assign("OPTIMIZED_LIBRARY", OPTIMIZED_LIBRARY, envir = globalenv())
  assign("mesma_lib", mesma_lib, envir = globalenv())

  compressed_templates_accessor <- precompute_compressed_templates(mesma_lib, "full")
  assign("compressed_templates_accessor", compressed_templates_accessor, envir = globalenv())
  assign(".COMPRESSED_TEMPLATES_ACCESSOR", compressed_templates_accessor, envir = globalenv())

  
  } # END else (standard endmember mode)
  # === END MODE BRANCH ===

# ========================================================================== 
# TRAINING DISABLED
# This repository/script no longer performs model *training*. Training is expected
# to be performed offline and any required artifacts (normalization params,
# libraries, models) should be provided on disk. The script will continue to
# perform inference and validation where possible.
# ==========================================================================

# Load validation locations for confusion matrix (held-out 20% from stratified split)
val_locations_file <- file.path(OUT_DIR, "validation_locations.csv")
validation_location_ids <- character(0)

if (file.exists(val_locations_file)) {
  val_locs_df <- tryCatch(read.csv(val_locations_file, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(val_locs_df) && "location_id" %in% names(val_locs_df)) {
    validation_location_ids <- unique(as.character(val_locs_df$location_id))
  }
} 

# Store validation location IDs globally for later confusion matrix computation
if (length(validation_location_ids) > 0) {
  assign("validation_location_ids", validation_location_ids, envir = globalenv())
} else {
  # no validation locations available
}

# NOTE: do NOT create or populate `df_tasks` here — inference code will build
# task tables from the inference input. Training/task-generation code was
# removed; keep this script focused on inference + validation only.

if (isTRUE(TESTING_MODE)) cat("[DEBUG] Reached line 6377 - about to define fit_one_task function\n")

# ============================================================================
# OOB FRACTION PREDICTION UNCERTAINTY
# Uses OOB validation residuals to add realistic prediction error to MC draws.
# For each true class, we store the distribution of prediction residuals
# (predicted_fraction - true_fraction) observed in OOB validation.
# During MC, we sample from these residuals to perturb the predicted fractions.
# ============================================================================

# Global storage for OOB fraction residuals per true class
.OOB_FRACTION_RESIDUALS <- NULL  # List: true_class -> matrix of residuals (n_samples x n_classes)

# Store OOB fraction residuals for use in MC uncertainty propagation
# Called after OOB validation with known true classes
#
# Parameters:
#   residuals_by_class: Named list where each element is a matrix of residuals
#                       Rows = OOB samples, Cols = predicted classes
#                       Each row is (predicted_fractions - true_fractions) for one sample
#   true_class_names: Names of the true classes (should match list names)
store_oob_fraction_residuals <- function(residuals_by_class) {
  if (!is.list(residuals_by_class) || length(residuals_by_class) == 0) {
    warning("[OOB_FRAC] Invalid residuals - must be a non-empty named list")
    return(invisible(FALSE))
  }

  # Validate structure
  for (cls in names(residuals_by_class)) {
    resid_mat <- residuals_by_class[[cls]]
    if (!is.matrix(resid_mat) && !is.data.frame(resid_mat)) {
      warning(sprintf("[OOB_FRAC] Residuals for class '%s' must be a matrix", cls))
      return(invisible(FALSE))
    }
  }

  assign(".OOB_FRACTION_RESIDUALS", residuals_by_class, envir = globalenv())

  if (isTRUE(DEBUG_UNCERTAINTY)) {
    cat("[OOB_FRAC] Stored OOB fraction residuals for MC uncertainty:\n")
    for (cls in names(residuals_by_class)) {
      n_samples <- nrow(residuals_by_class[[cls]])
      cat(sprintf("  %s: %d samples\n", cls, n_samples))
    }
  }

  return(invisible(TRUE))
}

# Sample a residual vector for a given dominant predicted class
# Returns a named vector of residuals to add to predicted fractions
sample_oob_residual <- function(dominant_class) {
  if (!exists(".OOB_FRACTION_RESIDUALS", envir = globalenv())) {
    return(NULL)
  }

  residuals_by_class <- get(".OOB_FRACTION_RESIDUALS", envir = globalenv())
  if (is.null(residuals_by_class)) return(NULL)

  # Match dominant class (case-insensitive)
  cls_lower <- tolower(dominant_class)
  available_classes <- tolower(names(residuals_by_class))
  match_idx <- which(available_classes == cls_lower)

  if (length(match_idx) == 0) {
    return(NULL)
  }

  resid_mat <- residuals_by_class[[match_idx]]
  if (is.null(resid_mat) || nrow(resid_mat) == 0) {
    return(NULL)
  }

  # Sample one row uniformly at random
  row_idx <- sample.int(nrow(resid_mat), 1)
  residual <- as.numeric(resid_mat[row_idx, ])
  names(residual) <- colnames(resid_mat)

  return(residual)
}

# ============================================================================
# CLASSIFICATION UNCERTAINTY: DIRICHLET PERTURBATION
# Propagates classification/misclassification uncertainty through bootstrap
# by resampling fractional covers from a Dirichlet distribution centered on
# the confusion matrix row for the true class.
# ============================================================================

# Global storage for confusion matrix and validation sample sizes
.CONFUSION_MATRIX <- NULL
.VALIDATION_SAMPLE_SIZES <- NULL

# Store the confusion matrix for use in Dirichlet perturbation
# Called after validation accuracy is computed
store_confusion_matrix <- function(conf_matrix, sample_sizes = NULL) {
  if (!is.matrix(conf_matrix) || nrow(conf_matrix) != ncol(conf_matrix)) {
    warning("[DIRICHLET] Invalid confusion matrix - must be square matrix")
    return(invisible(FALSE))
  }

  # Ensure row-normalized (each row sums to 1)
  row_sums <- rowSums(conf_matrix, na.rm = TRUE)
  conf_matrix_norm <- conf_matrix
  for (i in seq_len(nrow(conf_matrix_norm))) {
    if (row_sums[i] > 1e-9) {
      conf_matrix_norm[i, ] <- conf_matrix_norm[i, ] / row_sums[i]
    }
  }

  assign(".CONFUSION_MATRIX", conf_matrix_norm, envir = globalenv())

  if (!is.null(sample_sizes)) {
    assign(".VALIDATION_SAMPLE_SIZES", sample_sizes, envir = globalenv())
  }

  if (isTRUE(DEBUG_UNCERTAINTY)) {
    cat(sprintf("[DIRICHLET] Stored %dx%d confusion matrix for classification uncertainty\n",
                nrow(conf_matrix_norm), ncol(conf_matrix_norm)))
    cat(sprintf("[DIRICHLET] Classes: %s\n", paste(rownames(conf_matrix_norm), collapse = ", ")))
    if (!is.null(sample_sizes)) {
      cat(sprintf("[DIRICHLET] Validation sample sizes: %s\n",
                  paste(sprintf("%s=%d", names(sample_sizes), sample_sizes), collapse = ", ")))
    }
  }

  return(invisible(TRUE))
}

# Apply Dirichlet perturbation to fractional covers based on true class
# This simulates classification uncertainty by drawing from a Dirichlet
# distribution centered on the confusion matrix row for the true class
#
# Parameters:
#   fractions: Named numeric vector of fractional covers (must sum to ~1)
#   true_class: The true vegetation class for this location
#   conf_matrix: Row-normalized confusion matrix (optional, uses global if NULL)
#   sample_sizes: Named vector of validation sample sizes per class (optional)
#   concentration_scale: Multiplier for sample size to get Dirichlet alpha
#
# Returns:
#   Perturbed fractional covers (same names, still sum to 1)
#
apply_dirichlet_perturbation <- function(fractions, true_class,
                                          conf_matrix = NULL,
                                          sample_sizes = NULL,
                                          concentration_scale = DIRICHLET_CONCENTRATION_SCALE) {
  # Use global confusion matrix if not provided
 if (is.null(conf_matrix)) {
    if (exists(".CONFUSION_MATRIX", envir = globalenv())) {
      conf_matrix <- get(".CONFUSION_MATRIX", envir = globalenv())
    } else {
      # No confusion matrix available - return original fractions
      return(fractions)
    }
  }

  if (is.null(sample_sizes) && exists(".VALIDATION_SAMPLE_SIZES", envir = globalenv())) {
    sample_sizes <- get(".VALIDATION_SAMPLE_SIZES", envir = globalenv())
  }

  # Validate inputs
  if (is.null(conf_matrix) || !is.matrix(conf_matrix)) {
    return(fractions)
  }

  # Match true_class to confusion matrix row (case-insensitive)
  row_names <- tolower(rownames(conf_matrix))
  true_class_lower <- tolower(true_class)
  row_idx <- which(row_names == true_class_lower)

  if (length(row_idx) == 0) {
    # True class not in confusion matrix - return original
    return(fractions)
  }

  # Get the confusion row for the true class (this is our Dirichlet mean)
  mu <- as.numeric(conf_matrix[row_idx, ])
  class_names <- colnames(conf_matrix)

  # Ensure mu sums to 1 and has no zeros (add small epsilon for stability)
  mu[mu < 1e-6] <- 1e-6
  mu <- mu / sum(mu)

  # Determine concentration parameter alpha
  # Higher alpha = tighter concentration around mean = less noise
  if (!is.null(sample_sizes) && true_class_lower %in% tolower(names(sample_sizes))) {
    n_val <- sample_sizes[tolower(names(sample_sizes)) == true_class_lower]
    alpha <- concentration_scale * n_val
  } else {
    # Default: assume moderate validation sample size
    alpha <- concentration_scale * 20
  }

  # Ensure minimum alpha for numerical stability
  alpha <- max(alpha, 5)

  # Dirichlet parameter vector
  alpha_vec <- alpha * mu

  # Draw from Dirichlet distribution
  # Using the gamma distribution method: X_i ~ Gamma(alpha_i, 1), then normalize
  gamma_draws <- rgamma(length(alpha_vec), shape = alpha_vec, rate = 1)

  # Handle edge case where all draws are zero
  if (sum(gamma_draws) < 1e-10) {
    gamma_draws <- mu  # Fall back to mean
  }

  perturbed <- gamma_draws / sum(gamma_draws)

  # Map back to original fraction names
  # The confusion matrix classes may not exactly match the fraction names
  result <- fractions
  frac_names_lower <- tolower(names(fractions))

  for (i in seq_along(class_names)) {
    class_lower <- tolower(class_names[i])
    # Find matching fraction column (handle "frac_" prefix)
    match_idx <- which(frac_names_lower == class_lower |
                       frac_names_lower == paste0("frac_", class_lower))
    if (length(match_idx) > 0) {
      result[match_idx[1]] <- perturbed[i]
    }
  }

  # Re-normalize to ensure sum = 1 (in case of partial matching)
  if (sum(result) > 0) {
    result <- result / sum(result)
  }

  return(result)
}

# Batch apply Dirichlet perturbation to a data frame of coefficients
# Useful for bootstrap resampling of location-level results
#
# Parameters:
#   coef_df: Data frame with columns for fractional covers and a true_veg column
#   frac_cols: Names of fractional cover columns
#   true_class_col: Name of column containing true vegetation class
#
# Returns:
#   Data frame with perturbed fractional covers
#
apply_dirichlet_perturbation_batch <- function(coef_df, frac_cols, true_class_col = "true_veg") {
  if (!isTRUE(ENABLE_CLASSIFICATION_UNCERTAINTY)) {
    return(coef_df)
  }

  if (!true_class_col %in% names(coef_df)) {
    warning(sprintf("[DIRICHLET] True class column '%s' not found in data frame", true_class_col))
    return(coef_df)
  }

  # Check if confusion matrix is available
  if (!exists(".CONFUSION_MATRIX", envir = globalenv())) {
    if (isTRUE(DEBUG_UNCERTAINTY)) {
      cat("[DIRICHLET] No confusion matrix available - skipping perturbation\n")
    }
    return(coef_df)
  }

  result_df <- coef_df

  for (i in seq_len(nrow(coef_df))) {
    true_class <- as.character(coef_df[[true_class_col]][i])

    if (is.na(true_class) || true_class == "") {
      next
    }

    # Extract current fractions
    fracs <- as.numeric(coef_df[i, frac_cols])
    names(fracs) <- frac_cols

    # Apply perturbation
    perturbed <- apply_dirichlet_perturbation(fracs, true_class)

    # Store back
    result_df[i, frac_cols] <- perturbed
  }

  return(result_df)
}

# Simple i.i.d. residual bootstrap
# Resamples residuals with replacement (assumes approximate independence)
simple_residual_bootstrap <- function(residuals) {
  n <- length(residuals)
  if (n == 0) return(numeric(0))
  
  # Remove NA values for bootstrapping
  residuals_non_na <- residuals[!is.na(residuals)]
  n_non_na <- length(residuals_non_na)
  
  if (n_non_na == 0) return(residuals) # Return original if all are NA
  
  # Simple random sampling with replacement
  resampled_indices <- sample(1:n_non_na, n_non_na, replace = TRUE)
  
  # Create the bootstrapped residuals vector
  bootstrapped_residuals_non_na <- residuals_non_na[resampled_indices]
  
  # Place bootstrapped residuals back into the original structure with NAs
  residuals_boot <- residuals
  residuals_boot[!is.na(residuals)] <- bootstrapped_residuals_non_na
  
  return(residuals_boot)
}

  fit_one_task <- function(task_data) {
    # File-based debug log that won't be suppressed by sink()
    debug_log_file <- file.path(tempdir(), "fit_one_location_debug.log")
    write_debug <- function(msg) {
      tryCatch(write(paste0(Sys.time(), " [fit_one_task] ", msg), debug_log_file, append = TRUE), error = function(e) {})
    }
    
    if (is.null(task_data) || nrow(task_data) == 0) return(NULL)

    loc <- as.character(task_data$location_id[1])
    yr <- as.integer(task_data$pheno_year[1])
    
    write_debug(sprintf("Processing loc=%s yr=%d rows=%d", loc, yr, nrow(task_data)))

    # Extract lat/lon if available
    lat_val <- NA_real_
    lon_val <- NA_real_
    if ("lat" %in% names(task_data)) lat_val <- as.numeric(task_data$lat[1]) else if ("latitude" %in% names(task_data)) lat_val <- as.numeric(task_data$latitude[1])
    if ("lon" %in% names(task_data)) lon_val <- as.numeric(task_data$lon[1]) else if ("longitude" %in% names(task_data)) lon_val <- as.numeric(task_data$longitude[1])

    # Use MESMA parameters
    PARAMS <- MESMA_PARAMS
    
    # DEBUG: Verify PARAMS structure
    if (is.null(PARAMS) || is.null(PARAMS$indices) || is.null(PARAMS$means) || is.null(PARAMS$sds)) {
      cat(sprintf("[ERROR] loc=%s yr=%d: PARAMS structure is incomplete (PARAMS=%s, indices=%s, means=%s, sds=%s)\n",
                  loc, yr,
                  !is.null(PARAMS), 
                  !is.null(PARAMS$indices), 
                  !is.null(PARAMS$means), 
                  !is.null(PARAMS$sds)))
      write_debug(sprintf("loc=%s yr=%d: PARAMS incomplete - returning NULL", loc, yr))
      return(NULL)
    }
    
    # DEBUG: Log feature space
    write_debug(sprintf("loc=%s yr=%d: PARAMS$indices (%d features): %s", 
                       loc, yr, length(PARAMS$indices), paste(PARAMS$indices, collapse=", ")))
    if (isTRUE(TESTING_MODE)) {
      cat(sprintf("[FEATURE SPACE] loc=%s yr=%d: Using %d indices: %s\n",
                  loc, yr, length(PARAMS$indices), paste(head(PARAMS$indices, 5), collapse=", ")))
      cat(sprintf("[FEATURE SPACE] task_data columns: %s\n", paste(head(names(task_data), 10), collapse=", ")))
    }
    
    # CRITICAL: Validate that all required indices are present in task_data
    missing_indices <- setdiff(PARAMS$indices, names(task_data))
    if (length(missing_indices) > 0) {
      cat(sprintf("[ERROR] loc=%s yr=%d: Missing indices in inference data: %s\n",
                  loc, yr, paste(missing_indices, collapse=", ")))
      cat(sprintf("[ERROR] Available indices in task_data: %s\n",
                  paste(intersect(PARAMS$indices, names(task_data)), collapse=", ")))
      write_debug(sprintf("loc=%s yr=%d: Missing indices %s - returning NULL", 
                         loc, yr, paste(missing_indices, collapse=", ")))
      return(NULL)
    }
    
    # For inference/validation data, do NOT interpolate - use only actual observations
    # Use base_indices for building pentad matrix (raw indices that exist in the data)
    base_indices_for_build <- if (!is.null(PARAMS$base_indices)) PARAMS$base_indices else PARAMS$indices
    raw_mat <- build_pentad_matrix(task_data, base_indices_for_build, interpolate = FALSE)
    if (is.null(raw_mat)) {
      if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG] loc=%s yr=%d: build_pentad_matrix returned NULL\n", loc, yr))
      write_debug(sprintf("loc=%s yr=%d: build_pentad_matrix returned NULL", loc, yr))
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
      write_debug(sprintf("loc=%s yr=%d: n_valid=0 - returning NULL", loc, yr))
      return(NULL)
    } else {
      MIN_VALID_FRACTION <- 0.01
      if (exists("TESTING_MODE") && isTRUE(TESTING_MODE) && n_valid < (length(y_raw) * MIN_VALID_FRACTION)) {
        cat(sprintf("[WARN] loc=%s pheno_year=%d: low valid observations (%d < %.0f) - proceeding anyway\n",
                    loc, yr, n_valid, length(y_raw) * MIN_VALID_FRACTION))
      }
    }
    # ============================

    # --- New Logic: Estimate barren fraction using MSAVI linear model ---
    # Moved from end of function to allow early exit for pure barren pixels
    # 1. Get current MSAVI for the task

    # Require exact 'MSAVI' column only (no fallbacks). Compute summer median from raw MSAVI values.
    if (!"MSAVI" %in% names(task_data)) {
      available_cols <- paste(head(names(task_data), 20), collapse=", ")
      stop(sprintf("ERROR loc=%s yr=%d: 'MSAVI' column not found in task_data. This workflow requires raw 'MSAVI' values (no fallbacks). First 20 available columns: %s",
                   loc, yr, available_cols))
    }
    msavi_col <- "MSAVI"

    # Error if no MSAVI column found
    if (is.na(msavi_col)) {
      available_cols <- paste(head(names(task_data), 20), collapse=", ")
      stop(sprintf("ERROR loc=%s yr=%d: No MSAVI column found in task_data. Tried: MSAVI_median, MSAVI2_median, MSAVI, MSAVI2. First 20 available columns: %s",
                   loc, yr, available_cols))
    }

    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG MSAVI CHECK] loc=%s yr=%d: Using MSAVI column '%s', nrow=%d\n",
                loc, yr, msavi_col, nrow(task_data)))

    # Filter task_data for summer months (June-September)
    summer_task_data <- if ("month" %in% names(task_data)) {
      task_data[task_data$month %in% 6:9, ]
    } else {
      task_data
    }

    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG MSAVI] loc=%s yr=%d: summer_task_data has %d rows (filtered from %d)\n",
                loc, yr, nrow(summer_task_data), nrow(task_data)))

    # Error if no summer data available
    if (nrow(summer_task_data) == 0) {
      stop(sprintf("ERROR loc=%s yr=%d: No summer data (June-September) available for MSAVI extraction. Total rows: %d",
                   loc, yr, nrow(task_data)))
    }

    current_msavi <- median(summer_task_data[[msavi_col]], na.rm = TRUE)

    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG MSAVI] loc=%s yr=%d: Summer MSAVI median = %.4f (from %d values, %d valid)\n",
                loc, yr, current_msavi, length(summer_task_data[[msavi_col]]), sum(is.finite(summer_task_data[[msavi_col]]))))

    # Error if MSAVI is not finite
    if (!is.finite(current_msavi)) {
      stop(sprintf("ERROR loc=%s yr=%d: MSAVI is not finite (value=%s). Cannot calculate FVC without valid MSAVI data.",
                   loc, yr, as.character(current_msavi)))
    }

    # This avoids dependency on external calibration models.
    fvc_predicted <- current_msavi / (0.5 + current_msavi)
    # Log approximate MM-based mapping (note: not a calibrated FVC model)
    write_debug(sprintf("loc=%s yr=%d: MSAVI=%.4f -> FVC=%.4f (MM fallback)", loc, yr, current_msavi, fvc_predicted))

    # Clamp to [0, 1] for safety
    total_veg_cover <- pmin(pmax(fvc_predicted, 0), 1)
    barren_fraction <- 1 - total_veg_cover

    write_debug(sprintf("loc=%s yr=%d: MSAVI=%.4f -> FVC=%.4f (MM fallback)", loc, yr, current_msavi, total_veg_cover))

    # Verify barren_fraction is finite - error if not
    if (!is.finite(barren_fraction)) {
      stop(sprintf("ERROR loc=%s yr=%d: barren_fraction is not finite (value=%s) after MSAVI FVC prediction. This should never happen.",
                   loc, yr, as.character(barren_fraction)))
    }

    write_debug(sprintf("loc=%s yr=%d: Barren fraction from MSAVI linear model=%.4f, proceeding to unmixing", loc, yr, barren_fraction))

    # ===== MESMA UNMIXING =====
      # Get indices from params
      n_bins <- TEMPORAL_BUDGET
      indices <- PARAMS$indices
      base_indices <- if (!is.null(PARAMS$base_indices)) PARAMS$base_indices else indices
      n_base_idx <- length(base_indices)
      dual_mode <- isTRUE(PARAMS$dual_mode)
      l2_normalize <- isTRUE(PARAMS$l2_normalize)

      # Build feature vector based on representation mode
      # Also update valid_mask to match the new feature vector length
      if (dual_mode) {
        # DUAL MODE: concatenate raw and L2-normalized
        y_l2 <- l2_normalize_perindex(y_raw, n_base_idx, n_bins)
        y_work <- c(y_raw, y_l2)
        # Extend valid_mask: a position is valid if the corresponding raw position is valid
        valid_mask <- c(valid_mask, valid_mask)
      } else if (l2_normalize) {
        # L2 ONLY MODE
        y_work <- l2_normalize_perindex(y_raw, n_base_idx, n_bins)
        # valid_mask stays the same (same length as y_raw)
      } else {
        # RAW ONLY MODE
        y_work <- y_raw
        # valid_mask stays the same
      }

      # Z-score values
      n_idx <- length(indices)
      y_norm <- y_work
      for(k in seq_len(n_idx)) {
        idx_start <- (k-1)*n_bins + 1
        idx_end <- k*n_bins

        idx_name <- indices[k]
        param_idx <- which(names(PARAMS$means) == idx_name)
        if (length(param_idx) > 0) {
          mu <- PARAMS$means[param_idx]
          sigma <- PARAMS$sds[param_idx]
          if (!is.finite(sigma) || sigma < EPS_SIGMA) {
            if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_task] MESMA sigma for index %s is non-finite or too small (%.8f); using EPS_SIGMA=%.8e\n", idx_name, sigma, EPS_SIGMA))
            sigma <- EPS_SIGMA
          }
          y_norm[idx_start:idx_end] <- (y_work[idx_start:idx_end] - mu) / sigma
        }
      }
      y_norm[!is.finite(y_norm)] <- 0

      # Apply PCA-LDA transform to the z-scored observation so inference sits in the same
      # PCA-LDA-scored feature space as training. This is applied before masking so column
      # indices remain aligned between obs and library templates.
      if (!is.null(PARAMS$weights) && length(PARAMS$weights) == length(y_norm)) {
        y_norm_full_pca_lda <- tryCatch({ apply_pca_lda_transform(y_norm, PARAMS) }, error = function(e) { if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG] apply_pca_lda_transform failed for loc=%s yr=%d: %s\n", loc, yr, e$message)); y_norm })
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG PCA-LDA] Applied PCA-LDA transform to z-scored observation (loc=%s yr=%d).\n", loc, yr))
      } else {
        y_norm_full_pca_lda <- y_norm
      }

      # Mask the observation (z-scored only, no PCA-LDA weighting applied to data)
      y_norm_masked <- y_norm_full_pca_lda[valid_mask]

      # Extract masked PCA-LDA weights to pass to solver
      # The solver will apply weights directly to both endmembers and observations
      if (!is.null(PARAMS$weights) && length(PARAMS$weights) == length(y_norm_full_pca_lda)) {
        weights_masked <- PARAMS$weights[valid_mask]
      } else {
        weights_masked <- rep(1, length(y_norm_masked))
      }

      # Check if we have sufficient signal (using unweighted norm for validation)
      y_norm_val <- sqrt(sum(y_norm_masked^2, na.rm = TRUE))

      if (is.na(y_norm_val) || y_norm_val < 1e-9) {
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_task] loc=%s yr=%d: insufficient signal (norm=%.6g), skipping\n", loc, yr, y_norm_val))
        return(NULL)
      }

      # Use the masked, z-score normalized observation (NOT pre-weighted)
      # Weighting will be applied inside solve_weights_fcls via feature_weights parameter
      y_for_unmixing <- y_norm_masked

      # Perform MESMA using all endmembers (barren + all veg types)
      # Barren is kept as from MESMA - no replacement
      veg_kept <- names(mesma_lib)
      # SAFEGUARD: Ensure a 'barren' endmember exists
      if (!"barren" %in% veg_kept) {
        warning("No 'barren' endmember found in MESMA library - ensure barren endmembers are present.")
      }
      
      # DEBUG: Check library availability
      if (is.null(OPTIMIZED_LIBRARY) || length(OPTIMIZED_LIBRARY) == 0) {
        cat(sprintf("[ERROR] loc=%s yr=%d: OPTIMIZED_LIBRARY is NULL or empty\n", loc, yr))
        write_debug(sprintf("loc=%s yr=%d: OPTIMIZED_LIBRARY NULL or empty - returning NULL", loc, yr))
        return(NULL)
      }
      write_debug(sprintf("loc=%s yr=%d: OPTIMIZED_LIBRARY has %d veg types: %s", 
                         loc, yr, length(OPTIMIZED_LIBRARY), paste(names(OPTIMIZED_LIBRARY), collapse=", ")))
      if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG LIBRARY] loc=%s yr=%d: OPTIMIZED_LIBRARY contains %d vegetation types: %s\n", 
                  loc, yr, length(OPTIMIZED_LIBRARY), paste(names(OPTIMIZED_LIBRARY), collapse=", ")))
      
      # Pre-compute barren references for filtering (needed before variant loop)
      # Use ALL barren endmembers as references to support multiple barren prototypes
      barren_refs_eval_norm <- NULL
      if ("barren" %in% names(OPTIMIZED_LIBRARY)) {
        barren_lib <- OPTIMIZED_LIBRARY[["barren"]]
        if (!is.null(barren_lib) && !is.null(barren_lib$M_norm) && nrow(barren_lib$M_norm) > 0) {
          # Use all barren templates for filtering, apply valid_mask
          barren_refs_masked <- barren_lib$M_norm[, valid_mask, drop = FALSE]
          barren_refs_masked[!is.finite(barren_refs_masked)] <- 0
          # Normalize each barren endmember
          barren_refs_eval_norm <- t(apply(barren_refs_masked, 1, function(r) {
            nrm <- sqrt(sum(r^2))
            if (!is.finite(nrm) || nrm == 0) return(rep(0, length(r))) else return(r / nrm)
          }))
          # Remove zero rows
          keep_rows <- rowSums(abs(barren_refs_eval_norm)) > 0
          if (any(keep_rows)) {
            barren_refs_eval_norm <- barren_refs_eval_norm[keep_rows, , drop = FALSE]
            write_debug(sprintf("loc=%s yr=%d: %d barren references prepared for filtering", loc, yr, nrow(barren_refs_eval_norm)))
          } else {
            barren_refs_eval_norm <- NULL
          }
        }
      }
      
      write_debug(sprintf("loc=%s yr=%d: Starting variant loop for veg_kept: %s", loc, yr, paste(veg_kept, collapse=", ")))
      top_variants <- list()

      for(v in veg_kept) {
        lib <- OPTIMIZED_LIBRARY[[v]]
        if(is.null(lib)) {
          write_debug(sprintf("loc=%s yr=%d: veg=%s lib is NULL, skipping", loc, yr, v))
          next
        }
        
        write_debug(sprintf("loc=%s yr=%d: veg=%s lib has %d variants", loc, yr, v, 
                           if(!is.null(lib$M_norm)) nrow(lib$M_norm) else 0))

        # Build masked templates from raw library matrix, applying the SAME PCA-LDA transform
        # used in training (if available) so the observation and templates lie in identical
        # feature space before similarity ranking and unmixing.
        # NOTE: lib$M is already pruned (if PRUNE_ZERO_WEIGHT_FEATURES is enabled)
        # during precompute_optimized_library_weighted, so no need to prune again here
        M_full <- lib$M

        # Apply PCA-LDA transform to each variant row (if we have matching weights)
        if (!is.null(PARAMS$weights) && length(PARAMS$weights) == ncol(M_full)) {
          M_full_trans <- t(apply(M_full, 1, function(r) {
            r_num <- as.numeric(r)
            res <- tryCatch(apply_pca_lda_transform(r_num, PARAMS), error = function(e) { if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG] apply_pca_lda_transform failed for lib variant of veg=%s: %s\n", v, e$message)); r_num })
            res
          }))
        } else {
          M_full_trans <- M_full
        }

        # Mask to the valid observation columns
        lib_M_masked <- M_full_trans[, valid_mask, drop = FALSE]

        # Guard against NA/Inf in templates after masking
        lib_M_masked[!is.finite(lib_M_masked)] <- 0

        # Apply barren filtering to vegetation types (same logic as in evaluate_config and heatmap)
        keep_mask <- rep(TRUE, nrow(lib_M_masked))
        if (!is.null(barren_refs_eval_norm) && v != "barren" && nrow(lib_M_masked) > 0) {
          cosine_threshold <- BARREN_SIM_THRESHOLD

          # Normalize templates for cosine similarity
          lib_M_normalized <- t(apply(lib_M_masked, 1, function(r) {
            nrm <- sqrt(sum(r^2))
            if (!is.finite(nrm) || nrm == 0) return(rep(0, length(r))) else return(r / nrm)
          }))

          # Compute cosine similarity to barren references
          sim_mat <- lib_M_normalized %*% t(barren_refs_eval_norm)
          max_sims <- apply(sim_mat, 1, function(x) suppressWarnings(max(x, na.rm = TRUE)))
          keep_mask <- !(max_sims > cosine_threshold)

          # Always keep at least one variant
          if (!any(keep_mask)) {
            min_idx <- which.min(max_sims)
            keep_mask[min_idx] <- TRUE
          }

          n_filtered <- sum(!keep_mask)
          if (n_filtered > 0) {
            write_debug(sprintf("loc=%s yr=%d: veg=%s filtered %d/%d variants (barren similarity > %.2f)", 
                               loc, yr, v, n_filtered, length(keep_mask), cosine_threshold))
          }
        }

        # Apply filter
        lib_M_masked <- lib_M_masked[keep_mask, , drop = FALSE]
        lib_ids_kept <- lib$ids[keep_mask]

        if (nrow(lib_M_masked) == 0) {
          write_debug(sprintf("loc=%s yr=%d: veg=%s has 0 variants after barren filtering, skipping", loc, yr, v))
          next
        }

        # Renormalize each row after masking and filtering (unit L2 norm)
        if (nrow(lib_M_masked) == 1) {
          row_norm <- sqrt(sum(lib_M_masked^2, na.rm = TRUE))
          if (is.finite(row_norm) && row_norm >= 1e-9) lib_M_masked <- lib_M_masked / row_norm
          lib_M_norm_masked <- lib_M_masked
        } else {
          lib_M_norm_masked <- t(apply(lib_M_masked, 1, function(row) {
            row_norm <- sqrt(sum(row^2, na.rm = TRUE))
            if (!is.finite(row_norm) || row_norm < 1e-9) row else row / row_norm
          }))
        }

        # Compute similarities using weighted vectors (for ranking only)
        w_masked <- weights_masked
        if (is.null(w_masked) || length(w_masked) != length(y_for_unmixing)) w_masked <- rep(1, length(y_for_unmixing))
        y_for_sim <- y_for_unmixing * sqrt(pmax(w_masked, 0))
        denom_sim <- sqrt(sum(y_for_sim^2, na.rm = TRUE)); if (denom_sim < 1e-12) denom_sim <- 1
        y_sim_norm <- y_for_sim / denom_sim
        sims <- as.numeric(lib_M_norm_masked %*% y_sim_norm)
        # Include all variants for consideration (remove similarity-based top-K filtering)
        best_idx <- order(sims, decreasing=TRUE)
        write_debug(sprintf("loc=%s yr=%d: veg=%s computed similarities, keeping %d variants", loc, yr, v, length(best_idx)))
        
        top_variants[[v]] <- lapply(best_idx, function(i) {
          # Use masked template (only valid observations) for unmixing
          masked_vec <- lib_M_norm_masked[i, ]
          if (isTRUE(TESTING_MODE) && length(masked_vec) != length(y_for_unmixing)) {
            cat(sprintf("[WARN fit_one_task] loc=%s yr=%d veg=%s: template length mismatch (template=%d, y=%d)\n",
                        loc, yr, v, length(masked_vec), length(y_for_unmixing)))
          }
          list(vec = masked_vec, id = lib_ids_kept[i], similarity = sims[i])
        })
      }

      # Remove empty vegetation types
      empty_vegs <- names(top_variants)[sapply(top_variants, function(x) is.null(x) || length(x) == 0)]
      if (length(empty_vegs) > 0) {
        write_debug(sprintf("loc=%s yr=%d: Removing %d empty veg types: %s", loc, yr, length(empty_vegs), paste(empty_vegs, collapse=", ")))
        for (ev in empty_vegs) top_variants[[ev]] <- NULL
      }
      
      write_debug(sprintf("loc=%s yr=%d: After filtering, top_variants has %d veg types: %s", 
                         loc, yr, length(top_variants), paste(names(top_variants), collapse=", ")))

      if (length(top_variants) == 0) {
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_task] loc=%s yr=%d: no top_variants available after filtering empties, skipping\n", loc, yr))
        write_debug(sprintf("loc=%s yr=%d: no top_variants after filtering - returning NULL", loc, yr))
        return(NULL)
      }
      
      # Evaluate all combinations of endmembers
      # Pass unweighted observation; weighting happens inside solve_weights_fcls
      write_debug(sprintf("loc=%s yr=%d: Calling evaluate_all_combinations with y length=%d, %d veg types", 
                         loc, yr, length(y_for_unmixing), length(top_variants)))
      
      best_result <- tryCatch({
        evaluate_all_combinations(
          y_for_unmixing,
          top_variants,
          lambda = 0,
          feature_weights = weights_masked
        )
      }, error = function(e) {
        write_debug(sprintf("loc=%s yr=%d: evaluate_all_combinations ERROR: %s", loc, yr, e$message))
        cat(sprintf("[ERROR fit_one_task] evaluate_all_combinations failed for loc=%s year=%d: %s\n", as.character(loc), as.integer(yr), e$message))
        NULL
      })
      
      # BEGIN SAFETY BLOCK (Agent added)
      # tryCatch({
      
      write_debug(sprintf("loc=%s yr=%d: evaluate_all_combinations returned %s", 
                         loc, yr, if(is.null(best_result)) "NULL" else "result"))

      if (is.null(best_result)) {
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_task] loc=%s yr=%d: evaluate_all_combinations returned NULL\n", loc, yr))
        write_debug(sprintf("loc=%s yr=%d: evaluate_all_combinations returned NULL", loc, yr))
        return(NULL)
      }
      
      write_debug(sprintf("loc=%s yr=%d: best_result contents - ids=%s, w=%s, rmse=%s", 
                         loc, yr,
                         if(is.null(best_result$ids)) "NULL" else paste(names(best_result$ids), collapse=","),
                         if(is.null(best_result$w)) "NULL" else paste(round(best_result$w,3), collapse=","),
                         if(is.null(best_result$rmse)) "NULL" else round(best_result$rmse,4)))

      # Extract coefficients and create output
      chosen_ids <- best_result$ids
      coefs <- best_result$w
      rmse <- best_result$rmse
      residuals <- best_result$residuals
      E_best_masked <- if (!is.null(best_result$E_best)) best_result$E_best else NULL

      # DIAGNOSTIC: Check coefficient values immediately after extraction
      write_debug(sprintf("loc=%s yr=%d: [DIAGNOSTIC] coefs class=%s, length=%d, values=%s, any NA=%s",
                         loc, yr, class(coefs), length(coefs),
                         paste(round(coefs, 4), collapse=","),
                         any(is.na(coefs))))

      # CRITICAL DEBUG: Print MESMA coefficients to console for diagnosis
      cat(sprintf("[MESMA OUTPUT] loc=%s yr=%d: %s (rmse=%.4f)\n",
                  loc, yr,
                  paste(sprintf("%s=%.4f", names(coefs), coefs), collapse=", "),
                  rmse))

      
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

      # --- Monte Carlo error propagation through unmixing (optional) ---
      mc_coef_sd <- NULL
      mc_coef_q025 <- NULL
      mc_coef_q975 <- NULL

      run_monte_carlo <- isTRUE(ENABLE_UNCERTAINTY) && exists("ENABLE_MONTE_CARLO") && isTRUE(ENABLE_MONTE_CARLO)
      n_mc <- if (exists("MC_N_DRAWS") && is.finite(MC_N_DRAWS)) as.integer(MC_N_DRAWS) else 0L

      if (run_monte_carlo && n_mc >= 10 && !is.null(E_best_masked) && is.matrix(E_best_masked) && ncol(E_best_masked) == length(coefs)) {
        # Noise model: y is z-scored; use residual scale by default.
        mc_sigma <- rmse
        if (exists("MC_NOISE_SD") && is.finite(MC_NOISE_SD) && MC_NOISE_SD > 0) mc_sigma <- as.numeric(MC_NOISE_SD)
        if (exists("MC_NOISE_SCALE") && is.finite(MC_NOISE_SCALE) && MC_NOISE_SCALE > 0) mc_sigma <- mc_sigma * as.numeric(MC_NOISE_SCALE)
        if (!is.finite(mc_sigma) || mc_sigma < 0) mc_sigma <- 0

        # Optional: endmember bundles (sample endmembers from per-band variance estimated from the candidate pool)
        enable_bundles <- exists("ENABLE_ENDMEMBER_BUNDLES") && isTRUE(ENABLE_ENDMEMBER_BUNDLES)
        bundle_params <- NULL
        if (enable_bundles) {
          bundle_params <- list()
          verbose_bundle <- isTRUE(DEBUG_UNCERTAINTY) || isTRUE(TESTING_MODE)

          for (v in names(chosen_ids)) {
            cand_list <- top_variants[[v]]
            if (is.null(cand_list) || length(cand_list) < 3) next
            Mv <- do.call(rbind, lapply(cand_list, function(z) as.numeric(z$vec)))
            if (is.null(Mv) || nrow(Mv) < 3) next

            bundle_result <- compute_bundle_covariance(Mv, verbose = verbose_bundle)

            if (!is.null(bundle_result) && !is.null(bundle_result$A)) {
              bundle_params[[v]] <- list(
                mu = bundle_result$mu,
                A = bundle_result$A
              )
            }
          }
        }

        sample_from_bundle <- function(mu, A) {
          z <- stats::rnorm(length(mu))
          x <- as.numeric(mu + A %*% z)
          x[!is.finite(x)] <- 0
          nrm <- sqrt(sum(x^2))
          if (!is.finite(nrm) || nrm < 1e-9) return(mu)
          x / nrm
        }

        # Check if OOB fraction uncertainty is enabled and residuals are available
        enable_oob_frac_uncertainty <- exists("ENABLE_OOB_FRACTION_UNCERTAINTY") &&
                                        isTRUE(ENABLE_OOB_FRACTION_UNCERTAINTY) &&
                                        exists(".OOB_FRACTION_RESIDUALS", envir = globalenv())

        # Determine dominant class for OOB residual sampling (class with highest coefficient)
        dominant_class <- NULL
        if (enable_oob_frac_uncertainty && length(coefs) > 0) {
          dominant_class <- names(which.max(coefs))
        }

        w_draws <- matrix(NA_real_, nrow = length(coefs), ncol = n_mc)
        rownames(w_draws) <- names(coefs)

        for (b in seq_len(n_mc)) {
          y_mc <- y_for_unmixing
          if (length(mc_sigma) == 1 && is.finite(mc_sigma) && mc_sigma > 0) {
            y_mc <- y_mc + stats::rnorm(length(y_mc), mean = 0, sd = mc_sigma)
          }

          E_mc <- E_best_masked
          if (enable_bundles && !is.null(bundle_params) && length(bundle_params) > 0) {
            for (j in seq_along(names(chosen_ids))) {
              v <- names(chosen_ids)[j]
              bp <- bundle_params[[v]]
              if (!is.null(bp) && !is.null(bp$mu) && !is.null(bp$A)) {
                E_mc[, j] <- sample_from_bundle(bp$mu, bp$A)
              }
            }
          }

          res_mc <- solve_weights_fcls(E_mc, y_mc, feature_weights = weights_masked)
          if (!is.null(res_mc) && !is.null(res_mc$w) && length(res_mc$w) == length(coefs)) {
            w_mc <- as.numeric(res_mc$w)
            names(w_mc) <- names(coefs)

            # Apply OOB fraction residual perturbation if enabled
            if (enable_oob_frac_uncertainty && !is.null(dominant_class)) {
              oob_resid <- sample_oob_residual(dominant_class)
              if (!is.null(oob_resid)) {
                # Match residual names to coefficient names and add perturbation
                for (vname in names(w_mc)) {
                  vname_lower <- tolower(vname)
                  resid_names_lower <- tolower(names(oob_resid))
                  match_idx <- which(resid_names_lower == vname_lower)
                  if (length(match_idx) > 0) {
                    w_mc[vname] <- w_mc[vname] + oob_resid[match_idx[1]]
                  }
                }
                # Re-enforce sum-to-one and non-negativity constraints
                w_mc[w_mc < 0] <- 0
                w_sum <- sum(w_mc)
                if (w_sum > 1e-9) {
                  w_mc <- w_mc / w_sum
                }
              }
            }

            w_draws[, b] <- w_mc
          }
        }

        if (ncol(w_draws) > 1) {
          mc_coef_sd <- apply(w_draws, 1, stats::sd, na.rm = TRUE)
          mc_coef_q025 <- apply(w_draws, 1, stats::quantile, probs = 0.025, na.rm = TRUE)
          mc_coef_q975 <- apply(w_draws, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
        }
      }

      # Build coefficient dataframe with variant-level detail
      # The veg names are stored as names of chosen_ids (which is a named vector)
      # Values of chosen_ids are variant IDs like "populus_opt_1", but can be NA if veg type not used
      early_zero_applied <- FALSE

      veg_names <- names(chosen_ids)
      if (is.null(veg_names) || length(veg_names) == 0) {
        # Fallback: try to extract from variant IDs if names not available
        veg_names <- sapply(strsplit(as.character(chosen_ids), "_(v|opt)"), function(x) {
          if (length(x) > 0 && !is.na(x[1]) && nchar(x[1]) > 0) x[1] else NA_character_
        })
      }
      coef_df <- data.frame(
        location_id = loc,
        pheno_year = yr,
        lat = lat_val,
        lon = lon_val,
        Veg = veg_names,
        variant_id = as.character(chosen_ids),
        coef = coefs,
        rmse = rmse,
        coef_median = NA_real_,
        coef_025 = if (!is.null(mc_coef_q025)) as.numeric(mc_coef_q025) else NA,
        coef_975 = if (!is.null(mc_coef_q975)) as.numeric(mc_coef_q975) else NA,
        coef_sd = {
          base_sd <- as.numeric(coef_sd_vec)
          if (!is.null(mc_coef_sd)) {
            base_sd <- sqrt(pmax(base_sd, 0)^2 + pmax(as.numeric(mc_coef_sd), 0)^2)
          }
          base_sd
        },
        interval = NA,
        n_obs = nrow(task_data),
        inseparable_variant_flag = FALSE,
        inseparable_variant_details = NA_character_,
        stringsAsFactors = FALSE
      )

      # DIAGNOSTIC: Check coef_df immediately after creation
      write_debug(sprintf("loc=%s yr=%d: [DIAGNOSTIC] coef_df created with %d rows, coef class=%s, coef values=%s",
                         loc, yr, nrow(coef_df), class(coef_df$coef),
                         paste(sprintf("%s=%.4f", coef_df$Veg, coef_df$coef), collapse=", ")))

      # --- Rename high-similarity populus/tamarix variants to "woody_unknown" ---
      # Variants with cross-class similarity > 0.95 between populus and tamarix are indistinguishable
      if (exists("WOODY_UNKNOWN_VARIANTS", envir = globalenv())) {
        woody_unknown_list <- get("WOODY_UNKNOWN_VARIANTS", envir = globalenv())
        if (length(woody_unknown_list) > 0) {
          rename_mask <- coef_df$variant_id %in% woody_unknown_list
          rename_mask[is.na(rename_mask)] <- FALSE
          if (any(rename_mask)) {
            n_renamed <- sum(rename_mask)
            old_vegs <- coef_df$Veg[rename_mask]
            coef_df$Veg[rename_mask] <- "woody_unknown"
            write_debug(sprintf("loc=%s yr=%d: Renamed %d variants to 'woody_unknown' (originally: %s)",
                               loc, yr, n_renamed, paste(unique(old_vegs), collapse=", ")))
          }
        }
      }

      # --- Mark inseparable variants (drop instead of assign Veg = 'unknown') if detected in similarity tables ---
      # Wrapped in tryCatch to prevent non-critical metadata lookups from crashing the task
      tryCatch({
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
              mask[is.na(mask)] <- FALSE
              if (any(mask, na.rm = TRUE)) {
                coef_df$inseparable_variant_flag[mask] <- TRUE
                # Build detail strings: collect related vegs and other variant ids
                details <- vapply(coef_df$variant_id[mask], function(vid) {
                  rows <- apply(sim_tbl, 1, function(r) any(vid == as.character(r), na.rm = TRUE))
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
                    if (!(isTRUE(a_id %in% coef_df$variant_id) && isTRUE(b_id %in% coef_df$variant_id))) next

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

                    idx_drop <- which(coef_df$variant_id == drop_id)
                    if (length(idx_drop) > 0) {
                      coef_df$inseparable_variant_flag[idx_drop] <- TRUE
                      coef_df$inseparable_variant_details[idx_drop] <- paste0("High similarity (", sprintf("%.3f", sim_val), ") with barren variant ", other_id, " - kept")
                      if (isTRUE(TESTING_MODE)) cat(sprintf("[INFO] Flagged variant %s (similarity=%.3f > 0.95 with barren %s) but keeping coefficient\n", drop_id, sim_val, other_id))
                    }
                  }
                }
              }
            }
          }
      }, error = function(e) {
          if (isTRUE(TESTING_MODE)) cat(sprintf("[WARN fit_one_task] loc=%s yr=%d: Inseparable variant check failed: %s\n", loc, yr, e$message))
      })

      # Barren is kept as from MESMA - no replacement
      # Aggregate by vegetation type (sum coefficients for same veg type)
      # Must happen AFTER all filtering but BEFORE creating diagnostics/coef_agg usage
      write_debug(sprintf("loc=%s yr=%d: [DEBUG 3] Before aggregation: nrow(coef_df)=%d", loc, yr, nrow(coef_df)))
      if (nrow(coef_df) > 0) {
        # DIAGNOSTIC: Check coef values before aggregation
        write_debug(sprintf("loc=%s yr=%d: [DIAGNOSTIC] Before aggregate - coef_df$coef class=%s, any NA=%s, values=%s",
                           loc, yr, class(coef_df$coef), any(is.na(coef_df$coef)),
                           paste(sprintf("%s=%.4f", coef_df$Veg, coef_df$coef), collapse=", ")))
        coef_agg <- aggregate(coef ~ Veg, data = coef_df, FUN = sum)
        # DIAGNOSTIC: Check coef values after aggregation
        write_debug(sprintf("loc=%s yr=%d: [DIAGNOSTIC] After aggregate - coef_agg$coef class=%s, any NA=%s, values=%s",
                           loc, yr, class(coef_agg$coef), any(is.na(coef_agg$coef)),
                           paste(sprintf("%s=%.4f", coef_agg$Veg, coef_agg$coef), collapse=", ")))
      } else {
        # CRITICAL FIX: When coef_df is empty (e.g., MESMA picked 100% barren),
        # create coef_agg with vegetation types from library, initialized to 0
        # This allows PPI scaling logic to distribute vegetation cover among available types
        veg_types_nonbarren <- setdiff(veg_kept, "barren")
        if (length(veg_types_nonbarren) > 0) {
          coef_agg <- data.frame(
            Veg = veg_types_nonbarren,
            coef = rep(0, length(veg_types_nonbarren)),
            stringsAsFactors = FALSE
          )
        } else {
          coef_agg <- data.frame(Veg = character(0), coef = numeric(0), stringsAsFactors = FALSE)
        }
      }
      write_debug(sprintf("loc=%s yr=%d: [DEBUG 4] After aggregation: nrow(coef_agg)=%d", loc, yr, nrow(coef_agg)))

      # --- Scale vegetation fractions by index-derived total vegetation cover ---
      # MESMA provides relative proportions among vegetation types.
      # We multiply the vegetation fractions by total_veg_cover (from MSAVI) to get absolute fractions.
      # Barren fraction from MESMA is preserved - we do NOT replace it with index-derived estimate.
      # NO equal distribution fallback - if MESMA says 0 vegetation, we keep it as 0.

      veg_coefs_mask <- tolower(coef_agg$Veg) != "barren"
      veg_coefs_mask[is.na(veg_coefs_mask)] <- FALSE
      sum_original_veg_coefs <- sum(coef_agg$coef[veg_coefs_mask], na.rm = TRUE)

      # Scale vegetation fractions: multiply by total_veg_cover / sum_of_veg_fractions
      # This preserves MESMA's relative proportions while scaling to index-derived total veg cover
      scale_factor <- 1.0
      if (is.finite(sum_original_veg_coefs) && sum_original_veg_coefs > 1e-9) {
        scale_factor <- total_veg_cover / sum_original_veg_coefs
        coef_agg$coef[veg_coefs_mask] <- coef_agg$coef[veg_coefs_mask] * scale_factor
        write_debug(sprintf("loc=%s yr=%d: Scaled veg fractions by %.4f (total_veg_cover=%.4f, sum_mesma_veg=%.4f)",
                           loc, yr, scale_factor, total_veg_cover, sum_original_veg_coefs))
      } else {
        # MESMA returned 0 vegetation - keep it as 0, don't fabricate equal fractions
        write_debug(sprintf("loc=%s yr=%d: MESMA returned 0 vegetation (sum=%.6f), keeping as-is (no equal distribution fallback)",
                           loc, yr, sum_original_veg_coefs))
      }

      # Keep MESMA's barren fraction - do NOT replace with index-derived estimate
      # Just ensure barren row exists in coef_agg
      barren_idx <- which(tolower(coef_agg$Veg) == "barren")
      if (length(barren_idx) == 0) {
        # Add barren row from MESMA if it doesn't exist (use 1 - sum of scaled veg fractions)
        mesma_barren <- 1 - sum(coef_agg$coef[veg_coefs_mask], na.rm = TRUE)
        mesma_barren <- pmax(0, pmin(1, mesma_barren))
        barren_row <- data.frame(Veg = "barren", coef = mesma_barren)
        coef_agg <- rbind(coef_agg, barren_row)
        write_debug(sprintf("loc=%s yr=%d: Added barren row with fraction %.4f (1 - sum_veg)", loc, yr, mesma_barren))
      }
      # Note: If barren already exists in coef_agg, we keep MESMA's original value

      # --- Propagate scaling to coef_df (no equal distribution fallback) ---
      write_debug(sprintf("loc=%s yr=%d: [DEBUG 5] Before scaling: nrow(coef_df)=%d, sum_original_veg_coefs=%.4f", loc, yr, nrow(coef_df), sum_original_veg_coefs))

      # Wrap scaling in a defensive tryCatch to diagnose errors
      tryCatch({
        # Apply the same scale_factor to coef_df vegetation fractions
        # NO rebuilding coef_df from coef_agg - if MESMA returned empty, keep it empty
        if (nrow(coef_df) > 0 && is.finite(sum_original_veg_coefs) && sum_original_veg_coefs > 1e-9) {
          v_rows <- tolower(coef_df$Veg) != "barren"
          v_rows[is.na(v_rows)] <- FALSE
          coef_df$coef[v_rows] <- coef_df$coef[v_rows] * scale_factor
          if ("coef_sd" %in% names(coef_df)) coef_df$coef_sd[v_rows] <- coef_df$coef_sd[v_rows] * abs(scale_factor)
          if ("coef_025" %in% names(coef_df)) coef_df$coef_025[v_rows] <- coef_df$coef_025[v_rows] * scale_factor
          if ("coef_975" %in% names(coef_df)) coef_df$coef_975[v_rows] <- coef_df$coef_975[v_rows] * scale_factor
        }
        # If MESMA returned 0 vegetation, coef_df veg fractions stay at 0 - no fallback
        write_debug(sprintf("loc=%s yr=%d: [DEBUG 6] After scaling: nrow(coef_df)=%d", loc, yr, nrow(coef_df)))
      }, error = function(e) {
        # Detailed diagnostics to help track down TRUE/FALSE NA issues
        write_debug(sprintf("loc=%s yr=%d: ERROR in scaling block: %s", loc, yr, e$message))
        write_debug(sprintf("loc=%s yr=%d: DIAG: sum_original_veg_coefs=%s total_veg_cover=%s nrow(coef_df)=%d nrow(coef_agg)=%d", loc, yr, as.character(sum_original_veg_coefs), as.character(total_veg_cover), if(exists("coef_df")) nrow(coef_df) else -1, if(exists("coef_agg")) nrow(coef_agg) else -1))
        if (exists("coef_df")) {
          write_debug(sprintf("loc=%s yr=%d: DIAG: coef_df head Veg: %s", loc, yr, paste(head(as.character(coef_df$Veg), 5), collapse=",")))
          write_debug(sprintf("loc=%s yr=%d: DIAG: coef_df head variant_id: %s", loc, yr, paste(head(as.character(coef_df$variant_id), 5), collapse=",")))
        }
        stop(e)
      })
      
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
        if ("coef_sd" %in% names(barren_row_df)) barren_row_df$coef_sd <- 0
        if ("coef_025" %in% names(barren_row_df)) barren_row_df$coef_025 <- NA
        if ("coef_975" %in% names(barren_row_df)) barren_row_df$coef_975 <- NA
        coef_df <- rbind(coef_df, barren_row_df)
      } else {
        # No vegetation detected, use PPI barren fraction
        if (isTRUE(is.finite(barren_fraction) && barren_fraction > 0)) {
          # Return barren only result
          coef_df <- data.frame(
            location_id = loc,
            pheno_year = yr,
            lat = lat_val,
            lon = lon_val,
            Veg = "barren",
            variant_id = "barren_ppi",
            coef = barren_fraction,
            rmse = if(exists("rmse") && is.numeric(rmse)) rmse else 0,
            coef_025 = NA,
            coef_975 = NA,
            coef_sd = 0,
            interval = NA,
            n_obs = nrow(task_data),
            inseparable_variant_flag = FALSE,
            inseparable_variant_details = NA_character_,
            stringsAsFactors = FALSE
          )
        } else {
          # No vegetation and no barren - fail
          write_debug(sprintf("loc=%s yr=%d: No vegetation detected (100%% barren solution) - returning NULL", loc, yr))
          if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_task] loc=%s yr=%d: coef_df empty after barren removal\n", loc, yr))
          return(NULL)
        }
      }

      # Create diagnostics dataframe
      diag_df <- data.frame(
        location_id = loc,
        pheno_year = yr,
        stringsAsFactors = FALSE
      )
      # Add each vegetation type's total fraction
      for (v in coef_agg$Veg) {
        if (!is.na(v)) {
          v_mask <- coef_agg$Veg == v
          v_mask[is.na(v_mask)] <- FALSE
          diag_df[[paste0(v, "_fraction")]] <- coef_agg$coef[v_mask][1]
        }
      }
      diag_df$barren_fraction_ppi_based <- barren_fraction # Add PPI-based barren fraction

      # Reconstruct E_best for returning
      E_best_list <- list()
      for(v in names(chosen_ids)){
        vid <- chosen_ids[[v]]
        # Defensive: ensure we have candidate variants for this veg
        if (is.null(top_variants[[v]]) || length(top_variants[[v]]) == 0) {
          write_debug(sprintf("loc=%s yr=%d: No candidate variants present for veg=%s when reconstructing E_best", loc, yr, v))
          next
        }
        found <- FALSE
        for(variant in top_variants[[v]]){
          # Defensive equality check to avoid NA logicals causing crashes
          if (!is.na(variant$id) && !is.na(vid) && isTRUE(variant$id == vid)) {
            E_best_list[[v]] <- variant$vec
            found <- TRUE
            break
          }
        }
        # Fallback: if chosen variant not found (e.g., NA mismatch), use first candidate and log
        if (!found) {
          E_best_list[[v]] <- top_variants[[v]][[1]]$vec
          write_debug(sprintf("loc=%s yr=%d: Warning - chosen variant '%s' not found for veg=%s; using first candidate variant_id=%s", loc, yr, as.character(vid), v, as.character(top_variants[[v]][[1]]$id)))
        }
      }
      E_best <- do.call(cbind, E_best_list)

      write_debug(sprintf("loc=%s yr=%d: SUCCESS - preparing to return coef_df with %d rows", loc, yr, nrow(coef_df)))
      # DEBUG: Log coefficient values before returning
      if (nrow(coef_df) > 0) {
        coef_summary <- sprintf("coefs: %s", paste(sprintf("%s=%.4f", coef_df$Veg, coef_df$coef), collapse=", "))
        write_debug(sprintf("loc=%s yr=%d: %s", loc, yr, coef_summary))
        na_count <- sum(is.na(coef_df$coef))
        if (na_count > 0) {
          write_debug(sprintf("loc=%s yr=%d: WARNING - %d/%d coefficients are NA!", loc, yr, na_count, nrow(coef_df)))
        }
      }
      return(tryCatch({
        list(
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
        )
      }, error = function(e) {
        write_debug(sprintf("loc=%s yr=%d: ERROR in final return: %s", loc, yr, e$message))
        if (isTRUE(TESTING_MODE)) cat(sprintf("[ERROR fit_one_task] loc=%s yr=%d: final return failed: %s\n", loc, yr, e$message))
        stop(e)
      }))
      
      # }, error = function(e) {
      #    write_debug(sprintf("loc=%s yr=%d: ERROR in post-processing: %s", loc, yr, e$message))
      #    # Also print to stdout to be visible in terminal
      #    cat(sprintf("[CRITICAL ERROR fit_one_task] Post-processing failed for loc=%s yr=%d: %s\n", loc, yr, e$message))
      #    return(NULL)
      # })

      # ===== END MESMA UNMIXING =====
  # These variables are already set in the earlier conditional blocks above.
  # No need to reassign them here
  # Single-stage MESMA is used, so mesma_lib and OPTIMIZED_LIBRARY are set accordingly

  }

# New function: Process all years for a single location with multi-year bootstrap
  fit_one_location <- function(location_data) {
    # File-based debug log that won't be suppressed by sink()
    debug_log_file <- file.path(tempdir(), "fit_one_location_debug.log")
    write_debug <- function(msg) {
      tryCatch(write(paste0(Sys.time(), " - ", msg), debug_log_file, append = TRUE), error = function(e) {})
    }
    
    if (is.null(location_data) || nrow(location_data) == 0) {
      write_debug("Received NULL or empty location_data")
      if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_location] Received NULL or empty location_data\n"))
      return(NULL)
    }

    loc <- as.character(location_data$location_id[1])
    write_debug(sprintf("Processing location: %s, rows: %d", loc, nrow(location_data)))
    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_location] Processing location: %s, rows: %d\n", loc, nrow(location_data)))

    # Get all phenological years for this location
    years <- sort(unique(location_data$pheno_year))
    years <- years[!is.na(years)]
    
    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_location] Location %s has years: %s\n", loc, paste(years, collapse=", ")))

    if (length(years) == 0) {
      if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_location] Location %s has no valid years - returning NULL\n", loc))
      return(NULL)
    }

    # Process each year individually first (to get point estimates and chosen variants)
    year_results <- list()
    y_vecs_by_year <- list()
    chosen_ids_by_year <- list()
    w_hat_by_year <- list()

    # Wrap per-location processing in tryCatch so a single error does not discard valid year results
    out <- tryCatch({
      for (yr in years) {
        year_data <- location_data[location_data$pheno_year == yr, , drop = FALSE]
        write_debug(sprintf("loc=%s yr=%d: year_data has %d rows", loc, yr, nrow(year_data)))
        
        # === SPLINE MODE: Use spline-based unmixing if enabled ===
        if (isTRUE(USE_SPLINE_ENDMEMBERS) && exists("SPLINE_LIBRARY") && exists("SPLINE_PARAMS")) {
          res_yr <- tryCatch({
            fit_one_task_spline(year_data, SPLINE_LIBRARY, SPLINE_PARAMS$indices, SPLINE_PARAMS, loc, yr)
          }, error = function(e) {
            if (isTRUE(TESTING_MODE)) cat(sprintf("[ERROR fit_one_task_spline] loc=%s yr=%d: %s\n", loc, yr, e$message))
            NULL
          })
        } else {
          # Standard MESMA unmixing
          res_yr <- tryCatch({
            fit_one_task(year_data)
          }, error = function(e) {
            msg <- sprintf("loc=%s yr=%d: ERROR in fit_one_task: %s", loc, yr, e$message)
            write_debug(msg)
            message(paste("[ERROR fit_one_task CRASH]", msg))
            # Capture a compact call stack to pinpoint the failing if() condition
            calls <- tryCatch(sys.calls(), error = function(e2) NULL)
            if (!is.null(calls)) {
              tail_calls <- tail(calls, 12)
              call_str <- paste(vapply(tail_calls, function(x) paste(deparse(x), collapse = ""), character(1)), collapse = " | ")
              write_debug(sprintf("loc=%s yr=%d: call stack (tail): %s", loc, yr, call_str))
            }
            NULL
          })
        }
        # === END SPLINE MODE BRANCH ===

        write_debug(sprintf("loc=%s yr=%d: fit_one_task returned %s", loc, yr, if(is.null(res_yr)) "NULL" else "result"))
        if (!is.null(res_yr)) {
          tryCatch({
            year_results[[as.character(yr)]] <- res_yr
            write_debug(sprintf("loc=%s yr=%d: Added to year_results (now has %d entries)", loc, yr, length(year_results)))
          }, error = function(e) {
            if (isTRUE(TESTING_MODE)) cat(sprintf("[ERROR assigning year_results] loc=%s yr=%d: %s\n", loc, yr, e$message))
            stop(e)  # Re-throw to be caught by outer handler
          })

          # Store data needed for multi-year bootstrap - wrapped in tryCatch to isolate errors
          tryCatch({
            # For inference/validation data, do NOT interpolate - use only actual observations
            raw_mat_yr <- build_pentad_matrix(year_data, MESMA_PARAMS$indices, interpolate = FALSE)
            if (!is.null(raw_mat_yr)) {
              y_raw_yr <- as.numeric(raw_mat_yr)

              # Check dimensions match expectations
              n_bins <- TEMPORAL_BUDGET
              expected_length <- length(MESMA_PARAMS$indices) * n_bins
              if (length(y_raw_yr) != expected_length) {
                if (isTRUE(TESTING_MODE)) {
                  cat(sprintf("[WARN fit_one_location] loc=%s yr=%d: y_raw_yr length mismatch (got %d, expected %d), skipping multi-year processing for this year\n",
                    loc, yr, length(y_raw_yr), expected_length))
                }
              } else {
                # Apply same preprocessing as fit_one_task
                y_s1_yr <- y_raw_yr
                for(k in seq_along(MESMA_PARAMS$indices)) {
                  idx_start <- (k-1)*n_bins + 1
                  idx_end <- k*n_bins
                  mu <- MESMA_PARAMS$means[k]
                  sigma <- MESMA_PARAMS$sds[k]
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
          write_debug(sprintf("loc=%s yr=%d: fit_one_task returned NULL", loc, yr))
        }
      }

      write_debug(sprintf("loc=%s: After year loop, year_results has %d entries", loc, length(year_results)))
      # If we have results for multiple years and ENABLE_UNCERTAINTY, do multi-year bootstrap
      if (length(year_results) == 0) {
        write_debug(sprintf("loc=%s: year_results empty - returning NULL", loc))
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_location] loc=%s: no successful year results\n", loc))
        return(NULL)
      }

      if (isTRUE(ENABLE_UNCERTAINTY) && isTRUE(ENABLE_MULTI_YEAR_BOOTSTRAP) && length(years) >= 1 && length(y_vecs_by_year) == length(years)) {
        # Get comp_templates from first successful year result
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
            con_out <- tryCatch(file(sink_file, open = "wt"), error = function(e) NULL)
            con_msg <- tryCatch(file(paste0(sink_file, ".msg"), open = "wt"), error = function(e) NULL)
            tryCatch({
              if (!is.null(con_out) && inherits(con_out, "connection") && isOpen(con_out)) {
                sink(con_out, type = "output")
              }
              if (!is.null(con_msg) && inherits(con_msg, "connection") && isOpen(con_msg)) {
                sink(con_msg, type = "message")
              }

              for (b in seq_len(B_loc)) {
                
                # === SPLINE MODE: Bootstrap on per-observation residuals ===
                if (isTRUE(USE_SPLINE_ENDMEMBERS) && !is.null(res_yr$obs_residuals)) {
                  # Bootstrap residuals using simple i.i.d. resampling
                  residuals_boot <- simple_residual_bootstrap(res_yr$residuals)
                  
                  # Refit each observation with bootstrapped residuals added to the original observations
                  n_obs_yr <- length(res_yr$obs_residuals)
                  boot_coefs_obs <- matrix(NA, nrow = n_obs_yr, ncol = length(res_yr$mean_coefs))
                  colnames(boot_coefs_obs) <- names(res_yr$mean_coefs)
                  
                  resid_ptr <- 1  # Pointer into the concatenated bootstrapped residuals
                  doy_order <- order(res_yr$doys, na.last = TRUE)  # Same order as residuals were concatenated
                  
                  for (obs_idx in doy_order) {
                    obs_E <- res_yr$obs_E[[obs_idx]]
                    obs_y <- res_yr$obs_y[[obs_idx]]
                    
                    if (is.null(obs_E) || is.null(obs_y) || length(obs_y) == 0) next
                    
                    # Get the bootstrapped residuals for this observation
                    n_valid <- length(obs_y)
                    if (resid_ptr + n_valid - 1 > length(residuals_boot)) break
                    
                    boot_resid <- residuals_boot[resid_ptr:(resid_ptr + n_valid - 1)]
                    resid_ptr <- resid_ptr + n_valid
                    
                    # Add bootstrapped residuals to observed y to create perturbed observation
                    y_boot_obs <- obs_y + boot_resid
                    
                    # Re-solve FCLS with perturbed observations
                    w_boot <- tryCatch({
                      solve_weights_fcls_spline(obs_E, y_boot_obs)
                    }, error = function(e) NULL)
                    
                    if (!is.null(w_boot) && length(w_boot) == ncol(boot_coefs_obs)) {
                      boot_coefs_obs[obs_idx, ] <- w_boot
                    }
                  }
                  
                  # Aggregate bootstrapped coefficients
                  boot_mean_coefs <- colMeans(boot_coefs_obs, na.rm = TRUE)
                  if (sum(boot_mean_coefs, na.rm = TRUE) > 0) {
                    boot_mean_coefs <- boot_mean_coefs / sum(boot_mean_coefs, na.rm = TRUE)
                  }
                  
                  # Collect bootstrapped coefficients by vegetation class
                  for (ci in seq_along(names(boot_mean_coefs))) {
                    vg <- names(boot_mean_coefs)[ci]
                    if (is.na(vg)) next
                    val <- boot_mean_coefs[ci]
                    if (!is.finite(val)) next
                    if (is.null(boot_coef_by_veg[[vg]])) boot_coef_by_veg[[vg]] <- numeric(0)
                    boot_coef_by_veg[[vg]] <- c(boot_coef_by_veg[[vg]], val)
                  }
                  
                  # Compute bootstrapped RMSE for spline mode
                  # RMSE is the standard deviation of residuals after bootstrap perturbation
                  if (!is.null(res_yr$mean_rmse) && is.finite(res_yr$mean_rmse)) {
                    # Add variability to RMSE estimate based on bootstrapped residuals
                    boot_rmse_val <- sqrt(mean(residuals_boot^2, na.rm = TRUE))
                    if (is.finite(boot_rmse_val)) {
                      boot_rmse <- c(boot_rmse, boot_rmse_val)
                    }
                  }
                  
                  next  # Skip the standard MESMA bootstrap below
                }
                # === END SPLINE MODE ===

                # Bootstrap residuals using simple i.i.d. resampling
                residuals_boot <- simple_residual_bootstrap(res_yr$residuals)
                y_boot <- res_yr$y_hat + residuals_boot

                # OPTIMIZATION: Reuse the endmember matrix from the initial fit
                # Instead of searching all variant combinations again, just re-solve weights
                # with the same E_best matrix. This is ~100-1000x faster per bootstrap rep.
                E_boot <- res_yr$E_best
                if (is.null(E_boot) || !is.matrix(E_boot) || ncol(E_boot) < 1) {
                  next  # Skip if no valid endmember matrix
                }

                res_b <- tryCatch({
                  solve_weights_fcls(E_boot, y_boot, feature_weights = res_yr$weights_masked)
                }, error = function(e) NULL)

                if (is.null(res_b)) next

                coefs_b <- res_b$w
                # Use the original chosen variant IDs from the initial fit
                chosen_ids_b <- res_yr$coef_df$variant_id
                names(chosen_ids_b) <- res_yr$coef_df$Veg
                names(coefs_b) <- res_yr$coef_df$Veg
                
                dfb <- data.frame(
                    Veg = sapply(strsplit(chosen_ids_b, "_(v|opt)"), `[`, 1),
                    coef = coefs_b,
                    stringsAsFactors = FALSE
                )

                # collect coefficients
                for (i in seq_len(nrow(dfb))) {
                  vg <- as.character(dfb$Veg[i])
                  if (is.na(vg)) next
                  val <- as.numeric(dfb$coef[i])
                  if (!is.finite(val)) next
                  if (is.null(boot_coef_by_veg[[vg]])) boot_coef_by_veg[[vg]] <- numeric(0)
                  boot_coef_by_veg[[vg]] <- c(boot_coef_by_veg[[vg]], val)
                }

                # collect chosen variant ids
                for (i in seq_along(chosen_ids_b)) {
                  vg <- names(chosen_ids_b)[i]
                  if (is.na(vg) || !is.character(vg)) next
                  vid <- chosen_ids_b[[i]]
                  if (is.na(vid)) next
                  if (is.null(boot_variant_counts[[vg]])) boot_variant_counts[[vg]] <- list()
                  boot_variant_counts[[vg]][[as.character(vid)]] <- (boot_variant_counts[[vg]][[as.character(vid)]] %||% 0) + 1L
                }

                # rmse
                if (!is.null(res_b$rmse) && is.finite(as.numeric(res_b$rmse))) boot_rmse <- c(boot_rmse, as.numeric(res_b$rmse))
              }

            }, finally = {
              try(sink(type = "message"), silent = TRUE)
              try(sink(type = "output"), silent = TRUE)
              if (!is.null(con_out)) try(close(con_out), silent = TRUE)
              if (!is.null(con_msg)) try(close(con_msg), silent = TRUE)
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
                           coef_median = as.numeric(median(vals, na.rm = TRUE)),
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
                              # Store median bootstrap estimate for plotting/prediction preference
                              if ("coef_median" %in% names(coef_ci_df)) {
                                curr_df$coef_median[match_idx] <- coef_ci_df$coef_median[r_idx]
                              }
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

      write_debug(sprintf("Returning %d year results", length(year_results)))
      # Return results for all years
      return(year_results)

    }, error = function(e) {
      write_debug(sprintf("ERROR in fit_one_location: %s", e$message))
      cat(sprintf("[ERROR fit_one_location] loc=%s: %s\n", loc, e$message))
      if (length(year_results) > 0) {
        if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG fit_one_location] loc=%s: returning partial %d year(s) after error\n", loc, length(year_results)))
        return(year_results)
      }
      return(NULL)
    })

    write_debug(sprintf("Final return for %s: out is %s with %d items", loc, 
                        if(is.null(out)) "NULL" else "non-NULL",
                        if(is.null(out)) 0 else length(out)))
    out
  }

  # Training is disabled; do not short-circuit validation here. Initialize
  # validation containers and let the validation routines run normally below.
  validation_coefs <- data.frame()
  validation_results_list <- list()

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
    
    # Use spline mode if enabled
    if (isTRUE(USE_SPLINE_ENDMEMBERS) && exists("SPLINE_LIBRARY") && exists("SPLINE_PARAMS")) {
      if (isTRUE(TESTING_MODE)) cat('[DEBUG] Using SPLINE MODE\n')
      res <- tryCatch({ 
        fit_one_task_spline(task_data, SPLINE_LIBRARY, SPLINE_PARAMS$indices, SPLINE_PARAMS, loc, yr)
      }, error = function(e) { cat(sprintf('DEBUG UNMIX ERROR (spline): %s\n', e$message)); NULL })
    } else {
      res <- tryCatch({ fit_one_task(task_data) }, error = function(e) { cat(sprintf('DEBUG UNMIX ERROR: %s\n', e$message)); NULL })
    }
    
    if (!is.null(res)) {
      cat('--- RESULT SUMMARY ---\n')
      if (!is.null(res$diagnostics)) print(res$diagnostics)
      if (!is.null(res$coef_df)) print(head(res$coef_df, n = 20))
    }
    cat('=== END DEBUG UNMIX ===\n\n')
    res
  }
  cat("Starting main processing loop...\n")
  
  # === SPLINE MODE STATUS ===
  if (isTRUE(USE_SPLINE_ENDMEMBERS)) {
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════════════╗\n")
    cat("║  SPLINE ENDMEMBER MODE ACTIVE                                    ║\n")
    cat("║  Each class represented by smooth spline curves across DOY       ║\n")
    cat("║  Unmixing performed at observation-level with DOY-specific       ║\n")
    cat("║  endmember values from spline functions                          ║\n")
    cat("╚══════════════════════════════════════════════════════════════════╝\n")
    cat("\n")
  }
  # === END SPLINE MODE STATUS ===

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

aggregate_to_global_pattern <- function(all_coefs, method = "location_bootstrap") {
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")
    # Debug: print method to diagnose "Unknown method" errors
    if (exists("TESTING_MODE") && isTRUE(TESTING_MODE)) {
      cat(sprintf("[DEBUG] aggregate_to_global_pattern called with method='%s'\n", method))
    }
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
      add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
      ggplot2::geom_line(size = 1.2) +
      ggplot2::geom_point(size = 2) +
      ggplot2::coord_cartesian(ylim = c(0, NA))
    
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
      ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      ggplot2::scale_fill_brewer(palette = "Set2")
    p
  }

# Executing main processing steps directly (function removed). This section previously defined
# `main_processing_block <- function() { ... }`. To avoid large function compile stalls, the body
# is now executed at top-level during script run. The original nested helper functions remain in scope.

    # Assign df_tasks for the main processing loop (inference data, not training data)
    if (exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0) {
      df_tasks <- df_tasks_inference
    } else {
      df_tasks <- data.frame()  # Empty if no inference data
    }

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
  
  if (!isTRUE(QUIET_MODE)) {
    cat("Preparing locations for batched processing (multi-year bootstrap)...\n")

    # DEBUG: Check what columns df_tasks has
    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG VALIDATION DATA] df_tasks has %d rows, %d columns\n", nrow(df_tasks), ncol(df_tasks)))
    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG VALIDATION DATA] Column names: %s\n", paste(head(names(df_tasks), 50), collapse=", ")))
    if (isTRUE(TESTING_MODE)) cat(sprintf("[DEBUG VALIDATION DATA] Has PPI: %s, Has PPI_raw: %s\n", 
                "PPI" %in% names(df_tasks), "PPI_raw" %in% names(df_tasks)))
    if ("PPI" %in% names(df_tasks) && isTRUE(TESTING_MODE)) {
      cat(sprintf("[DEBUG VALIDATION DATA] PPI summary: min=%.4f, median=%.4f, max=%.4f, NA=%d\n",
                  min(df_tasks$PPI, na.rm=TRUE), median(df_tasks$PPI, na.rm=TRUE), 
                  max(df_tasks$PPI, na.rm=TRUE), sum(is.na(df_tasks$PPI))))
    }
  }

  # Group by LOCATION (not location-year)
  target_locations <- unique(df_tasks$location_id)

  n_locs_to_process <- length(target_locations)
  
  # Skip processing if no locations to process (training data processing disabled)
  if (n_locs_to_process == 0) {
    cat("[INFO] No training locations to process - skipping main processing loop\n")
    training_results_list <- list()
  } else {
  # Use centralized BATCH_SIZE setting from the top-level CONFIG (to change batch size, edit the USER-TUNABLE PARAMETERS block)
  loc_batches <- split(target_locations, ceiling(seq_along(target_locations) / BATCH_SIZE))
  n_batches <- length(loc_batches)
  pb_width <- min(40L, max(4L, n_batches))

  if (!isTRUE(QUIET_MODE)) cat(sprintf("Processing %d locations in %d batches (approx %d locations/batch)...\n",
              n_locs_to_process, length(loc_batches), BATCH_SIZE))

  # Collect full per-task results for downstream reporting (variant trajectories, diagnostics, uncertainty)
  training_results_list <- list()

  start_time <- Sys.time()

  for (i in seq_along(loc_batches)) {
    batch_locs <- loc_batches[[i]]
    batch_df <- df_tasks[df_tasks$location_id %in% batch_locs, ]
    batch_location_list <- split(batch_df, batch_df$location_id)
    # Suppress any verbose output from per-location processing
    sink_file <- tempfile()
    con_out <- tryCatch(file(sink_file, open = "wt"), error = function(e) NULL)
    con_msg <- tryCatch(file(paste0(sink_file, ".msg"), open = "wt"), error = function(e) NULL)
    tryCatch({
      if (!is.null(con_out) && inherits(con_out, "connection") && isOpen(con_out)) {
        sink(con_out, type = "output")
      }
      if (!is.null(con_msg) && inherits(con_msg, "connection") && isOpen(con_msg)) {
        sink(con_msg, type = "message")
      }
      batch_results <- .run_map(batch_location_list, fit_one_location, show_pb = FALSE)
    }, finally = {
      try(sink(type = "message"), silent = TRUE)
      try(sink(type = "output"), silent = TRUE)
      if (!is.null(con_out)) try(close(con_out), silent = TRUE)
      if (!is.null(con_msg)) try(close(con_msg), silent = TRUE)
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
                    list(coef_df = r_in$coef_df, diagnostics = r_in$diagnostics)
                  })(r)
      }
    }

  if (!isTRUE(QUIET_MODE)) cat(sprintf("\r  [Batch %d/%d complete]  ", i, n_batches))
  }
  cat("\n")
  
  end_time <- Sys.time()
  processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(
    "\nTraining processing finished in %.2f seconds (%.2f minutes)\n",
    processing_time, processing_time / 60
  ))
  } # End of n_locs_to_process > 0 conditional

  # Show diagnostic information from debug log
  debug_log_file <- file.path(tempdir(), "fit_one_location_debug.log")
  if (file.exists(debug_log_file) && !isTRUE(QUIET_MODE)) {
    cat(sprintf("[DEBUG LOG] fit_one_location debug output: %s\n", debug_log_file))
    log_lines <- readLines(debug_log_file, warn = FALSE)
    if (length(log_lines) > 0) {
      # Show summary statistics
      error_lines <- grep("ERROR", log_lines, value = TRUE)
      null_lines <- grep("returned NULL", log_lines, value = TRUE)
      
      cat(sprintf("[DEBUG LOG] Total log lines: %d\n", length(log_lines)))
      cat(sprintf("[DEBUG LOG] Error messages: %d\n", length(error_lines)))
      cat(sprintf("[DEBUG LOG] NULL returns: %d\n", length(null_lines)))
      
      # Show last 30 lines for context
      cat(sprintf("\n[DEBUG LOG] Last 30 lines:\n"))
      tail_lines <- tail(log_lines, 30)
      for (line in tail_lines) cat("  ", line, "\n")
    }
  }

  # ==========================================================================
  # EXTRACT AND VALIDATE RESULTS (use validation_coefs if available)
  # ==========================================================================
  
  cat("\n=== USING VALIDATION RESULTS FOR CONFUSION MATRIX ===\n")
  
  # validation_coefs already created from validation processing above
  # No need to extract from training_results_list since we skipped training processing
  
  if (is.null(validation_coefs) || nrow(validation_coefs) == 0) {
    cat("[WARNING] No validation coefficients were produced!\n")
    validation_coefs <- data.frame()  # Empty for safety
  } else {
    cat(sprintf("[VALIDATION] Using %d coefficient rows from %d locations for confusion matrix\n",
                nrow(validation_coefs), length(unique(validation_coefs$location_id))))
    
    # Show class distribution in results
    if ("Veg" %in% names(validation_coefs)) {
      result_classes <- table(tolower(validation_coefs$Veg))
      cat("[VALIDATION] Predicted class distribution:\n")
      print(result_classes)
    }
  }
  
  # ==========================================================================
  # INFERENCE PROCESSING (separate from validation)
  # ==========================================================================
  cat("\n=== STARTING INFERENCE PROCESSING ===\n")
  
  # Generate variant similarity heatmap before any inference processing
  cat("[INFERENCE] Generating variant similarity heatmap...\n")
  ensure_variant_similarity_heatmap(force = TRUE)
  
  # Load inference data from INFERENCE_CSV if not already loaded
  if (!isTRUE(SKIP_INFERENCE) && !isTRUE(TESTING_MODE)) {
    cat("[INFERENCE] Loading inference data from separate CSV file...\n")
    load_and_prepare_inference_data()
  }

  # Respect SKIP_INFERENCE: bypass inference loop entirely while still allowing
  # validation reporting and downstream outputs to proceed.
  if (isTRUE(SKIP_INFERENCE)) {
    cat("[INFERENCE] SKIP_INFERENCE=TRUE -> skipping inference processing.\n")
    inference_coefs <- NULL
  } else if (exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0) {
    # Prepare inference task list - use df_tasks_inference loaded from INFERENCE_CSV
    df_tasks_inference_proc <- df_tasks_inference
    
    if ("Veg" %in% names(df_tasks_inference_proc)) {
      df_tasks_inference_proc$Veg <- tolower(as.character(df_tasks_inference_proc$Veg))
    }
    
    if ("date" %in% names(df_tasks_inference_proc)) {
      df_tasks_inference_proc$date <- as.Date(df_tasks_inference_proc$date)
      if (!"pheno_year" %in% names(df_tasks_inference_proc)) {
        df_tasks_inference_proc$pheno_year <- assign_pheno_year(df_tasks_inference_proc$date)
      }
      df_tasks_inference_proc$doy <- pheno_doy(df_tasks_inference_proc$date)
      df_tasks_inference_proc$doy[df_tasks_inference_proc$doy < 1 | df_tasks_inference_proc$doy > 366] <- NA_integer_
    }
    
    # Get unique inference locations
    inference_locations <- df_tasks_inference_proc |>
      dplyr::filter(!is.na(.data$location_id) & trimws(.data$location_id) != "") |>
      dplyr::distinct(.data$location_id)
    inference_locations$location_id <- trimws(as.character(inference_locations$location_id))
    
    # Apply MAX_INFERENCE_LOCATIONS limit if set
    if (exists("MAX_INFERENCE_LOCATIONS") && nrow(inference_locations) > MAX_INFERENCE_LOCATIONS) {
      set.seed(123) # deterministic sampling for reproducibility
      sampled_loc_ids <- sample(inference_locations$location_id, MAX_INFERENCE_LOCATIONS, replace = FALSE)
      inference_locations <- inference_locations[inference_locations$location_id %in% sampled_loc_ids, , drop = FALSE]
      df_tasks_inference_proc <- df_tasks_inference_proc[df_tasks_inference_proc$location_id %in% sampled_loc_ids, ]
      cat(sprintf("[INFERENCE] Reduced to %d locations (random sample due to MAX_INFERENCE_LOCATIONS=%d)\n",
                  nrow(inference_locations), MAX_INFERENCE_LOCATIONS))
    }

    # Build inference_loc_years with defensive checks
    required_cols <- c("location_id", "pheno_year")
    if (!all(required_cols %in% names(df_tasks_inference_proc))) {
      missing_cols <- setdiff(required_cols, names(df_tasks_inference_proc))
      stop(sprintf("[INFERENCE ERROR] Missing required columns in df_tasks_inference_proc: %s", paste(missing_cols, collapse = ", ")))
    }

    tmp_inf <- df_tasks_inference_proc
    tmp_inf$location_id <- as.character(tmp_inf$location_id)
    tmp_inf$pheno_year <- as.integer(tmp_inf$pheno_year)
    tmp_inf$location_id <- trimws(tmp_inf$location_id)

    valid_mask <- !is.na(tmp_inf$location_id) & tmp_inf$location_id != "" & !is.na(tmp_inf$pheno_year) & tmp_inf$pheno_year > 0
    if (any(valid_mask)) {
      inference_loc_years <- unique(tmp_inf[valid_mask, c("location_id", "pheno_year")])
      inference_loc_years <- data.frame(location_id = inference_loc_years[,1], pheno_year = inference_loc_years[,2], stringsAsFactors = FALSE)
    } else {
      inference_loc_years <- data.frame(location_id = character(0), pheno_year = integer(0), stringsAsFactors = FALSE)
    }

    cat(sprintf("[INFERENCE] Processing %d locations (%d location-year pairs)\n",
                nrow(inference_locations), nrow(inference_loc_years)))
    
    # Process inference in batches
    inference_location_list <- as.character(inference_locations$location_id)
    inference_loc_batches <- split(inference_location_list, 
                                    ceiling(seq_along(inference_location_list) / BATCH_SIZE))
    
    cat(sprintf("[INFERENCE] Processing %d locations in %d batches (approx %d locations/batch)...\n",
                length(inference_location_list), length(inference_loc_batches), BATCH_SIZE))
    
    inference_results_list <- list()
    inference_start_time <- Sys.time()
    
    for (i in seq_along(inference_loc_batches)) {
      batch_locs <- inference_loc_batches[[i]]
      batch_df <- df_tasks_inference_proc[df_tasks_inference_proc$location_id %in% batch_locs, ]
      batch_location_list <- split(batch_df, batch_df$location_id)
      
      # Diagnostic: check batch data structure
      if (isTRUE(TESTING_MODE)) {
        cat(sprintf("  [Batch %d] Processing %d locations, batch_df has %d rows\n", 
                    i, length(batch_location_list), nrow(batch_df)))
        if (length(batch_location_list) > 0) {
          first_loc_name <- names(batch_location_list)[1]
          first_loc_data <- batch_location_list[[1]]
          cat(sprintf("  [Batch %d] First location '%s' has %d rows\n",
                      i, first_loc_name, nrow(first_loc_data)))
          cat(sprintf("  [Batch %d] All columns (%d): %s\n",
                      i, length(names(first_loc_data)), paste(names(first_loc_data), collapse=", ")))
          if ("pheno_year" %in% names(first_loc_data)) {
            years_sample <- head(sort(unique(first_loc_data$pheno_year)), 5)
            cat(sprintf("  [Batch %d] First location pheno_years (first 5): %s\n",
                        i, paste(years_sample, collapse=", ")))
          } else {
            cat(sprintf("  [Batch %d] ERROR: 'pheno_year' column missing!\n", i))
          }
          # Check for required index columns
          required_indices <- c("NDVI", "EVI", "PPI", "DVI")
          missing_indices <- setdiff(required_indices, names(first_loc_data))
          if (length(missing_indices) > 0) {
            cat(sprintf("  [Batch %d] WARNING: Missing required indices: %s\n",
                        i, paste(missing_indices, collapse=", ")))
          }
        }
        # Save debug log location before sinking
        debug_log <- file.path(tempdir(), "fit_one_location_debug.log")
        cat(sprintf("  [Batch %d] Debug log: %s\n", i, debug_log))
      } else {
        # Not in testing mode: still set the debug log path (used for later conditional printing),
        # but do not print verbose per-batch diagnostics to stdout.
        debug_log <- file.path(tempdir(), "fit_one_location_debug.log")
      }

      # Suppress any verbose output from per-location processing
      sink_file_inf <- tempfile()
      con_out_inf <- tryCatch(file(sink_file_inf, open = "wt"), error = function(e) NULL)
      con_msg_inf <- tryCatch(file(paste0(sink_file_inf, ".msg"), open = "wt"), error = function(e) NULL)
      tryCatch({
        if (!is.null(con_out_inf) && inherits(con_out_inf, "connection") && isOpen(con_out_inf)) {
          sink(con_out_inf, type = "output")
        }
        if (!is.null(con_msg_inf) && inherits(con_msg_inf, "connection") && isOpen(con_msg_inf)) {
          sink(con_msg_inf, type = "message")
        }
        batch_results <- .run_map(batch_location_list, fit_one_location, show_pb = FALSE)
      }, finally = {
        try(sink(type = "message"), silent = TRUE)
        try(sink(type = "output"), silent = TRUE)
        if (!is.null(con_out_inf)) try(close(con_out_inf), silent = TRUE)
        if (!is.null(con_msg_inf)) try(close(con_msg_inf), silent = TRUE)
        try(unlink(sink_file_inf), silent = TRUE)
        try(unlink(paste0(sink_file_inf, ".msg")), silent = TRUE)
      })
      
      # Store results
      n_null_results <- 0
      n_empty_years <- 0
      n_empty_coefs <- 0
      n_stored <- 0
      
      for (nm in names(batch_results)) {
        loc_result <- batch_results[[nm]]
        if (is.null(loc_result)) {
          n_null_results <- n_null_results + 1
          next
        }
        
        # fit_one_location returns a list of year results
        for (yr_char in names(loc_result)) {
          r <- loc_result[[yr_char]]
          if (is.null(r)) {
            n_empty_years <- n_empty_years + 1
            next
          }
          if (is.null(r$coef_df) || nrow(r$coef_df) == 0) {
            n_empty_coefs <- n_empty_coefs + 1
            next
          }
          
          res_key <- if ('pheno_year' %in% names(r$coef_df)) {
            paste(nm, unique(r$coef_df$pheno_year), sep = "_")
          } else {
            paste(nm, yr_char, sep = "_")
          }
          inference_results_list[[res_key]] <- list(coef_df = r$coef_df, diagnostics = r$diagnostics)
          n_stored <- n_stored + 1
        }
      }
      
      if (isTRUE(TESTING_MODE)) {
        cat(sprintf("  [Batch %d] Results: %d stored, %d null locations, %d null years, %d empty coefs\n",
                    i, n_stored, n_null_results, n_empty_years, n_empty_coefs))

        # Show debug log content if all locations failed
        if (n_null_results == length(batch_results) && file.exists(debug_log)) {
          cat(sprintf("  [Batch %d] All locations returned NULL - showing last 20 lines of debug log:\n", i))
          log_lines <- tryCatch(readLines(debug_log), error = function(e) character(0))
          if (length(log_lines) > 0) {
            tail_lines <- tail(log_lines, 20)
            for (line in tail_lines) cat("    ", line, "\n")
          }
        }

        cat(sprintf("\r  [Inference Batch %d/%d complete]  ", i, length(inference_loc_batches)))
      } else {
        # Quiet/production mode: emit a simple progress indicator without verbose batch diagnostics
        if (!isTRUE(QUIET_MODE)) cat(sprintf("\r  [Inference Batch %d/%d complete]  ", i, length(inference_loc_batches)))
      }
    }
    cat("\n")
    
    inference_end_time <- Sys.time()
    inference_processing_time <- as.numeric(difftime(inference_end_time, inference_start_time, units = "secs"))
    cat(sprintf("[INFERENCE] Processing finished in %.2f seconds (%.2f minutes)\n",
                inference_processing_time, inference_processing_time / 60))
    
    # Combine inference results
    inference_coefs <- do.call(rbind, lapply(inference_results_list, function(r) {
      if (!is.null(r$coef_df)) r$coef_df else NULL
    }))

    # Remove rows with NA or empty Veg values
    if (!is.null(inference_coefs) && nrow(inference_coefs) > 0) {
      valid_veg <- !is.na(inference_coefs$Veg) &
                   nchar(trimws(as.character(inference_coefs$Veg))) > 0 &
                   trimws(as.character(inference_coefs$Veg)) != "NA"
      na_veg_count <- sum(!valid_veg)
      if (na_veg_count > 0) {
        cat(sprintf("[INFERENCE] Removing %d rows with NA/empty Veg values\n", na_veg_count))
        inference_coefs <- inference_coefs[valid_veg, ]
      }
    }

    cat(sprintf("[INFERENCE] Processed %d inference coefficient rows\n",
                if(!is.null(inference_coefs)) nrow(inference_coefs) else 0))

    # --- INFERENCE: PPI-based barren estimation (always enabled) ---
    if (!is.null(inference_coefs) && nrow(inference_coefs) > 0 && exists("df_tasks_inference_proc") && ("PPI" %in% names(df_tasks_inference_proc) || "PPI_raw" %in% names(df_tasks_inference_proc))) {
      cat("[INFERENCE] Running PPI-based barren estimation (inference mode)\n")
      ppi_inf_full <- location_bootstrap_ppi(inference_coefs, df_tasks_inference_proc, B = BOOTSTRAP_B, seed = 123)
      if (is.null(ppi_inf_full) || nrow(ppi_inf_full) == 0) {
        cat("[INFERENCE] PPI inference aggregation returned no results (no matching loc-year PPI values).\n")
      } else {
        ppi_inf_veg <- ppi_inf_full[!tolower(trimws(ppi_inf_full$Veg)) %in% c("barren"), ]
        ppi_inf_barren <- ppi_inf_full[tolower(trimws(ppi_inf_full$Veg)) %in% c("barren"), ]

        if (!is.null(ppi_inf_veg) && nrow(ppi_inf_veg) > 0) {
          p_inf_ppi_ts <- ggplot(ppi_inf_veg, aes(x = year, y = global_coef, color = Veg, group = Veg)) +
            add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
            geom_line(linewidth = 1) +
            geom_point(show.legend = FALSE) +
            geom_ribbon(aes(ymin = coef_025, ymax = coef_975, fill = Veg), alpha = 0.15, color = NA) +
            labs(title = "PPI-Normalized Vegetation Fractions",
                 x = "Year", y = "Total Normalized Fraction", color = "Veg", fill = "Veg") +
            theme_minimal()
          ggsave(file.path(OUT_DIR, "inference_ppi_normalized_timeseries.png"), p_inf_ppi_ts, width = 8, height = 6)
          cat(sprintf("Saved inference PPI-normalized time series plot to: %s\n", file.path(OUT_DIR, "inference_ppi_normalized_timeseries.png")))

        }

        if (!is.null(ppi_inf_barren) && nrow(ppi_inf_barren) > 0) {
          p_inf_barren <- ggplot(ppi_inf_barren, aes(x = year, y = global_coef)) +
            add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
            geom_line(color = "saddlebrown", linewidth = 1) +
            geom_point(color = "saddlebrown") +
            geom_ribbon(aes(ymin = coef_025, ymax = coef_975), alpha = 0.15, fill = "saddlebrown", color = NA) +
            labs(title = "PPI-Based Barren Fraction", x = "Year", y = "Barren Fraction") +
            scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
            theme_minimal()
          ggsave(file.path(OUT_DIR, "inference_ppi_barren_cover.png"), p_inf_barren, width = 8, height = 6)
          cat(sprintf("Saved inference PPI-based barren cover plot to: %s\n", file.path(OUT_DIR, "inference_ppi_barren_cover.png")))


          # Herbs vs Woody (inference PPI)
          herbs_woody_inf <- ppi_inf_full[tolower(trimws(ppi_inf_full$Veg)) %in% c("herbs", "woody"), ]
          if (nrow(herbs_woody_inf) > 0) {
            herbs_woody_inf$Veg <- ifelse(tolower(herbs_woody_inf$Veg) == "herbs", "Herbs", "Woody")
            p_inf_ppi_herbs_woody <- ggplot(herbs_woody_inf, aes(x = year, y = global_coef, color = Veg, group = Veg)) +
              add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
              geom_line(linewidth = 1) +
              geom_point(show.legend = FALSE) +
              geom_ribbon(aes(ymin = coef_025, ymax = coef_975, fill = Veg), alpha = 0.15, color = NA) +
              scale_color_manual(values = c("Herbs" = "#2E8B57", "Woody" = "#8B4513")) +
              scale_fill_manual(values = c("Herbs" = "#2E8B57", "Woody" = "#8B4513")) +
              labs(title = "Inference PPI: Herbs vs Woody Vegetation",
                   x = "Year", y = "Total Normalized Fraction", color = "Type", fill = "Type") +
              theme_minimal()
            ggsave(file.path(OUT_DIR, "inference_ppi_herbs_vs_woody.png"), p_inf_ppi_herbs_woody, width = 8, height = 6)
            cat(sprintf("Saved inference PPI herbs vs woody plot to: %s\n", file.path(OUT_DIR, "inference_ppi_herbs_vs_woody.png")))


            # --- Stacked area (Proportion) and Woody/Herbs ratio (Inference PPI) ---
            df_wide <- tryCatch({
              tidyr::pivot_wider(herbs_woody_inf |> dplyr::select(year, Veg, global_coef), names_from = Veg, values_from = global_coef)
            }, error = function(e) NULL)
            if (!is.null(df_wide) && all(c("Herbs","Woody") %in% names(df_wide))) {
              df_wide$Herbs <- as.numeric(df_wide$Herbs)
              df_wide$Woody <- as.numeric(df_wide$Woody)
              df_wide$total <- rowSums(df_wide[, c("Herbs","Woody")], na.rm = TRUE)

              df_prop <- df_wide |> dplyr::filter(is.finite(total) & total > 0) |> dplyr::mutate(Herbs = Herbs/total, Woody = Woody/total) |> tidyr::pivot_longer(cols = c("Herbs","Woody"), names_to = "Veg", values_to = "prop")
              p_inf_ppi_stacked <- ggplot(df_prop, aes(x = year, y = prop, fill = Veg)) +
                geom_area() +
                scale_fill_manual(values = c("Herbs" = "#2E8B57", "Woody" = "#8B4513")) +
                labs(title = "Inference PPI: Herbs vs Woody (Proportion, stacked)", x = "Year", y = "Proportion", fill = "Type") +
                theme_minimal()
              ggsave(file.path(OUT_DIR, "inference_ppi_herbs_vs_woody_stacked.png"), p_inf_ppi_stacked, width = 8, height = 6)
              cat(sprintf("Saved inference PPI herbs vs woody stacked plot to: %s\n", file.path(OUT_DIR, "inference_ppi_herbs_vs_woody_stacked.png")))

              df_ratio <- df_wide |> dplyr::mutate(ratio = ifelse(is.finite(Herbs) & Herbs > 0, Woody / Herbs, NA_real_))
              p_inf_ppi_ratio <- ggplot(df_ratio, aes(x = year, y = ratio)) +
                geom_line(color = "#8B4513", linewidth = 1) +
                geom_point() +
                labs(title = "Inference PPI: Woody / Herbs Ratio", x = "Year", y = "Woody / Herbs") +
                theme_minimal()
              ggsave(file.path(OUT_DIR, "inference_ppi_woody_over_herbs.png"), p_inf_ppi_ratio, width = 8, height = 6)
              cat(sprintf("Saved inference PPI woody/herbs ratio plot to: %s\n", file.path(OUT_DIR, "inference_ppi_woody_over_herbs.png")))
            } else {
              cat("Cannot create inference PPI stacked/ratio plot: missing Herbs/Woody rows.\n")
            }

            # --- Tamarix vs Populus vs Woody_unknown stacked plot (Inference PPI) ---
            woody_types_inf <- ppi_inf_full[tolower(trimws(ppi_inf_full$Veg)) %in% c("tamarix", "populus", "woody_unknown"), ]
            if (nrow(woody_types_inf) > 0) {
              woody_types_inf$Veg <- dplyr::case_when(
                tolower(woody_types_inf$Veg) == "tamarix" ~ "Tamarix",
                tolower(woody_types_inf$Veg) == "populus" ~ "Populus",
                tolower(woody_types_inf$Veg) == "woody_unknown" ~ "Woody Unknown",
                TRUE ~ woody_types_inf$Veg
              )

              df_wide_woody <- tryCatch({
                tidyr::pivot_wider(woody_types_inf |> dplyr::select(year, Veg, global_coef), names_from = Veg, values_from = global_coef)
              }, error = function(e) NULL)

              if (!is.null(df_wide_woody)) {
                woody_cols <- intersect(c("Tamarix", "Populus", "Woody Unknown"), names(df_wide_woody))
                if (length(woody_cols) >= 2) {
                  for (col in woody_cols) {
                    df_wide_woody[[col]] <- as.numeric(df_wide_woody[[col]])
                    df_wide_woody[[col]][is.na(df_wide_woody[[col]])] <- 0
                  }
                  df_wide_woody$total <- rowSums(df_wide_woody[, woody_cols, drop = FALSE], na.rm = TRUE)

                  df_prop_woody <- df_wide_woody |>
                    dplyr::filter(is.finite(total) & total > 0) |>
                    dplyr::mutate(across(all_of(woody_cols), ~ . / total)) |>
                    tidyr::pivot_longer(cols = all_of(woody_cols), names_to = "Veg", values_to = "prop")

                  p_inf_ppi_woody_stacked <- ggplot(df_prop_woody, aes(x = year, y = prop, fill = Veg)) +
                    geom_area() +
                    scale_fill_manual(values = c("Tamarix" = "#CD853F", "Populus" = "#228B22", "Woody Unknown" = "#808080")) +
                    labs(title = "Inference PPI: Tamarix vs Populus vs Woody Unknown (Proportion, stacked)",
                         x = "Year", y = "Proportion", fill = "Type") +
                    theme_minimal()
                  ggsave(file.path(OUT_DIR, "inference_ppi_tamarix_populus_woody_unknown_stacked.png"), p_inf_ppi_woody_stacked, width = 8, height = 6)
                  cat(sprintf("Saved inference PPI Tamarix vs Populus vs Woody Unknown stacked plot to: %s\n", file.path(OUT_DIR, "inference_ppi_tamarix_populus_woody_unknown_stacked.png")))
                } else {
                  cat("Cannot create inference PPI Tamarix/Populus/Woody_unknown stacked plot: need at least 2 of these veg types.\n")
                }
              } else {
                cat("Cannot create inference PPI Tamarix/Populus/Woody_unknown stacked plot: pivot failed.\n")
              }
            } else {
              cat("No Tamarix/Populus/Woody_unknown data available for inference PPI woody types stacked plot.\n")
            }
          } else {
            cat("No herbs/woody data available for inference PPI herbs vs woody plot.\n")
          }
        } else {
          cat("[INFERENCE] No barren rows found in PPI inference aggregation results.\n")
        }
      }
    } else {
      cat("[INFERENCE] PPI data not available for inference; skipping PPI barren estimation.\n")
    }

    # --- INFERENCE: MSAVI-based aggregation ---
    if (!is.null(inference_coefs) && nrow(inference_coefs) > 0 && exists("df_tasks_inference_proc") && ("MSAVI" %in% names(df_tasks_inference_proc) || "MSAVI_raw" %in% names(df_tasks_inference_proc))) {
      cat("[INFERENCE] Running MSAVI-based aggregation (inference mode)\n")
      msavi_inf_full <- tryCatch({ location_bootstrap_msavi(inference_coefs, df_tasks_inference_proc, B = BOOTSTRAP_B, seed = 123) }, error = function(e) { cat(sprintf("[MSAVI INFERENCE] failed: %s\n", e$message)); NULL })
      if (!is.null(msavi_inf_full) && nrow(msavi_inf_full) > 0) {
        msavi_inf_veg <- msavi_inf_full[!tolower(trimws(msavi_inf_full$Veg)) %in% c("barren"), ]
        msavi_inf_barren <- msavi_inf_full[tolower(trimws(msavi_inf_full$Veg)) %in% c("barren"), ]

        if (!is.null(msavi_inf_veg) && nrow(msavi_inf_veg) > 0) {
          p_inf_msavi_ts <- ggplot(msavi_inf_veg, aes(x = year, y = global_coef, color = Veg, group = Veg)) +
            add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
            geom_line(linewidth = 1) +
            geom_point(show.legend = FALSE) +
            geom_ribbon(aes(ymin = coef_025, ymax = coef_975, fill = Veg), alpha = 0.15, color = NA) +
            labs(title = "MSAVI-Normalized Vegetation Fractions",
                 x = "Year", y = "Total Normalized Fraction", color = "Veg", fill = "Veg") +
            theme_minimal()
          ggsave(file.path(OUT_DIR, "inference_msavi_normalized_timeseries.png"), p_inf_msavi_ts, width = 8, height = 6)
          cat(sprintf("Saved inference MSAVI-normalized time series plot to: %s\n", file.path(OUT_DIR, "inference_msavi_normalized_timeseries.png")))

        }

        if (!is.null(msavi_inf_barren) && nrow(msavi_inf_barren) > 0) {
          p_inf_msavi_barren <- ggplot(msavi_inf_barren, aes(x = year, y = global_coef)) +
            add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
            geom_line(color = "saddlebrown", linewidth = 1) +
            geom_point(color = "saddlebrown") +
            geom_ribbon(aes(ymin = coef_025, ymax = coef_975), alpha = 0.15, fill = "saddlebrown", color = NA) +
            labs(title = "MSAVI-Based Barren Fraction", x = "Year", y = "Barren Fraction") +
            scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
            theme_minimal()
          ggsave(file.path(OUT_DIR, "inference_msavi_barren_cover.png"), p_inf_msavi_barren, width = 8, height = 6)
          cat(sprintf("Saved inference MSAVI-based barren cover plot to: %s\n", file.path(OUT_DIR, "inference_msavi_barren_cover.png")))


          # Herbs vs Woody (inference MSAVI)
          herbs_woody_msavi_inf <- msavi_inf_full[tolower(trimws(msavi_inf_full$Veg)) %in% c("herbs", "woody"), ]
          if (nrow(herbs_woody_msavi_inf) > 0) {
            herbs_woody_msavi_inf$Veg <- ifelse(tolower(herbs_woody_msavi_inf$Veg) == "herbs", "Herbs", "Woody")
            p_inf_msavi_herbs_woody <- ggplot(herbs_woody_msavi_inf, aes(x = year, y = global_coef, color = Veg, group = Veg)) +
              add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
              geom_line(linewidth = 1) +
              geom_point(show.legend = FALSE) +
              geom_ribbon(aes(ymin = coef_025, ymax = coef_975, fill = Veg), alpha = 0.15, color = NA) +
              scale_color_manual(values = c("Herbs" = "#2E8B57", "Woody" = "#8B4513")) +
              scale_fill_manual(values = c("Herbs" = "#2E8B57", "Woody" = "#8B4513")) +
              labs(title = "Inference MSAVI: Herbs vs Woody Vegetation",
                   x = "Year", y = "Total Normalized Fraction", color = "Type", fill = "Type") +
              theme_minimal()
            ggsave(file.path(OUT_DIR, "inference_msavi_herbs_vs_woody.png"), p_inf_msavi_herbs_woody, width = 8, height = 6)
            cat(sprintf("Saved inference MSAVI herbs vs woody plot to: %s\n", file.path(OUT_DIR, "inference_msavi_herbs_vs_woody.png")))


            # --- Stacked area (Proportion) and Woody/Herbs ratio (Inference MSAVI) ---
            df_wide <- tryCatch({
              tidyr::pivot_wider(herbs_woody_msavi_inf |> dplyr::select(year, Veg, global_coef), names_from = Veg, values_from = global_coef)
            }, error = function(e) NULL)
            if (!is.null(df_wide) && all(c("Herbs","Woody") %in% names(df_wide))) {
              df_wide$Herbs <- as.numeric(df_wide$Herbs)
              df_wide$Woody <- as.numeric(df_wide$Woody)
              df_wide$total <- rowSums(df_wide[, c("Herbs","Woody")], na.rm = TRUE)

              df_prop <- df_wide |> dplyr::filter(is.finite(total) & total > 0) |> dplyr::mutate(Herbs = Herbs/total, Woody = Woody/total) |> tidyr::pivot_longer(cols = c("Herbs","Woody"), names_to = "Veg", values_to = "prop")
              p_inf_msavi_stacked <- ggplot(df_prop, aes(x = year, y = prop, fill = Veg)) +
                geom_area() +
                scale_fill_manual(values = c("Herbs" = "#2E8B57", "Woody" = "#8B4513")) +
                labs(title = "Inference MSAVI: Herbs vs Woody (Proportion, stacked)", x = "Year", y = "Proportion", fill = "Type") +
                theme_minimal()
              ggsave(file.path(OUT_DIR, "inference_msavi_herbs_vs_woody_stacked.png"), p_inf_msavi_stacked, width = 8, height = 6)
              cat(sprintf("Saved inference MSAVI herbs vs woody stacked plot to: %s\n", file.path(OUT_DIR, "inference_msavi_herbs_vs_woody_stacked.png")))

              df_ratio <- df_wide |> dplyr::mutate(ratio = ifelse(is.finite(Herbs) & Herbs > 0, Woody / Herbs, NA_real_))
              p_inf_msavi_ratio <- ggplot(df_ratio, aes(x = year, y = ratio)) +
                geom_line(color = "#8B4513", linewidth = 1) +
                geom_point() +
                labs(title = "Inference MSAVI: Woody / Herbs Ratio", x = "Year", y = "Woody / Herbs") +
                theme_minimal()
              ggsave(file.path(OUT_DIR, "inference_msavi_woody_over_herbs.png"), p_inf_msavi_ratio, width = 8, height = 6)
              cat(sprintf("Saved inference MSAVI woody/herbs ratio plot to: %s\n", file.path(OUT_DIR, "inference_msavi_woody_over_herbs.png")))
            } else {
              cat("Cannot create inference MSAVI stacked/ratio plot: missing Herbs/Woody rows.\n")
            }
          } else {
            cat("No herbs/woody data available for inference MSAVI herbs vs woody plot.\n")
          }
        } else {
          cat("[INFERENCE] No barren rows found in MSAVI inference aggregation results.\n")
        }
      } else {
        cat("[INFERENCE] MSAVI inference aggregation returned no results.\n")
      }
    } else {
      cat("[INFERENCE] MSAVI data not available for inference; skipping MSAVI aggregation.\n")
    }

    # --- INFERENCE: NDVI-based aggregation ---
    if (!is.null(inference_coefs) && nrow(inference_coefs) > 0 && exists("df_tasks_inference_proc") && ("NDVI" %in% names(df_tasks_inference_proc) || "NDVI_raw" %in% names(df_tasks_inference_proc))) {
      cat("[INFERENCE] Running NDVI-based aggregation (inference mode)\n")
      ndvi_inf_full <- tryCatch({ location_bootstrap_ndvi(inference_coefs, df_tasks_inference_proc, B = BOOTSTRAP_B, seed = 123) }, error = function(e) { cat(sprintf("[NDVI INFERENCE] failed: %s\n", e$message)); NULL })
      if (!is.null(ndvi_inf_full) && nrow(ndvi_inf_full) > 0) {
        ndvi_inf_veg <- ndvi_inf_full[!tolower(trimws(ndvi_inf_full$Veg)) %in% c("barren"), ]
        ndvi_inf_barren <- ndvi_inf_full[tolower(trimws(ndvi_inf_full$Veg)) %in% c("barren"), ]

        if (!is.null(ndvi_inf_veg) && nrow(ndvi_inf_veg) > 0) {
          p_inf_ndvi_ts <- ggplot(ndvi_inf_veg, aes(x = year, y = global_coef, color = Veg, group = Veg)) +
            add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
            geom_line(linewidth = 1) +
            geom_point(show.legend = FALSE) +
            geom_ribbon(aes(ymin = coef_025, ymax = coef_975, fill = Veg), alpha = 0.15, color = NA) +
            labs(title = "NDVI-Normalized Vegetation Fractions",
                 x = "Year", y = "Total Normalized Fraction", color = "Veg", fill = "Veg") +
            theme_minimal()
          ggsave(file.path(OUT_DIR, "inference_ndvi_normalized_timeseries.png"), p_inf_ndvi_ts, width = 8, height = 6)
          cat(sprintf("Saved inference NDVI-normalized time series plot to: %s\n", file.path(OUT_DIR, "inference_ndvi_normalized_timeseries.png")))

        }

        if (!is.null(ndvi_inf_barren) && nrow(ndvi_inf_barren) > 0) {
          p_inf_ndvi_barren <- ggplot(ndvi_inf_barren, aes(x = year, y = global_coef)) +
            add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
            geom_line(color = "saddlebrown", linewidth = 1) +
            geom_point(color = "saddlebrown") +
            geom_ribbon(aes(ymin = coef_025, ymax = coef_975), alpha = 0.15, fill = "saddlebrown", color = NA) +
            labs(title = "NDVI-Based Barren Fraction", x = "Year", y = "Barren Fraction") +
            scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
            theme_minimal()
          ggsave(file.path(OUT_DIR, "inference_ndvi_barren_cover.png"), p_inf_ndvi_barren, width = 8, height = 6)
          cat(sprintf("Saved inference NDVI-based barren cover plot to: %s\n", file.path(OUT_DIR, "inference_ndvi_barren_cover.png")))


          # Herbs vs Woody (inference NDVI)
          herbs_woody_ndvi_inf <- ndvi_inf_full[tolower(trimws(ndvi_inf_full$Veg)) %in% c("herbs", "woody"), ]
          if (nrow(herbs_woody_ndvi_inf) > 0) {
            herbs_woody_ndvi_inf$Veg <- ifelse(tolower(herbs_woody_ndvi_inf$Veg) == "herbs", "Herbs", "Woody")
            p_inf_ndvi_herbs_woody <- ggplot(herbs_woody_ndvi_inf, aes(x = year, y = global_coef, color = Veg, group = Veg)) +
              add_excluded_years_shade(is_date = FALSE) + add_year_lines(is_date = FALSE) +
              geom_line(linewidth = 1) +
              geom_point(show.legend = FALSE) +
              geom_ribbon(aes(ymin = coef_025, ymax = coef_975, fill = Veg), alpha = 0.15, color = NA) +
              scale_color_manual(values = c("Herbs" = "#2E8B57", "Woody" = "#8B4513")) +
              scale_fill_manual(values = c("Herbs" = "#2E8B57", "Woody" = "#8B4513")) +
              labs(title = "Inference NDVI: Herbs vs Woody Vegetation",
                   x = "Year", y = "Total Normalized Fraction", color = "Type", fill = "Type") +
              theme_minimal()
            ggsave(file.path(OUT_DIR, "inference_ndvi_herbs_vs_woody.png"), p_inf_ndvi_herbs_woody, width = 8, height = 6)
            cat(sprintf("Saved inference NDVI herbs vs woody plot to: %s\n", file.path(OUT_DIR, "inference_ndvi_herbs_vs_woody.png")))


            # --- Stacked area (Proportion) and Woody/Herbs ratio (Inference NDVI) ---
            df_wide <- tryCatch({
              tidyr::pivot_wider(herbs_woody_ndvi_inf |> dplyr::select(year, Veg, global_coef), names_from = Veg, values_from = global_coef)
            }, error = function(e) NULL)
            if (!is.null(df_wide) && all(c("Herbs","Woody") %in% names(df_wide))) {
              df_wide$Herbs <- as.numeric(df_wide$Herbs)
              df_wide$Woody <- as.numeric(df_wide$Woody)
              df_wide$total <- rowSums(df_wide[, c("Herbs","Woody")], na.rm = TRUE)

              df_prop <- df_wide |> dplyr::filter(is.finite(total) & total > 0) |> dplyr::mutate(Herbs = Herbs/total, Woody = Woody/total) |> tidyr::pivot_longer(cols = c("Herbs","Woody"), names_to = "Veg", values_to = "prop")
              p_inf_ndvi_stacked <- ggplot(df_prop, aes(x = year, y = prop, fill = Veg)) +
                geom_area() +
                scale_fill_manual(values = c("Herbs" = "#2E8B57", "Woody" = "#8B4513")) +
                labs(title = "Inference NDVI: Herbs vs Woody (Proportion, stacked)", x = "Year", y = "Proportion", fill = "Type") +
                theme_minimal()
              ggsave(file.path(OUT_DIR, "inference_ndvi_herbs_vs_woody_stacked.png"), p_inf_ndvi_stacked, width = 8, height = 6)
              cat(sprintf("Saved inference NDVI herbs vs woody stacked plot to: %s\n", file.path(OUT_DIR, "inference_ndvi_herbs_vs_woody_stacked.png")))

              df_ratio <- df_wide |> dplyr::mutate(ratio = ifelse(is.finite(Herbs) & Herbs > 0, Woody / Herbs, NA_real_))
              p_inf_ndvi_ratio <- ggplot(df_ratio, aes(x = year, y = ratio)) +
                geom_line(color = "#8B4513", linewidth = 1) +
                geom_point() +
                labs(title = "Inference NDVI: Woody / Herbs Ratio", x = "Year", y = "Woody / Herbs") +
                theme_minimal()
              ggsave(file.path(OUT_DIR, "inference_ndvi_woody_over_herbs.png"), p_inf_ndvi_ratio, width = 8, height = 6)
              cat(sprintf("Saved inference NDVI woody/herbs ratio plot to: %s\n", file.path(OUT_DIR, "inference_ndvi_woody_over_herbs.png")))
            } else {
              cat("Cannot create inference NDVI stacked/ratio plot: missing Herbs/Woody rows.\n")
            }
          } else {
            cat("No herbs/woody data available for inference NDVI herbs vs woody plot.\n")
          }
        } else {
          cat("[INFERENCE] No barren rows found in NDVI inference aggregation results.\n")
        }
      } else {
        cat("[INFERENCE] NDVI inference aggregation returned no results.\n")
      }
    } else {
      cat("[INFERENCE] NDVI data not available for inference; skipping NDVI aggregation.\n")
    }

    # --- INFERENCE: No-index normalization aggregation ---
    if (!is.null(inference_coefs) && nrow(inference_coefs) > 0) {
      cat("[INFERENCE] Running no-index normalization aggregation (inference mode)\n")
      noindex_inf_full <- tryCatch({ location_bootstrap_noindex(inference_coefs, df_tasks_inference_proc, B = BOOTSTRAP_B, seed = 123) }, error = function(e) { cat(sprintf("[NOINDEX INFERENCE] failed: %s\n", e$message)); NULL })
      if (!is.null(noindex_inf_full) && nrow(noindex_inf_full) > 0) {
        noindex_inf_veg <- noindex_inf_full[!tolower(trimws(noindex_inf_full$Veg)) %in% c("barren"), ]
        noindex_inf_barren <- noindex_inf_full[tolower(trimws(noindex_inf_full$Veg)) %in% c("barren"), ]

        if (!is.null(noindex_inf_veg) && nrow(noindex_inf_veg) > 0) {
          p_inf_noindex_ts <- ggplot(noindex_inf_veg, aes(x = year, y = global_coef, color = Veg, group = Veg)) +
            add_year_lines(is_date = FALSE) +
            geom_line(linewidth = 1) +
            geom_point(show.legend = FALSE) +
            geom_ribbon(aes(ymin = coef_025, ymax = coef_975, fill = Veg), alpha = 0.15, color = NA) +
            labs(title = "No-Index Normalized Vegetation Fractions",
                 x = "Year", y = "Total Normalized Fraction", color = "Veg", fill = "Veg") +
            theme_minimal()
          ggsave(file.path(OUT_DIR, "inference_noindex_normalized_timeseries.png"), p_inf_noindex_ts, width = 8, height = 6)
          cat(sprintf("Saved inference no-index normalized time series plot to: %s\n", file.path(OUT_DIR, "inference_noindex_normalized_timeseries.png")))

        }

        if (!is.null(noindex_inf_barren) && nrow(noindex_inf_barren) > 0) {
          p_inf_noindex_barren <- ggplot(noindex_inf_barren, aes(x = year, y = global_coef)) +
            add_year_lines(is_date = FALSE) +
            geom_line(color = "saddlebrown", linewidth = 1) +
            geom_point(color = "saddlebrown") +
            geom_ribbon(aes(ymin = coef_025, ymax = coef_975), alpha = 0.15, fill = "saddlebrown", color = NA) +
            labs(title = "No-Index Based Barren Fraction", x = "Year", y = "Barren Fraction") +
            scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
            theme_minimal()
          ggsave(file.path(OUT_DIR, "inference_noindex_barren_cover.png"), p_inf_noindex_barren, width = 8, height = 6)
          cat(sprintf("Saved inference no-index based barren cover plot to: %s\n", file.path(OUT_DIR, "inference_noindex_barren_cover.png")))


          # Herbs vs Woody (inference NoIndex)
          herbs_woody_noindex_inf <- noindex_inf_full[tolower(trimws(noindex_inf_full$Veg)) %in% c("herbs", "woody"), ]
          if (nrow(herbs_woody_noindex_inf) > 0) {
            herbs_woody_noindex_inf$Veg <- ifelse(tolower(herbs_woody_noindex_inf$Veg) == "herbs", "Herbs", "Woody")
            p_inf_noindex_herbs_woody <- ggplot(herbs_woody_noindex_inf, aes(x = year, y = global_coef, color = Veg, group = Veg)) +
              add_year_lines(is_date = FALSE) +
              geom_line(linewidth = 1) +
              geom_point(show.legend = FALSE) +
              geom_ribbon(aes(ymin = coef_025, ymax = coef_975, fill = Veg), alpha = 0.15, color = NA) +
              scale_color_manual(values = c("Herbs" = "#2E8B57", "Woody" = "#8B4513")) +
              scale_fill_manual(values = c("Herbs" = "#2E8B57", "Woody" = "#8B4513")) +
              labs(title = "Inference No-Index: Herbs vs Woody Vegetation",
                   x = "Year", y = "Total Normalized Fraction", color = "Type", fill = "Type") +
              theme_minimal()
            ggsave(file.path(OUT_DIR, "inference_noindex_herbs_vs_woody.png"), p_inf_noindex_herbs_woody, width = 8, height = 6)
            cat(sprintf("Saved inference NoIndex herbs vs woody plot to: %s\n", file.path(OUT_DIR, "inference_noindex_herbs_vs_woody.png")))


            # --- Stacked area (Proportion) and Woody/Herbs ratio (Inference NoIndex) ---
            df_wide <- tryCatch({
              tidyr::pivot_wider(herbs_woody_noindex_inf |> dplyr::select(year, Veg, global_coef), names_from = Veg, values_from = global_coef)
            }, error = function(e) NULL)
            if (!is.null(df_wide) && all(c("Herbs","Woody") %in% names(df_wide))) {
              df_wide$Herbs <- as.numeric(df_wide$Herbs)
              df_wide$Woody <- as.numeric(df_wide$Woody)
              df_wide$total <- rowSums(df_wide[, c("Herbs","Woody")], na.rm = TRUE)

              df_prop <- df_wide |> dplyr::filter(is.finite(total) & total > 0) |> dplyr::mutate(Herbs = Herbs/total, Woody = Woody/total) |> tidyr::pivot_longer(cols = c("Herbs","Woody"), names_to = "Veg", values_to = "prop")
              p_inf_noindex_stacked <- ggplot(df_prop, aes(x = year, y = prop, fill = Veg)) +
                geom_area() +
                scale_fill_manual(values = c("Herbs" = "#2E8B57", "Woody" = "#8B4513")) +
                labs(title = "Inference NoIndex: Herbs vs Woody (Proportion, stacked)", x = "Year", y = "Proportion", fill = "Type") +
                theme_minimal()
              ggsave(file.path(OUT_DIR, "inference_noindex_herbs_vs_woody_stacked.png"), p_inf_noindex_stacked, width = 8, height = 6)
              cat(sprintf("Saved inference NoIndex herbs vs woody stacked plot to: %s\n", file.path(OUT_DIR, "inference_noindex_herbs_vs_woody_stacked.png")))

              df_ratio <- df_wide |> dplyr::mutate(ratio = ifelse(is.finite(Herbs) & Herbs > 0, Woody / Herbs, NA_real_))
              p_inf_noindex_ratio <- ggplot(df_ratio, aes(x = year, y = ratio)) +
                geom_line(color = "#8B4513", linewidth = 1) +
                geom_point() +
                labs(title = "Inference NoIndex: Woody / Herbs Ratio", x = "Year", y = "Woody / Herbs") +
                theme_minimal()
              ggsave(file.path(OUT_DIR, "inference_noindex_woody_over_herbs.png"), p_inf_noindex_ratio, width = 8, height = 6)
              cat(sprintf("Saved inference NoIndex woody/herbs ratio plot to: %s\n", file.path(OUT_DIR, "inference_noindex_woody_over_herbs.png")))
            } else {
              cat("Cannot create inference NoIndex stacked/ratio plot: missing Herbs/Woody rows.\n")
            }
          } else {
            cat("No herbs/woody data available for inference NoIndex herbs vs woody plot.\n")
          }
        } else {
          cat("[INFERENCE] No barren rows found in no-index inference aggregation results.\n")
        }
      } else {
        cat("[INFERENCE] No-index inference aggregation returned no results.\n")
      }
    } else {
      cat("[INFERENCE] No inference coefficients available; skipping no-index aggregation.\n")
    }

  } else {
    cat("[INFERENCE] No inference data to process\n")
    inference_coefs <- NULL
  }
  
  # Combine validation and inference results for downstream processing
  all_coefs <- rbind(validation_coefs, inference_coefs)
  cat(sprintf("[COMBINED] Total coefficient rows: %d (validation: %d, inference: %d)\n",
              if(!is.null(all_coefs)) nrow(all_coefs) else 0,
              if(!is.null(validation_coefs)) nrow(validation_coefs) else 0,
              if(!is.null(inference_coefs)) nrow(inference_coefs) else 0))

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

  # Training results CSV saving is disabled by default (user requested removal of training outputs)
  if (!is.null(all_coefs) && nrow(all_coefs) > 0) {
    cat("[INFO] Training results CSV saving skipped (training outputs removed by config)\n")
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
  n_year_results <- 0
  if (!is.null(all_coefs) && nrow(all_coefs) > 0) {
    n_year_results <- length(unique(paste(all_coefs$location_id, all_coefs$pheno_year, sep = "_")))
  }

  if (n_year_results > 0 && exists("processing_time")) {
    cat(sprintf("Average time per year-result: %.2f seconds\n", processing_time / n_year_results))
  } else {
    cat("Average time per year-result: N/A (0 results)\n")
  }

  # Ensure required library/templates exist for visualization
  ensure_library_and_templates()

  # Ensure library + templates exist
  ensure_library_and_templates()

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
    
    # --- 1. Identify Valid Locations and Get True Labels from HELD-OUT validation set ---
    labels_df <- NULL
    val_locations_file <- file.path(OUT_DIR, "validation_locations.csv")
    
    # PRIORITY 1: Use the held-out validation_locations.csv (20% stratified split created during training)
    if (file.exists(val_locations_file)) {
      vloc <- tryCatch(read.csv(val_locations_file, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(vloc) && "location_id" %in% names(vloc) && "Veg" %in% names(vloc)) {
        labels_df <- vloc %>% dplyr::select(location_id, Veg) %>% dplyr::distinct()
        cat(sprintf("[VALIDATION] Using held-out validation set from %s (%d locations)\n", 
                    val_locations_file, nrow(labels_df)))
        # Log class distribution in held-out set
        held_out_dist <- table(labels_df$Veg)
        cat(sprintf("[VALIDATION] Held-out set class distribution: %s\n",
                    paste(names(held_out_dist), "=", held_out_dist, collapse=", ")))
      }
    }
    
    # FALLBACK: Try df_tasks_inference or df_tasks if validation file not found
    if (is.null(labels_df)) {
      cat("[WARNING] validation_locations.csv not found, falling back to df_tasks labels (not a true held-out set)\n")
      if (exists("df_tasks_inference") && !is.null(df_tasks_inference) && "Veg" %in% names(df_tasks_inference)) {
         labels_df <- df_tasks_inference %>% dplyr::select(location_id, Veg) %>% dplyr::distinct()
      } else if (exists("df_tasks") && !is.null(df_tasks) && "Veg" %in% names(df_tasks)) {
         labels_df <- df_tasks %>% dplyr::select(location_id, Veg) %>% dplyr::distinct()
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

    # AGENT: Filter validation to TRAIN_YEARS only (validation uses same years as training)
    if ("pheno_year" %in% names(val_coefs_labeled) && exists("TRAIN_YEARS") && !is.null(TRAIN_YEARS)) {
      orig_n <- nrow(val_coefs_labeled)
      val_coefs_labeled <- val_coefs_labeled %>% dplyr::filter(pheno_year %in% TRAIN_YEARS)
      if (!isTRUE(QUIET_MODE)) cat(sprintf("[VALIDATION] Filtered to TRAIN_YEARS=%s: %d -> %d rows\n", 
                  paste(TRAIN_YEARS, collapse=", "), orig_n, nrow(val_coefs_labeled)))
    }
    
    # Check class distribution in validation set
    if (nrow(val_coefs_labeled) > 0 && "true_veg" %in% names(val_coefs_labeled)) {
      val_class_dist <- table(val_coefs_labeled$true_veg)
      cat(sprintf("[VALIDATION] Class distribution: %s\n", 
                  paste(names(val_class_dist), "=", val_class_dist, collapse=", ")))
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
       
       # AGENT: Aggregate variants to vegetation class level (e.g. "herbs_opt_1" -> "herbs")
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

    # --- 4. Row-Normalized Confusion Matrix (excluding barren) ---
    if (length(all_veg_cols) > 0 && nrow(val_coefs_wide) > 0) {
      # Filter to non-barren vegetation classes only
      matrix_veg_cols <- all_veg_cols[!grepl("barren", all_veg_cols, ignore.case = TRUE)]
      
      if (length(matrix_veg_cols) == 0) {
        cat("\n[NOTICE] No non-barren vegetation columns - skipping confusion matrix\n")
      } else {
        # Filter rows to non-barren true classes
        val_coefs_veg <- val_coefs_wide %>%
          dplyr::filter(tolower(true_veg) != "barren", !is.na(true_veg), true_veg != "")
        
        n_true_classes <- length(unique(val_coefs_veg$true_veg))
        
        if (n_true_classes < 2) {
          cat(sprintf("\n[NOTICE] Only %d non-barren class in validation - skipping confusion matrix\n", n_true_classes))
        } else {
          cat("\n=== CONFUSION MATRIX (Row-Normalized, Excluding Barren) ===\n")
          cat("Rows: True Class | Columns: Mean Predicted Fraction\n")
          cat("(Each row normalized to sum to 1.0)\n\n")
          
          # For each prediction, normalize vegetation fractions (excluding barren) to sum to 1
          # NOTE: When veg_sum is near zero, preserve original values instead of zeroing them
          val_coefs_norm <- val_coefs_veg %>%
            dplyr::mutate(
              veg_sum = rowSums(dplyr::across(dplyr::all_of(matrix_veg_cols)), na.rm = TRUE)
            ) %>%
            dplyr::mutate(
              dplyr::across(
                dplyr::all_of(matrix_veg_cols),
                ~ ifelse(veg_sum > 1e-9, .x / veg_sum, .x),
                .names = "{.col}_norm"
              )
            )
          
          # Column names after normalization
          norm_cols <- paste0(matrix_veg_cols, "_norm")
          
          # Average by true class
          avg_fractions <- val_coefs_norm %>%
            dplyr::group_by(true_veg) %>%
            dplyr::summarize(
              dplyr::across(dplyr::all_of(norm_cols), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop"
            )
          
          # Clean column names (remove frac_ prefix and _norm suffix)
          clean_names <- sub("_norm$", "", colnames(avg_fractions))
          clean_names <- sub("^frac_", "", clean_names)
          colnames(avg_fractions) <- clean_names
          
          # Convert to matrix
          mat_data <- as.matrix(avg_fractions[, -1])  # Remove true_veg column
          rownames(mat_data) <- avg_fractions$true_veg
          
          # Ensure square matrix with all classes in both dimensions
          all_classes <- sort(unique(c(rownames(mat_data), colnames(mat_data))))
          
          # Add missing columns
          for (cls in all_classes) {
            if (!cls %in% colnames(mat_data)) {
              new_col <- matrix(0, nrow = nrow(mat_data), ncol = 1)
              colnames(new_col) <- cls
              mat_data <- cbind(mat_data, new_col)
            }
          }
          
          # Add missing rows
          for (cls in all_classes) {
            if (!cls %in% rownames(mat_data)) {
              new_row <- matrix(0, nrow = 1, ncol = ncol(mat_data))
              rownames(new_row) <- cls
              mat_data <- rbind(mat_data, new_row)
            }
          }
          
          # Reorder to match all_classes
          mat_data <- mat_data[all_classes, all_classes, drop = FALSE]
          
          # Row-normalize (already done per-prediction, but verify)
          mat_rownorm <- mat_data
          row_sums <- rowSums(mat_rownorm, na.rm = TRUE)
          for (i in seq_len(nrow(mat_rownorm))) {
            if (row_sums[i] > 1e-9) {
              mat_rownorm[i, ] <- mat_rownorm[i, ] / row_sums[i]
            }
          }
          
          # Display matrix
          print(round(mat_rownorm, 3))
          
          # Compute diagonal accuracy
          diag_vals <- diag(mat_rownorm)
          mean_diagonal <- mean(diag_vals, na.rm = TRUE)
          
          cat(sprintf("\nMean diagonal (correctly predicted fraction): %.3f (%.1f%%)\n",
                      mean_diagonal, mean_diagonal * 100))
          cat(sprintf("Per-class accuracy:\n"))
          for (i in seq_along(all_classes)) {
            cat(sprintf("  %s: %.3f (%.1f%%)\n",
                        all_classes[i], diag_vals[i], diag_vals[i] * 100))
          }

          # Store confusion matrix for Dirichlet perturbation in bootstrap
          # Compute validation sample sizes per class
          val_sample_sizes <- table(val_coefs_veg$true_veg)
          val_sample_sizes <- setNames(as.numeric(val_sample_sizes), names(val_sample_sizes))

          store_confusion_matrix(mat_rownorm, sample_sizes = val_sample_sizes)
          cat(sprintf("\n[DIRICHLET] Confusion matrix stored for classification uncertainty propagation\n"))

          # Store OOB fraction residuals for MC uncertainty propagation
          # Residual = predicted_fraction - true_fraction (where true is 1 for own class, 0 for others)
          if (isTRUE(ENABLE_OOB_FRACTION_UNCERTAINTY)) {
            cat("\n[OOB_FRAC] Computing and storing OOB fraction residuals for MC uncertainty...\n")

            # Build residuals by true class
            residuals_by_class <- list()

            for (true_cls in all_classes) {
              # Get rows where true class matches
              cls_data <- val_coefs_norm[val_coefs_norm$true_veg == true_cls, , drop = FALSE]
              if (nrow(cls_data) == 0) next

              # Build true fraction vector: 1 for own class, 0 for others
              # Compute residuals: predicted - true
              n_samples <- nrow(cls_data)
              resid_mat <- matrix(0, nrow = n_samples, ncol = length(all_classes))
              colnames(resid_mat) <- all_classes

              for (i in seq_len(n_samples)) {
                for (pred_cls in all_classes) {
                  # Get predicted fraction (normalized)
                  pred_col <- paste0("frac_", pred_cls, "_norm")
                  pred_frac <- if (pred_col %in% names(cls_data)) cls_data[[pred_col]][i] else 0
                  if (!is.finite(pred_frac)) pred_frac <- 0

                  # True fraction: 1 if pred_cls == true_cls, else 0
                  true_frac <- if (pred_cls == true_cls) 1.0 else 0.0

                  # Residual = predicted - true
                  resid_mat[i, pred_cls] <- pred_frac - true_frac
                }
              }

              residuals_by_class[[true_cls]] <- resid_mat
            }

            if (length(residuals_by_class) > 0) {
              store_oob_fraction_residuals(residuals_by_class)
              cat(sprintf("[OOB_FRAC] Stored residuals for %d classes\n", length(residuals_by_class)))
            }
          }
        }
      }
    }

    # --- 5. Artificial Mix Logic ---
    if (length(all_veg_cols) > 0) {
      cat("DEBUG: Entering artificial mix section\n")
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
             mixed_coefs <- (sub_a[idx_a, all_veg_cols] + sub_b[idx_b, all_veg_cols]) / 2
             
             # Calculate relative fraction of Class A in the mix
             # Ideally should be 0.5 if A was 1.0 and B was 0.0 for A
             
             frac_col_a <- paste0("frac_", tolower(class_a))
             
             if (frac_col_a %in% names(mixed_coefs)) {
               # Rescale to total vegetation
               barren_col <- "frac_barren"
               veg_frac_cols_local <- all_veg_cols[all_veg_cols != barren_col]
               
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

    # --- 6. Compute Per-Class Prediction RMSE (for CI Inflation) ---
    per_class_rmse <- numeric(0)
    if (length(all_veg_cols) > 0 && nrow(val_coefs_wide) > 0) {
      cat("\n--- VALIDATION RMSE (for CI adjustment) ---\n")
      for (v_class in unique(c(val_coefs_wide$true_veg, sub("^frac_", "", all_veg_cols)))) {
        frac_col <- paste0("frac_", tolower(v_class))
        if(frac_col %in% names(val_coefs_wide)) {
           predicted <- val_coefs_wide[[frac_col]]
           # Truth is 1 if true_veg matches class, 0 otherwise
           truth <- ifelse(tolower(val_coefs_wide$true_veg) == tolower(v_class), 1, 0)
           
           # Filter out NAs
           valid_idx <- !is.na(predicted) & !is.na(truth)
           if (sum(valid_idx) > 2) {
             mse <- mean((predicted[valid_idx] - truth[valid_idx])^2)
             rmse_val <- sqrt(mse)
             per_class_rmse[tolower(v_class)] <- rmse_val
             cat(sprintf("  %s: RMSE = %.4f (N=%d)\n", v_class, rmse_val, sum(valid_idx)))
           }
        }
      }
    }

    cat("Validation accuracy computed", prefix, ".\n", sep = "")
    return(per_class_rmse)
  }

  # Compute accuracy on VALIDATION data only (held-out locations, TRAIN_YEARS)
  validation_rmse_adjustment <- NULL
  if (exists("validation_coefs") && !is.null(validation_coefs) && nrow(validation_coefs) > 0) {
    cat("\n=== VALIDATION ACCURACY (held-out locations, TRAIN_YEARS only) ===\n")
    val_coefs <- validation_coefs
    validation_rmse_adjustment <- report_validation_accuracy(val_coefs, "held-out validation set")
  } else {
    cat("[WARNING] No validation coefficients available for accuracy computation\n")
  }

  # ==========================================================================
  # ARTIFICIAL 50/50 MIX VALIDATION (Spectral-level mixing)
  # Create synthetic mixed spectra from validation samples and validate unmixing
  # ==========================================================================
  cat("\n=== ARTIFICIAL 50/50 MIX VALIDATION (Spectral-level) ===\n")

  artificial_mix_validation <- function() {
    # Check prerequisites
    if (!exists("df_validation") || is.null(df_validation) || nrow(df_validation) == 0) {
      cat("[SKIP] No df_validation available for artificial mixing\n")
      return(NULL)
    }

    if (!exists("SPLINE_PARAMS") || is.null(SPLINE_PARAMS) || !("indices" %in% names(SPLINE_PARAMS))) {
      cat("[SKIP] SPLINE_PARAMS$indices not available - cannot identify spectral columns\n")
      return(NULL)
    }

    spectral_cols <- SPLINE_PARAMS$indices
    if (length(spectral_cols) == 0) {
      cat("[SKIP] No spectral columns identified\n")
      return(NULL)
    }

    # Check that spectral columns exist in df_validation
    available_cols <- intersect(spectral_cols, names(df_validation))
    if (length(available_cols) < length(spectral_cols) * 0.5) {
      cat(sprintf("[SKIP] Too few spectral columns available in df_validation (%d/%d)\n",
                  length(available_cols), length(spectral_cols)))
      return(NULL)
    }
    spectral_cols <- available_cols

    # Get vegetation classes (excluding barren for mixing)
    df_val_with_veg <- df_validation
    if (!"Veg" %in% names(df_val_with_veg)) {
      cat("[SKIP] No Veg column in df_validation\n")
      return(NULL)
    }
    df_val_with_veg$Veg <- tolower(trimws(as.character(df_val_with_veg$Veg)))
    veg_classes <- unique(df_val_with_veg$Veg[df_val_with_veg$Veg != "barren"])

    if (length(veg_classes) < 2) {
      cat(sprintf("[SKIP] Need at least 2 vegetation classes for mixing, found: %d\n", length(veg_classes)))
      return(NULL)
    }

    cat(sprintf("[MIX] Using %d spectral columns and %d vegetation classes: %s\n",
                length(spectral_cols), length(veg_classes), paste(veg_classes, collapse=", ")))

    # Filter to TRAIN_YEARS if defined
    if (exists("TRAIN_YEARS") && !is.null(TRAIN_YEARS) && "pheno_year" %in% names(df_val_with_veg)) {
      df_val_with_veg <- df_val_with_veg[df_val_with_veg$pheno_year %in% TRAIN_YEARS, ]
      cat(sprintf("[MIX] Filtered to TRAIN_YEARS: %d rows\n", nrow(df_val_with_veg)))
    }

    # Aggregate spectra to location-year level (mean across observations)
    # This gives us one "representative spectrum" per location-year
    required_cols <- c("location_id", "pheno_year", "Veg", spectral_cols)
    if ("doy" %in% names(df_val_with_veg)) required_cols <- c(required_cols, "doy")
    if ("date" %in% names(df_val_with_veg)) required_cols <- c(required_cols, "date")
    if ("lat" %in% names(df_val_with_veg)) required_cols <- c(required_cols, "lat")
    if ("lon" %in% names(df_val_with_veg)) required_cols <- c(required_cols, "lon")

    df_val_subset <- df_val_with_veg[, intersect(required_cols, names(df_val_with_veg)), drop = FALSE]

    # Get unique location-years per class
    loc_year_by_class <- list()
    for (vc in veg_classes) {
      class_data <- df_val_subset[df_val_subset$Veg == vc, ]
      if (nrow(class_data) > 0) {
        # Get unique location-year combinations
        class_data$loc_year <- paste(class_data$location_id, class_data$pheno_year, sep = "_")
        loc_year_by_class[[vc]] <- unique(class_data$loc_year)
      }
    }

    # Generate all pairs of vegetation classes
    class_pairs <- utils::combn(veg_classes, 2, simplify = FALSE)

    mix_results <- list()
    n_mixes_total <- 0

    for (pair in class_pairs) {
      class_a <- pair[1]
      class_b <- pair[2]

      loc_years_a <- loc_year_by_class[[class_a]]
      loc_years_b <- loc_year_by_class[[class_b]]

      if (is.null(loc_years_a) || is.null(loc_years_b) ||
          length(loc_years_a) == 0 || length(loc_years_b) == 0) {
        cat(sprintf("[MIX] Skipping %s + %s: insufficient samples\n", class_a, class_b))
        next
      }

      # Create up to 50 mixes per pair
      n_mix <- min(length(loc_years_a), length(loc_years_b), 50)

      cat(sprintf("[MIX] Creating %d synthetic 50/50 mixes: %s + %s\n", n_mix, class_a, class_b))

      set.seed(42 + which(sapply(class_pairs, function(p) identical(p, pair))))
      sampled_a <- sample(loc_years_a, n_mix, replace = TRUE)
      sampled_b <- sample(loc_years_b, n_mix, replace = TRUE)

      mixed_data_list <- list()

      for (i in seq_len(n_mix)) {
        # Get data for each component
        loc_year_a <- sampled_a[i]
        loc_year_b <- sampled_b[i]

        parts_a <- strsplit(loc_year_a, "_")[[1]]
        parts_b <- strsplit(loc_year_b, "_")[[1]]

        loc_a <- paste(parts_a[-length(parts_a)], collapse = "_")
        yr_a <- as.integer(parts_a[length(parts_a)])
        loc_b <- paste(parts_b[-length(parts_b)], collapse = "_")
        yr_b <- as.integer(parts_b[length(parts_b)])

        data_a <- df_val_subset[df_val_subset$location_id == loc_a &
                                 df_val_subset$pheno_year == yr_a &
                                 df_val_subset$Veg == class_a, ]
        data_b <- df_val_subset[df_val_subset$location_id == loc_b &
                                 df_val_subset$pheno_year == yr_b &
                                 df_val_subset$Veg == class_b, ]

        if (nrow(data_a) == 0 || nrow(data_b) == 0) next

        # Match observations by DOY if possible, otherwise use mean spectra
        if ("doy" %in% names(data_a) && "doy" %in% names(data_b)) {
          # Find common or closest DOYs
          common_doys <- intersect(data_a$doy, data_b$doy)

          if (length(common_doys) >= 3) {
            # Use common DOYs
            data_a_matched <- data_a[data_a$doy %in% common_doys, ]
            data_b_matched <- data_b[data_b$doy %in% common_doys, ]

            # Ensure same number of rows by matching DOYs
            data_a_matched <- data_a_matched[order(data_a_matched$doy), ]
            data_b_matched <- data_b_matched[order(data_b_matched$doy), ]

            # Take first occurrence per DOY
            data_a_matched <- data_a_matched[!duplicated(data_a_matched$doy), ]
            data_b_matched <- data_b_matched[!duplicated(data_b_matched$doy), ]

            common_doys_final <- intersect(data_a_matched$doy, data_b_matched$doy)
            data_a_matched <- data_a_matched[data_a_matched$doy %in% common_doys_final, ]
            data_b_matched <- data_b_matched[data_b_matched$doy %in% common_doys_final, ]
          } else {
            # Not enough common DOYs, aggregate to mean
            data_a_matched <- data_a[1, , drop = FALSE]
            data_b_matched <- data_b[1, , drop = FALSE]
            for (sc in spectral_cols) {
              data_a_matched[[sc]] <- mean(data_a[[sc]], na.rm = TRUE)
              data_b_matched[[sc]] <- mean(data_b[[sc]], na.rm = TRUE)
            }
            data_a_matched$doy <- median(data_a$doy, na.rm = TRUE)
            data_b_matched$doy <- median(data_b$doy, na.rm = TRUE)
          }
        } else {
          # No DOY, use first row with mean values
          data_a_matched <- data_a[1, , drop = FALSE]
          data_b_matched <- data_b[1, , drop = FALSE]
          for (sc in spectral_cols) {
            data_a_matched[[sc]] <- mean(data_a[[sc]], na.rm = TRUE)
            data_b_matched[[sc]] <- mean(data_b[[sc]], na.rm = TRUE)
          }
        }

        # Create 50/50 mixed spectra
        n_obs <- min(nrow(data_a_matched), nrow(data_b_matched))
        if (n_obs == 0) next

        mixed_data <- data_a_matched[1:n_obs, , drop = FALSE]
        mixed_data$location_id <- sprintf("MIX_%s_%s_%d", class_a, class_b, i)
        mixed_data$pheno_year <- yr_a  # Use year from A
        mixed_data$Veg <- paste0("mix_", class_a, "_", class_b)
        mixed_data$true_frac_a <- 0.5
        mixed_data$true_frac_b <- 0.5
        mixed_data$class_a <- class_a
        mixed_data$class_b <- class_b

        # Mix spectral values: 0.5 * A + 0.5 * B
        for (sc in spectral_cols) {
          val_a <- data_a_matched[[sc]][1:n_obs]
          val_b <- data_b_matched[[sc]][1:n_obs]
          if (length(val_a) == n_obs && length(val_b) == n_obs) {
            mixed_data[[sc]] <- 0.5 * val_a + 0.5 * val_b
          }
        }

        mixed_data_list[[length(mixed_data_list) + 1]] <- mixed_data
      }

      if (length(mixed_data_list) > 0) {
        pair_key <- paste(class_a, class_b, sep = "_")
        mix_results[[pair_key]] <- list(
          class_a = class_a,
          class_b = class_b,
          mixed_data = do.call(rbind, mixed_data_list)
        )
        n_mixes_total <- n_mixes_total + nrow(mix_results[[pair_key]]$mixed_data)
      }
    }

    if (n_mixes_total == 0) {
      cat("[MIX] No valid mixes created\n")
      return(NULL)
    }

    cat(sprintf("[MIX] Created %d total artificial mix samples across %d class pairs\n",
                n_mixes_total, length(mix_results)))

    # Run unmixing on mixed data
    cat("[MIX] Running unmixing on artificial mixes...\n")

    all_mix_coefs <- list()

    for (pair_key in names(mix_results)) {
      mr <- mix_results[[pair_key]]
      mixed_df <- mr$mixed_data

      if (is.null(mixed_df) || nrow(mixed_df) == 0) next

      # Ensure required columns for unmixing
      if (!"date" %in% names(mixed_df) && "doy" %in% names(mixed_df) && "pheno_year" %in% names(mixed_df)) {
        # Approximate date from DOY and year
        mixed_df$date <- as.Date(paste0(mixed_df$pheno_year, "-01-01")) + mixed_df$doy - 1
      }

      # Run inference using run_inference_silent if available
      if (exists("run_inference_silent")) {
        inf_result <- tryCatch({
          run_inference_silent(mixed_df)
        }, error = function(e) {
          cat(sprintf("[MIX ERROR] %s: %s\n", pair_key, e$message))
          NULL
        })

        if (!is.null(inf_result) && !is.null(inf_result$all_coefs) && nrow(inf_result$all_coefs) > 0) {
          coefs <- inf_result$all_coefs
          coefs$class_a <- mr$class_a
          coefs$class_b <- mr$class_b
          coefs$true_frac_a <- 0.5
          coefs$true_frac_b <- 0.5
          all_mix_coefs[[pair_key]] <- coefs
        }
      }
    }

    if (length(all_mix_coefs) == 0) {
      cat("[MIX] Unmixing produced no results\n")
      return(NULL)
    }

    mix_coefs_combined <- do.call(rbind, all_mix_coefs)

    # Analyze results
    cat("\n--- ARTIFICIAL MIX VALIDATION RESULTS ---\n")
    cat(sprintf("Total unmixed mix samples: %d\n", length(unique(mix_coefs_combined$location_id))))

    # For each class pair, check if predicted fractions match expected 50/50
    for (pair_key in names(all_mix_coefs)) {
      coefs <- all_mix_coefs[[pair_key]]
      class_a <- coefs$class_a[1]
      class_b <- coefs$class_b[1]

      # Pivot to wide format
      coefs_wide <- coefs %>%
        dplyr::mutate(Veg_class = sub("(_opt_[0-9]+|_single|_ppi)$", "", tolower(Veg))) %>%
        dplyr::group_by(location_id, pheno_year, Veg_class) %>%
        dplyr::summarise(coef = sum(coef, na.rm = TRUE), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = Veg_class, values_from = coef, values_fill = 0)

      frac_col_a <- tolower(class_a)
      frac_col_b <- tolower(class_b)

      if (frac_col_a %in% names(coefs_wide) && frac_col_b %in% names(coefs_wide)) {
        pred_frac_a <- coefs_wide[[frac_col_a]]
        pred_frac_b <- coefs_wide[[frac_col_b]]

        # Calculate metrics
        mean_pred_a <- mean(pred_frac_a, na.rm = TRUE)
        mean_pred_b <- mean(pred_frac_b, na.rm = TRUE)
        sd_pred_a <- sd(pred_frac_a, na.rm = TRUE)
        sd_pred_b <- sd(pred_frac_b, na.rm = TRUE)

        # Expected is 0.5 for both
        rmse_a <- sqrt(mean((pred_frac_a - 0.5)^2, na.rm = TRUE))
        rmse_b <- sqrt(mean((pred_frac_b - 0.5)^2, na.rm = TRUE))

        cat(sprintf("\n  Mix %s + %s (N=%d):\n", class_a, class_b, nrow(coefs_wide)))
        cat(sprintf("    %s: predicted=%.3f ± %.3f (expected=0.500, RMSE=%.3f)\n",
                    class_a, mean_pred_a, sd_pred_a, rmse_a))
        cat(sprintf("    %s: predicted=%.3f ± %.3f (expected=0.500, RMSE=%.3f)\n",
                    class_b, mean_pred_b, sd_pred_b, rmse_b))
        cat(sprintf("    Sum: %.3f (expected=1.000)\n", mean_pred_a + mean_pred_b))
      }
    }

    return(mix_coefs_combined)
  }

  # Run artificial mix validation
  artificial_mix_coefs <- tryCatch({
    artificial_mix_validation()
  }, error = function(e) {
    cat(sprintf("[MIX ERROR] Artificial mix validation failed: %s\n", e$message))
    NULL
  })

  cat("=== END ARTIFICIAL MIX VALIDATION ===\n\n")

  # Results are now in validation_coefs and inference_coefs (combined as all_coefs)
  cat(sprintf("[RESULTS] Validation: %d rows, Inference: %d rows, Combined: %d rows\n",
              if(exists("validation_coefs") && !is.null(validation_coefs)) nrow(validation_coefs) else 0,
              if(exists("inference_coefs") && !is.null(inference_coefs)) nrow(inference_coefs) else 0,
              if(exists("all_coefs") && !is.null(all_coefs)) nrow(all_coefs) else 0))
              
  # --- Apply Validation RMSE Adjustment to CIs ---
  if (!is.null(validation_rmse_adjustment) && length(validation_rmse_adjustment) > 0) {
    cat("\n[UNCERTAINTY ADJUSTMENT] Incorporating validation RMSE into Confidence Intervals...\n")
    
    # Function to adjust CIs
    adjust_ci <- function(df, rmse_vec) {
      if (is.null(df) || nrow(df) == 0) return(df)
      changed_rows <- 0
      
      # Iterate over known classes in RMSE vector
      for (v_class in names(rmse_vec)) {
        rmse_val <- rmse_vec[[v_class]]
        # Find rows for this class (case-insensitive)
        idx <- which(tolower(trimws(df$Veg)) == tolower(v_class))
        
        if (length(idx) > 0) {
           # Extract current intervals
           current_lower <- df$coef_025[idx]
           current_upper <- df$coef_975[idx]
           current_coef <- df$coef[idx]
           
           # Approximate standard error including validation RMSE
           # Current Width ~ 2 * 1.96 * SE_fit
           # New SE ~ sqrt(SE_fit^2 + RMSE_val^2)
           
           # Recover implicit SE from CI width
           # If width is 0 (or NA), assume SE is 0 (or skip)
           width <- current_upper - current_lower
           width[is.na(width)] <- 0
           implicit_se <- width / (2 * 1.96)
           
           # Combine variances
           new_se <- sqrt(implicit_se^2 + rmse_val^2)
           
           # Construct new intervals centered on prediction
           # (We assume symmetric error distribution governed by new_se)
           new_lower <- current_coef - 1.96 * new_se
           new_upper <- current_coef + 1.96 * new_se
           
           # Clamp to [0, 1]
           new_lower <- pmax(0, pmin(1, new_lower))
           new_upper <- pmax(0, pmin(1, new_upper))
           
           # Update DF
           df$coef_025[idx] <- new_lower
           df$coef_975[idx] <- new_upper
           df$coef_sd[idx] <- new_se # Update SD if present
           
           changed_rows <- changed_rows + length(idx)
        }
      }
      return(df)
    }

    # Apply to all_coefs (which combines val & inference)
    if (!is.null(all_coefs) && nrow(all_coefs) > 0) {
       all_coefs <- adjust_ci(all_coefs, validation_rmse_adjustment)
    }
    # Also update constituent dataframes so they stay in sync if used later
    if (exists("inference_coefs") && !is.null(inference_coefs) && nrow(inference_coefs) > 0) {
       inference_coefs <- adjust_ci(inference_coefs, validation_rmse_adjustment)
    }
    if (exists("validation_coefs") && !is.null(validation_coefs) && nrow(validation_coefs) > 0) {
       validation_coefs <- adjust_ci(validation_coefs, validation_rmse_adjustment)
    }
    cat(sprintf("[UNCERTAINTY ADJUSTMENT] Adjusted CIs for ~%d rows across %d classes.\n", 
                nrow(all_coefs), length(validation_rmse_adjustment)))
  }

  # Fail fast if BOTH validation and inference result lists are missing/empty — indicates upstream processing issue
  if ((!exists("validation_coefs") || is.null(validation_coefs) || nrow(validation_coefs) == 0) &&
      (!exists("inference_coefs") || is.null(inference_coefs) || nrow(inference_coefs) == 0)) {
    stop(paste0("ERROR: No results collected (both validation and inference are empty). Upstream processing failed — check earlier logs and data filters."))
  }

  cat("Processing results...\n")

  # Diagnostics: report how many results were collected
  cat(sprintf("Validation results: %d rows\n", if(exists("validation_coefs") && !is.null(validation_coefs)) nrow(validation_coefs) else 0))
  cat(sprintf("Inference results: %d rows\n", if(exists("inference_coefs") && !is.null(inference_coefs)) nrow(inference_coefs) else 0))

  # Check barren fraction distribution in combined results
  if (!is.null(all_coefs) && nrow(all_coefs) > 0 && "Veg" %in% names(all_coefs)) {
    barren_rows <- all_coefs[tolower(all_coefs$Veg) == "barren", ]
    if (nrow(barren_rows) > 0 && "coef" %in% names(barren_rows)) {
      barren_one_count <- sum(barren_rows$coef == 1, na.rm = TRUE)
      barren_one_pct <- barren_one_count / nrow(all_coefs) * 100

      # ERROR if ALL samples are 100% barren
      if (barren_one_pct >= 99.9) {
        stop(sprintf(
          "[FATAL ERROR] All samples are 100%% barren (%.1f%%, %d/%d predictions)!\n\n" ,
          "This indicates a critical unmixing failure. Possible causes:\n",
          "1. PPI calibration issue: FVC models may not be properly calibrated\n",
          "2. Library mismatch: Endmember library may not match the inference data characteristics\n",
          "3. Barren shortcut threshold too aggressive: Try lowering BARREN_SHORTCUT_THRESHOLD\n",
          "4. Data preprocessing issue: Check that spectral indices are calculated correctly\n",
          "\nPlease review the debug log and check PPI values and library quality.\n"
        ), barren_one_pct, barren_one_count, nrow(all_coefs))
      }

      cat(sprintf("Barren fraction = 1 in %.1f%% of predictions (%d/%d) - within acceptable limits\n",
                  barren_one_pct, barren_one_count, nrow(all_coefs)))
    }
  }

  # Verify we have results to process
  if (is.null(all_coefs) || nrow(all_coefs) == 0) {
    cat("ERROR: No results to process!\n")
    cat("Most likely causes:\n")
    cat("1. Insufficient data quality: Tasks were filtered due to:\n")
    cat("   - Too few observations per location-year (<15 observations)\n")
    cat("   - Insufficient temporal coverage (<25 unique days of year)\n")
    cat("   - Too many missing pentad bins (>85% NA values)\n")
    cat("2. Data filtering issues or missing indices\n")
    cat("3. No valid validation/inference data available (no location-year pairs found)\n")
    cat("\nTo adjust data quality thresholds, modify MIN_OBS, MIN_UNIQUE_DOYS, and MIN_PENTAD_COVERAGE in fit_one_task function.\n")
    stop("No valid results to process")
  }

  # all_coefs is already built from combined validation_coefs and inference_coefs
  # Ensure required columns and data types
  cat(sprintf("Combined coefficients: %d rows\n", nrow(all_coefs)))
  
  required_coef_cols <- c("location_id", "pheno_year", "Veg", "coef", "rmse", "coef_025", "coef_975", "interval")
  missing <- setdiff(required_coef_cols, names(all_coefs))
  if (length(missing) > 0) {
    cat(sprintf("[NOTICE] Filling missing coefficient columns with NA: %s\n", paste(missing, collapse = ", ")))
    for (col in missing) all_coefs[[col]] <- NA
  }
  all_coefs$location_id <- as.character(all_coefs$location_id)
  all_coefs$pheno_year <- as.integer(all_coefs$pheno_year)
  # Ensure backward-compatible 'year' column exists for downstream code that expects it
  if (!"year" %in% names(all_coefs)) all_coefs$year <- all_coefs$pheno_year
  all_coefs$Veg <- as.character(all_coefs$Veg)
  all_coefs$coef <- as.numeric(all_coefs$coef)
  if ("rmse" %in% names(all_coefs)) all_coefs$rmse <- as.numeric(all_coefs$rmse)
  if ("coef_025" %in% names(all_coefs)) all_coefs$coef_025 <- as.numeric(all_coefs$coef_025)
  if ("coef_975" %in% names(all_coefs)) all_coefs$coef_975 <- as.numeric(all_coefs$coef_975)
  if ("interval" %in% names(all_coefs)) all_coefs$interval <- as.numeric(all_coefs$interval)
  all_coefs$location_id <- trimws(all_coefs$location_id)

  # Note: variant_trajectory and diagnostics are not currently tracked in the new processing flow
  all_variants_pca <- NULL
  all_diagnostics <- NULL

  q_dvi_data <- NULL

  # Variant similarity (this runs normally - not inside the if(FALSE) block)
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

    unique_locations <- unique(trimws(as.character(all_coefs$location_id)))
    unique_locations <- unique_locations[!is.na(unique_locations) & unique_locations != ""]
    if (length(unique_locations) == 0) {
      stop("No valid location IDs found in results")
    }

    true_veg_map <- gpts_map |> dplyr::select(location_id, true_veg = Veg)
    if (!is.character(all_coefs$location_id)) {
      all_coefs$location_id <- as.character(all_coefs$location_id)
    }
    if (!is.character(true_veg_map$location_id)) {
      true_veg_map$location_id <- as.character(true_veg_map$location_id)
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

      # Wrapper function to handle both spline and standard modes
      fit_task_wrapper <- function(task_data) {
        if (is.null(task_data) || nrow(task_data) == 0) return(NULL)
        loc <- as.character(task_data$location_id[1])
        yr <- as.integer(task_data$pheno_year[1])
        
        if (isTRUE(USE_SPLINE_ENDMEMBERS) && exists("SPLINE_LIBRARY") && exists("SPLINE_PARAMS")) {
          res <- tryCatch({
            fit_one_task_spline(task_data, SPLINE_LIBRARY, SPLINE_PARAMS$indices, SPLINE_PARAMS, loc, yr)
          }, error = function(e) NULL)
        } else {
          res <- fit_one_task(task_data)
        }
        return(res)
      }

      sink_file <- tempfile()
      con_out <- tryCatch(file(sink_file, open = "wt"), error = function(e) NULL)
      con_msg <- tryCatch(file(paste0(sink_file, ".msg"), open = "wt"), error = function(e) NULL)
      # Suppress output and messages while running heavy inference
      results <- NULL
      tryCatch({
        if (!is.null(con_out) && inherits(con_out, "connection") && isOpen(con_out)) {
          sink(con_out, type = "output")
        }
        if (!is.null(con_msg) && inherits(con_msg, "connection") && isOpen(con_msg)) {
          sink(con_msg, type = "message")
        }
        results <- .run_map(task_list, fit_task_wrapper, show_pb = FALSE)
      }, error = function(e) {
        warning(sprintf("run_inference_silent: underlying inference error: %s", e$message))
      }, finally = {
        # Attempt to restore sinks
        try(sink(type = "message"), silent = TRUE)
        try(sink(type = "output"), silent = TRUE)
        if (!is.null(con_out)) try(close(con_out), silent = TRUE)
        if (!is.null(con_msg)) try(close(con_msg), silent = TRUE)
        try(unlink(sink_file), silent = TRUE)
        try(unlink(paste0(sink_file, ".msg")), silent = TRUE)
      })

      all_coefs_local <- combine_results_from_list(results)
      list(results = results, all_coefs = all_coefs_local)
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

    # Collect uncertainty data for merging
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

    best_fit_summary <- do.call(rbind, lapply(unique_locations, function(loc) {
      yrs <- unique(all_coefs$pheno_year[all_coefs$location_id == loc])
      tv <- true_veg_map$true_veg[true_veg_map$location_id == loc]
      tv_val <- if (length(tv) > 0) tv[1] else NA_character_
      tv_lower <- tolower(tv_val)

      do.call(rbind, lapply(yrs, function(yr) {
        row <- all_coefs[all_coefs$location_id == loc & all_coefs$pheno_year == yr & tolower(all_coefs$Veg) == tv_lower, , drop = FALSE]
        pred <- NA_real_
        if (nrow(row) == 1) {
          if ("coef_median" %in% names(row) && is.finite(row$coef_median)) pred <- row$coef_median else pred <- row$coef
        }
        pred_abs <- pred
        rmse_val <- if (nrow(row) == 1 && "rmse" %in% names(row)) row$rmse else NA_real_
        # Prefer median bootstrap-based veg totals when available
        if ("coef_median" %in% names(all_coefs)) {
          sum_veg_coef <- sum(all_coefs$coef_median[all_coefs$location_id == loc & all_coefs$pheno_year == yr & tolower(all_coefs$Veg) != "barren"], na.rm = TRUE)
        } else {
          sum_veg_coef <- sum(all_coefs$coef[all_coefs$location_id == loc & all_coefs$pheno_year == yr & tolower(all_coefs$Veg) != "barren"], na.rm = TRUE)
        }
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

  timing_info$end_time <- Sys.time()
  total_time <- as.numeric(difftime(timing_info$end_time, timing_info$start_time, units = "secs"))

  cat(sprintf("\nTotal execution time: %.1f seconds (%.1f minutes)\n", total_time, total_time / 60))
  if (!is.null(timing_info$lib_construction_done) && !is.null(timing_info$moving_var_done)) {
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
    "Main processing: %.1f seconds\n",
    as.numeric(difftime(timing_info$end_time, timing_info$pca_computation_done, units = "secs"))
  ))

  cat("\nMESMA fitting completed successfully!\n")

# --- Modular helpers to re-run or run pieces of the main processing separately ---
process_batches <- function(BATCH_SIZE = 6, overwrite_wb = FALSE) { 
  cat("[INFO] process_batches: building batch list\n"); flush.console()
  target_locations <- unique(df_tasks$location_id)
  n_locs_to_process <- length(target_locations)
  loc_batches <- split(target_locations, ceiling(seq_along(target_locations) / BATCH_SIZE))
  results_local <- list()

  for (i in seq_along(loc_batches)) {
    batch_locs <- loc_batches[[i]]
    batch_df <- df_tasks[df_tasks$location_id %in% batch_locs, ]
    batch_location_list <- split(batch_df, batch_df$location_id)
    # Suppress any verbose output from per-location processing in helper
    sink_file <- tempfile()
    con_out <- tryCatch(file(sink_file, open = "wt"), error = function(e) NULL)
    con_msg <- tryCatch(file(paste0(sink_file, ".msg"), open = "wt"), error = function(e) NULL)
    tryCatch({
      if (!is.null(con_out) && inherits(con_out, "connection") && isOpen(con_out)) {
        sink(con_out, type = "output")
      }
      if (!is.null(con_msg) && inherits(con_msg, "connection") && isOpen(con_msg)) {
        sink(con_msg, type = "message")
      }
      batch_results <- .run_map(batch_location_list, fit_one_location, show_pb = FALSE)
    }, finally = {
      try(sink(type = "message"), silent = TRUE)
      try(sink(type = "output"), silent = TRUE)
      if (!is.null(con_out)) try(close(con_out), silent = TRUE)
      if (!is.null(con_msg)) try(close(con_msg), silent = TRUE)
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

      # Build location data (temporary, memory-efficient)
      loc_data <- do.call(rbind, lapply(loc_result, function(yr_res) yr_res$coef_df))
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

cat("\n=== Script execution finished ===\n")


