# fit_veg_spline.R
# Spline-based Spectral Mixture Analysis (Continuous Endmembers)
# Fits a smoothing spline to the full training dataset for each class (wrapping DOY)
# and performs NNLS unmixing against these continuous profiles.

# -----------------------------------------------------------------------------
# IMPORTS & SETUP
# -----------------------------------------------------------------------------
library(zoo)
library(dplyr)
library(sf)
library(ggplot2)
library(scales)
library(nnls)
library(terra)
library(magrittr)
library(mgcv) # For cyclic splines if needed, or we use smooth.spline with padding
library(lubridate)
library(tidyr) 

# --- CONFIGURATION ---
INPUT_CSV <- "C:\\Users\\yolan\\Downloads\\LS_S2_Harmonized_Timeseries.csv"
INFERENCE_CSV <- "C:\\Users\\yolan\\OneDrive\\Documenten\\UGENT\\Master\\masterproef\\GIS\\landsat_lower_inference.csv"
OUTPUT_DIR <- "phenology_results/veg_spline_fit"

TRAIN_YEARS <- 2024
ALLOWED_VEG <- c("populus", "tamarix", "phragmites") # Barren is added automatically
OPTIMAL_INDICES <- c("WDVI", "WDVI_BY_SWIR1", "GVI", "NIRv", "PSRI", "MSAVI2", "NDMI", "PPI", "EVI", "NDTI", "SATVI", "CIG", "BSI", "NBR", "TCW", "TCB", "NDSI")
FITTER_EXCLUDE_INDICES <- c("MSAVI", "NIRv", "OSAVI", "NDVI", "NDMI", "MSAVI2", "EVI", "SATVI")
RAW_BANDS <- c("blue", "green", "red", "nir", "swir1", "swir2")

ENABLE_OUTLIER_REMOVAL <- TRUE
OUTLIER_MAD_THRESHOLD <- 3.5
MIN_ENDMEMBER_SAMPLES <- 5

PARALLEL_ENABLE <- TRUE
PARALLEL_WORKERS <- 4

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------------------------------

# Source external helpers if available
if (file.exists("mesma_helpers.R")) source("mesma_helpers.R")
if (file.exists("ppi_helpers.R")) source("ppi_helpers.R")
if (file.exists("init_mesma.R")) try(source("init_mesma.R"), silent = TRUE)

# Define core helpers
assign_pheno_year <- function(d) {
  d <- as.Date(d)
  ifelse(is.na(d), NA_integer_, ifelse(lubridate::month(d) >= 3, lubridate::year(d), lubridate::year(d) - 1))
}

pheno_doy <- function(d) {
  d <- tryCatch(as.Date(d), error = function(e) NA)
  month <- lubridate::month(d)
  ifelse(is.na(d), NA_integer_,
    ifelse(month >= 3,
      as.integer(d - as.Date(paste0(lubridate::year(d), "-03-01"))) + 1L,
      as.integer(d - as.Date(paste0(lubridate::year(d) - 1, "-03-01"))) + 1L
    )
  )
}

normalize_band_names <- function(df) {
  bands <- RAW_BANDS
  if (is.null(df) || nrow(df) == 0) return(df)
  current_names <- names(df)
  for (b in bands) {
    candidates <- c(b, toupper(b), tools::toTitleCase(b), paste0('band_', b))
    for (cand in candidates) {
      if (cand %in% current_names && !(b %in% current_names)) {
        names(df)[names(df) == cand] <- b
        break
      }
    }
  }
  df
}

# --- Index Computation Helpers (Ported from fit_veg_mixture_mesma.R) ---

if (!exists("SOIL_LINE_SLOPE")) SOIL_LINE_SLOPE <- 1.0

compute_soil_line_slope <- function(input_df, min_samples = MIN_ENDMEMBER_SAMPLES, assign_global_dvi = TRUE) {
  if (!all(c('nir','red','Veg') %in% names(input_df))) {
    cat("[SOIL LINE] Required columns 'nir','red' or 'Veg' missing; using default SOIL_LINE_SLOPE=1.0\n")
    return(invisible(NA_real_))
  }
  bare_soil_df <- input_df[tolower(input_df$Veg) == 'barren' & is.finite(input_df$nir) & is.finite(input_df$red), , drop = FALSE]
  if (nrow(bare_soil_df) > min_samples) {
    soil_line_model <- tryCatch(lm(nir ~ red, data = bare_soil_df), error = function(e) NULL)
    if (!is.null(soil_line_model)) {
      slope <- as.numeric(coef(soil_line_model)[2])
      assign("SOIL_LINE_SLOPE", slope, envir = globalenv())
      cat(sprintf("[SOIL LINE] Calculated SOIL_LINE_SLOPE=%.4f from %d bare soil pixels\n", slope, nrow(bare_soil_df)))
      if (assign_global_dvi) {
        dvi_soil_calc <- mean(bare_soil_df$nir - bare_soil_df$red, na.rm = TRUE)
        if (is.finite(dvi_soil_calc)) {
          assign("GLOBAL_TRAINING_DVI_SOIL", dvi_soil_calc, envir = globalenv())
          cat(sprintf("[SOIL LINE] Set GLOBAL_TRAINING_DVI_SOIL=%.6f based on bare soil DVI\n", dvi_soil_calc))
        }
      }
      return(invisible(slope))
    }
  }
  cat("[SOIL LINE] Using default SOIL_LINE_SLOPE=1.0 (insufficient soil samples)\n")
  assign("SOIL_LINE_SLOPE", 1.0, envir = globalenv())
  return(invisible(1.0))
}

compute_indices_from_bands <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  eps <- 1e-9
  sl_slope <- if (exists("SOIL_LINE_SLOPE")) get("SOIL_LINE_SLOPE") else 1.0

  if (all(c('nir','red') %in% names(df))) df$DVI <- as.numeric(df$nir) - as.numeric(df$red)
  if (all(c('nir','red') %in% names(df))) df$WDVI <- as.numeric(df$nir) - sl_slope * as.numeric(df$red)
  if (all(c('nir','red') %in% names(df))) df$OSAVI <- (as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + as.numeric(df$red) + 0.16)
  if (all(c('red','green','blue') %in% names(df))) df$MCARI <- ((as.numeric(df$red) - as.numeric(df$green)) - 0.2*(as.numeric(df$red) - as.numeric(df$blue))) * (as.numeric(df$red) / (as.numeric(df$green) + eps))
  if (all(c('green','red') %in% names(df))) df$PRI <- (as.numeric(df$green) - as.numeric(df$red)) / (as.numeric(df$green) + as.numeric(df$red) + eps)
  if (all(c('nir','red') %in% names(df))) df$NIRv <- as.numeric(df$nir) * ((as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + as.numeric(df$red) + eps))
  if (all(c('red','blue','nir') %in% names(df))) df$PSRI <- (as.numeric(df$red) - as.numeric(df$blue)) / (as.numeric(df$nir) + eps)
  if (all(c('nir','swir2') %in% names(df))) df$NBR <- (as.numeric(df$nir) - as.numeric(df$swir2)) / (as.numeric(df$nir) + as.numeric(df$swir2) + eps)

  if (all(c('green','red','nir','swir1','swir2','blue') %in% names(df))) {
    df$TCB <- 0.3029 * as.numeric(df$blue) + 0.2786 * as.numeric(df$green) + 0.4733 * as.numeric(df$red) + 0.5599 * as.numeric(df$nir) + 0.508 * as.numeric(df$swir1) + 0.1872 * as.numeric(df$swir2)
    df$TCG <- -0.2941 * as.numeric(df$blue) - 0.243 * as.numeric(df$green) - 0.5424 * as.numeric(df$red) + 0.7276 * as.numeric(df$nir) + 0.0713 * as.numeric(df$swir1) - 0.1608 * as.numeric(df$swir2)
    df$TCW <- 0.1511 * as.numeric(df$blue) + 0.1973 * as.numeric(df$green) + 0.3283 * as.numeric(df$red) + 0.3407 * as.numeric(df$nir) - 0.7117 * as.numeric(df$swir1) - 0.4559 * as.numeric(df$swir2)
    df$GVI <- df$TCG
  }

  if (all(c('nir','red') %in% names(df))) df$MSAVI2 <- (2 * as.numeric(df$nir) + 1 - sqrt(pmax(0, (2 * as.numeric(df$nir) + 1)^2 - 8 * (as.numeric(df$nir) - as.numeric(df$red))))) / 2
  if (all(c('nir','swir1') %in% names(df))) df$NDMI <- (as.numeric(df$nir) - as.numeric(df$swir1)) / (as.numeric(df$nir) + as.numeric(df$swir1) + eps)
  if ('WDVI' %in% names(df) && 'swir1' %in% names(df)) df$WDVI_BY_SWIR1 <- as.numeric(df$WDVI) / (as.numeric(df$swir1) + eps)
  if (all(c('nir','red','blue') %in% names(df))) df$EVI <- 2.5 * ((as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + 6 * as.numeric(df$red) - 7.5 * as.numeric(df$blue) + 1 + eps))
  if (all(c('swir1','swir2') %in% names(df))) df$NDTI <- (as.numeric(df$swir1) - as.numeric(df$swir2)) / (as.numeric(df$swir1) + as.numeric(df$swir2) + eps)
  if (all(c('green','swir1') %in% names(df))) df$NDSI <- (as.numeric(df$green) - as.numeric(df$swir1)) / (as.numeric(df$green) + as.numeric(df$swir1) + eps)
  if (all(c('swir1','red','swir2') %in% names(df))) df$SATVI <- ((as.numeric(df$swir1) - as.numeric(df$red)) / (as.numeric(df$swir1) + as.numeric(df$red) + 0.5 + eps)) * 1.5 - (as.numeric(df$swir2) / 2)
  if (all(c('nir','green') %in% names(df))) df$CIG <- (as.numeric(df$nir) / (as.numeric(df$green) + eps)) - 1
  if (all(c('swir1','red','nir','blue') %in% names(df))) {
    term1 <- as.numeric(df$swir1) + as.numeric(df$red)
    term2 <- as.numeric(df$nir) + as.numeric(df$blue)
    df$BSI <- (term1 - term2) / (term1 + term2 + eps)
  }
  if ('NIRv' %in% names(df)) df$NIRv <- df$NIRv * 1.3

  if (all(c('nir','red') %in% names(df)) && !"NDVI" %in% names(df)) {
     df$NDVI <- (as.numeric(df$nir) - as.numeric(df$red)) / (as.numeric(df$nir) + as.numeric(df$red) + eps)
  }

  df
}

# --- Solvers & Normalization ---

solve_weights_fcls <- function(E, y) {
  if (is.null(E) || ncol(E) < 1) return(rep(0, 0))
  lambda <- 10 * mean(abs(y), na.rm = TRUE)
  if (lambda == 0) lambda <- 10
  E_aug <- rbind(E, rep(lambda, ncol(E)))
  y_aug <- c(y, lambda)
  res <- nnls::nnls(E_aug, y_aug)
  w <- res$x
  if (sum(w) > 0) w <- w / sum(w) else w <- rep(1/ncol(E), ncol(E))
  return(w)
}

remove_large_outliers <- function(df, candidates = OPTIMAL_INDICES, mad_thresh = OUTLIER_MAD_THRESHOLD) {
  if (!ENABLE_OUTLIER_REMOVAL) return(df)
  cat("[OUTLIER] Removing outliers using MAD threshold...\n")
  keep_mask <- rep(TRUE, nrow(df))
  candidates <- intersect(candidates, names(df))
  for (col in candidates) {
    vals <- df[[col]]
    if (all(is.na(vals))) next
    med <- median(vals, na.rm = TRUE)
    m <- mad(vals, na.rm = TRUE)
    if (is.na(m) || m == 0) next
    outliers <- abs(vals - med) > (mad_thresh * m)
    keep_mask[outliers & !is.na(outliers)] <- FALSE
  }
  n_removed <- sum(!keep_mask)
  cat(sprintf("[OUTLIER] Removed %d observations.\n", n_removed))
  df[keep_mask, , drop = FALSE]
}

normalize_data <- function(df, cols, params = NULL) {
  cols <- intersect(cols, names(df))
  if (is.null(params)) {
    params <- list(means = numeric(), sds = numeric())
    for (col in cols) {
      vals <- df[[col]]
      params$means[col] <- mean(vals, na.rm = TRUE)
      params$sds[col] <- sd(vals, na.rm = TRUE)
    }
  }
  for (col in cols) {
    m <- params$means[col]
    s <- params$sds[col]
    if (is.na(s) || s == 0) s <- 1
    df[[col]] <- (df[[col]] - m) / s
  }
  return(list(df = df, params = params))
}

# --- VALIDATION REPORTING ---
report_validation_accuracy <- function(val_coefs, label = "") {
    if (is.null(val_coefs) || nrow(val_coefs) == 0) {
      cat("[NOTICE] No validation coefficients found.")
      return(invisible(NULL))
    }
    prefix <- if (nzchar(label)) paste0(" (", label, ")") else ""
    labels_df <- NULL
    if ("Veg" %in% names(val_coefs)) {
       labels_df <- val_coefs %>% dplyr::select(location_id, Veg) %>% dplyr::distinct()
    } else if (file.exists("validation_locations.csv")) {
       vloc <- tryCatch(read.csv("validation_locations.csv"), error = function(e) NULL)
       if (!is.null(vloc) && "location_id" %in% names(vloc) && "Veg" %in% names(vloc)) {
          labels_df <- vloc %>% dplyr::select(location_id, Veg) %>% dplyr::distinct()
       }
    }
    if (is.null(labels_df)) {
       cat("[NOTICE] Could not find ground truth labels (Veg column). Cannot compute accuracy.")
       return(invisible(NULL))
    }
    labels_df$Veg <- tolower(trimws(as.character(labels_df$Veg)))
    labels_df <- labels_df %>% dplyr::rename(true_veg = Veg)
    if (!"true_veg" %in% names(val_coefs)) {
        val_coefs <- val_coefs %>% dplyr::left_join(labels_df, by = "location_id")
    }
    cat(sprintf("\n=== VALIDATION ACCURACY ON HELD-OUT SET%s ===\n", prefix))
    cat(sprintf("Validation set: %d locations\n", length(unique(val_coefs$location_id))))
    
    meta_cols <- c("location_id", "pheno_year", "doy", "date", "Veg", "true_veg", "RMSE")
    class_cols <- setdiff(names(val_coefs), c(meta_cols, names(labels_df), OPTIMAL_INDICES, "NDVI", "MSAVI2", "PPI"))
    class_cols <- class_cols[sapply(val_coefs[class_cols], is.numeric)]
    
    val_coefs_wide <- val_coefs
    for (cc in class_cols) {
        val_coefs_wide[[paste0("frac_", cc)]] <- val_coefs_wide[[cc]]
    }
    all_veg_cols <- paste0("frac_", class_cols)
    
    if (length(all_veg_cols) > 0) {
      cat("\n--- CONFUSION MATRIX (Winner-Takes-All) ---")
      frac_mat <- as.matrix(val_coefs_wide[, all_veg_cols])
      preds <- apply(frac_mat, 1, function(r) {
         idx <- which.max(r)
         if(length(idx) == 0) return(NA_character_)
         sub("^frac_", "", all_veg_cols[idx])
      })
      val_coefs_wide$predicted_class <- preds
      all_classes <- sort(unique(c(tolower(val_coefs_wide$true_veg), val_coefs_wide$predicted_class)))
      all_classes <- all_classes[!is.na(all_classes)]
      conf_mat <- table(
         True = factor(tolower(val_coefs_wide$true_veg), levels = all_classes),
         Predicted = factor(val_coefs_wide$predicted_class, levels = all_classes)
      )
      cat("\n[Validation Confusion Matrix - Counts (excluding barren)]\n")
      veg_cols <- setdiff(colnames(conf_mat), "barren")
      veg_rows <- setdiff(rownames(conf_mat), "barren")
      if (length(veg_cols) > 0 && length(veg_rows) > 0) print(conf_mat[veg_rows, veg_cols]) else print(conf_mat)
      
      conf_mat_pct <- sweep(conf_mat, 1, rowSums(conf_mat), "/")
      conf_mat_pct[is.nan(conf_mat_pct)] <- 0
      cat("\n[Validation Confusion Matrix - Row Normalized (Recall, excluding barren)]\n")
      if (length(veg_cols) > 0 && length(veg_rows) > 0) print(round(conf_mat_pct[veg_rows, veg_cols], 2)) else print(round(conf_mat_pct, 2))
      
      oa <- sum(diag(conf_mat)) / sum(conf_mat)
      cat(sprintf("\nOverall Validation Accuracy: %.2f%%\n", oa * 100))
    }
    cat("Validation accuracy computed", prefix, ".\n", sep = "")
}

# -----------------------------------------------------------------------------
# CORE LOGIC
# -----------------------------------------------------------------------------

build_spline_library <- function(df_train, indices, allowed_veg) {
  cat("\n=== Building Spline Library (One Curve Per Class) ===\n")
  classes <- unique(c("barren", allowed_veg))
  library_splines <- list()
  plot_data_list <- list()
  
  for (veg in classes) {
    cat(sprintf("Fitting splines for class: %s\n", veg))
    sub <- df_train[tolower(df_train$Veg) == veg, ]
    if (nrow(sub) < 10) {
      cat(sprintf("  [WARNING] Not enough samples (%d) for %s. Skipping.\n", nrow(sub), veg))
      next
    }
    class_splines <- list()
    for (idx in indices) {
      if (!idx %in% names(sub)) next
      d <- sub$doy
      v <- sub[[idx]]
      mask <- is.finite(d) & is.finite(v)
      d <- d[mask]
      v <- v[mask]
      if (length(d) < 5) {
        class_splines[[idx]] <- function(x) rep(mean(v, na.rm=TRUE), length(x))
        next
      }
      # Circular pad
      d_aug <- c(d - 365, d, d + 365)
      v_aug <- c(v, v, v)
      fit <- tryCatch({
        smooth.spline(x = d_aug, y = v_aug, spar = 0.7)
      }, error = function(e) {
        warning(paste("Spline fit failed for", veg, idx, ":", e$message))
        NULL
      })
      if (is.null(fit)) {
         class_splines[[idx]] <- function(x) rep(mean(v, na.rm=TRUE), length(x))
         next
      }
      spline_fun <- function(new_doy, fit_obj = fit) {
        predict(fit_obj, new_doy)$y
      }
      class_splines[[idx]] <- spline_fun
      pred_doy <- 1:365
      pred_val <- spline_fun(pred_doy)
      plot_data_list[[paste(veg, idx, sep="_")]] <- data.frame(
        Veg = veg, Index = idx, DOY = pred_doy, Value = pred_val
      )
    }
    library_splines[[veg]] <- class_splines
  }
  
  all_plots <- do.call(rbind, plot_data_list)
  if (!is.null(all_plots)) {
    p <- ggplot(all_plots, aes(x = DOY, y = Value, color = Veg)) +
      geom_line(size = 1) +
      facet_wrap(~Index, scales = "free_y") +
      theme_minimal() +
      labs(title = "Class Phenology Splines")
    ggsave(file.path(OUTPUT_DIR, "library_splines.png"), p, width = 12, height = 8)
  }
  return(library_splines)
}

fit_one_location_spline <- function(task_df, spline_lib, indices) {
  n_obs <- nrow(task_df)
  classes <- names(spline_lib)
  n_classes <- length(classes)
  coefs <- matrix(NA, nrow = n_obs, ncol = n_classes)
  colnames(coefs) <- classes
  rmses <- numeric(n_obs)
  
  for (i in seq_len(n_obs)) {
    doy <- task_df$doy[i]
    if (is.na(doy)) next
    E <- matrix(0, nrow = length(indices), ncol = n_classes)
    y <- numeric(length(indices))
    valid_idx_mask <- rep(TRUE, length(indices))
    
    for (k in seq_along(indices)) {
      idx <- indices[k]
      val <- task_df[[idx]][i]
      if (is.na(val)) {
        valid_idx_mask[k] <- FALSE
        next
      }
      y[k] <- val
      for (c_idx in seq_along(classes)) {
        cls <- classes[c_idx]
        fun <- spline_lib[[cls]][[idx]]
        if (!is.null(fun)) {
          E[k, c_idx] <- fun(doy)
        } else {
          valid_idx_mask[k] <- FALSE
        }
      }
    }
    if (sum(valid_idx_mask) < 2) next
    E_curr <- E[valid_idx_mask, , drop = FALSE]
    y_curr <- y[valid_idx_mask]
    w <- solve_weights_fcls(E_curr, y_curr)
    coefs[i, ] <- w
    pred <- E_curr %*% w
    rmses[i] <- sqrt(mean((y_curr - pred)^2))
  }
  res <- cbind(task_df, as.data.frame(coefs), RMSE = rmses)
  return(res)
}

# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

cat("Loading data...\n")
if (!file.exists(INPUT_CSV)) stop("Input CSV not found!")

df <- read.csv(INPUT_CSV)
df <- normalize_band_names(df)
df$date <- as.Date(if ("prediction_date" %in% names(df)) df$prediction_date else df$date)
df$pheno_year <- assign_pheno_year(df$date)
df$doy <- pheno_doy(df$date)
df <- df[!is.na(df$doy), ]

# Compute Indices & Soil Line
compute_soil_line_slope(df)
df <- compute_indices_from_bands(df)

# Check for PPI or add it
if (!"PPI" %in% names(df) && exists("auto_add_ppi_columns")) {
  cat("Calculating PPI...\n")
  ppi_df <- tryCatch(auto_add_ppi_columns(df), error = function(e) {
     cat(sprintf("[WARNING] auto_add_ppi_columns failed: %s. Using df without new PPI.\n", e$message))
     return(NULL)
  })
  if (!is.null(ppi_df)) df <- ppi_df
}

# Prepare Training Data
cat("Preparing Training Data...\n")
df_train_full <- df[df$pheno_year %in% TRAIN_YEARS, ]
if (nrow(df_train_full) == 0) stop(sprintf("No training data found for years: %s", paste(TRAIN_YEARS, collapse=", ")))

if ("Veg" %in% names(df_train_full)) {
    df_train_full$Veg <- tolower(trimws(df_train_full$Veg))
    df_train_full <- df_train_full[df_train_full$Veg %in% c("barren", ALLOWED_VEG), ]
} else {
    stop("Input CSV must contain a 'Veg' column for training!")
}

df_train_full <- remove_large_outliers(df_train_full)
if (nrow(df_train_full) == 0) stop("All training data removed by filtering!")

# Split Train/Test
cat("[SPLIT] Performing stratified 80/20 split...\n")
set.seed(42)
loc_veg_summary <- df_train_full %>%
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
  n_test <- ceiling(length(v_locs) * 0.20)
  if (n_test > 0) {
    selected <- sample(v_locs, n_test)
    test_locs_list[[i]] <- data.frame(location_id = selected, Veg = v, stringsAsFactors = FALSE)
  }
}
test_locs_df <- do.call(rbind, test_locs_list)

if (!is.null(test_locs_df) && nrow(test_locs_df) > 0) {
  write.csv(test_locs_df, file.path(OUTPUT_DIR, "validation_locations.csv"), row.names = FALSE)
  df_validation <- df_train_full[df_train_full$location_id %in% test_locs_df$location_id, ]
  df_train <- df_train_full[!df_train_full$location_id %in% test_locs_df$location_id, ]
  cat(sprintf("[SPLIT] Train: %d rows (%d locs), Val: %d rows (%d locs)\n", 
              nrow(df_train), length(unique(df_train$location_id)), 
              nrow(df_validation), length(unique(df_validation$location_id))))
} else {
  df_train <- df_train_full
  df_validation <- NULL
  cat("[SPLIT] Warning: No test locations selected.\n")
}

# Normalize
cat("Normalizing Training Data...\n")
avail_indices <- intersect(OPTIMAL_INDICES, names(df_train))
avail_indices <- setdiff(avail_indices, FITTER_EXCLUDE_INDICES)
if(length(avail_indices) < 2) stop("Not enough indices available!")

norm_res <- normalize_data(df_train, avail_indices)
df_train <- norm_res$df
norm_params <- norm_res$params

# Build Library
spline_lib <- build_spline_library(df_train, avail_indices, ALLOWED_VEG)

# Validation Unmixing
if (!is.null(df_validation) && nrow(df_validation) > 0) {
  cat("\n=== RUNNING VALIDATION ===\n")
  norm_val <- normalize_data(df_validation, avail_indices, params = norm_params)
  df_validation <- norm_val$df
  val_locs <- unique(df_validation$location_id)
  
  val_results_list <- lapply(val_locs, function(lid) {
     sub <- df_validation[df_validation$location_id == lid, ]
     fit_one_location_spline(sub, spline_lib, avail_indices)
  })
  val_results <- do.call(rbind, val_results_list)
  report_validation_accuracy(val_results, label = "Held-Out Set")
  write.csv(val_results, file.path(OUTPUT_DIR, "validation_results.csv"), row.names = FALSE)
}

# Inference
cat("\n=== RUNNING INFERENCE ===\n")
if (file.exists(INFERENCE_CSV)) {
  cat("Loading Inference CSV...\n")
  df_inf <- read.csv(INFERENCE_CSV)
  df_inf <- normalize_band_names(df_inf)
  df_inf$date <- as.Date(if ("prediction_date" %in% names(df_inf)) df_inf$prediction_date else df_inf$date)
  df_inf$pheno_year <- assign_pheno_year(df_inf$date)
  df_inf$doy <- pheno_doy(df_inf$date)
  df_inf <- compute_indices_from_bands(df_inf)
  if (!"PPI" %in% names(df_inf) && exists("auto_add_ppi_columns")) {
     ppi_inf <- tryCatch(auto_add_ppi_columns(df_inf), error = function(e) NULL)
     if (!is.null(ppi_inf)) df_inf <- ppi_inf
  }
  norm_inf <- normalize_data(df_inf, avail_indices, params = norm_params)
  df_inf <- norm_inf$df
} else {
  cat("Using all data for inference...\n")
  df_inf <- df
  norm_inf <- normalize_data(df_inf, avail_indices, params = norm_params)
  df_inf <- norm_inf$df
}

cat(sprintf("Starting Unmixing on %d locations...\n", length(unique(df_inf$location_id))))
location_ids <- unique(df_inf$location_id)

process_loc <- function(lid) {
  sub <- df_inf[df_inf$location_id == lid, ]
  fit_one_location_spline(sub, spline_lib, avail_indices)
}

if (PARALLEL_ENABLE) {
  library(future)
  library(future.apply)
  plan(multisession, workers = PARALLEL_WORKERS)
  results_list <- future_lapply(location_ids, process_loc)
} else {
  results_list <- lapply(location_ids, process_loc)
}

results_df <- do.call(rbind, results_list)
write.csv(results_df, file.path(OUTPUT_DIR, "spline_unmixing_results.csv"), row.names = FALSE)

agg_res <- results_df %>%
  group_by(pheno_year) %>%
  summarise(across(all_of(names(spline_lib)), mean, na.rm=TRUE))
print(agg_res)
write.csv(agg_res, file.path(OUTPUT_DIR, "global_annual_means.csv"), row.names = FALSE)

cat("Done.\n")
