suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source("mesma_config.R")
source("ppi_helpers.R")
source("mesma_helpers.R")

DIAG_MAX_LOCATIONS <- 1500L
DIAG_APPLY_OUTLIERS <- FALSE

get_interpolate_method <- function(val) {
  if (is.null(val)) return(NULL)
  if (is.logical(val)) {
    if (isTRUE(val)) return("linear")
    return("none")
  }
  match.arg(tolower(as.character(val)), c("linear", "whittaker", "none"))
}

normalize_band_names <- function(df, bands = c("blue", "green", "red", "nir", "swir1", "swir2")) {
  if (is.null(df)) stop("normalize_band_names: df is NULL")
  if (nrow(df) == 0) return(df)
  current_names <- names(df)
  for (b in bands) {
    candidates <- c(b, toupper(b), tools::toTitleCase(b), paste0("band_", b), toupper(paste0("band_", b)), paste0("Band_", b))
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

remove_large_outliers <- function(df, candidates = NULL, mad_thresh = OUTLIER_MAD_THRESHOLD) {
  if (!isTRUE(ENABLE_OUTLIER_REMOVAL)) return(df)

  interp_method <- NULL
  if (exists("MESMA_PARAMS") && !is.null(MESMA_PARAMS$interpolate_inference)) {
    interp_method <- MESMA_PARAMS$interpolate_inference
  } else if (exists("INTERPOLATE_INFERENCE")) {
    interp_method <- INTERPOLATE_INFERENCE
  }
  interp_method <- get_interpolate_method(interp_method)

  if (is.null(candidates)) {
    candidates <- intersect(unique(c(OPTIMAL_INDICES, RAW_BANDS)), names(df))
  } else {
    candidates <- intersect(candidates, names(df))
  }
  if (length(candidates) == 0 || !"location_id" %in% names(df)) return(df)
  if (!"pheno_year" %in% names(df) && "date" %in% names(df)) {
    df$pheno_year <- assign_pheno_year(df$date)
  }

  grp <- interaction(df$location_id, ifelse(is.na(df$pheno_year), "NA", as.character(df$pheno_year)), drop = TRUE)
  groups <- split(seq_len(nrow(df)), grp)
  removed_idx <- logical(nrow(df))

  if (identical(interp_method, "whittaker")) {
    for (rows in groups) {
      sub <- df[rows, , drop = FALSE]
      if (length(rows) < 5) { removed_idx[rows] <- TRUE; next }
      for (col in candidates) {
        if (!is.numeric(sub[[col]])) next
        v <- sub[[col]]
        finite_idx <- is.finite(v)
        if (sum(finite_idx) < 3) next
        medv <- stats::median(v[finite_idx], na.rm = TRUE)
        madv <- stats::mad(v[finite_idx], na.rm = TRUE)
        if (!is.finite(madv) || madv <= 0) next
        mask <- rep(FALSE, length(v))
        mask[finite_idx] <- abs(v[finite_idx] - medv) > mad_thresh * madv
        removed_idx[rows[mask]] <- TRUE
      }
    }
    return(df[!removed_idx, , drop = FALSE])
  }

  for (rows in groups) {
    sub <- df[rows, , drop = FALSE]
    if (length(rows) < 5) { removed_idx[rows] <- TRUE; next }
    out_mask <- rep(FALSE, nrow(sub))
    has_date <- "date" %in% names(sub) && any(!is.na(sub$date))
    if (!has_date || length(rows) < 10) next
    sub$doy <- as.numeric(format(sub$date, "%j"))
    for (col in candidates) {
      if (!is.numeric(sub[[col]])) next
      colv <- sub[[col]]
      finite_idx <- is.finite(colv) & is.finite(sub$doy)
      if (sum(finite_idx) < 5) next
      tryCatch({
        x <- sub$doy[finite_idx]; y <- colv[finite_idx]
        n_unique <- length(unique(x))
        fit1 <- stats::smooth.spline(x, y, df = min(OUTLIER_SPLINE_MAX_DF, length(x) / 2, n_unique - 1))
        pred1 <- predict(fit1, x)$y
        res1 <- y - pred1
        mad1 <- stats::mad(res1, na.rm = TRUE)
        if (!is.finite(mad1) || mad1 <= 1e-6) return(NULL)
        keep_mask <- abs(res1 - stats::median(res1, na.rm = TRUE)) <= (mad_thresh * 1.5 * mad1)
        if (sum(keep_mask) >= 5) {
          n_unique2 <- length(unique(x[keep_mask]))
          fit2 <- stats::smooth.spline(x[keep_mask], y[keep_mask], df = min(OUTLIER_SPLINE_MAX_DF, sum(keep_mask) / 2, n_unique2 - 1))
          pred_final <- predict(fit2, x)$y
        } else {
          pred_final <- pred1
        }
        residuals <- y - pred_final
        med_res <- stats::median(residuals, na.rm = TRUE)
        mad_res <- stats::mad(residuals, na.rm = TRUE)
        if (!is.finite(mad_res) || mad_res <= 0) return(NULL)
        this_mask <- rep(FALSE, length(colv))
        this_mask[finite_idx] <- abs(residuals - med_res) > mad_thresh * mad_res
        out_mask <- out_mask | this_mask
      }, error = function(e) NULL)
    }
    if (any(out_mask, na.rm = TRUE)) removed_idx[rows[which(out_mask)]] <- TRUE
  }
  df[!removed_idx, , drop = FALSE]
}

prep_inference_like_fit <- function(path, norm_params) {
  df_inf <- readr::read_csv(path, show_col_types = FALSE)
  df_inf <- normalize_band_names(df_inf)
  if ("...1" %in% names(df_inf) && !"location_id" %in% names(df_inf)) {
    if (is.character(df_inf$...1) || is.numeric(df_inf$...1)) names(df_inf)[names(df_inf) == "...1"] <- "location_id"
  }

  if ("prediction_date" %in% names(df_inf)) {
    df_inf$date <- as.Date(df_inf$prediction_date)
  } else if ("date" %in% names(df_inf)) {
    df_inf$date <- as.Date(df_inf$date)
  }

  if (!"location_id" %in% names(df_inf)) {
    if (all(c("lon", "lat") %in% names(df_inf))) {
      df_inf$location_id <- make_location_id(df_inf$lon, df_inf$lat)
    } else {
      stop("Inference data has no location_id and no lon/lat columns.")
    }
  }
  df_inf$location_id <- as.character(df_inf$location_id)

  loc_ids <- unique(stats::na.omit(df_inf$location_id))
  if (length(loc_ids) > DIAG_MAX_LOCATIONS) {
    set.seed(42)
    keep_locs <- sample(loc_ids, DIAG_MAX_LOCATIONS)
    df_inf <- df_inf[df_inf$location_id %in% keep_locs, , drop = FALSE]
    cat(sprintf("[DIAG] Sampled %d/%d inference locations for diagnostic speed.\n", DIAG_MAX_LOCATIONS, length(loc_ids)))
  }

  df_inf <- apply_oli_etm_bias_correction(df_inf, dataset_label = "inference", log_prefix = "[DIAG BIAS CORR]")
  df_inf <- compute_indices_from_bands(df_inf)

  if ("NDDI" %in% names(df_inf)) {
    df_inf <- df_inf[!(df_inf$NDDI > NDDI_DUST_THRESHOLD), , drop = FALSE]
    if (isTRUE(DIAG_APPLY_OUTLIERS)) {
      df_inf <- remove_large_outliers(df_inf)
    } else {
      cat("[DIAG] Skipping spline outlier removal for diagnostic speed.\n")
    }
  }

  if (!"pheno_year" %in% names(df_inf)) df_inf$pheno_year <- assign_pheno_year(df_inf$date)
  if (!"doy" %in% names(df_inf)) df_inf$doy <- pheno_doy(df_inf$date)
  if (!"zenith.angle" %in% names(df_inf)) df_inf$zenith.angle <- NA_real_

  if (!"PPI" %in% names(df_inf) || all(!is.finite(df_inf$PPI))) {
    dvi_soil_vec_inf <- compute_dvi_soil_per_location(df_inf, quantile_p = 0.10)
    df_inf <- add_ppi_columns(df_inf, dvi_soil = dvi_soil_vec_inf)
  }

  df_inf <- apply_stored_normalization(df_inf, norm_params, cols = names(norm_params$INDEX_SCALES), lat_default = 40.2)
  df_inf <- backup_and_normalize_ppi(df_inf, label = "Inference")
  df_inf
}

cat("Loading training outputs...\n")
train_df <- readRDS("preprocessed_data.rds")
norm_params <- readRDS("training_norm_params.rds")
feature_cols <- intersect(names(norm_params$INDEX_SCALES), names(train_df))
cat(sprintf("Using %d normalized features: %s\n", length(feature_cols), paste(feature_cols, collapse = ", ")))

train_df$Veg <- tolower(trimws(as.character(train_df$Veg)))
train_use <- train_df %>%
  filter(Veg %in% tolower(ALLOWED_VEG)) %>%
  select(all_of(c("Veg", feature_cols, "satellite", "pheno_year")))

train_centroids <- train_use %>%
  group_by(Veg) %>%
  summarise(across(all_of(feature_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

cat("Preparing inference data with fit-time preprocessing...\n")
inf_df <- prep_inference_like_fit(INFERENCE_CSV, norm_params)
inf_use <- inf_df %>% select(any_of(c(feature_cols, "satellite", "pheno_year")))

# Feature-shift summary
feature_summary <- bind_rows(
  tibble(feature = feature_cols,
         dataset = "training",
         mean = sapply(feature_cols, function(x) mean(train_use[[x]], na.rm = TRUE)),
         sd = sapply(feature_cols, function(x) sd(train_use[[x]], na.rm = TRUE))),
  tibble(feature = feature_cols,
         dataset = "inference",
         mean = sapply(feature_cols, function(x) mean(inf_use[[x]], na.rm = TRUE)),
         sd = sapply(feature_cols, function(x) sd(inf_use[[x]], na.rm = TRUE)))
)
write.csv(feature_summary, "inference_feature_summary.csv", row.names = FALSE)

shift_tbl <- tibble(
  feature = feature_cols,
  train_mean = sapply(feature_cols, function(x) mean(train_use[[x]], na.rm = TRUE)),
  inf_mean = sapply(feature_cols, function(x) mean(inf_use[[x]], na.rm = TRUE)),
  mean_shift = inf_mean - train_mean,
  train_sd = sapply(feature_cols, function(x) sd(train_use[[x]], na.rm = TRUE)),
  inf_sd = sapply(feature_cols, function(x) sd(inf_use[[x]], na.rm = TRUE))
) %>% arrange(desc(abs(mean_shift)))
write.csv(shift_tbl, "inference_feature_shift.csv", row.names = FALSE)

# Nearest-centroid assignment
cent_mat <- as.matrix(train_centroids[, feature_cols, drop = FALSE])
rownames(cent_mat) <- train_centroids$Veg
x_mat <- as.matrix(inf_use[, feature_cols, drop = FALSE])
x_mat[!is.finite(x_mat)] <- 0
cent_mat[!is.finite(cent_mat)] <- 0

nearest <- apply(x_mat, 1, function(x) {
  d <- rowSums((t(t(cent_mat) - x))^2)
  names(which.min(d))
})
inf_nearest <- tibble(nearest_class = nearest, satellite = inf_df$satellite, pheno_year = inf_df$pheno_year)
write.csv(inf_nearest, "inference_nearest_centroid.csv", row.names = FALSE)

nearest_counts <- inf_nearest %>% count(nearest_class, sort = TRUE)
nearest_by_sat <- inf_nearest %>% count(satellite, nearest_class, sort = TRUE)
write.csv(nearest_counts, "inference_nearest_centroid_counts.csv", row.names = FALSE)
write.csv(nearest_by_sat, "inference_nearest_centroid_by_satellite.csv", row.names = FALSE)

cat("\n=== Nearest centroid counts ===\n")
print(nearest_counts)
cat("\n=== Largest mean shifts (inference - training) ===\n")
print(head(shift_tbl, 10))
cat("\nWrote:\n")
cat(" - inference_feature_summary.csv\n")
cat(" - inference_feature_shift.csv\n")
cat(" - inference_nearest_centroid.csv\n")
cat(" - inference_nearest_centroid_counts.csv\n")
cat(" - inference_nearest_centroid_by_satellite.csv\n")
