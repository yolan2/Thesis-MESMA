# ==============================================================================
# HLS Visualization Script (CSV -> Plots)
# Author: Phenology Analysis Toolkit
# Date: 2025
# Purpose: Visualize per-location time series and trends from extracted CSV.
#          Includes per-location year-scale graph (each year different color)
#          and time-based outlier detection.
# ==============================================================================

OUTPUT_DIR <- "phenology_results"
INPUT_CSV <- file.path(OUTPUT_DIR, "hls_phenology_data.csv")
PLOTS_DIR <- file.path(OUTPUT_DIR, "plots")
if(!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)
MAX_LOC_PER_INDEX <- 30
# Limit of locations per vegetation-type plot (to keep plots readable)
MAX_LOC_PER_VEGTYPE <- 12

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(lubridate)
  library(mgcv)
  library(rlang)
})

# Helper: phenological year (March 1-based)
assign_pheno_year <- function(d) {
  d <- as.Date(d)
  ifelse(is.na(d), NA_integer_, ifelse(lubridate::month(d) >= 3, lubridate::year(d), lubridate::year(d) - 1))
}

# Stabilization constants for index calculations
EPS_IDX_DENOM <- 1e-6

# Silence NSE notes for R CMD check / linters
if(getRversion() >= "2.15.1") utils::globalVariables(c(
  "outlier_time", "doy", "lwr", "upr", "fit", "location_id", "frac_out", "coverage", "year", "Veg",
  "date", "ms_peak", "FVC_est", ".tmp_year"
))

cat("=== HLS Visualization ===\n")
if(!file.exists(INPUT_CSV)) stop("Input CSV not found: ", INPUT_CSV)

df <- suppressMessages(readr::read_csv(INPUT_CSV, show_col_types = FALSE))
if(nrow(df) == 0) stop("No rows in input CSV")

# If Veg/coverage are missing, attach from GeoJSON by rounded coordinates (consistent with fitter)
if(!("Veg" %in% names(df)) || all(is.na(df$Veg))) {
  if(requireNamespace("sf", quietly = TRUE)) {
    pp <- try(sf::st_read("C:/MAP/input/pure_pixels.geojson", quiet = TRUE), silent = TRUE)
    if(!inherits(pp, "try-error")) {
      pp <- pp[!sf::st_is_empty(pp), , drop = FALSE]
      pp <- sf::st_transform(pp, 4326)
      coords <- sf::st_coordinates(pp)
      pp$lon_round <- round(coords[,1], 4)
      pp$lat_round <- round(coords[,2], 4)
      if(!("location_id" %in% names(df))) {
        df <- df %>% dplyr::mutate(location_id = paste0(round(target_lat,4), "_", round(target_lon,4)))
      }
      if(all(c("target_lon","target_lat") %in% names(df))) {
        df$lon_round <- round(df$target_lon, 4)
        df$lat_round <- round(df$target_lat, 4)
        keep <- intersect(c("Veg","coverage","lon_round","lat_round"), names(pp))
        j <- dplyr::left_join(df, as.data.frame(pp)[, keep, drop = FALSE], by = c("lat_round","lon_round"))
        if(!("Veg" %in% names(df)) && ("Veg" %in% names(j))) df$Veg <- j$Veg
        if(!("coverage" %in% names(df)) && ("coverage" %in% names(j))) df$coverage <- j$coverage
        df <- df %>% dplyr::select(-dplyr::any_of(c("lat_round","lon_round")))
        cat("Attached Veg/coverage from GeoJSON (rounded coord join).\n")
      }
    }
  }
}

# Compute additional vegetation indices if raw bands are available
if(all(c("red_value", "green_value", "blue_value", "nir_value") %in% names(df))) {
  clamp_idx <- function(x, lo, hi) {
    x[!is.finite(x)] <- NA_real_
    x <- pmax(lo, pmin(hi, x))
    x
  }
  if(!"NDVI" %in% names(df)) df$NDVI <- (df$nir_value - df$red_value) / pmax(df$nir_value + df$red_value, EPS_IDX_DENOM)
  if(!"EVI"  %in% names(df)) df$EVI  <- 2.5 * (df$nir_value - df$red_value) / pmax(df$nir_value + 6 * df$red_value - 7.5 * df$blue_value + 1, EPS_IDX_DENOM)
  if(!"SAVI" %in% names(df)) df$SAVI <- ((df$nir_value - df$red_value) / pmax(df$nir_value + df$red_value + 0.5, EPS_IDX_DENOM)) * 1.5
  if(!"MSAVI"%in% names(df)) df$MSAVI<- (2 * df$nir_value + 1 - sqrt(pmax((2 * df$nir_value + 1)^2 - 8 * (df$nir_value - df$red_value), 0))) / 2
  if(!"GNDVI"%in% names(df)) df$GNDVI<- (df$nir_value - df$green_value) / pmax(df$nir_value + df$green_value, EPS_IDX_DENOM)
  if(!"ARI"  %in% names(df)) df$ARI  <- (1 / pmax(df$green_value, EPS_IDX_DENOM)) - (1 / pmax(df$red_value, EPS_IDX_DENOM))
  if(!"CRI"  %in% names(df)) df$CRI  <- (1 / pmax(df$blue_value, EPS_IDX_DENOM)) - (1 / pmax(df$green_value, EPS_IDX_DENOM))
  if(!"RGRI" %in% names(df)) df$RGRI <- df$red_value / pmax(df$green_value, EPS_IDX_DENOM)
  if(!"BGI"  %in% names(df)) df$BGI  <- df$blue_value / pmax(df$green_value, EPS_IDX_DENOM)
  # Clamp to reasonable ranges to avoid extreme outliers in plots
  df$NDVI <- clamp_idx(df$NDVI, -1, 1)
  df$EVI  <- clamp_idx(df$EVI,  -1, 5)
  df$SAVI <- clamp_idx(df$SAVI, -1, 1)
  df$MSAVI<- clamp_idx(df$MSAVI,-1, 1)
  df$GNDVI<- clamp_idx(df$GNDVI,-1, 1)
  df$ARI  <- clamp_idx(df$ARI,  -100, 100)
  df$CRI  <- clamp_idx(df$CRI,  -100, 100)
  df$RGRI <- clamp_idx(df$RGRI,  0, 10)
  df$BGI  <- clamp_idx(df$BGI,   0, 10)
}

# Estimate MSAVI soil/veg and compute FVC per location-year for unmixing
get_msavi_params <- function(df) {
  MSAVI_SOIL <- 0.05; MSAVI_VEG <- 0.54
  if("MSAVI" %in% names(df)) {
    barren_vals <- numeric(0)
    if("barren_location" %in% names(df)) {
      bl <- df$barren_location
      if(is.logical(bl)) sel <- bl else if(is.numeric(bl)) sel <- is.finite(bl) & bl > 0 else sel <- grepl("true|barren|bare", tolower(as.character(bl)))
      barren_vals <- df$MSAVI[sel & is.finite(df$MSAVI)]
    }
    if(length(barren_vals) == 0 && "Veg" %in% names(df)) {
      sel <- grepl("barr|bare", tolower(df$Veg))
      barren_vals <- df$MSAVI[sel & is.finite(df$MSAVI)]
    }
    if(length(barren_vals) >= 1) MSAVI_SOIL <- mean(barren_vals, na.rm = TRUE)
    ms <- df$MSAVI; ms <- ms[is.finite(ms)]
    if(length(ms) >= 30) { q98 <- as.numeric(stats::quantile(ms, 0.98, na.rm = TRUE, type = 7)); if(is.finite(q98)) MSAVI_VEG <- max(0.6, q98) } else if(length(ms) > 0) {
      mx <- max(ms, na.rm = TRUE); if(is.finite(mx)) MSAVI_VEG <- max(0.6, mx)
    }
  }
  if(!(is.finite(MSAVI_SOIL) && is.finite(MSAVI_VEG))) { MSAVI_SOIL <- 0.1; MSAVI_VEG <- 0.8 }
  if(MSAVI_VEG <= MSAVI_SOIL) MSAVI_VEG <- MSAVI_SOIL + 0.2
  list(soil = MSAVI_SOIL, veg = MSAVI_VEG)
}
params <- get_msavi_params(df)
MSAVI_SOIL <- params$soil; MSAVI_VEG <- params$veg
cat(sprintf("Using MSAVI_SOIL=%.3f, MSAVI_VEG=%.3f (visualize)\n", MSAVI_SOIL, MSAVI_VEG))
if("MSAVI" %in% names(df) && "date" %in% names(df)) {
  df <- df %>% dplyr::mutate(.tmp_year = assign_pheno_year(date))
  fvc_by <- df %>% dplyr::group_by(location_id, .tmp_year) %>% dplyr::summarise(ms_peak = suppressWarnings(max(MSAVI, na.rm = TRUE)), .groups='drop') %>%
    dplyr::mutate(FVC_est = {denom <- (MSAVI_VEG - MSAVI_SOIL); ifelse(is.finite(ms_peak) & denom > 1e-6, pmax(0, pmin(1, ((ms_peak - MSAVI_SOIL)/denom)^2)), NA_real_)})
  df <- df %>% dplyr::left_join(fvc_by, by = c("location_id",".tmp_year"))
  # Avoid duplicate year column: rename only if 'year' not already present
  if("year" %in% names(df)) {
    df <- df %>% dplyr::select(-.tmp_year)
  } else {
    df <- df %>% dplyr::rename(year = .tmp_year)
  }
} else {
  if(!"year" %in% names(df) && "date" %in% names(df)) df <- df %>% dplyr::mutate(year = lubridate::year(date))
  df$FVC_est <- NA_real_
}

# Collect locations dropped because they are mostly outliers
OUTLIER_LOCATIONS <- character(0)

# Choose index to visualize
choose_index <- function(data, preferred=c("NDVI","SAVI","EVI","MSAVI","GNDVI","ARI","CRI","RGRI","BGI","NDWI","NDWI_GREEN")) {
  for(cn in preferred) if(cn %in% names(data) && sum(is.finite(data[[cn]]), na.rm = TRUE) > 3) return(cn)
  if("nir_value" %in% names(data) && sum(is.finite(data$nir_value), na.rm = TRUE) > 3) return("nir_value")
  if("red_value" %in% names(data) && sum(is.finite(data$red_value), na.rm = TRUE) > 3) return("red_value")
  return(NULL)
}

# Time-based outlier detection per location using rolling median/IQR in time (doy)
flag_time_outliers <- function(vec, doy, k=7) {
  # k = window in days for local neighborhood
  if(length(vec) != length(doy)) return(rep(FALSE, length(vec)))
  ord <- order(doy)
  v <- vec[ord]; d <- doy[ord]
  flags <- rep(FALSE, length(v))
  for(i in seq_along(v)) {
    w <- which(abs(d - d[i]) <= k)
    local <- v[w]
    local <- local[is.finite(local)]
  # require slightly fewer local neighbors (relaxed from 5 -> 4)
  if(length(local) >= 4 && is.finite(v[i])) {
      q1 <- quantile(local, 0.25, na.rm = TRUE); q3 <- quantile(local, 0.75, na.rm = TRUE)
      iqr <- q3 - q1; lower <- q1 - 1.5*iqr; upper <- q3 + 1.5*iqr
      flags[i] <- (v[i] < lower || v[i] > upper)
    }
  }
  out <- rep(FALSE, length(vec)); out[ord] <- flags; out
}

# Plot per-location, year-colored time series and save
DO_GAMM_SMOOTH <- TRUE

plot_location_years <- function(data, loc_id, idx_col, out_dir = PLOTS_DIR) {
  d <- data %>% dplyr::filter(.data$location_id == loc_id)
  if(nrow(d) < 5) return(NULL)
  d <- d %>% dplyr::mutate(year = year(date))
  # Determine vegetation label for this location (most frequent Veg if available)
  veg_lab <- NA_character_
  if("Veg" %in% names(d) && any(!is.na(d$Veg))) {
    tb <- sort(table(d$Veg), decreasing = TRUE)
    if(length(tb) > 0) veg_lab <- as.character(names(tb)[1])
  }
  # Outliers by time window
  d$outlier_time <- flag_time_outliers(d[[idx_col]], yday(d$date), k=10)
  # Create a non-outlier subset for connecting lines so lines are not drawn across outlier points
  d_nonout <- d %>% dplyr::filter(!.data$outlier_time)
  title_txt <- paste0("Time Series at ", loc_id,
                      if(!is.na(veg_lab) && nzchar(veg_lab)) paste0(" — Veg: ", veg_lab) else "",
                      " (", idx_col, ")")
  p <- ggplot2::ggplot() +
    # lines only between non-outlier observations (preserve grouping by year)
    ggplot2::geom_line(data = d_nonout, ggplot2::aes(x = yday(date), y = .data[[idx_col]], color = as.factor(.data$year), group = .data$year), alpha = 0.7) +
    # show points for all observations, marking outliers specially
    ggplot2::geom_point(data = d, ggplot2::aes(x = yday(date), y = .data[[idx_col]], color = as.factor(.data$year), group = .data$year, shape = .data$outlier_time), alpha = 0.8) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 4, `FALSE` = 16)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title_txt, x = "Day of Year", y = idx_col, color = "Year", shape = "Time Outlier")

  if(DO_GAMM_SMOOTH) {
  sub <- d %>% dplyr::mutate(doy = yday(date)) %>% dplyr::filter(is.finite(.data[[idx_col]]) & !.data$outlier_time)
    if(nrow(sub) >= 10) {
      # Cyclic smooth over DOY; provide knots for periodicity
      try({
        knts <- list(doy = c(0.5, 366.5))
  fm <- stats::as.formula(paste0(idx_col, " ~ s(doy, bs='cc', k=20)"))
  # unmix response using FVC_est and an estimated barren seasonal curve if available
        dat <- sub
  if("FVC_est" %in% names(dat) && is.finite(mean(dat$FVC_est, na.rm = TRUE)) && "MSAVI" %in% names(df) && "Veg" %in% names(df)) {
          # crude barren seasonal: fit a GCV GAM on all rows labelled barren/bare (visual-only approximation)
          barr <- df %>% dplyr::filter(grepl("barr|bare", tolower(.data$Veg))) %>% dplyr::mutate(doy = yday(date))
          if(nrow(barr) >= 10 && idx_col %in% names(barr)) {
            g_barr <- try(mgcv::gam(stats::as.formula(paste0(idx_col, " ~ s(doy, bs='cc', k=20)")), data = barr, method = "REML"), silent = TRUE)
            if(!inherits(g_barr, "try-error")) {
              barr_mu <- as.numeric(predict(g_barr, newdata = data.frame(doy = dat$doy), se.fit = FALSE))
              covv <- pmax(0, pmin(1, dat$FVC_est))
              ok <- is.finite(covv) & covv > 0.1 & is.finite(barr_mu)
              adj <- dat[[idx_col]]
              adj[ok] <- (adj[ok] - (1 - covv[ok]) * barr_mu[ok]) / covv[ok]
              dat[[idx_col]] <- adj
            }
          }
        }
        g <- mgcv::gam(fm, data = dat, knots = knts, method = "REML")
        grid <- data.frame(doy = 1:365)
        pr <- predict(g, newdata = grid, se.fit = TRUE)
        grid$fit <- as.numeric(pr$fit)
        grid$upr <- grid$fit + 1.96*as.numeric(pr$se.fit)
        grid$lwr <- grid$fit - 1.96*as.numeric(pr$se.fit)
        p <- p +
          ggplot2::geom_ribbon(data = grid, ggplot2::aes(x = .data$doy, ymin = .data$lwr, ymax = .data$upr), inherit.aes = FALSE, fill = "grey70", alpha = 0.3) +
          ggplot2::geom_line(data = grid, ggplot2::aes(x = .data$doy, y = .data$fit), inherit.aes = FALSE, color = "black", linewidth = 1)
      }, silent = TRUE)
    }
  }
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(out_dir, paste0("ts_years_", gsub("[^A-Za-z0-9_]+","_", loc_id), "_", idx_col, ".png"))
  ggplot2::ggsave(out, plot = p, width = 9, height = 4.5, dpi = 150)
  cat("✓ Saved:", out, "\n")
  p
}

# Prepare location ids
if(!("location_id" %in% names(df))) {
  df <- df %>% dplyr::mutate(location_id = paste0(round(target_lat,4), "_", round(target_lon,4)))
}

# Explicitly set vegetation column to 'Veg'
VEG_COL <- "Veg"


## Plot multiple locations per vegetation type, multiple years on DOY x-axis
plot_vegtype_multiloc <- function(data, veg_val, idx_col, out_dir = PLOTS_DIR) {
  if(is.null(VEG_COL) || !(VEG_COL %in% names(data))) return(NULL)
  d <- data %>%
    dplyr::filter(.data[[VEG_COL]] == veg_val, is.finite(.data[[idx_col]]))
  if(nrow(d) < 10) return(NULL)

  # Keep top N locations by observation count for readability
  top_locs <- d %>%
    dplyr::count(.data$location_id, sort = TRUE) %>%
    utils::head(MAX_LOC_PER_VEGTYPE) %>%
    dplyr::pull(.data$location_id)

  d <- d %>% dplyr::filter(.data$location_id %in% top_locs) %>% dplyr::arrange(date)
  if(dplyr::n_distinct(d$location_id) < 2) return(NULL)

  # Ensure year column exists
  if(!"year" %in% names(d)) d <- d %>% dplyr::mutate(year = lubridate::year(date))

  # Per-location outlier detection (same logic as per-location plots)
  d <- d %>% dplyr::group_by(.data$location_id) %>%
    dplyr::mutate(outlier_time = flag_time_outliers(.data[[idx_col]], yday(date), k = 10)) %>%
    dplyr::ungroup()

  # Drop locations that are mostly outliers (relaxed: >50% of observations flagged)
  OUTLIER_LOC_FRAC <- 0.5
  loc_frac <- d %>% dplyr::group_by(.data$location_id) %>% dplyr::summarise(frac_out = mean(.data$outlier_time, na.rm = TRUE))
  bad_locs <- loc_frac %>% dplyr::filter(.data$frac_out > OUTLIER_LOC_FRAC) %>% dplyr::pull(.data$location_id)
  if(length(bad_locs) > 0) {
    # remember which locations were dropped globally
    OUTLIER_LOCATIONS <<- unique(c(OUTLIER_LOCATIONS, as.character(bad_locs)))
  d <- d %>% dplyr::filter(!.data$location_id %in% bad_locs)
  }

  # For plotting, remove individual outlier observations (do not show them)
  d_plot <- d %>% dplyr::filter(!.data$outlier_time)
  if(nrow(d_plot) < 10) return(NULL)

  # Base plot: lines drawn only between non-outlier observations; points show all observations and mark outliers
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = d_plot, ggplot2::aes(x = yday(date), y = .data[[idx_col]], color = .data$location_id, linetype = as.factor(.data$year), group = interaction(.data$location_id, .data$year)), alpha = 0.7) +
    ggplot2::geom_point(data = d, ggplot2::aes(x = yday(date), y = .data[[idx_col]], color = .data$location_id, linetype = as.factor(.data$year), group = interaction(.data$location_id, .data$year), shape = .data$outlier_time), size = 0.8, alpha = 0.8) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 4, `FALSE` = 16)) +
    ggplot2::scale_x_continuous(breaks = c(1,32,60,91,121,152,182,213,244,274,305,335,365), minor_breaks = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste0("By Vegetation Type: ", veg_val, " (", idx_col, ")"),
      x = "Day of Year", y = idx_col, color = "Location", linetype = "Year", shape = "Outlier"
    )

  # Add a single GAM smoothing across all selected locations (common seasonal fit)
  grid <- data.frame(doy = 1:365)
  sub_all <- d_plot %>% dplyr::mutate(doy = yday(date)) %>% dplyr::filter(is.finite(.data[[idx_col]]))
  if(nrow(sub_all) >= 10) {
    try({
      knts <- list(doy = c(0.5, 366.5))
      # Use coverage as a linear predictor when available; fallback to doy-only if coverage missing
      use_cov <- ("coverage" %in% names(sub_all)) && sum(is.finite(sub_all$coverage)) >= 10
  if(use_cov) {
        fm_all <- stats::as.formula(paste0(idx_col, " ~ s(doy, bs='cc', k=20) + coverage"))
        # choose coverage levels to illustrate linear effect (10th, 50th, 90th percentiles)
        cov_q <- as.numeric(quantile(sub_all$coverage, probs = c(0.1, 0.5, 0.9), na.rm = TRUE))
        preds <- lapply(cov_q, function(cv) {
          nd <- data.frame(doy = grid$doy, coverage = cv)
          pr <- predict(mgcv::gam(fm_all, data = sub_all, knots = knts, method = "REML"), newdata = nd, se.fit = TRUE)
          data.frame(doy = grid$doy, fit = as.numeric(pr$fit), upr = as.numeric(pr$fit) + 1.96*as.numeric(pr$se.fit), lwr = as.numeric(pr$fit) - 1.96*as.numeric(pr$se.fit), coverage = cv)
        })
        gf_all <- do.call(rbind, preds)
        # Plot each coverage level as separate line (shows linear variation)
        p <- p +
          ggplot2::geom_ribbon(data = gf_all, ggplot2::aes(x = .data$doy, ymin = .data$lwr, ymax = .data$upr, group = .data$coverage), inherit.aes = FALSE, fill = "grey70", alpha = 0.12, show.legend = FALSE) +
          ggplot2::geom_line(data = gf_all, ggplot2::aes(x = .data$doy, y = .data$fit, group = .data$coverage, linetype = as.factor(.data$coverage)), inherit.aes = FALSE, color = "black", linewidth = 1)
      } else {
        fm_all <- stats::as.formula(paste0(idx_col, " ~ s(doy, bs='cc', k=20)"))
        # unmix response with a crude barren seasonal as above
        dat_all <- sub_all
        if("FVC_est" %in% names(dat_all) && is.finite(mean(dat_all$FVC_est, na.rm = TRUE)) && "Veg" %in% names(df)) {
          barr <- df %>% dplyr::filter(grepl("barr|bare", tolower(.data$Veg))) %>% dplyr::mutate(doy = yday(date))
          if(nrow(barr) >= 10 && idx_col %in% names(barr)) {
            g_barr <- try(mgcv::gam(stats::as.formula(paste0(idx_col, " ~ s(doy, bs='cc', k=20)")), data = barr, method = "REML"), silent = TRUE)
            if(!inherits(g_barr, "try-error")) {
              barr_mu <- as.numeric(predict(g_barr, newdata = data.frame(doy = dat_all$doy), se.fit = FALSE))
              covv <- pmax(0, pmin(1, dat_all$FVC_est))
              ok <- is.finite(covv) & covv > 0.1 & is.finite(barr_mu)
              adj <- dat_all[[idx_col]]
              adj[ok] <- (adj[ok] - (1 - covv[ok]) * barr_mu[ok]) / covv[ok]
              dat_all[[idx_col]] <- adj
            }
          }
        }
        g_all <- mgcv::gam(fm_all, data = dat_all, knots = knts, method = "REML")
        pr_all <- predict(g_all, newdata = grid, se.fit = TRUE)
        gf_all <- data.frame(doy = grid$doy, fit = as.numeric(pr_all$fit), upr = as.numeric(pr_all$fit) + 1.96*as.numeric(pr_all$se.fit), lwr = as.numeric(pr_all$fit) - 1.96*as.numeric(pr_all$se.fit))
        p <- p +
          ggplot2::geom_ribbon(data = gf_all, ggplot2::aes(x = .data$doy, ymin = .data$lwr, ymax = .data$upr), inherit.aes = FALSE, fill = "grey70", alpha = 0.25) +
          ggplot2::geom_line(data = gf_all, ggplot2::aes(x = .data$doy, y = .data$fit), inherit.aes = FALSE, color = "black", linewidth = 1)
      }
    }, silent = TRUE)
  }

  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(out_dir, paste0("ts_by_veg_", gsub("[^A-Za-z0-9_]+","_", veg_val), "_", idx_col, ".png"))
  ggplot2::ggsave(out, plot = p, width = 10, height = 5, dpi = 150)
  cat("✓ Saved:", out, "\n")
  p
}

# Helper to choose best index by coverage
choose_best_index <- function(data, candidates, min_n = 5) {
  cand <- candidates[candidates %in% names(data)]
  if(length(cand) == 0) return(NULL)
  counts <- sapply(cand, function(cn) sum(is.finite(data[[cn]]), na.rm = TRUE))
  counts[is.na(counts)] <- -Inf
  best <- cand[which.max(counts)]
  if(length(best) == 0 || counts[which.max(counts)] < min_n) return(NULL)
  best
}

# Select indices: NDVI, MSAVI, plus new indices (GNDVI, ARI, CRI, RGRI, BGI) if available; also YI and NDRI/RGI

# Always include NDVI and MSAVI if present
sel <- c()
add_if <- function(col, min_n = 5) {
  if(col %in% names(df) && sum(is.finite(df[[col]]), na.rm = TRUE) >= min_n) assign("sel", c(get("sel", inherits = TRUE), col), inherits = TRUE)
}
add_if("NDVI"); add_if("MSAVI")
# New indices
add_if("GNDVI"); add_if("ARI"); add_if("CRI"); add_if("RGRI"); add_if("BGI")
# Existing extras
add_if("YI"); if("NDRI" %in% names(df) && sum(is.finite(df$NDRI), na.rm = TRUE) >= 5) sel <- c(sel, "NDRI") else add_if("RGI")
sel <- unique(sel)
if(length(sel) == 0) {
  idx <- choose_index(df)
  if(is.null(idx)) stop("No plottable index or band found in data")
  sel <- idx
}

for(idx in sel) {
  subdir <- file.path(PLOTS_DIR, idx)
  df_idx <- df %>% dplyr::filter(is.finite(.data[[idx]]))
  if(nrow(df_idx) == 0) next
  # Use all years present in the data
  locs <- df_idx %>% dplyr::count(location_id, sort = TRUE)
  locs <- utils::head(locs, MAX_LOC_PER_INDEX)
  for(li in seq_len(nrow(locs))) {
    plot_location_years(df_idx, locs$location_id[li], idx, out_dir = subdir)
  }

  # Per vegetation-type plots (multiple locations per plot)
  if(!is.null(VEG_COL) && (VEG_COL %in% names(df_idx))) {
    subdir_veg <- file.path(PLOTS_DIR, idx, "by_vegtype")
    veg_vals <- df_idx %>%
      dplyr::filter(!is.na(.data[[VEG_COL]])) %>%
      dplyr::pull(.data[[VEG_COL]]) %>%
      unique()
    for(v in veg_vals) {
      # For each vegetation type, plot only one year (most recent)
      plot_vegtype_multiloc(df_idx, v, idx, out_dir = subdir_veg)
    }
  # Combined phenology GAM across vegetation types for this index
  # Save all combined-index plots to a single directory for easier review
  combined_dir <- file.path(PLOTS_DIR, 'combined_veg_all_indices')
  if(!dir.exists(combined_dir)) dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
    # Fit one GAM per veg type (using coverage if available) and plot all fits together
    veg_gam_list <- list()
    grid <- data.frame(doy = 1:365)
    for(v in veg_vals) {
      dd0 <- df_idx %>% dplyr::filter(.data[[VEG_COL]] == v) %>% dplyr::mutate(doy = yday(date)) %>% dplyr::filter(is.finite(.data[[idx]]))
      n0 <- nrow(dd0)
      if(n0 == 0) next
      # flag and remove outlier observations
      dd <- dd0 %>% dplyr::mutate(outlier_time = flag_time_outliers(.data[[idx]], doy, k=10)) %>% dplyr::filter(!outlier_time)
      n1 <- nrow(dd)
      # report status
      cat(sprintf("Veg '%s': rows_before=%d, rows_after_outlier_filter=%d\n", as.character(v), n0, n1))
      # require at least 6 observations to fit (less strict)
      if(n1 < 6) {
        cat(sprintf("  Skipping '%s': not enough non-outlier rows (%d)\n", as.character(v), n1))
        next
      }
      # attempt GAM fit with visible errors
      tryCatch({
        use_cov <- ("coverage" %in% names(dd)) && sum(is.finite(dd$coverage)) >= 8
        if(use_cov) {
          fm <- stats::as.formula(paste0(idx, " ~ s(doy, bs='cc', k=20) + coverage"))
          newdata <- data.frame(doy = grid$doy, coverage = 1)
        } else {
          fm <- stats::as.formula(paste0(idx, " ~ s(doy, bs='cc', k=20)"))
          newdata <- data.frame(doy = grid$doy)
        }
        g <- mgcv::gam(fm, data = dd, method = "REML")
        pr <- predict(g, newdata = newdata, se.fit = TRUE)
        veg_gam_list[[as.character(v)]] <- data.frame(doy = grid$doy, fit = as.numeric(pr$fit), upr = as.numeric(pr$fit) + 1.96*as.numeric(pr$se.fit), lwr = as.numeric(pr$fit) - 1.96*as.numeric(pr$se.fit), Veg = as.character(v))
        cat(sprintf("  Fitted GAM for '%s' (n=%d)\n", as.character(v), n1))
      }, error = function(e) {
        cat(sprintf("  GAM error for '%s': %s\n", as.character(v), conditionMessage(e)))
      })
    }
    if(length(veg_gam_list) > 0) {
      vg_all <- do.call(rbind, veg_gam_list)
      p_comb <- ggplot2::ggplot(vg_all, ggplot2::aes(x = doy, y = fit, color = Veg, fill = Veg)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr), alpha = 0.2, color = NA) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste0("Combined phenology GAM - ", idx), subtitle = "Grouped by Vegetation Type", x = "Day of Year", y = idx, color = "Veg", fill = "Veg")
      outf <- file.path(combined_dir, paste0("combined_phenology_", idx, ".png"))
      ggplot2::ggsave(outf, plot = p_comb, width = 10, height = 5, dpi = 150)
      cat("✓ Saved combined phenology:", outf, "\n")
    }
  } else {
    cat("ℹ Skipping vegetation-type plots: no vegetation column found in this dataset.\n")
  }
}

if(length(OUTLIER_LOCATIONS) > 0) {
  cat("\n=== Outlier locations dropped (mostly outlier observations) ===\n")
  cat(paste(sort(unique(OUTLIER_LOCATIONS)), collapse = "\n"), "\n")
} else {
  cat("\n=== No outlier-heavy locations were dropped ===\n")
}

cat("=== Visualization complete ===\n")
