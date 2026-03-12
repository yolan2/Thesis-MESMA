# =============================================================================
# satellite_bias_check.R
# Compare Landsat 4/5/7 vs Landsat 8/9 observations at the same location and date.
# Goal: detect inter-satellite biases and assess cross-sensor consistency.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)   # panel layouts  (install if missing: install.packages("patchwork"))
})

# ── Config ────────────────────────────────────────────────────────────────────
# Primary training CSV path.  Additional files (e.g. inference sets) may be
# combined via INFERENCE_CSV defined in mesma_config.R below.
INPUT_CSV <- "C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_low (3).csv"
OUT_DIR   <- "satellite_bias_check"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Optionally source mesma_config to pick up any extra CSVs the user has set
if (file.exists("mesma_config.R")) {
  source("mesma_config.R")
}

# Compose list of CSVs to load for bias calculation; include training plus
# any inference or additional files specified in config variables.
BIAS_INPUTS <- unique(c(INPUT_CSV,
                        if (exists("INFERENCE_CSV") && nzchar(INFERENCE_CSV)) INFERENCE_CSV,
                        if (exists("ADDITIONAL_BIAS_CSVS")) ADDITIONAL_BIAS_CSVS))
# ensure non-null
BIAS_INPUTS <- BIAS_INPUTS[!is.na(BIAS_INPUTS) & nzchar(BIAS_INPUTS)]
cat(sprintf("[BIAS CHECK] Loading %d CSV(s) for bias computation:\n  %s\n",
            length(BIAS_INPUTS), paste(BIAS_INPUTS, collapse = "\n  ")))

BANDS <- c("Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2")

# -----------------------------------------------------------------------------
# additional spectral indices to compute and compare
INDICES <- c("NDVI", "EVI", "PPI", "MSAVI")

# helper to add indices to a data.table in place
# This now relies on the centralized ppi_helpers.R for PPI so that the
# calculation (M parameter, dvi_soil baseline) remains consistent with the
# MESMA fit pipeline.  We force a fixed M=0.584568 and a constant soil baseline
# for all sensors (DEFAULT_DVI_SOIL) to satisfy the user's requirement.
compute_indices <- function(dt) {
  # NDVI, EVI, MSAVI are simple algebraic formulas
  dt[, NDVI := (NIR - Red) / (NIR + Red)]
  dt[, EVI  := 2.5 * (NIR - Red) / (NIR + 6*Red - 7.5*Blue + 1)]
  dt[, MSAVI := (2 * NIR + 1 - sqrt((2 * NIR + 1)^2 - 8 * (NIR - Red))) / 2]

  # Prepare for helper-based PPI
  dt[, DVI := NIR - Red]   # required by add_ppi_columns

  # apply helper; it returns a data.frame so coerce back to data.table
  tmp <- add_ppi_columns(dt, dvi_soil = DEFAULT_DVI_SOIL)
  setDT(tmp)
  dt[, PPI := tmp$PPI]
}

# robust outlier filter used only for affine OLS coefficient estimation
# returns a logical mask for finite paired observations that remain after
# removing extreme values in x, y, and their difference using MAD thresholds.
robust_pair_mask <- function(x, y, mad_thresh = 4) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10) return(ok)

  x_ok <- x[ok]
  y_ok <- y[ok]
  d_ok <- y_ok - x_ok

  keep_local <- rep(TRUE, length(x_ok))
  for (vals in list(x_ok, y_ok, d_ok)) {
    medv <- stats::median(vals, na.rm = TRUE)
    madv <- stats::mad(vals, center = medv, constant = 1, na.rm = TRUE)
    if (is.finite(madv) && madv > 0) {
      keep_local <- keep_local & (abs(vals - medv) <= mad_thresh * madv)
    }
  }

  keep <- rep(FALSE, length(x))
  keep[which(ok)] <- keep_local
  keep
}

# all features we'll compute biases for
FEATURES <- c(BANDS, INDICES)

# Maximum days allowed between two observations to still call them "same date"
# (0 = exact same date only)
DATE_TOL_DAYS <- 0L

# ── Load data ─────────────────────────────────────────────────────────────────
# ensure PPI helpers available and constants fixed
if (file.exists("ppi_helpers.r")) source("ppi_helpers.r")
# enforce global M constant used by add_ppi_columns/ppi()
PPI_FULL_VEG_COVER <- 0.584568
# choose a constant DVI soil baseline for all sensors (OLI+ETM) per directive
# note: DEFAULT_DVI_SOIL is defined in ppi_helpers.R; fallback if missing
if (!exists("DEFAULT_DVI_SOIL")) DEFAULT_DVI_SOIL <- 0.0308

cat("Loading CSV(s) …\n")
if (length(BIAS_INPUTS) == 0) stop("No input CSVs specified for bias check")
# read and rbind; fill=TRUE to allow differing column sets
list_df <- lapply(BIAS_INPUTS, function(f) {
  if (!file.exists(f)) stop(sprintf("Input CSV not found: %s", f))
  dt <- fread(f)
  setnames(dt, names(dt), trimws(names(dt)))
  dt
})
df <- data.table::rbindlist(list_df, fill = TRUE)
cat(sprintf("  combined dataset contains %d rows\n", nrow(df)))

cat(sprintf("  %d rows | satellites: %s\n",
    nrow(df), paste(sort(unique(df$satellite)), collapse = ", ")))

# compute derived indices so we can compare them later
compute_indices(df)

# ── Identify co-located, same-date observations from different satellites ─────
# Strategy: self-join on (location_id, date) keeping only pairs where satellites differ.
cat("Finding co-temporal multi-satellite pairs …\n")

df[, date := as.IDate(date)]  # ensure Date type for arithmetic

# Separate the two sensor families
ls_old <- df[satellite == "LANDSAT_457"]
ls_new <- df[satellite == "LANDSAT_89"]

# Inner join on location_id × date
pairs <- ls_old[ls_new, on = .(location_id, date), nomatch = 0L, allow.cartesian = TRUE]
# After join: columns from ls_old get no suffix (x), ls_new get .i (y convention in data.table)
# Rename for clarity; handle all requested features
for (f in FEATURES) {
  setnames(pairs, f,        paste0(f, "_457"))
  setnames(pairs, paste0("i.", f), paste0(f, "_89"))
}
# Drop duplicated metadata columns added by join
drop_cols <- grep("^i\\.", names(pairs), value = TRUE)
pairs[, (drop_cols) := NULL]

cat(sprintf("  %d matched loc-date pairs\n", nrow(pairs)))

if (nrow(pairs) == 0L) {
  stop("No matching pairs found. Check satellite labels or DATE_TOL_DAYS setting.")
}

# ── Per-band bias statistics ──────────────────────────────────────────────────
cat("\n=== Per-feature bias (LANDSAT_457 minus LANDSAT_89) ===\n")

bias_stats <- rbindlist(lapply(FEATURES, function(b) {
  col457 <- paste0(b, "_457")
  col89  <- paste0(b, "_89")
  x   <- pairs[[col89]]
  y   <- pairs[[col457]]
  diff <- y - x

  # OLS affine fit: ETM+ = slope * OLI + intercept  (Roy et al. 2016 approach)
  # Fit only after robust outlier removal on the matched pairs.
  ok <- robust_pair_mask(x, y)
  ols_slope     <- NA_real_
  ols_intercept <- NA_real_
  ols_p_value   <- NA_real_
  ols_significant <- FALSE
  if (sum(ok) >= 10) {
    fit <- lm(y[ok] ~ x[ok])
    ols_slope     <- unname(coef(fit)[2])
    ols_intercept <- unname(coef(fit)[1])
    ols_p_value   <- summary(fit)$coefficients[2, 4]
    ols_significant <- !is.na(ols_p_value) && ols_p_value < 0.05
  }

  data.table(
    band          = b,
    n             = sum(!is.na(diff)),
    n_ols         = sum(ok, na.rm = TRUE),
    ols_slope     = ols_slope,
    ols_intercept = ols_intercept,
    ols_p_value   = ols_p_value,
    ols_significant = ols_significant,
    median_bias   = median(diff, na.rm = TRUE),
    sd_bias       = sd(diff, na.rm = TRUE),
    rmse          = sqrt(mean(diff^2, na.rm = TRUE)),
    mae           = mean(abs(diff), na.rm = TRUE),
    r2            = cor(x, y, use = "complete.obs")^2
  )
}))

bias_stats[, `:=`(
  ols_slope     = round(ols_slope,     6),
  ols_intercept = round(ols_intercept, 6),
  ols_p_value   = round(ols_p_value,   6),
  median_bias   = round(median_bias,   5),
  sd_bias       = round(sd_bias,       5),
  rmse          = round(rmse,          5),
  mae           = round(mae,           5),
  r2            = round(r2,            4)
)]
print(bias_stats)

fwrite(bias_stats, file.path(OUT_DIR, "bias_stats_features.csv"))
cat(sprintf("\nSaved: %s/bias_stats_features.csv\n", OUT_DIR))
cat("  Columns: ols_slope + ols_intercept + ols_p_value + ols_significant added (affine ETM+ = slope*OLI + intercept)\n")

# ── Melt to long format for plotting ─────────────────────────────────────────
# The long table is reused for several figures below.  ggplot's geom_violin
# silently drops any factor level that has fewer than two non-NA y-values
# (it prints a warning "Groups with fewer than two datapoints have been
# dropped").  that's exactly what the user reported: SWIR1/SWIR2 only
# produced a white box from the boxplot while the violin in the same slot
# vanished.  those bands had extremely low variance after whatever filtering
# was applied, so the density estimator collapsed and the polygon was omitted.
# We therefore create a separate column for plotting where we jitter groups
# with <2 unique values.  the jitter is infinitesimal and only affects the
# plot – the original `diff` is left untouched for statistics.

pairs_long <- melt(
  pairs,
  id.vars = c("location_id", "date", "year", "lat", "lon", "vegetation"),
  measure.vars = patterns("_457$", "_89$"),
  variable.name = "band_idx",
  value.name    = c("val_457", "val_89")
)
# convert to factor so ggplot respects the order defined in FEATURES; this
# guarantees that raw bands (first entries of FEATURES) appear left of all
# computed indices.
pairs_long[, band := factor(FEATURES[band_idx], levels = FEATURES)]
pairs_long[, diff := val_457 - val_89]

# create a plotting column that will be jittered if needed
pairs_long[, diff_plot := diff]
# identify bands with fewer than two unique non-NA diff values
short_bands <- pairs_long[!is.na(diff), .(uniq = uniqueN(diff)), by = band][uniq < 2, band]
if (length(short_bands) > 0) {
  warning("Band(s) with <2 unique differences (no violin shape will be drawn): ",
          paste(short_bands, collapse = ", "))
  # add tiny noise so geom_violin will at least draw a flat shape
  pairs_long[band %in% short_bands, diff_plot := diff +
               runif(.N, min = -1e-6, max = 1e-6)]
}

# ── Plot 1: Scatter plots (1:1 line) per band ─────────────────────────────────
cat("Plotting scatter plots …\n")

# Subsample for plotting speed (max 20k points per band)
MAX_PTS <- 20000L
set.seed(42)
pairs_samp <- pairs_long[, .SD[sample(.N, min(.N, MAX_PTS))], by = band]

scatter_plots <- lapply(FEATURES, function(b) {
  d <- pairs_samp[band == b]
  lim <- range(c(d$val_457, d$val_89), na.rm = TRUE)
  bias_val <- bias_stats[band == b, median_bias]
  r2_val   <- bias_stats[band == b, r2]
  ggplot(d, aes(x = val_457, y = val_89)) +
    geom_point(alpha = 0.15, size = 0.4, colour = "#2166ac") +
    geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.7) +
    geom_smooth(method = "lm", se = FALSE, colour = "#f4a582", linewidth = 0.7,
                linetype = "dashed") +
    coord_equal(xlim = lim, ylim = lim) +
    labs(
      title = b,
      subtitle = sprintf("mean bias = %+.4f  |  R² = %.4f", bias_val, r2_val),
      x = "LANDSAT 4/5/7",
      y = "LANDSAT 8/9"
    ) +
    theme_bw(base_size = 9)
})

scatter_panel <- wrap_plots(scatter_plots, ncol = 3) +
  plot_annotation(title = "Feature Values: L457 vs L89",
                  subtitle = sprintf("n = %d matched pairs  |  red = 1:1, orange = OLS fit", nrow(pairs)))

ggsave(file.path(OUT_DIR, "scatter_by_feature.png"),
       scatter_panel, width = 14, height = 9, dpi = 150)
cat(sprintf("Saved: %s/scatter_by_feature.png\n", OUT_DIR))

# ── Plot 2: Bias distribution (violin + boxplot) ──────────────────────────────
cat("Plotting bias distributions …\n")

# create a colour vector large enough for all features; the default
# Set2 palette only holds eight colours, so interpolate if we need more.
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  stop("Package 'RColorBrewer' is required for coloured plots")
}
fill_cols <- RColorBrewer::brewer.pal(min(length(FEATURES), 8), "Set2")
if (length(FEATURES) > length(fill_cols)) {
  fill_cols <- colorRampPalette(fill_cols)(length(FEATURES))
}
names(fill_cols) <- FEATURES

p_bias_dist <- ggplot(pairs_long, aes(x = band, y = diff_plot, fill = band)) +
  # vertical separator between raw bands and indices (after BANDS entries)
  geom_vline(xintercept = length(BANDS) + 0.5, colour = "grey40", linetype = "dotted",
             inherit.aes = FALSE) +
  # draw the box on a transparent background so it doesn't obscure a
  # violin that happens to lie entirely inside the IQR.  the box outline
  # is retained for reference.
  geom_boxplot(width = 0.15, outlier.shape = NA, colour = "grey10", fill = NA) +
  geom_violin(trim = TRUE, alpha = 0.6, colour = NA,
              width = 0.8, scale = "width", na.rm = TRUE) +
  geom_hline(yintercept = 0, colour = "red", linewidth = 0.7) +
  scale_fill_manual(values = fill_cols, guide = "none") +
  labs(
    title = "Bias by Feature (L457 - L89)",
    x = "Feature", y = "Difference"
  ) +
  ylim(-0.2, 0.2) +
  theme_bw(base_size = 11)

ggsave(file.path(OUT_DIR, "bias_distribution_by_feature.png"),
       p_bias_dist, width = 10, height = 6, dpi = 150)
cat(sprintf("Saved: %s/bias_distribution_by_feature.png\n", OUT_DIR))

# ── Plot 3: Mean bias through time ────────────────────────────────────────────
cat("Plotting temporal bias trend …\n")

time_bias <- pairs_long[, .(mean_bias = mean(diff, na.rm = TRUE),
                             sd_bias   = sd(diff,   na.rm = TRUE),
                             n         = .N),
                         by = .(year, band)]
time_bias[, se := sd_bias / sqrt(n)]

p_time <- ggplot(time_bias, aes(x = year, y = mean_bias, colour = band,
                                 group = band, fill = band)) +
  geom_ribbon(aes(ymin = mean_bias - se, ymax = mean_bias + se),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  facet_wrap(~band, scales = "free_y", ncol = 3) +
  scale_colour_brewer(palette = "Dark2", guide = "none") +
  scale_fill_brewer(palette = "Dark2", guide = "none") +
  labs(
    title = "Mean Bias by Year (L457 - L89)",
    x = "Year", y = "Mean difference ± SE"
  ) +
  theme_bw(base_size = 9)

ggsave(file.path(OUT_DIR, "bias_trend_by_year.png"),
       p_time, width = 14, height = 8, dpi = 150)
cat(sprintf("Saved: %s/bias_trend_by_year.png\n", OUT_DIR))

# ── Plot 4: Bias by vegetation class (if available) ───────────────────────────
if ("vegetation" %in% names(pairs_long) && any(!is.na(pairs_long$vegetation) & pairs_long$vegetation != "")) {
  cat("Plotting bias by vegetation class …\n")

  p_veg <- ggplot(pairs_long[vegetation != "" & !is.na(vegetation)],
                  aes(x = vegetation, y = diff, fill = vegetation)) +
    geom_boxplot(outlier.shape = NA, staplewidth = 0.4) +
    geom_hline(yintercept = 0, colour = "red", linewidth = 0.6) +
    facet_wrap(~band, scales = "free_y", ncol = 3) +
    scale_fill_brewer(palette = "Paired", guide = "none") +
    coord_cartesian(ylim = quantile(pairs_long$diff, c(0.01, 0.99), na.rm = TRUE)) +
    labs(
      title = "Bias by Vegetation (L457 - L89)",
      x = "Vegetation class", y = "Difference"
    ) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(OUT_DIR, "bias_by_vegetation.png"),
         p_veg, width = 14, height = 9, dpi = 150)
  cat(sprintf("Saved: %s/bias_by_vegetation.png\n", OUT_DIR))
}

# ── Plot 5: DOY (seasonal) bias ───────────────────────────────────────────────
cat("Plotting seasonal bias …\n")

if ("doy" %in% names(pairs)) {
  # propagate doy into long table (deduplicate key before merge to avoid cartesian expansion)
  doy_key <- unique(pairs[, .(location_id, date, doy)])
  pairs_long2 <- merge(pairs_long,
                       doy_key,
                       by = c("location_id", "date"), all.x = TRUE)
  pairs_long2[, doy_bin := cut(doy, breaks = seq(0, 365, by = 30),
                                labels = FALSE, include.lowest = TRUE)]

  doy_bias <- pairs_long2[, .(mean_bias = mean(diff, na.rm = TRUE)),
                           by = .(doy_bin, band)]

  p_doy <- ggplot(doy_bias, aes(x = doy_bin, y = mean_bias, colour = band,
                                 group = band)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    facet_wrap(~band, scales = "free_y", ncol = 3) +
    scale_colour_brewer(palette = "Dark2", guide = "none") +
    scale_x_continuous(breaks = 1:12,
                       labels = month.abb) +
    labs(
      title = "Monthly Bias (L457 - L89)",
      x = "Month", y = "Mean difference"
    ) +
    theme_bw(base_size = 9)

  ggsave(file.path(OUT_DIR, "bias_seasonal.png"),
         p_doy, width = 14, height = 8, dpi = 150)
  cat(sprintf("Saved: %s/bias_seasonal.png\n", OUT_DIR))
}

# ── Summary table: per-location average bias ─────────────────────────────────
cat("Writing per-location bias summary …\n")

loc_bias <- pairs_long[, .(mean_bias = mean(diff, na.rm = TRUE),
                            n_pairs   = .N),
                        by = .(location_id, lat, lon, band)]
fwrite(loc_bias, file.path(OUT_DIR, "per_location_bias.csv"))
cat(sprintf("Saved: %s/per_location_bias.csv\n", OUT_DIR))

# ── Done ──────────────────────────────────────────────────────────────────────
cat("\n=== satellite_bias_check.R complete ===\n")
cat(sprintf("All outputs written to: %s/\n", OUT_DIR))
