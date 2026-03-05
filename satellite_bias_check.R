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
INPUT_CSV <- "C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_low (3).csv"
OUT_DIR   <- "satellite_bias_check"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

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

cat("Loading CSV …\n")
df <- fread(INPUT_CSV)
setnames(df, names(df), trimws(names(df)))   # strip any accidental whitespace

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
  diff   <- pairs[[col457]] - pairs[[col89]]
  data.table(
    band   = b,
    n      = sum(!is.na(diff)),
    mean_bias = mean(diff, na.rm = TRUE),
    median_bias = median(diff, na.rm = TRUE),
    sd_bias = sd(diff, na.rm = TRUE),
    rmse   = sqrt(mean(diff^2, na.rm = TRUE)),
    mae    = mean(abs(diff), na.rm = TRUE),
    r2     = cor(pairs[[col457]], pairs[[col89]], use = "complete.obs")^2
  )
}))

bias_stats[, `:=`(
  mean_bias   = round(mean_bias,   5),
  median_bias = round(median_bias, 5),
  sd_bias     = round(sd_bias,     5),
  rmse        = round(rmse,        5),
  mae         = round(mae,         5),
  r2          = round(r2,          4)
)]
print(bias_stats)

fwrite(bias_stats, file.path(OUT_DIR, "bias_stats_features.csv"))
cat(sprintf("\nSaved: %s/bias_stats_features.csv\n", OUT_DIR))

# ── Melt to long format for plotting ─────────────────────────────────────────
pairs_long <- melt(
  pairs,
  id.vars = c("location_id", "date", "year", "lat", "lon", "vegetation"),
  measure.vars = patterns("_457$", "_89$"),
  variable.name = "band_idx",
  value.name    = c("val_457", "val_89")
)
pairs_long[, band := FEATURES[band_idx]]
pairs_long[, diff := val_457 - val_89]

# ── Plot 1: Scatter plots (1:1 line) per band ─────────────────────────────────
cat("Plotting scatter plots …\n")

# Subsample for plotting speed (max 20k points per band)
MAX_PTS <- 20000L
set.seed(42)
pairs_samp <- pairs_long[, .SD[sample(.N, min(.N, MAX_PTS))], by = band]

scatter_plots <- lapply(FEATURES, function(b) {
  d <- pairs_samp[band == b]
  lim <- range(c(d$val_457, d$val_89), na.rm = TRUE)
  bias_val <- bias_stats[band == b, mean_bias]
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
  plot_annotation(title = "Feature values: LANDSAT 4/5/7 vs LANDSAT 8/9 (same location & date)",
                  subtitle = sprintf("n = %d matched pairs  |  red = 1:1, orange = OLS fit", nrow(pairs)))

ggsave(file.path(OUT_DIR, "scatter_by_feature.png"),
       scatter_panel, width = 14, height = 9, dpi = 150)
cat(sprintf("Saved: %s/scatter_by_feature.png\n", OUT_DIR))

# ── Plot 2: Bias distribution (violin + boxplot) ──────────────────────────────
cat("Plotting bias distributions …\n")

p_bias_dist <- ggplot(pairs_long, aes(x = band, y = diff, fill = band)) +
  geom_violin(trim = TRUE, alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, colour = "grey10", fill = "white") +
  geom_hline(yintercept = 0, colour = "red", linewidth = 0.7) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(
    title = "Bias distribution per feature  (LANDSAT_457 − LANDSAT_89)",
    x = "Feature", y = "Difference"
  ) +
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
    title = "Mean bias per year  (LANDSAT_457 − LANDSAT_89)",
    x = "Year", y = "Mean difference ± SE"
  ) +
  theme_bw(base_size = 9)

ggsave(file.path(OUT_DIR, "bias_trend_by_year.png"),
       p_time, width = 14, height = 8, dpi = 150)
cat(sprintf("Saved: %s/bias_trend_by_year.png\n", OUT_DIR))

# ── Plot 4: Spatial bias maps ─────────────────────────────────────────────────
cat("Plotting spatial bias maps for all features …\n")

# compute mean bias per location for every band/feature
spatial_bias <- pairs_long[, .(mean_bias = mean(diff, na.rm = TRUE)),
                            by = .(location_id, lat, lon, band)]

# loop over each feature to create a separate map file
for (b in FEATURES) {
  sb <- spatial_bias[band == b]
  if (nrow(sb) == 0) next

  p_map <- ggplot(sb, aes(x = lon, y = lat, colour = mean_bias)) +
    geom_point(size = 1.2, alpha = 0.7) +
    scale_colour_gradient2(
      low = "#d6604d", mid = "white", high = "#4393c3",
      midpoint = 0, name = sprintf("%s bias\n(457−89)", b)
    ) +
    coord_quickmap() +
    labs(
      title = sprintf("Spatial distribution of %s bias  (LANDSAT_457 − LANDSAT_89)", b),
      x = "Longitude", y = "Latitude"
    ) +
    theme_bw(base_size = 11)

  fname <- sprintf("spatial_bias_%s.png", b)
  ggsave(file.path(OUT_DIR, fname),
         p_map, width = 10, height = 7, dpi = 150)
  cat(sprintf("Saved: %s/%s\n", OUT_DIR, fname))
}

# ── Plot 5: Bias by vegetation class (if available) ───────────────────────────
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
      title = "Bias by vegetation class  (LANDSAT_457 − LANDSAT_89)",
      x = "Vegetation class", y = "Difference"
    ) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(OUT_DIR, "bias_by_vegetation.png"),
         p_veg, width = 14, height = 9, dpi = 150)
  cat(sprintf("Saved: %s/bias_by_vegetation.png\n", OUT_DIR))
}

# ── Plot 6: DOY (seasonal) bias ───────────────────────────────────────────────
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
      title = "Seasonal (monthly) bias  (LANDSAT_457 − LANDSAT_89)",
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
