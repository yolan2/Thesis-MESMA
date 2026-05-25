# =============================================================================
# plot_mesma_maps.R
# Maps of MESMA inference output: vegetation fraction change 1986-2025
# =============================================================================
# Reads inference_results.csv and produces:
#   1. Per-class fraction maps (one panel per decade or selected years)
#   2. Change maps: fraction(end) - fraction(start) per class
#   3. Dominant vegetation map per year
# Points are aggregated into a raster grid (mean per cell) — no interpolation.
# =============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)

# --- Config ------------------------------------------------------------------
INFERENCE_CSV <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/R_MESMA/inference_results/inference_results_Landsat_Harmonized_Bands_1985_2025_low_3_.csv"
OUT_DIR       <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/R_MESMA/maps"
YEAR_START    <- 1986
YEAR_END      <- 2025
# Years to show in the faceted time-slice maps (NULL = all years)
SLICE_YEARS   <- sort(unique(c(seq(YEAR_START, YEAR_END, by = 5), YEAR_END)))
# Classes to map (NULL = all non-barren)
VEG_CLASSES   <- NULL   # e.g. c("populus", "tamarix", "herbs")
DPI           <- 150
# Grid cell size in degrees (~1 km at mid-latitudes ≈ 0.01°)
GRID_RES      <- 0.01
TREND_METRIC_COL <- "coef"      # PPI-normalized vegetation fraction column
REQUIRE_PPI_NORMALIZED_TREND <- TRUE
# -----------------------------------------------------------------------------

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

df_raw <- readr::read_csv(INFERENCE_CSV, show_col_types = FALSE)

resolve_trend_metric_col <- function(df, preferred_col = "coef") {
  candidates <- unique(c(preferred_col, "coef_ppi_norm", "coef_abs", "ppi_norm_coef", "coef"))
  for (col in candidates) {
    if (col %in% names(df)) {
      vals <- suppressWarnings(as.numeric(df[[col]]))
      if (any(is.finite(vals))) return(col)
    }
  }
  NA_character_
}

assert_not_unnormalized_metric <- function(df, metric_col) {
  unnormalized_cols <- c("coef_raw", "coef_rel", "rel_coef", "raw_coef", "unnormalized_coef")
  present_unorm <- intersect(unnormalized_cols, names(df))
  if (metric_col %in% present_unorm) {
    stop(sprintf("Selected trend metric '%s' is un-normalized. Use a PPI-normalized metric.", metric_col))
  }
}

trend_metric_col <- resolve_trend_metric_col(df_raw, preferred_col = TREND_METRIC_COL)
if (is.na(trend_metric_col)) stop("No usable PPI-normalized trend metric column found in inference CSV")
assert_not_unnormalized_metric(df_raw, trend_metric_col)
if (isTRUE(REQUIRE_PPI_NORMALIZED_TREND) && grepl("raw|rel|unnormal", trend_metric_col, ignore.case = TRUE)) {
  stop(sprintf("Resolved trend metric '%s' is not PPI-normalized.", trend_metric_col))
}
cat(sprintf("[MAP] Using PPI-normalized trend metric column: %s\n", trend_metric_col))

# Normalize class names
df_raw <- df_raw %>%
  mutate(
    Veg        = tolower(trimws(Veg)),
    pheno_year = as.integer(pheno_year),
    coef       = as.numeric(.data[[trend_metric_col]])
  ) %>%
  filter(pheno_year >= YEAR_START, pheno_year <= YEAR_END)

# Fill missing class fractions with 0 for each location and year
all_locs <- df_raw %>% distinct(location_id, lon, lat)
all_years <- df_raw %>% distinct(pheno_year)
all_classes <- df_raw %>% distinct(Veg)

grid_template <- expand_grid(all_locs, all_years, all_classes)

# Collapse variants: sum coef per location-year-class
df <- df_raw %>%
  group_by(location_id, pheno_year, lat, lon, Veg) %>%
  summarise(coef = sum(coef, na.rm = TRUE), .groups = "drop") %>%
  right_join(grid_template, by = c("location_id", "lon", "lat", "pheno_year", "Veg")) %>%
  mutate(coef = replace_na(coef, 0))

# Filter to requested veg classes (exclude barren by default)
if (!is.null(VEG_CLASSES)) {
  df <- df %>% filter(Veg %in% tolower(VEG_CLASSES))
} else {
  df <- df %>% filter(Veg != "barren")
}

veg_classes <- sort(unique(df$Veg))
cat(sprintf("[MAP] Classes: %s | Locations: %d | Years: %d-%d\n",
            paste(veg_classes, collapse = ", "),
            n_distinct(df$location_id), YEAR_START, YEAR_END))

# --- Raster aggregation helper -----------------------------------------------
# Snaps each point to the nearest grid cell centre and averages coef values.
# No spatial interpolation — empty cells stay NA.
snap_to_grid <- function(df, res = GRID_RES) {
  df %>%
    mutate(
      lon_cell = round(floor(lon / res) * res + res / 2, 6),
      lat_cell = round(floor(lat / res) * res + res / 2, 6)
    ) %>%
    group_by(lon_cell, lat_cell, pheno_year, Veg) %>%
    summarise(coef = mean(coef, na.rm = TRUE), .groups = "drop") %>%
    rename(lon = lon_cell, lat = lat_cell)
}

df_grid <- snap_to_grid(df)
cat(sprintf("[MAP] Grid resolution: %.3f deg | Grid cells: %d\n",
            GRID_RES, n_distinct(paste(df_grid$lon, df_grid$lat))))

# Colour palettes
dom_colors <- c(
  "herbs"   = "#9ACD32",
  "populus" = "#006400",
  "tamarix" = "#D95F02"
)
missing_cls <- setdiff(veg_classes, names(dom_colors))
if (length(missing_cls) > 0) {
  dom_colors[missing_cls] <- scales::hue_pal()(length(missing_cls))
}

# --- Reference locations -------------------------------------------------------
ref_locs <- data.frame(
  name = c("Taitema Lake", "Nuo'erma", "Da Xihaizi"),
  lon  = c(88.42063426100638, 88.128011, 87.559930),
  lat  = c(39.60191330443622, 40.282116, 40.558023)
)

# Layers added to every spatial map for context
map_context_layers <- list(
  geom_point(data = ref_locs, aes(x = lon, y = lat),
             inherit.aes = FALSE, shape = 21, size = 2.5,
             fill = "white", color = "black", stroke = 0.8),
  ggrepel::geom_text_repel(data = ref_locs, aes(x = lon, y = lat, label = name),
                            inherit.aes = FALSE, size = 2.5, fontface = "bold",
                            box.padding = 0.4, point.padding = 0.3,
                            segment.size = 0.3, color = "black",
                            bg.color = "white", bg.r = 0.15)
)

# Shared map theme
map_theme <- theme_minimal() +
  theme(
    panel.grid    = element_blank(),
    axis.title    = element_blank(),
    axis.text     = element_text(size = 7),
    legend.title  = element_text(size = 8),
    legend.text   = element_text(size = 7),
    strip.text    = element_text(size = 8),
    plot.title    = element_text(hjust = 0.5, size = 11),
    plot.subtitle = element_text(hjust = 0.5, size = 8)
  )

# =============================================================================
# 1. Time-slice raster maps (one facet per year, one file per class)
# =============================================================================
df_slices <- df_grid %>% filter(pheno_year %in% SLICE_YEARS)
slice_limits <- quantile(df_slices$coef, c(0.02, 0.98), na.rm = TRUE)

for (cls in veg_classes) {
  df_cls <- df_slices %>% filter(Veg == cls)
  if (nrow(df_cls) == 0) next

  p <- ggplot(df_cls, aes(x = lon, y = lat, fill = coef)) +
    geom_tile(width = GRID_RES, height = GRID_RES) +
    facet_wrap(~pheno_year, nrow = 2) +
    scale_fill_distiller(palette = "YlGn", direction = 1,
                         limits = slice_limits, oob = scales::squish,
                         name = "Fraction", na.value = "grey92") +
    map_context_layers +
    labs(title    = sprintf("%s \u2014 vegetation fraction", tools::toTitleCase(cls)),
         subtitle = sprintf("Grid: %.3f\u00b0 | Years: %s",
                            GRID_RES, paste(SLICE_YEARS, collapse = ", "))) +
    coord_equal() +
    map_theme

  fn <- file.path(OUT_DIR, sprintf("map_fraction_%s_slices.png", cls))
  ggsave(fn, p, width = 10, height = 5, dpi = DPI)
  cat(sprintf("[MAP] Saved: %s\n", fn))
}

# =============================================================================
# 1b. 2025 fraction maps (single-year, one file per class)
# =============================================================================
df_2025 <- df_grid %>% filter(pheno_year == 2025)
frac_2025_limits <- quantile(df_2025$coef, c(0.02, 0.98), na.rm = TRUE)

if (nrow(df_2025) > 0) {
  for (cls in veg_classes) {
    df_cls <- df_2025 %>% filter(Veg == cls)
    if (nrow(df_cls) == 0) next

    p <- ggplot(df_cls, aes(x = lon, y = lat, fill = coef)) +
      geom_tile(width = GRID_RES, height = GRID_RES) +
      scale_fill_distiller(palette = "YlGn", direction = 1,
                           limits = frac_2025_limits, oob = scales::squish,
                           name = "Fraction", na.value = "grey92") +
      map_context_layers +
      labs(title    = sprintf("%s \u2014 vegetation fraction 2025",
                              tools::toTitleCase(cls)),
           subtitle = sprintf("Grid: %.3f\u00b0", GRID_RES)) +
      coord_equal() +
      map_theme

    fn <- file.path(OUT_DIR, sprintf("map_fraction_%s_2025.png", cls))
    ggsave(fn, p, width = 7, height = 5, dpi = DPI)
    cat(sprintf("[MAP] Saved: %s\n", fn))
  }
} else {
  cat("[MAP] No 2025 data found — skipping 2025 fraction maps.\n")
}

# =============================================================================
# 2. Change raster maps: mean(last 3 years) - mean(first 3 years)
# =============================================================================
years_all <- sort(unique(df_grid$pheno_year))
yr_early  <- head(years_all, 3)
yr_late   <- tail(years_all, 3)

df_early <- df_grid %>%
  filter(pheno_year %in% yr_early) %>%
  group_by(lon, lat, Veg) %>%
  summarise(coef_early = mean(coef, na.rm = TRUE), .groups = "drop")

df_late <- df_grid %>%
  filter(pheno_year %in% yr_late) %>%
  group_by(lon, lat, Veg) %>%
  summarise(coef_late = mean(coef, na.rm = TRUE), .groups = "drop")

df_change <- inner_join(df_early, df_late, by = c("lon", "lat", "Veg")) %>%
  mutate(delta = coef_late - coef_early)

change_limit <- quantile(abs(df_change$delta), 0.98, na.rm = TRUE)


for (cls in veg_classes) {
  df_cls <- df_change %>% filter(Veg == cls)
  if (nrow(df_cls) == 0) next

  p <- ggplot(df_cls, aes(x = lon, y = lat, fill = delta)) +
    geom_tile(width = GRID_RES, height = GRID_RES) +
    scale_fill_distiller(palette = "RdBu", direction = 1,
                         limits = c(-change_limit, change_limit), oob = scales::squish,
                         name = "\u0394 fraction", na.value = "grey92") +
    map_context_layers +
    labs(title    = sprintf("%s \u2014 fraction change", tools::toTitleCase(cls)),
         subtitle = sprintf("Early: %s  \u2192  Late: %s",
                            paste(yr_early, collapse = "/"),
                            paste(yr_late,  collapse = "/"))) +
    coord_equal() +
    map_theme

  fn <- file.path(OUT_DIR, sprintf("map_change_%s.png", cls))
  ggsave(fn, p, width = 7, height = 5, dpi = DPI)
  cat(sprintf("[MAP] Saved: %s\n", fn))
}

# =============================================================================
# 3. Dominant vegetation raster map
# =============================================================================
df_dom <- df_grid %>%
  group_by(lon, lat, pheno_year) %>%
  slice_max(coef, n = 1, with_ties = FALSE) %>%
  ungroup()

df_dom_plot <- if (!is.null(SLICE_YEARS)) {
  df_dom %>% filter(pheno_year %in% SLICE_YEARS)
} else {
  df_dom
}

p_dom <- ggplot(df_dom_plot, aes(x = lon, y = lat, fill = Veg)) +
  geom_tile(width = GRID_RES, height = GRID_RES) +
  facet_wrap(~pheno_year, nrow = 2) +
  scale_fill_manual(values = dom_colors, name = "Dominant class") +
  map_context_layers +
  labs(title    = "Dominant vegetation class",
       subtitle = sprintf("Grid: %.3f\u00b0 | Years: %s",
                          GRID_RES, paste(SLICE_YEARS, collapse = ", "))) +
  coord_equal() +
  map_theme

fn_dom <- file.path(OUT_DIR, "map_dominant_veg.png")
ggsave(fn_dom, p_dom, width = 10, height = 5, dpi = DPI)
cat(sprintf("[MAP] Saved: %s\n", fn_dom))

# =============================================================================
# 4. Time-series trend: mean fraction per class over all years
# =============================================================================
# Using the grid approximation to get a true global average that is not
# weighted by unequal sampling density.
df_trend <- df_grid %>%
  group_by(pheno_year, Veg) %>%
  summarise(mean_frac = mean(coef, na.rm = TRUE),
            sd_frac   = sd(coef,   na.rm = TRUE),
            n         = n(), .groups = "drop") %>%
  mutate(se = sd_frac / sqrt(n))

trend_colors <- dom_colors[intersect(names(dom_colors), unique(df_trend$Veg))]

p_trend <- ggplot(df_trend,
                  aes(x = pheno_year, y = mean_frac, color = Veg, fill = Veg)) +
  geom_ribbon(aes(ymin = mean_frac - se, ymax = mean_frac + se),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = trend_colors, name = NULL) +
  scale_fill_manual(values  = trend_colors, name = NULL) +
  scale_x_continuous(breaks = seq(YEAR_START, YEAR_END, by = 5)) +
  labs(title = "Mean vegetation fraction over time", x = NULL, y = "Mean fraction") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")

fn_trend <- file.path(OUT_DIR, "trend_mean_fraction.png")
ggsave(fn_trend, p_trend, width = 8, height = 4, dpi = DPI)
cat(sprintf("[MAP] Saved: %s\n", fn_trend))

cat("[MAP] Done.\n")
