# gw_inference_correlation.R
# Correlate annual inference results with groundwater depths (Alagan, Yingsu)
# - inspects Excel sheets `Alagan` and `Yingsu` at runtime
# - implements delayed response (lag in years; default = 1)
# - outputs lag-wise correlations, regression stats and plots

# -------- USER CONFIG --------
inference_csv <- "inference_results/inference_results_Landsat_Harmonized_Bands_1985_2025_low_3_.csv"   # relative to project root
gw_xlsx <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/groundwater_depths (1).xlsx"
sheets_to_use <- c("Alagan", "Yingsu")
# Optional manual mapping: set sheet name -> inference `location_id` (numeric). If NULL, will use regional mean for that sheet.
# Example: sheet_loc_map <- list(Alagan = 123, Yingsu = 456)
sheet_loc_map <- list(Alagan = NULL, Yingsu = NULL)

veg_to_use <- "all"            # change to the veg/label you want to correlate (or NULL/'all' for all)
# Vegetation classes to exclude from plotting/analysis (e.g. 'barren')
exclude_vegs <- c('barren')
metric_col <- "coef"               # expected PPI-normalized column in `inference_csv`
require_ppi_normalized_trend <- TRUE  # fail fast when only known un-normalized metrics are available
# evaluate lags from min_lag_years through max_lag_years; negative values allowed
min_lag_years <- -2                # include additional lags down to -2 (GW leads by two years)
max_lag_years <- 3                   # evaluate lags 0..max_lag_years (default includes 1-year lag)
default_lag <- 1                     # used for plotting / example output
out_dir <- "gw_correlation"  # output CSVs + plots
# -----------------------------

# Required packages
pkgs <- c("readxl", "dplyr", "lubridate", "ggplot2", "broom", "tidyr", "writexl")
inst <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(inst)) install.packages(inst, repos = "https://cloud.r-project.org")

library(readxl); library(dplyr); library(lubridate); library(ggplot2)
library(broom); library(tidyr)

# define a lightweight ggplot2 theme used across MESMA scripts, with a
# fallback in case the function isn't provided elsewhere (e.g. when this
# script is run standalone.
if (!exists("theme_mesma", mode = "function")) {
  theme_mesma <- function(base_size = 9, base_family = "") {
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        axis.title = element_text(size = base_size),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
      )
  }
}

# helper: construct a vegetation color palette matching the logic used by
# fit_veg_mixture_mesma.R. This ensures correlation plots use the same base
# colours as the fitting workflow (responds to VEG_CALIBRATION_COLORS or
# falls back to RColorBrewer Set1).
build_veg_palette <- function(veg_levels) {
  pal <- NULL
  if (exists("VEG_CALIBRATION_COLORS", inherits = TRUE)) {
    supplied <- get("VEG_CALIBRATION_COLORS", inherits = TRUE)
    matched <- sapply(veg_levels, function(v) {
      nm <- names(supplied)
      im <- which(tolower(nm) == tolower(v))
      if (length(im) > 0) supplied[im[1]] else NA_character_
    }, USE.NAMES = FALSE)
    if (!all(is.na(matched))) pal <- setNames(matched, veg_levels)
  }
  if (is.null(pal) || any(is.na(pal))) {
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      stop("RColorBrewer required for default veg palette")
    }
    nveg <- max(3, length(veg_levels))
    brewer_cols <- RColorBrewer::brewer.pal(n = nveg, name = "Set1")
    brewer_cols <- brewer_cols[seq_len(length(veg_levels))]
    names(brewer_cols) <- veg_levels
    if (is.null(pal)) pal <- brewer_cols else {
      na_idx <- which(is.na(pal))
      if (length(na_idx) > 0) pal[na_idx] <- brewer_cols[na_idx]
    }
  }
  pal
}

# Resolve the metric column used for GW correlation, preferring explicit
# PPI-normalized columns when present.
resolve_inference_metric_col <- function(df, preferred_col = "coef") {
  candidates <- unique(c(
    preferred_col,
    "coef_ppi_norm", "coef_abs", "ppi_norm_coef", "coef"
  ))
  for (col in candidates) {
    if (col %in% names(df)) {
      vals <- suppressWarnings(as.numeric(df[[col]]))
      if (any(is.finite(vals))) return(col)
    }
  }
  NA_character_
}

# Guard against accidentally using relative/raw (un-normalized) trend columns.
assert_ppi_normalized_metric <- function(df, metric_col) {
  unnormalized_cols <- c("coef_raw", "coef_rel", "rel_coef", "raw_coef", "unnormalized_coef")
  present_unorm <- intersect(unnormalized_cols, names(df))
  if (length(present_unorm) > 0 && metric_col %in% present_unorm) {
    stop(sprintf(
      "Selected metric '%s' is an un-normalized trend column. Choose a PPI-normalized metric.",
      metric_col
    ))
  }

  # Variant rows can duplicate class estimates. Collapse to one value per
  # location-year-veg before trend/correlation stats.
  if (!all(c("location_id", "pheno_year", "Veg", metric_col) %in% names(df))) {
    stop("Inference data is missing required columns for metric validation")
  }
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# remove any pre-existing outputs for excluded vegetation classes (keeps output folder clean)
if(exists('exclude_vegs') && length(exclude_vegs) > 0) {
  patt <- paste0('(?i)(', paste0(exclude_vegs, collapse = '|'), ')')
  rem <- list.files(out_dir, pattern = patt, full.names = TRUE)
  if(length(rem) > 0) {
    file.remove(rem)
    message('Removed existing output files for excluded veg(s): ', paste(basename(rem), collapse = ', '))
  }
}

# Safe sheet read with informative messages (robust to locked/oneDrive files)
safe_read_sheet <- function(path, sheet) {
  # Try direct read first
  tryCatch(
    {
      df <- readxl::read_excel(path, sheet = sheet)
      message(sprintf("Read sheet '%s' -> %d rows, %d cols (direct)", sheet, nrow(df), ncol(df)))
      return(df)
    },
    error = function(e_primary) {
      message(sprintf("Direct read failed for '%s': %s", path, e_primary$message))
      # Attempt to copy to a temporary file and read the copy (helps when file is locked by Excel/OneDrive)
      tmp <- tempfile(fileext = ".xlsx")
      copied <- tryCatch({ file.copy(path, tmp, overwrite = TRUE) }, error = function(e) FALSE)
      if(!copied) {
        # look for a sibling copy (created externally) as a fallback
        path_dir <- dirname(path)
        # prefer explicit "copy" variants, otherwise pick any other groundwater_depths*.xlsx in the same folder
        candidates <- list.files(path_dir, pattern = "groundwater.*copy.*\\.xlsx$", full.names = TRUE, ignore.case = TRUE)
        if(length(candidates) == 0) {
          candidates <- list.files(path_dir, pattern = "groundwater.*depths.*\\.xlsx$", full.names = TRUE, ignore.case = TRUE)
          candidates <- setdiff(candidates, path)
        }
        alt_found <- if(length(candidates) > 0) candidates[1] else NULL
        if(!is.null(alt_found)) {
          message(sprintf("Could not create temp copy; using existing sibling copy: %s", alt_found))
          tryCatch({
            df2 <- readxl::read_excel(alt_found, sheet = sheet)
            message(sprintf("Read sheet '%s' -> %d rows, %d cols (from sibling copy)", sheet, nrow(df2), ncol(df2)))
            return(df2)
          }, error = function(e_alt) {
            stop(sprintf("Failed to read sheet '%s' from sibling copy '%s': %s", sheet, alt_found, e_alt$message))
          })
        }
        stop(sprintf("Could not copy locked file '%s' to temporary location. Close the file in Excel or provide an unlocked copy. (orig error: %s)", path, e_primary$message))
      }
      message(sprintf("Copied '%s' -> '%s' (temporary) and retrying read", path, tmp))
      tryCatch({
        df2 <- readxl::read_excel(tmp, sheet = sheet)
        message(sprintf("Read sheet '%s' -> %d rows, %d cols (from temp copy)", sheet, nrow(df2), ncol(df2)))
        return(df2)
      }, error = function(e2) {
        stop(sprintf("Failed to read sheet '%s' from temporary copy '%s': %s", sheet, tmp, e2$message))
      })
    }
  )
}

# Try to detect date/year and depth columns
detect_gw_columns <- function(df) {
  names_low <- tolower(names(df))
  date_col <- names(df)[grepl("date|day|time", names_low)][1]
  year_col <- names(df)[grepl("year|yr", names_low)][1]
  depth_col <- names(df)[grepl("depth|gwd|groundwater|waterlevel|wl|gw|level", names_low)][1]
  lat_col <- names(df)[grepl("lat|latitude", names_low)][1]
  lon_col <- names(df)[grepl("lon|long|longitude", names_low)][1]
  list(date_col = date_col, year_col = year_col, depth_col = depth_col, lat_col = lat_col, lon_col = lon_col)
}

# Convert groundwater sheet into year / mean_depth table
gw_to_yearly <- function(df) {
  cols <- detect_gw_columns(df)
  if(!is.null(cols$date_col)) {
    df <- df %>% mutate(.gw_date = as.Date(!!rlang::sym(cols$date_col)))
    df <- df %>% mutate(.gw_year = year(.data$.gw_date))
  } else if(!is.null(cols$year_col)) {
    df <- df %>% mutate(.gw_year = as.integer(!!rlang::sym(cols$year_col)))
  } else {
    stop('No date/year column detected in groundwater sheet; please provide a year or date column')
  }
  if(is.null(cols$depth_col)) stop('No depth column detected (looked for "depth|gwd|groundwater").')

  df2 <- df %>%
    filter(!is.na(.data$.gw_year)) %>%
    group_by(.gw_year) %>%
    summarize(gw_mean = mean(!!rlang::sym(cols$depth_col), na.rm = TRUE),
              gw_median = median(!!rlang::sym(cols$depth_col), na.rm = TRUE),
              n_obs = sum(!is.na(!!rlang::sym(cols$depth_col)))) %>%
    ungroup()
  df2
}

# Compute lagged correlations between annual inference and groundwater
compute_lagged_stats <- function(inf_df, gw_yearly_df, location_id, veg=NULL, metric_col = "coef", min_lag = -2, max_lag = 3) {
  inf_sub <- inf_df %>% filter(location_id == !!location_id)
  if(!is.null(veg)) inf_sub <- inf_sub %>% filter(Veg == veg)
  if(nrow(inf_sub) == 0) stop('No inference rows for chosen location_id / veg')
  if(!(metric_col %in% names(inf_sub))) stop(sprintf('Metric column "%s" not found in inference data', metric_col))

  inf_sub <- inf_sub %>% mutate(pheno_year = as.integer(pheno_year)) %>% select(location_id, pheno_year, !!rlang::sym(metric_col))
  names(inf_sub)[names(inf_sub)==metric_col] <- 'inf_metric'

  lags <- seq(min_lag, max_lag)
  out <- tibble(lag = lags, n = NA_integer_, pearson_r = NA_real_, p_value = NA_real_, lm_slope = NA_real_, lm_p = NA_real_, r2 = NA_real_)

  for(L in lags) {
    # assume groundwater at year Y-L influences pheno in year Y -> join on pheno_year == gw_year + L
    gw_lagged <- gw_yearly_df %>% mutate(pheno_year = .data$.gw_year + L)
    joined <- inner_join(inf_sub, gw_lagged, by = 'pheno_year')
    out$n[out$lag==L] <- nrow(joined)
    if(nrow(joined) >= 3) {
      ct <- cor.test(joined$inf_metric, joined$gw_mean, method = 'pearson')
      out$pearson_r[out$lag==L] <- ct$estimate
      out$p_value[out$lag==L] <- ct$p.value
      m <- lm(inf_metric ~ gw_mean, data = joined)
      s <- summary(m)
      out$lm_slope[out$lag==L] <- coef(m)[2]
      out$lm_p[out$lag==L] <- coef(s)[[4]][2]
      out$r2[out$lag==L] <- s$r.squared
    }
  }
  out
}

# Compute lagged stats when `inf_yearly_df` is already an annual time-series (pheno_year, inf_metric)
compute_lagged_stats_timeseries <- function(inf_yearly_df, gw_yearly_df, metric_col = 'inf_metric', min_lag = -2, max_lag = 3) {
  # inf_yearly_df: columns pheno_year, inf_metric
  if(!("pheno_year" %in% names(inf_yearly_df)) || !(metric_col %in% names(inf_yearly_df))) stop('inf_yearly_df must contain pheno_year and inf_metric')
  lags <- seq(min_lag, max_lag)
  out <- tibble(lag = lags, n = NA_integer_, pearson_r = NA_real_, p_value = NA_real_, lm_slope = NA_real_, lm_p = NA_real_, r2 = NA_real_)
  for(L in lags) {
    gw_lagged <- gw_yearly_df %>% mutate(pheno_year = .data$.gw_year + L)
    joined <- inner_join(inf_yearly_df %>% mutate(pheno_year = as.integer(pheno_year)), gw_lagged, by = 'pheno_year')
    out$n[out$lag==L] <- nrow(joined)
    if(nrow(joined) >= 3) {
      ct <- cor.test(joined[[metric_col]], joined$gw_mean, method = 'pearson')
      out$pearson_r[out$lag==L] <- ct$estimate
      out$p_value[out$lag==L] <- ct$p.value
      m <- lm(reformulate('gw_mean', response = metric_col), data = joined)
      s <- summary(m)
      out$lm_slope[out$lag==L] <- coef(m)[2]
      out$lm_p[out$lag==L] <- coef(s)[[4]][2]
      out$r2[out$lag==L] <- s$r.squared
    }
  }
  out
}

# ------------------- Trend helpers -------------------
# compute linear trend stats (slope, p-value, r2) for x/y columns in df
compute_trend_stats <- function(df, x_col, y_col) {
  if(!(x_col %in% names(df)) || !(y_col %in% names(df))) stop('missing columns for trend computation')
  df2 <- df %>% filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]]))
  if(nrow(df2) < 3) return(tibble(variable = y_col, slope = NA_real_, intercept = NA_real_, p_value = NA_real_, r2 = NA_real_, n = nrow(df2)))
  m <- lm(as.formula(paste(y_col, '~', x_col)), data = df2)
  s <- summary(m)
  tibble(variable = y_col,
         slope = unname(coef(m)[2]),
         intercept = unname(coef(m)[1]),
         p_value = coef(s)[2,4],
         r2 = s$r.squared,
         n = nrow(df2))
}

# Plot time-series trends for groundwater+inference (combined only)
plot_trends_for_sheet <- function(inf_yearly_df, gw_yearly_df, joined_df, sheet, id_tag, veg, out_dir, default_lag, veg_color = '#1f77b4', restrict01 = TRUE) {
  # use joined_df which already aligns years unlagged, scale groundwater to fit on same axis
  if(!is.null(joined_df) && nrow(joined_df) >= 3) {
    comb <- joined_df %>% arrange(pheno_year)
    
    # Scale groundwater to overlay on the same axis as the vegetation fraction
    min_inf <- min(comb$inf_metric, na.rm = TRUE); max_inf <- max(comb$inf_metric, na.rm = TRUE)
    min_gw  <- min(comb$gw_mean, na.rm = TRUE);  max_gw  <- max(comb$gw_mean, na.rm = TRUE)
    if(abs(max_gw - min_gw) < .Machine$double.eps) scale_factor <- 1 else scale_factor <- (max_inf - min_inf) / (max_gw - min_gw)
    
    comb <- comb %>% mutate(gw_scaled = (gw_mean - min_gw) * scale_factor + min_inf)

    # compute simple linear trend stats for annotations 
    tg_inf <- tryCatch(compute_trend_stats(comb %>% rename(year = pheno_year), 'year', 'inf_metric'), error = function(e) NULL)
    tg_gw  <- tryCatch(compute_trend_stats(comb %>% rename(year = pheno_year), 'year', 'gw_mean'), error = function(e) NULL)

    ann_txt <- c()
    if(!is.null(tg_inf)) ann_txt <- c(ann_txt, sprintf('inf slope=%.3g (p=%.2g)', tg_inf$slope, tg_inf$p_value))
    if(!is.null(tg_gw))  ann_txt <- c(ann_txt, sprintf('gw slope=%.3g (p=%.2g)', tg_gw$slope, tg_gw$p_value))

    y_pad <- (max_inf - min_inf) * 0.08
    p_comb <- ggplot(comb, aes(x = pheno_year)) +
      geom_line(aes(y = inf_metric, color = 'inf')) + geom_point(aes(y = inf_metric, color = 'inf')) +
      geom_line(aes(y = gw_scaled, color = 'gw'), linetype = 'dashed') + geom_point(aes(y = gw_scaled, color = 'gw')) +
      geom_smooth(aes(y = inf_metric), method = 'lm', se = FALSE, color = veg_color) +
      geom_smooth(aes(y = gw_scaled), method = 'lm', se = FALSE, color = '#d62728', linetype = 'dashed') +
      coord_cartesian(ylim = c(min_inf - y_pad, max_inf + y_pad * 4)) +
      scale_y_continuous(
        name = 'vegetation fraction',
        sec.axis = sec_axis(~ (. - min_inf)/scale_factor + min_gw, name = 'groundwater depth (m)')
      ) +
      scale_color_manual(name = '', values = c('inf' = veg_color, 'gw' = '#d62728'), labels = c('inference', 'gw')) +
      scale_x_continuous(limits = c(1984, NA)) +
      labs(x = 'year') +
      annotate('text', x = min(comb$pheno_year, na.rm = TRUE), y = max_inf + y_pad * 3.5, hjust = 0, label = paste(ann_txt, collapse = '\n'), size = 3) +
      theme_mesma()

    fc <- file.path(out_dir, sprintf('trend_combined_%s_%s_%s.png', sheet, id_tag, veg))
    ggsave(filename = fc, plot = p_comb, width = 7, height = 4)
    message('Saved combined trend plot: ', fc)
  } else {
    message('Skipping combined trend plot (insufficient joined data) for: ', sheet, ' / ', veg)
  }
}

# Plot dual-axis time series (veg fraction vs GW level)
# veg_color: colour to use for the vegetation series (should match fit script palette)
plot_dual_axis_time_series <- function(inf_yearly_df, gw_yearly_df, sheet, id_tag, veg, out_dir, veg_color = '#1f77b4') {
  if(is.null(inf_yearly_df) || nrow(inf_yearly_df) < 2 || is.null(gw_yearly_df) || nrow(gw_yearly_df) < 2) {
    message('Skipping dual-axis plot (insufficient data): ', sheet, ' / ', veg)
    return(invisible(NULL))
  }
  # align years
  gw_df <- gw_yearly_df %>% rename(gw_year = .gw_year)
  joined <- inner_join(inf_yearly_df, gw_df, by = c('pheno_year' = 'gw_year'))
  if(nrow(joined) < 2) { message('Not enough overlapping years for dual-axis plot'); return(invisible(NULL)) }

  min_inf <- min(joined$inf_metric, na.rm = TRUE); max_inf <- max(joined$inf_metric, na.rm = TRUE)
  min_gw <- min(joined$gw_mean, na.rm = TRUE); max_gw <- max(joined$gw_mean, na.rm = TRUE)
  # avoid zero division
  if(abs(max_gw - min_gw) < .Machine$double.eps) scale_factor <- 1 else scale_factor <- (max_inf - min_inf) / (max_gw - min_gw)
  # scaled groundwater for plotting on same primary axis
  joined <- joined %>% mutate(gw_scaled = (gw_mean - min_gw) * scale_factor + min_inf)

  p <- ggplot(joined, aes(x = pheno_year)) +
    geom_line(aes(y = inf_metric, color = 'veg'), linewidth = 0.9) + geom_point(aes(y = inf_metric, color = 'veg')) +
    geom_line(aes(y = gw_scaled, color = 'gw'), linewidth = 0.9, linetype = 'dashed') + geom_point(aes(y = gw_scaled, color = 'gw')) +
    scale_y_continuous(name = 'vegetation fraction', limits = c(0,1), sec.axis = sec_axis(~ (. - min_inf)/scale_factor + min_gw, name = 'groundwater depth (m)')) +
    scale_color_manual('', values = c('veg' = veg_color, 'gw' = '#d62728'), labels = c('vegetation', 'groundwater')) +
    scale_x_continuous(limits = c(1984, NA)) +
    labs(x = 'year') +
    theme_mesma()
  invisible(joined)
}

# Scatter GW (year t) vs vegetation (year t+lag) — default lag = 1
plot_scatter_lagged <- function(inf_yearly_df, gw_yearly_df, lag = 1, sheet, id_tag, veg, out_dir, veg_color = '#1f77b4', restrict01 = TRUE) {
  if(is.null(inf_yearly_df) || nrow(inf_yearly_df) < 2 || is.null(gw_yearly_df) || nrow(gw_yearly_df) < 2) {
    message('Skipping lagged scatter (insufficient data): ', sheet, ' / ', veg)
    return(invisible(NULL))
  }
  gw_lagged <- gw_yearly_df %>% mutate(pheno_year = .gw_year + lag)
  joined <- inner_join(inf_yearly_df, gw_lagged, by = 'pheno_year')
  if(nrow(joined) < 3) { message('Not enough paired years for lagged scatter (need >=3)'); return(invisible(NULL)) }

  ct <- cor.test(joined$gw_mean, joined$inf_metric, method = 'pearson')
  lm_mod <- lm(inf_metric ~ gw_mean, data = joined)
  p <- ggplot(joined, aes(x = gw_mean, y = inf_metric)) +
    geom_point(color = veg_color) + geom_smooth(method = 'lm', se = TRUE, color = veg_color) +
    { if(restrict01) scale_y_continuous(limits = c(0,1)) else scale_y_continuous() } +
    labs(x = 'GW mean (year t)', y = 'vegetation fraction (year t+lag)') +
    annotate('text', x = Inf, y = Inf, hjust = 1, vjust = 1, size = 3,
             label = sprintf('r=%.2f, p=%.2g', ct$estimate, ct$p.value)) +
    theme_mesma()
  # CSV suppressed by user request — lagged pair table not written to CSV
  message('(CSV suppressed) lagged pairs would have been (not written)')
  invisible(joined)
}

# ------------------- Run analysis -------------------
message('Reading inference results...')
inf <- readr::read_csv(inference_csv, show_col_types = FALSE)

metric_col_resolved <- resolve_inference_metric_col(inf, preferred_col = metric_col)
if (is.na(metric_col_resolved)) {
  stop("No usable PPI-normalized inference metric column found in inference CSV")
}
metric_col <- metric_col_resolved
assert_ppi_normalized_metric(inf, metric_col)
if (isTRUE(require_ppi_normalized_trend) && metric_col != "coef" && grepl("raw|rel|unnormal", metric_col, ignore.case = TRUE)) {
  stop(sprintf("Resolved metric '%s' is not PPI-normalized.", metric_col))
}
message(sprintf("Using inference metric column for GW trend/correlation: %s", metric_col))

# Use per-location-year-class medians to remove variant duplication before
# any trend/correlation aggregation.
if (all(c("location_id", "pheno_year", "lat", "lon", "Veg", metric_col) %in% names(inf))) {
  inf <- inf %>%
    mutate(.metric_value = suppressWarnings(as.numeric(.data[[metric_col]]))) %>%
    group_by(location_id, pheno_year, lat, lon, Veg) %>%
    summarize(.metric_value = median(.metric_value, na.rm = TRUE), .groups = "drop") %>%
    rename(!!metric_col := .metric_value)
}

# restrict satellite/inference data to years after 1984
if('pheno_year' %in% names(inf)) {
  inf <- inf %>% mutate(pheno_year = as.integer(pheno_year)) %>% filter(pheno_year >= 1984)
  message(sprintf('Filtered inference to pheno_year >= 1984; remaining rows: %d', nrow(inf)))
}
message(sprintf('Inference rows: %d, columns: %d', nrow(inf), ncol(inf)))

if(!is.null(veg_to_use)) {
  if(!('Veg' %in% names(inf))) stop('`Veg` column not found in inference CSV')
  # only warn when the user specified a concrete veg that doesn't exist;
  # the special value 'all' (case-insensitive) is handled later and should
  # not trigger a warning even though no row equals "all".
  if(!identical(tolower(veg_to_use), 'all') && !any(inf$Veg == veg_to_use)) {
    warning(sprintf('veg "%s" not found in inference data — results may be empty', veg_to_use))
  }
}

# For each groundwater sheet: inspect structure, convert to yearly, try to match, compute lags
# Each sheet is processed separately.
# Users who want to analyse both should include them individually in `sheets_to_use`.
final_reports <- list()
corr_list <- list()   # accumulate all lagged correlation tables for final Excel
for(sheet in sheets_to_use) {
  message('\n---- Sheet: ', sheet, ' ----')
  gw_raw <- tryCatch(safe_read_sheet(gw_xlsx, sheet), error = function(e) { message(e$message); return(NULL) })
  if(is.null(gw_raw)) next
  gw_yearly <- tryCatch({ gw_to_yearly(gw_raw) }, error = function(e) { message('Could not convert groundwater to yearly: ', e$message); return(NULL) })
  if(is.null(gw_yearly)) next

  message('Columns (first 12): ', paste(head(names(gw_raw), 12), collapse = ', '))
  print(utils::head(gw_raw, 6))

  # detect columns and convert to yearly
  cols <- detect_gw_columns(gw_raw)
  message('Detected - date_col: ', cols$date_col %||% 'NONE', '; year_col: ', cols$year_col %||% 'NONE', '; depth_col: ', cols$depth_col %||% 'NONE')

  gw_yearly <- tryCatch({ gw_to_yearly(gw_raw) }, error = function(e) { message('Could not convert groundwater to yearly: ', e$message); return(NULL) })
  if(is.null(gw_yearly)) next
  gw_yearly <- gw_yearly %>% rename(.gw_year = .data$.gw_year)  # keep consistent name inside function

  # Determine analysis mode: per-location (mapping) OR regional timeseries
  mapped_loc <- NULL
  if(!is.null(sheet_loc_map) && !is.null(sheet_loc_map[[sheet]])) mapped_loc <- sheet_loc_map[[sheet]]

  # decide chosen location id (if any) and id base tag
  chosen_loc <- NULL
  if(!is.null(mapped_loc)) {
    chosen_loc <- mapped_loc
    message(sprintf('Using manual mapping: location_id = %s for sheet %s', chosen_loc, sheet))
    id_base <- paste0('loc', chosen_loc)
  } else {
    id_base <- 'regional_mean'
    message(sprintf('No manual mapping for sheet %s — will use regional mean of inference metric', sheet))
  }

  # determine list of vegetation targets to evaluate
  if(is.null(veg_to_use) || identical(veg_to_use, 'all')) {
    veg_list <- sort(unique(inf$Veg))
  } else if(is.character(veg_to_use) && length(veg_to_use) >= 1) {
    veg_list <- veg_to_use
  } else {
    veg_list <- c(veg_to_use)
  }
  # exclude any vegetation classes configured in `exclude_vegs`
  if(exists('exclude_vegs') && length(exclude_vegs) > 0) {
    skip <- intersect(veg_list, exclude_vegs)
    if(length(skip) > 0) {
      message('Skipping excluded veg(s): ', paste(skip, collapse = ', '))
      veg_list <- setdiff(veg_list, skip)
    }
  }
  # build palette for all vegs we will loop over
  veg_palette <- build_veg_palette(veg_list)

  veg_results <- list()
  # prepare common label formatter for lag axis (numeric only)
  label_fun <- function(x) {
    as.character(x)
  }

  for(veg in veg_list) {
    message(sprintf('\nProcessing sheet=%s, veg=%s, id=%s', sheet, veg, id_base))
    # pick colour for this vegetation class (fallback to blue if somehow missing)
    veg_color <- if (!is.null(veg_palette) && veg %in% names(veg_palette)) veg_palette[[veg]] else '#1f77b4'
    if(!is.null(chosen_loc)) {
      stats_tbl <- compute_lagged_stats(inf, gw_yearly, chosen_loc, veg = veg, metric_col = metric_col,
                                       min_lag = min_lag_years, max_lag = max_lag_years)
      joined_def <- inner_join(inf %>% filter(location_id==chosen_loc, Veg==veg) %>% mutate(pheno_year = as.integer(pheno_year)) %>% select(pheno_year, inf_metric = !!rlang::sym(metric_col)), gw_yearly %>% mutate(pheno_year = .data$.gw_year + default_lag), by='pheno_year')
      inf_yearly_for_trend <- inf %>% filter(location_id==chosen_loc, Veg==veg) %>% group_by(pheno_year) %>% summarize(inf_metric = mean(!!rlang::sym(metric_col), na.rm = TRUE), n_sites = n()) %>% ungroup()
      id_tag <- paste0(id_base, '_', veg)
    } else {
      # regional aggregated timeseries for this veg
      inf_yearly_for_trend <- inf %>% filter(Veg == veg) %>% group_by(pheno_year) %>% summarize(inf_metric = mean(!!rlang::sym(metric_col), na.rm = TRUE), n_sites = n()) %>% ungroup()
      stats_tbl <- compute_lagged_stats_timeseries(inf_yearly_for_trend, gw_yearly, metric_col = 'inf_metric',
                                                     min_lag = min_lag_years, max_lag = max_lag_years)
      joined_def <- inner_join(inf_yearly_for_trend %>% mutate(pheno_year = as.integer(pheno_year)), gw_yearly %>% mutate(pheno_year = .data$.gw_year + default_lag), by='pheno_year')
      id_tag <- paste0(id_base, '_', veg)
    }

    stats_tbl <- stats_tbl %>%
      arrange(lag)
    # sanity check: ensure requested minimum lag is actually present
    if(exists('min_lag_years') && min(stats_tbl$lag, na.rm = TRUE) > min_lag_years) {
      warning(sprintf('Computed lags did not include min_lag_years (%s); check configuration', min_lag_years))
    }
    out_csv <- file.path(out_dir, sprintf('gw_correlation_%s_%s.csv', sheet, id_tag))
    # CSV suppressed by user request — lagged summary not written to CSV
    message('(CSV suppressed) lagged summary would have been: ', out_csv)

    # append to master correlations list for Excel export (include sheet/id/veg)
    corr_list[[length(corr_list) + 1]] <- stats_tbl %>% mutate(sheet = sheet, id = id_tag, veg = veg)

    # plot correlation vs lag (report lagged correlations numerically)
    p1 <- ggplot(stats_tbl, aes(x = lag, y = pearson_r)) +
        geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey60', linewidth = 0.4) +
        geom_line(color = veg_color, linewidth = 0.9) +
        geom_point(size = 3.5, color = veg_color) +
        # only label lags other than -1 and 0
        geom_text(
          data = subset(stats_tbl, !(lag %in% c(-1, 0))),
          aes(label = ifelse(n >= 3, sprintf('r=%.2f\np=%.3f', pearson_r, p_value), 'n<3')),
          vjust = -0.6, size = 2.7, lineheight = 0.9
        ) +
        scale_y_continuous(limits = c(-1, 1), expand = expansion(mult = c(0.05, 0.28))) +
        coord_cartesian(clip = 'off') +
        scale_x_continuous(breaks = stats_tbl$lag, labels = label_fun) +
        labs(x = 'Lag (years)', y = 'Pearson r',
             title = sprintf('%s — %s', sheet, veg)) +
        theme_mesma() +
        theme(
          plot.title  = element_text(size = 9, face = 'bold', margin = margin(b = 6)),
          plot.margin = margin(t = 16, r = 10, b = 8, l = 10),
          panel.grid.major.x = element_blank()
        )
    f1 <- file.path(out_dir, sprintf('r_vs_lag_%s_%s_%s.png', sheet, id_tag, veg))
    ggsave(filename = f1, plot = p1, width = 6, height = 4)
    message('Saved plot: ', f1)

    # create UNLAGGED joined timeseries for plotting and z-scored comparisons
    joined_unlagged <- if(nrow(inf_yearly_for_trend) > 0) inner_join(inf_yearly_for_trend %>% mutate(pheno_year = as.integer(pheno_year)), gw_yearly %>% mutate(pheno_year = .gw_year), by = 'pheno_year') else tibble()
    if(nrow(joined_unlagged) > 0) {
      jf_unlag <- file.path(out_dir, sprintf('joined_timeseries_%s_unlagged_%s_%s.csv', sheet, id_tag, veg))
      # CSV suppressed by user request — unlagged joined timeseries not written to CSV
      message('(CSV suppressed) joined timeseries (unlagged) would have been: ', jf_unlag)
    }
    # compute general (no-lag) correlation p-value for this veg
    if(nrow(joined_unlagged) >= 3) {
      ct_gen <- cor.test(joined_unlagged$inf_metric, joined_unlagged$gw_mean, method = 'pearson')
      p_gen <- ct_gen$p.value
    } else {
      p_gen <- NA_real_
    }

    # UNLAGGED scatter plot (GW year t vs veg year t) — plots must be unlagged per request
    tryCatch({
      if(nrow(joined_unlagged) >= 3) {
        plot_scatter_lagged(inf_yearly_for_trend, gw_yearly, lag = 0, sheet = sheet, id_tag = id_tag, veg = veg, out_dir = out_dir, veg_color = veg_color)
      } else {
        message('Not enough paired years for unlagged scatter (need >=3) for veg ', veg)
      }
    }, error = function(e) message('Unlagged scatter plotting skipped for ', sheet, ': ', e$message))

    # Combined z-scored trend plot (unlagged pairing)
    tryCatch({
      if(!grepl("Alagan_Yingsu", sheet)) {
        plot_trends_for_sheet(if(nrow(inf_yearly_for_trend)>0) inf_yearly_for_trend else NULL, gw_yearly, if(nrow(joined_unlagged)>0) joined_unlagged else NULL, sheet, id_tag, veg, out_dir, default_lag, veg_color = veg_color)
        message('Saved combined trend plot for: ', paste0(sheet, ' / ', id_tag, ' / ', veg))
      } else {
        message('Skipping trend plot for combined sheet: ', sheet)
      }
    }, error = function(e) message('Trend plotting skipped for ', sheet, ': ', e$message))

    # dual-axis time series: vegetation fraction (primary) and groundwater (secondary) — UNLAGGED
    tryCatch({
      plot_dual_axis_time_series(if(nrow(inf_yearly_for_trend)>0) inf_yearly_for_trend else NULL, gw_yearly, sheet, id_tag, veg, out_dir, veg_color = veg_color)
    }, error = function(e) message('Dual-axis plot skipped for ', sheet, ': ', e$message))

    veg_results[[veg]] <- list(stats = stats_tbl, joined_unlagged = joined_unlagged, p_overall = p_gen)
  }

  # after processing all vegs for this sheet, create combined overlay plots
  if(length(veg_results) > 0 && !grepl("Alagan_Yingsu", sheet)) {
    combined <- bind_rows(lapply(veg_results, function(x) x$stats), .id = 'veg')
    if(nrow(combined) > 0) {
      p_all <- ggplot(combined, aes(x = lag, y = pearson_r, color = veg, group = veg)) +
        geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey60', linewidth = 0.4) +
        geom_line(linewidth = 0.9) +
        geom_point(size = 3) +
        # suppress labels at lag -1 and 0
        geom_text(
          data = subset(combined, !(lag %in% c(-1, 0)) & !is.na(pearson_r)),
          aes(label = sprintf("r=%.2f\np=%.3f", pearson_r, p_value)),
          vjust = -0.6, hjust = 0.5, size = 2.6, lineheight = 0.9, show.legend = FALSE
        ) +
        scale_color_manual(values = veg_palette) +
        scale_y_continuous(limits = c(-1, 1), expand = expansion(mult = c(0.05, 0.28))) +
        scale_x_continuous(breaks = unique(combined$lag), labels = label_fun(unique(combined$lag))) +
        labs(x = 'Lag (years)', y = 'Pearson r', title = paste('All vegetation classes —', sheet)) +
        theme_mesma() +
        theme(
          plot.title       = element_text(size = 9, face = 'bold', margin = margin(b = 6)),
          plot.margin      = margin(t = 16, r = 10, b = 8, l = 10),
          panel.grid.major.x = element_blank(),
          legend.title     = element_blank()
        )
      fp <- file.path(out_dir, sprintf('r_vs_lag_%s_allvegs.png', sheet))
      ggsave(fp, plot = p_all, width = 6, height = 4)
      message('Saved combined correlation overview: ', fp)
    }
  }

  final_reports[[sheet]] <- list(chosen_location = if(!is.null(chosen_loc)) chosen_loc else id_base, veg_results = veg_results)
}

# Print concise summary
message('\nSummary:')
for(sheet in names(final_reports)) {
  fr <- final_reports[[sheet]]
  for(veg in names(fr$veg_results)) {
    tr <- fr$veg_results[[veg]]$stats
    unlag <- tr %>% filter(lag == 0) %>% slice(1)
    best <- tr %>% filter(!is.na(pearson_r)) %>% arrange(desc(abs(pearson_r))) %>% slice(1)
    u_r <- if(nrow(unlag) == 1) unlag$pearson_r else NA_real_
    u_p <- if(nrow(unlag) == 1) unlag$p_value else NA_real_
    u_n <- if(nrow(unlag) == 1) unlag$n else NA_integer_
    overall_p <- fr$veg_results[[veg]]$p_overall
    b_lag <- if(nrow(best) == 1) best$lag       else NA_integer_
    b_r   <- if(nrow(best) == 1) best$pearson_r else NA_real_
    b_p   <- if(nrow(best) == 1) best$p_value   else NA_real_
    b_n   <- if(nrow(best) == 1) best$n         else NA_integer_
    message(sprintf("Sheet %s / veg %s -> id=%s | unlagged (lag=0): r=%.3f p=%.3g n=%d; overall p=%.3g; best lag=%d: r=%.3f p=%.3g n=%d",
                    sheet, veg, fr$chosen_location, u_r, u_p, u_n, overall_p,
                    b_lag, b_r, b_p, b_n))
  }
}

# write ALL correlations (lagged and unlagged rows) to a single Excel file
if(length(corr_list) > 0) {
  all_corr <- bind_rows(corr_list) %>% select(sheet, id, veg, everything())
  xlsx_out <- file.path(out_dir, 'all_correlations.xlsx')
  writexl::write_xlsx(all_corr, xlsx_out)
  message('Saved Excel with all correlations (coefficients + p-values): ', xlsx_out)
} else {
  message('No correlation tables to write to Excel.')
}


# list generated files for user convenience
files <- list.files(out_dir, full.names = TRUE)
message('\nGenerated files:')
for(f in files) message(' - ', f)

message('\nOutputs written to: ', normalizePath(out_dir))
message('Done.')

# EOF