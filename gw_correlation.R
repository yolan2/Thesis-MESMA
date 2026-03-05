# gw_inference_correlation.R
# Correlate annual inference results with groundwater depths (Alagan, Yingsu)
# - inspects Excel sheets `Alagan` and `Yingsu` at runtime
# - implements delayed response (lag in years; default = 1)
# - outputs lag-wise correlations, regression stats and plots

# -------- USER CONFIG --------
inference_csv <- "inference_results/inference_results.csv"   # relative to project root
gw_xlsx <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/groundwater_depths (1).xlsx"
sheets_to_use <- c("Alagan", "Yingsu")
# Optional manual mapping: set sheet name -> inference `location_id` (numeric). If NULL, will use regional mean for that sheet.
# Example: sheet_loc_map <- list(Alagan = 123, Yingsu = 456)
sheet_loc_map <- list(Alagan = NULL, Yingsu = NULL)

veg_to_use <- "all"            # change to the veg/label you want to correlate (or NULL/'all' for all)
# Vegetation classes to exclude from plotting/analysis (e.g. 'barren')
exclude_vegs <- c('barren')
metric_col <- "coef"               # column in `inference_csv` to correlate (must be numeric)
max_lag_years <- 3                   # evaluate lags 0..max_lag_years (default includes 1-year lag)
default_lag <- 1                     # used for plotting / example output
out_dir <- "temp_results/gw_correlation"  # output CSVs + plots
# -----------------------------

# Required packages
pkgs <- c("readxl", "dplyr", "lubridate", "ggplot2", "broom", "tidyr", "writexl")
inst <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(inst)) install.packages(inst, repos = "https://cloud.r-project.org")

library(readxl); library(dplyr); library(lubridate); library(ggplot2)
library(broom); library(tidyr)

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
compute_lagged_stats <- function(inf_df, gw_yearly_df, location_id, veg=NULL, metric_col = "coef", max_lag = 3) {
  inf_sub <- inf_df %>% filter(location_id == !!location_id)
  if(!is.null(veg)) inf_sub <- inf_sub %>% filter(Veg == veg)
  if(nrow(inf_sub) == 0) stop('No inference rows for chosen location_id / veg')
  if(!(metric_col %in% names(inf_sub))) stop(sprintf('Metric column "%s" not found in inference data', metric_col))

  inf_sub <- inf_sub %>% mutate(pheno_year = as.integer(pheno_year)) %>% select(location_id, pheno_year, !!rlang::sym(metric_col))
  names(inf_sub)[names(inf_sub)==metric_col] <- 'inf_metric'

  out <- tibble(lag = 0:max_lag, n = NA_integer_, pearson_r = NA_real_, p_value = NA_real_, lm_slope = NA_real_, lm_p = NA_real_, r2 = NA_real_)

  for(L in 0:max_lag) {
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
compute_lagged_stats_timeseries <- function(inf_yearly_df, gw_yearly_df, metric_col = 'inf_metric', max_lag = 3) {
  # inf_yearly_df: columns pheno_year, inf_metric
  if(!("pheno_year" %in% names(inf_yearly_df)) || !(metric_col %in% names(inf_yearly_df))) stop('inf_yearly_df must contain pheno_year and inf_metric')
  out <- tibble(lag = 0:max_lag, n = NA_integer_, pearson_r = NA_real_, p_value = NA_real_, lm_slope = NA_real_, lm_p = NA_real_, r2 = NA_real_)
  for(L in 0:max_lag) {
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
# Note: single-variable 'mean' trend plots have been removed per request.
plot_trends_for_sheet <- function(inf_yearly_df, gw_yearly_df, joined_df, sheet, id_tag, veg, out_dir, default_lag) {
  # combined z-scored series (use joined_df which already aligns years unlagged)
  if(!is.null(joined_df) && nrow(joined_df) >= 3) {
    comb <- joined_df %>% arrange(pheno_year) %>% mutate(z_inf = as.numeric(scale(inf_metric)), z_gw = as.numeric(scale(gw_mean)))
    # compute simple linear trend stats for annotations (do not save separate 'mean' trend plots)
    tg_inf <- tryCatch(compute_trend_stats(comb %>% rename(year = pheno_year), 'year', 'inf_metric'), error = function(e) NULL)
    tg_gw  <- tryCatch(compute_trend_stats(comb %>% rename(year = pheno_year), 'year', 'gw_mean'), error = function(e) NULL)

    ann_txt <- c()
    if(!is.null(tg_inf)) ann_txt <- c(ann_txt, sprintf('inf slope=%.3g (p=%.2g)', tg_inf$slope, tg_inf$p_value))
    if(!is.null(tg_gw))  ann_txt <- c(ann_txt, sprintf('gw slope=%.3g (p=%.2g)', tg_gw$slope, tg_gw$p_value))

    p_comb <- ggplot(comb, aes(x = pheno_year)) +
      geom_line(aes(y = z_inf, color = 'inf')) + geom_point(aes(y = z_inf, color = 'inf')) +
      geom_line(aes(y = z_gw, color = 'gw')) + geom_point(aes(y = z_gw, color = 'gw')) +
      geom_smooth(aes(y = z_inf), method = 'lm', se = FALSE, color = '#1f77b4') +
      geom_smooth(aes(y = z_gw), method = 'lm', se = FALSE, color = '#d62728') +
      scale_color_manual(name = '', values = c('inf' = '#1f77b4', 'gw' = '#d62728'), labels = c('inference (z)', 'gw (z)')) +
      labs(x = 'year', y = 'z-score') +
      annotate('text', x = min(comb$pheno_year, na.rm = TRUE), y = max(c(comb$z_inf, comb$z_gw), na.rm = TRUE), hjust = 0, label = paste(ann_txt, collapse = '; '), size = 3) +
      theme_mesma()

    fc <- file.path(out_dir, sprintf('trend_combined_%s_%s_%s.png', sheet, id_tag, veg))
    ggsave(filename = fc, plot = p_comb, width = 7, height = 4)
    message('Saved combined trend plot: ', fc)
  } else {
    message('Skipping combined trend plot (insufficient joined data) for: ', sheet, ' / ', veg)
  }
}

# Plot dual-axis time series (veg fraction vs GW level)
plot_dual_axis_time_series <- function(inf_yearly_df, gw_yearly_df, sheet, id_tag, veg, out_dir) {
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
    scale_y_continuous(name = 'vegetation fraction', sec.axis = sec_axis(~ (. - min_inf)/scale_factor + min_gw, name = 'groundwater depth (m)')) +
    scale_color_manual('', values = c('veg' = '#1f77b4', 'gw' = '#d62728'), labels = c('vegetation', 'groundwater')) +
    labs(x = 'year') +
    theme_mesma()
  # CSV suppressed by user request — joined timeseries not written to CSV
  message('(CSV suppressed) joined timeseries would have been: ', jf)
  invisible(joined)
}

# Scatter GW (year t) vs vegetation (year t+lag) — default lag = 1
plot_scatter_lagged <- function(inf_yearly_df, gw_yearly_df, lag = 1, sheet, id_tag, veg, out_dir) {
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
    geom_point() + geom_smooth(method = 'lm', se = TRUE) +
    labs(x = 'GW mean (year t)', y = 'vegetation fraction (year t+lag)') +
    annotate('text', x = Inf, y = Inf, hjust = 1, vjust = 1, size = 3,
             label = sprintf('r=%.2f, p=%.2g, n=%d', ct$estimate, ct$p.value, nrow(joined))) +
    theme_mesma()
  # CSV suppressed by user request — lagged pair table not written to CSV
  message('(CSV suppressed) lagged pairs would have been (not written)')
  invisible(joined)
}

# ------------------- Run analysis -------------------
message('Reading inference results...')
inf <- readr::read_csv(inference_csv, show_col_types = FALSE)
message(sprintf('Inference rows: %d, columns: %d', nrow(inf), ncol(inf)))

if(!is.null(veg_to_use)) {
  if(!('Veg' %in% names(inf))) stop('`Veg` column not found in inference CSV')
  if(!any(inf$Veg == veg_to_use)) warning(sprintf('veg "%s" not found in inference data — results may be empty', veg_to_use))
}

# For each groundwater sheet: inspect structure, convert to yearly, try to match, compute lags
# If both Alagan and Yingsu are requested, treat them as a single combined groundwater series
if(all(c('Alagan','Yingsu') %in% sheets_to_use)) {
  message('Treating Alagan and Yingsu as the same location — combining their groundwater series')
  sheets_to_use <- unique(c(setdiff(sheets_to_use, c('Alagan','Yingsu')), 'Alagan_Yingsu'))
}
final_reports <- list()
corr_list <- list()   # accumulate all lagged correlation tables for final Excel
for(sheet in sheets_to_use) {
  message('\n---- Sheet: ', sheet, ' ----')
  if(sheet == 'Alagan_Yingsu') {
    # read both sheets and combine into an annual series
    gw1 <- tryCatch(safe_read_sheet(gw_xlsx, 'Alagan'), error = function(e) { message('Alagan read failed: ', e$message); NULL })
    gw2 <- tryCatch(safe_read_sheet(gw_xlsx, 'Yingsu'), error = function(e) { message('Yingsu read failed: ', e$message); NULL })
    if(is.null(gw1) && is.null(gw2)) { message('Neither Alagan nor Yingsu could be read; skipping'); next }
    gw1y <- if(!is.null(gw1)) tryCatch(gw_to_yearly(gw1), error = function(e) { message('gw_to_yearly failed for Alagan: ', e$message); NULL }) else NULL
    gw2y <- if(!is.null(gw2)) tryCatch(gw_to_yearly(gw2), error = function(e) { message('gw_to_yearly failed for Yingsu: ', e$message); NULL }) else NULL
    gw_yearly <- bind_rows(gw1y, gw2y) %>% group_by(.gw_year) %>% summarize(gw_mean = mean(gw_mean, na.rm = TRUE), gw_median = mean(gw_median, na.rm = TRUE), n_wells = n()) %>% ungroup()
    message(sprintf('Combined Alagan+Yingsu -> %d annual rows (years %s-%s)', nrow(gw_yearly), min(gw_yearly$.gw_year, na.rm=TRUE), max(gw_yearly$.gw_year, na.rm=TRUE)))
  } else {
    gw_raw <- tryCatch(safe_read_sheet(gw_xlsx, sheet), error = function(e) { message(e$message); return(NULL) })
    if(is.null(gw_raw)) next
    gw_yearly <- tryCatch({ gw_to_yearly(gw_raw) }, error = function(e) { message('Could not convert groundwater to yearly: ', e$message); return(NULL) })
    if(is.null(gw_yearly)) next
  }

  if(sheet != 'Alagan_Yingsu') {
    message('Columns (first 12): ', paste(head(names(gw_raw), 12), collapse = ', '))
    print(utils::head(gw_raw, 6))

    # detect columns and convert to yearly
    cols <- detect_gw_columns(gw_raw)
    message('Detected - date_col: ', cols$date_col %||% 'NONE', '; year_col: ', cols$year_col %||% 'NONE', '; depth_col: ', cols$depth_col %||% 'NONE')

    gw_yearly <- tryCatch({ gw_to_yearly(gw_raw) }, error = function(e) { message('Could not convert groundwater to yearly: ', e$message); return(NULL) })
    if(is.null(gw_yearly)) next
    gw_yearly <- gw_yearly %>% rename(.gw_year = .data$.gw_year)  # keep consistent name inside function
  } else {
    # combined Alagan+Yingsu: no single `gw_raw` available — set safe `cols` (no lat/lon)
    cols <- list(date_col = NULL, year_col = '.gw_year', depth_col = 'gw_mean', lat_col = NULL, lon_col = NULL)
  }

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

  veg_results <- list()
  for(veg in veg_list) {
    message(sprintf('\nProcessing sheet=%s, veg=%s, id=%s', sheet, veg, id_base))
    if(!is.null(chosen_loc)) {
      stats_tbl <- compute_lagged_stats(inf, gw_yearly, chosen_loc, veg = veg, metric_col = metric_col, max_lag = max_lag_years)
      joined_def <- inner_join(inf %>% filter(location_id==chosen_loc, Veg==veg) %>% mutate(pheno_year = as.integer(pheno_year)) %>% select(pheno_year, inf_metric = !!rlang::sym(metric_col)), gw_yearly %>% mutate(pheno_year = .data$.gw_year + default_lag), by='pheno_year')
      inf_yearly_for_trend <- inf %>% filter(location_id==chosen_loc, Veg==veg) %>% group_by(pheno_year) %>% summarize(inf_metric = mean(!!rlang::sym(metric_col), na.rm = TRUE), n_sites = n()) %>% ungroup()
      id_tag <- paste0(id_base, '_', veg)
    } else {
      # regional aggregated timeseries for this veg
      inf_yearly_for_trend <- inf %>% filter(Veg == veg) %>% group_by(pheno_year) %>% summarize(inf_metric = mean(!!rlang::sym(metric_col), na.rm = TRUE), n_sites = n()) %>% ungroup()
      stats_tbl <- compute_lagged_stats_timeseries(inf_yearly_for_trend, gw_yearly, metric_col = 'inf_metric', max_lag = max_lag_years)
      joined_def <- inner_join(inf_yearly_for_trend %>% mutate(pheno_year = as.integer(pheno_year)), gw_yearly %>% mutate(pheno_year = .data$.gw_year + default_lag), by='pheno_year')
      id_tag <- paste0(id_base, '_', veg)
    }

    stats_tbl <- stats_tbl %>% arrange(lag)
    out_csv <- file.path(out_dir, sprintf('gw_correlation_%s_%s.csv', sheet, id_tag))
    # CSV suppressed by user request — lagged summary not written to CSV
    message('(CSV suppressed) lagged summary would have been: ', out_csv)

    # append to master correlations list for Excel export (include sheet/id/veg)
    corr_list[[length(corr_list) + 1]] <- stats_tbl %>% mutate(sheet = sheet, id = id_tag, veg = veg)

    # plot correlation vs lag (report lagged correlations numerically)
    p1 <- ggplot(stats_tbl, aes(x = lag, y = pearson_r)) +
      geom_point(size = 3) + geom_line() +
      # place labels slightly closer to the point and give extra vertical room so they are not clipped
      geom_text(aes(label = ifelse(n >= 3, sprintf('r=%.2f\\n n=%d', pearson_r, n), 'n<3')), vjust = -0.5, size = 3) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
      coord_cartesian(clip = 'off') +
      labs(x = 'lag (years)', y = 'Pearson r') +
      theme_minimal() +
      theme(plot.margin = margin(t = 14, r = 8, b = 8, l = 8))
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

    # UNLAGGED scatter plot (GW year t vs veg year t) — plots must be unlagged per request
    tryCatch({
      if(nrow(joined_unlagged) >= 3) {
        plot_scatter_lagged(inf_yearly_for_trend, gw_yearly, lag = 0, sheet = sheet, id_tag = id_tag, veg = veg, out_dir = out_dir)
      } else {
        message('Not enough paired years for unlagged scatter (need >=3) for veg ', veg)
      }
    }, error = function(e) message('Unlagged scatter plotting skipped for ', sheet, ': ', e$message))

    # Combined z-scored trend plot (unlagged pairing)
    tryCatch({
      plot_trends_for_sheet(if(nrow(inf_yearly_for_trend)>0) inf_yearly_for_trend else NULL, gw_yearly, if(nrow(joined_unlagged)>0) joined_unlagged else NULL, sheet, id_tag, veg, out_dir, default_lag)
      message('Saved combined trend plot for: ', paste0(sheet, ' / ', id_tag, ' / ', veg))
    }, error = function(e) message('Trend plotting skipped for ', sheet, ': ', e$message))

    # dual-axis time series: vegetation fraction (primary) and groundwater (secondary) — UNLAGGED
    tryCatch({
      plot_dual_axis_time_series(if(nrow(inf_yearly_for_trend)>0) inf_yearly_for_trend else NULL, gw_yearly, sheet, id_tag, veg, out_dir)
    }, error = function(e) message('Dual-axis plot skipped for ', sheet, ': ', e$message))

    veg_results[[veg]] <- list(stats = stats_tbl, joined_unlagged = joined_unlagged) 

    # dual-axis time series: vegetation fraction (primary) and groundwater (secondary) — UNLAGGED
    tryCatch({
      plot_dual_axis_time_series(if(nrow(inf_yearly_for_trend)>0) inf_yearly_for_trend else NULL, gw_yearly, sheet, id_tag, veg, out_dir)
    }, error = function(e) message('Dual-axis plot skipped for ', sheet, ': ', e$message))

    veg_results[[veg]] <- list(stats = stats_tbl, joined_unlagged = joined_unlagged)
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
    message(sprintf("Sheet %s / veg %s -> id=%s | unlagged (lag=0): r=%.3f p=%.3g n=%d; best lag=%d: r=%.3f p=%.3g n=%d",
                    sheet, veg, fr$chosen_location, u_r, u_p, u_n,
                    best$lag, best$pearson_r, best$p_value, best$n))
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