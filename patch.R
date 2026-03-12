
apply_stored_normalization <- function(df, norm_params, cols = unique(c(OPTIMAL_INDICES, RAW_BANDS)), lat_default = 40.2) {
  cat('Applying stored normalization parameters to data...\n')

  if (!'zenith.angle' %in% names(df)) df\.angle <- NA_real_
  if (is.null(norm_params\) || length(norm_params\) == 0) {
    stop('[INDEX_SCALES] Missing stored normalization parameters (INDEX_SCALES)')
  }

  if (!'PPI' %in% names(df) || all(!is.finite(df\))) {
    if (all(c('nir', 'red') %in% names(df)) && !'DVI' %in% names(df)) {
      df\ <- as.numeric(df\) - as.numeric(df\)
    }
    if (exists('add_ppi_columns')) {
      dvi_soil_vec <- compute_dvi_soil_per_location(df)
      df <- add_ppi_columns(df, dvi_soil = dvi_soil_vec)
      cat('[PPI] Added PPI to data before applying stored normalization (per-location dvi_soil + per-location M).\n')
    }
  }
  if (!'PPI' %in% names(df) || all(!is.finite(df\))) {
    stop('[PPI] Missing or all PPI values non-finite after attempted auto-add; refusing to continue')
  }
  if (!'PPI_raw' %in% names(df)) df\ <- df\

  cat(sprintf('[apply_stored_normalization] Applying INDEX_SCALES to %d indices\n', length(norm_params\)))
  n_scaled <- 0L
  for (col in names(norm_params\)) {
    if (col %in% names(df)) {
      params <- norm_params\[[col]]
      if (!is.list(params) || !all(c('mean', 'sd') %in% names(params))) {
        stop(sprintf('[INDEX_SCALES] Invalid params for index \\'%s\\' (expected list(mean, sd))', col))
      }
      mu <- params\
      sigma <- params\
      if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0) {
        stop(sprintf('[INDEX_SCALES] Non-finite or non-positive mean/sd for index \\'%s\\'', col))
      }
      df[[col]] <- (df[[col]] - mu) / sigma
      n_scaled <- n_scaled + 1L
    }
  }
  cat(sprintf('[apply_stored_normalization] Z-scored %d feature columns\n', n_scaled))

  df
}

