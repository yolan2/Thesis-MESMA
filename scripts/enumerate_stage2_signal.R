MESMA_NO_AUTO_RUN <- TRUE
source('fit_veg_mixture_mesma.R')

if (!exists('df_tasks') || is.null(df_tasks) || nrow(df_tasks) == 0) {
  stop('df_tasks not found in workspace; cannot enumerate')
}

# Unique task rows (by location, pheno_year)
uniq_tasks <- unique(df_tasks[c('location_id', 'pheno_year')])
res_list <- vector('list', nrow(uniq_tasks))

for (i in seq_len(nrow(uniq_tasks))) {
  loc <- as.character(uniq_tasks$location_id[i])
  yr <- as.integer(uniq_tasks$pheno_year[i])
  task_df <- df_tasks[df_tasks$location_id == loc & df_tasks$pheno_year == yr, , drop = FALSE]
  diag <- tryCatch({ fit_one_task(task_df, compute_stage2_diagnostics = TRUE) }, error = function(e) { list(error = e$message, location_id = loc, pheno_year = yr) })
  res_list[[i]] <- diag
}

# Convert to data.frame (robust scalar extraction)
get_scalar <- function(x, name) {
  if (is.null(x) || is.null(x[[name]])) return(NA)
  v <- x[[name]]
  if (length(v) == 0) return(NA)
  return(v[[1]])
}

rows <- list()
for (x in res_list) {
  if (is.null(x)) next
  # Some error handlers return lists with only error/msg fields; handle robustly
  loc <- get_scalar(x, 'location_id')
  yr <- get_scalar(x, 'pheno_year')
  n_raw_bins <- get_scalar(x, 'n_raw_bins')
  n_valid_mask <- get_scalar(x, 'n_valid_mask')
  y_s2_norm_val <- get_scalar(x, 'y_s2_norm_val')
  y_s2_masked_mean <- get_scalar(x, 'y_s2_masked_mean')
  weights_fallback <- get_scalar(x, 'weights_s2_fallback_used')

  rows[[length(rows)+1]] <- data.frame(
    location_id = as.character(ifelse(is.na(loc), NA_character_, loc)),
    pheno_year = as.integer(ifelse(is.na(yr), NA_integer_, yr)),
    n_raw_bins = as.integer(ifelse(is.na(n_raw_bins), NA_integer_, n_raw_bins)),
    n_valid_mask = as.integer(ifelse(is.na(n_valid_mask), NA_integer_, n_valid_mask)),
    y_s2_norm_val = as.numeric(ifelse(is.na(y_s2_norm_val), NA_real_, y_s2_norm_val)),
    y_s2_masked_mean = as.numeric(ifelse(is.na(y_s2_masked_mean), NA_real_, y_s2_masked_mean)),
    weights_s2_fallback_used = as.logical(ifelse(is.na(weights_fallback), NA, weights_fallback)),
    stringsAsFactors = FALSE
  )
}

if (length(rows) == 0) stop('No diagnostics captured')
out_df <- do.call(rbind, rows)

out_path <- 'stage2_signal_summary.csv'
write.csv(out_df, out_path, row.names = FALSE)
cat(sprintf('Wrote summary to %s\n', out_path))
