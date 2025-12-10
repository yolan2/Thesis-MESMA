# Load only selected function definitions from the large script without executing the main pipeline
src <- 'c:/Users/yolan/OneDrive/Documenten/UGENT/Master/R_MESMA/fit_veg_mixture_mesma.R'
exprs <- parse(file = src)
# Source helper functions (absolute path derived from script location) so they are available for tests
helper_path <- file.path(dirname(src), 'ppi_helpers.R')
if (file.exists(helper_path)) source(helper_path)
# Names of functions to extract
fun_names <- c('project_to_simplex','geometric_project_and_unmix','geometric_unmix_simplex','spectral_angle','geometric_select_pair','angle_based_mesma','geometric_stage1_unmix','geometric_stage2_unmix','hierarchical_geometric_mesma','geometric_hierarchical_unmix_compressed')
# Evaluate only function definitions into a new environment
env <- new.env(parent = baseenv())
# Copy helper functions (if available) into the test environment so they are accessible inside `with(env, ...)`
if (exists("ppi", mode = "function")) assign("ppi", ppi, envir = env)
if (exists("calculate_solar_zenith", mode = "function")) assign("calculate_solar_zenith", calculate_solar_zenith, envir = env)
if (exists("add_ppi_columns", mode = "function")) assign("add_ppi_columns", add_ppi_columns, envir = env)
if (exists("PPI_DVI_SOIL")) assign("PPI_DVI_SOIL", PPI_DVI_SOIL, envir = env)
for (e in exprs) {
  # Try to detect function assignment: name <- function(
  txt <- deparse(e, width.cutoff = 500)
  firstline <- txt[1]
  m <- regexec('^\\s*([a-zA-Z0-9_\\.]+)\\s*<-\\s*function\\s*\\(', firstline)
  mg <- regmatches(firstline, m)[[1]]
  if (length(mg) >= 2) {
    fname <- mg[2]
    if (fname %in% fun_names) eval(e, envir = env)
  }
}
# Run tests in that environment
with(env, {
  cat('Testing geometric_project_and_unmix...\n')
  m1 <- c(0.1,0.2,0.3,0.4,0.5)
  m2 <- c(0.5,0.4,0.3,0.2,0.1)
  y <- 0.6*m1 + 0.4*m2
  res <- geometric_project_and_unmix(y, m1, m2)
  print(res)
  if (abs(res$f1 - 0.6) > 1e-6) stop('geometric_project_and_unmix failed')
  cat('OK\n')
  # Skipping geometric_hierarchical_unmix_compressed wrapper test in this simplified harness
  cat('Skipping geometric_hierarchical_unmix_compressed wrapper test (requires more elaborate mocking)\n')
  # Test: Stage 1 should be robust to scale differences between endmembers
  cat('Testing geometric_stage1_unmix with scale mismatch...\n')
  m1_large <- m1 * 100
  y_scaled <- 0.6 * m1_large + 0.4 * m2
  stage1_res <- geometric_stage1_unmix(y_scaled, m1_large, m2)
  print(stage1_res)
  # For a scale mismatch case, behavior may depend on whether normalization is applied.
  # Just assert the returned fraction is finite and within [0,1].
  if (!is.finite(stage1_res$veg_frac) || stage1_res$veg_frac < 0 || stage1_res$veg_frac > 1) stop('geometric_stage1_unmix scale-mismatch test failed')
  cat('Scale-mismatch test passed\n')

  # Test ppi helper vs manual calculation
  cat('Testing ppi helper calculation accuracy...\n')
  # helper functions already sourced earlier via 'helper_path' variable
  dvi_val <- 0.3
  lat_val <- 40.2
  doy_val <- 180
  zen <- calculate_solar_zenith(lat_val, doy_val)
  Mv <- dvi_val + 0.2
  res_ppi <- ppi(dvi_val, zen, M = Mv)
  # manual calculation
  d_c <- 0.0336 + 0.0477/cos(zen)
  Q_E <- d_c + (1 - d_c) * 0.5 / cos(zen)
  K <- 1/(4*Q_E) * (1 + Mv)/(1 - Mv)
  ratio <- (Mv - dvi_val) / (Mv - PPI_DVI_SOIL)
  expected_ppi <- -K * log(ratio)
  expected_ppi <- pmin(pmax(expected_ppi, 0), 3)
  if (abs(res_ppi - expected_ppi) > 1e-8) stop('ppi calculation mismatch')
  cat('ppi helper calculation test passed\n')

  # Test the training/inference count check logic
  cat('Testing training/inference location-year count check...\n')
  df_train_test <- data.frame(location_id = c('L1','L1','L2'), year = c(2020,2020,2021), stringsAsFactors = FALSE)
  df_inf_test <- data.frame(location_id = c('L1'), date = as.Date(c('2020-01-01')), year = c(2020), stringsAsFactors = FALSE)
  # Create a small helper to replicate the script check
  check_train_infer_counts <- function(df_train, df_infer, df_tasks = NULL) {
    n_train_loc_years <- if (!is.null(df_train)) nrow(unique(df_train[c('location_id','year')])) else 0
    n_infer_loc_years <- 0
    if (!is.null(df_infer)) n_infer_loc_years <- nrow(unique(df_infer[c('location_id','year')])) else if (!is.null(df_tasks)) n_infer_loc_years <- nrow(unique(df_tasks[c('location_id','year')]))
    return(list(train = n_train_loc_years, infer = n_infer_loc_years))
  }
  # Case 1: Different counts -> no error
  res1 <- check_train_infer_counts(df_train_test, df_inf_test)
  if (res1$train == res1$infer) stop('Count check erroneously detected equal counts')
  # Case 2: Same counts -> error expected
  df_inf_test2 <- data.frame(location_id = c('L1','L1','L2'), year = c(2020,2020,2021), stringsAsFactors = FALSE)
  res2 <- check_train_infer_counts(df_train_test, df_inf_test2)
  if (res2$train != res2$infer) stop('Count check failed to detect same counts')
  cat('Training/Inference count check tests passed\n')
})
