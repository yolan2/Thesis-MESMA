# Test: ensure inference sampling works and low observation counts are allowed in inference
source(file.path('..', 'fit_veg_mixture_mesma.R'))

# Create test df_inf with 3000 locations, some with few DOYs
nloc <- 3000
set.seed(123)
loc_ids <- paste0('L_', sprintf('%07d', seq_len(nloc)))
df_inf_test <- data.frame(
  location_id = rep(loc_ids, each = 1),
  date = Sys.Date() - sample(1:1000, nloc, replace = TRUE),
  PPI = runif(nloc),
  DVI = runif(nloc),
  DVI_max = runif(nloc) + 0.1,
  Veg = sample(c('phragmites','populus','tamarix','barren'), nloc, replace = TRUE),
  stringsAsFactors = FALSE
)

# Simulate the inference selection logic
inference_location_ids <- unique(df_inf_test$location_id)
cat(sprintf('Initial unique inference locations: %d\n', length(inference_location_ids)))

if (length(inference_location_ids) > INFERENCE_MAX_LOCATIONS) {
  set.seed(42)
  inference_location_ids <- sample(inference_location_ids, INFERENCE_MAX_LOCATIONS)
  cat(sprintf('After sampling, inference locations limited to %d\n', length(inference_location_ids)))
}

# Create minimal df_tasks and test_loc_years similarly
print(head(inference_location_ids, 5))
