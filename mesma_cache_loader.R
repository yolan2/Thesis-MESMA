# Helper functions to load and verify a cached MESMA model and
# to rebuild the compressed templates accessor for inference-only runs.

# Load and verify the MESMA cache saved by fit_veg_mixture_mesma.R
rload_mesma_cache <- function(cache_dir = file.path(OUT_DIR, "mesma_cache")) {
  cat("\n=== LOADING MESMA MODEL CACHE ===\n")

  if (!dir.exists(cache_dir)) {
    stop(sprintf("Cache directory not found: %s", cache_dir))
  }

  # Load and verify manifest
  manifest_file <- file.path(cache_dir, "manifest.rds")
  if (!file.exists(manifest_file)) {
    stop("Cache manifest not found")
  }

  manifest <- readRDS(manifest_file)
  cat(sprintf("Cache version: %s\n", manifest$version))
  cat(sprintf("Cache created: %s\n", manifest$created))

  # Verify checksums
  if (!is.null(manifest$files) && length(manifest$files) > 0) {
    for (f in manifest$files) {
      fpath <- file.path(cache_dir, f)
      if (!file.exists(fpath)) {
        stop(sprintf("Missing cache file: %s", f))
      }
      current_checksum <- tryCatch(tools::md5sum(fpath), error = function(e) stop(sprintf("Failed to compute checksum for %s: %s", fpath, e$message)))
      recorded <- manifest$checksum[[f]]
      if (!is.na(current_checksum) && !is.null(recorded) && current_checksum != recorded) {
        warning(sprintf("Checksum mismatch for %s", f))
      }
    }
  }

  cache <- list()

  # Load library
  lib_path <- file.path(cache_dir, "mesma_library.rds")
  if (!file.exists(lib_path)) stop("mesma_library.rds not found in cache")
  lib_data <- readRDS(lib_path)
  cache$lib <- lib_data$lib
  cache$mesma_lib <- lib_data$mesma_lib
  cache$lib_factor_pca <- lib_data$lib_factor_pca
  cache$lib_factor_lda <- lib_data$lib_factor_lda
  cache$veg_counts <- lib_data$veg_counts
  cache$avail <- lib_data$avail
  cache$ALLOWED_VEG <- lib_data$ALLOWED_VEG
  cache$BAND_SCALE <- lib_data$BAND_SCALE
  cache$COMPRESSED_STAGE1_LIB <- lib_data$COMPRESSED_STAGE1_LIB

  # Load projections
  proj_path <- file.path(cache_dir, "projection_matrices.rds")
  if (!file.exists(proj_path)) stop("projection_matrices.rds not found in cache")
  proj_data <- readRDS(proj_path)
  cache$GLOBAL_PCA <- proj_data$GLOBAL_PCA
  cache$GLOBAL_LDA <- proj_data$GLOBAL_LDA
  cache$pca_rank <- proj_data$pca_rank

  # Load compressed templates (optional)
  templates_path <- file.path(cache_dir, "compressed_templates.rds")
  if (file.exists(templates_path)) {
    template_data <- readRDS(templates_path)
    cache$compressed_templates <- template_data
  } else {
    cache$compressed_templates <- list()
  }

  # Load configuration
  cfg_path <- file.path(cache_dir, "config_params.rds")
  if (file.exists(cfg_path)) cache$config <- readRDS(cfg_path) else cache$config <- list()

  # Load metadata
  meta_path <- file.path(cache_dir, "training_metadata.rds")
  if (file.exists(meta_path)) cache$metadata <- readRDS(meta_path) else cache$metadata <- list()

  cat("Cache loaded successfully\n")
  if (!is.null(cache$metadata$n_training_samples)) {
    cat(sprintf("Trained on %d samples from %d locations\n",
                cache$metadata$n_training_samples,
                cache$metadata$n_locations_trained))
  }
  if (!is.null(cache$metadata$training_years)) {
    cat(sprintf("Training years: %s\n", paste(cache$metadata$training_years, collapse = ", ")))
  }
  if (!is.null(cache$metadata$vegetation_types)) {
    cat(sprintf("Vegetation types: %s\n", paste(cache$metadata$vegetation_types, collapse = ", ")))
  }

  return(cache)
}


# Function to rebuild the compressed templates accessor from cached data
rebuild_template_accessor <- function(template_data) {
  template_db <- new.env(hash = TRUE, parent = emptyenv())

  if (!is.null(template_data) && length(template_data) > 0) {
    for (veg in names(template_data)) {
      for (key in names(template_data[[veg]])) {
        assign(key, template_data[[veg]][[key]], envir = template_db)
      }
    }
  }

  accessor <- function(veg, variant_id, projection, grid_type = "medium") {
    key <- paste(veg, variant_id, projection, grid_type, sep = "|")
    if (exists(key, envir = template_db, inherits = FALSE)) {
      get(key, envir = template_db, inherits = FALSE)
    } else {
      NULL
    }
  }

  return(accessor)
}


# Convenience: load cache and register accessor into globalenv()
# Usage: register_mesma_cache(cache_dir)
register_mesma_cache <- function(cache_dir = file.path(OUT_DIR, "mesma_cache")) {
  c <- rload_mesma_cache(cache_dir)
  if (!is.null(c$compressed_templates) && length(c$compressed_templates) > 0) {
    accessor <- rebuild_template_accessor(c$compressed_templates)
    assign('.COMPRESSED_TEMPLATES_ACCESSOR', accessor, envir = globalenv())
    cat("Registered .COMPRESSED_TEMPLATES_ACCESSOR in global environment\n")
  } else {
    cat("No compressed templates found in cache; .COMPRESSED_TEMPLATES_ACCESSOR not registered\n")
  }
  # Also return the loaded cache
  invisible(c)
}
