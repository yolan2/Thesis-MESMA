path <- "phenology_results/veg_mixture_fit/mesma_cache/compressed_templates.rds"
if (!file.exists(path)) {
  cat("File not found:", path, "\n")
} else {
  data <- readRDS(path)
  cat("Structure of compressed_templates_accessor:\n")
  if (length(data) > 0) {
    veg <- names(data)[1]
    cat("First vegetation:", veg, "\n")
    if (length(data[[veg]]) > 0) {
      vid <- names(data[[veg]])[1]
      cat("First variant:", vid, "\n")
      cat("Keys in variant:", paste(names(data[[veg]][[vid]]), collapse=", "), "\n")
    } else {
      cat("Vegetation", veg, "is empty.\n")
    }
  } else {
    cat("Accessor is empty.\n")
  }
}
