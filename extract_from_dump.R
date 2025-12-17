#!/usr/bin/env Rscript
# Extract data from error dump using debugger

cat("=== LOADING ERROR DUMP ===\n")
load("mesma_error_dump.rda")

cat("\n=== AVAILABLE IN DUMP ===\n")
cat("Class:", class(mesma_error_dump), "\n")
cat("Length:", length(mesma_error_dump), "\n")

# The dump contains environment frames - we need to extract variables
# Let's look at the main processing block frame
cat("\n=== SEARCHING FOR INFERENCE DATA ===\n")

# Try to find the frame with df_tasks_inference
for (i in seq_along(mesma_error_dump)) {
  env <- mesma_error_dump[[i]]
  if (exists("df_tasks_inference", envir=env, inherits=FALSE)) {
    cat(sprintf("Found df_tasks_inference in frame %d\n", i))
    df_tasks_inference <- get("df_tasks_inference", envir=env)

    cat(sprintf("Rows: %d\n", nrow(df_tasks_inference)))
    cat(sprintf("Columns: %d\n", ncol(df_tasks_inference)))
    cat("\nFirst few column names:\n")
    print(head(colnames(df_tasks_inference), 30))

    cat("\n=== CHECKING INDICES COMPLETENESS ===\n")
    required_indices <- c("NBR", "NDTI", "NDVI", "EVI", "SAVI", "MSAVI", "kNDVI",
                          "NDMI", "NDWI", "TCB", "TCG", "TCW", "PPI")

    for (idx in required_indices) {
      if (idx %in% colnames(df_tasks_inference)) {
        n_finite <- sum(is.finite(df_tasks_inference[[idx]]))
        pct <- 100 * n_finite / nrow(df_tasks_inference)
        marker <- if (pct < 50) " ***LOW***" else ""
        cat(sprintf("  %-8s: %6.1f%% valid (%7d/%7d)%s\n",
                    idx, pct, n_finite, nrow(df_tasks_inference), marker))
      } else {
        cat(sprintf("  %-8s: MISSING ***\n", idx))
      }
    }

    break
  }
}

# Try to find training data
cat("\n=== SEARCHING FOR TRAINING DATA ===\n")
for (i in seq_along(mesma_error_dump)) {
  env <- mesma_error_dump[[i]]
  if (exists("df_train", envir=env, inherits=FALSE)) {
    cat(sprintf("Found df_train in frame %d\n", i))
    df_train <- get("df_train", envir=env)
    cat(sprintf("Rows: %d\n", nrow(df_train)))

    cat("\nTraining indices completeness:\n")
    for (idx in required_indices) {
      if (idx %in% colnames(df_train)) {
        n_finite <- sum(is.finite(df_train[[idx]]))
        pct <- 100 * n_finite / nrow(df_train)
        cat(sprintf("  %-8s: %6.1f%% valid\n", idx, pct))
      }
    }
    break
  }
}

# Try to find normalization params
cat("\n=== SEARCHING FOR NORMALIZATION PARAMS ===\n")
for (i in seq_along(mesma_error_dump)) {
  env <- mesma_error_dump[[i]]
  if (exists("TRAINING_NORM_PARAMS", envir=env, inherits=FALSE)) {
    cat(sprintf("Found TRAINING_NORM_PARAMS in frame %d\n", i))
    TRAINING_NORM_PARAMS <- get("TRAINING_NORM_PARAMS", envir=env)
    cat("Indices in normalization params:\n")
    print(names(TRAINING_NORM_PARAMS))
    break
  }
}

# Check for results
cat("\n=== SEARCHING FOR RESULTS ===\n")
for (i in seq_along(mesma_error_dump)) {
  env <- mesma_error_dump[[i]]
  if (exists("results", envir=env, inherits=FALSE)) {
    cat(sprintf("Found results in frame %d\n", i))
    results <- get("results", envir=env)
    cat("Results class:", class(results), "\n")
    cat("Results length:", length(results), "\n")
    if (length(results) > 0) {
      cat("First result:\n")
      print(str(results[[1]]))
    }
    break
  }
}
