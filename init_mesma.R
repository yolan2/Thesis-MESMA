setup_parallel_backend <- function(workers = NULL) {
  if (!requireNamespace("future", quietly = TRUE)) {
    warning("setup_parallel_backend: 'future' package not available; parallel execution disabled.")
    return(function() {})
  }

  # Capture current plan to allow restoration
  old_plan <- future::plan()

  # Determine the number of workers: respect MESMA_PARALLEL env var; otherwise use detectCores()-1
  if (is.null(workers)) {
    env <- Sys.getenv("MESMA_PARALLEL", "")
    if (nzchar(env)) {
      w <- suppressWarnings(as.integer(env))
      if (is.na(w) || w <= 0) w <- max(1L, parallel::detectCores() - 1L)
      workers <- w
    } else {
      workers <- max(1L, parallel::detectCores() - 1L)
    }
  }

  if (is.null(workers) || workers <= 1L) {
    cat("[PARALLEL] Parallel disabled (workers <= 1). Running sequentially.\n")
    return(function() { invisible(future::plan(old_plan)) })
  }

  future::plan(future::multisession, workers = workers)
  cat(sprintf("[PARALLEL] Set up 'future' multisession with %d workers\n", workers))

  cleanup <- function() {
    future::plan(old_plan)
    cat("[PARALLEL] Restored previous plan\n")
    invisible(NULL)
  }
  cleanup
}
