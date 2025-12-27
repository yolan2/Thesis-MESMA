#!/usr/bin/env Rscript
# Small wrapper to run the full MESMA fit or inference from the command line
# Usage examples:
#  Rscript run_fit_mesma.R --mode=fit --input=/path/to/training.csv --outdir=./phenology_results/veg_mixture_fit --workers=4
#  Rscript run_fit_mesma.R --mode=infer --inference_csv=/path/to/inference.csv --workers=2

args_raw <- commandArgs(trailingOnly = TRUE)

print_help <- function() {
  cat("run_fit_mesma.R - wrapper to run full MESMA fit or inference (sources fit_veg_mixture_mesma.R)\n")
  cat("\nOptions:\n")
  cat("  --mode=fit|infer        Run 'fit' (build library & save cache) or 'infer' (run inference only). Default: fit\n")
  cat("  --input=PATH            CSV path to training input (sets INPUT_CSV)\n")
  cat("  --inference_csv=PATH    CSV path for inference tasks (sets INFERENCE_CSV)\n")
  cat("  --outdir=PATH           Output directory (sets OUTPUT_DIR / OUT_DIR)\n")
  cat("  --workers=N             Number of parallel workers (sets MESMA_PARALLEL env var)\n")
  cat("  --no-parallel           Disable parallel processing (sets MESMA_PARALLEL=FALSE)\n")
  cat("  --help                  Print this help and exit\n")
  cat("\nExamples:\n  Rscript run_fit_mesma.R --mode=fit --input=data/train.csv --outdir=phenology_results/veg_mixture_fit --workers=4\n")
}

# Simple key=value parser
parse_args <- function(raw) {
  args <- list()
  for (a in raw) {
    if (a == "--help" || a == "-h") {
      args$help <- TRUE
      next
    }
    if (grepl("=", a)) {
      kv <- strsplit(a, "=", fixed = TRUE)[[1]]
      key <- sub('^--', '', kv[1])
      val <- paste(kv[-1], collapse = "=")
      args[[key]] <- val
    } else {
      # flags like --no-parallel
      key <- sub('^--', '', a)
      args[[key]] <- TRUE
    }
  }
  args
}

opts <- parse_args(args_raw)
if (!is.null(opts$help)) {
  print_help()
  quit(save = "no", status = 0)
}

# If script was launched via Rscript --file=..., set wd to script dir so relative paths work
suppressWarnings({
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- ca[grepl("--file=", ca)]
  if (length(file_arg) == 1) {
    spath <- sub("--file=", "", file_arg)
    if (nzchar(spath)) {
      script_dir <- dirname(normalizePath(spath))
      setwd(script_dir)
    }
  }
})

cat(sprintf("[RUNNER] Working dir: %s\n", getwd()))

if (!file.exists("fit_veg_mixture_mesma.R")) {
  stop("fit_veg_mixture_mesma.R not found in current working directory. Please run from project root or pass full paths.")
}

# Source main script
cat("[RUNNER] Sourcing fit_veg_mixture_mesma.R...\n")
source("fit_veg_mixture_mesma.R")

# Apply CLI overrides to key globals
if (!is.null(opts$input)) {
  INPUT_CSV <- opts$input
  cat(sprintf("[RUNNER] Overriding INPUT_CSV -> %s\n", INPUT_CSV))
}
if (!is.null(opts$inference_csv)) {
  INFERENCE_CSV <- opts$inference_csv
  cat(sprintf("[RUNNER] Overriding INFERENCE_CSV -> %s\n", INFERENCE_CSV))
}
if (!is.null(opts$outdir)) {
  OUTPUT_DIR <- opts$outdir
  OUT_DIR <- file.path(OUTPUT_DIR, "veg_mixture_fit")
  cat(sprintf("[RUNNER] Overriding OUTPUT_DIR -> %s\n", OUTPUT_DIR))
}
if (!is.null(opts$workers)) {
  Sys.setenv(MESMA_PARALLEL = opts$workers)
  cat(sprintf("[RUNNER] Set MESMA_PARALLEL=%s\n", opts$workers))
}
if (!is.null(opts$`no-parallel`) && isTRUE(opts$`no-parallel`)) {
  Sys.setenv(MESMA_PARALLEL = "FALSE")
  cat("[RUNNER] Parallel disabled (MESMA_PARALLEL=FALSE)\n")
}

mode <- if (!is.null(opts$mode)) tolower(opts$mode) else "fit"
cat(sprintf("[RUNNER] Mode: %s\n", mode))

run_fit <- function() {
  cat("[RUNNER] Running full fit (main_processing_block)...\n")
  if (!exists("main_processing_block") || !is.function(main_processing_block)) {
    stop("main_processing_block() not defined. The fit script may have changed.")
  }
  res <- tryCatch({
    main_processing_block()
    list(success = TRUE)
  }, error = function(e) {
    cat(sprintf("[ERROR] Fit failed: %s\n", e$message))
    list(success = FALSE, error = e$message)
  })
  if (!isTRUE(res$success)) quit(save = "no", status = 1)
  cat("[RUNNER] Fit completed successfully. Model cache should be saved by the script.\n")
}

run_infer <- function() {
  cat("[RUNNER] Running inference (process_inference if available)...\n")
  if (exists("process_inference") && is.function(process_inference)) {
    res <- tryCatch({ ri <- process_inference(); list(success = TRUE, result = ri) }, error = function(e) { cat(sprintf("[ERROR] process_inference failed: %s\n", e$message)); list(success = FALSE, error = e$message) })
    if (!isTRUE(res$success)) quit(save = "no", status = 1)
    cat("[RUNNER] Inference completed successfully.\n")
  } else if (exists("main_processing_block") && is.function(main_processing_block)) {
    cat("[RUNNER] process_inference not found; running main_processing_block() as fallback (may include fit+infer).\n")
    run_fit()
  } else {
    stop("No inference function found (process_inference or main_processing_block)")
  }
}

if (mode == "fit") {
  run_fit()
} else if (mode %in% c("infer", "inference")) {
  run_infer()
} else {
  stop("Unknown mode. Use --mode=fit or --mode=infer")
}

cat("[RUNNER] Done.\n")
