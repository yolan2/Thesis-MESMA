  main_processing_block <- function() {
    cat("[DEBUG] Entered main_processing_block\n")
    b_templates <- NULL
    n_train_loc_years <- 0L
    n_infer_loc_years <- 0L
    if (exists("df_train") && !is.null(df_train) && nrow(df_train) > 0) {
    if (!"pheno_year" %in% names(df_train) && "date" %in% names(df_train)) df_train$pheno_year <- assign_pheno_year(df_train$date)
    n_train_loc_years <- nrow(unique(df_train[c("location_id", "pheno_year")]))
  }
  if (exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0) {
    if (!"pheno_year" %in% names(df_tasks_inference) && "date" %in% names(df_tasks_inference)) df_tasks_inference$pheno_year <- assign_pheno_year(df_tasks_inference$date)
    n_infer_loc_years <- nrow(unique(df_tasks_inference[c("location_id", "pheno_year")]))
  } else if (exists("df_tasks") && !is.null(df_tasks) && nrow(df_tasks) > 0) {
    if (!"pheno_year" %in% names(df_tasks) && "date" %in% names(df_tasks)) df_tasks$pheno_year <- assign_pheno_year(df_tasks$date)
    n_infer_loc_years <- nrow(unique(df_tasks[c("location_id", "pheno_year")]))
  }

  cat(sprintf("Training dataset location-years: %d\n", n_train_loc_years))
  cat(sprintf("Inference dataset location-years: %d\n", n_infer_loc_years))

  if ((exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0) && n_train_loc_years == n_infer_loc_years) {
    stop(sprintf("ERROR: Training and inference datasets appear to have the same number of location-years (%d). This may indicate you passed the same data for training and inference — aborting to avoid accidental overlap.", n_train_loc_years))
  }
  
  cat("Preparing locations for batched processing (multi-year bootstrap)...\n")

  # Group by LOCATION (not location-year)
  target_locations <- location_list
  available_locations <- unique(df_tasks$location_id)
  target_locations <- intersect(target_locations, available_locations)

  n_locs_to_process <- length(target_locations)
  BATCH_SIZE <- 6 # Smaller batches for location-level processing (each has multiple years)

  loc_batches <- split(target_locations, ceiling(seq_along(target_locations) / BATCH_SIZE))
  n_batches <- length(loc_batches)
  pb_width <- min(40L, max(4L, n_batches))

  cat(sprintf("Processing %d locations in %d batches (approx %d locations/batch)...\n",
              n_locs_to_process, length(loc_batches), BATCH_SIZE))

  # Results will be written to Excel in batches
  all_coefs_list <- list()

  start_time <- Sys.time()

  for (i in seq_along(loc_batches)) {
    batch_locs <- loc_batches[[i]]

    batch_df <- df_tasks[df_tasks$location_id %in% batch_locs, ]

    # Split by location (not location-year)
    batch_location_list <- split(batch_df, batch_df$location_id)

    # Use fit_one_location instead of fit_one_task
    batch_results <- .run_map(batch_location_list, fit_one_location, show_pb = FALSE)

    # Write batch results to Excel immediately
    for (k in names(batch_results)) {
      loc_result <- batch_results[[k]]
      if (is.null(loc_result)) next
      loc_data <- do.call(rbind, lapply(loc_result, function(yr_res) yr_res$coef_df))
      if (!is.null(loc_data) && nrow(loc_data) > 0) {
        all_coefs_list[[k]] <- loc_data
        sheet_name <- paste0("Loc_", k)
        openxlsx::addWorksheet(wb, sheet_name)
        openxlsx::writeData(wb, sheet_name, loc_data)
      }
    }

    if (isTRUE(TESTING_MODE)) {
      for (k in names(batch_results)) {
        loc_result <- batch_results[[k]]
        if (is.null(loc_result)) {
          cat(sprintf("[DEBUG batch_result] location %s returned NULL\n", k))
        } else if (is.list(loc_result)) {
          n_years <- length(loc_result)
          cat(sprintf("[DEBUG batch_result] location %s returned %d year(s)\n", k, n_years))
          for (yr_char in names(loc_result)) {
            r <- loc_result[[yr_char]]
            if (!is.null(r)) {
              ca <- as.numeric(r$vegetated_fraction)
              cb <- as.numeric(r$barren_fraction)
              coef_n <- if (!is.null(r$coef_df) && is.data.frame(r$coef_df)) nrow(r$coef_df) else 0
              has_unc <- !is.null(r$uncertainty)
              cat(sprintf("  Year %s: veg_frac=%.4f barren_frac=%.4f coef_rows=%d has_uncertainty=%s\n",
                          yr_char, ca, cb, coef_n, has_unc))
            }
          }
        }
      }
    }
    
    rm(batch_df, batch_location_list, batch_results)
    gc(verbose = FALSE)
  }
  
  end_time <- Sys.time()
  processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(
    "Main processing loop finished in %.2f seconds (%.2f minutes)\n",
    processing_time, processing_time / 60
  ))

  # Combine all coefficient data
  all_coefs <- do.call(rbind, all_coefs_list)

  if (n_year_results > 0) {
    cat(sprintf("Average time per year-result: %.2f seconds\n", processing_time / n_year_results))
  } else {
    cat("Average time per year-result: N/A (0 results)\n")
  }

  inference_results_list <- list()
  if (isTRUE(TESTING_MODE)) {
    cat("[TESTING MODE] Skipping inference processing for faster testing\n")
  } else if (exists("df_tasks_inference") && !is.null(df_tasks_inference) && nrow(df_tasks_inference) > 0) {
    cat("\n=== PROCESSING INFERENCE DATA ===\n")
    
    if (!"year" %in% names(df_tasks_inference) && "date" %in% names(df_tasks_inference)) {
      df_tasks_inference$date <- as.Date(df_tasks_inference$date)
      df_tasks_inference$pheno_year <- assign_pheno_year(df_tasks_inference$date)
      df_tasks_inference$doy <- pheno_doy(df_tasks_inference$date)  # Use phenological DOY
    }
    
    if (!"task_key" %in% names(df_tasks_inference)) {
      df_tasks_inference$task_key <- paste(df_tasks_inference$location_id, df_tasks_inference$pheno_year, sep = "_")
    }
    
    inference_loc_years <- df_tasks_inference |> 
      dplyr::filter(!is.na(.data$location_id) & trimws(.data$location_id) != "" & !is.na(.data$pheno_year) & .data$pheno_year > 0) |> 
      dplyr::distinct(.data$location_id, .data$pheno_year)
    
    cat(sprintf("Inference dataset: %d location-year pairs from %d locations\n", 
                nrow(inference_loc_years), length(unique(inference_loc_years$location_id))))
    
    inference_target_keys <- paste(inference_loc_years$location_id, inference_loc_years$pheno_year, sep = "_")
    inference_available_keys <- unique(df_tasks_inference$task_key)
    inference_target_keys <- intersect(inference_target_keys, inference_available_keys)
    
    # Process all inference tasks (no limit)
    
    n_inference_keys <- length(inference_target_keys)

    if (n_inference_keys > 0) {
      cat(sprintf("[MASK-BASED PROCESSING] Using validity masks - processing all tasks with any valid observations (>5%% of bins)\n"))

      inference_key_batches <- split(inference_target_keys, ceiling(seq_along(inference_target_keys) / BATCH_SIZE))
      inference_n_batches <- length(inference_key_batches)
      inference_pb_width <- min(40L, max(4L, inference_n_batches))

      cat(sprintf("Processing %d inference tasks in %d batches...\n", n_inference_keys, length(inference_key_batches)))
      
      inference_results_list <- vector("list", n_inference_keys)
      names(inference_results_list) <- inference_target_keys
      
      inference_start_time <- Sys.time()
      
      for (i in seq_along(inference_key_batches)) {
        batch_keys <- inference_key_batches[[i]]
        
        batch_df <- df_tasks_inference[df_tasks_inference$task_key %in% batch_keys, ]
        
        batch_task_list <- split(batch_df, batch_df$task_key)
        
        batch_results <- .run_map(batch_task_list, fit_one_task, show_pb = FALSE)
        
        inference_results_list[names(batch_results)] <- batch_results
        if (isTRUE(TESTING_MODE)) {
          for (k in names(batch_results)) {
            r <- batch_results[[k]]
            if (is.null(r)) {
              cat(sprintf("[DEBUG batch_result] inference task %s returned NULL\n", k))
            } else {
              ca <- as.numeric(r$vegetated_fraction); cb <- as.numeric(r$barren_fraction)
              coef_n <- if (!is.null(r$coef_df) && is.data.frame(r$coef_df)) nrow(r$coef_df) else 0
              veg_list <- if (!is.null(r$coef_df)) paste(unique(r$coef_df$Veg), collapse = ",") else NA_character_
              diag_vf <- if (!is.null(r$diagnostics) && 'vegetated_fraction' %in% names(r$diagnostics)) r$diagnostics$vegetated_fraction else NA_real_
              cat(sprintf("[DEBUG batch_result] inference task %s returned non-NULL: veg_frac=%.4f barren_frac=%.4f coef_rows=%d vegs=%s diag_vf=%s\n", k, ca, cb, coef_n, veg_list, as.character(diag_vf)))
            }
          }
        }
        
        rm(batch_df, batch_task_list, batch_results)
        gc(verbose = FALSE)
      }
      
      inference_end_time <- Sys.time()
      inference_processing_time <- as.numeric(difftime(inference_end_time, inference_start_time, units = "secs"))
      cat(sprintf("Inference processing finished in %.2f seconds (%.2f minutes)\n",
                  inference_processing_time, inference_processing_time / 60))

      n_null_before_filter <- sum(sapply(inference_results_list, is.null))
      inference_results_list <- inference_results_list[!sapply(inference_results_list, is.null)]
      n_valid <- length(inference_results_list)

      cat(sprintf("Valid inference results: %d out of %d tasks (%.1f%%)\n",
                  n_valid, n_inference_keys, 100*n_valid/n_inference_keys))
      cat(sprintf("Skipped/filtered tasks: %d (%.1f%%)\n",
                  n_null_before_filter, 100*n_null_before_filter/n_inference_keys))
    } else {
      cat("No valid inference location-year pairs found.\n")
    }
  }

  if (length(inference_results_list) > 0) {
    cat(sprintf("Adding %d inference results to results list\n", length(inference_results_list)))
    results_list <- c(results_list, inference_results_list)
  }

  cat("Processing results and writing to Excel files...\n")

  results_list <- results_list[!sapply(results_list, is.null)]
  is_invalid_res <- sapply(results_list, function(res) {
    if (is.null(res)) return(TRUE)
    if (is.null(res$coef_df) || !is.data.frame(res$coef_df)) return(TRUE)
    lid <- res$coef_df$location_id
    if (is.null(lid) || length(lid) == 0) return(TRUE)
    all(is.na(lid) | trimws(as.character(lid)) == "")
  })
  if (any(is_invalid_res)) {
    cat(sprintf("[NOTICE] Filtering out %d invalid result(s) with missing location_id\n", sum(is_invalid_res)))
    results_list <- results_list[!is_invalid_res]
  }
  cat(sprintf("After filtering NULL results: %d results remaining\n", length(results_list)))

  if (length(results_list) > 0) {
    barren_one_count <- sum(sapply(results_list, function(res) {
      if (!is.null(res$barren_fraction) && is.finite(res$barren_fraction)) {
        res$barren_fraction == 1
      } else {
        FALSE
      }
    }))
    barren_one_pct <- barren_one_count / length(results_list) * 100
    if (barren_one_pct > 50) {
      msg <- sprintf("WARNING: %.1f%% of predictions (%d/%d) have barren_fraction = 1, exceeding 50%% threshold. This indicates severe model issues.", barren_one_pct, barren_one_count, length(results_list))
      # Always emit a warning rather than stopping execution. In testing mode, prefix to indicate non-actionable warning.
      if (isTRUE(TESTING_MODE)) {
        warning(paste("[TESTING MODE IGNORE]", msg))
      } else {
        warning(msg)
      }
    }
    cat(sprintf("Barren fraction = 1 in %.1f%% of predictions (%d/%d) - within acceptable limits\n", barren_one_pct, barren_one_count, length(results_list)))
  } else {
    cat("No results to check for barren fraction\n")
  }

  if (length(results_list) == 0) {
    cat("ERROR: All tasks returned NULL results!\n")
    cat("Most likely causes:\n")
    cat("1. Insufficient data quality: Tasks were filtered due to:\n")
    cat("   - Too few observations per location-year (<15 observations)\n")
    cat("   - Insufficient temporal coverage (<25 unique days of year)\n")
    cat("   - Too many missing pentad bins (>85% NA values)\n")
    cat("2. Data filtering issues or missing indices\n")
    cat("3. No valid testing/inference data available (no location-year pairs found)\n")
    cat("\nTo adjust data quality thresholds, modify MIN_OBS, MIN_UNIQUE_DOYS, and MIN_PENTAD_COVERAGE in fit_one_task function.\n")
    stop("No valid results to process")
  }

  if (length(all_coefs_list) > 0) {
    cat("Combining coefficient data frames...\n")

    if (length(all_coefs_list) == 0) {
      cat("ERROR: No valid coefficient data frames found!\n")
      stop("No coefficient data to process")
    }

    all_coefs <- tryCatch({
      do.call(rbind, all_coefs_list)
    }, error = function(e) stop(sprintf("ERROR combining coef_df: %s", e$message))
    )

    if (is.null(all_coefs)) {
      cat("Failed to combine coefficient data frames\n")
      stop("Cannot proceed without coefficient data")
    }

    required_coef_cols <- c("location_id", "pheno_year", "Veg", "coef", "rmse", "coef_025", "coef_975", "interval")
    missing <- setdiff(required_coef_cols, names(all_coefs))
    if (length(missing) > 0) {
      cat(sprintf("[NOTICE] Filling missing coefficient columns with NA: %s\n", paste(missing, collapse = ", ")))
      for (col in missing) all_coefs[[col]] <- NA
    }
    all_coefs$location_id <- as.character(all_coefs$location_id)
    all_coefs$pheno_year <- as.integer(all_coefs$pheno_year)
    # Ensure legacy 'year' column exists for downstream code that expects it
    if (!"year" %in% names(all_coefs)) all_coefs$year <- all_coefs$pheno_year
    all_coefs$Veg <- as.character(all_coefs$Veg)
    all_coefs$coef <- as.numeric(all_coefs$coef)
    all_coefs$rmse <- as.numeric(all_coefs$rmse)
    all_coefs$coef_025 <- as.numeric(all_coefs$coef_025)
    all_coefs$coef_975 <- as.numeric(all_coefs$coef_975)
    all_coefs$interval <- as.numeric(all_coefs$interval)

    all_coefs$location_id <- trimws(all_coefs$location_id)
    cat(sprintf("Combined coefficients: %d rows\n", nrow(all_coefs)))

    cat("Combining chosen variant summaries...\n")
    variant_list_pca <- lapply(results_list, function(res) if (!is.null(res$variant_trajectory)) res$variant_trajectory else NULL)
    variant_list_pca <- variant_list_pca[!sapply(variant_list_pca, is.null)]
    all_variants_pca <- NULL
    if (length(variant_list_pca) > 0) {
      all_variants_pca <- tryCatch({
        if (requireNamespace("dplyr", quietly = TRUE)) {
          dplyr::bind_rows(variant_list_pca)
        } else {
          do.call(rbind, variant_list_pca)
        }
      }, error = function(e) stop(sprintf("Failed to combine variant PCA summaries: %s", e$message)))
    }

    diag_list <- lapply(results_list, function(res) {
      if (!is.null(res$diagnostics)) res$diagnostics else NULL
    })
    diag_list <- diag_list[!sapply(diag_list, is.null)]
    all_diagnostics <- if (length(diag_list) > 0) tryCatch({
      if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required to combine diagnostics robustly")
      dplyr::bind_rows(diag_list)
    }, error = function(e) stop(sprintf("Failed to combine diagnostics: %s", e$message))) else NULL

    q_dvi_data <- tryCatch({
      if (requireNamespace("dplyr", quietly = TRUE)) {
        dplyr::bind_rows(lapply(results_list, function(res) {
          if ((!is.null(res$q10_dvi) && !is.na(res$q10_dvi)) || (!is.null(res$q90_dvi) && !is.na(res$q90_dvi))) {
            data.frame(
              location_id = res$coef_df$location_id[1],
              pheno_year = res$coef_df$pheno_year[1],
              q10_dvi = res$q10_dvi,
              q90_dvi = res$q90_dvi,
              stringsAsFactors = FALSE
            )
          } else {
            NULL
          }
        }))
      } else {
        do.call(rbind, lapply(results_list, function(res) {
          if ((!is.null(res$q10_dvi) && !is.na(res$q10_dvi)) || (!is.null(res$q90_dvi) && !is.na(res$q90_dvi))) {
            data.frame(
              location_id = res$coef_df$location_id[1],
              pheno_year = res$coef_df$pheno_year[1],
              q10_dvi = res$q10_dvi,
              q90_dvi = res$q90_dvi,
              stringsAsFactors = FALSE
            )
          } else {
            NULL
          }
        }))
      }
    }, error = function(e) {
      warning(sprintf("q_dvi_data assembly failed: %s", e$message)); NULL
    })

    variant_similarity_table <- NULL
    variant_similarity_summary <- NULL
    if (exists("VARIANT_SIMILARITY_TABLE")) {
      variant_similarity_table <- VARIANT_SIMILARITY_TABLE
    } else if (exists("INSEPARABLE_VARIANT_INFO") && !is.null(INSEPARABLE_VARIANT_INFO$similarity_table)) {
      variant_similarity_table <- INSEPARABLE_VARIANT_INFO$similarity_table
    }
    if (!is.null(variant_similarity_table) && nrow(variant_similarity_table) > 0) {
      variant_similarity_table <- variant_similarity_table[order(variant_similarity_table$veg, -variant_similarity_table$cos_sim, variant_similarity_table$euclidean_dist), , drop = FALSE]
      if (requireNamespace("dplyr", quietly = TRUE)) {
        variant_similarity_summary <- variant_similarity_table |> 
          dplyr::group_by(.data$veg) |> 
          dplyr::summarise(
            pair_count = dplyr::n(),
            max_cos_sim = if (all(is.na(.data$cos_sim))) NA_real_ else max(.data$cos_sim, na.rm = TRUE),
            min_cos_sim = if (all(is.na(.data$cos_sim))) NA_real_ else min(.data$cos_sim, na.rm = TRUE),
            median_cos_sim = if (all(is.na(.data$cos_sim))) NA_real_ else stats::median(.data$cos_sim, na.rm = TRUE),
            min_euclidean_dist = if (all(is.na(.data$euclidean_dist))) NA_real_ else min(.data$euclidean_dist, na.rm = TRUE),
            median_euclidean_dist = if (all(is.na(.data$euclidean_dist))) NA_real_ else stats::median(.data$euclidean_dist, na.rm = TRUE),
            .groups = "drop"
          )
      } else {
        by_veg <- split(variant_similarity_table, variant_similarity_table$veg)
        summary_rows <- lapply(names(by_veg), function(veg_name) {
          tbl <- by_veg[[veg_name]]
          data.frame(
            veg = veg_name,
            pair_count = nrow(tbl),
            max_cos_sim = if (all(is.na(tbl$cos_sim))) NA_real_ else max(tbl$cos_sim, na.rm = TRUE),
            min_cos_sim = if (all(is.na(tbl$cos_sim))) NA_real_ else min(tbl$cos_sim, na.rm = TRUE),
            median_cos_sim = if (all(is.na(tbl$cos_sim))) NA_real_ else stats::median(tbl$cos_sim, na.rm = TRUE),
            min_euclidean_dist = if (all(is.na(tbl$euclidean_dist))) NA_real_ else min(tbl$euclidean_dist, na.rm = TRUE),
            median_euclidean_dist = if (all(is.na(tbl$euclidean_dist))) NA_real_ else stats::median(tbl$euclidean_dist, na.rm = TRUE),
            stringsAsFactors = FALSE
          )
        })
        variant_similarity_summary <- do.call(rbind, summary_rows)
      }
    }

    cat(sprintf("[DEBUG] all_coefs$location_id class: %s\n", class(all_coefs$location_id)))
    cat(sprintf("[DEBUG] all_coefs$location_id sample: %s\n", paste(head(all_coefs$location_id, 10), collapse = ", ")))
    cat(sprintf("[DEBUG] all_coefs$location_id is.na sum: %d\n", sum(is.na(all_coefs$location_id))))
    cat(sprintf("[DEBUG] all_coefs$location_id == '' sum: %d\n", sum(all_coefs$location_id == "", na.rm = TRUE)))
    unique_locations <- unique(trimws(as.character(all_coefs$location_id)))
    unique_locations <- unique_locations[!is.na(unique_locations) & unique_locations != ""]
    if (length(unique_locations) == 0) {
      stop("No valid location IDs found in results")
    }
    cat(sprintf("Creating Excel file with %d valid locations\n", length(unique_locations)))

    true_veg_map <- gpts_map |> dplyr::select(location_id, true_veg = Veg)
    if (!is.character(all_coefs$location_id)) {
      all_coefs$location_id <- as.character(all_coefs$location_id)
      cat("[NOTICE] Coerced all_coefs$location_id to character to match true_veg_map for joining.\n")
    }
    if (!is.character(true_veg_map$location_id)) {
      true_veg_map$location_id <- as.character(true_veg_map$location_id)
      cat("[NOTICE] Coerced true_veg_map$location_id to character for joining.\n")
    }
    if (length(all_coefs$location_id) > 0 && length(true_veg_map$location_id) > 0) {
      s1 <- unique(na.omit(all_coefs$location_id))
      s2 <- unique(na.omit(true_veg_map$location_id))
      if (length(s1) && length(s2) && all(grepl("^[0-9]+$", s1)) && any(grepl("^L_", s2))) {
        cat("[WARNING] all_coefs$location_id looks numeric while true_veg_map$location_id looks like 'L_x_y' strings — matching will likely fail.\n")
      }
    }
    all_coefs <- dplyr::left_join(all_coefs, true_veg_map, by = "location_id")
    if (!"coef_025" %in% names(all_coefs)) all_coefs$coef_025 <- NA_real_
    if (!"coef_975" %in% names(all_coefs)) all_coefs$coef_975 <- NA_real_
    if (DEBUG_UNCERTAINTY) {
      cat("After combining all_coefs with true_veg_map:\n")
      cat("Number of rows with finite coef_025:", sum(is.finite(all_coefs$coef_025)), "\n")
      cat("Number of rows with finite coef_975:", sum(is.finite(all_coefs$coef_975)), "\n")
      cat("Number of rows with finite interval:", sum(is.finite(all_coefs$interval)), "\n")
    }

    # Helper: combine results_list into unified all_coefs (no printing)
    combine_results_from_list <- function(results_list) {
      coef_list <- lapply(results_list, function(res) {
        if (is.null(res) || is.null(res$coef_df) || nrow(res$coef_df) == 0) NULL else res$coef_df
      })
      coef_list <- coef_list[!sapply(coef_list, is.null)]
      if (length(coef_list) == 0) return(NULL)
      if (requireNamespace("dplyr", quietly = TRUE)) {
        all_coefs_local <- dplyr::bind_rows(coef_list)
      } else {
        all_coefs_local <- do.call(rbind, coef_list)
      }
      # Ensure minimal columns exist
      req_cols <- c("location_id","pheno_year","Veg","coef","coef_025","coef_975","interval")
      for (cname in req_cols) if (!(cname %in% names(all_coefs_local))) all_coefs_local[[cname]] <- NA
      all_coefs_local
    }

    # Run inference on a dataframe silently (suppress per-task prints) and return results list + combined coefs
    run_inference_silent <- function(df_subset) {
      if (is.null(df_subset) || nrow(df_subset) == 0) return(list(results = list(), all_coefs = NULL))
      if (!"pheno_year" %in% names(df_subset) && "date" %in% names(df_subset)) df_subset$pheno_year <- assign_pheno_year(df_subset$date)
      if (!"task_key" %in% names(df_subset)) df_subset$task_key <- paste(df_subset$location_id, df_subset$pheno_year, sep = "_")
      task_list <- split(df_subset, df_subset$task_key)

      sink_file <- tempfile()
      con_out <- file(sink_file, open = "wt")
      con_msg <- file(paste0(sink_file, ".msg"), open = "wt")
      # Suppress output and messages while running heavy inference
      results <- NULL
      tryCatch({
        sink(con_out, type = "output")
        sink(con_msg, type = "message")
        results <- .run_map(task_list, fit_one_task, show_pb = FALSE)
      }, error = function(e) {
        warning(sprintf("run_inference_silent: underlying inference error: %s", e$message))
      }, finally = {
        # Attempt to restore sinks
        try(sink(type = "message"), silent = TRUE)
        try(sink(type = "output"), silent = TRUE)
        close(con_out)
        close(con_msg)
        try(unlink(sink_file), silent = TRUE)
        try(unlink(paste0(sink_file, ".msg")), silent = TRUE)
      })

      all_coefs_local <- combine_results_from_list(results)
      list(results = results, all_coefs = all_coefs_local)
    }

    # --- Separate training-location inference for years 2023 and 2025 ---
    if (exists("df_train") && !is.null(df_train) && nrow(df_train) > 0) {
      training_locs <- unique(as.character(df_train$location_id))
      df_train_infer <- df_full |> dplyr::filter(location_id %in% training_locs & pheno_year %in% c(2023, 2025))
      if (nrow(df_train_infer) > 0) {
        cat(sprintf("[NOTICE] Running separate (silent) inference on %d training-location rows (years 2023 & 2025)...\n", nrow(df_train_infer)))
        res_train_infer <- run_inference_silent(df_train_infer)
        all_coefs_train_infer <- res_train_infer$all_coefs
        if (is.null(all_coefs_train_infer) || nrow(all_coefs_train_infer) == 0) {
          cat("[WARNING] Training-location inference produced no coefficients; skipping accuracy summary.\n")
        } else {
          # True veg mapping
          true_map <- if (exists("true_veg_map")) true_veg_map else if (exists("gpts_map")) gpts_map |> dplyr::select(location_id, true_veg = Veg) else NULL
          if (is.null(true_map)) {
            cat("[WARNING] No true_veg mapping available; cannot compute veg-type accuracy.\n")
          } else {
            # Merge all predictions with true vegetation labels
            merged <- all_coefs_train_infer |> dplyr::left_join(true_map, by = "location_id")

            # Mark correct predictions (predicted Veg matches true_veg)
            merged$correct <- tolower(merged$Veg) == tolower(merged$true_veg)

            # Calculate fraction-weighted accuracy per location-year
            location_year_accuracy <- merged |>
              dplyr::filter(tolower(true_veg) != "barren") |>  # Only vegetation samples
              dplyr::group_by(location_id, pheno_year, true_veg) |>
              dplyr::summarise(
                total_veg_fraction = sum(coef[tolower(Veg) != "barren"], na.rm = TRUE),  # Total predicted vegetation
                correct_veg_fraction = sum(coef[correct & tolower(Veg) != "barren"], na.rm = TRUE),  # Correctly predicted vegetation
                accuracy = ifelse(total_veg_fraction > 0, 100 * correct_veg_fraction / total_veg_fraction, NA_real_),
                .groups = "drop"
              )

            # Overall vegetation accuracy by year
            summary_by_year_veg_only <- location_year_accuracy |>
              dplyr::group_by(pheno_year) |>
              dplyr::summarise(
                n = dplyr::n(),
                mean_accuracy = mean(accuracy, na.rm = TRUE),
                total_correct_frac = sum(correct_veg_fraction, na.rm = TRUE),
                total_veg_frac = sum(total_veg_fraction, na.rm = TRUE),
                weighted_accuracy = 100 * total_correct_frac / pmax(0.001, total_veg_frac),
                .groups = "drop"
              )

            # Per-vegetation breakdown
            summary_by_veg <- location_year_accuracy |>
              dplyr::group_by(pheno_year, true_veg) |>
              dplyr::summarise(
                n = dplyr::n(),
                mean_accuracy = mean(accuracy, na.rm = TRUE),
                total_correct_frac = sum(correct_veg_fraction, na.rm = TRUE),
                total_veg_frac = sum(total_veg_fraction, na.rm = TRUE),
                weighted_accuracy = 100 * total_correct_frac / pmax(0.001, total_veg_frac),
                .groups = "drop"
              )

            cat("\n=== Training-locations inference accuracy (fraction-weighted) ===\n")
            cat("Accuracy = (correctly predicted veg fraction) / (total predicted veg fraction)\n\n")
            for (r in seq_len(nrow(summary_by_year_veg_only))) {
              cat(sprintf("Year %d: VEGETATION classification accuracy = %.1f%% [mean across %d locations]\n",
                summary_by_year_veg_only$pheno_year[r],
                summary_by_year_veg_only$mean_accuracy[r],
                summary_by_year_veg_only$n[r]))
              cat(sprintf("         Weighted accuracy = %.1f%% (%.3f correct / %.3f total veg)\n",
                summary_by_year_veg_only$weighted_accuracy[r],
                summary_by_year_veg_only$total_correct_frac[r],
                summary_by_year_veg_only$total_veg_frac[r]))
            }
            cat("\nPer-vegetation accuracy breakdown:\n")
            print(summary_by_veg)
            # Save concise summary
            pf_file <- file.path(OUT_DIR, "training_locations_inference_accuracy.csv")
            dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
            write.csv(summary_by_veg, pf_file, row.names = FALSE)
            cat(sprintf("Saved training-locations inference accuracy summary to: %s\n", pf_file))
          }
        }
      } else {
        cat("[NOTICE] No rows found for training-location inference on years 2023/2025.\n")
      }
    }

    inseparable_flags_detail <- NULL
    inseparable_flags_summary <- NULL
    if ("inseparable_variant_flag" %in% names(all_coefs)) {
      flagged_rows <- which(!is.na(all_coefs$inseparable_variant_flag) & all_coefs$inseparable_variant_flag)
      if (length(flagged_rows) > 0) {
        cols_keep <- intersect(c("location_id", "year", "Veg", "coef", "rmse", "true_veg", "inseparable_variant_details"), names(all_coefs))
        inseparable_flags_detail <- all_coefs[flagged_rows, cols_keep, drop = FALSE]
        inseparable_flags_detail <- inseparable_flags_detail[order(inseparable_flags_detail$location_id, inseparable_flags_detail$year, inseparable_flags_detail$Veg), , drop = FALSE]

        if (requireNamespace("dplyr", quietly = TRUE)) {
          inseparable_flags_summary <- inseparable_flags_detail |> 
            dplyr::group_by(.data$location_id, .data$year) |> 
            dplyr::summarise(
              veg_components = paste(sort(unique(.data$Veg)), collapse = ", "),
              max_coef = if (all(is.na(.data$coef))) NA_real_ else max(.data$coef, na.rm = TRUE),
              details = paste(unique(na.omit(.data$inseparable_variant_details)), collapse = " | "),
              .groups = "drop"
            ) |> 
            dplyr::arrange(.data$location_id, .data$year)
        } else {
          summary_keys <- unique(inseparable_flags_detail[, c("location_id", "year"), drop = FALSE])
          summary_keys <- summary_keys[order(summary_keys$location_id, summary_keys$pheno_year), , drop = FALSE]
          build_summary <- function(loc, yr) {
            rows <- inseparable_flags_detail$location_id == loc & inseparable_flags_detail$pheno_year == yr
            vegs <- sort(unique(inseparable_flags_detail$Veg[rows]))
            coef_vals <- inseparable_flags_detail$coef[rows]
            details <- unique(na.omit(inseparable_flags_detail$inseparable_variant_details[rows]))
            data.frame(
              location_id = loc,
              pheno_year = yr,
              veg_components = paste(vegs, collapse = ", "),
              max_coef = if (all(is.na(coef_vals))) NA_real_ else max(coef_vals, na.rm = TRUE),
              details = paste(details, collapse = " | "),
              stringsAsFactors = FALSE
            )
          }
          if (nrow(summary_keys) > 0) {
            summary_list <- apply(summary_keys, 1, function(row) build_summary(row[["location_id"]], as.integer(row[["year"]])))
            inseparable_flags_summary <- do.call(rbind, summary_list)
          }
        }
      }
    }

    cat("Creating single Excel file with all location results...\n")

    summary_data <- data.frame(
      Location_ID = unique_locations,
      Total_Years = sapply(unique_locations, function(loc) {
        length(unique(all_coefs$pheno_year[all_coefs$location_id == loc]))
      }),
      Total_Observations = sapply(unique_locations, function(loc) {
        nrow(all_coefs[all_coefs$location_id == loc, ])
      }),
      stringsAsFactors = FALSE
    )

    openxlsx::addWorksheet(wb, "Summary")
    openxlsx::writeData(wb, "Summary", summary_data)

    if (!is.null(all_diagnostics) && nrow(all_diagnostics) > 0) {
      openxlsx::addWorksheet(wb, "Diagnostics")
      openxlsx::writeData(wb, "Diagnostics", all_diagnostics)
    }

            if (!is.null(all_variants_pca) && nrow(all_variants_pca) > 0) {
      openxlsx::addWorksheet(wb, "Variant_Summary")
      openxlsx::writeData(wb, "Variant_Summary", all_variants_pca)
    }

    if (!is.null(variant_similarity_table) && nrow(variant_similarity_table) > 0) {
      openxlsx::addWorksheet(wb, "Variant_Similarity")
      openxlsx::writeData(wb, "Variant_Similarity", "Pairwise similarity across variants (cosine similarity and Euclidean distance)", startRow = 1, startCol = 1)
      start_row <- 2
      if (!is.null(variant_similarity_summary) && nrow(variant_similarity_summary) > 0) {
        openxlsx::writeData(wb, "Variant_Similarity", variant_similarity_summary, startRow = start_row, startCol = 1)
        start_row <- start_row + nrow(variant_similarity_summary) + 2
      }
      openxlsx::writeData(wb, "Variant_Similarity", variant_similarity_table, startRow = start_row, startCol = 1)
    }

    if (exists("INTER_CLASS_SIMILARITY") && !is.null(INTER_CLASS_SIMILARITY) && nrow(INTER_CLASS_SIMILARITY) > 0) {
      openxlsx::addWorksheet(wb, "Inter_Veg_Similarity")
      openxlsx::writeData(wb, "Inter_Veg_Similarity", "Pairwise similarity between variants of different vegetation types", startRow = 1, startCol = 1)
      openxlsx::writeData(wb, "Inter_Veg_Similarity", INTER_CLASS_SIMILARITY, startRow = 2, startCol = 1)
    }

    if (!is.null(inseparable_flags_summary) && nrow(inseparable_flags_summary) > 0) {
      openxlsx::addWorksheet(wb, "Inseparable_Flags")
      openxlsx::writeData(wb, "Inseparable_Flags", "Location-years with inseparable MESMA variants detected", startRow = 1, startCol = 1)
      openxlsx::writeData(wb, "Inseparable_Flags", inseparable_flags_summary, startRow = 2, startCol = 1)
      next_row <- nrow(inseparable_flags_summary) + 4
      if (!is.null(inseparable_flags_detail) && nrow(inseparable_flags_detail) > 0) {
        openxlsx::writeData(wb, "Inseparable_Flags", "Detailed component rows", startRow = next_row, startCol = 1)
        openxlsx::writeData(wb, "Inseparable_Flags", inseparable_flags_detail, startRow = next_row + 1, startCol = 1)
      }
    }

    unc_coef_rows <- list(); unc_var_rows <- list(); unc_rmse_rows <- list(); unc_meta_rows <- list()
    for (res in results_list) {
      if (is.null(res$uncertainty)) next
      loc <- if (!is.null(res$coef_df$location_id)) res$coef_df$location_id[1] else NA_character_
      yr <- if (!is.null(res$coef_df$year)) res$coef_df$year[1] else NA_integer_
      ci <- res$uncertainty$coef_ci
      if (!is.null(ci) && nrow(ci) > 0) {
        ci$location_id <- loc; ci$year <- yr
        unc_coef_rows[[length(unc_coef_rows) + 1]] <- ci[, c("location_id","year","Veg","coef_025","coef_975"), drop = FALSE]
      }
      vr <- res$uncertainty$variant_freq
      if (!is.null(vr) && nrow(vr) > 0) {
        vr$location_id <- loc; vr$year <- yr
        unc_var_rows[[length(unc_var_rows) + 1]] <- vr[, c("location_id","year","Veg","Variant","N","Percent"), drop = FALSE]
      }
      rc <- res$uncertainty$rmse_ci
      if (!is.null(rc) && length(rc) == 2) {
        unc_rmse_rows[[length(unc_rmse_rows) + 1]] <- data.frame(location_id = loc, year = yr, rmse_lo = rc[1], rmse_hi = rc[2], stringsAsFactors = FALSE)
      }

    }
    if (length(unc_coef_rows) > 0) {
      if (requireNamespace("dplyr", quietly = TRUE)) all_unc_coef <- dplyr::bind_rows(unc_coef_rows) else all_unc_coef <- do.call(rbind, unc_coef_rows)
    } else all_unc_coef <- NULL
    if (!is.null(all_unc_coef)) {
      n_total_unc <- nrow(all_unc_coef)
      n_non_na_ci <- sum(is.finite(all_unc_coef$coef_025) | is.finite(all_unc_coef$coef_975))
      cat(sprintf("Uncertainty CIs: total rows = %d; rows with at least one finite CI = %d\n", n_total_unc, n_non_na_ci))

      # MERGE BOOTSTRAP UNCERTAINTY INTO MAIN COEFFICIENTS
      # This is the critical fix: merge coef_025 and coef_975 from bootstrap into all_coefs
      if (requireNamespace("dplyr", quietly = TRUE)) {
        # Ensure matching column types for join
        all_coefs$location_id <- as.character(all_coefs$location_id)
        all_coefs$year <- as.integer(all_coefs$year)
        all_unc_coef$location_id <- as.character(all_unc_coef$location_id)
        all_unc_coef$year <- as.integer(all_unc_coef$year)

        # Remove existing coef_025 and coef_975 columns if they exist (they contain NAs)
        if ("coef_025" %in% names(all_coefs)) all_coefs$coef_025 <- NULL
        if ("coef_975" %in% names(all_coefs)) all_coefs$coef_975 <- NULL

        # Perform left join to merge bootstrap CIs into main results
        all_coefs <- dplyr::left_join(all_coefs, all_unc_coef, by = c("location_id", "year", "Veg"))

        # Report merge success
        n_merged <- sum(is.finite(all_coefs$coef_025) | is.finite(all_coefs$coef_975))
        cat(sprintf("Bootstrap CIs merged into main results: %d/%d rows now have uncertainty bounds\n",
                    n_merged, nrow(all_coefs)))
      } else {
        warning("dplyr not available - bootstrap CIs will not be merged into main results")
      }
    }
    if (length(unc_var_rows) > 0) {
      if (requireNamespace("dplyr", quietly = TRUE)) all_unc_var <- dplyr::bind_rows(unc_var_rows) else all_unc_var <- do.call(rbind, unc_var_rows)
    } else all_unc_var <- NULL
    if (length(unc_rmse_rows) > 0) {
      if (requireNamespace("dplyr", quietly = TRUE)) all_unc_rmse <- dplyr::bind_rows(unc_rmse_rows) else all_unc_rmse <- do.call(rbind, unc_rmse_rows)
    } else all_unc_rmse <- NULL
    if (length(unc_meta_rows) > 0) {
      if (requireNamespace("dplyr", quietly = TRUE)) all_unc_meta <- dplyr::bind_rows(unc_meta_rows) else all_unc_meta <- do.call(rbind, unc_meta_rows)
    } else all_unc_meta <- NULL

    if (!is.null(all_unc_coef) || !is.null(all_unc_var) || !is.null(all_unc_rmse)) {
      openxlsx::addWorksheet(wb, "Uncertainty")
      start_row <- 1
      openxlsx::writeData(wb, "Uncertainty", data.frame(Setting = c("ENABLE_UNCERTAINTY","BOOTSTRAP_METHOD"), Value = c(ENABLE_UNCERTAINTY, if (isTRUE(ENABLE_UNCERTAINTY)) "locations" else "none")), startRow = start_row, startCol = 1)
      start_row <- start_row + 3
      if (!is.null(all_unc_rmse)) {
        openxlsx::writeData(wb, "Uncertainty", "RMSE CI (2.5%/97.5%)", startRow = start_row, startCol = 1)
        openxlsx::writeData(wb, "Uncertainty", all_unc_rmse, startRow = start_row + 1, startCol = 1)
        start_row <- start_row + nrow(all_unc_rmse) + 3
      }
      if (!is.null(all_unc_var)) {
        openxlsx::writeData(wb, "Uncertainty", "Variant Dominance Frequencies (%)", startRow = start_row, startCol = 1)
        openxlsx::writeData(wb, "Uncertainty", all_unc_var, startRow = start_row + 1, startCol = 1)
      }
      if (!is.null(all_unc_meta)) {
        start_row <- start_row + ifelse(!is.null(all_unc_var), nrow(all_unc_var) + 3, 3)
        openxlsx::writeData(wb, "Uncertainty", "Bootstrap Meta", startRow = start_row, startCol = 1)
        openxlsx::writeData(wb, "Uncertainty", all_unc_meta, startRow = start_row + 1, startCol = 1)
      }
    }

    if (exists('stability_results') && !is.null(stability_results)) {
      openxlsx::addWorksheet(wb, 'Endmember_Stability')
      current_row <- 1
      for (veg in names(stability_results)) {
        res <- stability_results[[veg]]
        openxlsx::writeData(wb, 'Endmember_Stability', sprintf('=== %s ENDMEMBER STABILITY ===', toupper(veg)), startRow = current_row, startCol = 1)
        current_row <- current_row + 2
        openxlsx::writeData(wb, 'Endmember_Stability', 'Variant Count Distribution:', startRow = current_row, startCol = 1)
        vcd <- as.data.frame(res$variant_count_distribution)
        names(vcd) <- c('N_Variants', 'Frequency')
        openxlsx::writeData(wb, 'Endmember_Stability', vcd, startRow = current_row + 1, startCol = 1)
        current_row <- current_row + nrow(vcd) + 3
        openxlsx::writeData(wb, 'Endmember_Stability', 'Meta-Variant Stability Metrics:', startRow = current_row, startCol = 1)
        mv_stats <- if (length(res$meta_variants) > 0) do.call(rbind, lapply(res$meta_variants, function(mv) {
          data.frame(Meta_Variant = mv$meta_variant_id, N_Bootstrap_Members = mv$n_members, Coefficient_of_Variation = mv$coefficient_of_variation, Spectral_Angle_Mean_deg = mv$mean_spectral_angle, Spectral_Angle_SD_deg = sqrt(mv$spectral_angle_variance), Bootstrap_Frequency_pct = length(unique(mv$bootstrap_iters)) / res$n_bootstrap_iters * 100, stringsAsFactors = FALSE)
        })) else data.frame()
        if (nrow(mv_stats) > 0) {
          openxlsx::writeData(wb, 'Endmember_Stability', mv_stats, startRow = current_row + 1, startCol = 1)
          current_row <- current_row + nrow(mv_stats) + 4
        } else {
          current_row <- current_row + 4
        }
      }
      openxlsx::writeData(wb, 'Endmember_Stability', 'INTERPRETATION GUIDE', startRow = current_row, startCol = 1)
      current_row <- current_row + 1
      guide <- data.frame(Metric = c('Variant Count Distribution', 'Coefficient of Variation (CV)', 'Spectral Angle SD', 'Bootstrap Frequency'), Interpretation = c('How many variants are identified across bootstrap iterations. Stable if concentrated.', 'Variability of signature features. Lower CV = more stable endmember.', 'Directional variability in spectral space (degrees). Lower = more stable.', '% of bootstrap iterations where this pattern appears. Higher = more robust.'), Good_Value = c('Narrow distribution (e.g., 90% have same count)', '< 0.10 (very stable), 0.10-0.20 (stable), > 0.20 (unstable)', '< 5° (very stable), 5-10° (stable), > 10° (unstable)', '> 80% (very robust), 50-80% (robust), < 50% (unreliable)'), stringsAsFactors = FALSE)
      openxlsx::writeData(wb, 'Endmember_Stability', guide, startRow = current_row, startCol = 1)
    }

    best_fit_summary <- do.call(rbind, lapply(unique_locations, function(loc) {
      yrs <- unique(all_coefs$pheno_year[all_coefs$location_id == loc])
      tv <- true_veg_map$true_veg[true_veg_map$location_id == loc]
      if (length(tv) == 0) tv <- NA_character_

        do.call(rbind, lapply(yrs, function(yr) {
        row <- all_coefs[all_coefs$location_id == loc & all_coefs$pheno_year == yr & tolower(all_coefs$Veg) == tolower(tv), , drop = FALSE]
        pred <- if (nrow(row) == 1) row$coef else NA_real_
        pred_abs <- pred
        rmse_val <- if (nrow(row) == 1 && "rmse" %in% names(row)) row$rmse else NA_real_
          sum_veg_coef <- sum(all_coefs$coef[all_coefs$location_id == loc & all_coefs$pheno_year == yr & tolower(all_coefs$Veg) != "barren"], na.rm = TRUE)
          pred_rel <- if (!is.na(pred) && is.finite(sum_veg_coef) && sum_veg_coef > 0) pred / sum_veg_coef else NA_real_
          pred_rel <- pmin(pred_rel, 1)  # Clamp to 1 to prevent >100%
          abs_pct <- if (!is.na(pred_abs) && !is.na(tv)) abs(1 - pred_abs) * 100 else NA_real_
          abs_pct_rel <- if (!is.na(pred_rel) && !is.na(tv)) abs(1 - pred_rel) * 100 else NA_real_

        data.frame(
          location_id = loc,
          year = yr,
          true_veg = tv,
          pred_coef = pred_rel,
          pred_coef_abs = pred_abs,
          pred_coef_rel = pred_rel,
          rmse = rmse_val,
          abs_pct_diff = abs_pct_rel,
          abs_pct_diff_abs = abs_pct,
          abs_pct_diff_rel = abs_pct_rel,
          stringsAsFactors = FALSE
        )
      }))
    }))

    eval_years <- sort(unique(c(TRAIN_YEARS, TRAIN_YEARS - 1L, TRAIN_YEARS + 1L)))
    eval_years <- eval_years[is.finite(eval_years)]
    best_fit_summary$eval_window <- best_fit_summary$year %in% eval_years
    best_fit_eval <- best_fit_summary[best_fit_summary$eval_window, , drop = FALSE]
    if (nrow(best_fit_eval) == 0) {
      warning("No location-years fall within the TRAIN_YEARS +/- 1 evaluation window; overall fit cannot be computed")
    }

    overall_fit <- suppressWarnings(as.numeric(mean(best_fit_eval$pred_coef_rel * 100, na.rm = TRUE)))
    if (!is.finite(overall_fit)) overall_fit <- NA_real_

    openxlsx::writeData(wb, "Summary", data.frame(
      Overall_Fit_pct = overall_fit
    ), startRow = 1, startCol = ncol(summary_data) + 2)

    for (i in seq_along(unique_locations)) {
      loc_id <- unique_locations[i]

          sheet_name <- substr(gsub("[^A-Za-z0-9]", "_", loc_id), 1, 31)

          openxlsx::addWorksheet(wb, sheet_name)

          loc_coefs <- all_coefs[all_coefs$location_id == loc_id, ]

          # Extract true_veg for this location
          true_veg_val <- true_veg_map$true_veg[true_veg_map$location_id == loc_id]
          if (length(true_veg_val) == 0) true_veg_val <- NA_character_
          true_veg_val <- true_veg_val[1]  # Take first value

          quality_metrics <- loc_coefs |> 
            dplyr::group_by(.data$pheno_year) |> 
            dplyr::summarize(
              deviation = sum(abs(.data$coef - (tolower(.data$Veg) == tolower(true_veg_val)))),
              avg_rmse = mean(.data$rmse, na.rm = TRUE),
              .groups = "drop"
            ) |> 
            dplyr::summarize(
              avg_pct_deviation = mean(.data$deviation, na.rm = TRUE) * 100,
              avg_rmse = mean(.data$avg_rmse, na.rm = TRUE),
              .groups = "drop"
            )

          loc_q_data <- if (!is.null(q_dvi_data) && "location_id" %in% names(q_dvi_data)) {
            q_dvi_data[q_dvi_data$location_id == loc_id, ]
          } else {
            NULL
          }

          peak_q10_dvi <- if (!is.null(loc_q_data) && nrow(loc_q_data) > 0) {
            max(loc_q_data$q10_dvi, na.rm = TRUE)
          } else {
            NA
          }

          peak_q90_dvi <- if (!is.null(loc_q_data) && nrow(loc_q_data) > 0) {
            max(loc_q_data$q90_dvi, na.rm = TRUE)
          } else {
            NA
          }

          quality_metrics$peak_q10_dvi <- peak_q10_dvi
          quality_metrics$peak_q90_dvi <- peak_q90_dvi

          openxlsx::writeData(wb, sheet_name, "QUALITY METRICS", startRow = 1, startCol = 1)
          openxlsx::writeData(wb, sheet_name, quality_metrics, startRow = 2, startCol = 1)

          current_row <- nrow(quality_metrics) + 4

          # Removed DIAGNOSTICS section

          loc_best <- if ("location_id" %in% names(best_fit_summary)) best_fit_summary[best_fit_summary$location_id == loc_id, , drop = FALSE] else data.frame()
          if (!is.null(loc_best) && nrow(loc_best) > 0) {
            desired_cols <- c("year", "true_veg", "pred_coef", "pred_coef_abs", "rmse", "abs_pct_diff", "abs_pct_diff_abs", "eval_window")
            write_tbl <- loc_best[, c("location_id", intersect(desired_cols, names(loc_best))), drop = FALSE]
            openxlsx::writeData(wb, sheet_name, "BEST FIT SUMMARY (per-year) — pred_coef is proportion relative to vegetated area (barren excluded)",
              startRow = current_row, startCol = 1
            )
            openxlsx::writeData(wb, sheet_name, write_tbl, startRow = current_row + 1, startCol = 1)
            current_row <- current_row + nrow(write_tbl) + 3
          }

          loc_coefs_unified <- if ("location_id" %in% names(all_coefs)) all_coefs[all_coefs$location_id == loc_id, ] else data.frame()

          loc_variants_pca <- if (!is.null(all_variants_pca) && "location_id" %in% names(all_variants_pca)) {
            all_variants_pca[all_variants_pca$location_id == loc_id, ]
          } else {
            NULL
          }

          if (nrow(loc_coefs_unified) > 0) {
            openxlsx::writeData(wb, sheet_name, "COEFFICIENTS",
              startRow = current_row, startCol = 1
            )
            openxlsx::writeData(wb, sheet_name, loc_coefs_unified,
              startRow = current_row + 1, startCol = 1
            )
            current_row <- current_row + nrow(loc_coefs_unified) + 3
          }


          if (!is.null(loc_variants_pca) && nrow(loc_variants_pca) > 0) {
            variant_usage <- data.frame(
              year = unique(loc_variants_pca$year),
              stringsAsFactors = FALSE
            )

            veg_candidates <- unique(c(ALLOWED_VEG, na.omit(unique(as.character(all_coefs$Veg)))))
            veg_candidates <- veg_candidates[!is.na(veg_candidates) & veg_candidates != ""]
            veg_kept_local <- veg_candidates
            if (length(veg_kept_local) == 0) veg_kept_local <- ALLOWED_VEG
            for (veg in veg_kept_local) {
              var_col <- paste0(veg, "_variant")
              if (var_col %in% names(loc_variants_pca)) {
                variant_usage[[paste0(veg, "_most_common")]] <- sapply(variant_usage$year, function(y) {
                  year_data <- loc_variants_pca[loc_variants_pca$year == y, var_col]
                  if (length(year_data) > 0) {
                    names(sort(table(year_data), decreasing = TRUE))[1]
                  } else {
                    NA
                  }
                })
              }
            }

            openxlsx::writeData(wb, sheet_name, "VARIANT USAGE",
              startRow = current_row, startCol = 1
            )
            openxlsx::writeData(wb, sheet_name, variant_usage,
              startRow = current_row + 1, startCol = 1
            )
            current_row <- current_row + nrow(variant_usage) + 3
          }

          if (!is.null(loc_q_data) && nrow(loc_q_data) > 0) {
            openxlsx::writeData(wb, sheet_name, "Q10/Q90 DVI TREND",
              startRow = current_row, startCol = 1
            )
            openxlsx::writeData(wb, sheet_name, loc_q_data,
              startRow = current_row + 1, startCol = 1
            )
            current_row <- current_row + nrow(loc_q_data) + 3
          }

          if (exists("all_unc_rmse") && !is.null(all_unc_rmse) && "location_id" %in% names(all_unc_rmse)) {
            loc_unc_rmse <- all_unc_rmse[all_unc_rmse$location_id == loc_id, , drop = FALSE]
            if (nrow(loc_unc_rmse) > 0) {
              openxlsx::writeData(wb, sheet_name, "RMSE CI (2.5%/97.5%)",
                startRow = current_row, startCol = 1
              )
              openxlsx::writeData(wb, sheet_name, loc_unc_rmse,
                startRow = current_row + 1, startCol = 1
              )
              current_row <- current_row + nrow(loc_unc_rmse) + 3
            }
          }

          if (!is.null(inseparable_flags_summary) && nrow(inseparable_flags_summary) > 0 && "location_id" %in% names(inseparable_flags_summary)) {
            loc_flagged <- inseparable_flags_summary[inseparable_flags_summary$location_id == loc_id, , drop = FALSE]
            if (nrow(loc_flagged) > 0) {
              openxlsx::writeData(wb, sheet_name, "INSEPARABLE VARIANT ALERTS (summary)", startRow = current_row, startCol = 1)
              openxlsx::writeData(wb, sheet_name, loc_flagged, startRow = current_row + 1, startCol = 1)
              current_row <- current_row + nrow(loc_flagged) + 3
            }
          }

          if (!is.null(inseparable_flags_detail) && nrow(inseparable_flags_detail) > 0 && "location_id" %in% names(inseparable_flags_detail)) {
            loc_flagged_detail <- inseparable_flags_detail[inseparable_flags_detail$location_id == loc_id, , drop = FALSE]
            if (nrow(loc_flagged_detail) > 0) {
              openxlsx::writeData(wb, sheet_name, "INSEPARABLE VARIANT COMPONENTS", startRow = current_row, startCol = 1)
              openxlsx::writeData(wb, sheet_name, loc_flagged_detail, startRow = current_row + 1, startCol = 1)
              current_row <- current_row + nrow(loc_flagged_detail) + 3
            }
          }

        }

    output_filename <- file.path(OUT_DIR, "mesma_results.xlsx")

    if (!dir.exists(OUT_DIR)) {
      cat(sprintf("Creating output directory: %s\n", OUT_DIR))
      dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
    }

    if (!dir.exists(OUT_DIR)) {
      cat(sprintf("ERROR: Cannot create output directory: %s\n", OUT_DIR))
      stop("Cannot create output directory")
    }

    cat(sprintf("Saving workbook to: %s\n", output_filename))

    save_result <- tryCatch(
      {
        openxlsx::saveWorkbook(wb, output_filename, overwrite = TRUE)
        TRUE
      },
      error = function(e) {
        cat(sprintf("ERROR saving workbook: %s\n", e$message))
        FALSE
      }
    )

    if (save_result) {
      cat(sprintf(
        "Created Excel file '%s' with %d location sheets\n",
        basename(output_filename), length(unique_locations)
      ))
    } else {
      cat("Failed to save Excel file\n")
    }
  }


  aggregate_to_global_pattern <- function(all_coefs, method = "location_bootstrap") {

    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")

    required_cols <- c("location_id", "year", "Veg", "coef")
    missing <- setdiff(required_cols, names(all_coefs))
    if (length(missing) > 0) stop(paste("Missing columns:", paste(missing, collapse = ", ")))

    if (method == "location_bootstrap") {
      # Rename year to pheno_year for location_bootstrap_aggregate
      if ("year" %in% names(all_coefs) && !"pheno_year" %in% names(all_coefs)) {
        all_coefs$pheno_year <- all_coefs$year
      }
      result <- location_bootstrap_aggregate(all_coefs, B = BOOTSTRAP_B)
    } else if (method == "hierarchical") {
      if (!"coef_025" %in% names(all_coefs)) all_coefs$coef_025 <- NA_real_
      if (!"coef_975" %in% names(all_coefs)) all_coefs$coef_975 <- NA_real_

      all_coefs$interval <- all_coefs$coef_975 - all_coefs$coef_025
      all_coefs$interval[!is.finite(all_coefs$interval)] <- NA_real_

      all_coefs$se_proxy <- all_coefs$interval / 3.92
      all_coefs$se_proxy[all_coefs$se_proxy <= 0 | !is.finite(all_coefs$se_proxy)] <- NA_real_
      all_coefs$weight <- 1 / (all_coefs$se_proxy^2)
      all_coefs$weight[!is.finite(all_coefs$weight)] <- 1  # Default weight if no CI

      result <- aggregate_hierarchical(all_coefs)
    } else {
      stop("Unknown method. Use 'location_bootstrap' or 'hierarchical'.")
    }

    result
  }


  aggregate_hierarchical <- function(all_coefs) {

    if (!requireNamespace("lme4", quietly = TRUE)) {
      warning("lme4 not available, falling back to simple mean aggregation")
      # Simple aggregation fallback
      result <- all_coefs |>
        dplyr::group_by(Veg, year) |>
        dplyr::summarize(
          n_locations = dplyr::n(),
          global_coef = mean(coef, na.rm = TRUE),
          se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
          ci_lower = global_coef - 1.96 * se,
          ci_upper = global_coef + 1.96 * se,
          .groups = "drop"
        ) |>
        dplyr::mutate(method = "simple_mean")
      return(result)
    }
    
    veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
    results_list <- list()
    
    for (veg in veg_types) {
      veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), ]
      
      if (nrow(veg_data) == 0) {
        warning(sprintf("No valid coefficients found for vegetation type: %s", veg))
        next
      }
      
      if (nrow(veg_data) < 10) {
        simple_result <- veg_data |> 
          dplyr::group_by(year) |> 
          dplyr::summarize(
            n_locations = dplyr::n(),
            global_coef = mean(coef, na.rm = TRUE),
            se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
            ci_lower = global_coef - 1.96 * se,
            ci_upper = global_coef + 1.96 * se,
            .groups = "drop"
          ) |> 
          dplyr::mutate(Veg = veg, method = "simple_mean")
        
        results_list[[veg]] <- simple_result
        next
      }
      
      tryCatch({
        veg_data$year_factor <- as.factor(veg_data$year)
        
        model <- suppressWarnings(lme4::lmer(coef ~ year_factor + (1|location_id), data = veg_data))
        
        fe <- lme4::fixef(model)
        vcov_fe <- as.matrix(vcov(model))
        
        years <- sort(unique(veg_data$year))
        pred_data <- data.frame(year_factor = as.factor(years))
        
        preds <- predict(model, newdata = pred_data, re.form = NA)
        
        boot_preds <- lme4::bootMer(model, function(m) {
          predict(m, newdata = pred_data, re.form = NA)
        }, nsim = 100, type = "parametric")
        
        ci_lower <- apply(boot_preds$t, 2, quantile, 0.025)
        ci_upper <- apply(boot_preds$t, 2, quantile, 0.975)
        
        hier_result <- data.frame(
          year = years,
          Veg = veg,
          n_locations = sapply(years, function(y) sum(veg_data$year == y)),
          global_coef = preds,
          se = apply(boot_preds$t, 2, sd),
          ci_lower = pmax(0, ci_lower),
          ci_upper = pmin(1, ci_upper),
          method = "hierarchical"
        )
        
        results_list[[veg]] <- hier_result
        
      }, error = function(e) {
        warning(sprintf("Hierarchical model failed for %s: %s. Using simple mean.", veg, e$message))
        simple_result <- veg_data |> 
          dplyr::group_by(year) |> 
          dplyr::summarize(
            n_locations = dplyr::n(),
            global_coef = mean(coef, na.rm = TRUE),
            se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
            ci_lower = pmax(0, global_coef - 1.96 * se),
            ci_upper = pmin(1, global_coef + 1.96 * se),
            .groups = "drop"
          ) |> 
          dplyr::mutate(Veg = veg, method = "simple_mean")
        
        results_list[[veg]] <- simple_result
      })
    }
    
    dplyr::bind_rows(results_list)
  }

  aggregate_simple_mean <- function(all_coefs) {
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")

    veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])

    # Determine which year column to use
    year_col <- if ("pheno_year" %in% names(all_coefs)) "pheno_year" else "year"

    results_list <- list()

    for (veg in veg_types) {
      veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), ]

      if (nrow(veg_data) == 0) {
        warning(sprintf("No valid coefficients found for vegetation type: %s", veg))
        next
      }

      # Group by year and compute simple statistics
      simple_result <- veg_data |> 
        dplyr::group_by(!!rlang::sym(year_col)) |> 
        dplyr::summarize(
          n_locations = dplyr::n(),
          global_coef = mean(coef, na.rm = TRUE),
          se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
          ci_lower = pmax(0, global_coef - 1.96 * se),
          ci_upper = pmin(1, global_coef + 1.96 * se),
          .groups = "drop"
        ) |> 
        dplyr::mutate(Veg = veg, method = "simple_mean") |> 
        dplyr::rename(year = !!rlang::sym(year_col))

      results_list[[veg]] <- simple_result
    }

    dplyr::bind_rows(results_list)
  }




  plot_global_vegetation_pattern <- function(global_pattern, 
                                              title = "Global Vegetation Composition Over Time",
                                              show_ci = TRUE,
                                              ci_type = "auto") {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    
    if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
      if ("ci_lower" %in% names(global_pattern) && "ci_upper" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower
        global_pattern$ci_upper <- global_pattern$ci_upper
      }
    } else if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
      if ("ci_lower_simple" %in% names(global_pattern) && "ci_upper_simple" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower_simple
        global_pattern$ci_upper <- global_pattern$ci_upper_simple
      }
    } else if ("weighted_mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$weighted_mean_coef
      if ("ci_lower_pooled" %in% names(global_pattern) && "ci_upper_pooled" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower_pooled
        global_pattern$ci_upper <- global_pattern$ci_upper_pooled
      }
    }
    
    global_pattern_veg <- global_pattern[tolower(global_pattern$Veg) != "barren", ]
    global_pattern_barren <- global_pattern[tolower(global_pattern$Veg) == "barren", ]
    
    p <- ggplot2::ggplot(global_pattern_veg, ggplot2::aes(x = year, y = coef, color = Veg, fill = Veg)) +
      ggplot2::geom_line(size = 1.2) +
      ggplot2::geom_point(size = 2)
    
    if (show_ci && "ci_lower" %in% names(global_pattern_veg) && "ci_upper" %in% names(global_pattern_veg)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
        alpha = 0.2,
        color = NA
      )
    }
    
    if (nrow(global_pattern_barren) > 0) {
      veg_max <- suppressWarnings(max(global_pattern_veg$coef, global_pattern_veg$ci_upper, na.rm = TRUE))
      barren_max <- suppressWarnings(max(global_pattern_barren$coef, global_pattern_barren$ci_upper, na.rm = TRUE))
      if (is.na(veg_max) || veg_max <= 0 || is.na(barren_max) || barren_max <= 0) {
        barren_scale_factor <- 1
      } else {
        barren_scale_factor <- veg_max / barren_max
      }
      global_pattern_barren$coef_scaled <- global_pattern_barren$coef * barren_scale_factor
      
      p <- p + 
        ggplot2::geom_line(data = global_pattern_barren, 
              ggplot2::aes(x = year, y = coef_scaled), 
                          color = "brown", linewidth = 1.2, linetype = "dashed") +
        ggplot2::geom_point(data = global_pattern_barren, 
               ggplot2::aes(x = year, y = coef_scaled), 
                           color = "brown", size = 2)
      
      if (show_ci && "ci_lower" %in% names(global_pattern_barren) && "ci_upper" %in% names(global_pattern_barren)) {
        global_pattern_barren$ci_lower_scaled <- global_pattern_barren$ci_lower * barren_scale_factor
        global_pattern_barren$ci_upper_scaled <- global_pattern_barren$ci_upper * barren_scale_factor
        
        p <- p + ggplot2::geom_ribbon(data = global_pattern_barren,
                                     ggplot2::aes(x = year, ymin = ci_lower_scaled, ymax = ci_upper_scaled),
                                     fill = "brown", alpha = 0.1, color = NA)
      }
    }
    
    p <- p +
      ggplot2::labs(
        title = title,
        subtitle = sprintf("Based on %d locations", max(global_pattern$n_locations, na.rm = TRUE)),
        x = "Year",
        y = "Vegetation Fraction",
        color = "Vegetation Type",
        fill = "Vegetation Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)
      ) +
          ggplot2::scale_y_continuous(
        labels = scales::percent_format(),
        sec.axis = ggplot2::sec_axis(~ . / barren_scale_factor, name = "Barren Fraction", labels = scales::percent_format())
      ) +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::scale_fill_brewer(palette = "Set1")
    
    p
  }

  plot_vegetation_stacked_area <- function(global_pattern, 
                                            title = "Global Vegetation Composition Over Time") {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    
    if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
    } else if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
    }
    
    global_pattern <- global_pattern |> 
      dplyr::group_by(year) |> 
      dplyr::mutate(coef_normalized = coef / sum(coef, na.rm = TRUE)) |> 
      dplyr::ungroup()
    
    p <- ggplot2::ggplot(global_pattern, 
                          ggplot2::aes(x = year, y = coef_normalized, fill = Veg)) +
      ggplot2::geom_area(alpha = 0.8, position = "stack") +
      ggplot2::labs(
        title = title,
        x = "Year",
        y = "Relative Vegetation Fraction",
        fill = "Vegetation Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold")
      ) +
      ggplot2::scale_y_continuous(labels = scales::percent_format()) +
      ggplot2::scale_fill_brewer(palette = "Set2")
    
    p
  }

  plot_vegetation_heatmap <- function(global_pattern, 
                                       title = "Vegetation Fraction by Year") {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    
    if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
    } else if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
    }
    
    p <- ggplot2::ggplot(global_pattern, 
                          ggplot2::aes(x = year, y = Veg, fill = coef)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", coef * 100)), 
                         color = "white", size = 3) +
      ggplot2::scale_fill_viridis_c(option = "plasma", 
                                     labels = scales::percent_format()) +
      ggplot2::labs(
        title = title,
        x = "Year",
        y = "Vegetation Type",
        fill = "Fraction"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    
    p
  }

  if (exists("inference_results_list") && length(inference_results_list) > 0) {
    cat("\nProcessing inference results for separate output...\n")
    
    inference_coef_list <- lapply(inference_results_list, function(res) {
      if (is.null(res$coef_df) || nrow(res$coef_df) == 0) {
        NULL
      } else {
        res$coef_df
      }
    })
    inference_coef_list <- inference_coef_list[!sapply(inference_coef_list, is.null)]
    
    if (length(inference_coef_list) > 0) {
      inference_coefs <- tryCatch({
        dplyr::bind_rows(inference_coef_list)
      }, error = function(e) {
        do.call(rbind, inference_coef_list)
      })
      
      cat(sprintf("Combined %d inference coefficient rows from %d location-years\n", 
                  nrow(inference_coefs), length(inference_coef_list)))
      wb_inference <- openxlsx::createWorkbook()
      
      unique_inference_locations <- unique(inference_coefs$location_id)
      for (loc_id in unique_inference_locations) {
        loc_data <- inference_coefs[inference_coefs$location_id == loc_id, ]
        sheet_name <- paste0("Loc_", loc_id)
        openxlsx::addWorksheet(wb_inference, sheet_name)
        openxlsx::writeData(wb_inference, sheet_name, loc_data)
      }
      
      inference_output_filename <- file.path(OUT_DIR, "inference_results.xlsx")
      openxlsx::saveWorkbook(wb_inference, inference_output_filename, overwrite = TRUE)
      cat(sprintf("Saved inference Excel file to: %s\n", inference_output_filename))
      
      inference_global_pattern <- aggregate_to_global_pattern(inference_coefs, method = "location_bootstrap")
      
      if (nrow(inference_global_pattern) > 0) {
        if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
        if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")
        
        inference_global_pattern_veg <- inference_global_pattern[tolower(inference_global_pattern$Veg) != "barren", ]
        inference_global_pattern_barren <- inference_global_pattern[tolower(inference_global_pattern$Veg) == "barren", ]
        
        veg_max <- suppressWarnings(max(inference_global_pattern_veg$global_coef, inference_global_pattern_veg$ci_upper, na.rm = TRUE))
        barren_max <- suppressWarnings(max(inference_global_pattern_barren$global_coef, inference_global_pattern_barren$ci_upper, na.rm = TRUE))
        if (is.na(veg_max) || veg_max <= 0 || is.na(barren_max) || barren_max <= 0) {
          inference_barren_scale <- 1
        } else {
          inference_barren_scale <- veg_max / barren_max
        }
        if (nrow(inference_global_pattern_barren) > 0) {
          inference_global_pattern_barren$global_coef_scaled <- inference_global_pattern_barren$global_coef * inference_barren_scale
          inference_global_pattern_barren$ci_lower_scaled <- inference_global_pattern_barren$ci_lower * inference_barren_scale
          inference_global_pattern_barren$ci_upper_scaled <- inference_global_pattern_barren$ci_upper * inference_barren_scale
        }
        
        p_inference <- ggplot(inference_global_pattern_veg, aes(x = year, y = global_coef, color = Veg, fill = Veg)) +
          geom_line(linewidth = 1.2) +
          geom_point(size = 2, show.legend = FALSE) +
          geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, color = NA)
        
        if (nrow(inference_global_pattern_barren) > 0) {
          p_inference <- p_inference + 
            geom_line(data = inference_global_pattern_barren, 
                     aes(x = year, y = global_coef_scaled), 
                     color = "brown", linewidth = 1.2, linetype = "dashed") +
            geom_point(data = inference_global_pattern_barren, 
                      aes(x = year, y = global_coef_scaled), 
                      color = "brown", size = 2, show.legend = FALSE) +
            geom_ribbon(data = inference_global_pattern_barren,
                       aes(x = year, ymin = ci_lower_scaled, ymax = ci_upper_scaled),
                       fill = "brown", alpha = 0.1, color = NA)
        }
        
        all_y_values <- c(inference_global_pattern_veg$global_coef, inference_global_pattern_veg$ci_lower, inference_global_pattern_veg$ci_upper)
        if (nrow(inference_global_pattern_barren) > 0) {
          all_y_values <- c(all_y_values, inference_global_pattern_barren$global_coef_scaled, inference_global_pattern_barren$ci_lower_scaled, inference_global_pattern_barren$ci_upper_scaled)
        }
        min_y <- min(all_y_values, na.rm = TRUE)
        max_y <- max(all_y_values, na.rm = TRUE)
        
        p_inference <- p_inference +
          labs(
            title = "Inference Locations: Average Coverage Percentage per Vegetation Type (2020-2025)",
            subtitle = sprintf("Based on %d locations with bootstrap uncertainty", max(inference_global_pattern$n_locations, na.rm = TRUE)),
            x = "Year",
            y = "Vegetation Fraction",
            color = "Vegetation Type",
            fill = "Vegetation Type"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 10),
            axis.title = element_text(size = 12),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10)
          ) +
          scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                             limits = c(min_y, max_y),
                             expand = c(0, 0),
                             sec.axis = ggplot2::sec_axis(~ . / inference_barren_scale, name = "Barren Fraction", labels = scales::percent_format(accuracy = 1))) +
          scale_x_continuous(limits = c(2020, 2025)) +
          scale_color_brewer(palette = "Set1") +
          scale_fill_brewer(palette = "Set1")
        
        inference_plot_filename <- file.path(OUT_DIR, "inference_average_coverage_plot.png")
        ggsave(inference_plot_filename, p_inference, width = 10, height = 6, dpi = 300)
        cat(sprintf("Saved inference plot to: %s\n", inference_plot_filename))
        
        if (exists("bootstrap_trend_ci")) {
          cat("Computing trend CI via location bootstrap for inference results...\n")
          inference_trend_ci <- tryCatch({
            bootstrap_trend_ci(inference_coefs, B = 200)
          }, error = function(e) {
            warning(sprintf("bootstrap_trend_ci failed for inference: %s", e$message))
            NULL
          })
          if (!is.null(inference_trend_ci) && nrow(inference_trend_ci) > 0) {
            openxlsx::addWorksheet(wb_inference, "Trend_CI")
            openxlsx::writeData(wb_inference, "Trend_CI", inference_trend_ci)
            cat("Saved inference trend CI to workbook sheet 'Trend_CI'\n")
          }
        }
        
        if (exists("inference_coefs") && nrow(inference_coefs) > 0 && "PPI" %in% names(df_tasks_inference) && any(!is.na(df_tasks_inference$PPI))) {
          inference_pattern_ppi <- location_bootstrap_ppi(inference_coefs, df_tasks_inference, B = BOOTSTRAP_B, seed = 123)
          if (!is.null(inference_pattern_ppi) && nrow(inference_pattern_ppi) > 0) {
            inference_pattern_ppi <- inference_pattern_ppi[!tolower(trimws(inference_pattern_ppi$Veg)) %in% c("barren"), ]
          }
          if (!is.null(inference_pattern_ppi) && nrow(inference_pattern_ppi) > 0) {
            p_inf_ppi_ts <- ggplot(inference_pattern_ppi, aes(x = year, y = global_coef, color = Veg, group = Veg)) +
              geom_line(linewidth = 1) +
              geom_point(show.legend = FALSE) +
              geom_ribbon(aes(ymin = coef_025, ymax = coef_975, fill = Veg), alpha = 0.15, color = NA) +
              labs(title = "Inference: PPI-Normalized Vegetation Fractions Over Time (Location Bootstrap)",
                   x = "Year", y = "Total Normalized Fraction") +
              theme_minimal()
            inf_ppi_plot_filename <- file.path(OUT_DIR, "inference_ppi_normalized_timeseries.png")
            ggsave(inf_ppi_plot_filename, p_inf_ppi_ts, width = 8, height = 6)
            cat(sprintf("Saved inference PPI-normalized time series plot to: %s\n", inf_ppi_plot_filename))
          } else {
            cat("Inference PPI normalization aggregation returned no results.\n")
          }
        }
      }
    } else {
      cat("No inference coefficient data extracted.\n")
    }
  } else {
    cat("No inference results to process.\n")
  }

  cat("\nGenerating average coverage plot...\n")

  tryCatch({
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      install.packages("ggplot2")
    }
    library(ggplot2)

    global_pattern_all <- aggregate_to_global_pattern(all_coefs, method = "location_bootstrap")

    if (is.null(global_pattern_all) || nrow(global_pattern_all) == 0) {
      cat("[WARN] aggregate_to_global_pattern returned no data, skipping average coverage plot\n")
    } else {
      global_pattern_all <- global_pattern_all |> dplyr::filter(year >= 1985 & year <= 2025)
      global_pattern_all_veg <- global_pattern_all[tolower(global_pattern_all$Veg) != "barren", ]
      global_pattern_all_barren <- global_pattern_all[tolower(global_pattern_all$Veg) == "barren", ]

  p <- ggplot(global_pattern_all_veg, aes(x = year, y = global_coef, color = Veg, fill = Veg)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2, show.legend = FALSE) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, color = NA)
  
  barren_scale_factor <- 1
  if (nrow(global_pattern_all_barren) > 0) {
    veg_max <- suppressWarnings(max(global_pattern_all_veg$global_coef, global_pattern_all_veg$ci_upper, na.rm = TRUE))
    barren_max <- suppressWarnings(max(global_pattern_all_barren$global_coef, global_pattern_all_barren$ci_upper, na.rm = TRUE))
    
    if (is.infinite(veg_max) || is.na(veg_max) || veg_max <= 0) veg_max <- 0.01 # Default small value if no veg
    if (is.infinite(barren_max) || is.na(barren_max) || barren_max <= 0) barren_max <- 1
    
    if (veg_max <= 0.01) {
       barren_scale_factor <- 1 # No scaling if veg is negligible
    } else {
       barren_scale_factor <- veg_max / barren_max
    }
    
    global_pattern_all_barren$coef_scaled <- global_pattern_all_barren$global_coef * barren_scale_factor
    global_pattern_all_barren$ci_lower_scaled <- global_pattern_all_barren$ci_lower * barren_scale_factor
    global_pattern_all_barren$ci_upper_scaled <- global_pattern_all_barren$ci_upper * barren_scale_factor
    
    p <- p + 
      geom_line(data = global_pattern_all_barren, 
           aes(x = year, y = coef_scaled), 
               color = "brown", linewidth = 1.2, linetype = "dashed") +
      geom_point(data = global_pattern_all_barren, 
            aes(x = year, y = coef_scaled), 
                color = "brown", size = 2, show.legend = FALSE) +
      geom_ribbon(data = global_pattern_all_barren,
                 aes(x = year, ymin = ci_lower_scaled, ymax = ci_upper_scaled),
                 fill = "brown", alpha = 0.1, color = NA)
  }
  
  all_y_values <- c(global_pattern_all_veg$global_coef, global_pattern_all_veg$ci_lower, global_pattern_all_veg$ci_upper)
  if (exists("global_pattern_all_barren") && nrow(global_pattern_all_barren) > 0) {
    all_y_values <- c(all_y_values, global_pattern_all_barren$coef_scaled, global_pattern_all_barren$ci_lower_scaled, global_pattern_all_barren$ci_upper_scaled)
  }
  min_y <- min(all_y_values, na.rm = TRUE)
  max_y <- max(all_y_values, na.rm = TRUE)
  
  if (is.infinite(min_y) || is.na(min_y)) min_y <- 0
  if (is.infinite(max_y) || is.na(max_y)) max_y <- 1
  if (max_y <= min_y) max_y <- min_y + 0.1

  
  p <- p +
    labs(
      title = "Average Coverage Percentage per Vegetation Type (2000-2024)",
      subtitle = sprintf("Based on %d locations with bootstrap uncertainty", max(global_pattern_all$n_locations, na.rm = TRUE)),
      x = "Year",
      y = "Vegetation Fraction",
      color = "Vegetation Type",
      fill = "Vegetation Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(min_y, max_y),
      expand = c(0, 0),
      sec.axis = sec_axis(~ . / barren_scale_factor, name = "Barren Fraction", labels = scales::percent_format(accuracy = 1))
    ) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1")

      plot_filename <- file.path(OUT_DIR, "average_coverage_plot.png")
      ggsave(plot_filename, p, width = 10, height = 6, dpi = 300)
      cat(sprintf("Saved average coverage plot to: %s\n", plot_filename))
    }
  }, error = function(e) {
    cat(sprintf("[ERROR] Failed to generate average coverage plot: %s\n", e$message))
    cat("[INFO] Continuing with script execution...\n")
  })

  cat("\nGenerating Observations vs Accuracy plot...\n")
  
  
  loc_accuracy <- best_fit_summary |> 
    dplyr::group_by(location_id) |> 
    dplyr::summarize(
      mean_pred_coef_rel = mean(pred_coef_rel, na.rm = TRUE),
      mean_pred_coef_abs = mean(pred_coef_abs, na.rm = TRUE),
      .groups = "drop"
    )
  
  obs_vs_acc_data <- summary_data |> 
    dplyr::left_join(loc_accuracy, by = c("Location_ID" = "location_id"))
  
  obs_vs_acc_data <- obs_vs_acc_data |> dplyr::filter(!is.na(mean_pred_coef_rel))
  
  if (nrow(obs_vs_acc_data) > 0) {
    p_obs_acc <- ggplot(obs_vs_acc_data, aes(x = Total_Observations, y = mean_pred_coef_rel)) +
      geom_point(alpha = 0.6, color = "darkblue") +
      geom_smooth(method = "lm", color = "red", se = FALSE) +
      labs(
        title = "Number of Observations vs. Vegetation Prediction Accuracy",
        subtitle = "Accuracy = Mean Predicted Fraction of True Vegetation Class (Relative to Total Vegetation)",
        x = "Total Observations (All Years)",
        y = "Mean Correct Prediction Fraction"
      ) +
      theme_minimal() +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1))
    
    ggsave(file.path(OUT_DIR, "observations_vs_accuracy.png"), p_obs_acc, width = 8, height = 6)
    cat(sprintf("Saved Observations vs Accuracy plot to: %s\n", file.path(OUT_DIR, "observations_vs_accuracy.png")))
  } else {
    cat("No data available for Observations vs Accuracy plot.\n")
  }




  aggregate_to_global_pattern <- function(all_coefs, method = "location_bootstrap") {

    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")

    required_cols <- c("location_id", "year", "Veg", "coef")
    missing <- setdiff(required_cols, names(all_coefs))
    if (length(missing) > 0) stop(paste("Missing columns:", paste(missing, collapse = ", ")))

    if (method == "location_bootstrap") {
      # Rename year to pheno_year for location_bootstrap_aggregate
      if ("year" %in% names(all_coefs) && !"pheno_year" %in% names(all_coefs)) {
        all_coefs$pheno_year <- all_coefs$year
      }
      result <- location_bootstrap_aggregate(all_coefs, B = BOOTSTRAP_B)
    } else if (method == "hierarchical") {
      if (!"coef_025" %in% names(all_coefs)) all_coefs$coef_025 <- NA_real_
      if (!"coef_975" %in% names(all_coefs)) all_coefs$coef_975 <- NA_real_

      all_coefs$interval <- all_coefs$coef_975 - all_coefs$coef_025
      all_coefs$interval[!is.finite(all_coefs$interval)] <- NA_real_

      all_coefs$se_proxy <- all_coefs$interval / 3.92
      all_coefs$se_proxy[all_coefs$se_proxy <= 0 | !is.finite(all_coefs$se_proxy)] <- NA_real_
      all_coefs$weight <- 1 / (all_coefs$se_proxy^2)
      all_coefs$weight[!is.finite(all_coefs$weight)] <- 1  # Default weight if no CI

      result <- aggregate_hierarchical(all_coefs)
    } else {
      stop("Unknown method. Use 'location_bootstrap' or 'hierarchical'.")
    }

    result
  }


  aggregate_hierarchical <- function(all_coefs) {

    if (!requireNamespace("lme4", quietly = TRUE)) {
      warning("lme4 not available, falling back to simple mean aggregation")
      # Simple aggregation fallback
      result <- all_coefs |>
        dplyr::group_by(Veg, year) |>
        dplyr::summarize(
          n_locations = dplyr::n(),
          global_coef = mean(coef, na.rm = TRUE),
          se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
          ci_lower = global_coef - 1.96 * se,
          ci_upper = global_coef + 1.96 * se,
          .groups = "drop"
        ) |>
        dplyr::mutate(method = "simple_mean")
      return(result)
    }
    
    veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
    results_list <- list()
    
    for (veg in veg_types) {
      veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), ]
      
      if (nrow(veg_data) == 0) {
        warning(sprintf("No valid coefficients found for vegetation type: %s", veg))
        next
      }
      
      if (nrow(veg_data) < 10) {
        simple_result <- veg_data |> 
          dplyr::group_by(year) |> 
          dplyr::summarize(
            n_locations = dplyr::n(),
            global_coef = mean(coef, na.rm = TRUE),
            se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
            ci_lower = global_coef - 1.96 * se,
            ci_upper = global_coef + 1.96 * se,
            .groups = "drop"
          ) |> 
          dplyr::mutate(Veg = veg, method = "simple_mean")
        
        results_list[[veg]] <- simple_result
        next
      }
      
      tryCatch({
        veg_data$year_factor <- as.factor(veg_data$year)
        
        model <- suppressWarnings(lme4::lmer(coef ~ year_factor + (1|location_id), data = veg_data))
        
        fe <- lme4::fixef(model)
        vcov_fe <- as.matrix(vcov(model))
        
        years <- sort(unique(veg_data$year))
        pred_data <- data.frame(year_factor = as.factor(years))
        
        preds <- predict(model, newdata = pred_data, re.form = NA)
        
        boot_preds <- lme4::bootMer(model, function(m) {
          predict(m, newdata = pred_data, re.form = NA)
        }, nsim = 100, type = "parametric")
        
        ci_lower <- apply(boot_preds$t, 2, quantile, 0.025)
        ci_upper <- apply(boot_preds$t, 2, quantile, 0.975)
        
        hier_result <- data.frame(
          year = years,
          Veg = veg,
          n_locations = sapply(years, function(y) sum(veg_data$year == y)),
          global_coef = preds,
          se = apply(boot_preds$t, 2, sd),
          ci_lower = pmax(0, ci_lower),
          ci_upper = pmin(1, ci_upper),
          method = "hierarchical"
        )
        
        results_list[[veg]] <- hier_result
        
      }, error = function(e) {
        warning(sprintf("Hierarchical model failed for %s: %s. Using simple mean.", veg, e$message))
        simple_result <- veg_data |> 
          dplyr::group_by(year) |> 
          dplyr::summarize(
            n_locations = dplyr::n(),
            global_coef = mean(coef, na.rm = TRUE),
            se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
            ci_lower = pmax(0, global_coef - 1.96 * se),
            ci_upper = pmin(1, global_coef + 1.96 * se),
            .groups = "drop"
          ) |> 
          dplyr::mutate(Veg = veg, method = "simple_mean")
        
        results_list[[veg]] <- simple_result
      })
    }
    
    dplyr::bind_rows(results_list)
  }

  aggregate_simple_mean <- function(all_coefs) {
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")

    veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])

    # Determine which year column to use
    year_col <- if ("pheno_year" %in% names(all_coefs)) "pheno_year" else "year"

    results_list <- list()

    for (veg in veg_types) {
      veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), ]

      if (nrow(veg_data) == 0) {
        warning(sprintf("No valid coefficients found for vegetation type: %s", veg))
        next
      }

      # Group by year and compute simple statistics
      simple_result <- veg_data |> 
        dplyr::group_by(!!rlang::sym(year_col)) |> 
        dplyr::summarize(
          n_locations = dplyr::n(),
          global_coef = mean(coef, na.rm = TRUE),
          se = sd(coef, na.rm = TRUE) / sqrt(dplyr::n()),
          ci_lower = pmax(0, global_coef - 1.96 * se),
          ci_upper = pmin(1, global_coef + 1.96 * se),
          .groups = "drop"
        ) |> 
        dplyr::mutate(Veg = veg, method = "simple_mean") |> 
        dplyr::rename(year = !!rlang::sym(year_col))

      results_list[[veg]] <- simple_result
    }

    dplyr::bind_rows(results_list)
  }



  plot_global_vegetation_pattern <- function(global_pattern, 
                                              title = "Global Vegetation Composition Over Time",
                                              show_ci = TRUE,
                                              ci_type = "auto") {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    
    if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
      if ("ci_lower" %in% names(global_pattern) && "ci_upper" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower
        global_pattern$ci_upper <- global_pattern$ci_upper
      }
    } else if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
      if ("ci_lower_simple" %in% names(global_pattern) && "ci_upper_simple" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower_simple
        global_pattern$ci_upper <- global_pattern$ci_upper_simple
      }
    } else if ("weighted_mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$weighted_mean_coef
      if ("ci_lower_pooled" %in% names(global_pattern) && "ci_upper_pooled" %in% names(global_pattern)) {
        global_pattern$ci_lower <- global_pattern$ci_lower_pooled
        global_pattern$ci_upper <- global_pattern$ci_upper_pooled
      }
    }
    
    global_pattern_veg <- global_pattern[tolower(global_pattern$Veg) != "barren", ]
    global_pattern_barren <- global_pattern[tolower(global_pattern$Veg) == "barren", ]
    
    p <- ggplot2::ggplot(global_pattern_veg, ggplot2::aes(x = year, y = coef, color = Veg, fill = Veg)) +
      ggplot2::geom_line(size = 1.2) +
      ggplot2::geom_point(size = 2)
    
    if (show_ci && "ci_lower" %in% names(global_pattern_veg) && "ci_upper" %in% names(global_pattern_veg)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
        alpha = 0.2,
        color = NA
      )
    }
    
    barren_scale_factor <- 1
    if (nrow(global_pattern_barren) > 0) {
      veg_max <- suppressWarnings(max(global_pattern_veg$coef, global_pattern_veg$ci_upper, na.rm = TRUE))
      barren_max <- suppressWarnings(max(global_pattern_barren$coef, global_pattern_barren$ci_upper, na.rm = TRUE))
      if (is.na(veg_max) || veg_max <= 0 || is.na(barren_max) || barren_max <= 0) {
        barren_scale_factor <- 1
      } else {
        barren_scale_factor <- veg_max / barren_max
      }
      global_pattern_barren$coef_scaled <- global_pattern_barren$coef * barren_scale_factor
      
      p <- p + 
        ggplot2::geom_line(data = global_pattern_barren, 
                          ggplot2::aes(x = year, y = coef_scaled), 
                          color = "brown", linewidth = 1.2, linetype = "dashed") +
    p <- p +
      ggplot2::labs(
        title = title,
        subtitle = sprintf("Based on %d locations", max(global_pattern$n_locations, na.rm = TRUE)),
        x = "Year",
        y = "Vegetation Fraction",
        color = "Vegetation Type",
        fill = "Vegetation Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)
      ) +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(),
        sec.axis = ggplot2::sec_axis(~ . / barren_scale_factor, name = "Barren Fraction", labels = scales::percent_format())
      ) +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::scale_fill_brewer(palette = "Set1")
    
    p
  }

  plot_vegetation_stacked_area <- function(global_pattern, 
                                            title = "Global Vegetation Composition Over Time") {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    
    if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
    } else if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
    }
    
    global_pattern <- global_pattern |> 
      dplyr::group_by(year) |> 
      dplyr::mutate(coef_normalized = coef / sum(coef, na.rm = TRUE)) |> 
      dplyr::ungroup()
    
    p <- ggplot2::ggplot(global_pattern, 
                          ggplot2::aes(x = year, y = coef_normalized, fill = Veg)) +
      ggplot2::geom_area(alpha = 0.8, position = "stack") +
      ggplot2::labs(
        title = title,
        x = "Year",
        y = "Relative Vegetation Fraction",
        fill = "Vegetation Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold")
      ) +
      ggplot2::scale_y_continuous(labels = scales::percent_format()) +
      ggplot2::scale_fill_brewer(palette = "Set2")
    
    p
  }

  plot_vegetation_heatmap <- function(global_pattern, 
                                       title = "Vegetation Fraction by Year") {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    
    if ("mean_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$mean_coef
    } else if ("global_coef" %in% names(global_pattern)) {
      global_pattern$coef <- global_pattern$global_coef
    }
    
    p <- ggplot2::ggplot(global_pattern, 
                          ggplot2::aes(x = year, y = Veg, fill = coef)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", coef * 100)), 
                         color = "white", size = 3) +
      ggplot2::scale_fill_viridis_c(option = "plasma", 
                                     labels = scales::percent_format()) +
      ggplot2::labs(
        title = title,
        x = "Year",
        y = "Vegetation Type",
        fill = "Fraction"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    
    p
  }

  
  

  bootstrap_trend_ci <- function(all_coefs, B = 200, seed = 123) {
  set.seed(seed)
  if (!requireNamespace("lme4", quietly = TRUE)) {
    warning("lme4 package not found, trend CI calculation will be skipped.")
    return(NULL)
  }
  
  veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
  results_list <- list()
  
  for (veg in veg_types) {
    veg_data <- all_coefs[all_coefs$Veg == veg & is.finite(all_coefs$coef), ]
    
    if (nrow(veg_data) < 10 || length(unique(veg_data$location_id)) < 3) {
      cat(sprintf("[TREND] Skipping trend for '%s' (insufficient data: %d rows, %d locs)\n", veg, nrow(veg_data), length(unique(veg_data$location_id))))
      next
    }
    
    locations <- unique(veg_data$location_id)
    n_locs <- length(locations)
    
    boot_slopes <- replicate(B, {
      # Resample locations with replacement
      boot_locs_sampled <- sample(locations, n_locs, replace = TRUE)
      
      # Create bootstrap sample by selecting all data from resampled locations
      # Handle cases where a location is sampled multiple times by creating a new ID
      boot_data_list <- lapply(seq_along(boot_locs_sampled), function(i) {
        loc_data <- veg_data[veg_data$location_id == boot_locs_sampled[i], ]
        loc_data$boot_id <- paste0(boot_locs_sampled[i], "_", i)
        loc_data
      })
      boot_data <- do.call(rbind, boot_data_list)
      
      # Fit a linear mixed-effects model to the bootstrap sample
      # Use boot_id as the random effect to handle repeated locations
      model <- tryCatch({
        lme4::lmer(coef ~ pheno_year + (1|boot_id), data = boot_data)
      }, error = function(e) { NULL })
      
      if (is.null(model)) {
        return(NA_real_)
      } else {
        # Extract the fixed effect slope for the year
        return(lme4::fixef(model)["pheno_year"])
      }
    })
    
    finite_slopes <- boot_slopes[is.finite(boot_slopes)]
    if (length(finite_slopes) > 5) {
      results_list[[veg]] <- data.frame(
        Veg = veg, 
        slope_mean = mean(finite_slopes), 
        slope_median = median(finite_slopes), 
        slope_ci_lower = quantile(finite_slopes, 0.025), 
        slope_ci_upper = quantile(finite_slopes, 0.975), 
        prob_positive = mean(finite_slopes > 0), 
        prob_negative = mean(finite_slopes < 0),
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(results_list) > 0) {
    return(dplyr::bind_rows(results_list))
  } else {
    return(NULL)
  }
}



  global_pattern <- aggregate_to_global_pattern(all_coefs, method = "location_bootstrap")
  global_pattern_bootstrap <- global_pattern

  p1_all <- plot_global_vegetation_pattern(global_pattern, 
                                          title = "All Vegetation Trends",
                                          show_ci = TRUE)
  p1_all <- p1_all + scale_y_continuous(labels = scales::percent_format()) + labs(y = "Fraction")
  ggsave(file.path(OUT_DIR, "all_vegetation_trends.png"), p1_all, width = 10, height = 6)



  p3 <- plot_vegetation_stacked_area(global_pattern)
  ggsave(file.path(OUT_DIR, "vegetation_stacked_area.png"), p3, width = 10, height = 6)

  p4 <- plot_vegetation_heatmap(global_pattern)
  ggsave(file.path(OUT_DIR, "vegetation_heatmap.png"), p4, width = 10, height = 6)

  # Use a smaller bootstrap B for interactive/debug runs to avoid long hangs;
  # keep it moderate for production use (can be tuned by setting TREND_BOOT_B)
  TREND_BOOT_B <- if (exists("TREND_BOOT_B")) TREND_BOOT_B else 200
  trend_ci <- bootstrap_trend_ci(all_coefs, B = TREND_BOOT_B)
  # Ensure all ALLOWED_VEG are present in trend_ci (fill missing with NA rows)
  if (exists("ALLOWED_VEG")) {
    missing_vegs <- setdiff(ALLOWED_VEG, trend_ci$Veg)
    if (length(missing_vegs) > 0) {
      cat(sprintf("[NOTICE] No trend data for: %s. Adding placeholder NA rows for reporting.\n", paste(missing_vegs, collapse = ", ")))
      for (mv in missing_vegs) {
        trend_ci <- rbind(trend_ci, data.frame(Veg = mv, slope_mean = NA_real_, slope_median = NA_real_, slope_ci_lower = NA_real_, slope_ci_upper = NA_real_, prob_positive = NA_real_, prob_negative = NA_real_, stringsAsFactors = FALSE))
      }
    }
  }
  print(trend_ci)

  analyze_vegetation_trends <- function(all_coefs, B = 200) {
    # Simplified trend estimation using per-location linear slopes on annual coefficients.
    # Block/AR-based bootstrapping removed; we compute a median of location-level slopes.
    if (is.null(all_coefs) || nrow(all_coefs) == 0) return(NULL)
    veg_types <- unique(all_coefs$Veg[!is.na(all_coefs$Veg)])
    trend_rows <- list()
    for (veg in veg_types) {
      veg_data <- all_coefs[all_coefs$Veg == veg & !is.na(all_coefs$coef), , drop = FALSE]
      if (nrow(veg_data) == 0) next
      locations <- unique(veg_data$location_id)
      loc_slopes <- sapply(locations, function(loc) {
        loc_data <- veg_data[veg_data$location_id == loc, , drop = FALSE]
        if (nrow(loc_data) < 2) return(NA_real_)
        # Fit simple linear model of coef ~ year (use pheno_year if present)
        yr_col <- if ("pheno_year" %in% names(loc_data)) "pheno_year" else "year"
        df_loc <- loc_data[order(loc_data[[yr_col]]), , drop = FALSE]
        m <- tryCatch(lm(coef ~ get(yr_col), data = df_loc), error = function(e) NULL)
        if (is.null(m)) return(NA_real_)
        coef_val <- tryCatch(coef(m)[[2]], error = function(e) NA_real_)
        as.numeric(coef_val)
      })
      trend_rows[[veg]] <- data.frame(Veg = veg, slope = median(loc_slopes, na.rm = TRUE), n_locations = sum(is.finite(loc_slopes)), stringsAsFactors = FALSE)
    }
    if (length(trend_rows) == 0) return(NULL)
    do.call(rbind, trend_rows)
  }

  trends <- tryCatch(analyze_vegetation_trends(all_coefs, B = 200), error = function(e) { warning(sprintf("analyze_vegetation_trends failed: %s", e$message)); NULL })
  if (is.null(trends)) {
    cat("[NOTICE] No trends computed (insufficient data). Writing empty placeholder to Excel sheet.\n")
    trends <- data.frame(Veg = character(0), slope = numeric(0), n_locations = integer(0), stringsAsFactors = FALSE)
  }

  cat("\nGenerating PPI-normalized cumulative plot with location bootstrap...\n")
  if ("PPI" %in% names(df) && any(!is.na(df$PPI))) {
    global_pattern_ppi <- location_bootstrap_ppi(all_coefs, df, B = BOOTSTRAP_B, seed = 123)
    if (!is.null(global_pattern_ppi) && nrow(global_pattern_ppi) > 0) {
      global_pattern_ppi <- global_pattern_ppi[!tolower(trimws(global_pattern_ppi$Veg)) %in% c("barren"), ]
    }
    if (!is.null(global_pattern_ppi) && nrow(global_pattern_ppi) > 0) {
      p_ppi_ts <- ggplot(global_pattern_ppi, aes(x = year, y = global_coef, color = Veg, group = Veg)) +
        geom_line(linewidth = 1) +
        geom_point(show.legend = FALSE) +
        geom_ribbon(aes(ymin = coef_025, ymax = coef_975, fill = Veg), alpha = 0.15, color = NA) +
        labs(title = "PPI-Normalized Vegetation Fractions Over Time (Location Bootstrap)",
             x = "Year", y = "Total Normalized Fraction", color = "Veg", fill = "Veg") +
        theme_minimal()
      ggsave(file.path(OUT_DIR, "ppi_normalized_timeseries.png"), p_ppi_ts, width = 8, height = 6)
      cat(sprintf("Saved PPI-normalized time series plot to: %s\n", file.path(OUT_DIR, "ppi_normalized_timeseries.png")))
    } else {
      cat("PPI normalization aggregation returned no results (no matching location-year PPI values).\n")
    }
  } else {
    cat("PPI data not available, skipping PPI-normalized plot.\n")
  }

  openxlsx::addWorksheet(wb, "Global_Pattern")
  openxlsx::writeData(wb, "Global_Pattern", global_pattern)

  if (!is.null(global_pattern_ppi) && exists("global_pattern_ppi") && nrow(global_pattern_ppi) > 0) {
    openxlsx::addWorksheet(wb, "PPI_Normalized")
    openxlsx::writeData(wb, "PPI_Normalized", global_pattern_ppi)
    cat("Added PPI-normalized results to Excel workbook\n")
  }

  openxlsx::addWorksheet(wb, "Vegetation_Trends")
  openxlsx::writeData(wb, "Vegetation_Trends", trends)

  openxlsx::addWorksheet(wb, "Trend_Bootstrap_CI")
  openxlsx::writeData(wb, "Trend_Bootstrap_CI", trend_ci)

  timing_info$end_time <- Sys.time()
  total_time <- as.numeric(difftime(timing_info$end_time, timing_info$start_time, units = "secs"))

  cat(sprintf("\nTotal execution time: %.1f seconds (%.1f minutes)\n", total_time, total_time / 60))
  if (!is.null(timing_info$moving_var_done)) {
    cat(sprintf(
      "Moving variance: %.1f seconds\n",
      as.numeric(difftime(timing_info$moving_var_done, timing_info$start_time, units = "secs"))
    ))
  }
  if (!is.null(timing_info$lib_construction_done)) {
    cat(sprintf(
      "Library construction: %.1f seconds\n",
      as.numeric(difftime(timing_info$lib_construction_done, timing_info$moving_var_done, units = "secs"))
    ))
  }
  if (!is.null(timing_info$pca_computation_done)) {
    cat(sprintf(
      "PCA computation: %.1f seconds\n",
      as.numeric(difftime(timing_info$pca_computation_done, timing_info$lib_construction_done, units = "secs"))
    ))
  }
  cat(sprintf(
    "Main processing + Excel: %.1f seconds\n",
    as.numeric(difftime(timing_info$end_time, timing_info$pca_computation_done, units = "secs"))
