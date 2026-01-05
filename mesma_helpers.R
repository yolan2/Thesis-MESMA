
analyze_library_similarity <- function(mesma_lib, compressed_templates_accessor, grid_type = "full") {
  cat("\n=== INTER-CLASS VARIANT SIMILARITY ANALYSIS ===\n")

  # DIAGNOSTICS
  cat(sprintf("[DIAGNOSTIC] names(mesma_lib): %s\n", paste(names(mesma_lib), collapse=", ")))
  cat(sprintf("[DIAGNOSTIC] names(compressed_templates_accessor): %s\n", paste(names(compressed_templates_accessor), collapse=", ")))
  
  # Collect compressed 'full' templates in a safe manner
  all_ids <- character()
  all_vegs <- character()
  all_vecs <- list()
  all_nsamples <- numeric()

  diag_count <- 0
  for (veg in names(mesma_lib)) {
    if (is.null(mesma_lib[[veg]])) next
    # Try to access raw_lib_templates from global environment for fallback
    rt <- if (exists("raw_lib_templates", envir = globalenv()) && !is.null(get("raw_lib_templates", envir = globalenv())[[veg]])) get("raw_lib_templates", envir = globalenv())[[veg]] else NULL
    for (variant in filter_variants_by_min_samples(mesma_lib[[veg]], min_samples = MIN_ENDMEMBER_SAMPLES, veg = veg, raw_template = rt)) {
      vid <- if (!is.null(variant$variant_id)) variant$variant_id else if (!is.null(variant$id)) variant$id else NA_character_
      
      # If vid is still NA, try to construct one
      if (is.na(vid)) {
         # heuristic to generate a unique ID
         vid <- paste0(veg, "_", length(all_ids) + 1)
      }

      vec <- NULL
      
      # Robust lookup strategy
      # 1. Try exact match
      if (!is.null(compressed_templates_accessor[[veg]]) && !is.null(compressed_templates_accessor[[veg]][[vid]]) && !is.null(compressed_templates_accessor[[veg]][[vid]][["full"]])) {
         vec <- compressed_templates_accessor[[veg]][[vid]][["full"]]
      }
      
      # 2. Try case-insensitive veg match if failed
      if (is.null(vec)) {
         c_vegs <- names(compressed_templates_accessor)
         # Find veg key that matches case-insensitively
         v_match <- c_vegs[tolower(c_vegs) == tolower(veg)]
         if (length(v_match) > 0) {
            used_veg_key <- v_match[1]
            acc_veg <- compressed_templates_accessor[[used_veg_key]]
            
            # Try exact vid match in the resolved veg
            if (!is.null(acc_veg[[vid]]) && !is.null(acc_veg[[vid]][["full"]])) {
               vec <- acc_veg[[vid]][["full"]]
            } else {
               # 3. Try case-insensitive vid match
               c_vids <- names(acc_veg)
               vid_match <- c_vids[tolower(c_vids) == tolower(vid)]
               if (length(vid_match) > 0) {
                  vec <- acc_veg[[vid_match[1]]][["full"]]
               }
            }
         }
      }

      # 3. Fallback to using raw_mat from variant if available
      if (is.null(vec) && !is.null(variant$raw_mat)) {
         vec <- as.numeric(variant$raw_mat)
         vec[!is.finite(vec)] <- NA_real_
      }
      # 4. Fallback to using vec from variant if available
      if (is.null(vec) && !is.null(variant$vec)) {
         vec <- as.numeric(variant$vec)
         vec[!is.finite(vec)] <- NA_real_
      }

      # Diagnostic for the first few attempts
      if (diag_count < 5) {
         cat(sprintf("[DIAGNOSTIC] Checking veg='%s', vid='%s' -> Found? %s\n", veg, vid, !is.null(vec)))
         if (is.na(vid) || vid == "NA") {
             cat(sprintf("   [DEBUG] variant class: %s, names: %s\n", class(variant), paste(names(variant), collapse=", ")))
         }
         diag_count <- diag_count + 1
      }

      if (!is.null(vec) && length(vec) > 0 && !is.na(vid)) {
        all_vecs[[length(all_vecs) + 1]] <- as.numeric(vec)
        all_ids <- c(all_ids, vid)
        all_vegs <- c(all_vegs, veg)
        all_nsamples <- c(all_nsamples, if (!is.null(variant$n_samples)) variant$n_samples else 0)
      }
    }
  }

  cat(sprintf("[DIAGNOSTIC] Collected %d variant vectors for similarity analysis.\n", length(all_vecs)))

  # Print variant counts per vegetation type (without removing any)
  if (length(all_vecs) > 0) {
    veg_counts <- table(all_vegs)
    cat("[INFO] Variants per vegetation type:\n")
    print(veg_counts)
    
    if (min(veg_counts) < 2 && length(veg_counts) > 1) {
       cat("[INFO] Note: Some vegetation types have few variants. All variants will be retained for visualization.\n")
    }
  }

  if (length(all_vecs) < 2) {
    cat("[WARN] Not enough variant templates for similarity heatmap.\n")
    write.csv(data.frame(info = "insufficient_variants", n_variants = length(all_vecs)), file.path(OUT_DIR, "variant_similarity_matrix_info.csv"), row.names = FALSE)
    cat("==============================================\n\n")
    return(invisible(NULL))
  }

  # Normalize lengths by truncating to shortest vector if needed
  lens <- sapply(all_vecs, length)
  if (length(unique(lens)) > 1) {
    min_len <- min(lens)
    cat(sprintf("[WARN] Variant vectors have differing lengths; truncating to min length %d\n", min_len))
    all_vecs <- lapply(all_vecs, function(v) v[1:min_len])
  }

  n_v <- length(all_vecs)
  
  # Vectorized Similarity Calculation
  vec_mat <- do.call(rbind, all_vecs)
  vec_mat[!is.finite(vec_mat)] <- 0
  row_norms <- sqrt(rowSums(vec_mat^2))
  row_norms[row_norms == 0] <- 1
  vec_mat_norm <- vec_mat / row_norms
  sim_mat <- tcrossprod(vec_mat_norm)
  sim_mat[sim_mat > 1.01] <- 1.01
  sim_mat[sim_mat < -0.01] <- -0.01
  
  rownames(sim_mat) <- all_ids
  colnames(sim_mat) <- all_ids

  sim_df <- as.data.frame(as.table(sim_mat))
  colnames(sim_df) <- c("Var1", "Var2", "Similarity")
  
  # AGENT FIX: Ensure ordered_ids is defined
  ord_idx <- order(all_vegs, all_ids)
  ordered_ids <- all_ids[ord_idx]
  ordered_vegs <- all_vegs[ord_idx]
  
  # AGENT CHANGE: Add sample counts to weight tile thickness and allocate contiguous positions so there are no gaps between variant bars
  count_map <- setNames(all_nsamples, all_ids)

  # Determine ordered unique ids and counts
  counts_ordered <- as.numeric(count_map[ordered_ids])
  counts_ordered[is.na(counts_ordered)] <- 0
  total_counts <- sum(counts_ordered)
  if (!is.finite(total_counts) || total_counts <= 0) {
    n_ord <- length(ordered_ids)
    widths_prop <- rep(1 / n_ord, n_ord)
  } else {
    widths_prop <- counts_ordered / total_counts
    # Avoid zero widths (give tiny epsilon to zero-count variants so they remain visible)
    eps <- 1e-6
    widths_prop[widths_prop == 0] <- eps
    widths_prop <- widths_prop / sum(widths_prop)
  }

  # Compute continuous positions in [0,1]
  cumw <- cumsum(widths_prop)
  lefts <- c(0, cumw[-length(cumw)])
  rights <- cumw
  names(lefts) <- ordered_ids
  names(rights) <- ordered_ids

  # For Y axis we use reversed ordering (top-to-bottom)
  lefts_rev <- rev(lefts)
  rights_rev <- rev(rights)
  names(lefts_rev) <- rev(ordered_ids)
  names(rights_rev) <- rev(ordered_ids)

  # Map rectangle coordinates into sim_df
  sim_df$xmin <- lefts[as.character(sim_df$Var1)]
  sim_df$xmax <- rights[as.character(sim_df$Var1)]
  sim_df$ymin <- lefts_rev[as.character(sim_df$Var2)]
  sim_df$ymax <- rights_rev[as.character(sim_df$Var2)]

  # Compute centers for axis ticks/labels
  x_centers <- (lefts + rights) / 2
  y_centers <- (lefts_rev + rights_rev) / 2

  # Boundaries at vegetation type changes (in cumulative units)
  if (length(ordered_vegs) > 1) {
    change_idx <- which(ordered_vegs[-1] != ordered_vegs[-length(ordered_vegs)])
    vlines_x <- cumsum(widths_prop)[change_idx]
    hlines_y <- cumsum(widths_prop)[change_idx]
  } else {
    vlines_x <- numeric(0)
    hlines_y <- numeric(0)
  }

  # Build heatmap using fixed rectangles (no gaps between adjacent variants)
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  library(ggplot2)
  p_heat <- ggplot(sim_df) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Similarity), color = NA) +
    # Add thick black lines at vegetation boundaries
    geom_vline(xintercept = vlines_x, color = "black", size = 1.2) +
    geom_hline(yintercept = hlines_y, color = "black", size = 1.2) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.9, limits = c(-0.01, 1.01), na.value = "white") +
    scale_x_continuous(expand = c(0,0), breaks = x_centers, labels = ordered_ids) +
    scale_y_continuous(expand = c(0,0), breaks = y_centers, labels = rev(ordered_ids)) +
    coord_equal() +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6), axis.text.y = element_text(size = 6)) +
    labs(title = "Pairwise Cosine Similarity of All Variants", x = NULL, y = NULL)

  tryCatch({
    ggsave(file.path(OUT_DIR, "variant_similarity_heatmap.png"), p_heat, width = 12, height = 10)
    write.csv(sim_mat, file.path(OUT_DIR, "variant_similarity_matrix.csv"))
    cat(sprintf("Saved similarity heatmap to: %s\n", file.path(OUT_DIR, "variant_similarity_heatmap.png")))
  }, error = function(e) {
    cat(sprintf("[ERROR] Failed to generate similarity heatmap: %s\n", e$message))
  })

  # == AGENT FIX: Print #loc-years per variant ==
  cat("\n=== VARIANT SAMPLE SIZES (Loc-Years) ===\n")
  sample_stats <- data.frame(
    Vegetation = all_vegs,
    VariantID = all_ids,
    LocYears = all_nsamples
  )
  # Sort by vegetation then ID for display
  sample_stats <- sample_stats[order(sample_stats$Vegetation, sample_stats$VariantID), ]
  print(sample_stats, row.names = FALSE)
  write.csv(sample_stats, file.path(OUT_DIR, "variant_sample_sizes.csv"), row.names = FALSE)
  cat("=========================================\n\n")

  cat("==============================================\n\n")
  invisible(TRUE)
}


# Small test helper to exercise analyze_library_similarity without running the whole pipeline
# This can be invoked with: test_analyze_library_similarity()

test_analyze_library_similarity <- function() {
  cat("Running test_analyze_library_similarity()...\n")
  # Build minimal synthetic mesma_lib / compressed_templates_accessor
  pequena_lib <- list(
    vegA = list(list(variant_id = "vegA_v1", n_samples = 12)),
    vegB = list(list(variant_id = "vegB_v1", n_samples = 15))
  )
  compressed <- list(
    vegA = list(vegA_v1 = list(full = rep(1.0, 24))),
    vegB = list(vegB_v1 = list(full = rep(0.9, 24)))
  )

  # Use a temporary OUT_DIR to avoid clobbering real output
  tmp_out <- file.path(tempdir(), "test_variant_similarity")
  dir.create(tmp_out, showWarnings = FALSE, recursive = TRUE)
  old_out <- if (exists("OUT_DIR", inherits = FALSE)) get("OUT_DIR") else NULL
  assign("OUT_DIR", tmp_out, envir = .GlobalEnv)

  res <- tryCatch({
    analyze_library_similarity(pequena_lib, compressed)
    list(success = TRUE, dir = tmp_out, files = list.files(tmp_out))
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })

  if (!is.null(old_out)) assign("OUT_DIR", old_out, envir = .GlobalEnv) else rm("OUT_DIR", envir = .GlobalEnv)

  cat(sprintf("Test result: %s\n", if (isTRUE(res$success)) "OK" else paste0("FAILED: ", res$error)))
  if (!is.null(res$files)) cat(sprintf("Files created: %s\n", paste(res$files, collapse = ", ")))
  invisible(TRUE)
}


# Helper: ensure the variant similarity heatmap is generated once and early
ensure_global_dvi_soil_baseline <- function() {
  if (exists("GLOBAL_TRAINING_DVI_SOIL", envir = globalenv()) && !is.na(get("GLOBAL_TRAINING_DVI_SOIL", envir = globalenv())) && is.finite(get("GLOBAL_TRAINING_DVI_SOIL", envir = globalenv()))) {
    cat(sprintf("[PPI] GLOBAL_TRAINING_DVI_SOIL already set: %.6f\n", get("GLOBAL_TRAINING_DVI_SOIL", envir = globalenv())))
    return(invisible(get("GLOBAL_TRAINING_DVI_SOIL", envir = globalenv())))
  }

  cat("[PPI] GLOBAL_TRAINING_DVI_SOIL not found or invalid. Attempting to restore or recalculate...\n")
  dvi_soil_found <- FALSE
  dvi_soil_val <- NA_real_

  # 1. Attempt to restore from cache metadata or raw templates (existing logic)
  possible_cache_dirs <- unique(c(if (exists("OUT_DIR")) file.path(OUT_DIR, "mesma_cache") else NULL,
                                  "phenology_results/veg_mixture_fit/mesma_cache",
                                  file.path("mesma_cache")))
  for (cd in possible_cache_dirs) {
    if (!is.null(cd)) {
      # Try from training_metadata.rds
      if (!dvi_soil_found && file.exists(file.path(cd, "training_metadata.rds"))) {
        tryCatch({
          tm <- readRDS(file.path(cd, "training_metadata.rds"))
          if (!is.null(tm$dvi_soil) && is.finite(tm$dvi_soil)) {
            dvi_soil_val <- tm$dvi_soil
            assign("GLOBAL_TRAINING_DVI_SOIL", dvi_soil_val, envir = globalenv())
            cat(sprintf("[PPI] Restored GLOBAL_TRAINING_DVI_SOIL=%.6f from cache metadata\n", dvi_soil_val))
            dvi_soil_found <- TRUE
          }
        }, error = function(e) { cat(sprintf("[PPI ERROR] Failed to read dvi_soil from training_metadata.rds: %s\n", e$message)) })
      }
      
      # Try from raw_templates.rds if not yet found
      if (!dvi_soil_found && file.exists(file.path(cd, "raw_templates.rds"))) {
        tryCatch({
          raw_tpl <- readRDS(file.path(cd, "raw_templates.rds"))
          if (!is.null(raw_tpl[["barren"]])) {
            barren_raw <- raw_tpl[["barren"]]
            dvi_vals <- NULL
            if (is.data.frame(barren_raw) || is.matrix(barren_raw)) {
              if ("DVI" %in% colnames(barren_raw)) {
                dvi_vals <- barren_raw[, "DVI"]
              } else if ("nir" %in% colnames(barren_raw) && "red" %in% colnames(barren_raw)) {
                dvi_vals <- barren_raw[, "nir"] - barren_raw[, "red"]
              }
            }
            if (!is.null(dvi_vals) && length(dvi_vals) > 0 && any(is.finite(dvi_vals))) {
              dvi_soil_val <- median(dvi_vals, na.rm=TRUE)
              assign("GLOBAL_TRAINING_DVI_SOIL", dvi_soil_val, envir = globalenv())
              cat(sprintf("[PPI] Computed fallback GLOBAL_TRAINING_DVI_SOIL=%.6f from cached raw 'barren' templates\n", dvi_soil_val))
              dvi_soil_found <- TRUE
            }
          }
        }, error = function(e) { cat(sprintf("[PPI ERROR] Failed to compute dvi_soil from raw_templates.rds: %s\n", e$message)) })
      }
      if (dvi_soil_found) break # Stop searching if found
    }
  }

  # 2. If still not found, try to recalculate from df_train if available
  if (!dvi_soil_found && exists("df_train") && !is.null(df_train) && nrow(df_train) > 0) {
    cat("[PPI] No DVI soil baseline found in cache. Attempting to recalculate from df_train...\n")
    df_train_copy <- df_train # Work on a copy to avoid side effects
    
    # Ensure DVI and Veg columns for calculation
    if (!"DVI" %in% names(df_train_copy) && all(c("nir", "red") %in% names(df_train_copy))) {
      df_train_copy$DVI <- df_train_copy$nir - df_train_copy$red
    }
    
    if ("Veg" %in% names(df_train_copy) && "DVI" %in% names(df_train_copy)) {
      barren_dvi <- df_train_copy$DVI[is.finite(df_train_copy$DVI) & tolower(df_train_copy$Veg) == 'barren']
      if (length(barren_dvi) > 0) {
        dvi_soil_val <- as.numeric(median(barren_dvi, na.rm = TRUE))
        assign("GLOBAL_TRAINING_DVI_SOIL", dvi_soil_val, envir = globalenv())
        cat(sprintf("[PPI] Recalculated GLOBAL_TRAINING_DVI_SOIL=%.6f from df_train 'barren' observations\n", dvi_soil_val))
        dvi_soil_found <- TRUE
      } else {
        cat("[PPI WARNING] df_train exists but contains no 'barren' observations for recalculation.\n")
      }
    } else {
      cat("[PPI WARNING] df_train exists but missing 'Veg' or 'DVI' columns for recalculation.\n")
    }
  } else if (!dvi_soil_found) {
    cat("[PPI WARNING] df_train not found in memory, skipping recalculation from training data.\n")
  }

  if (!dvi_soil_found) {
    stop("[PPI ERROR] GLOBAL_TRAINING_DVI_SOIL could not be established after all attempts (cache or recalculation from df_train). Cannot proceed with PPI calculation.")
  }
  return(invisible(dvi_soil_val))
}

ensure_library_and_templates <- function(force = FALSE) {
  # Attempt to construct or expose `mesma_lib` and `compressed_templates_accessor`
  # so visualizations can run reliably. If these objects already exist, this function is a no-op.
  # If both objects already exist in globalenv, nothing to do
  if (!isTRUE(force) && exists("mesma_lib", envir = globalenv()) && (exists("compressed_templates_accessor", envir = globalenv()) || exists(".COMPRESSED_TEMPLATES_ACCESSOR", envir = globalenv()))) {
    return(invisible(TRUE))
  }

  cat("[INFO] Ensuring MESMA library and compressed templates are available for visualization...\n")

  # Attempt to use local variables if available in parent frame
  ms <- NULL; ct <- NULL
  if (exists("mesma_lib", envir = globalenv())) ms <- get("mesma_lib", envir = globalenv())
  if (exists("compressed_templates_accessor", envir = globalenv())) ct <- get("compressed_templates_accessor", envir = globalenv())
  if (is.null(ct) && exists(".COMPRESSED_TEMPLATES_ACCESSOR", envir = globalenv())) ct <- get(".COMPRESSED_TEMPLATES_ACCESSOR", envir = globalenv())
  if (is.null(ms) && exists("mesma_lib", envir = parent.frame())) ms <- get("mesma_lib", envir = parent.frame())
  if (is.null(ct) && exists("compressed_templates_accessor", envir = parent.frame())) ct <- get("compressed_templates_accessor", envir = parent.frame())
  if (is.null(ct) && exists(".COMPRESSED_TEMPLATES_ACCESSOR", envir = parent.frame())) ct <- get(".COMPRESSED_TEMPLATES_ACCESSOR", envir = parent.frame())

  # If still missing, attempt to build them from available inputs


  # If still missing, try to load from saved cache files
  if (is.null(ms)) {
    possible_cache_dirs <- unique(c(if (exists("OUT_DIR")) file.path(OUT_DIR, "mesma_cache") else NULL,
                                    "phenology_results/veg_mixture_fit/mesma_cache",
                                    file.path("mesma_cache")))
    cat(sprintf("[DEBUG] Checking possible cache dirs for 'mesma_library.rds': %s\n", paste(possible_cache_dirs, collapse = ", ")))
    for (cd in possible_cache_dirs) {
      if (is.null(cd)) next
      lib_path <- file.path(cd, "mesma_library.rds")
      cat(sprintf("[DEBUG] Checking %s -> exists=%s\n", normalizePath(cd, mustWork = FALSE), file.exists(lib_path)))
      if (file.exists(lib_path)) {
        tryCatch({
          ms_load <- readRDS(lib_path)
          if ("mesma_lib" %in% names(ms_load) && is.list(ms_load$mesma_lib)) {
             ms <- ms_load$mesma_lib
             cat("[INFO] Extracted 'mesma_lib' from cached wrapper object.\n")
          } else {
             ms <- ms_load
          }
          assign("mesma_lib", ms, envir = globalenv())
          cat(sprintf("[INFO] Loaded 'mesma_lib' from cache: %s\n", lib_path))
          break
        }, error = function(e) {
          cat(sprintf("[ERROR] Failed to read '%s': %s\n", lib_path, e$message))
        })
      }
    }
  }

  if (!is.null(ms) && is.null(ct)) {
    # Try to load compressed templates from cache dirs
    possible_cache_dirs <- unique(c(if (exists("OUT_DIR")) file.path(OUT_DIR, "mesma_cache") else NULL,
                                    "phenology_results/veg_mixture_fit/mesma_cache",
                                    file.path("mesma_cache")))
    cat(sprintf("[DEBUG] Checking possible cache dirs for 'compressed_templates.rds': %s\n", paste(possible_cache_dirs, collapse = ", ")))
    for (cd in possible_cache_dirs) {
      if (is.null(cd)) next
      tpl_path <- file.path(cd, "compressed_templates.rds")
      cat(sprintf("[DEBUG] Checking %s -> exists=%s\n", normalizePath(cd, mustWork = FALSE), file.exists(tpl_path)))
      if (file.exists(tpl_path)) {
        tryCatch({
          ct2 <- readRDS(tpl_path)
          
          # FIX: Check for flat keys and restructure if needed to ensure nested [[veg]][[vid]][[grid_type]] format
          is_flat <- FALSE
          if (length(ct2) > 0 && is.list(ct2[[1]]) && length(ct2[[1]]) > 0) {
             first_key <- names(ct2[[1]])[1]
             if (!is.null(first_key) && grepl("|", first_key, fixed=TRUE)) is_flat <- TRUE
          }
          
          if (is_flat) {
             cat("[INFO] Detected flat-key cache structure. Restructuring to nested format...\n")
             ct_nested <- list()
             for (veg in names(ct2)) {
                ct_nested[[veg]] <- list()
                for (key in names(ct2[[veg]])) {
                   parts <- strsplit(key, "|", fixed=TRUE)[[1]]
                   # key constructed as: paste(veg, variant$variant_id, grid_type, sep = "|")
                   # So parts[1]=veg, parts[last]=grid_type, middle=vid
                   if (length(parts) >= 3) {
                      gtype <- parts[length(parts)]
                      vid <- paste(parts[2:(length(parts)-1)], collapse="|")
                      
                      if (is.null(ct_nested[[veg]][[vid]])) ct_nested[[veg]][[vid]] <- list()
                      ct_nested[[veg]][[vid]][[gtype]] <- ct2[[veg]][[key]]
                   }
                }
             }
             ct2 <- ct_nested
          }
          
          assign("compressed_templates_accessor", ct2, envir = globalenv())
          assign(".COMPRESSED_TEMPLATES_ACCESSOR", ct2, envir = globalenv())
          ct <- ct2
          cat(sprintf("[INFO] Loaded 'compressed_templates_accessor' from cache: %s\n", tpl_path))
          break
        }, error = function(e) {
          cat(sprintf("[ERROR] Failed to read '%s': %s\n", tpl_path, e$message))
        })
      }
    }
  }

  if (is.null(ms) || is.null(ct)) {
    cat("[WARN] After attempts, 'mesma_lib' or 'compressed_templates_accessor' still not available. Some visualizations may be skipped.\n")
    return(invisible(FALSE))
  }

  invisible(TRUE)
}

ensure_variant_similarity_heatmap <- function(force = FALSE) {
  if (exists("VARIANT_SIMILARITY_HEATMAP_DONE", envir = globalenv()) && isTRUE(get("VARIANT_SIMILARITY_HEATMAP_DONE", envir = globalenv())) && !isTRUE(force)) {
    return(invisible(TRUE))
  }

  # Ensure that core objects are present (try to create them if missing)
  ensure_library_and_templates(force = force)

  ms <- NULL; ct <- NULL
  if (exists("mesma_lib", envir = globalenv())) ms <- get("mesma_lib", envir = globalenv())
  if (exists("compressed_templates_accessor", envir = globalenv())) ct <- get("compressed_templates_accessor", envir = globalenv())
  if (is.null(ct) && exists(".COMPRESSED_TEMPLATES_ACCESSOR", envir = globalenv())) ct <- get(".COMPRESSED_TEMPLATES_ACCESSOR", envir = globalenv())

  if (is.null(ms) || is.null(ct)) {
    cat("[WARN] Cannot generate variant similarity heatmap: 'mesma_lib' or 'compressed_templates_accessor' not available yet\n")
    return(invisible(FALSE))
  }

  tryCatch({
    analyze_library_similarity(ms, ct, grid_type = "full")
    assign("VARIANT_SIMILARITY_HEATMAP_DONE", TRUE, envir = globalenv())
    cat("[INFO] Variant similarity heatmap generated (ensure_variant_similarity_heatmap).\n")
  }, error = function(e) {
    cat(sprintf("[ERROR] Failed to generate variant similarity heatmap: %s\n", e$message))
  })
  invisible(TRUE)
}
