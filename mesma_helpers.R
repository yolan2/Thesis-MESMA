
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
      if (diag_count < 5 && exists("TESTING_MODE") && isTRUE(TESTING_MODE)) {
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

  # Apply barren similarity filter to vegetation variants (exclude barren-like vegetation)
  barren_threshold_exists <- exists("BARREN_SIM_THRESHOLD", envir = globalenv())
  cat(sprintf("[DEBUG] BARREN_SIM_THRESHOLD exists: %s\n", barren_threshold_exists))
  
  if (length(all_vecs) > 0 && barren_threshold_exists) {
    barren_threshold <- get("BARREN_SIM_THRESHOLD", envir = globalenv())
    cat(sprintf("[DEBUG] BARREN_SIM_THRESHOLD value: %.2f\n", barren_threshold))
    
    barren_indices <- which(tolower(all_vegs) == "barren")
    cat(sprintf("[DEBUG] Found %d barren variants at indices: %s\n", 
                length(barren_indices), paste(barren_indices, collapse=", ")))
    
    if (length(barren_indices) > 0) {
      # Normalize all vectors
      vec_mat_all <- do.call(rbind, all_vecs)
      vec_mat_all[!is.finite(vec_mat_all)] <- 0
      row_norms_all <- sqrt(rowSums(vec_mat_all^2))
      row_norms_all[row_norms_all == 0] <- 1
      vec_mat_norm_all <- vec_mat_all / row_norms_all
      
      # Get barren references
      barren_refs <- vec_mat_norm_all[barren_indices, , drop = FALSE]
      
      # Calculate cosine similarity for non-barren vegetation variants
      veg_indices <- which(tolower(all_vegs) != "barren")
      cat(sprintf("[DEBUG] Found %d vegetation variants to check\n", length(veg_indices)))
      
      if (length(veg_indices) > 0) {
        veg_vecs_norm <- vec_mat_norm_all[veg_indices, , drop = FALSE]
        sim_to_barren <- veg_vecs_norm %*% t(barren_refs)
        max_sim_to_barren <- apply(sim_to_barren, 1, function(x) suppressWarnings(max(x, na.rm = TRUE)))
        
        cat(sprintf("[DEBUG] Max similarities to barren: min=%.3f, max=%.3f, mean=%.3f\n",
                    min(max_sim_to_barren, na.rm=TRUE), 
                    max(max_sim_to_barren, na.rm=TRUE),
                    mean(max_sim_to_barren, na.rm=TRUE)))
        
        # Filter out vegetation variants too similar to barren
        keep_veg_mask <- max_sim_to_barren <= barren_threshold
        filtered_veg_indices <- veg_indices[keep_veg_mask]
        n_removed <- sum(!keep_veg_mask)
        
        if (n_removed > 0) {
          cat(sprintf("[BARREN FILTER] Removed %d vegetation variant(s) with cosine similarity > %.2f to barren\n", 
                      n_removed, barren_threshold))
          # Show which variants were removed
          removed_ids <- all_ids[veg_indices[!keep_veg_mask]]
          removed_sims <- max_sim_to_barren[!keep_veg_mask]
          for (j in seq_along(removed_ids)) {
            cat(sprintf("  - %s (similarity: %.3f)\n", removed_ids[j], removed_sims[j]))
          }
        } else {
          cat(sprintf("[BARREN FILTER] No vegetation variants exceeded similarity threshold of %.2f\n", barren_threshold))
        }
        
        # Keep all barren variants + filtered vegetation variants
        keep_indices <- c(barren_indices, filtered_veg_indices)
        all_vecs <- all_vecs[keep_indices]
        all_ids <- all_ids[keep_indices]
        all_vegs <- all_vegs[keep_indices]
        all_nsamples <- all_nsamples[keep_indices]
        
        cat(sprintf("[BARREN FILTER] Retained %d variants for heatmap after barren filtering\n", length(all_vecs)))
      }
    }
  } else if (length(all_vecs) > 0 && !barren_threshold_exists) {
    cat("[WARNING] BARREN_SIM_THRESHOLD not found in global environment - skipping barren filtering\n")
  }

  # Print variant counts per vegetation type (after barren filtering)
  if (length(all_vecs) > 0) {
    veg_counts <- table(all_vegs)
    cat("[INFO] Variants per vegetation type (after filtering):\n")
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

  # === STAGE 1: WOODY vs HERBS SIMILARITY HEATMAP (with Stage 1 LDA weights) ===
  # Generate a heatmap showing woody (populus+tamarix combined) vs herbs differentiation
  woody_classes <- if (exists("HIERARCHICAL_WOODY_CLASSES", envir = globalenv())) {
    get("HIERARCHICAL_WOODY_CLASSES", envir = globalenv())
  } else {
    c("populus", "tamarix")
  }

  # Get Stage 1 LDA weights if available
  stage1_weights <- NULL
  if (exists("HIERARCHICAL_LDA_PARAMS", envir = globalenv())) {
    hlp <- get("HIERARCHICAL_LDA_PARAMS", envir = globalenv())
    if (!is.null(hlp$stage1) && !is.null(hlp$stage1$weights)) {
      stage1_weights <- hlp$stage1$weights
      cat("[STAGE 1 HEATMAP] Using Stage 1 LDA weights for woody vs herbs similarity\n")
    }
  }

  # Identify woody and herbs indices
  woody_indices <- which(tolower(all_vegs) %in% tolower(woody_classes))
  herbs_indices <- which(tolower(all_vegs) == "herbs")
  stage1_indices <- c(woody_indices, herbs_indices)

  if (length(woody_indices) >= 1 && length(herbs_indices) >= 1) {
    cat("\n=== STAGE 1: WOODY vs HERBS VARIANT SIMILARITY ===\n")

    s1_vecs <- all_vecs[stage1_indices]
    s1_ids <- all_ids[stage1_indices]
    s1_vegs <- all_vegs[stage1_indices]
    s1_nsamples <- all_nsamples[stage1_indices]

    # Create stage1 class labels (woody vs herbs)
    s1_stage1_class <- ifelse(tolower(s1_vegs) %in% tolower(woody_classes), "woody", "herbs")

    # Build matrix and apply Stage 1 weights if available
    s1_vec_mat <- do.call(rbind, s1_vecs)
    s1_vec_mat[!is.finite(s1_vec_mat)] <- 0

    if (!is.null(stage1_weights) && length(stage1_weights) == ncol(s1_vec_mat)) {
      # Apply Stage 1 LDA weights
      s1_vec_mat_weighted <- sweep(s1_vec_mat, 2, sqrt(stage1_weights), "*")
      cat(sprintf("  Applied Stage 1 weights to %d features\n", length(stage1_weights)))
    } else {
      s1_vec_mat_weighted <- s1_vec_mat
      if (!is.null(stage1_weights)) {
        cat(sprintf("[WARNING] Stage 1 weights length (%d) != features (%d), using unweighted\n",
                    length(stage1_weights), ncol(s1_vec_mat)))
      }
    }

    # Normalize and compute similarity
    s1_row_norms <- sqrt(rowSums(s1_vec_mat_weighted^2))
    s1_row_norms[s1_row_norms == 0] <- 1
    s1_vec_mat_norm <- s1_vec_mat_weighted / s1_row_norms
    s1_sim_mat <- tcrossprod(s1_vec_mat_norm)
    s1_sim_mat[s1_sim_mat > 1.01] <- 1.01
    s1_sim_mat[s1_sim_mat < -0.01] <- -0.01

    rownames(s1_sim_mat) <- s1_ids
    colnames(s1_sim_mat) <- s1_ids

    s1_sim_df <- as.data.frame(as.table(s1_sim_mat))
    colnames(s1_sim_df) <- c("Var1", "Var2", "Similarity")

    # Order by stage1 class (herbs first, then woody)
    s1_ord_idx <- order(s1_stage1_class, s1_ids)
    s1_ordered_ids <- s1_ids[s1_ord_idx]
    s1_ordered_classes <- s1_stage1_class[s1_ord_idx]

    # Use sample counts for weighting
    s1_count_map <- setNames(s1_nsamples, s1_ids)
    s1_counts_ordered <- as.numeric(s1_count_map[s1_ordered_ids])
    s1_counts_ordered[is.na(s1_counts_ordered)] <- 0
    s1_total_counts <- sum(s1_counts_ordered)

    if (!is.finite(s1_total_counts) || s1_total_counts <= 0) {
      s1_widths_prop <- rep(1 / length(s1_ordered_ids), length(s1_ordered_ids))
    } else {
      s1_widths_prop <- s1_counts_ordered / s1_total_counts
      eps <- 1e-6
      s1_widths_prop[s1_widths_prop == 0] <- eps
      s1_widths_prop <- s1_widths_prop / sum(s1_widths_prop)
    }

    s1_cumw <- cumsum(s1_widths_prop)
    s1_lefts <- c(0, s1_cumw[-length(s1_cumw)])
    s1_rights <- s1_cumw
    names(s1_lefts) <- s1_ordered_ids
    names(s1_rights) <- s1_ordered_ids

    s1_lefts_rev <- rev(s1_lefts)
    s1_rights_rev <- rev(s1_rights)
    names(s1_lefts_rev) <- rev(s1_ordered_ids)
    names(s1_rights_rev) <- rev(s1_ordered_ids)

    s1_sim_df$xmin <- s1_lefts[as.character(s1_sim_df$Var1)]
    s1_sim_df$xmax <- s1_rights[as.character(s1_sim_df$Var1)]
    s1_sim_df$ymin <- s1_lefts_rev[as.character(s1_sim_df$Var2)]
    s1_sim_df$ymax <- s1_rights_rev[as.character(s1_sim_df$Var2)]

    s1_x_centers <- (s1_lefts + s1_rights) / 2
    s1_y_centers <- (s1_lefts_rev + s1_rights_rev) / 2

    # Boundaries at stage1 class changes
    if (length(s1_ordered_classes) > 1) {
      s1_change_idx <- which(s1_ordered_classes[-1] != s1_ordered_classes[-length(s1_ordered_classes)])
      s1_vlines_x <- cumsum(s1_widths_prop)[s1_change_idx]
      s1_hlines_y <- cumsum(s1_widths_prop)[s1_change_idx]
    } else {
      s1_vlines_x <- numeric(0)
      s1_hlines_y <- numeric(0)
    }

    # Build Stage 1 heatmap
    p_heat_s1 <- ggplot(s1_sim_df) +
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Similarity), color = NA) +
      geom_vline(xintercept = s1_vlines_x, color = "black", size = 1.2) +
      geom_hline(yintercept = s1_hlines_y, color = "black", size = 1.2) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.8, limits = c(-0.01, 1.01), na.value = "white") +
      scale_x_continuous(expand = c(0,0), breaks = s1_x_centers, labels = s1_ordered_ids) +
      scale_y_continuous(expand = c(0,0), breaks = s1_y_centers, labels = rev(s1_ordered_ids)) +
      coord_equal() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
            axis.text.y = element_text(size = 7)) +
      labs(title = "Stage 1: Woody vs Herbs Variant Similarity (LDA-weighted)",
           x = NULL, y = NULL)

    tryCatch({
      ggsave(file.path(OUT_DIR, "variant_similarity_stage1_woody_herbs.png"), p_heat_s1, width = 10, height = 8)
      write.csv(s1_sim_mat, file.path(OUT_DIR, "variant_similarity_stage1_woody_herbs.csv"))
      cat(sprintf("Saved Stage 1 (woody vs herbs) similarity heatmap to: %s\n",
                  file.path(OUT_DIR, "variant_similarity_stage1_woody_herbs.png")))

      # Print cross-class similarity statistics
      woody_idx_local <- which(s1_stage1_class == "woody")
      herbs_idx_local <- which(s1_stage1_class == "herbs")

      if (length(woody_idx_local) > 0 && length(herbs_idx_local) > 0) {
        cross_sims_s1 <- s1_sim_mat[woody_idx_local, herbs_idx_local, drop = FALSE]
        cat(sprintf("\n[STAGE 1: WOODY-HERBS CROSS-SIMILARITY STATS]\n"))
        cat(sprintf("  Min: %.3f, Max: %.3f, Mean: %.3f, Median: %.3f\n",
                    min(cross_sims_s1), max(cross_sims_s1), mean(cross_sims_s1), median(cross_sims_s1)))
      }
    }, error = function(e) {
      cat(sprintf("[ERROR] Failed to generate Stage 1 similarity heatmap: %s\n", e$message))
    })

    cat("================================================\n\n")
  } else {
    cat(sprintf("[INFO] Skipping Stage 1 heatmap: woody=%d, herbs=%d variants\n",
                length(woody_indices), length(herbs_indices)))
  }

  # === STAGE 2: POPULUS vs TAMARIX SIMILARITY HEATMAP (with Stage 2 LDA weights) ===
  # Generate a dedicated heatmap showing populus vs tamarix differentiation using Stage 2 weights

  # Get Stage 2 LDA weights if available
  stage2_weights <- NULL
  if (exists("HIERARCHICAL_LDA_PARAMS", envir = globalenv())) {
    hlp <- get("HIERARCHICAL_LDA_PARAMS", envir = globalenv())
    if (!is.null(hlp$stage2) && !is.null(hlp$stage2$weights)) {
      stage2_weights <- hlp$stage2$weights
      cat("[STAGE 2 HEATMAP] Using Stage 2 LDA weights for populus vs tamarix similarity\n")
    }
  }

  pt_indices <- which(tolower(all_vegs) %in% tolower(woody_classes))

  if (length(pt_indices) >= 2) {
    cat("\n=== STAGE 2: POPULUS vs TAMARIX VARIANT SIMILARITY ===\n")

    pt_vecs <- all_vecs[pt_indices]
    pt_ids <- all_ids[pt_indices]
    pt_vegs <- all_vegs[pt_indices]
    pt_nsamples <- all_nsamples[pt_indices]

    # Build matrix and apply Stage 2 weights if available
    pt_vec_mat <- do.call(rbind, pt_vecs)
    pt_vec_mat[!is.finite(pt_vec_mat)] <- 0

    if (!is.null(stage2_weights) && length(stage2_weights) == ncol(pt_vec_mat)) {
      # Apply Stage 2 LDA weights
      pt_vec_mat_weighted <- sweep(pt_vec_mat, 2, sqrt(stage2_weights), "*")
      cat(sprintf("  Applied Stage 2 weights to %d features\n", length(stage2_weights)))
    } else {
      pt_vec_mat_weighted <- pt_vec_mat
      if (!is.null(stage2_weights)) {
        cat(sprintf("[WARNING] Stage 2 weights length (%d) != features (%d), using unweighted\n",
                    length(stage2_weights), ncol(pt_vec_mat)))
      }
    }

    # Normalize and compute similarity
    pt_row_norms <- sqrt(rowSums(pt_vec_mat_weighted^2))
    pt_row_norms[pt_row_norms == 0] <- 1
    pt_vec_mat_norm <- pt_vec_mat_weighted / pt_row_norms
    pt_sim_mat <- tcrossprod(pt_vec_mat_norm)
    pt_sim_mat[pt_sim_mat > 1.01] <- 1.01
    pt_sim_mat[pt_sim_mat < -0.01] <- -0.01

    rownames(pt_sim_mat) <- pt_ids
    colnames(pt_sim_mat) <- pt_ids

    pt_sim_df <- as.data.frame(as.table(pt_sim_mat))
    colnames(pt_sim_df) <- c("Var1", "Var2", "Similarity")

    # Order by vegetation type
    pt_ord_idx <- order(pt_vegs, pt_ids)
    pt_ordered_ids <- pt_ids[pt_ord_idx]
    pt_ordered_vegs <- pt_vegs[pt_ord_idx]

    # Use sample counts for weighting
    pt_count_map <- setNames(pt_nsamples, pt_ids)
    pt_counts_ordered <- as.numeric(pt_count_map[pt_ordered_ids])
    pt_counts_ordered[is.na(pt_counts_ordered)] <- 0
    pt_total_counts <- sum(pt_counts_ordered)

    if (!is.finite(pt_total_counts) || pt_total_counts <= 0) {
      pt_n_ord <- length(pt_ordered_ids)
      pt_widths_prop <- rep(1 / pt_n_ord, pt_n_ord)
    } else {
      pt_widths_prop <- pt_counts_ordered / pt_total_counts
      eps <- 1e-6
      pt_widths_prop[pt_widths_prop == 0] <- eps
      pt_widths_prop <- pt_widths_prop / sum(pt_widths_prop)
    }

    pt_cumw <- cumsum(pt_widths_prop)
    pt_lefts <- c(0, pt_cumw[-length(pt_cumw)])
    pt_rights <- pt_cumw
    names(pt_lefts) <- pt_ordered_ids
    names(pt_rights) <- pt_ordered_ids

    pt_lefts_rev <- rev(pt_lefts)
    pt_rights_rev <- rev(pt_rights)
    names(pt_lefts_rev) <- rev(pt_ordered_ids)
    names(pt_rights_rev) <- rev(pt_ordered_ids)

    pt_sim_df$xmin <- pt_lefts[as.character(pt_sim_df$Var1)]
    pt_sim_df$xmax <- pt_rights[as.character(pt_sim_df$Var1)]
    pt_sim_df$ymin <- pt_lefts_rev[as.character(pt_sim_df$Var2)]
    pt_sim_df$ymax <- pt_rights_rev[as.character(pt_sim_df$Var2)]

    pt_x_centers <- (pt_lefts + pt_rights) / 2
    pt_y_centers <- (pt_lefts_rev + pt_rights_rev) / 2

    # Boundaries at vegetation type changes
    if (length(pt_ordered_vegs) > 1) {
      pt_change_idx <- which(pt_ordered_vegs[-1] != pt_ordered_vegs[-length(pt_ordered_vegs)])
      pt_vlines_x <- cumsum(pt_widths_prop)[pt_change_idx]
      pt_hlines_y <- cumsum(pt_widths_prop)[pt_change_idx]
    } else {
      pt_vlines_x <- numeric(0)
      pt_hlines_y <- numeric(0)
    }

    # Build Stage 2 heatmap
    p_heat_pt <- ggplot(pt_sim_df) +
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Similarity), color = NA) +
      geom_vline(xintercept = pt_vlines_x, color = "black", size = 1.2) +
      geom_hline(yintercept = pt_hlines_y, color = "black", size = 1.2) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.8, limits = c(-0.01, 1.01), na.value = "white") +
      scale_x_continuous(expand = c(0,0), breaks = pt_x_centers, labels = pt_ordered_ids) +
      scale_y_continuous(expand = c(0,0), breaks = pt_y_centers, labels = rev(pt_ordered_ids)) +
      coord_equal() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
            axis.text.y = element_text(size = 8)) +
      labs(title = "Stage 2: Populus vs Tamarix Variant Similarity (LDA-weighted)",
           x = NULL, y = NULL)

    tryCatch({
      ggsave(file.path(OUT_DIR, "variant_similarity_stage2_populus_tamarix.png"), p_heat_pt, width = 10, height = 8)
      write.csv(pt_sim_mat, file.path(OUT_DIR, "variant_similarity_stage2_populus_tamarix.csv"))
      cat(sprintf("Saved Stage 2 (populus vs tamarix) similarity heatmap to: %s\n",
                  file.path(OUT_DIR, "variant_similarity_stage2_populus_tamarix.png")))

      # Print cross-class similarity statistics
      populus_idx <- which(tolower(pt_vegs) == "populus")
      tamarix_idx <- which(tolower(pt_vegs) == "tamarix")

      if (length(populus_idx) > 0 && length(tamarix_idx) > 0) {
        cross_sims <- pt_sim_mat[populus_idx, tamarix_idx, drop = FALSE]
        cat(sprintf("\n[STAGE 2: POPULUS-TAMARIX CROSS-SIMILARITY STATS]\n"))
        cat(sprintf("  Min: %.3f, Max: %.3f, Mean: %.3f, Median: %.3f\n",
                    min(cross_sims), max(cross_sims), mean(cross_sims), median(cross_sims)))

        # Identify most similar pairs
        top_n <- min(5, length(cross_sims))
        cross_df <- data.frame(
          Populus = rep(rownames(cross_sims), ncol(cross_sims)),
          Tamarix = rep(colnames(cross_sims), each = nrow(cross_sims)),
          Similarity = as.vector(cross_sims)
        )
        cross_df <- cross_df[order(-cross_df$Similarity), ]
        cat(sprintf("  Top %d most similar populus-tamarix pairs:\n", top_n))
        for (i in seq_len(top_n)) {
          cat(sprintf("    %s <-> %s: %.3f\n",
                      cross_df$Populus[i], cross_df$Tamarix[i], cross_df$Similarity[i]))
        }

        # Store high-similarity (>0.9) populus-tamarix pairs for "woody_unknown" renaming
        WOODY_UNKNOWN_THRESHOLD <- 0.9
        high_sim_pairs <- cross_df[cross_df$Similarity > WOODY_UNKNOWN_THRESHOLD, , drop = FALSE]
        if (nrow(high_sim_pairs) > 0) {
          # Collect all variant IDs that should be renamed to "woody_unknown"
          woody_unknown_variants <- unique(c(as.character(high_sim_pairs$Populus),
                                              as.character(high_sim_pairs$Tamarix)))
          assign("WOODY_UNKNOWN_VARIANTS", woody_unknown_variants, envir = globalenv())
          cat(sprintf("\n[WOODY UNKNOWN] Identified %d variants with cross-class similarity > %.2f:\n",
                      length(woody_unknown_variants), WOODY_UNKNOWN_THRESHOLD))
          for (vid in woody_unknown_variants) {
            cat(sprintf("    - %s\n", vid))
          }
        } else {
          assign("WOODY_UNKNOWN_VARIANTS", character(0), envir = globalenv())
          cat(sprintf("\n[WOODY UNKNOWN] No variants exceeded cross-class similarity threshold of %.2f\n",
                      WOODY_UNKNOWN_THRESHOLD))
        }
      }
    }, error = function(e) {
      cat(sprintf("[ERROR] Failed to generate Stage 2 similarity heatmap: %s\n", e$message))
    })

    cat("================================================\n\n")
  } else {
    cat(sprintf("[INFO] Skipping Stage 2 heatmap: only %d variants found for woody classes\n", length(pt_indices)))
  }

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
  stop("[PPI ERROR] GLOBAL DVI soil baseline functionality has been removed. Compute per-location DJF medians or provide explicit dvi_soil to add_ppi_columns().")
}

# Helper: add a shaded rectangle covering excluded years (1992-1999) for ggplot2 time-series plots.
# Usage: + add_excluded_years_shade(is_date = TRUE)  # x axis is Date
#        + add_excluded_years_shade(is_date = FALSE) # x axis is numeric (year)
add_excluded_years_shade <- function(start_year = 1992, end_year = 1999, is_date = FALSE, fill = "grey70", alpha = 0.35) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  if (is_date) {
    xmin <- as.Date(paste0(start_year, "-01-01"))
    xmax <- as.Date(paste0(end_year, "-12-31"))
    return(ggplot2::annotate("rect", xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill, alpha = alpha))
  } else {
    return(ggplot2::annotate("rect", xmin = start_year, xmax = end_year, ymin = -Inf, ymax = Inf, fill = fill, alpha = alpha))
  }
}

# Helper: add vertical lines at specified years and label them with letters and the year in brackets (e.g., A (2007))
# Usage: + add_year_lines()  # numeric year x-axis
#        + add_year_lines(is_date = TRUE) # Date x-axis
add_year_lines <- function(years = c(2007, 2010, 2014), labels = NULL, is_date = FALSE, color = "black", linetype = "dashed", size = 0.6, text_size = 3, text_vjust = -0.5) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  if (is.null(labels)) {
    # Default to letters (A,B,C,...) and append the year in brackets
    base_labels <- LETTERS[seq_along(years)]
    labels <- paste0(base_labels, " (", years, ")")
  } else {
    # If user provided labels, still append years in brackets for clarity
    labels <- paste0(labels, " (", years, ")")
  }
  xs <- if (is_date) as.Date(paste0(years, "-01-01")) else as.numeric(years)
  layers <- unlist(lapply(seq_along(xs), function(i) list(
    ggplot2::geom_vline(xintercept = xs[i], color = color, linetype = linetype, size = size),
    ggplot2::annotate("text", x = xs[i], y = Inf, label = labels[i], vjust = text_vjust, size = text_size)
  )), recursive = FALSE)
  return(layers)
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


  # Cache loading removed - library should already exist in session

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
