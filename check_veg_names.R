# Check what vegetation names are actually in your data

cat("=== VEGETATION NAME MISMATCH DIAGNOSTIC ===\n\n")

if (exists("df") && "Veg" %in% names(df)) {
  cat("Unique vegetation names in dataframe:\n")
  veg_values <- unique(df$Veg[!is.na(df$Veg)])
  for (v in sort(veg_values)) {
    count <- sum(df$Veg == v, na.rm = TRUE)
    cat(sprintf("  '%s' : %d observations\n", v, count))
  }

  cat("\n")
  cat("Expected names (line 3541 filter):\n")
  cat("  'phragmites', 'populus', 'tamarix'\n\n")

  cat("ALLOWED_VEG (line 50):\n")
  if (exists("ALLOWED_VEG")) {
    cat(sprintf("  %s\n", paste(sprintf("'%s'", ALLOWED_VEG), collapse=", ")))
  } else {
    cat("  ALLOWED_VEG not defined yet\n")
  }

  cat("\n")
  cat("Names that WILL BE KEPT at line 3541:\n")
  filter_set <- c("phragmites", "populus", "tamarix")
  kept_names <- intersect(veg_values, filter_set)
  if (length(kept_names) > 0) {
    cat(sprintf("  %s\n", paste(sprintf("'%s'", kept_names), collapse=", ")))
  } else {
    cat("  NONE - ALL VEGETATION DATA WILL BE FILTERED OUT!\n")
  }

  cat("\n")
  cat("Names that WILL BE REMOVED at line 3541:\n")
  removed_names <- setdiff(veg_values, filter_set)
  removed_names <- removed_names[removed_names != "" & tolower(removed_names) != "barren"]
  if (length(removed_names) > 0) {
    for (v in removed_names) {
      count <- sum(df$Veg == v, na.rm = TRUE)
      cat(sprintf("  '%s' : %d observations LOST\n", v, count))
    }
  } else {
    cat("  None\n")
  }

  cat("\n=== DIAGNOSIS ===\n")
  if (length(kept_names) == 0) {
    cat("ERROR: NO vegetation names match the filter at line 3541!\n")
    cat("This means:\n")
    cat("  1. mesma_lib will be built with NO vegetation types\n")
    cat("  2. At inference, veg_kept will be empty\n")
    cat("  3. Stage 2 will NEVER run\n\n")
    cat("SOLUTION:\n")
    cat("  Either:\n")
    cat("  A) Rename your Veg column values to: 'phragmites', 'populus', or 'tamarix'\n")
    cat("  B) Change line 3541 to match YOUR vegetation names\n")
    cat("  C) Make line 3541 case-insensitive:\n")
    cat("     vegs <- vegs[tolower(vegs) %in% c('phragmites', 'populus', 'tamarix')]\n")
  } else {
    cat(sprintf("OK: %d vegetation type(s) will be kept\n", length(kept_names)))
    if (length(removed_names) > 0) {
      cat(sprintf("WARNING: %d vegetation type(s) will be removed\n", length(removed_names)))
    }
  }

  cat("\n")
  if (exists("lib")) {
    cat("Current lib contains:\n")
    cat(sprintf("  %s\n", paste(names(lib), collapse=", ")))
  }

  if (exists("mesma_lib")) {
    cat("\nCurrent mesma_lib contains:\n")
    cat(sprintf("  %s\n", paste(names(mesma_lib), collapse=", ")))
  }

} else {
  cat("ERROR: df not found or missing 'Veg' column\n")
}
