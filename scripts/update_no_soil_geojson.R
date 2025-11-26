#!/usr/bin/env Rscript
# Update geojson: set `no soil` = 0 for features with Veg == 'barren'
# Usage:
# Rscript scripts/update_no_soil_geojson.R "C:/path/to/file.geojson" [--inplace] [--backup-dir C:/path/to/backups]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript scripts/update_no_soil_geojson.R <input_geojson>")
  quit(status = 2)
}

infile <- args[1]
opt_inplace <- any(grepl("--inplace", args))
opt_backup_dir <- NULL
bk_arg_idx <- grep("--backup-dir", args)
if (length(bk_arg_idx)) {
  # expect next token to be path
  if (length(args) >= bk_arg_idx + 1) opt_backup_dir <- args[bk_arg_idx + 1]
}

if (!file.exists(infile)) {
  stop(sprintf("Input file not found: %s", infile))
}

if (!requireNamespace("sf", quietly = TRUE)) {
  stop("This script requires the 'sf' package. Install with: install.packages('sf')")
}

library(sf)

# Read geojson
cat(sprintf("Reading: %s\n", infile))
geo <- tryCatch({ sf::st_read(infile, quiet = TRUE) }, error = function(e) stop(e$message))

# Ensure property columns exist via names(geo)
nm <- names(geo)
cat(sprintf("Found %d properties (columns)\n", length(nm)))

## Look for a vegetation property among common names, prefer order:
## 'vegetation', 'veg', 'class' (case-insensitive). This covers GeoJSON files
## that use different field names. If none present, abort.
veg_col <- NULL
matched_cols <- names(geo)[tolower(names(geo)) %in% c("vegetation", "veg", "class")]
if (length(matched_cols) > 0) {
  veg_col <- matched_cols[1]
  cat(sprintf("Using vegetation column: %s\n", veg_col))
} else {
  stop("No vegetation column found in geojson properties (looked for: 'vegetation','veg','class'). Aborting.")
}

# Determine no soil column name(s); ensure both 'no soil' and 'no_soil' are present
no_soil_col <- NULL
if (any(tolower(nm) == "no soil")) {
  no_soil_col <- nm[which(tolower(nm) == "no soil")[1]]
} else {
  no_soil_col <- "no soil"
  geo[[no_soil_col]] <- NA_real_
}

# Also ensure an underscore variant exists so downstream tools/readers that
# sanitize property names still see the information (write both fields).
underscore_col <- NULL
if (any(tolower(nm) == "no_soil")) {
  underscore_col <- nm[which(tolower(nm) == "no_soil")[1]]
} else if (!("no_soil" %in% names(geo))) {
  underscore_col <- "no_soil"
  geo[[underscore_col]] <- as.numeric(geo[[no_soil_col]])
} else {
  underscore_col <- "no_soil"
}

# Count rows to modify
vals <- tolower(as.character(geo[[veg_col]]))
is_barren <- !is.na(vals) & vals == "barren"
count_barren <- sum(is_barren)
cat(sprintf("Detected %d features with Veg == 'barren'\n", count_barren))

if (count_barren == 0) {
  cat("No updates required. Exiting gracefully.\n")
  quit(status = 0)
}

# Create backup
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
basename_in <- basename(infile)
backup_dir <- if (!is.null(opt_backup_dir)) opt_backup_dir else dirname(infile)
if (!dir.exists(backup_dir)) dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)
backup_path <- file.path(backup_dir, paste0(basename_in, ".bak.", ts))
file.copy(infile, backup_path, overwrite = FALSE)
cat(sprintf("Created backup: %s\n", backup_path))

# Set both `no soil` and `no_soil` to numeric 0 for barren features
geo[[no_soil_col]][is_barren] <- 0
geo[[underscore_col]][is_barren] <- 0

# Write output
out_path <- if (opt_inplace) infile else file.path(dirname(infile), paste0("updated_", basename_in))
cat(sprintf("Writing updated geojson to: %s\n", out_path))
tryCatch({
  sf::st_write(geo, out_path, driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
  cat("Write complete.\n")
}, error = function(e) {
  stop(sprintf("Failed to write updated geojson: %s", e$message))
})

cat(sprintf("Updated %d features: set '%s' = 0 for Veg == 'barren'\n", count_barren, no_soil_col))

# Exit success
quit(status = 0)

# Post-processing: some GDAL/GEOJSON writers sanitize property names (spaces -> underscores).
# To ensure both forms are present in the saved GeoJSON we read the written file and
# add a mirrored property where one exists but the other is missing.
try({
  if (tolower(tools::file_ext(out_path)) %in% c("geojson", "json")) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      warning("jsonlite not installed; cannot ensure both 'no soil' and 'no_soil' keys are present in output file.")
    } else {
      j <- jsonlite::fromJSON(out_path, simplifyVector = FALSE)
      if (is.list(j) && !is.null(j$features) && length(j$features) > 0) {
        modified <- FALSE
        for (i in seq_along(j$features)) {
          props <- j$features[[i]]$properties
          if (is.null(props)) next
          has_spaced <- any(names(props) == 'no soil')
          has_under <- any(names(props) == 'no_soil')
          # If underscore exists but spaced name not present, duplicate
          if (!has_spaced && has_under) {
            j$features[[i]]$properties[['no soil']] <- props[['no_soil']]
            # Also add a sanitized variant with dot (used by sf on read)
            j$features[[i]]$properties[['no.soil']] <- props[['no_soil']]
            modified <- TRUE
          }
          # If spaced exists but underscore not present, duplicate and add dot variant
          if (!has_under && has_spaced) {
            j$features[[i]]$properties[['no_soil']] <- props[['no soil']]
            j$features[[i]]$properties[['no.soil']] <- props[['no soil']]
            modified <- TRUE
          }
          # Ensure dot variant exists if either of the other forms exist
          if (!has_spaced && !has_under && !is.null(props[['no.soil']])) {
            # nothing to do explicitly
          } else if (!has_spaced && has_under && is.null(props[['no.soil']])) {
            j$features[[i]]$properties[['no.soil']] <- props[['no_soil']]
            modified <- TRUE
          } else if (has_spaced && !has_under && is.null(props[['no.soil']])) {
            j$features[[i]]$properties[['no.soil']] <- props[['no soil']]
            modified <- TRUE
          }
        }
        if (modified) {
          # Write back (pretty for readability)
          json_text <- jsonlite::toJSON(j, auto_unbox = TRUE, pretty = TRUE)
          writeLines(json_text, out_path)
          cat(sprintf("Post-processed GeoJSON: ensured both 'no_soil' and 'no soil' exist in %s\n", out_path))
        }
      }
    }
  }
}, silent = TRUE)
