# Quick demo: how veg_ref resolution works in the MESMA script
# This demo shows the same fallback behavior implemented in compress_and_unmix_year
# It does NOT source the full main script to avoid execution of heavy initialization.

# Example 1: dly_year contains Veg
# Example 2: dly_year lacks Veg and gpts_map misses mapping -> fallback to GLOBAL_PCA
dly1 <- data.frame(location_id = "L_1_1", Veg = NA, date = as.Date("2020-04-01"))

# Simulate a location -> Veg mapping (gpts_map)
gpts_map <- data.frame(location_id = c("L_1_1", "L_3_3"), Veg = c("Grassland", "Shrub"), stringsAsFactors = FALSE)

resolve_veg_ref <- function(dly_year, gpts_map = NULL) {
  veg_ref <- if ("Veg" %in% names(dly_year)) unique(na.omit(as.character(dly_year$Veg)))[1] else NA_character_
  if (is.na(veg_ref) || !nzchar(as.character(veg_ref))) {
    # Try cross-reference using gpts_map
    if ("location_id" %in% names(dly_year) && !is.null(gpts_map) && is.data.frame(gpts_map) && "location_id" %in% names(gpts_map) && "Veg" %in% names(gpts_map)) {
      locs <- unique(na.omit(as.character(dly_year$location_id)))
      if (length(locs) > 0) {
        map_ids <- as.character(gpts_map$location_id)
        match_idx <- match(locs[1], map_ids)
        if (is.na(match_idx)) match_idx <- match(tolower(locs[1]), tolower(map_ids))
        if (!is.na(match_idx)) {
          mapped <- as.character(gpts_map$Veg[match_idx])
          if (!is.na(mapped) && nzchar(mapped)) return(tolower(mapped))
        }
      }
    }
    cat("[WARN] veg not defined for this row and cross-reference failed; returning NA\n")
    return(NA_character_)
  }
  tolower(as.character(veg_ref))
}

cat("Demo1 -> veg_ref (cross-ref): ", resolve_veg_ref(dly1, gpts_map = gpts_map), "\n")

# Example 2: dly_year lacks Veg and gpts_map misses mapping -> should warn
dly2 <- data.frame(location_id = "L_2_2", Veg = NA, date = as.Date("2020-04-02"))
cat("Demo2 -> veg_ref (cross-ref): ", resolve_veg_ref(dly2, gpts_map = gpts_map), "\n")

# Example 3: dly_year contains Veg
dly3 <- data.frame(location_id = "L_3_3", Veg = "Shrub", date = as.Date("2020-04-03"))
cat("Demo3 -> veg_ref (cross-ref): ", resolve_veg_ref(dly3, gpts_map = gpts_map), "\n")

# Instructions: Run this script from R using source('examples/quick_veg_ref_demo.R') to verify resolution logic.
