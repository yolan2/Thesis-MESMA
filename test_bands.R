
source("ppi_helpers.R")

normalize_band_names <- function(df, bands = c("blue", "green", "red", "nir", "swir1", "swir2")) {
  current_names <- names(df)
  for (b in bands) {
    candidates <- c(b, toupper(b), tools::toTitleCase(b), paste0("band_", b), toupper(paste0("band_", b)), paste0("Band_", b))
    for (cand in candidates) {
      if (cand %in% current_names && !(b %in% current_names)) {
        names(df)[names(df) == cand] <- b
        current_names <- names(df)
        break
      }
    }
  }
  df
}

df <- read.csv("C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_low (3).csv")[1:10, ]
df <- normalize_band_names(df)
df <- compute_indices_from_bands(df)
print(df[1, c("blue", "green", "red", "nir", "swir1", "swir2", "EVI", "PSRI", "NDMI", "NDTI", "MSI", "MSAVI")])

