
source("mesma_config.R")
source("mesma_helpers.R")

df <- read.csv("C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_low (3).csv")
df_pre <- df[df$year < 2015, ][1:100, ]
df_post <- df[df$year >= 2015, ][1:100, ]

for(d in list(df_pre, df_post)) {
  d <- normalize_band_names(d)
  d <- compute_indices_from_bands(d)
  
  cat("Mean EVI raw:", mean(d$EVI, na.rm=T), "\n")
}

