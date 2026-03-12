
df <- read.csv("C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_low (3).csv")
df$NDDI <- (df$Red - df$NIR)/(df$Red + df$NIR + 1e-9)
cat("pre-2015 total:", sum(df$year < 2015, na.rm=TRUE), "\n")
cat("pre-2015 < 0.4:", sum(df$year < 2015 & !is.na(df$NDDI) & df$NDDI <= 0.4, na.rm=TRUE), "\n")

