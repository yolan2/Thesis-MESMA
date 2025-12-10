# Script to create vegetation fraction trends graph and CSV for 2017-2025

library(ggplot2)
library(dplyr)
library(openxlsx)

# Define paths
OUTPUT_DIR <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/phenology_results"
OUT_DIR <- file.path(OUTPUT_DIR, "veg_mixture_fit")
RESULTS_FILE <- file.path(OUT_DIR, "mesma_results.xlsx")

# Load the results
if (!file.exists(RESULTS_FILE)) {
  stop("Results file not found: ", RESULTS_FILE)
}

# Get all sheet names
sheets <- getSheetNames(RESULTS_FILE)

# Exclude summary sheets, keep only location IDs (numeric)
location_sheets <- sheets[!sheets %in% c("Summary", "Diagnostics", "Variant_Summary", "Uncertainty", "Global_Pattern", "Vegetation_Trends", "Trend_Bootstrap_CI")]

# Read all location data
all_data <- list()
for (sheet in location_sheets) {
  tryCatch({
    data <- read.xlsx(RESULTS_FILE, sheet = sheet)
    # Find the data rows (skip quality metrics)
    data_rows <- which(data$QUALITY.METRICS == "location_id") + 1
    if (length(data_rows) > 0) {
      data <- data[data_rows:nrow(data), ]
        colnames(data) <- c("location_id", "year", "condition_number", "residual_sum_of_squares", "r_squared", "vegetated_fraction", "barren_fraction", NA, NA)
        # Accept either 'pheno_year' or 'year' column; normalize to 'pheno_year' for consistency
        if ("pheno_year" %in% names(data)) {
          time_col <- "pheno_year"
        } else if ("year" %in% names(data)) {
          time_col <- "year"
        } else {
          stop(sprintf("No 'year' or 'pheno_year' column found in sheet %s", sheet))
        }
        data <- data %>% select(location_id, all_of(time_col), vegetated_fraction, barren_fraction)
        # Rename the temporal column to 'pheno_year' if necessary
        if (time_col != "pheno_year") names(data)[names(data) == time_col] <- "pheno_year"
        data <- data %>% mutate(across(everything(), as.numeric))
      all_data[[sheet]] <- data
    }
  }, error = function(e) {
    cat("Error reading sheet", sheet, ":", e$message, "\n")
  })
}

# Combine all data
combined_data <- bind_rows(all_data, .id = "sheet")

# If the MESMA output 'year' column reflects phenological year (March 1 start), use it as 'pheno_year'
if (!"pheno_year" %in% names(combined_data)) combined_data$pheno_year <- combined_data$year

# Filter for years 2017-2025 and calculate vegetation fraction trends
veg_fraction_data <- combined_data %>%
  filter(pheno_year >= 2017, pheno_year <= 2025) %>%
  group_by(pheno_year) %>%
  summarize(
    mean_vegetated_fraction = mean(vegetated_fraction, na.rm = TRUE),
    sd_vegetated_fraction = sd(vegetated_fraction, na.rm = TRUE),
    n_locations = n(),
    .groups = "drop"
  )

# Create the plot
  p <- ggplot(veg_fraction_data, aes(x = pheno_year, y = mean_vegetated_fraction)) +
  geom_line(color = "darkgreen", size = 1.2) +
  geom_point(color = "darkgreen", size = 3) +
  geom_errorbar(aes(ymin = mean_vegetated_fraction - sd_vegetated_fraction,
                    ymax = mean_vegetated_fraction + sd_vegetated_fraction),
                width = 0.2, color = "darkgreen") +
  labs(
    title = "Vegetation Fraction Trends (2017-2025)",
    x = "Year",
    y = "Mean Vegetated Fraction"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  scale_x_continuous(breaks = 2017:2025)

# Save the plot
plot_file <- file.path(OUT_DIR, "vegetation_fraction_trends_2017_2025.png")
ggsave(plot_file, p, width = 10, height = 6, dpi = 300)

# Save the data to CSV
csv_file <- file.path(OUT_DIR, "vegetation_fraction_trends_2017_2025.csv")
write.csv(veg_fraction_data, csv_file, row.names = FALSE)

cat("Graph saved to:", plot_file, "\n")
cat("CSV data saved to:", csv_file, "\n")
cat("Data summary:\n")
print(veg_fraction_data)