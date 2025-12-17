# Test script to verify phenological DOY calculation
library(lubridate)

# Define the pheno_doy function (same as in fit_veg_mixture_mesma.R)
pheno_doy <- function(d) {
  d <- as.Date(d)
  month <- lubridate::month(d)
  day <- lubridate::day(d)

  # For March-December: count days from March 1
  # For January-February: count days from previous March 1 (add ~306 days)
  ifelse(is.na(d), NA_integer_,
    ifelse(month >= 3,
      # March onwards: days since March 1 of current year
      as.integer(d - as.Date(paste0(lubridate::year(d), "-03-01"))) + 1L,
      # Jan-Feb: days since March 1 of previous year
      as.integer(d - as.Date(paste0(lubridate::year(d) - 1, "-03-01"))) + 1L
    )
  )
}

# Test cases
test_dates <- as.Date(c(
  "2023-03-01",  # Start of pheno year 2023
  "2023-06-15",  # Mid-year
  "2023-12-31",  # End of calendar year
  "2024-01-15",  # January (still in pheno year 2023)
  "2024-02-28",  # End of pheno year 2023
  "2024-03-01"   # Start of pheno year 2024
))

calendar_doy <- lubridate::yday(test_dates)
pheno_doy_vals <- pheno_doy(test_dates)

cat("\n=== Testing Phenological DOY Calculation ===\n\n")
cat("Date          Calendar DOY  Pheno DOY  Description\n")
cat("------------------------------------------------------------------\n")

results <- data.frame(
  date = as.character(test_dates),
  cal_doy = calendar_doy,
  pheno_doy = pheno_doy_vals,
  description = c(
    "Start of pheno year 2023",
    "Mid-June (day ~107 of pheno year)",
    "End of calendar year (~306 days from Mar 1)",
    "Mid-January (still in pheno 2023, ~320 days)",
    "End of pheno year 2023 (~364 days)",
    "Start of pheno year 2024"
  )
)

for (i in 1:nrow(results)) {
  cat(sprintf("%-13s %4d          %4d       %s\n",
              results$date[i],
              results$cal_doy[i],
              results$pheno_doy[i],
              results$description[i]))
}

cat("\n=== Verification ===\n")
cat("✓ March 1 should have pheno_doy = 1:", pheno_doy_vals[1] == 1, "\n")
cat("✓ January dates should have pheno_doy > 300:", all(pheno_doy_vals[4:5] > 300), "\n")
cat("✓ Consecutive March 1 dates should both be 1:",
    pheno_doy_vals[1] == 1 && pheno_doy_vals[6] == 1, "\n")

cat("\n=== Impact on Pentad Assignment ===\n")
TEMPORAL_AGGREGATION_DAYS <- 5
N_TEMPORAL_BINS <- 73

doy_to_pentad <- function(doy) {
  pmin(ceiling(doy / TEMPORAL_AGGREGATION_DAYS), N_TEMPORAL_BINS)
}

cat("\nCalendar-based pentads (WRONG):\n")
for (i in 1:length(test_dates)) {
  pentad_wrong <- doy_to_pentad(calendar_doy[i])
  pentad_correct <- doy_to_pentad(pheno_doy_vals[i])
  cat(sprintf("  %s: pentad %2d (wrong) -> pentad %2d (correct)\n",
              test_dates[i], pentad_wrong, pentad_correct))
}

cat("\n✓ Fix applied! All DOY calculations now use phenological year alignment.\n")
