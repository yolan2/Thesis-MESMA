
library(data.table)

# Soil Endmember from evaluate_index_linearity.R
soil_spec <- c(
  blue  = 0.08,
  green = 0.12,
  red   = 0.18,
  nir   = 0.24,
  swir1 = 0.35,
  swir2 = 0.30
)

# Create a 1-row data.table
df <- as.data.table(t(soil_spec))

# Calculate indices
eps <- 1e-9
df[, `:=`(
    DVI   = nir - red,
    OSAVI = (nir - red) / (nir + red + 0.16),
    MCARI = ((red - green) - 0.2*(red - blue)) * (red / (green + eps)),
    NIRv  = (nir * ((nir - red) / (nir + red + eps))) * 1.3,
    PSRI  = (red - blue) / (nir + eps),
    NBR   = (nir - swir2) / (nir + swir2 + eps),
    TCW   = (swir1 - swir2) / (swir1 + swir2 + eps),
    NDVI   = (nir - red) / (nir + red + eps),
    MSAVI2 = (2 * nir + 1 - sqrt(pmax(0, (2 * nir + 1)^2 - 8 * (nir - red)))) / 2,
    NDMI   = (nir - swir1) / (nir + swir1 + eps),
    TCB    = 0.3029 * blue + 0.2786 * green + 0.4733 * red + 0.5599 * nir + 0.508 * swir1 + 0.1872 * swir2,
    GVI    = -0.2941 * blue - 0.243 * green - 0.5424 * red + 0.7276 * nir + 0.0713 * swir1 - 0.1608 * swir2
)]

# Apply Linearizations
df[, NIRv_lin := 2 * NIRv - NIRv^2]
df[, MSAVI2_lin := log(MSAVI2 + 1)]
df[, PSRI_lin := sign(PSRI) * sqrt(abs(PSRI))]
df[, TCW_lin := log(pmax(TCW + 1, 1e-6))]

# Add PPI check (simplified)
# Assume DVI_max = 0.5, Zenith = 0
ppi_check <- function(dvi, zenith=0, M=0.5, dvi.soil=0.09, G=0.5) {
  d_c <- 0.0336 + 0.0477/cos(zenith)
  Q_E <- d_c + (1 - d_c) * G / cos(zenith)
  K <- 1/(4*Q_E) * (1 + M)/(1 - M)
  res <- - K * log( (M - dvi) / (M - dvi.soil) )
  pmin(pmax(res, 0), 3)
}

df[, PPI := ppi_check(DVI)]

print(t(df))

# Print values for OPTIMAL_INDICES (linearized where applicable)
# "MCARI", "NIRv", "PSRI", "NBR", "TCW", "PPI", "TCB", "NDVI", "MSAVI2", "NDMI", "GVI"
cols <- c("MCARI", "NIRv_lin", "PSRI_lin", "NBR", "TCW_lin", "PPI", "TCB", "NDVI", "MSAVI2_lin", "NDMI", "GVI")
cat("\n=== SOIL BASELINE VALUES ===\n")
print(t(df[, ..cols]))
