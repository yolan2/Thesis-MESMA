# Demo: Show adaptive scaling of barren vs vegetation in plots
# Source the main functions file
source(file.path('..', 'fit_veg_mixture_mesma.R'))

library(ggplot2)

years <- seq(2017, 2024)
veg_types <- c('populus', 'phragmites', 'barren')

# create global patterns with veg small (0.05-0.2) and barren big (0.6-0.85)
set.seed(42)
rows <- expand.grid(year = years, Veg = veg_types, stringsAsFactors = FALSE)
rows$global_coef <- NA_real_
rows$ci_lower <- NA_real_
rows$ci_upper <- NA_real_
rows$n_locations <- 10

for (i in seq_len(nrow(rows))) {
  v <- rows$Veg[i]
  y <- rows$year[i]
  if (v == 'barren') {
    val <- 0.6 + 0.25 * sin((y - 2017) * 0.7 + i/5)
  } else {
    val <- 0.02 + 0.15 * runif(1)
  }
  rows$global_coef[i] <- val
  rows$ci_lower[i] <- pmax(0, val - 0.02)
  rows$ci_upper[i] <- pmin(1, val + 0.02)
}

p <- plot_global_vegetation_pattern(rows, title = 'Demo: Veg vs Barren adaptive scaling')

# Save demo plot
ggsave(file.path('..', 'demo_barren_veg_plot.png'), p, width = 10, height = 6, dpi = 300)
cat('Saved demo plot to: demo_barren_veg_plot.png\n')
