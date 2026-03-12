library(dplyr)
library(tidyr)
library(ggplot2)

# --- Load data ---
kon <- read.csv("C:/MAP/january_averages_results/preprocessed/Landsat_Harmonized_Bands_1985_2025_kon (1)/trend_mean_corrected.csv")
low <- read.csv("C:/MAP/january_averages_results/preprocessed/Landsat_Harmonized_Bands_1985_2025_low (3)/trend_mean_corrected.csv")

# --- Filter PPI only and join on pheno_year ---
joined <- inner_join(
  kon %>% filter(index == "PPI") %>% select(pheno_year, mean_val),
  low %>% filter(index == "PPI") %>% select(pheno_year, mean_val),
  by = "pheno_year",
  suffix = c("_kon", "_low")
)

cat("Matched years:", nrow(joined), "\n\n")

# --- Correlation ---
# One-sided test: H1 = negative correlation (alternative = "less")
ct <- cor.test(joined$mean_val_kon, joined$mean_val_low, method = "pearson", alternative = "less")
cat(sprintf("Pearson r        = %.4f\n", ct$estimate))
cat(sprintf("p-value (H1: r<0) = %.4g\n", ct$p.value))
cat(sprintf("Spearman r       = %.4f\n\n", cor(joined$mean_val_kon, joined$mean_val_low, method = "spearman")))
cat(if (ct$p.value < 0.05) "=> Significant negative correlation\n" else "=> No significant negative correlation\n")

# --- Scatter plot ---
p <- ggplot(joined, aes(x = mean_val_kon, y = mean_val_low)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("r = %.3f\np (r<0) = %.3g", ct$estimate, ct$p.value),
           hjust = -0.1, vjust = 1.3, size = 4) +
  labs(
    title = "PPI Correlation: Kon vs Low",
    x = "Kon PPI (mean_val)",
    y = "Low PPI (mean_val)"
  ) +
  theme_bw()

ggsave("correlate_kon_low_PPI.png", p, width = 7, height = 5, dpi = 150)
cat("Plot saved to correlate_kon_low_PPI.png\n")
