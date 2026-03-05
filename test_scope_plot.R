library(dplyr)
library(ggplot2)
INDICES_OF_INTEREST <- c("MSAVI","NDVI","PPI","OSAVI","NIRv","NBR","TCW","NDMI","TCB","GVI","EVI")
plot_global_averages <- function(summary_df, title_prefix, output_dir) {
  if (is.null(summary_df) || nrow(summary_df) == 0) return(invisible(NULL))
  idxs <- intersect(INDICES_OF_INTEREST,
                    gsub("^avg_", "", grep("^avg_", names(summary_df), value = TRUE)))
  if (length(idxs) == 0) return(invisible(NULL))

  long <- do.call(rbind, lapply(idxs, function(idx) {
    data.frame(index = idx,
               mean = summary_df[[paste0("avg_", idx)]],
               ci_lower = summary_df[[paste0(idx, "_ci_lower")]],
               ci_upper = summary_df[[paste0(idx, "_ci_upper")]],
               stringsAsFactors = FALSE)
  }))

  p <- ggplot(long, aes(x = index, y = mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
    theme_minimal() +
    labs(title = paste0(title_prefix, " global averages with 95% CI"),
         x = "Index", y = "Mean value")

  fn <- file.path(output_dir, paste0(gsub("\\s+", "_", title_prefix), "_global_CI.png"))
  ggsave(fn, plot = p, width = 10, height = 6)
  cat("Saved global CI plot to:", fn, "\n")
}

# read sample csv and attempt to plot
scope_dir <- "C:/MAP/january_averages_results/Landsat_Harmonized_Bands_1985_2025_train__3_/low"
df <- read.csv(file.path(scope_dir, "Winter_Global.csv"))
print(names(df))
plot_global_averages(df, 'scope_test', scope_dir)
