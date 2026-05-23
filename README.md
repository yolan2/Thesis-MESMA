# Vegetation Mixture Mapping via MESMA

## Overview
This repository contains an R-based processing pipeline for performing Multiple Endmember Spectral Mixture Analysis (MESMA) on remote sensing data (e.g., Landsat time series) to map vegetation dynamics. It also includes downstream analysis for correlating vegetation fractions with groundwater levels and checking satellite biases.

## Repository Structure

### ⚙️ Configuration & Helpers
*   **`mesma_config.R`**: Central configuration file defining paths, parameters, and global settings for the MESMA pipeline.
*   **`mesma_helpers.R`**: Core helper functions for running MESMA, handling spectral libraries, and calculating fractions.
*   **`ppi_helpers.r`**: Helper functions related to the Pixel Purity Index (PPI) or generic endmember extraction functions.

### 🚀 Processing Pipeline
*   **`preprocess_data.R`**: Handles the initial cleaning, formatting, and preparation of raw satellite time series data. Outputs processed components like `preprocessed_data.rds`.
*   **`january_averages.R`**: Specific script likely calculating winter baseline metrics or removing seasonal bias using January averages.
*   **`fit_veg_mixture_mesma.R`**: The main execution script that fits the vegetation mixture model using the configured MESMA algorithm.

### 📊 Validation & Downstream Analysis
*   **`satellite_bias_check.R`**: Analyzes biases or shifts between different satellite sensors (e.g., Landsat 7 vs. Landsat 8/9).
*   **`gw_correlation.R`**: Correlates the derived vegetation fractions from MESMA with groundwater (gw) data.

### 📈 Plotting & Visualization
*   **`plot_mesma_maps.R`**: Generates spatial maps of the resulting vegetation fractions.
*   **`plot_pca_lda_scatter.R`**: Creates scatter plots mapping the PCA/LDA feature space used to separate the different endmembers (vegetation classes).

---

## Data & Results Directories

*   **`downloads/`**: Raw satellite imagery and time-series data fetched for specific ROIs (e.g., populus, tamarix, phragmites, barren).
*   **`inference_results/`**: Outputs from models, feature shift summaries, and nearest centroid CSV files.
*   **`maps/` & `training_visualizations/`**: Generated plots, maps, and visual QA/QC outputs.
*   **`gw_correlation/` & `satellite_bias_check/`**: Dedicated output folders for specific downstream processing tasks.
*   **`temp_results/`**: Scratch space for intermediate pipeline operations.

---

## Detailed Documentation
For deep dives into specific methodologies used in this project, please refer to the markdown files in the `docs/` directory:

*   **[`docs/BOOTSTRAP_PIPELINE.md`](docs/BOOTSTRAP_PIPELINE.md)**: Details on the bootstrapping approach used in the analysis.
*   **[`docs/COSINE_SIMILARITY.md`](docs/COSINE_SIMILARITY.md)**: Explains the spectral matching (cosine similarity) technique.
*   **[`PCA_LDA_scatter_README.md`](PCA_LDA_scatter_README.md)**: Documentation around the Principal Component Analysis (PCA) and Linear Discriminant Analysis (LDA) transformations for endmember separability.

## Quick Start
1. Ensure your dependencies (`R` environment) are properly configured.
2. Review and adjust parameters in `mesma_config.R`.
3. Run `Rscript preprocess_data.R` to prepare the data.
4. Run `Rscript fit_veg_mixture_mesma.R` or other scripts to execute the MESMA modeling.
