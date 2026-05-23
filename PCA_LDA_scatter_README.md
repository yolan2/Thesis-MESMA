# PCA-LDA Scatter Plot Script

This repository includes a helper script (`plot_pca_lda_scatter.R`) that
visualizes training endmember samples in the PCA-LDA feature space.

## Purpose

- **Goal:** produce a 2‑D scatter plot of `loc-year` training samples
  projected onto the first two linear discriminants from an LDA trained
  on PCA-transformed features.
- **Use case:** inspect class separability and check that endmember
  training locations cluster according to vegetation type in the
  discriminant space.

## How it works

1. **Dependencies and configuration**
   - Loads tidyverse packages (`dplyr`, `ggplot2`) and `MASS` for LDA.
   - Sources `mesma_config.R` to inherit global constants such as
     `TRAIN_YEARS`, `OPTIMAL_INDICES`, etc.

2. **Data loading and feature selection**
   - Reads `preprocessed_data.rds` (output from preprocessing pipeline)
     and `training_norm_params.rds` for feature normalization scales.
   - Selects a set of feature indices (`avail`) using available norm
     parameters or default OPTIMAL_INDICES that exist as columns in the
     data frame.

3. **Filtering and labelling**
   - Restricts rows to the phenological training years (`TRAIN_YEARS`).
   - Computes `pheno_year` and `doy` if missing.
   - Filters to allowed vegetation classes plus `barren` and adds a
     `target_class` column (lower-cased vegetation label).

4. **Trace construction**
   - Splits the filtered data by `location_id` and `pheno_year` into
     individual loc-year traces.
   - For each trace with at least `MIN_PENTADS_PER_TRAIN_SAMPLE` rows,
     builds a pentad-averaged matrix (one row per pentad, columns =
     selected indices).  The final pentad is left `NA` per MESMA
     convention.
   - Fills missing interior pentads by linear interpolation.
   - Flattens the matrix to a vector and L2-normalizes it (if
     `ENABLE_LDA_L2_NORMALIZATION` is TRUE, which it is by default).

5. **Z-scoring**
   - After collecting all vectors into a matrix (`X_mat`), the script
     optionally z-scores each index×pentad column group (duplicate of
     training logic with `ENABLE_ZSCORE_AFTER_L2 = TRUE`).
   - Non-finite entries are then zeroed.

6. **PCA**
   - Performs PCA on the z-scored matrix.
   - Retains enough principal components to cover
     `PCA_VARIANCE_THRESHOLD` (default 0.95) of variance, but enforces
     LDA stability caps: at least (`n_classes - 1`) PCs, at most
     `min(n_min - 2, 20)` where `n_min` is the smallest class count.

7. **LDA**
   - Runs `MASS::lda` on the selected PCA scores using the vegetation
     class labels.
   - Suppresses collinearity warnings as done in training.
   - Computes the projection of all samples into the LDA space
     (`X_ld`) and records explained variance percentages.

8. **Plotting**
   - Constructs a data frame containing LD1/LD2 coordinates, class,
     location, year, and `loc_year` composite label.
   - Chooses a colour palette keyed to vegetation classes, adding any
     unexpected classes with rainbow colours.
   - Shapes points by year (2023 circle, 2024 triangle, etc.) for visual
     encoding of training time.
   - Draws points with optional 80% confidence ellipses per class.
   - Labels axes with percent variance and sets a subtitle summarizing
     sample counts.
   - Exports the plot to `training_lda_scatter.png` (9×6.5 in at 150 DPI).

9. **Diagnostics output**
   - After saving, prints sample counts per class and class-centroid
     coordinates along LD1 and LD2.

## Running

From the repository root run:

```powershell
Rscript plot_pca_lda_scatter.R
```

The generated image file will appear in the same directory.

---

This lightweight explanation can be adapted or expanded as needed for
documentation or thesis write‑up.