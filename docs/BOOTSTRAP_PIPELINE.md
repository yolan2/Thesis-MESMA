# Bootstrapping pipeline (MESMA)

This document describes the bootstrapping and uncertainty pipeline. Implementation lives in `mesma_bootstrap.R`.

---

## Overview

The pipeline layers several perturbation and resampling steps, applied in this order:

| Step | What it does | Uncertainty source |
|------|-------------|-------------------|
| **1. Cover-index bootstrap** (within-location) | Resamples cover-index observations and recomputes the median | Noise and timing in the cover index |
| **2. Location bootstrap** (between-location) | Resamples locations (or spatial blocks) with replacement | Sampling variability across locations |
| **3. Local residual bootstrap** (optional) | Resamples spectral residuals for per-pixel uncertainty | Per-location unmixing noise |

---


---

## 2. Cover-index bootstrap (within-location)

**Purpose:** Quantify the uncertainty in the vegetation cover index (PPI, NDVI, or MSAVI) for each location and year.

**How it works:** For each location–year, the pipeline collects all summer-season observations of the cover index. It then resamples these observations with replacement and computes the median each time. The median is used because it is robust to outliers.

$$
\tilde{x}_{\ell,y}^{*(b)} = \mathrm{median}\bigl(x_{\ell,y,i_1^*}, \dots, x_{\ell,y,i_n^*}\bigr)
$$

### Seasonal detrending

Before bootstrapping, the cover index can optionally be detrended to remove the effect of sampling timing within the summer season. A polynomial is fit to the day-of-year trend:

$$
x \approx \beta_0 + \beta_1 \cdot \mathrm{DOY} + \beta_2 \cdot \mathrm{DOY}^2 + \beta_3 \cdot \mathrm{DOY}^3
$$

The detrended value subtracts the fitted seasonal pattern (centered at its mean), so that early-summer and late-summer observations become more comparable.

---

## 3. Location bootstrap (between-location)

**Purpose:** Quantify sampling variability by resampling which locations contribute to the global average.

**How it works:** We draw locations (or spatial blocks) with replacement to create a bootstrap resample, then recompute the yearly global mean for each resample. Repeating this many times builds up a distribution of possible global averages.

### Simple location resampling

From the set of all locations, draw $N$ locations with replacement. Some locations will appear multiple times, others not at all. Compute the weighted mean for each year and class.

### Spatial block bootstrap

When locations are spatially clustered, resampling individual locations does not correctly account for spatial correlation. Instead, we group locations into spatial grid cells (blocks) and resample entire blocks.

The block size is estimated from the data by fitting an exponential variogram model. The variogram measures how dissimilar pairs of locations become as the distance between them increases:

$$
\gamma(h) \approx s \cdot \bigl(1 - e^{-h/r}\bigr)
$$

The range parameter $r$ (in km) captures the distance over which locations are correlated, and this becomes the block width.

If spatial coordinates are missing or the variogram fit fails, the pipeline falls back to simple i.i.d. location resampling.

### Aggregation within a bootstrap replicate

For cover-weighted summaries, each bootstrap replicate combines the resampled location weights with the cover values. The yearly global cover for class $k$ is:

$$
\bar{c}_{y,k}^{*(b)} = \frac{1}{N} \sum_{\ell} m_\ell^{*(b)} \cdot c_{\ell,y,k}^{*(b)}
$$

where $m_\ell^{*(b)}$ is the number of times location $\ell$ was drawn in that replicate.

### Derived categories

Bootstrap reporting is now done per-species (e.g. Populus, Tamarix).

---

## 4. Summarizing the bootstrap results

For each year and class, the pipeline stores all $B$ bootstrap replicates of the statistic of interest (either fraction or absolute cover). From these, it reports:

**Standard error** — the standard deviation of the bootstrap replicates:

$$
\widehat{\mathrm{SE}}(T_{y,k}) = \sqrt{\frac{1}{B-1}\sum_{b=1}^{B}\bigl(T_{y,k}^{*(b)} - \bar{T}_{y,k}^*\bigr)^2}
$$

**95% confidence interval** — the 2.5th and 97.5th percentiles of the bootstrap distribution:

$$
\bigl[Q_{0.025}(T_{y,k}^*),\ Q_{0.975}(T_{y,k}^*)\bigr]
$$

All reported fractions are clipped to the range [0, 1].

---

## 5. Local residual bootstrap (optional)

This is a separate, optional bootstrap at the individual location–year level (controlled by `ENABLE_MULTI_YEAR_BOOTSTRAP`). It is designed to quantify uncertainty in a single unmixing result, not in global aggregates.

Given a fitted spectrum and the residuals from that fit, resampled residuals are added back to the fitted spectrum to create a synthetic observation. The unmixing is then re-solved on this synthetic spectrum (keeping the same endmember set). Repeating this gives a distribution of coefficients for that one location–year.

$$
\mathbf{y}^{*(b)} = \hat{\mathbf{y}} + \mathbf{e}^{*(b)}
$$

This avoids repeating the expensive combinatorial endmember search — only the final coefficient solve is re-done.

---

## 6. Configuration reference

Key settings in `mesma_config.R`:

| Setting | What it controls |
|---------|-----------------|
| `BOOTSTRAP_B` | Number of bootstrap replicates |
| `ENABLE_SPATIAL_BLOCK_BOOTSTRAP` | Whether to use spatial block resampling |
| `BOOTSTRAP_BLOCK_KM` | Block size in km (overrides variogram estimate if set) |
| `ENABLE_UNCERTAINTY` | Master switch for all uncertainty computations |


---

## 7. Key takeaways

1. The bootstrap primarily quantifies **between-location** variability for aggregated patterns — it does not fix model misspecification.
2. **Spatial dependence matters**: block bootstrap and effective-sample-size adjustments are essential when locations are clustered.
3. **Cover-weighted summaries** add a second layer of randomness from bootstrapping the cover index itself. This is appropriate when the index is noisy and unevenly sampled across the season.

