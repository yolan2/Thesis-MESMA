# =============================================================================
# plot_pca_lda_scatter.R
# -----------------------------------------------------------------------------
# Produces a scatter plot of training endmembers (one point per loc-year)
# projected into PCA-LDA space, coloured by vegetation class.
#
# The pipeline mirrors the logic inside train_feature_pipeline():
#   1. Load preprocessed data and filter to TRAIN_YEARS
#   2. Build a temporal feature vector per loc-year (pentad averages, L2-norm,
#      z-score) — the same representation used for training
#   3. PCA (retain 95 % variance)  →  LDA on PCA scores
#   4. Project each loc-year sample onto LD1/LD2 and plot
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(MASS)
  library(nnls)
})

# ── Configuration constants ───────────────────────────────────────────────────
source("mesma_config.R")

# ── Data loading ──────────────────────────────────────────────────────────────
cat("Loading preprocessed data...\n")
df          <- readRDS("preprocessed_data.rds")
norm_params <- readRDS("training_norm_params.rds")
INDEX_SCALES <- norm_params$INDEX_SCALES

# Feature columns: prefer training norm indices, fall back to OPTIMAL_INDICES
if (!is.null(INDEX_SCALES) && length(INDEX_SCALES) > 0) {
  avail <- intersect(names(INDEX_SCALES), names(df))
} else {
  avail <- intersect(OPTIMAL_INDICES, names(df))
}
if (length(avail) == 0) stop("No feature indices found in preprocessed data.")
cat(sprintf("Features (%d): %s\n", length(avail), paste(avail, collapse = ", ")))

# ── Data filtering ────────────────────────────────────────────────────────────
# Compute phenological year if not present
if (!"pheno_year" %in% names(df)) {
  df$date <- as.Date(df$date)
  yr       <- as.integer(format(df$date, "%Y"))
  march1   <- as.Date(paste0(yr, "-03-01"))
  df$pheno_year <- ifelse(df$date >= march1, yr, yr - 1L)
}

df <- df[df$pheno_year %in% TRAIN_YEARS,  , drop = FALSE]
df <- df[!is.na(df$Veg) & nzchar(trimws(df$Veg)), , drop = FALSE]
df$target_class <- tolower(trimws(df$Veg))

allowed_all <- unique(c(tolower(ALLOWED_VEG), "barren"))
df <- df[df$target_class %in% allowed_all, , drop = FALSE]

if (nrow(df) == 0) stop("No training rows remain after filtering.")
cat(sprintf("Training rows: %d  |  classes: %s\n",
            nrow(df), paste(sort(unique(df$target_class)), collapse = ", ")))

# ── Phenological DOY ──────────────────────────────────────────────────────────
if (!"doy" %in% names(df) || any(is.na(df$doy))) {
  df$date  <- as.Date(df$date)
  yr       <- as.integer(format(df$date, "%Y"))
  march1   <- as.Date(paste0(yr, "-03-01"))
  pheno_s  <- ifelse(df$date >= march1,
                     as.integer(march1),
                     as.integer(as.Date(paste0(yr - 1L, "-03-01"))))
  df$doy   <- as.integer(df$date) - pheno_s + 1L
}

# ── Build loc-year feature vectors ────────────────────────────────────────────
N_PENT  <- TEMPORAL_BUDGET
PDAYS   <- TEMPORAL_AGGREGATION_DAYS
MIN_OBS <- MIN_PENTADS_PER_TRAIN_SAMPLE

df$pentad <- pmin(ceiling(df$doy / PDAYS), N_PENT)

traces <- split(df, list(df$location_id, df$pheno_year), drop = TRUE)
cat(sprintf("Total loc-year traces: %d\n", length(traces)))

X_list     <- list()
y_class    <- character(0)
y_location <- character(0)
y_year     <- character(0)

for (nm in names(traces)) {
  sub <- traces[[nm]]
  if (nrow(sub) < MIN_OBS) next

  # Pentad-averaged matrix  [N_PENT × n_features]
  mat <- matrix(NA_real_, nrow = N_PENT, ncol = length(avail))
  colnames(mat) <- avail

  for (p in seq_len(N_PENT - 1L)) {           # last pentad intentionally left NA
    rows <- sub[sub$pentad == p, , drop = FALSE]
    if (nrow(rows) == 0L) next
    for (j in seq_along(avail)) {
      v <- rows[[avail[j]]]
      v <- v[is.finite(v)]
      if (length(v) > 0L) mat[p, j] <- mean(v)
    }
  }

  # Linear interpolation for missing interior pentads (mirrors INTERPOLATE_INFERENCE="linear")
  for (j in seq_along(avail)) {
    y  <- mat[, j]
    ok <- which(is.finite(y))
    if (length(ok) < 2L) next
    mat[, j] <- approx(ok, y[ok], xout = seq_len(N_PENT), rule = 2L)$y
  }

  vec <- as.numeric(mat)

  # L2-normalize whole vector (ENABLE_LDA_L2_NORMALIZATION = TRUE)
  v_clean <- vec
  v_clean[!is.finite(v_clean)] <- 0
  nrm <- sqrt(sum(v_clean^2))
  if (is.finite(nrm) && nrm > 1e-9) vec <- vec / nrm

  X_list[[length(X_list) + 1L]] <- vec

  cls  <- names(sort(table(sub$target_class), decreasing = TRUE))[1]
  loc  <- as.character(sub$location_id[1])
  yr_s <- as.character(sub$pheno_year[1])

  y_class    <- c(y_class,    cls)
  y_location <- c(y_location, loc)
  y_year     <- c(y_year,     yr_s)
}

n_samples <- length(X_list)
cat(sprintf("Retained loc-year samples: %d\n", n_samples))
if (n_samples < 4L) stop("Too few samples for PCA-LDA.")

X_mat <- do.call(rbind, X_list)

# ── Z-score per feature×pentad column group (ENABLE_ZSCORE_AFTER_L2 = TRUE) ──
# Store per-index means and SDs so the collinearity guard can modify them
z_means <- numeric(length(avail))
z_sds   <- numeric(length(avail))
names(z_means) <- avail
names(z_sds)   <- avail

X_z <- X_mat
for (k in seq_along(avail)) {
  idx <- (k - 1L) * N_PENT + seq_len(N_PENT)
  v   <- X_mat[, idx]
  mu  <- mean(v, na.rm = TRUE)
  sg  <- sd(v,   na.rm = TRUE)
  if (is.na(sg) || sg < 1e-9) sg <- 1
  z_means[k] <- mu
  z_sds[k]   <- sg
  X_z[, idx] <- (v - mu) / sg
}
X_z[!is.finite(X_z)] <- 0

# ── Collinearity guard (mirrors train_feature_pipeline) ──────────────────────
# Detect classes trapped inside the convex hull of others and boost separating
# features by shrinking their z-score SDs.  This keeps the scatter plot
# consistent with the trained MESMA_PARAMS feature space.
coll_guard_enable <- exists("ENABLE_COLLINEARITY_GUARD") && isTRUE(ENABLE_COLLINEARITY_GUARD)

if (coll_guard_enable && length(unique(y_class)) >= 3) {
  coll_threshold <- if (exists("COLLINEARITY_THRESHOLD")) COLLINEARITY_THRESHOLD else 1.5
  coll_boost_max <- if (exists("COLLINEARITY_BOOST_FACTOR")) COLLINEARITY_BOOST_FACTOR else 3.0

  classes_coll <- unique(y_class)
  centroids_coll <- do.call(rbind, lapply(classes_coll, function(cls) {
    colMeans(X_z[y_class == cls, , drop = FALSE], na.rm = TRUE)
  }))
  rownames(centroids_coll) <- classes_coll

  coll_boost <- rep(1.0, length(avail))
  names(coll_boost) <- avail

  for (ci in seq_along(classes_coll)) {
    cls_name <- classes_coll[ci]
    target_centroid <- centroids_coll[ci, ]
    other_centroids <- centroids_coll[-ci, , drop = FALSE]

    E_coll <- t(other_centroids)
    delta_coll <- sqrt(mean(E_coll^2, na.rm = TRUE)) * 100
    if (!is.finite(delta_coll) || delta_coll < 1e-8) delta_coll <- 1.0

    E_aug_coll <- rbind(E_coll, delta_coll * rep(1, ncol(E_coll)))
    y_aug_coll <- c(target_centroid, delta_coll)

    res_coll <- nnls::nnls(E_aug_coll, y_aug_coll)
    alpha_coll <- res_coll$x
    s_coll <- sum(alpha_coll)
    if (s_coll > 0) alpha_coll <- alpha_coll / s_coll

    approx_coll <- as.numeric(E_coll %*% alpha_coll)
    residual_coll <- target_centroid - approx_coll

    sep_by_index <- numeric(length(avail))
    within_sd_by_index <- numeric(length(avail))
    for (k in seq_along(avail)) {
      col_s <- (k - 1L) * N_PENT + 1L
      col_e <- k * N_PENT
      sep_by_index[k] <- sqrt(mean(residual_coll[col_s:col_e]^2))
      cls_vals <- X_z[y_class == cls_name, col_s:col_e, drop = FALSE]
      within_sd_by_index[k] <- mean(apply(cls_vals, 2, sd, na.rm = TRUE), na.rm = TRUE)
    }
    within_sd_by_index[within_sd_by_index < 1e-10] <- 1e-10

    sep_distance <- sqrt(sum(residual_coll^2))
    class_spread <- sqrt(mean(apply(X_z[y_class == cls_name, , drop = FALSE], 2,
                                    var, na.rm = TRUE), na.rm = TRUE))
    if (!is.finite(class_spread) || class_spread < 1e-10) class_spread <- 1e-10
    relative_sep <- sep_distance / class_spread

    cat(sprintf("  [COLLINEARITY] Class '%s': rel_sep=%.4f (threshold=%.2f)\n",
                cls_name, relative_sep, coll_threshold))

    if (relative_sep < coll_threshold) {
      snr_by_index <- sep_by_index / within_sd_by_index
      max_snr <- max(snr_by_index, na.rm = TRUE)
      if (!is.finite(max_snr) || max_snr < 1e-10) max_snr <- 1.0
      snr_norm <- snr_by_index / max_snr
      index_boost <- 1.0 + (coll_boost_max - 1.0) * snr_norm
      coll_boost <- pmax(coll_boost, index_boost)
      cat(sprintf("  [COLLINEARITY] >>> '%s' trapped — boosting separating features\n", cls_name))
    }
  }

  if (any(coll_boost > 1.0 + 1e-6)) {
    cat(sprintf("  [COLLINEARITY] Boosting %d/%d indices (max=%.2f)\n",
                sum(coll_boost > 1.0 + 1e-6), length(avail), max(coll_boost)))
    z_sds <- z_sds / coll_boost
    for (k in seq_along(avail)) {
      idx <- (k - 1L) * N_PENT + seq_len(N_PENT)
      X_z[, idx] <- (X_mat[, idx] - z_means[k]) / z_sds[k]
    }
    X_z[!is.finite(X_z)] <- 0
    cat("  [COLLINEARITY] Z-scores recomputed with boosted SDs\n")
  }
}

# ── PCA ───────────────────────────────────────────────────────────────────────
vars      <- apply(X_z, 2, var)
keep_cols <- vars > 1e-9
X_pca_in  <- X_z[, keep_cols, drop = FALSE]
pca_res   <- prcomp(X_pca_in, center = FALSE, scale. = FALSE)

cum_var <- cumsum(pca_res$sdev^2) / sum(pca_res$sdev^2)
n_pcs   <- which(cum_var > PCA_VARIANCE_THRESHOLD)[1L]
if (is.na(n_pcs)) n_pcs <- ncol(pca_res$x)

# Guard for LDA stability
n_classes       <- length(unique(y_class))
class_counts    <- table(y_class)
n_min           <- min(class_counts)
min_pcs_for_lda <- max(1L, n_classes - 1L)
max_pcs_for_lda <- min(20L, max(min_pcs_for_lda, n_min - 2L))
n_pcs <- max(min_pcs_for_lda, min(n_pcs, max_pcs_for_lda))

cat(sprintf("PCA: using %d PCs (%.1f %% var retained)\n",
            n_pcs, 100 * cum_var[n_pcs]))

# ── LDA on PCA scores ─────────────────────────────────────────────────────────
lda_res <- withCallingHandlers(
  MASS::lda(pca_res$x[, 1:n_pcs, drop = FALSE],
            grouping = as.factor(y_class)),
  warning = function(w) invokeRestart("muffleWarning")
)
n_ld <- ncol(lda_res$scaling)
cat(sprintf("LDA: %d discriminant axis/axes\n", n_ld))

# Proportion of between-class variance explained by each axis
svd_sq  <- lda_res$svd^2
pct_var <- 100 * svd_sq / sum(svd_sq)

# ── Project into LDA space ────────────────────────────────────────────────────
X_ld <- pca_res$x[, 1:n_pcs, drop = FALSE] %*% lda_res$scaling  # [n × n_ld]
colnames(X_ld) <- paste0("LD", seq_len(n_ld))

# ── Plot data frame ───────────────────────────────────────────────────────────
plot_df <- data.frame(
  X_ld,
  class    = factor(y_class),
  location = y_location,
  year     = y_year,
  loc_year = paste(y_location, y_year, sep = "·"),
  stringsAsFactors = FALSE
)

# ── Colour palette ────────────────────────────────────────────────────────────
palette_base <- c(
  "populus"    = "#2ca02c",
  "tamarix"    = "#d62728",
  "phragmites" = "#1f77b4",
  "herbs"      = "#bcbd22",
  "barren"     = "#8c564b"
)
classes_present <- levels(plot_df$class)
extra_classes   <- setdiff(classes_present, names(palette_base))
if (length(extra_classes) > 0L) {
  extra_cols <- grDevices::rainbow(length(extra_classes), s = 0.8, v = 0.85)
  names(extra_cols) <- extra_classes
  palette_base <- c(palette_base, extra_cols)
}
pal <- palette_base[classes_present]   # keep only classes present in data

# ── Axis labels ───────────────────────────────────────────────────────────────
x_col <- "LD1"
y_col <- if (n_ld >= 2L) "LD2" else "LD1"

x_label <- sprintf("LD1  (%.1f %% between-class variance)", pct_var[1])
y_label <- if (n_ld >= 2L) {
  sprintf("LD2  (%.1f %% between-class variance)", pct_var[2])
} else {
  "LD1 (only one discriminant axis)"
}

subtitle_text <- sprintf(
  "1 point = 1 loc-year  |  training years: %s  |  n = %d samples across %d classes",
  paste(TRAIN_YEARS, collapse = ", "), nrow(plot_df), n_classes
)

# ── ggplot ────────────────────────────────────────────────────────────────────
p <- ggplot(plot_df,
            aes(x = LD1,
                y = if (n_ld >= 2L) LD2 else LD1,
                colour = class,
                shape  = year)) +
  geom_point(size = 3.2, alpha = 0.80) +
  {
    if (n_ld >= 2L && nrow(plot_df) >= 4L) {
      stat_ellipse(
        aes(fill = class),
        geom  = "polygon",
        alpha = 0.10,
        level = 0.80,
        show.legend = FALSE
      )
    }
  } +
  scale_colour_manual(values = pal, name = "Class") +
  scale_fill_manual(values = pal, name = "Class") +
  scale_shape_manual(
    values = c("2023" = 16, "2024" = 17, "2025" = 15,
               "2022" = 18, "2021" = 8,  "2020" = 3),
    name = "Year"
  ) +
  labs(
    title    = "Training Endmembers: PCA-LDA",
    subtitle = subtitle_text,
    x        = x_label,
    y        = y_label
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position   = "right",
    plot.subtitle     = element_text(size = 9, colour = "grey40"),
    panel.grid.minor  = element_blank()
  )

# ── Save + display ────────────────────────────────────────────────────────────
out_file <- "training_lda_scatter.png"
ggsave(out_file, p, width = 9, height = 6.5, dpi = 150)
cat(sprintf("Saved: %s\n", out_file))

# Show class distribution per axis
cat("\nSamples per class:\n")
print(sort(table(y_class), decreasing = TRUE))

cat("\nLD1 class centroids:\n")
print(round(tapply(X_ld[, 1], y_class, mean), 3))
if (n_ld >= 2L) {
  cat("LD2 class centroids:\n")
  print(round(tapply(X_ld[, 2], y_class, mean), 3))
}
