import sys

with open("fit_veg_mixture_mesma.R", "r", encoding="utf-8") as f:
    text = f.read()

old_mask = """        # Mask the observation (z-scored only, no PCA-LDA weighting applied to data)
        y_norm_masked <- y_norm_full_pca_lda[valid_mask]

        # Extract masked PCA-LDA weights to pass to solver
        # The solver will apply weights directly to both endmembers and observations
        if (!is.null(weights_for_mask) && length(weights_for_mask) == length(y_norm_full_pca_lda)) {
          weights_masked <- weights_for_mask[valid_mask]
        } else {
          weights_masked <- rep(1, length(y_norm_masked))
        }"""

new_mask = """        # PRESERVE FULL VECTOR so features align with unmasked library endmembers!
        y_norm_masked <- y_norm_full_pca_lda

        # Pass full weights to solver
        if (!is.null(weights_for_mask) && length(weights_for_mask) == length(y_norm_full_pca_lda)) {
          weights_masked <- weights_for_mask
        } else {
          weights_masked <- rep(1, length(y_norm_masked))
        }"""

text = text.replace(old_mask, new_mask)

with open("fit_veg_mixture_mesma.R", "w", encoding="utf-8") as f:
    f.write(text)

print("Patch 3 applied")
