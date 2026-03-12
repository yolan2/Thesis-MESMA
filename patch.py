import sys

with open("fit_veg_mixture_mesma.R", "r", encoding="utf-8") as f:
    text = f.read()

# 1. run_fcls
old_run_fcls = """    run_fcls <- function(E_fit, Y_fit, wts) {
      # delta: match spectral magnitude for sum-to-one enforcement
      delta <- sqrt(mean(E_fit^2)) * 100
      if (!is.finite(delta) || delta < 1e-8) delta <- 1.0
      aug_row <- delta * rep(1, n_endmembers)
      E_aug <- rbind(E_fit, aug_row)
      w_out <- matrix(0, nrow = n_samples, ncol = n_endmembers)
      for (i in seq_len(n_samples)) {
        y_aug <- c(Y_fit[, i], delta)
        res   <- nnls::nnls(E_aug, y_aug)
        w     <- res$x
        w[!is.finite(w)] <- 0
        s <- sum(w)
        if (s > 0) w <- w / s else w <- rep(1 / n_endmembers, n_endmembers)
        w_out[i, ] <- w
      }
      w_out
    }"""

new_run_fcls = """    run_fcls <- function(E_fit, Y_fit, wts) {
      # delta: match spectral magnitude for sum-to-one enforcement
      delta <- sqrt(mean(E_fit^2, na.rm=TRUE)) * 100
      if (!is.finite(delta) || delta < 1e-8) delta <- 1.0
      aug_row <- delta * rep(1, n_endmembers)
      
      w_out <- matrix(0, nrow = n_samples, ncol = n_endmembers)
      
      for (i in seq_len(n_samples)) {
        y_i <- Y_fit[, i]
        
        # Determine valid features (non-NA in both y and E)
        valid_y <- is.finite(y_i)
        valid_E <- apply(E_fit, 1, function(x) all(is.finite(x)))
        valid <- valid_y & valid_E
        
        E_i <- E_fit
        E_i[!is.finite(E_i)] <- 0 # Impute any residual before masking so nnls doesn't crash on NAs in unmasked parts
        
        if (!all(valid)) {
          y_i[!valid] <- 0
          E_i[!valid, ] <- 0
        }
        
        if (i == 1 && n_samples > 0) {
            cat(sprintf("[DEBUG batch_fcls] sum(valid)=%d / %d. NAs: y=%d, E=%d\\n", sum(valid), length(valid), sum(!valid_y), sum(!valid_E)), file="debug_log.txt", append=TRUE)
        }
        
        E_aug <- rbind(E_i, aug_row)
        y_aug <- c(y_i, delta)
        
        res   <- nnls::nnls(E_aug, y_aug)
        w     <- res$x
        w[!is.finite(w)] <- 0
        s <- sum(w)
        if (s > 0) w <- w / s else w <- rep(1 / n_endmembers, n_endmembers)
        w_out[i, ] <- w
      }
      w_out
    }"""

text = text.replace(old_run_fcls, new_run_fcls)

# 2. solve_weights_fcls
old_solve_weights = """  solve_weights_fcls <- function(E, y, feature_weights = NULL, max_iter = 500, tol = 1e-8) {
    if (is.null(E) || ncol(E) < 1) return(NULL)

    # Ensure numeric
    E_fit <- as.matrix(E)
    y_fit <- as.numeric(y)
    n_endmembers <- ncol(E_fit)
    n_bands <- nrow(E_fit)

    # Handle length mismatch
    if (length(y_fit) != n_bands) {
      if (length(y_fit) > n_bands) y_fit <- y_fit[1:n_bands] else y_fit <- c(y_fit, rep(0, n_bands - length(y_fit)))
    }

    # Impute non-finite values
    y_fit[!is.finite(y_fit)] <- 0
    E_fit[!is.finite(E_fit)] <- 0"""

new_solve_weights = """  solve_weights_fcls <- function(E, y, feature_weights = NULL, max_iter = 500, tol = 1e-8) {
    if (is.null(E) || ncol(E) < 1) return(NULL)

    # Ensure numeric
    E_fit <- as.matrix(E)
    y_fit <- as.numeric(y)
    n_endmembers <- ncol(E_fit)
    n_bands <- nrow(E_fit)

    # Handle length mismatch
    if (length(y_fit) != n_bands) {
      if (length(y_fit) > n_bands) y_fit <- y_fit[1:n_bands] else y_fit <- c(y_fit, rep(0, n_bands - length(y_fit)))
    }

    # Mask missing values instead of just zero-padding
    valid_y <- is.finite(y_fit)
    valid_E <- apply(E_fit, 1, function(x) all(is.finite(x)))
    valid <- valid_y & valid_E
    
    cat(sprintf("[DEBUG solve_weights] sum(valid)=%d / %d. NAs: y=%d, E=%d\\n", sum(valid), length(valid), sum(!valid_y), sum(!valid_E)), file="debug_log.txt", append=TRUE)

    E_fit[!is.finite(E_fit)] <- 0
    if (!all(valid)) {
      y_fit[!valid] <- 0
      E_fit[!valid, ] <- 0
    }
    y_fit[!is.finite(y_fit)] <- 0"""

text = text.replace(old_solve_weights, new_solve_weights)

# Save
with open("fit_veg_mixture_mesma.R", "w", encoding="utf-8") as f:
    f.write(text)

print("Patch applied to fit_veg_mixture_mesma.R")
