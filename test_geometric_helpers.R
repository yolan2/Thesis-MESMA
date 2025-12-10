project_to_simplex <- function(v) { n <- length(v); u <- sort(v, decreasing = TRUE); cssv <- cumsum(u); rho <- max(which(u - (cssv - 1) / seq_along(u) > 0)); theta <- (cssv[rho] - 1) / rho; w <- pmax(v - theta, 0); if (sum(w) <= 0) rep(1 / n, n) else w / sum(w) }

geometric_project_and_unmix <- function(y, m1, m2) {
  em_line <- m2 - m1
  em_norm_sq <- sum(em_line^2)
  if (em_norm_sq < 1e-10) return(list(f1 = 0.5, f2 = 0.5, y_proj = m1, residual = sqrt(sum((y - m1)^2))))
  y_minus_m1 <- y - m1
  t <- sum(y_minus_m1 * em_line) / em_norm_sq
  t_clamped <- max(0, min(1, t))
  y_proj <- m1 + t_clamped * em_line
  residual <- sqrt(sum((y - y_proj)^2))
  f2 <- t_clamped; f1 <- 1 - f2
  list(f1 = f1, f2 = f2, y_proj = y_proj, residual = residual, t = t)
}

geometric_unmix_simplex <- function(y, M) {
  N <- ncol(M)
  if (N == 1) return(list(f = 1, residual = sqrt(sum((y - M[, 1])^2)), y_proj = M[, 1]))
  if (N == 2) {
    res <- geometric_project_and_unmix(y, M[,1], M[,2])
    list(f = c(res$f1, res$f2), residual = res$residual, y_proj = res$y_proj)
  } else {
    MtM <- t(M) %*% M; Mty <- t(M) %*% y; ridge <- 1e-8 * diag(N)
    f_unconstrained <- tryCatch(solve(MtM + ridge, Mty), error = function(e) { if (requireNamespace('MASS', quietly=T)) MASS::ginv(MtM + ridge) %*% Mty else rep(1/N, N) })
    f <- project_to_simplex(as.numeric(f_unconstrained))
    y_proj <- as.numeric(M %*% f); residual <- sqrt(sum((y - y_proj)^2))
    list(f = f, residual = residual, y_proj = y_proj)
  }
}

spectral_angle <- function(a, b) { cs <- sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2))); acos(pmax(-1, pmin(1, cs))) }

geometric_select_pair <- function(y, M1, M2) {
  best_dist <- Inf; best_m1 <- NULL; best_m2 <- NULL
  for (i in seq_along(M1)) {
    m1 <- M1[[i]]$vec
    best_angle <- Inf; best_idx <- NULL
    for (j in seq_along(M2)) {
      m2 <- M2[[j]]$vec
      diff <- y - m1; diff2 <- m1 - m2; angle <- spectral_angle(diff, diff2)
      if (angle < best_angle) { best_angle <- angle; best_idx <- j }
    }
    if (!is.null(best_idx)) {
      m2 <- M2[[best_idx]]$vec
      E <- cbind(m1, m2); w <- tryCatch(solve(E, y), error = function(e) rep(0,2)); w <- pmax(0, w)
      if (sum(w) > 0) w <- w / sum(w)
      dist <- sqrt(sum((y - E %*% w)^2))
      if (dist < best_dist) { best_dist <- dist; best_m1 <- M1[[i]]; best_m2 <- M2[[best_idx]] }
    }
  }
  list(m1 = best_m1, m2 = best_m2, dist = best_dist)
}

angle_based_mesma <- function(y, library_list) {
  veg_names <- names(library_list)
  n_libs <- length(veg_names)
  if (n_libs == 0) return(NULL)
  if (n_libs == 1) {
    lib <- library_list[[1]]
    best_dist <- Inf
    best_idx <- 1
    for (i in seq_along(lib)) {
      dist <- sqrt(sum((y - lib[[i]]$vec)^2))
      if (dist < best_dist) { best_dist <- dist; best_idx <- i }
    }
    return(list(fractions = setNames(1.0, veg_names[1]), chosen = setNames(lib[[best_idx]]$id, veg_names[1]), residual = best_dist))
  }
  if (n_libs == 2) {
    r <- geometric_select_pair(y, library_list[[1]], library_list[[2]])
    if (is.null(r$m1) || is.null(r$m2)) return(NULL)
    u <- geometric_project_and_unmix(y, r$m1$vec, r$m2$vec)
    fractions <- c(u$f1, u$f2)
    names(fractions) <- veg_names
    chosen <- c(r$m1$id, r$m2$id); names(chosen) <- veg_names
    return(list(fractions = fractions, chosen = chosen, residual = u$residual))
  }
  NULL
}

# Tests
m1 <- c(0.1,0.2,0.3,0.4,0.5)
m2 <- c(0.5,0.4,0.3,0.2,0.1)
y <- 0.6*m1 + 0.4*m2
res <- geometric_project_and_unmix(y, m1, m2)
cat(sprintf('Project unmix result: f1=%.4f f2=%.4f residual=%.8f\n', res$f1, res$f2, res$residual))

if (abs(res$f1 - 0.6) > 1e-6) stop('geometric_project_and_unmix test failed')

libs <- list(veg1 = list(list(vec = m1, id = 'v1')), veg2 = list(list(vec = m2, id = 'v2')))
ans <- angle_based_mesma(y, libs)
cat('Angle-based result fractions:', paste(round(ans$fractions,4), collapse=', '), '\n')
if (abs(ans$fractions['veg1'] - 0.6) > 0.05) stop('angle_based_mesma test failed')

cat('All minimal tests passed\n')
