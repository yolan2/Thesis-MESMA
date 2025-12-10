# Test functions for geometric MESMA implementation
# Updated to match the detailed implementation

# Source the functions from the main script
source('fit_veg_mixture_mesma.R', local = TRUE)

test_geometric_mesma <- function() {
  cat("Testing geometric MESMA functions...\n")

  # Test cos_angle
  a <- c(1, 0, 0)
  b <- c(0, 1, 0)
  angle <- cos_angle(a, b)
  expected <- 0  # orthogonal vectors
  if (abs(angle) > 1e-6) stop("cos_angle test failed")

  # Test geometric_select_best_partner
  y <- c(0.3, 0.4, 0.3)
  M1 <- list(list(vec = c(0.1, 0.2, 0.3), id = 'm1'))
  M2 <- list(list(vec = c(0.5, 0.4, 0.3), id = 'm2'))
  res <- geometric_select_best_partner(y, M1, M2)
  if (is.null(res$m1) || is.null(res$m2)) stop("geometric_select_best_partner test failed")

  # Test solve_fclsu (requires quadprog)
  if (requireNamespace("quadprog", quietly = TRUE)) {
    E <- matrix(c(0.1, 0.2, 0.3, 0.5, 0.4, 0.3), nrow = 3)
    y <- c(0.3, 0.4, 0.3)
    res <- solve_fclsu(E, y)
    if (length(res$w) != 2 || any(res$w < 0) || abs(sum(res$w) - 1) > 1e-6) stop("solve_fclsu test failed")
  }

  # Test unmix_stage2_geometric
  veg_libs <- list(
    veg1 = list(list(vec = c(0.1, 0.2, 0.3), id = 'v1')),
    veg2 = list(list(vec = c(0.5, 0.4, 0.3), id = 'v2'))
  )
  res <- unmix_stage2_geometric(y, veg_libs)
  if (is.null(res)) stop("unmix_stage2_geometric test failed")

  # Test geometric_stage1_with_selection
  barren <- c(0.1, 0.1, 0.1)
  veg <- c(0.5, 0.5, 0.5)
  y_stage1 <- 0.3 * barren + 0.7 * veg
  res <- geometric_stage1_with_selection(y_stage1, barren, veg)
  if (abs(res$veg_frac - 0.7) > 0.1) stop("geometric_stage1_with_selection test failed")

  cat("All geometric MESMA tests passed!\n")
}

test_hierarchical_geometric_mesma <- function() {
  cat("Testing hierarchical geometric MESMA...\n")

  # Create test data
  barren <- c(0.1, 0.1, 0.1)
  veg <- c(0.5, 0.5, 0.5)
  y <- 0.2 * barren + 0.8 * veg

  veg_libs <- list(
    grass = list(list(vec = c(0.4, 0.5, 0.4), id = 'g1'), list(vec = c(0.6, 0.4, 0.5), id = 'g2')),
    shrub = list(list(vec = c(0.3, 0.6, 0.3), id = 's1'))
  )

  # Test hierarchical function
  res <- hierarchical_geometric_mesma(y, barren, veg, veg_libs)
  if (is.null(res$fractions)) stop("hierarchical_geometric_mesma test failed")

  # Check that fractions sum to 1
  total <- sum(res$fractions)
  if (abs(total - 1) > 1e-6) stop("Fractions do not sum to 1")

  cat("Hierarchical geometric MESMA test passed!\n")
}

# Run tests
test_geometric_mesma()
test_hierarchical_geometric_mesma()

cat("All tests completed successfully!\n")