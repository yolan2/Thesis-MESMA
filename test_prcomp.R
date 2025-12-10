
X <- matrix(rnorm(100), ncol=10)
X[, 1] <- 0 # Constant column
tryCatch({
  prcomp(X, center = FALSE, scale. = FALSE)
  print("Success with scale=FALSE")
}, error = function(e) {
  print(paste("Error with scale=FALSE:", e$message))
})

tryCatch({
  prcomp(X, center = FALSE, scale. = TRUE)
  print("Success with scale=TRUE")
}, error = function(e) {
  print(paste("Error with scale=TRUE:", e$message))
})
