# Test build_mesma_variants with synthetic data
src <- 'c:/Users/yolan/OneDrive/Documenten/UGENT/Master/R_MESMA/fit_veg_mixture_mesma.R'
source(src)

# Create synthetic reduced_data
n_traces <- 20
n_time <- 365
K <- 3
Z_sample <- matrix(0, nrow = n_time, ncol = K)
set.seed(1)

Z_list <- lapply(1:n_traces, function(i) {
  mat <- matrix(runif(n_time * K, -0.1, 0.9), nrow = n_time, ncol = K)
  colnames(mat) <- paste0('idx', 1:K)
  mat
})

# Build X_feat flattened
X_feat <- do.call(rbind, lapply(Z_list, function(m) as.numeric(m)))

reduced_data <- list(testveg = list(features = X_feat, Z_matrices = Z_list, trace_info = lapply(1:n_traces, function(i) list(location_id = i, year = 2020)), n_samples = n_traces))

raw_lib_templates <- list(testveg = list(T = Z_list[[1]]))

res <- build_mesma_variants(reduced_data, raw_lib_templates, min_cluster_size = 2)
print(str(res))
