
source('mesma_config.R')
source('mesma_helpers.R')
source('ppi_helpers.r')

B <- 10
# mock function
bootstrap_hierarchical_means <- function(df, metrics = c('MSAVI', 'NDVI', 'PPI'), group_col = 'location_id', B = 1000) {
    if (!group_col %in% names(df)) return(rep(NA, length(metrics) * 2))

    cols <- unique(c(group_col, metrics, 'lat', 'lon'))
    df <- df |> dplyr::select(dplyr::any_of(cols))

    df[[group_col]] <- as.character(df[[group_col]])
    ids <- unique(df[[group_col]])
    n_ids <- length(ids)
    if (n_ids < 2) return(rep(NA, length(metrics) * 2))

    coords_df <- NULL
    if (all(c('lat','lon') %in% names(df))) {
      coords_df <- unique(df[, c(group_col, 'lat', 'lon')])
      names(coords_df)[names(coords_df) == group_col] <- 'location_id'
    }

    data_map <- lapply(ids, function(id) {
      mat <- as.matrix(df[df[[group_col]] == id, metrics, drop = FALSE])        
      if (is.null(dim(mat))) mat <- matrix(mat, ncol = length(metrics))
      colnames(mat) <- metrics
      mat
    })
    names(data_map) <- ids

    boot_replicates <- replicate(B, {
      sel_ids <- sample(ids, n_ids, replace = TRUE)
      selected_mats <- data_map[sel_ids]

      cluster_means <- vapply(selected_mats, function(mat) {
        n_obs <- nrow(mat)
        if (is.na(n_obs) || n_obs <= 0) return(rep(NA_real_, length(metrics)))  
        rows <- sample.int(n_obs, n_obs, replace = TRUE)
        colMeans(mat[rows, , drop = FALSE], na.rm = TRUE)
      }, numeric(length(metrics)))
      if (is.null(dim(cluster_means))) cluster_means <- matrix(cluster_means, nrow = length(metrics))
      rowMeans(cluster_means, na.rm = TRUE)
    })

    if (is.null(dim(boot_replicates))) {
      boot_replicates <- matrix(boot_replicates, nrow = length(metrics))        
    }
    cis <- apply(boot_replicates, 1, function(x) stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))

    output <- as.vector(cis)
    return(output)
}

df <- data.frame(
  location_id = rep(c('a','b','c'), each=5),
  MSAVI_norm = runif(15)
)
print(bootstrap_hierarchical_means(df, metrics='MSAVI_norm', B=10))


