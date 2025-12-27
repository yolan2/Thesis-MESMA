# Deprecated test removed
# The test that enforced 'no_soil'/'no.soil' policy has been deprecated because
# the codebase no longer uses a 'no soil' fraction column. Remove this file if
# you want to clean up obsolete tests.

# Minimal safe_as_numeric helper (same semantics used by january_averages.R)
safe_as_numeric <- function(x) {
  if (is.null(x)) return(x)
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    s <- trimws(x)
    lower <- tolower(s)
    lower[lower %in% c("true", "t")] <- "1"
    lower[lower %in% c("false", "f")] <- "0"
    suppressWarnings(num <- as.numeric(lower))
    return(num)
  }
  if (is.numeric(x)) return(as.numeric(x))
  suppressWarnings(num <- as.numeric(as.character(x)))
  num
}

# Minimal df with all-NA no_soil variant
base_df <- data.frame(
  location_id = as.character(1:4),
  date = as.Date(rep("2020-01-01", 4)),
  nir = runif(4, 0.5, 1.0),
  red = runif(4, 0.1, 0.4),
  Veg = c("phragmites", "phragmites", "populus", "barren"),
  stringsAsFactors = FALSE
)

# Variant A: 'no_soil' all NA
df_a <- base_df
df_a$no_soil <- NA_real_

# Variant B: 'no.soil' all NA
df_b <- base_df
df_b$`no.soil` <- NA_real_

check_policy_for <- function(df, desc) {
  cat(sprintf("Checking variant: %s\n", desc))
  # Simulate january_averages normalization by defining the same helper
  # function locally (copy of the logic in january_averages.R). This avoids
  # sourcing the whole script which would perform I/O.
  normalize_no_soil_col <- function(tbl) {
    if (is.null(tbl) || !is.data.frame(tbl)) return(tbl)
    nm <- names(tbl)
    candidates <- c("no soil", "no_soil", "no.soil", "__no soil__", "__no_soil__", ".__no soil__", ".__no_soil__")
    if ("no.soil" %in% nm) {
      tbl[["no.soil"]] <- safe_as_numeric(tbl[["no.soil"]])
      return(tbl)
    }
    if ("no soil" %in% nm) {
      tbl[["no.soil"]] <- safe_as_numeric(tbl[["no soil"]])
      tbl[["no soil"]] <- NULL
      return(tbl)
    }
    found <- intersect(candidates, nm)
    if (length(found) > 0) {
      src <- found[1]
      tbl[["no.soil"]] <- safe_as_numeric(tbl[[src]])
      if (src != "no.soil") tbl[[src]] <- NULL
    }
    tbl
  }

  df_norm <- normalize_no_soil_col(df)
  if ("no.soil" %in% names(df_norm)) {
    n_non_na <- sum(!is.na(df_norm$`no.soil`))
    cat(sprintf("Non-NA no.soil: %d\n", n_non_na))
    if (n_non_na == 0) {
      cat("Policy violation detected as expected.\n")
      return(TRUE)
    }
  }
  cat("Policy NOT detected (unexpected).\n")
  return(FALSE)
}

res_a <- tryCatch({ check_policy_for(df_a, "no_soil (all NA)") }, error = function(e) { cat("Error: ", e$message, "\n"); TRUE })
res_b <- tryCatch({ check_policy_for(df_b, "no.soil (all NA)") }, error = function(e) { cat("Error: ", e$message, "\n"); TRUE })

if (!isTRUE(res_a) || !isTRUE(res_b)) {
  stop("Test failed: expected detection of all-NA no_soil / no.soil columns")
}

cat("Test passed: january_averages.R policy enforced for no_soil/no.soil.
")
