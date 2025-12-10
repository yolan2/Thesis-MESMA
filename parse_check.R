tryCatch({
  parse(file = "c:/Users/yolan/OneDrive/Documenten/UGENT/Master/R_MESMA/fit_veg_mixture_mesma.R")
  cat("PARSE_OK\n")
}, error = function(e) {
  cat("PARSE_ERROR:", e$message, "\n")
})
