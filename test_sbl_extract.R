library(stringr)
file <- 'c:/Users/yolan/OneDrive/Documenten/UGENT/Master/R_MESMA/fit_veg_mixture_mesma.R'
lines <- readLines(file)
start_line <- which(grepl('^\\s*solve_weights_ols\\s*<-\\s*function\\(', lines))
if (length(start_line) != 1) stop('ERROR: cannot find function start')
line_idx <- start_line
brace_count <- 0
found_start <- FALSE
end_line <- NULL
for (i in seq(start_line, length(lines))) {
  l <- lines[i]
  l_no_comment <- sub('#.*$', '', l)
  brace_count <- brace_count + stringr::str_count(l_no_comment, '\\{') - stringr::str_count(l_no_comment, '\\}')
  if (!found_start && grepl('\\{', l_no_comment)) found_start <- TRUE
  if (found_start && brace_count == 0) { end_line <- i; break }
}
if (is.null(end_line)) stop('ERROR: cannot find function end')
func_text <- paste(lines[start_line:end_line], collapse='\n')
# Evaluate function in new environment
env <- new.env()
eval(parse(text=func_text), envir = env)
if (!exists('solve_weights_ols', envir = env)) stop('ERROR: function not created')
# Test function
set.seed(1)
E <- matrix(c(1,0.2,0.2, 0.1,1,0.3, 0.2,0.1,1), nrow=3, ncol=3)
y <- as.numeric(E %*% c(0.5, 0.3, 0.2) + rnorm(3, 0, 1e-3))
cat('Testing SBL: dims E=', paste(dim(E), collapse='x'), ', length y=', length(y), ', var(y)=', var(y), '\n')
res <- tryCatch(env$solve_weights_sparse_bayesian(E, y, max_iter = 100, tol = 1e-8, enforce_constraints = TRUE, verbose = FALSE), error = function(e) {
  cat('SBL test raised error: ', e$message, '\n')
  NULL
})
if (!is.null(res)) {
  cat('SBL test result OK: w=', paste(round(res$w,4), collapse=', '), '\n')
} else {
  stop('SBL test failed')
}
