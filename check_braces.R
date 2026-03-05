lines <- readLines("january_averages.R")
cnt <- 0
for(i in seq_along(lines)){
  opens <- lengths(gregexpr("{", lines[i], fixed=TRUE))
  closes <- lengths(gregexpr("}", lines[i], fixed=TRUE))
  cnt <- cnt + opens - closes
  if(cnt < 0) cat("index", i, "negative count", cnt, "\n")
}
cat("final count", cnt, "\n")
