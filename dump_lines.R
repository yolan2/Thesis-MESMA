args <- commandArgs(trailingOnly=TRUE)
file <- args[1]
start <- as.integer(args[2])
end <- as.integer(args[3])
lines <- readLines(file)
for(i in seq(start, end)){
  cat(sprintf("%5d: %s\n", i, lines[i]))
}
