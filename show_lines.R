lines <- readLines('january_averages.R')
for(i in seq_along(lines)) {
  if (i >= 1760 && i <= 1780) cat(i, ':', lines[i], '\n')
}
