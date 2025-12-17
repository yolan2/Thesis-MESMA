lines <- readLines('fit_veg_mixture_mesma.R')
# Remove comments
lines <- sub('#.*', '', lines)
balance <- 0
for (i in 1:length(lines)) {
  line <- lines[i]
  open <- nchar(gsub('[^\\{]', '', line))
  close <- nchar(gsub('[^\\}]', '', line))
  balance <- balance + open - close
  if (balance < 0) {
    cat('Negative brace balance at line ', i, ': ', line, '\n')
    break
  }
}
cat('Final brace balance: ', balance, '\n')

# Check parentheses
balance_p <- 0
for (i in 1:length(lines)) {
  line <- lines[i]
  open_p <- nchar(gsub('[^\\(]', '', line))
  close_p <- nchar(gsub('[^\\)]', '', line))
  balance_p <- balance_p + open_p - close_p
  if (balance_p < 0) {
    cat('Negative paren balance at line ', i, ': ', line, '\n')
    break
  }
}
cat('Final paren balance: ', balance_p, '\n')

# Check quotes
total_quotes <- 0
for (i in 1:length(lines)) {
  line <- lines[i]
  quotes <- nchar(gsub('[^"]', '', line))
  total_quotes <- total_quotes + quotes
}
cat('Total quotes: ', total_quotes, '\n')

if (balance == 0 && balance_p == 0 && total_quotes %% 2 == 0) cat('All balanced\n') else cat('Imbalance detected\n')