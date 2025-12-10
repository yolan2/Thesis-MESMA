path <- 'c:/Users/yolan/OneDrive/Documenten/UGENT/Master/R_MESMA/fit_veg_mixture_mesma.R'
s <- readLines(path, warn=FALSE)
cnt <- 0
maxcnt <- 0
last_zero_line <- 0
for(i in seq_along(s)){
  line <- s[i]
  # Ignore full-line comments
  if (grepl('^\\s*#', line, perl = TRUE)) next
  nopen <- gregexpr('\\{', line)[[1]]
  nopen <- ifelse(nopen[1] == -1, 0, length(nopen))
  nclose <- gregexpr('\\}', line)[[1]]
  nclose <- ifelse(nclose[1] == -1, 0, length(nclose))
  delta <- nopen - nclose
  prev_cnt <- cnt
  cnt <- cnt + delta
  if(delta != 0){
    cat(sprintf("%5d delta=%3d cnt=%3d %s\n", i, delta, cnt, substr(line,1,140)))
  }
  if(cnt > maxcnt){ maxcnt <- cnt; cat(sprintf("NEW_MAX %d at line %d\n", maxcnt, i)) }
  if(cnt == 0) last_zero_line <- i
}
cat('FINAL COUNT', cnt, '\n')
cat('Last zero at line', last_zero_line, '\n')
cat('Lines after last zero (first 200 shown):\n')
for (i in (last_zero_line+1):min(length(s), last_zero_line + 200)) cat(sprintf("%5d %s\n", i, substr(s[i], 1, 150)))

cat('\nAnalyzing just the region after last zero for unmatched braces...\n')
cnt2 <- 0
for (i in seq(last_zero_line+1, length(s))) {
  line <- s[i]
  nopen <- gregexpr('\\{', line)[[1]]
  nopen <- ifelse(nopen[1] == -1, 0, length(nopen))
  nclose <- gregexpr('\\}', line)[[1]]
  nclose <- ifelse(nclose[1] == -1, 0, length(nclose))
  delta <- nopen - nclose
  cnt2 <- cnt2 + delta
  if(delta != 0) cat(sprintf("%5d delta=%3d cnt=%3d %s\n", i, delta, cnt2, substr(line,1,140)))
}
cat('FINAL COUNT REGION', cnt2, '\n')

# Build stack of brace positions to find unmatched openings
stack <- list()
for (i in seq(last_zero_line+1, length(s))) {
  line <- s[i]
  # find all { and } positions in order
  opens <- gregexpr('\\{', line)[[1]]
  closes <- gregexpr('\\}', line)[[1]]
  pos <- c()
  type <- c()
  if (opens[1] != -1) { pos <- c(pos, opens); type <- c(type, rep('{', length(opens))) }
  if (closes[1] != -1) { pos <- c(pos, closes); type <- c(type, rep('}', length(closes))) }
  if (length(pos) > 0) {
    ord <- order(pos)
    for (k in ord) {
      if (type[k] == '{') {
        stack <- c(stack, list(list(line=i, text=substr(line,1,140))))
      } else {
        if (length(stack) > 0) stack <- stack[-length(stack)]
        else {
          cat(sprintf('Unmatched closing brace at %d\n', i))
        }
      }
    }
  }
}
cat('Unmatched openings remaining:', length(stack), '\n')
if (length(stack) > 0) {
  for (j in seq_along(stack)) cat(sprintf('Unclosed { at line %d: %s\n', stack[[j]]$line, stack[[j]]$text))
}
