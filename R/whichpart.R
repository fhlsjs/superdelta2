## aux function that returns a list of the first p smallest x, which uses partial sorting.
## Note: In very rare cases, due to ties, sum(IDX[, i]) may be > 800.
whichpart <- function(x, p) {
  xp <- sort(x, partial=p)[p]
  return(x <= xp)
}