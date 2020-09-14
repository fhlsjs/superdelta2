normalize.global <- function(X){
  sweep(X, 2, colMeans(X))
}