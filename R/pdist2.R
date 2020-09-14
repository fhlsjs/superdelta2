########## A fast implementation of superdelta2 ##########
## Aux function that computes pairwise distance^2
pdist2 <- function(A) {
  m <- nrow(A)
  norm2 <- apply(A, 1, function(rvec) crossprod(rvec,rvec))
  ## tmp is the norm2 part, in matrix form
  tmp <- matrix(rep(norm2, m), nrow = m)
  tmp2 <- tmp +  t(tmp)
  return(tmp2 - 2 * tcrossprod(A))
}

