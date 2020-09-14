########## A fast implementation of superdelta2 ##########
#library(superdelta)
#library(genefilter)
## vector/matrix-based Newton-Raphson method.

## matrix-based approx. chisq distribution function
FastFChi2 <- function(x, df, R2max, ngrid=200){
  m <- nrow(x)
  dx <- R2max/ngrid; chi2grid <- seq(0, R2max, dx)
  chi2vec <- pchisq(chi2grid, df)
  ## use rule=2 for cases in which x is slightly larger than R2max due
  ## to numerical errors.
  return(matrix(approx(chi2grid, chi2vec, xout=x, rule=2)[["y"]], m))
}

Newton <- function(mu.trim, R, K, ngrid=200, tol=1e-3, niter=3){
  m <- length(mu.trim)
  k <- 0; err <- Inf; mu.now=mu.trim
  dx <- 2*R/ngrid; mygrid <- t(sapply(1:m, function(i) seq(-R[i], R[i], dx[i])))
  ## Calculate the Chisq-distribution just once
  Fchi2 <- FastFChi2(-sweep(mygrid^2, 1, R^2), df=K-1, R2max=max(R)^2)
  ## xFchi2, x2Fch2 are just the pre-calculated intermediate terms
  xFchi2 <- mygrid*Fchi2; x2Fchi2 <- mygrid*xFchi2; Bmu <- 1
  for (k in 1:niter){
    phi <- dnorm(sweep(mygrid, 1, mu.now))
    ## dx gets cancelled; so no need to multiply each term by dx
    Bmu <- rowSums(Fchi2*phi); Tmu <- rowSums(xFchi2*phi)
    Hmu <- rowSums(x2Fchi2*phi); Smu <- Tmu/Bmu
    dS <- (Hmu*Bmu - Tmu^2)/Bmu^2
    mu.next <- mu.now - (Smu-mu.trim)/dS
    mu.now <- mu.next; err <- mu.next-mu.now
  }
  return(cbind(mu=mu.now, err=err))
}

SphereTrim <- function(Rlist, CenterMat, trim, sigma, prop = 1.0, bias.correct = FALSE,  ...){
  m <- nrow(Rlist[[1]]); G <- length(Rlist)
  ## re-center
  Rlist2 <- lapply(1:G, function(g) sweep(Rlist[[g]], 1, CenterMat[,g]))
  if (prop < 1){
    idx1 <- sample(m, round(m*prop))
    trimmed.m <- round(length(idx1)*(1-trim), 0)
    Lx2 <- Reduce("+", Map(function(x) x[, idx1]^2, Rlist2))
    threshs <- rowOrderStats(Lx2, which = trimmed.m + 1)
  } else {
    trimmed.m <- round((m-1)*(1-trim), 0)
    Lx2 <- Reduce("+", Map(function(x) x^2, Rlist2))
    threshs <- rowOrderStats(Lx2, which = trimmed.m + 2)
    idx1 <- 1:m
  }
  IDX <- sweep(Lx2, 1, threshs, "<")
  Rtrim1 <- do.call(cbind, Map(function(x) {
    x2 <- x[,idx1]; rowSums(x2*IDX)
  }, Rlist2))/trimmed.m
  if (bias.correct) {
    ## bias correction
    mu2 <- sqrt(rowSums(Rtrim1^2))/sigma
    R.normalized <- sweep(Rtrim1, 1, mu2, "/")
    RR <- sqrt(threshs)/sigma
    ## BiasCorrection.results <- t(sapply(1:m, function(i) Newton1(mu2[i], RR[i], K=G-1, ...)))
    BiasCorrection.results <- Newton(mu2, RR, K = G-1, ...)
    R.est1 <- sweep(R.normalized, 1, BiasCorrection.results[,1], "*")
  } else {
    R.est1 <- Rtrim1
  }
  ## Adding the centers back. Rtrim1 is the trimmed estimator without
  ## out bias correction.
  return(list(Rtrim = R.est1 + CenterMat, Rtrim1 = Rtrim1 + CenterMat, threshs = threshs))
}

