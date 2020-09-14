## functions related to bias correction
## Sfunc <- function(R, mu, K, ngrid=1000){
##   n <- length(R);   Rmax <- max(R)
##   dx <- 2*Rmax/ngrid; mygrid <- seq(-Rmax, Rmax, dx)
##   x.minus.mu <- sapply(1:n, function(i) mygrid-mu[i])
##   Fchi2.grid <- pchisq(mygrid, K)
##   phi.grid <- dnorm(x.minus.mu)
##   h.grid <- sweep(phi.grid, 1, Fchi2.grid, "*")
##   S.top <- colSums(sweep(h.grid, 1, mygrid, "*"))
##   S.bottom <- 
## }
## K <- 3; R <- 2:4; mu <- round(runif(3), 2)


Sfunc1 <- function(mu, R, K, ngrid = 200){
  dx <- 2*R/ngrid; mygrid <- seq(-R, R, dx)
  phi.grid <- dnorm(mygrid-mu)
  Fchi2.grid <- pchisq(R^2-mygrid^2, K-1)
  h.grid <- Fchi2.grid*phi.grid
  Bmu <- sum(h.grid)*dx
  Tmu <- sum(h.grid*mygrid)*dx
  Hmu <- sum(h.grid*(mygrid^2))*dx
  return(c(Tmu = Tmu, Bmu = Bmu, Hmu = Hmu, S = Tmu/Bmu, dS = (Hmu*Bmu - Tmu^2)/Bmu^2))
}

Sderiv1 <- function(mu, R, K, ngrid = 200){
  dx <- 2*R/ngrid; mygrid <- seq(-R, R, dx)
  phi.grid <- dnorm(mygrid-mu)
  Fchi2.grid <- pchisq(R^2-mygrid^2, K-1)
  h.grid <- Fchi2.grid*phi.grid
  Bmu <- sum(h.grid)*dx
  Tmu <- sum(h.grid*mygrid)*dx
  Hmu <- sum(h.grid*(mygrid^2))*dx
  return((Hmu*Bmu - Tmu^2)/Bmu^2)
}

## ## This is the univariate version of the inverse of S, using
## ## Newton-Raphson.
## Newton1a <- function(mu.trim, R, K, ngrid = 200, tol = 1e-3, max.iter = 10){
##   k <- 0; err <- Inf; mu.now = mu.trim
##   ## While loop with implementation condition
##   while ((abs(err) > tol) & (k < max.iter)) {
##     rk <- Sfunc1(mu.now, R, K, ngrid = ngrid)
##     ## The update step
##     mu.next <- mu.now - (rk[["S"]] - mu.trim)/rk[["dS"]]
##     k <- k+1; err <- mu.next - mu.now
##     mu.now <- mu.next
##   }
##   if (abs(err) > tol) warning("Maximum number of iteration reached but the error is still greater than the tolerance level.")
##   return(list(mu = mu.now, niter = k, err = err))
## }

## Newton1 <- function(mu.trim, R, K, ngrid = 200, tol = 1e-3, max.iter = 10){
##   k <- 0; err <- Inf; mu.now = mu.trim
##   dx <- 2*R/ngrid; mygrid <- seq(-R, R, dx)
##   ## Calculate the Chisq-distribution just once
##   Fchi2 <- pchisq(R^2-mygrid^2, K-1)
##   ## xFchi2, x2Fchi2 are just the pre-calculated intermediate terms
##   xFchi2 <- mygrid*Fchi2; x2Fchi2 <- mygrid*xFchi2
##   while ((abs(err) > tol) & (k < max.iter)) {
##     phi <- dnorm(mygrid-mu.now)
##     ## dx gets cancelled; so no need to multiply each term by dx
##     Bmu <- sum(Fchi2*phi); Tmu <- sum(xFchi2*phi); Hmu <- sum(x2Fchi2*phi)
##     ## dS follows the rule of finding d(f(x)/g(x))
##     Smu <- Tmu/Bmu; dS <- (Hmu*Bmu - Tmu^2)/Bmu^2
##     mu.next <- mu.now - (Smu-mu.trim)/dS
##     k <- k+1; err <- mu.next-mu.now
##     mu.now <- mu.next
##   }
##   if (abs(err) > tol) warning("Maximum number of iteration reached but the error is still greater than the tolerance level.")
##   return(c(mu = mu.now, niter = k, err = err))
## }


## FastNewton <- function(mus.trim, RR, K, min.n.interpolation = 100, ngrid = 200, tol = 1e-3, max.iter = 10){
##   ## Step 1. Linear interpolation
##   Grid1 <- seq(0, max(RR), by = min(RR)/min.n.interpolation)
##   gammahat <- approxfun(Grid1, pchisq(Grid1, K-1))
##   ## Step 2. Individual Newton Raphson
##   mus <- sapply(mus.trim, function(x) Newton1(mu.trim, RR, K, gammahat, ngrid = ngrid, tol = tol, max.iter = max.iter))
##   return(mus)
## }



