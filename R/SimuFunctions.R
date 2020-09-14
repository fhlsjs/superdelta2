## based on R's rnbinom() function.
## nus=list(nu1, nu2, nu3): vectors of mean counts in the three groups;
## ns=c(n1, n2, n3): sample sizes in the three groups
## kappa, a are two shape parameters in the NBP model.
Sim1 <- function(nus, ns, kappa, a, l, u){
  mm <- sapply(nus, length); G <- length(mm)
  if (var(mm) !=0) {
    stop("nus must be a list of three vectors of mean counts with the same lengths.") } else {
      ngenes <- mm[1]
    }
  X <- NULL
  for (g in 1:G){
    gg <- nus[[g]]^(2-a)/kappa; pg <- gg/(nus[[g]] + gg)
    Xg <- t(sapply(1:ngenes, function(i) rnbinom(ns[1], size = gg[i], prob = pg[i])))
    X <- cbind(X, Xg)
  }
  ## add sample specific noise
  alphas <- runif(n1 + n2 + n3, min = l, max = u)
  Y <- round(sweep(X, 2, alphas, "*"))
  ## return the counts
  return(Y)
}

############################################################
##### SIM2
############################################################
Sim2 <- function(mus, ns, kappa, a, l, u){
  mm <- sapply(mus, length); G <- length(mm)
  if (var(mm) !=0) {
    stop("mus must be a list of three vectors of mean log-counts with the same lengths.") } else {
      ngenes <- mm[1]
    }
  Y <- NULL; alphas <- c()
  for (g in 1:G){
    ## generate alpha.j
    alphas.g <- runif(ns[g], min = l, max = u); alphas <- c(alphas, alphas.g)
    nus.g <- 2^(mus[[g]]) %*% t(alphas.g)
    ## turn nus.g, kappa, and a into gamma_g and p_g in NB model
    gg <- nus.g^(2-a)/kappa; pg <- gg/(nus.g + gg)
    Yg <- matrix(mapply(function(g, p) {rnbinom(n = 1, size = g, prob = p)}, gg, pg), nrow = ngenes)
    Y <- cbind(Y, Yg)
  }
  return(list(counts = Y, alphas = alphas))
}

# x <- array(1:9, c(3,3))
# x
# sapply(x, function(x) x * 10)
# sapply(c(1:3), function(x) x^2)
