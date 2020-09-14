## Ancillary function: Global normalization
## the main function of *weighted* F-tests and post hoc tests, with pre-trim.
## Only need to change lines 16-32, and lines 72-76
superdelta2 <- function(mydata, offset = 1.0, Grps, W = NULL, trim = .2, adjp.thresh = 0.05, prop = 1.0){
  mydata <- log2(mydata + offset)
  m <- nrow(mydata); N <- ncol(mydata); G <- length(unique(Grps))
  ## if null then equal weights, otherwise scale weights to summation = 1
  if (is.null(W)){
    W <- rep(1/N, N)
  } else{
    W <- W/sum(W)
  }
  if (length(W)!=N){
    stop("Length of input weights must match data!")
  }
  ## Generate a Grp factor based on Ns
  #Grps <- rep(1:G, times = Ns)
  ## Generate a Grp matrix (design matrix) to assist the calculation of mean/RSS
  Gmat <- model.matrix(~ 0 + as.factor(Grps))
  ## Per-group summation of weights
  W.g <- W%*%Gmat
  ## Below we follow equation (18) to estimate weighted within-group sum of squares
  GLBNor.data <- normalize.global(mydata)
  #GLBNor.center <- matrix(data = NA, nrow = m, ncol = G)
  #for (g in 1:G){
  #GLBNor.center[,g] <- apply(GLBNor.data[,(Nc[g]+1):Nc[g+1]], 1, w_mean, w = W[(Nc[g]+1):Nc[g+1]])
  #}
  ## GLBNor.center is m*G (within-group weighted mean of globally normalized expressions)
  GLBNor.center <- GLBNor.data %*% diag(W) %*% Gmat
  GLBNor.center <- sweep(GLBNor.center, 2, W.g, FUN = "/")
  GLBNor.res <- GLBNor.data - GLBNor.center %*% t(Gmat)
  ## Compute weighted within-group sum of squares (using \tilt{Res}_{i,g}), equation (21)
  GLBNor.res.sq <- GLBNor.res^2
  #GLBNor.WW <- matrix(data = NA, nrow = m, ncol = G)
  #for (g in 1:G){
  #GLBNor.WW[,g] <- sum(W[(Nc[g]+1):Nc[g+1]])*apply(GLBNor.res.sq[,(Nc[g]+1):Nc[g+1]], 1, w_mean, w = W[(Nc[g]+1):Nc[g+1]])
  #}
  ## GLBNor.WW is m*G (within-group residual sum of squares)
  GLBNor.WW <- GLBNor.res.sq %*% diag(W) %*% Gmat
  WWGRSS_hat <- apply(GLBNor.WW, 1, sum)
  Bottom <- WWGRSS_hat/(N-G)
  ## Next line is copied from the old version
  sigma2hat <- (m-1)/m*mean(Bottom); sigma <- sqrt(sigma2hat)
  ## centering data (subtract row means). Z is m*N
  #Z <- sweep(mydata, 1, rowMeans(mydata))
  ## 1. Per-group Z means (m*G matrix). Zbar is m*G
  #Zbar <- Z %*% Gmat %*% diag(1/Ns)
  ## Per-group centered Z. ResProj is N*N
  #ResProj <- (diag(N) - Gmat %*% diag(1/Ns) %*% t(Gmat))
  #Ztilde <- Z %*% ResProj #Ztilde is m*N
  #ZtildeBar <- colMeans(Ztilde)
  #WGRSS.N <- sweep(Ztilde, 2, ZtildeBar)^2 #WGRSS.N is m*N matrix
  #WGRSS <- m^2 / (m-1)^2 * rowSums(WGRSS.N)
  #Bottom <- WGRSS / (N-G)
  ## estimating sigma2
  #sigma2hat <- (m-1)/m * mean(Bottom); sigma <- sqrt(sigma2hat)
  ## Now let us compute individual R_ik,g.  Needed for trimmed mean
  ## calculation and pairwise t-tests. First TIME CONSUMING operation.
  ## Ybarbar is m*1 (overall weighted mean), equation (23)
  Ybarbar <- mydata %*% matrix(data = W, nrow = N, ncol = 1)
  ## Z is m*N (matrix of residuals), equation (23)
  Z <- mydata - Ybarbar %*% matrix(data = rep(1,N), nrow = 1)
  ## Ybar.g is m*G (per-group weighted mean), equation (23)
  Ybar.g <- mydata %*% diag(W) %*% Gmat
  Ybar.g <- sweep(Ybar.g, 2, W.g, FUN = "/")
  ## Zbar.g is m*G (matrix of per-group weighted mean minus overall weighted mean), equation (23)
  Zbar.g <- Ybar.g - Ybarbar %*% matrix(data = rep(1,G), nrow = 1)
  ## compute individual R_ik,g, equation (24)
  ## print(G)
  Rlist <- lapply(1:G, function(g) sqrt(W.g[g])*outer(Zbar.g[,g], Zbar.g[,g], "-"))
  ## BGRSSmat is the between group RSS
  Rlist.norm2 <- Reduce("+", Map(function(x) x^2, Rlist))
  ## Find the indices of BGRSS that are less than 80% quantile. Second TIME CONSUMING part.
  trimmed.m <- round((m-1)*(1-trim), 0)
  if (trim <= 0.03){
    warning("Your trimming factor may be too small. Please consider increasing the value of it.")
    threshs <- rep(Inf, m); IDX <- matrix(TRUE, m, m)
  } else if (trim >= 0.8) {
    stop("Your trimming factor is too large. Please consider decreasing the value of it.")
  } else {
    ## pre-trim
    pretrim <- SphereTrim(Rlist, matrix(0, m, G), trim, sigma, prop = .1)
    Rtrim0 <- pretrim$Rtrim; PreTrimRadii <- sqrt(pretrim$threshs)
    ## The main trimming, with bias correction
    maintrim <- SphereTrim(Rlist, Rtrim0, trim = trim, sigma, prop = prop, bias.correct = TRUE)
    Rtrim <- maintrim$Rtrim; TrimRadii <- sqrt(maintrim$threshs)
    ## for debugging use: Rtrim1 is the estimated centers after
    ## pre-trim but not bias-corrected.
    Rtrim1 <- maintrim$Rtrim1
  }
  ## The Top part of the F-statistics
  ## Top <- Reduce("+", Map(function(x) x^2, Rlist.trimmean))/(G-1)
  Top <- rowSums(Rtrim^2)/(G-1)
  ## The F-statistics
  Fstats <- Top/Bottom
  ################### All pairwise Tukey-style t-statistics ####################
  pairs <- combn(G, 2)
  pairs.psd <- matrix(sort(unique(Grps))[pairs], nrow = 2)
  pairs.names <- apply(pairs.psd, 2, function(x) paste(paste0("Group", x), collapse = "_vs_"))
  colnames(pairs) <- pairs.names
  ## MeanDiffs are the Top of paired t-tests
  MeanDiffs <- sapply(1:ncol(pairs), function(j) {
    a <- pairs[1, j]; b <- pairs[2, j]
    Rtrim[,a]/sqrt(W.g[a]) - Rtrim[,b]/sqrt(W.g[b])
  })
  colnames(MeanDiffs) <- pairs.names
  ## Unlike standard two-sample t-stats, Tukey-style t-statistics can use STD pooled from all samples. 
  ## Note that the constant that needs to be divided is the same as in the two-sample t-statistics.
  t.const <- sapply(1:ncol(pairs), function(j) {
    a <- pairs[1, j]; b <- pairs[2, j]
    sqrt(1/W.g[a] + 1/W.g[b])
  })
  Tukey.Bottom <- sqrt(Bottom) %*% t(t.const)
  tstats <- MeanDiffs/Tukey.Bottom
  ## if the original data have rownames (gene names), use them. Otherwise use row indices instead
  if (is.null(rownames(mydata))) { rownames(mydata) <- 1:m }
  names(Fstats) <- rownames(mydata)
  rownames(tstats) <- rownames(mydata)
  ## p-values and multiple testing adjustment (based on BH)
  pvalues.F <- pf(Fstats, G-1, N-G, lower.tail = FALSE)
  padj.F <- p.adjust(pvalues.F, method = "BH")
  ## two-sided t-test is F-test with df1 = 1
  pvalues.t <- pf(tstats^2, 1, N-G, lower.tail = FALSE)
  ## run BH proc only for those genes selected by F-test
  SigGenes <- which(padj.F < adjp.thresh)
  #if (length(SigGenes) != 0 ) {
  ## This line is very important in rearranging the adjusted p-values
    #padj.t <- matrix(p.adjust(pvalues.t[SigGenes,], method = "BH"), nrow = length(SigGenes))
    #rownames(padj.t) <- names(SigGenes)
    #colnames(padj.t) <- colnames(tstats)
  #}
  padj.t <- matrix(p.adjust(pvalues.t, method = "BH"), nrow = nrow(tstats))
  rownames(padj.t) <- rownames(tstats)
  colnames(padj.t) <- colnames(tstats)
  ## Output of Top and Bottom may not be necessary later on.
  return(list(SigGenes = SigGenes, Fstats = Fstats, pvalues.F = pvalues.F, padj.F = padj.F,
              Log2FC = MeanDiffs, tstats = tstats, pvalues.t = pvalues.t, padj.t = padj.t,
              WBGRMS = Top, WWGRMS = Bottom, sigma2hat = sigma2hat, Rtrim = Rtrim))
}