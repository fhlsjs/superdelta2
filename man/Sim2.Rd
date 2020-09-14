\name{Sim2}
\alias{Sim2}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
  %%  ~~function to do ... ~~
  Second type of simulation scheme designed in \code{superdelta2} package to generate random number of read count data from a Negative Binomial Poisson (NBP) model with 3 groups, given mean \bold{log}-counts parameter.
}

\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  This function is developed based on R's rnbinom() function.
  \code{mus = list(mu1, mu2, mu3)} are the vectors of mean log-counts in the three groups.
  \code{ns = c(n1, n2, n3)} are the sample sizes in the three groups.
  \code{kappa} and \code{a} are two shape parameters in the NBP model.
  \code{l} and \code{u} are the lower and upper bounds (uniform distribution) of simulated sample specific noise.
}

\usage{
  Sim2(mus, ns, kappa, a, l, u)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mus}{
  %%     ~~Describe \code{x} here~~
  \code{mus = list(mu1, mu2, mu3)} are the log-mean parameter vectors of the three groups respectively.
  }
  \item{ns}{\code{ns = c(n1, n2, n3)} are the sample sizes of the three groups respectively.}
  \item{kappa}{The first shape parameter.}
  \item{a}{The second shape parameter.}
  \item{l}{Lower bound of uniform distribution of simulated sample specific noise alpha_j.}
  \item{u}{Upper bound of uniform distribution of simulated sample specific noise alpha_j.}
  %%     ~~Describe \code{x} here~~
}

\details{
  An NBP distribution is an integer-valued distribution with three parameters, the location parameter \code{mu}, and two shape parameters \code{kappa} and \code{a}. The mean and variance of X ~ NBP(mu, kappa, a) are:
  \emph{E(X) = mu}
  \emph{Var(X) = mu + mu*kappa^a}
  Note that the NBP distribution describes a nonlinear relationship between the mean and variance of genes. Other technical details of the NBP distribution, such as the probability density function and its relationship with the negative binomial (NB) distribution are not covered here.
  %%  ~~ If necessary, more details than the description above ~~
}

\value{
%%  ~Describe the value returned
This function returns a list with the following objects:
  %%  If it is a LIST, use
  \item{counts }{A count matrix simulated by the scheme discussed above.}
  \item{alphas }{Simulated sample specific noise alpha_j added to each sample/library.}
  %% ...
}

\references{
  %% ~put references to the literature/web site here ~
}

\author{
  %%  ~~who you are~~
  Yuhang Liu, Xing Qiu, Jinfeng Zhang, and Zihan Cui
}

\note{
  %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  SIM1, SIM2, SIM3
  %% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
  ## The code block below generates sample data \code{SIM2}.
  ngenes <- 5000; n1 <- n2 <- n3 <- 50; ns <- c(n1,n2,n3)
  mu1 <- mu2 <- mu3 <- rep(log2(100), ngenes)
  mu2[1:600] <- log2(150); mu3[401:1000] <- log2(75)
  mus <- list(mu1, mu2, mu3)
  set.seed(2020)
  mod2 <- Sim2(mus, ns, kappa = 0.06, a = 2.2, l = 12, u = 30)
  SIM2 <- mod2$counts; alphas <- mod2$alphas
  
  ## The code block below generates sample data \code{SIM3}.
  ngenes <- 5000; n1 <- n2 <- n3 <- 50; ns <- c(n1,n2,n3)
  mu1 <- mu2 <- mu3 <- rep(log2(100), ngenes)
  mu2[1:600] <- log2(150); mu3[401:1000] <- log2(50)
  mus <- list(mu1, mu2, mu3)
  set.seed(2020)
  mod3 <- Sim2(mus, ns, kappa = 0.6, a = 2, l = 10, u = 20)
  SIM3 <- mod3$counts; alphas <- mod3$alphas
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}% use one of  RShowDoc("KEYWORDS")
\keyword{distribution}% __ONLY ONE__ keyword per line
