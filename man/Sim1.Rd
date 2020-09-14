\name{Sim1}
\alias{Sim1}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
  First type of simulation scheme designed in \code{superdelta2} package to generate random number of read count data from a Negative Binomial Poisson (NBP) model with 3 groups, given mean counts parameter.
  %%  ~~function to do ... ~~
}

\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  This function is developed based on R's rnbinom() function.
  \code{nus = list(nu1, nu2, nu3)} are the vectors of mean counts in the three groups.
  \code{ns = c(n1, n2, n3)} are the sample sizes in the three groups.
  \code{kappa} and \code{a} are two shape parameters in the NBP model.
  \code{l} and \code{u} are the lower and upper bounds (uniform distribution) of simulated sample specific noise.
}

\usage{
  Sim1(nus, ns, kappa, a, l, u)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nus}{
  %%     ~~Describe \code{x} here~~
  \code{nus = list(nu1, nu2, nu3)} are the mean parameter vectors of the three groups respectively.
  }
  \item{ns}{\code{ns = c(n1, n2, n3)} are the sample sizes of the three groups respectively.}
  \item{kappa}{The first shape parameter.}
  \item{a}{The second shape parameter.}
  \item{l}{Lower bound of uniform distribution of simulated sample specific noise alpha_j.}
  \item{u}{Upper bound of uniform distribution of simulated sample specific noise alpha_j.}
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
  %%  If it is a LIST, use
  This function returns a count matrix \code{Y} simulated by the scheme discussed above.
  %% \item{comp1 }{Description of 'comp1'}
  %% \item{comp2 }{Description of 'comp2'}
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
  ngenes <- 5000; n1 <- n2 <- n3 <- 50; ns <- c(n1,n2,n3)
  nu1 <- nu2 <- nu3 <- rep(100, ngenes)
  nu2[1:600] <- 150; nu3[401:1000] <- 75
  nus <- list(nu1, nu2, nu3)
  set.seed(2020)
  SIM1 <- Sim1(nus, ns, kappa = 0.06, a = 2.2, l = 12, u = 30)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen} % use one of  RShowDoc("KEYWORDS")
\keyword{distribution} % __ONLY ONE__ keyword per line
