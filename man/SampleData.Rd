\name{SampleData}
\alias{SampleData}
\alias{SIM1}
\alias{SIM2}
\alias{SIM3}
\docType{data}

\title{
  %%   ~~ data name/kind ... ~~
  Sample datasets to demonstrate performance of \code{superdelta2} method.
}

\description{
  %%  ~~ A concise (1-5 lines) description of the dataset. ~~
  This R workspace contains three simulated datasets based on negative binomial Poisson (NBP) distribution.
  Specifically, \code{SIM1} is obtained from function \code{Sim1}.
  \code{SIM2} and \code{SIM3} are obtained from function \code{Sim2}, with different parameters.
}

\usage{data(SampleData)}
\format{
  Three data matrices each with 150 samples (50 per group) and 5000 genes.
  \describe{
    \item{\code{SIM1}}{An integer-valued matrix simulated by \code{Sim1}. For more details see the help page of \code{Sim1}.}
    \item{\code{SIM2}}{An integer-valued matrix simulated by \code{Sim2}. For more details see the help page of \code{Sim2}.}
    \item{\code{SIM3}}{An integer-valued matrix simulated by \code{Sim2}. For more details see the help page of \code{Sim2}.}
  }
}

\details{
  %%  ~~ If necessary, more details than the __description__ above ~~
  \preformatted{
  }
}

\references{
  %%  ~~ possibly secondary sources and usages ~~
}

\examples{
  \code{
  ## Load the sample data
  data(SampleData)
  ## Number of genes and samples
  ngenes <- 5000; n1 <- n2 <- n3 <- 50
  Groups <- c(rep("A",n1), rep("B",n2), rep("C",n3))
  mod1 <- superdelta2(mydata = SIM1, offset = 1, Grps = Groups)
  mod2 <- superdelta2(mydata = SIM2, offset = 1, Grps = Groups)
  mod3 <- superdelta2(mydata = SIM3, offset = 1, Grps = Groups)
  }
}

\keyword{datasets}
