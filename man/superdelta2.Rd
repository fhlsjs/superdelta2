\name{superdelta2}
\alias{superdelta2}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
  %%  ~~function to do ... ~~
  Main function of superdelta2 package to implement differential expression analysis in a multiple (>=3) group setting, as well with Tukey style post-hoc pairwise comparison. Input arguments are flexible and adaptable.
}

\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  This function takes in the raw read count data matrix from RNA-Seq experiment, with rows corresponding to genes and 
  columns corresponding to samples/libraries, and conducts differential expression analysis based on the robust superdelta2 method. Whether or not to use the total number of mapped read counts per-library as sample weights (the \code{W} parameter) is up to user's decision.
}

\usage{
  superdelta2(mydata, offset = 1, Grps, W = NULL, trim = .2,
  adjp.thresh = 0.05, prop = 1.0)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{mydata}{A data matrix with rows being genes and columns being samples. Note that this argument 
  requires RAW read count of RNA-seq gene expression data, but not pre-processed data.}
  \item{offset}{A value to be added when taking logarithm to avoid log of zero. 
  Default value is 1. Commonly used values include 0.5, 1, 2, and 3.}
  \item{Grps}{A character vector of length equal to the total number of samples in the data matrix, indicating group labels. This object is used to create design matrix for further computation.}
  \item{W}{Sample weights to be used. Default is NULL, which means unweighted (equally weighted). Users may want to use per-sample total read counts as the weighting metric. This is shown to have marginal gain in some cases.}
  \item{trim}{Trimming proportion when removing the most extreme values of between-group sum of squares. Default is 0.2.}
  \item{adjp.thresh}{Cutoff of adjusted p-value to define significance. Default is 0.05.}
  \item{prop}{The proportion of reference genes used in the main (second
  round) spherical trimming for bias correction. The default value is
  1.0, which means using all genes for bias correction.  Using a small
  proportion can improve the computational speed, at the cost of
  slightly less accurate bias correction.}
}

\details{
  If not \code{NULL}, \code{W} must be a vector of positive per-sample
  weights of length equal to the number of columns (samples) in the data matrix,
  to specify per-sample weights. The program will automatically
  normalize the weights so that they sum up to 1. We recommend using the total
  number of mapped read counts as per-sample weights.
}

\value{
The \code{superdelta2} function will return a list with the following objects:
  \item{SigGenes}{A named vector of the significant gene IDs. Typically
the name of this vector is inherited from the rownames of the input data matrix.}
  \item{Fstats}{A vector of length equal to the number of genes in the data matrix indicating the super-delta F-statistic for each gene.}
  \item{pvalues.F}{Raw p-values obtained from the super-delta F-statistics.}
  \item{padj.F}{Benjamini adjusted p-values obtained from the super-delta F-statistics.}
  \item{Log2FC}{Log2 fold-change of post-hoc pairwise t-tests. Technically this is equal to the numerator of Tukey style t-statistics.}
  \item{tstats}{A matrix of row dimension equal to the number of genes in the input data matrix, providing the super-delta post-hoc Tukey t-statistics.}
  \item{pvalues.t}{Raw p-values obtained from the super-delta Tukey t-statistics.}
  \item{padj.t}{Benjamini adjusted p-values obtained from the super-delta Tukey t-statistics.}
  \item{WBGRMS}{Estimated weighted between-group mean square, namely the top part of F-statistics.}
  \item{WWGRMS}{Estimated weighted within-group mean square, namely the bottom part of F-statistics.}
  \item{sigma2hat}{Estimated sigma square (variance of the error terms).}
  \item{Rtrim}{Saved R object containing the information of the main trim.}
}

\references{
  %% ~put references to the literature/web site here ~
  \cite{Liu, Y., Zhang, J., & Qiu, X. (2017). Super-delta: a new differential gene expression analysis procedure with robust data normalization. BMC bioinformatics, 18(1), 582.}
}

\author{
  %%  ~~who you are~~
  Yuhang Liu, Xing Qiu, Jinfeng Zhang, and Zihan Cui
}

\note{
  %%  ~~further notes~~
  The \code{superdelta2} function implements a general-purpose differential gene expression analysis procedure that can be used in a multiple group (>=3) setting. It inherits the basic idea of \emph{superdelta} (see \bold{References} for details) with a robust internal normalization and an asymptotically unbiased estimator of the ``oracle'' between-group difference. It applies a robust trimming procedure to the estimated between-group sum of squares and takes advantage of all the within-group sum of squares, to end up with the \emph{superdelta} F-statistics. In addition, \code{superdelta2} is accompanied with a Tukey style pairwise comparison to also obtain post-hoc t-statistics.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  SIM1, SIM2, SIM3
}

\examples{
  ## Load the sample data
  data(SampleData)
  ## Number of genes and samples
  ngenes <- 5000; n1 <- n2 <- n3 <- 50; ns <- c(n1,n2,n3)
  Groups <- c(rep("A",n1), rep("B",n2), rep("C",n3))
  mod1 <- superdelta2(mydata = SIM1, offset = 1, Grps = Groups)
  mod2 <- superdelta2(mydata = SIM2, offset = 1, Grps = Groups)
  mod3 <- superdelta2(mydata = SIM3, offset = 1, Grps = Groups)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{robust}
\keyword{htest}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line

