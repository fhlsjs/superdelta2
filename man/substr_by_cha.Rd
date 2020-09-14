\name{substr_by_cha}
\alias{substr_by_cha}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Divide character strings given a certain character.
}

\description{
  This is a handy function to divide character strings given a certain character. This is especially useful
  when sometimes gene names show up in the form of "GENESYMBOL|ENTREZID"
}

\usage{
  substr_by_cha(cha_str, delim = "|", type = c("before", "after"))
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{cha_str}{A character string or a vector of character string.}
  \item{delim}{Delimiter character. Default is "|".}
  \item{type}{Return the partial character string before or after the delimiter.}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
  Return the desired substring (an individual or a vector).
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
  %% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
  substr_by_cha("ABCB1|5243", delim = "|", type = "before")
  ## [1] "ABCB1"
  substr_by_cha(c("ABCB1-5243", "ABL1-25", "ABL2-27"), delim = "-", type = "after")
  ## [1] "5243"  "25"  "27"
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}% use one of RShowDoc("KEYWORDS")
\keyword{misc}