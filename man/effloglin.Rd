\name{effloglin}
\alias{effloglin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting Log-Linear Models by Message Passing}
\description{
  Fit log-linear models to multidimensional contingency tables by
  Iterative Proportional Fitting. 
}
\usage{
effloglin(table, margin, fit = FALSE, eps = 0.01, iter=20, print = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{A contingency table}
  \item{margin}{A generating class for a hierarchical log--linear model}
  \item{fit}{If TRUE, the fitted values are returned.}
  \item{eps}{Convergence limit; see 'details' below.}
  \item{iter}{Maximum number of iterations allowed}
  \item{print}{If TRUE, iteration details are printed.}
}
\details{
  The function differs from
  \code{loglin} in that 1) data can be given in the form of a
  list of sufficient marginals and 2) the model is fitted only on the
  cliques of the triangulated interaction graph of the model. This means
  that the full table is not fitted, which means that \code{effloglin}
  is efficient (in terms of storage requirements). However
  \code{effloglin} is implemented entirely in R and is therefore slower
  than \code{loglin}.
}
\value{
  A list with compnents
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
%\references{}
\author{Søren Højsgaard, sorenh@agrsci.dk}
% \note{ ~~further notes~~ 

%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
\seealso{\code{\link{loglin}}}
\examples{
data(reinis)
glist <-list(c("smoke", "mental"), c("mental", "phys"), c("phys", "systol"
), c("systol", "smoke"))

stab <- lapply(glist, function(gg) tableMargin(reinis, gg))
fv3 <- effloglin(stab, glist, print=FALSE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

