\name{ggmfit}
\alias{ggmfit}
\alias{ggmfitr}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Iterative proportional fitting of graphical Gaussian model}
\description{
  Fit graphical Gaussian model by iterative proportional fitting.
}
\usage{
ggmfit( S, n.obs, glist, start=NULL, eps=1e-12, iter=1000, details=0, ...)
ggmfitr(S, n.obs, glist, start=NULL, eps=1e-12, iter=1000, details=0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Empirical covariance matrix}
  \item{n.obs}{Number of observations}
  \item{glist}{Generating class for model (a list)}
  \item{start}{Initial value for concentration matrix}
  \item{eps}{Convergence criterion}
  \item{iter}{Maximum number of iterations}
  \item{details}{Controlling the amount of output.}
  \item{...}{Optional arguments; currently not used}
}
\details{
  \code{ggmfit} is based on a C implementation.  
  \code{ggmfitr} is implemented purely in R (and is provided mainly as a
  benchmark for the C-version).
}
\value{
  A list with
  \item{lrt}{Likelihood ratio statistic (-2logL)}
  \item{df}{Degrees of freedom}
  \item{logL}{log likelihood}
  \item{K}{Estimated concentration matrix (inverse covariance matrix)}

}
%\references{FIXME}
\author{
  Søren Højsgaard, sorenh at mail dot dk
}
% \note{ ~~further notes~~ 

%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
\seealso{\code{\link{cmod}}, \code{\link{loglin}}}
\examples{
## Fitting "butterfly model" to mathmark data
## Notice that the output from the two fitting functions is not
## entirely identical.
data(math)
ddd <- cov.wt(math, method="ML")
glist <- list(c("al","st","an"), c("me","ve","al"))
ggmfit (ddd$cov, ddd$n.obs, glist)
ggmfitr(ddd$cov, ddd$n.obs, glist)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{models}% __ONLY ONE__ keyword per line
