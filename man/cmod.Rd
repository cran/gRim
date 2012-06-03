\name{cmod}
\alias{cmod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graphical Gaussian model}
\description{
  Specification of graphical Gaussian model. The 'c' in the name
  \code{cmod} refers to that it is a (graphical) model for 'c'ontinuous variables
}
\usage{
cmod(formula, data, marginal = NULL, fit = TRUE, details=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{Model specification in one of the following forms:
    1) a right-hand sided formula,
    2) as a list of generators,
    3) an undirected graph (represented either as a graphNEL object or as an adjacency
    matrix).
    Notice that there are certain model
    specification shortcuts, see Section 'details' below}
  \item{data}{Data in one of the following forms: 1) A dataframe or 2) a
  list with elements \code{cov} and \code{n.obs} (such as returned by
  the \code{cov.wt()} function.) }
  \item{marginal}{Should only a subset of the variables be used in
    connection with the model specification shortcuts}
  \item{fit}{Should the model be fitted.}
  \item{details}{Control the amount of output; for debugging purposes.}
}
\details{
  The independence model can be specified as \code{~.^1} and the
  saturated model as \code{~.^.}. The \code{marginal}
  argument can be used for specifying the independence or saturated
  models for only a subset of the variables.
}
\value{
  An object of class \code{cModel} (a list)
%   ~Describe the value returned
%   If it is a LIST, use
%   \item{comp1 }{Description of 'comp1'}
%   \item{comp2 }{Description of 'comp2'}
%   ...
}
%\references{ ~put references to the literature/web site here ~ }
\author{
  Søren Højsgaard, sorenh@math.aau.dk
}
% \note{ ~~further notes~~ 

%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
\seealso{
  \code{\link{dmod}}
  \code{\link{mmod}}
  \code{\link{ggmfit}}}
\examples{
## Graphical Gaussian model
data(carcass)
cm1<-cmod(~.^., carcass)

## Stepwise selection based on BIC
cm2<-backward(cm1,k=log(nrow(carcass)))

## Stepwise selection with fixed edges
cm3<-backward(cm1,k=log(nrow(carcass)),
 fixinMAT=matrix(c("LeanMeat","Meat11","Meat12","Meat13","LeanMeat","Fat11","Fat12","Fat13"),
 ncol=2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

