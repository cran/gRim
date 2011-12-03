\name{dmod}
\alias{dmod}
\alias{print.dModel}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Log--linear model}
\description{
  Specification of log--linear (graphical) model. The 'd' in the name
  \code{dmod} refers to that it is a (graphical) model for 'd'iscrete variables
}
\usage{
dmod(formula, data, marginal, interactions=NULL, fit = TRUE, details=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{Model specification in one of the following forms:
    1) a right-hand sided formula,
    2) as a list of generators,
    3) an undirected graph (represented either as a graphNEL object or as an adjacency
    matrix).
    Notice that there are certain model
    specification shortcuts, see Section 'details' below.}
  \item{data}{Either a table or a dataframe. In the latter case, the
    dataframe will be coerced to a table. See 'details' below.}
  \item{interactions}{A number given the highest order interactions in
    the model, see Section 'details' below.}
  \item{marginal}{Should only a subset of the variables be used in
    connection with the model specification shortcuts}
  \item{fit}{Should the model be fitted.}
  \item{details}{Control the amount of output; for debugging purposes.}
}

\details{
  The independence model can be specified as \code{~.^1} and the
  saturated model as \code{~.^.}.
  Setting e.g. \code{interactions=3} implies that there will be at most
  three factor interactions in the model. 
  
  Data can be specified as a table
  of counts or as a dataframe. If data is a dataframe then it will be
  converted to a table (using \code{xtabs()}). This means that if the
  dataframe contains numeric values then the you can get a very sparse
  and high dimensional table. When a dataframe contains numeric values
  it may be worthwhile to discretize data using the \code{cut()}
  function.

  The \code{marginal} argument can be used for specifying the
  independence or saturated models for only a subset of the
  variables. When \code{marginal} is given the corresponding marginal
  table of data is formed and used in the analysis (notice that this is
  different from the behaviour of \code{loglin()} which uses the full
  table. 
}
\value{
  An object of class \code{dModel}
%   ~Describe the value returned
%   If it is a LIST, use
%   \item{comp1 }{Description of 'comp1'}
%   \item{comp2 }{Description of 'comp2'}
%   ...
}
%\references{ ~put references to the literature/web site here ~ }
\author{
  Søren Højsgaard, sorenh at mail dot dk
}
% \note{ ~~further notes~~ 

%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
\seealso{
  \code{\link{cmod}}
  \code{\link{mmod}}
}
\examples{

## Graphical log-linear model
data(reinis)
dm1<-dmod(~.^., reinis)
dm2<-backward(dm1, k=2)
dm3<-backward(dm1, k=2, fixin=list(c("family","phys","systol")))
## At most 3-factor interactions
dm1<-dmod(~.^., data=reinis,interactions=3)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

