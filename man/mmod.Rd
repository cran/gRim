\name{mmod}
\alias{mmod}
\alias{coef.mModel}
\alias{coefficients.mModel}
\alias{print.mModel}
\alias{summary.mModel}
\alias{mmod_dimension}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Mixed interaction model. 
}
\description{
  A mixed interaction model is a model (often with conditional
  independence restrictions) for a combination of discrete and
  continuous variables. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mmod(formula, data, marginal = NULL, fit = TRUE, details = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{marginal}{
%%     ~~Describe \code{marginal} here~~
}
  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{details}{
%%     ~~Describe \code{details} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  An object of class \code{mModel} and the more general class
  \code{iModel}. 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
  Søren Højsgaard, \email{sorenh@math.aau.dk}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{dmod}}
  \code{\link{cmod}}
}
\examples{
### FIXME: To be written
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
