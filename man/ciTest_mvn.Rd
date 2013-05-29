\name{ciTest_mvn}
\alias{ciTest_mvn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Test for conditional independence in the multivariate normal distribution
}
\description{
    Test for conditional independence in the multivariate normal
    distribution.
}
\usage{
ciTest_mvn(x, set = NULL, statistic = "DEV", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    A list with elements \code{cov} and \code{n.obs} (such as returned
    from calling \code{cov.wt()} on a dataframe. See examples below.)
}
  \item{set}{
%%     ~~Describe \code{set} here~~
}
  \item{statistic}{
%%     ~~Describe \code{statistic} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
  \code{set} can be 1) a vector or 2) a right-hand sided formula in
  which variables are separated by '+'. In either case, it is tested if
  the first two variables in the \code{set} are conditionally
  independent given the remaining variables in \code{set}. 
  (Notice an abuse of the '+' operator in the right-hand sided formula:
  The order of the variables does matter.)
  
  If \code{set} is \code{NULL} then it is tested whether the first two
  variables are conditionally independent given the remaining variables.

  \code{x} must be a list with components \code{cov} and \code{n.obs}
  such as returned by  calling \code{cov.wt( , method='ML')} on a dataframe.

}
\value{
  An object of class 'citest' (which is a list).
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
  \code{\link{ciTest}}
  
  \code{\link{ciTest.table}}
  \code{\link{ciTest_table}}
  
  \code{\link{ciTest.data.frame}}
  \code{\link{ciTest_df}}
  
%  \code{\link{ciTest.list}}
%  \code{\link{ciTest_mvn}}

  \code{\link{chisq.test}}  

}
\examples{
data(carcass)
ciTest(cov.wt(carcass, method='ML'), set=~Fat11+Meat11+Fat12)
ciTest_mvn(cov.wt(carcass, method='ML'), set=~Fat11+Meat11+Fat12)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}

