\name{ciTest_table}
\alias{ciTest_table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Test for conditional independence in a contingency table
}
\description{
    Test for conditional independence in a contingency table
}
\usage{
ciTest_table(x, set = NULL, statistic = "dev", method = "chisq",
adjust.df = TRUE, slice.info = TRUE, L = 20, B = 200, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    A contingency table.
}
  \item{set}{
%%     ~~Describe \code{set} here~~
}
  \item{statistic}{
%%     ~~Describe \code{statistic} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{adjust.df}{
%%     ~~Describe \code{adjust.df} here~~
}
  \item{slice.info}{
%%     ~~Describe \code{slice.info} here~~
}
  \item{L}{
%%     ~~Describe \code{L} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{...}{
    Additional arguments.
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
  Søren Højsgaard, sorenh@math.aau.dk 
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{ciTest}}
  
%  \code{\link{ciTest.table}}
%  \code{\link{ciTest_table}}
  
  \code{\link{ciTest.data.frame}}
  \code{\link{ciTest_df}}
  
  \code{\link{ciTest.list}}
  \code{\link{ciTest_mvn}}

  \code{\link{chisq.test}}  

}
\examples{
data(reinis)
ciTest(reinis, set=~smo+phy+sys)
ciTest_table(reinis, set=~smo+phy+sys)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
