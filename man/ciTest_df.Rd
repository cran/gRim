\name{ciTest_df}
\alias{ciTest_df}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Test for conditional independence in a dataframe
}
\description{
  Test for conditional independence in a dataframe.
}
\usage{
ciTest_df(x, set = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    A dataframe.
}
  \item{set}{
%%     ~~Describe \code{set} here~~
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

  If \code{set} consists only of factors then \code{x[,set]} is
  converted to a contingency table and the test is made in this table
  using \code{ciTest_table()}.

  If \code{set} consists only of numeric values and integers then \code{x[,set]} is
  converted to a list with components \code{cov} and \code{n.obs} by
  calling \code{cov.wt(x[,set], method='ML')}. This list is then passed
  on to \code{ciTest_mvn()} which makes the test.
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
  Søren Højsgaard, sorenh at mail dot dk
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{ciTest}}
  
  \code{\link{ciTest.table}}
  \code{\link{ciTest_table}}
  
%  \code{\link{ciTest.data.frame}}
%  \code{\link{ciTest_df}}
  
  \code{\link{ciTest.list}}
  \code{\link{ciTest_mvn}}

  \code{\link{chisq.test}}  

}
\examples{
data(milkcomp1)
ciTest(milkcomp1, set=~tre+fat+pro)
ciTest_df(milkcomp1, set=~tre+fat+pro)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}

