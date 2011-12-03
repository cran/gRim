\name{CGstats}
\alias{CGstats}
\alias{CGstats.data.frame}
\alias{CGstats_internal}
\alias{print.CGstats}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Mean, covariance and counts for grouped data
}
\description{
  \code{CGstats} provides what corresponds to calling \code{cow.wt} on
  different strata of data where the strata are defined by the
  combinations of factors in data.   
}
\usage{
CGstats(object, varnames = NULL, homogeneous = TRUE, simplify = TRUE)
\method{CGstats}{data.frame}(object, varnames = NULL, homogeneous = TRUE, simplify = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{varnames}{
%%     ~~Describe \code{varnames} here~~
}
  \item{homogeneous}{
%%     ~~Describe \code{homogeneous} here~~
}
  \item{simplify}{
%%     ~~Describe \code{simplify} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list whose form depends on the type of input data and the varnames.
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
  \code{\link{cov.wt}}
}
\examples{

data(milkcomp)
% milkcomp <- subset(milkcomp, (treat %in% c("a","b")) & (lactime %in% c("t1", "t2")))
% milkcomp <- milkcomp[,-1]
% milkcomp$treat 		<- factor(milkcomp$treat)
% milkcomp$lactime 	<- factor(milkcomp$lactime)

CGstats(milkcomp)
CGstats(milkcomp, c(1,2))
CGstats(milkcomp, c("lactime","treat"))
CGstats(milkcomp, c(3,4))
CGstats(milkcomp, c("fat","protein"))

CGstats(milkcomp, c(2,3,4), simplify=FALSE)
CGstats(milkcomp, c(2,3,4), homogeneous=FALSE)
CGstats(milkcomp, c(2,3,4), simplify=FALSE, homogeneous=FALSE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
