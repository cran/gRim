\name{testdelete}
\alias{testdelete}
\alias{testdelete.iModel}
\alias{print.testdelete}

\alias{testdelete.mModel} %% FIXME: Make this unnecessary

\title{Test deletion of edge from an interaction model}
\description{
  Tests if an edge can be deleted from an interaction model. }
\usage{
testdelete(object, edge, k=2, details=1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A model; an object of class \code{iModel}.}
  \item{edge}{An edge in the model; either as a right-hand sided formula
  or as a vector}
  \item{k}{Penalty parameter used when calculating change in AIC}
  \item{details}{The amount of details to be printed; 0 surpresses
    all information}
  \item{\dots}{Further arguments to be passed on to the underlying
    functions for testing; that is to CItable and CImvn}

}
\details{
  If the  model is decomposable and the edge is contained in one clique only
  then the test is made in the marginal model given by that clique. In
  that case, if
  the model is a log-linear model then degrees of freedom are adjusted
  for sparsity

}
\value{
  A list.
}
%%\references{ ~put references to the literature/web site here ~ }
\author{
  Søren Højsgaard, \email{sorenh@math.aau.dk}
}
% \note{ ~~further notes~~ 

%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
\seealso{\code{\link{testadd}}}
\examples{
## ## ## testdelete
## ## ## 

## ## Discrete model
## ## 
data(reinis)
## A decomposable model
##
mf <- ~smoke:phys:mental+smoke:systol:mental
object <- dmod(mf, data=reinis)

testdelete(object,c("phys","mental"))
testdelete(object,c("smoke","mental"))
#testdelete(object,c("systol","phys"))


## A non-decomposable model
##
mf <- ~smoke:phys+phys:mental+smoke:systol+systol:mental
object <- dmod(mf, data=reinis)

testdelete(object,c("phys","mental"))
#testdelete(object,c("systol","phys"))
#testdelete(object,c("smoke","mental"))


## ## Continuous model
## ## 
data(math)
## A decomposable model
##
mf <- ~me:ve:al+me:al:an
object <- cmod(mf, data=math)

testdelete(object,c("ve","al"))
testdelete(object,c("me","al"))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{htest}% __ONLY ONE__ keyword per line
