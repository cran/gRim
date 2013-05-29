\name{testadd}
\alias{testadd}
\alias{testadd.iModel}
\alias{print.testadd}

\alias{testadd.mModel} %% FIXME: Make this unnecessary

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test addition of edge to graphical model}
\description{
  Performs a test of addition of an edge to a graphical model (an
  \code{iModel} object). 
}
\usage{
testadd(object, edge, k=2, details=1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A model; an object of class \code{iModel}.}
  \item{edge}{An edge; either as a vector or as a right hand sided formula.}
  \item{k}{Penalty parameter used when calculating change in AIC}
  \item{details}{The amount of details to be printed; 0 surpresses all information}
  \item{\dots}{Further arguments to be passed on to the underlying
    functions for testing; that is to CItable and CImvn}
}
\details{
  Let M0 be the model and e=\{u,v\} be an edge and let M1 be the model
  obtained by adding e to M0. If M1 is decomposable AND e is contained
  in one clique C only of M1 then the test is carried out in the
  C-marginal model. In this case, and if the model is a log-linear model
  then the degrees of freedom is adjusted for sparsity.
}
\value{
  A list
}
%\references{ ~put references to the literature/web site here ~ }
\author{
  Søren Højsgaard, \email{sorenh@math.aau.dk}
}
% \note{ ~~further notes~~ 

%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
\seealso{\code{\link{testdelete}}}
\examples{
## ## ## testadd
## ## ## 

## ## Discrete model
## ## 
data(reinis)
## A decomposable model
##
mf <- ~smoke:phys:mental+smoke:systol:mental
object <- dmod(mf, data=reinis)
testadd(object,c("systol","phys"))


## A non-decomposable model
##
mf <- ~smoke:phys+phys:mental+smoke:systol+systol:mental
object <- dmod(mf, data=reinis)
testadd(object,c("phys","systol"))


## ## Continuous model
## ## 
data(math)
## A decomposable model
##
mf <- ~me:ve:al+al:an
object <- cmod(mf, data=math)
testadd(object,c("me","an"))

## A non-decomposable model
##
mf <- ~me:ve+ve:al+al:an+an:me
object <- cmod(mf, data=math)
testadd(object,c("me","al"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{htest}% __ONLY ONE__ keyword per line
