\name{modify_glist}
\alias{modify_glist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Modify generating class for a graphical/hierarchical model
}
\description{
    Modify generating class for a graphical/hierarchical model by 1)
    adding edges, 2) deleting edges, 3) adding terms and 4) deleting terms.
}
\usage{
modify_glist(glist, items, details = 0)
}
\arguments{
  \item{glist}{
    A list of vectors where each vector is a generator of the model.
}
\item{items}{
  A list with edges / terms to be added and deleted. See section
  'details' below.
}
  \item{details}{
    Control the amount of output (for debugging purposes).
}
}
\details{
  
  The \code{items} is a list with named entries as
  \code{list(add.edge=, drop.edge=, add.term=, drop.term=)}

  Not all entries need to be in the list. The corresponding actions are
  carried out in the order in which they appear in the list.

  See section 'examples' below for examples. 
  
  Notice that the operations do not in general commute: Adding an edge
  which is already in a generating class and then removing the edge
  again does not give the original generating class.
}
\value{
  A generating class for the modified model. The elements of the list
  are character vectors. 
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
  \code{\link{cmod}}
  \code{\link{dmod}}
  \code{\link{mmod}}
}
\examples{
glist <- list(c(1,2,3),c(2,3,4))

## Add edges
modify_glist(glist, items=list(add.edge=c(1,4)))
modify_glist(glist, items=list(add.edge=~1:4))

## Add terms
modify_glist(glist, items=list(add.term=c(1,4)))
modify_glist(glist, items=list(add.term=~1:4))

## Notice: Only the first term is added as the second is already 
## in the model.
modify_glist(glist, items=list(add.term=list(c(1,4),c(1,3))))
modify_glist(glist, items=list(add.term=~1:4+1:3))

## Notice: Operations are carried out in the order given in the
## items list and hence we get different results: 
modify_glist(glist, items=list(drop.edge=c(1,4), add.edge=c(1,4)))
modify_glist(glist, items=list(add.edge=c(1,4), drop.edge=c(1,4)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

