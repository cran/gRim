\name{getEdges}
\alias{getEdges}
\alias{getEdges.list}
\alias{getEdges.graphNEL}
\alias{getEdges.matrix}

\alias{getInEdges}
\alias{getOutEdges}

\alias{getInEdgesMAT}
\alias{getOutEdgesMAT}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Find edges in a graph and edges not in an undirected graph.
}
\description{
  Returns the edges of a graph (or edges not in a graph) where the graph
  can be either a graphNEL object, a list of generators or an adjacency matrix.   
}
\usage{
getEdges(object, type = "unrestricted", ingraph=TRUE, discrete=NULL, ...)
\method{getEdges}{list}(object, type = "unrestricted", ingraph=TRUE, discrete=NULL, ...)
\method{getEdges}{graphNEL}(object, type = "unrestricted", ingraph=TRUE, discrete=NULL, ...)
\method{getEdges}{matrix}(object, type = "unrestricted", ingraph=TRUE, discrete=NULL, ...)

getInEdges(object, type = "unrestricted", discrete=NULL, ...)
getOutEdges(object, type = "unrestricted", discrete=NULL, ...)
getInEdgesMAT(adjmat, type = "unrestricted", discrete=NULL, ...)
getOutEdgesMAT(adjmat, type = "unrestricted", discrete=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    An object representing a graph; either a generator list,
    a graphNEL object or an adjacency matrix.
  }
  \item{type}{
    Either "unrestricted" or "decomposable"
  }
  \item{ingraph}{
    If TRUE the result is the edges in the graph; if FALSE the result is
    the edges not in the graph.
  }
  \item{discrete}{
    This argument is relevant only if \code{object} specifies a marked
    graph in which some vertices represent discrete variables and some
    represent continuous variables. 
  }
  
  \item{\dots}{
    Additional arguments; currently not used.
  }
  \item{adjmat}{
    An adjacency matrix
  }
}

\details{
  When \code{ingraph=TRUE}:
  If type="decomposable" then \code{getEdges()} returns those edges e for
  which the graph with e removed is decomposable.
  
  When \code{ingraph=FALSE}:
  Likewise, if type="decomposable" then \code{getEdges()} returns those
  edges e for which the graph with e added is decomposable.

  The functions \code{getInEdges()} and \code{getInEdges()} are just
  wrappers for calls to \code{getEdges()}.
  
  The workhorses are \code{getInEdgesMAT()} and \code{getOutEdgesMAT()} and these work on
  adjacency matrices.

  Regarding the argument \code{discrete}, please see the documentation
  of \code{\link[gRbase]{mcsmarked}}.
}
\value{
  A p * 2 matrix with edges.
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
  Søren Højsgaard, sorenh@math.aau.dk
}
\note{
  These functions work on undirected graphs. The behaviour is
  undocumented for directed graphs.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{edgeList}},
  \code{\link{nonEdgeList}}
  \code{\link{mcsmarked}}
}

\examples{
gg     <- ug(~a:b:d+a:c:d+c:e)
glist  <- maxClique(gg)[[1]]
adjmat <- as.adjMAT(gg)

#### On a glist
getEdges(glist)
getEdges(glist,type="decomposable")
# Deleting (a,d) would create a 4-cycle

getEdges(glist, ingraph=FALSE)
getEdges(glist,type="decomposable", ingraph=FALSE)
# Adding (e,b) would create a 4-cycle

#### On a graphNEL
getEdges(gg)
getEdges(gg,type="decomposable")
# Deleting (a,d) would create a 4-cycle

getEdges(gg, ingraph=FALSE)
getEdges(gg,type="decomposable", ingraph=FALSE)
# Adding (e,b) would create a 4-cycle

#### On an adjacency matrix
getEdges(adjmat)
getEdges(adjmat,type="decomposable")
# Deleting (a,d) would create a 4-cycle

getEdges(adjmat, ingraph=FALSE)
getEdges(adjmat,type="decomposable", ingraph=FALSE)
# Adding (e,b) would create a 4-cycle


## Marked graphs; vertices a,b are discrete; c,d are continuous
UG <- ug(~a:b:c+b:c:d)
disc <- c("a","b")
getEdges(UG)
getEdges(UG, discrete=disc)
## Above: same results; there are 5 edges in the graph

getEdges(UG, type="decomposable")
## Above: 4 edges can be removed and will give a decomposable graph
##(only removing the edge (b,c) would give a non-decomposable model)

getEdges(UG, type="decomposable", discrete=c("a","b"))
## Above: 3 edges can be removed and will give a strongly decomposable
## graph. Removing (b,c) would create a 4--cycle and removing (a,b)
## would create a forbidden path; a path with only continuous vertices
## between two discrete vertices.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
