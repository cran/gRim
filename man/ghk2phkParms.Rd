\name{ghk2phkParms}
\alias{ghk2phkParms}
\alias{ghk2pmsParms}
\alias{phk2ghkParms}
\alias{phk2pmsParms}
\alias{pms2ghkParms}
\alias{pms2phkParms}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Conversion between different parametrizations of mixed interaction models.
}
\description{
  Functions to convert between canonical parametrization (g,h,K), moment
  parametrization (p,m,S) and mixed parametrization (p,h,K).
}
\usage{
ghk2phkParms(parms)
ghk2pmsParms(parms)
phk2ghkParms(parms)
phk2pmsParms(parms)
pms2ghkParms(parms)
pms2phkParms(parms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parms}{Parameters of a mixed interaction model
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Parameters of a mixed interaction model.
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
