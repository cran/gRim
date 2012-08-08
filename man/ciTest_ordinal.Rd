\name{ciTest_ordinal}
\alias{ciTest_ordinal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to compute Monte Carlo and asymptotic tests of conditional independence for ordinal and/or nominal variables.
}
\description{
The function computes tests of independence of two variables, say u and v, given
a set of variables, say S. The deviance, Wilcoxon, Kruskal-Wallis and Jonkheere-Terpstra tests are supported. 
Asymptotic and Monte Carlo p-values are computed.
}
\usage{
ciTest_ordinal(x, set = NULL, statistic = "dev", N = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    A dataframe or table.
  }
  \item{set}{
    The variable set (u,v,S), given either as an integer vector of the column numbers of a dataframe or dimension numbers
    of a table, or as a character vector with the corresponding variable or dimension names.
  }
  \item{statistic}{
    Either "deviance", "wilcoxon", "kruskal" or "jt".
  }
  \item{N}{
    The number of Monte Carlo samples. If N<=0 then Monte Carlo p-values are not computed.
  }
  \item{\dots}{
    Additional arguments, currently not used
  }
}
\details{
  The deviance test is appropriate when u and v are nominal;
  Wilcoxon, when u is binary and v is ordinal;
  Kruskal-Wallis, when u is nominal and v is ordinal;
  Jonckheere-Terpstra, when both u and v are ordinal.
}
\value{
A list including the test statistic, the asymptotic p-value and, when computed, the Monte Carlo p-value.  
\item{P}{Asymptotic p-value}
\item{montecarlo.P}{Monte Carlo p-value}
}
\references{
See Edwards D. (2000), "Introduction to Graphical Modelling", 2nd ed., Springer-Verlag, pp. 130-153.
}
\author{
Flaminia Musella, David Edwards, Søren Højsgaard, \email{sorenh@math.aau.dk}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{ciTest_table}},   \code{\link{ciTest}}
}
\examples{
library(gRim)
data(dumping)
ciTest_ordinal(dumping, c(2,1,3),stat="jt",N=1000)
ciTest_ordinal(dumping, c("Operation","Symptom","Centre"),stat="jt",N=1000)
ciTest_ordinal(dumping, ~Operation+Symptom+Centre,stat="jt",N=1000)

data(reinis)
ciTest_ordinal(reinis, c(1,3,4:6),N=1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
