
logLik.iModel <- function(object,...)
  structure(object$fitinfo$logL, df=object$fitinfo$dimension["df"], class="logLik")

## Returns (df, AIC=-2logL + k df), so the objective is to mimimize this quantity
##
extractAIC.iModel <- function(fit, scale, k = 2, ...){
  c(fit$fitinfo$df, fit$fitinfo$lrt - 2*fit$fitinfo$dimension["df"])
}

summary.iModel <- function(object, ...){
  glist <- object$glist

  isg   <- object$isGraphical
  isd   <- object$isDecomposable
  cq    <- maxClique(ugList(glist))$maxCliques
  ans   <- structure(list(glist=glist, isGraphical=isg, isDecomposable=isd, cliques=cq),
                     class="iModelsummary")
  ans
}

print.iModelsummary <- function(x,...){
  cat(sprintf("is graphical=%s; is decomposable=%s\n", x$isGraphical, x$isDecomposable))
  cat("generators (glist):\n")
  str(x$glist, give.head=FALSE, comp.str=" ", no.list=TRUE)
  cat("EXPERIMENTAL: components: ", names(x),"\n")
  invisible(x)
}

.extractFIT <- function(object,...){
  c(object[[1]], object$df)
}


formula.iModel <- function(x,...){
	list2rhsFormula(x$glist)
}

terms.iModel <- function(x,...){
	x$glist
}


