##############################################################
##
## Fit graphical Gaussian model
##
## FIXME: Move to gRbase package
##
##############################################################

#' @title Iterative proportional fitting of graphical Gaussian model
#' 
#' @description Fit graphical Gaussian model by iterative proportional fitting.
#' 
#' @details \code{ggmfit} is based on a C implementation.  \code{ggmfitr} is
#'     implemented purely in R (and is provided mainly as a benchmark for the
#'     C-version).
#' 
#' @aliases ggmfit ggmfitr
#' @param S Empirical covariance matrix
#' @param n.obs Number of observations
#' @param glist Generating class for model (a list)
#' @param start Initial value for concentration matrix
#' @param eps Convergence criterion
#' @param iter Maximum number of iterations
#' @param details Controlling the amount of output.
#' @param ... Optional arguments; currently not used
#' @return A list with \item{lrt}{Likelihood ratio statistic (-2logL)}
#'     \item{df}{Degrees of freedom} \item{logL}{log likelihood}
#'     \item{K}{Estimated concentration matrix (inverse covariance matrix)}
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{cmod}}, \code{\link{loglin}}
#' @keywords multivariate models
#' @examples
#' 
#' ## Fitting "butterfly model" to mathmark data
#' ## Notice that the output from the two fitting functions is not
#' ## entirely identical.
#' data(math)
#' ddd <- cov.wt(math, method="ML")
#' glist <- list(c("al", "st", "an"), c("me", "ve", "al"))
#' ggmfit (ddd$cov, ddd$n.obs, glist)
#' ggmfitr(ddd$cov, ddd$n.obs, glist)
#' 
#' @export ggmfit
ggmfit <- function(S, n.obs, glist, start=NULL, 
                   eps=1e-12, iter=1000, details=0, ...)
{

  #cat("ggmfit:\n"); print(S); print(n.obs); print(glist)
  
  glist.save <- glist.num <- glist
  data.vn    <- colnames(S)

  ## The used variables
  usevar <- unique.default(unlist(glist))       

  ## Check that the used variables are in S
  zzz <- match(usevar, data.vn)
  if (any(is.na(zzz)))
    stop("Variables ", usevar[is.na(zzz)], " not in data\n")

  ## Get variables in the right order: The order in the S
  usevar  <-  data.vn[sort(zzz)]

  ## Possibly consider only submatrix of S
  S <- S[usevar, usevar, drop=FALSE]
  data.vn <- colnames(S)

  vn <- seq_along(data.vn)
  nvar <- length(vn)
  
  ## Numerical (indices) representation of glist
  glist.num <- lapply(glist, match, data.vn)

  glen    <- sapply(glist.num, length)
  ng      <- length(glist.num)

  clist.num   <- lapply(glist.num, function(x) vn[-x])  
  clen    <- sapply(clist.num,length)

  gg <- as.integer(unlist(glist.num)-1)# print(gg); 
  cc <- as.integer(unlist(clist.num)-1)# print(cc)    

  if (is.null(start)){
    start <- diag(1/diag(S))   #print(start)
  }

  ## dyn.load("ggmfit.dll")
  xxx<-.C("Cggmfit", S=S, n=as.integer(n.obs), K=start, nvar=nvar, ngen=ng, 
          glen=glen, glist=gg, clen=clen, clist=cc, 
          logL=numeric(1), eps=as.numeric(eps),
          iter=as.integer(iter), converged=as.integer(1),
          details=as.integer(details),
          PACKAGE="gRim")
  ## dyn.unload("ggmfit.dll")

  xxx             <- xxx[c("logL", "K", "iter")]  

  dimnames(xxx$K) <- dimnames(S)
  detK  <- det(xxx$K)
  dev   <- -n.obs * log(det(S %*% xxx$K))            ## deviance to the saturated model  
  df    <-  sum(xxx$K==0) / 2

  ans  <- list(dev=dev, df=df, detK=detK, nvar=nvar,S=S,n.obs=n.obs)
  ans   <- c(ans, xxx)
  
  return(ans)  
}





#' @export  
ggmfitr <- function(S, n.obs, glist, start=NULL, 
                    eps=1e-12, iter=1000, details=0, ...)
{

    ##
## Calculate logL for N(0,\Sigma) model.
##
## Sigma = Covariance matrix parameter
## K     = Sigma inverse
## S     = sample covariance matrix
## n     = sample size
##
ell <- function(Sigma, S, n){

    shdet <- function(Sigma){
        prod(eigen(Sigma)[[1]])
    }
    p <- dim(S)[1]
    const <- -n * p/2 * log(2 * pi)
    const - n/2 * log(shdet(Sigma)) - n/2 * sum(diag( solve(Sigma) %*% S )) 
}

ellK <- function (K, S, n)
{
    value <- (n/2) * (log(det(K)) - sum(rowSums(K * S)))
    value
}


    
  if (is.null(start)){
    K     <- diag(1/diag(S))
  } else {
    K     <- start
  }

  dimnames(K)<-dimnames(S)
  vn <- colnames(S); #print(vn)

  x <- lapply(glist, match, vn)
  
  varIndex=1:nrow(K)
  itcount=0

  if (length(x)){
    my.complement <- function(C) return(setdiff(varIndex,C))
    x.complements <- lapply(x, my.complement)
                                        #print("x"); print(x)
                                        #print("x.comp");print(x.complements)
    
    if(length(x.complements[[1]])==0){
      return(list(K=solve(S)))
    }
    logLvec <- NULL
    repeat {    

      for(j in 1:length(x)){
        C     <- x[[j]]
        notC  <- x.complements[[j]]
        #print(C); print(S[C,C,drop=FALSE])
        K[C,C] <- solve( S[C,C,drop=FALSE] ) +
          K[C,notC,drop=FALSE]%*%solve(K[notC,notC,drop=FALSE])%*%K[notC,C,drop=FALSE]
        #print(K)
      }
      logL <- ellK(K,S,n.obs)
      logLvec <- c(logLvec, logL)
      itcount <- itcount+1
      if (itcount>1){
        if (logL - prevlogL < eps){  
          converged=TRUE
          break
        }
      } else {
        if (itcount==iter){
          converged=FALSE
          break
        } 
      }
      prevlogL <- logL
    }    
  }

  df <- sum(K[upper.tri(K)] == 0)
  #ans <- list(K=K, logL=logL, converged=converged, itcount=itcount)
  ans <- list(dev=-2*logL, df=df, logL=logL, K=K, S=S,n.obs=n.obs,
              itcount=itcount, converged=converged,logLvec=logLvec)
  return(ans)
}

