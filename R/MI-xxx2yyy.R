#####################################################################
####
#### Conversion between different parametrizations of mixed models
####
#####################################################################

pms2ghkParms <- function(parms){
  parms <- unclass(parms)  
  switch(parms[['gentype']],
         "discrete"={
           res    <- list(g=log(parms[["p"]]),h=NULL, K=NULL, gentype="discrete")
         },
         "mixed"={
           KK     <- .cholsolve(parms[["Sigma"]])           
           Q      <- nrow(KK)
           #.logdetSig <- -log(det(KK))
           .logdetSig <- -c(determinant.matrix(KK)[['modulus']])

           mu     <- parms[["mu"]]
           hh     <- KK %*% mu # h = Sigma.inv %*% mu
           quad   <- colSumsPrim(hh * mu)
           gg     <- log(parms[["p"]]) + (- .logdetSig - Q*log(2*pi) - quad) / 2
           res    <- list(g=gg, h=hh, K=KK, gentype="mixed")
         },
         "continuous"={
           KK     <- .cholsolve(parms[["Sigma"]])
           Q      <- nrow(KK)
           #detSig <- 1/det(KK) ##FIXME: brug det.matrix isf.
           detSig <- 1/c(determinant.matrix(KK, logarithm=FALSE)[['modulus']])
           mu     <- parms[["mu"]]           
           hh     <- KK %*% mu # h = Sigma.inv %*% mu
           quad   <- colSumsPrim(hh * mu)
           gg     <- (- log(detSig) - Q*log(2*pi) - quad) / 2
           res    <- list(g=gg, h=hh, K=KK,gentype="continuous")
         })
  
  val <- c(res, parms[-(1:4)])
  class(val)<- c("ghk","MIparms")
  return(val)
}

pms2phkParms <- function(parms){
  parms <- unclass(parms)  
  switch(parms[['gentype']],
         "discrete"={
           res    <- list(p=parms[["p"]],h=NULL, K=NULL, gentype="discrete")
         },
         "mixed"={
           KK     <- .cholsolve(parms[["Sigma"]])           
           hh     <- KK %*% parms[["mu"]] # h = Sigma.inv %*% mu
           res    <-list(p=parms[['p']], h=hh, K=KK, gentype="mixed")
         },
         "continuous"={
           KK     <- .cholsolve(parms[["Sigma"]])
           hh     <- KK %*% parms[["mu"]] # h = Sigma.inv %*% mu           
           res    <- list(p=parms[['p']], h=hh, K=KK,gentype="continuous")
         })
  
  val <- c(res, parms[-(1:4)])
  class(val)<- c("phk","MIparms")
  return(val)
}


phk2ghkParms <- function(parms){
  parms <- unclass(parms)    
  switch(parms[['gentype']],
         "discrete"={
           res    <- list(g=log(parms[["p"]]),h=NULL, K=NULL, gentype="discrete")
         },
         "mixed"={
           KK     <- parms[["K"]]           
           Q      <- nrow(KK)
           detSig <- 1/det(KK)
           hh     <- parms[["h"]]
           mu     <- .cholsolve(KK) %*% hh # K.inv %*% h
           quad   <- colSumsPrim(hh * mu)

           gg     <- log(parms[["p"]]) + (- log(detSig) - Q*log(2*pi) - quad) / 2
           res    <-list(g=gg, h=hh, K=KK, gentype="mixed")
         },
         "continuous"={
           KK     <- parms[["K"]]           
           Q      <- nrow(KK)
           detSig <- 1/det(KK)
           mu     <- .cholsolve(KK) %*% hh # K.inv %*% h
           quad   <- colSumsPrim(hh * mu)

           gg     <- (- log(detSig) - Q*log(2*pi) - quad) / 2
           res    <-list(g=gg, h=hh, K=KK, gentype="continuous")
         })
  
  val <- c(res, parms[-(1:4)])
  class(val)<- c("ghk","MIparms")
  return(val)
}

ghk2pmsParms<-function(parms){

  parms <- unclass(parms)
  switch(parms[['gentype']],
         "discrete"={
           zzz <- parms[['g']]
           pp  <- exp(zzz-mean.default(zzz))
           res <- list(p=pp/sum(pp), mu=NULL, Sigma=NULL, gentype="discrete")
         },
         "mixed"={
           Sigma     <- .cholsolve(parms[['K']])
           hh        <- parms[['h']]
           mu        <- Sigma %*% hh        # Kinv %*% h
           g.quad    <- parms[['g']] + colSumsPrim(hh * mu)/2
           pp        <- exp( g.quad - mean.default(g.quad) )           
           res       <- list(p=pp/sum(pp), mu=mu, Sigma=Sigma, gentype="mixed")  
         },
         "continuous"={
           Sigma     <- .cholsolve(parms[['K']])
           mu        <- Sigma %*% parms[['h']]
           res       <- list(p=1, mu=mu, Sigma=Sigma, gentype="continuous")
         })

  val <- c(res, parms[-(1:4)])
  class(val)<- c("pms","MIparms")
  return(val)
}



phk2pmsParms<-function(parms){
  parms <- unclass(parms)  
  switch(parms[['gentype']],
         "discrete"={
           res <- list(p=parms[['g']], mu=NULL, Sigma=NULL, gentype="discrete")
         },
         "mixed"={
           Sigma     <- .cholsolve(parms[['K']])
           mu        <- Sigma %*% parms[['h']] # Kinv %*% h
           res       <- list(p=parms[['p']], mu=mu, Sigma=Sigma, gentype="mixed")  
         },
         "continuous"={
           Sigma     <- .cholsolve(parms[['K']])
           mu        <- Sigma %*% parms[['h']]
           res       <- list(p=1, mu=mu, Sigma=Sigma, gentype="continuous")
         })

  val <- c(res, parms[-(1:4)])
  class(val)<- c("pms","MIparms")
  return(val)
}

ghk2phkParms<-function(parms){
  parms <- unclass(parms)  
  switch(parms[['gentype']],
         "discrete"={
           zzz <- parms[['g']]
           pp  <- exp(zzz-mean.default(zzz)) 
           res <- list(p=pp/sum(pp), mu=NULL, Sigma=NULL, gentype="discrete")
         },
         "mixed"={
           Sigma     <- .cholsolve(parms[['K']])
           hh        <- parms[['h']]
           mu        <- Sigma %*% hh # Kinv %*% h
           g.quad    <- parms[['g']] + colSumsPrim(hh * mu)/2
           pp        <- exp( g.quad - mean.default(g.quad))           
           res       <- list(p=pp/sum(pp), h=hh, K=parms[['K']], gentype="mixed")  
         },
         "continuous"={
           Sigma     <- .cholsolve(parms[['K']])
           mu        <- Sigma %*% parms[['h']]
           res       <- list(p=1, h=parms[['h']], K=parms[['K']], gentype="continuous")
         })

  val <- c(res, parms[-(1:4)])
  class(val)<- c("phk","MIparms")
  return(val)
}








### Normalizes ghK representation
###

.normalize.ghkParms <- function(parms){

  K.idx <- 3
  h         <- parms[['h']]
  #mu        <- solve.default(parms[['K']], h)
  mu        <- .cholsolve(parms[[K.idx]]) %*% h
  logdetK   <- .logdet(parms[[K.idx]])

  Q         <- nrow(parms[[K.idx]])

  quad   <- colSumsPrim(h * mu)
  zzz    <- parms[['g']] + quad / 2
  ppp    <- exp( zzz - mean.default(zzz))
  normconst     <- sum(ppp)
  pppn    <- ppp / normconst

  g.new <- log(pppn) + (logdetK - Q*log(2*pi) - quad)/2

  parms[['g']] <- g.new
  parms
}




CGstats2mmodParms <- function(parms, type="ghk"){
  type <- match.arg(type, c("ghk","pms"))

  ans <- switch(type,
                "pms"={moment2pmsParms(parms)},
                "ghk"={pms2ghkParms(moment2pmsParms(parms))}
                )                
  return(ans)
}

moment2pmsParms <- function(SS){
  parms <- list(p=SS$n.obs/sum(SS$n.obs), mu=SS$center, Sigma=SS$cov, n.total=sum(SS$n.obs), gentype="mixed")
  #ans   <- c(list(parms=parms), SS[-(1:3)])
  ans   <- c(parms, SS[-(1:3)])
  class(ans)<-c("pms","MIparms")
  ans
}
