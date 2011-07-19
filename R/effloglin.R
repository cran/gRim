##
## Implementation of efficient IPS algorith for log-linear models
## Based on fitting model using a grain structure
##
## Argument names are chosen so as to match those of loglin()
##

## FIXME: cpt2() is a bad name, and it is crap that it needs a gmData object...

effloglin <- function(table, margin, fit=FALSE, eps=0.01, iter=20, print=TRUE){
  
  ## uu  <- ugListMAT(margin)
  uu  <- ugList(margin, result="matrix")
  vn  <- colnames(uu)
  tri <- triangulateMAT(uu)
  rip <- ripMAT(tri)

  cliq     <- rip$cliques
  len.cliq <- length(cliq)
  itcount  <- 1

  ## "host clique" for each generator
  ##
  ghost <- rep(NA, length(margin))
  seqcliq <- seq_along(cliq)
  for (kk in 1:length(margin)){
    ##cat("kk:", kk,"\n")
    gg <- margin[[kk]]
    for (ii in seqcliq){
      ##cat ("ii", ii, "\n")
      zz <- charmatch(gg, cliq[[ii]])
      if (!any(is.na(zz))){
        ghost[kk] <- ii
        break
      } 
    }
  }
    
  if (is.array(table)){
    Nobs <- sum(table)
    stlist  <- lapply(margin, function(xx)
                      {
                        tableMargin(table, xx)
                      }
                      )    
  } else {
    Nobs <- sum(table[[1]])
    stlist <- table
  }

  zzz       <- unlist(lapply(stlist, dimnames), recursive=FALSE)
  vl        <- zzz[uniquePrim(names(zzz))]
  pot.list  <- lapply(cliq, function(cq) parray(cq, levels=vl[cq], values=1, normalize="all"))

  cat("effloglin\n")
  print(as.data.frame.table(pot.list[[1]]))
  
  ## ## Potential list over cliques
  ## Clique marginals
  prob.list  <- propagateLS(pot.list, rip, initialize=TRUE)        
  
  logL <- 0
  zzz <- vector("numeric", length(margin))
  repeat{
    cat(sprintf("---------- iteration: %i -----------\n", itcount))
    for (ss in seq_along(margin)){
      gg  <- margin[[ss]]
      st  <- stlist[[ss]]
      cq  <- cliq[[ghost[ss]]]
      cq.idx  <- ghost[ss]      
      cpot    <- prob.list[[cq.idx]]
      ##adjust  <- tableOp(st, tableMargin(cpot, gg)*Nobs, "/")
      
      tm <- tableMargin(cpot, gg)*Nobs
      adjust  <- st / tm
      zzz[ss] <- max(abs(log(adjust)))
                                        #zzz[ss] <- max(abs(st-tm))
      logL <- logL + sum(st * log(adjust))
                                        #pot.list[[cq.idx]] <- tableOp(pot.list[[cq.idx]], adjust, "*")
      pot.list[[cq.idx]] <- tableOp2(pot.list[[cq.idx]], adjust, `*`)
      prob.list  <- propagateLS(pot.list, rip, initialize=TRUE)
    }
    
    if (print)
      cat("max deviation (obs-fitted):", max(zzz), "\n")
    if (max(zzz)<eps || itcount>=iter)
      break()

    itcount <- itcount + 1
  }
   
  vl <- unlist(lapply(stlist, dimnames), recursive=FALSE)[vn]
  nlev <- unlistPrim(lapply(vl, length))
  
  gn <- lapply(margin, charmatch, vn)
  nparm <- .loglinGenDim(gn, nlev)
  df    <- prod(nlev) - 1 - nparm
  
  ans <- list(potlist=pot.list, margin=margin, vn=vn, rip=rip, ghost=ghost, stlist=stlist, logL=logL, nparm=nparm, df=df)
  
  ## Create full joint:
  ##
  if (fit){
    pjoint <- prob.list[[1]]
    if (length(prob.list)>1){
      for (ii in 2:length(prob.list)){
        pjoint <- tableOp(pjoint, tableOp(prob.list[[ii]], tableMargin(prob.list[[ii]], rip$sep[[ii]]), "/"),"*")
      }
    }
    pjoint <- tablePerm(pjoint, vn)*Nobs
    
    ans <- c(ans, list(fit=pjoint))
  }  
  ## class(ans) <- "effloglin"
  return(ans)
}






## ## List of clique marginals
## cqlist  <- lapply(rip$cliques, function(x) cpt2(x, gmData=gmd)/Nobs)

## ## List of separator marginals
## splist  <- lapply(rip$separators[-1], function(x) cpt2(x, gmData=gmd)/Nobs)

## prev.logL <- sum(cqlist[[1]] * log(prob.list[[1]]))
## if (length(cliq)>1){
##   for (ii in 2:length(cliq))
##     prev.logL <- prev.logL + sum(cqlist[[ii]] * log(prob.list[[ii]])) - sum(splist[[ii-1]]* log(tableMargin(prob.list[[ii]], rip$sep[[ii]])))
## }

## logL <- sum(cqlist[[1]] * log(prob.list[[1]]))
## if (length(cliq)>1){
##   for (ii in 2:length(cliq))
##     logL <- logL + sum(cqlist[[ii]] * log(prob.list[[ii]])) - sum(splist[[ii-1]]* log(tableMargin(prob.list[[ii]], rip$sep[[ii]])))
## }
## d.logL <- logL-prev.logL
## cat(sprintf("logL: %8.5f diff logL: %7.6f\n", logL,d.logL))


## function (table, margin, start = rep(1, length(table)), fit = FALSE, 
##     eps = 0.1, iter = 20, param = FALSE, print = TRUE) 
## {
##     rfit <- fit
##     dtab <- dim(table)
##     nvar <- length(dtab)
##     ncon <- length(margin)
##     conf <- matrix(0, nrow = nvar, ncol = ncon)
##     nmar <- 0
##     varnames <- names(dimnames(table))
##     for (k in seq_along(margin)) {
##         tmp <- margin[[k]]
##         if (is.character(tmp)) {
##             tmp <- match(tmp, varnames)
##             margin[[k]] <- tmp
##         }
##         conf[1:length(tmp), k] <- tmp
##         nmar <- nmar + prod(dtab[tmp])
##     }
##     ntab <- length(table)
##     if (length(start) != ntab) 
##         stop("'start' and 'table' must be same length")
##     storage.mode(conf) <- "integer"
##     z <- .C("loglin", as.integer(nvar), as.integer(dtab), as.integer(ncon), 
##         conf, as.integer(ntab), as.double(table), fit = as.double(start), 
##         locmar = integer(ncon), as.integer(nmar), marginals = double(nmar), 
##         as.integer(ntab), u = double(ntab), as.double(eps), as.integer(iter), 
##         dev = double(iter), nlast = integer(1), ifault = integer(1), 
##         PACKAGE = "stats")
##     switch(z$ifault, stop("this should not happen"), stop("this should not happen"), 
##         warning("algorithm did not converge"), stop("incorrect specification of 'table' or 'start'"))
##     if (print) 
##         cat(z$nlast, "iterations: deviation", z$dev[z$nlast], 
##             "\n")
##     fit <- z$fit
##     attributes(fit) <- attributes(table)
##     observed <- as.vector(table[start > 0])
##     expected <- as.vector(fit[start > 0])
##     pearson <- sum((observed - expected)^2/expected)
##     observed <- as.vector(table[table * fit > 0])
##     expected <- as.vector(fit[table * fit > 0])
##     lrt <- 2 * sum(observed * log(observed/expected))
##     subsets <- function(x) {
##         y <- list(vector(mode(x), length = 0))
##         for (i in seq_along(x)) {
##             y <- c(y, lapply(y, "c", x[i]))
##         }
##         y[-1]
##     }
##     df <- rep.int(0, 2^nvar)
##     for (k in seq_along(margin)) {
##         terms <- subsets(margin[[k]])
##         for (j in seq_along(terms)) df[sum(2^(terms[[j]] - 1))] <- prod(dtab[terms[[j]]] - 
##             1)
##     }
##     if (!is.null(varnames) && all(nzchar(varnames))) {
##         for (k in seq_along(margin)) margin[[k]] <- varnames[margin[[k]]]
##     }
##     else {
##         varnames <- as.character(1:ntab)
##     }
##     y <- list(lrt = lrt, pearson = pearson, df = ntab - sum(df) - 
##         1, margin = margin)
##     if (rfit) 
##         y$fit <- fit
##     if (param) {
##         fit <- log(fit)
##         terms <- seq_along(df)[df > 0]
##         parlen <- length(terms) + 1
##         parval <- list(parlen)
##         parnam <- character(parlen)
##         parval[[1]] <- mean(fit)
##         parnam[1] <- "(Intercept)"
##         fit <- fit - parval[[1]]
##         dyadic <- NULL
##         while (any(terms > 0)) {
##             dyadic <- cbind(dyadic, terms%%2)
##             terms <- terms%/%2
##         }
##         dyadic <- dyadic[order(rowSums(dyadic)), ]
##         for (i in 2:parlen) {
##             vars <- which(dyadic[i - 1, ] > 0)
##             parval[[i]] <- apply(fit, vars, mean)
##             parnam[i] <- paste(varnames[vars], collapse = ".")
##             fit <- sweep(fit, vars, parval[[i]], check.margin = FALSE)
##         }
##         names(parval) <- parnam
##         y$param <- parval
##     }
##     return(y)
## }
