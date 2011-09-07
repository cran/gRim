pFormula <- function (formula, varNames, marginal=NULL, interactions=NULL, v.sep = ":", g.sep = "+", ignore.power.value=FALSE) 
{

  check.listformula <- function(list.formula, used.var){
    if (any(is.na(pmatch(unlistPrim(list.formula), used.var, 
                         duplicates.ok = TRUE)))) 
      stop("An invalid variable specification has been found\n")
    list.formula <- lapply(list.formula, function(x) {
      ii <- pmatch(x, used.var)
      used.var[ii]
    })

    modnames <- uniquePrim(unlist(list.formula))
    if(any(is.na(match(modnames, used.var))))
      stop("Variables in model not contained in the variable set. Perhaps a problem with 'marginal'?")
    
    list.formula
  }

  set.interactions <- function(list.formula, interactions){

    zz <- lapply(list.formula, function(ss){
      if (length(ss)<=interactions){
        list(ss)
      } else {
        combn(ss, interactions, simplify=FALSE)
      }
    })
    zz <- removeRedundant(unlist(zz, recursive=FALSE))
    zz
  }
  
  list.formula2str.formula <- function(list.formula){
    tmp <- lapply(list.formula, paste, collapse = v.sep)
    str.formula <- paste(unlistPrim(tmp), collapse = g.sep, sep = "")
    str.formula
  }
    
  clformula <- class(formula)

  if (length(marginal)>0){
    used.var <- marginal
  } else {
    used.var <- varNames
  }
  
  switch(clformula,
         "formula"={
           ##cat("A formula is given\n")
           pow <- .extract.power(formula)
           if (is.numeric(pow)){
             ##cat("pow          : ", pow, "\n")
             ## A power formula
             if (ignore.power.value){
               if (pow>1)
                 pow <- -1
             }
             if (pow == -1){
               ## cat("The saturated model\n")
               list.formula <- list(used.var)
             } else {
               pow <- min(c(pow, length(used.var)))
               list.formula <- combn(used.var, pow, simplify=FALSE)
             }
             if (!is.null(interactions))
               list.formula <- set.interactions(list.formula, interactions)
           } else {
             ## A "proper formula"
             tmp <- as.formula(paste("~", pow))
             list.formula <- rhsFormula2list(tmp)
             list.formula <- check.listformula(list.formula, used.var)
             if (!is.null(interactions))
               list.formula <- set.interactions(list.formula, interactions)             
           }
         },
         "list"={
           list.formula <- check.listformula(formula, used.var)      
           if (!is.null(interactions))
             list.formula <- set.interactions(list.formula, interactions)
         },
         "graphNEL"={
           list.formula <- maxCliqueMAT(as.adjMAT(formula))[[1]]
           list.formula <- check.listformula(list.formula, used.var)
           if (!is.null(interactions))
             list.formula <- set.interactions(list.formula, interactions)
         },
         "matrix"={
           list.formula <- maxCliqueMAT(formula)[[1]]
           list.formula <- check.listformula(list.formula, used.var)
           if (!is.null(interactions))
             list.formula <- set.interactions(list.formula, interactions)

         })

  formula      <- list2rhsFormula(list.formula)           
  ##            cat("str.formula  :", str.formula, "\n")
  ##            cat("formula      : ");  print(formula)
  ##            cat("list.formula :\n"); print(list.formula)      
    
  value <- list(glist = list.formula,
                formula = formula,
                ## str.formula = str.formula,                 
                varNames = uniquePrim(unlist(list.formula)))
  return(value)
}

.extract.power <- function (fff) 
{
    mimf <- paste(as.formula(fff))[2]
    mimf.split <- unlist(strsplit(mimf, ""))
    if (length(grep("[:alpha:]", mimf)) > 0) {
        pow <- mimf
    }
    else {
        has.hat <- match("^", mimf.split)
        sub <- unlist(strsplit(mimf, "\\^"))
        if (!is.na(has.hat)) {
            pow <- ifelse(sub[2] == ".", -1, as.numeric(sub[2]))
        }
        else {
            pow <- length(unlist(strsplit(sub, "\\.")))
        }
    }
    return(pow)
}
