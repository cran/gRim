
print.citest <- function(x,...){
  if (length(x$varNames) > 2){
    cat("Testing", x$varNames[1], "_|_", x$varNames[2], "|",x$varNames[-(1:2)],"\n")
  } else {
    cat("Testing", x$varNames[1], "_|_", x$varNames[2], "\n")
  }
  cat(sprintf("Statistic (%s): %8.3f df: %s p-value: %6.4f method: %s\n",
              x$statname, x$statistic, x$df, x$p.value, x$method))
}

summary.citest <- function(object,...){
  str(object,max.level=1)
  return(invisible(object))
}
