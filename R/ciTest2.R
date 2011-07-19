ciTest <- function(x,set=NULL,...){
  UseMethod("ciTest")
}

ciTest.table <- function(x, set=NULL, ...){
  ciTest_table(x, set, ...)
}

ciTest.list <- function(x, set=NULL, ...){
  ciTest_mvn(x, set, ...)
}

ciTest.data.frame <- function(x, set=NULL, ...){
  ciTest_df(x, set, ...)
}


