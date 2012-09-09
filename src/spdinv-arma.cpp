/* File: spdinv-arma.cpp */
#include <Rcpp.h>
#include <RcppArmadillo.h>

RcppExport SEXP C_spdinv_arma ( SEXP X_ ){
  arma::mat X    = Rcpp::as<arma::mat>(X_);
  arma::mat Xinv = arma::inv( arma::sympd(X) );
  return(Rcpp::wrap(Xinv));
}

