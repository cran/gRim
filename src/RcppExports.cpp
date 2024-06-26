// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/gRim.h"
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fit2way_
NumericVector fit2way_(const NumericVector& tab1, const NumericVector& tab2, const CharacterVector& R, const CharacterVector& vn);
RcppExport SEXP _gRim_fit2way_(SEXP tab1SEXP, SEXP tab2SEXP, SEXP RSEXP, SEXP vnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type tab1(tab1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type tab2(tab2SEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type vn(vnSEXP);
    rcpp_result_gen = Rcpp::wrap(fit2way_(tab1, tab2, R, vn));
    return rcpp_result_gen;
END_RCPP
}
// inv_qr_
mat inv_qr_(mat& X);
RcppExport SEXP _gRim_inv_qr_(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(inv_qr_(X));
    return rcpp_result_gen;
END_RCPP
}
// conips_ggm_
List conips_ggm_(arma::mat& S, List& elst, umat& emat, int& nobs, arma::mat K, int& maxit, double& eps, int& convcrit, int& print, List& aux);
RcppExport SEXP _gRim_conips_ggm_(SEXP SSEXP, SEXP elstSEXP, SEXP ematSEXP, SEXP nobsSEXP, SEXP KSEXP, SEXP maxitSEXP, SEXP epsSEXP, SEXP convcritSEXP, SEXP printSEXP, SEXP auxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< List& >::type elst(elstSEXP);
    Rcpp::traits::input_parameter< umat& >::type emat(ematSEXP);
    Rcpp::traits::input_parameter< int& >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K(KSEXP);
    Rcpp::traits::input_parameter< int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int& >::type convcrit(convcritSEXP);
    Rcpp::traits::input_parameter< int& >::type print(printSEXP);
    Rcpp::traits::input_parameter< List& >::type aux(auxSEXP);
    rcpp_result_gen = Rcpp::wrap(conips_ggm_(S, elst, emat, nobs, K, maxit, eps, convcrit, print, aux));
    return rcpp_result_gen;
END_RCPP
}
// covips_ggm_
List covips_ggm_(mat& S, List& elst, umat& emat, int& nobs, mat& K, int& maxit, double& eps, int& convcrit, int& print, List& aux);
RcppExport SEXP _gRim_covips_ggm_(SEXP SSEXP, SEXP elstSEXP, SEXP ematSEXP, SEXP nobsSEXP, SEXP KSEXP, SEXP maxitSEXP, SEXP epsSEXP, SEXP convcritSEXP, SEXP printSEXP, SEXP auxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< List& >::type elst(elstSEXP);
    Rcpp::traits::input_parameter< umat& >::type emat(ematSEXP);
    Rcpp::traits::input_parameter< int& >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< mat& >::type K(KSEXP);
    Rcpp::traits::input_parameter< int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int& >::type convcrit(convcritSEXP);
    Rcpp::traits::input_parameter< int& >::type print(printSEXP);
    Rcpp::traits::input_parameter< List& >::type aux(auxSEXP);
    rcpp_result_gen = Rcpp::wrap(covips_ggm_(S, elst, emat, nobs, K, maxit, eps, convcrit, print, aux));
    return rcpp_result_gen;
END_RCPP
}
// ncd_ggm_
List ncd_ggm_(mat& S, List& elst, umat& emat, int& nobs, mat K, int maxit, double& eps, int& convcrit, int print, List& aux);
RcppExport SEXP _gRim_ncd_ggm_(SEXP SSEXP, SEXP elstSEXP, SEXP ematSEXP, SEXP nobsSEXP, SEXP KSEXP, SEXP maxitSEXP, SEXP epsSEXP, SEXP convcritSEXP, SEXP printSEXP, SEXP auxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< List& >::type elst(elstSEXP);
    Rcpp::traits::input_parameter< umat& >::type emat(ematSEXP);
    Rcpp::traits::input_parameter< int& >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< mat >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int& >::type convcrit(convcritSEXP);
    Rcpp::traits::input_parameter< int >::type print(printSEXP);
    Rcpp::traits::input_parameter< List& >::type aux(auxSEXP);
    rcpp_result_gen = Rcpp::wrap(ncd_ggm_(S, elst, emat, nobs, K, maxit, eps, convcrit, print, aux));
    return rcpp_result_gen;
END_RCPP
}
// clone_
SEXP clone_(SEXP& x);
RcppExport SEXP _gRim_clone_(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(clone_(x));
    return rcpp_result_gen;
END_RCPP
}
// ggm_logL_
double ggm_logL_(mat& S, mat& K, int nobs);
RcppExport SEXP _gRim_ggm_logL_(SEXP SSEXP, SEXP KSEXP, SEXP nobsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< mat& >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    rcpp_result_gen = Rcpp::wrap(ggm_logL_(S, K, nobs));
    return rcpp_result_gen;
END_RCPP
}
// parm_ghk2pms_
List parm_ghk2pms_(List parms);
RcppExport SEXP _gRim_parm_ghk2pms_(SEXP parmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type parms(parmsSEXP);
    rcpp_result_gen = Rcpp::wrap(parm_ghk2pms_(parms));
    return rcpp_result_gen;
END_RCPP
}
// C_pms2ghk
RcppExport SEXP C_pms2ghk(SEXP parms_);
RcppExport SEXP _gRim_C_pms2ghk(SEXP parms_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type parms_(parms_SEXP);
    rcpp_result_gen = Rcpp::wrap(C_pms2ghk(parms_));
    return rcpp_result_gen;
END_RCPP
}
// C_ghk2pms
RcppExport SEXP C_ghk2pms(SEXP parms_);
RcppExport SEXP _gRim_C_ghk2pms(SEXP parms_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type parms_(parms_SEXP);
    rcpp_result_gen = Rcpp::wrap(C_ghk2pms(parms_));
    return rcpp_result_gen;
END_RCPP
}
// parm_normalize_ghk_
List parm_normalize_ghk_(List parms);
RcppExport SEXP _gRim_parm_normalize_ghk_(SEXP parmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type parms(parmsSEXP);
    rcpp_result_gen = Rcpp::wrap(parm_normalize_ghk_(parms));
    return rcpp_result_gen;
END_RCPP
}

// validate (ensure exported C++ functions exist before calling them)
static int _gRim_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _gRim_RcppExport_registerCCallable() { 
    R_RegisterCCallable("gRim", "_gRim_RcppExport_validate", (DL_FUNC)_gRim_RcppExport_validate);
    return R_NilValue;
}
