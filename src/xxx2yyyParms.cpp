#include <RcppArmadillo.h>
#include <string>

RcppExport SEXP C_pms2ghk ( SEXP parms_ ){
  using namespace arma;
  using namespace Rcpp;
  Rcpp::List parms(parms_);
  std::string gentype = as<std::string>(parms["gentype"]); 
  if (strcmp(gentype.c_str(),"discrete")==0){
    //Rcpp::Rcout << "generator: discrete" << std::endl;
    vec DD = as<vec>(parms["p"]);
    int ii, ndd = DD.n_elem;
    NumericVector gg1 = parms["p"];
    NumericVector gg2 = clone(gg1);
    for (ii=0; ii<ndd; ii++){ gg2(ii) = log(DD(ii)); }

    return List::create(Named("g", wrap(gg2)),  Named("h", R_NilValue),
			Named("K", R_NilValue), Named("gentype", "discrete"));
  } else {
      mat QQ = as<mat>(parms["Sigma"]);
      mat LL = as<mat>(parms["mu"]);
      mat QQ_out  = inv( sympd(QQ) );
      mat LL_out  = QQ_out * LL;
      int ii, ndd = LL_out.n_cols, Q =QQ.n_rows;
      double cst  = log(det(QQ_out)) - Q*log(2*datum::pi);

      if (strcmp(gentype.c_str(),"mixed")==0){
	//Rcpp::Rcout << "generator: mixed" << std::endl;      
	vec dd = vec(Q);
	vec DD = as<arma::vec>(parms["p"]);
	vec gg  = vec(ndd);

	for (ii=0; ii<ndd; ii++){
	  dd = LL.col(ii);
	  gg(ii)   = log(DD(ii)) + (cst - as_scalar(dd.t() * QQ_out * dd)) / 2;
	}
	NumericVector gg1 = parms["p"];
	NumericVector gg2 = clone(gg1);
	for (ii=0; ii<ndd; ii++){ gg2(ii) = gg(ii); }
	return List::create(Named("g", wrap(gg2)),    Named("h", wrap(LL_out)),
			    Named("K", wrap(QQ_out)), Named("gentype", gentype.c_str()));	  	
      } else { 
	//Rcpp::Rcout << "generator: continuous" << std::endl;      	
	double gg2;
	gg2 = (cst - as_scalar(LL.t() * QQ_out * LL))/2;
	return List::create(Named("g", wrap(gg2)),    Named("h", wrap(LL_out)),
			    Named("K", wrap(QQ_out)), Named("gentype", gentype.c_str()));	  
      }
  }
}



RcppExport SEXP C_ghk2pms ( SEXP parms_ ){
  using namespace arma;
  using namespace Rcpp;
  Rcpp::List parms(parms_);
  std::string gentype = as<std::string>(parms["gentype"]); 
  if (strcmp(gentype.c_str(),"discrete")==0){
    //Rcpp::Rcout << "generator: discrete" << std::endl;
    vec DD = as<arma::vec>(parms["g"]);
    int ii, ndd = DD.n_elem;
    rowvec DD_out = rowvec(ndd);
    double mm = mean(DD);
    for (ii=0; ii<ndd; ii++) {DD_out(ii) = exp(DD(ii)-mm);}
    double ss = sum(DD_out);
    for (ii=0; ii<ndd; ii++) {DD_out(ii) = DD_out(ii)/ss;}
    
    NumericVector gg1 = parms["g"];
    NumericVector gg2 = clone(gg1);
    for (ii=0; ii<ndd; ii++){ gg2(ii) = DD_out(ii); }

    return List::create(Named("p", wrap(gg2)),      Named("mu", R_NilValue),
			Named("Sigma", R_NilValue), Named("gentype", "discrete"));
  } else {
    mat QQ = as<arma::mat>(parms["K"]);
    mat LL = as<arma::mat>(parms["h"]);
    mat QQ_out = inv( sympd(QQ) );
    mat LL_out = QQ_out * LL;
    int ii, ndd = LL_out.n_cols, Q =QQ_out.n_rows;
    rowvec DD_out  = rowvec(ndd);

    if (strcmp(gentype.c_str(),"mixed")==0){
      //Rcpp::Rcout << "generator: mixed" << std::endl;      
      vec dd = vec(Q);
      arma::vec DD = as<arma::vec>(parms["g"]);

      vec quad = vec(ndd);
      
      for (ii=0; ii<ndd; ii++){
	dd   = LL.col(ii);
	quad(ii) = DD(ii) + as_scalar(dd.t() * QQ_out * dd)/2;
      }
      double mm = mean(quad);
      for (ii=0; ii<ndd; ii++){ DD_out(ii) = exp(quad(ii)-mm); }
      double ss = sum(DD_out);
      for (ii=0; ii<ndd; ii++) {DD_out(ii) = DD_out(ii)/ss;}

      NumericVector gg1 = parms["g"];
      NumericVector gg2 = clone(gg1);
      for (ii=0; ii<ndd; ii++){ gg2(ii) = DD_out(ii); }
      return Rcpp::List::create(Rcpp::Named("p", wrap(gg2)),
				Rcpp::Named("mu", wrap(LL_out)),
				Rcpp::Named("Sigma", wrap(QQ_out)),
				Rcpp::Named("gentype", gentype.c_str())	);	
    } else {
      //Rcpp::Rcout << "generator: continuous" << std::endl;      	
      return Rcpp::List::create(Rcpp::Named("p", 1),
				Rcpp::Named("mu", wrap(LL_out)),
				Rcpp::Named("Sigma", wrap(QQ_out)),
				Rcpp::Named("gentype", gentype.c_str())	);	
    } 
  }
}









