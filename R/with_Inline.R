
.src.colmult2 <- 
'/* Calculates the product of matrices X and Y */
  int nprot=0;
  PROTECT(v = coerceVector(v, REALSXP)); nprot++;
  PROTECT(M = coerceVector(M, REALSXP)); nprot++;	
  double *vptr;    vptr = REAL(v);
  double *mptr;    mptr = REAL(M);

  int nrM = INTEGER(GET_DIM(M))[0];
  int ncM = INTEGER(GET_DIM(M))[1];
  
  SEXP ans; PROTECT(ans = allocMatrix(REALSXP, nrM, ncM)); nprot++;
  double *ansptr; ansptr = REAL(ans);

  int ii,jj;
  
  for (jj=0; jj<ncM; jj++){
    for (ii=0; ii<nrM; ii++){
		ansptr[ii+nrM*jj] = vptr[jj] * mptr[ii+nrM*jj];
	}
  }
  
  UNPROTECT(nprot);
  return ans;
'

.colmult2 <- cfunction(signature(v = "vector", M="matrix"),
                       body = .src.colmult2, convention=".Call") 


.src.vMMt2 <- 
'/* Calculates the product of matrices X and Y */
  int nprot=0;
  PROTECT(v = coerceVector(v, REALSXP)); nprot++;
  PROTECT(M = coerceVector(M, REALSXP)); nprot++;	
  double *vptr;    vptr = REAL(v);
  double *mptr;    mptr = REAL(M);

  int nrM = INTEGER(GET_DIM(M))[0];
  int ncM = INTEGER(GET_DIM(M))[1];
  
  SEXP ans; PROTECT(ans = allocMatrix(REALSXP, nrM, nrM)); nprot++;
  double *ansptr; ansptr = REAL(ans);

  double sum;
  int ii,jj,kk;
  
  for (ii=0; ii<nrM; ii++){
    for (jj=ii; jj<nrM; jj++){
		sum=0;
		for (kk=0; kk<ncM; kk++){
		   sum = sum + mptr[ii+nrM*kk]*mptr[jj+nrM*kk]*vptr[kk];
		}
		ansptr[ii+nrM*jj] = sum;
		ansptr[jj+nrM*ii] = sum;
	}
  }
  
  UNPROTECT(nprot);
  return ans;
'

.vMMt2 <- cfunction(signature(v = "vector", M="matrix"),
                    body = .src.vMMt2, convention=".Call") 
