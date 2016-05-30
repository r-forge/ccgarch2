#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h> /* for dgemv etc. */
#include <R_ext/Lapack.h> /* for dpotrf etc. */

SEXP dcc_ll2(SEXP Rin, SEXP epsi){
  int i, j, ione = 1, info1, info2, nobs = Rf_nrows(Rin), 
    cordim = Rf_ncols(Rin), ndim= Rf_ncols(epsi), IPIV[ndim]; 
  double  one = 1.0, zero = 0.0, *rdetR, 
    *rR1, *rvecR, *rI, *rzz, *rz1, *rz2, *rz_row, *rtmp1, *rtmp2, *rll;
  
  SEXP ll, R1, I, z1, z2, z_row, detR, tmp1, tmp2, vecR, eps;
  
  eps = PROTECT(duplicate(epsi));
  vecR = PROTECT(duplicate(Rin));
  
  R1 = PROTECT(allocMatrix(REALSXP, ndim, ndim));
  I = PROTECT(allocMatrix(REALSXP, ndim, ndim));
  ll = PROTECT(allocVector(REALSXP, nobs));
  z2 = PROTECT(allocVector(REALSXP, ndim));
  z1 = PROTECT(allocVector(REALSXP, ndim));
  z_row = PROTECT(allocVector(REALSXP, ndim));
  detR = PROTECT(allocVector(REALSXP, ione));
  tmp1 = PROTECT(allocVector(REALSXP, ione));
  tmp2 = PROTECT(allocVector(REALSXP, ione));
  
  rI = REAL(I);
  rvecR = REAL(vecR);
  rR1 = REAL(R1);
  rzz = REAL(eps);
  rz1 = REAL(z1);
  rz2 = REAL(z2);
  rz_row = REAL(z_row);
  rll = REAL(ll);
  rdetR = REAL(detR);
  rtmp1 = REAL(tmp1);
  rtmp2 = REAL(tmp2);
  
  /* initializing vectors */
  for(i =0; i < ndim; i++){
    rz1[i] = 0.0;
    rz2[i] = 0.0;
    rz_row[i] = 0.0;
  }
  
  /* main loop */
  for(i=0; i< nobs; i++){
    rll[i] = 0.0;       /* initializing */
  for(j = 0; j < cordim; j++){
    rI[j] = 0.0;                    /* initializing the working matrix */
  rR1[j] = 0.0;
  rR1[j] = rvecR[i + j*nobs];     /* R1 at time t */
  }
  
  /* z at time t */
  for(j = 0; j < ndim; j++){
    rz1[j] = rzz[i + j*nobs];
  }
  F77_CALL(dcopy)(&ndim, rz1, &ione, rz2 ,&ione);    /* copy rz1 to rz2 */
  
  /* LU decomposition of rR1 */
  F77_CALL(dgetrf)(&ndim, &ndim, rR1, &ndim, IPIV, &info1);
  
  rdetR[0] = 0.0;   /* initializing */
  rtmp1[0] = 0.0;   /* initializing */
  rtmp2[0] = 0.0;   /* initializing */
  /* log-determinant of rR1. This is possible because rR1 is always positive definite. */
  for (j = 0; j < ndim; j++) {
    rdetR[0] += log(fabs(rR1[j*(ndim + 1)]));    /* fabs() is a basic function in C returning an absolute value of an argument. */
  }
  
  /* inverse of rR1 */
  F77_CALL(dgetri)(&ndim, rR1, &ndim, IPIV, rI, &ndim, &info2);      /* now rR1 is the inverse of rR1. rI is a working matrix */
  
  F77_CALL(dgemv)("N", &ndim, &ndim, &one, rR1, &ndim, rz1, &ione, &zero, rz_row, &ione);      /* rR1%*%rz1 = rz_row */
  rtmp1[0] = F77_CALL(ddot)(&ndim, rz_row, &ione, rz2, &ione);    /* inner product of rz_row and rz2 to make rz1%*%rR1%*%rz1*/
  rtmp2[0] = F77_CALL(ddot)(&ndim, rz2, &ione, rz2, &ione);       /* inner product of rz2 to make sum(rz2^2) */
  rll[i] =  0.5*(rdetR[0] + rtmp1[0] - rtmp2[0]);
  }
  UNPROTECT(11);
  return(ll); 
}

