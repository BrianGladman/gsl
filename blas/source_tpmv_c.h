/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

#include "matrix_access.h"

  size_t i, j;
  size_t ix;
  const int nounit = ( Diag == CblasNonUnit );

  if(TransA == CblasNoTrans) {
    /* form  x:= A*x */

    if(Uplo == CblasUpper) {
      for(i=0; i<N; i++) {
        BASE_TYPE atmp_r = REAL(Ap, 1, M_PACKEDTRUP_INDEX(N, 0, i, i));
	BASE_TYPE atmp_i = IMAG(Ap, 1, M_PACKEDTRUP_INDEX(N, 0, i, i));
	BASE_TYPE temp_r;
	BASE_TYPE temp_i;
	if(nounit) {
	  temp_r = atmp_r*REAL(X,incX,i) - atmp_i*IMAG(X,incX,i);
	  temp_i = atmp_r*IMAG(X,incX,i) + atmp_i*REAL(X,incX,i);
	}
	else {
	  temp_r = REAL(X,incX,i);
	  temp_i = IMAG(X,incX,i);
	}
        for(j=i+1; j<N; j++) {
	  atmp_r = REAL(Ap, 1, M_PACKEDTRUP_INDEX(N, 0, i, j));
	  atmp_i = IMAG(Ap, 1, M_PACKEDTRUP_INDEX(N, 0, i, j));
          temp_r += atmp_r * REAL(X,incX,j) - atmp_i * IMAG(X,incX,j);
	  temp_i += atmp_r * IMAG(X,incX,j) + atmp_i * REAL(X,incX,j);
        }
	REAL(X,incX,i) = temp_r;
	IMAG(X,incX,i) = temp_i;
      }
    }
    else {
      ix = N-1;
      for(i=0; i<N; i++) {        
        BASE_TYPE atmp_r = REAL(Ap, 1, M_PACKEDTRLO_INDEX(N, 0, N-1-i, N-1-i));
	BASE_TYPE atmp_i = IMAG(Ap, 1, M_PACKEDTRLO_INDEX(N, 0, N-1-i, N-1-i));
	BASE_TYPE temp_r;
	BASE_TYPE temp_i;
	if(nounit) {
	  temp_r = atmp_r*REAL(X,incX,ix) - atmp_i*IMAG(X,incX,ix);
	  temp_i = atmp_r*IMAG(X,incX,ix) + atmp_i*REAL(X,incX,ix);
	}
	else {
	  temp_r = REAL(X,incX,ix);
	  temp_i = IMAG(X,incX,ix);
	}
	if(ix > 0) {
	  const size_t j_max = ix-1;
	  for(j=0; j<=j_max; j++) {
	    atmp_r  = REAL(Ap, 1, M_PACKEDTRLO_INDEX(N, 0, N-1-i, j));
	    atmp_i  = IMAG(Ap, 1, M_PACKEDTRLO_INDEX(N, 0, N-1-i, j));
            temp_r += atmp_r * REAL(X,incX,j) - atmp_i * IMAG(X,incX,j);
	    temp_i += atmp_r * IMAG(X,incX,j) + atmp_i * REAL(X,incX,j);
          }
	}
	REAL(X,incX,ix) = temp_r;
	IMAG(X,incX,ix) = temp_i;
	ix--;
      }
    }

  }
  else {
    /* form  x := A'*x */

    if(Uplo == CblasUpper) {
      ix = N-1;
      for(i=0; i<N; i++) {        
        BASE_TYPE atmp_r = REAL(Ap, 1, M_PACKEDTRUP_INDEX(N, 0, N-1-i, N-1-i));
	BASE_TYPE atmp_i = IMAG(Ap, 1, M_PACKEDTRUP_INDEX(N, 0, N-1-i, N-1-i));
	BASE_TYPE temp_r;
	BASE_TYPE temp_i;
	if(nounit) {
	  temp_r = atmp_r*REAL(X,incX,ix) - atmp_i*IMAG(X,incX,ix);
	  temp_i = atmp_r*IMAG(X,incX,ix) + atmp_i*REAL(X,incX,ix);
	}
	else {
	  temp_r = REAL(X,incX,ix);
	  temp_i = IMAG(X,incX,ix);
	}
	if(ix > 0) {
	  const size_t j_max = ix-1;
	  for(j=0; j<=j_max; j++) {
	    atmp_r  = REAL(Ap, 1, M_PACKEDTRUP_INDEX(N, 0, j, N-1-i));
	    atmp_i  = IMAG(Ap, 1, M_PACKEDTRUP_INDEX(N, 0, j, N-1-i));
            temp_r += atmp_r * REAL(X,incX,j) - atmp_i * IMAG(X,incX,j);
	    temp_i += atmp_r * IMAG(X,incX,j) + atmp_i * REAL(X,incX,j);
          }
	}
	REAL(X,incX,ix) = temp_r;
	IMAG(X,incX,ix) = temp_i;
	ix--;
      }
    }
    else {
      for(i=0; i<N; i++) {
        BASE_TYPE atmp_r = REAL(Ap, 1, M_PACKEDTRLO_INDEX(N, 0, i, i));
	BASE_TYPE atmp_i = IMAG(Ap, 1, M_PACKEDTRLO_INDEX(N, 0, i, i));
	BASE_TYPE temp_r;
	BASE_TYPE temp_i;
	if(nounit) {
	  temp_r = atmp_r*REAL(X,incX,i) - atmp_i*IMAG(X,incX,i);
	  temp_i = atmp_r*IMAG(X,incX,i) + atmp_i*REAL(X,incX,i);
	}
	else {
	  temp_r = REAL(X,incX,i);
	  temp_i = IMAG(X,incX,i);
	}
        for(j=i+1; j<N; j++) {
	  atmp_r = REAL(Ap, 1, M_PACKEDTRLO_INDEX(N, 0, j, i));
	  atmp_i = IMAG(Ap, 1, M_PACKEDTRLO_INDEX(N, 0, j, i));
          temp_r += atmp_r * REAL(X,incX,j) - atmp_i * IMAG(X,incX,j);
	  temp_i += atmp_r * IMAG(X,incX,j) + atmp_i * REAL(X,incX,j);
        }
	REAL(X,incX,i) = temp_r;
	IMAG(X,incX,i) = temp_i;
      }
    }
  }
