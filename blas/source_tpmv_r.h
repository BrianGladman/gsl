
#include "matrix_access.h"

  size_t i, j;
  size_t ix, jx;
  const int nounit = ( Diag == CblasNonUnit );

  if(TransA == CblasNoTrans) {
    /* form  x:= A*x */

    if(Uplo == CblasUpper) {
      ix = 0;
      for(i=0; i<N; i++) {
        BASE_TYPE atmp = M_PACKEDTRUP_ACCESS(Ap, N, 0, i, i);
        BASE_TYPE temp = ( nounit ? X[ix] * atmp : X[ix] );
	jx = (i+1)*incX;
        for(j=i+1; j<N; j++) {
	  atmp  = M_PACKEDTRUP_ACCESS(Ap, N, 0, i, j);
          temp += atmp * X[jx];
	  jx += incX;
        }
        X[ix] = temp;
	ix += incX;
      }
    }
    else {
      ix = (N-1)*incX;
      for(i=0; i<N; i++) {
        BASE_TYPE atmp = M_PACKEDTRLO_ACCESS(Ap, N, 0, N-1-i, N-1-i);
        BASE_TYPE temp = ( nounit ? X[ix] * atmp : X[ix] );
	if(i < N-1) {
	  const size_t j_max = N-2-i;
	  size_t jx = 0;
	  for(j=0; j<=j_max; j++) {
	    atmp  = M_PACKEDTRLO_ACCESS(Ap, N, 0, N-1-i, j);
            temp += atmp * X[jx];
	    jx += incX;
          }
	}
	X[ix] = temp;
	ix -= incX;
      }
    }

  }
  else {
    /* form  x := A'*x */

    if(Uplo == CblasUpper) {
      ix = (N-1)*incX;
      for(i=0; i<N; i++) {
        BASE_TYPE atmp = M_PACKEDTRUP_ACCESS(Ap, N, 0, N-1-i, N-1-i);
        BASE_TYPE temp = ( nounit ? X[ix] * atmp : X[ix] );
	if(i < N-1) {
	  const size_t j_max = N-2-i;
	  size_t jx = 0;
	  for(j=0; j<=j_max; j++) {
	    atmp  = M_PACKEDTRUP_ACCESS(Ap, N, 0, j, N-1-i);
            temp += atmp * X[jx];
	    jx += incX;
          }
	}
	X[ix] = temp;
	ix -= incX;
      }
    }
    else {
      ix = 0;
      for(i=0; i<N; i++) {
        BASE_TYPE atmp = M_PACKEDTRLO_ACCESS(Ap, N, 0, i, i);
        BASE_TYPE temp = ( nounit ? X[ix] * atmp : X[ix] );
	jx = (i+1)*incX;
        for(j=i+1; j<N; j++) {
	  atmp  = M_PACKEDTRLO_ACCESS(Ap, N, 0, j, i);
          temp += atmp * X[jx];
	  jx += incX;
        }
        X[ix] = temp;
	ix += incX;
      }
    }
  }
