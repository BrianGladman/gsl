
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
      for(i=N-1; i+1>=1; i--) {
        size_t k = (i*(i+3))/2;
        BASE_TYPE temp = ( nounit ? X[i*incX] * Ap[k] : X[i*incX] );
	const size_t j_min = ( i>1 ? i-1 : 0 );
	k -= N-i;
	for(j=j_min; j+1>=1; j--) {
          temp += Ap[k] * X[j * incX];
          k -= N-j;
        }
	X[i*incX] = temp;
      }
    }
    else {
      ix = 0;
      for(i=0; i<N; i++) {
        size_t k = (i*(i+3))/2;
        BASE_TYPE temp = ( nounit ? X[ix] * Ap[k] : X[ix] );
	k += i+1;
	jx = (i+1)*incX;
        for(j=i+1; j<N; j++) {
          temp += Ap[k] * X[jx];
          k += j+1;
	  jx += incX;
        }
        X[ix] = temp;
	ix += incX;
      }
    }
  }
