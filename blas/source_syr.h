/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

  size_t i, j;
  size_t ix, jx;

  if(Uplo == CblasUpper) {
    jx = 0;
    for(j=0; j<N; j++) {
      BASE_TYPE tmp = alpha * X[jx];
      ix = jx;
      for(i=j; i<N; i++) {
        A[lda*j + i] += X[ix] * tmp;
        ix += incX;
      }
      jx += incX;
    }
  }
  else {
    jx = 0;
    for(j=0; j<N; j++) {
      BASE_TYPE tmp = alpha * X[jx];
      ix = 0;
      for(i=0; i<=j; i++) {
        A[lda*j + i] += X[ix]*tmp;
        ix += incX;
      }
      jx += incX;
    }
  }
