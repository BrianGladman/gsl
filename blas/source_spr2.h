/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

  size_t i, j;
  size_t k = 0;

  if(Uplo == CblasUpper) {
    for(i=0; i<N; i++) {
      for(j=i; j<N; j++) {
        A[k] += alpha*(X[incX*i]*Y[incY*j] + X[incX*j]*Y[incY*i]);
	k++;
      }
    }
  }
  else {
    for(i=0; i<N; i++) {
      for(j=0; j<=i; j++) {
        A[k] += alpha*(X[incX*i]*Y[incY*j] + X[incX*j]*Y[incY*i]);
	k++;
      }
    }
  }
