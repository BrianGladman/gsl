/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

  BASE_TYPE max = 0.0;
  CBLAS_INDEX n;
  CBLAS_INDEX result;
  for(n=0; n<N; n++) {
    const BASE_TYPE a = fabs(REAL(X, incX, n)) + fabs(IMAG(X, incX, n));
    if(a > max) {
      max = a;
      result = n;
    }
  }
  return result;
