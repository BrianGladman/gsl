/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

  BASE_TYPE max = 0.0;
  size_t i = 0;
  CBLAS_INDEX n;
  CBLAS_INDEX result;
  for(n=0; n<N; n++) {
    if(fabs(X[i]) > max) {
      max = fabs(X[i]);
      result = n;
    }
    i += incX;
  }
  return result;
