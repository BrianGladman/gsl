/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

  BASE_TYPE r = 0.0;
  size_t n;
  for(n=0; n<N; n++) {
    r += fabs(REAL(X, incX, n)) + fabs(IMAG(X, incX, n));
  }
  return r;
