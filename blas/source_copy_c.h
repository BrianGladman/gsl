/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

  size_t n;
  for(n=0; n<N; n++) {
    REAL(Y, incY, n) = REAL(X, incX, n);
    IMAG(Y, incY, n) = IMAG(X, incX, n);
  }
