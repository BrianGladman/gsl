/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

  BASE_TYPE rr = 0.0;
  BASE_TYPE ri = 0.0;
  size_t n;
  for(n=0; n<N; n++) {
    rr += REAL(X, incX, n)*REAL(Y, incY, n) - CONJ_SIGN * IMAG(X, incX, n)*IMAG(Y, incY, n);
    ri += REAL(X, incX, n)*IMAG(Y, incY, n) + CONJ_SIGN * IMAG(X, incX, n)*REAL(Y, incY, n);
  }
  REAL0(result) = rr;
  IMAG0(result) = ri;
