/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

  if (fabs(REAL0(alpha)) + fabs(IMAG0(alpha)) != 0.0 ) {
    size_t i;
    for(i=0; i<N; i++) {
      const BASE_TYPE tmpR = REAL0(alpha) * REAL(X, incX, i) - IMAG0(alpha) * IMAG(X, incX, i);
      const BASE_TYPE tmpI = REAL0(alpha) * IMAG(X, incX, i) + IMAG0(alpha) * REAL(X, incX, i);
      REAL(Y, incY, i) += tmpR;
      IMAG(Y, incY, i) += tmpI;
    }
  }
