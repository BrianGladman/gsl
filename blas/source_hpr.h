/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

  size_t i, j, k;
  size_t kk = 0;
  const BASE_TYPE conj = -1.0;

  if(Uplo == CblasUpper) {
    for(j=0; j<N; j++) {
      const BASE_TYPE tmpR = alpha * REAL(X, incX, j);
      const BASE_TYPE tmpI = alpha * IMAG(X, incX, j);
      i = j;
      for(k=kk; k<kk+N-j; k++) {
        REAL(Ap, 1, k) += REAL(X, incX, i) * tmpR - conj * IMAG(X, incX, i) * tmpI;
	IMAG(Ap, 1, k) += REAL(X, incX, i) * tmpI + conj * IMAG(X, incX, i) * tmpR;
	i++;
      }
      kk += N - j;
    }
  }
  else {
    for(j=0; j<N; j++) {
      const BASE_TYPE tmpR = alpha * REAL(X, incX, j);
      const BASE_TYPE tmpI = alpha * IMAG(X, incX, j);
      i = 0;
      for(k=kk; k<kk+j+1; k++) {
        REAL(Ap, 1, k) += REAL(X, incX, i) * tmpR - conj * IMAG(X, incX, i) * tmpI;
	IMAG(Ap, 1, k) += REAL(X, incX, i) * tmpI + conj * IMAG(X, incX, i) * tmpR;
        i++;
      }
      kk += j+1;
    }
  }
