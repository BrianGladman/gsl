/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

  size_t i, j, k;
  size_t kk = 0;
  const BASE_TYPE conj = -1.0;
  const BASE_TYPE aR = REAL0(alpha);
  const BASE_TYPE aI = IMAG0(alpha);
  const BASE_TYPE bR = REAL0(beta);
  const BASE_TYPE bI = IMAG0(beta);

  for(i=0; i<N; i++) {
    const BASE_TYPE yR = REAL(Y, incY, i);
    const BASE_TYPE yI = IMAG(Y, incY, i);
    REAL(Y, incY, i) = yR * bR - yI * bI;
    IMAG(Y, incY, i) = yR * bI + yI * bR;
  }

  if(Uplo == CblasUpper) {
    for(j=0; j<N; j++) {
      BASE_TYPE tmp1R = aR * REAL(X, incX, j) - aI * IMAG(X, incX, j);
      BASE_TYPE tmp1I = aR * IMAG(X, incX, j) + aI * REAL(X, incX, j);
      BASE_TYPE tmp2R = 0.0;
      BASE_TYPE tmp2I = 0.0;
      REAL(Y, incY, j) += tmp1R * REAL(Ap, 1, kk) - tmp1I * IMAG(Ap, 1, kk);
      IMAG(Y, incY, j) += tmp1R * IMAG(Ap, 1, kk) + tmp1I * REAL(Ap, 1, kk);
      i = j;
      for(k=kk+1; k<kk+N-j; k++) {
        const BASE_TYPE apkR = REAL(Ap, 1, k);
        const BASE_TYPE apkI = IMAG(Ap, 1, k);
        i++;
	REAL(Y, incY, i) += tmp1R * apkR - tmp1I * conj * apkI;
	IMAG(Y, incY, i) += tmp1I * apkR + tmp1R * conj * apkI;
	tmp2R += apkR * REAL(X, incX, i) - apkI * IMAG(X, incX, i);
	tmp2I += apkR * IMAG(X, incX, i) + apkI * REAL(X, incX, i);
      }
      REAL(Y, incY, j) += aR * tmp2R - aI * tmp2I;
      IMAG(Y, incY, j) += aR * tmp2I + aI * tmp2R;
      kk += N - j;
    }
  }
  else {
    for(j=0; j<N; j++) {
      BASE_TYPE tmp1R = aR * REAL(X, incX, j) - aI * IMAG(X, incX, j);
      BASE_TYPE tmp1I = aR * IMAG(X, incX, j) + aI * REAL(X, incX, j);
      BASE_TYPE tmp2R = 0.0;
      BASE_TYPE tmp2I = 0.0;
      i = 0;
      for(k=kk; k<kk+j; k++) {
        const BASE_TYPE apkR = REAL(Ap, 1, k);
	const BASE_TYPE apkI = IMAG(Ap, 1, k);
	REAL(Y, incY, i) += tmp1R * apkR - tmp1I * conj * apkI;
	IMAG(Y, incY, i) += tmp1I * apkR + tmp1R * conj * apkI;
	tmp2R += apkR * REAL(X, incX, i) - apkI * IMAG(X, incX, i);
	tmp2I += apkR * IMAG(X, incX, i) + apkI * REAL(X, incX, i);
        i++;
      }
      REAL(Y, incY, j) += tmp1R * REAL(Ap, 1, kk + j) - tmp1I * IMAG(Ap, 1, kk + j);
      IMAG(Y, incY, j) += tmp1R * IMAG(Ap, 1, kk + j) + tmp1I * REAL(Ap, 1, kk + j);
      REAL(Y, incY, j) += aR*tmp2R - aI*tmp2I;
      IMAG(Y, incY, j) += aI*tmp2R + aR*tmp2I;
      kk += j+1;
    }
  }
