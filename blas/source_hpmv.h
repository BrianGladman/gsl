
  size_t i, j, k;
  size_t kk = 0;
  const BASE_TYPE conj = -1.0;

  for(i=0; i<N; i++) {
    BASE_TYPE tmpR = REAL(Y, incY, i) * REAL0(beta) - IMAG(Y, incY, i) * IMAG0(beta);
    BASE_TYPE tmpI = REAL(Y, incY, i) * IMAG0(beta) + IMAG(Y, incY, i) * REAL0(beta);
    REAL(Y, incY, i) = tmpR;
    IMAG(Y, incY, i) = tmpI;
  }
  
  if(Uplo == CblasUpper) {
    for(j=0; j<N; j++) {
      BASE_TYPE tmp1R = REAL0(alpha) * REAL(X, incX, j) - IMAG0(alpha) * IMAG(X, incX, j);
      BASE_TYPE tmp1I = REAL0(alpha) * IMAG(X, incX, j) + IMAG0(alpha) * REAL(X, incX, j);
      BASE_TYPE tmp2R = 0.0;
      BASE_TYPE tmp2I = 0.0;
      i = 0;
      for(k=kk; k<kk+j-1; k++) {
        BASE_TYPE apkR = REAL(Ap, 1, k);
	BASE_TYPE apkI = IMAG(Ap, 1, k);
	REAL(Y, incY, i) += tmp1R * apkR - tmp1I * apkI;
	IMAG(Y, incY, i) += tmp1R * apkI + tmp1I * apkR;
	tmp2R += apkR * REAL(X, incX, i) - conj * apkI * IMAG(X, incX, i);
	tmp2I += apkR * IMAG(X, incX, i) + conj * apkI * REAL(X, incX, i);
        i++;
      }
      REAL(Y, incY, j) += tmp1R * REAL(Ap, 1, kk + j - 1) - tmp1I * IMAG(Ap, 1, kk + j - 1);
      IMAG(Y, incY, j) += tmp1R * IMAG(Ap, 1, kk + j - 1) + tmp1I * REAL(Ap, 1, kk + j - 1);
      kk += j;
    }
  }
  else {
    for(j=0; j<N; j++) {
      BASE_TYPE tmp1R = REAL0(alpha) * REAL(X, incX, j) - IMAG0(alpha) * IMAG(X, incX, j);
      BASE_TYPE tmp1I = REAL0(alpha) * IMAG(X, incX, j) + IMAG0(alpha) * REAL(X, incX, j);
      BASE_TYPE tmp2R = 0.0;
      BASE_TYPE tmp2I = 0.0;
      REAL(Y, incY, j) += tmp1R * REAL(Ap, 1, kk) - tmp1I * IMAG(Ap, 1, kk);
      IMAG(Y, incY, j) += tmp1R * IMAG(Ap, 1, kk) + tmp1I * REAL(Ap, 1, kk);
      i = j;
      for(k=kk+1; k<=kk+N-j; k++) {
        BASE_TYPE apkR = REAL(Ap, 1, k);
        BASE_TYPE apkI = IMAG(Ap, 1, k);
        i++;
	REAL(Y, incY, i) += tmp1R * apkR - tmp1I * apkI;
	IMAG(Y, incY, i) += tmp1R * apkI + tmp1I * apkR;
	tmp2R += apkR * REAL(X, incX, i) - apkI * IMAG(X, incX, i);
	tmp2I += apkR * IMAG(X, incX, i) + apkI * REAL(X, incX, i);
      }
      REAL(Y, incY, j) += REAL0(alpha) * tmp2R - IMAG0(alpha) * tmp2I;
      IMAG(Y, incY, j) += REAL0(alpha) * tmp2I + IMAG0(alpha) * tmp2R;
      kk += N - j + 1;
    }
  }
