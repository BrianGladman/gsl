
  size_t i, j, k;
  size_t kk = 0;

  if(Uplo == CblasUpper) {
    for(j=0; j<N; j++) {
      BASE_TYPE tmpR = alpha * REAL(X, incX, j);
      BASE_TYPE tmpI = alpha * IMAG(X, incX, j);
      i = 0;
      for(k=kk; k<kk+j; k++) {
        REAL(Ap, 1, k) += REAL(X, incX, i) * tmpR - IMAG(X, incX, i) * tmpI;
	IMAG(Ap, 1, k) += REAL(X, incX, i) * tmpI + IMAG(X, incX, i) * tmpR;
        i++;
      }
      kk += j;
    }
  }
  else {
    for(j=0; j<N; j++) {
      BASE_TYPE tmpR = alpha * REAL(X, incX, j);
      BASE_TYPE tmpI = alpha * IMAG(X, incX, j);
      i = j;
      for(k=kk; k<=kk+N-j; k++) {
        REAL(Ap, 1, k) += REAL(X, incX, i) * tmpR - IMAG(X, incX, i) * tmpI;
	IMAG(Ap, 1, k) += REAL(X, incX, i) * tmpI + IMAG(X, incX, i) * tmpR;
      }
      kk += N - j + 1;
    }
  }
