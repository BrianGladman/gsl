
  size_t i, j;
  const BASE_TYPE conj = -1.0;

  if(Uplo == CblasUpper) {
    for(j=0; j<N; j++) {
      BASE_TYPE tmpR = alpha * REAL(X, incX, j);
      BASE_TYPE tmpI = alpha * IMAG(X, incX, j);
      for(i=j; i<N; i++) {
        REAL(A, 1, lda*j + i) += REAL(X, incX, i) * tmpR - conj * IMAG(X, incX, i) * tmpI;
	IMAG(A, 1, lda*j + i) += REAL(X, incX, i) * tmpI + conj * IMAG(X, incX, i) * tmpR;
      }
    }
  }
  else {
    for(j=0; j<N; j++) {
      BASE_TYPE tmpR = alpha * REAL(X, incX, j);
      BASE_TYPE tmpI = alpha * IMAG(X, incX, j);
      for(i=0; i<=j; i++) {
        REAL(A, 1, lda*j + i) += REAL(X, incX, i) * tmpR - conj * IMAG(X, incX, i) * tmpI;
	IMAG(A, 1, lda*j + i) += REAL(X, incX, i) * tmpI + conj * IMAG(X, incX, i) * tmpR;
      }
    }
  }
