
  size_t i, j;
  const BASE_TYPE conj = -1.0;
  
  if(Uplo == CblasUpper) {
    for(j=0; j<N; j++) {
      const BASE_TYPE tmp1R = REAL0(alpha) * REAL(Y, incY, j) - conj * IMAG0(alpha) * IMAG(Y, incY, j);
      const BASE_TYPE tmp1I = REAL0(alpha) * IMAG(Y, incY, j) + conj * IMAG0(alpha) * REAL(Y, incY, j);
      const BASE_TYPE tmp2R = REAL0(alpha) * REAL(X, incX, j) - IMAG0(alpha) * IMAG(X, incX, j);
      const BASE_TYPE tmp2I = REAL0(alpha) * IMAG(X, incX, j) + IMAG0(alpha) * REAL(X, incX, j);
      for(i=j; i<N; i++) {
        REAL(A, 1, lda*j + i) += REAL(X, incX, i) * tmp1R - conj * IMAG(X, incX, i) * tmp1I;
	IMAG(A, 1, lda*j + i) += REAL(X, incX, i) * tmp1I + conj * IMAG(X, incX, i) * tmp1R;
        REAL(A, 1, lda*j + i) += REAL(Y, incY, i) * tmp2R - conj * IMAG(Y, incY, i) * tmp2I;
	IMAG(A, 1, lda*j + i) += REAL(Y, incY, i) * tmp2I + conj * IMAG(Y, incY, i) * tmp2R;
      }
    }
  }
  else {
    for(j=0; j<N; j++) {
      const BASE_TYPE tmp1R = REAL0(alpha) * REAL(Y, incY, j) - conj * IMAG0(alpha) * IMAG(Y, incY, j);
      const BASE_TYPE tmp1I = REAL0(alpha) * IMAG(Y, incY, j) + conj * IMAG0(alpha) * REAL(Y, incY, j);
      const BASE_TYPE tmp2R = REAL0(alpha) * REAL(X, incX, j) - IMAG0(alpha) * IMAG(X, incX, j);
      const BASE_TYPE tmp2I = REAL0(alpha) * IMAG(X, incX, j) + IMAG0(alpha) * REAL(X, incX, j);
      for(i=0; i<=j; i++) {
        REAL(A, 1, lda*j + i) += REAL(X, incX, i) * tmp1R - conj * IMAG(X, incX, i) * tmp1I;
	IMAG(A, 1, lda*j + i) += REAL(X, incX, i) * tmp1I + conj * IMAG(X, incX, i) * tmp1R;
        REAL(A, 1, lda*j + i) += REAL(Y, incY, i) * tmp2R - conj * IMAG(Y, incY, i) * tmp2I;
	IMAG(A, 1, lda*j + i) += REAL(Y, incY, i) * tmp2I + conj * IMAG(Y, incY, i) * tmp2R;
      }
    }
  }
