
  size_t i, j;
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
      REAL(Y, incY, j) += tmp1R * REAL(A, 1, lda*j + j) - tmp1I * IMAG(A, 1, lda*j + j);
      REAL(Y, incY, j) += tmp1R * IMAG(A, 1, lda*j + j) + tmp1I * REAL(A, 1, lda*j + j);
      for(i=j+1; i<GSL_MIN(N,j+K+1); i++) {
	REAL(Y, incY, i) += tmp1R * REAL(A, 1, lda*j + i) - tmp1I * IMAG(A, 1, lda*j + i);
	IMAG(Y, incY, i) += tmp1R * IMAG(A, 1, lda*j + i) + tmp1I * REAL(A, 1, lda*j + i);
	tmp2R += REAL(A, 1, lda*j + i) * REAL(X, incX, i) - conj * IMAG(A, 1, lda*j + i) * IMAG(X, incX, i);
	tmp2R += REAL(A, 1, lda*j + i) * IMAG(X, incX, i) + conj * IMAG(A, 1, lda*j + i) * REAL(X, incX, i);
      }
      REAL(Y, incY, j) += REAL0(alpha) * tmp2R - IMAG0(alpha) * tmp2I;
      IMAG(Y, incY, j) += REAL0(alpha) * tmp2I + IMAG0(alpha) * tmp2R;
    }
  }
  else {
    for(j=0; j<N; j++) {
      BASE_TYPE tmp1R = REAL0(alpha) * REAL(X, incX, j) - IMAG0(alpha) * IMAG(X, incX, j);
      BASE_TYPE tmp1I = REAL0(alpha) * IMAG(X, incX, j) + IMAG0(alpha) * REAL(X, incX, j);
      BASE_TYPE tmp2R = 0.0;
      BASE_TYPE tmp2I = 0.0;
      for(i=GSL_MAX(0,j-K); i<j; i++) {
        REAL(Y, incY, i) += tmp1R * REAL(A, 1, lda*j + i) - tmp1I * IMAG(A, 1, lda*j + i);
	IMAG(Y, incY, i) += tmp1R * IMAG(A, 1, lda*j + i) + tmp1I * REAL(A, 1, lda*j + i);
	tmp2R += REAL(A, 1, lda*j + i) * REAL(X, incX, i) - conj * IMAG(A, 1, lda*j + i) * IMAG(X, incX, i);
	tmp2I += REAL(A, 1, lda*j + i) * IMAG(X, incX, i) + conj * IMAG(A, 1, lda*j + i) * REAL(X, incX, i);
      }
      REAL(Y, incY, j) += tmp1R * REAL(A, 1, lda*j + j) - tmp1I * IMAG(A, 1, lda*j + j);
      REAL(Y, incY, j) += REAL0(alpha) * tmp2R - IMAG0(alpha) * tmp2I;
      IMAG(Y, incY, j) += tmp1R * IMAG(A, 1, lda*j + j) + tmp1I * REAL(A, 1, lda*j + j);
      IMAG(Y, incY, j) += REAL0(alpha) * tmp2I + IMAG0(alpha) * tmp2R;
    }
  }
