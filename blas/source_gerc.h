
  size_t i, j;
  const BASE_TYPE conj = -1.0;

  for(j=0; j<M; j++) {
    const BASE_TYPE tmpR = REAL0(alpha) * REAL(X, incX, j) - IMAG0(alpha) * IMAG(X, incX, j);
    const BASE_TYPE tmpI = REAL0(alpha) * IMAG(X, incX, j) + IMAG0(alpha) * REAL(X, incX, j);
    for(i=0; i<N; i++) {
      REAL(A, 1, lda*j + i) += REAL(Y, incY, i) * tmpR - conj * IMAG(Y, incY, i) * tmpI;
      IMAG(A, 1, lda*j + i) += REAL(Y, incY, i) * tmpI + conj * IMAG(Y, incY, i) * tmpR;
    }
  }
