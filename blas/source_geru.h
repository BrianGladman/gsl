
  size_t i, j;

  for(j=0; j<M; j++) {
    BASE_TYPE tmpR = REAL0(alpha) * REAL(Y, incY, j) - IMAG0(alpha) * IMAG(Y, incY, j);
    BASE_TYPE tmpI = REAL0(alpha) * IMAG(Y, incY, j) + IMAG0(alpha) * REAL(Y, incY, j);
    for(i=0; i<N; i++) {
      REAL(A, 1, lda*j + i) += REAL(X, incX, i) * tmpR - IMAG(X, incX, i) * tmpI;
      REAL(A, 1, lda*j + i) += REAL(X, incX, i) * tmpI + IMAG(X, incX, i) * tmpR;
    }
  }
