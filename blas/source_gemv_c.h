
  size_t i, j;
  size_t lenX, lenY;

  if(TransA == CblasNoTrans) {
    lenX = N;
    lenY = M;
  }
  else {
    lenX = M;
    lenY = N;
  }

  for(i=0; i<lenY; i++) {
    BASE_TYPE tmpR = REAL(Y, incY, i) * REAL0(beta) - IMAG(Y, incY, i) * IMAG0(beta);
    BASE_TYPE tmpI = REAL(Y, incY, i) * IMAG0(beta) - IMAG(Y, incY, i) * REAL0(beta);
    REAL(Y, incY, i) = tmpR;
    IMAG(Y, incY, i) = tmpI;
  }

  if(TransA == CblasNoTrans) {
    for(i=0; i<lenY; i++) {
      BASE_TYPE dotR = 0.0;
      BASE_TYPE dotI = 0.0;
      for(j=0; j<lenX; j++) {
        dotR += REAL(X, incX, j) * REAL(A, 1, lda*i + j) - IMAG(X, incX, j) * IMAG(A, 1, lda*i + j);
	dotI += REAL(X, incX, j) * IMAG(A, 1, lda*i + j) + IMAG(X, incX, j) * REAL(A, 1, lda*i + j);
      }
      REAL(Y, incY, i) += REAL0(alpha) * dotR - IMAG0(alpha) * dotI;
      IMAG(Y, incY, i) += REAL0(alpha) * dotI + IMAG0(alpha) * dotR;
    }
  }
  else {
    for(j=0; j<lenX; j++) {
      BASE_TYPE tmpR = REAL0(alpha) * REAL(X, incX, j) - IMAG0(alpha) * IMAG(X, incX, j);
      BASE_TYPE tmpI = REAL0(alpha) * IMAG(X, incX, j) + IMAG0(alpha) * REAL(X, incX, j);
      for(i=0; i<lenY; i++) {
        REAL(Y, incY, i) += REAL(A, 1, lda*j + i) * tmpR - IMAG(A, 1, lda*j + i) * tmpI;
	IMAG(Y, incY, i) += REAL(A, 1, lda*j + i) * tmpI + IMAG(A, 1, lda*j + i) * tmpR;
      }
    }
  }
