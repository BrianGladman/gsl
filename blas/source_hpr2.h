
  size_t i, j;

  if(Uplo == CblasUpper) {
    size_t k = 0;
    for(i=0; i<N; i++) {
      for(j=i; j<N; j++) {
        BASE_TYPE tmp1R = REAL(X, incX, i) * REAL(Y, incY, j) - IMAG(X, incX, i) * IMAG(Y, incY, j);
	BASE_TYPE tmp1I = REAL(X, incX, i) * IMAG(Y, incY, j) + IMAG(X, incX, i) * REAL(Y, incY, j);
	BASE_TYPE tmp2R = REAL(X, incX, j) * REAL(Y, incY, i) - IMAG(X, incX, j) * IMAG(Y, incY, i);
	BASE_TYPE tmp2I = REAL(X, incX, j) * IMAG(Y, incY, i) + IMAG(X, incX, j) * REAL(Y, incY, i);
	BASE_TYPE tmpR  = REAL0(alpha) * tmp1R - IMAG0(alpha) * tmp1I
	                + REAL0(alpha) * tmp2R + IMAG0(alpha) * tmp2I;
	BASE_TYPE tmpI  = REAL0(alpha) * tmp1I + IMAG0(alpha) * tmp1R
	                + REAL0(alpha) * tmp2I - IMAG0(alpha) * tmp2R;
        REAL(Ap, 1, k) += tmpR;
	IMAG(Ap, 1, k) += tmpI;
	k++;
      }
    }
  }
  else {
    size_t k = 0;
    for(i=0; i<N; i++) {
      for(j=0; j<=i; j++) {
        BASE_TYPE tmp1R = REAL(X, incX, i) * REAL(Y, incY, j) - IMAG(X, incX, i) * IMAG(Y, incY, j);
	BASE_TYPE tmp1I = REAL(X, incX, i) * IMAG(Y, incY, j) + IMAG(X, incX, i) * REAL(Y, incY, j);
	BASE_TYPE tmp2R = REAL(X, incX, j) * REAL(Y, incY, i) - IMAG(X, incX, j) * IMAG(Y, incY, i);
	BASE_TYPE tmp2I = REAL(X, incX, j) * IMAG(Y, incY, i) + IMAG(X, incX, j) * REAL(Y, incY, i);
	BASE_TYPE tmpR  = REAL0(alpha) * tmp1R - IMAG0(alpha) * tmp1I
	                + REAL0(alpha) * tmp2R + IMAG0(alpha) * tmp2I;
	BASE_TYPE tmpI  = REAL0(alpha) * tmp1I + IMAG0(alpha) * tmp1R
	                + REAL0(alpha) * tmp2I - IMAG0(alpha) * tmp2R;
        REAL(Ap, 1, k) += tmpR;
	IMAG(Ap, 1, k) -= tmpI; /* -= for lower triangle due to conjugate */
	k++;
      }
    }
  }
