
  BASE_TYPE rr = 0.0;
  BASE_TYPE ri = 0.0;
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    rr += REAL(X, incX, i)*REAL(Y, incY, j) - CONJ_SIGN * IMAG(X, incX, i)*IMAG(Y, incY, j);
    ri += REAL(X, incX, i)*IMAG(Y, incY, j) + CONJ_SIGN * IMAG(X, incX, i)*REAL(Y, incY, j);
    i += incX;
    j += incY;
  }
  REAL0(result) = rr;
  IMAG0(result) = ri;
