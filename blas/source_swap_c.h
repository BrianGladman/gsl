
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    BASE_TYPE tmpr = REAL(X, incX, i);
    BASE_TYPE tmpi = IMAG(X, incX, i);
    REAL(X, incX, i) = REAL(Y, incY, j);
    IMAG(X, incX, i) = IMAG(Y, incY, j);
    REAL(Y, incY, j) = tmpr;
    IMAG(Y, incY, j) = tmpi;
    i += incX;
    j += incY;
  }
