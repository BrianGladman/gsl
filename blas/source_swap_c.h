
  size_t n;
  for(n=0; n<N; n++) {
    BASE_TYPE tmpr = REAL(X, incX, n);
    BASE_TYPE tmpi = IMAG(X, incX, n);
    REAL(X, incX, n) = REAL(Y, incY, n);
    IMAG(X, incX, n) = IMAG(Y, incY, n);
    REAL(Y, incY, n) = tmpr;
    IMAG(Y, incY, n) = tmpi;
  }
