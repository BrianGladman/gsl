
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    REAL(Y, incY, j) = REAL(X, incX, i);
    IMAG(Y, incY, j) = IMAG(X, incX, i);
    i += incX;
    j += incY;
  }
