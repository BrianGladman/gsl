
  BASE_TYPE r = 0.0;
  size_t n;
  size_t i;
  for(n=0; n<N; n++) {
    r += fabs(REAL(X, incX, i)) + fabs(IMAG(X, incX, i));
    i += incX;
  }
  return r;
