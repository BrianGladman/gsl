
  BASE_TYPE r = 0.0;
  size_t n;
  size_t i = 0;
  for(n=0; n<N; n++) {
    r += fabs(X[i]);
    i += incX;
  }
  return r;
