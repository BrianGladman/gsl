
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    BASE_TYPE tmp = X[i];
    X[i] = Y[j];
    Y[j] = tmp;
    i += incX;
    j += incY;
  }
