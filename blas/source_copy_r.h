
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    Y[j] = X[i];
    i += incX;
    j += incY;
  }
