
  size_t i;
  size_t ix = 0;
  size_t iy = 0;
  for(i=0; i<N; i++) {
    BASE_TYPE x = X[ix];
    BASE_TYPE y = Y[iy];
    X[ix] =  c*x + s*y;
    Y[iy] = -s*x + c*y;
    ix += incX;
    iy += incY;
  }
