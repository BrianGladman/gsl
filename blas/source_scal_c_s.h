
  size_t n;
  size_t ix = 0;
  for(n=0; n<N; n++) {
    REAL(X,incX,ix) *= alpha;
    IMAG(X,incX,ix) *= alpha;
    ix += incX;
  }
