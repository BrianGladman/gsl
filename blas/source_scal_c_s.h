
  size_t n;
  for(n=0; n<N; n++) {
    REAL(X,incX,n) *= alpha;
    IMAG(X,incX,n) *= alpha;
  }
