
  size_t i, j;
  size_t jy = 0;

  for(j=0; j<M; j++) {
    const BASE_TYPE tmp = alpha * Y[jy];
    size_t ix = 0;
    for(i=0; i<N; i++) {
      A[lda*j + i] += X[ix]*tmp;
      ix += incX;
    }
    jy += incY;
  }
