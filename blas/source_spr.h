
  size_t i, j, k;
  size_t ix, iy, jx, jy, kk;
  
  kk = 0;

  if(Uplo == CblasUpper) {
    jx = 0;
    for(j=0; j<N; j++) {
      BASE_TYPE tmp = alpha * X[jx];
      ix = 0;
      for(k=kk; k<kk+j; k++) {
        Ap[k] += X[ix] * tmp;
        ix += incX;
      }
      jx += incX;
      kk += j;
    }
  }
  else {
    jx = 0;
    for(j=0; j<N; j++) {
      BASE_TYPE tmp = alpha * X[jx];
      ix = jx;
      for(k=kk; k<=kk+N-j; k++) {
        Ap[k] += X[ix]*tmp;
      }
      jx += incX;
      kk += N - j + 1;
    }
  }
