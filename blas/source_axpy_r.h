
  size_t i;

  if(incX == 1 && incY == 1) {
    size_t m = N % 4;
    for(i=0; i<m; i++) {
      Y[i] += alpha * X[i];
    }
    for(i=m; i<N-3; i += 4) {
      Y[i]   += alpha*X[i];
      Y[i+1] += alpha*X[i+1];
      Y[i+2] += alpha*X[i+2];
      Y[i+3] += alpha*X[i+3];
    }
  }
  else {
    size_t ix = 0;
    size_t iy = 0;
    
    for(i=0; i<N; i++) {
      Y[iy] += alpha*X[ix];
      ix += incX;
      iy += incY;
    }
  }
