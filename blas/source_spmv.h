
  size_t i, j, k;
  size_t ix, iy, jx, jy;
  size_t kk = 0;

  iy=0;
  for(i=0; i<N; i++) {
    Y[iy] *= beta;
    iy += incY;
  }
  
  if(Uplo == CblasUpper) {
    jx = 0;
    jy = 0;
    for(j=0; j<N; j++) {
      BASE_TYPE tmp1 = alpha*X[jx];
      BASE_TYPE tmp2 = 0.0;
      ix = 0;
      iy = 0;
      for(k=kk; k<kk+j-1; k++) {
        BASE_TYPE apk = Ap[k];
        Y[iy] += tmp1 * apk;
        tmp2  += apk * X[ix];
        ix += incX;
        iy += incY;
      }
      Y[jy] += tmp1*Ap[kk + j - 1] + alpha*tmp2;
      jx += incX;
      jy += incY;
      kk += j;
    }
  }
  else {
    jx = 0;
    jy = 0;
    for(j=0; j<N; j++) {
      BASE_TYPE tmp1 = alpha*X[jx];
      BASE_TYPE tmp2 = 0.0;
      Y[jy] += tmp1*Ap[kk];
      ix = jx;
      iy = jy;
      for(k=kk+1; k<=kk+N-j; k++) {
        BASE_TYPE apk = Ap[k];
        ix += incX;
	iy += incY;
	Y[iy] += tmp1 * apk;
	tmp2  += apk * X[ix];
      }
      Y[jy] += alpha*tmp2;
      jx += incX;
      jy += incY;
      kk += N - j + 1;
    }
  }
