
  size_t i, j;
  size_t ix, iy;
  size_t lenX, lenY;
  
  if(TransA == CblasNoTrans) {
    lenX = N;
    lenY = M;
  }
  else {
    lenX = M;
    lenY = N;
  }

  iy = 0;
  for(i=0; i<lenY; i++) {
    Y[iy] *= beta;
    iy += incY;
  }

  if(TransA == CblasNoTrans) {
    iy = 0;
    for(i=0; i<lenY; i++) {
      BASE_TYPE dot = 0.0;
      ix = 0;
      for(j=0; j<lenX; j++) {
        dot += X[ix]*A[lda*i + j];
        ix += incX;
      }
      Y[iy] += alpha * dot;
      iy += incY;
    }
  }
  else {
    ix = 0;
    for(j=0; j<lenX; j++) {
      BASE_TYPE tmp = alpha * X[ix];
      iy = 0;
      for(i=0; i<lenY; i++) {
        Y[iy] += A[lda*j + i] * tmp;
	iy += incY;
      }
      ix += incX;
    }
  }
