
  size_t i, j;
  size_t iy;
  size_t lenX, lenY;

  if(alpha == 0.0 && beta == 1.0) return;  

  if(TransA == CblasNoTrans) {
    lenX = N;
    lenY = M;
  }
  else {
    lenX = M;
    lenY = N;
  }

  /* form  y := beta*y */
  if(beta != 1.0) {
    iy = 0;
    for(i=0; i<lenY; i++) {
      Y[iy] *= beta;
      iy += incY;
    }
  }

  if(alpha == 0.0) return;

  if(TransA == CblasNoTrans) {
    /* form  y := alpha*A*x + y */
    iy = 0;
    for(i=0; i<lenY; i++) {
      BASE_TYPE temp = 0.0;
      for(j=GSL_MAX(0,i-KL); j<GSL_MIN(lenX, i+KU+1); j++) {
        temp += X[incX * j] * A[lda*i + j];
      }
      Y[iy] += alpha * temp;
      iy += incY;
    }
  }
  else {
    /* form  y := alpha*A'*x + y */
    for(j=0; j<lenX; j++) {
      BASE_TYPE temp = alpha * X[incX * j];
      for(i = GSL_MAX(0, j-KU); i<GSL_MIN(lenY, j+KL+1); i++) {
        Y[incY * i] += temp * A[lda*j + i];
      }
    }
  }
