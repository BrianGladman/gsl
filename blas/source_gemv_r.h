#define KU GSL_MAX(M,N)
#define KL GSL_MAX(M,N)
#include "source_gbmv_r.h"
#undef KU
#undef KL

#if 0
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

  if(beta != 1.0) {
    iy = 0;
    for(i=0; i<lenY; i++) {
      Y[iy] *= beta;
      iy += incY;
    }
  }

  if(alpha == 0.0) return;

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
#endif
