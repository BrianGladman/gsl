/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

  size_t i, j;
  size_t ix, iy;
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
      const size_t j_min = ( KL > i ? 0 : i-KL );
      ix = j_min * incX;
      for(j=j_min; j<GSL_MIN(lenX, i+KU+1); j++) {
        temp += X[ix] * A[lda*i + j];
	ix += incX;
      }
      Y[iy] += alpha * temp;
      iy += incY;
    }
  }
  else {
    /* form  y := alpha*A'*x + y */
    ix = 0;
    for(j=0; j<lenX; j++) {
      const BASE_TYPE temp = alpha * X[ix];
      if(temp != 0.0) {
        const size_t i_min = ( KU > j ? 0 : j-KU );
	iy = i_min * incY;
        for(i=i_min; i<GSL_MIN(lenY, j+KL+1); i++) {
          Y[iy] += temp * A[lda*j + i];
	  iy += incY;
        }
      }
      ix += incX;
    }
  }
