/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

  BASE_TYPE scale = 0.0;
  BASE_TYPE ssq   = 1.0;
  size_t n;
  for(n=0; n<N; n++) {
    const BASE_TYPE xi = REAL(X,incX,n);
    const BASE_TYPE yi = IMAG(X,incX,n);
    if (xi != 0) {
      const BASE_TYPE axi = fabs(xi);

      if(scale < axi) {
        ssq   = 1.0 + ssq*(scale/axi)*(scale/axi);
        scale = axi;
      }
      else {
        ssq += (axi/scale)*(axi/scale);
      }
    }
    if (yi != 0) {
      const BASE_TYPE ayi = fabs(yi);

      if(scale < ayi) {
        ssq   = 1.0 + ssq*(scale/ayi)*(scale/ayi);
        scale = ayi;
      }
      else {
        ssq += (ayi/scale)*(ayi/scale);
      }
    }
  }
  return scale * sqrt(ssq);
