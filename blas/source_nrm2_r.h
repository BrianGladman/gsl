/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

  BASE_TYPE scale = 0.0;
  BASE_TYPE ssq   = 1.0;
  size_t n;
  size_t i = 0;
  for(n=0; n<N; n++) {
    BASE_TYPE axi = fabs(X[i]);
    if(scale < axi) {
      ssq   = 1.0 + ssq*(scale/axi)*(scale/axi);
      scale = axi;
    }
    else {
      ssq += (axi/scale)*(axi/scale);
    }
    i += incX;
  }
  return scale * sqrt(ssq);
