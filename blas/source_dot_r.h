/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

  ACC_TYPE r = INIT_VAL;
  size_t n;
  size_t i = 0;
  size_t j = 0;
  for(n=0; n<N; n++) {
    r += X[i]*Y[j];
    i += incX;
    j += incY;
  }
  return r;
