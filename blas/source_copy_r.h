/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

  size_t n;
  size_t i = 0;
  size_t j = 0;
  for(n=0; n<N; n++) {
    Y[j] = X[i];
    i += incX;
    j += incY;
  }
