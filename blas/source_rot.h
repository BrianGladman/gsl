/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

  size_t i;
  size_t ix = 0;
  size_t iy = 0;
  for(i=0; i<N; i++) {
    const BASE_TYPE x = X[ix];
    const BASE_TYPE y = Y[iy];
    X[ix] =  c*x + s*y;
    Y[iy] = -s*x + c*y;
    ix += incX;
    iy += incY;
  }
