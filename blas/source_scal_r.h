/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

  size_t n;
  size_t ix = 0;
  for(n=0; n<N; n++) {
    X[ix] *= alpha;
    ix += incX;
  }
