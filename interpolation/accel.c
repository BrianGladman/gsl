/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include "gsl_interp.h"


gsl_interp_accel *
gsl_interp_accel_new(int heuristic, int cache_size)
{
  gsl_interp_accel * a = (gsl_interp_accel *) malloc(sizeof(gsl_interp_accel));
  if(a != 0) {
    a->cache_size = (cache_size > 0 ? cache_size : 1);
    a->cache_lo   = 0;
    a->cache_hi   = a->cache_lo + a->cache_size;
    a->heuristic  = heuristic;
    a->miss_count = 0;
    a->hit_count  = 0;
  }
  return a;
}


void
gsl_interp_accel_free(gsl_interp_accel * a)
{
  if(a != 0) free(a);
}
