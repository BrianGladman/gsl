/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdlib.h>
#include "bsearch.h"
#include "gsl_interp.h"


gsl_interp_accel *
gsl_interp_accel_new(void)
{
  gsl_interp_accel * a = (gsl_interp_accel *) malloc(sizeof(gsl_interp_accel));
  if(a != 0) {
    a->cache = 0;
    a->hit_count = 0;
    a->miss_count = 0;
  }
  return a;
}


size_t
gsl_interp_accel_find(gsl_interp_accel * a, const double xa[], size_t len, double x)
{
  size_t x_index = a->cache;
 
  if(x < xa[x_index]) {
    a->miss_count++;
    a->cache = interp_bsearch(xa, x, 0, x_index);
  }
  else if(x > xa[x_index + 1]) {
    a->miss_count++;
    a->cache = interp_bsearch(xa, x, x_index, len-1);
  }
  else {
    a->hit_count++;
  }
  
  return a->cache;
}


void
gsl_interp_accel_free(gsl_interp_accel * a)
{
  if(a != 0) {
    free(a);
  }
}
