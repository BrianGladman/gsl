/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include "gsl_interp.h"


gsl_interp_iter *
gsl_interp_iter_new(int type, int cache_size)
{
  gsl_interp_iter * it = (gsl_interp_iter *) malloc(sizeof(gsl_interp_iter));
  if(it != 0) {
    it->cache_size = (cache_size > 0 ? cache_size : 1);
    it->cache_lo = 0;
    it->cache_hi = it->cache_lo + it->cache_size;
    it->type = type;
    it->miss_count = 0;
    it->hit_count = 0;
  }
  return it;
}


void
gsl_interp_iter_free(gsl_interp_iter * iter)
{
  if(iter != 0) free(iter);
}
