/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_errno.h>
#include "bsearch.h"
#include "gsl_interp.h"

#define locMax(a, b)  ((a) > (b) ? (a) : (b))
#define locMin(a, b)  ((a) < (b) ? (a) : (b))


gsl_interp_linear *
gsl_interp_linear_new(const double x_array[], const double y_array[], int size)
{
  if(size <= 1)
    return 0;
  else {
    gsl_interp_linear * interp = (gsl_interp_linear *) malloc(sizeof(gsl_interp_linear));
    if(interp != 0) {
      interp->xmin = x_array[0];
      interp->xmax = x_array[size-1];
      interp->size = size;
    }
    return interp;
  }
}


void
gsl_interp_linear_free(gsl_interp_linear * interp)
{
  if(interp != 0) free(interp);
}


int
gsl_interp_linear_eval_direct_impl(const gsl_interp_linear * interp,
                                   const double x_array[], const double y_array[],
                                   double x,
                                   double * y
                                   )
{
  if(x > interp->xmax) {
    *y = y_array[interp->size - 1];
    return GSL_EDOM;
  }
  else if(x < interp->xmin) {
    *y = y_array[0];
    return GSL_EDOM;
  }
  else {
    int    index = bsearch(x_array, x, 0, interp->size - 1);
    double x_lo = x_array[index];
    double y_lo = y_array[index];
    double x_hi = x_array[index + 1];
    double y_hi = y_array[index + 1];
    double dx = x_hi - x_lo;
    if(dx > 0.0) {
      *y = y_lo + (x - x_lo)/dx * (y_hi - y_lo);
      return GSL_SUCCESS;
    }
    else {
      *y = 0.0;
      return GSL_FAILURE;
    }
  }
}


int
gsl_interp_linear_eval_impl(const gsl_interp_linear * interp,
                            gsl_interp_iter * iter,
                            const double x_array[], const double y_array[],
		            double x,
			    double * y
		            )
{
  if(x < interp->xmin) {
    *y = y_array[0];
    return GSL_EDOM;
  }
  else if(x > interp->xmax) {
    *y = y_array[interp->size - 1];
    return GSL_EDOM;
  }
  else {
    int     index_lo = iter->cache_lo;
    int     index_hi = iter->cache_hi;
    int     index;
    double  x_lo;
    double  x_hi;
    double  y_lo;
    double  y_hi;
    double  dx;


    /* find correct bin; check for iterator cache hits */

    if(x < x_array[index_lo]) {
      /* cache miss: lo side */

      if(iter->type == GSL_INTERP_LOCALSTEP) {
        /* reverse step heuristic */
	index_lo = locMax(0, index_lo - iter->cache_size);
	index_hi = locMin(interp->size - 1, index_lo + iter->cache_size);
	if(x >= x_array[index_lo] && x <= x_array[index_hi]) {
	  ++iter->hit_count;
	  index = CHECK_BSEARCH(x_array, x, index_lo, index_hi);
	}
	else {
	  ++iter->miss_count;
	  index = bsearch(x_array, x, 0, index_lo);
	}
      }
      else {
        ++iter->miss_count;
	index = bsearch(x_array, x, 0, index_lo);
      }
    }
    else if(x > x_array[index_hi]) {
      /* cache miss: hi side */

      if(iter->type == GSL_INTERP_LOCALSTEP) {
        /* forward step heuristic */
	index_hi = locMin(interp->size - 1, index_hi + iter->cache_size);
	index_lo = locMax(0, index_hi - iter->cache_size);
	if(x >= x_array[index_lo] && x <= x_array[index_hi]) {
	  ++iter->hit_count;
	  index = CHECK_BSEARCH(x_array, x, index_lo, index_hi);
	}
	else {
	  ++iter->miss_count;
	  index = bsearch(x_array, x, index_hi, interp->size - 1);
	}
      }
      else {
        ++iter->miss_count;
	index = bsearch(x_array, x, index_hi, interp->size - 1);
      }
    }
    else {
      /* cache hit */
      ++iter->hit_count;
      index = CHECK_BSEARCH(x_array, x, index_lo, index_hi);
    }

      
    /* adjust iterator cache */
    index_lo = locMax(0, index - iter->cache_size/2);
    index_hi = locMin(interp->size - 1, index_lo + iter->cache_size);

    /* evaluate */
    x_lo = x_array[index];
    x_hi = x_array[index + 1];
    y_lo = y_array[index];
    y_hi = y_array[index + 1];
    dx = x_hi - x_lo;
    if(dx > 0.0) {
      *y = y_lo + (x - x_lo)/dx * (y_hi - y_lo);
      return GSL_SUCCESS;
    }
    else {
      *y = 0.0;
      return GSL_FAILURE;
    }
  }
}
