/* falsepos.c -- false position root finding algorithm */
/* $Id$ */


/* config headers */
#include <config.h>

/* standard headers */
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* gsl headers */
#include <gsl_errno.h>
#include <gsl_roots.h>

/* roots headers */
#include "roots.h"

/* Note that sometimes this function checks against 2 * epsilon instead of
   epsilon. This is because if the interval is twice as wide as we want, we
   can return the midpoint of the interval and that will be within the desired
   accuracy. */
/* See the documentation for more information about this function. */
int
gsl_root_falsepos(double * root, double (* f)(double), double * lower_bound,
                  double * upper_bound, double rel_epsilon, double abs_epsilon,
                  unsigned int max_iterations, double max_deltay)
{
  /* Validate input. */
  {
    int status;
    
    /* Validate arguments. */
    if (_gsl_root_validate_bfp_args(root, f, lower_bound, upper_bound,
                                    rel_epsilon, abs_epsilon, max_iterations,
                                    max_deltay) == GSL_FAILURE)
      return GSL_FAILURE; /* GSL_ERROR was already called. */
    
    /* Make sure we have a root. */
    if (_gsl_root_ivt_guar(f, *lower_bound, *upper_bound) == GSL_FAILURE)
      return GSL_FAILURE; /* GSL_ERROR was already called. */

    /* Check if the user unwittingly provided the root. If _g_r_s_u returns
       non-zero, it has taken care of things and we should return what it
       returned. */
    if ((status = _gsl_root_silly_user(root, f, *lower_bound, *upper_bound,
                                       rel_epsilon, abs_epsilon, max_deltay)))
      return status;
  }

  /* Doh! It looks like we'll have to do actual work. */
  {
    double splitpoint, fl, fu, fs, old_lower_bound, old_upper_bound;
    unsigned int iterations;
  
    /* Evaluate the function under search at *lower_bound and *upper_bound. */
    _BARF_FPCALL(f, *lower_bound, fl);
    _BARF_FPCALL(f, *upper_bound, fu);
  
    for (iterations = 0; iterations < max_iterations; iterations++) {
      /* Draw a line between f(*lower_bound) and f(*upper_bound) and note
         where that crosses the X axis; that's where we will split the
         interval. */
      _BARF_ZERO(fl - fu);
      splitpoint = *upper_bound - (fu * (*lower_bound - *upper_bound)
                                   / (fl - fu));
      _BARF_FPCALL(f, splitpoint, fs);

      /* If the split point is the root exactly, we're done. */
      if (fs == 0.0) {
        *root = splitpoint;
        return GSL_SUCCESS;
      }
      
      /* Make note of the soon-to-be-old lower and upper bounds. */
      old_lower_bound = *lower_bound;
      old_upper_bound = *upper_bound;

      /* Split the interval and discard the part which doesn't contain the
         root. */
      if ((fl > 0.0 && fs < 0.0) || (fl < 0.0 && fs > 0.0)) {
        *upper_bound = splitpoint;
        fu = fs;
      }
      /* If the root is not between *lower_bound and the split point, it is
         guaranteed to be between the split point and *upper_bound, because
         there is definitely a root between *lower_bound and *upper_bound, and
         f(*lower_bound), f(*upper_bound), and f(splitpoint) are all known to
         be nonzero. */
      else {
        *lower_bound = splitpoint;
        fl = fs;
      }
      
      /* Now, let's check if we're finished. FIXME.1: this test needs help! */
      _BARF_TOLS(*lower_bound, *upper_bound, rel_epsilon, abs_epsilon);
      _BARF_DELTAY(fl, fu, max_deltay);
      /* If the interval has collapsed to within twice what we need, the root
         is its midpoint. */
      if (_WITHIN_TOL(*lower_bound, *upper_bound, 2 * rel_epsilon,
                      2 * abs_epsilon)) {
        *root = (*lower_bound + *upper_bound) / 2.0;
        return GSL_SUCCESS;
      }
      /* If the lower bound stayed the same and the upper bound moved less
         than epsilon, the root is *upper_bound. */
      if (old_lower_bound == *lower_bound
          && _WITHIN_TOL(old_upper_bound, *upper_bound, rel_epsilon,
                         abs_epsilon)) {
        *root = *upper_bound;
        return GSL_SUCCESS;
      }
      /* If the upper bound stayed the same and the lower bound moved less
         than epsilon, the root is *lower_bound. */
      if (old_upper_bound == *upper_bound
          && _WITHIN_TOL(old_lower_bound, *lower_bound, rel_epsilon,
                         abs_epsilon)) {
        *root = *lower_bound;
        return GSL_SUCCESS;
      }
    }
    
    /* Uh oh, ran out of iterations. */
    GSL_ERROR("exceeded maximum number of iterations", GSL_ETIMEOUT);
  }
}
