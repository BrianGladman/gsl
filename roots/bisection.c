/* bisection.c -- bisection root finding algorithm */
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

/* Note that this function checks against 2 * epsilon, not epsilon. This is
   because if the interval is twice as wide as we want, we can return the
   midpoint of the interval and that will be within the desired accuracy. */
/* See the documentation for more information about this function. */
int
gsl_root_bisection(double * root, double (* f)(double), double * lower_bound,
                   double * upper_bound, double rel_epsilon,
                   double abs_epsilon, unsigned int max_iterations,
                   double max_deltay)
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
    double midpoint, fl, fu, fm;
    unsigned int iterations;
  
    /* Evaluate the function under search at *lower_bound and *upper_bound. */
    _BARF_FPCALL(f, *lower_bound, fl);
    _BARF_FPCALL(f, *upper_bound, fu);
  
    for (iterations = 0; iterations < max_iterations; iterations++) {
      /* Find the midpoint of the interval; this is where we will split it. */
      midpoint = (*upper_bound + *lower_bound) / 2.0;
      _BARF_FPCALL(f, midpoint, fm);
      
      /* If the midpoint is the root exactly, we're done. */
      if (fm == 0.0) {
        *root = midpoint;
        return GSL_SUCCESS;
      }
      
      /* Discard the half of the interval which doesn't contain the root. */
      if ((fl > 0.0 && fm < 0.0) || (fl < 0.0 && fm > 0.0)) {
        *upper_bound = midpoint;
        fu = fm;
      }
      /* If the root is not between *lower_bound and midpoint, it is guaranteed
         to be between midpoint and *upper_bound. There is definitely a root
         between *lower_bound and *upper_bound, and f(*lower_bound),
         f(*upper_bound), and f(midpoint) are all non-zero. */
      else {
        *lower_bound = midpoint;
        fl = fm;
      }
      
      /* Now, let's check if we're finished. */
      _BARF_TOLS(*lower_bound, *upper_bound, 2 * rel_epsilon, 2 * abs_epsilon);
      _BARF_DELTAY(fl, fu, max_deltay);
      if (_WITHIN_TOL(*lower_bound, *upper_bound, 2 * rel_epsilon,
                      2 * abs_epsilon)) {
        *root = (*lower_bound + *upper_bound) / 2.0;
        return GSL_SUCCESS;
      }
    }

    /* Uh oh, ran out of iterations. */
    GSL_ERROR("exceeded maximum number of iterations", GSL_ETIMEOUT);
  }
}
