/* bisection.c -- bisection root finding algorithm */

#include <config.h>

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <gsl_errno.h>
#include <gsl_roots.h>

#include "roots.h"

/* Note that this function checks against 2 * epsilon, not epsilon.
   This is because if the interval is twice as wide as we want, we can
   return the midpoint of the interval and that will be within the
   desired accuracy. */

int
gsl_root_bisection (double *root, double (*f) (double), 
		    double *lower_bound, double *upper_bound, 
		    double rel_epsilon, double abs_epsilon, 
		    unsigned int max_iterations)
{
  unsigned int iterations;
  double midpoint, fl, fu, fm;

  if (*lower_bound >= *upper_bound)
    GSL_ERROR ("lower bound larger than upper_bound", GSL_EINVAL);
 
  if (rel_epsilon < 0.0 || abs_epsilon < 0.0)
    GSL_ERROR ("relative or absolute tolerance negative", GSL_EBADTOL);

 if (rel_epsilon < DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER)
    GSL_ERROR ("relative tolerance too small", GSL_EBADTOL);

  SAFE_FUNC_CALL (f, *lower_bound, fl);

  if (fl == 0.0)
    {
      *root = *lower_bound;
      return GSL_SUCCESS;
    }

  SAFE_FUNC_CALL (f, *upper_bound, fu);

  if (fu == 0.0)
    {
      *root = *upper_bound;
      return GSL_SUCCESS;
    }

  if ((fl < 0 && fu < 0.0) || (fl > 0 && fu > 0))
    {
      GSL_ERROR ("endpoints do not straddle y=0", GSL_EINVAL);
    }

  if (WITHIN_TOL (*lower_bound, *upper_bound, 2 * rel_epsilon, 2 * abs_epsilon))
    {
      *root = (*lower_bound + *upper_bound) / 2.0;
      return GSL_SUCCESS;
    }
  
  for (iterations = 0; iterations < max_iterations; iterations++)
    {
      midpoint = (*upper_bound + *lower_bound) / 2.0;

      SAFE_FUNC_CALL (f, midpoint, fm);
      
      if (fm == 0.0)
	{
	  *root = midpoint;
	  return GSL_SUCCESS;
	}
      
      /* Discard the half of the interval which doesn't contain the root. */

      if ((fl > 0.0 && fm < 0.0) || (fl < 0.0 && fm > 0.0))
	{
	  *upper_bound = midpoint;
	  fu = fm;
	}
      else
	{
	  *lower_bound = midpoint;
	  fl = fm;
	}
      

      if (WITHIN_TOL (*lower_bound, *upper_bound, 2 * rel_epsilon, 2 * abs_epsilon))
	{
	  *root = (*lower_bound + *upper_bound) / 2.0;
	  return GSL_SUCCESS;
	}
      else
	{
	  CHECK_TOL (*lower_bound, *upper_bound, 2 * rel_epsilon, 2 * abs_epsilon);
	}
    }
  
  GSL_ERROR ("exceeded maximum number of iterations", GSL_EMAXITER);
  
}
