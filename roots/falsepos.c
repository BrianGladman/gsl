/* falsepos.c -- false position root finding algorithm */

#include <config.h>

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <gsl_math.h>
#include <gsl_errno.h>
#include <gsl_roots.h>

#include "roots.h"

/* Note that sometimes this function checks against 2 * epsilon instead of
   epsilon. This is because if the interval is twice as wide as we want, we
   can return the midpoint of the interval and that will be within the desired
   accuracy. */

int
gsl_root_falsepos (double *root, double (*f) (double), 
		   double *lower_bound, double *upper_bound, 
		   double rel_epsilon, double abs_epsilon,
		   unsigned int max_iterations)
{
  unsigned int iterations;
  enum { UPPER, LOWER } moved = 0;
  double splitpoint, fl, fu, fs, old_lower_bound, old_upper_bound;
  
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
      *upper_bound = *lower_bound;
      return GSL_SUCCESS;
    }

  SAFE_FUNC_CALL (f, *upper_bound, fu);

  if (fu == 0.0)
    {
      *root = *upper_bound;
      *lower_bound = *upper_bound;
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

      /* Draw a line between f(*lower_bound) and f(*upper_bound) and
	 note where that crosses the X axis; that's where we will
	 split the interval. */

      double dx = *lower_bound - *upper_bound ;

      splitpoint = *upper_bound - (fu * dx / (fl - fu));

      SAFE_FUNC_CALL (f, splitpoint, fs);
      
      /* If the split point is the root exactly, we're done. */

      if (fs == 0.0)
	{
	  *root = splitpoint;
	  *lower_bound = splitpoint;
	  *upper_bound = splitpoint;
	  return GSL_SUCCESS;
	}
      
      /* Make note of the soon-to-be-old lower and upper bounds. */

      old_lower_bound = *lower_bound;
      old_upper_bound = *upper_bound;
      
      /* Split the interval and discard the part which doesn't contain
	 the root. */

      if ((fl > 0.0 && fs < 0.0) || (fl < 0.0 && fs > 0.0))
	{
	  *upper_bound = splitpoint;
	  fu = fs;
	  moved = UPPER;
	}
      else
	{
	  *lower_bound = splitpoint;
	  fl = fs;
	  moved = LOWER;
	}
      
      /* FIXME: this test needs help! */

      CHECK_TOL (*lower_bound, *upper_bound, rel_epsilon, abs_epsilon);

      /* If the interval has collapsed to within twice what we need,
	 the root is its midpoint. */

      if (WITHIN_TOL (*lower_bound, *upper_bound, 
		      2 * rel_epsilon, 2 * abs_epsilon))
	{
	  *root = (*lower_bound + *upper_bound) / 2.0;
	  return GSL_SUCCESS;
	}

      /* If the lower bound stayed the same and the upper bound moved less
	 than epsilon, the root is *upper_bound. */

      if (moved == UPPER && WITHIN_TOL (old_upper_bound, *upper_bound, 
					rel_epsilon, abs_epsilon))
	{
	  *root = *upper_bound;
	  return GSL_SUCCESS;
	}

      /* If the upper bound stayed the same and the lower bound moved less
	 than epsilon, the root is *lower_bound. */

      if (moved == LOWER && WITHIN_TOL (old_lower_bound, *lower_bound, 
					rel_epsilon, abs_epsilon))
	{
	  *root = *lower_bound;
	  return GSL_SUCCESS;
	}
    }
  
  GSL_ERROR ("exceeded maximum number of iterations", GSL_EMAXITER);

}
