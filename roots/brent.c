/* brent.c -- brent root finding algorithm 

   The original paper is,

   "An algorithm with guaranteed convergence for finding a zero of a
   function", R. P. Brent, The Computer Journal, Volume 14, Number 4,
   pages 422-425.

*/

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

int
gsl_root_brent (double *root, double (*f) (double),
		double *lower_bound, double *upper_bound,
		double rel_epsilon, double abs_epsilon,
		unsigned int max_iterations)
{
  unsigned int iterations;
  double fu, fl, a, b, c, d, e, tol, fa, fb, fc, m;

  if (*lower_bound >= *upper_bound)
    GSL_ERROR ("lower bound larger than upper_bound", GSL_EINVAL);

  if (rel_epsilon < 0.0 || abs_epsilon < 0.0)
    GSL_ERROR ("relative or absolute tolerance negative", GSL_EBADTOL);

  if (rel_epsilon < GSL_DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER)
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

  a = *lower_bound;
  fa = fl;

  b = *upper_bound;
  fb = fu;

  c = b;
  fc = fb;

  d = b - a;
  e = b - a;

  for (iterations = 0; iterations < max_iterations; iterations++)
    {
      int ac_equal = 0;

      if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0))
	{
	  ac_equal = 1;
	  c = a;
	  fc = fa;
	  d = b - a;
	  e = b - a;
	}

      if (fabs (fc) < fabs (fb))
	{
	  ac_equal = 1;
	  a = b;
	  b = c;
	  c = a;
	  fa = fb;
	  fb = fc;
	  fc = fa;
	}

      tol = 0.5 * GSL_MAX (rel_epsilon * fabs (b), abs_epsilon);
      m = 0.5 * (c - b);

      if (fb == 0)
	{
	  *root = b;
	  *lower_bound = b;
	  *upper_bound = b;

	  return GSL_SUCCESS;
	}
      
      if (fabs (m) <= tol)
	{
	  *root = b;
	  *lower_bound = b;
	  *upper_bound = c;

	  return GSL_SUCCESS;
	}

      if (fabs (e) < tol || fabs (fa) <= fabs (fb))
	{
	  d = m;		/* use bisection */
	  e = m;
	}
      else
	{
	  double p, q, r;	/* use inverse cubic interpolation */
	  double s = fb / fa;

	  if (ac_equal)
	    {
	      p = 2 * m * s;
	      q = 1 - s;
	    }
	  else
	    {
	      q = fa / fc;
	      r = fb / fc;
	      p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
	      q = (q - 1) * (r - 1) * (s - 1);
	    }

	  if (p > 0)
	    {
	      q = -q;
	    }
	  else
	    {
	      p = -p;
	    }

	  if (2 * p < GSL_MIN (3 * m * q - fabs (tol * q), fabs (e * q)))
	    {
	      e = d;
	      d = p / q;
	    }
	  else
	    {
	      /* interpolation failed, fall back to bisection */

	      d = m;
	      e = m;
	    }
	}

      a = b;
      fa = fb;

      if (fabs (d) > tol)
	{
	  b += d;
	}
      else
	{
	  b += (m > 0 ? +tol : -tol);
	}

      SAFE_FUNC_CALL (f, b, fb);

      /* Update the best estimate of the root and bounds on each
         iteration */

      *root = b;

      if (b < c)
	{
	  *lower_bound = b;
	  *upper_bound = c;
	}
      else
	{
	  *lower_bound = c;
	  *upper_bound = b;
	}
    }

  GSL_ERROR ("exceeded maximum number of iterations", GSL_EMAXITER);

}
