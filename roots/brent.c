/* brent.c -- brent root finding algorithm */

#include <config.h>

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

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
  
  for (iterations = 0; iterations < max_iterations; iterations++)
    {
      if ((fb < 0 && fc < 0.0) || (fb > 0 && fc > 0))
	{
	  c = a;
	  fc = fa;
	  d = b - a;
	  e = d ;
	}
      
      if (fabs(fc) < fabs(fb))
	{
	  a = b ;
	  b = c ; 
	  c = a ;
	  fa = fb ;
	  fb = fc ;
	  fc = fa ; 
	}

      tol = 2 * rel * fabs(b) + 0.5 * abs
      err = 0.5 * (c - b);

      if (fabs(err) <= tol || fb == 0)
	{
	  *root = b;
	  return GSL_SUCCESS;
	}

      if (fabs(e) >= tol && fabs(fa) > fabs(fb))
	{
	  s = fb / fa;
	  if (a == c)
	    {
	      p = 2 * xm * s;
	      q = 1 - s;
	    }
	  else
	    {
	      q = fa / fc;
	      r = fb / fc;
	      p = s * (2 * xm * q * (q - r) - (b - a) * (r - 1)) ;
	      q = (q - 1) * (r - 1) * (s - 1);
	    }
	  if (p <= 0)
	    {
	      p = -p ;
	    }
	  else
	    {
	      q = -q ;
	    }
	  
	  s = e ;

	  if (2 * p < MIN(3 * xm * q - fabs(tol * q), fabs(s * q)))
	    {
	      e = d ;
	      d = p/q ;
	    }
	  else
	    {
	      d = xm ;
	      e = d ;
	    }
	}

      a = b ;
      fa = fb ;

      if (fabs(d) > tol)
	{
	  b += d;
	}
      else
	{
	  b += (xm < 0 ? -tol1 : +tol1) ;
	}

      SAFE_FUNC(f, *b, fb);
	
    }
  
  GSL_ERROR ("exceeded maximum number of iterations", GSL_EMAXITER);
  
}
