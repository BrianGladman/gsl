/* newton.c -- Newton's Method root finding algorithm */

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
gsl_root_newton (double *root,
		 double (*f) (double),
		 double (*df) (double),
		 void (*fdf) (double, double *, double *),
		 double *guess, 
		 double rel_epsilon, double abs_epsilon,
		 unsigned int max_iterations)
{
  double new_guess, fg, dfg, fnew, dfnew;
  unsigned int iterations;
    
  if (rel_epsilon < 0.0 || abs_epsilon < 0.0)
    GSL_ERROR ("relative or absolute tolerance negative", GSL_EBADTOL);

  if (rel_epsilon < DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER)
    GSL_ERROR ("relative tolerance too small", GSL_EBADTOL);
  
  SAFE_FUNC_CALL (f, *guess, fg);
  
  if (fg == 0.0)
    {
      *root = *guess;
      return GSL_SUCCESS;
    }
  
  SAFE_FUNC_CALL (df, *guess, dfg) ;
  
  for (iterations = 0; iterations < max_iterations; iterations++)
    {
      
      /* Draw a line tangent to f at guess and note where that crosses the X
	 axis; that's our new guess. */

      if (dfg == 0)
	{
	  GSL_ERROR("derivative is zero", GSL_EZERODIV);
	}

      new_guess = *guess - (fg / dfg);
      
      if (fdf) 
	{
	  (*fdf) (new_guess, &fnew, &dfnew);
	}
      else
	{
	  fnew = (*f) (new_guess);
	  dfnew = (*df) (new_guess) ;
	}
      
      if (!GSL_ISREAL (fnew))
	{
	  GSL_ERROR ("function not continuous", GSL_EBADFUNC);
	}
      
      if (fnew == 0.0)
	{
	  *root = new_guess;
	  return GSL_SUCCESS;
	}
      
      /* Now, let's check if we're finished. FIXME.2: In certain fairly
	 unusual cases, this test can cause the root finder to stop early. */

      if (WITHIN_TOL (*guess, new_guess, rel_epsilon, abs_epsilon))
	{
	  *root = new_guess;
	  return GSL_SUCCESS;
	}
      else
	{
	  CHECK_TOL (*guess, new_guess, rel_epsilon, abs_epsilon);
	}

      if (!GSL_ISREAL (dfnew))
	{
	  GSL_ERROR ("function not differentiable", GSL_EBADFUNC);
	}
      
      /* Rotate the guesses. */

      *guess = new_guess;
      fg = fnew;
      dfg = dfnew;
    }
  
  GSL_ERROR ("exceeded maximum number of iterations", GSL_EMAXITER);
}
