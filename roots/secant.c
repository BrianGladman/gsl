/* secant.c -- secant method root finding algorithm */

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
gsl_root_secant (double *root, 
		 double (*f) (double), 
		 double *guess1, double *guess2, 
		 double rel_epsilon, double abs_epsilon,
		 unsigned int max_iterations)
{
  unsigned int iterations;
  double new_guess, fg1, fg2, fnew;

  if (rel_epsilon < 0.0 || abs_epsilon < 0.0)
    GSL_ERROR ("relative or absolute tolerance negative", GSL_EBADTOL);

  if (rel_epsilon < GSL_DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER)
    GSL_ERROR ("relative tolerance too small", GSL_EBADTOL);

  SAFE_FUNC_CALL (f, *guess1, fg1);

  if (fg1 == 0.0)
    {
      *root = *guess1;
      return GSL_SUCCESS;
    }

  SAFE_FUNC_CALL (f, *guess2, fg2);

  if (fg2 == 0.0)
    {
      *root = *guess2;
      return GSL_SUCCESS;
    }

  for (iterations = 0; iterations < max_iterations; iterations++)
    {
      /* Draw a line through f(*guess1) and f(*guess2) and note where that
	 crosses the X axis; that's our new guess. */

      if (fg1 == fg2) 
	{
	  GSL_ERROR("secant approximation to derivative is 0", GSL_EZERODIV);
	}
							       
      new_guess = *guess2 - (fg2 * (*guess1 - *guess2) / (fg1 - fg2));
      
      SAFE_FUNC_CALL (f, new_guess, fnew);
      
      /* If the new point is the root exactly, we're done. */

      if (fnew == 0.0)
	{
	  *root = new_guess;
	  return GSL_SUCCESS;
	}
      
      /* Rotate the guesses. */
      *guess2 = *guess1;
      fg2 = fg1;
      *guess1 = new_guess;
      fg1 = fnew;

      if (WITHIN_TOL (*guess1, *guess2, rel_epsilon, abs_epsilon))
	{
	  *root = *guess1;
	  return GSL_SUCCESS;
	}
      else
	{
	  CHECK_TOL (*guess1, *guess2, rel_epsilon, abs_epsilon);
	}
    }
  
  GSL_ERROR ("exceeded maximum number of iterations", GSL_EMAXITER);

}
