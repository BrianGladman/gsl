/* newton.c -- Newton's Method root finding algorithm */
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

/* See the documentation for information about this function. */
int
gsl_root_newton(double * root,
                double (* fdf)(double *, double *, double, int, int),
                double * guess, double rel_epsilon, double abs_epsilon,
                unsigned int max_iterations, double max_step_size)
{
  /* Validate input. */
  {
    double y, yprime;
    
    /* Validate arguments. */  /* FIXME: can't use the same check for fdf */

    if (_gsl_root_validate_sn_args(root, fdf, guess, guess, rel_epsilon,
                                   abs_epsilon, max_iterations, max_step_size))
      return GSL_FAILURE; /* GSL_ERROR was already called. */
    
    /* Did the user unwittingly provide the root? */
    _BARF_FDFPCALL(fdf, &y, &yprime, *guess, WANTED, !WANTED);
    if (y == 0.0) {
      *root = *guess;
      return GSL_SUCCESS;
    }
  }

  /* Doh! It looks like we'll have to do actual work. */
  {
    double new_guess, fg, dfg, fnew, dfnew;
    unsigned int iterations;
  
    _BARF_FDFPCALL(fdf, &fg, &dfg, *guess, WANTED, WANTED);
  
    for (iterations = 0; iterations < max_iterations; iterations++) {
      /* Draw a line tangent to f at guess and note where that crosses the X
         axis; that's our new guess. */
      _BARF_ZERO(dfg);

      new_guess = *guess - (fg / dfg);

      (*fdf)(&fnew, &dfnew, new_guess, 1, 1); 

      if (!GSL_ISREAL(fnew))
	{
	  GSL_ERROR("function not continuous", GSL_EBADFUNC);
	}

      /* If the new point is the root exactly, we're done. */
      if (fnew == 0.0) {
        *root = new_guess;
        return GSL_SUCCESS;
      }

      /* Make sure we're staying on Earth. */
      _BARF_DELTAX(*guess, new_guess, max_step_size);

      /* Now, let's check if we're finished. FIXME.2: In certain fairly
         unusual cases, this test can cause the root finder to stop early. */
      _BARF_TOLS(*guess, new_guess, rel_epsilon, abs_epsilon);
      if (_WITHIN_TOL(*guess, new_guess, rel_epsilon, abs_epsilon)) {
        *root = new_guess;
        return GSL_SUCCESS;
      }

      if (!GSL_ISREAL(dfnew))
	{
	  GSL_ERROR("function not differentiable", GSL_EBADFUNC);
	}

      /* Rotate the guesses. */
      *guess = new_guess;
      fg = fnew;
      dfg = dfnew;
    }
    
    /* Uh oh, ran out of iterations. */
    GSL_ERROR("exceeded maximum number of iterations", GSL_EMAXITER);
  }
}
