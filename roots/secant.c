/* secant.c -- secant method root finding algorithm */
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

/* See the documentation for more information about this function. */
int
gsl_root_secant(double * root, double (* f)(double), double * guess1,
                double * guess2, double rel_epsilon, double abs_epsilon,
                unsigned int max_iterations, double max_step_size)
{
  /* Validate input. */
  {
    int status;
    
    /* Validate arguments. */
    if (_gsl_root_validate_sm_args(root, f, guess1, guess2, rel_epsilon,
                                   abs_epsilon, max_iterations, max_step_size))
      return GSL_FAILURE; /* GSL_ERROR was already called. */
    
    /* Check if the user unwittingly provided the root. If _g_r_s_u returns
       non-zero, it has taken care of things and we should return what it
       returned. */
    if ((status = _gsl_root_silly_user(root, f, *guess1, *guess2, rel_epsilon,
                                       abs_epsilon, _IGNORE_DELTAY)))
      return status;
  }

  /* Doh! It looks like we'll have to do actual work. */
  {
    double new_guess, fg1, fg2, fnew;
    unsigned int iterations;
  
    _BARF_FPCALL(f, *guess1, fg1);
    _BARF_FPCALL(f, *guess2, fg2);
  
    for (iterations = 0; iterations < max_iterations; iterations++) {
      /* Draw a line through f(*guess1) and f(*guess2) and note where that
         crosses the X axis; that's our new guess. */
      _BARF_ZERO(fg1 - fg2);
      new_guess = *guess2 - (fg2 * (*guess1 - *guess2) / (fg1 - fg2));
      _BARF_FPCALL(f, new_guess, fnew);

      /* If the new point is the root exactly, we're done. */
      if (fnew == 0.0) {
        *root = new_guess;
        return GSL_SUCCESS;
      }

      /* Make sure we're staying on Earth. */
      _BARF_DELTAX(*guess1, new_guess, max_step_size);

      /* Rotate the guesses. */
      *guess2 = *guess1;
      fg2 = fg1;
      *guess1 = new_guess;
      fg1 = fnew;
      
      /* Now, let's check if we're finished. FIXME.1: this test needs help! */
      _BARF_TOLS(*guess1, *guess2, rel_epsilon, abs_epsilon);
      if (_WITHIN_TOL(*guess1, *guess2, rel_epsilon, abs_epsilon)) {
        *root = *guess1;
        return GSL_SUCCESS;
      }
    }
    
    /* Uh oh, ran out of iterations. */
    GSL_ERROR("exceeded maximum number of iterations", GSL_ETIMEOUT);
  }
}
