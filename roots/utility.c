/* utility.c -- various root finding utility routines */
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


/* Validate arguments common to gsl_root_bisection and gsl_root_falsepos.
   Return GSL_SUCCESS if all arguments are okay, complain appropriately (i.e.
   call GSL_ERROR and return GSL_FAILURE) otherwise. */
int
_gsl_root_validate_bfp_args(void * root, void * f, double * lower_bound,
                            double * upper_bound, double rel_epsilon,
                            double abs_epsilon, unsigned int max_iterations,
                            double max_deltay)
{
  /* Is the maximum delta-y too small? */
  if (max_deltay < GSL_ROOT_MIN_MAX_DELTAY)
    GSL_ERROR("invalid argument(s)", GSL_EINVAL);

  /* The rest of the arguments are common. */
  if (_gsl_root_validate_args_p(root, f, lower_bound, upper_bound,
                                rel_epsilon, abs_epsilon, max_iterations))
    return GSL_SUCCESS;
  else
    GSL_ERROR("invalid argument(s)", GSL_EINVAL);
}

/* Validate the arguments common to all four low level functions. Return false
   if any of the following conditions do not hold, otherwise return true.

   * No pointer arguments are null.
   * The maximum number of iterations is non-zero.
   * Relative and absolute error are non-negative.
   * The relative error is not too small.
   * The upper bound is larger than the lower bound. */
int
_gsl_root_validate_args_p(void * root, void * f, double * lower_bound,
                          double * upper_bound, double rel_epsilon,
                          double abs_epsilon, unsigned int max_iterations)
{
  /* Are any pointers null? */
  if ((root == NULL) || (f == NULL) || (lower_bound == NULL)
      || (upper_bound == NULL))
    return 0;
  /* Did the user tell us to do no iterations? */
  if (max_iterations == 0)
    return 0;
  /* Did the user try to pawn a negative tolerance off on us? */
  if (rel_epsilon < 0.0 || abs_epsilon < 0.0)
    return 0;
  /* Is the relative error too small? */
  if (rel_epsilon < DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER)
    return 0;
  /* Did the user give a lower bound that not less than the upper bound? */
  if (*lower_bound >= *upper_bound)
    return 0;

  /* All is well. */
  return 1;
}

/* Verify that the supplied interval is guaranteed by the Intermediate Value
   Theorem to contain a root and complain appropriately if it is not. (It
   might actually be a discontinuity, but we check for that elsewhere.) Return
   GSL_SUCCESS if all is well, otherwise, call GSL_ERROR and return
   GSL_FAILURE. */
int
_gsl_root_ivt_guar(double (* f)(double), double lower_bound,
                     double upper_bound)
{
    double fl, fu;

    _BARF_FPCALL(f, lower_bound, fl);
    _BARF_FPCALL(f, upper_bound, fu);

    if (fl * fu > 0.0)
      GSL_ERROR("interval not guaranteed to contain a root", GSL_EINVAL);
    else
      return GSL_SUCCESS;
}

/* Check if the user has the root but doesn't know it. If lower_bound or
   upper_bound is a root of f, or the interval [upper_bound, lower_bound] is
   within tolerance, return 1 and set *root appropriately. Otherwise, return
   0. On error, call GSL_ERROR and return GSL_FAILURE. */
int
_gsl_root_silly_user(double * root, double (* f)(double), double lower_bound,
                     double upper_bound, double rel_epsilon,
                       double abs_epsilon, double max_deltay)
{
  double fl, fu;
  
  /* Is lower_bound the root? */
  _BARF_FPCALL(f, lower_bound, fl);
  if (fl == 0.0) {
    *root = lower_bound;
    return 1;
  }

  /* Is upper_bound the root? */
  _BARF_FPCALL(f, upper_bound, fu);
  if (fu == 0.0) {
    *root = upper_bound;
    return 1;
  }
  
  /* Are lower_bound and upper_bound within tolerance? */
  _BARF_TOLS(lower_bound, upper_bound, 2 * rel_epsilon, 2 * abs_epsilon);
  _BARF_DELTAY(fl, fu, max_deltay);
  if (_WITHIN_TOL(lower_bound, upper_bound, 2 * rel_epsilon,
                  2 * abs_epsilon)) {
    *root = (lower_bound + upper_bound) / 2.0;
    return 1;
  }

  /* No? Bummer. */
  return 0;
}
