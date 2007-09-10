/* roots/utility.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Reid Priedhorsky, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* utility.c -- various root finding utility routines */

/* config headers */
#include <config.h>

/* standard headers */
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* gsl headers */
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

/* roots headers */
#include "roots.h"


/* Validate arguments common to gsl_root_bisection and gsl_root_falsepos.
   Return GSL_SUCCESS if all arguments are okay, complain appropriately (i.e.
   call GSL_ERROR and return GSL_FAILURE) otherwise. */
int
_gsl_root_validate_bfp_args (double *root, double (*f) (double),
                             double *lower_bound,
                             double *upper_bound, double epsrel,
                             double epsabs, unsigned int max_iterations,
                             double max_deltay)
{
  /* Is the maximum delta-y too small? */
  if (max_deltay < GSL_ROOT_MIN_MAX_DELTAY)
    GSL_ERROR ("maximum delta-y negative, zero, or too small", GSL_EBADTOL);

  /* Did the user give a lower bound that not less than the upper bound? */
  if (*lower_bound >= *upper_bound)
    GSL_ERROR ("lower bound larger than upper_bound", GSL_EINVAL);

  /* The rest of the arguments are common. */
  return _gsl_root_validate_args (root, f, lower_bound, upper_bound,
                                  epsrel, epsabs, max_iterations);
}

/* Validate arguments commond to gsl_root_secant and gsl_root_newtons. Return
   GSL_SUCCESS if all arguments are okay, complain appropriately (i.e. call
   GSL_ERROR and return GSL_FAILURE) otherwise. */
int
_gsl_root_validate_sn_args (double *root, double (*f) (double),
                            double *guess1,
                            double *guess2, double epsrel,
                            double epsabs, unsigned int max_iterations,
                            double max_step_size)
{
  /* Is the maximum step size ridiculous? */
  if (max_step_size <= 0.0)
    GSL_ERROR ("maximum step size <= 0", GSL_EBADTOL);

  /* The rest of the arguments are common. */
  return _gsl_root_validate_args (root, f, guess1, guess2, epsrel,
                                  epsabs, max_iterations);
}

/* Validate the arguments common to all four low level functions. Return
   GSL_SUCCESS if all of the following arguments hold true, call GSL_ERROR and
   return GSL_FAILURE otherwise.

   * No pointer arguments are null.
   * The maximum number of iterations is non-zero.
   * Relative and absolute error are non-negative.
   * The relative error is not too small. */
int
_gsl_root_validate_args (double *root, double (*f) (double), double *where1,
                         double *where2, double epsrel,
                         double epsabs, unsigned int max_iterations)
{
  /* Are any pointers null? */
  if ((root == NULL) || (f == NULL) || (where1 == NULL)
      || (where2 == NULL))
    GSL_ERROR ("pointer argument null", GSL_EINVAL);
  /* Did the user tell us to do no iterations? */
  if (max_iterations == 0)
    GSL_ERROR ("maximum iterations 0", GSL_EINVAL);
  /* Did the user try to pawn a negative tolerance off on us? */
  if (epsrel < 0.0 || epsabs < 0.0)
    GSL_ERROR ("relative or absolute tolerance negative", GSL_EBADTOL);
  /* Is the relative error too small? */
  if (epsrel < GSL_DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER)
    GSL_ERROR ("relative tolerance too small", GSL_EBADTOL);

  /* All is well. */
  return GSL_SUCCESS;
}

/* Verify that the supplied interval is guaranteed by the Intermediate Value
   Theorem to contain a root and complain appropriately if it is not. (It
   might actually be a discontinuity, but we check for that elsewhere.) Return
   GSL_SUCCESS if all is well, otherwise, call GSL_ERROR and return
   GSL_FAILURE. */
int
_gsl_root_ivt_guar (double (*f) (double), double lower_bound,
                    double upper_bound)
{
  double fl, fu;

  _BARF_FPCALL (f, lower_bound, fl);
  _BARF_FPCALL (f, upper_bound, fu);

  if (fl * fu > 0.0)
    {
      GSL_ERROR ("interval not guaranteed to contain a root", GSL_EINVAL);
    }
  else
    {
      return GSL_SUCCESS;
    }
}

/* Check if the user has the root but doesn't know it. If lower_bound or
   upper_bound is a root of f, or the interval [upper_bound, lower_bound] is
   within tolerance, return 1 and set *root appropriately. Otherwise, return
   0. On error, call GSL_ERROR and return GSL_FAILURE. Only worry about
   max_deltay if it is greater than 0 (which implies that you should validate
   arguments _before_ calling this function). */
int
_gsl_root_silly_user (double *root, double (*f) (double), double lower_bound,
                      double upper_bound, double epsrel,
                      double epsabs, double max_deltay)
{
  double fl, fu;

  /* Is lower_bound the root? */
  _BARF_FPCALL (f, lower_bound, fl);
  if (fl == 0.0)
    {
      *root = lower_bound;
      return 1;
    }

  /* Is upper_bound the root? */
  _BARF_FPCALL (f, upper_bound, fu);
  if (fu == 0.0)
    {
      *root = upper_bound;
      return 1;
    }

  /* Are lower_bound and upper_bound within tolerance? */
  _BARF_TOLS (lower_bound, upper_bound, 2 * epsrel, 2 * epsabs);
  if (max_deltay > 0.0)
    _BARF_DELTAY (fl, fu, max_deltay);
  if (_WITHIN_TOL (lower_bound, upper_bound, 2 * epsrel,
                   2 * epsabs))
    {
      *root = (lower_bound + upper_bound) / 2.0;
      return 1;
    }

  /* No? Bummer. */
  return 0;
}
