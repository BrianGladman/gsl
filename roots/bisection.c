/* bisection.c -- bisection root finding algorithm */
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

/* Note that when determining if we're finished, this function checks against
   2 * epsilon, not epsilon. This is because if the interval is twice as wide
   as we want, we can simply return the midpoint of the interval, and that
   will be within the desired accuracy. */
/* See the documentation for more information about this function. */
int
gsl_root_bisection(double * root, double (* f)(double), double * lower_bound,
                   double * upper_bound, double epsilon,
                   unsigned int max_iterations)
{
  /* First, let's make sure all the arguments are okay. */
  /* Are any pointers null? */
  if ((root == NULL) || (f == NULL) || (lower_bound == NULL)
      || (upper_bound == NULL))
    GSL_ERROR("pointer argument null", GSL_EINVAL);
  /* Did the user tell us to do no iterations? */
  if (max_iterations == 0)
    GSL_ERROR("requested 0 iterations", GSL_EINVAL);
  /* Did the user ask for more accuracy than we can give or a negative
     epsilon? */
  if (epsilon < GSL_ROOT_EPSILON_BUFFER * DBL_EPSILON)
    GSL_ERROR("requested accuracy unattainable", GSL_EINVAL);
  /* Did the user give a lower bound that is greater than the upper bound? */
  if (*lower_bound >= *upper_bound)
    GSL_ERROR("lower bound not less than upper bound", GSL_EINVAL);
  /* Did the user give an interval which might not contain a root? */
  { 
    double fl, fu;
    _GSL_ROOT_FPCALL(f, *lower_bound, fl);
    _GSL_ROOT_FPCALL(f, *upper_bound, fu);
    if (fl * fu > 0.0)
      GSL_ERROR("interval not guaranteed to contain a root", GSL_EINVAL);
  }

  /* Check if the silly user has the answer but doesn't know it. */
  {
    double fl, fu;

    _GSL_ROOT_FPCALL(f, *lower_bound, fl);
    if (fl == 0.0) {
      *root = *lower_bound;
      return GSL_SUCCESS;
    }

    _GSL_ROOT_FPCALL(f, *upper_bound, fu);
    if (fu == 0.0) {
      *root = *upper_bound;
      return GSL_SUCCESS;
    }

    if (_GSL_ROOT_ERR(*lower_bound, *upper_bound) <= 2.0 * epsilon) {
      *root = (*upper_bound + *lower_bound) / 2.0;
      return GSL_SUCCESS;
    }
  }

  /* Doh! It looks like we'll have to do actual work. */
  {
    double midpoint, fl, fu, fm;
    unsigned int iterations;
  
    /* Evaluate the function under search at lower_bound and upper_bound. */
    _GSL_ROOT_FPCALL(f, *lower_bound, fl);
    _GSL_ROOT_FPCALL(f, *upper_bound, fu);
  
    for (iterations = 0; iterations < max_iterations; iterations++) {
      /* Chop the interval in half. */
      midpoint = (*upper_bound + *lower_bound) / 2.0;
      _GSL_ROOT_FPCALL(f, midpoint, fm);
      
      /* If the midpoint is the root exactly, we're done. */
      if (fm == 0.0) {
        *root = midpoint;
        return GSL_SUCCESS;
      }
      
      /* Discard the half of the interval which doesn't contain the root. */
      if (fl * fm < 0.0) {
        *upper_bound = midpoint;
        fu = fm;
      }
      /* If the root is not between *lower_bound and midpoint, it is guaranteed
         to be between midpoint and *upper_bound. There is definitely a root
         between *lower_bound and *upper_bound, and f(*lower_bound),
         f(*upper_bound), and f(midpoint) are all non-zero. */
      else {
        *lower_bound = midpoint;
        fl = fm;
      }
      
      /* Now, let's check if we're finished. */
      if (_GSL_ROOT_ERR(*lower_bound, *upper_bound) < 2.0 * epsilon) {
        *root = (*upper_bound + *lower_bound) / 2.0;
        return GSL_SUCCESS;
      }
    }

    /* Uh oh, ran out of iterations. */
    GSL_ERROR("exceeded maximum number of iterations", GSL_ETIMEOUT);
  }
}
