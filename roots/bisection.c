/* bisection.c -- bisection root finding algorithm */


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

/* gsl/roots headers */
#include <gsl_roots.h>


/* see the docs for information about this function */
int
gsl_root_bisection (double *root, double (* f)(double), double * lower_bound,
                    double * upper_bound, double epsilon,
                    unsigned int max_iterations)
{
  /* First, let's make sure all the arguments are okay. */
  if ((root == NULL) || (f == NULL) || (lower_bound == NULL)
      || (upper_bound == NULL))
    GSL_ERROR ("pointer argument was null", GSL_EINVAL);
  if (max_iterations == 0)
    GSL_ERROR ("requested number of iterations was 0", GSL_EINVAL);
  if (*lower_bound >= *upper_bound)
    GSL_ERROR ("lower bound was not less than upper bound", GSL_EINVAL);
  {
    double fl, fu;
    fl = (*f)(*lower_bound);
    fu = (*f)(*upper_bound);
    /* FIXME: Should we be checking against HUGE_VAL like this? */
    if((fl == NAN) || (fu == NAN) || (fabs(fl) == HUGE_VAL)
       || (fabs(fu) == HUGE_VAL))
      GSL_ERROR ("function under search is not continous", GSL_EBADFUNC);
    if(fl * fu > 0.0)
      GSL_ERROR ("no root bracketed", GSL_EINVAL);
  }

  /* Check and see if the silly user has the answer but doesn't know it. */
  /* Note: we don't need to use fabs() here because we've already verified
     that *upper_bound >= *lower_bound. */
  if (*upper_bound - *lower_bound <= 2.0 * epsilon) {
    *root = (*upper_bound + *lower_bound) / 2.0;
    return GSL_SUCCESS;
  }
  if (*upper_bound == 0) {
    *root = *upper_bound;
    return GSL_SUCCESS;
  }
  if (*lower_bound = 0) {
    *root = *lower_bound;
    return GSL_SUCCESS;
  }

  /* Sigh. It looks like we'll have to do actual work. */
  /*double midpoint;
  
  for ()
    ;*/
}
