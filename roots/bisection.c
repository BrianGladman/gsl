
/* bisection.c -- bisection root finding algorithm */

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

typedef struct
  {
    double f_lower, f_upper;
  }
bisection_state_t;

static int bisection_init (void * vstate, gsl_function * f, double * root, gsl_interval * x);
static int bisection_iterate (void * vstate, gsl_function * f, double * root, gsl_interval * x);

static int
bisection_init (void * vstate, gsl_function * f, double * root, gsl_interval * x)
{
  bisection_state_t * state = (bisection_state_t *) vstate;

  double x_lower = x->lower ;
  double x_upper = x->upper ;

  double f_lower, f_upper ;

  *root = 0.5 * (x_lower + x_upper) ;

  SAFE_FUNC_CALL (f, x_lower, &f_lower);
  SAFE_FUNC_CALL (f, x_upper, &f_upper);
  
  state->f_lower = f_lower;
  state->f_upper = f_upper;

  if ((f_lower < 0.0 && f_upper < 0.0) || (f_lower > 0.0 && f_upper > 0.0))
    {
      GSL_ERROR ("endpoints do not straddle y=0", GSL_EINVAL);
    }

  return GSL_SUCCESS;

}

static int
bisection_iterate (void * vstate, gsl_function * f, double * root, gsl_interval * x)
{
  bisection_state_t * state = (bisection_state_t *) vstate;

  double x_bisect, f_bisect;

  const double x_lower = x->lower ;
  const double x_upper = x->upper ;

  const double f_lower = state->f_lower; 
  const double f_upper = state->f_upper;

  if (f_lower == 0.0)
    {
      *root = x_lower ;
      x->upper = x_lower;
      return GSL_SUCCESS;
    }
  
  if (f_upper == 0.0)
    {
      *root = x_upper ;
      x->lower = x_upper;
      return GSL_SUCCESS;
    }
  
  x_bisect = (x_lower + x_upper) / 2.0;
  
  SAFE_FUNC_CALL (f, x_bisect, &f_bisect);
      
  if (f_bisect == 0.0)
    {
      *root = x_bisect;
      x->lower = x_bisect;
      x->upper = x_bisect;
      return GSL_SUCCESS;
    }
      
  /* Discard the half of the interval which doesn't contain the root. */
  
  if ((f_lower > 0.0 && f_bisect < 0.0) || (f_lower < 0.0 && f_bisect > 0.0))
    {
      *root = 0.5 * (x_lower + x_bisect) ;
      x->upper = x_bisect;
      state->f_upper = f_bisect;
    }
  else
    {
      *root = 0.5 * (x_bisect + x_upper) ;
      x->lower = x_bisect;
      state->f_lower = f_bisect;
    }

  return GSL_SUCCESS;
}


static const gsl_root_fsolver_type bisection_type =
{"bisection",				/* name */
 sizeof (bisection_state_t),
 &bisection_init,
 &bisection_iterate};

const gsl_root_fsolver_type  * gsl_root_fsolver_bisection = &bisection_type;
