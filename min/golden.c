
/* goldensection.c -- goldensection minimum finding algorithm */

#include <config.h>

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>

#include "min.h"

typedef struct
  {
    double dummy;
  }
goldensection_state_t;

static int goldensection_init (void * vstate, gsl_function * f, double minimum, double f_minimum, gsl_interval x, double f_lower, double f_upper);
static int goldensection_iterate (void * vstate, gsl_function * f, double * minimum, double * f_minimum, gsl_interval * x, double * f_lower, double * f_upper);

static int
goldensection_init (void * vstate, gsl_function * f, double minimum, double f_minimum, gsl_interval x, double f_lower, double f_upper)
{
  goldensection_state_t * state = (goldensection_state_t *) vstate;

  /* no initialization require, prevent warnings about unused variables */

  state = 0;
  f = 0;
  minimum = 0;
  f_minimum = 0;
  x.lower = 0;
  f_lower = 0;
  x.upper = 0;
  f_upper = 0;

  return GSL_SUCCESS;
}

static int
goldensection_iterate (void * vstate, gsl_function * f, double * minimum, double * f_minimum, gsl_interval * x, double * f_lower, double * f_upper)
{
  goldensection_state_t * state = (goldensection_state_t *) vstate;

  const double x_minimum = *minimum ;
  const double x_lower = x->lower ;
  const double x_upper = x->upper ;

  const double f_min = *f_minimum;

  const double golden = 0.318966 ; /* golden = (3 - sqrt(5))/2 */
  
  const double w_lower = (x_minimum - x_lower);
  const double w_upper = (x_upper - x_minimum);

  double x_new, f_new;

  state = 0 ; /* avoid warning about unused parameters */
  
  x_new = x_minimum + golden * ((w_upper > w_lower) ? w_upper : -w_lower) ;

  SAFE_FUNC_CALL (f, x_new, &f_new);

  if (f_new < f_min)
    {
      *minimum = x_new ;
      *f_minimum = f_new ;
      return GSL_SUCCESS;
    }
  else if (x_new < x_minimum && f_new > f_min)
    {
      x->lower = x_new ;
      *f_lower = f_new ;
      return GSL_SUCCESS;
    }
  else if (x_new > x_minimum && f_new > f_min)
    {
      x->upper = x_new ;
      *f_upper = f_new ;
      return GSL_SUCCESS;
    }
  else
    {
      return GSL_FAILURE;
    }
}


static const gsl_min_fminimizer_type goldensection_type =
{"goldensection",				/* name */
 sizeof (goldensection_state_t),
 &goldensection_init,
 &goldensection_iterate};

const gsl_min_fminimizer_type  * gsl_min_fminimizer_goldensection = &goldensection_type;
