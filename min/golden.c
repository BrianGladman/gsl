
/* goldensection.c -- goldensection minimum finding algorithm */

#include <config.h>

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <gsl_math.h>
#include <gsl_errno.h>
#include <gsl_min.h>

#include "min.h"

typedef struct
  {
    double f_lower, f_upper, f_minimum;
  }
goldensection_state_t;

int goldensection_init (void * vstate, gsl_function * f, double * minimum, gsl_interval * x);
int goldensection_iterate (void * vstate, gsl_function * f, double * minimum, gsl_interval * x);

int
goldensection_init (void * vstate, gsl_function * f, double * minimum, gsl_interval * x)
{
  goldensection_state_t * state = (goldensection_state_t *) vstate;

  double x_lower = x->lower ;
  double x_upper = x->upper ;

  double f_lower, f_upper, f_minimum ;

  SAFE_FUNC_CALL (f, x_lower, &f_lower);
  SAFE_FUNC_CALL (f, x_upper, &f_upper);
  SAFE_FUNC_CALL (f, *minimum, &f_minimum);
  
  state->f_lower = f_lower;
  state->f_upper = f_upper;
  state->f_minimum = f_minimum;

  if (f_minimum >= f_lower || f_minimum >= f_upper)
    {
      GSL_ERROR ("endpoints do not enclose a minimum", GSL_EINVAL);
    }

  return GSL_SUCCESS;

}

int
goldensection_iterate (void * vstate, gsl_function * f, double * minimum, gsl_interval * x)
{
  goldensection_state_t * state = (goldensection_state_t *) vstate;

  const double x_minimum = *minimum ;
  const double x_lower = x->lower ;
  const double x_upper = x->upper ;

  const double f_minimum = state->f_minimum;

  const double golden = 0.318966 ; /* golden = (3 - sqrt(5))/2 */
  
  const double w_lower = (x_minimum - x_lower);
  const double w_upper = (x_upper - x_minimum);

  double x_new, f_new;
  
  x_new = x_minimum + golden * ((w_upper > w_lower) ? w_upper : -w_lower) ;

  SAFE_FUNC_CALL (f, x_new, &f_new);

  if (f_new < f_minimum)
    {
      *minimum = x_new ;
      state->f_minimum = f_new ;
      return GSL_SUCCESS;
    }
  else if (x_new < x_minimum && f_new > f_minimum)
    {
      x->lower = x_new ;
      state->f_lower = f_new ;
      return GSL_SUCCESS;
    }
  else if (x_new > x_minimum && f_new > f_minimum)
    {
      x->upper = x_new ;
      state->f_upper = f_new ;
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
