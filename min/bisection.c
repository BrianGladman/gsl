
/* bisection.c -- bisection minimum finding algorithm */

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
    double f_lower, f_upper, f_middle;
  }
bisection_state_t;

int bisection_init (void * vstate, gsl_function * f, double * minimum, gsl_interval * x);
int bisection_iterate (void * vstate, gsl_function * f, double * minimum, gsl_interval * x);

int
bisection_init (void * vstate, gsl_function * f, double * minimum, gsl_interval * x)
{
  bisection_state_t * state = (bisection_state_t *) vstate;

  double x_lower = x->lower ;
  double x_upper = x->upper ;

  double f_lower, f_upper, f_middle ;

  SAFE_FUNC_CALL (f, x_lower, &f_lower);
  SAFE_FUNC_CALL (f, x_upper, &f_upper);
  SAFE_FUNC_CALL (f, *minimum, &f_middle);
  
  state->f_lower = f_lower;
  state->f_upper = f_upper;
  state->f_middle = f_middle;

  if (f_middle >= f_lower || f_middle >= f_upper)
    {
      GSL_ERROR ("endpoints do not enclose a minimum", GSL_EINVAL);
    }

  return GSL_SUCCESS;

}

int
bisection_iterate (void * vstate, gsl_function * f, double * minimum, gsl_interval * x)
{
  bisection_state_t * state = (bisection_state_t *) vstate;

  const double x_lower = x->lower ;
  const double x_upper = x->upper ;

  const double f_middle = state->f_middle;

  const double golden = 0.318966 ;
  
  const double w_lower = (*minimum - x_lower);
  const double w_upper = (x_upper - *minimum);

  if (w_lower > w_upper) 
    {
      double x_new = *minimum - golden * w_lower ;

      double f_new;
      SAFE_FUNC_CALL (f, x_new, &f_new);

      if (f_new < f_middle)
        {
          *minimum = x_new ;
          state->f_middle = f_new ;
        }
      else 
        {
          x->lower = x_new ;
          state->f_lower = f_new ;
        }
    } 
  else
    {
      double x_new = *minimum + golden * w_upper ;

      double f_new;
      SAFE_FUNC_CALL (f, x_new, &f_new);

      if (f_new < f_middle)
        {
          *minimum = x_new ;
          state->f_middle = f_new ;
        }
      else 
        {
          x->upper = x_new ;
          state->f_upper = f_new ;
        }
    }

  return GSL_SUCCESS;
}


static const gsl_min_fsolver_type bisection_type =
{"bisection",				/* name */
 sizeof (bisection_state_t),
 &bisection_init,
 &bisection_iterate};

const gsl_min_fsolver_type  * gsl_min_fsolver_bisection = &bisection_type;
