/* secant.c -- secant root finding algorithm 

   The secant algorithm is a variant of the Newton algorithm with the
   derivative term replaced by a numerical estimate from the last two
   function evaluations.

   x[i+1] = x[i] - f(x[i]) / f'_est

   where f'_est = (f(x[i]) - f(x[i-1])) / (x[i] - x[i-1])

   The exact derivative is used for the initial value of f'_est.

*/

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
    double f;
    double df;
  }
secant_state_t;

static int secant_init (void * vstate, gsl_function_fdf * fdf, double * root);
static int secant_iterate (void * vstate, gsl_function_fdf * fdf, double * root);

static int
secant_init (void * vstate, gsl_function_fdf * fdf, double * root)
{
  secant_state_t * state = (secant_state_t *) vstate;

  const double x = *root;

  GSL_FN_FDF_EVAL_F_DF (fdf, x, &(state->f), &(state->df));
  
  return GSL_SUCCESS;

}

static int
secant_iterate (void * vstate, gsl_function_fdf * fdf, double * root)
{
  secant_state_t * state = (secant_state_t *) vstate;
  
  const double x = *root ;
  const double f = state->f;
  const double df = state->df;

  double x_new, f_new, df_new;

  if (state->df == 0.0)
    {
      GSL_ERROR("derivative is zero", GSL_EZERODIV);
    }

  x_new = x - (f / df);

  f_new = GSL_FN_FDF_EVAL_F(fdf, x_new) ;
  df_new = (f_new - f) / (x_new - x) ;

  *root = x_new ;

  state->f = f_new ;
  state->df = df_new ;

  if (!GSL_IS_REAL (f_new))
    {
      GSL_ERROR ("function not continuous", GSL_EBADFUNC);
    }

  if (!GSL_IS_REAL (df_new))
    {
      GSL_ERROR ("function not differentiable", GSL_EBADFUNC);
    }
      
  return GSL_SUCCESS;
}


static const gsl_root_fdfsolver_type secant_type =
{"secant",				/* name */
 sizeof (secant_state_t),
 &secant_init,
 &secant_iterate};

const gsl_root_fdfsolver_type  * gsl_root_fdfsolver_secant = &secant_type;
