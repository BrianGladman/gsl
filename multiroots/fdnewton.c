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
    double f, df;
  }
newton_state_t;

int newton_init (void * vstate, gsl_function_fdf * fdf, gsl_vector * x);
int newton_iterate (void * vstate, gsl_function_fdf * fdf, gsl_vector * x);

int
newton_init (void * vstate, gsl_function_fdf * fdf, gsl_vector * x)
{
  newton_state_t * state = (newton_state_t *) vstate;

  const double x = *root ;

  state->f = GSL_FN_FDF_EVAL_F (fdf, x);
  state->df = GSL_FN_FDF_EVAL_DF (fdf, x) ;

  return GSL_SUCCESS;
}

int
newton_iterate (void * vstate, gsl_function_fdf * fdf, double * root)
{
  newton_state_t * state = (newton_state_t *) vstate;
  
  double root_new, f_new, df_new;

  if (state->df == 0.0)
    {
      GSL_ERROR("derivative is zero", GSL_EZERODIV);
    }

  root_new = *root - (state->f / state->df);

  *root = root_new ;
  
  GSL_FN_FDF_EVAL_F_DF(fdf, root_new, &f_new, &df_new);

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


static const gsl_root_fdfsolver_type newton_type =
{"newton",				/* name */
 sizeof (newton_state_t),
 &newton_init,
 &newton_iterate};

const gsl_root_fdfsolver_type  * gsl_root_fdfsolver_newton = &newton_type;
