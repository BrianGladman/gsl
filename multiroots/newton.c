#include <config.h>

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <gsl_math.h>
#include <gsl_errno.h>
#include <gsl_multiroots.h>

typedef struct
  {
    double f, df;
  }
newton_state_t;

int newton_init (void * vstate, gsl_multiroot_function_fdf * fdf, gsl_vector * x, gsl_vector * f, gsl_matrix * J, gsl_vector * dx);
int newton_iterate (void * vstate, gsl_multiroot_function_fdf * fdf, gsl_vector * x, gsl_vector * f, gsl_matrix * J, gsl_vector * dx);

int
newton_init (void * vstate, gsl_multiroot_function_fdf * fdf, gsl_vector * x, gsl_vector * f, gsl_matrix * J, gsl_vector * dx)
{
  newton_state_t * state = (newton_state_t *) vstate;

  return GSL_SUCCESS;
}

int
newton_iterate (void * vstate, gsl_multiroot_function_fdf * fdf, gsl_vector * x, gsl_vector * f, gsl_matrix * J, gsl_vector * dx)
{
  newton_state_t * state = (newton_state_t *) vstate;
  
      
  return GSL_SUCCESS;
}


static const gsl_multiroot_fdfsolver_type newton_type =
{"newton",				/* name */
 sizeof (newton_state_t),
 &newton_init,
 &newton_iterate};

const gsl_multiroot_fdfsolver_type  * gsl_multiroot_fdfsolver_newton = &newton_type;
