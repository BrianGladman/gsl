/* steepest_descent.c -- the dumb steepest descent algorithm */

#include <gsl_multimin.h>
#include <gsl_blas_types.h>
#include <gsl_blas.h>

typedef struct
  {
    double dummy;
  }
steepest_descent_state_t;

int 
steepest_descent_alloc(void *vstate, size_t n)
{
  return GSL_SUCCESS;
}

int 
steepest_descent_restart(void *vstate)
{
  return GSL_SUCCESS;
}

void
steepest_descent_free(void *vstate)
{
  /* nothing */
}

int 
steepest_descent_direction(void *state,gsl_multimin_fdf_history *h ,gsl_vector * dir) 
{
  size_t i;
  
  for(i = 0; i<dir->size; i++) 
    {
      gsl_vector_set(dir,i,-gsl_vector_get(h->g,i));
    }
  return GSL_SUCCESS;
}

static const gsl_multimin_fdf_minimizer_type steepest_descent_type =
{"steepest_descent",			/* name */
 sizeof (steepest_descent_state_t),
 &steepest_descent_alloc,
 &steepest_descent_restart,
 &steepest_descent_direction,
 &steepest_descent_free};

const gsl_multimin_fdf_minimizer_type *gsl_multimin_fdf_minimizer_steepest_descent = &steepest_descent_type;
