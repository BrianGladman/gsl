/* conjugate.c -- Conjugate gradients */

#include <gsl_multimin.h>
#include <gsl_blas_types.h>
#include <gsl_blas.h>
#include <gsl_machine.h>

typedef struct
  {
    int first_move;
    gsl_vector *last_direction;
    int (* gamma)(const gsl_multimin_fdf_history *h,double *value);
    /*    gsl_vector *delta_gradient;*/
  }
conjugate_state_t;

int
gamma_polak_ribiere(const gsl_multimin_fdf_history *h,double *value)
{
  double scalar;
  double gamma;
  size_t i;
  double tmp;
  
  gsl_blas_ddot(h->g1,h->g1,&scalar);
  /* FIXME: good idea?*/
  if (scalar < (GSL_DBL_EPSILON*GSL_DBL_EPSILON)) 
    {
      return GSL_FAILURE;
    }
  else 
    {
      gamma = 0;
      for(i = 0; i<h->g->size; i++) 
	{
	  tmp = gsl_vector_get(h->g,i);
	  gamma += tmp * (tmp - gsl_vector_get(h->g1,i));
	}
      gamma /= scalar;
      *value = gamma;
      return GSL_SUCCESS;
    }
}

int
gamma_fletcher_reeves(const gsl_multimin_fdf_history *h,double *value)
{
  double scalar;
  double gamma;
  
  gsl_blas_ddot(h->g1,h->g1,&scalar);
  /* FIXME: good idea?*/
  if (scalar < (GSL_DBL_EPSILON*GSL_DBL_EPSILON)) 
    {
      return GSL_FAILURE;
    }
  else 
    {
      gsl_blas_ddot(h->g,h->g,&gamma);
      gamma /= scalar;
      *value = gamma;
      return GSL_SUCCESS;
    }
}

int
conjugate_alloc(void *vstate, size_t n)
{
  conjugate_state_t *state = (conjugate_state_t *)vstate;

  state->first_move = 1;
  state->last_direction = gsl_vector_calloc(n);
  if (state->last_direction == 0) 
    {
      GSL_ERROR_RETURN ("failed to allocate PR conjugate gradient internal struct",
			GSL_ENOMEM, 0);     
    }
  /*      state->delta_gradient = gsl_vector_calloc(n);
	  if (state->delta_gradient == 0) 
	  {
	  gsl_vector_free(state->last_direction);
	  GSL_ERROR_RETURN ("failed to allocate PR conjugate gradient internal struct",
	  GSL_ENOMEM, 0);     
	  }*/
  return GSL_SUCCESS;
}

int
conjugate_alloc_pr(void *vstate, size_t n)
{
  conjugate_state_t *state = (conjugate_state_t *)vstate;

  state->gamma = gamma_polak_ribiere;
  return conjugate_alloc(state,n);
}

int
conjugate_alloc_fr(void *vstate, size_t n)
{
  conjugate_state_t *state = (conjugate_state_t *)vstate;

  state->gamma = gamma_fletcher_reeves;
  return conjugate_alloc(state,n);
}

int
conjugate_restart(void *vstate)
{
  conjugate_state_t *state = (conjugate_state_t *)vstate;

  state->first_move = 1;
  return GSL_SUCCESS;
}

void
conjugate_free(void *vstate)
{
  conjugate_state_t *state = (conjugate_state_t *)vstate;

  gsl_vector_free(state->last_direction);
  /*gsl_vector_free(state->delta_gradient);*/
}

int 
conjugate_direction(void *vstate,gsl_multimin_fdf_history *h ,gsl_vector * dir) 
{
  conjugate_state_t *state = (conjugate_state_t *)vstate;
  double gamma;
  double tmp;
  int status;
  size_t i;

  if (state->first_move) 
    {
      for(i = 0; i<h->g->size; i++) 
	{
	  tmp = -gsl_vector_get(h->g,i);
	  gsl_vector_set(state->last_direction,i,tmp);
	  gsl_vector_set(dir,i,tmp);
	}
      /*      gsl_vector_copy(state->last_direction,h->g);
      gsl_blas_dscal(-1,state->last_direction);
      gsl_vector_copy(dir,state->last_direction);*/
      state->first_move = 0;
      return GSL_SUCCESS;
    }
  else
    {
      status = (state->gamma)(h,&gamma);
      if (status == GSL_FAILURE) 
	{
	  /* auto-restarting, the gradient is to small */
	  for(i = 0; i<h->g->size; i++) 
	    {
	      gsl_vector_set(state->last_direction,i,-gsl_vector_get(h->g,i));
	    }
	}
      else 
	{
	  gsl_blas_dscal(gamma,state->last_direction); 
	  for(i = 0; i<h->g->size; i++) 
	    {
	      tmp = gsl_vector_get(state->last_direction,i)-gsl_vector_get(h->g,i);
	      gsl_vector_set(state->last_direction,i,tmp);
	    }
	  /*	  gsl_vector_copy(state->delta_gradient,h->g);
	  gsl_blas_daxpy(-1,h->g1,state->delta_gradient);
	  gsl_blas_ddot(h->g,state->delta_gradient,&gamma);
	  gamma /= scalar;
	  gsl_blas_dscal(gamma,state->last_direction); 
	  gsl_blas_daxpy(-1,h->g,state->last_direction); */
	}
      gsl_vector_copy(dir,state->last_direction);
      return GSL_SUCCESS;
    }
}

static const gsl_multimin_fdf_minimizer_type conjugate_type_pr =
{"conjugate_pr",			/* name */
 sizeof (conjugate_state_t),
 &conjugate_alloc_pr,
 &conjugate_restart,
 &conjugate_direction,
 &conjugate_free};

const gsl_multimin_fdf_minimizer_type *gsl_multimin_fdf_minimizer_conjugate_pr = &conjugate_type_pr;

static const gsl_multimin_fdf_minimizer_type conjugate_type_fr =
{"conjugate_fr",			/* name */
 sizeof (conjugate_state_t),
 &conjugate_alloc_fr,
 &conjugate_restart,
 &conjugate_direction,
 &conjugate_free};

const gsl_multimin_fdf_minimizer_type *gsl_multimin_fdf_minimizer_conjugate_fr = &conjugate_type_fr;
