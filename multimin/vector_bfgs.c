/* multimin/vector_bfgs.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Fabrice Rossi
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* vector_bfgs.c -- Limited memory Broyden-Fletcher-Goldfarb-Shanno conjugate gradient method */

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas_types.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_machine.h>

typedef struct
  {
    int first_move;
    gsl_vector *p;
    gsl_vector *v;
  }
vector_bfgs_state_t;

int 
vector_bfgs_direction(void *vstate,gsl_multimin_fdf_history *h ,
		      gsl_vector * dir) 
{
  vector_bfgs_state_t *state = (vector_bfgs_state_t *)vstate;
  double tmp;
  size_t i;
  double p_dot_v,inv_p_dot_v;
  double p_dot_g;
  double v_dot_g;
  double v_dot_v;
  double B,A;
  if (state->first_move) 
    {
      for(i = 0; i<h->g->size; i++) 
	{
	  tmp = -gsl_vector_get(h->g,i);
	  gsl_vector_set(dir,i,tmp);
	}
      state->first_move = 0;
      return GSL_SUCCESS;
    }
  else
    {
      for(i = 0; i<h->x->size; i++) 
	{
	  gsl_vector_set(state->p,i,
			 gsl_vector_get(h->x,i) - gsl_vector_get(h->x1,i));
	  tmp = gsl_vector_get(h->g,i);
	  gsl_vector_set(state->v,i,
			 tmp - gsl_vector_get(h->g1,i));
	  gsl_vector_set(dir,i, -tmp);
	}
      gsl_blas_ddot(state->p,state->v,&p_dot_v);
      /*      printf("delta gradient "); gsl_vector_fprintf (stdout,state->v, "%g");
	      printf("delta position "); gsl_vector_fprintf (stdout,state->p, "%g");*/
      if (fabs(p_dot_v) > (GSL_DBL_EPSILON*GSL_DBL_EPSILON)) 
	{
	  /* the other case induces an auto restart */
	  inv_p_dot_v = 1.0/p_dot_v;
	  gsl_blas_ddot(state->p,h->g,&p_dot_g);
	  B = inv_p_dot_v * p_dot_g;

	  gsl_blas_ddot(state->v,h->g,&v_dot_g);
	  gsl_blas_ddot(state->v,state->v,&v_dot_v);
	  A = - (1 + v_dot_v * inv_p_dot_v) * B + v_dot_g * inv_p_dot_v;

	  /*	  printf ("1/pv=%g, pg=%g, vg=%g, vv=%g, A=%g, B=%g\n",inv_p_dot_v,
		  p_dot_g,v_dot_g,v_dot_v,A,B);*/

	  gsl_blas_daxpy(A,state->p,dir);
	  gsl_blas_daxpy(B,state->v,dir);
	}
      /*      printf("dir "); gsl_vector_fprintf (stdout,dir, "%g");*/
 
      return GSL_SUCCESS;
    }
}

int
vector_bfgs_alloc(void *vstate, size_t n)
{
  vector_bfgs_state_t *state = (vector_bfgs_state_t *)vstate;

  state->first_move = 1;
  state->v = gsl_vector_calloc(n);
  if (state->v == 0) 
    {
      GSL_ERROR_VAL ("failed to allocate PR conjugate gradient internal struct",
			GSL_ENOMEM, 0);     
    }
  state->p = gsl_vector_calloc(n);
  if (state->p == 0) 
    {
      gsl_vector_free(state->v);
      GSL_ERROR_VAL ("failed to allocate PR conjugate gradient internal struct",
			GSL_ENOMEM, 0);     
    }
  return GSL_SUCCESS;
}

int
vector_bfgs_restart(void *vstate)
{
  vector_bfgs_state_t *state = (vector_bfgs_state_t *)vstate;

  state->first_move = 1;
  return GSL_SUCCESS;
}

void
vector_bfgs_free(void *vstate)
{
  vector_bfgs_state_t *state = (vector_bfgs_state_t *)vstate;

  gsl_vector_free(state->v);
  gsl_vector_free(state->p);
}

static const gsl_multimin_fdf_minimizer_type vector_bfgs_type =
{"vector_bfgs",			/* name */
 sizeof (vector_bfgs_state_t),
 &vector_bfgs_alloc,
 &vector_bfgs_restart,
 &vector_bfgs_direction,
 &vector_bfgs_free};

const gsl_multimin_fdf_minimizer_type *gsl_multimin_fdf_minimizer_vector_bfgs = &vector_bfgs_type;

