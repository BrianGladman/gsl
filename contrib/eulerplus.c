/* eulerplus.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * Copyright (C) 2004 Alfonso Acosta & Daniel Rodríguez
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Improved Euler Method */

/* Authors:  Alfonso Acosta & Daniel Rodríguez
 */

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "odeiv_util.h"

typedef struct
{
  double *k1;
  double *k2;
  double *ytmp;
}
eulerplus_state_t;

static void *
eulerplus_alloc (size_t dim)
{
  eulerplus_state_t *state =
    (eulerplus_state_t *) malloc (sizeof (eulerplus_state_t));

  if (state == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for eulerplus_state",
		      GSL_ENOMEM);
    }

  state->k1 = (double *) malloc (dim * sizeof (double));

  if (state->k1 == 0)
    {
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for k1", GSL_ENOMEM);
    }

  state->k2 = (double *) malloc (dim * sizeof (double));

  if (state->k2 == 0)
    {
      free (state->k1);    
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for k2", GSL_ENOMEM);
    }

  
  state->ytmp = (double *) malloc (dim * sizeof (double));

  if (state->ytmp == 0)
    {
      free (state->k2);
      free (state->k1);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for ytmp", GSL_ENOMEM);
    }

  return state;
}


static int
eulerplus_apply (void *vstate,
           size_t dim,
           double t,
           double h,
           double y[],
           double yerr[],
           const double dydt_in[],
           double dydt_out[], 
           const gsl_odeiv_system * sys)
{
  eulerplus_state_t *state = (eulerplus_state_t *) vstate;

  double temp;
  size_t i;
  int status = 0;

  double *const k1 = state->k1;
  double *const k2 = state->k2;
  double *const ytmp = state->ytmp;


  if (dydt_in != NULL)
    {
      DBL_MEMCPY (k1, dydt_in, dim);
    }
  else
    {
      int s = GSL_ODEIV_FN_EVAL (sys, t, y, k1);
      GSL_STATUS_UPDATE (&status, s);
    }
  /* aplicamos euler */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] = y[i] + h * k1[i] ;
    }
  
  {
      int s = GSL_ODEIV_FN_EVAL (sys, t+h, ytmp, k2);
      GSL_STATUS_UPDATE (&status, s);	       
  }

  for (i = 0; i < dim; i++)
   {  
      temp = (k1[i]+k2[i]) / 2; /* la pendiente en x + h */ 
      if (dydt_out != NULL)
          dydt_out[i] = temp;
      yerr[i] = h * temp;
      y[i] += h*temp;
    }


  return status;
}

static int
eulerplus_reset (void *vstate, size_t dim)
{
  eulerplus_state_t *state = (eulerplus_state_t *) vstate;

  DBL_ZERO_MEMSET (state->k1, dim);
  DBL_ZERO_MEMSET (state->k2, dim);
  DBL_ZERO_MEMSET (state->ytmp, dim);

  return GSL_SUCCESS;
}

static unsigned int
eulerplus_order (void *vstate)
{
  eulerplus_state_t *state = (eulerplus_state_t *) vstate;
  state = 0; /* prevent warnings about unused parameters */
  return 2;
}

static void
eulerplus_free (void *vstate)
{
  eulerplus_state_t *state = (eulerplus_state_t *) vstate;
  free (state->k1);
  free (state->k2);
  free (state->ytmp);
  free (state);
}

static const gsl_odeiv_step_type eulerplus_type = { "eulerplus",    /* name */
  1,                            /* can use dydt_in */
  0,                            /* gives exact dydt_out */
  &eulerplus_alloc,
  &eulerplus_apply,
  &eulerplus_reset,
  &eulerplus_order,
  &eulerplus_free
};

const gsl_odeiv_step_type *gsl_odeiv_step_eulerplus = &eulerplus_type;
