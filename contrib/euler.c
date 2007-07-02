/* euler.c
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

/* Euler Method*/

/* Authors:  Alfonso Acosta & Daniel Rodríguez
 */

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "odeiv_util.h"

typedef struct
{
  double *k;
}
euler_state_t;

static void *
euler_alloc (size_t dim)
{
  euler_state_t *state = (euler_state_t *) malloc (sizeof (euler_state_t));

  if (state == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for euler_state", GSL_ENOMEM);
    }

  state->k = (double *) malloc (dim * sizeof (double));

  if (state->k == 0)
    {
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for k", GSL_ENOMEM);
    }

  return state;
}


static int
euler_apply (void *vstate,
           size_t dim,
           double t,
           double h,
           double y[],
           double yerr[],
           const double dydt_in[],
           double dydt_out[], 
           const gsl_odeiv_system * sys)
{
  euler_state_t *state = (euler_state_t *) vstate;
  
  double temp;
  size_t i;
  int status = 0;

  double *const k = state->k;


  if (dydt_in != NULL)
    {
      DBL_MEMCPY (k, dydt_in, dim);
    }
  else
    {
      int s = GSL_ODEIV_FN_EVAL (sys, t, y, k);
      GSL_STATUS_UPDATE (&status, s);
    }
  

  for (i = 0; i < dim; i++)
    {
      temp = y[i]; /* save y[i] */ 
      y[i] = h * k[i];   
      yerr[i] = h * y[i];
      y[i] += temp;
      if (dydt_out != NULL)
          dydt_out[i] = k[i];
    }

  return status;
}

static int
euler_reset (void *vstate, size_t dim)
{
  euler_state_t *state = (euler_state_t *) vstate;

  DBL_ZERO_MEMSET (state->k, dim);

  return GSL_SUCCESS;
}

static unsigned int
euler_order (void *vstate)
{
  euler_state_t *state = (euler_state_t *) vstate;
  state = 0; /* prevent warnings about unused parameters */
  return 1;
}

static void
euler_free (void *vstate)
{
  euler_state_t *state = (euler_state_t *) vstate;
  free (state->k);
  free (state);
}

static const gsl_odeiv_step_type euler_type = { "euler",    /* name */
  1,                            /* can use dydt_in */
  0,                            /* gives exact dydt_out */
  &euler_alloc,
  &euler_apply,
  &euler_reset,
  &euler_order,
  &euler_free
};

const gsl_odeiv_step_type *gsl_odeiv_step_euler = &euler_type;
