/* abm.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * Copyright (C) 2004 Daniel Rodríguez
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

/* 4 step Adams-Bashford-Moulton */

/* Authors: Daniel Rodríguez
 */

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "odeiv_util.h"

#define FALSE 0
#define TRUE !FALSE

typedef struct
{
  double valueT;
  double * valueY;
  double * valueF;
} node_t;

typedef struct
{
  node_t * TY;
  int initialized;
} abm_state_t;

static void *
abm_alloc (size_t dim)
{
  abm_state_t *state = (abm_state_t *) malloc (sizeof (abm_state_t));

  if (state == NULL) {
    GSL_ERROR_NULL ("failed to allocate space for abm_state", GSL_ENOMEM);
  }
  state->initialized = 3;

  state->TY = (node_t*) malloc (4 * sizeof (node_t));
  if (state->TY == NULL) {
    free(state);
    GSL_ERROR_NULL("failed to allocate space for valueY", GSL_ENOMEM);
  }
  state->TY[0].valueY = (double*) malloc (dim * sizeof (double));
  if (state->TY[0].valueY == NULL) {
    free(state);
    GSL_ERROR_NULL("failed to allocate space for valueY", GSL_ENOMEM);
  }
  state->TY[1].valueY = (double*) malloc (dim * sizeof (double));
  if (state->TY[1].valueY == NULL) {
    free(state->TY[0].valueY);
    free(state);
    GSL_ERROR_NULL("failed to allocate space for valueY", GSL_ENOMEM);
  }
  state->TY[2].valueY = (double*) malloc (dim * sizeof (double));
  if (state->TY[2].valueY == NULL) {
    free(state->TY[1].valueY);
    free(state->TY[0].valueY);
    free(state);
    GSL_ERROR_NULL("failed to allocate space for valueY", GSL_ENOMEM);
  }
  state->TY[3].valueY = (double*) malloc (dim * sizeof (double));
  if (state->TY[3].valueY == NULL) {
    free(state->TY[2].valueY);
    free(state->TY[1].valueY);
    free(state->TY[0].valueY);
    free(state);
    GSL_ERROR_NULL("failed to allocate space for valueY", GSL_ENOMEM);
  }

  state->TY[0].valueF = (double*) malloc (dim * sizeof (double));
  if (state->TY[0].valueF == NULL) {
    free(state->TY[3].valueY);
    free(state->TY[2].valueY);
    free(state->TY[1].valueY);
    free(state->TY[0].valueY);
    free(state);
    GSL_ERROR_NULL("failed to allocate space for valueY", GSL_ENOMEM);
  }
  state->TY[1].valueF = (double*) malloc (dim * sizeof (double));
  if (state->TY[1].valueF == NULL) {
    free(state->TY[0].valueF);
    free(state->TY[3].valueY);
    free(state->TY[2].valueY);
    free(state->TY[1].valueY);
    free(state->TY[0].valueY);
    free(state);
    GSL_ERROR_NULL("failed to allocate space for valueY", GSL_ENOMEM);
  }
  state->TY[2].valueF = (double*) malloc (dim * sizeof (double));
  if (state->TY[2].valueF == NULL) {
    free(state->TY[1].valueF);
    free(state->TY[0].valueF);
    free(state->TY[3].valueY);
    free(state->TY[2].valueY);
    free(state->TY[1].valueY);
    free(state->TY[0].valueY);
    free(state);
    GSL_ERROR_NULL("failed to allocate space for valueY", GSL_ENOMEM);
  }
  state->TY[3].valueF = (double*) malloc (dim * sizeof (double));
  if (state->TY[3].valueF == NULL) {
    free(state->TY[2].valueF);
    free(state->TY[1].valueF);
    free(state->TY[0].valueF);
    free(state->TY[3].valueY);
    free(state->TY[2].valueY);
    free(state->TY[1].valueY);
    free(state->TY[0].valueY);
    free(state);
    GSL_ERROR_NULL("failed to allocate space for valueY", GSL_ENOMEM);
  }

  return state;
}


static int
abm_apply (void *vstate,
           size_t dim,
           double t,
           double h,
           double y[],
           double yerr[],
           const double dydt_in[],
           double dydt_out[], 
           const gsl_odeiv_system * sys)
{
  abm_state_t *state = (abm_state_t *) vstate;

  int * initialized = &(state->initialized);

  double temp, proxt;
  size_t i;
  int status = 0;
  int s;

  double * p = (double*) malloc(dim * sizeof(double));
  double * Ftmp = (double*) malloc(dim * sizeof(double));
  node_t *TY = state->TY;

  if (*initialized) { /* obtain the first 3 steps with RK4 */
    TY[3-*initialized].valueT = t;
    DBL_MEMCPY (TY[3-*initialized].valueY, y, dim);
    s = GSL_ODEIV_FN_EVAL
      (sys, t, TY[3-*initialized].valueY, TY[3-*initialized].valueF);
    GSL_STATUS_UPDATE (&status, s);
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rk4;
    gsl_odeiv_step * step = gsl_odeiv_step_alloc(T, dim);
    gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc(dim);
    int s2 = gsl_odeiv_evolve_apply(e, NULL, step, sys, &t, t+4*h, &h, y);
    GSL_STATUS_UPDATE (&status, s2);
    DBL_MEMCPY(TY[3-*initialized+1].valueY, y, dim);
    TY[3-*initialized+1].valueT = t;
    s = GSL_ODEIV_FN_EVAL
      (sys, t, TY[3-*initialized+1].valueY, TY[3-*initialized+1].valueF);
    GSL_STATUS_UPDATE (&status, s);
    for (i=0; i<dim; i++) {
      yerr[i] = h * TY[3-*initialized+1].valueF[i];
      if (dydt_out != NULL)
	dydt_out[i] = TY[3-*initialized+1].valueF[i];
    }
    (*initialized)--;
  } else {
    for (i=0; i<dim; i++) {
      p[i] =
	TY[3].valueY[i] + (h/24)*(-9*TY[0].valueF[i]+37*TY[1].valueF[i]-
				  59*TY[2].valueF[i]+55*TY[3].valueF[i]);
    }
    proxt = TY[3].valueT + h;
    GSL_ODEIV_FN_EVAL (sys, proxt, p, Ftmp);
    GSL_STATUS_UPDATE (&status, s);
    for (i=0; i<dim; i++) {
      temp = (h/24)*(TY[1].valueF[i]-5*TY[2].valueF[i]+
		     19*TY[3].valueF[i]+9*Ftmp[i]);
      p[i] = TY[3].valueY[i] + temp;
      yerr[i] = h * temp;
      if (dydt_out != NULL)
	dydt_out[i] = temp;
    }
    DBL_MEMCPY(TY[0].valueY, TY[1].valueY, dim);
    DBL_MEMCPY(TY[1].valueY, TY[2].valueY, dim);
    DBL_MEMCPY(TY[2].valueY, TY[3].valueY, dim);
    DBL_MEMCPY(TY[3].valueY, p, dim);
    TY[0].valueT = TY[1].valueT;
    TY[1].valueT = TY[2].valueT;
    TY[2].valueT = TY[3].valueT;
    TY[3].valueT = proxt;
    DBL_MEMCPY(TY[0].valueF, TY[1].valueF, dim);
    DBL_MEMCPY(TY[1].valueF, TY[2].valueF, dim);
    DBL_MEMCPY(TY[2].valueF, TY[3].valueF, dim);
    GSL_ODEIV_FN_EVAL (sys, proxt, TY[3].valueY, Ftmp);
    GSL_STATUS_UPDATE (&status, s);
    DBL_MEMCPY(TY[3].valueF, Ftmp, dim);
    
    DBL_MEMCPY(y, TY[3].valueY, dim);
  }
  return status;
}

static int
abm_reset (void *vstate, size_t dim)
{
  abm_state_t *state = (abm_state_t *) vstate;

  DBL_ZERO_MEMSET (state->TY[0].valueY, dim);
  DBL_ZERO_MEMSET (state->TY[1].valueY, dim);
  DBL_ZERO_MEMSET (state->TY[2].valueY, dim);
  DBL_ZERO_MEMSET (state->TY[3].valueY, dim);
  DBL_ZERO_MEMSET (state->TY[0].valueF, dim);
  DBL_ZERO_MEMSET (state->TY[1].valueF, dim);
  DBL_ZERO_MEMSET (state->TY[2].valueF, dim);
  DBL_ZERO_MEMSET (state->TY[3].valueF, dim);
  state->TY[0].valueT = 0;
  state->TY[1].valueT = 0;
  state->TY[2].valueT = 0;
  state->TY[3].valueT = 0;
  state->initialized = 3;

  return GSL_SUCCESS;
}

static unsigned int
abm_order (void *vstate)
{
  abm_state_t *state = (abm_state_t *) vstate;
  state = 0; /* prevent warnings about unused parameters */
  return 4;
}

static void
abm_free (void *vstate)
{
  abm_state_t *state = (abm_state_t *) vstate;
  free (state->TY[0].valueY);
  free (state->TY[1].valueY);
  free (state->TY[2].valueY);
  free (state->TY[3].valueY);
  free (state->TY[0].valueF);
  free (state->TY[1].valueF);
  free (state->TY[2].valueF);
  free (state->TY[3].valueF);
  free (state->TY);
  free (state);
}

static const gsl_odeiv_step_type abm_type = { "abm",    /* name */
  0,                            /* can use dydt_in */
  0,                            /* gives exact dydt_out */
  &abm_alloc,
  &abm_apply,
  &abm_reset,
  &abm_order,
  &abm_free
};

const gsl_odeiv_step_type *gsl_odeiv_step_abm = &abm_type;

