/* multifit_nlinear/cgst.c
 * 
 * Copyright (C) 2016 Patrick Alken
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

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>

/*
 * This module contains an implementation of the Steihaug-Toint
 * conjugate gradient algorithm for nonlinear optimization problems.
 * This implementation closely follows the following works:
 *
 * [1] T. Steihaug, The conjugate gradient method and trust regions
 *     in large scale optimization, SIAM J. Num. Anal., 20(3) 1983.
 */

typedef struct
{
  size_t n;                  /* number of observations */
  size_t p;                  /* number of parameters */
  gsl_vector *z;             /* Gauss-Newton step, size p */
  gsl_vector *r;             /* steepest descent step, size p */
  gsl_vector *d;             /* steepest descent step, size p */
  gsl_vector *fvv;           /* D_v^2 f(x), size n */
  gsl_vector *workp;         /* workspace, length p */
  gsl_vector *workn;         /* workspace, length n */
  double normg;              /* || D g || */

  void *update_state;        /* workspace for parameter update method */
  void *solver_state;        /* workspace for linear solver */

  double cgtol;              /* tolerance for CG solution */
  size_t cgmaxit;            /* maximum CG iterations */

  /* tunable parameters */
  gsl_multifit_nlinear_parameters params;
} cgst_state_t;

#include "common.c"

static void * cgst_alloc (const void * params, const size_t n, const size_t p);
static void cgst_free(void *vstate);
static int cgst_init(const void *vtrust_state, void *vstate);
static int cgst_preloop(const void * vtrust_state, void * vstate);
static int cgst_step(const void * vtrust_state, const double delta,
                     gsl_vector * dx, void * vstate);
static int cgst_rcond(const gsl_matrix *J, double *rcond, void *vstate);
static double cgst_avratio(void *vstate);
static double cgst_calc_tau(const gsl_vector * p, const gsl_vector * d,
                            const double delta, cgst_state_t * state);

static void *
cgst_alloc (const void * params, const size_t n, const size_t p)
{
  const gsl_multifit_nlinear_parameters *mparams = (const gsl_multifit_nlinear_parameters *) params;
  cgst_state_t *state;
  
  state = calloc(1, sizeof(cgst_state_t));
  if (state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate st state", GSL_ENOMEM);
    }

  state->z = gsl_vector_alloc(p);
  if (state->z == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for z", GSL_ENOMEM);
    }

  state->r = gsl_vector_alloc(p);
  if (state->r == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for r", GSL_ENOMEM);
    }

  state->d = gsl_vector_alloc(p);
  if (state->d == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for d", GSL_ENOMEM);
    }

  state->workp = gsl_vector_alloc(p);
  if (state->workp == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for workp", GSL_ENOMEM);
    }

  state->workn = gsl_vector_alloc(n);
  if (state->workn == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for workn", GSL_ENOMEM);
    }

  state->fvv = gsl_vector_alloc(n);
  if (state->fvv == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for fvv", GSL_ENOMEM);
    }

  state->update_state = (mparams->update->alloc)();
  if (state->update_state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for update state", GSL_ENOMEM);
    }

  state->solver_state = (mparams->solver->alloc)(n, p);
  if (state->solver_state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for solver state", GSL_ENOMEM);
    }

  state->n = n;
  state->p = p;
  state->params = *mparams;

  return state;
}

static void
cgst_free(void *vstate)
{
  cgst_state_t *state = (cgst_state_t *) vstate;
  const gsl_multifit_nlinear_parameters *params = &(state->params);

  if (state->z)
    gsl_vector_free(state->z);

  if (state->r)
    gsl_vector_free(state->r);

  if (state->d)
    gsl_vector_free(state->d);

  if (state->workp)
    gsl_vector_free(state->workp);

  if (state->workn)
    gsl_vector_free(state->workn);

  if (state->fvv)
    gsl_vector_free(state->fvv);

  if (state->update_state)
    (params->update->free)(state->update_state);

  if (state->solver_state)
    (params->solver->free)(state->solver_state);

  free(state);
}

/*
cgst_init()
  Initialize st solver

Inputs: vtrust_state - trust state
        vstate       - workspace

Return: success/error
*/

static int
cgst_init(const void *vtrust_state, void *vstate)
{
  cgst_state_t *state = (cgst_state_t *) vstate;

  /* set default parameters */
  state->cgmaxit = 50;
  state->cgtol = 1.0e-6;

  (void)vtrust_state;

  return GSL_SUCCESS;
}

/*
cgst_preloop()
  Initialize st method prior to iteration loop.

Notes: on output,
*/

static int
cgst_preloop(const void * vtrust_state, void * vstate)
{
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  cgst_state_t *state = (cgst_state_t *) vstate;
  const size_t p = state->p;
  size_t i;

  /* Step 1 of [1], section 2 */

  for (i = 0; i < p; ++i)
    {
      double gi = gsl_vector_get(trust_state->g, i);

      gsl_vector_set(state->z, i, 0.0);
      gsl_vector_set(state->r, i, -gi);
      gsl_vector_set(state->d, i, -gi);
    }

  /* compute || g || */
  state->normg = gsl_blas_dnrm2(trust_state->g);

  return GSL_SUCCESS;
}

/*
cgst_step()
  Calculate a new step vector

Return:
GSL_SUCCESS if CG solution found
GSL_EMAXITER if no solution found
*/

static int
cgst_step(const void * vtrust_state, const double delta,
          gsl_vector * dx, void * vstate)
{
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  cgst_state_t *state = (cgst_state_t *) vstate;
  const gsl_matrix * J = trust_state->J;
  double alpha, beta, u;
  double norm_Jd, norm_r, norm_rp1;
  size_t i;

  for (i = 0; i < state->cgmaxit; ++i)
    {
      /* compute || J d_i || */
      gsl_blas_dgemv(CblasNoTrans, 1.0, J, state->d, 0.0, state->workn);
      norm_Jd = gsl_blas_dnrm2(state->workn);

      /* Step 2 of [1], section 2 */
      if (norm_Jd == 0.0)
        {
          double tau = cgst_calc_tau(state->z, state->d, delta, state);

          /* dx = z_i + tau*d_i */
          scaled_addition(1.0, state->z, tau, state->d, dx);

          return GSL_SUCCESS;
        }

      /* Step 3 of [1], section 2 */

      norm_r = gsl_blas_dnrm2(state->r);
      u = norm_r / norm_Jd;
      alpha = u * u;

      /* workp <= z_{i+1} = z_i + alpha_i*d_i */
      scaled_addition(1.0, state->z, alpha, state->d, state->workp);
      u = gsl_blas_dnrm2(state->workp);
      if (u >= delta)
        {
          double tau = cgst_calc_tau(state->z, state->d, delta, state);

          /* dx = z_i + tau*d_i */
          scaled_addition(1.0, state->z, tau, state->d, dx);

          return GSL_SUCCESS;
        }

      /* store z_{i+1} */
      gsl_vector_memcpy(state->z, state->workp);

      /* Step 4 of [1], section 2 */

      /* compute: workp <= alpha B d_i = alpha J^T J d_i,
       * where J d_i is already stored in workn */
      gsl_blas_dgemv(CblasTrans, alpha, J, state->workn, 0.0, state->workp);

      /* r_{i+1} = r_i - alpha*B*d_i */
      gsl_vector_sub(state->r, state->workp);
      norm_rp1 = gsl_blas_dnrm2(state->r);

      u = norm_rp1 / state->normg;
      if (u < state->cgtol)
        {
          gsl_vector_memcpy(dx, state->z);
          return GSL_SUCCESS;
        }

      /* Step 5 of [1], section 2 */

      /* compute u = ||r_{i+1}|| / ||r_i|| */
      u = norm_rp1 / norm_r;
      beta = u * u;

      /* compute: d_{i+1} = r_{i+1} + beta*d_i */
      scaled_addition(1.0, state->r, beta, state->d, state->d);
    }

  /* failed to find CG solution */
  return GSL_EMAXITER;
}

static int
cgst_rcond(const gsl_matrix *J, double *rcond, void *vstate)
{
  int status;
  cgst_state_t *state = (cgst_state_t *) vstate;

  *rcond = 0.0; /* XXX */

  return GSL_SUCCESS;
}

static double
cgst_avratio(void *vstate)
{
  cgst_state_t *state = (cgst_state_t *) vstate;
  return 0.0;
}

/*
cgst_calc_tau()
  Compute tau > 0 such that:

|| p + tau*d || = delta
*/

static double
cgst_calc_tau(const gsl_vector * p, const gsl_vector * d,
              const double delta, cgst_state_t * state)
{
  double norm_p, norm_d, u;
  double t1, t2, tau;

  norm_p = gsl_blas_dnrm2(p);
  norm_d = gsl_blas_dnrm2(d);

  /* compute (p, d) */
  gsl_blas_ddot(p, d, &u);

  t1 = u / (norm_d * norm_d);
  t2 = t1*u + delta*delta - norm_p*norm_p;
  tau = -t1 + sqrt(t2) / norm_d;

  return tau;
}

static const gsl_multifit_nlinear_trs cgst_type =
{
  "steihaug-toint",
  cgst_alloc,
  cgst_init,
  cgst_preloop,
  cgst_step,
  NULL,
  cgst_rcond,
  cgst_free
};

const gsl_multifit_nlinear_trs *gsl_multifit_nlinear_trs_cgst = &cgst_type;
