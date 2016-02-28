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
  gsl_vector *diag;
  double normg;              /* || D g || */
  double pred_red;           /* predicted reduction */

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
static double cgst_calc_rho(const gsl_multifit_nlinear_trust_state * trust_state,
                            const gsl_vector * dx, const gsl_vector * f_trial,
                            cgst_state_t * state);
static double cgst_calc_tau(const gsl_vector * p, const gsl_vector * d,
                            const gsl_vector * diag, const double delta,
                            cgst_state_t * state);

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

  state->diag = gsl_vector_alloc(p);
  if (state->diag == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for diag", GSL_ENOMEM);
    }
  gsl_vector_set_all(state->diag, 1.0);

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

  if (state->diag)
    gsl_vector_free(state->diag);

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
  int status;
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  cgst_state_t *state = (cgst_state_t *) vstate;

  /* set default parameters */
  state->cgmaxit = 50;
  state->cgtol = 1.0e-6;

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
  int status;
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  cgst_state_t *state = (cgst_state_t *) vstate;
  const gsl_multifit_nlinear_parameters *params = trust_state->params;
  const size_t p = state->p;
  size_t i;

  /* Step 1 of [1], section 2 */

  for (i = 0; i < p; ++i)
    {
      double gi = gsl_vector_get(trust_state->g, i);
      double Di = gsl_vector_get(state->diag, i);

      gsl_vector_set(state->z, i, 0.0);
      gsl_vector_set(state->r, i, -gi);
      gsl_vector_set(state->d, i, -gi / (Di * Di));
    }

  /* compute || D g || */
  state->normg = scaled_norm(state->diag, trust_state->g);

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
  const gsl_multifit_nlinear_parameters *params = trust_state->params;
  const gsl_matrix * J = trust_state->J;
  const gsl_vector * diag = state->diag;
  const size_t p = state->p;
  double alpha, beta, u;
  double norm_Jd, norm_Dinvr, norm_Drp1;
  size_t i, k;

  for (i = 0; i < state->cgmaxit; ++i)
    {
      /* compute || J d_i || */
      gsl_blas_dgemv(CblasNoTrans, 1.0, J, state->d, 0.0, state->workn);
      norm_Jd = gsl_blas_dnrm2(state->workn);

      /* Step 2 of [1], section 2 */
      if (norm_Jd == 0.0)
        {
          double tau = cgst_calc_tau(state->z, state->d, diag, delta, state);

          /* dx = z_i + tau*d_i */
          scaled_addition(1.0, state->z, tau, state->d, dx);

          return GSL_SUCCESS;
        }

      /* Step 3 of [1], section 2 */

      norm_Dinvr = invscaled_norm(diag, state->r);
      u = norm_Dinvr / norm_Jd;
      alpha = u * u;

      /* workp <= z_{i+1} = z_i + alpha_i*d_i */
      scaled_addition(1.0, state->z, alpha, state->d, state->workp);
      u = scaled_norm(diag, state->workp);
      if (u >= delta)
        {
          double tau = cgst_calc_tau(state->z, state->d, diag, delta, state);

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
      norm_Drp1 = scaled_norm(diag, state->r);

      u = norm_Drp1 / state->normg;
      if (u < state->cgtol)
        {
          gsl_vector_memcpy(dx, state->z);
          return GSL_SUCCESS;
        }

      /* Step 5 of [1], section 2 */

      /* compute u = ||D^{-1} r_{i+1}|| / ||D^{-1} r_i|| */
      u = invscaled_norm(diag, state->r) / norm_Dinvr;
      beta = u * u;

      /* compute: d_{i+1} = D^{-2} r_{i+1} + beta*d_i */
      for (k = 0; k < p; ++k)
        {
          double Dk = gsl_vector_get(diag, k);
          double dk = gsl_vector_get(state->d, k);
          double rk = gsl_vector_get(state->r, k);

          gsl_vector_set(state->d, k, rk / (Dk * Dk) + beta*dk);
        }
    }

  /* failed to find CG solution */
  return GSL_EMAXITER;
}

/*
cgst_calc_rho()
  Calculate ratio of actual reduction to predicted
reduction, given by Eq 4.4 of More, 1978.

Inputs: trust_state - trust state
        dx          - proposed step, size p
        f_trial     - f(x + dx)
        state       - workspace

Notes:
1) On input, state->pred_red must contain predicted reduction
*/

static double
cgst_calc_rho(const gsl_multifit_nlinear_trust_state * trust_state,
              const gsl_vector * dx, const gsl_vector * f_trial,
              cgst_state_t * state)
{
  const double normf = gsl_blas_dnrm2(trust_state->f);
  const double normf_trial = gsl_blas_dnrm2(f_trial);
  double rho;
  double actual_reduction;
  double pred_reduction;
  double u, norm_Jp;

  /* if ||f(x+dx)|| > ||f(x)|| reject step immediately */
  if (normf_trial >= normf)
    return -1.0;

  /* compute numerator of rho (actual reduction) */
  u = normf_trial / normf;
  actual_reduction = 1.0 - u*u;

  /* compute denominator of rho (predicted reduction) */

  /* compute J*dx */
  gsl_blas_dgemv(CblasNoTrans, 1.0, trust_state->J, dx, 0.0, state->workn);
  norm_Jp = gsl_blas_dnrm2(state->workn);

  u = norm_Jp / normf;
  pred_reduction = -u * u;

  gsl_blas_ddot(dx, trust_state->g, &u);
  pred_reduction -= 2.0 * u / (normf * normf);

  if (pred_reduction > 0.0)
    rho = actual_reduction / pred_reduction;
  else
    rho = -1.0;

  return rho;
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

|| D (p + tau*d) || = delta
*/

static double
cgst_calc_tau(const gsl_vector * p, const gsl_vector * d,
              const gsl_vector * diag, const double delta,
              cgst_state_t * state)
{
  const size_t N = p->size;
  gsl_vector_view Dp = gsl_vector_subvector(state->workp, 0, N);
  gsl_vector_view Dd = gsl_vector_subvector(state->workn, 0, N);
  double norm_Dp, norm_Dd, u;
  double t1, t2, tau;
  size_t i;

  for (i = 0; i < N; ++i)
    {
      double Di = gsl_vector_get(diag, i);
      double pi = gsl_vector_get(p, i);
      double di = gsl_vector_get(d, i);

      gsl_vector_set(&Dp.vector, i, Di * pi);
      gsl_vector_set(&Dd.vector, i, Di * di);
    }

  norm_Dp = gsl_blas_dnrm2(&Dp.vector);
  norm_Dd = gsl_blas_dnrm2(&Dd.vector);

  /* compute (Dp, Dd) */
  gsl_blas_ddot(&Dp.vector, &Dd.vector, &u);

  t1 = u / (norm_Dd * norm_Dd);
  t2 = t1*u + delta*delta - norm_Dp*norm_Dp;
  tau = -t1 + sqrt(t2) / norm_Dd;

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
