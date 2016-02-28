/* multifit_nlinear/lm.c
 * 
 * Copyright (C) 2014, 2015, 2016 Patrick Alken
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
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_eigen.h>

/*
 * This module contains an implementation of the Levenberg-Marquardt
 * algorithm for nonlinear optimization problems. This implementation
 * closely follows the following works:
 *
 * [1] H. B. Nielsen, K. Madsen, Introduction to Optimization and
 *     Data Fitting, Informatics and Mathematical Modeling,
 *     Technical University of Denmark (DTU), 2010.
 *
 * [2] J. J. More, The Levenberg-Marquardt Algorithm: Implementation
 *     and Theory, Lecture Notes in Mathematics, v630, 1978.
 */

typedef struct
{
  size_t n;                  /* number of observations */
  size_t p;                  /* number of parameters */
  gsl_vector *fvv;           /* D_v^2 f(x), size n */
  gsl_vector *vel;           /* geodesic velocity (standard LM step), size p */
  gsl_vector *acc;           /* geodesic acceleration, size p */
  gsl_vector *workp;         /* workspace, length p */
  gsl_vector *workn;         /* workspace, length n */
  double mu;                 /* LM parameter mu */

  void *update_state;        /* workspace for parameter update method */
  void *solver_state;        /* workspace for linear solver */

  gsl_matrix *JTJ;           /* J^T J for rcond calculation */
  gsl_eigen_symm_workspace *eigen_p;

  double avratio;            /* current |a| / |v| */

  /* tunable parameters */
  gsl_multifit_nlinear_parameters params;
} lm_state_t;

#include "common.c"

static void * lm_alloc (const void * params, const size_t n, const size_t p);
static void lm_free(void *vstate);
static int lm_init(const void *vtrust_state, void *vstate);
static int lm_preloop(const void * vtrust_state, void * vstate);
static int lm_step(const void * vtrust_state, const double delta,
                   gsl_vector * dx, void * vstate);
static int lm_check_step(const void * vtrust_state, const gsl_vector * dx,
                         const gsl_vector * f_trial, const double rho, void * vstate);
static int lm_rcond(const gsl_matrix *J, double *rcond, void *vstate);
static double lm_avratio(void *vstate);

static void *
lm_alloc (const void * params, const size_t n, const size_t p)
{
  const gsl_multifit_nlinear_parameters *mparams = (const gsl_multifit_nlinear_parameters *) params;
  lm_state_t *state;
  
  state = calloc(1, sizeof(lm_state_t));
  if (state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate lm state", GSL_ENOMEM);
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

  state->vel = gsl_vector_alloc(p);
  if (state->vel == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for vel", GSL_ENOMEM);
    }

  state->acc = gsl_vector_alloc(p);
  if (state->acc == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for acc", GSL_ENOMEM);
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

  state->JTJ = gsl_matrix_alloc(p, p);
  if (state->JTJ == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for JTJ", GSL_ENOMEM);
    }

  state->eigen_p = gsl_eigen_symm_alloc(p);
  if (state->eigen_p == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for eigen workspace", GSL_ENOMEM);
    }

  state->n = n;
  state->p = p;
  state->params = *mparams;

  return state;
}

static void
lm_free(void *vstate)
{
  lm_state_t *state = (lm_state_t *) vstate;
  const gsl_multifit_nlinear_parameters *params = &(state->params);

  if (state->workp)
    gsl_vector_free(state->workp);

  if (state->workn)
    gsl_vector_free(state->workn);

  if (state->fvv)
    gsl_vector_free(state->fvv);

  if (state->vel)
    gsl_vector_free(state->vel);

  if (state->acc)
    gsl_vector_free(state->acc);

  if (state->update_state)
    (params->update->free)(state->update_state);

  if (state->solver_state)
    (params->solver->free)(state->solver_state);

  if (state->JTJ)
    gsl_matrix_free(state->JTJ);

  if (state->eigen_p)
    gsl_eigen_symm_free(state->eigen_p);

  free(state);
}

/*
lm_init()
  Initialize LM solver

Inputs: vtrust_state - trust state
        vstate       - workspace

Return: success/error
*/

static int
lm_init(const void *vtrust_state, void *vstate)
{
  int status;
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  lm_state_t *state = (lm_state_t *) vstate;
  const gsl_multifit_nlinear_parameters *params = trust_state->params;

  /* initialize LM parameter mu */
  status = (params->update->init)(trust_state->J,
                                  trust_state->x,
                                  &(state->mu),
                                  state->update_state);
  if (status)
    return status;

  gsl_vector_set_zero(state->vel);
  gsl_vector_set_zero(state->acc);

  /* set default parameters */
  state->avratio = 0.0;

  return GSL_SUCCESS;
}

/*
lm_preloop()
  Initialize LM method for new Jacobian matrix
*/

static int
lm_preloop(const void * vtrust_state, void * vstate)
{
  int status;
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  lm_state_t *state = (lm_state_t *) vstate;
  const gsl_multifit_nlinear_parameters *params = trust_state->params;

  /* initialize linear least squares solver */
  status = (params->solver->init)(trust_state->J, state->solver_state);
  if (status)
    return status;

  return GSL_SUCCESS;
}

/*
lm_step()
  Calculate a new step vector by solving the linear
least squares system:

[      J     ] v = - [ f ]
[ sqrt(mu) D ]       [ 0 ]
*/

static int
lm_step(const void * vtrust_state, const double delta,
        gsl_vector * dx, void * vstate)
{
  int status;
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  lm_state_t *state = (lm_state_t *) vstate;
  const gsl_multifit_nlinear_parameters *params = trust_state->params;
  const size_t p = state->p;
  size_t i;

  /* prepare the linear solver with current LM parameter mu */
  status = (params->solver->presolve)(state->mu, state->solver_state);
  if (status)
    return status;

  /*
   * solve: [     J      ] v = - [ f ]
   *        [ sqrt(mu)*I ]       [ 0 ]
   */
  status = (params->solver->solve)(trust_state->f,
                                   trust_state->g,
                                   trust_state->J,
                                   state->vel,
                                   state->solver_state);
  if (status)
    return status;

  /* compute unscaled velocity, v := D^{-1} v */
  for (i = 0; i < p; ++i)
    {
      double di = gsl_vector_get(trust_state->diag, i);
      double *vi = gsl_vector_ptr(state->vel, i);
      *vi /= di;
    }

  if (params->accel)
    {
      /* compute geodesic acceleration */
      status = gsl_multifit_nlinear_eval_fvv(params->h_fvv,
                                             trust_state->x,
                                             state->vel,
                                             trust_state->f,
                                             trust_state->J,
                                             trust_state->sqrt_wts,
                                             trust_state->fdf,
                                             state->fvv,
                                             state->workp);
      if (status)
        return status;

      /*
       * solve: [     J      ] a = - [ fvv ]
       *        [ sqrt(mu)*D ]       [  0  ]
       */
      status = (params->solver->solve)(state->fvv,
                                       NULL,
                                       trust_state->J,
                                       state->acc,
                                       state->solver_state);
      if (status)
        return status;
    }

  /* compute unscaled acceleration, a := D^{-1} a */
  for (i = 0; i < p; ++i)
    {
      double di = gsl_vector_get(trust_state->diag, i);
      double *ai = gsl_vector_ptr(state->acc, i);
      *ai /= di;
    }

  /* compute (scaled) step dx_scaled = D*(v + 1/2 a) */
  for (i = 0; i < p; ++i)
    {
      double vi = gsl_vector_get(state->vel, i);
      double ai = gsl_vector_get(state->acc, i);
      double di = gsl_vector_get(trust_state->diag, i);
      double dxi = vi + 0.5 * ai;

      gsl_vector_set(dx, i, dxi * di);
    }

  return GSL_SUCCESS;
}

/*
lm_check_step()
  Test whether a new step should be accepted, and
update mu parameter accordingly

Inputs: vtrust_state - trust state
        dx           - step vector, size p
        f_trial      - proposed residual vector f(x + dx)
        rho          - (output)
        vstate       - workspace

Return:
GSL_SUCCESS to accept step
GSL_FAILURE to reject step

Notes:
1) state->mu is updated according to whether step
is accepted or rejected
*/

static int
lm_check_step(const void * vtrust_state, const gsl_vector * dx,
              const gsl_vector * f_trial, const double rho, void * vstate)
{
  int status;
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  lm_state_t *state = (lm_state_t *) vstate;
  const gsl_multifit_nlinear_parameters *params = trust_state->params;

  if (rho > 0.0)
    status = GSL_SUCCESS;
  else
    status = GSL_FAILURE;

  /* if using geodesic acceleration, check that |a|/|v| < alpha */
  if (params->accel)
    {
      double anorm = gsl_blas_dnrm2(state->acc);
      double vnorm = gsl_blas_dnrm2(state->vel);

      /* store |a| / |v| */
      state->avratio = anorm / vnorm;

      /* reject step if acceleration is too large compared to velocity */
      if (state->avratio > params->avmax)
        status = GSL_FAILURE;
    }

  /* update state->mu */
  if (status == GSL_SUCCESS)
    {
      /* step accepted, decrease mu */
      int s = (params->update->accept)(rho, &(state->mu), state->update_state);
      if (s)
        return s;
    }
  else
    {
      /* step rejected, increase mu */
      int s = (params->update->reject)(&(state->mu), state->update_state);
      if (s)
        return s;
    }

  return status;
}

static int
lm_rcond(const gsl_matrix *J, double *rcond, void *vstate)
{
  int status;
  lm_state_t *state = (lm_state_t *) vstate;
  gsl_vector *eval = state->workp;
  double eval_min, eval_max;

  /* compute J^T J */
  gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, J, 0.0, state->JTJ);

  /* compute eigenvalues of J^T J */
  status = gsl_eigen_symm(state->JTJ, eval, state->eigen_p);
  if (status)
    return status;

  gsl_vector_minmax(eval, &eval_min, &eval_max);

  if (eval_max > 0.0 && eval_min > 0.0)
    {
      *rcond = sqrt(eval_min / eval_max);
    }
  else
    {
      /* compute eigenvalues are not accurate; possibly due
       * to rounding errors in forming J^T J */
      *rcond = 0.0;
    }

  return GSL_SUCCESS;
}

static double
lm_avratio(void *vstate)
{
  lm_state_t *state = (lm_state_t *) vstate;
  return state->avratio;
}

static const gsl_multifit_nlinear_trs lm_type =
{
  "levenberg-marquardt",
  lm_alloc,
  lm_init,
  lm_preloop,
  lm_step,
  lm_check_step,
  lm_rcond,
  lm_free
};

const gsl_multifit_nlinear_trs *gsl_multifit_nlinear_trs_lm = &lm_type;
