/* multifit_nlinear/dogleg.c
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
 * This module contains an implementation of the Powell dogleg
 * algorithm for nonlinear optimization problems. This implementation
 * closely follows the following works:
 *
 * [1] H. B. Nielsen, K. Madsen, Introduction to Optimization and
 *     Data Fitting, Informatics and Mathematical Modeling,
 *     Technical University of Denmark (DTU), 2010.
 */

typedef struct
{
  size_t n;                  /* number of observations */
  size_t p;                  /* number of parameters */
  gsl_vector *dx_gn;         /* Gauss-Newton step, size p */
  gsl_vector *dx_sd;         /* steepest descent step, size p */
  double norm_gn;            /* || dx_gn || */
  double norm_sd;            /* || dx_sd || */
  gsl_vector *fvv;           /* D_v^2 f(x), size n */
  gsl_vector *vel;           /* geodesic velocity (standard step), size p */
  gsl_vector *acc;           /* geodesic acceleration, size p */
  gsl_vector *workp;         /* workspace, length p */
  gsl_vector *workn;         /* workspace, length n */
  double alpha;              /* || g ||^2 / || J*g ||^2 */
  double pred_red;           /* predicted reduction */

  void *update_state;        /* workspace for parameter update method */
  void *solver_state;        /* workspace for linear solver */

  double avratio;            /* current |a| / |v| */

  /* tunable parameters */
  gsl_multifit_nlinear_parameters params;
} dogleg_state_t;

#include "common.c"

static void * dogleg_alloc (const void * params, const size_t n, const size_t p);
static void dogleg_free(void *vstate);
static int dogleg_init(const void *vtrust_state, void *vstate);
static int dogleg_preloop(const void * vtrust_state, void * vstate);
static int dogleg_step(const void * vtrust_state, const double delta,
                       gsl_vector * dx, void * vstate);
static int dogleg_rcond(const gsl_matrix *J, double *rcond, void *vstate);
static double dogleg_avratio(void *vstate);
static double dogleg_calc_rho(const gsl_multifit_nlinear_trust_state * trust_state,
                              const gsl_vector * dx, const gsl_vector * f_trial,
                              dogleg_state_t * state);

static void *
dogleg_alloc (const void * params, const size_t n, const size_t p)
{
  const gsl_multifit_nlinear_parameters *mparams = (const gsl_multifit_nlinear_parameters *) params;
  dogleg_state_t *state;
  
  state = calloc(1, sizeof(dogleg_state_t));
  if (state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate dogleg state", GSL_ENOMEM);
    }

  state->dx_gn = gsl_vector_alloc(p);
  if (state->dx_gn == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for dx_gn", GSL_ENOMEM);
    }

  state->dx_sd = gsl_vector_alloc(p);
  if (state->dx_sd == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for dx_sd", GSL_ENOMEM);
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

  state->n = n;
  state->p = p;
  state->params = *mparams;

  return state;
}

static void
dogleg_free(void *vstate)
{
  dogleg_state_t *state = (dogleg_state_t *) vstate;
  const gsl_multifit_nlinear_parameters *params = &(state->params);

  if (state->dx_gn)
    gsl_vector_free(state->dx_gn);

  if (state->dx_sd)
    gsl_vector_free(state->dx_sd);

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

  free(state);
}

/*
dogleg_init()
  Initialize dogleg solver

Inputs: vtrust_state - trust state
        vstate       - workspace

Return: success/error
*/

static int
dogleg_init(const void *vtrust_state, void *vstate)
{
  int status;
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  dogleg_state_t *state = (dogleg_state_t *) vstate;

  gsl_vector_set_zero(state->vel);
  gsl_vector_set_zero(state->acc);

  /* set default parameters */
  state->avratio = 0.0;

  return GSL_SUCCESS;
}

/*
dogleg_preloop()
  Initialize dogleg method prior to iteration loop.
This involves computing the Gauss-Newton step and
steepest descent step

Notes: on output,
1) state->dx_gn contains Gauss-Newton step
2) state->dx_sd contains steepest descent step
3) state->alpha = ||g||^2 / ||Jg||^2
*/

static int
dogleg_preloop(const void * vtrust_state, void * vstate)
{
  int status;
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  dogleg_state_t *state = (dogleg_state_t *) vstate;
  const gsl_multifit_nlinear_parameters *params = trust_state->params;
  double norm_g,  /* || g || */
         norm_Jg; /* || J*g || */
  double u;

  /* initialize linear least squares solver */
  status = (params->solver->init)(trust_state->J, state->solver_state);
  if (status)
    return status;

  /* prepare the linear solver to compute Gauss-Newton step */
  status = (params->solver->presolve)(0.0, trust_state->diag, state->solver_state);
  if (status)
    return status;

  /* solve: J dx_gn = -f for Gauss-Newton step */
  status = (params->solver->solve)(trust_state->f,
                                   trust_state->g,
                                   trust_state->J,
                                   state->dx_gn,
                                   state->solver_state);
  if (status)
    return status;

  /* now calculate the steepest descent step */

  /* compute: workn = J*g */
  gsl_blas_dgemv(CblasNoTrans, 1.0, trust_state->J, trust_state->g, 0.0, state->workn);

  /* compute |g| and |Jg| */
  norm_g = gsl_blas_dnrm2(trust_state->g);
  norm_Jg = gsl_blas_dnrm2(state->workn);

  /* alpha = |g|^2 / |Jg|^2 */
  u = norm_g / norm_Jg;
  state->alpha = u * u;

  /* dx_sd = -alpha * g */
  gsl_vector_memcpy(state->dx_sd, trust_state->g);
  gsl_vector_scale(state->dx_sd, -state->alpha);

  /* store norms */
  state->norm_gn = gsl_blas_dnrm2(state->dx_gn);
  state->norm_sd = gsl_blas_dnrm2(state->dx_sd);

  return GSL_SUCCESS;
}

/*
dogleg_step()
  Calculate a new step vector

Notes: on output,
1) state->pred_red contains the predicted reduction
*/

static int
dogleg_step(const void * vtrust_state, const double delta,
            gsl_vector * dx, void * vstate)
{
  int status;
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  dogleg_state_t *state = (dogleg_state_t *) vstate;
  const gsl_multifit_nlinear_parameters *params = trust_state->params;
  const double alpha = state->alpha;

  if (state->norm_gn < 0.0)
    {
      /* failed to find Gauss-Newton step due to singular
       * Jacobian; use steepest descent direction */
      gsl_vector_memcpy(dx, state->dx_sd);

      /* truncate step if outside trust region */
      if (state->norm_sd >= delta)
        gsl_vector_scale(dx, delta / state->norm_sd);

      return GSL_SUCCESS;
    }

  if (state->norm_gn > delta)
    {
      if (state->norm_sd >= delta)
        {
          /* both Gauss-Newton and steepest descent steps are outside
           * trust region; truncate steepest descent step to trust
           * region boundary */
          gsl_vector_memcpy(dx, state->dx_sd);
          gsl_vector_scale(dx, delta / state->norm_sd);
          state->pred_red = delta * (2.0*state->norm_sd - delta) / (2.0 * alpha);
        }
      else
        {
          /* Gauss-Newton step is outside trust region, but steepest
           * descent is inside; use dogleg step */

          double normf = gsl_blas_dnrm2(trust_state->f);
          double normg = gsl_blas_dnrm2(trust_state->g);
          double c, beta, dterm, diff_sq;

          /* compute: workp = dx_gn - dx_sd */
          gsl_vector_memcpy(state->workp, state->dx_gn);
          gsl_vector_sub(state->workp, state->dx_sd);

          /* compute: diff_sq = || dx_gn - dx_sd ||^2 */
          gsl_blas_ddot(state->workp, state->workp, &diff_sq);

          /* c = dx_sd . (dx_gn - dx_sd) */
          gsl_blas_ddot(state->dx_sd, state->workp, &c);

          /* dterm = delta^2 - || dx_sd ||^2 */
          dterm = (delta + state->norm_sd) * (delta - state->norm_sd);

          if (c > 0.0)
            {
              beta = dterm / (c + sqrt(c*c + diff_sq*dterm));
            }
          else
            {
              beta = (sqrt(c*c + diff_sq*dterm) - c) / diff_sq;
            }

          /* dx = dx_sd + beta*(dx_gn - dx_sd) */
          gsl_vector_memcpy(dx, state->dx_sd);
          gsl_blas_daxpy(beta, state->workp, dx);

          state->pred_red = 0.5*alpha*(1.0 - beta)*(1.0 - beta)*normg*normg +
                            0.5*beta*(2.0 - beta)*normf*normf;
        }
    }
  else
    {
      double normf = gsl_blas_dnrm2(trust_state->f);

      /* Gauss-Newton step is inside trust region, use it as final step */
      gsl_vector_memcpy(dx, state->dx_gn);
      state->pred_red = 0.5 * normf * normf;
    }

  return GSL_SUCCESS;
}

/*
dogleg_calc_rho()
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
dogleg_calc_rho(const gsl_multifit_nlinear_trust_state * trust_state,
                const gsl_vector * dx, const gsl_vector * f_trial,
                dogleg_state_t * state)
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
dogleg_rcond(const gsl_matrix *J, double *rcond, void *vstate)
{
  int status;
  dogleg_state_t *state = (dogleg_state_t *) vstate;

  *rcond = 0.0; /* XXX */

  return GSL_SUCCESS;
}

static double
dogleg_avratio(void *vstate)
{
  dogleg_state_t *state = (dogleg_state_t *) vstate;
  return state->avratio;
}

static const gsl_multifit_nlinear_method dogleg_type =
{
  "dogleg",
  dogleg_alloc,
  dogleg_init,
  dogleg_preloop,
  dogleg_step,
  NULL,
  dogleg_rcond,
  dogleg_free
};

const gsl_multifit_nlinear_method *gsl_multifit_nlinear_method_dogleg = &dogleg_type;
