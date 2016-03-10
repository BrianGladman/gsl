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
 *
 * [2] J. E. Dennis and H. H. W. Mei, Two new unconstrained optimization
 *     algorithms which use function and gradient values, J. Opt. Theory and
 *     Appl., 28(4), 1979.
 */

typedef struct
{
  size_t n;                  /* number of observations */
  size_t p;                  /* number of parameters */
  gsl_vector *dx_gn;         /* Gauss-Newton step, size p */
  gsl_vector *dx_sd;         /* steepest descent step, size p */
  double norm_gn;            /* || dx_gn || */
  double norm_sd;            /* || dx_sd || */
  double norm_g;             /* || g || */
  double norm_Jg;            /* || J g || */
  gsl_vector *workp;         /* workspace, length p */
  gsl_vector *workn;         /* workspace, length n */

  void *update_state;        /* workspace for parameter update method */
  void *solver_state;        /* workspace for linear solver */

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
static int dogleg_double_step(const void * vtrust_state, const double delta,
                              gsl_vector * dx, void * vstate);
static int dogleg_preduction(const void * vtrust_state, const gsl_vector * dx,
                             double * pred, void * vstate);
static double dogleg_beta(const double t, const double delta, dogleg_state_t *state);

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
  (void)vtrust_state;
  (void)vstate;

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
3) state->norm_g contains || g ||
4) state->norm_Jg contains || J g ||
*/

static int
dogleg_preloop(const void * vtrust_state, void * vstate)
{
  int status;
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  dogleg_state_t *state = (dogleg_state_t *) vstate;
  const gsl_multifit_nlinear_parameters *params = trust_state->params;
  double u;
  double alpha; /* ||g||^2 / ||Jg||^2 */

  /* initialize linear least squares solver */
  status = (params->solver->init)(trust_state->J, state->solver_state);
  if (status)
    return status;

  /* prepare the linear solver to compute Gauss-Newton step */
  status = (params->solver->presolve)(0.0, state->solver_state);
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
  state->norm_g = gsl_blas_dnrm2(trust_state->g);
  state->norm_Jg = gsl_blas_dnrm2(state->workn);

  /* alpha = |g|^2 / |Jg|^2 */
  u = state->norm_g / state->norm_Jg;
  alpha = u * u;

  /* dx_sd = -alpha * g */
  gsl_vector_memcpy(state->dx_sd, trust_state->g);
  gsl_vector_scale(state->dx_sd, -alpha);

  /* store norms */
  state->norm_gn = gsl_blas_dnrm2(state->dx_gn);
  state->norm_sd = gsl_blas_dnrm2(state->dx_sd);

  return GSL_SUCCESS;
}

/*
dogleg_step()
  Calculate a new step vector
*/

static int
dogleg_step(const void * vtrust_state, const double delta,
            gsl_vector * dx, void * vstate)
{
  dogleg_state_t *state = (dogleg_state_t *) vstate;

  if (state->norm_gn > delta)
    {
      if (state->norm_sd >= delta)
        {
          /* both Gauss-Newton and steepest descent steps are outside
           * trust region; truncate steepest descent step to trust
           * region boundary */
          gsl_vector_memcpy(dx, state->dx_sd);
          gsl_vector_scale(dx, delta / state->norm_sd);
        }
      else
        {
          /* Gauss-Newton step is outside trust region, but steepest
           * descent is inside; use dogleg step */

          double beta = dogleg_beta(1.0, delta, state);

          /* dx = dx_sd + beta*(dx_gn - dx_sd) */
          scaled_addition(beta, state->workp, 1.0, state->dx_sd, dx);
        }
    }
  else
    {
      /* Gauss-Newton step is inside trust region, use it as final step */
      gsl_vector_memcpy(dx, state->dx_gn);
    }

  (void)vtrust_state;

  return GSL_SUCCESS;
}

/*
dogleg_double_step()
  Calculate a new step with double dogleg method. Based on
section 3 of [2]
*/

static int
dogleg_double_step(const void * vtrust_state, const double delta,
                   gsl_vector * dx, void * vstate)
{
  const double alpha_fac = 0.8; /* recommended value from Dennis and Mei */
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  dogleg_state_t *state = (dogleg_state_t *) vstate;

  if (state->norm_gn > delta)
    {
      double t, u, v, c;

      /* compute: u = ||g||^2 / ||J*g||^2 */
      v = state->norm_g / state->norm_Jg;
      u = v * v;

      /* compute: v = g^T dx_gn */
      gsl_blas_ddot(trust_state->g, state->dx_gn, &v);

      /* compute: c = ||g||^4 / (||J*g||^2 * |g^T dx_gn|) */
      c = u * (state->norm_g / fabs(v)) * state->norm_g;

      /* compute: t = 1 - alpha_fac*(1-c) */
      t = 1.0 - alpha_fac*(1.0 - c);

      if (t * state->norm_gn <= delta)
        {
          /* set dx = (delta / ||dx_gn||) dx_gn */
          gsl_vector_memcpy(dx, state->dx_gn);
          gsl_vector_scale(dx, delta / state->norm_gn);
        }
      else if (state->norm_sd >= delta)
        {
          /* both Gauss-Newton and steepest descent steps are outside
           * trust region; truncate steepest descent step to trust
           * region boundary */
          gsl_vector_memcpy(dx, state->dx_sd);
          gsl_vector_scale(dx, delta / state->norm_sd);
        }
      else
        {
          /* Cauchy point is inside, Gauss-Newton is outside trust region;
           * use double dogleg step */
          double beta = dogleg_beta(t, delta, state);

          /* dx = dx_sd + beta*(t*dx_gn - dx_sd) */
          scaled_addition(beta, state->workp, 1.0, state->dx_sd, dx);
        }
    }
  else
    {
      /* Gauss-Newton step is inside trust region, use it as final step */
      gsl_vector_memcpy(dx, state->dx_gn);
    }

  return GSL_SUCCESS;
}

static int
dogleg_preduction(const void * vtrust_state, const gsl_vector * dx,
                  double * pred, void * vstate)
{
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  dogleg_state_t *state = (dogleg_state_t *) vstate;

  *pred = quadratic_preduction(trust_state->f, trust_state->J, dx, state->workn);

  return GSL_SUCCESS;
}

/*
dogleg_beta()
  This function finds beta in [0,1] such that the step

dx = dx_sd + beta*(t*dx_gn - dx_sd)

has norm

||dx|| = delta

Inputs: t     - amount of Gauss-Newton step to use for dogleg
                (= 1 for classical dogleg, <= 1 for double dogleg)
        delta - trust region radius
        state - workspace

Notes:
1) On output, state->workp contains t*dx_gn - dx_sd
*/

static double
dogleg_beta(const double t, const double delta, dogleg_state_t *state)
{
  double dot1, dot2;
  double dterm;
  double beta;

  /* compute: workp = t*dx_gn - dx_sd */
  scaled_addition(t, state->dx_gn, -1.0, state->dx_sd, state->workp);

  /* compute: dot1 = dx_sd . (t*dx_gn - dx_sd) */
  gsl_blas_ddot(state->dx_sd, state->workp, &dot1);

  /* compute: dot2 = || t*dx_gn - dx_sd ||^2 */
  gsl_blas_ddot(state->workp, state->workp, &dot2);

  /* dterm = delta^2 - || dx_sd ||^2 */
  dterm = (delta + state->norm_sd) * (delta - state->norm_sd);

  if (dot1 > 0.0)
    {
      beta = dterm / (dot1 + sqrt(dot1*dot1 + dot2*dterm));
    }
  else
    {
      beta = (sqrt(dot1*dot1 + dot2*dterm) - dot1) / dot2;
    }

  return beta;
}

static const gsl_multifit_nlinear_trs dogleg_type =
{
  "dogleg",
  dogleg_alloc,
  dogleg_init,
  dogleg_preloop,
  dogleg_step,
  dogleg_preduction,
  dogleg_free
};

const gsl_multifit_nlinear_trs *gsl_multifit_nlinear_trs_dogleg = &dogleg_type;

static const gsl_multifit_nlinear_trs ddogleg_type =
{
  "double-dogleg",
  dogleg_alloc,
  dogleg_init,
  dogleg_preloop,
  dogleg_double_step,
  dogleg_preduction,
  dogleg_free
};

const gsl_multifit_nlinear_trs *gsl_multifit_nlinear_trs_ddogleg = &ddogleg_type;
