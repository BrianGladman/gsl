/* multifit_nlinear/trust.c
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
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_eigen.h>

#include "common.c"

/*
 * This module contains a high level driver for a general trust
 * region nonlinear least squares solver. This container handles
 * the computation of all of the quantities relevant to all trust
 * region methods, including:
 *
 * residual vector: f_k = f(x_k)
 * Jacobian matrix: J_k = J(x_k)
 * gradient vector: g_k = J_k^T f_k
 * scaling matrix:  D_k
 */

typedef struct
{
  size_t n;                  /* number of observations */
  size_t p;                  /* number of parameters */
  gsl_vector *diag;          /* D = diag(J^T J) */
  gsl_vector *x_trial;       /* trial parameter vector */
  gsl_vector *f_trial;       /* trial function vector */
  gsl_vector *workp;         /* workspace, length p */
  gsl_vector *workn;         /* workspace, length n */

  void *method_state;        /* workspace for trust region method */

  double avratio;            /* current |a| / |v| */

  /* tunable parameters */
  gsl_multifit_nlinear_parameters params;
} trust_state_t;

static void * trust_alloc (const gsl_multifit_nlinear_parameters * params,
                           const size_t n, const size_t p);
static void trust_free(void *vstate);
static int trust_init(void *vstate, const gsl_vector * swts,
                      gsl_multifit_nlinear_fdf *fdf, const gsl_vector *x,
                      gsl_vector *f, gsl_matrix *J, gsl_vector *g);
static int trust_iterate(void *vstate, const gsl_vector *swts,
                         gsl_multifit_nlinear_fdf *fdf,
                         gsl_vector *x, gsl_vector *f, gsl_matrix *J,
                         gsl_vector *g, gsl_vector *dx);
static int trust_rcond(const gsl_matrix *J, double *rcond, void *vstate);
static double trust_avratio(void *vstate);

static void *
trust_alloc (const gsl_multifit_nlinear_parameters * params,
             const size_t n, const size_t p)
{
  trust_state_t *state;
  
  state = calloc(1, sizeof(trust_state_t));
  if (state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate lm state", GSL_ENOMEM);
    }

  state->diag = gsl_vector_alloc(p);
  if (state->diag == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for diag", GSL_ENOMEM);
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

  state->x_trial = gsl_vector_alloc(p);
  if (state->x_trial == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for x_trial", GSL_ENOMEM);
    }

  state->f_trial = gsl_vector_alloc(n);
  if (state->f_trial == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for f_trial", GSL_ENOMEM);
    }

  state->method_state = (params->method->alloc)(params, n, p);
  if (state->method_state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for method state", GSL_ENOMEM);
    }

  state->n = n;
  state->p = p;
  state->params = *params;

  return state;
}

static void
trust_free(void *vstate)
{
  trust_state_t *state = (trust_state_t *) vstate;
  const gsl_multifit_nlinear_parameters *params = &(state->params);

  if (state->diag)
    gsl_vector_free(state->diag);

  if (state->workp)
    gsl_vector_free(state->workp);

  if (state->workn)
    gsl_vector_free(state->workn);

  if (state->x_trial)
    gsl_vector_free(state->x_trial);

  if (state->f_trial)
    gsl_vector_free(state->f_trial);

  if (state->method_state)
    (params->method->free)(state->method_state);

  free(state);
}

/*
trust_init()
  Initialize trust region solver

Inputs: vstate - workspace
        swts   - sqrt(W) vector
        fdf    - user callback functions
        x      - initial parameter values
        f      - (output) f(x) vector
        J      - (output) J(x) matrix
        g      - (output) J(x)' f(x) vector

Return: success/error
*/

static int
trust_init(void *vstate, const gsl_vector *swts,
           gsl_multifit_nlinear_fdf *fdf, const gsl_vector *x,
           gsl_vector *f, gsl_matrix *J, gsl_vector *g)
{
  int status;
  trust_state_t *state = (trust_state_t *) vstate;
  const gsl_multifit_nlinear_parameters *params = &(state->params);

  /* evaluate function and Jacobian at x and apply weight transform */
  status = gsl_multifit_nlinear_eval_f(fdf, x, swts, f);
  if (status)
   return status;

  status = gsl_multifit_nlinear_eval_df(x, f, swts, params->h_df,
                                        params->fdtype, fdf, J, state->workn);
  if (status)
    return status;

  /* compute g = J^T f */
  gsl_blas_dgemv(CblasTrans, 1.0, J, f, 0.0, g);

  /* initialize diagonal scaling matrix D */
  (params->scale->init)(J, state->diag);

  /* initialize trust region method solver */
  status = (params->method->init)(x, J, state->diag, state->method_state);
  if (status)
    return status;

  /* set default parameters */
  state->avratio = 0.0;

  return GSL_SUCCESS;
}

/*
trust_iterate()
  This function performs 1 iteration of the trust region algorithm.
It calls a user-specified method for computing the next step
(LM or dogleg), then tests if the computed step is acceptable.

Args: vstate - trust workspace
      swts   - data weights (NULL if unweighted)
      fdf    - function and Jacobian pointers
      x      - on input, current parameter vector
               on output, new parameter vector x + dx
      f      - on input, f(x)
               on output, f(x + dx)
      J      - on input, J(x)
               on output, J(x + dx)
      g      - on input, g(x) = J(x)' f(x)
               on output, g(x + dx) = J(x + dx)' f(x + dx)
      dx     - (output only) parameter step vector

Return:
1) GSL_SUCCESS if we found a step which reduces the cost
function

2) GSL_ENOPROG if 15 successive attempts were to made to
find a good step without success
*/

static int
trust_iterate(void *vstate, const gsl_vector *swts,
              gsl_multifit_nlinear_fdf *fdf, gsl_vector *x,
              gsl_vector *f, gsl_matrix *J, gsl_vector *g,
              gsl_vector *dx)
{
  int status;
  trust_state_t *state = (trust_state_t *) vstate;
  const gsl_multifit_nlinear_parameters *params = &(state->params);
  gsl_vector *x_trial = state->x_trial;       /* trial x + dx */
  gsl_vector *f_trial = state->f_trial;       /* trial f(x + dx) */
  gsl_vector *diag = state->diag;             /* diag(D) */
  double rho;                                 /* ratio dF/dL */
  int foundstep = 0;                          /* found step dx */
  int bad_steps = 0;                          /* consecutive rejected steps */

  /* initialize trust region method with this Jacobian */
  status = (params->method->init_J)(J, state->method_state);
  if (status)
    return status;

  /* loop until we find an acceptable step dx */
  while (!foundstep)
    {
      /* calculate new step */
      status = (params->method->step)(x, f, g, J, diag, swts, fdf,
                                      dx, state->method_state);
      if (status)
        return status;

      /* compute x_trial = x + dx */
      trial_step(x, dx, x_trial);

      /* compute f(x + dx) */
      status = gsl_multifit_nlinear_eval_f(fdf, x_trial, swts, f_trial);
      if (status)
        return status;

      /* check if new step should be accepted */
      status = (params->method->check_step)(f, f_trial, g, diag,
                                            &rho, state->method_state);

      if (status == GSL_SUCCESS)
        {
          /* step was accepted */

          /* compute J <- J(x + dx) */
          status = gsl_multifit_nlinear_eval_df(x_trial, f_trial, swts,
                                                params->h_df, params->fdtype,
                                                fdf, J, state->workn);
          if (status)
            return status;

          /* update x <- x + dx */
          gsl_vector_memcpy(x, x_trial);

          /* update f <- f(x + dx) */
          gsl_vector_memcpy(f, f_trial);

          /* compute new g = J^T f */
          gsl_blas_dgemv(CblasTrans, 1.0, J, f, 0.0, g);

          /* update scaling matrix D */
          (params->scale->update)(J, diag);

          foundstep = 1;
          bad_steps = 0;
        }
      else
        {
          /* if more than 15 consecutive rejected steps, report no progress */
          if (++bad_steps > 15)
            return GSL_ENOPROG;
        }
    }

  return GSL_SUCCESS;
} /* trust_iterate() */

static int
trust_rcond(const gsl_matrix *J, double *rcond, void *vstate)
{
  int status;
  trust_state_t *state = (trust_state_t *) vstate;
  const gsl_multifit_nlinear_parameters *params = &(state->params);

  status = (params->method->rcond) (J, rcond, state->method_state);

  return status;
}

static double
trust_avratio(void *vstate)
{
  trust_state_t *state = (trust_state_t *) vstate;
  return state->avratio;
}

static const gsl_multifit_nlinear_type trust_type =
{
  "trust-region",
  trust_alloc,
  trust_init,
  trust_iterate,
  trust_rcond,
  trust_avratio,
  trust_free
};

const gsl_multifit_nlinear_type *gsl_multifit_nlinear_trust = &trust_type;
