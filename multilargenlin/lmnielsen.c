/* multilargenlin/lmnielsen.c
 * 
 * Copyright (C) 2015 Patrick Alken
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
#include <gsl/gsl_multilarge.h>
#include <gsl/gsl_multilarge_nlin.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>

/*
 * This module contains an implementation of the Levenberg-Marquardt
 * algorithm for nonlinear optimization problems. This implementation
 * closely follows the following work:
 *
 * [1] H. B. Nielsen, K. Madsen, Introduction to Optimization and
 *     Data Fitting, Informatics and Mathematical Modeling,
 *     Technical University of Denmark (DTU), 2010.
 */

typedef struct
{
  size_t p;                  /* number of model parameters */
  gsl_vector *sdiag;         /* sqrt(diag(J^T J)), size p */
  gsl_vector *g;             /* gradient vector, g = J^T f, size p */
  gsl_vector *x_trial;       /* trial parameter vector, size p */
  long nu;                   /* nu */
  double sqrt_mu;            /* square root of LM damping parameter mu */
  double sqrt_tau;           /* initial scale factor for mu */
  double normf;              /* || f(x) || */
  double normf_trial;        /* || f(x + dx) || */
  int eval_J;                /* 1 if we are currently accumulating full (J,f) */
  int init;                  /* 0 if not yet initialized */

  gsl_multilarge_linear_workspace *linear_workspace_p;
} lmn_state_t;

#define GSL_LM_SQRT_ONE_THIRD    (0.577350269189626)

static void *lmn_alloc (const gsl_multilarge_linear_type * T,
                        const size_t p);
static void *lmn_alloc_normal (const size_t p);
static void *lmn_alloc_tsqr (const size_t p);
static void lmn_free(void *vstate);
static int lmn_init(const gsl_vector * x, gsl_multilarge_function_fdf * fdf,
                    void * work, void * vstate);
static int lmn_accumulate(gsl_matrix * J, gsl_vector * f, void * vstate);
static int lmn_iterate(gsl_vector * x, gsl_vector * dx,
                       gsl_multilarge_function_fdf * fdf,
                       void * fdf_work, void * vstate);
static int lmn_eval(const int eval_J, const gsl_vector * x,
                    gsl_multilarge_function_fdf * fdf,
                    void * work, lmn_state_t * state);
static void lmn_trial_step(const gsl_vector * x, const gsl_vector * dx,
                           gsl_vector * x_trial);
static double lmn_calc_dL(const double sqrt_mu, const gsl_vector *dx,
                          const gsl_vector *g);

static void *
lmn_alloc (const gsl_multilarge_linear_type * T, const size_t p)
{
  lmn_state_t *state;

  if (p == 0)
    {
      GSL_ERROR_NULL ("p must be a positive integer",
                      GSL_EINVAL);
    }
  
  state = calloc(1, sizeof(lmn_state_t));
  if (state == 0)
    {
      GSL_ERROR_NULL ("failed to allocate lm state", GSL_ENOMEM);
    }

  state->sdiag = gsl_vector_alloc(p);
  if (state->sdiag == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for sdiag", GSL_ENOMEM);
    }

  state->g = gsl_vector_alloc(p);
  if (state->g == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for g", GSL_ENOMEM);
    }

  state->x_trial = gsl_vector_alloc(p);
  if (state->x_trial == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for x_trial", GSL_ENOMEM);
    }

  state->linear_workspace_p = gsl_multilarge_linear_alloc(T, p);
  if (state->linear_workspace_p == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for multilarge workspace",
                      GSL_ENOMEM);
    }

  state->p = p;
  state->sqrt_tau = sqrt(1.0e-3);
  state->sqrt_mu = 0.0;
  state->nu = 2;
  state->normf = 0.0;
  state->normf_trial = 0.0;
  state->eval_J = 0;
  state->init = 0;

  return state;
}

static void *
lmn_alloc_normal (const size_t p)
{
  return lmn_alloc(gsl_multilarge_linear_normal, p);
}

static void *
lmn_alloc_tsqr (const size_t p)
{
  return lmn_alloc(gsl_multilarge_linear_tsqr, p);
}

static void
lmn_free(void *vstate)
{
  lmn_state_t *state = (lmn_state_t *) vstate;

  if (state->sdiag)
    gsl_vector_free(state->sdiag);

  if (state->g)
    gsl_vector_free(state->g);

  if (state->x_trial)
    gsl_vector_free(state->x_trial);

  if (state->linear_workspace_p)
    gsl_multilarge_linear_free(state->linear_workspace_p);

  free(state);
}

/*
lmn_init()
  Initialize LM solver

Inputs: x      - initial parameter vector
        fdf    - user-supplied callback function
        work   - workspace to provide to fdf
        vstate - local workspace
*/

static int
lmn_init(const gsl_vector * x, gsl_multilarge_function_fdf * fdf,
         void * work, void * vstate)
{
  int status;
  lmn_state_t *state = (lmn_state_t *) vstate;

  /* set init to 0 to compute diag(J^T J) in lmn_accumulate() */
  state->init = 0;
  gsl_vector_set_zero(state->sdiag);

  /* accumulate initial Jacobian and residual vector into system */
  status = lmn_eval(1, x, fdf, work, state);
  if (status)
    return status;

  state->init = 1;

  /* initialize mu = tau * max [ (J^T J)_{ii} ] using square roots */
  state->sqrt_mu = state->sqrt_tau * gsl_vector_max(state->sdiag);

  return GSL_SUCCESS;
}

/*
lmn_accumulate()
  Accumulate block of rows (J,f) into linear least squares
system

Inputs: J      - Jacobian block, n-by-p
        f      - residual vector, n-by-1
        vstate - workspace

Notes:
1) state->normf is updated to track || f ||
*/

static int
lmn_accumulate(gsl_matrix * J, gsl_vector * f, void * vstate)
{
  int status = GSL_SUCCESS;
  lmn_state_t *state = (lmn_state_t *) vstate;

  if (state->eval_J)
    {
      /*
       * we are accumulating full (J,f) system for a new step size
       * calculation
       */

      if (state->init == 0)
        {
          /* initialization step: compute diag(J^T J) for initial mu estimate */
          size_t i;

          /* diag(J^T J) = sum_i diag(J_i^T J_i) */
          for (i = 0; i < state->p; ++i)
            {
              gsl_vector_view c = gsl_matrix_column(J, i);
              double val = gsl_blas_dnrm2(&c.vector);
              double *di = gsl_vector_ptr(state->sdiag, i);
              *di += val;
            }
        }

      /* update g += J^T f */
      gsl_blas_dgemv(CblasTrans, 1.0, J, f, 1.0, state->g);

      /* update normf */
      state->normf = gsl_hypot(state->normf, gsl_blas_dnrm2(f));

      /* scale f <- -f for step size computation */
      gsl_vector_scale(f, -1.0);

      /* accumulate J and -f into large linear system */
      status = gsl_multilarge_linear_accumulate(J, f, state->linear_workspace_p);
    }
  else
    {
      /* we are only accumulating f to test a potential new step size */

      /* update normf_trial */
      state->normf_trial = gsl_hypot(state->normf_trial, gsl_blas_dnrm2(f));
    }

  return status;
}

/*
lmn_iterate()
  Performe one iteration of LM solver with Nielsen
updating criteria for damping parameter.

Inputs: x        - (input/output)
                   on input, current parameter vector
                   on output, new vector x <- x + dx
        dx       - (output) new step size dx, size p
        fdf      - user-supplied (J,f) function
        fdf_work - workspace argument for fdf
        vstate   - workspace
*/

static int
lmn_iterate(gsl_vector * x, gsl_vector * dx,
            gsl_multilarge_function_fdf * fdf,
            void * fdf_work, void * vstate)
{
  int status;
  lmn_state_t *state = (lmn_state_t *) vstate;
  gsl_vector *x_trial = state->x_trial; /* trial x + dx */
  gsl_vector *g = state->g;             /* gradient J^T f */
  int foundstep = 0;                    /* found step dx */
  double dF;                            /* F(x) - F(x + dx) */
  double dL;                            /* L(0) - L(dx) */
  gsl_multilarge_linear_workspace * linear_p = state->linear_workspace_p;
  double rnorm;
  double snorm;

  /* loop until we find an acceptable step dx */
  while (!foundstep)
    {
      /*
       * solve:
       * [      J     ] dx = [ -f ]
       * [ sqrt(mu)*I ]      [  0 ]
       */
      status = gsl_multilarge_linear_solve(state->sqrt_mu, dx, &rnorm, &snorm, linear_p);
      if (status)
        return status;

      /* compute x_trial = x + dx */
      lmn_trial_step(x, dx, x_trial);

      /* accumulate f(x+dx) without Jacobian */
      status = lmn_eval(0, x_trial, fdf, fdf_work, state);
      if (status)
        return status;

      /* compute F(x) - F(x+dx) = 0.5*(||f||^2 - ||f_trial||^2) */
      dF = 0.5 * (state->normf + state->normf_trial) *
                 (state->normf - state->normf_trial);

      /* compute L(0) - L(dx) = 0.5 * dx^T (mu*dx - g) */
      dL = lmn_calc_dL(state->sqrt_mu, dx, g);

      /* check that rho = dF/dL > 0 */
      if ((dL > 0.0) && (dF >= 0.0))
        {
          /* reduction in error, step acceptable */

          double b;

          /* update LM parameter mu */
          b = 2.0 * (dF / dL) - 1.0;
          b = 1.0 - b*b*b;
          if (b > 0.0)
            state->sqrt_mu *= GSL_MAX(GSL_LM_SQRT_ONE_THIRD, sqrt(b));
          else
            state->sqrt_mu *= GSL_LM_SQRT_ONE_THIRD;

          state->nu = 2;

          /* compute new J(x+dx) and f(x+dx) */
          status = lmn_eval(1, x_trial, fdf, fdf_work, state);
          if (status)
            return status;

          /* update x <- x + dx */
          gsl_vector_memcpy(x, x_trial);

          foundstep = 1;
        }
      else
        {
          long nu2;

          /* step did not reduce error, reject step */
          state->sqrt_mu *= sqrt((double) state->nu);
          nu2 = state->nu << 1; /* 2*nu */
          if (nu2 <= state->nu)
            {
              /*
               * nu has wrapped around / overflown, reset mu and nu
               * to original values and break to force another iteration
               */
              state->nu = 2;
              state->sqrt_mu = state->sqrt_tau * gsl_vector_max(state->sdiag);
              break;
            }

          state->nu = nu2;
        }
    }

  return GSL_SUCCESS;
}

static gsl_vector *
lmn_gradient(void * vstate)
{
  lmn_state_t *state = (lmn_state_t *) vstate;
  return state->g;
}

static double
lmn_normf(void * vstate)
{
  lmn_state_t *state = (lmn_state_t *) vstate;
  return state->normf;
}

/*
lmn_eval()
  Evaluate user-supplied Jacobian and residual vector
function. During iteration, it may be necessary to
compute F(x+dx) = 1/2 ||f(x+dx)||^2 without accumulating
the Jacobian. In this case the flag 'eval_J' determines
whether we need the full (J,f) accumulation.

Inputs: eval_J   - if 1, request full (J,f) system from fdf
                   if 0, request only f from fdf
        x        - parameter vector for evaluation of (J,f)
        fdf      - user-supplied function
        fdf_work - workspace argument to fdf
        state    - local workspace
*/

static int
lmn_eval(const int eval_J, const gsl_vector * x,
         gsl_multilarge_function_fdf * fdf,
         void * fdf_work, lmn_state_t * state)
{
  int status;

  if (eval_J)
    {
      /* reset parameters for new (J,f) accumulation */
      gsl_multilarge_linear_reset(state->linear_workspace_p);
      state->normf = 0.0;
      gsl_vector_set_zero(state->g);
    }

  state->normf_trial = 0.0;
  state->eval_J = eval_J;
  
  status = (*(fdf->fdf)) (eval_J, x, fdf->params, fdf_work);

  return status;
}

/* compute x_trial = x + dx */
static void
lmn_trial_step(const gsl_vector * x, const gsl_vector * dx,
               gsl_vector * x_trial)
{
  size_t i, N = x->size;

  for (i = 0; i < N; i++)
    {
      double dxi = gsl_vector_get (dx, i);
      double xi = gsl_vector_get (x, i);
      gsl_vector_set (x_trial, i, xi + dxi);
    }
}

/*
lmn_calc_dL()
  Compute dL = L(0) - L(dx) = 1/2 dx^T (mu * dx - g)
*/

static double
lmn_calc_dL(const double sqrt_mu, const gsl_vector *dx, const gsl_vector *g)
{
  const size_t p = dx->size;
  const double mu = sqrt_mu * sqrt_mu;
  size_t i;
  double dL = 0.0;

  for (i = 0; i < p; ++i)
    {
      double dxi = gsl_vector_get(dx, i);
      double gi = gsl_vector_get(g, i);

      dL += dxi * (mu * dxi - gi);
    }

  dL *= 0.5;

  return dL;
}

static const gsl_multilarge_nlinear_type lmnormal_type =
{
  "lmnormal",
  lmn_alloc_normal,
  lmn_init,
  lmn_accumulate,
  lmn_iterate,
  lmn_gradient,
  lmn_normf,
  lmn_free
};

static const gsl_multilarge_nlinear_type lmtsqr_type =
{
  "lmtsqr",
  lmn_alloc_tsqr,
  lmn_init,
  lmn_accumulate,
  lmn_iterate,
  lmn_gradient,
  lmn_normf,
  lmn_free
};

const gsl_multilarge_nlinear_type *gsl_multilarge_nlinear_lmnormal = &lmnormal_type;
const gsl_multilarge_nlinear_type *gsl_multilarge_nlinear_lmtsqr = &lmtsqr_type;
