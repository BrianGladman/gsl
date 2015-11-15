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

#include "oct.h"

/*
 * This module contains an implementation of the Levenberg-Marquardt
 * algorithm for nonlinear optimization problems. This implementation
 * closely follows the following works:
 *
 * [1] H. B. Nielsen, K. Madsen, Introduction to Optimization and
 *     Data Fitting, Informatics and Mathematical Modeling,
 *     Technical University of Denmark (DTU), 2010.
 */

typedef struct
{
  gsl_matrix *A;             /* J^T J */
  gsl_matrix *A_copy;        /* copy of J^T J */
  gsl_matrix *J;             /* Jacobian J(x) */
  gsl_vector *diag;          /* D = diag(J^T J) */
  gsl_vector *rhs;           /* rhs vector = -g = -J^T f */
  gsl_vector *x_trial;       /* trial parameter vector */
  gsl_vector *f_trial;       /* trial function vector */
  gsl_vector *work;          /* workspace length p */
  long nu;                   /* nu */
  double mu;                 /* LM damping parameter mu */
  double tau;                /* initial scale factor for mu */
  double normf;              /* || f(x) || */
  double normf_trial;        /* || f(x + dx) || */
  int eval_J;                /* 1 if we are currently accumulating full (J,f) */

  gsl_multilarge_linear_workspace *linear_workspace_p;
} lmn_state_t;

#define LM_ONE_THIRD         (0.333333333333333)

static void *lmn_alloc (const gsl_multilarge_linear_type * T,
                        const size_t p);
static void *lmn_alloc_normal (const size_t p);
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

  state->A = gsl_matrix_alloc(p, p);
  if (state->A == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for A", GSL_ENOMEM);
    }

  state->diag = gsl_vector_alloc(p);
  if (state->diag == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for diag", GSL_ENOMEM);
    }

  state->rhs = gsl_vector_alloc(p);
  if (state->rhs == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for rhs", GSL_ENOMEM);
    }

  state->work = gsl_vector_alloc(p);
  if (state->work == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for work", GSL_ENOMEM);
    }

  state->A_copy = gsl_matrix_alloc(p, p);
  if (state->A_copy == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for A_copy", GSL_ENOMEM);
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

  state->tau = 1.0e-3;
  state->mu = state->tau;
  state->nu = 2;
  state->normf = 0.0;
  state->normf_trial = 0.0;

  return state;
}

static void *
lmn_alloc_normal (const size_t p)
{
  return lmn_alloc(gsl_multilarge_linear_normal, p);
}

static void
lmn_free(void *vstate)
{
  lmn_state_t *state = (lmn_state_t *) vstate;

  if (state->A)
    gsl_matrix_free(state->A);

  if (state->J)
    gsl_matrix_free(state->J);

  if (state->diag)
    gsl_vector_free(state->diag);

  if (state->rhs)
    gsl_vector_free(state->rhs);

  if (state->work)
    gsl_vector_free(state->work);

  if (state->A_copy)
    gsl_matrix_free(state->A_copy);

  if (state->x_trial)
    gsl_vector_free(state->x_trial);

  if (state->f_trial)
    gsl_vector_free(state->f_trial);

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

  /* accumulate initial Jacobian and residual vector into system */
  status = lmn_eval(1, x, fdf, work, state);
  if (status)
    return status;

  fprintf(stderr, "normf = %.12e\n", state->normf);

  /* XXX */
  state->mu = 607079115151.04102;

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
  int status;
  lmn_state_t *state = (lmn_state_t *) vstate;

  if (state->eval_J)
    {
      /*
       * we are accumulating full (J,f) system for a new step size
       * calculation
       */

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
  int foundstep = 0;                    /* found step dx */
  double dF;                            /* F(x) - F(x + dx) */
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
      fprintf(stderr, "mu = %.12e\n", state->mu);
      status = gsl_multilarge_linear_solve(sqrt(state->mu), dx, &rnorm, &snorm, linear_p);
      if (status)
        return status;

      printv_octave(dx, "dx");

      lmn_trial_step(x, dx, x_trial);

      /* accumulate f(x+dx) without Jacobian */
      status = lmn_eval(0, x_trial, fdf, fdf_work, state);
      if (status)
        return status;

      /* compute F(x) - F(x+dx) = 0.5*(||f||^2 - ||f_trial||^2) */
      dF = 0.5 * (state->normf + state->normf_trial) *
                 (state->normf - state->normf_trial);

      exit(1);
    }

  return GSL_SUCCESS;
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
      state->normf = 0.0;
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

static const gsl_multilarge_nlinear_type lmnormal_type =
{
  "lmnormal",
  &gsl_multilarge_linear_normal,
  lmn_alloc_normal,
  lmn_init,
  lmn_accumulate,
  lmn_iterate,
  lmn_free
};

const gsl_multilarge_nlinear_type *gsl_multilarge_nlinear_lmnormal = &lmnormal_type;
