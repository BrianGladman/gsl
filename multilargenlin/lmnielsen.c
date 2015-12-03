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
  gsl_matrix *A;             /* J^T J, size p-by-p */
  gsl_matrix *work_A;        /* temporary J^T J matrix, size p-by-p */
  gsl_vector *rhs;           /* rhs vector, -J^T f, size p */
  gsl_vector *x_trial;       /* trial parameter vector, size p */
  gsl_vector *diag;          /* diagonal elements of D^T D matrix added to normal equations, size p */
  gsl_vector *D;             /* scale factors for Cholesky decomposition, size p */
  long nu;                   /* nu */
  double mu;                 /* LM damping parameter mu */
  double tau;                /* initial scale factor for mu */
  double normf;              /* || f(x) || */
  double normf_trial;        /* || f(x + dx) || */
  int scale;                 /* use scale matrix in each LM step */
} lmn_state_t;

#define GSL_LM_SQRT_ONE_THIRD    (0.577350269189626)
#define GSL_LM_ONE_THIRD         (0.333333333333333)

static void *lmn_alloc (const int scale, const size_t p);
static void lmn_free(void *vstate);
static int lmn_init(const gsl_vector * x, gsl_multilarge_function_fdf * fdf,
                    void * vstate);
static int lmn_iterate(gsl_vector * x, gsl_vector * dx,
                       gsl_multilarge_function_fdf * fdf, void * vstate);
static gsl_vector *lmn_gradient(void * vstate);
static double lmn_normf(void * vstate);
static int lmn_rcond(double *rcond, void * vstate);
static int lmn_eval(const int evaldf, const gsl_vector * x,
                    gsl_multilarge_function_fdf * fdf,
                    lmn_state_t * state);
static int lmn_solve(const double mu, const gsl_matrix * A, const gsl_vector * b,
                     gsl_vector * x, lmn_state_t * state);
static int lmn_solve_cholesky(gsl_matrix * A, const gsl_vector * b, gsl_vector * x,
                              lmn_state_t * state);
static void lmn_trial_step(const gsl_vector * x, const gsl_vector * dx,
                           gsl_vector * x_trial);
static double lmn_calc_rho(const double mu, const gsl_vector * dx,
                           const gsl_vector * minus_g, lmn_state_t * state);

static void *
lmn_alloc (const int scale, const size_t p)
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

  state->work_A = gsl_matrix_alloc(p, p);
  if (state->work_A == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for temporary A workspace", GSL_ENOMEM);
    }

  state->rhs = gsl_vector_alloc(p);
  if (state->rhs == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for rhs", GSL_ENOMEM);
    }

  state->x_trial = gsl_vector_alloc(p);
  if (state->x_trial == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for x_trial", GSL_ENOMEM);
    }

  state->diag = gsl_vector_alloc(p);
  if (state->diag == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for diag", GSL_ENOMEM);
    }

  state->D = gsl_vector_alloc(p);
  if (state->D == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for D", GSL_ENOMEM);
    }

  state->p = p;
  state->tau = 1.0e-3;
  state->mu = 0.0;
  state->nu = 2;
  state->normf = 0.0;
  state->normf_trial = 0.0;
  state->scale = scale;

  return state;
}

static void *
lmn_alloc_noscale (const size_t p)
{
  return lmn_alloc(0, p);
}

static void *
lmn_alloc_scale (const size_t p)
{
  return lmn_alloc(1, p);
}

static void
lmn_free(void *vstate)
{
  lmn_state_t *state = (lmn_state_t *) vstate;

  if (state->A)
    gsl_matrix_free(state->A);

  if (state->work_A)
    gsl_matrix_free(state->work_A);

  if (state->rhs)
    gsl_vector_free(state->rhs);

  if (state->x_trial)
    gsl_vector_free(state->x_trial);

  if (state->diag)
    gsl_vector_free(state->diag);

  if (state->D)
    gsl_vector_free(state->D);

  free(state);
}

/*
lmn_init()
  Initialize LM solver

Inputs: x      - initial parameter vector
        fdf    - user-supplied callback function
        vstate - local workspace
*/

static int
lmn_init(const gsl_vector * x, gsl_multilarge_function_fdf * fdf,
         void * vstate)
{
  int status;
  lmn_state_t *state = (lmn_state_t *) vstate;
  gsl_vector_view diag = gsl_matrix_diagonal(state->A);

  /* compute initial J^T J, J^T f, and ||f|| */
  status = lmn_eval(1, x, fdf, state);
  if (status)
    return status;

  if (state->scale)
    {
      /* set mu = tau and use D^T D = diag(A) */
      gsl_vector_memcpy(state->diag, &diag.vector);
      state->mu = state->tau;
    }
  else
    {
      /* mu = tau * max [ A_ii ] with D = I */
      state->mu = state->tau * gsl_vector_max(&diag.vector);
    }

  return GSL_SUCCESS;
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
        vstate   - workspace
*/

static int
lmn_iterate(gsl_vector * x, gsl_vector * dx,
            gsl_multilarge_function_fdf * fdf,
            void * vstate)
{
  int status;
  lmn_state_t *state = (lmn_state_t *) vstate;
  gsl_vector *x_trial = state->x_trial; /* trial x + dx */
  gsl_vector *rhs = state->rhs;         /* negative gradient -J^T f */
  int foundstep = 0;                    /* found step dx */
  double rho;                           /* ratio [ F(x) - F(x + dx) ] / [ L(0) - L(dx) ] */

  /* loop until we find an acceptable step dx */
  while (!foundstep)
    {
      /* solve: (J^T J + mu D^T D) dx = - J^T f */
      status = lmn_solve(state->mu, state->A, state->rhs, dx, state);
      if (status)
        return status;

      /* compute x_trial = x + dx */
      lmn_trial_step(x, dx, x_trial);

      /* compute ||f(x+dx)|| */
      status = lmn_eval(0, x_trial, fdf, state);
      if (status)
        return status;

      rho = lmn_calc_rho(state->mu, dx, rhs, state);

      /* check that rho > 0 */
      if (rho > 0.0)
        {
          /* reduction in error, step acceptable */

          double b;

          /* update LM parameter mu */
          b = 2.0 * rho - 1.0;
          b = 1.0 - b*b*b;
          if (b > 0.0)
            state->mu *= GSL_MAX(GSL_LM_ONE_THIRD, b);
          else
            state->mu *= GSL_LM_ONE_THIRD;

          state->nu = 2;

          /* compute new JTJ(x+dx), JTf(x+dx) and ||f(x+dx)|| */
          status = lmn_eval(1, x_trial, fdf, state);
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
          state->mu *= (double) state->nu;
          nu2 = state->nu << 1; /* 2*nu */
          if (nu2 <= state->nu)
            {
              gsl_vector_view diag = gsl_matrix_diagonal(state->A);

              /*
               * nu has wrapped around / overflown, reset mu and nu
               * to original values and break to force another iteration
               */
              state->nu = 2;

              if (state->scale)
                state->mu = state->tau;
              else
                state->mu = state->tau * gsl_vector_max(&diag.vector);

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
  return state->rhs;
}

static double
lmn_normf(void * vstate)
{
  lmn_state_t *state = (lmn_state_t *) vstate;
  return state->normf;
}

static int
lmn_rcond(double *rcond, void * vstate)
{
  /*lmn_state_t *state = (lmn_state_t *) vstate;*/
  return 0.0; /* XXX */
}

/*
lmn_eval()
  Evaluate user-supplied Jacobian and residual vector
function. During iteration, it may be necessary to
compute F(x+dx) = 1/2 ||f(x+dx)||^2 without accumulating
the Jacobian. In this case the flag 'evaldf' determines
whether we need the full (J,f) accumulation.

Inputs: evaldf   - if 1, request full (J,f) system from fdf
                   if 0, request only f from fdf
        x        - parameter vector for evaluation of (J,f)
        fdf      - user-supplied function
        state    - local workspace

Notes:
1) If evaldf == 1, on output,
state->A contains J^T J
state->rhs contains -J^T f
state->normf contains ||f(x)||

2) If evaldf == 0, on output,
state->normf_trial contains ||f(x)||
*/

static int
lmn_eval(const int evaldf, const gsl_vector * x,
         gsl_multilarge_function_fdf * fdf,
         lmn_state_t * state)
{
  int status;

  if (evaldf)
    {
      GSL_MULTILARGE_EVAL_FDF(fdf, x, state->A, state->rhs, &(state->normf), status);
      gsl_vector_scale(state->rhs, -1.0);
    }
  else
    GSL_MULTILARGE_EVAL_FDF(fdf, x, NULL, NULL, &(state->normf_trial), status);

  return status;
}

static int
lmn_solve(const double mu, const gsl_matrix * A, const gsl_vector * b,
          gsl_vector * x, lmn_state_t * state)
{
  int status;
  gsl_vector_view diag = gsl_matrix_diagonal(state->work_A);

  gsl_matrix_tricpy('L', 1, state->work_A, A);
  gsl_vector_add_constant(&diag.vector, mu);

  status = lmn_solve_cholesky(state->work_A, b, x, state);

  return status;
}

static int
lmn_solve_cholesky(gsl_matrix * A, const gsl_vector * b, gsl_vector * x,
                   lmn_state_t * state)
{
  int status;

  /* compute Cholesky decomposition of A^T A with scaling */
  status = gsl_linalg_cholesky_decomp2(A, state->D);
  if (status)
    return status;

  status = gsl_linalg_cholesky_solve2(A, state->D, b, x);
  if (status)
    return status;

  return GSL_SUCCESS;
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

static double
lmn_calc_rho(const double mu, const gsl_vector * dx, const gsl_vector * minus_g,
             lmn_state_t * state)
{
  double rho;
  double actual_reduction;
  double pred_reduction;
  double u, norm_dx;
  size_t i;

  /* if ||f(x+dx)|| > ||f(x)|| reject step immediately */
  if (state->normf_trial > state->normf)
    return 0.0;

  /* compute numerator of rho */
  u = state->normf_trial / state->normf;
  actual_reduction = 1.0 - u*u;

  /* compute denominator of rho */
  norm_dx = gsl_blas_dnrm2(dx); /* || D dx || */
  u = norm_dx / state->normf;
  pred_reduction = mu * u * u;

  for (i = 0; i < dx->size; ++i)
    {
      double dxi = gsl_vector_get(dx, i);
      double mgi = gsl_vector_get(minus_g, i);

      pred_reduction += (dxi / state->normf) * (mgi / state->normf);
    }

  if (pred_reduction > 0.0)
    rho = actual_reduction / pred_reduction;
  else
    rho = 0.0;

  return rho;
}

static const gsl_multilarge_nlinear_type lmnielsen_type =
{
  "lmnielsen",
  lmn_alloc_noscale,
  lmn_init,
  lmn_iterate,
  lmn_gradient,
  lmn_normf,
  lmn_rcond,
  lmn_free
};

const gsl_multilarge_nlinear_type *gsl_multilarge_nlinear_lmnielsen = &lmnielsen_type;
