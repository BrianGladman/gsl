/* multilargenlin/lm.c
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
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_multilarge_nlin.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>

#include "oct.h"
#define DEBUG   0

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
  gsl_vector *dx;            /* step vector, size p */
  gsl_vector *x_trial;       /* trial parameter vector, size p */
  gsl_vector *DTD;           /* diagonal elements of D^T D matrix added to normal equations, size p */
  gsl_vector *tau_qr;        /* Householder scalars for QR decomposition, size p */
  gsl_vector *workp;         /* temporary workspace, size p */
  long nu;                   /* nu */
  double lambda;             /* LM damping parameter */
  double tau;                /* initial scale factor for lambda */
  double normf;              /* || f(x) || */
  double normf_trial;        /* || f(x + dx) || */
  int scale;                 /* use scale matrix in each LM step */

  gsl_eigen_symm_workspace *eigen_p;
} lm_state_t;

#define GSL_LM_ONE_THIRD         (0.333333333333333)

static void *lm_alloc (const int scale, const size_t p);
static void *lm_alloc_noscale (const size_t p);
static void *lm_alloc_scale (const size_t p);
static void lm_free(void *vstate);
static int lm_init(const gsl_vector * x, gsl_multilarge_function_fdf * fdf,
                   void * vstate);
static int lm_iterate(gsl_vector * x, gsl_vector * dx,
                      gsl_multilarge_function_fdf * fdf, void * vstate);
static gsl_vector *lm_gradient(void * vstate);
static double lm_normf(void * vstate);
static int lm_rcond(double *rcond, void * vstate);
static int lm_eval(const int evaldf, const gsl_vector * x,
                   gsl_multilarge_function_fdf * fdf,
                   lm_state_t * state);
static int lm_solve(gsl_vector * dx, lm_state_t * state);
static int lm_solve_regularize(lm_state_t * state);
static int lm_solve_cholesky(gsl_matrix * A, const gsl_vector * b, gsl_vector * x,
                             lm_state_t * state);
static int lm_solve_QR(gsl_matrix * A, const gsl_vector * b,
                       gsl_vector * x, lm_state_t *state);
static void lm_trial_step(const gsl_vector * x, const gsl_vector * dx,
                          gsl_vector * x_trial);
static double lm_calc_rho(const double lambda, const gsl_vector * dx,
                          const gsl_vector * minus_g, lm_state_t * state);
static int lm_update_diag(lm_state_t * state);

static void *
lm_alloc (const int scale, const size_t p)
{
  lm_state_t *state;

  if (p == 0)
    {
      GSL_ERROR_NULL ("p must be a positive integer",
                      GSL_EINVAL);
    }
  
  state = calloc(1, sizeof(lm_state_t));
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

  state->dx = gsl_vector_alloc(p);
  if (state->dx == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for dx", GSL_ENOMEM);
    }

  state->x_trial = gsl_vector_alloc(p);
  if (state->x_trial == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for x_trial", GSL_ENOMEM);
    }

  state->DTD = gsl_vector_alloc(p);
  if (state->DTD == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for DTD", GSL_ENOMEM);
    }

  state->tau_qr = gsl_vector_alloc(p);
  if (state->tau_qr == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for Householder vector", GSL_ENOMEM);
    }

  state->workp = gsl_vector_alloc(p);
  if (state->workp == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for temporary p workspace", GSL_ENOMEM);
    }

  state->eigen_p = gsl_eigen_symm_alloc(p);
  if (state->eigen_p == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for eigen workspace", GSL_ENOMEM);
    }

  state->p = p;
  state->tau = 1.0e-3;
  state->lambda = 0.0;
  state->nu = 2;
  state->normf = 0.0;
  state->normf_trial = 0.0;
  state->scale = scale;

  return state;
}

static void *
lm_alloc_noscale (const size_t p)
{
  return lm_alloc(0, p);
}

static void *
lm_alloc_scale (const size_t p)
{
  return lm_alloc(1, p);
}

static void
lm_free(void *vstate)
{
  lm_state_t *state = (lm_state_t *) vstate;

  if (state->A)
    gsl_matrix_free(state->A);

  if (state->work_A)
    gsl_matrix_free(state->work_A);

  if (state->rhs)
    gsl_vector_free(state->rhs);

  if (state->dx)
    gsl_vector_free(state->dx);

  if (state->x_trial)
    gsl_vector_free(state->x_trial);

  if (state->DTD)
    gsl_vector_free(state->DTD);

  if (state->tau_qr)
    gsl_vector_free(state->tau_qr);

  if (state->workp)
    gsl_vector_free(state->workp);

  if (state->eigen_p)
    gsl_eigen_symm_free(state->eigen_p);

  free(state);
}

/*
lm_init()
  Initialize LM solver

Inputs: x      - initial parameter vector
        fdf    - user-supplied callback function
        vstate - local workspace
*/

static int
lm_init(const gsl_vector * x, gsl_multilarge_function_fdf * fdf,
        void * vstate)
{
  int status;
  lm_state_t *state = (lm_state_t *) vstate;

  /* compute initial J^T J, J^T f, and ||f|| */
  status = lm_eval(1, x, fdf, state);
  if (status)
    return status;

  state->lambda = state->tau;

  if (!state->scale)
    {
      /* lambda = tau * max [ A_ii ] with DTD = I */
      gsl_vector_view diag = gsl_matrix_diagonal(state->A);
      state->lambda *= gsl_vector_max(&diag.vector);
    }

  /* initialize scaling matrix DTD */
  gsl_vector_set_zero(state->DTD);
  lm_update_diag(state);

  return GSL_SUCCESS;
}

/*
lm_iterate()
  Performe one iteration of LM solver with Nielsen
updating criteria for damping parameter.

Inputs: x        - (input/output)
                   on input, current parameter vector
                   on output, new vector x <- x + dx
        dx       - (output) new step size dx, size p
        fdf      - user-supplied (J,f) function
        vstate   - workspace

Notes:
1) After each iteration,
state->A contains J(x)^T J(x)
state->rhs contains -J(x)^T f(x)
state->normf contains ||f(x)||
*/

static int
lm_iterate(gsl_vector * x, gsl_vector * dx,
           gsl_multilarge_function_fdf * fdf,
           void * vstate)
{
  int status;
  lm_state_t *state = (lm_state_t *) vstate;
  gsl_vector *x_trial = state->x_trial; /* trial x + dx */
  gsl_vector *rhs = state->rhs;         /* negative gradient -J^T f */
  int foundstep = 0;                    /* found step dx */
  double rho;                           /* ratio [ F(x) - F(x + dx) ] / [ L(0) - L(dx) ] */
  int bad_steps = 0;                    /* number of consecutive rejected steps */

  /* loop until we find an acceptable step dx */
  while (!foundstep)
    {
      /* solve: (J^T J + lambda*D) dx = - J^T f */
      status = lm_solve(state->dx, state);
      if (status)
        return status;

      /* compute x_trial = x + dx */
      lm_trial_step(x, state->dx, x_trial);

      /* compute ||f(x+dx)|| */
      status = lm_eval(0, x_trial, fdf, state);
      if (status)
        return status;

      rho = lm_calc_rho(state->lambda, state->dx, rhs, state);

#if DEBUG
      /*XXX*/
      fprintf(stderr, "rho = %.12e, ||dx|| = %.12e, ||x|| = %.12e, ||x+dx|| = %.12e, lambda = %.12e\n",
               rho, gsl_blas_dnrm2(state->dx), gsl_blas_dnrm2(x), gsl_blas_dnrm2(x_trial), state->lambda);
#endif

      /* check that rho > 0 */
      if (rho > 0.0)
        {
          /* reduction in error, step acceptable */

          double b;

          /* update LM parameter */
          b = 2.0 * rho - 1.0;
          b = 1.0 - b*b*b;
          state->lambda *= GSL_MAX(GSL_LM_ONE_THIRD, b);
          state->nu = 2;

          /* compute new JTJ(x+dx), JTf(x+dx) and ||f(x+dx)|| */
          status = lm_eval(1, x_trial, fdf, state);
          if (status)
            return status;

          /* update x <- x + dx */
          gsl_vector_memcpy(x, x_trial);

          /* store step dx */
          gsl_vector_memcpy(dx, state->dx);

          /* update diag(D^T D) if using scaling */
          lm_update_diag(state);

          foundstep = 1;
          bad_steps = 0;
        }
      else
        {
          long nu2;

          /* step did not reduce error, reject step */

          /* if more than 15 consecutive rejected steps, report no progress */
          if (++bad_steps > 15)
            return GSL_ENOPROG;

          state->lambda *= (double) state->nu;

          nu2 = state->nu << 1; /* 2*nu */
          if (nu2 <= state->nu)
            {
              gsl_vector_view diag = gsl_matrix_diagonal(state->A);

              /*
               * nu has wrapped around / overflown, reset lambda and nu
               * to original values and break to force another iteration
               */
              state->nu = 2;

              if (state->scale)
                state->lambda = state->tau;
              else
                state->lambda = state->tau * gsl_vector_max(&diag.vector);

              break;
            }

          state->nu = nu2;
        }
    }

  return GSL_SUCCESS;
}

static gsl_vector *
lm_gradient(void * vstate)
{
  lm_state_t *state = (lm_state_t *) vstate;
  return state->rhs;
}

static double
lm_normf(void * vstate)
{
  lm_state_t *state = (lm_state_t *) vstate;
  return state->normf;
}

static int
lm_rcond(double *rcond, void * vstate)
{
  lm_state_t *state = (lm_state_t *) vstate;
  int status;
  double eval_min, eval_max;
  gsl_vector *eval = state->workp; /* temporary workspace for eigenvalues */

  /* compute eigenvalues of A = J^T J */

  /* copy lower triangle of A to temporary workspace */
  gsl_matrix_tricpy('L', 1, state->work_A, state->A);

  /* compute eigenvalues of A */
  status = gsl_eigen_symm(state->work_A, eval, state->eigen_p);
  if (status)
    return status;

  gsl_vector_minmax(eval, &eval_min, &eval_max);

#if 0
  /*XXX*/
  {
    gsl_vector_scale(eval, 1.0 / eval_max);
    gsl_sort_vector(eval);
    printv_octave(eval, "eval");
  }
#endif

  if (eval_max > 0.0 && eval_min > 0.0)
    {
      *rcond = sqrt(eval_min / eval_max);
    }
  else
    {
      /*
       * computed eigenvalues are not accurate, probably due
       * to rounding errors in forming J^T J
       */
      *rcond = 0.0;
    }

  return GSL_SUCCESS;
}

/*
lm_eval()
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
lm_eval(const int evaldf, const gsl_vector * x,
        gsl_multilarge_function_fdf * fdf,
        lm_state_t * state)
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
lm_solve(gsl_vector * dx, lm_state_t * state)
{
  int status;
  gsl_error_handler_t *err_handler;

  /* compute: work_A = A + lambda*D */
  status = lm_solve_regularize(state);
  if (status)
    return status;

  /* turn off error handler in case Cholesky fails */
  err_handler = gsl_set_error_handler_off();

  /* solve: (A + lambda*D) x = b using Cholesky decomposition */
  status = lm_solve_cholesky(state->work_A, state->rhs, dx, state);

  /* restore error handler */
  gsl_set_error_handler(err_handler);

  if (status)
    {
      /* Cholesky failed, restore matrix and use QR */
      lm_solve_regularize(state);

      status = lm_solve_QR(state->work_A, state->rhs, dx, state);
      if (status)
        return status;
    }

  return GSL_SUCCESS;
}

static int
lm_solve_regularize(lm_state_t * state)
{
  /* copy lower triangle of A to work_A */
  gsl_matrix_tricpy('L', 1, state->work_A, state->A);

  if (state->scale)
    {
      /* add lambda*D to A; state->DTD contains diag(D^T D) */
      size_t i;

      for (i = 0; i < state->p; ++i)
        {
          double Di = gsl_vector_get(state->DTD, i);
          double *Aii = gsl_matrix_ptr(state->work_A, i, i);
          *Aii += state->lambda * Di;
        }
    }
  else
    {
      /* add lambda*I to A */
      gsl_vector_view diag = gsl_matrix_diagonal(state->work_A);
      gsl_vector_add_constant(&diag.vector, state->lambda);
    }

  return GSL_SUCCESS;
}

static int
lm_solve_cholesky(gsl_matrix * A, const gsl_vector * b, gsl_vector * x,
                  lm_state_t * state)
{
  int status;

  /* compute Cholesky decomposition of A^T A with scaling */
  status = gsl_linalg_cholesky_decomp2(A, state->workp);
  if (status)
    return status;

  status = gsl_linalg_cholesky_solve2(A, state->workp, b, x);
  if (status)
    return status;

  return GSL_SUCCESS;
}

static int
lm_solve_QR(gsl_matrix * A, const gsl_vector * b,
            gsl_vector * x, lm_state_t *state)
{
  int status;
  gsl_vector *D = state->workp;
  gsl_vector *tau = state->tau_qr;

  /* scale: A <- diag(D) A diag(D) to try to reduce cond(A) */
  status = gsl_linalg_cholesky_scale(A, D);
  if (status)
    return status;

  /* copy lower triangle of A to upper */
  gsl_matrix_transpose_tricpy('L', 0, A, A);

  status = gsl_linalg_QR_decomp(A, tau);
  if (status)
    return status;

  /* scale rhs vector: b <- diag(D) b */
  gsl_vector_memcpy(x, b);
  gsl_vector_mul(x, D);

  /* solve system */
  status = gsl_linalg_QR_svx(A, tau, x);
  if (status)
    return status;

  /* undo scaling */
  gsl_vector_mul(x, D);

  return GSL_SUCCESS;
}

/* compute x_trial = x + dx */
static void
lm_trial_step(const gsl_vector * x, const gsl_vector * dx,
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
lm_calc_rho()
  Calculate ratio of actual reduction to predicted
reduction, given by Eq 4.4 of More, 1978.
*/

static double
lm_calc_rho(const double lambda, const gsl_vector * dx,
            const gsl_vector * minus_g, lm_state_t * state)
{
  double rho;
  double actual_reduction;
  double pred_reduction;
  double u, norm_dx;
  size_t i;

  /* if ||f(x+dx)|| > ||f(x)|| reject step immediately */
  if (state->normf_trial >= state->normf)
    return 0.0;

  /* compute numerator of rho */
  u = state->normf_trial / state->normf;
  actual_reduction = 1.0 - u*u;

  /* compute || D dx || */
  if (state->scale)
    {
      for (i = 0; i < state->p; ++i)
        {
          double dxi = gsl_vector_get(dx, i);
          double Di = gsl_vector_get(state->DTD, i);

          gsl_vector_set(state->workp, i, dxi * sqrt(Di));
        }

      norm_dx = gsl_blas_dnrm2(state->workp);
    }
  else
    norm_dx = gsl_blas_dnrm2(dx); /* D = I */

  /* compute denominator of rho */
  u = norm_dx / state->normf;
  pred_reduction = lambda * u * u;

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

/*
lm_update_diag()
  After a new step x has been accepted, update the
diag vector, representing the elements of diag(D^T D). This
is done according to Eq. 6.3 of More, 1978.
*/

static int
lm_update_diag(lm_state_t * state)
{
  if (state->scale)
    {
      size_t i;

      for (i = 0; i < state->p; ++i)
        {
          double *diagi = gsl_vector_ptr(state->DTD, i);
          double Aii = gsl_matrix_get(state->A, i, i);

          if (Aii == 0.0)
            Aii = 1.0;

          *diagi = GSL_MAX(*diagi, Aii);
        }
    }

  return GSL_SUCCESS;
}

static const gsl_multilarge_nlinear_type lm_type =
{
  "lm_nielsen",
  lm_alloc_noscale,
  lm_init,
  lm_iterate,
  lm_gradient,
  lm_normf,
  lm_rcond,
  lm_free
};

static const gsl_multilarge_nlinear_type lmns_type =
{
  "lm_nielsen_scaled",
  lm_alloc_scale,
  lm_init,
  lm_iterate,
  lm_gradient,
  lm_normf,
  lm_rcond,
  lm_free
};

const gsl_multilarge_nlinear_type *gsl_multilarge_nlinear_lm = &lm_type;
const gsl_multilarge_nlinear_type *gsl_multilarge_nlinear_lms = &lmns_type;
