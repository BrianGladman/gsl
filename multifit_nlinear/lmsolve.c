/* multifit_nlinear/lmsolve.c
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

/*
 * This module handles the solution of the linear least squares
 * system:
 *
 * [    J     ] v = - [ f ]
 * [ lambda*D ]       [ 0 ]
 *
 * using either the normal equations or QR approach. The solvers are
 * organized into 4 sequential steps:
 *
 * 1. init1: initialize solver for a given f(x) and J(x), independent of lambda
 * 2. init2: further solver initialization once lambda is selected
 * 3. solve_vel: solve above linear system for geodesic velocity (this is the
 *               standard LM step)
 * 4. solve_acc: solve a similar linear system for the geodesic acceleration:
 *
 * [    J     ] a = - [ fvv ]
 * [ lambda*D ]       [  0  ]
 *
 * only the right hand side is different in this case (instead of f(x) we
 * use the second directional derivative in the velocity direction).
 */

#include "qrsolv.c"

static int normal_init1(const gsl_vector * f, const gsl_matrix * J,
                        const gsl_vector * g, void * vstate);
static int normal_init2(const double lambda, void * vstate);
static int normal_solve_vel(gsl_vector *v, void *vstate);
static int normal_solve_acc(const gsl_matrix *J, const gsl_vector *fvv,
                            gsl_vector *a, void *vstate);
static int normal_solve(const gsl_vector * b, gsl_vector *x, lm_state_t *state);
static int normal_regularize(const double lambda,
                             const gsl_vector * diag, gsl_matrix * A);
static int normal_solve_QR(gsl_matrix * A, const gsl_vector * b,
                           gsl_vector * x, lm_state_t *state);

static int qr_init(const gsl_vector * f, const gsl_matrix * J,
                   void * vstate);
static int qr_solve(const double lambda, const gsl_vector *g,
                    gsl_vector *dx, void *vstate);

/* compute A = J^T J */
static int
normal_init1(const gsl_vector * f, const gsl_matrix * J,
             const gsl_vector * g, void * vstate)
{
  lm_state_t *state = (lm_state_t *) vstate;

  (void)f; /* avoid unused parameter warning */

  /* prepare rhs vector = -g = -J^T f */
  gsl_vector_memcpy(state->rhs_vel, g);
  gsl_vector_scale(state->rhs_vel, -1.0);

  /* compute A = J^T J */
  gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, J, 0.0, state->A);

  return GSL_SUCCESS;
}

/*
normal_init2()
  Compute the Cholesky decomposition of J^T J + lambda * D^T D

Inputs: lambda - LM parameter
        vstate - workspace

Notes:
1) On output, state->work_JTJ contains the Cholesky decomposition of
J^T J + lambda D^T D

2) On output, state->workp contains scale factors needed for a
solution of the Cholesky system
*/

static int
normal_init2(const double lambda, void * vstate)
{
  lm_state_t *state = (lm_state_t *) vstate;
  gsl_matrix *JTJ = state->work_JTJ;
  gsl_error_handler_t *err_handler;
  int status;

  /* copy lower triangle of A to workspace */
  gsl_matrix_tricpy('L', 1, JTJ, state->A);

  /* augment normal equations with LM term: A -> A + lambda D^T D */
  status = normal_regularize(lambda, state->diag, JTJ);
  if (status)
    return status;

  /* turn off error handler in case Cholesky fails */
  err_handler = gsl_set_error_handler_off();

  status = gsl_linalg_cholesky_decomp2(JTJ, state->workp);

  /* restore error handler */
  gsl_set_error_handler(err_handler);

  if (status == GSL_SUCCESS)
    {
      state->chol = 1;
    }
  else
    {
      state->chol = 0;
    }

  return GSL_SUCCESS;
}

/* compute velocity by solving (J^T J + lambda D^T D) v = -J^T f */
static int
normal_solve_vel(gsl_vector *v, void *vstate)
{
  lm_state_t *state = (lm_state_t *) vstate;
  int status;

  status = normal_solve(state->rhs_vel, v, state);

  return status;
}

/* compute acceleration by solving (J^T J + lambda D^T D) a = -J^T fvv */
static int
normal_solve_acc(const gsl_matrix *J, const gsl_vector *fvv,
                 gsl_vector *a, void *vstate)
{
  lm_state_t *state = (lm_state_t *) vstate;
  int status;

  /* compute rhs = -J^T fvv */
  gsl_blas_dgemv(CblasTrans, -1.0, J, fvv, 0.0, state->rhs_acc);

  status = normal_solve(state->rhs_acc, a, state);

  return status;
}

/* solve: (J^T J + lambda D^T D) x = b */
static int
normal_solve(const gsl_vector * b, gsl_vector *x, lm_state_t *state)
{
  int status;
  gsl_matrix *JTJ = state->work_JTJ;

  if (state->chol == 1)
    {
      /* we have a Cholesky factorization of J^T J + lambda D^T D */
      status = gsl_linalg_cholesky_solve2(JTJ, state->workp, b, x);
      if (status)
        return status;
    }
  else
    {
    }

#if 0
  if (status)
    {
      /* Cholesky failed, restore matrix and use QR */
      gsl_matrix_tricpy('L', 1, JTJ, state->A);
      normal_regularize(lambda, state->diag, JTJ);

      status = normal_solve_QR(JTJ, g, dx, state);
      if (status)
        return status;
    }

  /* reverse step to go downhill */
  gsl_vector_scale(dx, -1.0);
#endif

  return GSL_SUCCESS;
}

/* A <- A + lambda D^T D */
static int
normal_regularize(const double lambda, const gsl_vector * diag,
                  gsl_matrix * A)
{
  const size_t p = diag->size;
  size_t i;

  for (i = 0; i < p; ++i)
    {
      double di = gsl_vector_get(diag, i);
      double *Aii = gsl_matrix_ptr(A, i, i);

      *Aii += lambda * di * di;
    }

  return GSL_SUCCESS;
}

static int
normal_solve_QR(gsl_matrix * A, const gsl_vector * b,
                gsl_vector * x, lm_state_t *state)
{
  int status;
  gsl_vector *D = state->workp;
  gsl_vector *tau = state->tau;

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

/* compute J = Q R PT and qtf = Q^T f */
static int
qr_init(const gsl_vector * f, const gsl_matrix * J, void * vstate)
{
  lm_state_t *state = (lm_state_t *) vstate;
  int signum;

  gsl_matrix_memcpy(state->R, J);
  gsl_linalg_QRPT_decomp(state->R, state->tau, state->perm,
                         &signum, state->workp);

  gsl_vector_memcpy(state->qtf, f);
  gsl_linalg_QR_QTvec(state->R, state->tau, state->qtf);

  return GSL_SUCCESS;
}

static int
qr_solve(const double lambda, const gsl_vector *g, gsl_vector *dx, void *vstate)
{
  lm_state_t *state = (lm_state_t *) vstate;
  const double sqrt_lambda = sqrt(lambda);
  int status;

  /*
   * solve:
   *
   * [       J       ] dx = [ f ]
   * [ sqrt(lamba) D ]      [ 0 ]
   *
   * using QRPT factorization of J
   */
  status = qrsolv(state->R, state->perm, sqrt_lambda, state->diag,
                  state->qtf, dx, state->workp, state->workn);

  /* reverse step to go downhill */
  gsl_vector_scale(dx, -1.0);

  (void)g; /* avoid unused parameter warning */

  return status;
}
