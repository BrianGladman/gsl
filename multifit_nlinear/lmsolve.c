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
 * [    J     ] dx = - [ f ]
 * [ lambda*D ]        [ 0 ]
 *
 * using either the normal equations or QR approach
 */

static int normal_init(const gsl_matrix * J, void * vstate);
static int normal_solve(const double lambda, gsl_vector *dx,
                        void *vstate);
static int normal_regularize(const double lambda,
                             const gsl_vector * diag, gsl_matrix * A);
static int normal_solve_cholesky(gsl_matrix * A, const gsl_vector * b,
                                 gsl_vector * x, lm_state_t * state);
static int normal_solve_QR(gsl_matrix * A, const gsl_vector * b,
                           gsl_vector * x, lm_state_t *state);

/* compute A = J^T J */
static int
normal_init(const gsl_matrix * J, void * vstate)
{
  lm_state_t *state = (lm_state_t *) vstate;
  return gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, J, 0.0, state->A);
}

/* compute step dx by solving (J^T J + lambda D^T D) dx = -J^T f */
static int
normal_solve(const double lambda, gsl_vector *dx, void *vstate)
{
  lm_state_t *state = (lm_state_t *) vstate;
  int status;
  gsl_matrix *JTJ = state->work_JTJ;
  gsl_error_handler_t *err_handler;

  /* copy lower triangle of A to workspace */
  gsl_matrix_tricpy('L', 1, JTJ, state->A);

  /* augment normal equations with LM term: A -> A + lambda D^T D */
  status = normal_regularize(lambda, state->diag, JTJ);
  if (status)
    return status;

  /* turn off error handler in case Cholesky fails */
  err_handler = gsl_set_error_handler_off();

  status = normal_solve_cholesky(JTJ, state->rhs, dx, state);

  /* restore error handler */
  gsl_set_error_handler(err_handler);

  if (status)
    {
      /* Cholesky failed, restore matrix and use QR */
      gsl_matrix_tricpy('L', 1, JTJ, state->A);
      normal_regularize(lambda, state->diag, JTJ);

      status = normal_solve_QR(JTJ, state->rhs, dx, state);
      if (status)
        return status;
    }

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

/* solve: A x = b using Cholesky decomposition of A */
static int
normal_solve_cholesky(gsl_matrix * A, const gsl_vector * b,
                      gsl_vector * x, lm_state_t * state)
{
  int status;

  /* compute Cholesky decomposition of A with scaling */
  status = gsl_linalg_cholesky_decomp2(A, state->workp);
  if (status)
    return status;

  status = gsl_linalg_cholesky_solve2(A, state->workp, b, x);
  if (status)
    return status;

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
