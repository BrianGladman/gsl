/* multilarge_nlinear/lmsolve.c
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

static int lm_solve_decomp(const double mu, const gsl_matrix * JTJ,
                           const gsl_vector * DTD, gsl_matrix * A,
                           lm_state_t * state);
static int lm_solve(const gsl_matrix * A, const gsl_vector * b,
                    gsl_vector * x, lm_state_t * state);
static int lm_solve_regularize(const double mu, const gsl_matrix * JTJ,
                               const gsl_vector * DTD, gsl_matrix * A);

/*
lm_solve_decomp()
  Compute Cholesky decomposition of

  J^T J + mu D^T D

If Cholesky fails, a QR decomposition is used

Inputs: mu    - LM parameter
        JTJ   - J^T J matrix (lower half)
        DTD   - diag(D^T D)
        A     - (output) where to store Cholesky (or QR) decomposition
        state - workspace

Return: success/error

Notes:
1) If Cholesky succeeds, on output:
  a) state->chol is set to 1
  b) state->chol_D contains Cholesky scale factors

2) If Cholesky fails, on output:
  a) state->chol is set to 0
  b) state->chol_D contains Cholesky scale factors
  c) state->tau contains QR Householder scalars
*/

static int
lm_solve_decomp(const double mu, const gsl_matrix * JTJ,
                const gsl_vector * DTD, gsl_matrix * A,
                lm_state_t * state)
{
  int status;
  gsl_error_handler_t *err_handler;

  /* compute: A = J^T J + mu D^T D */
  status = lm_solve_regularize(mu, JTJ, DTD, A);
  if (status)
    return status;

  /* turn off error handler in case Cholesky fails */
  err_handler = gsl_set_error_handler_off();

  /* compute Cholesky decomposition of A with scaling */
  status = gsl_linalg_cholesky_decomp2(A, state->chol_D);

  /* restore error handler */
  gsl_set_error_handler(err_handler);

  if (status)
    {
      state->chol = 0;

      /* re-form A matrix */
      lm_solve_regularize(mu, JTJ, DTD, A);

      /* scale: A <- diag(D) A diag(D) to try to reduce cond(A) */
      status = gsl_linalg_cholesky_scale(A, state->chol_D);
      if (status)
        return status;

      /* copy lower triangle of A to upper */
      gsl_matrix_transpose_tricpy('L', 0, A, A);

      /* perform QR decomposition */
      status = gsl_linalg_QR_decomp(A, state->tau);
      if (status)
        return status;
    }
  else
    {
      state->chol = 1;
    }

  return GSL_SUCCESS;
}

/* solve A x = -b */
static int
lm_solve(const gsl_matrix * A, const gsl_vector * b,
         gsl_vector * x, lm_state_t * state)
{
  int status;
  gsl_vector *D = state->chol_D;

  if (state->chol == 1)
    {
      /* solve: A x = b using Cholesky decomposition of A */
      status = gsl_linalg_cholesky_solve2(A, D, b, x);
      if (status)
        return status;
    }
  else
    {
      /* solve: A x = b using QR decomposition of A */

      /* scale rhs vector: b <- diag(D) b */
      gsl_vector_memcpy(x, b);
      gsl_vector_mul(x, D);

      /* solve system */
      status = gsl_linalg_QR_svx(A, state->tau, x);
      if (status)
        return status;

      /* undo scaling */
      gsl_vector_mul(x, D);
    }

  /* reverse step for downhill direction */
  gsl_vector_scale(x, -1.0);

  return GSL_SUCCESS;
}

/* compute: A = J^T J + mu D^T D */
static int
lm_solve_regularize(const double mu, const gsl_matrix * JTJ,
                    const gsl_vector * DTD, gsl_matrix * A)
{
  size_t i;

  /* copy lower triangle of J^T J to A */
  gsl_matrix_tricpy('L', 1, A, JTJ);

  for (i = 0; i < DTD->size; ++i)
    {
      double Di = gsl_vector_get(DTD, i);
      double *Aii = gsl_matrix_ptr(A, i, i);
      *Aii += mu * Di;
    }

  return GSL_SUCCESS;
}
