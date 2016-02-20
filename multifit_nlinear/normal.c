/* multifit_nlinear/normal.c
 * 
 * Copyright (C) 2015, 2016 Patrick Alken
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
 * This module calculates the solution of the linear least squares
 * system:
 *
 * [ J^T J + mu D^T D ] v = -J^T f, for geodesic velocity
 *
 * and
 *
 * [ J^T J + mu D^T D ] a = -J^T fvv, for geodesic acceleration
 */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>

typedef struct
{
  gsl_matrix *A;             /* J^T J */
  gsl_matrix *work_JTJ;      /* copy of J^T J */
  gsl_vector *rhs_vel;       /* -J^T f */
  gsl_vector *rhs_acc;       /* -J^T fvv */
  gsl_vector *tau;           /* Householder scalars for QR, size p */
  gsl_vector *workp;         /* workspace, size p */
  int chol;                  /* Cholesky factorization successful */
} normal_state_t;

static void *normal_alloc (const size_t n, const size_t p);
static int normal_init(const gsl_matrix * J, void * vstate);
static int normal_presolve(const double mu, const gsl_vector * diag, void * vstate);
static int normal_solve_rhs(const gsl_vector * b, gsl_vector *x, normal_state_t *state);
static int normal_regularize(const double mu,
                             const gsl_vector * diag, gsl_matrix * A);

static void *
normal_alloc (const size_t n, const size_t p)
{
  normal_state_t *state;

  (void)n;
  
  state = calloc(1, sizeof(normal_state_t));
  if (state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate normal state", GSL_ENOMEM);
    }

  state->A = gsl_matrix_alloc(p, p);
  if (state->A == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for A", GSL_ENOMEM);
    }

  state->work_JTJ = gsl_matrix_alloc(p, p);
  if (state->work_JTJ == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for JTJ workspace",
                      GSL_ENOMEM);
    }

  state->rhs_vel = gsl_vector_alloc(p);
  if (state->rhs_vel == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for rhs_vel", GSL_ENOMEM);
    }

  state->rhs_acc = gsl_vector_alloc(p);
  if (state->rhs_acc == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for rhs_acc", GSL_ENOMEM);
    }

  state->tau = gsl_vector_alloc(p);
  if (state->tau == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for tau", GSL_ENOMEM);
    }

  state->workp = gsl_vector_alloc(p);
  if (state->workp == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for workp", GSL_ENOMEM);
    }

  return state;
}

static void
normal_free(void *vstate)
{
  normal_state_t *state = (normal_state_t *) vstate;

  if (state->A)
    gsl_matrix_free(state->A);

  if (state->work_JTJ)
    gsl_matrix_free(state->work_JTJ);

  if (state->rhs_vel)
    gsl_vector_free(state->rhs_vel);

  if (state->rhs_acc)
    gsl_vector_free(state->rhs_acc);

  if (state->tau)
    gsl_vector_free(state->tau);

  if (state->workp)
    gsl_vector_free(state->workp);

  free(state);
}

/* compute A = J^T J */
static int
normal_init(const gsl_matrix * J, void * vstate)
{
  normal_state_t *state = (normal_state_t *) vstate;

  /* compute A = J^T J */
  gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, J, 0.0, state->A);

  return GSL_SUCCESS;
}

/*
normal_presolve()
  Compute the Cholesky decomposition of J^T J + mu * D^T D

Inputs: mu     - LM parameter
        diag   - diag(D)
        vstate - workspace

Notes:
1) On output, state->work_JTJ contains the Cholesky decomposition of
J^T J + mu D^T D

2) On output, state->workp contains scale factors needed for a
solution

3) If Cholesky fails, QR decompsition is used and on output, state->tau
contains Householder scalars
*/

static int
normal_presolve(const double mu, const gsl_vector * diag, void * vstate)
{
  normal_state_t *state = (normal_state_t *) vstate;
  gsl_matrix *JTJ = state->work_JTJ;
  gsl_vector *chol_D = state->workp;
  gsl_error_handler_t *err_handler;
  int status;

  /* copy lower triangle of A to workspace */
  gsl_matrix_tricpy('L', 1, JTJ, state->A);

  /* augment normal equations with LM term: A -> A + mu D^T D */
  status = normal_regularize(mu, diag, JTJ);
  if (status)
    return status;

  /* turn off error handler in case Cholesky fails */
  err_handler = gsl_set_error_handler_off();

  status = gsl_linalg_cholesky_decomp2(JTJ, chol_D);

  /* restore error handler */
  gsl_set_error_handler(err_handler);

  if (status == GSL_SUCCESS)
    {
      state->chol = 1;
    }
  else
    {
      state->chol = 0;

      /* Cholesky failed, restore matrix and use QR decomposition */
      gsl_matrix_tricpy('L', 1, JTJ, state->A);
      normal_regularize(mu, diag, JTJ);

      /* scale: JTJ <- diag(D) JTJ diag(D) to try to reduce cond(JTJ) */
      status = gsl_linalg_cholesky_scale(JTJ, chol_D);
      if (status)
        return status;

      /* copy lower triangle of A to upper */
      gsl_matrix_transpose_tricpy('L', 0, JTJ, JTJ);

      status = gsl_linalg_QR_decomp(JTJ, state->tau);
      if (status)
        return status;
    }

  return GSL_SUCCESS;
}

/*
normal_solve()
  Compute (J^T J + mu D^T D) x = -J^T f

Inputs: f      - right hand side vector f
        g      - J^T f (can be NULL)
        x      - (output) solution vector
        vstate - normal workspace
*/
static int
normal_solve(const gsl_vector * f, const gsl_vector * g,
             gsl_vector *x, void *vstate)
{
  normal_state_t *state = (normal_state_t *) vstate;
  int status;

  if (g != NULL)
    {
      status = normal_solve_rhs(g, x, state);
      if (status)
        return status;

      /* reverse step to go downhill */
      gsl_vector_scale(x, -1.0);
    }
  else
    {
#if 0
      /* compute rhs = -J^T f */
      gsl_blas_dgemv(CblasTrans, -1.0, J, f, 0.0, state->rhs_acc);

      status = normal_solve(state->rhs_acc, x, state);
      if (status)
        return status;
#endif
    }

  return GSL_SUCCESS;
}

/* solve: (J^T J + mu D^T D) x = b */
static int
normal_solve_rhs(const gsl_vector * b, gsl_vector *x, normal_state_t *state)
{
  int status;
  gsl_matrix *JTJ = state->work_JTJ;
  gsl_vector *chol_D = state->workp;

  if (state->chol == 1)
    {
      /* we have a Cholesky factorization of J^T J + mu D^T D */
      status = gsl_linalg_cholesky_solve2(JTJ, chol_D, b, x);
      if (status)
        return status;
    }
  else
    {
      /* scale rhs vector: b <- diag(D) b */
      gsl_vector_memcpy(x, b);
      gsl_vector_mul(x, chol_D);

      /* solve system */
      status = gsl_linalg_QR_svx(JTJ, state->tau, x);
      if (status)
        return status;

      /* undo scaling */
      gsl_vector_mul(x, chol_D);
    }

  return GSL_SUCCESS;
}

/* A <- A + mu D^T D */
static int
normal_regularize(const double mu, const gsl_vector * diag,
                  gsl_matrix * A)
{
  const size_t p = diag->size;
  size_t i;

  for (i = 0; i < p; ++i)
    {
      double di = gsl_vector_get(diag, i);
      double *Aii = gsl_matrix_ptr(A, i, i);

      *Aii += mu * di * di;
    }

  return GSL_SUCCESS;
}

static const gsl_multifit_nlinear_solver normal_type =
{
  "normal",
  normal_alloc,
  normal_init,
  normal_presolve,
  normal_solve,
  normal_free
};

const gsl_multifit_nlinear_solver *gsl_multifit_nlinear_solver_normal = &normal_type;
