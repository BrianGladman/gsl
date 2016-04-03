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
 * [ J^T J + mu*I ] v = -J^T f
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

typedef struct
{
  gsl_matrix *JTJ;           /* J^T J */
  gsl_matrix *work_JTJ;      /* copy of J^T J */
  gsl_vector *rhs;           /* -J^T f */
  gsl_permutation *perm;     /* permutation matrix for modified Cholesky */
} normal_state_t;

static void *normal_alloc (const size_t n, const size_t p);
static int normal_init(const gsl_matrix * J, void * vstate);
static int normal_presolve(const double mu, void * vstate);
static int normal_solve(const gsl_vector * f, const gsl_vector * g,
                        const gsl_matrix * J, gsl_vector *x, void *vstate);
static int normal_solve_rhs(const gsl_vector * b, gsl_vector *x, normal_state_t *state);
static int normal_regularize(const double mu, gsl_matrix * A);

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

  state->JTJ = gsl_matrix_alloc(p, p);
  if (state->JTJ == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for JTJ", GSL_ENOMEM);
    }

  state->work_JTJ = gsl_matrix_alloc(p, p);
  if (state->work_JTJ == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for JTJ workspace",
                      GSL_ENOMEM);
    }

  state->rhs = gsl_vector_alloc(p);
  if (state->rhs == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for rhs", GSL_ENOMEM);
    }

  state->perm = gsl_permutation_alloc(p);
  if (state->perm == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for perm", GSL_ENOMEM);
    }

  return state;
}

static void
normal_free(void *vstate)
{
  normal_state_t *state = (normal_state_t *) vstate;

  if (state->JTJ)
    gsl_matrix_free(state->JTJ);

  if (state->work_JTJ)
    gsl_matrix_free(state->work_JTJ);

  if (state->rhs)
    gsl_vector_free(state->rhs);

  if (state->perm)
    gsl_permutation_free(state->perm);

  free(state);
}

static int
normal_init(const gsl_matrix * J, void * vstate)
{
  normal_state_t *state = (normal_state_t *) vstate;

  /* compute J^T J */
  gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, J, 0.0, state->JTJ);

  return GSL_SUCCESS;
}

/*
normal_presolve()
  Compute the modified Cholesky decomposition of J^T J + mu * I.
Modified Cholesky is used in case mu = 0 and there are rounding
errors in forming J^T J which could lead to an indefinite matrix.

Inputs: mu     - LM parameter
        vstate - workspace

Notes:
1) On output, state->work_JTJ contains the Cholesky decomposition of
J^T J + mu*I
*/

static int
normal_presolve(const double mu, void * vstate)
{
  normal_state_t *state = (normal_state_t *) vstate;
  gsl_matrix *JTJ = state->work_JTJ;
  int status;

  /* copy lower triangle of A to workspace */
  gsl_matrix_tricpy('L', 1, JTJ, state->JTJ);

  /* augment normal equations: A -> A + mu*I */
  status = normal_regularize(mu, JTJ);
  if (status)
    return status;

  /* compute modified Cholesky decomposition */
  status = gsl_linalg_mcholesky_decomp(JTJ, state->perm, NULL);
  if (status)
    return status;

  return GSL_SUCCESS;
}

/*
normal_solve()
  Compute (J^T J + mu D^T D) x = -J^T f

Inputs: f      - right hand side vector f
        g      - J^T f (can be NULL)
        J      - Jacobian matrix
        x      - (output) solution vector
        vstate - normal workspace

Notes:
1) If g is not NULL, it is assumed to equal J^T f
   If g is NULL, J^T f is computed prior to solving
*/

static int
normal_solve(const gsl_vector * f, const gsl_vector * g,
             const gsl_matrix * J, gsl_vector *x, void *vstate)
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
      /* compute rhs = -J^T f */
      gsl_blas_dgemv(CblasTrans, -1.0, J, f, 0.0, state->rhs);

      status = normal_solve_rhs(state->rhs, x, state);
      if (status)
        return status;
    }

  return GSL_SUCCESS;
}

/* solve: (J^T J + mu D^T D) x = b */
static int
normal_solve_rhs(const gsl_vector * b, gsl_vector *x, normal_state_t *state)
{
  int status;
  gsl_matrix *JTJ = state->work_JTJ;

  status = gsl_linalg_mcholesky_solve(JTJ, state->perm, b, x);
  if (status)
    return status;

  return GSL_SUCCESS;
}

/* A <- A + mu*I */
static int
normal_regularize(const double mu, gsl_matrix * A)
{
  if (mu != 0.0)
    {
      gsl_vector_view d = gsl_matrix_diagonal(A);
      gsl_vector_add_constant(&d.vector, mu);
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
