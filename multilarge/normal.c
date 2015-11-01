/* normal.c
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multilarge.h>

typedef struct
{
  size_t nmax;          /* maximum rows to add at once */
  size_t p;             /* number of columns of LS matrix */
  gsl_matrix *ATA;      /* A^T A, p-by-p */
  gsl_vector *ATb;      /* A^T b, p-by-1 */
  double bTb;           /* b^T b */
  gsl_matrix *work_ATA; /* temporary workspace, p-by-p */
} normal_state_t;

static void *normal_alloc(const size_t nmax, const size_t p);
static void normal_free(void *vstate);
static int normal_reset(void *vstate);
static int normal_accumulate(const gsl_matrix * A,
                             const gsl_vector * b,
                             void * vstate);
static int normal_solve(const double lambda, gsl_vector * x,
                        double * rnorm, double * snorm,
                        void * vstate);

/*
normal_alloc()
  Allocate workspace for solving large linear least squares
problems using the normal equations approach

Inputs: nmax - maximum number of rows to accumulate at once
        p    - number of columns of LS matrix

Return: pointer to workspace
*/

static void *
normal_alloc(const size_t nmax, const size_t p)
{
  normal_state_t *state;

  if (nmax == 0)
    {
      GSL_ERROR_NULL("nmax must be a positive integer",
                     GSL_EINVAL);
    }

  if (p == 0)
    {
      GSL_ERROR_NULL("p must be a positive integer",
                     GSL_EINVAL);
    }

  state = calloc(1, sizeof(normal_state_t));
  if (!state)
    {
      GSL_ERROR_NULL("failed to allocate normal state", GSL_ENOMEM);
    }

  state->nmax = nmax;
  state->p = p;

  state->ATA = gsl_matrix_alloc(p, p);
  if (state->ATA == NULL)
    {
      normal_free(state);
      GSL_ERROR_NULL("failed to allocate ATA matrix", GSL_ENOMEM);
    }

  state->work_ATA = gsl_matrix_alloc(p, p);
  if (state->work_ATA == NULL)
    {
      normal_free(state);
      GSL_ERROR_NULL("failed to allocate temporary ATA matrix", GSL_ENOMEM);
    }

  state->ATb = gsl_vector_alloc(p);
  if (state->ATb == NULL)
    {
      normal_free(state);
      GSL_ERROR_NULL("failed to allocate ATb vector", GSL_ENOMEM);
    }

  return state;
}

static void
normal_free(void *vstate)
{
  normal_state_t *state = (normal_state_t *) vstate;

  if (state->ATA)
    gsl_matrix_free(state->ATA);

  if (state->work_ATA)
    gsl_matrix_free(state->work_ATA);

  if (state->ATb)
    gsl_vector_free(state->ATb);

  free(state);
}

static int
normal_reset(void *vstate)
{
  normal_state_t *state = (normal_state_t *) vstate;

  gsl_matrix_set_zero(state->ATA);
  gsl_vector_set_zero(state->ATb);
  state->bTb = 0.0;

  return GSL_SUCCESS;
}

/*
normal_accumulate()
  Add a new block of rows to the normal equations system

Inputs: A      - new block of rows, n-by-p
        b      - new rhs vector n-by-1
        vstate - workspace

Return: success/error
*/

static int
normal_accumulate(const gsl_matrix * A, const gsl_vector * b,
                  void * vstate)
{
  normal_state_t *state = (normal_state_t *) vstate;
  const size_t n = A->size1;

  if (A->size2 != state->p)
    {
      GSL_ERROR("columns of A do not match workspace", GSL_EBADLEN);
    }
  else if (n != b->size)
    {
      GSL_ERROR("A and b have different numbers of rows", GSL_EBADLEN);
    }
  else
    {
      int s;
      double bnorm;

      /* ATA += A^T A, using only the lower half of the matrix */
      s = gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, A, 1.0, state->ATA);
      if (s)
        return s;

      /* ATb += A^T b */
      s = gsl_blas_dgemv(CblasTrans, 1.0, A, b, 1.0, state->ATb);
      if (s)
        return s;

      /* bTb += b^T b */
      bnorm = gsl_blas_dnrm2(b);
      state->bTb += bnorm * bnorm;

      return GSL_SUCCESS;
    }
}

/*
normal_solve()
  Solve normal equations system:

(A^T A + \lambda^2 I) x = A^T b

using Cholesky decomposition

Inputs: lambda - regularization parameter
        x      - (output) solution vector p-by-1
        rnorm  - (output) residual norm ||b - A x||
        snorm  - (output) solution norm ||x||
        vstate - workspace

Return: success/error
*/

static int
normal_solve(const double lambda, gsl_vector * x,
             double * rnorm, double * snorm,
             void * vstate)
{
  normal_state_t *state = (normal_state_t *) vstate;

  if (x->size != state->p)
    {
      GSL_ERROR("solution vector does not match workspace", GSL_EBADLEN);
    }
  else
    {
      int status;
      gsl_vector_view d = gsl_matrix_diagonal(state->work_ATA);

      /* copy ATA matrix to temporary workspace and regularize */
      gsl_matrix_memcpy(state->work_ATA, state->ATA);
      gsl_vector_add_constant(&d.vector, lambda * lambda);

      /* compute Cholesky decomposition of A^T A */
      status = gsl_linalg_cholesky_decomp(state->work_ATA);
      if (status)
        return status;

      /* solve system (A^T A) x = A^T b */
      status = gsl_linalg_cholesky_solve(state->work_ATA, state->ATb, x);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}

static const gsl_multilarge_linear_type normal_type =
{
  "normal",
  normal_alloc,
  normal_reset,
  normal_accumulate,
  normal_solve,
  normal_free
};

const gsl_multilarge_linear_type * gsl_multilarge_linear_normal =
  &normal_type;
