/* L D L^T Decomposition
 *
 * Copyright (C) 2016 Patrick Alken
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 3, or (at your option) any
 * later version.
 *
 * This source is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 */

/*
 * L D L^T decomposition of a symmetrix positive definite matrix.
 *
 * This algorithm does:
 *   P A P' = L D L'
 * with
 *   L  := unit lower left triangle matrix
 *   D  := diagonal matrix
 *   L' := the transposed form of L.
 *   P  := permutation matrix
 *
 */

#include <config.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>

static int pcholesky_swap_rowcol(gsl_matrix * m, const size_t i, const size_t j);

/*
gsl_linalg_pcholesky_decomp()
  Perform Pivoted Cholesky LDLT decomposition of a symmetric positive
semidefinite matrix

Inputs: A - (input) symmetric, positive semidefinite matrix,
                    stored in lower triangle
            (output) lower triangle contains L; diagonal contains D

Return: success/error

Notes:
1) Based on algorithm 4.2.2 (Outer Product LDLT with Pivoting) of
Golub and Van Loan, Matrix Computations (4th ed).
*/

int
gsl_linalg_pcholesky_decomp (gsl_matrix * A, gsl_permutation * p)
{
  const size_t N = A->size1;

  if (N != A->size2)
    {
      GSL_ERROR("LDLT decomposition requires square matrix", GSL_ENOTSQR);
    }
  else if (p->size != N)
    {
      GSL_ERROR ("permutation length must match matrix size", GSL_EBADLEN);
    }
  else
    {
      gsl_vector_view diag = gsl_matrix_diagonal(A);
      size_t k;

      gsl_permutation_init(p);

      for (k = 0; k < N; ++k)
        {
          gsl_vector_view w;
          size_t j;

          /* compute j = max_idx { A_kk, ..., A_nn } */
          w = gsl_vector_subvector(&diag.vector, k, N - k);
          j = gsl_vector_max_index(&w.vector) + k;
          gsl_permutation_swap(p, k, j);

          pcholesky_swap_rowcol(A, k, j);

          if (k < N - 1)
            {
              double alpha = gsl_matrix_get(A, k, k);
              double alphainv = 1.0 / alpha;

              /* v = A(k+1:n, k) */
              gsl_vector_view v = gsl_matrix_subcolumn(A, k, k + 1, N - k - 1);

              /* m = A(k+1:n, k+1:n) */
              gsl_matrix_view m = gsl_matrix_submatrix(A, k + 1, k + 1, N - k - 1, N - k - 1);

              /* m = m - v v^T / alpha */
              gsl_blas_dsyr(CblasLower, -alphainv, &v.vector, &m.matrix);

              /* v = v / alpha */
              gsl_vector_scale(&v.vector, alphainv);
            }
        }

      /* copy the transposed (strictly) lower triangle to the upper triangle */
      gsl_matrix_transpose_tricpy('L', 0, A, A);
      
      return GSL_SUCCESS;
    }
}

int
gsl_linalg_pcholesky_solve(const gsl_matrix * LDLT,
                           const gsl_permutation * p,
                           const gsl_vector * b,
                           gsl_vector * x)
{
  if (LDLT->size1 != LDLT->size2)
    {
      GSL_ERROR ("LDLT matrix must be square", GSL_ENOTSQR);
    }
  else if (LDLT->size1 != p->size)
    {
      GSL_ERROR ("matrix size must match permutation size", GSL_EBADLEN);
    }
  else if (LDLT->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LDLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      int status;

      gsl_vector_memcpy (x, b);

      status = gsl_linalg_pcholesky_svx (LDLT, p, x);
      
      return status;
    }
}

int
gsl_linalg_pcholesky_svx(const gsl_matrix * LDLT,
                         const gsl_permutation * p,
                         gsl_vector * x)
{
  if (LDLT->size1 != LDLT->size2)
    {
      GSL_ERROR ("LDLT matrix must be square", GSL_ENOTSQR);
    }
  else if (LDLT->size1 != p->size)
    {
      GSL_ERROR ("matrix size must match permutation size", GSL_EBADLEN);
    }
  else if (LDLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      gsl_vector_const_view D = gsl_matrix_const_diagonal(LDLT);

      /* x := P b */
      gsl_permute_vector(p, x);

      /* solve: L w = P b */
      gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasUnit, LDLT, x);

      /* solve: D y = w */
      gsl_vector_div(x, &D.vector);

      /* solve: L^T z = y */
      gsl_blas_dtrsv(CblasLower, CblasTrans, CblasUnit, LDLT, x);

      /* compute: x = P^T z */
      gsl_permute_vector_inverse(p, x);

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_pcholesky_decomp2(gsl_matrix * A, gsl_permutation * p,
                             gsl_vector * S)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (M != N)
    {
      GSL_ERROR("cholesky decomposition requires square matrix", GSL_ENOTSQR);
    }
  else if (N != p->size)
    {
      GSL_ERROR ("matrix size must match permutation size", GSL_EBADLEN);
    }
  else if (N != S->size)
    {
      GSL_ERROR("S must have length N", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* compute scaling factors to reduce cond(A) */
      status = gsl_linalg_cholesky_scale(A, S);
      if (status)
        return status;

      /* apply scaling factors */
      status = gsl_linalg_cholesky_scale_apply(A, S);
      if (status)
        return status;

      /* compute Cholesky decomposition of diag(S) A diag(S) */
      status = gsl_linalg_pcholesky_decomp(A, p);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_pcholesky_solve2(const gsl_matrix * LDLT,
                            const gsl_permutation * p,
                            const gsl_vector * S,
                            const gsl_vector * b,
                            gsl_vector * x)
{
  if (LDLT->size1 != LDLT->size2)
    {
      GSL_ERROR ("LDLT matrix must be square", GSL_ENOTSQR);
    }
  else if (LDLT->size1 != p->size)
    {
      GSL_ERROR ("matrix size must match permutation size", GSL_EBADLEN);
    }
  else if (LDLT->size1 != S->size)
    {
      GSL_ERROR ("matrix size must match S", GSL_EBADLEN);
    }
  else if (LDLT->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LDLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      int status;

      gsl_vector_memcpy (x, b);

      status = gsl_linalg_pcholesky_svx2 (LDLT, p, S, x);
      
      return status;
    }
}

int
gsl_linalg_pcholesky_svx2(const gsl_matrix * LDLT,
                          const gsl_permutation * p,
                          const gsl_vector * S,
                          gsl_vector * x)
{
  if (LDLT->size1 != LDLT->size2)
    {
      GSL_ERROR ("LDLT matrix must be square", GSL_ENOTSQR);
    }
  else if (LDLT->size1 != p->size)
    {
      GSL_ERROR ("matrix size must match permutation size", GSL_EBADLEN);
    }
  else if (LDLT->size1 != S->size)
    {
      GSL_ERROR ("matrix size must match S", GSL_EBADLEN);
    }
  else if (LDLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* x := S b */
      gsl_vector_mul(x, S);

      /* solve: A~ x~ = b~, with A~ = S A S, b~ = S b */
      status = gsl_linalg_pcholesky_svx(LDLT, p, x);
      if (status)
        return status;

      /* compute: x = S x~ */
      gsl_vector_mul(x, S);

      return GSL_SUCCESS;
    }
}

/*
pcholesky_swap_rowcol()
  Swap rows and columns i and j of matrix m, assuming only the
lower triangle of m is up to date.
*/

static int
pcholesky_swap_rowcol(gsl_matrix * m, const size_t i, const size_t j)
{
  if (i != j)
    {
      const size_t N = m->size1;
      size_t k;

      /* fill in column i above diagonal using lower triangle */
      for (k = 0; k < i; ++k)
        {
          double mki = gsl_matrix_get(m, i, k);
          gsl_matrix_set(m, k, i, mki);
        }

      /* fill in row i above diagonal using lower triangle */
      for (k = i + 1; k < N; ++k)
        {
          double mik = gsl_matrix_get(m, k, i);
          gsl_matrix_set(m, i, k, mik);
        }

      /* fill in column j above diagonal using lower triangle */
      for (k = 0; k < j; ++k)
        {
          double mkj = gsl_matrix_get(m, j, k);
          gsl_matrix_set(m, k, j, mkj);
        }

      /* fill in row j above diagonal using lower triangle */
      for (k = j + 1; k < N; ++k)
        {
          double mjk = gsl_matrix_get(m, k, j);
          gsl_matrix_set(m, j, k, mjk);
        }

      gsl_matrix_swap_rows(m, i, j);
      gsl_matrix_swap_columns(m, i, j);
    }

  return GSL_SUCCESS;
}
