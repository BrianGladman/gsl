/* Modified Cholesky Decomposition
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

#include <config.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>

/*
 * This module contains routines related to the Modified Cholesky
 * Decomposition, which factors a symmetric indefinite matrix A as
 *
 * P (A + E) P^T = L D L^T
 *
 * where:
 * P: permutation matrix
 * E: small, non-negative diagonal matrix
 * L: unit lower triangular matrix
 * D: strictly positive diagonal matrix
 *
 * These routines follow these works closely:
 *
 * [1] P. E. Gill, W. Murray, M. H. Wright, Practical Optimization,
 *     Academic Press, 1981.
 *
 * [2] Dennis and Schnabel, Numerical Methods for Unconstrained Optimization
 *     and Nonlinear Equations, SIAM, 1996
 */

static size_t mcholesky_maxabs(const gsl_vector * v, double *maxabs);
static int mcholesky_swap_rowcol(gsl_matrix * m, const size_t i, const size_t j);

/*
gsl_linalg_mcholesky_decomp()
  Perform Pivoted Modified Cholesky LDLT decomposition of a symmetric positive
indefinite matrix:

P (A + E) P^T = L D L^T

Inputs: A - (input) symmetric, positive indefinite matrix,
                    stored in lower triangle
            (output) lower triangle contains L; diagonal contains D
        p - (output) permutation matrix P
        E - (output) perturbation matrix E

Return: success/error

Notes:
1) Based on algorithm 4.2.2 (Outer Product LDLT with Pivoting) of
Golub and Van Loan, Matrix Computations (4th ed), with modifications
described in [1] and [2]

2) E can be set to NULL if not required
*/

int
gsl_linalg_mcholesky_decomp (gsl_matrix * A, gsl_permutation * p,
                             gsl_vector * E)
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
      const double delta = GSL_DBL_EPSILON;
      double beta;
      double gamma = 0.0;
      double xi = 0.0;
      gsl_vector_view diag = gsl_matrix_diagonal(A);
      size_t i, j;

      gsl_permutation_init(p);

      /* compute:
       * gamma = max | A_{ii} |
       * xi = max_{i \ne j} | A_{ij} |
       */
      for (i = 0; i < N; ++i)
        {
          double aii = gsl_matrix_get(A, i, i);

          gamma = GSL_MAX(gamma, fabs(aii));

          for (j = 0; j < i; ++j)
            {
              double aij = gsl_matrix_get(A, i, j);
              xi = GSL_MAX(xi, fabs(aij));
            }
        }

      /* compute:
       * beta = sqrt[ max { gamma, xi/nu, eps } ]
       * with: nu = max{ sqrt(N^2 - 1), 1 }
       */
      if (N == 1)
        {
          beta = GSL_MAX(GSL_MAX(gamma, xi), GSL_DBL_EPSILON);
        }
      else
        {
          double nu = sqrt(N*N - 1.0);
          beta = GSL_MAX(GSL_MAX(gamma, xi / nu), GSL_DBL_EPSILON);
        }

      beta = sqrt(beta);

      for (j = 0; j < N; ++j)
        {
          double ajj, thetaj, u, alpha, alphainv;
          gsl_vector_view w;
          size_t q;

          /* compute q = max_idx { A_jj, ..., A_nn } */
          w = gsl_vector_subvector(&diag.vector, j, N - j);
          q = mcholesky_maxabs(&w.vector, NULL) + j;

          gsl_permutation_swap(p, q, j);
          mcholesky_swap_rowcol(A, q, j);

          /* theta_j = max_{j+1 <= i <= n} |A_{ij}| */
          if (j < N - 1)
            {
              w = gsl_matrix_subcolumn(A, j, j + 1, N - j - 1);
              mcholesky_maxabs(&w.vector, &thetaj);
            }
          else
            {
              thetaj = 0.0;
            }

          u = thetaj / beta;

          /* compute alpha = d_j */
          ajj = gsl_matrix_get(A, j, j);
          alpha = GSL_MAX(GSL_MAX(delta, fabs(ajj)), u * u);
          alphainv = 1.0 / alpha;

          if (j < N - 1)
            {
              /* v = A(j+1:n, j) */
              gsl_vector_view v = gsl_matrix_subcolumn(A, j, j + 1, N - j - 1);

              /* m = A(j+1:n, j+1:n) */
              gsl_matrix_view m = gsl_matrix_submatrix(A, j + 1, j + 1, N - j - 1, N - j - 1);

              /* m = m - v v^T / alpha */
              gsl_blas_dsyr(CblasLower, -alphainv, &v.vector, &m.matrix);

              /* v = v / alpha */
              gsl_vector_scale(&v.vector, alphainv);

            }

          if (E)
            gsl_vector_set(E, j, alpha - ajj);

          gsl_matrix_set(A, j, j, alpha);
        }

      if (E)
        {
          /* we currently have: P A P^T + E = L D L^T, permute E
           * so that we have: P (A + E) P^T = L D L^T */
          gsl_permute_vector_inverse(p, E);
        }

      /* copy the transposed (strictly) lower triangle to the upper triangle */
      gsl_matrix_transpose_tricpy('L', 0, A, A);
      
      return GSL_SUCCESS;
    }
}

int
gsl_linalg_mcholesky_solve(const gsl_matrix * LDLT,
                           const gsl_permutation * p,
                           const gsl_vector * b,
                           gsl_vector * x)
{
  int status = gsl_linalg_pcholesky_solve(LDLT, p, b, x);
  return status;
}

int
gsl_linalg_mcholesky_svx(const gsl_matrix * LDLT,
                         const gsl_permutation * p,
                         gsl_vector * x)
{
  int status = gsl_linalg_pcholesky_svx(LDLT, p, x);
  return status;
}

int
gsl_linalg_mcholesky_decomp2(gsl_matrix * A, gsl_permutation * p,
                             gsl_vector * E, gsl_vector * S)
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
      status = gsl_linalg_mcholesky_decomp(A, p, E);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_mcholesky_solve2(const gsl_matrix * LDLT,
                            const gsl_permutation * p,
                            const gsl_vector * S,
                            const gsl_vector * b,
                            gsl_vector * x)
{
  int status = gsl_linalg_pcholesky_solve2(LDLT, p, S, b, x);
  return status;
}

/*
mcholesky_maxabs()
  Compute:

val = max_i |v_i|

Inputs: v      - vector
        maxabs - (output) max abs value

Return: index corresponding to max_i |v_i|
*/

static size_t
mcholesky_maxabs(const gsl_vector * v, double *maxabs)
{
  const size_t n = v->size;
  size_t i;
  size_t idx = 0;
  double max = gsl_vector_get(v, idx);

  for (i = 1; i < n; ++i)
    {
      double vi = gsl_vector_get(v, i);
      double absvi = fabs(vi);

      if (absvi > max)
        {
          max = absvi;
          idx = i;
        }
    }

  if (maxabs)
    *maxabs = max;

  return idx;
}

/*
mcholesky_swap_rowcol()
  Swap rows and columns i and j of matrix m, assuming only the
lower triangle of m is up to date.
*/

static int
mcholesky_swap_rowcol(gsl_matrix * m, const size_t i, const size_t j)
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
