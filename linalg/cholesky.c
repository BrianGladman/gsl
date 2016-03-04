/* Cholesky Decomposition
 *
 * Copyright (C) 2000 Thomas Walter
 * Copyright (C) 2000, 2001, 2002, 2003, 2005, 2007 Brian Gough, Gerard Jungman
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
 *
 * 03 May 2000: Modified for GSL by Brian Gough
 * 29 Jul 2005: Additions by Gerard Jungman
 * 04 Mar 2016: Change Cholesky algorithm to gaxpy version by Patrick Alken
 */

/*
 * Cholesky decomposition of a symmetrix positive definite matrix.
 * This is useful to solve the matrix arising in
 *    periodic cubic splines
 *    approximating splines
 *
 * This algorithm does:
 *   A = L * L'
 * with
 *   L  := lower left triangle matrix
 *   L' := the transposed form of L.
 *
 */

#include <config.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

/*
gsl_linalg_cholesky_decomp()
  Perform Cholesky decomposition of a symmetric positive
definite matrix

Inputs: A - (input) symmetric, positive definite matrix
            (output) lower triangle contains Cholesky factor

Return: success/error

Notes:
1) Based on algorithm 4.2.1 (Gaxpy Cholesky) of Golub and
Van Loan, Matrix Computations (4th ed).
*/

int
gsl_linalg_cholesky_decomp (gsl_matrix * A)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (M != N)
    {
      GSL_ERROR("cholesky decomposition requires square matrix", GSL_ENOTSQR);
    }
  else
    {
      size_t j;

      for (j = 0; j < N; ++j)
        {
          double ajj;
          gsl_vector_view v = gsl_matrix_subcolumn(A, j, j, N - j); /* A(j:n,j) */

          if (j > 0)
            {
              gsl_vector_view w = gsl_matrix_subrow(A, j, 0, j);           /* A(j,1:j-1)^T */
              gsl_matrix_view m = gsl_matrix_submatrix(A, j, 0, N - j, j); /* A(j:n,1:j-1) */

              gsl_blas_dgemv(CblasNoTrans, -1.0, &m.matrix, &w.vector, 1.0, &v.vector);
            }

          ajj = gsl_matrix_get(A, j, j);

          if (ajj <= 0.0)
            {
              GSL_ERROR("matrix is not positive definite", GSL_EDOM);
            }

          ajj = sqrt(ajj);
          gsl_vector_scale(&v.vector, 1.0 / ajj);
        }

      /* now copy the transposed lower triangle to the upper triangle,
       * the diagonal is common. */
      gsl_matrix_transpose_tricpy('L', 0, A, A);
      
      return GSL_SUCCESS;
    }
}

int
gsl_linalg_cholesky_solve (const gsl_matrix * LLT,
                           const gsl_vector * b,
                           gsl_vector * x)
{
  if (LLT->size1 != LLT->size2)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else if (LLT->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      /* Copy x <- b */

      gsl_vector_memcpy (x, b);

      /* Solve for c using forward-substitution, L c = b */

      gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasNonUnit, LLT, x);

      /* Perform back-substitution, U x = c */

      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, LLT, x);


      return GSL_SUCCESS;
    }
}

int
gsl_linalg_cholesky_svx (const gsl_matrix * LLT,
                         gsl_vector * x)
{
  if (LLT->size1 != LLT->size2)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else if (LLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      /* Solve for c using forward-substitution, L c = b */

      gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasNonUnit, LLT, x);

      /* Perform back-substitution, U x = c */

      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, LLT, x);

      return GSL_SUCCESS;
    }
}

/*
gsl_linalg_cholesky_invert()
  Compute the inverse of a symmetric positive definite matrix in
Cholesky form.

Inputs: LLT - matrix in cholesky form on input
              A^{-1} = L^{-t} L^{-1} on output

Return: success or error
*/

int
gsl_linalg_cholesky_invert(gsl_matrix * LLT)
{
  if (LLT->size1 != LLT->size2)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else
    {
      size_t N = LLT->size1;
      size_t i, j;
      double sum;
      gsl_vector_view v1, v2;

      /* invert the lower triangle of LLT */
      for (i = 0; i < N; ++i)
        {
          double ajj;

          j = N - i - 1;

          gsl_matrix_set(LLT, j, j, 1.0 / gsl_matrix_get(LLT, j, j));
          ajj = -gsl_matrix_get(LLT, j, j);

          if (j < N - 1)
            {
              gsl_matrix_view m;
              
              m = gsl_matrix_submatrix(LLT, j + 1, j + 1,
                                       N - j - 1, N - j - 1);
              v1 = gsl_matrix_subcolumn(LLT, j, j + 1, N - j - 1);

              gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit,
                             &m.matrix, &v1.vector);

              gsl_blas_dscal(ajj, &v1.vector);
            }
        } /* for (i = 0; i < N; ++i) */

      /*
       * The lower triangle of LLT now contains L^{-1}. Now compute
       * A^{-1} = L^{-t} L^{-1}
       *
       * The (ij) element of A^{-1} is column i of L^{-1} dotted into
       * column j of L^{-1}
       */

      for (i = 0; i < N; ++i)
        {
          for (j = i + 1; j < N; ++j)
            {
              v1 = gsl_matrix_subcolumn(LLT, i, j, N - j);
              v2 = gsl_matrix_subcolumn(LLT, j, j, N - j);

              /* compute Ainv_{ij} = sum_k Linv_{ki} Linv_{kj} */
              gsl_blas_ddot(&v1.vector, &v2.vector, &sum);

              /* store in upper triangle */
              gsl_matrix_set(LLT, i, j, sum);
            }

          /* now compute the diagonal element */
          v1 = gsl_matrix_subcolumn(LLT, i, i, N - i);
          gsl_blas_ddot(&v1.vector, &v1.vector, &sum);
          gsl_matrix_set(LLT, i, i, sum);
        }

      /* copy the transposed upper triangle to the lower triangle */

      for (j = 1; j < N; j++)
        {
          for (i = 0; i < j; i++)
            {
              double A_ij = gsl_matrix_get (LLT, i, j);
              gsl_matrix_set (LLT, j, i, A_ij);
            }
        } 

      return GSL_SUCCESS;
    }
} /* gsl_linalg_cholesky_invert() */

int
gsl_linalg_cholesky_decomp_unit(gsl_matrix * A, gsl_vector * D)
{
  const size_t N = A->size1;
  size_t i, j;

  /* initial Cholesky */
  int stat_chol = gsl_linalg_cholesky_decomp(A);

  if(stat_chol == GSL_SUCCESS)
  {
    /* calculate D from diagonal part of initial Cholesky */
    for(i = 0; i < N; ++i)
    {
      const double C_ii = gsl_matrix_get(A, i, i);
      gsl_vector_set(D, i, C_ii*C_ii);
    }

    /* multiply initial Cholesky by 1/sqrt(D) on the right */
    for(i = 0; i < N; ++i)
    {
      for(j = 0; j < N; ++j)
      {
        gsl_matrix_set(A, i, j, gsl_matrix_get(A, i, j) / sqrt(gsl_vector_get(D, j)));
      }
    }

    /* Because the initial Cholesky contained both L and transpose(L),
       the result of the multiplication is not symmetric anymore;
       but the lower triangle _is_ correct. Therefore we reflect
       it to the upper triangle and declare victory.
       */
    for(i = 0; i < N; ++i)
      for(j = i + 1; j < N; ++j)
        gsl_matrix_set(A, i, j, gsl_matrix_get(A, j, i));
  }

  return stat_chol;
}

/*
gsl_linalg_cholesky_scale()
  This function computes scale factors diag(D), such that

diag(D) A diag(D)

has a condition number within a factor N of the matrix
with the smallest condition number over all possible
diagonal scalings. See Corollary 7.6 of:

N. J. Higham, Accuracy and Stability of Numerical Algorithms (2nd Edition),
SIAM, 2002.

Inputs: A     - (input/output)
                on input, symmetric positive definite matrix
                on output, diag(D) * A * diag(D)
        D     - (output) scale factors, D_i = 1 / sqrt(A_ii)
*/

int
gsl_linalg_cholesky_scale(gsl_matrix * A, gsl_vector * D)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (M != N)
    {
      GSL_ERROR("A is not a square matrix", GSL_ENOTSQR);
    }
  else if (N != D->size)
    {
      GSL_ERROR("D must have length N", GSL_EBADLEN);
    }
  else
    {
      gsl_vector_view diag = gsl_matrix_diagonal(A);
      double amin, amax;
      size_t i, j;

      /* compute minimum and maximum elements of diag(A) */
      gsl_vector_minmax(&diag.vector, &amin, &amax);

      if (amin <= 0.0)
        {
          GSL_ERROR ("matrix must be positive definite", GSL_EDOM);
        }

      /* compute D_i = 1/sqrt(A_{ii}) */
      for (i = 0; i < N; ++i)
        {
          double Aii = gsl_matrix_get(A, i, i);
          gsl_vector_set(D, i, 1.0 / sqrt(Aii));
        }

      /* compute: A <- diag(D) A diag(D) using lower triangle */
      for (j = 0; j < N; ++j)
        {
          double dj = gsl_vector_get(D, j);

          for (i = j; i < N; ++i)
            {
              double di = gsl_vector_get(D, i);
              double *Aij = gsl_matrix_ptr(A, i, j);
              *Aij *= di * dj;
            }
        }

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_cholesky_decomp2(gsl_matrix * A, gsl_vector * D)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (M != N)
    {
      GSL_ERROR("cholesky decomposition requires square matrix", GSL_ENOTSQR);
    }
  else if (N != D->size)
    {
      GSL_ERROR("D must have length N", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* compute scaling factors to reduce cond(A) */
      status = gsl_linalg_cholesky_scale(A, D);
      if (status)
        return status;

      /* compute Cholesky decomposition of diag(D) A diag(D) */
      status = gsl_linalg_cholesky_decomp(A);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_cholesky_svx2 (const gsl_matrix * LLT,
                          const gsl_vector * D,
                          gsl_vector * x)
{
  if (LLT->size1 != LLT->size2)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else if (LLT->size2 != D->size)
    {
      GSL_ERROR ("matrix size must match D", GSL_EBADLEN);
    }
  else if (LLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      /* b~ = diag(D) b */
      gsl_vector_mul(x, D);

      /* Solve for c using forward-substitution, L c = b~ */
      gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasNonUnit, LLT, x);

      /* Perform back-substitution, U x~ = c */
      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, LLT, x);

      /* compute original solution vector x = D x~ */
      gsl_vector_mul(x, D);

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_cholesky_solve2 (const gsl_matrix * LLT,
                            const gsl_vector * D,
                            const gsl_vector * b,
                            gsl_vector * x)
{
  if (LLT->size1 != LLT->size2)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else if (LLT->size1 != D->size)
    {
      GSL_ERROR ("matrix size must match D size", GSL_EBADLEN);
    }
  else if (LLT->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* Copy x <- b */
      gsl_vector_memcpy (x, b);

      status = gsl_linalg_cholesky_svx2(LLT, D, x);

      return status;
    }
}
