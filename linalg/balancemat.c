/* linalg/balance.c
 * 
 * Copyright (C) 2006 Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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

/* Balance a general matrix by scaling the rows and columns, so the
 * new row and column norms are the same order of magnitude.
 *
 * B =  D^-1 A D
 *
 * where D is a diagonal matrix
 * 
 * This is necessary for the unsymmetric eigenvalue problem since the
 * calculation can become numerically unstable for unbalanced
 * matrices.  
 *
 * See Golub & Van Loan, "Matrix Computations" (3rd ed), Section 7.5.7
 * and Wilkinson & Reinsch, "Handbook for Automatic Computation", II/11 p320.
 */

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_linalg.h>

int
gsl_linalg_balance_matrix (gsl_matrix * A, gsl_vector * D)
{
  double row_norm, col_norm;
  int not_converged;
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (M != N)
    {
      GSL_ERROR ("matrix for balancing must be square", GSL_EINVAL);
    }

  if (D->size != N)
    {
      GSL_ERROR ("length of D must match dimension of A", GSL_EINVAL);
    }

  gsl_vector_set_all (D, 1.0);

  not_converged = 1;

  while (not_converged)
    {
      size_t i, j;
      double g, f, s;

      not_converged = 0;

      for (i = 0; i < N; i++)
        {
          row_norm = 0.0;
          col_norm = 0.0;

          for (j = 0; j < N; ++j)
            {
              if (j != i)
                {
                  col_norm += fabs (gsl_matrix_get (A, j, i));
                  row_norm += fabs (gsl_matrix_get (A, i, j));
                }
            }

          if ((col_norm == 0.0) || (row_norm == 0.0) 
              || !gsl_finite (col_norm) || !gsl_finite (row_norm))
            {
              continue;
            }

          g = row_norm / 2.0;
          f = 1.0;
          s = col_norm + row_norm;

          /* FIXME: we could use frexp() here */

          /*
           * find the integer power of the machine radix which
           * comes closest to balancing the matrix
           */
          while (col_norm < g)
            {
              f *= 2.0;
              col_norm *= 4.0;
            }

          g = row_norm * 2.0;

          while (col_norm > g)
            {
              f /= 2.0;
              col_norm /= 4.0;
            }

          if ((row_norm + col_norm) < 0.95 * s * f)
            {
              not_converged = 1;

              gsl_vector_set (D, i, f);

              /* apply similarity transformation */

              if (f != 1.0)
                {
                  gsl_vector_view A_row_i = gsl_matrix_row (A, i);
                  gsl_vector_view A_col_i = gsl_matrix_column (A, i);
                  g = 1.0 / f;
                  gsl_blas_dscal (g, &A_row_i.vector);
                  gsl_blas_dscal (f, &A_col_i.vector);
                }
            }
        }
    }

  return GSL_SUCCESS;
}

/*
gsl_linalg_balance_accum()
  Accumulate a balancing transformation into a matrix.
This is used during the computation of Schur vectors since the
Schur vectors computed are the vectors for the balanced matrix.
We must at some point accumulate the balancing transformation into
the Schur vector matrix to get the vectors for the original matrix.

A -> D A

where D is the diagonal matrix

Inputs: A - matrix to transform
        D - vector containing diagonal elements of D
*/

int
gsl_linalg_balance_accum(gsl_matrix *A, gsl_vector *D)
{
  const size_t N = A->size1;

  if (N != D->size)
    {
      GSL_ERROR ("vector must match matrix size", GSL_EBADLEN);
    }
  else
    {
      size_t i;
      double s;
      gsl_vector_view r;

      for (i = 0; i < N; ++i)
        {
          s = gsl_vector_get(D, i);
          r = gsl_matrix_row(A, i);

          gsl_blas_dscal(s, &r.vector);
        }

      return GSL_SUCCESS;
    }
} /* gsl_linalg_balance_accum() */
