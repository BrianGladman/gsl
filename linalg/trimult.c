/* linalg/trimult.c
 *
 * Copyright (C) 2019 Patrick Alken
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * This module contains code to compute L^T L where L is a lower triangular matrix
 */

#include <config.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

static int triangular_multiply_L2(CBLAS_UPLO_t Uplo, gsl_matrix * T);
static int triangular_multiply_L3(CBLAS_UPLO_t Uplo, gsl_matrix * T);

#define CROSSOVER_TRIMULT       24

int
gsl_linalg_tri_LTL(gsl_matrix * L)
{
  return triangular_multiply_L3(CblasLower, L);
}

/*
triangular_multiply_L2()
  Compute L^T L or U U^T

Inputs: Uplo - CblasUpper or CblasLower
        T    - on output the upper (or lower) part of T
               is replaced by L^T L or U U^T

Return: success/error

Notes:
1) Based on LAPACK routine DLAUU2 using Level 2 BLAS
*/

static int
triangular_multiply_L2(CBLAS_UPLO_t Uplo, gsl_matrix * T)
{
  const size_t N = T->size1;

  if (N != T->size2)
    {
      GSL_ERROR ("matrix must be square", GSL_ENOTSQR);
    }
  else
    {
      gsl_vector_view v1, v2;
      size_t i;

      if (Uplo == CblasUpper)
        {
        }
      else
        {
          for (i = 0; i < N; ++i)
            {
              double Tii = gsl_matrix_get(T, i, i);

              if (i < N - 1)
                {
                  double tmp;

                  v1 = gsl_matrix_subcolumn(T, i, i, N - i);
                  gsl_blas_ddot(&v1.vector, &v1.vector, &tmp);
                  gsl_matrix_set(T, i, i, tmp);

                  if (i > 0)
                    {
                      gsl_matrix_view m = gsl_matrix_submatrix(T, i + 1, 0, N - i - 1, i);

                      v1 = gsl_matrix_subcolumn(T, i, i + 1, N - i - 1);
                      v2 = gsl_matrix_subrow(T, i, 0, i);

                      gsl_blas_dgemv(CblasTrans, 1.0, &m.matrix, &v1.vector, Tii, &v2.vector);
                    }
                }
              else
                {
                  v1 = gsl_matrix_row(T, N - 1);
                  gsl_blas_dscal(Tii, &v1.vector);
                }
            }
        }

      return GSL_SUCCESS;
    }
}

/*
triangular_multiply_L3()
  Compute L^T L or U U^T

Inputs: Uplo - CblasUpper or CblasLower
        T    - on output the upper (or lower) part of T
               is replaced by L^T L or U U^T

Return: success/error

Notes:
1) Based on ReLAPACK routine DLAUUM using Level 3 BLAS
*/

static int
triangular_multiply_L3(CBLAS_UPLO_t Uplo, gsl_matrix * T)
{
  const size_t N = T->size1;

  if (N != T->size2)
    {
      GSL_ERROR ("matrix must be square", GSL_ENOTSQR);
    }
  else if (N <= CROSSOVER_TRIMULT)
    {
      return triangular_multiply_L2(Uplo, T);
    }
  else
    {
      /* partition matrix:
       *
       * T11 T12
       * T21 T22
       *
       * where T11 is N1-by-N1
       */
      int status;
      const size_t N1 = GSL_LINALG_SPLIT(N);
      const size_t N2 = N - N1;
      gsl_matrix_view T11 = gsl_matrix_submatrix(T, 0, 0, N1, N1);
      gsl_matrix_view T12 = gsl_matrix_submatrix(T, 0, N1, N1, N2);
      gsl_matrix_view T21 = gsl_matrix_submatrix(T, N1, 0, N2, N1);
      gsl_matrix_view T22 = gsl_matrix_submatrix(T, N1, N1, N2, N2);

      /* recursion on T11 */
      status = triangular_multiply_L3(Uplo, &T11.matrix);
      if (status)
        return status;

      if (Uplo == CblasLower)
        {
          /* T11 += T21^T T21 */
          gsl_blas_dsyrk(Uplo, CblasTrans, 1.0, &T21.matrix, 1.0, &T11.matrix);

          /* T21 = T22^T * T21 */
          gsl_blas_dtrmm(CblasLeft, Uplo, CblasTrans, CblasNonUnit, 1.0, &T22.matrix, &T21.matrix);
        }
      else
        {
          /* T11 += T12 T12^T */
          gsl_blas_dsyrk(Uplo, CblasNoTrans, 1.0, &T12.matrix, 1.0, &T11.matrix);

          /* T12 = T12 * T22^T */
          gsl_blas_dtrmm(CblasRight, Uplo, CblasTrans, CblasNonUnit, 1.0, &T22.matrix, &T12.matrix);
        }

      /* recursion on T22 */
      status = triangular_multiply_L3(Uplo, &T22.matrix);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}
