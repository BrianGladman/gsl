/* spgmres.c
 * 
 * Copyright (C) 2014 Patrick Alken
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
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sparse.h>

gsl_splinalg_gmres_workspace *
gsl_splinalg_gmres_alloc(const size_t n)
{
  gsl_splinalg_gmres_workspace *w;

  w = calloc(1, sizeof(gsl_splinalg_gmres_workspace));
  if (!w)
    {
      GSL_ERROR_NULL("failed to allocate gmres workspace", GSL_ENOMEM);
    }

  w->n = n;
  w->m = GSL_MIN(n, 10);
  w->max_iter = 100;

  w->r = gsl_vector_alloc(n);
  if (!w->r)
    {
      gsl_splinalg_gmres_free(w);
      GSL_ERROR_NULL("failed to allocate r vector", GSL_ENOMEM);
    }

  w->V = gsl_matrix_alloc(n, w->m + 1);
  if (!w->V)
    {
      gsl_splinalg_gmres_free(w);
      GSL_ERROR_NULL("failed to allocate V matrix", GSL_ENOMEM);
    }

  w->H = gsl_matrix_alloc(w->m + 1, w->m);
  if (!w->H)
    {
      gsl_splinalg_gmres_free(w);
      GSL_ERROR_NULL("failed to allocate H matrix", GSL_ENOMEM);
    }

  return w;
} /* gsl_splinalg_gmres_alloc() */

void
gsl_splinalg_gmres_free(gsl_splinalg_gmres_workspace *w)
{
  if (w->r)
    gsl_vector_free(w->r);

  if (w->V)
    gsl_matrix_free(w->V);

  if (w->H)
    gsl_matrix_free(w->H);

  free(w);
} /* gsl_splinalg_gmres_free() */

int
gsl_splinalg_gmres_solve(const gsl_spmatrix *A, const gsl_vector *b,
                         gsl_vector *x,
                         gsl_splinalg_gmres_workspace *w)
{
  int s;

  /* initial guess x = 0 */
  gsl_vector_set_zero(x);

  s = gsl_splinalg_gmres_solve_x(A, b, x, w);

  return s;
} /* gsl_splinalg_gmres_solve() */

/*
gsl_splinalg_gmres_solve_x()
  Solve A*x = b using GMRES algorithm

Inputs: A - sparse square matrix
        b - right hand side vector
        x - (input/output) on input, initial estimate x_0;
            on output, solution vector
        w - workspace

Notes:
1) Based on algorithm 6.9 of

Saad, Iterative methods for sparse linear systems (2nd edition), SIAM
*/

int
gsl_splinalg_gmres_solve_x(const gsl_spmatrix *A, const gsl_vector *b,
                           gsl_vector *x,
                           gsl_splinalg_gmres_workspace *w)
{
  const size_t N = A->size1;

  if (N != A->size2)
    {
      GSL_ERROR("matrix must be square", GSL_ENOTSQR);
    }
  else if (N != b->size)
    {
      GSL_ERROR("matrix does not match right hand side", GSL_EBADLEN);
    }
  else if (N != x->size)
    {
      GSL_ERROR("matrix does not match solution vector", GSL_EBADLEN);
    }
  else if (N != w->n)
    {
      GSL_ERROR("matrix does not match workspace", GSL_EBADLEN);
    }
  else
    {
      const size_t M = w->m;
      size_t iter;          /* iteration number */
      size_t i, j;          /* looping */
      double beta;          /* || r_0 || = || b - A*x_0 || */
      gsl_vector_view v1 = gsl_matrix_column(w->V, 0); /* v_1 */
      gsl_matrix *H = w->H; /* Hessenberg matrix */
      double h;             /* element of H */

      for (iter = 0; iter < w->max_iter; ++iter)
        {
          /* Step 1: r = b - A*x, beta = ||r||, v1 = r / ||r|| */
          gsl_vector_memcpy(&v1.vector, b);
          gsl_spblas_dgemv(-1.0, A, x, 1.0, &v1.vector);
          beta = gsl_blas_dnrm2(&v1.vector);
          gsl_vector_scale(&v1.vector, 1.0 / beta);

          fprintf(stderr, "gsl: iter = %zu, residual = %e\n", iter, beta);

          gsl_matrix_set_zero(H);

          /* Step 2 */
          for (j = 0; j < M; ++j)
            {
              gsl_vector_view vj = gsl_matrix_column(w->V, j);
              gsl_vector_view wj = gsl_matrix_column(w->V, j + 1);

              /* Step 3: w_j = A*v_j */
              gsl_spblas_dgemv(1.0, A, &vj.vector, 0.0, &wj.vector);

              /* Step 4 */
              for (i = 0; i < j + 1; ++i)
                {
                  gsl_vector_view vi = gsl_matrix_column(w->V, i);

                  /* Step 5: h_{ij} = (w_j, v_i) */
                  gsl_blas_ddot(&wj.vector, &vi.vector, &h);
                  gsl_matrix_set(H, i, j, h);

                  /* Step 6: w_j -> w_j - h_{ij} v_i */
                  gsl_blas_daxpy(-h, &vi.vector, &wj.vector);
                }

              /* Step 7/8: h_{j+1,j} = || w_j || */
              h = gsl_blas_dnrm2(&wj.vector);
              gsl_matrix_set(H, j + 1, j, h);

              /* Step 9: v_{j+1} = w_j / h_{j+1,j} */
              gsl_vector_scale(&wj.vector, 1.0 / h);
            }
        }

      return GSL_SUCCESS;
    }
} /* gsl_splinalg_gmres_solve() */
