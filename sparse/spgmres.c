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
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sparse.h>
#include <gsl/gsl_multifit.h>

#include "../linalg/givens.c"
#include "../linalg/apply_givens.c"

gsl_splinalg_gmres_workspace *
gsl_splinalg_gmres_alloc(const size_t n)
{
  gsl_splinalg_gmres_workspace *w;

  w = calloc(1, sizeof(gsl_splinalg_gmres_workspace));
  if (!w)
    {
      GSL_ERROR_NULL("failed to allocate gmres workspace", GSL_ENOMEM);
    }

  w->max_iter = GSL_MIN(n, 10);
  w->n = n;
  w->m = w->max_iter;

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
  const double tol = 1.0e-6;

  /* initial guess x = 0 */
  gsl_vector_set_zero(x);

  s = gsl_splinalg_gmres_solve_x(A, b, tol, x, w);

  return s;
} /* gsl_splinalg_gmres_solve() */

/*
gsl_splinalg_gmres_solve_x()
  Solve A*x = b using GMRES algorithm

Inputs: A   - sparse square matrix
        b   - right hand side vector
        tol - stopping tolerance
        x   - (input/output) on input, initial estimate x_0;
              on output, solution vector
        w   - workspace

Notes:
1) Based on algorithm 6.10 of

Saad, Iterative methods for sparse linear systems (2nd edition), SIAM
*/

int
gsl_splinalg_gmres_solve_x(const gsl_spmatrix *A, const gsl_vector *b,
                           const double tol, gsl_vector *x,
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
      const size_t maxit = M;
      const double normb = gsl_blas_dnrm2(b); /* ||b|| */
      const double reltol = tol * normb;      /* tol*||b|| */
      size_t m, k;
      double beta;          /* || r_0 || = || b - A*x_0 || */
      double tau;           /* householder scalar */
      gsl_matrix *H;        /* Hessenberg matrix [h_1 h_2 ... h_m] */
      gsl_vector *tauv = gsl_vector_alloc(maxit + 1);
      gsl_vector *y = gsl_vector_alloc(maxit); /* y_m */
      gsl_vector *w = gsl_vector_calloc(maxit + 1); /* w */
      gsl_matrix_view Rm;   /* R_m = [h_1 h_2 ... h_m] */
      gsl_vector_view ym, wm;

      double *cc = malloc(maxit * sizeof(double));
      double *ss = malloc(maxit * sizeof(double));

      gsl_vector *z = gsl_vector_alloc(N);

      /*
       * The Hessenberg matrix will have the following structure:
       *
       * H = [ ||r_0|| | v_1 v_2 ... v_m     ]
       *     [   u_1   | u_2 u_3 ... u_{m+1} ]
       *
       * where v_j are the orthonormal vectors spanning the Krylov
       * subpsace of length j + 1 and u_{j+1} are the householder
       * vectors of length n - j - 1.
       * In fact, u_{j+1} has length n - j since u_{j+1}[0] = 1,
       * but this 1 is not stored.
       */
      H = gsl_matrix_calloc(N, maxit + 1);

      /* Step 1a: compute z = b - A*x_0 */
      gsl_vector_memcpy(z, b);
      gsl_spblas_dgemv(-1.0, A, x, 1.0, z);

      /* Step 1b */
      {
        /* hessenberg vector h_0 */
        gsl_vector_view h0 = gsl_matrix_column(H, 0);

        gsl_vector_memcpy(&h0.vector, z);
        tau = gsl_linalg_householder_transform(&h0.vector);

        /* store tau_1 */
        gsl_vector_set(tauv, 0, tau);

        /* initialize w */
        gsl_vector_set(w, 0, gsl_vector_get(&h0.vector, 0));
      }

      for (m = 1; m <= maxit; ++m)
        {
          size_t j = m - 1; /* C indexing */
          double c, s;

          /* v_m */
          gsl_vector_view vm = gsl_matrix_column(H, m);

          /* v_m(m:end) */
          gsl_vector_view vv = gsl_vector_subvector(&vm.vector, j, N - j);

          /* householder vector u_m for projection P_m */
          gsl_vector_view um = gsl_matrix_subcolumn(H, j, j, N - j);

          /* householder vector u_{m+1} for projection P_{m+1} */
          gsl_vector_view ump1 = gsl_matrix_subcolumn(H, m, m, N - m);

          /* Step 2a: form v_m = P_m e_m = e_m - tau_m w_m */
          gsl_vector_set_zero(&vm.vector);
          gsl_vector_memcpy(&vv.vector, &um.vector);
          tau = gsl_vector_get(tauv, j); /* tau_m */
          gsl_vector_scale(&vv.vector, -tau);
          gsl_vector_set(&vv.vector, 0, 1.0 - tau);

          /* Step 2a: v_m <- P_1 P_2 ... P_{m-1} v_m */
          for (k = j; k > 0 && k--; )
            {
              gsl_vector_view uk =
                gsl_matrix_subcolumn(H, k, k, N - k);
              gsl_vector_view vk =
                gsl_vector_subvector(&vm.vector, k, N - k);
              tau = gsl_vector_get(tauv, k);
              gsl_linalg_householder_hv(tau, &uk.vector, &vk.vector);
            }

          /* Step 2a: v_m <- A*v_m */
          gsl_spblas_dgemv(1.0, A, &vm.vector, 0.0, z);
          gsl_vector_memcpy(&vm.vector, z);

          /* Step 2a: v_m <- P_m ... P_1 v_m */
          for (k = 0; k <= j; ++k)
            {
              gsl_vector_view uk = gsl_matrix_subcolumn(H, k, k, N - k);
              gsl_vector_view vk = gsl_vector_subvector(&vm.vector, k, N - k);
              tau = gsl_vector_get(tauv, k);
              gsl_linalg_householder_hv(tau, &uk.vector, &vk.vector);
            }

          /* Steps 2c,2d: find P_{m+1} and set v_m <- P_{m+1} v_m */
          tau = gsl_linalg_householder_transform(&ump1.vector);
          gsl_vector_set(tauv, j + 1, tau);

          /* Step 2e: v_m <- J_{m-1} ... J_1 v_m */
          for (k = 0; k < j; ++k)
            apply_givens_vec(&vm.vector, k, k + 1, cc[k], ss[k]);

          /* Step 2g: find givens rotation J_m for v_m(m:m+1) */
          create_givens(gsl_vector_get(&vm.vector, j),
                        gsl_vector_get(&vm.vector, j + 1),
                        &c, &s);

          /* store givens rotation for later use */
          cc[j] = c;
          ss[j] = s;

          /* Step 2h: v_m <- J_m v_m */
          apply_givens_vec(&vm.vector, j, j + 1, c, s);

          /* Step 2h: w <- J_m w */
          apply_givens_vec(w, j, j + 1, c, s);

          /*
           * Step 2i: R_m = [ R_{m-1}, v_m ] - already taken care
           * of due to our memory storage scheme
           */
        }

      /* Step 3a: solve triangular system R_m y_m = w */
      m--;
      Rm = gsl_matrix_submatrix(H, 0, 1, m, m);
      wm = gsl_vector_subvector(w, 0, m);
      ym = gsl_vector_subvector(y, 0, m);
      gsl_vector_memcpy(&ym.vector, &wm.vector);
      gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit,
                     &Rm.matrix, &ym.vector);

      /*
       * Step 3b: update solution vector x; the loop below
       * uses a different but equivalent formulation from
       * Saad, algorithm 6.10, step 14
       */
      gsl_vector_set_zero(z);
      for (k = m; k > 0 && k--; )
        {
          double ymk = gsl_vector_get(&ym.vector, k);
          gsl_vector_view uk = gsl_matrix_subcolumn(H, k, k, N - k);
          gsl_vector_view zk = gsl_vector_subvector(z, k, N - k);

          /* z <- n_k e_k + z */
          gsl_vector_set(z, k, gsl_vector_get(z, k) + ymk);

          /* z <- P_k z */
          tau = gsl_vector_get(tauv, k);
          gsl_linalg_householder_hv(tau, &uk.vector, &zk.vector);
        }

      /* x <- x + z */
      gsl_vector_add(x, z);

      gsl_vector_fprintf(stdout, x, "%.12e");
      exit(1);

#if 0
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
#endif

      return GSL_SUCCESS;
    }
} /* gsl_splinalg_gmres_solve() */
