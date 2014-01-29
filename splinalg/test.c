/* test.c
 * 
 * Copyright (C) 2012-2014 Patrick Alken
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
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_splinalg.h>

/*
create_random_sparse()
  Create a random sparse matrix with approximately
M*N*density non-zero entries

Inputs: M       - number of rows
        N       - number of columns
        density - sparse density \in [0,1]
                  0 = no non-zero entries
                  1 = all m*n entries are filled
        r       - random number generator

Return: pointer to sparse matrix in triplet format (must be freed by caller)

Notes:
1) non-zero matrix entries are uniformly distributed in [0,1]
*/

static gsl_spmatrix *
create_random_sparse(const size_t M, const size_t N, const double density,
                     const gsl_rng *r)
{
  gsl_spmatrix *m = gsl_spmatrix_alloc(M, N);
  size_t nnzwanted = (size_t) round(M * N * GSL_MIN(density, 1.0));
  size_t n = 0;

  while (n <= nnzwanted)
    {
      /* generate a random row and column */
      size_t i = gsl_rng_uniform(r) * M;
      size_t j = gsl_rng_uniform(r) * N;
      double x;

      /* check if this position is already filled */
      if (gsl_spmatrix_get(m, i, j) != 0.0)
        continue;

      /* generate random m_{ij} and add it */
      x = gsl_rng_uniform(r);
      gsl_spmatrix_set(m, i, j, x);
      ++n;
    }

  return m;
} /* create_random_sparse() */

static void
create_random_vector(gsl_vector *v, const gsl_rng *r)
{
  size_t i;

  for (i = 0; i < v->size; ++i)
    {
      double x = gsl_rng_uniform(r);
      gsl_vector_set(v, i, x);
    }
} /* create_random_vector() */

int
test_vectors(gsl_vector *observed, gsl_vector *expected, const double tol,
             const char *str)
{
  int s = 0;
  size_t N = observed->size;
  size_t i;

  for (i = 0; i < N; ++i)
    {
      double x_obs = gsl_vector_get(observed, i);
      double x_exp = gsl_vector_get(expected, i);

      gsl_test_rel(x_obs, x_exp, tol, "N=%zu i=%zu %s", N, i, str);
    }

  return s;
} /* test_vectors() */

/*
test_toeplitz()
  Solve u''(x) = -pi^2 sin(pi*x), u(x) = sin(pi*x)
  epsrel is the relative error threshold with the exact solution
*/
static void
test_toeplitz(const size_t N, const double epsrel)
{
  const size_t n = N - 2;                     /* subtract 2 to exclude boundaries */
  const double h = 1.0 / (N - 1.0);           /* grid spacing */
  const double tol = 1.0e-9;
  gsl_spmatrix *T = gsl_spmatrix_alloc(n ,n); /* triplet format */
  gsl_vector *f = gsl_vector_alloc(n);        /* right hand side vector */
  gsl_vector *u = gsl_vector_calloc(n);       /* solution vector, u0 = 0 */
  gsl_splinalg_gmres_workspace *w = gsl_splinalg_gmres_alloc(n);
  size_t i;
  int status;

  /* construct the sparse matrix for the finite difference equation */

  /* first row of matrix */
  gsl_spmatrix_set(T, 0, 0, -2.0);
  gsl_spmatrix_set(T, 0, 1, 1.0);

  /* loop over interior grid points */
  for (i = 1; i < n - 1; ++i)
    {
      gsl_spmatrix_set(T, i, i + 1, 1.0);
      gsl_spmatrix_set(T, i, i, -2.0);
      gsl_spmatrix_set(T, i, i - 1, 1.0);
    }

  /* last row of matrix */
  gsl_spmatrix_set(T, n - 1, n - 1, -2.0);
  gsl_spmatrix_set(T, n - 1, n - 2, 1.0);

  /* scale by h^2 */
  gsl_spmatrix_scale(T, 1.0 / (h * h));

  /* construct right hand side vector */
  for (i = 0; i < n; ++i)
    {
      double xi = (i + 1) * h;
      double fi = -M_PI * M_PI * sin(M_PI * xi);
      gsl_vector_set(f, i, fi);
    }

  /* solve the system */
  status = gsl_splinalg_gmres_solve_x(T, f, tol, u, w);
  gsl_test(status, "toeplitz status s=%d N=%zu", status, N);

  /* check solution against analytic */
  for (i = 0; i < n; ++i)
    {
      double xi = (i + 1) * h;
      double u_gsl = gsl_vector_get(u, i);
      double u_exact = sin(M_PI * xi);

      gsl_test_rel(u_gsl, u_exact, epsrel, "toeplitz N=%zu i=%zu", N, i);
    }

  /* check that the residual satisfies ||r|| <= tol*||f|| */
  {
    gsl_vector *r = gsl_vector_alloc(n);
    double normr, normf;

    gsl_vector_memcpy(r, f);
    gsl_spblas_dgemv(-1.0, T, u, 1.0, r);

    normr = gsl_blas_dnrm2(r);
    normf = gsl_blas_dnrm2(f);

    status = (normr <= tol*normf) != 1;
    gsl_test(status, "toeplitz residual N=%zu normr=%.12e normf=%.12e",
             N, normr, normf);

    gsl_vector_free(r);
  }

  gsl_splinalg_gmres_free(w);
  gsl_spmatrix_free(T);
  gsl_vector_free(f);
  gsl_vector_free(u);
} /* test_toeplitz() */

static void
test_random(const size_t N, const gsl_rng *r)
{
  const double tol = 1.0e-8;
  int status;
  gsl_spmatrix *A = create_random_sparse(N, N, 0.3, r);
  gsl_vector *b = gsl_vector_alloc(N);
  gsl_vector *x = gsl_vector_calloc(N);
  gsl_vector *res = gsl_vector_alloc(N);
  gsl_splinalg_gmres_workspace *w = gsl_splinalg_gmres_alloc(N);

  create_random_vector(b, r);

  do
    {
      status = gsl_splinalg_gmres_solve_x(A, b, tol, x, w);
      printf("%.12e\n", w->normr);
    }
  while (status == GSL_CONTINUE);
  gsl_test(status, "random status s=%d N=%zu", status, N);

  /* compute r = ||b - A*x|| */
  gsl_vector_memcpy(res, b);
  gsl_spblas_dgemv(-1.0, A, x, 1.0, res);

  gsl_spmatrix_free(A);
  gsl_vector_free(b);
  gsl_vector_free(x);
  gsl_vector_free(res);
  gsl_splinalg_gmres_free(w);
} /* test_random() */

int
main()
{
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  test_toeplitz(7, 1.0e-1);
  test_toeplitz(543, 1.0e-5);
  test_toeplitz(1000, 1.0e-6);
  test_toeplitz(5000, 1.0e-7);

  test_random(100, r);

  gsl_rng_free(r);

  exit (gsl_test_summary());
} /* main() */
