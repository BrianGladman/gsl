/* multilarge/test.c
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
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multilarge.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_ieee_utils.h>

static void test_random_matrix_orth(gsl_matrix *m, const gsl_rng *r);
static void test_random_matrix_ill(gsl_matrix *m, const gsl_rng *r);
static void test_random_vector(gsl_vector *v, const gsl_rng *r,
                               const double lower, const double upper);
static void test_random_matrix(gsl_matrix *m, const gsl_rng *r,
                               const double lower, const double upper);
static void test_random_vector_noise(const gsl_rng *r, gsl_vector *y);
static void test_compare_vectors(const double tol, const gsl_vector * a,
                                 const gsl_vector * b, const char * desc);
static void test_multifit_solve(const double lambda, const gsl_matrix * X,
                                const gsl_vector * y, gsl_vector * c);
static void test_multilarge_solve(const gsl_multilarge_linear_type * T, const double lambda,
                                  const gsl_matrix * X, const gsl_vector * y, gsl_vector * c);

#include "test_normal.c"
#include "test_tsqr.c"

/* generate random square orthogonal matrix via QR decomposition */
static void
test_random_matrix_orth(gsl_matrix *m, const gsl_rng *r)
{
  const size_t M = m->size1;
  gsl_matrix *A = gsl_matrix_alloc(M, M);
  gsl_vector *tau = gsl_vector_alloc(M);
  gsl_matrix *R = gsl_matrix_alloc(M, M);

  test_random_matrix(A, r, -1.0, 1.0);
  gsl_linalg_QR_decomp(A, tau);
  gsl_linalg_QR_unpack(A, tau, m, R);

  gsl_matrix_free(A);
  gsl_matrix_free(R);
  gsl_vector_free(tau);
}

/* construct ill-conditioned matrix via SVD */
static void
test_random_matrix_ill(gsl_matrix *m, const gsl_rng *r)
{
  const size_t M = m->size1;
  const size_t N = m->size2;
  gsl_matrix *U = gsl_matrix_alloc(M, M);
  gsl_matrix *V = gsl_matrix_alloc(N, N);
  gsl_vector *S = gsl_vector_alloc(N);
  gsl_matrix_view Uv = gsl_matrix_submatrix(U, 0, 0, M, N);
  const double smin = 16.0 * GSL_DBL_EPSILON;
  const double smax = 10.0;
  const double ratio = pow(smin / smax, 1.0 / (N - 1.0));
  double s;
  size_t j;

  test_random_matrix_orth(U, r);
  test_random_matrix_orth(V, r);

  /* compute U * S */

  s = smax;
  for (j = 0; j < N; ++j)
    {
      gsl_vector_view uj = gsl_matrix_column(U, j);

      gsl_vector_scale(&uj.vector, s);
      s *= ratio;
    }

  /* compute m = (U * S) * V' */
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &Uv.matrix, V, 0.0, m);

  gsl_matrix_free(U);
  gsl_matrix_free(V);
  gsl_vector_free(S);
}

static void
test_random_vector(gsl_vector *v, const gsl_rng *r,
                   const double lower, const double upper)
{
  size_t i;
  size_t N = v->size;

  for (i = 0; i < N; ++i)
    {
      gsl_vector_set(v, i,
                     gsl_rng_uniform(r) * (upper - lower) + lower);
    }
}

static void
test_random_matrix(gsl_matrix *m, const gsl_rng *r,
                   const double lower, const double upper)
{
  size_t i, j;
  size_t M = m->size1;
  size_t N = m->size2;

  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
      {
        gsl_matrix_set(m, i, j,
                       gsl_rng_uniform(r) * (upper - lower) + lower);
      }
    }
}

static void
test_random_vector_noise(const gsl_rng *r, gsl_vector *y)
{
  size_t i;

  for (i = 0; i < y->size; ++i)
    {
      double *ptr = gsl_vector_ptr(y, i);
      *ptr += 1.0e-3 * gsl_rng_uniform(r);
    }
}

static void
test_compare_vectors(const double tol, const gsl_vector * a,
                     const gsl_vector * b, const char * desc)
{
  size_t i;

  for (i = 0; i < a->size; ++i)
    {
      double ai = gsl_vector_get(a, i);
      double bi = gsl_vector_get(b, i);

      gsl_test_rel(bi, ai, tol, "%s i=%zu", desc, i);
    }
}

/* solve least squares system with multifit SVD */
static void
test_multifit_solve(const double lambda, const gsl_matrix * X,
                    const gsl_vector * y, gsl_vector * c)
{
  const size_t n = X->size1;
  const size_t p = X->size2;
  gsl_multifit_linear_workspace *w =
    gsl_multifit_linear_alloc(n, p);
  double rnorm, snorm;

  gsl_multifit_linear_svd(X, w);
  gsl_multifit_linear_solve(lambda, X, y, c, &rnorm, &snorm, w);

  gsl_multifit_linear_free(w);
}

/* solve least squares system with multilarge */
static void
test_multilarge_solve(const gsl_multilarge_linear_type * T, const double lambda,
                      const gsl_matrix * X, const gsl_vector * y, gsl_vector * c)
{
  const size_t n = X->size1;
  const size_t p = X->size2;
  const size_t nblock = 5;
  const size_t nrows = n / nblock; /* number of rows per block */
  gsl_multilarge_linear_workspace *w =
    gsl_multilarge_linear_alloc(T, nrows, p);
  size_t rowidx = 0;
  double rnorm, snorm;

  while (rowidx < n)
    {
      size_t nleft = n - rowidx;
      size_t nr = GSL_MIN(nrows, nleft);
      gsl_matrix_const_view Xv = gsl_matrix_const_submatrix(X, rowidx, 0, nr, p);
      gsl_vector_const_view yv = gsl_vector_const_subvector(y, rowidx, nr);

      gsl_multilarge_linear_accumulate(&Xv.matrix, &yv.vector, w);

      rowidx += nr;
    }

  gsl_multilarge_linear_solve(lambda, c, &rnorm, &snorm, w);

  gsl_multilarge_linear_free(w);
}

int
main (void)
{
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  gsl_ieee_env_setup();

  test_normal(r);
  test_tsqr(r);

  gsl_rng_free(r);

  exit (gsl_test_summary ());
}
