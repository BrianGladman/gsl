/* linalg/test_cholesky.c
 *
 * Copyright (C) 2016 Patrick Alken
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
#include <gsl/gsl_test.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_permutation.h>

static int create_random_vector(gsl_vector * v, gsl_rng * r);
static int create_posdef_matrix(gsl_matrix * m, gsl_rng * r);
static int create_hilbert_matrix(gsl_matrix * m);

static int test_cholesky_decomp_eps(const int scale, const gsl_matrix * m,
                                    const double eps, const char * desc);
static int test_cholesky_decomp(gsl_rng * r);
static int test_pcholesky_decomp(gsl_rng * r);
int test_pcholesky_solve_eps(const int scale, const gsl_matrix * m, const gsl_vector * rhs,
                             const gsl_vector * sol, const double eps,
                             const char * desc);
static int test_pcholesky_solve(gsl_rng * r);

static int
create_random_vector(gsl_vector * v, gsl_rng * r)
{
  const size_t N = v->size;
  size_t i;

  for (i = 0; i < N; ++i)
    {
      double vi = gsl_rng_uniform(r);
      gsl_vector_set(v, i, vi);
    }

  return GSL_SUCCESS;
}

static int
create_symm_matrix(gsl_matrix * m, gsl_rng * r)
{
  const size_t N = m->size1;
  size_t i, j;

  for (i = 0; i < N; ++i)
    {
      for (j = 0; j <= i; ++j)
        {
          double mij = gsl_rng_uniform(r);
          gsl_matrix_set(m, i, j, mij);
        }
    }

  /* copy lower triangle to upper */
  gsl_matrix_transpose_tricpy('L', 0, m, m);

  return GSL_SUCCESS;
}

static int
create_posdef_matrix(gsl_matrix * m, gsl_rng * r)
{
  const size_t N = m->size1;
  const double alpha = 10.0 * N;
  size_t i;

  /* The idea is to make a symmetric diagonally dominant
   * matrix. Make a symmetric matrix and add alpha*I to
   * its diagonal */

  create_symm_matrix(m, r);

  for (i = 0; i < N; ++i)
    {
      double mii = gsl_matrix_get(m, i, i);
      gsl_matrix_set(m, i, i, mii + alpha);
    }

  return GSL_SUCCESS;
}

static int
create_hilbert_matrix(gsl_matrix * m)
{
  const size_t N = m->size1;
  size_t i, j;

  for (i = 0; i < N; i++)
    {
      for (j = 0; j < N; j++)
        {
          gsl_matrix_set(m, i, j, 1.0/(i+j+1.0));
        }
    }

  return GSL_SUCCESS;
}

static int
test_cholesky_decomp_eps(const int scale, const gsl_matrix * m,
                         const double eps, const char * desc)
{
  int s = 0;
  size_t i, j, M = m->size1, N = m->size2;

  gsl_matrix * V  = gsl_matrix_alloc(M, N);
  gsl_matrix * A  = gsl_matrix_alloc(M, N);
  gsl_matrix * L  = gsl_matrix_alloc(M, N);
  gsl_matrix * LT = gsl_matrix_alloc(N, N);
  gsl_vector * S = gsl_vector_alloc(N);

  gsl_matrix_memcpy(V, m);

  if (scale)
    s += gsl_linalg_cholesky_decomp2(V, S);
  else
    s += gsl_linalg_cholesky_decomp(V);
  
  /* compute L and LT */
  for (i = 0; i < N ; i++)
    {
      for (j = 0; j < N; j++)
        {
          double Vij = gsl_matrix_get(V, i, j);
          gsl_matrix_set (L, i, j, i >= j ? Vij : 0);
          gsl_matrix_set (LT, i, j, i <= j ? Vij : 0);
        }
    }

  if (scale)
    {
      /* L <- S^{-1} L, LT <- LT S^{-1} */
      for (i = 0; i < N; ++i)
        {
          double Si = gsl_vector_get(S, i);
          gsl_vector_view v = gsl_matrix_row(L, i);
          gsl_vector_view w = gsl_matrix_column(LT, i);

          gsl_vector_scale(&v.vector, 1.0 / Si);
          gsl_vector_scale(&w.vector, 1.0 / Si);
        }
    }
            
  /* compute A = L LT */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, L, LT, 0.0, A);

  for (i = 0; i < M; i++)
    {
      for (j = 0; j < N; j++)
        {
          double Aij = gsl_matrix_get(A, i, j);
          double mij = gsl_matrix_get(m, i, j);

          gsl_test_rel(Aij, mij, eps,
                       "%s: (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, M, N, i, j, Aij, mij);
        }
    }

  gsl_matrix_free(V);
  gsl_matrix_free(A);
  gsl_matrix_free(L);
  gsl_matrix_free(LT);
  gsl_vector_free(S);

  return s;
}

static int
test_cholesky_decomp(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 50;
  size_t N;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix * m = gsl_matrix_alloc(N, N);

      create_posdef_matrix(m, r);
      test_cholesky_decomp_eps(0, m, N * GSL_DBL_EPSILON, "cholesky_decomp unscaled random");
      test_cholesky_decomp_eps(1, m, 20.0 * N * GSL_DBL_EPSILON, "cholesky_decomp scaled random");

      if (N <= 12)
        {
          create_hilbert_matrix(m);
          test_cholesky_decomp_eps(0, m, GSL_DBL_EPSILON, "cholesky_decomp unscaled hilbert");
          test_cholesky_decomp_eps(1, m, N * GSL_DBL_EPSILON, "cholesky_decomp scaled hilbert");
        }

      gsl_matrix_free(m);
    }

  return s;
}

static int
test_mcholesky_decomp_eps(const int scale, const gsl_matrix * m,
                          const double eps, const char * desc)
{
  int s = 0;
  size_t i, j, M = m->size1, N = m->size2;

  gsl_matrix * V  = gsl_matrix_alloc(M, N);
  gsl_matrix * A  = gsl_matrix_alloc(M, N);
  gsl_matrix * L  = gsl_matrix_alloc(M, N);
  gsl_matrix * LT = gsl_matrix_alloc(N, N);
  gsl_vector * D = gsl_vector_alloc(N);
  gsl_vector * E = gsl_vector_alloc(N);

  gsl_matrix_memcpy(V, m);

#if 0
  if (scale)
    s += gsl_linalg_cholesky_decomp2(V, D);
  else
#endif
    s += gsl_linalg_mcholesky_decomp(V, E);
  
  /* compute L and LT */
  for (i = 0; i < N ; i++)
    {
      for (j = 0; j < N; j++)
        {
          double Vij = gsl_matrix_get(V, i, j);
          gsl_matrix_set (L, i, j, i >= j ? Vij : 0);
          gsl_matrix_set (LT, i, j, i <= j ? Vij : 0);
        }
    }

  if (scale)
    {
      /* L <- D^{-1} L, LT <- LT D^{-1} */
      for (i = 0; i < N; ++i)
        {
          double Di = gsl_vector_get(D, i);
          gsl_vector_view v = gsl_matrix_row(L, i);
          gsl_vector_view w = gsl_matrix_column(LT, i);

          gsl_vector_scale(&v.vector, 1.0 / Di);
          gsl_vector_scale(&w.vector, 1.0 / Di);
        }
    }
            
  /* compute A = L LT */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, L, LT, 0.0, A);

  for (i = 0; i < M; i++)
    {
      double Eii = gsl_vector_get(E, i);

      for (j = 0; j < N; j++)
        {
          double Aij = gsl_matrix_get(A, i, j);
          double mij = gsl_matrix_get(m, i, j);

          /* add diag(E) to m */
          if (i == j)
            mij += Eii;

          gsl_test_rel(Aij, mij, eps,
                       "%s: (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, M, N, i, j, Aij, mij);
        }
    }

  gsl_matrix_free(V);
  gsl_matrix_free(A);
  gsl_matrix_free(L);
  gsl_matrix_free(LT);
  gsl_vector_free(D);
  gsl_vector_free(E);

  return s;
}

static int
test_mcholesky_decomp(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 50;
  size_t N;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix * m = gsl_matrix_alloc(N, N);

      create_posdef_matrix(m, r);
      test_mcholesky_decomp_eps(0, m, 10.0 * N * GSL_DBL_EPSILON, "mcholesky_decomp unscaled random posdef");

      create_symm_matrix(m, r);
      test_mcholesky_decomp_eps(0, m, 10.0 * N * GSL_DBL_EPSILON, "mcholesky_decomp unscaled random symm");

      if (N <= 12)
        {
          create_hilbert_matrix(m);
          test_mcholesky_decomp_eps(0, m, GSL_DBL_EPSILON, "mcholesky_decomp unscaled hilbert");
        }

      gsl_matrix_free(m);
    }

  return s;
}

static int
test_pcholesky_decomp_eps(const int scale, const gsl_matrix * m,
                          const double eps, const char * desc)
{
  int s = 0;
  size_t i, j, M = m->size1, N = m->size2;

  gsl_matrix * V  = gsl_matrix_alloc(M, N);
  gsl_matrix * A  = gsl_matrix_alloc(M, N);
  gsl_matrix * L  = gsl_matrix_alloc(M, N);
  gsl_matrix * LT = gsl_matrix_alloc(N, N);
  gsl_vector * S = gsl_vector_alloc(N);
  gsl_permutation * perm = gsl_permutation_alloc(N);

  gsl_vector_view D = gsl_matrix_diagonal(V);

  gsl_matrix_memcpy(V, m);

  if (scale)
    s += gsl_linalg_pcholesky_decomp2(V, perm, S);
  else
    s += gsl_linalg_pcholesky_decomp(V, perm);

  /* check that D is decreasing */
  s = 0;
  for (i = 1; i < N; ++i)
    {
      double dprev = gsl_vector_get(&D.vector, i - 1);
      double di = gsl_vector_get(&D.vector, i);

      if (di > dprev)
        s = 1;
    }

  gsl_test(s, "%s: (%zu,%zu): D is not decreasing",
           desc, M, N);
  
  /* compute L and LT */
  for (i = 0; i < N ; i++)
    {
      for (j = 0; j < N; j++)
        {
          double Vij = gsl_matrix_get(V, i, j);

          if (i == j)
            {
              gsl_matrix_set(L, i, j, 1.0);
              gsl_matrix_set(LT, i, j, 1.0);
            }
          else if (i > j)
            {
              gsl_matrix_set(L, i, j, Vij);
              gsl_matrix_set(LT, i, j, 0.0);
            }
          else
            {
              gsl_matrix_set(L, i, j, 0.0);
              gsl_matrix_set(LT, i, j, Vij);
            }
        }
    }

  /* compute (L sqrt(D)) and (sqrt(D) LT) */
  for (i = 0; i < N; ++i)
    {
      gsl_vector_view v = gsl_matrix_column(L, i);
      gsl_vector_view w = gsl_matrix_row(LT, i);
      double di = gsl_vector_get(&D.vector, i);

      gsl_vector_scale(&v.vector, sqrt(di));
      gsl_vector_scale(&w.vector, sqrt(di));
    }

  /* compute A = L D LT */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, L, LT, 0.0, A);

  /* compute V = P S M S P^T */

  gsl_matrix_memcpy(V, m);

  /* compute S M S */
  if (scale)
    {
      gsl_linalg_cholesky_scale_apply(V, S);
      gsl_matrix_transpose_tricpy('L', 0, V, V);
    }

  /* compute M P^T */
  for (i = 0; i < N; ++i)
    {
      gsl_vector_view v = gsl_matrix_row(V, i);
      gsl_permute_vector(perm, &v.vector);
    }

  /* compute P M P^T */
  for (i = 0; i < N; ++i)
    {
      gsl_vector_view v = gsl_matrix_column(V, i);
      gsl_permute_vector(perm, &v.vector);
    }

  for (i = 0; i < M; i++)
    {
      for (j = 0; j < N; j++)
        {
          double Aij = gsl_matrix_get(A, i, j); /* L D L^T */
          double Bij = gsl_matrix_get(V, i, j); /* P M P^T */

          gsl_test_rel(Aij, Bij, eps,
                       "%s: (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, M, N, i, j, Aij, Bij);
        }
    }

  gsl_matrix_free(V);
  gsl_matrix_free(A);
  gsl_matrix_free(L);
  gsl_matrix_free(LT);
  gsl_vector_free(S);
  gsl_permutation_free(perm);

  return s;
}

static int
test_pcholesky_decomp(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 50;
  size_t N;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix * m = gsl_matrix_alloc(N, N);

      create_posdef_matrix(m, r);
      test_pcholesky_decomp_eps(0, m, 20.0 * N * GSL_DBL_EPSILON, "pcholesky_decomp unscaled random");
      test_pcholesky_decomp_eps(1, m, 20.0 * N * GSL_DBL_EPSILON, "pcholesky_decomp scaled random");

      if (N <= 12)
        {
          create_hilbert_matrix(m);
          test_pcholesky_decomp_eps(0, m, N * GSL_DBL_EPSILON, "pcholesky_decomp unscaled hilbert");
          test_pcholesky_decomp_eps(1, m, N * GSL_DBL_EPSILON, "pcholesky_decomp scaled hilbert");
        }

      gsl_matrix_free(m);
    }

  return s;
}

int
test_pcholesky_solve_eps(const int scale, const gsl_matrix * m, const gsl_vector * rhs,
                         const gsl_vector * sol, const double eps,
                         const char * desc)
{
  int s = 0;
  size_t i, N = m->size1;
  gsl_matrix * u  = gsl_matrix_alloc(N, N);
  gsl_vector * x = gsl_vector_calloc(N);
  gsl_vector * S = gsl_vector_alloc(N);
  gsl_permutation * perm = gsl_permutation_alloc(N);

  gsl_matrix_memcpy(u, m);

  if (scale)
    {
      s += gsl_linalg_pcholesky_decomp2(u, perm, S);
      s += gsl_linalg_pcholesky_solve2(u, perm, S, rhs, x);
    }
  else
    {
      s += gsl_linalg_pcholesky_decomp(u, perm);
      s += gsl_linalg_pcholesky_solve(u, perm, rhs, x);
    }

  for (i = 0; i < N; i++)
    {
      double xi = gsl_vector_get(x, i);
      double yi = gsl_vector_get(sol, i);

      gsl_test_rel(xi, yi, eps,
                   "%s: %3lu[%lu]: %22.18g   %22.18g\n",
                   desc, N, i, xi, yi);
    }

  gsl_vector_free(x);
  gsl_vector_free(S);
  gsl_matrix_free(u);
  gsl_permutation_free(perm);

  return s;
}

static int
test_pcholesky_solve(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 50;
  size_t N;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix * m = gsl_matrix_alloc(N, N);
      gsl_vector * rhs = gsl_vector_alloc(N);
      gsl_vector * sol = gsl_vector_alloc(N);

      create_posdef_matrix(m, r);
      create_random_vector(sol, r);
      gsl_blas_dsymv(CblasLower, 1.0, m, sol, 0.0, rhs);

      test_pcholesky_solve_eps(0, m, rhs, sol, 16.0 * N * GSL_DBL_EPSILON, "pcholesky_solve unscaled random");
      test_pcholesky_solve_eps(1, m, rhs, sol, 16.0 * N * GSL_DBL_EPSILON, "pcholesky_solve scaled random");

      if (N <= 4)
        {
          create_hilbert_matrix(m);
          gsl_blas_dsymv(CblasLower, 1.0, m, sol, 0.0, rhs);
          test_pcholesky_solve_eps(0, m, rhs, sol, 512.0 * N * GSL_DBL_EPSILON, "pcholesky_solve unscaled hilbert");
          test_pcholesky_solve_eps(1, m, rhs, sol, 1024.0 * N * GSL_DBL_EPSILON, "pcholesky_solve scaled hilbert");
        }

      gsl_matrix_free(m);
      gsl_vector_free(rhs);
      gsl_vector_free(sol);
    }

  return s;
}

int
main(void)
{
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  gsl_ieee_env_setup ();

#if 0
  gsl_test(test_cholesky_decomp(r),       "Cholesky Decomposition");
  gsl_test(test_mcholesky_decomp(r),      "Modified Cholesky Decomposition");
#endif
  gsl_test(test_pcholesky_decomp(r),      "Pivoted Cholesky Decomposition");
  gsl_test(test_pcholesky_solve(r),       "Pivoted Cholesky Solve");

  gsl_rng_free(r);

  exit (gsl_test_summary());
}
