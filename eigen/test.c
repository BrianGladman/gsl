/* eigen/test.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2006, 2007 Gerard Jungman, Patrick Alken
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

/* Author:  G. Jungman
 */

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sort.h>

gsl_matrix *
create_hilbert_matrix(int size)
{
  int i, j;
  gsl_matrix * m = gsl_matrix_alloc(size, size);
  for(i=0; i<size; i++) {
    for(j=0; j<size; j++) {
      gsl_matrix_set(m, i, j, 1.0/(i+j+1.0));
    }
  }
  return m;
}

gsl_matrix *
create_random_symm_matrix(int size)
{
  int i, j;
  unsigned long k = 1;
  gsl_matrix * m = gsl_matrix_alloc(size, size);
  for(i=0; i<size; i++) {
    for(j=i; j<size; j++) {
      double x;
      k = (69069 * k + 1) & 0xffffffffUL;
      x = k / 4294967296.0;
      gsl_matrix_set(m, i, j, x);
      gsl_matrix_set(m, j, i, x);
    }
  }
  return m;
}

gsl_matrix_complex *
create_random_herm_matrix(int size)
{
  int i, j;
  unsigned long k = 1;
  gsl_matrix_complex * m = gsl_matrix_complex_alloc(size, size);
  for(i=0; i<size; i++) {
    for(j=i; j<size; j++) {
      gsl_complex z;
      k = (69069 * k + 1) & 0xffffffffUL;
      GSL_REAL(z) = k / 4294967296.0;
      k = (69069 * k + 1) & 0xffffffffUL;
      GSL_IMAG(z) = (i == j) ? 0 : k / 4294967296.0;
      gsl_matrix_complex_set(m, i, j, z);
      gsl_matrix_complex_set(m, j, i, gsl_complex_conjugate(z));
    }
  }
  return m;
}

void
create_random_nonsymm_matrix(gsl_matrix *m, gsl_rng *r, int lower,
                             int upper)
{
  size_t i, j;

  for (i = 0; i < m->size1; ++i)
    {
      for (j = 0; j < m->size2; ++j)
      {
        gsl_matrix_set(m,
                       i,
                       j,
                       gsl_rng_uniform(r) * (upper - lower) + lower);
      }
    }
} /* create_random_nonsymm_matrix() */

void
test_eigen_results (size_t N, const gsl_matrix * m, const gsl_vector * eval, 
                    const gsl_matrix * evec, const char * desc,
                    const char * desc2)
{
  size_t i,j;

  gsl_vector * x = gsl_vector_alloc(N);
  gsl_vector * y = gsl_vector_alloc(N);

  /* check eigenvalues */

  for (i = 0; i < N; i++)
    {
      double ei = gsl_vector_get (eval, i);
      gsl_vector_const_view vi = gsl_matrix_const_column(evec, i);
      gsl_vector_memcpy(x, &vi.vector);
      /* compute y = m x (should = lambda v) */
      gsl_blas_dgemv (CblasNoTrans, 1.0, m, x, 0.0, y);
      for (j = 0; j < N; j++)
        {
          double xj = gsl_vector_get (x, j);
          double yj = gsl_vector_get (y, j);
          gsl_test_rel(yj, ei * xj,  1e8 * GSL_DBL_EPSILON, 
                       "%s, eigenvalue(%d,%d), %s", desc, i, j, desc2);
        }
    }

  /* check eigenvectors are orthonormal */

  for (i = 0; i < N; i++)
    {
      gsl_vector_const_view vi = gsl_matrix_const_column(evec, i);
      double nrm_v = gsl_blas_dnrm2(&vi.vector);
      gsl_test_rel (nrm_v, 1.0, N * GSL_DBL_EPSILON, "%s, normalized(%d), %s", 
                    desc, i, desc2);
    }

  for (i = 0; i < N; i++)
    {
      gsl_vector_const_view vi = gsl_matrix_const_column(evec, i);
      for (j = i + 1; j < N; j++)
        {
          gsl_vector_const_view vj = gsl_matrix_const_column(evec, j);
          double vivj;
          gsl_blas_ddot (&vi.vector, &vj.vector, &vivj);
          gsl_test_abs (vivj, 0.0, N * GSL_DBL_EPSILON, 
                        "%s, orthogonal(%d,%d), %s", desc, i, j, desc2);
        }
    }

  gsl_vector_free(x);
  gsl_vector_free(y);
}


void
test_eigenvalues (size_t N, const gsl_vector *eval, const gsl_vector * eval2, 
                  const char * desc, const char * desc2)
{
  size_t i;
  for (i = 0; i < N; i++)
    {
      double ei = gsl_vector_get (eval, i);
      double e2i = gsl_vector_get (eval2, i);
      gsl_test_rel(ei, e2i, 100*N*GSL_DBL_EPSILON, 
                   "%s, direct eigenvalue(%d), %s",
                   desc, i, desc2);
    }
}

void
test_eigenvalues_complex (const gsl_vector_complex * eval,
                          const gsl_vector_complex * eval2, 
                          const char * desc, const char * desc2)
{
  const size_t N = eval->size;
  size_t i;

  for (i = 0; i < N; i++)
    {
      gsl_complex ei = gsl_vector_complex_get (eval, i);
      gsl_complex e2i = gsl_vector_complex_get (eval2, i);
      gsl_test_rel(GSL_REAL(ei), GSL_REAL(e2i), 10*N*GSL_DBL_EPSILON, 
                   "%s, direct eigenvalue(%d) real, %s",
                   desc, i, desc2);
      gsl_test_rel(GSL_IMAG(ei), GSL_IMAG(e2i), 10*N*GSL_DBL_EPSILON, 
                   "%s, direct eigenvalue(%d) imag, %s",
                   desc, i, desc2);
    }
}

void
test_eigen_complex_results (size_t N, const gsl_matrix_complex * m, 
                            const gsl_vector * eval, 
                            const gsl_matrix_complex * evec, 
                            const char * desc,
                            const char * desc2)
{
  size_t i,j;

  gsl_vector_complex * x = gsl_vector_complex_alloc(N);
  gsl_vector_complex * y = gsl_vector_complex_alloc(N);

  /* check eigenvalues */

  for (i = 0; i < N; i++)
    {
      double ei = gsl_vector_get (eval, i);
      gsl_vector_complex_const_view vi = gsl_matrix_complex_const_column(evec, i);
      gsl_vector_complex_memcpy(x, &vi.vector);
      /* compute y = m x (should = lambda v) */
      gsl_blas_zgemv (CblasNoTrans, GSL_COMPLEX_ONE, m, x, 
                      GSL_COMPLEX_ZERO, y);
      for (j = 0; j < N; j++)
        {
          gsl_complex xj = gsl_vector_complex_get (x, j);
          gsl_complex yj = gsl_vector_complex_get (y, j);
          gsl_test_rel(GSL_REAL(yj), ei * GSL_REAL(xj), 1e8*GSL_DBL_EPSILON, 
                       "%s, eigenvalue(%d,%d), real, %s", desc, i, j, desc2);
          gsl_test_rel(GSL_IMAG(yj), ei * GSL_IMAG(xj), 1e8*GSL_DBL_EPSILON, 
                       "%s, eigenvalue(%d,%d), imag, %s", desc, i, j, desc2);
        }
    }

  /* check eigenvectors are orthonormal */

  for (i = 0; i < N; i++)
    {
      gsl_vector_complex_const_view vi = gsl_matrix_complex_const_column(evec, i);
      double nrm_v = gsl_blas_dznrm2(&vi.vector);
      gsl_test_rel (nrm_v, 1.0, N * GSL_DBL_EPSILON, "%s, normalized(%d), %s", 
                    desc, i, desc2);
    }

  for (i = 0; i < N; i++)
    {
      gsl_vector_complex_const_view vi = gsl_matrix_complex_const_column(evec, i);
      for (j = i + 1; j < N; j++)
        {
          gsl_vector_complex_const_view vj 
            = gsl_matrix_complex_const_column(evec, j);
          gsl_complex vivj;
          gsl_blas_zdotc (&vi.vector, &vj.vector, &vivj);
          gsl_test_abs (gsl_complex_abs(vivj), 0.0, N * GSL_DBL_EPSILON, 
                        "%s, orthogonal(%d,%d), %s", desc, i, j, desc2);
        }
    }

  gsl_vector_complex_free(x);
  gsl_vector_complex_free(y);
}

void
test_eigen_symm(const char * desc, const gsl_matrix * m)
{
  size_t N = m->size1;

  gsl_matrix * A = gsl_matrix_alloc(N, N);
  gsl_matrix * evec = gsl_matrix_alloc(N, N);
  gsl_vector * eval = gsl_vector_alloc(N);
  gsl_vector * eval2 = gsl_vector_alloc(N);

  gsl_eigen_symm_workspace * w1 = gsl_eigen_symm_alloc (N);
  gsl_eigen_symmv_workspace * w2 = gsl_eigen_symmv_alloc (N);

  gsl_matrix_memcpy(A, m);
  gsl_eigen_symmv(A, eval, evec, w2);
  test_eigen_results (N, m, eval, evec, desc, "unsorted");

  gsl_matrix_memcpy(A, m);
  gsl_eigen_symm(A, eval2, w1);
  test_eigenvalues (N, eval, eval2, desc, "unsorted");

  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);
  test_eigen_results (N, m, eval, evec, desc, "val/asc");

  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);
  test_eigen_results (N, m, eval, evec, desc, "val/desc");

  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
  test_eigen_results (N, m, eval, evec, desc, "abs/asc");

  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
  test_eigen_results (N, m, eval, evec, desc, "abs/desc");

  gsl_eigen_symm_free (w1);
  gsl_eigen_symmv_free (w2);

  gsl_matrix_free(A);
  gsl_matrix_free(evec);
  gsl_vector_free(eval);
  gsl_vector_free(eval2);
}


void
test_eigen_herm(const char * desc, const gsl_matrix_complex * m)
{
  size_t N = m->size1;

  gsl_matrix_complex * A = gsl_matrix_complex_alloc(N, N);
  gsl_matrix_complex * evec = gsl_matrix_complex_alloc(N, N);
  gsl_vector * eval = gsl_vector_alloc(N);
  gsl_vector * eval2 = gsl_vector_alloc(N);

  gsl_eigen_herm_workspace * w1 = gsl_eigen_herm_alloc (N);
  gsl_eigen_hermv_workspace * w2 = gsl_eigen_hermv_alloc (N);

  gsl_matrix_complex_memcpy(A, m);
  gsl_eigen_hermv(A, eval, evec, w2);
  test_eigen_complex_results (N, m, eval, evec, desc, "unsorted");

  gsl_matrix_complex_memcpy(A, m);
  gsl_eigen_herm(A, eval2, w1);
  test_eigenvalues (N, eval, eval2, desc, "unsorted");

  gsl_eigen_hermv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);
  test_eigen_complex_results (N, m, eval, evec, desc, "val/asc");

  gsl_eigen_hermv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);
  test_eigen_complex_results (N, m, eval, evec, desc, "val/desc");

  gsl_eigen_hermv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
  test_eigen_complex_results (N, m, eval, evec, desc, "abs/asc");

  gsl_eigen_hermv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
  test_eigen_complex_results (N, m, eval, evec, desc, "abs/desc");

  gsl_eigen_herm_free (w1);
  gsl_eigen_hermv_free (w2);

  gsl_matrix_complex_free(A);
  gsl_matrix_complex_free(evec);
  gsl_vector_free(eval);
  gsl_vector_free(eval2);
}

void
test_eigen_jacobi(const char * desc, const gsl_matrix * m)
{
  size_t N = m->size1;
  unsigned int nrot;

  gsl_matrix * A = gsl_matrix_alloc(N, N);
  gsl_matrix * evec = gsl_matrix_alloc(N, N);
  gsl_vector * eval = gsl_vector_alloc(N);

  gsl_matrix_memcpy(A, m);
  gsl_eigen_jacobi(A, eval, evec, 1000, &nrot);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);

  test_eigen_results (N, m, eval, evec, desc, "");

  gsl_matrix_free(A);
  gsl_matrix_free(evec);
  gsl_vector_free(eval);
}


int test_invert_jacobi(void)
{
  int s = 0;
  int i, j;
  gsl_matrix * hminv = gsl_matrix_alloc(10, 10);
  gsl_matrix * id    = gsl_matrix_alloc(10, 10);

  /* 10x10 Hilbert matrix */
  gsl_matrix * hm = create_hilbert_matrix(10);
  gsl_eigen_invert_jacobi(hm, hminv, 1000);

  /* gsl_linalg_matmult(hm, hminv, id); */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, hm, hminv, 0.0, id);

  for(i=0; i<10; i++) {
    for(j=0; j<10; j++) {
      double delta_ij = ( i == j ? 1.0 : 0.0 );
      double id_ij    = gsl_matrix_get(id, i, j);
      
      int rs = ( fabs(id_ij - delta_ij) > 5.0e-3 );
      s += rs;
      gsl_test_abs(id_ij, delta_ij, 5e-3, "invert hilbert(10) %d,%d", i,j);
    }
  }

  gsl_test (s, "gsl_eigen_jacobi_invert hilbert(10)");

  gsl_matrix_free(hm);
  gsl_matrix_free(hminv);
  gsl_matrix_free(id);

  return s;
}

/******************************************
 * nonsymm test code                      *
 ******************************************/

void
test_eigen_nonsymm_results (const gsl_matrix * m, 
                            const gsl_vector_complex * eval, 
                            const gsl_matrix_complex * evec, 
                            size_t count,
                            const char * desc,
                            const char * desc2)
{
  size_t i,j;
  size_t N = m->size1;

  gsl_vector_complex * x = gsl_vector_complex_alloc(N);
  gsl_vector_complex * y = gsl_vector_complex_alloc(N);
  gsl_matrix_complex * A = gsl_matrix_complex_alloc(N, N);

  /* we need a complex matrix for the blas routines, so copy m into A */
  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          gsl_complex z;
          GSL_SET_COMPLEX(&z, gsl_matrix_get(m, i, j), 0.0);
          gsl_matrix_complex_set(A, i, j, z);
        }
    }

  for (i = 0; i < N; i++)
    {
      gsl_complex ei = gsl_vector_complex_get (eval, i);
      gsl_vector_complex_const_view vi = gsl_matrix_complex_const_column(evec, i);
      double norm = gsl_blas_dznrm2(&vi.vector);

      /* check that eigenvector is normalized */
      gsl_test_rel(norm, 1.0, N * GSL_DBL_EPSILON,
                   "nonsymm(N=%u,cnt=%u), %s, normalized(%d), %s", N, count, i, desc, desc2);

      gsl_vector_complex_memcpy(x, &vi.vector);

      /* compute y = m x (should = lambda v) */
      gsl_blas_zgemv (CblasNoTrans, GSL_COMPLEX_ONE, A, x, 
                      GSL_COMPLEX_ZERO, y);

      /* compute x = lambda v */
      gsl_blas_zscal(ei, x);

      /* now test if y = x */
      for (j = 0; j < N; j++)
        {
          gsl_complex xj = gsl_vector_complex_get (x, j);
          gsl_complex yj = gsl_vector_complex_get (y, j);

          /* use abs here in case the values are close to 0 */
          gsl_test_abs(GSL_REAL(yj), GSL_REAL(xj), 1e8*GSL_DBL_EPSILON, 
                       "nonsymm(N=%u,cnt=%u), %s, eigenvalue(%d,%d), real, %s", N, count, desc, i, j, desc2);
          gsl_test_abs(GSL_IMAG(yj), GSL_IMAG(xj), 1e8*GSL_DBL_EPSILON, 
                       "nonsymm(N=%u,cnt=%u), %s, eigenvalue(%d,%d), imag, %s", N, count, desc, i, j, desc2);
        }
    }

  gsl_matrix_complex_free(A);
  gsl_vector_complex_free(x);
  gsl_vector_complex_free(y);
}

/* test if A Z = Q S */
void
test_eigen_schur(const gsl_matrix * A, const gsl_matrix * S,
                 const gsl_matrix * Q, const gsl_matrix * Z,
                 size_t count, const char * desc,
                 const char * desc2)
{
  const size_t N = A->size1;
  size_t i, j;

  gsl_matrix * T1 = gsl_matrix_alloc(N, N);
  gsl_matrix * T2 = gsl_matrix_alloc(N, N);

  /* compute T1 = A Z */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, Z, 0.0, T1);

  /* compute T2 = Q S */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Q, S, 0.0, T2);

  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          double x = gsl_matrix_get(T1, i, j);
          double y = gsl_matrix_get(T2, i, j);

          gsl_test_rel(x, y, 1.0e8 * GSL_DBL_EPSILON,
                       "%s(N=%u,cnt=%u), %s, schur(%d,%d)", desc, N, count, desc2, i, j);
        }
    }

  gsl_matrix_free (T1);
  gsl_matrix_free (T2);
}

void
test_eigen_nonsymm_matrix(const gsl_matrix * m, size_t count,
                          const char * desc,
                          gsl_eigen_nonsymmv_workspace *w)
{
  const size_t N = m->size1;
  gsl_matrix * A = gsl_matrix_alloc(N, N);
  gsl_matrix * Z = gsl_matrix_alloc(N, N);
  gsl_matrix_complex * evec = gsl_matrix_complex_alloc(N, N);
  gsl_vector_complex * eval = gsl_vector_complex_alloc(N);

  /*
   * calculate eigenvalues and eigenvectors - it is sufficient to
   * test gsl_eigen_nonsymmv() since that function calls
   * gsl_eigen_nonsymm() for the eigenvalues
   */ 
  gsl_matrix_memcpy(A, m);
  gsl_eigen_nonsymmv(A, eval, evec, w);
  test_eigen_nonsymm_results (m, eval, evec, count, desc, "unsorted");

  /* test sort routines */
  gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
  test_eigen_nonsymm_results (m, eval, evec, count, desc, "abs/asc");

  gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
  test_eigen_nonsymm_results (m, eval, evec, count, desc, "abs/desc");

  /* test Schur vectors */
  gsl_matrix_memcpy(A, m);
  gsl_eigen_nonsymmv_Z(A, eval, evec, Z, w);
  gsl_linalg_hessenberg_set_zero(A);
  test_eigen_schur(m, A, Z, Z, count, "nonsymm", desc);

  gsl_matrix_free(A);
  gsl_matrix_free(Z);
  gsl_matrix_complex_free(evec);
  gsl_vector_complex_free(eval);
}

void
test_eigen_nonsymm(void)
{
  size_t N_max = 50;
  size_t n, i;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  for (n = 1; n <= N_max; ++n)
    {
      gsl_matrix * m = gsl_matrix_alloc(n, n);
      gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc(n);

      for (i = 0; i < 5; ++i)
        {
          create_random_nonsymm_matrix(m, r, -10, 10);

          gsl_eigen_nonsymm_params(1, 0, w->nonsymm_workspace_p);
          test_eigen_nonsymm_matrix(m, i, "random, unbalanced", w);

          gsl_eigen_nonsymm_params(1, 1, w->nonsymm_workspace_p);
          test_eigen_nonsymm_matrix(m, i, "random, balanced", w);
        }

      gsl_matrix_free(m);
      gsl_eigen_nonsymmv_free(w);
    }

  gsl_rng_free(r);

  {
    double dat1[] = { 0, 1, 1, 1,
                      1, 1, 1, 1,
                      0, 0, 0, 0,
                      0, 0, 0, 0 };
    double dat2[] = { 1, 1, 0, 1,
                      1, 1, 1, 1,
                      1, 1, 1, 1,
                      0, 1, 0, 0 };
    gsl_matrix_view v;
    gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc(4);
    
    v = gsl_matrix_view_array (dat1, 4, 4);
    test_eigen_nonsymm_matrix(&v.matrix, 0, "integer", w);

    v = gsl_matrix_view_array (dat2, 4, 4);
    test_eigen_nonsymm_matrix(&v.matrix, 1, "integer", w);

    gsl_eigen_nonsymmv_free(w);
  }
} /* test_eigen_nonsymm() */

/******************************************
 * gensymm test code                      *
 ******************************************/

void
create_random_sym_matrix(gsl_matrix *m, gsl_rng *r, int lower, int upper)
{
  size_t i, j;

  for (i = 0; i < m->size1; ++i)
    {
      for (j = i; j < m->size2; ++j)
      {
        double x = gsl_rng_uniform(r) * (upper - lower) + lower;
        gsl_matrix_set(m, i, j, x);
        gsl_matrix_set(m, j, i, x);
      }
    }
} /* create_random_sym_matrix() */

/* with r \in (0,1) if m_{ij} = r^{|i - j|} then m is positive definite */
void
create_random_posdef_matrix(gsl_matrix *m, gsl_rng *r)
{
  size_t i, j;
  double x = gsl_rng_uniform(r);

  for (i = 0; i < m->size1; ++i)
    {
      for (j = i; j < m->size2; ++j)
      {
        double a = pow(x, (double) (j - i));

        gsl_matrix_set(m, i, j, a);
        gsl_matrix_set(m, j, i, a);
      }
    }
} /* create_random_posdef_matrix() */

void
test_eigen_gensymm_results (const gsl_matrix * A, 
                            const gsl_matrix * B,
                            const gsl_vector * eval, 
                            const gsl_matrix * evec, 
                            size_t count,
                            const char * desc,
                            const char * desc2)
{
  const size_t N = A->size1;
  size_t i, j;

  gsl_vector * x = gsl_vector_alloc(N);
  gsl_vector * y = gsl_vector_alloc(N);
  gsl_vector * z = gsl_vector_alloc(N);

  /* check A v = lambda B v */
  for (i = 0; i < N; i++)
    {
      double ei = gsl_vector_get (eval, i);
      gsl_vector_const_view vi = gsl_matrix_const_column(evec, i);
      double norm = gsl_blas_dnrm2(&vi.vector);

      /* check that eigenvector is normalized */
      gsl_test_rel(norm, 1.0, N * GSL_DBL_EPSILON,
                   "gensymm(N=%u,cnt=%u), %s, normalized(%d), %s", N, count,
                   i, desc, desc2);

      gsl_vector_memcpy(z, &vi.vector);

      /* compute y = A z */
      gsl_blas_dgemv (CblasNoTrans, 1.0, A, z, 0.0, y);

      /* compute x = B z */
      gsl_blas_dgemv (CblasNoTrans, 1.0, B, z, 0.0, x);

      /* compute x = lambda B z */
      gsl_blas_dscal(ei, x);

      /* now test if y = x */
      for (j = 0; j < N; j++)
        {
          double xj = gsl_vector_get (x, j);
          double yj = gsl_vector_get (y, j);

          gsl_test_rel(yj, xj, 1e9 * GSL_DBL_EPSILON, 
                       "gensymm(N=%u,cnt=%u), %s, eigenvalue(%d,%d), real, %s", N, count, desc, i, j, desc2);
        }
    }

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(z);
}

void
test_eigen_gensymm(void)
{
  size_t N_max = 50;
  size_t n, i;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  int s;

  for (n = 1; n <= N_max; ++n)
    {
      gsl_matrix * A = gsl_matrix_alloc(n, n);
      gsl_matrix * B = gsl_matrix_alloc(n, n);
      gsl_matrix * ma = gsl_matrix_alloc(n, n);
      gsl_matrix * mb = gsl_matrix_alloc(n, n);
      gsl_vector * eval = gsl_vector_alloc(n);
      gsl_vector * evalv = gsl_vector_alloc(n);
      gsl_vector * x = gsl_vector_alloc(n);
      gsl_vector * y = gsl_vector_alloc(n);
      gsl_matrix * evec = gsl_matrix_alloc(n, n);
      gsl_eigen_gensymm_workspace * w = gsl_eigen_gensymm_alloc(n);
      gsl_eigen_gensymmv_workspace * wv = gsl_eigen_gensymmv_alloc(n);

      for (i = 0; i < 10; ++i)
        {
          create_random_sym_matrix(A, r, -10, 10);
          create_random_posdef_matrix(B, r);

          gsl_matrix_memcpy(ma, A);
          gsl_matrix_memcpy(mb, B);

          gsl_eigen_gensymmv(ma, mb, evalv, evec, wv);
          test_eigen_gensymm_results(A, B, evalv, evec, i, "random", "unsorted");

          gsl_matrix_memcpy(ma, A);
          gsl_matrix_memcpy(mb, B);

          gsl_eigen_gensymm(ma, mb, eval, w);

          /* eval and evalv have to be sorted? not sure why */
          gsl_vector_memcpy(x, eval);
          gsl_vector_memcpy(y, evalv);
          gsl_sort_vector(x);
          gsl_sort_vector(y);
          test_eigenvalues(n, y, x, "gensymm, random", "unsorted");

          gsl_eigen_gensymmv_sort(evalv, evec, GSL_EIGEN_SORT_VAL_ASC);
          test_eigen_gensymm_results(A, B, evalv, evec, i, "random", "val/asc");

          gsl_eigen_gensymmv_sort(evalv, evec, GSL_EIGEN_SORT_VAL_DESC);
          test_eigen_gensymm_results(A, B, evalv, evec, i, "random", "val/desc");

          gsl_eigen_gensymmv_sort(evalv, evec, GSL_EIGEN_SORT_ABS_ASC);
          test_eigen_gensymm_results(A, B, evalv, evec, i, "random", "abs/asc");
          gsl_eigen_gensymmv_sort(evalv, evec, GSL_EIGEN_SORT_ABS_DESC);
          test_eigen_gensymm_results(A, B, evalv, evec, i, "random", "abs/desc");
        }

      gsl_matrix_free(A);
      gsl_matrix_free(B);
      gsl_matrix_free(ma);
      gsl_matrix_free(mb);
      gsl_vector_free(eval);
      gsl_vector_free(evalv);
      gsl_vector_free(x);
      gsl_vector_free(y);
      gsl_matrix_free(evec);
      gsl_eigen_gensymm_free(w);
      gsl_eigen_gensymmv_free(wv);
    }

  gsl_rng_free(r);
} /* test_eigen_gensymm() */

/******************************************
 * gen test code                          *
 ******************************************/

typedef struct
{
  gsl_matrix *A;
  gsl_matrix *B;
  gsl_vector_complex *alpha;
  gsl_vector *beta;
  gsl_vector_complex *alphav;
  gsl_vector *betav;
  gsl_vector_complex *eval;
  gsl_vector_complex *evalv;
  gsl_vector *x;
  gsl_vector *y;
  gsl_matrix *Q;
  gsl_matrix *Z;
  gsl_matrix_complex *evec;
  gsl_eigen_gen_workspace *gen_p;
  gsl_eigen_genv_workspace *genv_p;
} test_eigen_gen_workspace;

test_eigen_gen_workspace *
test_eigen_gen_alloc(const size_t n)
{
  test_eigen_gen_workspace *w;

  w = calloc(1, sizeof(test_eigen_gen_workspace));

  w->A = gsl_matrix_alloc(n, n);
  w->B = gsl_matrix_alloc(n, n);
  w->alpha = gsl_vector_complex_alloc(n);
  w->beta = gsl_vector_alloc(n);
  w->alphav = gsl_vector_complex_alloc(n);
  w->betav = gsl_vector_alloc(n);
  w->eval = gsl_vector_complex_alloc(n);
  w->evalv = gsl_vector_complex_alloc(n);
  w->x = gsl_vector_alloc(n);
  w->y = gsl_vector_alloc(n);
  w->Q = gsl_matrix_alloc(n, n);
  w->Z = gsl_matrix_alloc(n, n);
  w->evec = gsl_matrix_complex_alloc(n, n);
  w->gen_p = gsl_eigen_gen_alloc(n);
  w->genv_p = gsl_eigen_genv_alloc(n);

  return (w);
} /* test_eigen_gen_alloc() */

void
test_eigen_gen_free(test_eigen_gen_workspace *w)
{
  gsl_matrix_free(w->A);
  gsl_matrix_free(w->B);
  gsl_vector_complex_free(w->alpha);
  gsl_vector_free(w->beta);
  gsl_vector_complex_free(w->alphav);
  gsl_vector_free(w->betav);
  gsl_vector_complex_free(w->eval);
  gsl_vector_complex_free(w->evalv);
  gsl_vector_free(w->x);
  gsl_vector_free(w->y);
  gsl_matrix_free(w->Q);
  gsl_matrix_free(w->Z);
  gsl_matrix_complex_free(w->evec);
  gsl_eigen_gen_free(w->gen_p);
  gsl_eigen_genv_free(w->genv_p);
} /* test_eigen_gen_free() */

void
test_eigen_gen_results (const gsl_matrix * A, const gsl_matrix * B,
                        const gsl_vector_complex * alpha, 
                        const gsl_vector * beta,
                        const gsl_matrix_complex * evec, 
                        size_t count, const char * desc,
                        const char * desc2)
{
  const size_t N = A->size1;
  size_t i, j;
  gsl_matrix_complex *ma, *mb;
  gsl_vector_complex *x, *y;
  gsl_complex z_one, z_zero;

  ma = gsl_matrix_complex_alloc(N, N);
  mb = gsl_matrix_complex_alloc(N, N);
  y = gsl_vector_complex_alloc(N);
  x = gsl_vector_complex_alloc(N);

  /* ma <- A, mb <- B */
  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          gsl_complex z;

          GSL_SET_COMPLEX(&z, gsl_matrix_get(A, i, j), 0.0);
          gsl_matrix_complex_set(ma, i, j, z);

          GSL_SET_COMPLEX(&z, gsl_matrix_get(B, i, j), 0.0);
          gsl_matrix_complex_set(mb, i, j, z);
        }
    }

  GSL_SET_COMPLEX(&z_one, 1.0, 0.0);
  GSL_SET_COMPLEX(&z_zero, 0.0, 0.0);

  /* check eigenvalues */
  for (i = 0; i < N; ++i)
    {
      gsl_vector_complex_const_view vi =
        gsl_matrix_complex_const_column(evec, i);
      gsl_complex ai = gsl_vector_complex_get(alpha, i);
      double bi = gsl_vector_get(beta, i);

      /* compute x = alpha * B * v */
      gsl_blas_zgemv(CblasNoTrans, z_one, mb, &vi.vector, z_zero, x);
      gsl_blas_zscal(ai, x);

      /* compute y = beta * A v */
      gsl_blas_zgemv(CblasNoTrans, z_one, ma, &vi.vector, z_zero, y);
      gsl_blas_zdscal(bi, y);

      /* now test if y = x */
      for (j = 0; j < N; ++j)
        {
          gsl_complex xj = gsl_vector_complex_get(x, j);
          gsl_complex yj = gsl_vector_complex_get(y, j);

          gsl_test_abs(GSL_REAL(yj), GSL_REAL(xj), 1e8*GSL_DBL_EPSILON, 
                       "gen(N=%u,cnt=%u), %s, eigenvalue(%d,%d), real, %s",
                       N, count, desc, i, j, desc2);
          gsl_test_abs(GSL_IMAG(yj), GSL_IMAG(xj), 1e8*GSL_DBL_EPSILON, 
                       "gen(N=%u,cnt=%u), %s, eigenvalue(%d,%d), real, %s",
                       N, count, desc, i, j, desc2);
        }
    }

  gsl_matrix_complex_free(ma);
  gsl_matrix_complex_free(mb);
  gsl_vector_complex_free(y);
  gsl_vector_complex_free(x);
} /* test_eigen_gen_results() */

void
test_eigen_gen_pencil(const gsl_matrix * A, const gsl_matrix * B,
                      size_t count, const char * desc, int test_schur,
                      test_eigen_gen_workspace *w)
{
  const size_t N = A->size1;
  size_t i;

  gsl_matrix_memcpy(w->A, A);
  gsl_matrix_memcpy(w->B, B);

  if (test_schur)
    {
      gsl_eigen_genv_QZ(w->A, w->B, w->alphav, w->betav, w->evec, w->Q, w->Z, w->genv_p);
      test_eigen_schur(A, w->A, w->Q, w->Z, count, "genv/A", desc);
      test_eigen_schur(B, w->B, w->Q, w->Z, count, "genv/B", desc);
    }
  else
    gsl_eigen_genv(w->A, w->B, w->alphav, w->betav, w->evec, w->genv_p);

  test_eigen_gen_results(A, B, w->alphav, w->betav, w->evec, count, "random", "unsorted");

  gsl_matrix_memcpy(w->A, A);
  gsl_matrix_memcpy(w->B, B);

  if (test_schur)
    {
      gsl_eigen_gen_params(1, 1, 0, w->gen_p);
      gsl_eigen_gen_QZ(w->A, w->B, w->alpha, w->beta, w->Q, w->Z, w->gen_p);
      test_eigen_schur(A, w->A, w->Q, w->Z, count, "gen/A", desc);
      test_eigen_schur(B, w->B, w->Q, w->Z, count, "gen/B", desc);
    }
  else
    {
      gsl_eigen_gen_params(0, 0, 0, w->gen_p);
      gsl_eigen_gen(w->A, w->B, w->alpha, w->beta, w->gen_p);
    }

  /* compute eval = alpha / beta values */
  for (i = 0; i < N; ++i)
    {
      gsl_complex z, ai;
      double bi;

      ai = gsl_vector_complex_get(w->alpha, i);
      bi = gsl_vector_get(w->beta, i);
      GSL_SET_COMPLEX(&z, GSL_REAL(ai) / bi, GSL_IMAG(ai) / bi);
      gsl_vector_complex_set(w->eval, i, z);

      ai = gsl_vector_complex_get(w->alphav, i);
      bi = gsl_vector_get(w->betav, i);
      GSL_SET_COMPLEX(&z, GSL_REAL(ai) / bi, GSL_IMAG(ai) / bi);
      gsl_vector_complex_set(w->evalv, i, z);
    }

  /* sort eval and evalv and test them */
  gsl_eigen_nonsymmv_sort(w->eval, NULL, GSL_EIGEN_SORT_ABS_ASC);
  gsl_eigen_nonsymmv_sort(w->evalv, NULL, GSL_EIGEN_SORT_ABS_ASC);
  test_eigenvalues_complex(w->evalv, w->eval, "gen", "random");

  gsl_eigen_genv_sort(w->alphav, w->betav, w->evec, GSL_EIGEN_SORT_ABS_ASC);
  test_eigen_gen_results(A, B, w->alphav, w->betav, w->evec, count, "random", "abs/asc");
  gsl_eigen_genv_sort(w->alphav, w->betav, w->evec, GSL_EIGEN_SORT_ABS_DESC);
  test_eigen_gen_results(A, B, w->alphav, w->betav, w->evec, count, "random", "abs/desc");
} /* test_eigen_gen_pencil() */

void
test_eigen_gen(void)
{
  size_t N_max = 20;
  size_t n, i;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  int s;

  for (n = 1; n <= N_max; ++n)
    {
      gsl_matrix * A = gsl_matrix_alloc(n, n);
      gsl_matrix * B = gsl_matrix_alloc(n, n);
      test_eigen_gen_workspace * w = test_eigen_gen_alloc(n);

      for (i = 0; i < 5; ++i)
        {
          create_random_nonsymm_matrix(A, r, -10, 10);
          create_random_nonsymm_matrix(B, r, -10, 10);

          test_eigen_gen_pencil(A, B, i, "random", 0, w);
          test_eigen_gen_pencil(A, B, i, "random", 1, w);
        }

      gsl_matrix_free(A);
      gsl_matrix_free(B);
      test_eigen_gen_free(w);
    }

  gsl_rng_free(r);
} /* test_eigen_gen() */

int main()
{
  gsl_ieee_env_setup ();
  gsl_rng_env_setup ();

  {
    double r[] =  { 0,  0, -1,  0, 
                    0,  1,  0,  1,
                    -1,  0,  0,  0,
                    0,  1,  0,  0 };
    gsl_matrix_view s4 = gsl_matrix_view_array (r, 4, 4);

    test_eigen_symm("symm(4)", &s4.matrix);
  }
  
  {
    double c[] =  { 0,0,  0,0, -1,0,  0,0, 
                    0,0,  1,0,  0,0,  1,0,
                    -1,0,  0,0,  0,0,  0,0,
                    0,0,  1,0,  0,0,  0,0 };

    gsl_matrix_complex_view h4 = gsl_matrix_complex_view_array (c, 4, 4);

    test_eigen_herm("herm(4)", &h4.matrix);
  }

  {
    double r[] =  { 1,  0,  0,  0, 
                    0,  2,  0,  0,
                    0,  0,  3,  0,
                    0,  0,  0,  4 };
    gsl_matrix_view s4 = gsl_matrix_view_array (r, 4, 4);

    test_eigen_symm("symm(4) diag", &s4.matrix);
  }
  
  {
    double c[] =  { 1,0,  0,0, 0,0,  0,0, 
                    0,0,  2,0, 0,0,  0,0,
                    0,0,  0,0, 3,0,  0,0,
                    0,0,  0,0, 0,0,  4,0 };

    gsl_matrix_complex_view h4 = gsl_matrix_complex_view_array (c, 4, 4);

    test_eigen_herm("herm(4) diag", &h4.matrix);
  }

  {
    gsl_matrix *rs10 = create_random_symm_matrix (10);
    test_eigen_symm("symm(10)", rs10);
    gsl_matrix_free (rs10);
  }

  {
    gsl_matrix_complex *rh10 = create_random_herm_matrix (10);
    test_eigen_herm("herm(10)", rh10);
    gsl_matrix_complex_free (rh10);
  }

  test_eigen_nonsymm();
  test_eigen_gensymm();
  test_eigen_gen();

#if 0 /* Deprecated functions */
  {
    gsl_matrix *h5 = create_hilbert_matrix (5);
    test_eigen_jacobi("hilbert(5)", h5); 
    test_invert_jacobi(); 
    gsl_matrix_free (h5); 
  }
#endif

  exit (gsl_test_summary());
}
