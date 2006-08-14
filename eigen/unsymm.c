/* eigen/unsymm.c
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

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_matrix.h>

/*
 * This module computes the eigenvalues of a real unsymmetric
 * matrix, using the QR decomposition.
 *
 * See Golub & Van Loan, "Matrix Computations" (3rd ed),
 * algorithm 7.5.2
 */

inline static size_t schur_decomp(gsl_matrix * H, size_t top,
                                  size_t bot, gsl_vector_complex * eval,
                                  size_t evidx,
                                  size_t * nit,
                                  gsl_eigen_unsymm_workspace * w);
inline static size_t zero_subdiag_small_elements(gsl_matrix * A);
static inline int francis_qrstep(gsl_matrix * H,
                                 gsl_eigen_unsymm_workspace * w,
                                 double s, double t);
static void get_2b2_eigenvalues(gsl_matrix * A, gsl_complex * e1,
                                gsl_complex * e2);

/*
gsl_eigen_unsymm_alloc()

Allocate a workspace for solving the unsymmetric
eigenvalue/eigenvector problem

Inputs: n - size of matrix

Return: pointer to workspace
*/

gsl_eigen_unsymm_workspace *
gsl_eigen_unsymm_alloc(const size_t n)
{
  gsl_eigen_unsymm_workspace *w;

  if (n == 0)
    {
      GSL_ERROR_NULL ("matrix dimension must be positive integer",
                      GSL_EINVAL);
    }

  w = ((gsl_eigen_unsymm_workspace *)
       malloc (sizeof (gsl_eigen_unsymm_workspace)));

  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace", GSL_ENOMEM);
    }

  w->size = n;
  w->max_iterations = 30 * n;
  w->work = gsl_vector_alloc(n);

  return (w);
}

void
gsl_eigen_unsymm_free (gsl_eigen_unsymm_workspace * w)
{
  gsl_vector_free(w->work);
  free(w);
}

/*
gsl_eigen_unsymm()

Solve the unsymmetric eigenvalue problem

A x = \lambda x

for the eigenvalues \lambda using algorithm 7.5.2 of
Golub & Van Loan, "Matrix Computations" (3rd ed)

Inputs: A    - matrix
        eval - where to store eigenvalues
        w    - workspace

Notes: On output, A contains the upper block triangular Schur
       decomposition.
*/

int
gsl_eigen_unsymm (gsl_matrix * A, gsl_vector_complex * eval,
                  gsl_eigen_unsymm_workspace * w)
{
  /* check matrix and vector sizes */

  if (A->size1 != A->size2)
    {
      GSL_ERROR ("matrix must be square to compute eigenvalues", GSL_ENOTSQR);
    }
  else if (eval->size != A->size1)
    {
      GSL_ERROR ("eigenvalue vector must match matrix size", GSL_EBADLEN);
    }
  else
    {
      size_t N;
      gsl_complex lambda1, lambda2; /* eigenvalues */
      size_t nit;

      N = A->size1;

      /* special cases */
      if (N == 1)
      {
        GSL_SET_COMPLEX(&lambda1, gsl_matrix_get(A, 0, 0), 0.0);
        gsl_vector_complex_set(eval, 0, lambda1);
        return GSL_SUCCESS;
      }

      if (N == 2)
      {
        /*
         * The 2x2 case is special since the matrix is already
         * in upper quasi-triangular form so no Schur decomposition
         * is necessary
         */
        get_2b2_eigenvalues(A, &lambda1, &lambda2);
        gsl_vector_complex_set(eval, 0, lambda1);
        gsl_vector_complex_set(eval, 1, lambda2);
        return GSL_SUCCESS;
      }

      /* balance the matrix */
      gsl_linalg_balance_matrix(A,w->work); 

      /* compute the Hessenberg reduction of A */
      gsl_linalg_hessenberg(A,w->work); 

      nit = 0;

      /*
       * compute Schur decomposition of A and store eigenvalues
       * into eval
       */
      schur_decomp(A, 0, N - 1, eval, 0, &nit, w);

      if (nit > w->max_iterations)
        {
          GSL_ERROR("maximum iterations exceeded", GSL_EMAXITER);
        }

      return GSL_SUCCESS;
    }
} /* gsl_eigen_unsymm() */

/********************************************
 *           INTERNAL ROUTINES              *
 ********************************************/

/*
schur_decomp()
  Compute the Schur decomposition of the submatrix of H
starting from (top, top) to (bot, bot)

Inputs: H     - hessenberg matrix
        top   - top index
        bot   - bottom index
        eval  - where to store eigenvalues
        evidx - index into eval
        nit   - running total of number of QR iterations since
                we last found an eigenvalue
                (must be initialized before calling this function)
        w     - workspace

Return: number of eigenvalues found
*/

inline static size_t
schur_decomp(gsl_matrix * H, size_t top, size_t bot,
             gsl_vector_complex * eval, size_t evidx,
             size_t * nit,
             gsl_eigen_unsymm_workspace * w)
{
  gsl_matrix_view m;
  size_t N;    /* size of matrix */
  double s, t; /* shifts */
  size_t nev;  /* number of eigenvalues found so far */
  size_t q;
  gsl_complex lambda1, /* eigenvalues */
              lambda2;

  N = bot - top + 1;

  if (N == 1)
    {
      GSL_SET_COMPLEX(&lambda1, gsl_matrix_get(H, 0, 0), 0.0);
      gsl_vector_complex_set(eval, evidx, lambda1);
      *nit = 0;
      return 1;
    }
  else if (N == 2)
    {
      get_2b2_eigenvalues(H, &lambda1, &lambda2);
      gsl_vector_complex_set(eval, evidx, lambda1);
      gsl_vector_complex_set(eval, evidx + 1, lambda2);
      *nit = 0;
      return 2;
    }

  m = gsl_matrix_submatrix(H, top, top, N, N);

  nev = 0;
  while ((N > 2) && ((*nit)++ < w->max_iterations))
    {
      if ((*nit == 10) || (*nit == 20))
        {
          /*
           * We have gone 10 or 20 iterations without finding
           * a new eigenvalue, try a new choice of shifts.
           * See Numerical Recipes in C, eq 11.6.27
           */
          t = fabs(gsl_matrix_get(&m.matrix, N - 1, N - 2)) +
              fabs(gsl_matrix_get(&m.matrix, N - 2, N - 3));
          s = 1.5 * t;
          t *= t;
        }
      else
        {
          /* s = a_1 + a_2 = m_{mm} + m_{nn} where m = n - 1 */
          s = gsl_matrix_get(&m.matrix, N - 2, N - 2) +
              gsl_matrix_get(&m.matrix, N - 1, N - 1);

          /* t = a_1 * a_2 = m_{mm} * m_{nn} - m_{mn} * m_{nm} */
          t = (gsl_matrix_get(&m.matrix, N - 2, N - 2) *
               gsl_matrix_get(&m.matrix, N - 1, N - 1)) -
              (gsl_matrix_get(&m.matrix, N - 2, N - 1) *
               gsl_matrix_get(&m.matrix, N - 1, N - 2));
        }

      francis_qrstep(&m.matrix, w, s, t);
      q = zero_subdiag_small_elements(&m.matrix);

      if (q == 0)
        {
          /* no small subdiagonal element found */
          continue;
        }

      if (q == (N - 1))
        {
          /*
           * the last subdiagonal element of the matrix is 0 -
           * m_{NN} is a real eigenvalue
           */
          GSL_SET_COMPLEX(&lambda1,
                          gsl_matrix_get(&m.matrix, q, q), 0.0);
          gsl_vector_complex_set(eval, evidx + nev++, lambda1);
          *nit = 0;

          --N;
          m = gsl_matrix_submatrix(&m.matrix, 0, 0, N, N);
        }
      else if (q == (N - 2))
        {
          gsl_matrix_view v;

          /*
           * The bottom right 2x2 block of m is an eigenvalue
           * system
           */

          v = gsl_matrix_submatrix(&m.matrix, q, q, 2, 2);
          get_2b2_eigenvalues(&v.matrix, &lambda1, &lambda2);
          gsl_vector_complex_set(eval, evidx + nev++, lambda1);
          gsl_vector_complex_set(eval, evidx + nev++, lambda2);
          *nit = 0;

          N -= 2;
          m = gsl_matrix_submatrix(&m.matrix, 0, 0, N, N);
        }
      else if (q == 1)
        {
          /* the first matrix element is an eigenvalue */
          GSL_SET_COMPLEX(&lambda1,
                          gsl_matrix_get(&m.matrix, 0, 0), 0.0);
          gsl_vector_complex_set(eval, evidx + nev++, lambda1);
          *nit = 0;

          --N;
          m = gsl_matrix_submatrix(&m.matrix, 1, 1, N, N);
        }
      else if (q == 2)
        {
          gsl_matrix_view v;

          /* the upper left 2x2 block is an eigenvalue system */

          v = gsl_matrix_submatrix(&m.matrix, 0, 0, 2, 2);
          get_2b2_eigenvalues(&v.matrix, &lambda1, &lambda2);
          gsl_vector_complex_set(eval, evidx + nev++, lambda1);
          gsl_vector_complex_set(eval, evidx + nev++, lambda2);
          *nit = 0;

          N -= 2;
          m = gsl_matrix_submatrix(&m.matrix, 2, 2, N, N);
        }
      else
        {
          /*
           * There is a zero element on the subdiagonal somewhere
           * in the middle of the matrix - we can now operate
           * separately on the two submatrices split by this
           * element. q is the row index of the zero element.
           */

          /* operate on lower right (N - q)x(N - q) block first */
          nev += schur_decomp(&m.matrix,
                              q,
                              N - 1,
                              eval,
                              evidx + nev,
                              nit,
                              w);

          /* operate on upper left qxq block */
          nev += schur_decomp(&m.matrix,
                              0,
                              q - 1,
                              eval,
                              evidx + nev,
                              nit,
                              w);
          N = 0;
        }
    }

  if (N == 1)
    {
      GSL_SET_COMPLEX(&lambda1, gsl_matrix_get(&m.matrix, 0, 0), 0.0);
      gsl_vector_complex_set(eval, evidx + nev++, lambda1);
      *nit = 0;
    }
  else if (N == 2)
    {
      get_2b2_eigenvalues(&m.matrix, &lambda1, &lambda2);
      gsl_vector_complex_set(eval, evidx + nev++, lambda1);
      gsl_vector_complex_set(eval, evidx + nev++, lambda2);
      *nit = 0;
    }

  return (nev);
}

/*
zero_subdiag_small_elements()
  Sets to zero all elements on the subdiaganal of a matrix A
which satisfy

|A_{i,i-1}| <= eps * (|A_{i,i}| + |A_{i-1,i-1}|)

Inputs: A - matrix (must be at least 3x3)

Return: row index of small subdiagonal element or 0 if not found
*/

inline static size_t
zero_subdiag_small_elements(gsl_matrix * A)
{
  const size_t N = A->size1;
  size_t i;
  double dpel = gsl_matrix_get(A, N - 2, N - 2);

  for (i = N - 1; i > 0; --i)
    {
      double sel = gsl_matrix_get(A, i, i - 1);
      double del = gsl_matrix_get(A, i, i);

      if ((sel == 0.0) ||
          (fabs(sel) < GSL_DBL_EPSILON * (fabs(del) + fabs(dpel))))
        {
          gsl_matrix_set(A, i, i - 1, 0.0);
          return (i);
        }

      dpel = del;
    }

  return (0);
}

/*
francis_qrstep()
  Perform a Francis QR step.

See Golub & Van Loan, "Matrix Computations" (3rd ed),
algorithm 7.5.1

Inputs: H - unreduced upper Hessenberg matrix
        w - workspace
        s - sum of current shifts (a_1 + a_2)
        t - product of shifts (a_1 * a_2)
*/

static inline int
francis_qrstep(gsl_matrix * H, gsl_eigen_unsymm_workspace * w,
               double s, double t)
{
  const size_t N = H->size1;
  double x, y, z;
  double scale;
  size_t i;
  gsl_matrix_view m;
  double tau_i;
  size_t q, r;

  x = gsl_matrix_get(H, 0, 0) * gsl_matrix_get(H, 0, 0) +
      gsl_matrix_get(H, 0, 1) * gsl_matrix_get(H, 1, 0) -
      s*gsl_matrix_get(H, 0, 0) + t;
  y = gsl_matrix_get(H, 1, 0) *
      (gsl_matrix_get(H, 0, 0) + gsl_matrix_get(H, 1, 1) - s);
  z = gsl_matrix_get(H, 1, 0) * gsl_matrix_get(H, 2, 1);

  scale = fabs(x) + fabs(y) + fabs(z);
  if (scale != 0.0)
    {
      /* scale to prevent overflow or underflow */
      x /= scale;
      y /= scale;
      z /= scale;
    }

  for (i = 0; i < N - 2; ++i)
    {
      double dat[3];
      gsl_vector_view v3 = gsl_vector_view_array(dat, 3);
      dat[0] = x; dat[1] = y; dat[2] = z;

      tau_i = gsl_linalg_householder_transform(&v3.vector);

      if (tau_i != 0.0)
        {
          /* q = max(1, i - 1) */
          q = (1 > ((int)i - 1)) ? 0 : (i - 1);

          /* apply left householder matrix (I - tau_i v v') to H */
          m = gsl_matrix_submatrix(H, i, q, 3, N - q);
          gsl_linalg_householder_hm(tau_i, &v3.vector, &m.matrix);

          /* r = min(i + 3, N - 1) */
          r = ((i + 3) < (N - 1)) ? (i + 3) : (N - 1);

          /* apply right householder matrix (I - tau_i v v') to H */
          m = gsl_matrix_submatrix(H, 0, i, r + 1, 3);
          gsl_linalg_householder_mh(tau_i, &v3.vector, &m.matrix);
        }

      x = gsl_matrix_get(H, i + 1, i);
      y = gsl_matrix_get(H, i + 2, i);
      if (i < (N - 3))
        {
          z = gsl_matrix_get(H, i + 3, i);
        }

      scale = fabs(x) + fabs(y) + fabs(z);
      if (scale != 0.0)
        {
          /* scale to prevent overflow or underflow */
          x /= scale;
          y /= scale;
          z /= scale;
        }
    }

  {
    double dat[2];
    gsl_vector_view v2 = gsl_vector_view_array(dat, 2);
    dat[0] = x; dat[1] = y;
    
    tau_i = gsl_linalg_householder_transform(&v2.vector);

    m = gsl_matrix_submatrix(H, N - 2, N - 3, 2, 3);
    gsl_linalg_householder_hm(tau_i, &v2.vector, &m.matrix);

    m = gsl_matrix_submatrix(H, 0, N - 2, N, 2);
    gsl_linalg_householder_mh(tau_i, &v2.vector, &m.matrix);
  }

  return GSL_SUCCESS;
}

/*
get_2b2_eigenvalues()
  Compute the eigenvalues of a 2x2 real matrix

Inputs: A  - 2x2 matrix
        e1 - where to store eigenvalue 1
        e2 - where to store eigenvalue 2
*/

static void
get_2b2_eigenvalues(gsl_matrix * A, gsl_complex * e1,
                    gsl_complex * e2)
{
  double discr;      /* discriminant of characteristic poly */
  double a, b, c, d; /* matrix values */

  a = gsl_matrix_get(A, 0, 0);
  b = gsl_matrix_get(A, 0, 1);
  c = gsl_matrix_get(A, 1, 0);
  d = gsl_matrix_get(A, 1, 1);

  discr = (a + d)*(a + d) - 4.0*(a*d - b*c);
  if (discr < 0.0)
    {
      GSL_SET_REAL(e1, 0.5*(a + d));
      GSL_SET_REAL(e2, 0.5*(a + d));

      GSL_SET_IMAG(e1, 0.5*sqrt(-discr));
      GSL_SET_IMAG(e2, -0.5*sqrt(-discr));
    }
  else
    {
      GSL_SET_REAL(e1, 0.5*(a + d + sqrt(discr)));
      GSL_SET_REAL(e2, 0.5*(a + d - sqrt(discr)));

      GSL_SET_IMAG(e1, 0.0);
      GSL_SET_IMAG(e2, 0.0);
    }
}
