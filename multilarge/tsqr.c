/* tsqr.c
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

/*
 * This module implements the sequential TSQR algorithm
 * described in
 *
 * [1] Demmel, J., Grigori, L., Hoemmen, M. F., and Langou, J.
 *     "Communication-optimal parallel and sequential QR and LU factorizations",
 *     UCB Technical Report No. UCB/EECS-2008-89, 2008.
 *
 * The algorithm operates on a tall least squares system:
 *
 * [ A_1 ] x = [ b_1 ]
 * [ A_2 ]     [ b_2 ]
 * [ ... ]     [ ... ]
 * [ A_k ]     [ b_k ]
 *
 * as follows:
 *
 * 1. Initialize
 *    a. [Q_1,R_1] = qr(A_1)
 *    b. z_1 = Q_1^T b_1
 * 2. Loop i = 2:k
 *    a. [Q_i,R_i] = qr( [ R_{i-1} ; A_i ] )
 *    b. z_i = Q_i^T [ z_{i-1} ; b_i ]
 * 3. Output:
 *    a. R = R_k
 *    b. Q^T b = z_k
 *
 * Step 2(a) is optimized to take advantage
 * of the sparse structure of the matrix
 */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multilarge.h>
#include <gsl/gsl_multifit.h>

typedef struct
{
  size_t nmax;          /* maximum rows to add at once */
  size_t p;             /* number of columns of LS matrix */
  int init;             /* QR system has been initialized */
  int svd;              /* SVD of R has been computed */
  double normb;         /* || b || for computing residual norm */

  gsl_vector *tau;      /* Householder scalars, p-by-1 */
  gsl_matrix *R;        /* [ R ; A_i ], size (nmax + p)-by-p */
  gsl_vector *QTb;      /* [ Q^T b ; b_i ], size (nmax + p)-by-1 */

  gsl_multifit_linear_workspace *multifit_workspace_p;
} tsqr_state_t;

static void *tsqr_alloc(const size_t nmax, const size_t p);
static void tsqr_free(void *vstate);
static int tsqr_reset(void *vstate);
static int tsqr_accumulate(const gsl_matrix * A,
                             const gsl_vector * b,
                             void * vstate);
static int tsqr_solve(const double lambda, gsl_vector * x,
                        double * rnorm, double * snorm,
                        void * vstate);
static int tsqr_rcond(double * rcond, void * vstate);
static int tsqr_lcurve(gsl_vector * reg_param, gsl_vector * rho,
                       gsl_vector * eta, void * vstate);
static int tsqr_svd(tsqr_state_t * state);
static int tsqr_zero_R(gsl_matrix *R);
static double tsqr_householder_transform (const size_t N, const size_t j,
                                          gsl_vector * v);
static int tsqr_householder_hv (const size_t N, const size_t colidx, const double tau,
                                const gsl_vector * v, gsl_vector * w);
static int tsqr_QR_decomp (gsl_matrix * A, gsl_vector * tau);

/*
tsqr_alloc()
  Allocate workspace for solving large linear least squares
problems using the TSQR approach

Inputs: nmax - maximum number of rows to accumulate at once
        p    - number of columns of LS matrix

Return: pointer to workspace
*/

static void *
tsqr_alloc(const size_t nmax, const size_t p)
{
  tsqr_state_t *state;

  if (nmax == 0)
    {
      GSL_ERROR_NULL("nmax must be a positive integer",
                     GSL_EINVAL);
    }

  if (p == 0)
    {
      GSL_ERROR_NULL("p must be a positive integer",
                     GSL_EINVAL);
    }

  state = calloc(1, sizeof(tsqr_state_t));
  if (!state)
    {
      GSL_ERROR_NULL("failed to allocate tsqr state", GSL_ENOMEM);
    }

  state->nmax = nmax;
  state->p = p;
  state->init = 0;
  state->svd = 0;
  state->normb = 0.0;

  state->R = gsl_matrix_alloc(nmax + p, p);
  if (state->R == NULL)
    {
      tsqr_free(state);
      GSL_ERROR_NULL("failed to allocate R matrix", GSL_ENOMEM);
    }

  state->QTb = gsl_vector_alloc(nmax + p);
  if (state->QTb == NULL)
    {
      tsqr_free(state);
      GSL_ERROR_NULL("failed to allocate QTb vector", GSL_ENOMEM);
    }

  state->tau = gsl_vector_alloc(p);
  if (state->tau == NULL)
    {
      tsqr_free(state);
      GSL_ERROR_NULL("failed to allocate tau vector", GSL_ENOMEM);
    }

  state->multifit_workspace_p = gsl_multifit_linear_alloc(p, p);
  if (state->multifit_workspace_p == NULL)
    {
      tsqr_free(state);
      GSL_ERROR_NULL("failed to allocate multifit workspace", GSL_ENOMEM);
    }

  return state;
}

static void
tsqr_free(void *vstate)
{
  tsqr_state_t *state = (tsqr_state_t *) vstate;

  if (state->R)
    gsl_matrix_free(state->R);

  if (state->QTb)
    gsl_vector_free(state->QTb);

  if (state->tau)
    gsl_vector_free(state->tau);

  if (state->multifit_workspace_p)
    gsl_multifit_linear_free(state->multifit_workspace_p);

  free(state);
}

static int
tsqr_reset(void *vstate)
{
  tsqr_state_t *state = (tsqr_state_t *) vstate;

  gsl_matrix_set_zero(state->R);
  gsl_vector_set_zero(state->QTb);
  state->init = 0;
  state->svd = 0;
  state->normb = 0.0;

  return GSL_SUCCESS;
}

/*
tsqr_accumulate()
  Add a new block of rows to the QR system

Inputs: A      - new block of rows, n-by-p
        b      - new rhs vector n-by-1
        vstate - workspace

Return: success/error

Notes:
1) On output, the upper triangular portion of state->R(1:p,1:p)
contains current R matrix

2) state->QTb(1:p) contains current Q^T b vector
*/

static int
tsqr_accumulate(const gsl_matrix * A, const gsl_vector * b,
                void * vstate)
{
  tsqr_state_t *state = (tsqr_state_t *) vstate;
  const size_t n = A->size1;
  const size_t p = A->size2;

  if (p != state->p)
    {
      GSL_ERROR("columns of A do not match workspace", GSL_EBADLEN);
    }
  else if (n != b->size)
    {
      GSL_ERROR("A and b have different numbers of rows", GSL_EBADLEN);
    }
  else if (state->init == 0)
    {
      int status;
      gsl_vector_view tau = gsl_vector_subvector(state->tau, 0, GSL_MIN(n, p));
      gsl_matrix_view R = gsl_matrix_submatrix(state->R, 0, 0, n, p);
      gsl_vector_view QTb = gsl_vector_subvector(state->QTb, 0, n);

      /* this is the first matrix block A_1, compute its (dense) QR decomposition */

      /* copy A into the upper portion of state->R, so that R = [ A ; 0 ] */
      gsl_matrix_memcpy(&R.matrix, A);

      /* compute QR decomposition of A */
      status = gsl_linalg_QR_decomp(&R.matrix, &tau.vector);
      if (status)
        return status;

      /* compute Q^T b */
      gsl_vector_memcpy(&QTb.vector, b);
      gsl_linalg_QR_QTvec(&R.matrix, &tau.vector, &QTb.vector);

      /* compute ||b|| */
      state->normb = gsl_blas_dnrm2(b);

      state->init = 1;

      return GSL_SUCCESS;
    }
  else
    {
      int status;
      const size_t npp = n + p;
      gsl_vector_view tau = gsl_vector_subvector(state->tau, 0, p);
      gsl_matrix_view R = gsl_matrix_submatrix(state->R, 0, 0, npp, p);
      gsl_vector_view QTb = gsl_vector_subvector(state->QTb, 0, npp);
      gsl_matrix_view Ai = gsl_matrix_submatrix(state->R, p, 0, n, p);
      gsl_vector_view bi = gsl_vector_subvector(state->QTb, p, n);

      /* form state->R = [ R_{i-1} ; A_i ] */
      gsl_matrix_memcpy(&Ai.matrix, A);

      /* compute QR decomposition of [ R_{i-1] ; A_i ], accounting for
       * sparse structure */
      status = tsqr_QR_decomp(&R.matrix, &tau.vector);
      if (status)
        return status;

      {
        size_t i;

        /* compute Q^T [ QTb_{i - 1}; b_i ], accounting for the sparse
         * structure of the Householder reflectors */
        gsl_vector_memcpy(&bi.vector, b);
        for (i = 0; i < p; i++)
          {
            gsl_vector_const_view h = gsl_matrix_const_subcolumn (&R.matrix, i, i, npp - i);
            gsl_vector_view w = gsl_vector_subvector (&QTb.vector, i, npp - i);
            double ti = gsl_vector_get (&tau.vector, i);
            tsqr_householder_hv (p, i, ti, &(h.vector), &(w.vector));
          }
      }

      /* update ||b|| */
      state->normb = gsl_hypot(state->normb, gsl_blas_dnrm2(b));

      return GSL_SUCCESS;
    }
}

/*
tsqr_solve()
  Solve the least squares system:

chi^2 = || QTb - R x ||^2 + lambda^2 || x ||^2

using the SVD of R

Inputs: lambda - regularization parameter
        x      - (output) solution vector p-by-1
        rnorm  - (output) residual norm ||b - A x||
        snorm  - (output) solution norm ||x||
        vstate - workspace

Return: success/error
*/

static int
tsqr_solve(const double lambda, gsl_vector * x,
           double * rnorm, double * snorm,
           void * vstate)
{
  tsqr_state_t *state = (tsqr_state_t *) vstate;
  const size_t p = x->size;

  if (p != state->p)
    {
      GSL_ERROR("solution vector does not match workspace", GSL_EBADLEN);
    }
  else
    {
      int status;
      gsl_matrix_view R = gsl_matrix_submatrix(state->R, 0, 0, p, p);
      gsl_vector_view QTb = gsl_vector_subvector(state->QTb, 0, p);

      /* compute SVD of R if not already computed */
      if (state->svd == 0)
        {
          status = tsqr_svd(state);
          if (status)
            return status;
        }

      status = gsl_multifit_linear_solve(lambda, &R.matrix, &QTb.vector, x, rnorm, snorm,
                                         state->multifit_workspace_p);
      if (status)
        return status;

      /*
       * Since we're solving a reduced square system above, we need
       * to account for the full residual vector:
       *
       * rnorm = || [ Q1^T b - R x ; Q2^T b ] ||
       *
       * where Q1 is the thin Q factor of X, and Q2
       * are the remaining columns of Q. But:
       *
       * || Q2^T b ||^2 = ||b||^2 - ||Q1^T b||^2
       * 
       * so add this into the rnorm calculation
       */
      {
        double norm_Q1Tb = gsl_blas_dnrm2(&QTb.vector);
        double ratio = norm_Q1Tb / state->normb;
        double diff = 1.0 - ratio*ratio;

        if (diff > GSL_DBL_EPSILON)
          {
            double norm_Q2Tb = state->normb * sqrt(diff);
            *rnorm = gsl_hypot(*rnorm, norm_Q2Tb);
          }
      }

      return GSL_SUCCESS;
    }
}

/*
tsqr_lcurve()
  Compute L-curve of least squares system

Inputs: reg_param - (output) vector of regularization parameters
        rho       - (output) vector of residual norms
        eta       - (output) vector of solution norms
        vstate    - workspace

Return: success/error
*/

static int
tsqr_lcurve(gsl_vector * reg_param, gsl_vector * rho,
            gsl_vector * eta, void * vstate)
{
  tsqr_state_t *state = (tsqr_state_t *) vstate;
  const size_t p = state->p;
  gsl_vector_view QTb = gsl_vector_subvector(state->QTb, 0, p);
  int status;

  /* compute SVD of R if not already computed */
  if (state->svd == 0)
    {
      status = tsqr_svd(state);
      if (status)
        return status;
    }

  status = gsl_multifit_linear_lcurve(&QTb.vector, reg_param, rho, eta,
                                      state->multifit_workspace_p);

  return status;
}

static int
tsqr_rcond(double * rcond, void * vstate)
{
  tsqr_state_t *state = (tsqr_state_t *) vstate;

  /* compute SVD of R if not already computed */
  if (state->svd == 0)
    {
      int status = tsqr_svd(state);
      if (status)
        {
          GSL_ERROR("error computing SVD of R", GSL_EINVAL);
        }
    }

  *rcond = gsl_multifit_linear_rcond(state->multifit_workspace_p);

  return GSL_SUCCESS;
}

/*
tsqr_svd()
  Compute the SVD of the upper triangular
R factor. This allows us to compute the upper/lower
bounds on the regularization parameter and compute
the matrix reciprocal condition number.

Inputs: state - workspace

Return: success/error
*/

static int
tsqr_svd(tsqr_state_t * state)
{
  const size_t p = state->p;
  gsl_matrix_view R = gsl_matrix_submatrix(state->R, 0, 0, p, p);
  int status;

  /* XXX: zero R below the diagonal for SVD routine - this could be
   * optimized by developing SVD routine for triangular matrices */
  tsqr_zero_R(state->R);
  status = gsl_multifit_linear_svd(&R.matrix, state->multifit_workspace_p);

  state->svd = 1;

  return status;
}

/* zero everything below the diagonal */
static int
tsqr_zero_R(gsl_matrix *R)
{
  const size_t n = R->size1;
  const size_t p = R->size2;
  size_t j;

  for (j = 0; j < p; ++j)
    {
      gsl_vector_view v = gsl_matrix_subcolumn(R, j, j + 1, n - j - 1);
      gsl_vector_set_zero(&v.vector);
    }

  return GSL_SUCCESS;
}

/*
tsqr_householder_transform()
  This routine is an optimized version of
gsl_linalg_householder_transform(), designed for the QR
decomposition of M-by-N matrices of the form:

T = [ R ]
    [ A ]

where R is N-by-N upper triangular, and A is (M-N)-by-N dense.
This routine computes a householder transformation (tau,v) of a 
x so that P x = [ I - tau*v*v' ] x annihilates x(1:n-1). x will
be a subcolumn of the matrix T, and so its structure will be:

x = [ * ] <- 1 nonzero value
    [ 0 ] <- N - j - 1 zeros, where j is column of matrix in [0,N-1]
    [ * ] <- M-N nonzero values for the dense part A

Inputs: N      - number of columns in matrix T
        colidx - column number in [0, N-1]
        v      - on input, x vector
                 on output, householder vector v
*/

static double
tsqr_householder_transform (const size_t N, const size_t colidx, gsl_vector * v)
{
  /* replace v[0:M-1] with a householder vector (v[0:M-1]) and
     coefficient tau that annihilate v[1:M-1] */

  const size_t M = v->size ;

  if (M == 1)
    {
      return 0.0; /* tau = 0 */
    }
  else
    { 
      double alpha, beta, tau ;
      
      /* A portion of vector */
      gsl_vector_view x = gsl_vector_subvector (v, N - colidx, M - (N - colidx));

      /* compute xnorm = || v[1:M-1] ||, ignoring zero part of vector */
      double xnorm = gsl_blas_dnrm2(&x.vector);

      if (xnorm == 0) 
        {
          return 0.0; /* tau = 0 */
        }
      
      alpha = gsl_vector_get (v, 0) ;
      beta = - (alpha >= 0.0 ? +1.0 : -1.0) * hypot(alpha, xnorm) ;
      tau = (beta - alpha) / beta ;
      
      {
        double s = (alpha - beta);
        
        if (fabs(s) > GSL_DBL_MIN) 
          {
            gsl_blas_dscal (1.0 / s, &x.vector);
            gsl_vector_set (v, 0, beta) ;
          }
        else
          {
            gsl_blas_dscal (GSL_DBL_EPSILON / s, &x.vector);
            gsl_blas_dscal (1.0 / GSL_DBL_EPSILON, &x.vector);
            gsl_vector_set (v, 0, beta) ;
          }
      }
      
      return tau;
    }
}

static int
tsqr_householder_hv (const size_t N, const size_t colidx, const double tau,
                     const gsl_vector * v, gsl_vector * w)
{
  /* applies a householder transformation v to vector w */
  const size_t M = v->size;
 
  if (tau == 0)
    return GSL_SUCCESS ;

  {
    /* compute d = v'w */

    double w0 = gsl_vector_get(w,0);
    double d1, d;

    gsl_vector_const_view v1 = gsl_vector_const_subvector(v, N - colidx, M - (N - colidx));
    gsl_vector_view w1 = gsl_vector_subvector(w, N - colidx, M - (N - colidx));

    /* compute d1 = v(2:n)'w(2:n) */
    gsl_blas_ddot (&v1.vector, &w1.vector, &d1);

    /* compute d = v'w = w(1) + d1 since v(1) = 1 */
    d = w0 + d1;

    /* compute w = w - tau (v) (v'w) */

    gsl_vector_set (w, 0, w0 - tau * d);
    gsl_blas_daxpy (-tau * d, &v1.vector, &w1.vector);
  }
  
  return GSL_SUCCESS;
}

static int
tsqr_householder_hm (const size_t N, const size_t colidx, const double tau, const gsl_vector * v, gsl_matrix * A)
{
  /* applies a householder transformation v,tau to matrix m */

  if (tau == 0.0)
    {
      return GSL_SUCCESS;
    }
  else
    {
      gsl_vector_const_view v1 = gsl_vector_const_subvector (v, N - colidx, v->size - (N - colidx));
      gsl_matrix_view A1 = gsl_matrix_submatrix (A, N - colidx, 0, A->size1 - (N - colidx), A->size2);
      size_t j;

      for (j = 0; j < A->size2; j++)
        {
          double A0j = gsl_matrix_get (A, 0, j);
          double wj;
          gsl_vector_view A1j = gsl_matrix_column(&A1.matrix, j);

          gsl_blas_ddot (&A1j.vector, &v1.vector, &wj);
          wj += A0j;

          gsl_matrix_set (A, 0, j, A0j - tau *  wj);

          gsl_blas_daxpy (-tau * wj, &v1.vector, &A1j.vector);
        }

      return GSL_SUCCESS;
    }
}

/* Factorise a general M x N matrix A into
 *  
 *   A = Q R
 *
 * where Q is orthogonal (M x M) and R is upper triangular (M x N).
 *
 * Q is stored as a packed set of Householder transformations in the
 * strict lower triangular part of the input matrix.
 *
 * R is stored in the diagonal and upper triangle of the input matrix.
 *
 * The full matrix for Q can be obtained as the product
 *
 *       Q = Q_k .. Q_2 Q_1
 *
 * where k = MIN(M,N) and
 *
 *       Q_i = (I - tau_i * v_i * v_i')
 *
 * and where v_i is a Householder vector
 *
 *       v_i = [1, m(i+1,i), m(i+2,i), ... , m(M,i)]
 *
 * This storage scheme is the same as in LAPACK.  */

static int
tsqr_QR_decomp (gsl_matrix * A, gsl_vector * tau)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (tau->size != GSL_MIN (M, N))
    {
      GSL_ERROR ("size of tau must be MIN(M,N)", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      for (i = 0; i < GSL_MIN (M, N); i++)
        {
          /* Compute the Householder transformation to reduce the j-th
             column of the matrix to a multiple of the j-th unit vector,
             taking into account the sparse structure of A */

          gsl_vector_view c = gsl_matrix_subcolumn (A, i, i, M - i);
          double tau_i = tsqr_householder_transform(N, i, &c.vector);

          gsl_vector_set (tau, i, tau_i);

          /* Apply the transformation to the remaining columns and
             update the norms */

          if (i + 1 < N)
            {
              gsl_matrix_view m = gsl_matrix_submatrix (A, i, i + 1, M - i, N - (i + 1));
              tsqr_householder_hm (N, i, tau_i, &(c.vector), &(m.matrix));
            }
        }

      return GSL_SUCCESS;
    }
}

static const gsl_multilarge_linear_type tsqr_type =
{
  "tsqr",
  tsqr_alloc,
  tsqr_reset,
  tsqr_accumulate,
  tsqr_solve,
  tsqr_rcond,
  tsqr_lcurve,
  tsqr_free
};

const gsl_multilarge_linear_type * gsl_multilarge_linear_tsqr =
  &tsqr_type;
