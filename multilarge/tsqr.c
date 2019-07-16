/* tsqr.c
 * 
 * Copyright (C) 2015, 2019 Patrick Alken
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
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multilarge.h>

typedef struct
{
  size_t p;             /* number of columns of LS matrix */
  int nblocks;          /* number of blocks processed */
  double rnorm;         /* || b - A x || residual norm */
  int svd;              /* SVD of R factor has been computed */

  gsl_matrix *T;        /* block reflector matrix, p-by-p */
  gsl_matrix *R;        /* R factor, p-by-p */
  gsl_vector *QTb;      /* [ Q^T b ; b_i ], size p-by-1 */
  gsl_vector *work;     /* workspace, size p */
  gsl_vector *work3;    /* workspace, size 3*p */

  gsl_multifit_linear_workspace * multifit_workspace_p;
} tsqr_state_t;

static void *tsqr_alloc(const size_t p);
static void tsqr_free(void *vstate);
static int tsqr_reset(void *vstate);
static int tsqr_accumulate(gsl_matrix * A, gsl_vector * b,
                           void * vstate);
static int tsqr_solve(const double lambda, gsl_vector * x,
                      double * rnorm, double * snorm,
                      void * vstate);
static int tsqr_rcond(double * rcond, void * vstate);
static int tsqr_lcurve(gsl_vector * reg_param, gsl_vector * rho,
                       gsl_vector * eta, void * vstate);
static double tsqr_householder_transform (double *v0, gsl_vector * v);
static int tsqr_QR_decomp (gsl_matrix * S, gsl_matrix * A, gsl_matrix * T);
static int tsqr_svd(tsqr_state_t * state);

/*
tsqr_alloc()
  Allocate workspace for solving large linear least squares
problems using the TSQR approach

Inputs: p    - number of columns of LS matrix

Return: pointer to workspace
*/

static void *
tsqr_alloc(const size_t p)
{
  tsqr_state_t *state;

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

  state->p = p;
  state->nblocks = 0;
  state->rnorm = 0.0;

  state->R = gsl_matrix_alloc(p, p);
  if (state->R == NULL)
    {
      tsqr_free(state);
      GSL_ERROR_NULL("failed to allocate R matrix", GSL_ENOMEM);
    }

  state->QTb = gsl_vector_alloc(p);
  if (state->QTb == NULL)
    {
      tsqr_free(state);
      GSL_ERROR_NULL("failed to allocate QTb vector", GSL_ENOMEM);
    }

  state->T = gsl_matrix_alloc(p, p);
  if (state->T == NULL)
    {
      tsqr_free(state);
      GSL_ERROR_NULL("failed to allocate T matrix", GSL_ENOMEM);
    }

  state->work = gsl_vector_alloc(p);
  if (state->work == NULL)
    {
      tsqr_free(state);
      GSL_ERROR_NULL("failed to allocate workspace vector", GSL_ENOMEM);
    }

  state->work3 = gsl_vector_alloc(3 * p);
  if (state->work3 == NULL)
    {
      tsqr_free(state);
      GSL_ERROR_NULL("failed to allocate work3 vector", GSL_ENOMEM);
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

  if (state->T)
    gsl_matrix_free(state->T);

  if (state->work)
    gsl_vector_free(state->work);

  if (state->work3)
    gsl_vector_free(state->work3);

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
  state->nblocks = 0;
  state->rnorm = 0.0;
  state->svd = 0;

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

3) A and b are destroyed
*/

static int
tsqr_accumulate(gsl_matrix * A, gsl_vector * b, void * vstate)
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
  else if (state->nblocks == 0 && n < p)
    {
      GSL_ERROR ("n must be >= p", GSL_EBADLEN);
    }
  else if (state->nblocks == 0)
    {
      int status;
      gsl_matrix_view R = gsl_matrix_submatrix(A, 0, 0, p, p);
      gsl_vector_view QTb = gsl_vector_subvector(state->QTb, 0, p);
      gsl_vector_view b1 = gsl_vector_subvector(b, 0, p);

      /* this is the first matrix block A_1, compute its (dense) QR decomposition */

      /* compute QR decomposition of A */
      status = gsl_linalg_QR_decomp_r(A, state->T);
      if (status)
        return status;

      /* store upper triangular R factor in state->R */
      gsl_matrix_tricpy('U', 1, state->R, &R.matrix);

      /* compute Q^T b and keep the first p elements */
      gsl_linalg_QR_QTvec_r(A, state->T, b, state->work);
      gsl_vector_memcpy(&QTb.vector, &b1.vector);

      if (n > p)
        {
          gsl_vector_view b2 = gsl_vector_subvector(b, p, n - p);
          state->rnorm = gsl_blas_dnrm2(&b2.vector);
        }
      else
        state->rnorm = 0.0;

      state->nblocks = 1;

      return GSL_SUCCESS;
    }
  else
    {
      int status;

      /* compute QR decomposition of [ R_{i-1} ; A_i ], accounting for
       * sparse structure */
      status = tsqr_QR_decomp(state->R, A, state->T);
      if (status)
        return status;

      /*
       * Compute:
       *
       * Q^T [ QTb_{i-1} ] = [ QTb_{i-1} - w ]
       *     [    b_i    ]   [  b_i - V~ w   ]
       *
       * where:
       * 
       * w = T^T (QTb_{i-1} + V~^T b_i)
       *
       *       p
       * V = [ I  ] p
       *     [ V~ ] n
       */
      gsl_vector_memcpy(state->work, state->QTb);
      gsl_blas_dgemv(CblasTrans, 1.0, A, b, 1.0, state->work);                     /* w := w + V~^T b */
      gsl_blas_dtrmv(CblasUpper, CblasTrans, CblasNonUnit, state->T, state->work); /* w := T^T w */
      gsl_vector_sub(state->QTb, state->work);                                     /* QTb := QTb - w */

      /* update residual norm */
      gsl_blas_dgemv(CblasNoTrans, -1.0, A, state->work, 1.0, b);                  /* b := b - V~ w */
      state->rnorm = gsl_hypot(state->rnorm, gsl_blas_dnrm2(b));

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

  if (x->size != state->p)
    {
      GSL_ERROR ("solution vector does not match workspace", GSL_EBADLEN);
    }
  else if (lambda < 0.0)
    {
      GSL_ERROR ("regularization parameter should be non-negative", GSL_EINVAL);
    }
  else
    {
      if (lambda == 0.0)
        {
          /* solve: R x = Q^T b */
          gsl_vector_memcpy(x, state->QTb);
          gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, state->R, x);
          *rnorm = state->rnorm;
          *snorm = gsl_blas_dnrm2(x);
        }
      else
        {
          int status;

          /* compute SVD of R if not already computed */
          if (state->svd == 0)
            {
              status = tsqr_svd(state);
              if (status)
                return status;
            }

          status = gsl_multifit_linear_solve(lambda, state->R, state->QTb, x, rnorm, snorm,
                                             state->multifit_workspace_p);
          if (status)
            return status;

          *rnorm = gsl_hypot(*rnorm, state->rnorm);
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
  int status;
  size_t i;

  /* compute SVD of R if not already computed */
  if (state->svd == 0)
    {
      status = tsqr_svd(state);
      if (status)
        return status;
    }

  status = gsl_multifit_linear_lcurve(state->QTb, reg_param, rho, eta,
                                      state->multifit_workspace_p);

  /* now add contribution to rnorm from Q2 factor */
  for (i = 0; i < rho->size; ++i)
    {
      double *rhoi = gsl_vector_ptr(rho, i);
      *rhoi = gsl_hypot(*rhoi, state->rnorm);
    }

  return status;
}

static int
tsqr_rcond(double * rcond, void * vstate)
{
  tsqr_state_t *state = (tsqr_state_t *) vstate;
  return gsl_linalg_tri_rcond(CblasUpper, state->R, rcond, state->work3);
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

x = [ x0 ] <- 1 nonzero value for the diagonal element of R
    [ 0  ] <- N - j - 1 zeros, where j is column of matrix in [0,N-1]
    [ x  ] <- M-N nonzero values for the dense part A

Inputs: v0 - pointer to diagonal element of R
             on input, v0 = x0;
        v  - on input, x vector
             on output, householder vector v
*/

static double
tsqr_householder_transform (double *v0, gsl_vector * v)
{
  /* replace v[0:M-1] with a householder vector (v[0:M-1]) and
     coefficient tau that annihilate v[1:M-1] */

  double alpha, beta, tau ;
  
  /* compute xnorm = || [ 0 ; v ] ||, ignoring zero part of vector */
  double xnorm = gsl_blas_dnrm2(v);

  if (xnorm == 0) 
    {
      return 0.0; /* tau = 0 */
    }

  alpha = *v0;
  beta = - GSL_SIGN(alpha) * hypot(alpha, xnorm) ;
  tau = (beta - alpha) / beta ;
  
  {
    double s = (alpha - beta);
    
    if (fabs(s) > GSL_DBL_MIN) 
      {
        gsl_blas_dscal (1.0 / s, v);
        *v0 = beta;
      }
    else
      {
        gsl_blas_dscal (GSL_DBL_EPSILON / s, v);
        gsl_blas_dscal (1.0 / GSL_DBL_EPSILON, v);
        *v0 = beta;
      }
  }
  
  return tau;
}

/*
tsqr_QR_decomp()
  Compute the QR decomposition of the matrix

  [ S ] = Q [ R ]
  [ A ]     [ 0 ]

where S is N-by-N upper triangular and A is M-by-N dense.

Inputs: S   - on input, upper triangular N-by-N matrix
              on output, R factor in upper triangle
        A   - dense M-by-N matrix
        T   - (output) block reflector matrix, N-by-N

Notes:
1) Based on the Elmroth/Gustavson algorithm, taking into account the
sparse structure of the S matrix
*/

static int
tsqr_QR_decomp (gsl_matrix * S, gsl_matrix * A, gsl_matrix * T)
{
  const size_t M = A->size1;
  const size_t N = S->size1;

  if (N != S->size2)
    {
      GSL_ERROR ("S matrix must be square", GSL_ENOTSQR);
    }
  else if (N != A->size2)
    {
      GSL_ERROR ("S and A have different number of columns", GSL_EBADLEN);
    }
  else if (T->size1 != N || T->size2 != N)
    {
      GSL_ERROR ("T matrix has wrong dimensions", GSL_EBADLEN);
    }
  else if (N == 1)
    {
      /* base case, compute Householder transform for single column matrix */
      double * T00 = gsl_matrix_ptr(T, 0, 0);
      double * S00 = gsl_matrix_ptr(S, 0, 0);
      gsl_vector_view v = gsl_matrix_column(A, 0);
      *T00 = tsqr_householder_transform(S00, &v.vector);
      return GSL_SUCCESS;
    }
  else
    {
      /*
       * partition matrices:
       *
       *       N1  N2              N1  N2
       * N1 [ S11 S12 ] and  N1 [ T11 T12 ]
       * N2 [  0  S22 ]      N2 [  0  T22 ]
       * M  [  A1  A2 ]
       */
      int status;
      const size_t N1 = N / 2;
      const size_t N2 = N - N1;

      gsl_matrix_view S11 = gsl_matrix_submatrix(S, 0, 0, N1, N1);
      gsl_matrix_view S12 = gsl_matrix_submatrix(S, 0, N1, N1, N2);
      gsl_matrix_view S22 = gsl_matrix_submatrix(S, N1, N1, N2, N2);

      gsl_matrix_view A1 = gsl_matrix_submatrix(A, 0, 0, M, N1);
      gsl_matrix_view A2 = gsl_matrix_submatrix(A, 0, N1, M, N2);

      gsl_matrix_view T11 = gsl_matrix_submatrix(T, 0, 0, N1, N1);
      gsl_matrix_view T12 = gsl_matrix_submatrix(T, 0, N1, N1, N2);
      gsl_matrix_view T22 = gsl_matrix_submatrix(T, N1, N1, N2, N2);

      /*
       * Eq. 2: recursively factor
       *
       *       N1           N1
       * N1 [ S11 ] = Q1 [ R11 ] N1
       * N2 [  0  ]      [  0  ] N2
       * M  [  A1 ]      [  0  ] M
       */
      status = tsqr_QR_decomp(&S11.matrix, &A1.matrix, &T11.matrix);
      if (status)
        return status;

      /*
       * Eq. 3:
       *
       *    N2                 N2           N2
       * N1 [ R12  ] = Q1^T [ S12 ] = [   S12 - W  ] N1
       * N2 [ S22~ ]        [ S22 ]   [     S22    ] N2
       * M  [  A2~ ]        [  A2 ]   [ A2 - V1~ W ] M
       *
       * where W = T11^T ( S12 + V1~^T A2 ), using T12 as temporary storage, and
       *
       *        N1
       * V1 = [  I  ] N1
       *      [  0  ] N2
       *      [ V1~ ] M
       */
      gsl_matrix_memcpy(&T12.matrix, &S12.matrix);                                                    /* W := S12 */
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &A1.matrix, &A2.matrix, 1.0, &T12.matrix);        /* W := S12 + V1~^T A2 */
      gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, 1.0, &T11.matrix, &T12.matrix); /* W := T11^T W */
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, &A1.matrix, &T12.matrix, 1.0, &A2.matrix);     /* A2 := A2 - V1~ W */
      gsl_matrix_sub(&S12.matrix, &T12.matrix);                                                       /* R12 := S12 - W */

      /*
       * Eq. 4: recursively factor
       *
       * [ S22~ ] = Q2~ [ R22 ]
       * [  A2~ ]       [  0  ]
       */
      status = tsqr_QR_decomp(&S22.matrix, &A2.matrix, &T22.matrix);
      if (status)
        return status;

      /*
       * Eq. 13: update T12 := -T11 * V1^T * V2 * T22
       *
       * where:
       *
       *        N1                N2
       * V1 = [  I  ] N1   V2 = [  0  ] N1
       *      [  0  ] N2        [  I  ] N2
       *      [ V1~ ] M         [ V2~ ] M
       *
       * Note: V1^T V2 = V1~^T V2~
       */

      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &A1.matrix, &A2.matrix, 0.0, &T12.matrix);           /* T12 := V1~^T * V2~ */
      gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, -1.0, &T11.matrix, &T12.matrix); /* T12 := -T11 * T12 */
      gsl_blas_dtrmm(CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, &T22.matrix, &T12.matrix); /* T12 := T12 * T22 */

      return GSL_SUCCESS;
    }
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
  int status;

  status = gsl_multifit_linear_svd(state->R, state->multifit_workspace_p);
  if (status)
    {
      GSL_ERROR("error computing SVD of R", status);
    }

  state->svd = 1;

  return GSL_SUCCESS;
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

const gsl_multilarge_linear_type * gsl_multilarge_linear_tsqr = &tsqr_type;
