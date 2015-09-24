/* multifit/multiridge.c
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
 * Reference:
 *
 * [1] P. C. Hansen & D. P. O'Leary, "The use of the L-curve in
 * the regularization of discrete ill-posed problems",  SIAM J. Sci.
 * Comput. 14 (1993), pp. 1487-1503.
 */

#include <config.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "linear_common.c"

int
gsl_multifit_linear_ridge_solve (const double lambda,
                                 const gsl_vector * y,
                                 gsl_vector * c,
                                 gsl_matrix * cov,
                                 double *rnorm,
                                 double *snorm,
                                 gsl_multifit_linear_workspace * work)
{
  size_t rank;
  int status;
  double chisq;

  status = multifit_linear_solve(y, GSL_DBL_EPSILON, lambda, &rank, c,
                                 cov, rnorm, snorm, &chisq, work);

  return status;
} /* gsl_multifit_linear_ridge_solve() */

/*
gsl_multifit_linear_ridge_svd()
  Using regularization matrix L = \lambda G, with
G = diag(g_1,g_2,...,g_p), transform to Tikhonov standard form:

X~ = X G^{-1}
c~ = G c

and compute SVD of X~

Inputs: g    - Tikhonov matrix as a vector of diagonal elements
        X    - least squares matrix
        work - workspace

Return: success/error

Notes:
1) X~ is computed as well as its SVD which is stored in work
*/

int
gsl_multifit_linear_ridge_svd (const gsl_vector * g,
                               const gsl_matrix * X,
                               gsl_multifit_linear_workspace * work)
{
  const size_t n = work->n;
  const size_t p = work->p;

  if (p != g->size)
    {
      GSL_ERROR("g vector does not match workspace", GSL_EBADLEN);
    }
  else if (n != X->size1 || p != X->size2)
    {
      GSL_ERROR("X matrix does not match workspace", GSL_EBADLEN);
    }
  else
    {
      int status = GSL_SUCCESS;
      size_t j;

      /* construct X~ = X * G^{-1} matrix using work->A */
      for (j = 0; j < p; ++j)
        {
          gsl_vector_const_view Xj = gsl_matrix_const_column(X, j);
          gsl_vector_view Aj = gsl_matrix_column(work->A, j);
          double gj = gsl_vector_get(g, j);

          if (gj == 0.0)
            {
              GSL_ERROR("G matrix is singular", GSL_EDOM);
            }

          gsl_vector_memcpy(&Aj.vector, &Xj.vector);
          gsl_vector_scale(&Aj.vector, 1.0 / gj);
        }

      /* compute SVD of X~; do not balance for regularized problems */
      status = multifit_linear_svd (work->A, 0, work);

      return status;
    }
}

/*
gsl_multifit_linear_ridge_trans()
  Backtransform regularized solution vector using matrix
G = diag(g)
*/

int
gsl_multifit_linear_ridge_trans (const gsl_vector * g,
                                 gsl_vector * c,
                                 gsl_multifit_linear_workspace * work)
{
  const size_t p = work->p;

  if (p != g->size)
    {
      GSL_ERROR("g vector does not match workspace", GSL_EBADLEN);
    }
  else if (p != c->size)
    {
      GSL_ERROR("c vector does not match workspace", GSL_EBADLEN);
    }
  else
    {
      /* compute true solution vector c = G^{-1} c~ */
      gsl_vector_div(c, g);

      return GSL_SUCCESS;
    }
}

/*
gsl_multifit_linear_ridge_lcurve()
  Calculate L-curve using regularization parameters estimated
from singular values of least squares matrix

Inputs: y         - right hand side vector
        reg_param - (output) vector of regularization parameters
                    derived from singular values
        rho       - (output) vector of residual norms ||y - X c||
        eta       - (output) vector of solution norms ||lambda c||
        work      - workspace

Return: success/error
*/

int
gsl_multifit_linear_ridge_lcurve (const gsl_vector * y,
                                  gsl_vector * reg_param,
                                  gsl_vector * rho, gsl_vector * eta,
                                  gsl_multifit_linear_workspace * work)
{
  const size_t N = rho->size; /* number of points on L-curve */

  if (N < 3)
    {
      GSL_ERROR ("at least 3 points are needed for L-curve analysis",
                 GSL_EBADLEN);
    }
  else if (N != eta->size)
    {
      GSL_ERROR ("size of rho and eta vectors do not match",
                 GSL_EBADLEN);
    }
  else if (reg_param->size != eta->size)
    {
      GSL_ERROR ("size of reg_param and eta vectors do not match",
                 GSL_EBADLEN);
    }
  else
    {
      int status = GSL_SUCCESS;
      const size_t n = work->n;
      const size_t p = work->p;

      /* smallest regularization parameter */
      const double smin_ratio = 16.0 * GSL_DBL_EPSILON;

      double ratio, tmp;
      size_t i, j;

      gsl_matrix *A = work->A;
      gsl_vector *S = work->S;
      gsl_vector *xt = work->xt;
      gsl_vector_view workp = gsl_matrix_column(work->QSI, 0);
      gsl_vector *workp2 = work->D; /* D isn't used for regularized problems */

      const double s1 = gsl_vector_get(S, 0);
      const double sp = gsl_vector_get(S, p - 1);

      double dr; /* residual error from projection */
      double normy = gsl_blas_dnrm2(y);
      double normUTy;

      /* compute projection xt = U^T y */
      gsl_blas_dgemv (CblasTrans, 1.0, A, y, 0.0, xt);

      normUTy = gsl_blas_dnrm2(xt);
      dr = normy*normy - normUTy*normUTy;

      tmp = GSL_MAX(sp, s1*smin_ratio);
      gsl_vector_set(reg_param, N - 1, tmp);

      /* ratio so that reg_param(1) = s(1) */
      ratio = pow(s1 / tmp, 1.0 / (N - 1.0));

      /* calculate the regularization parameters */
      for (i = N - 1; i > 0 && i--; )
        {
          double rp1 = gsl_vector_get(reg_param, i + 1);
          gsl_vector_set(reg_param, i, ratio * rp1);
        }

      for (i = 0; i < N; ++i)
        {
          double lambda = gsl_vector_get(reg_param, i);
          double lambda_sq = lambda * lambda;

          for (j = 0; j < p; ++j)
            {
              double sj = gsl_vector_get(S, j);
              double xtj = gsl_vector_get(xt, j);
              double f = sj / (sj*sj + lambda_sq);

              gsl_vector_set(&workp.vector, j, f * xtj);
              gsl_vector_set(workp2, j, (1.0 - sj*f) * xtj);
            }

          gsl_vector_set(eta, i, gsl_blas_dnrm2(&workp.vector));
          gsl_vector_set(rho, i, gsl_blas_dnrm2(workp2));
        }

      if (n > p && dr > 0.0)
        {
          /* add correction to residual norm (see eqs 6-7 of [1]) */
          for (i = 0; i < N; ++i)
            {
              double rhoi = gsl_vector_get(rho, i);
              double *ptr = gsl_vector_ptr(rho, i);

              *ptr = sqrt(rhoi*rhoi + dr);
            }
        }

      /* restore D to identity matrix */
      gsl_vector_set_all(work->D, 1.0);

      return status;
    }
} /* gsl_multifit_linear_ridge_lcurve() */

/*
gsl_multifit_linear_ridge_lcorner()
  Determine point on L-curve of maximum curvature. For each
set of 3 points on the L-curve, the circle which passes through
the 3 points is computed. The radius of the circle is then used
as an estimate of the curvature at the middle point. The point
with maximum curvature is then selected.

Inputs: rho - vector of residual norms ||A x - b||
        eta - vector of solution norms ||L x||
        idx - (output) index i such that
              (log(rho(i)),log(eta(i)) is the point of
              maximum curvature

Return: success/error
*/

int
gsl_multifit_linear_ridge_lcorner(const gsl_vector *rho,
                                  const gsl_vector *eta,
                                  size_t *idx)
{
  const size_t n = rho->size;

  if (n < 3)
    {
      GSL_ERROR ("at least 3 points are needed for L-curve analysis",
                 GSL_EBADLEN);
    }
  else if (n != eta->size)
    {
      GSL_ERROR ("size of rho and eta vectors do not match",
                 GSL_EBADLEN);
    }
  else
    {
      int s = GSL_SUCCESS;
      size_t i;
      double x1, y1;      /* first point of triangle on L-curve */
      double x2, y2;      /* second point of triangle on L-curve */
      double rmin = -1.0; /* minimum radius of curvature */

      /* initial values */
      x1 = log(gsl_vector_get(rho, 0));
      y1 = log(gsl_vector_get(eta, 0));

      x2 = log(gsl_vector_get(rho, 1));
      y2 = log(gsl_vector_get(eta, 1));

      for (i = 1; i < n - 1; ++i)
        {
          /*
           * The points (x1,y1), (x2,y2), (x3,y3) are the previous,
           * current, and next point on the L-curve. We will find
           * the circle which fits these 3 points and take its radius
           * as an estimate of the curvature at this point.
           */
          double x3 = log(gsl_vector_get(rho, i + 1));
          double y3 = log(gsl_vector_get(eta, i + 1));

          double x21 = x2 - x1;
          double y21 = y2 - y1;
          double x31 = x3 - x1;
          double y31 = y3 - y1;
          double h21 = x21*x21 + y21*y21;
          double h31 = x31*x31 + y31*y31;
          double d = fabs(2.0 * (x21*y31 - x31*y21));
          double r = sqrt(h21*h31*((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))) / d;

          /* if d =~ 0 then there are nearly colinear points */
          if (gsl_finite(r))
            {
              /* check for smallest radius of curvature */
              if (r < rmin || rmin < 0.0)
                {
                  rmin = r;
                  *idx = i;
                }
            }

          /* update previous/current log values */
          x1 = x2;
          y1 = y2;
          x2 = x3;
          y2 = y3;
        }

      /* check if a minimum radius was found */
      if (rmin < 0.0)
        {
          /* possibly co-linear points */
          GSL_ERROR("failed to find minimum radius", GSL_EINVAL);
        }

      return s;
    }
} /* gsl_multifit_linear_ridge_lcorner() */
