/* multifit/multilinear.c
 * 
 * Copyright (C) 2000, 2007, 2010 Brian Gough
 * Copyright (C) 2013, 2015 Patrick Alken
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "linear_common.c"

int
gsl_multifit_linear (const gsl_matrix * X,
                     const gsl_vector * y,
                     gsl_vector * c,
                     gsl_matrix * cov,
                     double *chisq, gsl_multifit_linear_workspace * work)
{
  size_t rank;
  int status;
  double rnorm, snorm;

  status = gsl_multifit_linear_bsvd (X, work);
  if (status)
    return status;

  status = multifit_linear_solve (y, GSL_DBL_EPSILON, 0.0, &rank,
                                  c, cov, &rnorm, &snorm, chisq, work);

  return status;
}

/*
gsl_multifit_linear_svd()
  Perform SVD decomposition of the matrix X and store in work without
balancing
*/

int
gsl_multifit_linear_svd (const gsl_matrix * X,
                         gsl_multifit_linear_workspace * work)
{
  /* do not balance by default */
  int status = multifit_linear_svd(X, 0, work);

  return status;
} /* gsl_multifit_linear_svd() */

/*
gsl_multifit_linear_bsvd()
  Perform SVD decomposition of the matrix X and store in work with
balancing
*/

int
gsl_multifit_linear_bsvd (const gsl_matrix * X,
                          gsl_multifit_linear_workspace * work)
{
  int status = multifit_linear_svd(X, 1, work);

  return status;
} /* gsl_multifit_linear_bsvd() */

/* Estimation of values for given x */

int
gsl_multifit_linear_est (const gsl_vector * x,
                         const gsl_vector * c,
                         const gsl_matrix * cov, double *y, double *y_err)
{

  if (x->size != c->size)
    {
      GSL_ERROR ("number of parameters c does not match number of observations x",
         GSL_EBADLEN);
    }
  else if (cov->size1 != cov->size2)
    {
      GSL_ERROR ("covariance matrix is not square", GSL_ENOTSQR);
    }
  else if (c->size != cov->size1)
    {
      GSL_ERROR ("number of parameters c does not match size of covariance matrix cov",
         GSL_EBADLEN);
    }
  else
    {
      size_t i, j;
      double var = 0;
      
      gsl_blas_ddot(x, c, y);       /* y = x.c */

      /* var = x' cov x */

      for (i = 0; i < x->size; i++)
        {
          const double xi = gsl_vector_get (x, i);
          var += xi * xi * gsl_matrix_get (cov, i, i);

          for (j = 0; j < i; j++)
            {
              const double xj = gsl_vector_get (x, j);
              var += 2 * xi * xj * gsl_matrix_get (cov, i, j);
            }
        }

      *y_err = sqrt (var);

      return GSL_SUCCESS;
    }
}

/*
gsl_multifit_linear_residuals()
  Compute vector of residuals from fit

Inputs: X - design matrix
        y - rhs vector
        c - fit coefficients
        r - (output) where to store residuals
*/

int
gsl_multifit_linear_residuals (const gsl_matrix *X, const gsl_vector *y,
                               const gsl_vector *c, gsl_vector *r)
{
  if (X->size1 != y->size)
    {
      GSL_ERROR
        ("number of observations in y does not match rows of matrix X",
         GSL_EBADLEN);
    }
  else if (X->size2 != c->size)
    {
      GSL_ERROR ("number of parameters c does not match columns of matrix X",
                 GSL_EBADLEN);
    }
  else if (y->size != r->size)
    {
      GSL_ERROR ("number of observations in y does not match number of residuals",
                 GSL_EBADLEN);
    }
  else
    {
      /* r = y - X c */
      gsl_vector_memcpy(r, y);
      gsl_blas_dgemv(CblasNoTrans, -1.0, X, c, 1.0, r);

      return GSL_SUCCESS;
    }
} /* gsl_multifit_linear_residuals() */
