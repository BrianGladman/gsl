/* randist/mvgauss.c
 * 
 * Copyright (C) 2016 Timothée Flutre, Patrick Alken
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
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>

static int multivar_mean (const double data[], size_t d, size_t tda, size_t n,
                          double mean[], size_t stride);
static int multivar_vcov (const double data[], size_t d, size_t tda, size_t n,
                          double vcov[], size_t tda2);

/* Generate a random vector from a multivariate Gaussian distribution using
 * the Cholesky decomposition of the variance-covariance matrix, following
 * "Computational Statistics" from Gentle (2009), section 7.4.
 *
 * mu      mean vector (dimension d)
 * L       matrix resulting from the Cholesky decomposition of
 *         variance-covariance matrix Sigma = L L^T (dimension d x d)
 * result  output vector (dimension d)
 */
int
gsl_ran_multivariate_gaussian (const gsl_rng * r,
                               const gsl_vector * mu,
                               const gsl_matrix * L,
                               gsl_vector * result)
{
  const size_t M = L->size1;
  const size_t N = L->size2;

  if (M != N)
    {
      GSL_ERROR("requires square matrix", GSL_ENOTSQR);
    }
  else if (mu->size != M)
    {
      GSL_ERROR("incompatible dimension of mean vector with variance-covariance matrix", GSL_EBADLEN);
    }
  else if (result->size != M)
    {
      GSL_ERROR("incompatible dimension of result vector", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      for (i = 0; i < M; ++i)
        gsl_vector_set(result, i, gsl_ran_ugaussian(r));

      gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, L, result);
      gsl_vector_add(result, mu);

      return GSL_SUCCESS;
    }
}

/* Compute the log of the probability density function at a given quantile
 * vector for a multivariate Gaussian distribution using the Cholesky
 * decomposition of the variance-covariance matrix.
 *
 * x       vector of quantiles (dimension d)
 * mu      mean vector (dimension d)
 * L       matrix resulting from the Cholesky decomposition of
 *         variance-covariance matrix Sigma = L L^T (dimension d x d)
 * result  output of the density (dimension 1)
 * work    vector used for intermediate computations (dimension d)
 */
int
gsl_ran_multivariate_gaussian_log_pdf (const gsl_vector * x,
                                       const gsl_vector * mu,
                                       const gsl_matrix * L,
                                       double * result,
                                       gsl_vector * work)
{
  const size_t M = L->size1;
  const size_t N = L->size2;

  if (M != N)
    {
      GSL_ERROR("requires square matrix", GSL_ENOTSQR);
    }
  else if (mu->size != M)
    {
      GSL_ERROR("incompatible dimension of mean vector with variance-covariance matrix", GSL_EBADLEN);
    }
  else if (x->size != M)
    {
      GSL_ERROR("incompatible dimension of quantile vector", GSL_EBADLEN);
    }
  else if (work->size != M)
    {
      GSL_ERROR("incompatible dimension of work vector", GSL_EBADLEN);
    }
  else
    {
      size_t i;
      double quadForm;        /* (x - mu)' Sigma^{-1} (x - mu) */
      double logSqrtDetSigma; /* log [ sqrt(|Sigma|) ] */

      /* compute: work = x - mu */
      for (i = 0; i < M; ++i)
        {
          double xi = gsl_vector_get(x, i);
          double mui = gsl_vector_get(mu, i);
          gsl_vector_set(work, i, xi - mui);
        }

      /* compute: work = L^{-1} * (x - mu) */
      gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, L, work);

      /* compute: quadForm = (x - mu)' Sigma^{-1} (x - mu) */
      gsl_blas_ddot(work, work, &quadForm);

      /* compute: log [ sqrt(|Sigma|) ] = sum_i log L_{ii} */
      logSqrtDetSigma = 0.0;
      for (i = 0; i < M; ++i)
        {
          double Lii = gsl_matrix_get(L, i, i);
          logSqrtDetSigma += log(Lii);
        }

      *result = -0.5*quadForm - logSqrtDetSigma - 0.5*M*log(2.0*M_PI);

      return GSL_SUCCESS;
    }
}

int
gsl_ran_multivariate_gaussian_pdf (const gsl_vector * x,
                                   const gsl_vector * mu,
                                   const gsl_matrix * L,
                                   double * result,
                                   gsl_vector * work)
{
  double logpdf;
  int status = gsl_ran_multivariate_gaussian_log_pdf(x, mu, L, &logpdf, work);

  if (status == GSL_SUCCESS)
    *result = exp(logpdf);

  return status;
}

/* Compute the maximum-likelihood estimate of the mean vector of samples
 * from a multivariate Gaussian distribution.
 */
int
gsl_ran_multivariate_gaussian_mean (const gsl_matrix * samples,
                                    gsl_vector * mean){
  return multivar_mean (samples->data, samples->size2, samples->tda,
                        samples->size1,
                        mean->data, mean->stride);
}

/* Compute the maximum-likelihood estimate of the variance-covariance matrix
 * of samples from a multivariate Gaussian distribution.
 */
int
gsl_ran_multivariate_gaussian_vcov (const gsl_matrix * samples,
                                    gsl_matrix * vcov){
  return multivar_vcov (samples->data, samples->size2, samples->tda,
                        samples->size1,
                        vcov->data, vcov->tda);
}

/* Example from R (GPL): http://www.r-project.org/
 * (samples <- matrix(c(4.348817, 2.995049, -3.793431, 4.711934, 1.190864, -1.357363), nrow=3, ncol=2))
 * colMeans(samples) # 1.183478 1.515145
 */
static int
multivar_mean (const double data[], size_t d, size_t tda, size_t n,
               double mean[], size_t stride)
{
  size_t j = 0;

  for (j = 0; j < d; ++j)
    {
      mean[j * stride] = gsl_stats_mean(&(data[j]), tda, n);
    }

  return GSL_SUCCESS;
}

/* Example from R (GPL): http://www.r-project.org/
 * (samples <- matrix(c(4.348817, 2.995049, -3.793431, 4.711934, 1.190864, -1.357363), nrow=3, ncol=2))
 * cov(samples) # 19.03539 11.91384 \n 11.91384  9.28796
 */
static int
multivar_vcov (const double data[], size_t d, size_t tda, size_t n,
               double vcov[], size_t tda2)
{
  size_t j1 = 0, j2 = 0;

  for (j1 = 0; j1 < d; ++j1)
    {
      vcov[j1 * tda2 + j1] = gsl_stats_variance(&(data[j1]), tda, n);
      for (j2 = j1 + 1; j2 < d; ++j2)
        {
          vcov[j1 * tda2 + j2] = gsl_stats_covariance(&(data[j1]), tda,
                                                      &(data[j2]), tda, n);
          vcov[j2 * tda2 + j1] = vcov[j1 * tda2 + j2];
        }
      }
  
  return GSL_SUCCESS;
}
