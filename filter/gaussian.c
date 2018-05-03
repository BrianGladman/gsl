/* filter/gaussian.c
 *
 * Gaussian smoothing filters
 * 
 * Copyright (C) 2018 Patrick Alken
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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_filter.h>

static int gausswin(const double sigma, const size_t order, const size_t n, double * window);

/*
gsl_filter_gaussian_alloc()
  Allocate a workspace for Gaussian filtering.

Inputs: K - number of samples in window; if even, it is rounded up to
            the next odd, to have a symmetric window

Return: pointer to workspace
*/

gsl_filter_gaussian_workspace *
gsl_filter_gaussian_alloc(const size_t K)
{
  gsl_filter_gaussian_workspace *w;

  w = calloc(1, sizeof(gsl_filter_gaussian_workspace));
  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace", GSL_ENOMEM);
    }

  w->H = K / 2;
  w->K = 2 * w->H + 1;

  w->kernel = malloc(w->K * sizeof(double));
  if (w->kernel == 0)
    {
      gsl_filter_gaussian_free(w);
      GSL_ERROR_NULL ("failed to allocate space for kernel", GSL_ENOMEM);
      return NULL;
    }

  return w;
}

void
gsl_filter_gaussian_free(gsl_filter_gaussian_workspace * w)
{
  if (w->kernel)
    free(w->kernel);

  free(w);
}

/*
gsl_filter_gaussian()
  Apply a Gaussian filter to an input vector:

G_{sigma}(x) = exp [ -x^2 / (2 sigma^2) ]

Inputs: sigma - standard deviation of Gaussian
        order - derivative order of Gaussian
        x     - input vector, size n
        y     - (output) filtered vector, size n
        w     - workspace
*/

int
gsl_filter_gaussian(const double sigma, const size_t order, const gsl_vector * x, gsl_vector * y, gsl_filter_gaussian_workspace * w)
{
  if (x->size != y->size)
    {
      GSL_ERROR("input and output vectors must have same length", GSL_EBADLEN);
    }
  else if (sigma <= 0.0)
    {
      GSL_ERROR("sigma must be positive", GSL_EDOM);
    }
  else if (order > 2)
    {
      GSL_ERROR("order must be <= 2", GSL_EDOM);
    }
  else
    {
      const int n = (int) x->size; /* input vector length */
      const int H = (int) w->H;    /* kernel radius */
      int i;

      /* construct Gaussian kernel of length K */
      gausswin(sigma, order, w->K, w->kernel);

      for (i = 0; i < n; ++i)
        {
          int idx1 = GSL_MAX(i - H, 0);
          int idx2 = GSL_MIN(i + H, n - 1);
          double sum = 0.0;
          int j;

          for (j = idx1; j <= idx2; ++j)
            {
              double xj = gsl_vector_get(x, j);
              sum += xj * w->kernel[H + i - j];
            }

          gsl_vector_set(y, i, sum);
        }

      return GSL_SUCCESS;
    }
}

static int
gausswin(const double sigma, const size_t order, const size_t n, double * window)
{
  const double m = n - 1.0;
  const double mhalf = m / 2.0;
  const double variance = sigma * sigma;
  const double alpha = 0.5 / variance;
  double sum = 0.0;
  double norm = 0.0; /* normalization */
  size_t i;

  if (order == 0)
    {
      for (i = 0; i < n; ++i)
        {
          double xi =( (double)i - mhalf) / mhalf;
          double gi = exp(-alpha * xi * xi);

          window[i] = gi;
          sum += window[i];
        }

      norm = 1.0 / sum;
    }
  else if (order == 1)
    {
      for (i = 0; i < n; ++i)
        {
          double xi =( (double)i - mhalf) / mhalf;
          double gi = exp(-alpha * xi * xi);

          window[i] = -xi * gi / variance;
          sum += gi;
        }

      norm = 1.0 / (sum * mhalf);
    }
  else if (order == 2)
    {
    }
  else
    {
      GSL_ERROR("order must be <= 2", GSL_EDOM);
    }

  /* normalize */
  for (i = 0; i < n; ++i)
    window[i] *= norm;

  return GSL_SUCCESS;
}
