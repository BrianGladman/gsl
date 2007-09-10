/* statistics/wcovar_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Jim Davies, Brian Gough
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

static double 
FUNCTION(compute,wcovariance) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n, const double wmean);

static double
FUNCTION(compute,factor) (const BASE w[], const size_t wstride, const size_t n);

static double
FUNCTION(compute,wcovariance) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n, const double wmean)
{
  /* takes a dataset and finds the weighted covariance */

  long double wcovariance = 0 ;
  long double W = 0;

  size_t i;

  /* find the sum of the squares */
  for (i = 0; i < n; i++)
    {
      BASE wi = w[i * wstride];

      if (wi > 0) {
        const long double delta = (data[i * stride] - wmean);
        W += wi ;
        wcovariance += (delta * delta - wcovariance) * (wi / W);
      }
    }

  return wcovariance ;
}

static double
FUNCTION(compute,factor) (const BASE w[], const size_t wstride, const size_t n)
{
  /* Find the factor ``N/(N-1)'' which multiplies the raw std dev */

  long double a = 0 ;
  long double b = 0;
  long double factor;

  size_t i;

  /* find the sum of the squares */
  for (i = 0; i < n; i++)
    {
      BASE wi = w[i * wstride];

      if (wi > 0)
        {
          a += wi ;
          b += wi * wi ;
        }
    }

  factor = (a*a) / ((a*a) - b);

  return factor ;
}

double 
FUNCTION(gsl_stats,wcovariance_m) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n, const double wmean)
{
  const double covariance = FUNCTION(compute,wcovariance) (w, wstride, data, stride, n, wmean);
  const double scale = FUNCTION(compute,factor)(w, wstride, n);
  
  return scale * covariance;
}

double 
FUNCTION(gsl_stats,wcovariance) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n)
{
  const double wmean = FUNCTION(gsl_stats,wmean) (w, wstride, data, stride, n);
  return FUNCTION(gsl_stats,wcovariance_m)(w, wstride, data, stride, n, wmean);
}
