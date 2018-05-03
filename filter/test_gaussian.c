/* filter/test_gaussian.c
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

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

static void
test_gaussian_fixed(void)
{
  const size_t K = 7;
  const size_t n = 20;
  const double sigma = 1.0;
  const double expected_y[] = { 1.352898375626185,  2.106288519406315,  3.000000000000000,  4.000000000000000,
                                5.000000000000000,  6.000000000000000,  7.000000000000000,  8.000000000000000,
                                9.000000000000000,  10.000000000000000, 11.000000000000000, 12.000000000000000,
                                13.000000000000000, 14.000000000000000, 15.000000000000000, 16.000000000000000,
                                17.000000000000000, 15.767941092467399, 13.714904499976400, 10.987123123526164 };
  gsl_vector_const_view ev = gsl_vector_const_view_array(expected_y, n);
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_filter_gaussian_workspace *w = gsl_filter_gaussian_alloc(K);
  size_t i;

  for (i = 0; i < n; ++i)
    gsl_vector_set(x, i, i + 1.0);

  gsl_filter_gaussian(sigma, 0, x, y, w);

  compare_vectors(1.0e-12, y, &ev.vector, "test_gaussian_fixed");

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_filter_gaussian_free(w);
}

static void
test_gaussian_deriv(const double sigma, const size_t n, const size_t K)
{
  const double f_low = 1.0;
  const double f_high = 50.0;
  const double alpha = 2.0 * M_PI / (n - 1.0);
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *dx = gsl_vector_alloc(n);
  gsl_vector *y1 = gsl_vector_alloc(n);
  gsl_vector *y2 = gsl_vector_alloc(n);
  gsl_filter_gaussian_workspace *w = gsl_filter_gaussian_alloc(K);
  size_t i;

  /* make input signal composed of two sine waves at different frequencies */
  for (i = 0; i < n; ++i)
    {
      double xi = sin(alpha * f_low * i) + sin(alpha * f_high * i);
      double dxi = alpha * f_low * cos(alpha * f_low * i) +
                   alpha * f_high * cos(alpha * f_high * i);

      gsl_vector_set(x, i, xi);
      gsl_vector_set(dx, i, dxi);
    }

  /* compute y1 = G * dx(t)/dt */
  gsl_filter_gaussian(sigma, 0, dx, y1, w);

  /* compute y2 = dG/dt * x(t) */
  gsl_filter_gaussian(sigma, 1, x, y2, w);

  for (i = 0; i < n; ++i)
    {
      double ti = (double) i / (n - 1.0);

      printf("%f %.12e %.12e %.12e %.12e\n",
             ti,
             gsl_vector_get(x, i),
             gsl_vector_get(dx, i),
             gsl_vector_get(y1, i),
             gsl_vector_get(y2, i));
    }

  gsl_vector_free(x);
  gsl_vector_free(dx);
  gsl_vector_free(y1);
  gsl_vector_free(y2);
  gsl_filter_gaussian_free(w);
}

static void
test_gaussian(void)
{
  test_gaussian_fixed();

  test_gaussian_deriv(0.25, 1000, 211);
}
