/* movstat/test_qqr.c
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
#include <gsl/gsl_movstat.h>

/* compute moving QQR by explicitely constructing window and computing QQR */
int
slow_movqqr(const gsl_movstat_end_t etype, const double q, const gsl_vector * x,
            gsl_vector * y, const int H, const int J)
{
  const int n = (int) x->size;
  const int K = H + J + 1;
  double *window = malloc(K * sizeof(double));
  int i;

  for (i = 0; i < n; ++i)
    {
      int wsize = test_window(etype, i, H, J, x, window);
      double quant1, quant2;

      gsl_sort(window, 1, wsize);
      quant1 = gsl_stats_quantile_from_sorted_data(window, 1, wsize, q);
      quant2 = gsl_stats_quantile_from_sorted_data(window, 1, wsize, 1.0 - q);

      gsl_vector_set(y, i, quant2 - quant1);
    }

  free(window);

  return GSL_SUCCESS;
}

static void
test_qqr_proc(const double tol, const double q, const size_t n, const size_t H, const size_t J,
              const gsl_movstat_end_t etype, gsl_rng * rng_p)
{
  gsl_movstat_workspace * w = gsl_movstat_alloc2(H, J);
  gsl_vector * x = gsl_vector_alloc(n);
  gsl_vector * y = gsl_vector_alloc(n);
  gsl_vector * z = gsl_vector_alloc(n);
  char buf[2048];

  random_vector(x, rng_p);

  /* y = QQR(x) with slow brute force method */
  slow_movqqr(etype, q, x, y, H, J);

  /* y = QQR(x) with fast method */
  gsl_movstat_qqr(etype, x, q, z, w);

  /* test y = z */
  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u QQR random", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  /* z = QQR(x) in-place */
  gsl_vector_memcpy(z, x);
  gsl_movstat_qqr(etype, z, q, z, w);

  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u QQR random in-place", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  gsl_movstat_free(w);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(z);
}

static void
test_qqr(gsl_rng * rng_p)
{
  const double eps = 1.0e-10;

  test_qqr_proc(eps, 0.1, 100, 0, 0, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_qqr_proc(eps, 0.1, 1000, 3, 3, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_qqr_proc(eps, 0.2, 1000, 0, 5, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_qqr_proc(eps, 0.3, 1000, 5, 0, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_qqr_proc(eps, 0.25, 2000, 10, 5, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_qqr_proc(eps, 0.4, 2000, 5, 10, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_qqr_proc(eps, 0.25, 20, 50, 50, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_qqr_proc(eps, 0.25, 20, 10, 50, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_qqr_proc(eps, 0.25, 20, 50, 10, GSL_MOVSTAT_END_PADZERO, rng_p);

  test_qqr_proc(eps, 0.1, 100, 0, 0, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_qqr_proc(eps, 0.1, 1000, 3, 3, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_qqr_proc(eps, 0.2, 1000, 0, 5, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_qqr_proc(eps, 0.3, 1000, 5, 0, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_qqr_proc(eps, 0.25, 2000, 10, 5, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_qqr_proc(eps, 0.4, 2000, 5, 10, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_qqr_proc(eps, 0.25, 20, 50, 50, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_qqr_proc(eps, 0.25, 20, 10, 50, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_qqr_proc(eps, 0.25, 20, 50, 10, GSL_MOVSTAT_END_PADVALUE, rng_p);

  test_qqr_proc(eps, 0.1, 100, 0, 0, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_qqr_proc(eps, 0.1, 1000, 3, 3, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_qqr_proc(eps, 0.2, 1000, 0, 5, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_qqr_proc(eps, 0.3, 1000, 5, 0, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_qqr_proc(eps, 0.25, 2000, 10, 5, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_qqr_proc(eps, 0.4, 2000, 5, 10, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_qqr_proc(eps, 0.25, 20, 50, 50, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_qqr_proc(eps, 0.25, 20, 10, 50, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_qqr_proc(eps, 0.25, 20, 50, 10, GSL_MOVSTAT_END_TRUNCATE, rng_p);
}
